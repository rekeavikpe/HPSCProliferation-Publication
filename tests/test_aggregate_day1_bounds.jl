using Test, CSV, TOML
include("../revision/AggregateDayOneBounds.jl")
using .AggregateDayOneBounds
using .AggregateDayOneBounds.AggregateDayOne
using .AggregateDayOne.MixedExposure, .MixedExposure.RevisionCommon
using .RevisionCommon.Workflow, .Workflow.Calibration, .Calibration.Proliferation
const B=AggregateDayOneBounds

length(ARGS)==3 || error("Supply REFERENCE_FIT SPLIT_FIT COMPLETED_D1_RESULTS")
fit,split,source=abspath.(ARGS)
before=B.snapshot(source)
B.verify_parent(source,fit,split)
prep,p0=source_contexts(fit,split)
train=aggregate_dayone(prep.train)
p,bank,paths=B.saved_seeds(source,train,3)

@testset "Read-only completed fits and isolated outputs" begin
    @test isfile(joinpath(source,"COMPLETE"))
    @test length(bank)==3
    @test isapprox(B.actual(train,p),545.3364247994402;atol=1e-6,rtol=0)
    for out in (source,joinpath(source,"checks"),dirname(source),fit,split)
        @test_throws ErrorException B.safe_output(out,[source,fit,split])
    end
    mktempdir() do dir
        @test B.safe_output(joinpath(dir,"new"),[source])==joinpath(realpath(dir),"new")
        symlink(source,joinpath(dir,"alias"))
        @test_throws ErrorException B.safe_output(joinpath(dir,"alias","new"),[source])
    end
end

@testset "Only approved upper limits change" begin
    baseline=deepcopy(train.config)
    for multiplier in B.MULTIPLIERS
        ctx=B.bounded_context(train,multiplier)
        @test length(free_names(ctx,6))==11
        @test ctx.observations===train.observations
        @test ctx.y===train.y && ctx.sigma===train.sigma
        @test ctx.conditions===train.conditions
        @test prediction(ctx,6,p)==prediction(train,6,p)
        for (name,bounds) in train.config["bounds"]
            @test ctx.config["bounds"][name][1]==bounds[1]
            @test ctx.config["bounds"][name][2]==bounds[2]*(Symbol(name) in B.TARGETS ? multiplier : 1)
        end
        @test train.config==baseline
    end
    for target in B.TARGETS
        ctx=B.bounded_context(train,100;fixed=Dict(string(target)=>getproperty(p,target)))
        @test length(free_names(ctx,6))==10
        @test !(target in free_names(ctx,6))
        @test only(keys(ctx.config["fixed"]))==string(target)
        @test all(n in free_names(ctx,6) for n in parameter_names(6) if n!=target)
        @test all(isapprox.(collect(values(decode(ctx,6,encode(ctx,6,p)))),collect(values(p));rtol=1e-13))
        @test ctx.y==train.y && ctx.sigma==train.sigma
    end
    @test_throws ErrorException B.bounded_context(train,2)
    @test_throws ErrorException B.bounded_context(train,100;fixed=Dict("mu_max"=>.1))
end

@testset "Every planned fixed-value check has feasible training-only seeds" begin
    for target in B.TARGETS, value in getproperty(B.GRID,target)
        ctx=B.bounded_context(train,100;fixed=Dict(string(target)=>value))
        seeds=B.unique_seeds(ctx,bank;limit=2,feasible=true)
        @test !isempty(seeds)
        @test all(isfinite(B.actual(ctx,q)) for q in seeds)
        @test all(isapprox(getproperty(q,target),value;rtol=1e-12) for q in seeds)
        @test :mu_max in free_names(ctx,6)
    end
end

@testset "Original training and selection separation and scoring" begin
    @test length(prep.split.training)==19 && length(prep.split.selection)==12
    @test isempty(intersect(prep.split.training,prep.split.selection))
    held=aggregate_dayone(prep.held)
    @test length(train.y)==111 && length(held.y)==111
    @test all(r->r.state!="aggregate"||r.day>1,held.observations)
    @test isapprox(AggregateDayOne.score_row(held,p,"saved","selection40").bic,1197.0109964438373;atol=1e-6,rtol=0)
    @test isapprox(AggregateDayOne.score_row(comparison_context(prep.held),p0,"D0","selection40").bic,1344.8598464335823;atol=1e-6,rtol=0)
    predictions=prediction(train,6,p)
    altered=merge(prep.held,(y=prep.held.y .* 2,))
    @test prediction(train,6,p)==predictions
    @test length(MixedExposure.evaluations(prep))==8
end

@testset "Checkpoint resume, seed protection and changed-input rejection" begin
    mktempdir() do dir
        ctx=B.bounded_context(train,1)
        result=B.optimize_one(ctx,p,joinpath(dir,"start");lm_iterations=1,polish_iterations=2)
        @test result.nll<=B.actual(ctx,p)+1e-7
        hashes=B.snapshot(dir)
        resumed=B.optimize_one(ctx,p,joinpath(dir,"start");lm_iterations=1,polish_iterations=2)
        @test resumed.nll==result.nll
        @test B.snapshot(dir)==hashes
        @test_throws ErrorException B.optimize_one(ctx,p,joinpath(dir,"start");lm_iterations=2,polish_iterations=2)
        @test B.snapshot(dir)==hashes
        fresh=joinpath(dir,"locked")
        B.lock_run(fresh,Dict("test"=>true),source,[joinpath(source,"manifest.toml")])
        @test_throws ErrorException B.lock_run(fresh,Dict("test"=>false),source,[joinpath(source,"manifest.toml")])
    end
    @test B.snapshot(source)==before
end
