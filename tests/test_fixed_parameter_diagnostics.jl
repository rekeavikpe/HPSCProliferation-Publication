using Test
include("../profiles/run_fixed_parameter_diagnostics.jl")

@testset "Pointwise profiles preserve nuisance optimization and isolate jobs" begin
    # Curved ridge: reoptimizing y leaves a known one-dimensional profile for x.
    loglik=z->-((z[1]-.5)^2/.02+(z[2]-(.2+.6z[1]))^2/.002)
    reference=[.5,.5]
    coordinates=collect(range(.1,.9;length=17))
    results=Vector{Any}(undef,length(coordinates))
    starts=[[.1],[.9],[.5]]
    Threads.@threads for j in eachindex(coordinates)
        results[j]=optimize_diagnostic_point(loglik,reference,0.,1,coordinates[j],deepcopy(starts);maxiters=1000)
    end
    for (x,r) in zip(coordinates,results)
        @test isapprox(r.relative_loglik,-(x-.5)^2/.02;atol=1e-7)
        @test isapprox(r.z[2],.2+.6x;atol=1e-6)
        @test r.stopping_met
    end
    @test starts==[[.1],[.9],[.5]]
    @test reference==[.5,.5]
    @test length(unique(point_id.(coordinates)))==length(coordinates)
end

@testset "Confidence limits distinguish crossings from search bounds" begin
    x=collect(range(0.,1.;length=201));y=-((x.-.5).^2)./.02
    ci=diagnostic_limits(x,y,.5)
    @test ci.lower_crossed && ci.upper_crossed
    @test isapprox(ci.lower,.5-sqrt(-DIAGNOSTIC_THRESHOLD*.02);atol=1e-5)
    @test isapprox(ci.upper,.5+sqrt(-DIAGNOSTIC_THRESHOLD*.02);atol=1e-5)
    leftopen=diagnostic_limits(x,-0.5 .* (x.-.5).^2,.5)
    @test !leftopen.lower_crossed && !leftopen.upper_crossed
    @test leftopen.lower==0. && leftopen.upper==1.
    boundmle=diagnostic_limits(x,-4 .* (1 .- x),1.)
    @test boundmle.lower_crossed && !boundmle.upper_crossed
    @test diagnostic_grid([0.,.1,.5,1-4e-11,1-1e-11,1.],1.)==[0.,.1,.5,1.]
    @test isempty(diagnostic_curve_issues(y,fill(true,length(y))))
    @test !isempty(diagnostic_curve_issues([0.,-1.,-.1,-2.],fill(true,4)))
    @test !isempty(diagnostic_curve_issues(y,[false;fill(true,length(y)-1)]))
end

if !isempty(ARGS)
    @testset "Real-data design and source protection" begin
        prep=diagnostic_setup(ARGS[1],ARGS[2])
        @test length(prep.names)==11
        @test isempty(prep.ctx.config["fixed"])
        @test :K_G in prep.names && :K_Lac in prep.names && :Y_GX in prep.names
        @test length(prep.ctx.config["calibration_cohort"]["series_ids"])==19
        @test length(prep.ctx.y)==116
        @test isapprox(prep.nll,prep.saved_nll;atol=1e-6)
        @test_throws ErrorException safe_output(ARGS[1],[ARGS[1],ARGS[2]])
        @test_throws ErrorException safe_output(joinpath(ARGS[2],"overwrite"),[ARGS[1],ARGS[2]])
        for name in DIAGNOSTIC_NAMES
            i=findfirst(==(Symbol(name)),prep.names)
            bank=load_seed_bank(ARGS[2],prep.names,i,prep.ctx.config["bounds"],prep.z)
            coordinate=name=="Y_GX" ? bank.coordinates[108] : first(bank.coordinates)
            starts=diagnostic_starts(bank,coordinate)
            @test all(s->length(s)==10,starts)
            @test starts!==bank.vectors
        end
    end
end

if !isempty(ARGS)
    @testset "Saved-point export, numerical dip detection and safe grid extension" begin
        prep=diagnostic_setup(ARGS[1],ARGS[2])
        mktempdir() do out
            for name in DIAGNOSTIC_NAMES
                source=joinpath(ARGS[2],name)
                rows=collect(CSV.File(joinpath(source,"profile.csv")))
                nuisance=collect(CSV.File(joinpath(source,"nuisance_parameters.csv")))
                for (j,r) in enumerate(rows)
                    dest=point_directory(out,name,r.normalized_parameter);mkpath(dest)
                    values=Dict(Symbol(v.parameter)=>Float64(v.value) for v in nuisance if v.point==j)
                    p=NamedTuple{Tuple(prep.names)}(Tuple(values[n] for n in prep.names))
                    save_parameters(joinpath(dest,"parameters.csv"),p)
                    CSV.write(joinpath(dest,"point.csv"),[(parameter=name,parameter_value=r.parameter_value,
                        normalized_parameter=r.normalized_parameter,profile_loglik=r.profile_loglik,
                        stopping_met=true,retcode="Success",starts=1,best_seed=-r.profile_loglik)])
                    write(joinpath(dest,"COMPLETE"),"Synthetic checkpoint fixture from saved curves\n")
                end
                @test length(completed_points(out,name))==length(rows)
                @test isempty(extension_jobs(prep,out,name))
                result=redirect_stdout(devnull) do
                    export_diagnostic(prep,out,name)
                end
                @test isfile(joinpath(out,name,"profile.csv"))
                @test isfile(joinpath(out,name,"nuisance_parameters.csv"))
                if name=="Y_GX"
                    @test result.audit.lower_crossed==false && result.audit.upper_crossed==true
                    @test result.audit.numerical_review_required
                    @test occursin("sharp dips",result.audit.status)
                    @test isnan(result.audit.lower)
                    # Removing only the far-right fixture points must request a
                    # further grid batch, using the last completed nuisance fit.
                    for r in completed_points(out,name)
                        r.normalized_parameter>.8 && rm(r.directory;recursive=true)
                    end
                    jobs=extension_jobs(prep,out,name)
                    @test !isempty(jobs)
                    @test all(j->j.coordinate>last(completed_points(out,name)).normalized_parameter,jobs)
                    @test all(j->length(only(j.extra))==10,jobs)
                else
                    @test result.audit.lower_crossed==true && result.audit.upper_crossed==false
                    @test !result.audit.numerical_review_required
                    @test isnan(result.audit.upper)
                end
            end
        end
    end
end
