using Test, CSV, Statistics
include("../revision/AggregateDayOne.jl")
using .AggregateDayOne
using .AggregateDayOne.MixedExposure
using .MixedExposure.RevisionCommon
using .RevisionCommon.Workflow, .Workflow.Calibration, .Calibration.Proliferation

length(ARGS)==2 || error("Supply REFERENCE_FIT SPLIT_FIT")
prep,p = MixedExposure.source_contexts(ARGS...)

@testset "Preserved cohorts, P50 inputs and D1 aggregate initialization" begin
    @test length(prep.split.training)==19
    @test length(prep.split.selection)==12
    @test isempty(intersect(prep.split.training,prep.split.selection))
    for (cohort,raw,settings) in MixedExposure.evaluations(prep)
        base = replace_exposure(raw;edr_percentile=50,settings)
        ctx = aggregate_dayone(base)
        ref = comparison_context(base)
        @test ctx.y==ref.y
        @test ctx.sigma==ref.sigma
        @test ctx.observations==ref.observations
        @test ctx.config==base.config
        @test length(free_names(ctx,6))==11
        @test all(r->r.state!="aggregate" || r.day>1,ctx.observations)
        @test count(r->r.state!="aggregate",ctx.observations)==count(r->r.state!="aggregate",base.observations)
        mask = [r.state!="aggregate" || r.day>1 for r in base.observations]
        @test prediction(ref,6,p) ≈ prediction(base,6,p)[mask] rtol=1e-11
        @test all(isfinite,prediction(ctx,6,p))
        for (c,q) in ctx.conditions
            source = base.conditions[c]
            measured = only(r.mean for r in base.observations if r.condition==c && r.day==1 && r.state=="aggregate")
            @test q.aggregate_day1==measured
            @test q.u0==source.u0
            @test q.u0[4]==15.
            @test q.times==source.times
            @test q.feeds==source.feeds
            @test q.glucose==source.glucose
            for size in (15.,200.,300.,1000.)
                @test collect(q.lookup(size)) ≈ collect(source.lookup(size)) rtol=1e-12
                @test collect(q.lookup(size)) ≈ collect(raw.conditions[c].lookup(size)) rtol=1e-12
            end
            @test collect(q.lookup(15.))==collect(q.lookup(200.))
        end
        @test length(base.y)-length(ctx.y)==length(base.conditions)
    end
    @test length(aggregate_dayone(prep.train).y)==length(prep.train.y)-5
    @test length(aggregate_dayone(prep.held).y)==length(prep.held.y)-5
    # Changing selection D1 observations cannot change training initial conditions or residuals.
    training_before = prediction(aggregate_dayone(prep.train),6,p)
    obs = [r.state=="aggregate" && r.day==1 ? merge(r,(mean=r.mean+50.,)) : r for r in prep.held.observations]
    modified = merge(prep.held,(observations=obs,y=[r.mean for r in obs]))
    held = aggregate_dayone(modified)
    @test held.y==aggregate_dayone(prep.held).y
    @test all(held.conditions[c].aggregate_day1==aggregate_dayone(prep.held).conditions[c].aggregate_day1+50 for c in keys(held.conditions))
    @test prediction(aggregate_dayone(prep.train),6,p)==training_before
    missing = merge(prep.train,(observations=filter(r->!(r.day==1 && r.state=="aggregate"),prep.train.observations),))
    @test_throws ErrorException aggregate_dayone(missing)
end

@testset "Three-state first-day equations preserve the original cell and metabolite laws" begin
    for u in ([8e4,20.,0.],[1.5e5,10.,5.],[0.,0.,-1.]), t in (0.,12.,23.9), inhibition in (true,false)
        full=zeros(4);small=zeros(3)
        model_rhs!(full,vcat(u,15.),p,t,6,.01,1.;lactate_inhibition=inhibition)
        AggregateDayOne.first_day_rhs!(small,u,p,t,.01;lactate_inhibition=inhibition)
        @test small ≈ full[1:3] rtol=1e-14
    end
end

@testset "Continuity, absolute death clock and feeding" begin
    source = prep.train.conditions[1]
    times = [0.,12.,24.,36.,48.,72.,96.,120.]
    feeds = filter(f->f.time_h<=120,source.feeds)
    q = merge(source,(times,feeds,aggregate_day1=180.,aggregate_day1_n=3))
    sim = simulate_aggregate_dayone(6,p,q)
    # Independent four-state oracle: the old RHS with a zero fourth derivative before D1.
    function frozen_rhs!(du,u,unused,t)
        tau,edr=q.lookup(15.)
        model_rhs!(du,u,p,t,6,tau,edr)
        du[4]=0.
    end
    early=Proliferation.solve(Proliferation.ODEProblem(frozen_rhs!,q.u0,(0.,24.)),Proliferation.Tsit5();reltol=1e-12,abstol=1e-13)
    @test sim.predictions[1:3,1:3] ≈ hcat([early(t)[1:3] for t in times[1:3]]...) rtol=1e-9
    @test all(==(15.),sim.predictions[4,1:2])
    @test sim.initialization.before[4]==15.
    @test sim.predictions[4,3]==180.
    @test sim.initialization.after[1:3]==sim.initialization.before[1:3]
    @test [e.time_h for e in sim.events]==[f.time_h for f in feeds]
    # Independent post-D1 oracle uses absolute time, including the death term and feeds.
    function late_rhs!(du,u,unused,t)
        tau,edr=q.lookup(u[4])
        model_rhs!(du,u,p,t,6,tau,edr)
    end
    carry=vcat(early(24.)[1:3],q.aggregate_day1)
    start=24.
    for stop in sort(unique(vcat([f.time_h for f in feeds],last(times))))
        sol=Proliferation.solve(Proliferation.ODEProblem(late_rhs!,carry,(start,stop)),Proliferation.Tsit5();reltol=1e-12,abstol=1e-13)
        for (j,t) in enumerate(times)
            start<t<=stop || continue
            @test sim.predictions[:,j] ≈ sol(t) rtol=1e-8 atol=1e-7
        end
        carry=sol(stop)
        feed=findfirst(f->f.time_h==stop,feeds)
        feed!==nothing && (carry=exchange(carry,feeds[feed].fraction,q.glucose))
        start=stop
    end
    # A different D1 aggregate measurement must not propagate backwards into any earlier state.
    other = simulate_aggregate_dayone(6,p,merge(q,(aggregate_day1=300.,)))
    @test other.predictions[1:3,1:3]==sim.predictions[1:3,1:3]
    @test other.predictions[:,1:2]==sim.predictions[:,1:2]
    @test other.predictions[4,3]==300.
    @test !isapprox(other.predictions[4,end],sim.predictions[4,end];rtol=1e-4)
    # Aggregate ODE parameters cannot affect the first day while the size input is fixed.
    otherp=merge(p,(alpha=2p.alpha,K_eps=2p.K_eps))
    otherparam=simulate_aggregate_dayone(6,otherp,q)
    @test otherparam.predictions[:,1:3]==sim.predictions[:,1:3]
    @test_throws ErrorException simulate_aggregate_dayone(6,p,merge(q,(u0=[q.u0[1:3];20.],)))
    # Exercise a D1 feeding event even though the main experiments feed on later days.
    qfeed = merge(q,(feeds=[(time_h=24.,fraction=.5),(time_h=48.,fraction=.8)],))
    fed = simulate_aggregate_dayone(6,p,qfeed)
    expected = exchange(fed.initialization.before,.5,q.glucose)
    expected[4]=q.aggregate_day1
    @test fed.initialization.after ≈ expected rtol=1e-12
    @test fed.predictions[1:3,3]==fed.initialization.before[1:3]
    @test [e.time_h for e in fed.events]==[24.,48.]
    for e in fed.events
        @test e.after[1]==e.before[1]
        @test e.after[4]==e.before[4]
        @test e.after[2] ≈ (1-e.fraction)*e.before[2]+e.fraction*q.glucose
        @test e.after[3] ≈ (1-e.fraction)*e.before[3]
    end
end

@testset "No pre-D1 aggregate curve or scored initial aggregate measurement" begin
    base=prep.held
    ctx=aggregate_dayone(base)
    points=NamedTuple[];curves=NamedTuple[];initials=NamedTuple[]
    AggregateDayOne.export_cohort!(points,curves,initials,base,ctx,p,AggregateDayOne.DAYONE,"selection40")
    agg=filter(r->r.state=="aggregate",curves)
    @test minimum(r.day for r in agg)==1.
    @test minimum(r.day for r in curves if r.state=="cells")==0.
    @test all(r->r.role=="initial_condition",filter(r->r.state=="aggregate" && r.day==1,points))
    @test all(r->r.role=="scored",filter(r->r.state!="aggregate",points))
    @test all(r->r.before==15.,filter(r->r.state=="aggregate",initials))
end

@testset "Matched selection score and unchanged original result" begin
    original_nll=gaussian_nll(prep.held.y,prediction(prep.held,6,p),prep.held.sigma)
    @test round(bic(original_nll,11,length(prep.held.y));digits=2)==1458.73
    ref=comparison_context(prep.held)
    ctx=aggregate_dayone(prep.held)
    @test length(ref.y)==length(ctx.y)
    @test length(free_names(ref,6))==length(free_names(ctx,6))==11
    @test all(r->r.state!="aggregate" || r.day>1,ref.observations)
    @test all(r->r.day==1 && r.state!="aggregate",filter(r->r.day==1,ctx.observations))
end
