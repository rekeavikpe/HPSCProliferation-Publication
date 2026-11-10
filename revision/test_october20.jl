using Test
include("October20.jl")
length(ARGS)==2||error("Usage: test_october20.jl FIT SPLIT")
fit,split=ARGS;prep,external=source_contexts(fit,split)
p=read_parameters(joinpath(fit,"model5/selected_parameters.csv"),6)
@testset "October 20 analysis design" begin
    @test length(prep.split.training)==19
    @test length(prep.split.selection)==12
    @test isempty(intersect(prep.split.training,prep.split.selection))
    @test [legacy_model(i) for i in 0:6]==[0,1,3,4,5,6,7]
    d=dayone(prep.train)
    @test all(r.day>1 for r in d.observations)
    @test length(d.y)==length(d.sigma)==length(d.observations)
    @test all(q.times[1]>24 for q in values(d.conditions))
    # Exact continuation test detects incorrectly restarted death clock or feed timing.
    qs=copy(prep.train.conditions)
    for (c,q) in qs
        sim=simulate(6,p,q.u0,q.times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,reltol=1e-12,abstol=1e-13)
        qs[c]=merge(d.conditions[c],(u0=sim.predictions[:,findfirst(==(24.),q.times)],))
    end
    continuation=DayOneContext(d.config,d.observations,qs,d.y,d.sigma)
    ix=findall(r->r.day>1,prep.train.observations)
    @test maximum(abs.(prediction(continuation,6,p).-prediction(prep.train,6,p)[ix])./d.sigma)<1e-5
    for s in EIGHT
        ctx=dayone(route(prep.train,s));m=s=="baseline" ? 0 : 6
        @test length(ctx.y)==length(d.y)
        @test ctx.sigma==d.sigma
        @test all(isfinite,prediction(ctx,m,m==0 ? read_parameters(joinpath(fit,"model0/selected_parameters.csv"),0) : p))
    end
    cfg=deepcopy(prep.cfg);cfg["fixed"]=Dict(string(n)=>getproperty(p,n) for n in keys(p) if !(n in REQUESTED5))
    for line in ("585A1","SCTi003-A")
        ctx=cell_context(line,cfg)
        @test Set(free_names(ctx,6))==Set(REQUESTED5)
        @test length(ctx.config["fixed"])==6
        @test all(isfinite,prediction(ctx,6,p))
    end
    zero=read_parameters(joinpath(fit,"model0/selected_parameters.csv"),0)
    @test prediction(prep.train,0,zero)==prediction(exposure_context(prep.train,"lag_p75"),0,zero)
end
