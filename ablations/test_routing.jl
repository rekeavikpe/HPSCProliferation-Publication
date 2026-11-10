using Test, TOML
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("RoutedExposure.jl")
using .RoutedExposure
cfg=TOML.parsefile(joinpath(ROOT,"config/calibration_lac100_kg100.toml"))
ctx=context(cfg=cfg);p=decode(ctx,6,fill(.5,11))
@testset "Exposure routing, limits, and parameter accounting" begin
    @test length(free_names(ctx,0))==9
    for (scenario,(cell,agg)) in ROUTES
        scenario=="baseline"&&continue
        routed=routed_context(ctx,scenario)
        expected=9+(cell!="none")+(agg!="none")
        @test length(free_names(routed,6))==expected
        @test ctx.y==routed.y && ctx.sigma==routed.sigma
        for (c,q) in ctx.conditions, diameter in (15.,200.,350.,900.)
            tau,edr=q.lookup(diameter)
            expectedcell=cell=="shear" ? tau : cell=="edr" ? edr : 0.
            expectedagg=agg=="shear" ? tau : agg=="edr" ? edr : 0.
            @test routed.conditions[c].lookup(diameter)==(expectedcell,expectedagg)
        end
        @test routed.config["bounds"]["K_tau"]==cfg["bounds"][cell=="edr" ? "K_eps" : "K_tau"]
        @test routed.config["bounds"]["K_eps"]==cfg["bounds"][agg=="shear" ? "K_tau" : "K_eps"]
    end
    @test prediction(ctx,6,p)≈prediction(routed_context(ctx,"full"),6,p)
    p0=parameters(0,[getproperty(p,n) for n in parameter_names(0)])
    @test prediction(ctx,0,p0)≈prediction(routed_context(ctx,"baseline"),6,p)
    for u in ([2e5,15.,3.,150.],[8e5,4.,10.,350.])
        baseline=zeros(4);cell=zeros(4);agg=zeros(4)
        model_rhs!(baseline,u,p,36.,6,0.,0.)
        model_rhs!(cell,u,p,36.,6,.003,0.)
        model_rhs!(agg,u,p,36.,6,0.,.003)
        @test agg[1:3]≈baseline[1:3]
        @test agg[4]<baseline[4]
        @test cell[1]<baseline[1]
        @test cell[4]<baseline[4] # Shared growth retains an indirect aggregate response.
    end
end

@testset "Shared response parameter" begin
    for scenario in ("shear_both","edr_both")
        shared=routed_context(ctx,scenario;shared=true)
        separate=routed_context(ctx,scenario)
        @test length(free_names(shared,6))==10
        @test :K_eps ∉ free_names(shared,6)
        for fraction in (.2,.5,.8)
            q=decode(shared,6,fill(fraction,10))
            @test q.K_tau==q.K_eps
            @test collect(values(decode(shared,6,encode(shared,6,q))))≈collect(values(q))
            @test prediction(shared,6,q)≈prediction(separate,6,q)
        end
        cfg=deepcopy(ctx.config)
        cfg["fixed"]=Dict("K_G"=>p.K_G,"K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX)
        conditional=routed_context(context(cfg=cfg),scenario;shared=true)
        @test length(free_names(conditional,6))==7
        q=decode(conditional,6,fill(.5,7))
        @test q.K_tau==q.K_eps
        @test q.K_G==p.K_G && q.K_Lac==p.K_Lac && q.Y_GX==p.Y_GX
        @test response_labels(scenario;shared=true)["K_tau"] in ("K_{\\tau}","K_{\\varepsilon}")
    end
    @test length(free_names(routed_context(ctx,"full";shared=true),6))==11
end
