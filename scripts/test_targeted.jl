using Test, CSV
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration
@testset "Targeted-check pipeline and failed-start accounting" begin
    fitdir=joinpath(ROOT,"reference/reproduction")
    ctx=checked_context(fitdir);m=legacy_model(5)
    p=read_parameters(joinpath(fitdir,"model5/selected_parameters.csv"),m)
    cfg=deepcopy(ctx.config);cfg["lactate_inhibition"]=false;cfg["fixed"]["K_Lac"]=1.
    inactive=context(cfg=cfg)
    @test length(free_names(inactive,m))==10
    @test :K_Lac ∉ free_names(inactive,m)
    a=local_starts(inactive,m,p;n=2,seed=7)
    @test a==local_starts(inactive,m,p;n=2,seed=7)
    mktempdir() do dir
        jobs=[(ctx=inactive,model=m,z=z,dir=joinpath(dir,"run$i")) for (i,z) in enumerate(a)]
        run_jobs(jobs;lm_iterations=2,polish_iterations=5)
        @test all(isfile(joinpath(j.dir,"COMPLETE")) for j in jobs)
        mkpath(joinpath(dir,"run_failed"));write(joinpath(dir,"run_failed/FAILED"),"intentional test fixture")
        result=save_choice(inactive,m,dir;estimated_parameters=10)
        @test isfinite(result.nll)
        @test result.k==10
        @test result.failed_starts==1
        @test result.nll<=result.minimum_nll+1e-6
    end
end
