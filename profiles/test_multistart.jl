using Test
include("PackageWorkflow.jl")
using .PackageWorkflow
@testset "Multistart profiling and local CI brackets" begin
    p=make_problem((x,d)->-100sum(abs2,x.-0.5),[0.5,0.5],[:a,:b])
    sol=fit_mle(p)
    seeds=prob->begin
        @assert hasproperty(prob.f.f,Symbol("θₙ"))
        [[.1],[.9]]
    end
    prof,records=profile_parameter(p,sol,1;resolution=12,parallel=false,full_range=true,seed_provider=seeds)
    @test minimum(prof.parameter_values[1])==0
    @test maximum(prof.parameter_values[1])==1
    @test abs(prof.confidence_intervals[1].lower-(.5-sqrt(1.920729410347062/100)))<.01
    @test abs(prof.confidence_intervals[1].upper-(.5+sqrt(1.920729410347062/100)))<.01
    @test all(r.starts==3 for r in records)
    @test all(r.objective<=r.best_seed+1e-6 for r in records)
end
