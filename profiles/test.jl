using Test, ProfileLikelihood, SciMLBase
include("PackageWorkflow.jl")
using .PackageWorkflow
@testset "Original package workflow" begin
    ll=(z,_)->-((z[1]-.4)^2+(z[2]-.6)^2)/(2*.05^2)
    prob=make_problem(ll,[.45,.55],[:a,:b]); sol=fit_mle(prob;maxiters=100)
    @test maximum(abs.(sol.mle.-[.4,.6]))<1e-4
    prof,records=profile_parameter(prob,sol,1;resolution=20,maxiters=100,parallel=true)
    @test !isempty(records)
    @test all(r->r.stopping_met,records)
    @test isapprox(prof.confidence_intervals[1].lower,.4-1.95996398454*.05;atol=.002)
    @test isapprox(prof.confidence_intervals[1].upper,.4+1.95996398454*.05;atol=.002)
    @test all(isapprox(v,-(x-.4)^2/(2*.05^2);atol=1e-5) for (x,v) in zip(prof.parameter_values[1],prof.profile_values[1]))
    cache=[-eps(),1+eps()]
    inset!(cache,[0.,0.],[1.,1.])
    @test all(0 .< cache .< 1)
end
