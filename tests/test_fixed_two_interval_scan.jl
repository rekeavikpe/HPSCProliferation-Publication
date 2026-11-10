using Test, ProfileLikelihood
include("../profiles/PackageWorkflow.jl")
using .PackageWorkflow

@testset "Threshold scan retains the confidence interval on a known likelihood" begin
    # Correlated coordinates require nuisance reoptimization. The profile of
    # x has known SD 0.05 and a closed 95% confidence interval inside [0,1].
    loglik(z, _) = -0.5 * (((z[1]-0.5)/0.05)^2 + ((z[2]-0.4-0.7*(z[1]-0.5))/0.1)^2)
    prob=make_problem(loglik, [0.5,0.4], [:x,:y])
    sol=fit_mle(prob;maxiters=100)
    seeds = _ -> [copy(sol.mle[2:2])]
    short, calls=profile_parameter(prob,sol,1;resolution=200,maxiters=1000,parallel=false,full_range=false,seed_provider=seeds)
    full, _=profile_parameter(prob,sol,1;resolution=200,maxiters=1000,parallel=false,full_range=true,seed_provider=seeds)
    delta=sqrt(-2*ProfileLikelihood.get_chisq_threshold(0.95))*0.05
    ci=short.confidence_intervals[1]
    @test isapprox(ci.lower,0.5-delta;atol=2e-6)
    @test isapprox(ci.upper,0.5+delta;atol=2e-6)
    @test isapprox(ci.lower,full.confidence_intervals[1].lower;atol=2e-6)
    @test isapprox(ci.upper,full.confidence_intervals[1].upper;atol=2e-6)
    @test length(short.parameter_values[1]) < length(full.parameter_values[1])
    @test minimum(short.parameter_values[1]) > 0
    @test maximum(short.parameter_values[1]) < 1
    @test first(short.profile_values[1]) < ProfileLikelihood.get_chisq_threshold(0.95)
    @test last(short.profile_values[1]) < ProfileLikelihood.get_chisq_threshold(0.95)
    @test all(r -> r.starts==2, calls)
    @test all(r -> r.stopping_met, calls)
    @test maximum(abs(short.other_mles[1][j][1]-(0.4+0.7*(short.parameter_values[1][j]-0.5))) for j in eachindex(short.parameter_values[1])) < 1e-5
end
