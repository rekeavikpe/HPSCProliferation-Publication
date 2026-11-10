using Test
include(joinpath(@__DIR__,"..","src","BoundedLeastSquares.jl"))
using .BoundedLeastSquares
@testset "Bounded least squares: known solutions" begin
    for (target,expected) in [([0.3,0.7],[0.3,0.7]),([-0.5,1.5],[0.0,1.0]),([0.0,0.8],[0.0,0.8])]
        x,r=bounded_lm(z->z.-target,[0.8,0.2])
        @test maximum(abs.(x.-expected))<1e-4
        @test r.converged
    end
    x,r=bounded_lm(z->[10*(z[2]-z[1]^2),1-z[1]],[0.2,0.7])
    @test maximum(abs.(x.-[1,1]))<1e-3
    @test all(0 .<= x .<= 1)
end
