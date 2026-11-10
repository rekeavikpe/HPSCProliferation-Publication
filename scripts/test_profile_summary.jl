using Test
include("../src/ProfileSummary.jl")
using .ProfileSummary
row(x,y,ok=true)=(fixed_value=x,twice_delta_nll=y,stopping_met=ok)
@testset "Profile interval safeguards" begin
    rows=[row(1.,6.),row(10.,0.),row(100.,6.)]
    s=summarize_profile(rows,10.)
    @test 1<s.lower.estimate<10
    @test 10<s.upper.estimate<100
    @test !s.better_minimum
    @test ismissing(summarize_profile([row(1.,0.),row(10.,0.)],1.).upper.estimate)
    @test ismissing(summarize_profile([row(1.,6.,false),row(10.,0.)],10.).lower.estimate)
    @test summarize_profile([row(1.,-.1),row(10.,0.)],10.).better_minimum
end
