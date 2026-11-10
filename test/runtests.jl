using Test
include("../src/Proliferation.jl")
using .Proliferation

base = [6.174043e7, .1349, 9.634, 30.40495, 2.932e-6, .443, 3.62007e-6, 62.45, .1116]
function fixture(model)
    extra = model == 0 ? [] : model == 6 ? [.00579, .004024] : model == 7 ? [.00579, .01, .004024, .5] : [.00579, .01, .004024]
    parameters(model, vcat(base, extra))
end

@testset "Equation limits and conservation" begin
    u = [80000., 14., 1., 100.]
    d0, d = zeros(4), zeros(4)
    model_rhs!(d0, u, fixture(0), 24., 0)
    for m in UNIQUE_MODELS
        model_rhs!(d, u, fixture(m), 24., m, 0., 0.)
        @test d ≈ d0
    end
    p = fixture(6)
    model_rhs!(d, u, p, 24., 6, .004, .0003)
    growth = d[1]/u[1]+p.K_d*exp(-2)
    @test d[2]+p.Y_GLac*d[3] ≈ -p.Y_GX*growth*u[1]
    @test d[4] ≈ growth*(p.alpha-.0003/p.K_eps*u[4])
    for tau in (0., .001, .1), edr in (0., .001)
        model_rhs!(d0,u,fixture(1),24.,1,tau,edr)
        model_rhs!(d,u,fixture(2),24.,2,tau,edr)
        @test d == d0
    end
    model_rhs!(d0,u,fixture(4),24.,4,0.,.001)
    model_rhs!(d,u,fixture(4),24.,4,.01,.001)
    @test d[1] < d0[1]
    @test d[2:4] == d0[2:4]
end

@testset "Exact exchange and observation order" begin
    p = parameters(0, [1e8, 0., 10., 30., 1e-6, .5, 1e-6, 50., 0.])
    u = [80000., 4., 10., 100.]
    feeds = [(time_h=48., fraction=.5), (time_h=72., fraction=.8)]
    result = simulate(0,p,u,[0.,24.,48.,60.,72.,96.],feeds; fresh_glucose=14.)
    @test result.predictions[:,3] ≈ u
    @test result.predictions[:,4] ≈ [80000.,9.,5.,100.]
    @test result.predictions[:,5] ≈ [80000.,9.,5.,100.]
    @test result.predictions[:,6] ≈ [80000.,13.,1.,100.]
    @test length(result.events) == 2
    @test exchange(u,0.,14.) == u
    @test exchange(u,1.,14.) == [80000.,14.,0.,100.]
    @test_throws ArgumentError simulate(0,p,u,[48.,96.],vcat(feeds[1:1],feeds[1:1]);fresh_glucose=14.)
    @test_throws ArgumentError exchange(u,1.01,14.)
end

@testset "Exact inactive lactate inhibition limit" begin
    for m in UNIQUE_MODELS
        p=fixture(m);u=[80000.,14.,0.,100.];a=zeros(4);b=zeros(4)
        model_rhs!(a,u,p,24.,m,.004,.0003;lactate_inhibition=false)
        u[3]=100.
        model_rhs!(b,u,p,24.,m,.004,.0003;lactate_inhibition=false)
        @test a==b
        huge=merge(p,(K_Lac=1e14,))
        model_rhs!(b,u,huge,24.,m,.004,.0003)
        @test isapprox(a,b;rtol=1e-10)
        model_rhs!(b,u,p,24.,m,.004,.0003)
        @test b[1]<a[1]
    end
end

@testset "Analytic death trajectory and no clock reset" begin
    values = [1e8,.1,10.,30.,1e-6,.5,1e-6,50.,0.]
    p = parameters(0,values)
    times = [0.,12.,24.,48.,72.]
    result = simulate(0,p,[80000.,14.,0.,15.],times,[(time_h=24.,fraction=.5)];fresh_glucose=14.)
    expected = [80000*exp(-.1*12*(1-exp(-t/12))) for t in times]
    @test result.predictions[1,:] ≈ expected rtol=1e-8
end

@testset "Exposure and score" begin
    exposure = linear_exposure([200.,400.,600.],[.002,.004,.006],[.0001,.0003,.0005])
    @test collect(exposure(300.)) ≈ [.003,.0002]
    @test exposure(15.) == exposure(200.)
    @test exposure(1000.) == exposure(600.)
    @test_throws ArgumentError linear_exposure([200.,200.],[1.,2.],[1.,2.])
    @test gaussian_nll([1.,2.],[1.,2.],[1.,1.]) ≈ log(2pi)
    @test gaussian_nll([1.],[3.],[2.]) ≈ log(2)+log(2pi)/2+.5
    @test_throws ArgumentError gaussian_nll([1.],[1.],[0.])
    @test bic(10.,3,20) ≈ 20+3log(20)
end
