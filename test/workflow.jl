using Test, CSV
include("../src/Workflow.jl")
include("../src/ReplicateSplits.jl")
using .Workflow, .ReplicateSplits
@testset "Whole-culture random splitting and leakage checks" begin
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    s=random_split(series)
    @test s==random_split(series)
    @test s!=random_split(series;seed=20261004)
    @test isempty(intersect(s.training,s.selection))
    @test length(s.training)+length(s.selection)==31
    @test all(r->r.legacy_condition!=4,s.assignments)
    for c in (1,2,3,5,6)
        @test any(r->r.legacy_condition==c&&r.role=="calibration",s.assignments)
        @test any(r->r.legacy_condition==c&&r.role=="selection",s.assignments)
    end
    a=context(ids=s.training);b=context(ids=s.selection)
    @test length(a.y)==length(b.y)==116
    @test a.y!=b.y
    @test all(isfinite,a.sigma)
    @test a.conditions[1].u0[1]==20000
    @test legacy_model(5)==6
    @test legacy_model(6)==7
end
@testset "Predetermined split feasibility" begin
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    for seed in 20261003:20261012
        s=random_split(series;seed)
        a=context(ids=s.training)
        b=context(ids=s.selection,fallback_glucose=Dict(c=>q.glucose for (c,q) in a.conditions))
        @test all(isfinite,a.y)
        @test all(isfinite,b.y)
        @test all(>(0),a.sigma)
    end
end
