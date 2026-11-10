using Test, CSV, TOML
include("SplitSelection.jl")
using .SplitSelection
using .SplitSelection.Workflow
const cfgfile=joinpath(ROOT,"config/calibration_lac100_kg100.toml")
p=prepare(cfgfile);q=prepare(cfgfile)
@testset "Whole-culture split and training-only weights" begin
    @test p.split.assignments==q.split.assignments
    @test length(p.split.training)==19
    @test length(p.split.selection)==12
    @test isempty(intersect(p.split.training,p.split.selection))
    @test length(unique(vcat(p.split.training,p.split.selection)))==31
    for c in (1,2,3,5,6),role in ("calibration","selection")
        @test count(r->r.legacy_condition==c&&r.role==role,p.split.assignments)>0
    end
    obskey(r)=(r.condition,r.day,r.state)
    scales=Dict(obskey(r)=>s for (r,s) in zip(p.train.observations,p.train.sigma))
    @test all(s==scales[obskey(r)] for (r,s) in zip(p.held.observations,p.held.sigma))
    @test all(isfinite,p.held.sigma)&&all(>(0),p.held.sigma)
    @test all(r->r.n_vessels>=2,p.train.observations)
    raw=context(ids=p.split.selection,cfg=p.cfg)
    @test raw.y==p.held.y
    @test raw.sigma!=p.held.sigma
    @test all(p.held.conditions[c].u0==raw.conditions[c].u0 for c in keys(raw.conditions))
    mktempdir() do dir
        @test !checkpoint(p.train,6,dir)
        cfg=deepcopy(p.cfg);cfg["fixed"]["K_Lac"]=100.
        path=joinpath(dir,"fixed.toml");open(io->TOML.print(io,cfg),path,"w")
        @test_throws ErrorException prepare(path)
    end
end
