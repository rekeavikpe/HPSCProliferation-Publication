using Test, CSV
include("CrossValidation.jl")
using .CrossValidation, .CrossValidation.SplitSelection.Workflow
cfg=joinpath(ROOT,"config/calibration_lac100_kg100.toml")
ps=plans(); allps=plans(design="all126")
@testset "Condition and culture holdouts" begin
    @test length(ps)==15
    @test length(allps)==132
    @test ps==plans()
    large=ps[1].selection
    partial=filter(p->p.experiment=="partial500",ps)
    @test all(id->count(p->id in p.training,partial)==4,large)
    @test all(id->count(p->id in p.selection,partial)==5,large)
    @test length(unique(p.training for p in partial))==9
    for p in ps
        @test isempty(intersect(p.training,p.selection))
        prep=prepare_fold(cfg,p)
        @test all(isfinite,prep.held.sigma)&&all(>(0),prep.held.sigma)
        @test isempty(prep.cfg["fixed"])
        if p.experiment=="loco"
            @test length(prep.train.conditions)==4
            @test length(prep.held.conditions)==1
            @test isempty(intersect(keys(prep.train.conditions),keys(prep.held.conditions)))
            @test all(==("training_state_rms"),prep.weight_source)
        elseif p.experiment=="partial500"
            @test length(p.training)==35 && length(p.selection)==5
            @test length(prep.train.conditions)==6
        end
        # Held-out means/SDs must not affect scoring weights.
        altered=[merge(r,(mean=100r.mean,sd=1000r.sd)) for r in prep.held.observations]
        @test held_weights(prep.train,altered)[1]==prep.held.sigma
        # Perfect predictions have zero MRPE; noise remains descriptive replicate SD.
        @test all(r->r.mrpe_percent==0,metrics(prep.held,prep.held.y))
    end
end
