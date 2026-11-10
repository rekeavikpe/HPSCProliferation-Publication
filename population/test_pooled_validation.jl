using Test, Statistics
include("PooledValidation.jl")
using .PooledValidation, .PooledValidation.SplitSelection.Workflow
cfg=joinpath(ROOT,"config/calibration_lac100_kg100.toml")
ps=plans(); inv=inventory()
@testset "Pooled folds and training-only exclusions" begin
    @test length(inv)==46
    @test ps==plans()
    @test sort([length(p.selection) for p in ps if p.experiment=="pooled_full"])==[15,15,16]
    full=filter(p->p.experiment=="pooled_full",ps)
    @test sort(vcat([p.selection for p in full]...))==sort([r.series_id for r in inv])
    for p in ps
        @test isempty(intersect(p.training,p.selection))
        prep=prepare_fold(cfg,p)
        @test all(>(0),prep.held.sigma)
        @test all(isfinite,prep.held.y)
        altered=[merge(r,(mean=100r.mean,sd=1000r.sd)) for r in prep.held.observations]
        @test held_weights(prep.train,altered)[1]==prep.held.sigma
        if startswith(p.experiment,"pooled")
            @test length(prep.held.conditions)==8
            @test length(prep.train.conditions)==8
            for c in (7,8)
                @test [r.fraction for r in prep.train.conditions[c].feeds]==[.5, c==7 ? .2 : .5,c==7 ? .2 : .5,c==7 ? .2 : .5]
            end
        end
    end
    for p in full
        trim=only(filter(q->q.name==replace(p.name,"_full"=>"_trimmed"),ps))
        @test trim.selection==p.selection
        decisions=exclusions(p.training)
        @test all(r->r.series_id in p.training,decisions)
        @test all(r->r.condition<=6 || !startswith(r.decision,"excluded"),decisions)
        @test Set(setdiff(p.training,trim.training))==Set(r.series_id for r in decisions if startswith(r.decision,"excluded"))
    end
    # With only original data, pooled context must reproduce the established workflow.
    p=only(filter(p->p.name=="reference100",ps));prep=prepare_fold(cfg,p)
    old=context(ids=p.training,cfg=prep.cfg)
    @test old.y==prep.train.y
    @test old.sigma==prep.train.sigma
    for c in keys(old.conditions)
        @test old.conditions[c].u0==prep.train.conditions[c].u0
        @test old.conditions[c].feeds==prep.train.conditions[c].feeds
    end
end
