using CSV, Test, TOML, SHA
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
source=joinpath(ROOT,"reference/figure_data")
fitdir=joinpath(ROOT,"reference/reproduction")
curves=collect(CSV.File(joinpath(source,"training_curves.csv")))
ctx=context();winner=TOML.parsefile(joinpath(source,"source.toml"))["selected_model"]
@testset "Plotted sampling values and exchange events" begin
    for id in unique([0,winner])
        model=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
        predicted=prediction(ctx,model,p)
        for (i,r) in enumerate(ctx.observations)
            row=only(filter(x->x.model==id&&x.condition==r.condition&&x.day==r.day&&x.state==r.state&&x.side=="before_or_between_exchange",curves))
            @test isapprox(row.value,predicted[i];rtol=1e-9,atol=1e-9)
        end
        for (c,q) in ctx.conditions,event in q.feeds,state in STATES
            pair=filter(r->r.model==id&&r.condition==c&&r.day==event.time_h/24&&r.state==state,curves)
            @test length(pair)==2
            before=only(filter(r->r.side=="before_or_between_exchange",pair)).value
            after=only(filter(r->r.side=="after_exchange",pair)).value
            expected=state=="glucose" ? (1-event.fraction)*before+event.fraction*q.glucose : state=="lactate" ? (1-event.fraction)*before : before
            @test isapprox(after,expected;rtol=1e-12,atol=1e-12)
        end
    end
end
