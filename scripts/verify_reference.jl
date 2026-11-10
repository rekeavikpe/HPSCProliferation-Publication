using Test, CSV, TOML, SHA
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
ctx=context()
reference=joinpath(ROOT,"reference/scale_transfer")
checks=NamedTuple[]
@testset "Reproduce archived training scores and scale predictions" begin
    @test length(ctx.y)==116
    @test sort(collect(keys(ctx.conditions)))==[1,2,3,5,6]
    for public_id in 0:6
        m=legacy_model(public_id)
        dir=joinpath(reference,"legacy_model$m")
        expected=only(CSV.File(joinpath(dir,"summary.csv")))
        p=read_parameters(joinpath(dir,"selected_parameters.csv"),m)
        nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
        @test isapprox(nll,expected.nll;atol=1e-6,rtol=0)
        @test isapprox(bic(nll,length(p),116),expected.bic;atol=2e-6,rtol=0)
        push!(checks,(model=public_id,nll=nll,reference_nll=expected.nll,absolute_difference=abs(nll-expected.nll)))
    end
    testctx=context(conditions=(4,))
    for id in (0,5)
        m=legacy_model(id)
        p=read_parameters(joinpath(reference,"legacy_model$m/selected_parameters.csv"),m)
        pred=prediction(testctx,m,p)
        saved=filter(r->r.publication_model==id,collect(CSV.File(joinpath(reference,"test_predictions.csv"))))
        for (i,r) in enumerate(testctx.observations)
            expected=only(filter(s->s.day==r.day && s.state==r.state,saved))
            @test isapprox(pred[i],expected.predicted;rtol=1e-8,atol=1e-6)
        end
    end
end
out=joinpath(ROOT,"results/verification");mkpath(out)
CSV.write(joinpath(out,"reference_scores.csv"),checks)
println("Reference reproduction passed; this verifies saved parameters, not a new global optimum.")
