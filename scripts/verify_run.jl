using CSV, Test, TOML
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
length(ARGS)==1 || error("Usage: verify_run.jl CALIBRATION_DIRECTORY")
dir=abspath(ARGS[1]);isfile(joinpath(dir,"COMPLETE"))||error("Run incomplete")
ctx=context(cfg=TOML.parsefile(joinpath(dir,"calibration.toml")))
@testset "Saved fresh calibration outputs" begin
    for id in 0:6
        m=legacy_model(id);folder=joinpath(dir,"model$id")
        selected=only(CSV.File(joinpath(folder,"summary.csv")))
        p=read_parameters(joinpath(folder,"selected_parameters.csv"),m)
        pred=prediction(ctx,m,p)
        nll=gaussian_nll(ctx.y,pred,ctx.sigma)
        @test isapprox(nll,selected.nll;atol=1e-6,rtol=0)
        @test isapprox(bic(nll,length(p),length(ctx.y)),selected.bic;atol=2e-6,rtol=0)
        tight=prediction(ctx,m,p;reltol=1e-13,abstol=1e-14)
        @test maximum(abs.((tight.-pred)./ctx.sigma))<1e-5
        stages=Float64[]
        for start in filter(f->startswith(f,"run"),readdir(folder))
            for row in CSV.File(joinpath(folder,start,"stages.csv"))
                q=read_parameters(joinpath(folder,start,row.stage*"_parameters.csv"),m)
                actual=gaussian_nll(ctx.y,prediction(ctx,m,q),ctx.sigma)
                @test isapprox(actual,row.nll;atol=1e-6,rtol=0)
                push!(stages,actual)
            end
        end
        @test selected.nll<=minimum(stages)+1e-6
    end
end
