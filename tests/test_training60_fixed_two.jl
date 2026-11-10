using Test, SHA
include("../profiles/run_fixed_two.jl")
length(ARGS)==1 || error("Supply the saved 19-culture reference fit")
fit=abspath(only(ARGS))
inputs=[joinpath(fit,p) for p in ("calibration.toml","split.csv","observations.csv","model5/selected_parameters.csv","model5/summary.csv")]
before=sha256.(read.(inputs))
p=fixed_two_context(fit)
@testset "Only K_Lac and Y_GX fixed, K_G free in nine-parameter analysis" begin
    @test Set(keys(p.fixed))==Set(["K_Lac","Y_GX"])
    @test p.fixed["K_Lac"]==p.p.K_Lac
    @test p.fixed["Y_GX"]==p.p.Y_GX
    @test length(p.names)==9
    @test :K_G in p.names
    @test p.raw.config["fixed"]==Dict()
    @test p.ctx.config["bounds"]==p.raw.config["bounds"]
    @test p.ctx.observations==p.raw.observations
    @test p.ctx.y==p.raw.y && p.ctx.sigma==p.raw.sigma
    @test length(p.ctx.y)==116
    @test p.ctx.config["exposure_percentile"]==50
    @test Set(keys(p.ctx.conditions))==Set([1,2,3,5,6])
    @test length(p.ctx.config["calibration_cohort"]["series_ids"])==19
    @test length(p.ctx.config["calibration_cohort"]["selection_ids"])==12
    @test isapprox(gaussian_nll(p.ctx.y,prediction(p.ctx,p.model,p.p),p.ctx.sigma),p.saved_nll;atol=1e-6,rtol=0)
    z=encode(p.ctx,p.model,p.p)
    q=decode(p.ctx,p.model,z)
    @test q.K_Lac==p.fixed["K_Lac"] && q.Y_GX==p.fixed["Y_GX"]
    z[findfirst(==(:K_G),p.names)]+=.001
    q=decode(p.ctx,p.model,z)
    @test q.K_G!=p.p.K_G
    @test q.K_Lac==p.fixed["K_Lac"] && q.Y_GX==p.fixed["Y_GX"]
    @test prediction(p.ctx,p.model,q)!=prediction(p.ctx,p.model,p.p)
    for target in p.names
        @test length(setdiff(p.names,[target]))==8
        target==:K_G || @test :K_G in setdiff(p.names,[target])
    end
end
@testset "Source preservation and protected paths" begin
    @test sha256.(read.(inputs))==before
    @test_throws ErrorException safe_output(fit,[fit])
    @test_throws ErrorException safe_output(joinpath(fit,"new_profiles"),[fit])
    @test_throws ErrorException safe_output(dirname(fit),[fit])
    mktempdir() do dir
        symlink(fit,joinpath(dir,"alias"))
        @test_throws ErrorException safe_output(joinpath(dir,"alias","new"),[fit])
        @test safe_output(joinpath(dir,"profiles"),[fit])==joinpath(realpath(dir),"profiles")
    end
end
