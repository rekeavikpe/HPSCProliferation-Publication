using Test, TOML, CSV
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("../ablations/RoutedExposure.jl")
using .RoutedExposure
cfg=TOML.parsefile(joinpath(only(ARGS),"calibration.toml"))
@testset "Training-only downstream cohort" begin
 train=context(;cfg)
 @test length(cfg["calibration_cohort"]["series_ids"])==19
 @test length(cfg["calibration_cohort"]["selection_ids"])==12
 @test isempty(intersect(cfg["calibration_cohort"]["series_ids"],cfg["calibration_cohort"]["selection_ids"]))
 allcfg=deepcopy(cfg);delete!(allcfg,"calibration_cohort")
 @test train.y != context(;cfg=allcfg).y
 @test context(;cfg,conditions=(4,)).y==context(;cfg=allcfg,conditions=(4,)).y
 held=context(;cfg,ids=cfg["calibration_cohort"]["selection_ids"],fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
 @test held.y != train.y
 p=read_parameters(joinpath(only(ARGS),"model5/selected_parameters.csv"),6)
 fixed=deepcopy(cfg);fixed["fixed"]=Dict("K_G"=>p.K_G,"K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX)
 ctx=context(;cfg=fixed)
 @test ctx.y==train.y
 @test ctx.sigma==train.sigma
 @test length(free_names(ctx,6))==8
 @test prediction(ctx,6,p)==prediction(train,6,p)
 for scenario in ("shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both")
  routed=routed_context(ctx,scenario;shared=false)
  @test routed.y==train.y
  @test routed.sigma==train.sigma
 end
end
