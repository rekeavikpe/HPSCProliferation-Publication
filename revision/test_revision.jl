using Test
include("run_david_analyses.jl")
@testset "Revision design and scoring" begin
 cfg=TOML.parsefile(joinpath(ROOT,"config/calibration_lac100_kg100.toml"))
 prep=SplitSelection.prepare(joinpath(ROOT,"config/calibration_lac100_kg100.toml"))
 @test length(prep.split.training)==19
 @test length(prep.split.selection)==12
 @test isempty(intersect(prep.split.training,prep.split.selection))
 for scenario in EXPOSURES
  a=exposure_context(prep.train,scenario);b=exposure_context(prep.held,scenario)
  @test a.y==prep.train.y && a.sigma==prep.train.sigma
  @test b.y==prep.held.y && b.sigma==prep.held.sigma
  for (c,q) in a.conditions
   @test all(isfinite,q.lookup(15.))
   @test all(>=(0),q.lookup(300.))
   startswith(scenario,"eul_")&&(@test q.lookup(15.)==q.lookup(800.))
   scenario=="lag_p50"&&(@test q.lookup(300.)==prep.train.conditions[c].lookup(300.))
  end
 end
 @test length(BIOLOGICAL7)==7 && length(METABOLIC5)==5
 @test Set(METABOLIC5)==Set((:K_G,:K_Lac,:K_oG,:Y_GX,:Y_GLac))
 for line in ("585A1","SCTi003-A")
  ctx=cell_context(line,cfg)
  @test length(ctx.y)==(line=="585A1" ? 16 : 20)
  @test only(errors(ctx,ctx.y)[end:end]).mrpe_percent==0.
  @test only(errors(ctx,1.1ctx.y)[end:end]).mrpe_percent≈10.
  @test all(>(0),ctx.sigma)
  @test first(ctx.conditions[1].feeds).fraction==(line=="585A1" ? .8 : .5)
 end
 for s in ("shear_both","edr_both")
  @test length(free_names(RoutedExposure.routed_context(prep.train,s;shared=false),6))==11
  @test length(free_names(RoutedExposure.routed_context(prep.train,s;shared=true),6))==10
 end
end
