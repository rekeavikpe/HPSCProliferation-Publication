using CSV,TOML
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
fit,split=ARGS
cfg=TOML.parsefile(joinpath(fit,"calibration.toml"));ctx=context(cfg=cfg)
isempty(cfg["fixed"])||error("Reference must be all-free")
source=collect(CSV.File(joinpath(fit,"observations.csv")))
[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in source]==[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in ctx.observations]||error("Reference observations changed")
for id in (0,5)
 m=legacy_model(id);p=read_parameters(joinpath(fit,"model$id/selected_parameters.csv"),m)
 nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
 saved=only(filter(r->r.model==id,collect(CSV.File(joinpath(fit,"ranking.csv")))))
 isapprox(nll,saved.nll;atol=1e-6,rtol=0)||error("Reference likelihood changed for model$id")
end
println("Saved all-100 mL observations and selected-fit likelihoods verified.")
