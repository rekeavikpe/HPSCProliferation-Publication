using CSV, TOML, Statistics
include(joinpath(@__DIR__,"..","population","SplitSelection.jl"))
using .SplitSelection, .SplitSelection.Workflow, .SplitSelection.Workflow.Calibration, .SplitSelection.Workflow.Calibration.Proliferation
function verify()
 fit=joinpath(@__DIR__,"cache/training60/reference_fit")
 cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
 prep=SplitSelection.prepare(joinpath(fit,"calibration.toml")); train=prep.train;held=prep.held
 @assert Set(prep.split.training)==Set(cfg["calibration_cohort"]["series_ids"])
 @assert Set(prep.split.selection)==Set(cfg["calibration_cohort"]["selection_ids"])
 @assert length(cfg["calibration_cohort"]["series_ids"])==19
 @assert length(cfg["calibration_cohort"]["selection_ids"])==12
 @assert isempty(intersect(cfg["calibration_cohort"]["series_ids"],cfg["calibration_cohort"]["selection_ids"]))
 @assert length(train.y)==length(held.y)==116
 @assert cfg["exposure_percentile"]==50
 scores=collect(CSV.File(joinpath(fit,"selection_scores.csv")))
 mapping=Dict(0=>0,1=>1,2=>3,3=>4,4=>5,5=>6,6=>7)
 results=NamedTuple[];point_rows=NamedTuple[]
 for id in 0:6
  model=legacy_model(id);p=read_parameters(joinpath(fit,"model$id/selected_parameters.csv"),model)
  for (cohort,ctx) in (("training",train),("selection",held))
   yp=prediction(ctx,model,p);nll=gaussian_nll(ctx.y,yp,ctx.sigma)
   expected=cohort=="training" ? only(CSV.File(joinpath(fit,"model$id/summary.csv"))) : only(filter(r->r.storage_model==id,scores))
   @assert isapprox(nll,expected.nll;atol=1e-6,rtol=0) (id,cohort,nll,expected.nll)
   score=bic(nll,expected.k,length(yp));@assert isapprox(score,expected.bic;atol=2e-6,rtol=0)
   for state in (STATES...,"all")
    ix=[i for i in eachindex(yp) if (state=="all"||ctx.observations[i].state==state)&&ctx.y[i]!=0]
    push!(results,(manuscript_model=mapping[id],storage_model=id,cohort,state,n=length(ix),nll,bic=score,mrpe_percent=100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix]))))
   end
   for (r,y,s) in zip(ctx.observations,yp,ctx.sigma)
    push!(point_rows,(manuscript_model=mapping[id],storage_model=id,cohort,condition=r.condition,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,sigma=s,predicted=y))
   end
   println("Verified ",cohort," Model ",mapping[id]," NLL=",nll," BIC=",score);flush(stdout)
  end
 end
 ctx=context(cfg=cfg,conditions=(4,))
 for id in [0,5]
  model=legacy_model(id);p=read_parameters(joinpath(fit,"model$id/selected_parameters.csv"),model);yp=prediction(ctx,model,p)
  saved=filter(r->r.publication_model==id,collect(CSV.File(joinpath(fit,"figure_data/test_predictions.csv"))))
  @assert length(saved)==length(yp)
  for (i,r) in enumerate(ctx.observations)
   sr=only(filter(s->s.condition==r.condition && s.day==r.day && s.state==r.state,saved))
   @assert isapprox(sr.observed,r.mean;rtol=1e-10)
   @assert isapprox(sr.predicted,yp[i];rtol=1e-7,atol=1e-6)
  end
  println("Verified 500 mL Model ",mapping[id]," predictions: ",length(yp));flush(stdout)
 end
 CSV.write(joinpath(@__DIR__,"verification/primary_scores.csv"),results)
 CSV.write(joinpath(@__DIR__,"verification/primary_predictions.csv"),point_rows)
 println("PRIMARY VERIFICATION PASSED")
end
verify()
