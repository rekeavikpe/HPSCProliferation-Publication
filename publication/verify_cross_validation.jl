using CSV, TOML, Statistics
include(joinpath(@__DIR__,"..","population","CrossValidation.jl"))
using .CrossValidation, .CrossValidation.SplitSelection.Workflow, .CrossValidation.SplitSelection.Workflow.Calibration, .CrossValidation.SplitSelection.Workflow.Calibration.Proliferation
function verify_cv()
 root=joinpath(@__DIR__,"cache/cross_validation"); results=NamedTuple[]
 for plan in filter(p->p.experiment=="loco",plans())
  folder=joinpath(root,plan.name);prep=prepare_fold(joinpath(folder,"calibration.toml"),plan)
  model=legacy_model(5);p=read_parameters(joinpath(folder,"model5/selected_parameters.csv"),model)
  nll=gaussian_nll(prep.train.y,prediction(prep.train,model,p),prep.train.sigma)
  @assert isapprox(nll,only(CSV.File(joinpath(folder,"model5/summary.csv"))).nll;atol=1e-6,rtol=0)
  yp=prediction(prep.held,model,p)
  points=filter(r->r.storage_model==5 && r.cohort=="selection",collect(CSV.File(joinpath(folder,"predictions.csv"))))
  @assert length(points)==length(yp)
  for (i,r) in enumerate(prep.held.observations)
   saved=only(filter(x->x.condition==r.condition && x.day==r.day && x.state==r.state,points))
   @assert isapprox(saved.predicted,yp[i];rtol=1e-7,atol=1e-6)
  end
  for r in metrics(prep.held,yp)
   push!(results,merge((fold=plan.name,manuscript_condition=only(unique(x.condition for x in prep.held.observations)),),r))
  end
  println("Verified ",plan.name," (paper Condition ",only(unique(r.condition for r in prep.held.observations)),")");flush(stdout)
 end
 CSV.write(joinpath(@__DIR__,"verification/loco_metrics.csv"),results)
 println("LOCO VERIFICATION PASSED")
end
verify_cv()
