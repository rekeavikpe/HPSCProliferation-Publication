include(joinpath(@__DIR__,"..","revision","October20.jl"))
function verify_additional()
 cache=joinpath(@__DIR__,"cache")
 fit=joinpath(cache,"training60/reference_fit")
 prep=SplitSelection.prepare(joinpath(fit,"calibration.toml"))
 scores=collect(CSV.File(joinpath(cache,"october20/p75/selection_scores.csv")))
 held=exposure_context(prep.held,"lag_p75")
 for m in [0,1,3,4,5,6,7]
  p=read_parameters(joinpath(cache,"october20/p75/model$m/selected_parameters.csv"),m)
  r=score_row(held,m,p,"Model $m","selection40")
  old=only(filter(x->x.tag==r.tag,scores))
  @assert isapprox(r.nll,old.nll;atol=1e-6,rtol=0)
  @assert isapprox(r.bic,old.bic;atol=2e-6,rtol=0)
  println("Verified P75 ",r.tag," BIC=",r.bic);flush(stdout)
 end
 result=NamedTuple[]
 for line in ["585A1","SCTi003-A"]
  for (label,folder) in [("Recalibration 1",joinpath(cache,"october20/cells",line,"round1")),("Recalibration 2",joinpath(cache,"ipsc",line,"fixed_K_d_Y_GLac"))]
   ctx=cell_context(line,TOML.parsefile(joinpath(folder,"calibration.toml")))
   p=read_parameters(joinpath(folder,"reported_parameters.csv"),6)
   yp=prediction(ctx,6,p);nll=gaussian_nll(ctx.y,yp,ctx.sigma)
   saved=only(CSV.File(joinpath(folder,"summary.csv")))
   @assert isapprox(nll,saved.nll;atol=1e-6,rtol=0)
   err=100mean(abs.((yp.-ctx.y)./ctx.y))
   push!(result,(line,strategy=label,nll,mrpe_percent=err))
   println("Verified ",line," ",label," MRPE=",err);flush(stdout)
  end
 end
 CSV.write(joinpath(@__DIR__,"verification/ipsc_refit_scores.csv"),result)
 println("P75 AND iPSC FIT VERIFICATION PASSED")
end
verify_additional()
