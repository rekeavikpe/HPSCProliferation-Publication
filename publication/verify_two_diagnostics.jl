include(joinpath(@__DIR__,"..","profiles","run_fixed_parameter_diagnostics.jl"))
function review()
 root=joinpath(@__DIR__,"cache/two_diagnostics")
 ctx=context(cfg=TOML.parsefile(joinpath(root,"calibration.toml")))
 model=legacy_model(5)
 summary=only(CSV.File(joinpath(root,"mle_summary.csv")))
 p=read_parameters(joinpath(root,"mle_parameters.csv"),model)
 nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
 @assert isapprox(nll,summary.nll;atol=1e-6,rtol=0)
 @assert isempty(ctx.config["fixed"]) && length(ctx.config["calibration_cohort"]["series_ids"])==19
 println("Reference reproduced: ",nll);flush(stdout)
 details=NamedTuple[]
 for name in ["K_Lac","Y_GX"]
  points=collect(CSV.File(joinpath(root,name,"profile.csv")))
  nuisance=collect(CSV.File(joinpath(root,name,"nuisance_parameters.csv")))
  for (i,row) in enumerate(points)
   entries=filter(r->r.point==i,nuisance);@assert length(entries)==11
   values=Dict(Symbol(r.parameter)=>Float64(r.value) for r in entries)
   @assert isapprox(values[Symbol(name)],row.parameter_value;rtol=1e-12)
   candidate=parameters(model,[values[k] for k in parameter_names(model)])
   actual=summary.nll-gaussian_nll(ctx.y,prediction(ctx,model,candidate),ctx.sigma)
   @assert abs(actual-row.profile_loglik)<1e-5
   push!(details,(parameter=name,point=i,saved=Float64(row.profile_loglik),replayed=actual,error=actual-row.profile_loglik))
  end
  issues=diagnostic_curve_issues([r.profile_loglik for r in points],[r.stopping_met===true for r in points])
  @assert isempty(issues)
  mle=only(CSV.File(joinpath(root,name,"interval.csv"))).mle
  lo,hi=ctx.config["bounds"][name]
  mz=log(mle/lo)/log(hi/lo)
  lim=diagnostic_limits([r.normalized_parameter for r in points],[r.profile_loglik for r in points],mz)
  println(name," points=",length(points)," issues=",issues," lower_crossed=",lim.lower_crossed," upper_crossed=",lim.upper_crossed," physical_limits=",exp(log(lo)+lim.lower*log(hi/lo))," / ",exp(log(lo)+lim.upper*log(hi/lo)));flush(stdout)
 end
 CSV.write(joinpath(@__DIR__,"verification/two_diagnostic_replay.csv"),details)
 println("All ",length(details)," saved likelihoods reproduced. Maximum difference: ",maximum(abs(r.error) for r in details))
end
review()
