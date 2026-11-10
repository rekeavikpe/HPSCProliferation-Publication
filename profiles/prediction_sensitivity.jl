using CSV, TOML, Statistics
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)==2||error("Usage: prediction_sensitivity.jl ACTIVE_FIT INACTIVE_PROFILES")
active,out=abspath.(ARGS)
isfile(joinpath(out,"COMPLETE"))||error("Inactive profiles are incomplete")
meta=TOML.parsefile(joinpath(out,"settings.toml"))
get(meta,"smoke_test",true)&&error("Smoke profiles cannot support this analysis")
cfg=TOML.parsefile(joinpath(out,"calibration.toml"));acfg=TOML.parsefile(joinpath(active,"calibration.toml"))
!cfg["lactate_inhibition"]&&get(acfg,"lactate_inhibition",true)||error("Expected inactive and active formulations")
cfg["fixed"]["Y_GX"]==acfg["fixed"]["Y_GX"]||error("Yield fixation differs")
TOML.parsefile(joinpath(active,"settings.toml"))["source_fit"]==meta["source_fit"]||error("Source calibrations differ")
model=legacy_model(5);ctx=context(cfg=cfg)
p=read_parameters(joinpath(out,"mle_parameters.csv"),model);ap=read_parameters(joinpath(active,"mle_parameters.csv"),model)
mle=only(CSV.File(joinpath(out,"mle_summary.csv")))
mle.stopping_met||error("MLE stopping criterion was not met")
actual=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
isapprox(actual,mle.nll;atol=1e-6,rtol=0)||error("MLE no longer reproduces")
ensemble=[p];audit=NamedTuple[];threshold=-1.920729410347062
for name in meta["parameters"]
    points=collect(CSV.File(joinpath(out,name,"profile.csv")))
    nuisance=collect(CSV.File(joinpath(out,name,"nuisance_parameters.csv")))
    maximum(r.profile_loglik for r in points)>1e-5&&error("Profile found a better MLE; refit before constructing envelope: $name")
    for (j,r) in enumerate(points)
        accepted=r.stopping_met===true && r.profile_loglik>=threshold
        push!(audit,(parameter=name,point=j,profile_loglik=r.profile_loglik,stopping_met=r.stopping_met,accepted))
        accepted||continue
        values=Dict(n=>getproperty(p,n) for n in keys(p))
        for v in nuisance
            v.point==j&&(values[Symbol(v.parameter)]=v.value)
        end
        q=parameters(model,[values[n] for n in parameter_names(model)])
        ll=actual-gaussian_nll(ctx.y,prediction(ctx,model,q),ctx.sigma)
        isapprox(ll,r.profile_loglik;atol=1e-5,rtol=0)||error("Saved profile vector failed reevaluation")
        push!(ensemble,q)
    end
end
CSV.write(joinpath(out,"envelope_sample_audit.csv"),audit)
function trajectory(ctx,model,p,c)
    q=ctx.conditions[c];times=sort(unique(vcat(collect(0.:3.:last(q.times)),q.times,[f.time_h for f in q.feeds])))
    result=simulate(model,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,lactate_inhibition=get(ctx.config,"lactate_inhibition",true),reltol=1e-12,abstol=1e-13)
    entries=[(time_h=t,phase=0,values=result.predictions[:,i]) for (i,t) in enumerate(times)]
    append!(entries,[(time_h=e.time_h,phase=1,values=e.after) for e in result.events])
    sort!(entries;by=r->(r.time_h,r.phase))
end
trajectories=NamedTuple[];errors=NamedTuple[];observations=NamedTuple[]
for (cohort,conditions) in (("100mL calibration",(1,2,3,5,6)),("500mL retrospective",(4,)))
    ic=context(cfg=cfg,conditions=conditions);ac=context(cfg=acfg,conditions=conditions)
    append!(observations,[(;cohort,r...) for r in ic.observations])
    pred=prediction(ic,model,p);original=prediction(ac,model,ap)
    for state in STATES
        keep=findall(r->r.state==state,ic.observations);nonzero=filter(i->ic.y[i]!=0,keep)
        push!(errors,(cohort,state,active_mrpe=100mean(abs.((original[nonzero].-ic.y[nonzero])./ic.y[nonzero])),inactive_mrpe=100mean(abs.((pred[nonzero].-ic.y[nonzero])./ic.y[nonzero]))))
    end
    for c in conditions
        center=trajectory(ic,model,p,c);before=trajectory(ac,model,ap,c)
        samples=Vector{Matrix{Float64}}(undef,length(ensemble))
        Threads.@threads :greedy for j in eachindex(ensemble)
            samples[j]=hcat([r.values for r in trajectory(ic,model,ensemble[j],c)]...)
        end
        lower=reduce((a,b)->min.(a,b),samples);upper=reduce((a,b)->max.(a,b),samples)
        for (i,r) in enumerate(center),(k,state) in enumerate(STATES)
            push!(trajectories,(cohort,condition=c,state,time_h=r.time_h,phase=r.phase,inactive_mle=r.values[k],active_mle=before[i].values[k],profile_min=lower[k,i],profile_max=upper[k,i]))
        end
    end
end
CSV.write(joinpath(out,"sensitivity_observations.csv"),observations)
CSV.write(joinpath(out,"prediction_sensitivity.csv"),trajectories)
CSV.write(joinpath(out,"prediction_comparison.csv"),errors)
CSV.write(joinpath(out,"fit_comparison.csv"),[(formulation=n,nll=only(CSV.File(joinpath(d,"mle_summary.csv"))).nll,free_parameters=length(free_names(context(cfg=c),model))) for (n,d,c) in (("active inhibition; yield fixed",active,acfg),("inactive inhibition; yield fixed",out,cfg))])
write(joinpath(out,"prediction_sensitivity_notes.txt"),"Conditional profile-path sensitivity envelope from $(length(ensemble)) parameter vectors, including the MLE. Includes converged saved profile points within the nominal one-parameter likelihood-ratio threshold. Not a simultaneous 95% prediction interval: it does not cover the complete joint parameter region, uncertainty in fixed Y_GX, observation noise, or model discrepancy. No 500-mL observations were used to select parameters or envelope samples. The 500-mL evaluation is retrospective. Active/inactive fit comparison reports NLL and free dimensions only: Y_GX was estimated using these data, so counting it as known for ordinary BIC would be inappropriate.\n")
println("Prediction sensitivity saved for $(length(ensemble)) supported vectors: $out")
