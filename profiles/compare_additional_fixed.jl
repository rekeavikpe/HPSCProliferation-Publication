using CSV, TOML, Statistics
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)>=2||error("Usage: compare_additional_fixed.jl OUTPUT CASE_DIRECTORY...")
out=abspath(ARGS[1]);mkpath(out);model=legacy_model(5)
checks=NamedTuple[];fits=NamedTuple[];errors=NamedTuple[]
for folder in abspath.(ARGS[2:end])
    isfile(joinpath(folder,"COMPLETE"))||error("Incomplete suite: $folder")
    meta=TOML.parsefile(joinpath(folder,"settings.toml"));cfg=TOML.parsefile(joinpath(folder,"calibration.toml"))
    p=read_parameters(joinpath(folder,"mle_parameters.csv"),model);ctx=context(cfg=cfg)
    fit=only(CSV.File(joinpath(folder,"mle_summary.csv")))
    nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
    isapprox(nll,fit.nll;atol=1e-6,rtol=0)||error("MLE failed verification")
    for (name,value) in cfg["fixed"]
        getproperty(p,Symbol(name))==value||error("Fixed value changed")
    end
    label=basename(folder)
    for name in meta["parameters"]
        points=collect(CSV.File(joinpath(folder,name,"profile.csv")))
        ci=only(CSV.File(joinpath(folder,name,"interval.csv")))
        lower=any(r.parameter_value<ci.mle&&r.profile_loglik<=-1.920729410347062&&r.stopping_met===true for r in points)
        upper=any(r.parameter_value>ci.mle&&r.profile_loglik<=-1.920729410347062&&r.stopping_met===true for r in points)
        nuisance=collect(CSV.File(joinpath(folder,name,"nuisance_parameters.csv")))
        touching=String[]
        for r in nuisance
            r.parameter==name&&continue
            points[r.point].profile_loglik>=-1.920729410347062||continue
            lo,hi=cfg["bounds"][r.parameter]
            z=log(r.value/lo)/log(hi/lo)
            (z<1e-4||z>1-1e-4)&&push!(touching,r.parameter)
        end
        push!(checks,(case=label,parameter=name,mle=ci.mle,lower_crossed=lower,upper_crossed=upper,package_lower=ci.package_lower,package_upper=ci.package_upper,unconfirmed_points=ci.unconfirmed_points,better_minimum=ci.better_minimum,supported_nuisance_boundaries=join(sort(unique(touching)),";")))
    end
    push!(fits,(case=label,nll,free_parameters=length(free_names(ctx,model)),stopping_met=fit.stopping_met,additional_fixed_parameter=get(meta,"additional_fixed_parameter",""),additional_fixed_value=get(cfg["fixed"],get(meta,"additional_fixed_parameter",""),NaN)))
    for (cohort,conditions) in (("100mL calibration",(1,2,3,5,6)),("500mL retrospective",(4,)))
        c=context(cfg=cfg,conditions=conditions);y=prediction(c,model,p)
        for state in STATES
            ix=findall(r->r.state==state&&r.mean!=0,c.observations)
            push!(errors,(case=label,cohort,state,mrpe=100mean(abs.((y[ix].-c.y[ix])./c.y[ix]))))
        end
    end
end
CSV.write(joinpath(out,"profile_checks.csv"),checks)
CSV.write(joinpath(out,"fit_comparison.csv"),fits)
CSV.write(joinpath(out,"prediction_errors.csv"),errors)
println("Verified conditional fits, crossings, nuisance-bound diagnostics and predictions: $out")
