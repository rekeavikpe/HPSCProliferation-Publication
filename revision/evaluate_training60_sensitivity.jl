using CSV, TOML, Statistics
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("../ablations/RoutedExposure.jl")
using .RoutedExposure
length(ARGS)==3 || error("Usage: evaluate_training60_sensitivity.jl FIT ABLATION_RESULTS ORIGINAL_SPLIT")
fit,out,splitfit=abspath.(ARGS);cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
train=context(;cfg)
held=context(;ids=cfg["calibration_cohort"]["selection_ids"],cfg,fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
key(r)=(r.condition,r.day,r.state)
scales=Dict(key(r)=>s for (r,s) in zip(train.observations,train.sigma))
held=merge(held,(sigma=[scales[key(r)] for r in held.observations],))
shared=TOML.parsefile(joinpath(out,"manifest.toml"))["settings"]["shared_response_parameter"]
rank=NamedTuple[];metrics=NamedTuple[];runs=NamedTuple[]
for row in CSV.File(joinpath(out,"ranking.csv"))
    scenario=String(row.scenario);model=scenario=="baseline" ? 0 : 6
    ctx=scenario in ("baseline","full") ? held : routed_context(held,scenario;shared)
    p=read_parameters(joinpath(out,scenario,"selected_parameters.csv"),model)
    yp=prediction(ctx,model,p);nll=gaussian_nll(ctx.y,yp,ctx.sigma)
    push!(rank,(scenario,n=length(yp),k=row.k,nll,bic=bic(nll,row.k,length(yp)),stopping_met=row.stopping_met))
    for state in (STATES...,"all")
        ix=findall(i->(state=="all" || ctx.observations[i].state==state)&&ctx.y[i]!=0,eachindex(yp))
        push!(metrics,(scenario,state,mrpe_percent=100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix]))))
    end
    if scenario in ("baseline","full")
        for r in CSV.File(joinpath(fit,"selection_optimizer_bic.csv"))
            r.storage_model==(scenario=="baseline" ? 0 : 5) || continue
            push!(runs,(scenario,run=r.run,bic=r.bic,stopping_met=r.stopping_met))
        end
    else
        for name in filter(n->occursin(r"^run\d+$",n),readdir(joinpath(out,scenario)))
            dir=joinpath(out,scenario,name);isfile(joinpath(dir,"COMPLETE"))||error("Incomplete ablation run")
            rp=read_parameters(joinpath(dir,"selected_parameters.csv"),model)
            rr=only(CSV.File(joinpath(dir,"summary.csv")))
            v=gaussian_nll(ctx.y,prediction(ctx,model,rp),ctx.sigma)
            push!(runs,(scenario,run=parse(Int,name[4:end]),bic=bic(v,row.k,length(ctx.y)),stopping_met=rr.stopping_met))
        end
    end
end
CSV.write(joinpath(out,"selection_ranking.csv"),sort(rank;by=r->r.bic))
CSV.write(joinpath(out,"selection_state_errors.csv"),metrics)
CSV.write(joinpath(out,"selection_optimizer_bic.csv"),runs)
println("Eight formulations evaluated on the unchanged 40% selection cultures")
