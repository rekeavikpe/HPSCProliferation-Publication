using CSV, TOML, Statistics, SHA
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
function exposure_context(base,scenario)
    inputs=collect(CSV.File(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    conditions=copy(base.conditions)
    for (c,q) in base.conditions
        info=first(filter(r->r.condition==c,series))
        rows=sort(filter(r->r.scenario==scenario&&r.volume_ml==info.volume_ml&&r.rpm==info.rpm,inputs);by=r->r.size_um)
        lookup=if startswith(scenario,"lag_")
            @assert [r.size_um for r in rows]==[200,400,600,800,1000]
            linear_exposure([r.size_um for r in rows],[r.shear_pa for r in rows],[r.edr_m2_s3 for r in rows])
        else
            r=only(rows);let tau=r.shear_pa,edr=r.edr_m2_s3;size->(tau,edr);end
        end
        conditions[c]=merge(q,(;lookup))
    end
    merge(base,(;conditions))
end

length(ARGS)==3 || error("Usage: evaluate_training60_exposure.jl FIT SAVED_EXPOSURE_SELECTION OUTPUT")
fit,source,out=abspath.(ARGS);cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
isfile(joinpath(source,"COMPLETE")) || error("Exposure-selection source is incomplete")
saved=collect(CSV.File(joinpath(source,"split.csv")))
for (role,ids) in (("calibration",cfg["calibration_cohort"]["series_ids"]),("selection",cfg["calibration_cohort"]["selection_ids"]))
    Set(String(r.series_id) for r in saved if r.role==role)==Set(ids) || error("Exposure split mismatch")
end
base=context(;cfg);ext=context(;cfg,conditions=(4,))
held=context(;cfg,ids=cfg["calibration_cohort"]["selection_ids"],fallback_glucose=Dict(c=>q.glucose for (c,q) in base.conditions))
key(r)=(r.condition,r.day,r.state);scales=Dict(key(r)=>s for (r,s) in zip(base.observations,base.sigma))
held=merge(held,(sigma=[scales[key(r)] for r in held.observations],))
mkpath(out);scores=NamedTuple[];metrics=NamedTuple[];hashes=NamedTuple[]
for row in CSV.File(joinpath(source,"ranking.csv"))
    scenario=String(row.scenario);file=joinpath(source,scenario,"selected_parameters.csv")
    p=read_parameters(file,6);ctx=exposure_context(base,scenario)
    actual=gaussian_nll(ctx.y,prediction(ctx,6,p),ctx.sigma)
    isapprox(actual,row.training_nll;atol=1e-6,rtol=0) || error("Exposure training likelihood mismatch: $scenario")
    sel=exposure_context(held,scenario);nll=gaussian_nll(sel.y,prediction(sel,6,p),sel.sigma)
    isapprox(nll,row.selection_nll;atol=1e-6,rtol=0)||error("Exposure selection likelihood mismatch")
    push!(scores,(scenario,n=length(sel.y),k=row.k,nll,bic=bic(nll,row.k,length(sel.y))))
    push!(hashes,(scenario,source=file,sha256=bytes2hex(sha256(read(file)))))
    for (cohort,b) in (("selection40",held),("500 mL prediction",ext))
        c=exposure_context(b,scenario);yp=prediction(c,6,p)
        for state in (STATES...,"all")
            ix=findall(i->(state=="all"||c.observations[i].state==state)&&c.y[i]!=0,eachindex(yp))
            ni=filter(i->c.observations[i].n_vessels>=2,ix)
            push!(metrics,(scenario,cohort,state,mrpe_percent=100mean(abs.((yp[ix].-c.y[ix])./c.y[ix])),experimental_noise_percent=isempty(ni) ? NaN : 100mean([c.observations[i].sd/abs(c.y[i]) for i in ni])))
        end
    end
end
CSV.write(joinpath(out,"selection_scores.csv"),sort(scores;by=r->r.bic))
CSV.write(joinpath(out,"state_errors.csv"),metrics)
CSV.write(joinpath(out,"parameter_sources.csv"),hashes)
write(joinpath(out,"COMPLETE"),"Saved 60% exposure fits verified and evaluated without refitting\n")
