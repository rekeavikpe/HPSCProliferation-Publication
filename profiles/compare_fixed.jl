using CSV, TOML, Statistics
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)==1||error("Usage: compare_fixed.jl CONDITIONAL_PROFILE_DIRECTORY")
out=abspath(ARGS[1]);meta=TOML.parsefile(joinpath(out,"settings.toml"));fitdir=meta["source_fit"];id=meta["model"];model=legacy_model(id)
original=checked_context(fitdir);cfg=TOML.parsefile(joinpath(out,"calibration.toml"))
old=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
new=read_parameters(joinpath(out,"mle_parameters.csv"),model)
for (name,value) in meta["fixed_parameters"]
    @assert getproperty(new,Symbol(name))==value
end
@assert length(free_names(context(cfg=cfg),model))==9
rows=NamedTuple[]
for (cohort,conditions) in (("100mL calibration",(1,2,3,5,6)),("500mL retrospective",(4,)))
    before=context(cfg=original.config,conditions=conditions);after=context(cfg=cfg,conditions=conditions)
    a=prediction(before,model,old);b=prediction(after,model,new)
    for state in STATES
        keep=findall(r->r.state==state,after.observations)
        push!(rows,(cohort,state,original_mrpe=100mean(abs.((a[keep].-after.y[keep])./after.y[keep])),conditional_mrpe=100mean(abs.((b[keep].-after.y[keep])./after.y[keep])),max_difference_over_sigma=maximum(abs.((b[keep].-a[keep])./after.sigma[keep]))))
    end
end
CSV.write(joinpath(out,"prediction_comparison.csv"),rows)
println("Verified both fixed values, nine free parameters, and exported prediction comparisons.")
