using CSV, TOML, Statistics
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)>=2||error("Usage: lactate_limit.jl FIT_DIRECTORY OUTPUT_DIRECTORY [STARTS=4]")
fitdir=abspath(ARGS[1]);out=abspath(ARGS[2]);n=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 4
ctx=checked_context(fitdir);m=legacy_model(5)
prepare_output(out,fitdir;settings=Dict("analysis"=>"Model 5 exact inactive lactate inhibition sensitivity","starts"=>n,"seed"=>20261015,"lm_iterations"=>600,"polish_iterations"=>12000))
p=read_parameters(joinpath(fitdir,"model5/selected_parameters.csv"),m)
cfg=deepcopy(ctx.config);cfg["lactate_inhibition"]=false;cfg["fixed"]["K_Lac"]=1.0
open(joinpath(out,"calibration.toml"),"w") do io;TOML.print(io,cfg);end
inactive=context(cfg=cfg)
jobs=[(ctx=inactive,model=m,z=z,dir=joinpath(out,"model5/run$j")) for (j,z) in enumerate(local_starts(inactive,m,p;n=n,seed=20261015))]
run_jobs(jobs)
summary=save_choice(inactive,m,joinpath(out,"model5");estimated_parameters=10)
old=only(CSV.File(joinpath(fitdir,"model5/summary.csv")))
CSV.write(joinpath(out,"comparison.csv"),[(formulation="bounded lactate inhibition",nll=old.nll,bic=old.bic,k=11,stopping_met=old.stopping_met),(formulation="exact inactive inhibition",nll=summary.nll,bic=summary.bic,k=10,stopping_met=summary.stopping_met)])
q=read_parameters(joinpath(out,"model5/selected_parameters.csv"),m)
records=NamedTuple[]
for conditions in ((1,2,3,5,6),(4,))
    original=context(conditions=conditions,cfg=ctx.config);modified=context(conditions=conditions,cfg=cfg)
    a=prediction(original,m,p);b=prediction(modified,m,q)
    for state in STATES
        ix=findall(r->r.state==state,original.observations);y=original.y[ix];nz=filter(i->original.y[i]!=0,ix)
        push!(records,(cohort=length(conditions)==1 ? "500mL retrospective" : "100mL calibration",state=state,original_mrpe=100mean(abs.((a[nz].-original.y[nz])./original.y[nz])),inactive_mrpe=100mean(abs.((b[nz].-original.y[nz])./original.y[nz])),max_difference_over_sigma=maximum(abs.((b[ix].-a[ix])./original.sigma[ix]))))
    end
end
CSV.write(joinpath(out,"prediction_comparison.csv"),records)
metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete"
open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
write(joinpath(out,"COMPLETE"),"sensitivity completed; original formulation retained\n")
println(summary)
