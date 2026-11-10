using CSV, TOML
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration
length(ARGS)>=2||error("Usage: refine_selected.jl FIT_DIRECTORY OUTPUT_DIRECTORY [STARTS=4]")
fitdir=abspath(ARGS[1]);out=abspath(ARGS[2]);n=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 4
n>=2||error("At least two starts required")
ctx=checked_context(fitdir)
prepare_output(out,fitdir;settings=Dict("analysis"=>"targeted models 1 and 2 refinement","starts"=>n,"lm_iterations"=>600,"polish_iterations"=>12000,"seed"=>20261005))
jobs=NamedTuple[]
for id in (1,2)
    m=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),m)
    folder=joinpath(out,"model$id");original=joinpath(folder,"run0");mkpath(original)
    old=only(CSV.File(joinpath(fitdir,"model$id/summary.csv")))
    save_parameters(joinpath(original,"original_parameters.csv"),p)
    CSV.write(joinpath(original,"stages.csv"),[(stage="original",nll=old.nll,stopping_met=old.stopping_met)])
    write(joinpath(original,"COMPLETE"),"saved original fit, not a new optimization\n")
    for (j,z) in enumerate(local_starts(ctx,m,p;n=n,seed=20261005+id))
        push!(jobs,(ctx=ctx,model=m,z=z,dir=joinpath(folder,"run$j")))
    end
end
run_jobs(jobs)
rows=[save_choice(ctx,legacy_model(id),joinpath(out,"model$id")) for id in (1,2)]
CSV.write(joinpath(out,"summary.csv"),rows)
metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete"
open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
write(joinpath(out,"COMPLETE"),"completed; source calibration and selection lock unchanged\n")
println(rows)
