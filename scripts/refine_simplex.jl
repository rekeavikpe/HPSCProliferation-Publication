using CSV, TOML, LinearAlgebra
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration
length(ARGS)==2||error("Usage: refine_simplex.jl FIT_DIRECTORY OUTPUT_DIRECTORY")
fitdir=abspath(ARGS[1]);out=abspath(ARGS[2]);ctx=checked_context(fitdir)
prepare_output(out,fitdir;settings=Dict("analysis"=>"small-simplex restarts for Models 1 and 2","simplex_sizes"=>[1e-5,1e-6],"lm_iterations"=>0,"polish_iterations"=>12000,"stopping_tolerance"=>1e-8))
BLAS.set_num_threads(1)
jobs=[(id=id,run=j,scale=scale) for id in (1,2) for (j,scale) in enumerate((1e-5,1e-6))]
Threads.@threads :greedy for j in jobs
    m=legacy_model(j.id);p=read_parameters(joinpath(fitdir,"model$(j.id)/selected_parameters.csv"),m)
    fit_start(ctx,m,encode(ctx,m,p),joinpath(out,"model$(j.id)/run$(j.run)");lm_iterations=0,polish_iterations=12000,tight=true,simplex_scale=j.scale)
end
rows=[save_choice(ctx,legacy_model(id),joinpath(out,"model$id")) for id in (1,2)]
CSV.write(joinpath(out,"summary.csv"),rows)
metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete"
open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
write(joinpath(out,"COMPLETE"),"complete; small-simplex convergence is a local numerical diagnostic, not proof of a global optimum\n")
println(rows)
