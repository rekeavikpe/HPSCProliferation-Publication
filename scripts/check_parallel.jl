using CSV, Test, LinearAlgebra
include("../src/Workflow.jl")
using .Workflow
using .Workflow.Calibration: encode
BLAS.set_num_threads(1)
ctx=context()
out=joinpath(ROOT,"results/parallel_check_$(Threads.nthreads())")
mkpath(out)
tasks=[]
for id in 0:6
    m=legacy_model(id)
    p=read_parameters(joinpath(ROOT,"reference/scale_transfer/legacy_model$m/selected_parameters.csv"),m)
    for run in 1:6
        push!(tasks,(model=m,z=encode(ctx,m,p),dir=joinpath(out,"model$id/run$run")))
    end
end
elapsed=@elapsed Threads.@threads :dynamic for i in eachindex(tasks)
    task=tasks[i]
    fit_start(ctx,task.model,task.z,task.dir;lm_iterations=1,polish_iterations=1)
end
CSV.write(joinpath(out,"scores.csv"),[(task=i,nll=last(CSV.File(joinpath(t.dir,"stages.csv"))).nll) for (i,t) in enumerate(tasks)])
println("42 independent optimizer smoke jobs completed with $(Threads.nthreads()) threads in $elapsed seconds. Not full fits.")
