using CSV, TOML, SHA, LinearAlgebra
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("HybridRefinement.jl")
using .HybridRefinement
length(ARGS)>=2 || error("Usage: retest_refinement.jl OUT POPULATION_RESULTS...")
out=abspath(ARGS[1]);sources=abspath.(ARGS[2:end]);ispath(out)&&error("Choose new output")
cfg=TOML.parsefile(joinpath(first(sources),"calibration.toml"))
all(s->TOML.parsefile(joinpath(s,"calibration.toml"))==cfg,sources)||error("Configs differ")
iterations=parse(Int,get(ENV,"HPSC_HYBRID_ITERATIONS","3000"))
models=parse.(Int,split(get(ENV,"HPSC_HYBRID_MODELS","0,1,2,3,4,5,6"),','))
ctx=context(cfg=cfg);BLAS.set_num_threads(1);mkpath(out)
cp(joinpath(first(sources),"calibration.toml"),joinpath(out,"calibration.toml"))
CSV.write(joinpath(out,"observations.csv"),ctx.observations)
jobs=NamedTuple[]
for source in sources, folder in filter(n->startswith(n,"model"),readdir(source))
    id=parse(Int,replace(folder,"model"=>""));id in models || continue
    for run in filter(n->startswith(n,"run"),readdir(joinpath(source,folder)))
        push!(jobs,(id=id,run=run,source=joinpath(source,folder,run)))
    end
end
open(joinpath(out,"run.toml"),"w") do io
    TOML.print(io,Dict("source_runs"=>sources,"method"=>"original direct NelderMead after saved DE candidates",
        "iterations"=>iterations,"models"=>models,"refiner_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"HybridRefinement.jl"))))))
end
Threads.@threads :greedy for job in jobs
    m=legacy_model(job.id);dir=joinpath(out,"model$(job.id)",job.run)
    p=read_parameters(joinpath(job.source,"population_parameters.csv"),m)
    stage=only(refine_original(ctx,m,encode(ctx,m,p),dir;iterations))
    popnll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
    old=only(CSV.File(joinpath(job.source,"summary.csv")))
    CSV.write(joinpath(dir,"summary.csv"),[(model=job.id,run=job.run,bic=bic(stage.nll,length(p),length(ctx.y)),
        stopping_met=stage.stopping_met,population_bic=bic(popnll,length(p),length(ctx.y)),previous_hybrid_bic=old.bic)])
end
rows=[NamedTuple(only(CSV.File(joinpath(out,"model$(j.id)",j.run,"summary.csv")))) for j in jobs]
CSV.write(joinpath(out,"runs.csv"),rows)
write(joinpath(out,"COMPLETE"),"all refinements completed\n")
println(rows)
