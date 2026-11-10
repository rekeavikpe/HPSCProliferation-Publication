using CSV, TOML, SHA
include("../src/Workflow.jl")
using .Workflow

length(ARGS)>=2 || error("Usage: finalize.jl OUTPUT_CALIBRATION POPULATION_RESULTS...")
out=abspath(ARGS[1]);inputs=abspath.(ARGS[2:end])
ispath(out) && error("Choose a new output directory")
cfg=TOML.parsefile(joinpath(first(inputs),"calibration.toml"))
obs=read(joinpath(first(inputs),"observations.csv"))
metadata=[TOML.parsefile(joinpath(d,"run.toml")) for d in inputs]
all(m->m["status"]=="complete",metadata) || error("Population runs are incomplete or have failures")
all(d->TOML.parsefile(joinpath(d,"calibration.toml"))==cfg,inputs) || error("Configurations differ")
all(d->read(joinpath(d,"observations.csv"))==obs,inputs) || error("Observations differ")
repeats=first(metadata)["repeats"]
all(m->m["repeats"]==repeats,metadata) || error("Repeat counts differ")
all(m->m["max_evaluations"]==first(metadata)["max_evaluations"],metadata) || error("Search budgets differ")
ids=vcat([Int.(m["models"]) for m in metadata]...)
sort(ids)==collect(0:6) || error("Supply each model 0–6 exactly once")
ctx=context(cfg=cfg)
mkpath(out)
cp(joinpath(first(inputs),"calibration.toml"),joinpath(out,"calibration.toml"))
cp(joinpath(first(inputs),"observations.csv"),joinpath(out,"observations.csv"))
for (i,(dir,meta)) in enumerate(zip(inputs,metadata))
    cp(joinpath(dir,"source_hashes.csv"),joinpath(out,"source_hashes_input$i.csv"))
    for id in meta["models"]
        dest=joinpath(out,"model$id");cp(joinpath(dir,"model$id"),dest)
        for run in 1:repeats
            path=joinpath(dest,"run$run")
            for row in CSV.File(joinpath(path,"local","stages.csv"))
                stage=String(row.stage)
                cp(joinpath(path,"local",stage*"_parameters.csv"),joinpath(path,stage*"_parameters.csv"))
            end
        end
    end
end
rows=[assemble(ctx,legacy_model(id),joinpath(out,"model$id")) for id in 0:6]
sort!(rows;by=r->r.bic)
CSV.write(joinpath(out,"ranking.csv"),rows)
winner=first(rows)
CSV.write(joinpath(out,"selection_lock.csv"),[(model=winner.model,bic=winner.bic,
    parameter_sha256=bytes2hex(sha256(read(joinpath(out,"model$(winner.model)/selected_parameters.csv")))))])
open(joinpath(out,"run.toml"),"w") do io
    TOML.print(io,Dict("status"=>"complete","starts"=>repeats,"source_runs"=>inputs,
        "method"=>"population search followed by local refinement","max_evaluations"=>first(metadata)["max_evaluations"],
        "note"=>"No initial runs removed. Budget exhaustion is not convergence. Inspect stopping flags."))
end
write(joinpath(out,"COMPLETE"),"assembled population comparison; original fits unchanged\n")
println(rows)
