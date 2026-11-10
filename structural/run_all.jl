using CSV, TOML, SHA
root=normpath(joinpath(@__DIR__,".."))
out=isempty(ARGS) ? joinpath(root,"results/structural") : abspath(ARGS[1])
ispath(out)&&error("Choose a new output directory: $out")
mkpath(out)
julia=Base.julia_cmd()
# Separate processes avoid the symbolic package's shared logging state.
for id in 0:6
    run(`$julia --startup-file=no --project=$(@__DIR__) $(joinpath(@__DIR__,"run.jl")) $id pooled 0 global $out`)
    run(`$julia --startup-file=no --project=$(@__DIR__) $(joinpath(@__DIR__,"initial_state_rank.jl")) $id $out`)
end
summary=NamedTuple[]
for id in 0:6
    rows=collect(CSV.File(joinpath(out,"pooled_segment0_global/model$id/parameters.csv")))
    rank=TOML.parsefile(joinpath(out,"nominal_initial_state/model$id/result.toml"))
    push!(summary,(model=id,parameters=length(rows),global_generic_segment=count(r->r.classification=="globally",rows),local_only_generic_segment=join([r.parameter for r in rows if r.classification=="locally"],";"),nominal_initial_rank=rank["rank"]))
end
CSV.write(joinpath(out,"summary.csv"),summary)
files=vcat([joinpath(@__DIR__,f) for f in readdir(@__DIR__) if endswith(f,".jl")||endswith(f,".toml")],[joinpath(root,"data/processed/lagrangian_exposure.csv"),joinpath(root,"data/processed/measurements.csv")])
CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,root),sha256=bytes2hex(sha256(read(f)))) for f in files])
write(joinpath(out,"COMPLETE"),"complete; read docs/structural_identifiability.md for scope\n")
