using CSV, TOML, SHA
length(ARGS)==3 || error("Usage: export_split_combined.jl SPLIT_RESULTS ALL_100ML_FIGURE_DATA OUTPUT_DATA")
splitdir,full,out=abspath.(ARGS)
isfile(joinpath(splitdir,"COMPLETE")) || error("Split comparison is incomplete")
protocol=TOML.parsefile(joinpath(splitdir,"protocol.toml"))
sort(Int.(protocol["storage_models"]))==collect(0:6) || error("Expected all seven fitted models")
lock=only(CSV.File(joinpath(splitdir,"selection_lock.csv")))
meta=TOML.parsefile(joinpath(full,"source.toml"))
Int(lock.storage_model)==Int(meta["selected_model"]) || error("The split selected publication Model $(lock.publication_model), but the supplied 500-mL predictions use a different model. Supply the matching all-100-mL fitted model's prediction export before combining.")
fullfit=dirname(full)
a=TOML.parsefile(joinpath(splitdir,"calibration.toml"))
b=TOML.parsefile(joinpath(fullfit,"calibration.toml"))
for key in ("bounds","fixed","likelihood","exposure_percentile","initial_lactate_mM","initial_aggregate_um","initial_glucose","fresh_glucose")
    a[key]==b[key] || error("Split and all-100-mL settings differ: $key")
end
hashes=collect(CSV.File(joinpath(full,"parameter_hashes.csv")))
for id in unique([0,Int(lock.storage_model)])
    row=only(filter(r->r.model==id,hashes))
    bytes2hex(sha256(read(joinpath(fullfit,"model$id","selected_parameters.csv"))))==row.sha256 || error("All-data parameters differ from the prediction export")
end
ranking=collect(CSV.File(joinpath(splitdir,"training_ranking.csv")))
runs=NamedTuple[]
for id in 0:6
    fit=only(filter(r->r.model==id,ranking))
    for run in 1:protocol["starts"]
        dir=joinpath(splitdir,"model$id","run$run")
        isfile(joinpath(dir,"COMPLETE")) || error("Missing completed checkpoint: $dir")
        stages=collect(CSV.File(joinpath(dir,"stages.csv")))
        best=argmin(r->r.nll,stages)
        push!(runs,(model=id,run="run$run",bic=2best.nll+fit.k*log(fit.n),stopping_met=best.stopping_met,stage=best.stage))
    end
end
for file in ("test_trajectories.csv","test_predictions.csv")
    read(joinpath(full,file))==read(joinpath(fullfit,"scale_transfer",file)) || error("Prediction export differs from saved scale-transfer output: $file")
end
mkpath(out)
CSV.write(joinpath(out,"optimizer_bic.csv"),runs)
cp(joinpath(splitdir,"training_ranking.csv"),joinpath(out,"training_ranking.csv");force=true)
for file in ("test_trajectories.csv","test_predictions.csv","parameter_hashes.csv")
    cp(joinpath(full,file),joinpath(out,file);force=true)
end
for file in ("selection_lock.csv","selection_ranking.csv","split_counts.csv")
    cp(joinpath(splitdir,file),joinpath(out,file);force=true)
end
newmeta=Dict("selected_model"=>Int(lock.storage_model),"initial_starts_per_model"=>protocol["starts"],
    "source_run"=>basename(splitdir),"prediction_source_run"=>meta["source_run"],
    "bic_scope"=>"Calibration BIC from the 60% culture partition. The model used in the prediction panels was selected using predictive NLL on the held-out 40%, not training BIC. The 500-mL panels reuse its verified, previously completed all-100-mL calibration under matching settings; they do not use the 60%-fit parameters.")
open(io->TOML.print(io,newmeta),joinpath(out,"source.toml"),"w")
println("Combined inputs prepared. Left: split calibration BIC. Right: publication Model $(lock.publication_model), existing matching all-100-mL fit.")
