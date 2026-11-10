using CSV, TOML, SHA, Statistics
include("../population/SplitSelection.jl")
using .SplitSelection, .SplitSelection.Workflow, .SplitSelection.Workflow.Calibration
using .SplitSelection.Workflow.Calibration.Proliferation
using .SplitSelection.PublicationModels
length(ARGS)==2 || error("Usage: prepare_training60.jl SAVED_SPLIT OUTPUT")
source,out=abspath.(ARGS)
source!=out || error("Preserve source results")
prep=prepare(joinpath(source,"calibration.toml"))
saved=collect(CSV.File(joinpath(source,"split.csv")))
for (role,ids) in (("calibration",prep.split.training),("selection",prep.split.selection))
    Set(String(r.series_id) for r in saved if r.role==role)==Set(ids) || error("Saved culture split differs")
end
length(prep.split.training)==19 && length(prep.split.selection)==12 || error("Unexpected culture counts")
for (name,ctx) in (("calibration",prep.train),("selection",prep.held))
    obs=collect(CSV.File(joinpath(source,name*"_observations.csv")))
    [(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels,r.sigma) for r in obs]==[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels,s) for (r,s) in zip(ctx.observations,ctx.sigma)] || error("Saved $name observations or weights differ")
end
files=[joinpath(source,"calibration.toml"),joinpath(source,"split.csv"),joinpath(source,"protocol.toml")]
for id in 0:6
    append!(files,[joinpath(source,"model$id",f) for f in ("selected_parameters.csv","summary.csv","residuals.csv")])
end
for id in 0:6, (dir,_,names) in walkdir(joinpath(source,"model$id")), name in names
    (name=="stages.csv" || endswith(name,"_parameters.csv")) && push!(files,joinpath(dir,name))
end
append!(files,[joinpath(source,"calibration_observations.csv"),joinpath(source,"selection_observations.csv")])
manifest=Dict("source"=>source,"hashes"=>Dict(relpath(f,source)=>bytes2hex(sha256(read(f))) for f in files),"policy"=>"Estimate on 19 cultures only. Select on 12 cultures. No all-100-mL refit.")
mkpath(out)
lockfile=joinpath(out,"source_manifest.toml")
isfile(lockfile) && TOML.parsefile(lockfile)!=manifest && error("Source changed. Use a new output folder")
open(io->TOML.print(io,manifest),lockfile,"w")
cfg=deepcopy(prep.cfg)
cfg["calibration_cohort"]=Dict("series_ids"=>prep.split.training,"selection_ids"=>prep.split.selection,"source"=>source)
open(io->TOML.print(io,cfg),joinpath(out,"calibration.toml"),"w")
@assert context(;cfg).y==prep.train.y
@assert context(;cfg).sigma==prep.train.sigma
@assert context(;cfg,conditions=(4,)).y==context(;cfg=prep.cfg,conditions=(4,)).y
CSV.write(joinpath(out,"observations.csv"),prep.train.observations)
cp(joinpath(source,"split.csv"),joinpath(out,"split.csv");force=true)
protocol=TOML.parsefile(joinpath(source,"protocol.toml"))
open(io->TOML.print(io,Dict("starts"=>protocol["starts"],"parameter_source"=>"60% training cultures")),joinpath(out,"run.toml"),"w")
ranks=NamedTuple[];scores=NamedTuple[];metrics=NamedTuple[];runs=NamedTuple[]
for id in 0:6
    m=legacy_model(id);dir=joinpath(source,"model$id");p=read_parameters(joinpath(dir,"selected_parameters.csv"),m)
    row=NamedTuple(only(CSV.File(joinpath(dir,"summary.csv"))))
    actual=gaussian_nll(prep.train.y,prediction(prep.train,m,p),prep.train.sigma)
    isapprox(actual,row.nll;atol=1e-6,rtol=0) || error("Saved training fit differs: model$id")
    dest=joinpath(out,"model$id");mkpath(dest)
    for f in ("selected_parameters.csv","summary.csv","residuals.csv")
        cp(joinpath(dir,f),joinpath(dest,f);force=true)
    end
    push!(ranks,merge(row,(model=id,)))
    yp=prediction(prep.held,m,p);nll=gaussian_nll(prep.held.y,yp,prep.held.sigma)
    push!(scores,(model=publication_model(id),storage_model=id,n=length(yp),k=row.k,nll,bic=bic(nll,row.k,length(yp))))
    for state in (STATES...,"all")
        ix=findall(i->(state=="all" || prep.held.observations[i].state==state) && prep.held.y[i]!=0,eachindex(yp))
        push!(metrics,(model=publication_model(id),state,mrpe_percent=100mean(abs.((yp[ix].-prep.held.y[ix])./prep.held.y[ix]))))
    end
    for name in filter(n->occursin(r"^run\d+$",n),readdir(dir))
        run=parse(Int,name[4:end]);stagesfile=joinpath(dir,name,"stages.csv")
        isfile(stagesfile) || continue
        mkpath(joinpath(dest,name));cp(stagesfile,joinpath(dest,name,"stages.csv");force=true)
        run in 1:protocol["starts"] || continue
        stage=first(sort(collect(CSV.File(stagesfile));by=r->r.nll))
        paramfile=joinpath(dir,name,string(stage.stage)*"_parameters.csv")
        isfile(paramfile) || error("Missing saved run parameters: $paramfile")
        rp=read_parameters(paramfile,m)
        check=gaussian_nll(prep.train.y,prediction(prep.train,m,rp),prep.train.sigma)
        isapprox(check,stage.nll;atol=1e-6,rtol=0) || error("Saved run likelihood differs: $paramfile")
        v=gaussian_nll(prep.held.y,prediction(prep.held,m,rp),prep.held.sigma)
        push!(runs,(model=publication_model(id),storage_model=id,run,bic=bic(v,row.k,length(yp)),stopping_met=stage.stopping_met))
    end
end
if get(ENV,"HPSC_ALLOW_SUMMARY_ONLY","0")!="1"
    for id in 0:6
        count(r->r.storage_model==id,runs)==protocol["starts"] || error("Missing optimizer checkpoints for model$id. Use full server results, not a writing bundle.")
    end
end
sort!(scores;by=r->r.bic)
first(scores).storage_model==5 || error("Selection changed from Model 6. Review selection_scores before downstream analyses")
CSV.write(joinpath(out,"ranking.csv"),sort(ranks;by=r->r.bic))
CSV.write(joinpath(out,"selection_scores.csv"),scores)
CSV.write(joinpath(out,"selection_state_errors.csv"),metrics)
!isempty(runs) && CSV.write(joinpath(out,"selection_optimizer_bic.csv"),runs)
CSV.write(joinpath(out,"selection_lock.csv"),[(model=5,publication_model=6,criterion="40% selection BIC",parameter_sha256=bytes2hex(sha256(read(joinpath(out,"model5/selected_parameters.csv")))) )])
write(joinpath(out,"COMPLETE"),"Verified 19 training and 12 selection cultures. Copied parameters without refitting.\n")
println("Verified training-only parameter source: $out")
