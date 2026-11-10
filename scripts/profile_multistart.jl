using CSV, TOML, Random, SHA
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)>=4||error("Usage: profile_multistart.jl FIT_DIRECTORY OUTPUT_DIRECTORY MODEL PARAMETER [GRID_POINTS=21]")
fitdir=abspath(ARGS[1]);out=abspath(ARGS[2]);id=parse(Int,ARGS[3]);name=Symbol(ARGS[4]);npoints=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 21
npoints>=3||error("At least three grid points required")
fullctx=checked_context(fitdir);m=legacy_model(id)
name in parameter_names(m)||error("Parameter absent from model")
p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),m)
lo,hi=fullctx.config["bounds"][string(name)]
grid=sort(unique(vcat(exp.(range(log(lo),log(hi);length=npoints)),getproperty(p,name))))
resume=length(ARGS)>=6&&ARGS[6]=="--resume"
if resume
    metadata=TOML.parsefile(joinpath(out,"run.toml"))
    metadata["model"]==id&&metadata["parameter"]==string(name)&&metadata["source_fit"]==fitdir||error("Resume configuration differs")
    savedgrid=collect(CSV.File(joinpath(out,"grid.csv")))
    [r.fixed_value for r in savedgrid]==grid||error("Resume grid differs")
    for r in CSV.File(joinpath(out,"source_hashes.csv"))
        path=normpath(joinpath(ROOT,r.path))
        if basename(path) in ("Proliferation.jl","Workflow.jl","Calibration.jl","BoundedLeastSquares.jl","Manifest.toml","calibration.toml","observations.csv","ranking.csv","selected_parameters.csv")
            bytes2hex(sha256(read(path)))==r.sha256||error("Scientific source changed: $(r.path)")
        end
    end
    files=[@__FILE__,joinpath(ROOT,"src/TargetedChecks.jl")]
    CSV.write(joinpath(out,"resume_code_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in files])
else
    prepare_output(out,fitdir;settings=Dict("analysis"=>"multistart profile with bidirectional neighboring starts","model"=>id,"parameter"=>string(name),"seed"=>20261025,"grid_points"=>length(grid),"lm_iterations"=>600,"polish_iterations"=>12000))
end
CSV.write(joinpath(out,"grid.csv"),[(point=i,fixed_value=v) for (i,v) in enumerate(grid)])
previous=joinpath(fitdir,"profiles/model$id/$(name)")
previous_rows=isfile(joinpath(previous,"profile.csv")) ? collect(CSV.File(joinpath(previous,"profile.csv"))) : []
contexts=Any[];jobs=NamedTuple[]
for (i,value) in enumerate(grid)
    cfg=deepcopy(fullctx.config);cfg["fixed"][string(name)]=value
    ctx=context(cfg=cfg);push!(contexts,ctx)
    seeds=local_starts(ctx,m,p;n=3,seed=20261025+100id+i)
    if !isempty(previous_rows)
        nearest=argmin([abs(log(r.fixed_value/value)) for r in previous_rows])
        row=previous_rows[nearest]
        old=read_parameters(joinpath(previous,"point$nearest/$(row.stage)_parameters.csv"),m)
        seeds[2]=encode(ctx,m,old)
        seeds[3]=clamp.(seeds[2] .+ .03randn(MersenneTwister(20261025+i),length(seeds[2])),0.,1.)
    end
    for (j,z) in enumerate(seeds)
        dir=joinpath(out,"point$i/run$j")
        resume&&isfile(joinpath(dir,"COMPLETE"))&&continue
        isfile(joinpath(dir,"FAILED"))&&mv(joinpath(dir,"FAILED"),joinpath(dir,"FAILED_previous_$(time_ns())"))
        push!(jobs,(ctx=ctx,model=m,z=z,dir=dir))
    end
end
run_jobs(jobs)
firstpass=[read_parameters(collect_best(joinpath(out,"point$i")).choice.path,m) for i in eachindex(grid)]
neighbors=NamedTuple[]
for i in eachindex(grid)
    for (j,neighbor) in enumerate(unique([max(1,i-1),min(length(grid),i+1)]))
        neighbor==i&&continue
        dir=joinpath(out,"point$i/run$(3+j)")
        resume&&isfile(joinpath(dir,"COMPLETE"))&&continue
        isfile(joinpath(dir,"FAILED"))&&mv(joinpath(dir,"FAILED"),joinpath(dir,"FAILED_previous_$(time_ns())"))
        push!(neighbors,(ctx=contexts[i],model=m,z=encode(contexts[i],m,firstpass[neighbor]),dir=dir))
    end
end
run_jobs(neighbors)
baseline=only(CSV.File(joinpath(fitdir,"model$id/summary.csv"))).nll
rows=NamedTuple[]
for (i,value) in enumerate(grid)
    result=collect_best(joinpath(out,"point$i"));best=result.choice
    q=read_parameters(best.path,m)
    actual=gaussian_nll(contexts[i].y,prediction(contexts[i],m,q),contexts[i].sigma)
    isapprox(actual,best.nll;atol=1e-6,rtol=0)||error("Objective verification failed")
    save_parameters(joinpath(out,"point$i/selected_parameters.csv"),q)
    runs=unique(e.run for e in result.entries)
    agreeing=count(run->minimum(e.nll for e in result.entries if e.run==run)<=result.minimum_nll+1e-5,runs)
    push!(rows,(point=i,parameter=string(name),fixed_value=value,nll=actual,twice_delta_nll=2(actual-baseline),stopping_met=best.stopping_met,agreeing_starts=agreeing,total_starts=length(runs),run=best.run,stage=best.stage))
end
CSV.write(joinpath(out,"profile.csv"),rows)
metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete"
open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
write(joinpath(out,"COMPLETE"),"complete; agreement is not a guarantee of the global nuisance optimum; no automatic confidence intervals\n")
println("Completed profile: $out")
