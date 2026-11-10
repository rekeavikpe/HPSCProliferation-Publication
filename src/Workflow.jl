module Workflow
using CSV, TOML, Random, Statistics, LinearAlgebra, Optim, SHA, Dates
include("Calibration.jl")
using .Calibration, .Calibration.Proliferation
include("BoundedLeastSquares.jl")
using .BoundedLeastSquares
export context, prediction, fit_start, assemble, initial_starts, run_calibration, save_parameters, read_parameters, ROOT, STATES, legacy_model, config
const ROOT = normpath(joinpath(@__DIR__, ".."))
const STATES = Calibration.STATES
legacy_model(id) = UNIQUE_MODELS[id+1]
config() = TOML.parsefile(joinpath(ROOT,"config/calibration.toml"))

function context(; ids=nothing, conditions=(1,2,3,5,6), cfg=config(), fallback_glucose=Dict())
    series = collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    measurements = collect(CSV.File(joinpath(ROOT,"data/processed/measurements.csv")))
    # A saved training-only adapter carries culture IDs through profiles and refits.
    # Explicit IDs override this policy. The external 500 mL cohort is unaffected.
    if ids === nothing && haskey(cfg,"calibration_cohort") && all(c -> c in (1,2,3,5,6), conditions)
        ids = cfg["calibration_cohort"]["series_ids"]
        isempty(ids) && error("Training culture list is empty")
    end
    keep = ids === nothing ? Set(String(r.series_id) for r in series if r.condition in conditions) : Set(ids)
    rows = filter(r -> String(r.series_id) in keep, measurements)
    obs = NamedTuple[]
    for c in conditions, d in 1:6, state in STATES
        vals = Float64[r.value for r in rows if r.condition==c && r.day==d && r.state==state && r.status=="measured"]
        isempty(vals) && continue
        push!(obs,(condition=c,day=d,state=state,mean=mean(vals),sd=length(vals)>1 ? std(vals) : 0.,n_vessels=length(vals)))
    end
    feeds = collect(CSV.File(joinpath(ROOT,"data/processed/feeding.csv")))
    exposure = collect(CSV.File(joinpath(ROOT,"data/processed/lagrangian_exposure.csv")))
    setups = Dict()
    for c in conditions
        info = first(filter(r -> r.condition==c,series))
        glucose_values = Float64[r.value for r in rows if r.condition==c && r.day==0 && r.state=="glucose" && r.status=="measured"]
        isempty(glucose_values) && !haskey(fallback_glucose,c) && error("No measured day-zero glucose for condition $c")
        glucose = isempty(glucose_values) ? fallback_glucose[c] : mean(glucose_values)
        grid = sort(filter(r -> r.volume_ml==info.volume_ml && r.rpm==info.rpm && r.percentile==cfg["exposure_percentile"],exposure);by=r->r.size_um)
        lookup = linear_exposure([r.size_um for r in grid],[r.shear_pa for r in grid],[r.edr_m2_s3 for r in grid])
        times = sort(unique(Float64[24r.day for r in obs if r.condition==c]))
        isempty(times) && error("No observations for condition $c")
        events = [(time_h=Float64(r.time_h),fraction=r.exchange_fraction) for r in feeds if r.condition==c && r.time_h<=last(times)]
        setups[c] = (;u0=[Float64(info.seeding_cells_ml),glucose,cfg["initial_lactate_mM"],cfg["initial_aggregate_um"]],glucose,lookup,times,feeds=events)
    end
    y = Float64[r.mean for r in obs]
    sigma = Float64[max(r.sd,r.state in ("glucose","lactate") ? cfg["likelihood"]["metabolite_floor_mM"] : cfg["likelihood"]["relative_floor"]*abs(r.mean)) for r in obs]
    (;config=cfg,observations=obs,conditions=setups,y,sigma)
end

function prediction(ctx,model,p;reltol=1e-12,abstol=1e-13)
    values = zeros(length(ctx.y))
    for c in sort(collect(keys(ctx.conditions)))
        q = ctx.conditions[c]
        sim = simulate(model,p,q.u0,q.times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,reltol,abstol,lactate_inhibition=get(ctx.config,"lactate_inhibition",true))
        all(x->isfinite(x)&&x>=0,sim.predictions) || error("Invalid ODE prediction")
        for (i,r) in enumerate(ctx.observations)
            r.condition==c || continue
            values[i] = sim.predictions[findfirst(==(r.state),STATES),findfirst(==(24r.day),q.times)]
        end
    end
    values
end

save_parameters(path,p) = CSV.write(path,[(parameter=string(n),value=getproperty(p,n)) for n in keys(p)])
function read_parameters(path,model)
    values=Dict(Symbol(r.parameter)=>Float64(r.value) for r in CSV.File(path))
    parameters(model,[values[n] for n in parameter_names(model)])
end
function residual_function(ctx,model;reltol=1e-11,abstol=1e-12)
    z -> try
        (prediction(ctx,model,decode(ctx,model,z);reltol,abstol).-ctx.y)./ctx.sigma
    catch e
        e isa InterruptException && rethrow()
        fill(1e6,length(ctx.y))
    end
end
function initial_starts(ctx,model;nstarts=6,seed=20261003)
    rng=MersenneTwister(seed+model)
    residual=residual_function(ctx,model)
    starts=Vector{Float64}[]
    for run in 1:nstarts
        z=run==1 ? fill(.5,length(free_names(ctx,model))) : clamp.(.5 .+ .18randn(rng,length(free_names(ctx,model))),.02,.98)
        attempts=0
        while sum(abs2,residual(z))/2>1e10
            attempts+=1; attempts<=100 || error("Cannot find feasible initial parameters")
            z=clamp.(.5 .+ .15randn(rng,length(z)),.02,.98)
        end
        push!(starts,z)
    end
    starts
end
function fit_start(ctx,model,z0,dir;lm_iterations=300,polish_iterations=3000,tight=false,simplex_scale=.0002)
    mkpath(dir)
    residual=residual_function(ctx,model;reltol=tight ? 1e-12 : 1e-11,abstol=tight ? 1e-13 : 1e-12)
    objective(z)=all(x->isfinite(x)&&0<=x<=1,z) ? sum(abs2,residual(z))/2 : 1e100
    CSV.write(joinpath(dir,"start.csv"),[(parameter=string(n),normalized_log_value=z) for (n,z) in zip(free_names(ctx,model),z0)])
    t=time()
    z,lm=bounded_lm(residual,z0;maxiter=lm_iterations,h=1e-5)
    fit=optimize(objective,z,NelderMead(initial_simplex=Optim.AffineSimplexer(a=simplex_scale,b=0.)),Optim.Options(iterations=polish_iterations,g_abstol=1e-8))
    stages=NamedTuple[]
    for (stage,zz,stop) in (("least_squares",z,lm.converged),("polished",Optim.minimizer(fit),Optim.converged(fit)))
        p=decode(ctx,model,zz)
        nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
        save_parameters(joinpath(dir,stage*"_parameters.csv"),p)
        push!(stages,(stage=stage,nll=nll,stopping_met=stop))
    end
    CSV.write(joinpath(dir,"stages.csv"),stages)
    open(joinpath(dir,"diagnostics.toml"),"w") do io
        TOML.print(io,Dict("elapsed_seconds"=>time()-t,"lm_reason"=>lm.reason,"lm_iterations"=>lm.iterations,"lm_projected_gradient"=>lm.projected_gradient,"polish_iterations"=>Optim.iterations(fit),"simplex_scale"=>simplex_scale))
    end
    write(joinpath(dir,"COMPLETE"),"complete\n")
    stages
end
function assemble(ctx,model,dir)
    entries=NamedTuple[]
    for run in sort(filter(f->startswith(f,"run"),readdir(dir)))
        isfile(joinpath(dir,run,"COMPLETE")) || error("Incomplete start: $run")
        for r in CSV.File(joinpath(dir,run,"stages.csv"))
            push!(entries,(run=run,stage=String(r.stage),nll=r.nll,stopping_met=r.stopping_met))
        end
    end
    best=minimum(r.nll for r in entries)
    eligible=filter(r->r.stopping_met && r.nll<=best+1e-6,entries)
    chosen=first(sort(isempty(eligible) ? entries : eligible;by=r->r.nll))
    p=read_parameters(joinpath(dir,chosen.run,chosen.stage*"_parameters.csv"),model)
    save_parameters(joinpath(dir,"selected_parameters.csv"),p)
    pred=prediction(ctx,model,p)
    nll=gaussian_nll(ctx.y,pred,ctx.sigma)
    z=encode(ctx,model,p)
    bounds=join([string(n)*":"*(v<1e-4 ? "lower" : "upper") for (n,v) in zip(free_names(ctx,model),z) if v<1e-4||v>1-1e-4],";")
    # Values fitted on this dataset and subsequently fixed are not externally known.
    k=length(parameter_names(model))
    summary=(model=findfirst(==(model),UNIQUE_MODELS)-1,legacy_model=model,n=length(ctx.y),k=k,nll=nll,bic=bic(nll,k,length(ctx.y)),stopping_met=chosen.stopping_met,boundaries=bounds,run=chosen.run,stage=chosen.stage)
    CSV.write(joinpath(dir,"summary.csv"),[summary])
    CSV.write(joinpath(dir,"residuals.csv"),[(condition=r.condition,day=r.day,state=r.state,observed=ctx.y[i],predicted=pred[i],sigma=ctx.sigma[i]) for (i,r) in enumerate(ctx.observations)])
    summary
end
function run_calibration(out;ctx=context(),nstarts=6,seed=20261003,lm_iterations=300,polish_iterations=3000,refine=true)
    BLAS.set_num_threads(1)
    ispath(out) && error("Output already exists: use a new run directory")
    mkpath(out)
    open(joinpath(out,"run.toml"),"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads(),"seed"=>seed,"starts"=>nstarts,"lm_iterations"=>lm_iterations,"polish_iterations"=>polish_iterations,"refine"=>refine,"created"=>string(now()),"status"=>"running"))
    end
    effective=deepcopy(ctx.config)
    effective["analysis_id"]=get(ctx.config,"analysis_id","100mL-all-free-p50")
    merge!(effective["optimization"],Dict("seed"=>seed,"starts"=>nstarts,"local_iterations"=>lm_iterations,"polish_iterations"=>polish_iterations,"note"=>"Bounds-based seeded starts; own-optimum and training-only anchor refinements when enabled; no global optimum guarantee."))
    open(joinpath(out,"calibration.toml"),"w") do io; TOML.print(io,effective); end
    CSV.write(joinpath(out,"observations.csv"),ctx.observations)
    sourcefiles=[joinpath(ROOT,f) for f in ("Project.toml","Manifest.toml")]
    for sub in ("src","scripts","config","data/processed")
        append!(sourcefiles,[joinpath(ROOT,sub,f) for f in readdir(joinpath(ROOT,sub)) if isfile(joinpath(ROOT,sub,f))])
    end
    CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in sort(sourcefiles)])
    tasks=[(model=m,run=r,z=z) for m in UNIQUE_MODELS for (r,z) in enumerate(initial_starts(ctx,m;nstarts,seed))]
    Threads.@threads :dynamic for i in eachindex(tasks)
        task=tasks[i]
        dir=joinpath(out,"model$(findfirst(==(task.model),UNIQUE_MODELS)-1)","run$(task.run)")
        fit_start(ctx,task.model,task.z,dir;lm_iterations,polish_iterations)
        @info "Completed calibration start" model=findfirst(==(task.model),UNIQUE_MODELS)-1 start=task.run
    end
    summaries=[assemble(ctx,m,joinpath(out,"model$(i-1)")) for (i,m) in enumerate(UNIQUE_MODELS)]
    if refine
        anchor=read_parameters(joinpath(out,"model5/selected_parameters.csv"),6)
        tasks=NamedTuple[]
        for (i,m) in enumerate(UNIQUE_MODELS)
            dir=joinpath(out,"model$(i-1)")
            own=read_parameters(joinpath(dir,"selected_parameters.csv"),m)
            vals=[hasproperty(anchor,n) ? getproperty(anchor,n) : sqrt(prod(ctx.config["bounds"][string(n)])) for n in parameter_names(m)]
            for (r,p) in ((nstarts+1,own),(nstarts+2,parameters(m,vals)))
                push!(tasks,(model=m,run=r,z=encode(ctx,m,p),dir=dir))
            end
        end
        Threads.@threads :dynamic for i in eachindex(tasks)
            task=tasks[i]
            fit_start(ctx,task.model,task.z,joinpath(task.dir,"run$(task.run)");lm_iterations,polish_iterations=max(polish_iterations,5000),tight=true)
            @info "Completed refinement" model=findfirst(==(task.model),UNIQUE_MODELS)-1 start=task.run
        end
        summaries=[assemble(ctx,m,joinpath(out,"model$(i-1)")) for (i,m) in enumerate(UNIQUE_MODELS)]
    end
    sort!(summaries;by=r->r.bic)
    CSV.write(joinpath(out,"ranking.csv"),summaries)
    CSV.write(joinpath(out,"selection_lock.csv"),[(model=first(summaries).model,bic=first(summaries).bic,parameter_sha256=bytes2hex(sha256(read(joinpath(out,"model$(first(summaries).model)/selected_parameters.csv")))))])
    metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete";metadata["finished"]=string(now())
    open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
    write(joinpath(out,"COMPLETE"),"calibration complete; inspect convergence before interpretation\n")
    summaries
end
end
