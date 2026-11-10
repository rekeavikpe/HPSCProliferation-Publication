module TargetedChecks
using CSV, TOML, Random, SHA, LinearAlgebra
include("Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
export Workflow, checked_context, prepare_output, collect_best, save_choice, local_starts, run_jobs

function checked_context(fitdir)
    cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml"))
    ctx=context(cfg=cfg)
    saved=collect(CSV.File(joinpath(fitdir,"observations.csv")))
    length(saved)==length(ctx.observations)||error("Cohort differs from saved fit")
    for (a,b) in zip(saved,ctx.observations)
        all(getproperty(a,k)==getproperty(b,k) for k in (:condition,:day,:state,:n_vessels))||error("Observation identities differ")
        isapprox(a.mean,b.mean;rtol=1e-12)&&isapprox(a.sd,b.sd;rtol=1e-12)||error("Observation values differ")
    end
    isempty(cfg["fixed"])||error("Expected all-free calibration")
    ctx
end

function prepare_output(out,fitdir;settings=Dict())
    ispath(out)&&error("Choose a new output directory: $out")
    mkpath(out)
    settings=merge(Dict("source_fit"=>abspath(fitdir),"julia"=>string(VERSION),"threads"=>Threads.nthreads(),"status"=>"running"),settings)
    open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,settings);end
    files=[joinpath(ROOT,f) for f in ("src/Proliferation.jl","src/Workflow.jl","src/Calibration.jl","src/BoundedLeastSquares.jl","src/TargetedChecks.jl","Manifest.toml")]
    append!(files,[joinpath(ROOT,"scripts",f) for f in readdir(joinpath(ROOT,"scripts")) if endswith(f,".jl")])
    append!(files,[joinpath(fitdir,"calibration.toml"),joinpath(fitdir,"observations.csv"),joinpath(fitdir,"ranking.csv")])
    append!(files,[joinpath(fitdir,"model$id/selected_parameters.csv") for id in 0:6])
    CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in files])
end

function local_starts(ctx,m,p;n=4,seed=20261005)
    rng=MersenneTwister(seed);z=encode(ctx,m,p)
    residual=Workflow.residual_function(ctx,m)
    starts=Vector{Float64}[]
    for i in 1:n
        candidate=i==1 ? copy(z) : clamp.(z .+ (i==2 ? .003 : i==3 ? .03 : .1).*randn(rng,length(z)),0.,1.)
        for attempt in 1:20
            sum(abs2,residual(candidate))/2<1e10 && break
            candidate=(candidate+z)/2
        end
        sum(abs2,residual(candidate))/2<1e10||error("No feasible local start")
        push!(starts,candidate)
    end
    starts
end

function run_jobs(jobs;lm_iterations=600,polish_iterations=12000)
    BLAS.set_num_threads(1)
    Threads.@threads :greedy for i in eachindex(jobs)
        j=jobs[i]
        try
            fit_start(j.ctx,j.model,j.z,j.dir;lm_iterations,polish_iterations,tight=true)
        catch e
            e isa InterruptException&&rethrow()
            mkpath(j.dir)
            write(joinpath(j.dir,"FAILED"),sprint(showerror,e)*"\n")
            @warn "Targeted start failed; retained for audit" folder=j.dir exception=sprint(showerror,e)
            continue
        end
        @info "Targeted start complete" folder=j.dir
    end
end

function collect_best(folder)
    entries=NamedTuple[]
    for name in sort(filter(n->startswith(n,"run"),readdir(folder)))
        dir=joinpath(folder,name)
        isfile(joinpath(dir,"FAILED"))&&continue
        isfile(joinpath(dir,"COMPLETE"))||error("Incomplete start: $dir")
        for r in CSV.File(joinpath(dir,"stages.csv"))
            push!(entries,(run=name,stage=String(r.stage),nll=r.nll,stopping_met=r.stopping_met,path=joinpath(dir,r.stage*"_parameters.csv")))
        end
    end
    isempty(entries)&&error("No valid completed starts in $folder")
    best=minimum(r.nll for r in entries)
    converged=filter(r->r.stopping_met&&r.nll<=best+1e-6,entries)
    choice=first(sort(isempty(converged) ? entries : converged;by=r->r.nll))
    (choice=choice,minimum_nll=best,entries=entries)
end

function save_choice(ctx,m,folder;estimated_parameters=length(parameter_names(m)))
    result=collect_best(folder);r=result.choice;p=read_parameters(r.path,m)
    actual=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
    isapprox(actual,r.nll;atol=1e-6,rtol=0)||error("Saved objective mismatch")
    save_parameters(joinpath(folder,"selected_parameters.csv"),p)
    z=encode(ctx,m,p)
    bounds=join([string(n)*":"*(v<1e-4 ? "lower" : "upper") for (n,v) in zip(free_names(ctx,m),z) if v<1e-4||v>1-1e-4],";")
    failed=count(n->isfile(joinpath(folder,n,"FAILED")),readdir(folder))
    row=(model=findfirst(==(m),UNIQUE_MODELS)-1,nll=actual,bic=bic(actual,estimated_parameters,length(ctx.y)),k=estimated_parameters,stopping_met=r.stopping_met,boundaries=bounds,run=r.run,stage=r.stage,minimum_nll=result.minimum_nll,failed_starts=failed)
    CSV.write(joinpath(folder,"summary.csv"),[row])
    row
end
end
