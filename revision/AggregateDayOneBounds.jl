module AggregateDayOneBounds

include("AggregateDayOne.jl")
using .AggregateDayOne, .AggregateDayOne.MixedExposure
using .MixedExposure.RevisionCommon
using .RevisionCommon.Workflow, .Workflow.Calibration, .Calibration.Proliferation
using CSV, TOML, SHA, Dates, Printf, LinearAlgebra, Statistics

const MODEL = 6
const TARGETS = (:K_G, :K_eps)
const MULTIPLIERS = (1, 10, 100)
const ORIGINAL_UPPER = (K_G=100.0, K_eps=5.0)
const GRID = (K_G=[1.,3.,10.,30.,100.,300.,1000.,3000.,10000.],
              K_eps=[1e-5,.001,.01,.1,1.,5.,15.,50.,150.,500.])

sha(path) = bytes2hex(sha256(read(path)))
actual(ctx,p) = gaussian_nll(ctx.y,prediction(ctx,MODEL,p),ctx.sigma)
function atomic_text(path,text)
    mkpath(dirname(path)); temp=path*".tmp"
    write(temp,text); mv(temp,path;force=true)
end
function atomic_csv(path,rows)
    mkpath(dirname(path)); temp=path*".tmp"
    CSV.write(temp,rows); mv(temp,path;force=true)
end
function canonical(path)
    path=abspath(path)
    ispath(path) && return realpath(path)
    joinpath(canonical(dirname(path)),basename(path))
end
inside(path,parent) = path==parent || startswith(path,parent*"/")
function safe_output(out,inputs)
    dest=canonical(out)
    for source in inputs
        src=canonical(source)
        (inside(dest,src)||inside(src,dest)) && error("Output overlaps a protected input: $source")
    end
    dest
end
function snapshot(root)
    Dict(relpath(joinpath(dir,f),root)=>sha(joinpath(dir,f))
         for (dir,_,files) in walkdir(root) for f in files)
end

"Change only the two upper search limits. Fixed-value scans hold one target, never the other ten."
function bounded_context(train,multiplier;fixed=Dict{String,Float64}())
    multiplier in MULTIPLIERS || error("Unplanned bound multiplier")
    cfg=deepcopy(train.config)
    isempty(cfg["fixed"]) || error("Expected eleven free parameters in the source model")
    for target in TARGETS
        bounds=cfg["bounds"][string(target)]
        bounds[2]==getproperty(ORIGINAL_UPPER,target) || error("Original bounds changed")
        cfg["bounds"][string(target)]=[Float64(bounds[1]),multiplier*Float64(bounds[2])]
    end
    all(Symbol(n) in TARGETS for n in keys(fixed)) || error("Unexpected fixed parameter")
    for (name,value) in fixed
        lo,hi=cfg["bounds"][name]
        lo<=value<=hi || error("Fixed value outside diagnostic range")
    end
    cfg["fixed"]=copy(fixed)
    AggregateDayOneContext(cfg,train.observations,train.conditions,train.y,train.sigma)
end

function verify_parent(source,fit,splitfit)
    isfile(joinpath(source,"COMPLETE")) || error("The original aggregate-D1 analysis is incomplete")
    manifest=TOML.parsefile(joinpath(source,"manifest.toml"))
    manifest["settings"]["analysis"]=="model6_p50p50_aggregate_only_day1_fixed15_v2" || error("Wrong D1 design")
    manifest["settings"]["starts"]==100 || error("Expected the completed 100-start analysis")
    for (key,hash) in manifest["hashes"]
        path=if startswith(key,"results/training60_final/reference_fit/")
            joinpath(fit,relpath(key,"results/training60_final/reference_fit"))
        elseif startswith(key,"results/split6040_refined_hybrid_lac100_kg100/")
            joinpath(splitfit,relpath(key,"results/split6040_refined_hybrid_lac100_kg100"))
        else
            joinpath(ROOT,key)
        end
        isfile(path) && sha(path)==hash || error("Original analysis input differs: $key")
    end
    manifest
end

"Seeds are ranked exclusively by saved training likelihood, then numerically verified."
function saved_seeds(source,train,count)
    fit=joinpath(source,"aggregate_day1_fit")
    entries=NamedTuple[]
    for run in 1:100
        dir=joinpath(fit,"run$run")
        isfile(joinpath(dir,"COMPLETE")) || error("Missing original checkpoint: run$run")
        for row in CSV.File(joinpath(dir,"stages.csv"))
            push!(entries,(path=joinpath(dir,string(row.stage)*"_parameters.csv"),nll=Float64(row.nll)))
        end
    end
    sort!(entries;by=r->r.nll)
    selected=read_parameters(joinpath(fit,"selected_parameters.csv"),MODEL)
    nll=actual(train,selected)
    recorded=only(CSV.File(joinpath(fit,"summary.csv"))).nll
    isapprox(nll,recorded;atol=1e-6,rtol=0) || error("Saved D1 likelihood is not reproducible")
    isapprox(nll,first(entries).nll;atol=1e-6,rtol=0) || error("Saved D1 fit is not the best saved training stage")
    seeds=[selected]; paths=[joinpath(fit,"selected_parameters.csv")]
    for row in entries
        p=read_parameters(row.path,MODEL)
        any(maximum(abs.(log.(collect(values(p))./collect(values(q)))))<1e-8 for q in seeds) && continue
        isapprox(actual(train,p),row.nll;atol=1e-6,rtol=0) || error("Saved seed likelihood differs: $(row.path)")
        push!(seeds,p);push!(paths,row.path)
        length(seeds)>=count && break
    end
    selected,seeds,paths
end

boundary_string(ctx,p) = join([string(name)*":"*(z<=1e-5 ? "lower" : "upper")
    for (name,z) in zip(free_names(ctx,MODEL),encode(ctx,MODEL,p)) if z<=1e-5||z>=1-1e-5]," | ")

function unique_seeds(ctx,seeds;limit=typemax(Int),feasible=false)
    result=typeof(first(seeds))[]
    for p in seeds
        options=[p]
        if haskey(ctx.config["fixed"],"K_G")
            # A seed only: preserve its glucose growth factor at the median initial glucose.
            # mu_max remains free during every nuisance optimization.
            g=median([q.glucose for q in values(ctx.conditions)])
            mu=p.mu_max*(g+ctx.config["fixed"]["K_G"])/(g+p.K_G)
            lo,hi=ctx.config["bounds"]["mu_max"]
            push!(options,merge(p,(mu_max=clamp(mu,lo,hi),)))
        end
        for option in options
            # decode applies the scan's fixed value and clamps starting values to its box.
            q=decode(ctx,MODEL,encode(ctx,MODEL,option))
            any(maximum(abs.(encode(ctx,MODEL,q).-encode(ctx,MODEL,r)))<1e-9 for r in result) && continue
            if feasible
                good=try isfinite(actual(ctx,q)) catch e
                    e isa InterruptException && rethrow()
                    false
                end
                good || continue
            end
            push!(result,q)
            length(result)>=limit && return result
        end
    end
    result
end

"Each checkpoint retains its seed if optimization is worse, and records true optimizer status."
function optimize_one(ctx,seed,dir;lm_iterations,polish_iterations)
    mkpath(dir)
    input=Dict("seed"=>Dict(string(k)=>v for (k,v) in pairs(seed)),
               "bounds"=>ctx.config["bounds"],"fixed"=>ctx.config["fixed"],
               "lm_iterations"=>lm_iterations,"polish_iterations"=>polish_iterations)
    spec=joinpath(dir,"input.toml")
    if isfile(spec)
        TOML.parsefile(spec)==input || error("Checkpoint inputs changed: $dir")
    else
        atomic_text(spec,sprint(io->TOML.print(io,input)))
    end
    if isfile(joinpath(dir,"COMPLETE"))
        p=read_parameters(joinpath(dir,"selected_parameters.csv"),MODEL)
        row=only(CSV.File(joinpath(dir,"summary.csv")))
        isapprox(actual(ctx,p),row.nll;atol=1e-6,rtol=0) || error("Invalid diagnostic checkpoint: $dir")
        return (;p,nll=Float64(row.nll),stopping_met=Bool(row.stopping_met),stage=String(row.stage),dir)
    end
    z=encode(ctx,MODEL,seed);pseed=decode(ctx,MODEL,z)
    seed_nll=actual(ctx,pseed)
    save_parameters(joinpath(dir,"seed_parameters.csv"),pseed)
    stages=Workflow.fit_start(ctx,MODEL,z,joinpath(dir,"optimization");lm_iterations,
        polish_iterations,tight=true,simplex_scale=.01)
    candidates=[(p=pseed,nll=seed_nll,stopping_met=false,stage="seed")]
    for row in stages
        p=read_parameters(joinpath(dir,"optimization",row.stage*"_parameters.csv"),MODEL)
        push!(candidates,(;p,nll=Float64(row.nll),stopping_met=Bool(row.stopping_met),stage=row.stage))
    end
    minimum_nll=minimum(r.nll for r in candidates)
    # Prefer a converged result only if equivalent at numerical precision.
    converged=filter(r->r.stopping_met && r.nll<=minimum_nll+1e-7,candidates)
    best=first(sort(isempty(converged) ? candidates : converged;by=r->r.nll))
    best.nll<=seed_nll+1e-7 || error("Optimization degraded its starting solution")
    save_parameters(joinpath(dir,"selected_parameters.csv"),best.p)
    atomic_csv(joinpath(dir,"summary.csv"),[(nll=best.nll,seed_nll,improvement=seed_nll-best.nll,
        stopping_met=best.stopping_met,stage=best.stage,boundaries=boundary_string(ctx,best.p))])
    atomic_text(joinpath(dir,"COMPLETE"),"Numerical checkpoint complete. Inspect optimization status.\n")
    (;best...,dir)
end

function optimize_bank(ctx,seeds,dir;lm_iterations,polish_iterations,progress)
    seeds=unique_seeds(ctx,seeds;feasible=true)
    isempty(seeds) && error("No feasible nuisance starting point in $dir. Inspect the recorded design, do not infer non-identifiability.")
    results=Vector{Any}(undef,length(seeds))
    Threads.@threads :greedy for i in eachindex(seeds)
        results[i]=optimize_one(ctx,seeds[i],joinpath(dir,"start$(lpad(i,2,'0'))");lm_iterations,polish_iterations)
        progress("DONE $(relpath(dir)) start $i/$(length(seeds)) NLL=$(round(results[i].nll;digits=6))")
    end
    best=first(sort(results;by=r->r.nll))
    save_parameters(joinpath(dir,"selected_parameters.csv"),best.p)
    atomic_csv(joinpath(dir,"summary.csv"),[(nll=best.nll,stopping_met=best.stopping_met,
        selected_start=basename(best.dir),starts=length(results),
        nll_spread=maximum(r.nll for r in results)-minimum(r.nll for r in results),
        boundaries=boundary_string(ctx,best.p))])
    atomic_text(joinpath(dir,"COMPLETE"),"Optimization bank complete.\n")
    best
end

function lock_run(out,settings,source,extra)
    files=vcat([@__FILE__,joinpath(@__DIR__,"run_aggregate_day1_bounds.jl"),
        joinpath(ROOT,"scripts/run_aggregate_day1_bounds.sh"),joinpath(ROOT,"tests/test_aggregate_day1_bounds.jl")],extra)
    hashes=Dict(abspath(f)=>sha(f) for f in unique(files))
    manifest=Dict("settings"=>settings,"hashes"=>hashes,"original_result_hashes"=>snapshot(source))
    path=joinpath(out,"manifest.toml")
    if isfile(path)
        TOML.parsefile(path)==manifest || error("Resume inputs/settings changed. Existing results preserved. Use a new output directory.")
    else
        isdir(out) && !isempty(readdir(out)) && error("Refusing to use an unrecognized nonempty output directory")
        atomic_text(path,sprint(io->TOML.print(io,manifest)))
    end
    manifest
end

function run_checks(fit,splitfit,source,out;smoke=false)
    fit,splitfit,source,out=abspath.((fit,splitfit,source,out))
    safe_output(out,[fit,splitfit,source,joinpath(ROOT,"results/mixed_exposure_p50shear_p75edr")])
    parent=verify_parent(source,fit,splitfit)
    prep,p0=MixedExposure.source_contexts(fit,splitfit)
    train=aggregate_dayone(prep.train)
    seed_count=smoke ? 2 : parse(Int,get(ENV,"HPSC_AGGBOUND_SEEDS","8"))
    profile_count=smoke ? 2 : parse(Int,get(ENV,"HPSC_AGGBOUND_SCAN_SEEDS","4"))
    lm=smoke ? 1 : parse(Int,get(ENV,"HPSC_AGGBOUND_LM_ITERATIONS","200"))
    polish=smoke ? 2 : parse(Int,get(ENV,"HPSC_AGGBOUND_POLISH_ITERATIONS","3000"))
    seed_count>=2 && profile_count>=2 && lm>0 && polish>0 || error("Invalid diagnostic budgets")
    saved,bank,seedpaths=saved_seeds(source,train,seed_count)
    # Hash every parent input in its current local/server location, plus the saved result tree.
    files=String[]
    for key in keys(parent["hashes"])
        path=startswith(key,"results/training60_final/reference_fit/") ? joinpath(fit,relpath(key,"results/training60_final/reference_fit")) :
             startswith(key,"results/split6040_refined_hybrid_lac100_kg100/") ? joinpath(splitfit,relpath(key,"results/split6040_refined_hybrid_lac100_kg100")) : joinpath(ROOT,key)
        push!(files,path)
    end
    settings=Dict("analysis"=>"aggregate_D1_boundary_diagnostics_v1","julia"=>string(VERSION),
        "source"=>source,"reference_fit"=>fit,"split_fit"=>splitfit,"smoke"=>smoke,
        "saved_seeds"=>seed_count,"scan_seeds"=>profile_count,"lm_iterations"=>lm,"polish_iterations"=>polish,
        "multipliers"=>collect(MULTIPLIERS),"grid"=>Dict(string(k)=>v for (k,v) in pairs(GRID)),
        "selection_used_for_parameter_fitting"=>false,"intervals"=>"diagnostic grid only, no confidence limits")
    manifest=lock_run(out,settings,source,files)
    BLAS.set_num_threads(1)
    loglock=ReentrantLock()
    function progress(s)
        lock(loglock) do
            line=string(Dates.now()," ",s)
            open(io->println(io,line),joinpath(out,"progress.log"),"a")
            println(line);flush(stdout)
        end
    end
    progress("START boundary checks. Original 100-start results are read-only.")
    atomic_csv(joinpath(out,"seed_sources.csv"),[(rank=i,path=seedpaths[i],nll=actual(train,p)) for (i,p) in enumerate(bank)])
    atomic_csv(joinpath(out,"split.csv"),prep.split.assignments)
    saved_stop=Bool(only(CSV.File(joinpath(source,"aggregate_day1_fit/summary.csv"))).stopping_met)
    candidates=Any[(label="D1 saved fit",p=saved,ctx=train,nll=actual(train,saved),stopping_met=saved_stop)]
    previous=saved
    for multiplier in MULTIPLIERS
        ctx=bounded_context(train,multiplier)
        seeds=vcat([previous],bank)
        if multiplier>1 && !smoke
            # Deterministic starts inside the newly opened range, in addition to all retained seeds.
            for names in ((:K_G,),(:K_eps,),TARGETS)
                vals=Dict(n=>getproperty(previous,n) for n in parameter_names(MODEL))
                for n in names; vals[n]=getproperty(ORIGINAL_UPPER,n)*sqrt(multiplier); end
                push!(seeds,parameters(MODEL,[vals[n] for n in parameter_names(MODEL)]))
            end
        end
        progress("START bounds x$multiplier")
        result=optimize_bank(ctx,seeds,joinpath(out,"bounds_x$multiplier");lm_iterations=lm,polish_iterations=polish,progress)
        push!(candidates,(label="D1 bounds x$multiplier",p=result.p,ctx,nll=result.nll,stopping_met=result.stopping_met))
        previous=result.p
        progress("FINISHED bounds x$multiplier training NLL=$(result.nll)")
    end
    wide=bounded_context(train,100)
    anchor=first(sort(candidates;by=r->r.nll)).p
    # Lock the anchor before either scan, so parallel scheduling cannot change starting values.
    save_parameters(joinpath(out,"scan_anchor_parameters.csv"),anchor)
    grids=Dict{Symbol,Vector{Float64}}()
    for n in TARGETS
        grid=copy(smoke ? [getproperty(ORIGINAL_UPPER,n)] : getproperty(GRID,n))
        any(isapprox(x,getproperty(anchor,n);rtol=1e-10,atol=0) for x in grid) || push!(grid,getproperty(anchor,n))
        grids[n]=sort(grid)
    end
    scanresults=Dict{Symbol,Any}(); scanlock=ReentrantLock()
    Threads.@threads :greedy for target in collect(TARGETS)
        grid=grids[target];center=argmin(abs.(log.(grid./getproperty(anchor,target))))
        order=vcat([center],collect(center-1:-1:1),collect(center+1:length(grid)))
        completed=Dict{Int,Any}()
        for index in order
            value=grid[index]
            ctx=bounded_context(train,100;fixed=Dict(string(target)=>value))
            neighboring=isempty(completed) ? [anchor] : [completed[first(sort(collect(keys(completed));by=i->abs(log(grid[i]/value))))].p]
            seeds=unique_seeds(ctx,vcat([anchor],neighboring,bank);limit=profile_count,feasible=true)
            dest=joinpath(out,"scans",string(target),"point$(lpad(index,2,'0'))")
            progress("START scan $target=$value ($index/$(length(grid)))")
            result=optimize_bank(ctx,seeds,dest;lm_iterations=lm,polish_iterations=polish,progress)
            completed[index]=result
            atomic_csv(joinpath(dest,"point.csv"),[(parameter=string(target),value,nll=result.nll,
                stopping_met=result.stopping_met,nuisance_parameters=10,boundaries=boundary_string(ctx,result.p))])
            progress("FINISHED scan $target=$value training NLL=$(result.nll)")
        end
        lock(scanlock) do; scanresults[target]=(;grid,completed); end
    end
    scan_candidates=[(label="scan $(target) point $i",p=r.p,ctx=wide,nll=r.nll,stopping_met=r.stopping_met)
        for target in TARGETS for (i,r) in scanresults[target].completed]
    allcandidates=vcat(candidates,scan_candidates)
    bestscan=first(sort(allcandidates;by=r->r.nll))
    progress("START final unrestricted refinement")
    final=optimize_bank(wide,[bestscan.p,anchor],joinpath(out,"final_refinement");lm_iterations=lm,polish_iterations=polish,progress)
    push!(allcandidates,(label="final unrestricted refinement",p=final.p,ctx=wide,nll=final.nll,stopping_met=final.stopping_met))
    best=first(sort(allcandidates;by=r->r.nll))
    save_parameters(joinpath(out,"best_training_parameters.csv"),best.p)
    atomic_csv(joinpath(out,"training_candidates.csv"),[(label=r.label,nll=r.nll,stopping_met=r.stopping_met,
        K_G=r.p.K_G,K_eps=r.p.K_eps,boundaries=boundary_string(r.ctx,r.p)) for r in allcandidates])
    scanrows=NamedTuple[]
    for target in TARGETS
        grid,completed=scanresults[target]
        for i in eachindex(grid)
            r=completed[i]
            ctx=bounded_context(train,100;fixed=Dict(string(target)=>grid[i]))
            push!(scanrows,(parameter=string(target),value=grid[i],nll=r.nll,delta_nll=r.nll-best.nll,
                stopping_met=r.stopping_met,nuisance_boundaries=boundary_string(ctx,r.p)))
        end
    end
    atomic_csv(joinpath(out,"likelihood_scans.csv"),scanrows)
    progress("START evaluation with fixed training estimates")
    evaluations=vcat([(label="Current D0 model",p=p0)],[(label=r.label,p=r.p) for r in candidates],
                     [(label="Best D1 diagnostic candidate",p=best.p)])
    scores=NamedTuple[]; metrics=NamedTuple[]; points=NamedTuple[]; curves=NamedTuple[]; initials=NamedTuple[]
    for (cohort,raw,settings) in MixedExposure.evaluations(prep)
        base=replace_exposure(raw;edr_percentile=50,settings)
        for r in evaluations
            ctx=r.label=="Current D0 model" ? comparison_context(base) : aggregate_dayone(base)
            append!(metrics,AggregateDayOne.metrics(ctx,prediction(ctx,MODEL,r.p),r.label,cohort))
            cohort in ("training60","selection40","validation500") && push!(scores,AggregateDayOne.score_row(ctx,r.p,r.label,cohort))
            AggregateDayOne.export_cohort!(points,curves,initials,base,ctx,r.p,r.label,cohort)
        end
    end
    atomic_csv(joinpath(out,"scores.csv"),scores)
    atomic_csv(joinpath(out,"selection_scores.csv"),filter(r->r.cohort=="selection40",scores))
    atomic_csv(joinpath(out,"state_errors.csv"),metrics)
    atomic_csv(joinpath(out,"predictions.csv"),points)
    atomic_csv(joinpath(out,"curves.csv"),curves)
    atomic_csv(joinpath(out,"initial_states.csv"),initials)
    snapshot(source)==manifest["original_result_hashes"] || error("Original result files changed while diagnostics ran")
    flags=(best_training_nll=best.nll,saved_training_nll=actual(train,saved),
        improvement=actual(train,saved)-best.nll,best_source=best.label,
        best_stopping_met=best.stopping_met,boundaries=boundary_string(wide,best.p),
        scan_points_unconfirmed=count(r->!r.stopping_met,scanrows),
        scan_found_better_than_free_fits=minimum(r.nll for r in scan_candidates)<minimum(r.nll for r in candidates)-1e-5,
        original_results_unchanged=true,confidence_intervals_computed=false,smoke_test=smoke)
    atomic_csv(joinpath(out,"audit.csv"),[flags])
    report(out,scores,metrics,scanrows,flags)
    progress("FINISHED evaluation. Review audit.csv and index.html before adopting any fit.")
    atomic_text(joinpath(out,smoke ? "SMOKE_COMPLETE" : "COMPLETE"),"Diagnostics completed. Original results preserved. No documents changed. No confidence intervals computed.\n")
    flags
end

function report(out,scores,metrics,scans,flags)
    open(joinpath(out,"index.html"),"w") do io
        print(io,"<!doctype html><meta charset='utf-8'><title>Aggregate D1 boundary checks</title><style>body{font:17px Arial;margin:35px auto;max-width:1400px}table{border-collapse:collapse;margin:20px 0}td,th{border:1px solid #bbb;padding:7px;text-align:right}td:first-child{text-align:left}</style><h1>Aggregate day 1: parameter-bound checks</h1><p>Training-only refinements. Same 19/12 culture split, P50/P50 inputs, fixed 15 µm pre-D1 input and cohort-specific measured D1 aggregate initialization. Aggregate D1 observations are excluded in all comparisons. Other-state observations are unchanged. External parameters are not refitted.</p><p>Upper limits for K_G and K_eps are expanded by factors of 10 and 100 for numerical sensitivity only. Other bounds are unchanged. These are not new biological limits. Each one-dimensional likelihood scan refits the other ten parameters. Sparse scans are diagnostic and do not supply confidence intervals or certify identifiability.</p>")
        print(io,"<h2>Audit</h2><table>")
        for (k,v) in pairs(flags);print(io,"<tr><td>$k</td><td>$v</td></tr>");end
        print(io,"</table><h2>Selection scores (111 matched observations, 11 parameters)</h2><table><tr><th>Fit</th><th>NLL</th><th>BIC</th></tr>")
        for r in scores
            r.cohort=="selection40" || continue
            @printf(io,"<tr><td>%s</td><td>%.4f</td><td>%.4f</td></tr>",r.formulation,r.nll,r.bic)
        end
        print(io,"</table><h2>MRPE (%) by dataset and state</h2><p>Overall errors pool scored observations. iPSC results are direct transfers, not cell-line recalibrations. All-100 mL results overlap training and selection.</p><table><tr><th>Fit</th><th>Dataset</th><th>State</th><th>MRPE</th></tr>")
        for r in metrics
            r.condition==0 || continue
            @printf(io,"<tr><td>%s</td><td>%s</td><td>%s</td><td>%.4f</td></tr>",r.formulation,r.cohort,r.state,r.mrpe_percent)
        end
        print(io,"</table><h2>Fixed-value likelihood diagnostics</h2><table><tr><th>Parameter</th><th>Value</th><th>ΔNLL from best training candidate</th><th>Optimization confirmed</th><th>Other bounds reached</th></tr>")
        for r in scans
            @printf(io,"<tr><td>%s</td><td>%.6g</td><td>%.6f</td><td>%s</td><td>%s</td></tr>",r.parameter,r.value,r.delta_nll,string(r.stopping_met),r.nuisance_boundaries)
        end
        print(io,"</table><p>Completion means computation finished. Inspect convergence, other active bounds and whether the curve is still changing at the largest tested value. Existing manuscript intervals do not apply to these new fits.</p>")
    end
end

end
