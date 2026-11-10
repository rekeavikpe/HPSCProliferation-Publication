using CSV, TOML, SHA, Serialization, ProfileLikelihood, SciMLBase
include("PackageWorkflow.jl")
using .PackageWorkflow
include("FixedTwoSeeds.jl")
using .FixedTwoSeeds
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation

# Model 6 has storage ID 5. This is a separate uncertainty analysis of the
# published P50/P50 training fit, not a new candidate or a day-one fit.
function fixed_two_context(fitdir)
    raw=checked_context(fitdir)
    cohort=raw.config["calibration_cohort"]
    training=String.(cohort["series_ids"]); selection=String.(cohort["selection_ids"])
    length(training)==length(unique(training))==19 || error("Expected 19 unique training cultures")
    length(selection)==length(unique(selection))==12 || error("Expected 12 unique selection cultures")
    isempty(intersect(training,selection)) || error("Training and selection overlap")
    split=collect(CSV.File(joinpath(fitdir,"split.csv")))
    Set(String(r.series_id) for r in split if r.role=="calibration")==Set(training) || error("Training IDs differ from split")
    Set(String(r.series_id) for r in split if r.role=="selection")==Set(selection) || error("Selection IDs differ from split")
    raw.config["exposure_percentile"]==50 || error("Expected P50 inputs")
    get(raw.config,"lactate_inhibition",true) || error("Lactate inhibition must remain active")
    isempty(get(raw.config,"parameter_ties",Dict())) || error("Unexpected parameter ties")
    model=legacy_model(5)
    p=read_parameters(joinpath(fitdir,"model5/selected_parameters.csv"),model)
    saved_nll=Float64(only(CSV.File(joinpath(fitdir,"model5/summary.csv"))).nll)
    source_nll=gaussian_nll(raw.y,prediction(raw,model,p),raw.sigma)
    isapprox(source_nll,saved_nll;atol=1e-6,rtol=0) || error("Saved training likelihood differs")
    fixed=Dict("K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX)
    cfg=deepcopy(raw.config); merge!(cfg["fixed"],fixed)
    ctx=context(cfg=cfg)
    ctx.observations==raw.observations && ctx.y==raw.y && ctx.sigma==raw.sigma || error("Profile cohort or weights changed")
    ctx.config["bounds"]==raw.config["bounds"] || error("Bounds changed")
    names=collect(free_names(ctx,model))
    length(names)==9 && :K_G in names || error("Expected nine free parameters including K_G")
    Set(string.(names))==setdiff(Set(string.(parameter_names(model))),Set(keys(fixed))) || error("Wrong free parameters")
    (;raw,ctx,model,p,saved_nll,fixed,names)
end

function safe_output(out,inputs)
    canonical(p)=ispath(p) ? realpath(p) : joinpath(canonical(dirname(abspath(p))),basename(p))
    target=canonical(abspath(out))
    for p in inputs
        source=canonical(abspath(p))
        (target==source || startswith(target,source*"/") || startswith(source,target*"/")) && error("Output overlaps protected source: $p")
    end
    target
end

# The old parallel seed bug did not change the likelihood, but it could choose
# the wrong starting guesses. Re-evaluate every saved point and re-optimize the
# points bracketing each confidence limit before accepting an imported profile.
function recheck_saved_profile(dest,names,index,bounds,sol,loglik,bank;maxiters=10000)
    rows=collect(CSV.File(joinpath(dest,"profile.csv")))
    saved=collect(CSV.File(joinpath(dest,"nuisance_parameters.csv")))
    coordinates=Float64.(getproperty.(rows,:normalized_parameter))
    issorted(coordinates) || error("Recovered profile is not sorted: $dest")
    others=setdiff(eachindex(names),[index])
    vectors=Vector{Float64}[]
    for point in eachindex(rows)
        entries=filter(v->v.point==point,saved)
        length(entries)==length(names) && Set(String(v.parameter) for v in entries)==Set(string.(names)) || error("Incomplete recovered nuisance vector")
        values=Dict(String(v.parameter)=>Float64(v.value) for v in entries)
        z=[log(values[string(n)]/bounds[string(n)][1])/log(bounds[string(n)][2]/bounds[string(n)][1]) for n in names]
        isapprox(z[index],coordinates[point];atol=1e-10,rtol=0) || error("Recovered coordinate mismatch")
        actual=loglik(z,nothing)-sol.maximum
        isfinite(actual) && isapprox(actual,Float64(rows[point].profile_loglik);atol=1e-5,rtol=0) || error("Recovered likelihood mismatch: $(names[index]) point $point")
        push!(vectors,z)
    end
    threshold=ProfileLikelihood.get_chisq_threshold(.95)
    center=argmin(abs.(coordinates.-sol.mle[index]))
    check=Set([center])
    for direction in (collect(center:-1:1),collect(center:length(rows)))
        found=false
        for (a,b) in zip(direction[1:end-1],direction[2:end])
            if rows[a].profile_loglik>=threshold && rows[b].profile_loglik<=threshold
                union!(check,[a,b]);found=true;break
            end
        end
        found || union!(check,direction[max(1,length(direction)-1):end])
    end
    # Do not carry an unconfirmed optimizer result close to the interval.
    union!(check,[j for j in eachindex(rows) if rows[j].stopping_met!==true && rows[j].profile_loglik>=-3])
    diagnostics=NamedTuple[]
    needs_profile=false
    for point in sort(collect(check))
        fixed_coordinate=coordinates[point]
        objective=function(u,_)
            z=copy(sol.mle);z[index]=fixed_coordinate;z[others].=u
            sol.maximum-loglik(z,nothing)
        end
        f=PackageWorkflow.Optimization.OptimizationFunction(objective,PackageWorkflow.Optimization.AutoFiniteDiff())
        problem=PackageWorkflow.Optimization.OptimizationProblem(f,copy(vectors[point][others]);lb=zeros(length(others)),ub=ones(length(others)))
        starts=FixedTwoSeeds.seed_vectors(bank,fixed_coordinate)
        alg=AuditedBOBYQA(NamedTuple[],ReentrantLock(),_ -> deepcopy(starts))
        trial=SciMLBase.solve(problem,alg;maxiters,abstol=1e-8)
        improvement=-Float64(trial.objective)-Float64(rows[point].profile_loglik)
        stopping_met=SciMLBase.successful_retcode(trial.retcode)
        needs_profile |= improvement>1e-5 || !stopping_met
        push!(diagnostics,(point,parameter=string(names[index]),normalized_parameter=fixed_coordinate,saved_loglik=Float64(rows[point].profile_loglik),rechecked_loglik=-Float64(trial.objective),improvement,stopping_met,retcode=string(trial.retcode)))
        CSV.write(joinpath(dest,"recovery_recheck.csv"),diagnostics)
        println("Rechecked $(names[index]) point $point: likelihood improvement $improvement");flush(stdout)
        -trial.objective>1e-5 && error("Recheck found a better reference likelihood; inspect before reporting intervals")
    end
    if needs_profile
        archive=joinpath(dest,"before_seed_recheck")
        mkpath(archive)
        for file in ("profile.csv","nuisance_parameters.csv","optimizer_calls.csv","interval.csv","COMPLETE")
            if isfile(joinpath(archive,file))
                sha256(read(joinpath(dest,file)))==sha256(read(joinpath(archive,file))) || error("Recovery archive differs: $file")
            else
                cp(joinpath(dest,file),joinpath(archive,file))
            end
        end
        write(joinpath(dest,"RECOMPUTE_REQUIRED"),"Isolated-seed recheck requires a new profile. Original files preserved in before_seed_recheck.\n")
        rm(joinpath(dest,"COMPLETE");force=true)
        rm(joinpath(dest,"RECOVERED_NEEDS_CHECK"))
        println("Recomputing $(names[index]) after seed recheck; saved results preserved in $archive")
        return false
    end
    write(joinpath(dest,"RECOVERY_CHECKED"),"All saved likelihoods verified. Confidence-limit bracket points and any unconfirmed points near the interval were reoptimized using isolated seeds. Review the curve and audit before inference.\n")
    rm(joinpath(dest,"RECOVERED_NEEDS_CHECK"))
    println("REUSED profile $(names[index]): saved likelihoods and interval brackets checked");flush(stdout)
    true
end

function main()
    length(ARGS)>=3 || error("Usage: run_fixed_two.jl REFERENCE_FIT OUTPUT SEED_PROFILE_DIR [RESOLUTION=200] [PARAMETER=all] [--smoke]")
    fitdir,out,seed_dir=abspath.(ARGS[1:3])
    safe_output(out,[fitdir,seed_dir])
    resolution=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 200
    requested=length(ARGS)>=5 ? ARGS[5] : "all"; smoke="--smoke" in ARGS
    resolution>=3 || error("Resolution must be at least 3")
    prep=fixed_two_context(fitdir)
    raw,ctx,model,p,saved_nll,fixed=prep.raw,prep.ctx,prep.model,prep.p,prep.saved_nll,prep.fixed
    id=5
    seedcfg=TOML.parsefile(joinpath(seed_dir,"calibration.toml"))
    seedcfg["bounds"]==raw.config["bounds"] || error("Seed profile bounds differ")
    seedcfg["calibration_cohort"]["series_ids"]==raw.config["calibration_cohort"]["series_ids"] || error("Seed profiles used a different training cohort")
    inactive=false;mle_only=false;conditional_parent="";extra_fixed=""
    names=collect(free_names(ctx,model))
    indices=requested=="all" ? collect(eachindex(names)) : [something(findfirst(==(Symbol(n)),names)) for n in split(requested,",")]
    full_range=get(ENV,"HPSC_PROFILE_FULL_RANGE","0")=="1"
    branch_seeds=false
    start=encode(ctx,model,p)
    function loglik(z,_)
        all(x->isfinite(x)&&-1e-7<=x<=1+1e-7,z)||return -1e100
        try
            -gaussian_nll(ctx.y,prediction(ctx,model,decode(ctx,model,z)),ctx.sigma)
        catch e
            e isa InterruptException && rethrow()
            -1e100
        end
    end
    initial_nll=-loglik(start,nothing)
    settings=Dict("source_fit"=>fitdir,"conditional_parent"=>conditional_parent,"additional_fixed_parameter"=>extra_fixed,"lactate_inhibition"=>!inactive,"mle_only"=>mle_only,"model"=>id,"resolution_per_direction"=>resolution,"parameters"=>string.(names[indices]),"smoke_test"=>smoke,"engine"=>"ProfileLikelihood.jl","package_version"=>string(pkgversion(ProfileLikelihood)),"coordinates"=>"normalized log parameters","mle_optimizer"=>"Fminbox/BFGS/BackTracking","profile_optimizer"=>"NLopt.LN_BOBYQA","next_initial_estimate_method"=>"prev","fixed_parameters"=>fixed,"K_G_bounds"=>ctx.config["bounds"]["K_G"],"K_Lac_bounds"=>ctx.config["bounds"]["K_Lac"])
    settings["publication_model"]=6
    settings["analysis"]= "Nine profiles with K_Lac and Y_GX fixed at saved training estimates"
    settings["training_cultures"]=19
    settings["selection_cultures_not_used"]=12
    settings["full_range_repair"]=full_range
    settings["multistart_seed_directory"]=seed_dir
    settings["cross_branch_seeds"]=branch_seeds
    settings["profile_scope"]="fixed-two"
    files=[joinpath(ROOT,f) for f in ("src/Proliferation.jl","src/Workflow.jl","src/Calibration.jl","src/TargetedChecks.jl","profiles/PackageWorkflow.jl","profiles/run_fixed_two.jl","profiles/FixedTwoSeeds.jl","profiles/Project.toml","profiles/Manifest.toml","data/processed/measurements.csv","data/processed/series.csv","data/processed/feeding.csv","data/processed/lagrangian_exposure.csv")]
    append!(files,[joinpath(fitdir,"calibration.toml"),joinpath(fitdir,"split.csv"),joinpath(fitdir,"observations.csv"),joinpath(fitdir,"model$id/selected_parameters.csv"),joinpath(fitdir,"model$id/summary.csv")])
    !isempty(conditional_parent)&&append!(files,[joinpath(conditional_parent,f) for f in ("settings.toml","calibration.toml","mle_parameters.csv","mle_summary.csv")])
    if !isempty(seed_dir)
        push!(files,joinpath(seed_dir,"calibration.toml"))
        for i in indices, f in ("profile.csv","nuisance_parameters.csv")
            push!(files,joinpath(seed_dir,string(names[i]),f))
        end
    end
    hashes=[(path=abspath(f),sha256=bytes2hex(sha256(read(f)))) for f in files]
    if isdir(out)
        TOML.parsefile(joinpath(out,"settings.toml"))==settings||error("Resume settings differ")
        saved_hashes=collect(CSV.File(joinpath(out,"source_hashes.csv")))
        length(saved_hashes)==length(hashes)&&all(a.path==b.path&&a.sha256==b.sha256 for (a,b) in zip(saved_hashes,hashes))||error("Resume source hashes differ")
    else
        mkpath(out); open(io->TOML.print(io,settings),joinpath(out,"settings.toml"),"w")
        CSV.write(joinpath(out,"source_hashes.csv"),hashes)
    end
    open(io->TOML.print(io,ctx.config),joinpath(out,"calibration.toml"),"w")
    CSV.write(joinpath(out,"fixed_parameters.csv"),[(parameter=n,value=fixed[n],source=inactive&&n=="K_Lac" ? "unused placeholder; inhibition exactly inactive" : n==extra_fixed ? "conditional parent MLE: "*conditional_parent : "saved 100-mL calibration") for n in sort(collect(keys(fixed)))] ; header=["parameter","value","source"])
    prob=make_problem(loglik,start,names)
    checkpoint=joinpath(out,"mle_solution.bin")
    if isfile(checkpoint)
        sol=deserialize(checkpoint)
    elseif isfile(joinpath(out,"recovery.json"))
        # Rebuild the likelihood with current code rather than deserialize a
        # closure from the defective runner. Preserve the saved numerical MLE.
        restored=read_parameters(joinpath(out,"mle_parameters.csv"),model)
        restored_z=encode(ctx,model,restored)
        summary=only(CSV.File(joinpath(out,"mle_summary.csv")))
        restored_maximum=-Float64(summary.nll)
        isapprox(loglik(restored_z,nothing),restored_maximum;atol=1e-6,rtol=0) || error("Recovered reference likelihood differs")
        retcode=getproperty(SciMLBase.ReturnCode,Symbol(summary.retcode))
        alg=PackageWorkflow.OptimizationOptimJL.Fminbox(PackageWorkflow.Optim.BFGS(linesearch=PackageWorkflow.LineSearches.BackTracking()))
        sol=ProfileLikelihood.LikelihoodSolution{length(names),typeof(restored_z),typeof(prob),typeof(restored_maximum),typeof(retcode),typeof(alg)}(restored_z,prob,alg,restored_maximum,retcode)
        serialize(checkpoint,sol)
        println("Restored saved numerical MLE; no reference refit");flush(stdout)
    else
        println("Model 6: nine-parameter profile reference optimization on the 19 training cultures");flush(stdout)
        sol=fit_mle(prob;maxiters=smoke ? 100 : 5000)
        nll=-loglik(sol.mle,nothing)
        isfinite(nll)&&nll<=initial_nll+1e-5||error("Package MLE is worse than its starting point; inspect before profiling: $nll vs $initial_nll")
        isapprox(nll,-sol.maximum;atol=1e-6,rtol=0)||error("Package MLE objective mismatch")
        serialize(checkpoint,sol)
        save_parameters(joinpath(out,"mle_parameters.csv"),decode(ctx,model,sol.mle))
        CSV.write(joinpath(out,"mle_summary.csv"),[(model=id,saved_nll,initial_nll,nll,free_parameters=length(names),retcode=string(sol.retcode),stopping_met=SciMLBase.successful_retcode(sol.retcode),smoke_test=smoke)])
    end
    if mle_only
        write(joinpath(out,"MLE_COMPLETE"),"MLE-only conditional sensitivity; not a profile suite\n")
        println("Completed conditional MLE: $out")
        return
    end
    # Resume uses the package's saved MLE; newly constructed likelihood is checked above.
    isapprox(-loglik(sol.mle,nothing),-sol.maximum;atol=1e-6,rtol=0) || error("Resumed MLE differs")
    prob=sol.problem
    physical(i,z)=exp(log(ctx.config["bounds"][string(names[i])][1])+z*log(ctx.config["bounds"][string(names[i])][2]/ctx.config["bounds"][string(names[i])][1]))
    Threads.@threads :greedy for i in indices
        dest=joinpath(out,string(names[i]))
        isfile(joinpath(dest,"COMPLETE")) && !isfile(joinpath(dest,"RECOVERED_NEEDS_CHECK")) && !isfile(joinpath(dest,"RECOMPUTE_REQUIRED")) && continue
        bank=load_seed_bank(seed_dir,names,i,ctx.config["bounds"],sol.mle)
        if isfile(joinpath(dest,"RECOVERED_NEEDS_CHECK")) && !isfile(joinpath(dest,"RECOMPUTE_REQUIRED"))
            println("Checking recovered profile $(names[i])");flush(stdout)
            recheck_saved_profile(dest,names,i,ctx.config["bounds"],sol,loglik,bank;maxiters=smoke ? 200 : 10000) && continue
        end
        mkpath(dest); println("Profiling $(names[i]) with ProfileLikelihood.jl");flush(stdout)
        prof,records=profile_parameter(prob,sol,i;resolution,maxiters=smoke ? 200 : 10000,parallel=true,full_range,seed_provider=FixedTwoSeeds.seed_provider(bank))
        rows=NamedTuple[]; nuisance=NamedTuple[]
        for j in eachindex(prof.parameter_values[i])
            z=copy(sol.mle);z[i]=prof.parameter_values[i][j]
            others=prof.other_mles[i][j];z[setdiff(eachindex(z),[i])].=others
            actual=loglik(z,nothing)-sol.maximum
            isapprox(actual,prof.profile_values[i][j];atol=1e-5,rtol=0)||error("Profile objective verification failed")
            matches=filter(r->length(r.u)==length(others)&&maximum(abs.(r.u.-others))<1e-12&&abs(r.objective+actual)<1e-6,records)
            at_mle=maximum(abs.(z.-sol.mle))<1e-12
            stop=isempty(matches) ? (at_mle ? SciMLBase.successful_retcode(sol.retcode) : missing) : last(matches).stopping_met
            code=isempty(matches) ? (at_mle ? string(sol.retcode) : "not recorded") : last(matches).retcode
            push!(rows,(parameter=string(names[i]),parameter_value=physical(i,z[i]),profile_loglik=actual,normalized_parameter=z[i],stopping_met=stop,retcode=code))
            append!(nuisance,[(point=j,parameter=string(names[k]),value=physical(k,z[k])) for k in eachindex(z)])
        end
        CSV.write(joinpath(dest,"profile.csv"),rows);CSV.write(joinpath(dest,"nuisance_parameters.csv"),nuisance)
        CSV.write(joinpath(dest,"optimizer_calls.csv"),[(call=j,objective=r.objective,retcode=r.retcode,stopping_met=r.stopping_met,starts=r.starts,best_seed=r.best_seed,start_objective_spread=r.start_objective_spread) for (j,r) in enumerate(records)])
        ci=prof.confidence_intervals[i]
        CSV.write(joinpath(dest,"interval.csv"),[(parameter=string(names[i]),mle=physical(i,sol.mle[i]),package_lower=physical(i,ci.lower),package_upper=physical(i,ci.upper),lower_hits_bound=ci.lower<=1e-7,upper_hits_bound=ci.upper>=1-1e-7,unconfirmed_points=count(r->r.stopping_met!==true,rows),better_minimum=maximum(r.profile_loglik for r in rows)>1e-5)])
        write(joinpath(dest,"COMPLETE"),"ProfileLikelihood.jl profile complete; review interval qualifications.\n")
        rm(joinpath(dest,"RECOMPUTE_REQUIRED");force=true)
        rm(joinpath(dest,"RECOVERED_NEEDS_CHECK");force=true)
        println("FINISHED profile $(names[i])"); flush(stdout)
    end
    CSV.write(joinpath(out,"profile_confidence_intervals.csv"),vcat([collect(CSV.File(joinpath(out,string(names[i]),"interval.csv"))) for i in indices]...))
    write(joinpath(out,"COMPLETE"),smoke ? "SMOKE TEST ONLY; not final inference\n" : "ProfileLikelihood.jl suite complete; review flags before inference\n")
    println("Completed package workflow: $out")
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
