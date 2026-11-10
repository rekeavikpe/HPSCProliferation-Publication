using CSV, TOML, SHA, Serialization, ProfileLikelihood, SciMLBase
include("PackageWorkflow.jl")
using .PackageWorkflow
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation

function main()
    length(ARGS)>=2||error("Usage: profiles/run.jl FIT OUT [MODEL=5] [RESOLUTION=100] [all|PARAMETER] [--smoke] [--fix-lactate-yield | --fix-yield | --fix-three] [--inactive-lactate] [--mle-only] [--kg-upper=40] [--lac-upper=100] [--conditional-from=DIR --fix-at-mle=PARAMETER]")
    fitdir,out=abspath.(ARGS[1:2]); id=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 5
    resolution=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 100
    requested=length(ARGS)>=5 ? ARGS[5] : "all"; smoke="--smoke" in ARGS
    resolution>=3||error("Resolution must be at least 3")
    ctx=checked_context(fitdir); model=legacy_model(id)
    p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
    saved_nll=only(CSV.File(joinpath(fitdir,"model$id/summary.csv"))).nll
    source_nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
    isapprox(source_nll,saved_nll;atol=1e-6,rtol=0)||error("Saved calibration likelihood differs")
    inactive="--inactive-lactate" in ARGS
    mle_only="--mle-only" in ARGS
    inactive && "--fix-lactate-yield" in ARGS && error("Use --inactive-lactate --fix-yield; inactive K_Lac is not an estimated fixed value")
    fixed_three="--fix-three" in ARGS
    fixed_three && any(a->a in ARGS,("--inactive-lactate","--fix-yield","--fix-lactate-yield")) && error("--fix-three retains active lactate inhibition; use it without other fixation flags")
    fixed=Dict{String,Float64}()
    if fixed_three
        id==5||error("The fixed-three conditional analysis is specified for Model 5")
        get(ctx.config,"lactate_inhibition",true)||error("Expected source calibration with active lactate inhibition")
        fixed=Dict("K_G"=>p.K_G,"K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX)
        cfg=deepcopy(ctx.config);cfg["lactate_inhibition"]=true
        merge!(cfg["fixed"],fixed);ctx=context(cfg=cfg)
    end
    if "--fix-lactate-yield" in ARGS
        id==5||error("The conditional two-parameter analysis is specified for Model 5")
        fixed=Dict("K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX)
        cfg=deepcopy(ctx.config);merge!(cfg["fixed"],fixed);ctx=context(cfg=cfg)
    end
    if "--fix-yield" in ARGS || inactive
        id==5||error("This sensitivity analysis is specified for Model 5")
        cfg=deepcopy(ctx.config)
        if "--fix-yield" in ARGS
            fixed["Y_GX"]=p.Y_GX
        end
        if inactive
            cfg["lactate_inhibition"]=false
            fixed["K_Lac"]=1.0 # Unused placeholder: inhibition is exactly one in the equations.
        end
        merge!(cfg["fixed"],fixed);ctx=context(cfg=cfg)
    end
    parent_options=filter(a->startswith(a,"--conditional-from="),ARGS)
    fixed_options=filter(a->startswith(a,"--fix-at-mle="),ARGS)
    length(parent_options)<=1&&length(fixed_options)<=1||error("Specify one parent and one extra fixed parameter")
    conditional_parent="";extra_fixed=""
    if !isempty(parent_options)||!isempty(fixed_options)
        length(parent_options)==length(fixed_options)==1||error("Use --conditional-from and --fix-at-mle together")
        conditional_parent=abspath(split(only(parent_options),"=";limit=2)[2])
        extra_fixed=split(only(fixed_options),"=";limit=2)[2]
        isfile(joinpath(conditional_parent,"COMPLETE"))||error("Parent profile suite is incomplete")
        pm=TOML.parsefile(joinpath(conditional_parent,"settings.toml"))
        pm["source_fit"]==fitdir&&pm["model"]==id&&!pm["smoke_test"]||error("Conditional parent source or model differs")
        TOML.parsefile(joinpath(conditional_parent,"calibration.toml"))==ctx.config||error("Conditional parent configuration differs")
        Symbol(extra_fixed) in free_names(ctx,model)||error("Additional fixed parameter must be free in parent")
        p=read_parameters(joinpath(conditional_parent,"mle_parameters.csv"),model)
        parent_summary=only(CSV.File(joinpath(conditional_parent,"mle_summary.csv")))
        parent_summary.stopping_met||error("Parent MLE did not converge")
        isapprox(gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma),parent_summary.nll;atol=1e-6,rtol=0)||error("Parent MLE does not reproduce")
        fixed[extra_fixed]=getproperty(p,Symbol(extra_fixed))
        cfg=deepcopy(ctx.config);merge!(cfg["fixed"],fixed);ctx=context(cfg=cfg)
    end
    kg_options=filter(a->startswith(a,"--kg-upper="),ARGS)
    length(kg_options)<=1||error("Specify K_G upper bound only once")
    if !isempty(kg_options)
        upper=parse(Float64,split(only(kg_options),"=";limit=2)[2])
        lower=ctx.config["bounds"]["K_G"][1]
        isfinite(upper)&&lower<upper||error("Invalid K_G upper bound")
        lower<=p.K_G<=upper||error("Saved MLE is outside the requested K_G range; a separate initialization is required")
        cfg=deepcopy(ctx.config);cfg["bounds"]["K_G"]=[lower,upper];ctx=context(cfg=cfg)
    end
    lac_options=filter(a->startswith(a,"--lac-upper="),ARGS)
    length(lac_options)<=1||error("Specify K_Lac upper bound only once")
    if !isempty(lac_options)
        isempty(fixed)&&isempty(conditional_parent)||error("A new lactate cap must first be calibrated with parameters free")
        upper=parse(Float64,split(only(lac_options),"=";limit=2)[2]);lower=ctx.config["bounds"]["K_Lac"][1]
        isfinite(upper)&&lower<upper||error("Invalid K_Lac upper bound")
        cfg=deepcopy(ctx.config);cfg["bounds"]["K_Lac"]=[lower,upper];ctx=context(cfg=cfg)
        p=merge(p,(K_Lac=clamp(p.K_Lac,lower,upper),))
    end
    names=collect(free_names(ctx,model))
    indices=requested=="all" ? collect(eachindex(names)) : [something(findfirst(==(Symbol(requested)),names))]
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
    files=[joinpath(ROOT,f) for f in ("src/Proliferation.jl","src/Workflow.jl","src/Calibration.jl","src/TargetedChecks.jl","profiles/PackageWorkflow.jl","profiles/run.jl","profiles/Project.toml","profiles/Manifest.toml","data/processed/measurements.csv","data/processed/series.csv","data/processed/feeding.csv","data/processed/lagrangian_exposure.csv")]
    append!(files,[joinpath(fitdir,"calibration.toml"),joinpath(fitdir,"observations.csv"),joinpath(fitdir,"model$id/selected_parameters.csv"),joinpath(fitdir,"model$id/summary.csv")])
    !isempty(conditional_parent)&&append!(files,[joinpath(conditional_parent,f) for f in ("settings.toml","calibration.toml","mle_parameters.csv","mle_summary.csv")])
    hashes=[(path=abspath(f),sha256=bytes2hex(sha256(read(f)))) for f in files]
    if isdir(out)
        TOML.parsefile(joinpath(out,"settings.toml"))==settings||error("Resume settings differ")
        old=collect(CSV.File(joinpath(out,"source_hashes.csv")))
        length(old)==length(hashes)&&all(a.path==b.path&&a.sha256==b.sha256 for (a,b) in zip(old,hashes))||error("Resume source hashes differ")
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
    else
        println("ProfileLikelihood.jl MLE: Model $id");flush(stdout)
        sol=fit_mle(prob;maxiters=smoke ? 30 : 5000)
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
    prob=sol.problem
    physical(i,z)=exp(log(ctx.config["bounds"][string(names[i])][1])+z*log(ctx.config["bounds"][string(names[i])][2]/ctx.config["bounds"][string(names[i])][1]))
    Threads.@threads :greedy for i in indices
        dest=joinpath(out,string(names[i]));isfile(joinpath(dest,"COMPLETE"))&&continue
        mkpath(dest); println("Profiling $(names[i]) with ProfileLikelihood.jl");flush(stdout)
        prof,records=profile_parameter(prob,sol,i;resolution,maxiters=smoke ? 200 : 10000,parallel=true)
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
        CSV.write(joinpath(dest,"optimizer_calls.csv"),[(call=j,objective=r.objective,retcode=r.retcode,stopping_met=r.stopping_met) for (j,r) in enumerate(records)])
        ci=prof.confidence_intervals[i]
        CSV.write(joinpath(dest,"interval.csv"),[(parameter=string(names[i]),mle=physical(i,sol.mle[i]),package_lower=physical(i,ci.lower),package_upper=physical(i,ci.upper),lower_hits_bound=ci.lower<=1e-7,upper_hits_bound=ci.upper>=1-1e-7,unconfirmed_points=count(r->r.stopping_met!==true,rows),better_minimum=maximum(r.profile_loglik for r in rows)>1e-5)])
        write(joinpath(dest,"COMPLETE"),"ProfileLikelihood.jl profile complete; review interval qualifications.\n")
    end
    CSV.write(joinpath(out,"profile_confidence_intervals.csv"),vcat([collect(CSV.File(joinpath(out,string(names[i]),"interval.csv"))) for i in indices]...))
    write(joinpath(out,"COMPLETE"),smoke ? "SMOKE TEST ONLY; not final inference\n" : "ProfileLikelihood.jl suite complete; review flags before inference\n")
    println("Completed package workflow: $out")
end
main()
