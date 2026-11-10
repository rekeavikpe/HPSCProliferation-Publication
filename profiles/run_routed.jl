using CSV, TOML, SHA, Serialization, ProfileLikelihood, SciMLBase
include("PackageWorkflow.jl")
using .PackageWorkflow
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
include("../ablations/RoutedExposure.jl")
using .RoutedExposure

function main()
    length(ARGS)>=3||error("Usage: run_routed.jl SCENARIO_FIT OUTPUT all-active|conditional [RESOLUTION=100] [PARAMETER=all] [--smoke]")
    fitdir,out=abspath.(ARGS[1:2]);mode=ARGS[3]
    mode in ("all-active","conditional")||error("Unknown profile scope")
    resolution=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 100
    requested=length(ARGS)>=5 ? ARGS[5] : "all";smoke="--smoke" in ARGS
    resolution>=3||error("Resolution must be at least 3")
    routing=TOML.parsefile(joinpath(fitdir,"routing.toml"));scenario=routing["scenario"];shared=get(routing,"shared_response_parameter",false)
    id=scenario=="baseline" ? 0 : 5;model=legacy_model(id)
    raw=checked_context(fitdir)
    ctx=scenario=="baseline" ? raw : routed_context(raw,scenario;shared)
    p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
    saved_nll=only(CSV.File(joinpath(fitdir,"model$id/summary.csv"))).nll
    source_nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
    isapprox(source_nll,saved_nll;atol=1e-6,rtol=0)||error("Routed calibration likelihood differs")
    fixed=mode=="conditional" ? Dict("K_G"=>p.K_G,"K_Lac"=>p.K_Lac,"Y_GX"=>p.Y_GX) : Dict{String,Float64}()
    if !isempty(fixed)
        cfg=deepcopy(raw.config);merge!(cfg["fixed"],fixed)
        raw=context(cfg=cfg);ctx=scenario=="baseline" ? raw : routed_context(raw,scenario;shared)
    end
    inactive=false;mle_only=false;conditional_parent="";extra_fixed=""
    names=collect(free_names(ctx,model))
    indices=requested=="all" ? collect(eachindex(names)) : [something(findfirst(==(Symbol(n)),names)) for n in split(requested,",")]
    full_range=get(ENV,"HPSC_PROFILE_FULL_RANGE","0")=="1"
    seed_dir=get(ENV,"HPSC_PROFILE_SEED_DIR","")
    branch_seeds=get(ENV,"HPSC_PROFILE_BRANCH_SEEDS","0")=="1"
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
    settings["routing_scenario"]=scenario
    settings["full_range_repair"]=full_range
    settings["multistart_seed_directory"]=seed_dir
    settings["cross_branch_seeds"]=branch_seeds
    settings["parameter_latex"]=response_labels(scenario;shared)
    settings["profile_scope"]=mode
    files=[joinpath(ROOT,f) for f in ("src/Proliferation.jl","src/Workflow.jl","src/Calibration.jl","src/TargetedChecks.jl","profiles/PackageWorkflow.jl","profiles/run_routed.jl","ablations/RoutedExposure.jl","profiles/Project.toml","profiles/Manifest.toml","data/processed/measurements.csv","data/processed/series.csv","data/processed/feeding.csv","data/processed/lagrangian_exposure.csv")]
    append!(files,[joinpath(fitdir,"calibration.toml"),joinpath(fitdir,"routing.toml"),joinpath(fitdir,"observations.csv"),joinpath(fitdir,"model$id/selected_parameters.csv"),joinpath(fitdir,"model$id/summary.csv")])
    !isempty(conditional_parent)&&append!(files,[joinpath(conditional_parent,f) for f in ("settings.toml","calibration.toml","mle_parameters.csv","mle_summary.csv")])
    if !isempty(seed_dir)
        for i in indices, f in ("profile.csv","nuisance_parameters.csv")
            push!(files,joinpath(seed_dir,string(names[i]),f))
        end
    end
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
    prob=sol.problem
    physical(i,z)=exp(log(ctx.config["bounds"][string(names[i])][1])+z*log(ctx.config["bounds"][string(names[i])][2]/ctx.config["bounds"][string(names[i])][1]))
    Threads.@threads :greedy for i in indices
        dest=joinpath(out,string(names[i]));isfile(joinpath(dest,"COMPLETE"))&&continue
        mkpath(dest); println("Profiling $(names[i]) with ProfileLikelihood.jl");flush(stdout)
        seed_provider=nothing
        if !isempty(seed_dir)
            old=collect(CSV.File(joinpath(seed_dir,string(names[i]),"profile.csv")))
            bank=collect(CSV.File(joinpath(seed_dir,string(names[i]),"nuisance_parameters.csv")))
            other_indices=setdiff(eachindex(names),[i])
            vectors=Dict{Int,Vector{Float64}}()
            for j in eachindex(old)
                values=Dict(String(v.parameter)=>Float64(v.value) for v in bank if v.point==j)
                vectors[j]=[clamp(log(values[string(names[k])]/ctx.config["bounds"][string(names[k])][1])/log(ctx.config["bounds"][string(names[k])][2]/ctx.config["bounds"][string(names[k])][1]),0.,1.) for k in other_indices]
            end
            seed_provider=function(fixed_prob)
                hasproperty(fixed_prob.f.f,Symbol("θₙ"))||error("ProfileLikelihood fixed-coordinate interface changed")
                fixed=Float64(getproperty(fixed_prob.f.f,Symbol("θₙ")))
                nearest=sortperm([abs(Float64(v.normalized_parameter)-fixed) for v in old])[1:min(3,length(old))]
                if branch_seeds
                    # Include both sides of a branch change, not only the nearest points.
                    plausible=findall(v->Float64(v.profile_loglik)>-10,old)
                    anchors=isempty(plausible) ? Int[] : plausible[unique(round.(Int,range(1,length(plausible);length=min(16,length(plausible)))))]
                    nearest=unique(vcat(nearest,anchors))
                end
                [copy(sol.mle[other_indices]),[copy(vectors[j]) for j in nearest]...]
            end
        end
        prof,records=profile_parameter(prob,sol,i;resolution,maxiters=smoke ? 200 : 10000,parallel=true,full_range,seed_provider)
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
    end
    CSV.write(joinpath(out,"profile_confidence_intervals.csv"),vcat([collect(CSV.File(joinpath(out,string(names[i]),"interval.csv"))) for i in indices]...))
    write(joinpath(out,"COMPLETE"),smoke ? "SMOKE TEST ONLY; not final inference\n" : "ProfileLikelihood.jl suite complete; review flags before inference\n")
    println("Completed package workflow: $out")
end
main()
