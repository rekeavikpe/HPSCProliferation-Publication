include(joinpath(@__DIR__,"..","profiles","run_fixed_two.jl"))

function review()
    base=joinpath(@__DIR__,"verification")
    fit=joinpath(@__DIR__,"cache/training60/reference_fit")
    root=joinpath(@__DIR__,"cache/nine_profiles")
    prep=fixed_two_context(fit)
    ctx,model=prep.ctx,prep.model
    summary=only(CSV.File(joinpath(root,"mle_summary.csv")))
    estimate=read_parameters(joinpath(root,"mle_parameters.csv"),model)
    nll=gaussian_nll(ctx.y,prediction(ctx,model,estimate),ctx.sigma)
    @assert isapprox(nll,summary.nll;atol=1e-6,rtol=0)
    println("Reference NLL reproduced: ", nll, "; training cultures: ", length(ctx.config["calibration_cohort"]["series_ids"]), "; observations: ", length(ctx.y));flush(stdout)
    details=NamedTuple[]
    for name in ["K_Agg","K_d","K_G","K_oG","Y_GLac","alpha","mu_max","K_tau","K_eps"]
        points=collect(CSV.File(joinpath(root,name,"profile.csv")))
        nuisance=collect(CSV.File(joinpath(root,name,"nuisance_parameters.csv")))
        errors=Float64[]
        for (i,row) in enumerate(points)
            entries=filter(r->r.point==i,nuisance)
            @assert length(entries)==9
            values=Dict(Symbol(r.parameter)=>Float64(r.value) for r in entries)
            for (key,val) in prep.fixed
                values[Symbol(key)]=val
            end
            p=parameters(model,[values[k] for k in parameter_names(model)])
            actual=summary.nll-gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
            err=actual-row.profile_loglik
            push!(errors,abs(err))
            push!(details,(parameter=name,point=i,parameter_value=Float64(row.parameter_value),saved_loglik=Float64(row.profile_loglik),replayed_loglik=actual,error=err))
        end
        println(name," reproduced ",length(points)," points; max |difference| = ",maximum(errors));flush(stdout)
        @assert maximum(errors)<1e-5
    end
    CSV.write(joinpath(base,"nine_profile_replay.csv"),details)

end
review()
