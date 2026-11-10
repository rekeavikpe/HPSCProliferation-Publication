module PackageWorkflow
using ProfileLikelihood, Optimization, OptimizationOptimJL, OptimizationNLopt, NLopt, Optim, LineSearches, FiniteDiff, SciMLBase
export AuditedBOBYQA, make_problem, fit_mle, profile_parameter, inset!

inset!(x,lo,hi) = (x .= clamp.(x,lo .+ 1e-9 .* (hi.-lo),hi .- 1e-9 .* (hi.-lo)); x)

# Retain the bound-safe :prev restart used in the original all-11 script.
function ProfileLikelihood.set_next_initial_estimate!(sub_cache, param_vals, other_mles, prob, theta_n; next_initial_estimate_method=Val(:prev))
    next_initial_estimate_method==Val(:prev)||error("This compatibility fix supports :prev only")
    sub_cache .= other_mles[end]
    inset!(sub_cache,prob.lb,prob.ub)
    nothing
end

struct AuditedBOBYQA
    records::Vector{NamedTuple}
    lock::ReentrantLock
    seeds::Function
end
AuditedBOBYQA()=AuditedBOBYQA(NamedTuple[],ReentrantLock(), prob->Vector{Float64}[])
function SciMLBase.solve(prob::SciMLBase.OptimizationProblem, alg::AuditedBOBYQA; kwargs...)
    starts=Vector{Float64}[copy(prob.u0)]
    append!(starts,alg.seeds(prob))
    solutions=Any[]
    best_seed=minimum(prob.f(u,prob.p) for u in starts)
    for start in starts
        inset!(start,prob.lb,prob.ub)
        trial=Optimization.solve(SciMLBase.remake(prob;u0=start),NLopt.LN_BOBYQA();kwargs...)
        isfinite(trial.objective)&&push!(solutions,trial)
    end
    isempty(solutions)&&error("All nuisance optimizations failed")
    sol=solutions[argmin([x.objective for x in solutions])]
    # Different local algorithm checks the best derivative-free solution near the CI.
    if length(starts)>1 && sol.objective<20
        trial=Optimization.solve(SciMLBase.remake(prob;u0=copy(sol.u)),
            OptimizationOptimJL.Fminbox(Optim.BFGS(linesearch=LineSearches.BackTracking()));
            maxiters=1000,g_tol=1e-7)
        isfinite(trial.objective)&&trial.objective<sol.objective&&(sol=trial)
    end
    sol.objective<=best_seed+1e-6||error("Nuisance fit is worse than a supplied starting solution")
    lock(alg.lock) do
        push!(alg.records,(u=copy(sol.u),objective=Float64(sol.objective),retcode=string(sol.retcode),stopping_met=SciMLBase.successful_retcode(sol.retcode),starts=length(starts),best_seed=Float64(best_seed),start_objective_spread=maximum(x.objective for x in solutions)-minimum(x.objective for x in solutions)))
        if get(ENV,"HPSC_PROFILE_POINT_PROGRESS","0")=="1" && length(alg.records)%20==0
            println("Completed nuisance optimizations: $(length(alg.records)); latest relative objective: $(sol.objective)")
            flush(stdout)
        end
    end
    sol
end
function make_problem(loglik,start,names)
    z=Float64.(start); inset!(z,zeros(length(z)),ones(length(z)))
    LikelihoodProblem(loglik,z;data=nothing,syms=names,
        f_kwargs=(adtype=Optimization.AutoFiniteDiff(),),prob_kwargs=(lb=zeros(length(z)),ub=ones(length(z))))
end
function fit_mle(prob;maxiters=5000)
    mle(prob,OptimizationOptimJL.Fminbox(Optim.BFGS(linesearch=LineSearches.BackTracking()));maxiters,g_tol=1e-6)
end
function profile_parameter(prob,sol,index;resolution=100,maxiters=10000,parallel=true,full_range=false,seed_provider=nothing)
    ranges=construct_profile_ranges(sol,zeros(length(sol.mle)),ones(length(sol.mle)),fill(resolution,length(sol.mle)))
    alg=seed_provider===nothing ? AuditedBOBYQA() : AuditedBOBYQA(NamedTuple[],ReentrantLock(),seed_provider)
    prof=ProfileLikelihood.profile(prob,sol,[index];param_ranges=ranges,parallel,conf_level=.95,
        alg,next_initial_estimate_method=:prev,progress=false,maxiters,abstol=1e-8,
        min_steps=full_range ? 0 : min(10,resolution),
        threshold=full_range ? -Inf : ProfileLikelihood.get_chisq_threshold(.95),
        confidence_interval_method=full_range ? :extrema : :spline)
    if full_range
        x=prof.parameter_values[index];y=prof.profile_values[index]
        threshold=ProfileLikelihood.get_chisq_threshold(.95)
        mid=argmin(abs.(x.-sol.mle[index]))
        function endpoint(order,bound)
            for (a,b) in zip(order[1:end-1],order[2:end])
                if y[a]>=threshold && y[b]<=threshold
                    lo,hi=minmax(x[a],x[b]); f=z->prof.splines[index](z)-threshold
                    flo=f(lo)
                    for _ in 1:60
                        m=(lo+hi)/2
                        if (f(m)>=0)==(flo>=0);lo=m;else;hi=m;end
                    end
                    return (lo+hi)/2
                end
            end
            bound
        end
        lower=endpoint(collect(mid:-1:1),first(x))
        upper=endpoint(collect(mid:length(x)),last(x))
        prof.confidence_intervals[index]=ProfileLikelihood.ConfidenceInterval(lower,upper,.95)
    end
    prof,alg.records
end
end
