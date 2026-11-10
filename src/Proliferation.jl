module Proliferation

using SciMLBase: ODEProblem, solve, successful_retcode
using OrdinaryDiffEqTsit5: Tsit5

export model_rhs!, parameter_names, parameters, simulate, exchange,
       linear_exposure, gaussian_nll, bic, UNIQUE_MODELS

const UNIQUE_MODELS = (0, 1, 3, 4, 5, 6, 7)
const COMMON = (:K_Agg, :K_d, :K_G, :K_Lac, :K_oG, :Y_GLac, :Y_GX, :alpha, :mu_max)

function parameter_names(model)
    model in 0:7 || throw(ArgumentError("Model must be in 0:7"))
    model == 0 && return COMMON
    model == 6 && return (COMMON..., :K_tau, :K_eps)
    model == 7 && return (COMMON..., :K_tau, :delta_max, :K_eps, :gamma)
    return (COMMON..., :K_tau, :delta_max, :K_eps)
end

function parameters(model, values)
    names = parameter_names(model)
    length(values) == length(names) || throw(ArgumentError("Parameter count mismatch"))
    all(x -> isfinite(x) && x >= 0, values) || throw(ArgumentError("Invalid parameter value"))
    p = NamedTuple{names}(Tuple(Float64.(values)))
    for name in (:K_Agg, :K_G, :K_Lac, :K_tau, :K_eps)
        hasproperty(p, name) && getproperty(p, name) <= 0 && throw(ArgumentError("$name must be positive"))
    end
    return p
end

"State order: cells/mL, glucose mM, lactate mM, aggregate diameter µm; time in hours."
function model_rhs!(du, u, p, t, model, tau=0.0, edr=0.0; death_time_h=12.0, lactate_inhibition=true)
    # Retain the archived evaluation floors; they are not observation replacements.
    X, G, Lac, Agg = max(u[1], 1.0), max(u[2], .001), max(u[3], 0.0), max(u[4], 1.0)
    inhibition = lactate_inhibition ? p.K_Lac/(Lac+p.K_Lac) : 1.0
    base = p.mu_max * G/(G+p.K_G) * inhibition * p.K_Agg/(Agg^3+p.K_Agg)
    shear = if model == 0
        1.0
    elseif model == 3
        x = tau/p.K_tau
        1/(1+x+x^2/2)
    elseif model == 5
        p.K_tau^2/(tau^2+p.K_tau^2)
    elseif model == 7
        k = p.K_tau*(1+p.gamma*edr/(edr+p.K_eps))
        k/(tau+k)
    else
        p.K_tau/(tau+p.K_tau) # Model 2 is exactly the Model 1 hypothesis.
    end
    growth = base*shear
    metabolic_growth = model == 4 ? base : growth
    lactate_rate = p.K_oG*metabolic_growth*G
    glucose_rate = p.Y_GX*metabolic_growth+p.Y_GLac*lactate_rate
    du[1] = (growth-p.K_d*exp(-t/death_time_h))*X
    du[2] = -glucose_rate*X
    du[3] = lactate_rate*X
    du[4] = if model == 0
        p.alpha*base
    elseif model == 6
        growth*(p.alpha-edr/p.K_eps*Agg)
    else
        p.alpha*metabolic_growth-p.delta_max*edr/(edr+p.K_eps)*Agg
    end
    return nothing
end

"Medium exchange retains cells and aggregates and mixes extracellular metabolites."
function exchange(u, fraction, fresh_glucose; fresh_lactate=0.0)
    isfinite(fraction) && 0 <= fraction <= 1 || throw(ArgumentError("Invalid exchange fraction"))
    all(x -> isfinite(x) && x >= 0, (fresh_glucose, fresh_lactate)) || throw(ArgumentError("Invalid fresh medium"))
    v = copy(u)
    v[2] = (1-fraction)*u[2]+fraction*fresh_glucose
    v[3] = (1-fraction)*u[3]+fraction*fresh_lactate
    return v
end

"Linear size interpolation with explicit nearest-endpoint behavior outside the CFD grid."
function linear_exposure(sizes, shear, edr)
    length(sizes) == length(shear) == length(edr) >= 2 || throw(ArgumentError("Invalid exposure grid"))
    issorted(sizes) && allunique(sizes) || throw(ArgumentError("Sizes must increase strictly"))
    all(x -> isfinite(x) && x >= 0, vcat(sizes, shear, edr)) || throw(ArgumentError("Invalid exposure values"))
    return function(size)
        isfinite(size) || throw(ArgumentError("Nonfinite aggregate size"))
        x = clamp(size, first(sizes), last(sizes))
        i = clamp(searchsortedlast(sizes, x), 1, length(sizes)-1)
        weight = (x-sizes[i])/(sizes[i+1]-sizes[i])
        return ((1-weight)*shear[i]+weight*shear[i+1], (1-weight)*edr[i]+weight*edr[i+1])
    end
end

"Integrate between exchange events; record observations at the pre-exchange endpoint."
function simulate(model, p, u0, sample_times, feeds; fresh_glucose,
                  exposure=agg -> (0.0, 0.0), death_time_h=12.0, lactate_inhibition=true,
                  reltol=1e-10, abstol=1e-10, algorithm=Tsit5())
    keys(p) == parameter_names(model) || throw(ArgumentError("Parameter names/order mismatch"))
    parameters(model, collect(values(p)))
    length(u0) == 4 && all(x -> isfinite(x) && x >= 0, u0) || throw(ArgumentError("Invalid initial state"))
    isfinite(death_time_h) && death_time_h > 0 || throw(ArgumentError("Invalid death time"))
    !isempty(sample_times) && issorted(sample_times) && allunique(sample_times) &&
        all(x -> isfinite(x) && x >= 0, sample_times) || throw(ArgumentError("Invalid sample times"))
    times = [f.time_h for f in feeds]
    issorted(times) && allunique(times) && all(t -> 0 < t <= last(sample_times), times) || throw(ArgumentError("Invalid feed times"))
    for f in feeds
        exchange(Float64.(u0), f.fraction, fresh_glucose)
    end
    function rhs!(du, u, unused, t)
        tau, edr = model == 0 ? (0.0, 0.0) : exposure(u[4])
        all(x -> isfinite(x) && x >= 0, (tau, edr)) || throw(ArgumentError("Invalid exposure"))
        model_rhs!(du, u, p, t, model, tau, edr; death_time_h, lactate_inhibition)
    end
    predicted = Matrix{Float64}(undef, 4, length(sample_times))
    u, start = Float64.(u0), 0.0
    events = NamedTuple[]
    for stop in sort(unique(vcat(Float64.(times), Float64(last(sample_times)))))
        if stop > start
            solution = solve(ODEProblem(rhs!, u, (start, stop)), algorithm;
                             reltol, abstol, dense=true, maxiters=1_000_000)
            successful_retcode(solution) || error("Integration failed: $(solution.retcode)")
            all(s -> all(isfinite, s), solution.u) || error("Nonfinite solution")
            for (j, t) in enumerate(sample_times)
                if (start == 0 && t == 0) || start < t <= stop
                    predicted[:, j] = solution(t)
                end
            end
            u = copy(solution.u[end])
        else
            predicted[:, 1] = u
        end
        index = findfirst(==(stop), times)
        if index !== nothing
            before = copy(u)
            u = exchange(u, feeds[index].fraction, fresh_glucose)
            push!(events, (time_h=stop, fraction=feeds[index].fraction, before, after=copy(u)))
        end
        start = stop
    end
    return (predictions=predicted, events=events)
end

"Fixed-sigma Gaussian score. Missing observations are filtered by the caller, never imputed."
function gaussian_nll(observed, predicted, sigma)
    length(observed) == length(predicted) == length(sigma) > 0 || throw(ArgumentError("Score length mismatch"))
    all(isfinite, observed) && all(isfinite, predicted) || throw(ArgumentError("Nonfinite observation/prediction"))
    all(s -> isfinite(s) && s > 0, sigma) || throw(ArgumentError("Sigma must be positive"))
    return sum(log(s)+log(2pi)/2+((y-f)/s)^2/2 for (y,f,s) in zip(observed,predicted,sigma))
end

function bic(nll, parameter_count::Int, observation_count::Int)
    isfinite(nll) && parameter_count >= 0 && observation_count > 0 || throw(ArgumentError("Invalid BIC input"))
    return 2nll+parameter_count*log(observation_count)
end

end
