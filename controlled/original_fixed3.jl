# =============================================================================
# PROFILE LIKELIHOOD ANALYSIS — MODEL 6B
#
# Model 6B: MM Shear + Growth-Coupled EDR
#   Agg' = mu * (alpha - (eps/K_eps) * Agg)
#
#   Free  (8): K_Agg, K_d, K_G, K_oG, Y_GLac, alpha, mu_max, K_tau, K_eps
#   Fixed (2): K_Lac = 30.40495     (right-unidentifiable, P50 fit)
#              Y_GX  = 3.6201e-06   (left-unidentifiable,  P50 fit)
#
# Starting values taken from lowest-BIC calibration run.
# Numerical settings: PresetTimeCallback, abstol/reltol=1e-8,
#                     AutoFiniteDiff, LN_BOBYQA, :prev
# =============================================================================

println("="^60)
println("JULIA THREADING VERIFICATION")
println("="^60)
println("Julia threads available: $(Threads.nthreads())")
if Threads.nthreads() < 16
    println("WARNING: Only $(Threads.nthreads()) threads available!")
else
    println("SUCCESS: $(Threads.nthreads()) threads ready!")
end
thread_test = zeros(Int, Threads.nthreads())
Threads.@threads for i in 1:Threads.nthreads()
    thread_test[Threads.threadid()] += 1
end
println("Active threads: $(count(x -> x > 0, thread_test))/$(Threads.nthreads())")
println("="^60)

# =============================================================================
# PACKAGES
# =============================================================================

using DrWatson
@quickactivate "HPSCProliferation"

using CSV, DataFrames, DifferentialEquations, DiffEqCallbacks, ProfileLikelihood
using SciMLBase: ReturnCode
using Optimization, NLopt, Statistics, Optim, LineSearches
using OptimizationOptimJL, OptimizationNLopt, CairoMakie, LaTeXStrings, Printf, Dates
using LinearAlgebra
using XLSX

println("All packages loaded successfully!")

# =============================================================================
# FIXED PARAMETERS — taken from lowest-BIC calibration MLE
# =============================================================================

const FIXED_K_Lac = 30.40495      # P50 size-dependent fit (right-unidentifiable)
const FIXED_Y_GX  = 3.620070e-06  # P50 size-dependent fit (left-unidentifiable)
const FIXED_K_Agg = 6.174043e7    # calibrated (BIC-selected) value; flat profile above this, fixed

println("\nFixed parameters (not profiled):")
println(@sprintf("  K_Lac = %.6g   (right-unidentifiable)", FIXED_K_Lac))
println(@sprintf("  Y_GX  = %.6g   (left-unidentifiable)",  FIXED_Y_GX))

# =============================================================================
# EXPERIMENTAL CONDITIONS
# =============================================================================

const NUM_CONDITIONS = 6

const INITIAL_CONDITIONS_RAW = [
    [20000.0,  14.0, 0.0, 15.0],   # Condition 1: 2e4/40RPM/100mL
    [80000.0,  14.0, 0.0, 15.0],   # Condition 2: 8e4/40RPM/100mL
    [150000.0, 14.0, 0.0, 15.0],   # Condition 3: 15e4/40RPM/100mL
    [80000.0,  15.0, 0.0, 15.0],   # Condition 4: 8e4/30RPM/500mL
    [80000.0,  14.0, 0.0, 15.0],   # Condition 5: 8e4/30RPM/100mL
    [80000.0,  14.0, 0.0, 15.0]    # Condition 6: 8e4/60RPM/100mL
]

const TAU_VALUES = [
    0.003676928,   # Condition 1
    0.003676928,   # Condition 2
    0.003676928,   # Condition 3
    0.002713068,   # Condition 4
    0.002835895,   # Condition 5
    0.005445635    # Condition 6
]

const EPSILON_VALUES = [
    0.000335119,   # Condition 1
    0.000335119,   # Condition 2
    0.000335119,   # Condition 3
    0.000299712,   # Condition 4
    0.000173578,   # Condition 5
    0.000929992    # Condition 6
]

# =============================================================================
# SIZE-DEPENDENT LAGRANGIAN INPUTS (P50)
# tau and eps are interpolated from the Lagrangian size-series at the current
# aggregate size, replacing the fixed per-condition TAU_VALUES/EPSILON_VALUES.
# (TAU_VALUES/EPSILON_VALUES above are retained but no longer used.)
# =============================================================================
const LAG_FILE  = datadir("cfd_raw", "Particles_REBUILT_ALLNONAN_fromTXT.xlsx")
const STAT      = "Median"            # P50. Use "75th Percentile" for P75.
const SIZES     = Float64[200, 400, 600, 800, 1000]
const COND_VRPM = [(100,40),(100,40),(100,40),(500,30),(100,30),(100,60)]

function load_lagrangian_series(stat::String)
    lag = DataFrame(XLSX.readtable(LAG_FILE, "Particles_REBUILT_ALLNONAN"))
    shear = Vector{Vector{Float64}}(undef, NUM_CONDITIONS)
    edr   = Vector{Vector{Float64}}(undef, NUM_CONDITIONS)
    for (c, (V, R)) in enumerate(COND_VRPM)
        sh = Float64[]; ed = Float64[]
        for s in SIZES
            row = lag[(lag.Volume_mL .== V) .& (lag.RPM .== R) .& (lag.Particle_Size_um .== s), :]
            push!(sh, row[1, "Shear $stat"]); push!(ed, row[1, "EDR $stat"])
        end
        shear[c] = sh; edr[c] = ed
    end
    return shear, edr
end

const _LAG        = load_lagrangian_series(STAT)
const SHEAR_SERIES = _LAG[1]
const EDR_SERIES   = _LAG[2]


# Linear interpolation. Below the lowest CFD size (200um) we linearly
# extrapolate using the first segment's slope rather than clamping, since the
# experimental aggregates span ~42-356um and clamping below 200um removed all
# size-dependence in that range. Floored at a small positive value to keep
# tau/eps physical. Above the largest CFD size we still clamp to the endpoint.
@inline function clamp_interp(x, xs::Vector{Float64}, ys::Vector{Float64})
    if x <= xs[1]
        slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
        return max(ys[1] + slope * (x - xs[1]), 1e-8)
    end
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] + t * (ys[i+1] - ys[i])
end

const MAX_DAYS = [6, 6, 6, 5, 5, 6]

const REFRESHMENT_TIMES = [
    [4.02] .* 24,
    [2.02, 3.02, 4.02, 5.02] .* 24,
    [2.02, 3.02, 4.02, 5.02] .* 24,
    [2.02, 3.02, 4.02] .* 24,
    [2.02, 3.02, 4.02] .* 24,
    [2.02, 3.02, 4.02, 5.02] .* 24
]

const MEDIA_CHANGE_SCHEDULE = [
    [0.5],
    [0.5, 0.8, 0.8, 0.8],
    [0.5, 0.8, 0.8, 0.8],
    [0.5, 0.8, 0.8],
    [0.5, 0.8, 0.8],
    [0.5, 0.8, 0.8, 0.8]
]

# =============================================================================
# MEDIA CHANGE CALLBACK
# =============================================================================

function create_media_change_callback(condition_idx::Int)
    refreshment_times = REFRESHMENT_TIMES[condition_idx]
    change_schedule   = MEDIA_CHANGE_SCHEDULE[condition_idx]
    initial_glucose   = INITIAL_CONDITIONS_RAW[condition_idx][2]
    initial_lactate   = INITIAL_CONDITIONS_RAW[condition_idx][3]

    function affect_func!(integrator)
        for (i, time) in enumerate(refreshment_times)
            if abs(integrator.t - time) < 1e-6
                change_percent  = change_schedule[i]
                integrator.u[2] = change_percent * initial_glucose +
                                   (1 - change_percent) * integrator.u[2]
                integrator.u[3] = change_percent * initial_lactate +
                                   (1 - change_percent) * integrator.u[3]
                break
            end
        end
    end

    return PresetTimeCallback(refreshment_times, affect_func!)
end

# =============================================================================
# MODEL DEFINITION
# =============================================================================

# Model 6B — MM Shear + Growth-Coupled EDR  (9 free parameters)
# theta = [K_d, K_G, K_oG, Y_GLac, alpha, mu_max, K_tau, K_eps]  (K_Agg fixed)
# full_params = [theta..., tau_val, eps_val]
# K_Lac and Y_GX read from module-level constants FIXED_K_Lac / FIXED_Y_GX
# Agg' = mu * (alpha - (eps/K_eps) * Agg)
function make_model6B(sh::Vector{Float64}, ed::Vector{Float64})
    return function (du, u, p, t)
        X, G, Lac, Agg = u
        K_d, K_G, K_oG, Y_GLac, alpha, mu_max, K_tau, K_eps = p

        K_Agg = FIXED_K_Agg
        K_Lac = FIXED_K_Lac
        Y_GX  = FIXED_Y_GX

        tau = clamp_interp(Agg, SIZES, sh)   # size-dependent shear (live)
        eps = clamp_interp(Agg, SIZES, ed)   # size-dependent EDR (live)

        kappa = K_d * exp(-t / 12.0)
        shear = K_tau / (tau + K_tau)
        mu    = mu_max * (G/(G+K_G)) * (K_Lac/(Lac+K_Lac)) *
                (K_Agg/(Agg^3+K_Agg)) * shear
        lam   = K_oG * mu * G
        sigma = Y_GX * mu + Y_GLac * lam

        du[1] = (mu - kappa) * X
        du[2] = -sigma * X
        du[3] = lam * X
        du[4] = mu * (alpha - (eps / K_eps) * Agg)
    end
end

# =============================================================================
# DATA LOADING
# =============================================================================

function get_condition_data(df::DataFrame, condition_set::Int)
    condition_data = filter(row -> row.ConditionSet == condition_set, df)
    grouped        = groupby(condition_data, Symbol("t (days)"))

    times        = Float64[]
    cells_mean   = Float64[]; cells_std   = Float64[]
    glucose_mean = Float64[]; glucose_std = Float64[]
    lactate_mean = Float64[]; lactate_std = Float64[]
    agg_mean     = Float64[]; agg_std     = Float64[]

    for group in grouped
        push!(times, first(group[!, Symbol("t (days)")]))

        cells = [Float64(x) for x in group.ExperimentalCells_per_mL if !ismissing(x)]
        push!(cells_mean, isempty(cells) ? NaN : mean(cells))
        push!(cells_std,  length(cells) > 1 ? std(cells) :
                          (isempty(cells) ? NaN : mean(cells) * 0.1))

        glucose = [Float64(x) for x in group.ExperimentalGlucose_mM if !ismissing(x)]
        push!(glucose_mean, isempty(glucose) ? NaN : mean(glucose))
        push!(glucose_std,  length(glucose) > 1 ? std(glucose) : 0.5)

        lactate = [Float64(x) for x in group.ExperimentalLactate_mM if !ismissing(x)]
        push!(lactate_mean, isempty(lactate) ? NaN : mean(lactate))
        push!(lactate_std,  length(lactate) > 1 ? std(lactate) : 0.5)

        agg = [Float64(x) for x in group.AggregateSize_um if !ismissing(x)]
        push!(agg_mean, isempty(agg) ? NaN : mean(agg))
        push!(agg_std,  length(agg) > 1 ? std(agg) :
                        (isempty(agg) ? NaN : mean(agg) * 0.1))
    end

    sort_idx = sortperm(times)
    return (
        times     = times[sort_idx],
        cells     = (mean = cells_mean[sort_idx],   std = cells_std[sort_idx]),
        glucose   = (mean = glucose_mean[sort_idx], std = glucose_std[sort_idx]),
        lactate   = (mean = lactate_mean[sort_idx], std = lactate_std[sort_idx]),
        aggregate = (mean = agg_mean[sort_idx],     std = agg_std[sort_idx])
    )
end

function cache_experimental_data(df::DataFrame)
    return Dict{Int,Any}(i => get_condition_data(df, i) for i in 1:NUM_CONDITIONS)
end

const EXP_DATA_CACHE = Ref{Dict{Int,Any}}(Dict{Int,Any}())

# =============================================================================
# SIMULATION
# =============================================================================

function simulate_condition(params::AbstractVector, condition_idx::Int,
                             save_times::Vector{Float64})
    T = eltype(params)

    u0          = T.(INITIAL_CONDITIONS_RAW[condition_idx])
    rhs         = make_model6B(SHEAR_SERIES[condition_idx], EDR_SERIES[condition_idx])
    full_params = collect(params)                 # theta (9 free params); tau,eps now size-dependent
    tspan       = (0.0, Float64(MAX_DAYS[condition_idx]) * 24.0)

    prob             = ODEProblem(rhs, u0, tspan, full_params)
    callback         = create_media_change_callback(condition_idx)
    valid_save_times = filter(t -> 0 <= t <= tspan[2], save_times)

    try
        sol = solve(prob, Vern7(),
                    callback = callback,
                    tstops   = REFRESHMENT_TIMES[condition_idx],
                    abstol   = 1e-8,
                    reltol   = 1e-8,
                    saveat   = valid_save_times,
                    maxiters = 10_000)

        if sol.retcode == :Success || sol.retcode == SciMLBase.ReturnCode.Success
            return sol
        else
            return nothing
        end
    catch
        return nothing
    end
end

# =============================================================================
# LOG-LIKELIHOOD
# =============================================================================

function loglik(theta::AbstractVector, data)
    T  = eltype(theta)
    ll = zero(T)

    for cond_idx in 1:NUM_CONDITIONS
        exp_data   = EXP_DATA_CACHE[][cond_idx]
        save_times = exp_data.times .* 24.0
        sol        = simulate_condition(theta, cond_idx, save_times)

        if sol === nothing
            return T(-1e10)
        end

        for (i, t_day) in enumerate(exp_data.times)
            t_idx = argmin(abs.(sol.t .- t_day * 24.0))

            if !isnan(exp_data.cells.mean[i]) && exp_data.cells.std[i] > 0
                sigma    = max(exp_data.cells.std[i], exp_data.cells.mean[i] * 0.05)
                residual = (sol[1, t_idx] - exp_data.cells.mean[i]) / sigma
                ll      -= T(0.5) * residual^2 + T(log(sigma))
            end

            if !isnan(exp_data.glucose.mean[i]) && exp_data.glucose.std[i] > 0
                sigma    = max(exp_data.glucose.std[i], 0.5)
                residual = (sol[2, t_idx] - exp_data.glucose.mean[i]) / sigma
                ll      -= T(0.5) * residual^2 + T(log(sigma))
            end

            if !isnan(exp_data.lactate.mean[i]) && exp_data.lactate.std[i] > 0
                sigma    = max(exp_data.lactate.std[i], 0.5)
                residual = (sol[3, t_idx] - exp_data.lactate.mean[i]) / sigma
                ll      -= T(0.5) * residual^2 + T(log(sigma))
            end

            if !isnan(exp_data.aggregate.mean[i]) && exp_data.aggregate.std[i] > 0
                sigma    = max(exp_data.aggregate.std[i],
                               exp_data.aggregate.mean[i] * 0.05)
                residual = (sol[4, t_idx] - exp_data.aggregate.mean[i]) / sigma
                ll      -= T(0.5) * residual^2 + T(log(sigma))
            end
        end
    end

    return ll
end

# =============================================================================
# PARAMETER DEFINITIONS — Model 6B (9 free parameters)
# Starting values from lowest-BIC calibration run.
# K_Lac and Y_GX removed (fixed above).
# Bounds are wide and physically motivated — NOT set to match prior CI estimates.
# =============================================================================

const START = [
    0.136049,       # K_d
    16.081846,      # K_G
    3.060877e-06,   # K_oG      (K_Lac omitted — fixed)
    0.417865,       # Y_GLac    (Y_GX omitted  — fixed)
    62.303804,      # alpha
    0.127683,       # mu_max
    0.00934151,     # K_tau
    0.0111309       # K_eps
]


const PARAM_NAMES = [
    "K_d", "K_G", "K_oG", "Y_GLac",
    "alpha", "mu_max", "K_tau", "K_eps"
]

const PARAM_LABELS = [
    L"K_d", L"K_G", L"K_{oG}", L"Y_{GLac}",
    L"\alpha",  L"\mu_{\max}", L"K_{\tau}", L"K_{\varepsilon}"
]

# Wide bounds — give the inner optimizer room to move nuisance parameters
# freely at every profile point. Never set these to match CI estimates.
const LB = [
    0.001,   # K_d
    0.1,     # K_G
    1e-8,    # K_oG
    0.01,    # Y_GLac
    1.0,     # alpha
    0.001,   # mu_max
    1e-5,    # K_tau
    1e-5     # K_eps
]

const UB = [
    1.0,     # K_d
    40.0,    # K_G
    1e-4,    # K_oG
    2.0,     # Y_GLac
    500.0,   # alpha
    2.0,     # mu_max
    0.1,     # K_tau
    0.05     # K_eps
]

const N_PARAMS = length(START)

# =============================================================================
# EXECUTION
# =============================================================================

println("\n" * "="^80)
println("PROFILE LIKELIHOOD — MODEL 6B")
println("  MM Shear + Growth-Coupled EDR  (8 free, 3 fixed)")
println("  Fixed: K_Lac = $(FIXED_K_Lac), Y_GX = $(FIXED_Y_GX)")
println("="^80)

data_file = datadir("exp_raw", "Combined100mLand500mLDataAllRPM.csv")
if !isfile(data_file)
    error("Data file not found: $data_file\n   Current dir: $(pwd())")
end
println("Found: $(abspath(data_file))")

isdir("plots") || mkdir("plots")
current_date  = Dates.format(now(), "U_d_yyyy")
output_folder = joinpath("plots", "Model6B_ProfileLikelihood_$(current_date)")
isdir(output_folder) || mkdir(output_folder)
println("Output folder: $output_folder")

df = CSV.read(data_file, DataFrame)
println("\nLoaded: $(nrow(df)) rows")
println("Caching experimental data...")
EXP_DATA_CACHE[] = cache_experimental_data(df)
for i in 1:NUM_CONDITIONS
    println("  Condition $i: $(length(EXP_DATA_CACHE[][i].times)) timepoints")
end

# ---- Warm-up ----
println("\nWarm-up evaluation...")
_ = loglik(START, nothing)
times_eval = [@elapsed loglik(START, nothing) for _ in 1:3]
println("Avg likelihood eval: $(round(mean(times_eval)*1000, digits=1)) ms")

# ---- Parameter table ----
println("\nStarting parameters and bounds:")
println("="^90)
println(@sprintf("%-12s %18s %18s %18s", "Parameter", "Lower", "Start", "Upper"))
println("="^90)
for i in 1:N_PARAMS
    println(@sprintf("%-12s %18.4g %18.4g %18.4g",
                     PARAM_NAMES[i], LB[i], START[i], UB[i]))
end
println("="^90)

# ---- MLE ----
println("\nRunning MLE...")

likeliProb = LikelihoodProblem(
    loglik, START;
    data        = nothing,
    syms        = PARAM_NAMES,
    f_kwargs    = (adtype = Optimization.AutoFiniteDiff(),),
    prob_kwargs = (lb = LB, ub = UB),
)

mle_start = time()

# Robust MLE. The previous Fminbox(BFGS)+AutoFiniteDiff stalled because K_Agg
# enters as K_Agg/(Agg^3+K_Agg); near START its finite-difference gradient is
# tiny, so BFGS never stepped and left K_Agg pinned at START (lopsided profile).
# LN_BOBYQA is derivative-free (the same algorithm used for the profiles) and
# does not depend on that gradient. We run it from START plus several
# perturbations and keep the best, so the optimiser reliably reaches the joint
# optimum and the profiles come out properly shaped.
function run_mle(p0)
    prob = LikelihoodProblem(
        loglik, p0;
        data        = nothing,
        syms        = PARAM_NAMES,
        f_kwargs    = (adtype = Optimization.AutoFiniteDiff(),),
        prob_kwargs = (lb = LB, ub = UB),
    )
    sol = mle(prob, NLopt.LN_BOBYQA();
              maxiters = 20_000, xtol_rel = 1e-10, ftol_rel = 1e-12)
    return sol, loglik(sol.mle, nothing)
end

best_sol, best_ll = run_mle(START)
for f in (0.5, 1.5, 2.0)                # mild multistart (K_Agg now fixed)
    p = clamp.(collect(START) .* f, LB, UB)
    sol_i, ll_i = run_mle(p)
    if ll_i > best_ll
        global best_ll, best_sol = ll_i, sol_i
    end
end

likeliSol = best_sol
mle_time  = time() - mle_start

mle_params = likeliSol.mle
ll_start   = loglik(START,      nothing)
ll_mle     = loglik(mle_params, nothing)

println("MLE completed in $(round(mle_time, digits=2)) seconds!")
println("\nLog-likelihood:")
println("  Starting: $(round(ll_start, digits=4))")
println("  MLE:      $(round(ll_mle,   digits=4))")
println("  Delta:    $(round(ll_mle - ll_start, digits=4))")

println("\nMLE Results:")
println("="^80)
println(@sprintf("%-12s %15s %15s %15s", "Parameter", "Start", "MLE", "% Diff"))
println("="^80)
for i in 1:N_PARAMS
    diff_pct = (mle_params[i] - START[i]) / abs(START[i]) * 100
    println(@sprintf("%-12s %15.5g %15.5g %+14.1f%%",
                     PARAM_NAMES[i], START[i], mle_params[i], diff_pct))
end
println("="^80)

# ---- Boundary check ----
println("\nBoundary check:")
println("="^80)
boundary_params = String[]
for i in 1:N_PARAMS
    dist_lb     = (mle_params[i] - LB[i]) / (UB[i] - LB[i])
    dist_ub     = (UB[i] - mle_params[i]) / (UB[i] - LB[i])
    at_boundary = dist_lb < 0.01 || dist_ub < 0.01
    at_boundary && push!(boundary_params, PARAM_NAMES[i])
    status = at_boundary ? "BOUNDARY" : "ok"
    println(@sprintf("  %-12s: LB dist=%.4f, UB dist=%.4f  %s",
                     PARAM_NAMES[i], dist_lb, dist_ub, status))
end
println("="^80)
if isempty(boundary_params)
    println("  No boundary solutions!")
else
    println("  Boundary parameters: $(join(boundary_params, ", "))")
end

# Rebuild the likelihood problem anchored at the MLE so the profile ranges and
# starting point are consistent with the optimum found above.
likeliProb = LikelihoodProblem(
    loglik, collect(mle_params);
    data        = nothing,
    syms        = PARAM_NAMES,
    f_kwargs    = (adtype = Optimization.AutoFiniteDiff(),),
    prob_kwargs = (lb = LB, ub = UB),
)

# ---- Profile likelihood ----
println("\nComputing profile likelihood ($N_PARAMS parameters, 100 points each)...")

resolution   = 100
resolutions  = fill(resolution, N_PARAMS)
param_ranges = construct_profile_ranges(likeliSol, LB, UB, resolutions)

# Focus the profile range for parameters whose likelihood is flat far from the
# optimum (K_Agg and K_G saturate because Agg^3 stays ~1e7 while K_Agg ~1e8, and
# glucose stays near saturation for K_G). The 95% CI still closes inside these
# windows; restricting the plotted/profiled range to the informative region
# keeps the parabola from being crushed against a flat tail. Ranges are built
# as two half-ranges (MLE -> lower, MLE -> upper) as ProfileLikelihood expects.
function focused_range(mle_val, lo, hi, n)
    half = n ÷ 2
    left  = LinRange(mle_val, lo, half)
    right = LinRange(mle_val, hi, n - half)
    return (left, right)
end

let
    iG = findfirst(==("K_G"), PARAM_NAMES)
    if iG !== nothing
        param_ranges[iG] = focused_range(mle_params[iG], 2.0, 30.0, resolution)
    end
end

profile_start = time()
likeliProf    = ProfileLikelihood.profile(likeliProb, likeliSol;
    param_ranges                 = param_ranges,
    parallel                     = true,
    conf_level                   = 0.95,
    alg                          = NLopt.LN_BOBYQA(),
    next_initial_estimate_method = :prev,
    progress                     = true
)
profile_time = time() - profile_start
println("Profile completed in $(round(profile_time/60, digits=1)) minutes!")

# ---- CI table ----
println("\n95% Confidence Intervals:")
println("="^105)
println(@sprintf("%-12s %15s %15s %15s %15s  %s",
                 "Parameter", "Lower CI", "MLE", "Upper CI", "Width %", "Status"))
println("="^105)

well_identified   = String[]
poorly_identified = String[]

for i in 1:N_PARAMS
    ci        = likeliProf.confidence_intervals[i]
    width_pct = (ci.upper - ci.lower) / abs(mle_params[i]) * 100
    at_lb     = abs(ci.lower - LB[i]) / (abs(LB[i]) + 1e-30) < 0.01
    at_ub     = abs(ci.upper - UB[i]) / (abs(UB[i]) + 1e-30) < 0.01

    status = if at_lb && at_ub;    "Both bounds"
           elseif at_lb;           "Left-unidentifiable"
           elseif at_ub;           "Right-unidentifiable"
           elseif width_pct > 200; "Wide CI"
           else;                   "Identified"
           end

    if width_pct < 200 && !at_lb && !at_ub
        push!(well_identified, PARAM_NAMES[i])
    else
        push!(poorly_identified, PARAM_NAMES[i])
    end

    println(@sprintf("%-12s %15.5g %15.5g %15.5g %14.1f%%  %s",
                     PARAM_NAMES[i], ci.lower, mle_params[i],
                     ci.upper, width_pct, status))
end
println("="^105)

println("\nIDENTIFIABILITY SUMMARY — Model 6B")
println("  Identified     : $(isempty(well_identified)   ? "none" : join(well_identified,   ", "))")
println("  Non-identified : $(isempty(poorly_identified) ? "none" : join(poorly_identified, ", "))")

println("\nMLE parameters with 95% CI:")
println("  " * "-"^70)
for i in 1:N_PARAMS
    ci = likeliProf.confidence_intervals[i]
    @printf("  %-12s = %12.5g   95%% CI: [%10.5g, %10.5g]\n",
            PARAM_NAMES[i], mle_params[i], ci.lower, ci.upper)
end
println("  Fixed (not profiled):")
@printf("  %-12s = %12.5g\n", "K_Lac", FIXED_K_Lac)
@printf("  %-12s = %12.5g\n", "Y_GX",  FIXED_Y_GX)

# ---- Plot — 8 params → 2×4 grid ----
figProfile = plot_profiles(likeliProf;
    latex_names = PARAM_LABELS,
    show_mles   = true,
    shade_ci    = true,
    nrow        = 2,
    ncol        = 4,
    axis_kwargs = (width=200, height=140))

CairoMakie.resize_to_layout!(figProfile)

timestamp        = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
profile_filename = joinpath(output_folder, "ProfileLikelihood_Model6B_$(timestamp).png")
save(profile_filename, figProfile, px_per_unit=3)
println("\nProfile plot saved: $profile_filename")

println("\n" * "="^80)
println("COMPLETE")
println("Output: $output_folder")
println("="^80)




