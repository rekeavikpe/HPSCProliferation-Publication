module AggregateDayOne

include("MixedExposure.jl")
using .MixedExposure
using .MixedExposure.RevisionCommon
using .RevisionCommon.Workflow, .Workflow.Calibration, .Calibration.Proliferation
using CSV, Dates, Printf, TOML, Statistics
using .Proliferation: ODEProblem, solve, successful_retcode, Tsit5

export AggregateDayOneContext, aggregate_dayone, comparison_context,
       simulate_aggregate_dayone, run_aggregate_dayone

const MODEL = 6
const REFERENCE = "Original D0 initialization"
const DAYONE = "Aggregate initialized at D1"
const RESET_H = 24.0
const FIRST_DAY_SIZE_UM = 15.0

"The aggregate equation starts at day 1. Other states start at day 0 with a fixed size input."
struct AggregateDayOneContext
    config
    observations
    conditions
    y
    sigma
end

eligible(r) = r.state != "aggregate" || r.day > 1

"Apply the same scoring mask to the original model without changing its trajectory."
function comparison_context(base)
    ix = findall(eligible, base.observations)
    merge(base, (observations=base.observations[ix], y=base.y[ix], sigma=base.sigma[ix]))
end

function aggregate_dayone(base)
    qs = Dict{Int,Any}()
    for (c,q) in base.conditions
        rows = filter(r -> r.condition == c && r.state == "aggregate" && r.day == 1,
                      base.observations)
        length(rows) == 1 || error("Expected one day 1 aggregate mean for condition $c")
        r = only(rows)
        isfinite(r.mean) && r.mean > 0 || error("Invalid day 1 aggregate size for condition $c")
        q.u0[4] == FIRST_DAY_SIZE_UM || error("Expected the approved 15 µm first-day size input")
        last(q.times) > RESET_H || error("No post-day-1 observations for condition $c")
        qs[c] = merge(q, (aggregate_day1=Float64(r.mean), aggregate_day1_n=r.n_vessels))
    end
    scored = comparison_context(base)
    AggregateDayOneContext(base.config, scored.observations, qs, scored.y, scored.sigma)
end

"Original Model 6 cell and metabolite equations with a fixed size input, without an aggregate ODE."
function first_day_rhs!(du,u,p,t,tau; lactate_inhibition=true)
    X,G,Lac = max(u[1],1.0),max(u[2],.001),max(u[3],0.0)
    inhibition = lactate_inhibition ? p.K_Lac/(Lac+p.K_Lac) : 1.0
    growth = p.mu_max*G/(G+p.K_G)*inhibition*p.K_Agg/(FIRST_DAY_SIZE_UM^3+p.K_Agg)*p.K_tau/(tau+p.K_tau)
    lactate_rate = p.K_oG*growth*G
    glucose_rate = p.Y_GX*growth+p.Y_GLac*lactate_rate
    du[1] = (growth-p.K_d*exp(-t/12))*X
    du[2] = -glucose_rate*X
    du[3] = lactate_rate*X
    nothing
end

"Integrate only X, G and Lac through D1. The fourth returned row records the fixed input, not a size prediction."
function simulate_first_day(p,q,times,feeds;reltol,abstol,lactate_inhibition)
    tau,edr = q.lookup(FIRST_DAY_SIZE_UM)
    all(x -> isfinite(x) && x >= 0,(tau,edr)) || error("Invalid first-day exposure")
    feed_times = [f.time_h for f in feeds]
    issorted(feed_times) && allunique(feed_times) && all(t->0<t<=RESET_H,feed_times) || error("Invalid first-day feed times")
    for f in feeds
        exchange(Float64.(q.u0),f.fraction,q.glucose)
    end
    rhs!(du,u,unused,t) = first_day_rhs!(du,u,p,t,tau;lactate_inhibition)
    predicted = fill(FIRST_DAY_SIZE_UM,4,length(times))
    u,start = Float64.(q.u0[1:3]),0.0
    events = NamedTuple[]
    for stop in sort(unique(vcat(Float64.(feed_times),RESET_H)))
        solution = solve(ODEProblem(rhs!,u,(start,stop)),Tsit5();reltol,abstol,dense=true,maxiters=1_000_000)
        successful_retcode(solution) || error("First-day integration failed: $(solution.retcode)")
        all(s->all(isfinite,s),solution.u) || error("Nonfinite first-day solution")
        for (j,t) in enumerate(times)
            if (start==0 && t==0) || start<t<=stop
                predicted[1:3,j] = solution(t)
            end
        end
        u = copy(solution.u[end])
        index = findfirst(==(stop),feed_times)
        if index !== nothing
            before = vcat(u,FIRST_DAY_SIZE_UM)
            after = exchange(before,feeds[index].fraction,q.glucose)
            push!(events,(time_h=stop,fraction=feeds[index].fraction,before,after))
            u = after[1:3]
        end
        start = stop
    end
    (;predictions=predicted,events)
end

"Hold the size input at 15 µm before D1, then start the aggregate ODE at its measured D1 size."
function simulate_aggregate_dayone(m,p,q,times=q.times; reltol=1e-12,abstol=1e-13,
                                  lactate_inhibition=true)
    m == MODEL || error("The aggregate-only day 1 analysis is defined for Model 6")
    keys(p) == parameter_names(m) || error("Parameter names/order mismatch")
    parameters(m,collect(values(p)))
    length(q.u0)==4 && all(x->isfinite(x) && x>=0,q.u0) && q.u0[4]==FIRST_DAY_SIZE_UM ||
        error("Invalid initial states or first-day size input")
    !isempty(times) && issorted(times) && allunique(times) &&
        all(t -> isfinite(t) && t >= 0, times) && last(times) >= RESET_H ||
        error("Invalid sample times for aggregate day 1 initialization")
    isfinite(q.aggregate_day1) && q.aggregate_day1 > 0 || error("Invalid day 1 size")
    early_times = sort(unique(vcat(filter(<=(RESET_H), Float64.(times)), RESET_H)))
    early_feeds = filter(f -> f.time_h <= RESET_H, q.feeds)
    early = simulate_first_day(p,q,early_times,early_feeds;reltol,abstol,lactate_inhibition)
    # At a feeding endpoint observations are pre-exchange, but continuation is post-exchange.
    before = copy(early.predictions[:,end])
    carry = copy(before)
    at_reset = filter(e -> e.time_h == RESET_H, early.events)
    !isempty(at_reset) && (carry = copy(only(at_reset).after))
    carry[4] = q.aggregate_day1
    predicted = Matrix{Float64}(undef,4,length(times))
    for (j,t) in enumerate(times)
        t > RESET_H && continue
        predicted[:,j] = early.predictions[:,findfirst(==(t),early_times)]
        t == RESET_H && (predicted[4,j] = q.aggregate_day1)
    end
    events = copy(early.events)
    late_ix = findall(>(RESET_H),times)
    if !isempty(late_ix)
        # Kd*exp(-(s+24)/12) preserves the death clock when the solver's time is shifted.
        shifted = merge(p,(K_d=p.K_d*exp(-RESET_H/12),))
        late_feeds = [(time_h=f.time_h-RESET_H,fraction=f.fraction)
                      for f in q.feeds if RESET_H < f.time_h <= last(times)]
        late = simulate(m,shifted,carry,Float64.(times[late_ix]).-RESET_H,late_feeds;
                        fresh_glucose=q.glucose,exposure=q.lookup,reltol,abstol,lactate_inhibition)
        predicted[:,late_ix] = late.predictions
        append!(events,[merge(e,(time_h=e.time_h+RESET_H,)) for e in late.events])
    end
    all(x -> isfinite(x) && x >= 0, predicted) || error("Invalid aggregate day 1 prediction")
    (;predictions=predicted,events,initialization=(time_h=RESET_H,before,after=copy(carry)))
end

function Workflow.prediction(ctx::AggregateDayOneContext,m,p;reltol=1e-12,abstol=1e-13)
    yp = zeros(length(ctx.y))
    for (c,q) in ctx.conditions
        sim = simulate_aggregate_dayone(m,p,q;reltol,abstol,
                                       lactate_inhibition=get(ctx.config,"lactate_inhibition",true))
        for (i,r) in enumerate(ctx.observations)
            r.condition == c || continue
            yp[i] = sim.predictions[findfirst(==(r.state),STATES),findfirst(==(24r.day),q.times)]
        end
    end
    yp
end

function score_row(ctx,p,formulation,cohort)
    nll = gaussian_nll(ctx.y,prediction(ctx,MODEL,p),ctx.sigma)
    k = length(free_names(ctx,MODEL))
    (;formulation,cohort,n=length(ctx.y),k,nll,bic=bic(nll,k,length(ctx.y)))
end

function metrics(ctx,yp,formulation,cohort)
    rows = NamedTuple[]
    for c in (0,sort(collect(keys(ctx.conditions)))...),s in (STATES...,"all")
        ix = findall(r -> (c==0 || r.condition==c) && (s=="all" || r.state==s) && r.mean!=0,
                     ctx.observations)
        isempty(ix) && continue
        push!(rows,(formulation,cohort,condition=c,state=s,n=length(ix),
                    mrpe_percent=100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix]))))
    end
    rows
end

function export_cohort!(points,curves,initials,base,ctx,p,formulation,cohort)
    for (c,q) in ctx.conditions
        times = sort(unique(vcat(collect(0.:.5:maximum(q.times)),q.times,RESET_H)))
        reset = ctx isa AggregateDayOneContext
        sim = reset ? simulate_aggregate_dayone(MODEL,p,q,times;
                            lactate_inhibition=get(ctx.config,"lactate_inhibition",true)) :
                      simulate(MODEL,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,
                            exposure=q.lookup,reltol=1e-12,abstol=1e-13,
                            lactate_inhibition=get(ctx.config,"lactate_inhibition",true))
        for (j,t) in enumerate(times),(k,s) in enumerate(STATES)
            # The proposed aggregate trajectory begins at D1. No pre-D1 size curve is exported.
            reset && s=="aggregate" && t<RESET_H && continue
            side = reset && s=="aggregate" && t==RESET_H ? "initial" : "pre"
            push!(curves,(formulation,cohort,condition=c,day=t/24,state=s,side,value=sim.predictions[k,j]))
            for e in sim.events
                e.time_h == t || continue
                value = reset && s=="aggregate" && t==RESET_H ? q.aggregate_day1 : e.after[k]
                push!(curves,(formulation,cohort,condition=c,day=t/24,state=s,side="post",value))
            end
        end
        for r in base.observations
            r.condition == c || continue
            j = findfirst(==(24r.day),times)
            predicted = sim.predictions[findfirst(==(r.state),STATES),j]
            role = eligible(r) ? "scored" : (reset ? "initial_condition" : "excluded_for_matched_comparison")
            push!(points,(formulation,cohort,condition=c,day=r.day,state=r.state,
                          observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,predicted,role))
        end
        if reset
            before,after = sim.initialization.before,sim.initialization.after
            for (i,s) in enumerate(STATES)
                push!(initials,(cohort,condition=c,state=s,day=1,before=before[i],after=after[i],
                                source=i==4 ? "measured day 1 mean from this cohort" : "model trajectory from day 0",
                                n_day1_aggregate=q.aggregate_day1_n))
            end
        end
    end
end

function write_report(out,scores,comparison)
    open(joinpath(out,"index.html"),"w") do io
        print(io,"<!doctype html><meta charset='utf-8'><title>Aggregate day 1 initialization</title><style>body{font:17px Arial;max-width:1150px;margin:35px auto}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{padding:9px;border-bottom:1px solid #ccc;text-align:right}td:first-child,th:first-child{text-align:left}</style><h1>Model 6 with aggregate size initialized at day 1</h1><p>Both inputs are Lagrangian P50. The new fit uses the original 19 training cultures and all 11 Model 6 parameters. Cells, glucose and lactate retain their day 0 initial conditions. Before day 1, only the cell, glucose and lactate equations are integrated. Their size-dependent terms use a fixed 15 µm input. No aggregate-size equation is integrated in that interval. At day 1, the aggregate equation starts from the measured size mean within that cohort, while the other three predicted states continue without measurement resets. The culture clock and feeding schedule are preserved.</p><p><strong>First-day assumption:</strong> the 15 µm value is a fixed model input, not a measured aggregate diameter or a prediction of aggregation. The existing CFD interpolation uses its 200 µm endpoint for sizes below the available grid. That rule is retained for P50 shear and EDR. Aggregate growth limitation remains in the cell-growth equation with the fixed 15 µm input.</p><p>The day 1 aggregate measurement is an initial condition, not a scored prediction. The original D0 model is re-evaluated on exactly the same remaining observations for comparison. The new selection BIC must be compared with the matched original score here, not with the earlier score that included day 1 aggregate errors. Predictions at 500 mL and for the other experiments use their own measured day 1 size without fitting parameters to those datasets.</p><h2>Selection BIC on matched observations</h2><table><tr><th>Formulation</th><th>Observations</th><th>Parameters</th><th>NLL</th><th>BIC</th></tr>")
        for r in scores
            r.cohort == "selection40" || continue
            @printf(io,"<tr><td>%s</td><td>%d</td><td>%d</td><td>%.5f</td><td>%.5f</td></tr>",r.formulation,r.n,r.k,r.nll,r.bic)
        end
        print(io,"</table><h2>MRPE on matched observations</h2><p>Negative differences favour the D1 aggregate fit. Aggregate errors use days after day 1. The other three states retain all previously scored days. Overall MRPE pools eligible condition-day-state observations. The all-100 mL comparison is descriptive. The iPSC comparisons are direct transfers, not cell-line recalibrations.</p><table><tr><th>Subset</th><th>State</th><th>Original D0 (%)</th><th>Aggregate D1 (%)</th><th>Difference (percentage points)</th></tr>")
        for r in comparison
            r.condition == 0 || continue
            @printf(io,"<tr><td>%s</td><td>%s</td><td>%.3f</td><td>%.3f</td><td>%+.3f</td></tr>",r.cohort,r.state,r.reference_mrpe,r.aggregate_day1_mrpe,r.difference_pp)
        end
        print(io,"</table><p>Day 1 sizes and the retained first-three-state predictions are recorded in initial_states.csv. Unmodified original-model scores are supplied separately in reference_original_scores.csv for traceability.</p>")
    end
end

function run_aggregate_dayone(fit,splitfit,out)
    fit,splitfit,out = abspath.((fit,splitfit,out))
    starts = parse(Int,get(ENV,"HPSC_AGGDAY1_STARTS","100"))
    budget = parse(Int,get(ENV,"HPSC_AGGDAY1_EVALUATIONS","20000"))
    iterations = parse(Int,get(ENV,"HPSC_AGGDAY1_ITERATIONS","3000"))
    starts>0 && budget>=100 && iterations>0 || error("Invalid budgets")
    prep,p0 = MixedExposure.source_contexts(fit,splitfit)
    extra = [@__FILE__,joinpath(@__DIR__,"run_aggregate_day1.jl"),joinpath(@__DIR__,"MixedExposure.jl"),
             joinpath(fit,"calibration.toml"),joinpath(fit,"model5/selected_parameters.csv"),
             joinpath(splitfit,"split.csv"),joinpath(splitfit,"model5/selected_parameters.csv"),
             joinpath(splitfit,"model5/summary.csv")]
    settings = Dict("analysis"=>"model6_p50p50_aggregate_only_day1_fixed15_v2","starts"=>starts,
                    "evaluations"=>budget,"iterations"=>iterations,"julia"=>string(VERSION),
                    "first_day"=>"three-state ODE with fixed 15 um size input; no aggregate ODE before 24 h",
                    "first_day_size_um"=>FIRST_DAY_SIZE_UM,
                    "below_CFD_grid"=>"retain nearest endpoint (200 um) P50 exposure",
                    "day1_change"=>"aggregate diameter only, from cohort-specific measured mean",
                    "scoring"=>"exclude aggregate day<=1 for both comparison arms",
                    "selection_used_for_parameter_fitting"=>false)
    lock_inputs(out,settings,extra)
    progress(s) = open(io -> println(io,Dates.now()," ",s),joinpath(out,"progress.log"),"a")
    CSV.write(joinpath(out,"split.csv"),prep.split.assignments)
    save_parameters(joinpath(out,"reference_parameters.csv"),p0)
    train = aggregate_dayone(prep.train)
    progress("START aggregate_day1_fit")
    p,_ = fit_context(train,MODEL,joinpath(out,"aggregate_day1_fit");starts,budget,iterations)
    progress("FINISHED aggregate_day1_fit")
    progress("START evaluation")
    scores = NamedTuple[]; original_scores = NamedTuple[]; allmetrics = NamedTuple[]
    points = NamedTuple[]; curves = NamedTuple[]; initials = NamedTuple[]
    for (cohort,raw,settings) in MixedExposure.evaluations(prep)
        base = MixedExposure.replace_exposure(raw;edr_percentile=50,settings)
        for (tag,ctx,par) in ((REFERENCE,comparison_context(base),p0),
                              (DAYONE,aggregate_dayone(base),p))
            append!(allmetrics,metrics(ctx,prediction(ctx,MODEL,par),tag,cohort))
            if cohort in ("training60","selection40","validation500")
                push!(scores,score_row(ctx,par,tag,cohort))
            end
            export_cohort!(points,curves,initials,base,ctx,par,tag,cohort)
        end
        cohort in ("training60","selection40","validation500") &&
            push!(original_scores,score_row(base,p0,REFERENCE,cohort))
    end
    key(r) = (r.cohort,r.condition,r.state)
    reference = Dict(key(r)=>r for r in allmetrics if r.formulation==REFERENCE)
    comparison = [(cohort=r.cohort,condition=r.condition,state=r.state,n=r.n,
                   reference_mrpe=reference[key(r)].mrpe_percent,aggregate_day1_mrpe=r.mrpe_percent,
                   difference_pp=r.mrpe_percent-reference[key(r)].mrpe_percent)
                  for r in allmetrics if r.formulation==DAYONE]
    CSV.write(joinpath(out,"scores.csv"),scores)
    CSV.write(joinpath(out,"selection_scores.csv"),filter(r->r.cohort=="selection40",scores))
    CSV.write(joinpath(out,"reference_original_scores.csv"),original_scores)
    CSV.write(joinpath(out,"state_errors.csv"),allmetrics)
    CSV.write(joinpath(out,"comparison.csv"),comparison)
    CSV.write(joinpath(out,"predictions.csv"),points)
    CSV.write(joinpath(out,"curves.csv"),curves)
    CSV.write(joinpath(out,"initial_states.csv"),initials)
    write_report(out,scores,comparison)
    progress("FINISHED evaluation")
    write(joinpath(out,"COMPLETE"),"Fit and matched-observation evaluations completed.\n")
    println("Aggregate day 1 comparison saved: ",joinpath(out,"index.html"))
end

end
