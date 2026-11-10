module MixedExposure
include("run_david_analyses.jl")
using Dates, Printf
export mixed_lookup, replace_exposure, media_context, source_contexts, evaluations, run_mixed

const MODEL = 6 # publication and mechanistic Model 6, stored as model5 in original fits
const BASELINE = "P50 shear / P50 EDR"
const MIXED = "P50 shear / P75 EDR"

function mixed_lookup(volume, rpm; edr_percentile=75)
    rows=collect(CSV.File(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))
    function grid(percentile, speed)
        g=sort(filter(r->r.scenario=="lag_p$percentile" && r.volume_ml==volume && r.rpm==speed,rows);by=r->r.size_um)
        [r.size_um for r in g]==[200,400,600,800,1000] || error("Incomplete exposure grid: $volume mL, $speed RPM, P$percentile")
        g
    end
    speeds=sort(unique(r.rpm for r in rows if r.scenario=="lag_p50" && r.volume_ml==volume))
    minimum(speeds)<=rpm<=maximum(speeds) || error("RPM outside available CFD range")
    lo=maximum(filter(<=(rpm),speeds));hi=minimum(filter(>=(rpm),speeds))
    w=lo==hi ? 0. : (rpm-lo)/(hi-lo)
    a,b=grid(50,lo),grid(50,hi)
    c,d=grid(edr_percentile,lo),grid(edr_percentile,hi)
    linear_exposure([r.size_um for r in a],
        [(1-w)*x.shear_pa+w*y.shear_pa for (x,y) in zip(a,b)],
        [(1-w)*x.edr_m2_s3+w*y.edr_m2_s3 for (x,y) in zip(c,d)])
end

function replace_exposure(ctx; edr_percentile=75, settings=nothing)
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    qs=copy(ctx.conditions)
    for (c,q) in ctx.conditions
        volume,rpm=if settings===nothing
            r=first(filter(r->r.condition==c,series));(r.volume_ml,r.rpm)
        else
            settings[c]
        end
        qs[c]=merge(q,(lookup=mixed_lookup(volume,rpm;edr_percentile),))
    end
    merge(ctx,(conditions=qs,))
end

function source_contexts(fit,splitfit)
    cfgfile=joinpath(fit,"calibration.toml");prep=SplitSelection.prepare(cfgfile)
    cfg=prep.cfg
    cfg["exposure_percentile"]==50 || error("Expected original P50 reference")
    length(prep.split.training)==19 && length(prep.split.selection)==12 || error("Split counts changed")
    Set(cfg["calibration_cohort"]["series_ids"])==Set(prep.split.training) || error("Reference is not training-only")
    Set(cfg["calibration_cohort"]["selection_ids"])==Set(prep.split.selection) || error("Selection culture list differs")
    saved=collect(CSV.File(joinpath(splitfit,"split.csv")))
    Set(String(r.series_id) for r in saved if r.role=="calibration")==Set(prep.split.training) || error("Training split differs")
    Set(String(r.series_id) for r in saved if r.role=="selection")==Set(prep.split.selection) || error("Selection split differs")
    cfg["bounds"]["K_G"][2]==100 && cfg["bounds"]["K_Lac"][2]==100 || error("Approved bounds differ")
    length(free_names(prep.train,MODEL))==11 || error("Expected all 11 Model 6 parameters")
    p=read_parameters(joinpath(fit,"model5/selected_parameters.csv"),MODEL)
    p==read_parameters(joinpath(splitfit,"model5/selected_parameters.csv"),MODEL) || error("Reference parameters differ from training optimum")
    nll=gaussian_nll(prep.train.y,prediction(prep.train,MODEL,p),prep.train.sigma)
    recorded=only(CSV.File(joinpath(splitfit,"model5/summary.csv"))).nll
    isapprox(nll,recorded;atol=1e-6,rtol=0) || error("Reference training likelihood differs")
    prep,p
end

function media_context(cfg,protocol)
    raw=filter(r->r.protocol_percent==protocol && r.status=="measured",collect(CSV.File(joinpath(ROOT,"data/processed/media_exchange_measurements.csv"))))
    initial(s)=mean(Float64[r.value for r in raw if r.day==0 && r.state==s])
    obs=NamedTuple[]
    for day in 1:6,s in STATES
        vals=Float64[r.value for r in raw if r.day==day && r.state==s]
        isempty(vals) && continue
        push!(obs,(condition=protocol,day,state=s,mean=mean(vals),sd=length(vals)>1 ? std(vals) : 0.,n_vessels=length(vals)))
    end
    glucose=initial("glucose")
    q=(u0=[initial("cells"),glucose,initial("lactate"),Float64(cfg["initial_aggregate_um"])],glucose,
       lookup=mixed_lookup(100,40;edr_percentile=50),times=collect(24.:24.:144.),
       feeds=[(time_h=24.0d,fraction=d==2 ? .5 : protocol/100) for d in 2:5])
    (;config=cfg,conditions=Dict(protocol=>q),observations=obs,y=Float64[r.mean for r in obs],sigma=ones(length(obs)))
end

function evaluations(prep)
    cfg=prep.cfg
    allids=vcat(prep.split.training,prep.split.selection)
    cohorts=Tuple{String,Any,Any}[("training60",prep.train,nothing),("selection40",prep.held,nothing),
        ("all100_descriptive",context(;ids=allids,cfg),nothing),
        ("validation500",context(;conditions=(4,),cfg),nothing)]
    for protocol in (20,50)
        push!(cohorts,("media$protocol",media_context(cfg,protocol),Dict(protocol=>(100,40))))
    end
    raw=collect(CSV.File(joinpath(ROOT,"data/processed/cell_line_measurements.csv")))
    for line in ("585A1","SCTi003-A")
        rpm=only(unique(r.rpm for r in raw if r.cell_line==line && r.status=="measured"))
        push!(cohorts,("ipsc_$(line)_direct_transfer",cell_context(line,cfg),Dict(1=>(100,rpm))))
    end
    cohorts
end

function export_evaluation!(out,cohorts,tag,p,percentile)
    metrics=NamedTuple[];points=NamedTuple[];curves=NamedTuple[]
    for (cohort,base,settings) in cohorts
        ctx=replace_exposure(base;edr_percentile=percentile,settings)
        yp=prediction(ctx,MODEL,p)
        for c in (0,sort(collect(keys(ctx.conditions)))...),s in (STATES...,"all")
            ix=findall(r->(c==0||r.condition==c)&&(s=="all"||r.state==s)&&r.mean!=0,ctx.observations)
            isempty(ix) && continue
            push!(metrics,(formulation=tag,cohort,condition=c,state=s,n=length(ix),mrpe_percent=100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix]))))
        end
        append!(points,[(formulation=tag,cohort,condition=r.condition,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,predicted=yp[i]) for (i,r) in enumerate(ctx.observations)])
        append!(curves,[merge((formulation=tag,cohort),r) for r in export_curves(ctx,MODEL,p)])
    end
    slug=percentile==50 ? "reference" : "mixed"
    CSV.write(joinpath(out,"$(slug)_predictions.csv"),points)
    CSV.write(joinpath(out,"$(slug)_curves.csv"),curves)
    metrics
end

function write_report(out,comparison,scores)
    open(joinpath(out,"index.html"),"w") do io
        print(io,"<!doctype html><meta charset='utf-8'><title>Mixed exposure Model 6 comparison</title><style>body{font:17px Arial;max-width:1150px;margin:35px auto;color:#17202a}table{border-collapse:collapse;width:100%;margin:20px 0}th,td{padding:9px;border-bottom:1px solid #ccc;text-align:right}th:first-child,td:first-child{text-align:left}</style><h1>Model 6 exposure comparison</h1><p>Both formulations use the same 19 training cultures. The current P50/P50 training estimate is reused. The P50-shear/P75-EDR estimate minimizes the training likelihood across the new searches. Selection data and external experiments do not contribute to fitting.</p><h2>Selection BIC</h2><table><tr><th>Formulation</th><th>Parameters</th><th>Selection NLL</th><th>Selection BIC</th></tr>")
        for r in scores
            @printf(io,"<tr><td>%s</td><td>%d</td><td>%.5f</td><td>%.5f</td></tr>",r.formulation,r.k,r.selection_nll,r.selection_bic)
        end
        print(io,"</table><h2>MRPE by data subset and state</h2><p>Negative differences favour mixed exposure. Overall MRPE pools eligible condition-day-state observations, as in the current analysis. The all-100 mL row is descriptive and overlaps training and selection. iPSC rows are direct predictions with the ESI-017 parameters and are not the separately recalibrated cell-line comparisons.</p><table><tr><th>Subset</th><th>State</th><th>P50/P50 (%)</th><th>P50/P75 (%)</th><th>Difference (percentage points)</th></tr>")
        for r in comparison
            r.condition==0 || continue
            @printf(io,"<tr><td>%s</td><td>%s</td><td>%.3f</td><td>%.3f</td><td>%+.3f</td></tr>",r.cohort,r.state,r.reference_mrpe,r.mixed_mrpe,r.difference_pp)
        end
        print(io,"</table><p>The CSV files retain full numerical precision and include errors by condition. Measured day-zero initial conditions and feeding schedules are unchanged.</p>")
    end
end

function run_mixed(fit,splitfit,out)
    fit,splitfit,out=abspath.((fit,splitfit,out))
    starts=parse(Int,get(ENV,"HPSC_MIXED_STARTS","100"))
    budget=parse(Int,get(ENV,"HPSC_MIXED_EVALUATIONS","20000"))
    iterations=parse(Int,get(ENV,"HPSC_MIXED_ITERATIONS","3000"))
    starts>0 && budget>=100 && iterations>0 || error("Invalid budgets")
    prep,p0=source_contexts(fit,splitfit)
    extra=[@__FILE__,joinpath(@__DIR__,"run_mixed_exposure.jl"),joinpath(fit,"calibration.toml"),joinpath(fit,"model5/selected_parameters.csv"),joinpath(splitfit,"split.csv"),joinpath(splitfit,"model5/selected_parameters.csv"),joinpath(splitfit,"model5/summary.csv")]
    lock_inputs(out,Dict("analysis"=>"model6_p50shear_p75edr_v1","starts"=>starts,"evaluations"=>budget,"iterations"=>iterations,"julia"=>string(VERSION),"selection_used_for_fitting"=>false),extra)
    progress(s)=open(io->println(io,Dates.now()," ",s),joinpath(out,"progress.log"),"a")
    CSV.write(joinpath(out,"split.csv"),prep.split.assignments)
    save_parameters(joinpath(out,"reference_parameters.csv"),p0)
    progress("START mixed_fit")
    train=replace_exposure(prep.train)
    p,summary=fit_context(train,MODEL,joinpath(out,"mixed_fit");starts,budget,iterations)
    progress("FINISHED mixed_fit")
    progress("START evaluation")
    scores=NamedTuple[]
    for (tag,par,percentile) in ((BASELINE,p0,50),(MIXED,p,75))
        held=replace_exposure(prep.held;edr_percentile=percentile)
        nll=gaussian_nll(held.y,prediction(held,MODEL,par),held.sigma)
        k=length(free_names(held,MODEL))
        push!(scores,(formulation=tag,n=length(held.y),k,selection_nll=nll,selection_bic=bic(nll,k,length(held.y))))
    end
    cohorts=evaluations(prep)
    reference=export_evaluation!(out,cohorts,BASELINE,p0,50)
    mixed=export_evaluation!(out,cohorts,MIXED,p,75)
    key(r)=(r.cohort,r.condition,r.state)
    lookup=Dict(key(r)=>r for r in reference)
    comparison=[(cohort=r.cohort,condition=r.condition,state=r.state,n=r.n,reference_mrpe=lookup[key(r)].mrpe_percent,mixed_mrpe=r.mrpe_percent,difference_pp=r.mrpe_percent-lookup[key(r)].mrpe_percent) for r in mixed]
    CSV.write(joinpath(out,"selection_scores.csv"),scores)
    CSV.write(joinpath(out,"state_errors.csv"),vcat(reference,mixed))
    CSV.write(joinpath(out,"comparison.csv"),comparison)
    write_report(out,comparison,scores)
    progress("FINISHED evaluation")
    write(joinpath(out,"COMPLETE"),"All mixed-exposure fits and frozen-parameter comparisons completed.\n")
    println("Comparison saved: ",joinpath(out,"index.html"))
end
end
