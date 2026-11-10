module PooledValidation
using CSV, TOML, Random, Statistics
include("CrossValidation.jl")
using .CrossValidation: held_weights, metrics
using .CrossValidation.SplitSelection
using .SplitSelection.Workflow
using .SplitSelection.Workflow.Calibration.Proliferation
export plans, prepare_fold, held_weights, metrics, SplitSelection, pooled_context, inventory, exclusions

function inventory()
    base=[(series_id=String(r.series_id),condition=Int(r.condition),volume_ml=Int(r.volume_ml),rpm=Int(r.rpm),seeding=Float64(r.seeding_cells_ml),label=r.volume_ml==500 ? "500mL" : "100mL_condition$(Dict(1=>1,2=>2,3=>3,5=>4,6=>5)[r.condition])") for r in CSV.File(joinpath(ROOT,"data/processed/series.csv"))]
    media=collect(CSV.File(joinpath(ROOT,"data/processed/media_exchange_measurements.csv")))
    for p in (20,50), id in sort(unique(String(r.replicate) for r in media if r.protocol_percent==p))
        push!(base,(series_id="media$(p)_$id",condition=p==20 ? 7 : 8,volume_ml=100,rpm=40,seeding=80000.,label="media_$p"))
    end
    base
end
function measurements()
    rows=[(series_id=String(r.series_id),condition=Int(r.condition),day=Float64(r.day),state=String(r.state),value=ismissing(r.value) ? NaN : Float64(r.value),status=String(r.status)) for r in CSV.File(joinpath(ROOT,"data/processed/measurements.csv"))]
    for r in CSV.File(joinpath(ROOT,"data/processed/media_exchange_measurements.csv"))
        push!(rows,(series_id="media$(r.protocol_percent)_$(r.replicate)",condition=r.protocol_percent==20 ? 7 : 8,day=Float64(r.day),state=String(r.state),value=ismissing(r.value) ? NaN : Float64(r.value),status=String(r.status)))
    end
    rows
end
function pooled_context(ids,cfg;fallback_glucose=Dict())
    info=filter(r->r.series_id in ids,inventory()); rows=filter(r->r.series_id in ids && r.status=="measured",measurements())
    obs=NamedTuple[]; setups=Dict()
    exposure=collect(CSV.File(joinpath(ROOT,"data/processed/lagrangian_exposure.csv")))
    feeds=collect(CSV.File(joinpath(ROOT,"data/processed/feeding.csv")))
    for c in sort(unique(r.condition for r in info))
        sub=filter(r->r.condition==c,rows); q=first(filter(r->r.condition==c,info))
        for day in sort(unique(r.day for r in sub if r.day>0)), state in STATES
            vals=[r.value for r in sub if r.day==day && r.state==state]
            isempty(vals)&&continue
            push!(obs,(condition=c,day,state,mean=mean(vals),sd=length(vals)>1 ? std(vals) : 0.,n_vessels=length(vals)))
        end
        initial(state,default)=begin
            v=[r.value for r in sub if r.day==0 && r.state==state]; isempty(v) ? default : mean(v)
        end
        glucose=initial("glucose",get(fallback_glucose,c,NaN)); isfinite(glucose)||error("Missing initial glucose for $c")
        grid=sort(filter(r->r.volume_ml==q.volume_ml&&r.rpm==q.rpm&&r.percentile==cfg["exposure_percentile"],exposure);by=r->r.size_um)
        lookup=linear_exposure([r.size_um for r in grid],[r.shear_pa for r in grid],[r.edr_m2_s3 for r in grid])
        times=sort(unique(24r.day for r in obs if r.condition==c))
        events=c<=6 ? [(time_h=Float64(r.time_h),fraction=r.exchange_fraction) for r in feeds if r.condition==c&&r.time_h<=last(times)] : [(time_h=24.0d,fraction=d==2 ? .5 : c==7 ? .2 : .5) for d in 2:5 if 24d<=last(times)]
        u0=[c<=6 ? q.seeding : initial("cells",q.seeding),glucose,c<=6 ? cfg["initial_lactate_mM"] : initial("lactate",cfg["initial_lactate_mM"]),cfg["initial_aggregate_um"]]
        setups[c]=(;u0,glucose,lookup,times,feeds=events)
    end
    y=Float64[r.mean for r in obs]
    sigma=Float64[max(r.sd,r.state in ("glucose","lactate") ? cfg["likelihood"]["metabolite_floor_mM"] : cfg["likelihood"]["relative_floor"]*abs(r.mean)) for r in obs]
    (;config=cfg,observations=obs,conditions=setups,y,sigma)
end
function exclusions(ids)
    inv=inventory(); rows=measurements(); result=NamedTuple[]
    for c in sort(unique(r.condition for r in inv if r.series_id in ids))
        group=sort([r.series_id for r in inv if r.condition==c&&r.series_id in ids])
        terminal=maximum(r.day for r in rows if r.condition==c&&r.state=="cells"&&r.status=="measured")
        values=Dict(r.series_id=>r.value for r in rows if r.series_id in group&&r.state=="cells"&&r.day==terminal&&r.status=="measured")
        eligible=length(group)>=4 && length(values)==length(group)
        ranked=eligible ? sort(group;by=id->(values[id],id)) : String[]
        for id in group
            reason=!eligible ? "retained_insufficient_replicates_or_terminal_measurements" : id==first(ranked) ? "excluded_lowest" : id==last(ranked) ? "excluded_highest" : "retained"
            push!(result,(series_id=id,condition=c,terminal_day=terminal,terminal_cells=get(values,id,NaN),decision=reason))
        end
    end
    result
end
function plans(;design="pooled3",seed=20261014)
    design=="pooled3"||error("Only pooled3 supported")
    inv=inventory(); allids=sort([r.series_id for r in inv]); rng=MersenneTwister(seed)
    folds=[String[] for _ in 1:3]
    offset=0
    for c in sort(unique(r.condition for r in inv))
        group=shuffle(rng,sort([r.series_id for r in inv if r.condition==c]))
        for (i,id) in enumerate(group); push!(folds[mod1(i+offset,3)],id); end
        offset=mod(offset+length(group),3)
    end
    result=NamedTuple[]
    for i in 1:3
        held=sort(folds[i]); training=setdiff(allids,held)
        push!(result,(name="pooled$(i)_full",experiment="pooled_full",training,selection=held))
        dropped=[r.series_id for r in exclusions(training) if startswith(r.decision,"excluded")]
        push!(result,(name="pooled$(i)_trimmed",experiment="pooled_trimmed",training=setdiff(training,dropped),selection=held))
    end
    small=sort([r.series_id for r in inv if r.condition in (1,2,3,5,6)])
    large=sort([r.series_id for r in inv if r.condition==4])
    push!(result,(name="reference100",experiment="reference",training=small,selection=large))
    dropped=[r.series_id for r in exclusions(small) if startswith(r.decision,"excluded")]
    push!(result,(name="trimmed100",experiment="trimmed100",training=setdiff(small,dropped),selection=large))
    result
end
function prepare_fold(cfgfile,plan)
    cfg=TOML.parsefile(cfgfile); isempty(cfg["fixed"])||error("Requires all-free calibration")
    isempty(intersect(plan.training,plan.selection))||error("Replicate leakage")
    train=pooled_context(plan.training,cfg)
    held=pooled_context(plan.selection,cfg;fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
    sigma,weight_source=held_weights(train,held.observations);held=merge(held,(;sigma))
    assignments=[(series_id=r.series_id,source_label=r.label,legacy_condition=r.condition,publication_condition=r.label,role=r.series_id in plan.training ? "calibration" : "selection") for r in inventory() if r.series_id in union(plan.training,plan.selection)]
    split=(training=plan.training,selection=plan.selection,assignments=assignments)
    (;cfg,split,train,held,weight_source)
end
end
