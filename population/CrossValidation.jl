module CrossValidation
using CSV, TOML, Random, Statistics
include("SplitSelection.jl")
using .SplitSelection
using .SplitSelection.Workflow
export plans, prepare_fold, held_weights, metrics, SplitSelection
const CONDITIONS=(1,2,3,5,6)
const PUB=Dict(1=>1,2=>2,3=>3,5=>4,6=>5,4=>6)

function plans(;design="balanced9",seed=20261007)
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    small=sort([String(r.series_id) for r in series if r.volume_ml==100])
    large=sort([String(r.series_id) for r in series if r.volume_ml==500])
    length(small)==31 && length(large)==9 || error("Culture counts changed; review split design")
    result=NamedTuple[]
    push!(result,(name="reference100",experiment="reference",training=small,selection=large))
    for c in CONDITIONS
        held=sort([String(r.series_id) for r in series if r.condition==c])
        push!(result,(name="loco_condition$(PUB[c])",experiment="loco",training=setdiff(small,held),selection=held))
    end
    subsets=Vector{String}[]
    if design=="balanced9"
        shuffled=shuffle(MersenneTwister(seed),large)
        for i in 1:9
            push!(subsets,sort([shuffled[mod1(i+j,9)] for j in 0:3]))
        end
    elseif design=="all126"
        for a in 1:6,b in a+1:7,c in b+1:8,d in c+1:9
            push!(subsets,[large[a],large[b],large[c],large[d]])
        end
    else
        error("Design must be balanced9 or all126")
    end
    for (i,part) in enumerate(subsets)
        push!(result,(name="partial500_$(lpad(i,3,'0'))",experiment="partial500",training=vcat(small,part),selection=setdiff(large,part)))
    end
    result
end

function held_weights(train,observations)
    key(r)=(r.condition,r.day,r.state)
    exact=Dict(key(r)=>s for (r,s) in zip(train.observations,train.sigma))
    # An unseen condition has no matching calibration SD. Use a fixed training-only
    # RMS scale per state, rather than estimating a weight from held-out outcomes.
    pooled=Dict(state=>sqrt(mean(abs2,[s for (r,s) in zip(train.observations,train.sigma) if r.state==state])) for state in STATES)
    ([get(exact,key(r),pooled[r.state]) for r in observations],
     [haskey(exact,key(r)) ? "training_condition_day_state" : "training_state_rms" for r in observations])
end

function prepare_fold(cfgfile,plan)
    cfg=TOML.parsefile(cfgfile)
    isempty(cfg["fixed"]) || error("All fold parameters must be estimated from training cultures")
    isempty(get(cfg,"parameter_ties",Dict())) || error("Unexpected parameter ties")
    isempty(intersect(plan.training,plan.selection)) || error("Culture leakage")
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    cond(ids)=sort(unique(Int(r.condition) for r in series if String(r.series_id) in ids))
    train=context(ids=plan.training,conditions=cond(plan.training),cfg=cfg)
    held=context(ids=plan.selection,conditions=cond(plan.selection),cfg=cfg,
        fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
    sigma,weight_source=held_weights(train,held.observations)
    held=merge(held,(sigma=sigma,))
    assignments=[(series_id=String(r.series_id),source_label=String(r.source_label),legacy_condition=Int(r.condition),
        publication_condition=r.volume_ml==500 ? "500 mL" : "100 mL condition $(PUB[r.condition])",
        role=String(r.series_id) in plan.training ? "calibration" : "selection") for r in series if String(r.series_id) in union(plan.training,plan.selection)]
    split=(training=plan.training,selection=plan.selection,assignments=assignments)
    (;cfg,split,train,held,weight_source)
end

function metrics(ctx,yp)
    rows=NamedTuple[]
    for state in (STATES...,"all")
        ix=findall(r->(state=="all" || r.state==state)&&r.day>0&&r.mean!=0,ctx.observations)
        ni=filter(i->ctx.observations[i].n_vessels>=2,ix)
        push!(rows,(state,n=length(ix),n_noise=length(ni),
            mrpe_percent=isempty(ix) ? NaN : 100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix])),
            experimental_noise_percent=isempty(ni) ? NaN : 100mean([ctx.observations[i].sd/abs(ctx.y[i]) for i in ni])))
    end
    rows
end
end
