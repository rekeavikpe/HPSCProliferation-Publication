module Calibration

using CSV, TOML, Statistics
include("Proliferation.jl")
using .Proliferation
export load_calibration, predict, residuals, decode, encode, free_names, score, Proliferation

const STATES = ("cells", "glucose", "lactate", "aggregate")

function load_calibration(root, config_path)
    config = TOML.parsefile(config_path)
    cohort = config["cohort"]
    cohort in ("available_measurements", "legacy_cohort") || error("Unknown cohort")
    observations = collect(CSV.File(joinpath(root,"data/processed/summary_$(cohort).csv")))
    sort!(observations; by=r -> (r.condition,r.day,findfirst(==(r.state),STATES)))
    measurements = collect(CSV.File(joinpath(root,"data/processed/measurements.csv")))
    series = collect(CSV.File(joinpath(root,"data/processed/series.csv")))
    feeding = collect(CSV.File(joinpath(root,"data/processed/feeding.csv")))
    exposure_rows = collect(CSV.File(joinpath(root,"data/processed/lagrangian_exposure.csv")))
    conditions = Dict()
    for c in 1:6
        info = first(filter(r -> r.condition == c, series))
        glucose = mean([r.value for r in measurements if r.condition == c && r.day == 0 && r.state == "glucose" && r.status == "measured"])
        grid = sort(filter(r -> r.volume_ml == info.volume_ml && r.rpm == info.rpm && r.percentile == config["exposure_percentile"],exposure_rows);by=r->r.size_um)
        lookup = linear_exposure([r.size_um for r in grid],[r.shear_pa for r in grid],[r.edr_m2_s3 for r in grid])
        times = sort(unique(Float64[24r.day for r in observations if r.condition == c]))
        feeds = [(time_h=Float64(r.time_h),fraction=r.exchange_fraction) for r in feeding if r.condition == c]
        u0 = [Float64(info.seeding_cells_ml),glucose,config["initial_lactate_mM"],config["initial_aggregate_um"]]
        conditions[c] = (;u0,glucose,lookup,times,feeds)
    end
    y = Float64[r.mean for r in observations]
    sigma = Float64[max(r.sd,r.state in ("glucose","lactate") ? config["likelihood"]["metabolite_floor_mM"] : config["likelihood"]["relative_floor"]*abs(r.mean)) for r in observations]
    all(s->isfinite(s)&&s>0,sigma) || error("Invalid uncertainty")
    return (;config,observations,conditions,y,sigma)
end

free_names(ctx,model) = filter(n -> !haskey(ctx.config["fixed"],string(n)) && !haskey(get(ctx.config,"parameter_ties",Dict()),string(n)), parameter_names(model))

function decode(ctx,model,z)
    names=free_names(ctx,model)
    length(z)==length(names) || error("Free parameter count mismatch")
    values=Dict{Symbol,Float64}(Symbol(k)=>v for (k,v) in ctx.config["fixed"])
    for (name,x) in zip(names,z)
        lo,hi=ctx.config["bounds"][string(name)]
        values[name]=exp(log(lo)+x*log(hi/lo))
    end
    for (target,source) in get(ctx.config,"parameter_ties",Dict())
        haskey(values,Symbol(target)) && error("Tied parameter cannot also be fixed or free")
        values[Symbol(target)]=values[Symbol(source)]
    end
    return parameters(model,[values[n] for n in parameter_names(model)])
end

function encode(ctx,model,p)
    return [begin
        lo,hi=ctx.config["bounds"][string(n)]
        clamp(log(getproperty(p,n)/lo)/log(hi/lo),0.,1.)
    end for n in free_names(ctx,model)]
end

function predict(ctx,model,p;reltol=ctx.config["solver"]["reltol"],abstol=ctx.config["solver"]["abstol"])
    predictions=zeros(length(ctx.y))
    for c in 1:6
        condition=ctx.conditions[c]
        result=simulate(model,p,condition.u0,condition.times,condition.feeds;
            fresh_glucose=condition.glucose,exposure=condition.lookup,reltol,abstol)
        all(x->isfinite(x)&&x>=0,result.predictions) || error("Invalid predicted state")
        for (i,row) in enumerate(ctx.observations)
            row.condition==c || continue
            predictions[i]=result.predictions[findfirst(==(row.state),STATES),findfirst(==(24row.day),condition.times)]
        end
    end
    predictions
end

residuals(ctx,model,p) = (predict(ctx,model,p).-ctx.y)./ctx.sigma
score(ctx,model,p) = gaussian_nll(ctx.y,predict(ctx,model,p),ctx.sigma)

end
