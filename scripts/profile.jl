using CSV, TOML, LinearAlgebra
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
length(ARGS)>=3 || error("Usage: profile.jl CALIBRATION_DIRECTORY PUBLIC_MODEL PARAMETER [GRID_POINTS]")
fitdir=abspath(ARGS[1]);id=parse(Int,ARGS[2]);name=Symbol(ARGS[3]);npoints=length(ARGS)==4 ? parse(Int,ARGS[4]) : 21
npoints>=3 || error("At least three grid points required")
m=legacy_model(id);cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml"))
isempty(cfg["fixed"]) || error("This profile command expects the all-free fit")
name in parameter_names(m) || error("Parameter not in this model")
fullctx=context(cfg=cfg)
savedobs=collect(CSV.File(joinpath(fitdir,"observations.csv")))
length(savedobs)==length(fullctx.y) && all(isapprox(r.mean,fullctx.y[i];rtol=1e-12) for (i,r) in enumerate(savedobs)) || error("Profile must use the same full 100-mL cohort as its fit; split fits are not supported here")
p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),m)
lo,hi=cfg["bounds"][string(name)]
grid=sort(unique(vcat(exp.(range(log(lo),log(hi);length=npoints)),getproperty(p,name))))
out=joinpath(fitdir,"profiles/model$id/$(name)");ispath(out)&&error("Profile already exists");mkpath(out)
BLAS.set_num_threads(1)
Threads.@threads :dynamic for i in eachindex(grid)
    localcfg=deepcopy(cfg);localcfg["fixed"][string(name)]=grid[i]
    ctx=context(cfg=localcfg)
    z=encode(ctx,m,p)
    fit_start(ctx,m,z,joinpath(out,"point$i");polish_iterations=5000,tight=true)
end
baseline=only(CSV.File(joinpath(fitdir,"model$id/summary.csv"))).nll
rows=NamedTuple[]
for i in eachindex(grid)
    stages=sort(collect(CSV.File(joinpath(out,"point$i/stages.csv")));by=r->r.nll)
    best=first(stages)
    push!(rows,(parameter=string(name),fixed_value=grid[i],nll=best.nll,twice_delta_nll=2(best.nll-baseline),stopping_met=best.stopping_met,stage=String(best.stage)))
end
CSV.write(joinpath(out,"profile.csv"),rows)
println("Profile saved. Inspect stopping flags and improved minima before reporting intervals; flat/bound-truncated profiles are not finite confidence intervals.")
