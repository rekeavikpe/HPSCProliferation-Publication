using CSV, SHA, TOML
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
fitdir=isempty(ARGS) ? joinpath(ROOT,"results/reproduction_calibration") : abspath(ARGS[1])
out=length(ARGS)<2 ? joinpath(ROOT,"reference/figure_data") : abspath(ARGS[2])
mkpath(out)
lock=only(CSV.File(joinpath(fitdir,"selection_lock.csv")));winner=Int(lock.model)
ctx=context(cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml")));points=NamedTuple[];curves=NamedTuple[];starts=NamedTuple[]
for id in unique([0,winner])
    model=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
    for (c,q) in sort(collect(ctx.conditions);by=first)
        times=collect(0.:.25:last(q.times))
        sim=simulate(model,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,reltol=1e-12,abstol=1e-13)
        for (i,t) in enumerate(times)
            for (k,state) in enumerate(STATES)
                push!(curves,(model=id,condition=c,day=t/24,state=state,side="before_or_between_exchange",value=sim.predictions[k,i]))
            end
            event=findfirst(e->e.time_h==t,sim.events)
            if event!==nothing
                for (k,state) in enumerate(STATES)
                    push!(curves,(model=id,condition=c,day=t/24,state=state,side="after_exchange",value=sim.events[event].after[k]))
                end
            end
        end
    end
end
CSV.write(joinpath(out,"training_points.csv"),ctx.observations)
CSV.write(joinpath(out,"training_curves.csv"),curves)
for id in 0:6
    summary=only(CSV.File(joinpath(fitdir,"model$id/summary.csv")))
    for name in sort(filter(n->startswith(n,"run"),readdir(joinpath(fitdir,"model$id"))))
        stages=sort(collect(CSV.File(joinpath(fitdir,"model$id",name,"stages.csv")));by=r->r.nll)
        r=first(stages)
        push!(starts,(model=id,run=name,bic=bic(r.nll,summary.k,summary.n),stopping_met=r.stopping_met,stage=r.stage))
    end
end
CSV.write(joinpath(out,"optimizer_bic.csv"),starts)
for name in ("test_predictions.csv","test_trajectories.csv","test_state_errors.csv")
    cp(joinpath(fitdir,"scale_transfer",name),joinpath(out,name);force=true)
end
cp(joinpath(fitdir,"ranking.csv"),joinpath(out,"training_ranking.csv");force=true)
CSV.write(joinpath(out,"parameter_hashes.csv"),[(model=id,sha256=bytes2hex(sha256(read(joinpath(fitdir,"model$id/selected_parameters.csv"))))) for id in 0:6])
open(joinpath(out,"source.toml"),"w") do io
    TOML.print(io,Dict("source_run"=>basename(fitdir),"selected_model"=>winner,"initial_starts_per_model"=>TOML.parsefile(joinpath(fitdir,"run.toml"))["starts"],"meaning"=>"Saved calibration outputs; consult source convergence and analysis metadata"))
end
println("Exported figure data to $out")
