using CSV, SHA, Statistics, TOML
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
length(ARGS)==1 || error("Usage: predict_scale.jl CALIBRATION_DIRECTORY")
fitdir=abspath(ARGS[1]); lock=only(CSV.File(joinpath(fitdir,"selection_lock.csv")))
winner=Int(lock.model)
parameterfile=joinpath(fitdir,"model$winner/selected_parameters.csv")
bytes2hex(sha256(read(parameterfile)))==lock.parameter_sha256 || error("Selected parameters changed after selection lock")
ctx=context(conditions=(4,),cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml")));q=ctx.conditions[4]
out=joinpath(fitdir,"scale_transfer");ispath(out)&&error("Prediction output already exists");mkpath(out)
points=NamedTuple[];curves=NamedTuple[];errors=NamedTuple[]
for id in unique([winner,0])
    m=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),m)
    pred=prediction(ctx,m,p)
    for (i,r) in enumerate(ctx.observations)
        push!(points,(publication_model=id,condition=4,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,predicted=pred[i]))
    end
    times=collect(0.:.25:last(q.times))
    sim=simulate(m,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,reltol=1e-12,abstol=1e-13)
    for (i,t) in enumerate(times)
        for (k,state) in enumerate(STATES)
            push!(curves,(publication_model=id,day=t/24,state=state,side="before_or_between_exchange",value=sim.predictions[k,i]))
        end
        event=findfirst(e->e.time_h==t,sim.events)
        if event!==nothing
            for (k,state) in enumerate(STATES)
                push!(curves,(publication_model=id,day=t/24,state=state,side="after_exchange",value=sim.events[event].after[k]))
            end
        end
    end
    for state in STATES
        inds=findall(r->r.state==state,ctx.observations);y=ctx.y[inds];yp=pred[inds]
        nz=findall(!iszero,y)
        push!(errors,(model=id,state=state,n=length(y),mrpe_percent=100mean(abs.((yp[nz].-y[nz])./y[nz])),rmse=sqrt(mean((yp.-y).^2)),bias=mean(yp.-y)))
    end
end
CSV.write(joinpath(out,"test_predictions.csv"),points)
CSV.write(joinpath(out,"test_trajectories.csv"),curves)
CSV.write(joinpath(out,"test_state_errors.csv"),errors)
println("500-mL predictions saved; no post-day-zero outcomes used for fitting.")
