using CSV, SHA, Statistics, TOML, Test
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation

length(ARGS)==2 || error("Usage: predict_media_exchange.jl CALIBRATION_DIRECTORY OUTPUT_DIRECTORY")
fitdir,out=abspath.(ARGS)
ispath(out) && error("Output already exists: $out")
lock=only(CSV.File(joinpath(fitdir,"selection_lock.csv")))
Int(lock.model)==5 || error("This evaluation was specified for selected Model 5")
parameterfile=joinpath(fitdir,"model5/selected_parameters.csv")
bytes2hex(sha256(read(parameterfile)))==lock.parameter_sha256 || error("Selected parameters changed after selection")
cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml"))
get(cfg,"lactate_inhibition",true) || error("Primary evaluation requires active lactate inhibition")
rawfile=joinpath(ROOT,"data/processed/media_exchange_measurements.csv")
exposurefile=joinpath(ROOT,"data/processed/lagrangian_exposure.csv")
raw=collect(CSV.File(rawfile))
grid=sort(filter(r->r.volume_ml==100 && r.rpm==40 && r.percentile==cfg["exposure_percentile"],collect(CSV.File(exposurefile)));by=r->r.size_um)
lookup=linear_exposure([r.size_um for r in grid],[r.shear_pa for r in grid],[r.edr_m2_s3 for r in grid])
points=NamedTuple[]; curves=NamedTuple[]; errors=NamedTuple[]; events=NamedTuple[]; setups=NamedTuple[]
@testset "Frozen media-exchange evaluation" begin
    @test sort(unique(r.replicate for r in raw if r.protocol_percent==20)) == ["L85","L87","L89"]
    @test sort(unique(r.replicate for r in raw if r.protocol_percent==50)) == ["L84","L86","L88"]
    for protocol in (20,50)
        measured=filter(r->r.protocol_percent==protocol && r.status=="measured",raw)
        initial(state)=mean(Float64[r.value for r in measured if r.day==0 && r.state==state])
        glucose=initial("glucose")
        u0=[initial("cells"),glucose,initial("lactate"),Float64(cfg["initial_aggregate_um"])]
        feeds=[(time_h=24.0day,fraction=day==2 ? .5 : protocol/100) for day in 2:5]
        push!(setups,(protocol_percent=protocol,volume_ml=100,rpm=40,initial_cells_ml=u0[1],initial_glucose_mM=glucose,initial_lactate_mM=u0[3],initial_aggregate_um=u0[4],fresh_glucose_mM=glucose))
        times=collect(0.:.25:144.)
        for id in (5,0)
            model=legacy_model(id)
            p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
            sim=simulate(model,p,u0,times,feeds;fresh_glucose=glucose,exposure=lookup,reltol=1e-11,abstol=1e-12)
            tight=simulate(model,p,u0,collect(24.:24.:144.),feeds;fresh_glucose=glucose,exposure=lookup,reltol=1e-12,abstol=1e-13)
            @test all(isfinite,sim.predictions) && minimum(sim.predictions)>=0
            @test length(sim.events)==4
            for (j,t) in enumerate(times)
                for (k,state) in enumerate(STATES)
                    push!(curves,(model=id,protocol_percent=protocol,day=t/24,state=state,side="before_or_between_exchange",value=sim.predictions[k,j]))
                end
                e=findfirst(e->e.time_h==t,sim.events)
                e===nothing && continue
                event=sim.events[e]
                @test event.before ≈ sim.predictions[:,j]
                @test event.before[[1,4]]==event.after[[1,4]]
                @test event.after[2] ≈ (1-event.fraction)*event.before[2]+event.fraction*glucose
                @test event.after[3] ≈ (1-event.fraction)*event.before[3]
                for (k,state) in enumerate(STATES)
                    push!(curves,(model=id,protocol_percent=protocol,day=t/24,state=state,side="after_exchange",value=event.after[k]))
                    push!(events,(model=id,protocol_percent=protocol,day=t/24,state=state,fraction=event.fraction,before=event.before[k],after=event.after[k]))
                end
            end
            for (k,state) in enumerate(STATES), day in 1:6
                vals=Float64[r.value for r in measured if r.day==day && r.state==state]
                isempty(vals) && continue
                mu=mean(vals); sd=length(vals)>1 ? std(vals) : missing
                pred=sim.predictions[k,findfirst(==(24.0day),times)]
                scale=max(abs(mu),state in ("glucose","lactate") ? .5 : 1.)
                @test abs(pred-tight.predictions[k,day])/scale < 1e-5
                push!(points,(model=id,protocol_percent=protocol,day=day,state=state,observed=mu,sd=sd,n_vessels=length(vals),predicted=pred))
            end
            for state in STATES
                rows=filter(r->r.model==id && r.protocol_percent==protocol && r.state==state,points)
                differences=[r.predicted-r.observed for r in rows]
                relative=[abs((r.predicted-r.observed)/r.observed) for r in rows if r.observed!=0]
                push!(errors,(model=id,protocol_percent=protocol,state=state,n=length(rows),n_mrpe=length(relative),mrpe_percent=isempty(relative) ? missing : 100mean(relative),rmse=sqrt(mean(abs2,differences)),bias=mean(differences)))
            end
        end
    end
end
mkpath(out)
CSV.write(joinpath(out,"predictions.csv"),points)
CSV.write(joinpath(out,"trajectories.csv"),curves)
CSV.write(joinpath(out,"state_errors.csv"),errors)
CSV.write(joinpath(out,"exchange_events.csv"),events)
CSV.write(joinpath(out,"initial_conditions.csv"),setups)
cp(rawfile,joinpath(out,"measurements.csv"))
cp(joinpath(ROOT,"data/processed/media_exchange_source.json"),joinpath(out,"measurement_source.json"))
files=[rawfile,exposurefile,joinpath(fitdir,"selection_lock.csv"),joinpath(fitdir,"calibration.toml"),parameterfile,joinpath(fitdir,"model0/selected_parameters.csv"),@__FILE__,joinpath(ROOT,"src/Proliferation.jl"),joinpath(ROOT,"src/Workflow.jl"),joinpath(ROOT,"tools/import_media_exchange.py")]
CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(file,ROOT),sha256=bytes2hex(sha256(read(file)))) for file in files])
open(joinpath(out,"analysis.toml"),"w") do io
    TOML.print(io,Dict("source_fit"=>fitdir,"selected_model"=>5,"parameters_reestimated"=>false,
        "initial_glucose"=>"Protocol mean of measured day-zero glucose; not fitted",
        "fresh_glucose"=>"Assumed equal to protocol initial glucose; fresh medium not separately measured",
        "initial_aggregate"=>"Inherited calibration assumption; day-zero workbook zeros excluded",
        "feeding"=>"50% day 2; assigned 20% or 50% days 3,4,5; no harvest-day exchange",
        "observations"=>"Days 1-6, before exchange, mean and sample SD of labeled bioreactors",
        "mrpe"=>"100 * mean(abs(prediction - observed_mean)/abs(observed_mean)); zeros excluded and counted",
        "interpretation"=>"Frozen-parameter evaluation of previously inspected protocol data; no predictive intervals estimated"))
end
write(joinpath(out,"COMPLETE"),"Frozen predictions and numerical/event checks completed.\n")
println("Saved media-exchange predictions to $out")
