using CSV, SHA, Statistics, TOML, Test
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
length(ARGS)==2 || error("Usage: predict_cell_lines.jl CALIBRATION_DIRECTORY OUTPUT_DIRECTORY")
fitdir,out=abspath.(ARGS);ispath(out) && error("Output already exists")
lock=only(CSV.File(joinpath(fitdir,"selection_lock.csv")))
Int(lock.model)==5 || error("Expected selected Model 5")
parameterfile=joinpath(fitdir,"model5/selected_parameters.csv")
bytes2hex(sha256(read(parameterfile)))==lock.parameter_sha256 || error("Selection hash mismatch")
cfg=TOML.parsefile(joinpath(fitdir,"calibration.toml"))
get(cfg,"lactate_inhibition",true) || error("Expected active lactate inhibition")
rawfile=joinpath(ROOT,"data/processed/cell_line_measurements.csv")
exposurefile=joinpath(ROOT,"data/processed/lagrangian_exposure.csv")
raw=collect(CSV.File(rawfile)); cfd=collect(CSV.File(exposurefile))
points=NamedTuple[];curves=NamedTuple[];grids=NamedTuple[];setups=NamedTuple[]
@testset "Frozen cell-line transfer" begin
    for line in ("585A1","SCTi003-A")
        rows=filter(r->r.cell_line==line && r.status=="measured",raw)
        @test length(unique(r.replicate for r in rows))==3
        # Verify the archived import flag, not the biological treatment.
        # 585A1 experimentally contained DS; see data/processed/README.md.
        @test all(r->r.dextran_sulfate=="absent",rows)
        rpm=only(unique(r.rpm for r in rows)); days=Int(maximum(r.day for r in rows))
        grid40=sort(filter(r->r.volume_ml==100 && r.rpm==40 && r.percentile==cfg["exposure_percentile"],cfd);by=r->r.size_um)
        grid60=sort(filter(r->r.volume_ml==100 && r.rpm==60 && r.percentile==cfg["exposure_percentile"],cfd);by=r->r.size_um)
        @test [r.size_um for r in grid40]==[r.size_um for r in grid60]
        weight=(rpm-40)/20
        shear=[(1-weight)*a.shear_pa+weight*b.shear_pa for (a,b) in zip(grid40,grid60)]
        edr=[(1-weight)*a.edr_m2_s3+weight*b.edr_m2_s3 for (a,b) in zip(grid40,grid60)]
        sizes=[r.size_um for r in grid40];lookup=linear_exposure(sizes,shear,edr)
        for i in eachindex(sizes)
            push!(grids,(cell_line=line,rpm=rpm,size_um=sizes[i],shear_pa=shear[i],edr_m2_s3=edr[i],rpm_interpolated=rpm==50))
        end
        initial(state)=mean(Float64[r.value for r in rows if r.day==0 && r.state==state])
        glucose=initial("glucose");u0=[initial("cells"),glucose,initial("lactate"),Float64(cfg["initial_aggregate_um"])]
        feeds=[(time_h=24.0d,fraction=line=="SCTi003-A" && d==2 ? .5 : .8) for d in 2:days-1]
        push!(setups,(cell_line=line,initial_cells_ml=u0[1],initial_glucose_mM=glucose,initial_lactate_mM=u0[3],initial_aggregate_um=u0[4],fresh_glucose_mM=glucose))
        times=collect(0.:.25:24days)
        for id in (5,0)
            model=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),model)
            sim=simulate(model,p,u0,times,feeds;fresh_glucose=glucose,exposure=lookup,reltol=1e-11,abstol=1e-12)
            tight=simulate(model,p,u0,collect(24.:24.:24days),feeds;fresh_glucose=glucose,exposure=lookup,reltol=1e-12,abstol=1e-13)
            @test all(isfinite,sim.predictions) && minimum(sim.predictions)>=0
            for (j,t) in enumerate(times)
                for (k,state) in enumerate(STATES)
                    push!(curves,(model=id,cell_line=line,day=t/24,state=state,side="before_or_between_exchange",value=sim.predictions[k,j]))
                end
                e=findfirst(e->e.time_h==t,sim.events)
                e===nothing && continue
                @test sim.events[e].before ≈ sim.predictions[:,j]
                @test sim.events[e].after ≈ exchange(sim.events[e].before,sim.events[e].fraction,glucose)
                for (k,state) in enumerate(STATES)
                    push!(curves,(model=id,cell_line=line,day=t/24,state=state,side="after_exchange",value=sim.events[e].after[k]))
                end
            end
            for (k,state) in enumerate(STATES),day in 1:days
                vals=Float64[r.value for r in rows if r.day==day && r.state==state]
                isempty(vals) && continue
                observed=mean(vals);pred=sim.predictions[k,findfirst(==(24.0day),times)]
                @test abs(pred-tight.predictions[k,day])/max(abs(observed),.5)<1e-5
                push!(points,(cohort="$line frozen transfer",model=id,cell_line=line,day=day,state=state,observed=observed,sd=length(vals)>1 ? std(vals) : missing,n_vessels=length(vals),predicted=pred))
            end
        end
    end
end
mkpath(out)
CSV.write(joinpath(out,"predictions.csv"),points)
CSV.write(joinpath(out,"trajectories.csv"),curves)
CSV.write(joinpath(out,"exposure_grid.csv"),grids)
CSV.write(joinpath(out,"initial_conditions.csv"),setups)
cp(rawfile,joinpath(out,"measurements.csv"))
cp(joinpath(ROOT,"data/processed/cell_line_source.json"),joinpath(out,"measurement_source.json"))
files=[rawfile,exposurefile,parameterfile,joinpath(fitdir,"model0/selected_parameters.csv"),joinpath(fitdir,"selection_lock.csv"),joinpath(fitdir,"calibration.toml"),@__FILE__,joinpath(ROOT,"src/Proliferation.jl"),joinpath(ROOT,"tools/import_cell_lines.py")]
CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(file,ROOT),sha256=bytes2hex(sha256(read(file)))) for file in files])
open(joinpath(out,"analysis.toml"),"w") do io
    TOML.print(io,Dict("source_fit"=>fitdir,"parameters_reestimated"=>false,"cohort"=>"585A1 C9/C13/C21; SCTi003-A SCT3/SCT5/SCT7",
    "treatment_annotation"=>"Author confirmed dextran sulfate in 585A1. The archived no-DS import flag is superseded; there is no explicit DS term in the model. See publication/culture_metadata_corrections.csv.",
    "exposure"=>"585A1 50 rpm: linear interpolation between 40/60 rpm CFD at each size, following original script; SCTi003-A direct 40 rpm; size endpoints clamped",
    "initial_glucose"=>"Measured cell-line mean on day zero, replacing original rounded 22.5 mM",
    "fresh_glucose"=>"Assumed equal to measured initial mean, not independently assayed",
    "initial_aggregate"=>"Inherited calibration assumption; absent day-zero aggregate measurements not imputed into observations",
    "sampling"=>"Pre-exchange; 585A1 80% days 2/3, SCTi003-A 50% day 2 and 80% days 3/4",
    "scope"=>"Frozen ESI-017 parameter transfer only. Separate cell-line recalibration remains a distinct analysis, not validation."))
end
write(joinpath(out,"COMPLETE"),"Frozen cell-line prediction checks passed.\n")
