using CSV, TOML, SHA, Statistics, Test
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation

function exposure_context(base,scenario,inputs,series)
    conditions=copy(base.conditions)
    for (c,q) in base.conditions
        info=first(filter(r->r.condition==c,series))
        rows=filter(r->r.scenario==scenario && r.volume_ml==info.volume_ml && r.rpm==info.rpm,inputs)
        lookup=if startswith(scenario,"lag_")
            sort!(rows;by=r->r.size_um)
            [r.size_um for r in rows]==[200,400,600,800,1000]||error("Incomplete CFD size grid for $scenario, condition $c")
            linear_exposure([r.size_um for r in rows],[r.shear_pa for r in rows],[r.edr_m2_s3 for r in rows])
        else
            row=only(rows);let tau=row.shear_pa,edr=row.edr_m2_s3; size->(tau,edr);end
        end
        conditions[c]=merge(q,(;lookup))
    end
    merge(base,(;conditions))
end

function main()
    length(ARGS)==2||error("Usage: evaluate_exposure_transfer.jl SAVED_EXPOSURE_FITS OUTPUT")
    fit,out=abspath.(ARGS)
    isfile(joinpath(fit,"COMPLETE"))||error("Source exposure comparison incomplete")
    fit!=out||error("Use a separate evaluation directory")
    cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
    isempty(cfg["fixed"])||error("Expected original all-active exposure fits")
    inputs=collect(CSV.File(joinpath(fit,"exposure_inputs.csv")))
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    ranking=collect(CSV.File(joinpath(fit,"ranking.csv")))
    scenarios=Set(["lag_p50","lag_mean","lag_p25","lag_p75","lag_p90","lag_max","eul_va","eul_p50","eul_p75","eul_p90","eul_max"])
    Set(String(r.scenario) for r in ranking)==scenarios||error("Expected all eleven exposure summaries")
    train=context(;cfg);external=context(;cfg,conditions=(4,))
    obs=collect(CSV.File(joinpath(fit,"observations.csv")))
    @assert [(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in obs]==[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in train.observations]
    files=vcat([@__FILE__,joinpath(fit,"calibration.toml"),joinpath(fit,"ranking.csv"),joinpath(fit,"observations.csv"),joinpath(fit,"exposure_inputs.csv")],
        [joinpath(fit,r.scenario,"selected_parameters.csv") for r in ranking],
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),
        filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    hashes=Dict(abspath(f)=>bytes2hex(sha256(read(f))) for f in files)
    manifest=Dict("source_fit"=>fit,"hashes"=>hashes,"scope"=>"Saved all-100-mL exposure fits; 500-mL retrospective prediction without refitting", "order"=>"Original training-BIC order; no 500-mL selection", "model"=>6)
    if isdir(out)
        TOML.parsefile(joinpath(out,"manifest.toml"))==manifest||error("Evaluation inputs changed; choose a new output directory")
    else
        mkpath(out);open(io->TOML.print(io,manifest),joinpath(out,"manifest.toml"),"w")
    end
    points=NamedTuple[];metrics=NamedTuple[];checks=NamedTuple[]
    @testset "Saved exposure fits and scale-transfer checks" begin
        for r in ranking
            scenario=String(r.scenario)
            p=read_parameters(joinpath(fit,scenario,"selected_parameters.csv"),6)
            for (cohort,base) in (("100 mL calibration",train),("500 mL prediction",external))
                ctx=exposure_context(base,scenario,inputs,series)
                pred=prediction(ctx,6,p);tight=prediction(ctx,6,p;reltol=1e-13,abstol=1e-14)
                change=maximum(abs.(pred-tight)./ctx.sigma)
                @test change<1e-4
                nll=gaussian_nll(ctx.y,pred,ctx.sigma)
                if cohort=="100 mL calibration"
                    @test isapprox(nll,r.nll;atol=1e-6,rtol=0)
                    @test isapprox(bic(nll,r.k,length(ctx.y)),r.bic;atol=2e-6,rtol=0)
                end
                push!(checks,(scenario,cohort,n=length(ctx.y),evaluation_nll=nll,training_bic=r.bic,training_stopping_met=r.stopping_met,max_solver_change_sigma=change))
                for (i,row) in enumerate(ctx.observations)
                    push!(points,(scenario,cohort,condition=row.condition,day=row.day,state=row.state,observed=row.mean,sd=row.sd,n_vessels=row.n_vessels,sigma=ctx.sigma[i],predicted=pred[i]))
                end
                for state in (STATES...,"all")
                    ix=findall(i->state=="all"||ctx.observations[i].state==state,eachindex(ctx.y))
                    nz=filter(i->ctx.y[i]!=0,ix);noiseix=filter(i->ctx.observations[i].n_vessels>=2,nz)
                    push!(metrics,(scenario,cohort,state,n=length(ix),n_mrpe=length(nz),n_noise=length(noiseix),
                        mrpe_percent=isempty(nz) ? missing : 100mean(abs.((pred[nz].-ctx.y[nz])./ctx.y[nz])),
                        experimental_noise_percent=isempty(noiseix) ? missing : 100mean([ctx.observations[i].sd/abs(ctx.y[i]) for i in noiseix]),
                        rmse=state=="all" ? missing : sqrt(mean(abs2,pred[ix].-ctx.y[ix]))))
                end
            end
        end
    end
    CSV.write(joinpath(out,"state_errors.csv"),metrics)
    CSV.write(joinpath(out,"predictions.csv"),points)
    CSV.write(joinpath(out,"verification.csv"),checks)
    cp(joinpath(fit,"ranking.csv"),joinpath(out,"training_ranking.csv");force=true)
    write(joinpath(out,"metric_definitions.txt"),"MRPE = 100 × mean(abs(prediction-observed)/abs(observed)) across nonzero condition-day means. Experimental Noise (%) = 100 × mean(sample SD/abs(observed mean)) for means supported by at least two cultures. Overall pools eligible condition-day-state entries, rather than averaging four state percentages. SD is not SEM; these descriptive metrics do not model repeated-day covariance. Exposure order is inherited from training BIC; 500 mL does not select or refit parameters. The prediction evaluation is retrospective because 500 mL outcomes have already been inspected. Selected-fit convergence flags are retained in verification.csv.\n")
    write(joinpath(out,"COMPLETE"),"Saved fits verified and evaluated without optimization\n")
    println("Exposure scale-transfer evaluation saved: $out")
end
main()
