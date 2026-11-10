using CSV, TOML, SHA, Random, Statistics, LinearAlgebra, Test
import BlackBoxOptim as BBO
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("HybridRefinement.jl")
using .HybridRefinement
include("ExposureResume.jl")
const SCENARIOS=("lag_p50","lag_mean","lag_p25","lag_p75","lag_p90","lag_max","eul_va","eul_p50","eul_p75","eul_p90","eul_max")
function main(; scenarios=SCENARIOS, source_scenario=identity, lookup_factory=nothing,
    inputs_note="Paired marginal shear and EDR summaries, not necessarily simultaneous joint exposures; Eulerian constant, Lagrangian linear size interpolation with endpoint clamping", extra_sources=String[])
    2<=length(ARGS)<=5 || error("Usage: exposure_sensitivity.jl FINAL_CALIBRATION OUTPUT [STARTS=100] [EVALUATIONS=20000] [--resume]")
    resume=length(ARGS)==5 && ARGS[5]=="--resume"
    length(ARGS)==5 && !resume && error("Expected --resume")
    fit,out=abspath.(ARGS[1:2]);starts=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 100; budget=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 20000
    starts>0 && budget>=100 || error("Invalid search budget")
    ispath(out) && !resume && error("Output already exists; use --resume to preserve and reuse completed fits")
    resume && !isdir(out) && error("Resume directory does not exist")
    cfgfile=joinpath(fit,"calibration.toml");cfg=TOML.parsefile(cfgfile)
    isempty(cfg["fixed"]) || error("Expected all-free model selection configuration")
    get(cfg,"lactate_inhibition",true) || error("Expected active lactate inhibition")
    cfg["bounds"]["K_Lac"][2]==100 && cfg["bounds"]["K_G"][2]==100 || error("Expected final 100 mM bounds")
    base=context(cfg=cfg)
    inputsfile=joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")
    inputs=collect(CSV.File(inputsfile));series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    contexts=Dict()
    @testset "Exposure inputs and fixed training design" begin
        @test sort(collect(keys(base.conditions)))==[1,2,3,5,6]
        @test length(base.y)==116
        sourceobs=collect(CSV.File(joinpath(fit,"observations.csv")))
        @test [(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in sourceobs]==[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in base.observations]
        for scenario in scenarios
            conditions=copy(base.conditions)
            for (c,q) in base.conditions
                info=first(filter(r->r.condition==c,series))
                rows=filter(r->r.scenario==source_scenario(scenario) && r.volume_ml==info.volume_ml && r.rpm==info.rpm,inputs)
                lookup=if lookup_factory !== nothing
                    sort!(rows;by=r->r.size_um)
                    lookup_factory(rows,scenario)
                elseif startswith(scenario,"lag_")
                    sort!(rows;by=r->r.size_um)
                    @test [r.size_um for r in rows]==[200,400,600,800,1000]
                    linear_exposure([r.size_um for r in rows],[r.shear_pa for r in rows],[r.edr_m2_s3 for r in rows])
                else
                    r=only(rows);let tau=r.shear_pa,edr=r.edr_m2_s3; size->(tau,edr);end
                end
                if scenario=="lag_p50"
                    for size in (15.,200.,300.,500.,1000.,1200.);@test lookup(size)==q.lookup(size);end
                end
                conditions[c]=merge(q,(;lookup))
            end
            contexts[scenario]=merge(base,(;conditions))
            @test contexts[scenario].y==base.y && contexts[scenario].sigma==base.sigma
        end
    end
    mkpath(out);BLAS.set_num_threads(1)
    reference=only(filter(r->r.model==5,collect(CSV.File(joinpath(fit,"ranking.csv")))))

    settings=Dict("source_fit"=>fit,"model"=>5,"legacy_model"=>6,"n"=>116,"k"=>11,"scenarios"=>collect(scenarios),"starts_per_summary"=>starts,"population_size"=>50,"max_evaluations"=>budget,"local_iterations"=>3000,"seed_policy"=>"20261004 + 5000 + run; same seed schedule as final Model 5 calibration for every summary","threads"=>Threads.nthreads(),"julia"=>string(VERSION),"method"=>"adaptive_de_rand_1_bin_radiuslimited followed by original direct NelderMead","inference"=>"Sensitivity conditional on Model 5; no reselection among Models 0-6","ranking"=>"Minimum achieved training BIC; prefer converged stage only within 1e-6 NLL of minimum","likelihood"=>"Same fixed-sigma Gaussian working likelihood as source calibration; unmodelled repeated-day dependence retained","evaluation"=>"100-mL training only; no 500-mL outcomes used","inputs"=>inputs_note)

    files=[cfgfile,inputsfile,@__FILE__,joinpath(@__DIR__,"HybridRefinement.jl"),joinpath(@__DIR__,"Manifest.toml"),joinpath(fit,"ranking.csv"),joinpath(fit,"observations.csv")]
    append!(files,extra_sources)
    append!(files,filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)))
    append!(files,filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    if resume
        verify_resume_inputs(out,settings,files)
    else
        cp(cfgfile,joinpath(out,"calibration.toml"));cp(inputsfile,joinpath(out,"exposure_inputs.csv"))
        cp(joinpath(ROOT,"data/processed/exposure_summary_sources.json"),joinpath(out,"exposure_sources.json"))
        CSV.write(joinpath(out,"observations.csv"),base.observations)
        CSV.write(joinpath(out,"reference_model5.csv"),[NamedTuple(reference)])
        open(io->TOML.print(io,settings),joinpath(out,"analysis.toml"),"w")
        CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in unique(files)])
    end
    jobs=[(scenario=s,run=i) for s in scenarios for i in 1:starts]
    if resume
        pending=NamedTuple[]; audit=NamedTuple[]
        for job in jobs
            dir=joinpath(out,job.scenario,"run$(job.run)")
            reusable,reason=verify_saved_fit(dir,job,contexts[job.scenario])
            push!(audit,(scenario=job.scenario,run=job.run,reused=reusable,reason=reason))
            if reusable
                isfile(joinpath(dir,"COMPLETE")) || write(joinpath(dir,"COMPLETE"),"Recovered from verified saved summary and parameters\n")
            else
                push!(pending,job)
            end
        end
        CSV.write(joinpath(out,"resume_audit.csv"),audit)
        @info "Resume inventory" verified_completed=length(jobs)-length(pending) remaining=length(pending) total=length(jobs)
        # Preserve interrupted/invalid directories before retrying only those searches.
        for job in pending
            dir=joinpath(out,job.scenario,"run$(job.run)")
            if isdir(dir)
                archive=joinpath(out,"interrupted_attempts",job.scenario);mkpath(archive)
                mv(dir,joinpath(archive,"run$(job.run)_$(time_ns())"))
            end
        end
        for f in ("COMPLETE","failed_starts.csv")
            isfile(joinpath(out,f)) && rm(joinpath(out,f))
        end
        open(io->TOML.print(io,Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads(),"reused"=>length(jobs)-length(pending),"remaining"=>length(pending),"driver_sha256"=>bytes2hex(sha256(read(@__FILE__))))),joinpath(out,"resume_environment.toml"),"w")
    else
        pending=jobs
    end
    Threads.@threads :greedy for job in pending
        dir=joinpath(out,job.scenario,"run$(job.run)");mkpath(dir)
        try
            ctx=contexts[job.scenario];Random.seed!(20261004+5000+job.run)
            residual=Workflow.residual_function(ctx,6);objective=z->sum(abs2,residual(z))/2
            search=BBO.bboptimize(objective;SearchRange=(0.,1.),NumDimensions=11,Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=50,MaxFuncEvals=budget,MaxSteps=0,MaxStepsWithoutProgress=budget+100,MaxTime=0.,TraceMode=:silent)
            z=copy(BBO.best_candidate(search));p=decode(ctx,6,z)
            nll=gaussian_nll(ctx.y,prediction(ctx,6,p),ctx.sigma)
            save_parameters(joinpath(dir,"population_parameters.csv"),p)
            localstages=refine_original(ctx,6,z,joinpath(dir,"local");iterations=3000)
            stages=vcat([(stage="population",nll=nll,stopping_met=false)],localstages)
            minimum_nll=minimum(r.nll for r in stages)
            eligible=filter(r->r.stopping_met && r.nll<=minimum_nll+1e-6,stages)
            best=argmin(r->r.nll,isempty(eligible) ? stages : eligible)
            selected=best.stage=="population" ? p : read_parameters(joinpath(dir,"local/polished_parameters.csv"),6)
            actual=gaussian_nll(ctx.y,prediction(ctx,6,selected),ctx.sigma)
            isapprox(actual,best.nll;atol=1e-6,rtol=0) || error("Saved likelihood mismatch")
            coords=encode(ctx,6,selected)
            boundaries=join([string(n)*":"*(x<1e-4 ? "lower" : "upper") for (n,x) in zip(free_names(ctx,6),coords) if x<1e-4 || x>1-1e-4],";")
            save_parameters(joinpath(dir,"selected_parameters.csv"),selected)
            CSV.write(joinpath(dir,"stages.csv"),stages)
            CSV.write(joinpath(dir,"summary.csv"),[(scenario=job.scenario,run=job.run,nll=actual,bic=bic(actual,11,116),stopping_met=best.stopping_met,stage=best.stage,boundaries=boundaries,evaluations=BBO.f_calls(search),population_stop=string(BBO.stop_reason(search)))])
            write(joinpath(dir,"COMPLETE"),"Hybrid fit saved\n")
            @info "Exposure fit complete" scenario=job.scenario run=job.run bic=bic(actual,11,116)
        catch e
            e isa InterruptException && rethrow();write(joinpath(dir,"FAILED"),sprint(showerror,e))
        end
    end
    runs=NamedTuple[];failed=NamedTuple[]
    for job in jobs
        dir=joinpath(out,job.scenario,"run$(job.run)")
        if isfile(joinpath(dir,"FAILED"));push!(failed,(scenario=job.scenario,run=job.run,error=read(joinpath(dir,"FAILED"),String)))
        else;push!(runs,NamedTuple(only(CSV.File(joinpath(dir,"summary.csv")))));end
    end
    CSV.write(joinpath(out,"runs.csv"),runs)
    if !isempty(failed);CSV.write(joinpath(out,"failed_starts.csv"),failed);error("Failed starts retained; ranking not finalized");end
    ranking=NamedTuple[];points=NamedTuple[];errors=NamedTuple[]
    @testset "Selected exposure fits" begin
        for scenario in scenarios
            ctx=contexts[scenario];rows=filter(r->r.scenario==scenario,runs)
            minimum_nll=minimum(r.nll for r in rows)
            eligible=filter(r->r.stopping_met && r.nll<=minimum_nll+1e-6,rows)
            best=argmin(r->r.nll,isempty(eligible) ? rows : eligible)
            p=read_parameters(joinpath(out,scenario,"run$(best.run)/selected_parameters.csv"),6)
            pred=prediction(ctx,6,p);tight=prediction(ctx,6,p;reltol=1e-13,abstol=1e-14)
            solver_change=maximum(abs.(pred-tight)./ctx.sigma)
            @test solver_change<1e-4
            @test isapprox(gaussian_nll(ctx.y,pred,ctx.sigma),best.nll;atol=1e-6,rtol=0)
            save_parameters(joinpath(out,scenario,"selected_parameters.csv"),p)
            push!(ranking,merge(best,(n=116,k=11,minimum_nll=minimum_nll,converged_starts=count(r->r.stopping_met,rows),total_starts=starts,max_solver_change_sigma=solver_change,delta_reference_p50=best.bic-reference.bic,)))
            for (i,r) in enumerate(ctx.observations)
                push!(points,(scenario=scenario,condition=r.condition,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,sigma=ctx.sigma[i],predicted=pred[i]))
            end
            for state in (STATES...,"all")
                ix=findall(r->state=="all" || r.state==state,ctx.observations); nz=filter(i->ctx.y[i]!=0,ix);d=pred[ix].-ctx.y[ix]
                push!(errors,(scenario=scenario,state=state,n=length(ix),n_mrpe=length(nz),mrpe_percent=100mean(abs.((pred[nz].-ctx.y[nz])./ctx.y[nz])),experimental_noise_cv_percent=100mean([ctx.observations[i].sd/abs(ctx.y[i]) for i in nz if ctx.observations[i].n_vessels>=2]),rmse=state=="all" ? missing : sqrt(mean(abs2,d)),bias=state=="all" ? missing : mean(d)))
            end
        end
    end
    sort!(ranking;by=r->r.bic)
    CSV.write(joinpath(out,"ranking.csv"),[merge(r,(delta_bic=r.bic-first(ranking).bic,)) for r in ranking])
    CSV.write(joinpath(out,"predictions.csv"),points);CSV.write(joinpath(out,"state_errors.csv"),errors)
    write(joinpath(out,"COMPLETE"),"Training-only exposure comparison verified; check convergence and bounds\n")
end
if abspath(PROGRAM_FILE)==@__FILE__; main(); end
