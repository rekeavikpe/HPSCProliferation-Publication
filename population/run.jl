using CSV, TOML, Random, Statistics, SHA, LinearAlgebra
import BlackBoxOptim as BBO
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("HybridRefinement.jl")
using .HybridRefinement

function main()
    4 <= length(ARGS) <= 7 || error("Usage: run.jl CONFIG OUT REPEATS MAX_EVALUATIONS [MODEL_IDS=0,5] [REFERENCE_BIC_CSV] [--original-hybrid]")
    original_hybrid=length(ARGS)==7
    original_hybrid && ARGS[7]!="--original-hybrid" && error("Unknown refinement option")
    cfgfile,out=abspath.(ARGS[1:2])
    repeats,budget=parse.(Int,ARGS[3:4])
    ids=length(ARGS)>=5 ? parse.(Int,split(ARGS[5],',')) : [0,5]
    repeats>0 && budget>=100 || error("Positive repetitions and at least 100 evaluations required")
    all(id->id in 0:6,ids) && allunique(ids) || error("Invalid model IDs")
    ispath(out) && error("Choose a new output directory")
    cfg=TOML.parsefile(cfgfile);isempty(cfg["fixed"]) || error("Expected all-free calibration")
    ctx=context(cfg=cfg);BLAS.set_num_threads(1)
    mkpath(out)
    effective=deepcopy(cfg)
    if original_hybrid
        merge!(effective["optimization"],Dict("starts"=>repeats,"seed"=>20261004,"local_iterations"=>0,
            "polish_iterations"=>3000,"method"=>"differential evolution then direct NelderMead",
            "note"=>"Original hybrid structure; current likelihood, bounds and normalized-log coordinates retained."))
    end
    open(io->TOML.print(io,effective),joinpath(out,"calibration.toml"),"w")
    CSV.write(joinpath(out,"observations.csv"),ctx.observations)
    reference=Dict{Int,Float64}()
    if length(ARGS)>=6
        ref=abspath(ARGS[6]);cp(ref,joinpath(out,"reference_bic.csv"))
        for row in CSV.File(ref)
            id=Int(row.model);reference[id]=min(get(reference,id,Inf),row.bic)
        end
    end
    settings=Dict("status"=>"running","models"=>ids,"repeats"=>repeats,"max_evaluations"=>budget,
        "population_size"=>50,"seed"=>20261004,"method"=>"adaptive_de_rand_1_bin_radiuslimited",
        "coordinates"=>"normalized log parameters, same as current local calibration",
        "local_iterations"=>original_hybrid ? 0 : 300,"polish_iterations"=>3000,"threads"=>Threads.nthreads(),
        "local_method"=>original_hybrid ? "original direct NelderMead with clipped coordinates" : "least squares then penalty NelderMead",
        "julia"=>string(VERSION),"BlackBoxOptim"=>string(pkgversion(BBO)),
        "reference_note"=>"Reference BIC is descriptive only; never used to initialize or stop searches. Reference cohort/config must match.")
    write_settings()=open(io->TOML.print(io,settings),joinpath(out,"run.toml"),"w")
    write_settings()
    files=vcat([cfgfile,@__FILE__,joinpath(@__DIR__,"HybridRefinement.jl"),joinpath(@__DIR__,"Project.toml"),joinpath(@__DIR__,"Manifest.toml")],
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),
        filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in files])
    jobs=[(id=id,run=r) for id in ids for r in 1:repeats]
    Threads.@threads :greedy for j in jobs
        dir=joinpath(out,"model$(j.id)","run$(j.run)");mkpath(dir)
        try
            model=legacy_model(j.id);seed=20261004+1000j.id+j.run
            Random.seed!(seed)
            residual=Workflow.residual_function(ctx,model)
            objective=z->sum(abs2,residual(z))/2
            start=time()
            result=BBO.bboptimize(objective;SearchRange=(0.,1.),NumDimensions=length(free_names(ctx,model)),
                Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=50,
                MaxFuncEvals=budget,MaxSteps=0,MaxStepsWithoutProgress=budget+100,
                MaxTime=0.,TraceMode=:silent)
            z=BBO.best_candidate(result);p=decode(ctx,model,z)
            nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
            save_parameters(joinpath(dir,"population_parameters.csv"),p)
            open(joinpath(dir,"population.toml"),"w") do io
                TOML.print(io,Dict("seed"=>seed,"evaluations"=>BBO.f_calls(result),
                    "stop_reason"=>BBO.stop_reason(result),"elapsed_seconds"=>time()-start,"nll"=>nll))
            end
            stages=original_hybrid ? refine_original(ctx,model,z,joinpath(dir,"local");iterations=3000) :
                fit_start(ctx,model,z,joinpath(dir,"local");lm_iterations=300,polish_iterations=3000)
            candidates=vcat([(stage="population",nll=nll,stopping_met=false)],stages)
            CSV.write(joinpath(dir,"stages.csv"),candidates)
            best=argmin(r->r.nll,candidates)
            path=best.stage=="population" ? joinpath(dir,"population_parameters.csv") : joinpath(dir,"local",best.stage*"_parameters.csv")
            chosen=read_parameters(path,model)
            actual=gaussian_nll(ctx.y,prediction(ctx,model,chosen),ctx.sigma)
            isapprox(actual,best.nll;atol=1e-6,rtol=0) || error("Saved likelihood mismatch")
            save_parameters(joinpath(dir,"selected_parameters.csv"),chosen)
            k=length(free_names(ctx,model));score=bic(actual,k,length(ctx.y))
            CSV.write(joinpath(dir,"summary.csv"),[(model=j.id,run=j.run,n=length(ctx.y),k=k,
                population_bic=bic(nll,k,length(ctx.y)),bic=score,stopping_met=best.stopping_met,
                stage=best.stage,delta_reference=score-get(reference,j.id,NaN),evaluations=BBO.f_calls(result))])
            write(joinpath(dir,"COMPLETE"),"population and local stages verified\n")
            @info "Population comparison complete" model=j.id run=j.run bic=score
        catch e
            e isa InterruptException && rethrow()
            write(joinpath(dir,"FAILED"),sprint(showerror,e)*"\n")
            @warn "Run failed and retained" model=j.id run=j.run exception=sprint(showerror,e)
        end
    end
    rows=NamedTuple[];failures=NamedTuple[]
    for j in jobs
        dir=joinpath(out,"model$(j.id)","run$(j.run)")
        if isfile(joinpath(dir,"FAILED"))
            push!(failures,(model=j.id,run=j.run,error=read(joinpath(dir,"FAILED"),String)))
        else
            push!(rows,NamedTuple(only(CSV.File(joinpath(dir,"summary.csv")))))
        end
    end
    CSV.write(joinpath(out,"runs.csv"),rows)
    isempty(failures) || CSV.write(joinpath(out,"failures.csv"),failures)
    summary=[begin
        rr=filter(r->r.model==id,rows)
        (model=id,requested=repeats,completed=length(rr),failed=count(r->r.model==id,failures),
         best_bic=isempty(rr) ? NaN : minimum(r.bic for r in rr),
         median_bic=isempty(rr) ? NaN : median([r.bic for r in rr]),
         worst_bic=isempty(rr) ? NaN : maximum(r.bic for r in rr),
         within_2_of_reference=count(r->abs(r.delta_reference)<=2,rr),
         better_than_reference_by_2=count(r->r.delta_reference < -2,rr),
         stopping_met=count(r->r.stopping_met,rr))
    end for id in ids]
    CSV.write(joinpath(out,"summary.csv"),summary)
    settings["status"]=isempty(failures) ? "complete" : "complete_with_failures";write_settings()
    write(joinpath(out,"COMPLETE"),settings["status"]*"\n")
    println(summary)
end
main()
