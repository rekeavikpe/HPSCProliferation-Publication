using CSV, TOML, SHA, Random, LinearAlgebra
import BlackBoxOptim as BBO
include("SplitSelection.jl")
using .SplitSelection
using .SplitSelection.Workflow, .SplitSelection.Workflow.Calibration, .SplitSelection.Workflow.Calibration.Proliferation
using .SplitSelection.PublicationModels

function main()
    length(ARGS)==2 || error("Usage: run_split_selection.jl CONFIG OUTPUT")
    cfgfile,out=abspath.(ARGS)
    starts=parse(Int,get(ENV,"HPSC_SPLIT_STARTS","100"))
    budget=parse(Int,get(ENV,"HPSC_SPLIT_EVALUATIONS","20000"))
    iterations=parse(Int,get(ENV,"HPSC_SPLIT_LOCAL_ITERATIONS","3000"))
    ids=parse.(Int,split(get(ENV,"HPSC_SPLIT_STORAGE_IDS","0,1,2,3,4,5,6"),','))
    allunique(ids)&&all(in(0:6),ids) || error("Invalid storage IDs")
    starts>0&&budget>=100&&iterations>0 || error("Invalid budgets")
    prep=prepare(cfgfile);ctx=prep.train;BLAS.set_num_threads(1)
    protocol=Dict("split_seed"=>20261003,"optimizer_seed"=>20261005,"starts"=>starts,
        "max_evaluations"=>budget,"local_iterations"=>iterations,"storage_models"=>ids,
        "primary_selection"=>"held-out negative log-likelihood; no selection refitting",
        "tie_rule"=>"within 1e-6 NLL prefer fewer parameters, then lower training BIC",
        "weights"=>"calibration condition-day-state SD with configured floors; reused unchanged for selection",
        "initial_conditions"=>"partition-specific measured day-zero glucose; calibration fallback only if unavailable; fresh glucose assumed equal to initial glucose",
        "scope"=>"retrospective internal selection; culture-level random split; final refit and external evaluations are separate",
        "julia"=>string(VERSION),"BlackBoxOptim"=>string(pkgversion(BBO)))
    files=vcat([cfgfile,@__FILE__,joinpath(@__DIR__,"SplitSelection.jl"),joinpath(@__DIR__,"HybridRefinement.jl"),joinpath(@__DIR__,"Manifest.toml")],
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),
        filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    hashes=[(path=relpath(p,ROOT),sha256=bytes2hex(sha256(read(p)))) for p in sort(unique(files))]
    buffer=IOBuffer();CSV.write(buffer,hashes);fingerprint=String(take!(buffer))
    if isdir(out)
        TOML.parsefile(joinpath(out,"protocol.toml"))==protocol || error("Resume protocol differs; use a new output directory")
        read(joinpath(out,"source_hashes.csv"),String)==fingerprint || error("Source changed; preserve results and use a new output directory")
        saved=collect(CSV.File(joinpath(out,"split.csv")))
        [(r.series_id,r.role) for r in saved]==[(r.series_id,r.role) for r in prep.split.assignments] || error("Split manifest changed")
    else
        mkpath(out)
        open(io->TOML.print(io,protocol),joinpath(out,"protocol.toml"),"w")
        write(joinpath(out,"source_hashes.csv"),fingerprint)
        cp(cfgfile,joinpath(out,"calibration.toml"))
        CSV.write(joinpath(out,"split.csv"),prep.split.assignments)
    end
    for (name,q) in (("calibration",ctx),("selection",prep.held))
        CSV.write(joinpath(out,name*"_observations.csv"),[merge(r,(sigma=s,)) for (r,s) in zip(q.observations,q.sigma)])
    end
    counts=[(condition=c,calibration_cultures=count(r->r.legacy_condition==c&&r.role=="calibration",prep.split.assignments),selection_cultures=count(r->r.legacy_condition==c&&r.role=="selection",prep.split.assignments)) for c in (1,2,3,5,6)]
    CSV.write(joinpath(out,"split_counts.csv"),counts)
    get(ENV,"HPSC_SPLIT_PREPARE_ONLY","0")=="1" && return println("Prepared locked split: $(length(prep.split.training)) calibration / $(length(prep.split.selection)) selection cultures")
    jobs=[(id=id,run=r) for id in ids for r in 1:starts]
    Threads.@threads :greedy for job in jobs
        dir=joinpath(out,"model$(job.id)","run$(job.run)");m=legacy_model(job.id)
        checkpoint(ctx,m,dir) && continue
        mkpath(dir)
        try
            Random.seed!(20261005+1000job.id+job.run)
            residual=Workflow.residual_function(ctx,m)
            objective=z->sum(abs2,residual(z))/2
            pop=BBO.bboptimize(objective;SearchRange=(0.,1.),NumDimensions=length(free_names(ctx,m)),Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=50,MaxFuncEvals=budget,MaxSteps=0,MaxStepsWithoutProgress=budget+100,MaxTime=0.,TraceMode=:silent)
            z=BBO.best_candidate(pop);p=decode(ctx,m,z)
            nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
            save_parameters(joinpath(dir,"population_parameters.csv"),p)
            stages=HybridRefinement.refine_original(ctx,m,z,joinpath(dir,"local");iterations)
            cp(joinpath(dir,"local","polished_parameters.csv"),joinpath(dir,"polished_parameters.csv");force=true)
            CSV.write(joinpath(dir,"stages.csv"),vcat([(stage="population",nll,stopping_met=false)],stages))
            open(joinpath(dir,"population.toml"),"w") do io
                TOML.print(io,Dict("evaluations"=>BBO.f_calls(pop),"stop_reason"=>string(BBO.stop_reason(pop))))
            end
            # Verify both saved candidates before marking this run reusable.
            for r in CSV.File(joinpath(dir,"stages.csv"))
                q=read_parameters(joinpath(dir,String(r.stage)*"_parameters.csv"),m)
                isapprox(gaussian_nll(ctx.y,prediction(ctx,m,q),ctx.sigma),r.nll;atol=1e-6,rtol=0) || error("Saved objective mismatch")
            end
            rm(joinpath(dir,"FAILED");force=true)
            write(joinpath(dir,"COMPLETE.tmp"),"Both stages verified; calibration partition only\n")
            mv(joinpath(dir,"COMPLETE.tmp"),joinpath(dir,"COMPLETE");force=true)
            @info "Split calibration checkpoint" model=publication_model(job.id) run=job.run
        catch e
            e isa InterruptException && rethrow()
            write(joinpath(dir,"FAILED"),sprint(showerror,e))
            @warn "Run retained for retry" model=job.id run=job.run exception=sprint(showerror,e)
        end
    end
    all(job->isfile(joinpath(out,"model$(job.id)","run$(job.run)","COMPLETE")),jobs) || error("Some runs incomplete; repeat the same command to resume")
    chosen=evaluate(prep,out,ids)
    write(joinpath(out,"COMPLETE"),"Split comparison complete; no all-data refit or external evaluation performed\n")
    println("Held-out selection: publication Model $(chosen.publication_model). Inspect convergence and selection_ranking.csv before final refit.")
end
main()
