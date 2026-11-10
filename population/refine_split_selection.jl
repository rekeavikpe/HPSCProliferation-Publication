using CSV, TOML, SHA, Random, LinearAlgebra
include("SplitSelection.jl")
using .SplitSelection
using .SplitSelection.Workflow, .SplitSelection.Workflow.Calibration, .SplitSelection.Workflow.Calibration.Proliferation
using .SplitSelection.PublicationModels

function main()
    length(ARGS)==2 || error("Usage: refine_split_selection.jl ORIGINAL_SPLIT OUTPUT")
    source,out=abspath.(ARGS)
    (out==source || startswith(out,source*"/")) && error("Use a separate output directory")
    isfile(joinpath(source,"COMPLETE")) || error("Original split is incomplete")
    original=TOML.parsefile(joinpath(source,"protocol.toml"))
    sort(Int.(original["storage_models"]))==collect(0:6) || error("Expected seven models")
    prep=prepare(joinpath(source,"calibration.toml");seed=original["split_seed"])
    ctx=prep.train;BLAS.set_num_threads(1)
    assignments=collect(CSV.File(joinpath(source,"split.csv")))
    [(r.series_id,r.role) for r in assignments]==[(r.series_id,r.role) for r in prep.split.assignments] || error("Split assignment changed")
    for (cohort,c) in (("calibration",ctx),("selection",prep.held))
        saved=collect(CSV.File(joinpath(source,cohort*"_observations.csv")))
        length(saved)==length(c.observations) || error("Observation count changed")
        for (a,b,sigma) in zip(saved,c.observations,c.sigma)
            (a.condition,a.day,a.state,a.n_vessels)==(b.condition,b.day,b.state,b.n_vessels) || error("Observation identity changed")
            all(isapprox(x,y;atol=1e-12,rtol=1e-12) for (x,y) in ((a.mean,b.mean),(a.sd,b.sd),(a.sigma,sigma))) || error("Observation or weight changed")
        end
    end
    iterations=parse(Int,get(ENV,"HPSC_REFINE_ITERATIONS","12000"));iterations>0 || error("Invalid iteration limit")
    targets=[6,1,4,5]
    jobs=NamedTuple[]
    for pub in targets
        id=storage_model(pub);m=legacy_model(id);dir=joinpath(source,"model$id")
        bestp=read_parameters(joinpath(dir,"selected_parameters.csv"),m);zbest=encode(ctx,m,bestp)
        candidates=NamedTuple[]
        for r in 1:original["starts"]
            path=joinpath(dir,"run$r");checkpoint(ctx,m,path) || error("Incomplete source checkpoint")
            for stage in CSV.File(joinpath(path,"stages.csv"))
                push!(candidates,(nll=stage.nll,path=joinpath(path,String(stage.stage)*"_parameters.csv")))
            end
        end
        sort!(candidates;by=r->r.nll)
        # Two additional saved population starts; choose them by calibration score only.
        anchor=only(CSV.File(joinpath(dir,"summary.csv")))
        distinct=filter(r->basename(dirname(r.path))!=anchor.run,candidates)
        alternatives=NamedTuple[];seen=Set{String}()
        for row in distinct
            run=basename(dirname(row.path));run in seen && continue
            push!(seen,run);push!(alternatives,row);length(alternatives)==2 && break
        end
        for j in 1:6
            if j==1
                z=copy(zbest);kind="saved selected fit"
            elseif j<=3 && length(alternatives)>=j-1
                q=read_parameters(alternatives[j-1].path,m);z=encode(ctx,m,q);kind=relpath(alternatives[j-1].path,source)
            else
                rng=MersenneTwister(20261006+100id+j)
                z=clamp.(zbest.+.01randn(rng,length(zbest)),0.,1.);kind="seeded 1% normalized-log perturbation"
            end
            push!(jobs,(id=id,pub=pub,run=original["starts"]+j,z=z,origin=kind))
        end
    end
    files=vcat([@__FILE__,joinpath(@__DIR__,"SplitSelection.jl"),joinpath(@__DIR__,"HybridRefinement.jl"),joinpath(@__DIR__,"Manifest.toml")],
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),
        filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    for (dir,_,names) in walkdir(source)
        for name in names
            (endswith(name,".csv") || endswith(name,".toml")) && push!(files,joinpath(dir,name))
        end
    end
    buf=IOBuffer();CSV.write(buf,[(path=relpath(p,ROOT),sha256=bytes2hex(sha256(read(p)))) for p in sort(unique(files))]);fingerprint=String(take!(buf))
    settings=Dict("source"=>source,"publication_models"=>targets,"starts_per_target"=>6,"local_iterations"=>iterations,
        "optimizer"=>"Original NelderMead with clipped normalized-log coordinates; unchanged g_abstol=1e-8",
        "seed"=>20261006,"selection_policy"=>"Unchanged; held-out outcomes are not used for optimization or choosing training candidates")
    if isdir(out)
        TOML.parsefile(joinpath(out,"refinement_protocol.toml"))==settings || error("Resume settings changed")
        read(joinpath(out,"refinement_source_hashes.csv"),String)==fingerprint || error("Source changed; use a new output directory")
    else
        mkpath(out)
        open(io->TOML.print(io,settings),joinpath(out,"refinement_protocol.toml"),"w")
        write(joinpath(out,"refinement_source_hashes.csv"),fingerprint)
    end
    # Copy immutable original candidates; restart safely if copying was interrupted.
    if !isfile(joinpath(out,"BASE_COPIED"))
        for id in 0:6
            cp(joinpath(source,"model$id"),joinpath(out,"model$id");force=true)
        end
        for name in ("protocol.toml","calibration.toml","split.csv","split_counts.csv","calibration_observations.csv","selection_observations.csv")
            cp(joinpath(source,name),joinpath(out,name);force=true)
        end
        write(joinpath(out,"BASE_COPIED"),"Original candidates preserved\n")
    end
    Threads.@threads :greedy for job in jobs
        dir=joinpath(out,"model$(job.id)","run$(job.run)");m=legacy_model(job.id)
        checkpoint(ctx,m,dir) && continue
        mkpath(dir)
        try
            p=decode(ctx,m,job.z);nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
            save_parameters(joinpath(dir,"initial_parameters.csv"),p)
            stages=HybridRefinement.refine_original(ctx,m,job.z,joinpath(dir,"local");iterations)
            cp(joinpath(dir,"local","polished_parameters.csv"),joinpath(dir,"polished_parameters.csv");force=true)
            CSV.write(joinpath(dir,"stages.csv"),vcat([(stage="initial",nll=nll,stopping_met=false)],stages))
            write(joinpath(dir,"origin.txt"),job.origin*"\n")
            for row in CSV.File(joinpath(dir,"stages.csv"))
                p=read_parameters(joinpath(dir,String(row.stage)*"_parameters.csv"),m)
                isapprox(gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma),row.nll;atol=1e-6,rtol=0) || error("Saved objective mismatch")
            end
            rm(joinpath(dir,"FAILED");force=true)
            write(joinpath(dir,"COMPLETE.tmp"),"Verified targeted calibration-only refinement\n")
            mv(joinpath(dir,"COMPLETE.tmp"),joinpath(dir,"COMPLETE");force=true)
            @info "Targeted refinement complete" model=job.pub run=job.run
        catch e
            e isa InterruptException && rethrow()
            write(joinpath(dir,"FAILED"),sprint(showerror,e))
            @warn "Refinement retained for retry" model=job.pub exception=sprint(showerror,e)
        end
    end
    all(j->isfile(joinpath(out,"model$(j.id)","run$(j.run)","COMPLETE")),jobs) || error("Some refinements incomplete; rerun the same command")
    chosen=evaluate(prep,out,collect(0:6))
    old=collect(CSV.File(joinpath(source,"selection_ranking.csv")));new=collect(CSV.File(joinpath(out,"selection_ranking.csv")))
    comparison=map(0:6) do pub
        a=only(filter(r->r.publication_model==pub,old));b=only(filter(r->r.publication_model==pub,new))
        (publication_model=pub,old_training_nll=a.training_nll,new_training_nll=b.training_nll,
         old_training_bic=a.training_bic,new_training_bic=b.training_bic,
         old_stopping_met=a.training_stopping_met,new_stopping_met=b.training_stopping_met,
         old_selection_nll=a.selection_nll,new_selection_nll=b.selection_nll)
    end
    CSV.write(joinpath(out,"refinement_comparison.csv"),comparison)
    write(joinpath(out,"COMPLETE"),"Targeted refinement and held-out rescoring complete; original output untouched\n")
    println("Selected publication Model $(chosen.publication_model). Inspect refinement_comparison.csv; convergence is not guaranteed.")
end
main()
