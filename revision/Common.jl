module RevisionCommon
using CSV, TOML, SHA, Random, Statistics, LinearAlgebra
import BlackBoxOptim as BBO
include("../population/PooledValidation.jl")
using .PooledValidation
using .PooledValidation.SplitSelection
using .SplitSelection.Workflow, .SplitSelection.Workflow.Calibration, .SplitSelection.Workflow.Calibration.Proliferation
include("../ablations/RoutedExposure.jl")
using .RoutedExposure
export Workflow, PooledValidation, SplitSelection, RoutedExposure, fit_context, lock_inputs, errors, export_curves
function lock_inputs(out,settings,extra=String[])
    files=vcat([@__FILE__,joinpath(@__DIR__,"run_david_analyses.jl"),joinpath(@__DIR__,"CellContexts.jl"),joinpath(ROOT,"population/HybridRefinement.jl"),joinpath(ROOT,"population/SplitSelection.jl"),joinpath(ROOT,"population/PooledValidation.jl"),joinpath(ROOT,"population/Manifest.toml")],extra,
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    hashes=Dict(relpath(p,ROOT)=>bytes2hex(sha256(read(p))) for p in unique(files))
    manifest=Dict("settings"=>settings,"hashes"=>hashes)
    if isfile(joinpath(out,"manifest.toml"))
        TOML.parsefile(joinpath(out,"manifest.toml"))==manifest||error("Changed source/settings in $out. Preserve results and use a new output root.")
    else
        mkpath(out);open(io->TOML.print(io,manifest),joinpath(out,"manifest.toml"),"w")
    end
end
function fit_context(ctx,m,out;starts=100,budget=20000,iterations=3000,anchor=nothing)
    BLAS.set_num_threads(1);mkpath(out)
    Threads.@threads :greedy for run in 1:starts
        dir=joinpath(out,"run$run")
        SplitSelection.checkpoint(ctx,m,dir)&&continue
        mkpath(dir);Random.seed!(20261014+1000m+run)
        residual=Workflow.residual_function(ctx,m)
        pop=BBO.bboptimize(z->sum(abs2,residual(z))/2;SearchRange=(0.,1.),NumDimensions=length(free_names(ctx,m)),Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=50,MaxFuncEvals=budget,MaxSteps=0,MaxStepsWithoutProgress=budget+100,MaxTime=0.,TraceMode=:silent)
        z=copy(BBO.best_candidate(pop));p=decode(ctx,m,z)
        nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma)
        save_parameters(joinpath(dir,"population_parameters.csv"),p)
        stages=SplitSelection.HybridRefinement.refine_original(ctx,m,z,joinpath(dir,"local");iterations)
        cp(joinpath(dir,"local/polished_parameters.csv"),joinpath(dir,"polished_parameters.csv");force=true)
        allstages=vcat([(stage="population",nll,stopping_met=false)],stages)
        # Refining a source estimate is supplementary to, not a replacement for, independent searches.
        if anchor!==nothing && run==1
            z0=encode(ctx,m,anchor)
            extra=SplitSelection.HybridRefinement.refine_original(ctx,m,z0,joinpath(dir,"anchor");iterations)
            cp(joinpath(dir,"anchor/polished_parameters.csv"),joinpath(dir,"anchor_parameters.csv");force=true)
            push!(allstages,(stage="anchor",nll=only(extra).nll,stopping_met=only(extra).stopping_met))
        end
        CSV.write(joinpath(dir,"stages.csv"),allstages)
        write(joinpath(dir,"COMPLETE"),"saved stages\n")
        SplitSelection.checkpoint(ctx,m,dir)||error("Checkpoint verification failed")
        println("Completed $out start $run/$starts");flush(stdout)
    end
    # assemble counts all model parameters by default. Correct the descriptive BIC for fixed-parameter sensitivity fits.
    summary=Workflow.assemble(ctx,m,out)
    k=length(free_names(ctx,m));summary=merge(summary,(k=k,bic=bic(summary.nll,k,length(ctx.y))))
    CSV.write(joinpath(out,"summary.csv"),[summary])
    read_parameters(joinpath(out,"selected_parameters.csv"),m),summary
end
function errors(ctx,yp)
    [begin
        ix=findall(r->(s=="all"||r.state==s)&&r.day>0&&r.mean!=0,ctx.observations)
        noise=filter(i->ctx.observations[i].n_vessels>=2,ix)
        (state=s,n=length(ix),mrpe_percent=isempty(ix) ? NaN : 100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix])),experimental_noise_percent=isempty(noise) ? NaN : 100mean(ctx.observations[i].sd/abs(ctx.y[i]) for i in noise))
    end for s in (STATES...,"all")]
end
function export_curves(ctx,m,p)
    rows=NamedTuple[]
    for (c,q) in ctx.conditions
        times=sort(unique(vcat(collect(0.:.25:maximum(q.times)),q.times)))
        sim=simulate(m,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup)
        for (j,t) in enumerate(times)
            for (k,s) in enumerate(STATES);push!(rows,(condition=c,day=t/24,state=s,side="pre",value=sim.predictions[k,j]));end
            for e in sim.events
                e.time_h==t||continue
                for (k,s) in enumerate(STATES);push!(rows,(condition=c,day=t/24,state=s,side="post",value=e.after[k]));end
            end
        end
    end
    rows
end
end
