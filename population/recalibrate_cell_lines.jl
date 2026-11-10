using CSV, TOML, SHA, Random, Statistics, LinearAlgebra, Optim, Test
import BlackBoxOptim as BBO
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
const FREE=(:mu_max,:alpha,:K_G,:K_oG,:Y_GLac,:K_d,:K_Agg)
const BOUNDS=((.05,.5),(10.,200.),(1.,50.),(1e-7,1e-4),(.1,1.5),(.01,.5),(1e6,5e8))
const SCALES=(5e5,5.,5.,100.)
const WEIGHTS=(1.,1.,1.,1.5)
function main()
    3<=length(ARGS)<=5 || error("Usage: recalibrate_cell_lines.jl FIT FROZEN OUT [STARTS=20] [EVALUATIONS=20000]")
    fit,frozen,out=abspath.(ARGS[1:3]); starts=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 20;budget=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 20000
    starts>0 && budget>=100 || error("Invalid search budget")
    ispath(out) && error("Choose a new output directory")
    isfile(joinpath(frozen,"COMPLETE")) || error("Incomplete frozen predictions")
    source=joinpath(fit,"model5/selected_parameters.csv"); digest=bytes2hex(sha256(read(source)))
    lock=only(CSV.File(joinpath(fit,"selection_lock.csv")))
    lock.model==5 && lock.parameter_sha256==digest || error("Model 5 selection changed")
    hashes=collect(CSV.File(joinpath(frozen,"source_hashes.csv")))
    only(filter(r->endswith(r.path,"model5/selected_parameters.csv"),hashes)).sha256==digest || error("Frozen parameters do not match calibration")
    parent=read_parameters(source,6)
    obs=filter(r->r.model==5,collect(CSV.File(joinpath(frozen,"predictions.csv"))))
    grids=collect(CSV.File(joinpath(frozen,"exposure_grid.csv")))
    initials=collect(CSV.File(joinpath(frozen,"initial_conditions.csv")))
    setups=Dict()
    for line in ("585A1","SCTi003-A")
        rows=filter(r->r.cell_line==line,obs);days=Int(maximum(r.day for r in rows))
        q=only(filter(r->r.cell_line==line,initials)); g=sort(filter(r->r.cell_line==line,grids);by=r->r.size_um)
        lookup=linear_exposure([r.size_um for r in g],[r.shear_pa for r in g],[r.edr_m2_s3 for r in g])
        feeds=[(time_h=24.0d,fraction=line=="SCTi003-A" && d==2 ? .5 : .8) for d in 2:days-1]
        setups[line]=(;rows,days,lookup,feeds,u0=[q.initial_cells_ml,q.initial_glucose_mM,q.initial_lactate_mM,q.initial_aggregate_um],glucose=q.fresh_glucose_mM)
    end
    function decode(z)
        vals=[begin lo,hi=BOUNDS[i]; x=clamp(z[i],0.,1.); i in (4,7) ? exp(log(lo)+x*log(hi/lo)) : lo+x*(hi-lo) end for i in 1:7]
        merge(parent,NamedTuple{FREE}(Tuple(vals)))
    end
    function simulate_line(q,p;dense=false,tight=false)
        simulate(6,p,q.u0,dense ? collect(0.:.25:24q.days) : collect(24.:24.:24q.days),q.feeds;fresh_glucose=q.glucose,exposure=q.lookup,reltol=tight ? 1e-12 : 1e-11,abstol=tight ? 1e-13 : 1e-12)
    end
    function score(q,p)
        sim=simulate_line(q,p)
        all(isfinite,sim.predictions) && minimum(sim.predictions)>=0 || error("Invalid trajectory")
        sum(begin k=findfirst(==(r.state),STATES);WEIGHTS[k]*((sim.predictions[k,Int(r.day)]-r.observed)/SCALES[k])^2 end for r in q.rows)
    end
    mkpath(out);BLAS.set_num_threads(1)
    settings=Dict("source_fit"=>fit,"source_frozen"=>frozen,"parameters_reestimated"=>true,"starts"=>starts,"population_size"=>30,"max_evaluations"=>budget,"local_iterations"=>3000,"seed"=>20261005,"threads"=>Threads.nthreads(),"julia"=>string(VERSION),"free_parameters"=>collect(string.(FREE)),"fixed_parameters"=>["K_Lac","Y_GX","K_tau","K_eps"],"objective"=>"Original state-scaled weighted SSR: scales [500000,5,5,100], weights [1,1,1,1.5]; not the Gaussian calibration likelihood","coordinates"=>"Original linear coordinates except logarithmic K_oG and K_Agg; clamp to bounds","method"=>"adaptive_de_rand_1_bin_radiuslimited then NelderMead; seeded evaluation budget replaces original wall-clock budget","scope"=>"Cell-line-specific in-sample recalibration, not independent validation","initial_conditions"=>"Copied from frozen evaluation; measured day-zero glucose, fresh glucose assumed equal; no-DS cohorts")
    open(io->TOML.print(io,settings),joinpath(out,"analysis.toml"),"w")
    CSV.write(joinpath(out,"parameter_policy.csv"),[(parameter=string(n),free=n in FREE,source_value=getproperty(parent,n),lower=n in FREE ? BOUNDS[findfirst(==(n),FREE)][1] : missing,upper=n in FREE ? BOUNDS[findfirst(==(n),FREE)][2] : missing) for n in keys(parent)])
    files=[source,joinpath(fit,"model0/selected_parameters.csv"),joinpath(fit,"selection_lock.csv"),@__FILE__,joinpath(ROOT,"src/Proliferation.jl"),joinpath(@__DIR__,"Manifest.toml")]
    append!(files,[joinpath(frozen,n) for n in ("predictions.csv","exposure_grid.csv","initial_conditions.csv","source_hashes.csv","analysis.toml")])
    CSV.write(joinpath(out,"source_hashes.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in files])
    jobs=[(line=line,run=i,seed=20261005+1000j+i) for (j,line) in enumerate(("585A1","SCTi003-A")) for i in 1:starts]
    Threads.@threads :greedy for job in jobs
        dir=joinpath(out,job.line,"run$(job.run)");mkpath(dir);q=setups[job.line]
        objective=z->try score(q,decode(z)) catch e; e isa InterruptException && rethrow(); 1e100 end
        try
            Random.seed!(job.seed)
            global_result=BBO.bboptimize(objective;SearchRange=(0.,1.),NumDimensions=7,Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=30,MaxFuncEvals=budget,MaxSteps=0,MaxStepsWithoutProgress=budget+100,MaxTime=0.,TraceMode=:silent)
            zg=copy(BBO.best_candidate(global_result));pg=decode(zg);sg=score(q,pg)
            local_result=optimize(objective,zg,NelderMead(),Optim.Options(iterations=3000,g_abstol=1e-8))
            pl=decode(Optim.minimizer(local_result));sl=score(q,pl)
            use_local=sl<=sg+1e-10; chosen=use_local ? pl : pg; loss=use_local ? sl : sg
            save_parameters(joinpath(dir,"population_parameters.csv"),pg);save_parameters(joinpath(dir,"polished_parameters.csv"),pl);save_parameters(joinpath(dir,"selected_parameters.csv"),chosen)
            bounds=join([string(n)*":"*(abs(getproperty(chosen,n)-BOUNDS[i][1])<=1e-6*(BOUNDS[i][2]-BOUNDS[i][1]) ? "lower" : "upper") for (i,n) in enumerate(FREE) if min(abs(getproperty(chosen,n)-BOUNDS[i][1]),abs(getproperty(chosen,n)-BOUNDS[i][2]))<=1e-6*(BOUNDS[i][2]-BOUNDS[i][1])],";")
            CSV.write(joinpath(dir,"summary.csv"),[(cell_line=job.line,run=job.run,seed=job.seed,ssr=loss,population_ssr=sg,polished_ssr=sl,stage=use_local ? "polished" : "population",stopping_met=use_local && Optim.converged(local_result),iterations=Optim.iterations(local_result),evaluations=BBO.f_calls(global_result),boundaries=bounds)])
            write(joinpath(dir,"COMPLETE"),"Hybrid search completed\n")
        catch e
            e isa InterruptException && rethrow();write(joinpath(dir,"FAILED"),sprint(showerror,e))
        end
    end
    summaries=NamedTuple[]
    for job in jobs
        dir=joinpath(out,job.line,"run$(job.run)")
        isfile(joinpath(dir,"FAILED")) && error("Failed start retained in $dir; evaluation not finalized")
        push!(summaries,NamedTuple(only(CSV.File(joinpath(dir,"summary.csv")))))
    end
    CSV.write(joinpath(out,"runs.csv"),summaries)
    points=NamedTuple[NamedTuple(r) for r in CSV.File(joinpath(frozen,"predictions.csv"))]
    curves=NamedTuple[];ranking=NamedTuple[]
    @testset "Cell-line recalibration verification" begin
        for line in ("585A1","SCTi003-A")
            rows=filter(r->r.cell_line==line,summaries); minimum_ssr=minimum(r.ssr for r in rows)
            eligible=filter(r->r.stopping_met && r.ssr<=minimum_ssr+1e-8,rows)
            best=argmin(r->r.ssr,isempty(eligible) ? rows : eligible)
            p=read_parameters(joinpath(out,line,"run$(best.run)","selected_parameters.csv"),6);q=setups[line]
            @test score(q,p) ≈ best.ssr atol=1e-8
            for n in (:K_Lac,:Y_GX,:K_tau,:K_eps); @test getproperty(p,n)==getproperty(parent,n);end
            save_parameters(joinpath(out,line,"selected_parameters.csv"),p)
            push!(ranking,merge(best,(minimum_ssr=minimum_ssr,frozen_ssr=score(q,parent),)))
            sim=simulate_line(q,p;dense=true);tight=simulate_line(q,p;tight=true)
            for r in q.rows
                k=findfirst(==(r.state),STATES);pred=sim.predictions[k,Int(24r.day/.25)+1]
                @test abs(pred-tight.predictions[k,Int(r.day)])/max(abs(r.observed),.5)<1e-5
                push!(points,merge(NamedTuple(r),(cohort="$line recalibration (in sample, no DS)",predicted=pred)))
            end
            for (j,t) in enumerate(0.:.25:24q.days)
                for (k,state) in enumerate(STATES); push!(curves,(model=5,cell_line=line,day=t/24,state=state,side="before_or_between_exchange",value=sim.predictions[k,j]));end
                e=findfirst(e->e.time_h==t,sim.events)
                if e!==nothing
                    @test sim.events[e].before ≈ sim.predictions[:,j]
                    for (k,state) in enumerate(STATES);push!(curves,(model=5,cell_line=line,day=t/24,state=state,side="after_exchange",value=sim.events[e].after[k]));end
                end
            end
        end
    end
    CSV.write(joinpath(out,"summary.csv"),ranking)
    CSV.write(joinpath(out,"predictions.csv"),points);CSV.write(joinpath(out,"trajectories.csv"),curves)
    write(joinpath(out,"COMPLETE"),"Selected fits verified; consult stopping flags and bounds in summary.csv\n")
end
main()
