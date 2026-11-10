using CSV, TOML, SHA, Random, Statistics, LinearAlgebra, Test
import BlackBoxOptim as BBO
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
include("HybridRefinement.jl")
using .HybridRefinement

include("../ablations/RoutedExposure.jl")
using .RoutedExposure
const VARIANTS=("shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both")
const SHARED=get(ENV,"HPSC_ABLATION_SHARED","true")=="true"
ablated(ctx,s)=routed_context(ctx,s;shared=SHARED)

function main()
    length(ARGS)==2 || error("Usage: run_routed_ablations.jl ALL100_FIT OUTPUT")
    fit,out=abspath.(ARGS);cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
    isempty(cfg["fixed"]) || error("Expected all-free source fit")
    cfg["bounds"]["K_G"][2]==100 && cfg["bounds"]["K_Lac"][2]==100 || error("Expected final bounds")
    starts=parse(Int,get(ENV,"HPSC_ABLATION_STARTS","100"))
    budget=parse(Int,get(ENV,"HPSC_ABLATION_EVALUATIONS","20000"))
    iterations=parse(Int,get(ENV,"HPSC_ABLATION_ITERATIONS","3000"))
    starts>0 && budget>=100 && iterations>0 || error("Invalid budget")
    isfile(joinpath(fit,"figure_data/optimizer_bic.csv"))||error("Export source figure_data first")
    base=context(cfg=cfg);external=context(cfg=cfg,conditions=(4,))
    sourceobs=collect(CSV.File(joinpath(fit,"observations.csv")))
    @assert [(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in sourceobs]==[(r.condition,r.day,r.state,r.mean,r.sd,r.n_vessels) for r in base.observations]
    contexts=Dict(s=>ablated(base,s) for s in VARIANTS)
    pfull=read_parameters(joinpath(fit,"model5/selected_parameters.csv"),6)
    pzero=read_parameters(joinpath(fit,"model0/selected_parameters.csv"),0)
    @testset "Ablation limits and active counts" begin
        p=merge(pzero,(K_tau=.001,K_eps=.001))
        for u in ([2e5,15.,3.,150.],[8e5,4.,10.,350.])
            a=zeros(4);b=zeros(4)
            model_rhs!(a,u,pzero,36.,0);model_rhs!(b,u,p,36.,6,0.,0.)
            @test a≈b
        end
        for ctx in values(contexts)
            @test length(free_names(ctx,6)) in (10,11)
            @test ctx.y==base.y && ctx.sigma==base.sigma
        end
        for (c,q) in base.conditions,d in (15.,200.,350.,900.)
            tau,edr=q.lookup(d)
            @test contexts["shear_cell"].conditions[c].lookup(d)==(tau,0.)
            @test contexts["edr_aggregate"].conditions[c].lookup(d)==(0.,edr)
        end
    end
    settings=Dict("starts"=>starts,"evaluations"=>budget,"iterations"=>iterations,"source_fit"=>fit,"seed"=>20261006,
        "scope"=>haskey(cfg,"calibration_cohort") ? "60% training cultures only. 500 mL predictions without refitting." : "All-100-mL mechanistic sensitivity. 500-mL predictions without refitting.",
        "inactive_parameters"=>"Inactive exposure-channel placeholders are excluded; same-exposure both-channel variants tie K_eps to K_tau when shared_response_parameter=true.",
        "shared_response_parameter"=>SHARED,
        "routing"=>Dict(s=>collect(ROUTES[s]) for s in VARIANTS))
    files=vcat([@__FILE__,joinpath(@__DIR__,"HybridRefinement.jl"),joinpath(ROOT,"ablations/RoutedExposure.jl"),joinpath(@__DIR__,"Manifest.toml"),joinpath(fit,"calibration.toml"),joinpath(fit,"ranking.csv"),joinpath(fit,"observations.csv"),joinpath(fit,"figure_data/optimizer_bic.csv"),joinpath(fit,"figure_data/source.toml")],
        [joinpath(fit,"model$i/selected_parameters.csv") for i in (0,5)],
        filter(f->endswith(f,".jl"),readdir(joinpath(ROOT,"src");join=true)),
        filter(f->endswith(f,".csv"),readdir(joinpath(ROOT,"data/processed");join=true)))
    hashes=Dict(relpath(f,ROOT)=>bytes2hex(sha256(read(f))) for f in files)
    manifest=Dict("settings"=>settings,"hashes"=>hashes)
    if isdir(out)
        TOML.parsefile(joinpath(out,"manifest.toml"))==manifest || error("Inputs or budgets changed; choose a new output directory")
    else
        mkpath(out);open(io->TOML.print(io,manifest),joinpath(out,"manifest.toml"),"w")
    end
    BLAS.set_num_threads(1)
    jobs=[(scenario=s,run=i) for s in VARIANTS for i in 1:starts]
    Threads.@threads :greedy for job in jobs
        ctx=contexts[job.scenario];dir=joinpath(out,job.scenario,"run$(job.run)")
        if isfile(joinpath(dir,"COMPLETE"))
            p=read_parameters(joinpath(dir,"selected_parameters.csv"),6)
            saved=only(CSV.File(joinpath(dir,"summary.csv")))
            @assert isapprox(gaussian_nll(ctx.y,prediction(ctx,6,p),ctx.sigma),saved.nll;atol=1e-6,rtol=0)
            continue
        end
        if isdir(dir)
            archive=joinpath(out,"interrupted",job.scenario);mkpath(archive)
            mv(dir,joinpath(archive,"run$(job.run)_$(time_ns())"))
        end
        mkpath(dir);Random.seed!(20261006+job.run+1000findfirst(==(job.scenario),VARIANTS))
        residual=Workflow.residual_function(ctx,6)
        search=BBO.bboptimize(z->sum(abs2,residual(z))/2;SearchRange=(0.,1.),NumDimensions=length(free_names(ctx,6)),
            Method=:adaptive_de_rand_1_bin_radiuslimited,PopulationSize=50,MaxFuncEvals=budget,
            MaxSteps=0,MaxStepsWithoutProgress=budget+100,MaxTime=0.,TraceMode=:silent)
        z=copy(BBO.best_candidate(search));p=decode(ctx,6,z)
        nll=gaussian_nll(ctx.y,prediction(ctx,6,p),ctx.sigma)
        save_parameters(joinpath(dir,"population_parameters.csv"),p)
        stages=vcat([(stage="population",nll=nll,stopping_met=false)],refine_original(ctx,6,z,joinpath(dir,"local");iterations))
        bestscore=minimum(x.nll for x in stages);eligible=filter(x->x.stopping_met && x.nll<=bestscore+1e-6,stages)
        best=argmin(x->x.nll,isempty(eligible) ? stages : eligible)
        selected=best.stage=="population" ? p : read_parameters(joinpath(dir,"local/polished_parameters.csv"),6)
        save_parameters(joinpath(dir,"selected_parameters.csv"),selected)
        CSV.write(joinpath(dir,"stages.csv"),stages)
        bounds=join(["$n:"*(x<1e-4 ? "lower" : "upper") for (n,x) in zip(free_names(ctx,6),encode(ctx,6,selected)) if x<1e-4 || x>1-1e-4],";")
        CSV.write(joinpath(dir,"summary.csv"),[(scenario=job.scenario,run=job.run,nll=best.nll,bic=bic(best.nll,length(free_names(ctx,6)),length(ctx.y)),stopping_met=best.stopping_met,boundaries=bounds)])
        write(joinpath(dir,"COMPLETE.tmp"),"Verified saved stages\n");mv(joinpath(dir,"COMPLETE.tmp"),joinpath(dir,"COMPLETE"))
        @info "Ablation complete" scenario=job.scenario run=job.run total=starts
    end
    runs=[NamedTuple(only(CSV.File(joinpath(out,j.scenario,"run$(j.run)/summary.csv")))) for j in jobs]
    CSV.write(joinpath(out,"runs.csv"),runs)
    ranking=NamedTuple[];points=NamedTuple[];errors=NamedTuple[]
    for scenario in ("baseline",VARIANTS...,"full")
        model=scenario=="baseline" ? 0 : 6
        if haskey(contexts,scenario)
            ctx=contexts[scenario];group=filter(x->x.scenario==scenario,runs);v=minimum(x.nll for x in group)
            eligible=filter(x->x.stopping_met && x.nll<=v+1e-6,group)
            best=argmin(x->x.nll,isempty(eligible) ? group : eligible)
            p=read_parameters(joinpath(out,scenario,"run$(best.run)/selected_parameters.csv"),6)
            status=best.stopping_met;origin="new ablation";boundaries=best.boundaries
        else
            ctx=base;p=scenario=="baseline" ? pzero : pfull
            row=only(filter(x->x.model==(scenario=="baseline" ? 0 : 5),collect(CSV.File(joinpath(fit,"ranking.csv")))))
            status=row.stopping_met;origin=haskey(cfg,"calibration_cohort") ? "reused 60% training fit" : "reused all-100-mL fit";boundaries=string(row.boundaries)
        end
        k=length(free_names(ctx,model));pred=prediction(ctx,model,p);nll=gaussian_nll(ctx.y,pred,ctx.sigma)
        tight=prediction(ctx,model,p;reltol=1e-13,abstol=1e-14)
        @assert maximum(abs.(tight.-pred)./ctx.sigma)<1e-4
        push!(ranking,(scenario=scenario,n=length(ctx.y),k=k,nll=nll,bic=bic(nll,k,length(ctx.y)),stopping_met=status,boundaries=boundaries,source=origin))
        mkpath(joinpath(out,scenario));save_parameters(joinpath(out,scenario,"selected_parameters.csv"),p)
        adapter=joinpath(out,"profile_inputs",scenario);mkpath(joinpath(adapter,"model$(scenario=="baseline" ? 0 : 5)"))
        open(io->TOML.print(io,cfg),joinpath(adapter,"calibration.toml"),"w")
        open(io->TOML.print(io,Dict("scenario"=>scenario,"cell_input"=>ROUTES[scenario][1],"aggregate_input"=>ROUTES[scenario][2],"shared_response_parameter"=>SHARED,"labels"=>response_labels(scenario;shared=SHARED))),joinpath(adapter,"routing.toml"),"w")
        CSV.write(joinpath(adapter,"observations.csv"),base.observations)
        modeldir=joinpath(adapter,"model$(scenario=="baseline" ? 0 : 5)")
        save_parameters(joinpath(modeldir,"selected_parameters.csv"),p)
        CSV.write(joinpath(modeldir,"summary.csv"),[last(ranking)])
        for (cohort,c) in (("100mL calibration",ctx),("500mL retrospective prediction",haskey(contexts,scenario) ? ablated(external,scenario) : external))
            yhat=prediction(c,model,p)
            for (i,r) in enumerate(c.observations)
                push!(points,(scenario=scenario,cohort=cohort,condition=r.condition,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,predicted=yhat[i]))
            end
            for state in (STATES...,"all")
                ix=findall(i->(state=="all" || c.observations[i].state==state) && c.y[i]!=0,eachindex(c.y))
                push!(errors,(scenario=scenario,cohort=cohort,state=state,mrpe_percent=100mean(abs.((yhat[ix].-c.y[ix])./c.y[ix])),experimental_noise_percent=100mean([c.observations[i].sd/abs(c.y[i]) for i in ix if c.observations[i].n_vessels>=2])))
            end
        end
    end
    sort!(ranking;by=x->x.bic);CSV.write(joinpath(out,"ranking.csv"),ranking)
    CSV.write(joinpath(out,"predictions.csv"),points);CSV.write(joinpath(out,"state_errors.csv"),errors)
    # Package all optimizer distributions with the comparison; never fabricate reused boxes.
    source_bic=joinpath(fit,"figure_data/optimizer_bic.csv")
    isfile(source_bic)||error("Export source figure_data before running this workflow")
    metadata=TOML.parsefile(joinpath(fit,"figure_data/source.toml"))
    allruns=[(scenario=r.scenario,run=r.run,bic=r.bic,stopping_met=r.stopping_met) for r in runs]
    for r in CSV.File(source_bic)
        r.model in (0,5)||continue
        num=parse(Int,replace(r.run,"run"=>""));num<=metadata["initial_starts_per_model"]||continue
        push!(allruns,(scenario=r.model==0 ? "baseline" : "full",run=num,bic=Float64(r.bic),stopping_met=Bool(r.stopping_met)))
    end
    CSV.write(joinpath(out,"optimizer_bic.csv"),allruns)
    write(joinpath(out,"COMPLETE"),"Eight-version routing sensitivity completed with the saved source cohort\n")
end
main()
