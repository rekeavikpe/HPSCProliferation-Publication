# David's October 20 requests. Original publication model/condition numbering.
include("run_david_analyses.jl")
const REQUESTED5=(:mu_max,:K_d,:K_oG,:Y_GLac,:alpha)
const EIGHT=("baseline","shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both","full")

struct DayOneContext
    config; observations; conditions; y; sigma
end
function dayone(ctx)
    qs=Dict()
    for (c,q) in ctx.conditions
        u=[only(r.mean for r in ctx.observations if r.condition==c&&r.day==1&&r.state==s) for s in STATES]
        qs[c]=merge(q,(u0=u,times=filter(>(24),q.times),feeds=filter(f->f.time_h>24,q.feeds)))
    end
    ix=findall(r->r.day>1,ctx.observations)
    DayOneContext(ctx.config,ctx.observations[ix],qs,ctx.y[ix],ctx.sigma[ix])
end
function Workflow.prediction(ctx::DayOneContext,m,p;reltol=1e-12,abstol=1e-13)
    yp=zeros(length(ctx.y))
    # exp(-(t+24)/12) = exp(-24/12)*exp(-t/12): preserve absolute culture age.
    pp=merge(p,(K_d=p.K_d*exp(-2),))
    for (c,q) in ctx.conditions
        sim=simulate(m,pp,q.u0,q.times.-24,[(time_h=f.time_h-24,fraction=f.fraction) for f in q.feeds];fresh_glucose=q.glucose,exposure=q.lookup,reltol,abstol)
        all(x->isfinite(x)&&x>=0,sim.predictions)||error("Invalid day-one prediction")
        for (i,r) in enumerate(ctx.observations)
            r.condition==c||continue
            yp[i]=sim.predictions[findfirst(==(r.state),STATES),findfirst(==(24r.day),q.times)]
        end
    end
    yp
end
route(ctx,s)=s in ("baseline","full") ? ctx : RoutedExposure.routed_context(ctx,s;shared=false)
function report_metrics(ctx,m,p;tag="",cohort="")
    yp=prediction(ctx,m,p)
    [(;tag,cohort,condition=c,state=s,n=length(ix),mrpe_percent=100mean(abs.((yp[ix].-ctx.y[ix])./ctx.y[ix])))
      for c in (0,sort(collect(keys(ctx.conditions)))...),s in (STATES...,"all")
      for ix in [findall(r->(c==0||r.condition==c)&&(s=="all"||r.state==s)&&r.mean!=0,ctx.observations)] if !isempty(ix)]
end
function score_row(ctx,m,p,tag,cohort)
    nll=gaussian_nll(ctx.y,prediction(ctx,m,p),ctx.sigma);k=length(free_names(ctx,m))
    (;tag,cohort,n=length(ctx.y),k,nll,bic=bic(nll,k,length(ctx.y)))
end
function source_contexts(fit,split)
    prep=SplitSelection.prepare(joinpath(fit,"calibration.toml"));verify_split(prep,split)
    cfg=TOML.parsefile(joinpath(fit,"calibration.toml"))
    @assert Set(cfg["calibration_cohort"]["series_ids"])==Set(prep.split.training)
    @assert length(prep.split.training)==19&&length(prep.split.selection)==12
    for id in (0,5)
        m=legacy_model(id)
        reference=read_parameters(joinpath(fit,"model$id/selected_parameters.csv"),m)
        original=read_parameters(joinpath(split,"model$id/selected_parameters.csv"),m)
        all(isapprox(getproperty(reference,n),getproperty(original,n);rtol=1e-12,atol=0) for n in keys(original))||error("Reference fit differs from the selected training parameters")
    end
    # Same external likelihood convention as previous 500 mL assessments.
    external=context(cfg=cfg,conditions=(4,))
    prep,external
end
function fit_new(ctx,m,dir)
    fit_context(ctx,m,dir;starts=parse(Int,get(ENV,"HPSC_REVIEW_STARTS","100")),budget=parse(Int,get(ENV,"HPSC_REVIEW_EVALUATIONS","20000")),iterations=parse(Int,get(ENV,"HPSC_REVIEW_ITERATIONS","3000")))
end
function run_october(stage,fit,split,out)
    prep,external=source_contexts(fit,split);mkpath(out)
    files=[@__FILE__,joinpath(fit,"calibration.toml"),joinpath(split,"split.csv")]
    append!(files,[joinpath(split,"model$i/selected_parameters.csv") for i in 0:6])
    lock_inputs(out,Dict("stage"=>stage,"starts"=>get(ENV,"HPSC_REVIEW_STARTS","100"),"budget"=>get(ENV,"HPSC_REVIEW_EVALUATIONS","20000"),"iterations"=>get(ENV,"HPSC_REVIEW_ITERATIONS","3000"),"day1"=>"all four measured day1 states, separate partitions, score days>1, absolute death clock","numbering"=>"original R1"),files)
    CSV.write(joinpath(out,"split.csv"),prep.split.assignments)
    rows=NamedTuple[];metrics=NamedTuple[]
    if stage=="p75"
        for id in 0:6
            m=legacy_model(id);train=exposure_context(prep.train,"lag_p75");held=exposure_context(prep.held,"lag_p75")
            dir=joinpath(out,"model$m");mkpath(dir)
            if m==0
                p=read_parameters(joinpath(split,"model0/selected_parameters.csv"),0)
                save_parameters(joinpath(dir,"selected_parameters.csv"),p)
            elseif m==6 && isfile(joinpath(get(ENV,"HPSC_P75_SOURCE","results/david_revision_20261014/exposure_selection"),"lag_p75/selected_parameters.csv"))
                source=get(ENV,"HPSC_P75_SOURCE","results/david_revision_20261014/exposure_selection")
                manifest=TOML.parsefile(joinpath(source,"manifest.toml"))
                manifest["settings"]["starts"]==parse(Int,get(ENV,"HPSC_REVIEW_STARTS","100"))||error("P75 source search budget differs")
                saved=collect(CSV.File(joinpath(source,"split.csv")))
                Set(String(r.series_id) for r in saved if r.role=="calibration")==Set(prep.split.training)||error("P75 source split differs")
                Set(String(r.series_id) for r in saved if r.role=="selection")==Set(prep.split.selection)||error("P75 selection split differs")
                p=read_parameters(joinpath(source,"lag_p75/selected_parameters.csv"),6)
                old=only(CSV.File(joinpath(source,"lag_p75/summary.csv")))
                isapprox(gaussian_nll(train.y,prediction(train,6,p),train.sigma),old.nll;atol=1e-6,rtol=0)||error("P75 saved likelihood differs")
                save_parameters(joinpath(dir,"selected_parameters.csv"),p)
                write(joinpath(dir,"reused_source.txt"),source)
            else
                p,_=fit_new(train,m,dir)
            end
            push!(rows,score_row(held,m,p,"Model $m","selection40"))
            append!(metrics,report_metrics(held,m,p;tag="Model $m",cohort="selection40"))
            CSV.write(joinpath(out,"selection_scores.csv"),sort(rows;by=r->r.bic));CSV.write(joinpath(out,"errors.csv"),metrics)
        end
    elseif stage=="day1"
        for s in EIGHT
            m=s=="baseline" ? 0 : 6
            train=dayone(route(prep.train,s));held=dayone(route(prep.held,s))
            p,_=fit_new(train,m,joinpath(out,s))
            append!(rows,[score_row(c,m,p,s,label) for (label,c) in (("selection40",held),("500mL_day1_initialized",dayone(route(external,s))))])
            append!(metrics,report_metrics(held,m,p;tag=s,cohort="selection40"))
            CSV.write(joinpath(out,"scores.csv"),rows);CSV.write(joinpath(out,"errors.csv"),metrics)
            CSV.write(joinpath(out,s,"initial_states.csv"),[(cohort=label,condition=c,state=STATES[i],value=q.u0[i]) for (label,ctx) in (("training60",train),("selection40",held)) for (c,q) in ctx.conditions for i in 1:4])
        end
    elseif stage=="reuse"
        p=read_parameters(joinpath(fit,"model5/selected_parameters.csv"),6)
        for (label,c) in (("training60",prep.train),("selection40",prep.held),("500mL",external))
            append!(metrics,report_metrics(c,6,p;tag="Model 6",cohort=label))
        end
        CSV.write(joinpath(out,"model6_errors.csv"),metrics)
        CSV.write(joinpath(out,"candidate_selection_optimizer_bic.csv"),[(model=legacy_model(r.storage_model),run=r.run,bic=r.bic) for r in CSV.File(joinpath(fit,"selection_optimizer_bic.csv"))])
        CSV.write(joinpath(out,"candidate_selection_scores.csv"),[score_row(prep.held,legacy_model(id),read_parameters(joinpath(split,"model$id/selected_parameters.csv"),legacy_model(id)),"Model $(legacy_model(id))","selection40") for id in 0:6])
        abl=joinpath(dirname(fit),"separate_ablations")
        TOML.parsefile(joinpath(abl,"manifest.toml"))["settings"]["shared_response_parameter"]==false||error("Expected separate-response ablations")
        checks=NamedTuple[]
        for s in EIGHT
            m=s=="baseline" ? 0 : 6;pa=read_parameters(joinpath(abl,s,"selected_parameters.csv"),m)
            ctx=route(prep.train,s);nll=gaussian_nll(ctx.y,prediction(ctx,m,pa),ctx.sigma)
            saved=only(r for r in CSV.File(joinpath(abl,"ranking.csv")) if r.scenario==s)
            isapprox(nll,saved.nll;atol=1e-6,rtol=0)||error("Ablation source is not current training fit: $s")
            for (label,c) in (("selection40",route(prep.held,s)),("500mL",route(external,s)))
                push!(rows,score_row(c,m,pa,s,label));append!(metrics,report_metrics(c,m,pa;tag=s,cohort=label))
            end
            for n in free_names(ctx,m)
                lo,hi=ctx.config["bounds"][string(n)];v=getproperty(pa,n)
                push!(checks,(scenario=s,parameter=string(n),value=v,lower=lo,upper=hi,normalized=log(v/lo)/log(hi/lo)))
            end
        end
        CSV.write(joinpath(out,"ablation_scores.csv"),rows);CSV.write(joinpath(out,"all_errors.csv"),metrics);CSV.write(joinpath(out,"ablation_parameters.csv"),checks)
        cp(joinpath(abl,"selection_optimizer_bic.csv"),joinpath(out,"selection_optimizer_bic.csv");force=true)
        # Predicted fold expansion and actual growth/aggregate response along the existing trajectories.
        fold=NamedTuple[];response=NamedTuple[]
        for (label,ctx) in (("training60",prep.train),("selection40",prep.held),("all100_descriptive",context(ids=vcat(prep.split.training,prep.split.selection),cfg=prep.cfg)))
            yp=prediction(ctx,6,p)
            for (c,q) in ctx.conditions
                ix=findall(r->r.condition==c&&r.state=="cells",ctx.observations)
                a=only(i for i in ix if ctx.observations[i].day==1)
                for endpoint in ("day5","final_day")
                    b=endpoint=="day5" ? only(i for i in ix if ctx.observations[i].day==5) : argmax(i->ctx.observations[i].day,ix)
                    push!(fold,(cohort=label,condition=c,endpoint,final_day=ctx.observations[b].day,observed_ratio=ctx.y[b]/ctx.y[a],predicted_ratio=yp[b]/yp[a]))
                end
            end
        end
        for s in EIGHT
            m=s=="baseline" ? 0 : 6;pa=read_parameters(joinpath(abl,s,"selected_parameters.csv"),m)
            for (c,q) in route(external,s).conditions
                sim=simulate(m,pa,q.u0,q.times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup)
                for (j,t) in enumerate(q.times)
                    d=sim.predictions[4,j];a,b=q.lookup(d)
                    push!(response,(scenario=s,condition=c,time_h=t,aggregate_um=d,cell_input=a,aggregate_input=b,cell_factor=m==0 ? 1. : pa.K_tau/(pa.K_tau+a),aggregate_loss_fraction=m==0 ? 0. : b*d/(pa.K_eps*pa.alpha)))
                end
            end
        end
        CSV.write(joinpath(out,"fold_predictions.csv"),fold);CSV.write(joinpath(out,"ablation_response_terms.csv"),response)
        shear=NamedTuple[]
        for (c,q) in prep.train.conditions
            times=collect(0.:1.:maximum(q.times));sim=simulate(6,p,q.u0,times,q.feeds;fresh_glucose=q.glucose,exposure=q.lookup)
            for (j,t) in enumerate(times)
                tau,_=q.lookup(sim.predictions[4,j]);push!(shear,(condition=c,time_h=t,tau_pa=tau,inhibition_percent=100tau/(p.K_tau+tau)))
            end
        end
        CSV.write(joinpath(out,"shear_inhibition.csv"),shear)
        CSV.write(joinpath(out,"shear_inhibition_summary.csv"),[(condition=c,minimum_percent=minimum(r.inhibition_percent for r in shear if c==0||r.condition==c),maximum_percent=maximum(r.inhibition_percent for r in shear if c==0||r.condition==c)) for c in (0,1,2,3,5,6)])
        # Existing CFD summaries only, not an eleven-summary refit.
        summaries=collect(CSV.File(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))
        conditions=Set([(100,30),(100,40),(100,60),(500,30)])
        table=filter(r->(r.volume_ml,r.rpm) in conditions && (startswith(r.scenario,"eul_")||r.size_um in (200,400)),summaries)
        CSV.write(joinpath(out,"exposure_distributions.csv"),table)
        CSV.write(joinpath(out,"exposure_rpm_sweep.csv"),filter(r->r.volume_ml==100&&(r.scenario=="eul_va"||r.scenario=="lag_p50"&&r.size_um==200),summaries))
    else
        error("Unknown stage $stage")
    end
    write(joinpath(out,"COMPLETE"),"Computed. Review scores and diagnostics before wording conclusions.\n")
end
