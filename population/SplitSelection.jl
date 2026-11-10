module SplitSelection
using CSV, TOML, SHA, Statistics
include("../src/Workflow.jl")
include("../src/ReplicateSplits.jl")
include("../src/PublicationModels.jl")
using .Workflow, .Workflow.Calibration, .Workflow.Calibration.Proliferation
using .ReplicateSplits, .PublicationModels
include("HybridRefinement.jl")
using .HybridRefinement
export prepare, evaluate, checkpoint, Workflow, HybridRefinement

function prepare(cfgfile; seed=20261003)
    cfg=TOML.parsefile(cfgfile)
    isempty(cfg["fixed"]) || error("Split calibration requires all parameters free; do not reuse all-data estimates")
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    split=random_split(series;seed)
    isempty(intersect(split.training,split.selection)) || error("Culture leakage")
    train=context(;ids=split.training,cfg)
    held=context(;ids=split.selection,cfg,fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
    key(r)=(r.condition,r.day,r.state)
    scales=Dict(key(r)=>s for (r,s) in zip(train.observations,train.sigma))
    all(r->haskey(scales,key(r)),held.observations) || error("Selection observation has no calibration noise weight; inspect split coverage")
    held=merge(held,(sigma=[scales[key(r)] for r in held.observations],))
    # Held-out SDs are descriptive only. Scoring uses calibration SDs and floors.
    (;cfg,split,train,held)
end

function checkpoint(ctx,model,dir)
    isfile(joinpath(dir,"COMPLETE")) || return false
    for r in CSV.File(joinpath(dir,"stages.csv"))
        p=read_parameters(joinpath(dir,String(r.stage)*"_parameters.csv"),model)
        actual=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
        isapprox(actual,r.nll;atol=1e-6,rtol=0) || error("Corrupt checkpoint: $dir")
    end
    true
end

function evaluate(prep,out,ids)
    train,held=prep.train,prep.held
    training=NamedTuple[];scores=NamedTuple[];errors=NamedTuple[];points=NamedTuple[]
    for id in ids
        m=legacy_model(id);dir=joinpath(out,"model$id")
        fit=Workflow.assemble(train,m,dir);push!(training,fit)
        p=read_parameters(joinpath(dir,"selected_parameters.csv"),m)
        pred=prediction(held,m,p)
        nll=gaussian_nll(held.y,pred,held.sigma)
        push!(scores,(publication_model=publication_model(id),storage_model=id,legacy_model=m,
            training_n=length(train.y),training_nll=fit.nll,training_bic=fit.bic,k=fit.k,
            selection_n=length(held.y),selection_nll=nll,selection_mean_nll=nll/length(held.y),
            training_stopping_met=fit.stopping_met,training_run=fit.run,training_stage=fit.stage))
        for (cohort,ctx,yp) in (("calibration",train,prediction(train,m,p)),("selection",held,pred))
            for (i,r) in enumerate(ctx.observations)
                push!(points,(publication_model=publication_model(id),storage_model=id,cohort,condition=r.condition,day=r.day,state=r.state,observed=r.mean,sd=r.sd,n_vessels=r.n_vessels,sigma=ctx.sigma[i],predicted=yp[i]))
            end
            for state in (STATES...,"all")
                ix=findall(r->state=="all" || r.state==state,ctx.observations);nz=filter(i->ctx.y[i]!=0,ix)
                push!(errors,(publication_model=publication_model(id),storage_model=id,cohort,state,n=length(ix),n_mrpe=length(nz),mrpe_percent=isempty(nz) ? NaN : 100mean(abs.((yp[nz].-ctx.y[nz])./ctx.y[nz])),scaled_mse=mean(abs2,(yp[ix].-ctx.y[ix])./ctx.sigma[ix])))
            end
        end
    end
    sort!(training;by=r->r.bic);sort!(scores;by=r->r.selection_nll)
    CSV.write(joinpath(out,"training_ranking.csv"),training)
    CSV.write(joinpath(out,"selection_ranking.csv"),scores)
    CSV.write(joinpath(out,"state_errors.csv"),errors)
    CSV.write(joinpath(out,"predictions.csv"),points)
    # Near-equal predictive scores: prefer fewer fitted parameters, then training BIC.
    best=first(scores).selection_nll
    eligible=filter(r->r.selection_nll<=best+1e-6,scores)
    chosen=first(sort(eligible;by=r->(r.k,r.training_bic,r.publication_model)))
    CSV.write(joinpath(out,"selection_lock.csv"),[merge(chosen,(criterion="held-out NLL; within 1e-6 prefer fewer parameters then training BIC",parameter_sha256=bytes2hex(sha256(read(joinpath(out,"model$(chosen.storage_model)","selected_parameters.csv"))))))])
    chosen
end
end
