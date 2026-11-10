using Statistics
function replicate_report(cfgfile,out,ps)
    records=NamedTuple[]; points=NamedTuple[]
    inv=Dict(r.series_id=>r for r in inventory())
    for plan in ps
        prep=prepare_fold(cfgfile,plan)
        for id in (0,5)
            m=legacy_model(id);p=read_parameters(joinpath(out,plan.name,"model$id/selected_parameters.csv"),m)
            for sid in plan.selection
                ctx=pooled_context([sid],prep.cfg;fallback_glucose=Dict(c=>q.glucose for (c,q) in prep.train.conditions))
                yp=prediction(ctx,m,p)
                for (r,v) in zip(ctx.observations,yp)
                    push!(points,(fold=plan.name,experiment=plan.experiment,model=publication_model(id),series_id=sid,group=inv[sid].label,day=r.day,state=r.state,observed=r.mean,predicted=v))
                end
                for metric in metrics(ctx,yp)
                    push!(records,(fold=plan.name,experiment=plan.experiment,model=publication_model(id),series_id=sid,group=inv[sid].label,state=metric.state,n=metric.n,mrpe_percent=metric.mrpe_percent))
                end
            end
        end
    end
    CSV.write(joinpath(out,"replicate_prediction_errors.csv"),records)
    CSV.write(joinpath(out,"replicate_predictions.csv"),points)
    summary=NamedTuple[]
    for ex in unique(r.experiment for r in records), m in (0,6), state in (STATES...,"all"), group in vcat(["all_groups"],sort(unique(r.group for r in records)))
        rs=filter(r->r.experiment==ex&&r.model==m&&r.state==state&&(group=="all_groups"||r.group==group)&&isfinite(r.mrpe_percent),records)
        isempty(rs)&&continue
        v=[r.mrpe_percent for r in rs]
        push!(summary,(experiment=ex,model=m,group,state,n_replicates=length(v),mean_mrpe_percent=mean(v),sd_mrpe_percentage_points=length(v)>1 ? std(v) : NaN,variance_mrpe_percentage_points_squared=length(v)>1 ? var(v) : NaN))
    end
    CSV.write(joinpath(out,"replicate_error_summary.csv"),summary)
    folds=NamedTuple[]
    for ex in unique(r.experiment for r in records),m in (0,6),state in (STATES...,"all")
        rs=filter(r->r.experiment==ex&&r.model==m&&r.state==state&&isfinite(r.mrpe_percent),records)
        vals=[mean(r.mrpe_percent for r in rs if r.fold==f) for f in unique(r.fold for r in rs)]
        isempty(vals)&&continue
        push!(folds,(experiment=ex,model=m,state,n_folds=length(vals),mean_fold_mrpe_percent=mean(vals),variance_across_fold_means=length(vals)>1 ? var(vals) : NaN))
    end
    CSV.write(joinpath(out,"fold_error_summary.csv"),folds)
    println("Saved replicate_error_summary.csv and fold_error_summary.csv. Variances are descriptive, not independent-fold uncertainty estimates.")
end
