using CSV, Statistics, Distributions
root=normpath(joinpath(@__DIR__,".."));out=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(root,"results/david_revision/rpm_statistics");mkpath(out)
series=collect(CSV.File(joinpath(root,"data/processed/series.csv")))
measured=filter(r->r.status=="measured",collect(CSV.File(joinpath(root,"data/processed/measurements.csv"))))
feeds=collect(CSV.File(joinpath(root,"data/processed/feeding.csv")))
# The RPM comparison is restricted to the common 100 mL / 80,000 cells mL^-1 design.
matched=filter(r->r.volume_ml==100&&r.seeding_cells_ml==80000,series)
@assert Set(r.rpm for r in matched)==Set([30,40,60])
units=NamedTuple[]
for info in matched,state in ("cells","glucose","lactate","aggregate"),endpoint in ("mean_days1_4","day5")
    days=endpoint=="day5" ? [5] : collect(1:4)
    vals=filter(r->r.series_id==info.series_id&&r.state==state&&r.day in days,measured)
    Set(r.day for r in vals)==Set(days)||continue
    @assert length(vals)==length(days)
    push!(units,(series_id=String(info.series_id),legacy_condition=info.condition,publication_condition=info.condition==5 ? 4 : info.condition==6 ? 5 : info.condition,rpm=info.rpm,state,endpoint,value=mean(r.value for r in vals)))
end
CSV.write(joinpath(out,"culture_endpoints.csv"),units)
rows=NamedTuple[]
for state in ("cells","glucose","lactate","aggregate"),endpoint in ("mean_days1_4","day5"),(a,b) in ((30,40),(30,60),(40,60))
    x=[r.value for r in units if r.rpm==a&&r.state==state&&r.endpoint==endpoint];y=[r.value for r in units if r.rpm==b&&r.state==state&&r.endpoint==endpoint]
    length(x)>=2&&length(y)>=2||error("Insufficient independent cultures")
    vx=var(x)/length(x);vy=var(y)/length(y);se=sqrt(vx+vy);df=(vx+vy)^2/(vx^2/(length(x)-1)+vy^2/(length(y)-1));t=(mean(x)-mean(y))/se
    p=2ccdf(TDist(df),abs(t));ci=quantile(TDist(df),.975)*se
    # Samples precede exchange, so events at the endpoint do not affect that observation.
    time=endpoint=="day5" ? 120 : 96
    signature(rpm)=[(r.time_h,r.exchange_fraction,r.additional_handling) for r in feeds if r.condition==first(filter(z->z.rpm==rpm,matched)).condition&&r.time_h<time]
    equal=signature(a)==signature(b)
    push!(rows,(state,endpoint,rpm_a=a,rpm_b=b,n_a=length(x),n_b=length(y),mean_a=mean(x),mean_b=mean(y),sd_a=std(x),sd_b=std(y),mean_difference=mean(x)-mean(y),difference_ci_lower=mean(x)-mean(y)-ci,difference_ci_upper=mean(x)-mean(y)+ci,welch_df=df,p_value=p,protocol_matched=equal,note=equal ? "Recorded feeding and handling matched before endpoint. Unrecorded batch effects remain possible." : "Recorded handling differs. This is an operating-condition comparison, not an isolated RPM effect."))
end
# All 24 reported contrasts belong to one declared family. Intervals above are unadjusted.
ord=sortperm([r.p_value for r in rows]);adj=zeros(length(rows));running=0.
for (j,i) in enumerate(ord);global running=max(running,min(1.,(length(rows)-j+1)*rows[i].p_value));adj[i]=running;end
CSV.write(joinpath(out,"welch_holm.csv"),[merge(r,(p_holm=adj[i],family_size=length(rows))) for (i,r) in enumerate(rows)])
write(joinpath(out,"README.txt"),"Experimental cultures, not daily observations, are the independent units. Days 1–4 require all four measurements. Day 5 is separate. Two-sided Welch tests and Holm adjustment over all 24 reported contrasts. CI columns are unadjusted 95% mean-difference intervals. No causal RPM claim for protocols with different handling. Missing endpoints are not imputed.\n")
write(joinpath(out,"COMPLETE"),"24 culture-level contrasts\n")
println("RPM statistics saved: $out")
