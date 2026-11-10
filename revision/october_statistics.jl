using CSV,Statistics,Distributions
root=normpath(joinpath(@__DIR__,".."));out=only(ARGS);mkpath(out)
series=collect(CSV.File(joinpath(root,"data/processed/series.csv")))
raw=filter(r->r.status=="measured",collect(CSV.File(joinpath(root,"data/processed/measurements.csv"))))
function ranks(x)
    order=sortperm(x);rr=zeros(length(x));ties=Int[];j=1
    while j<=length(x)
        k=j;while k<length(x)&&x[order[k+1]]==x[order[j]];k+=1;end
        rr[order[j:k]].=(j+k)/2;push!(ties,k-j+1);j=k+1
    end
    rr,ties
end
function kw(groups)
    x=vcat(groups...);rr,ties=ranks(x);n=length(x);offset=0;h=0.
    for g in groups
        ix=offset+1:offset+length(g);h+=sum(rr[ix])^2/length(g);offset+=length(g)
    end
    correction=1-sum(t^3-t for t in ties)/(n^3-n)
    stat=(12h/(n*(n+1))-3(n+1))/correction
    stat,ccdf(Chisq(length(groups)-1),max(0.,stat))
end
units=NamedTuple[]
for s in series
    s.volume_ml in (100,500)&&s.seeding_cells_ml==80000||continue
    for endpoint in ("aggregate_day5","cells_fold_day1_day5")
        state=startswith(endpoint,"aggregate") ? "aggregate" : "cells"
        a=[Float64(r.value) for r in raw if r.series_id==s.series_id&&r.day==1&&r.state==state]
        b=[Float64(r.value) for r in raw if r.series_id==s.series_id&&r.day==5&&r.state==state]
        length(b)==1&&(state=="aggregate"||length(a)==1&&only(a)>0)||continue
        value=state=="aggregate" ? only(b) : only(b)/only(a)
        push!(units,(series_id=s.series_id,condition=s.condition,volume_ml=s.volume_ml,rpm=s.rpm,endpoint,value))
    end
end
CSV.write(joinpath(out,"culture_values.csv"),units)
CSV.write(joinpath(out,"culture_eligibility.csv"),[(series_id=s.series_id,condition=s.condition,rpm=s.rpm,volume_ml=s.volume_ml,endpoint,included=any(r->r.series_id==s.series_id&&r.endpoint==endpoint,units)) for s in series if s.seeding_cells_ml==80000 for endpoint in ("aggregate_day5","cells_fold_day1_day5")])
tests=NamedTuple[];summaries=NamedTuple[]
for endpoint in ("aggregate_day5","cells_fold_day1_day5")
    groups=[[r.value for r in units if r.endpoint==endpoint&&r.volume_ml==100&&r.rpm==rpm] for rpm in (30,40,60)]
    all(length(g)>=2 for g in groups)||error("Insufficient cultures for $endpoint")
    h,p=kw(groups);push!(tests,(endpoint,comparison="100 mL: 30 vs 40 vs 60 RPM",method="Kruskal-Wallis, tie corrected, chi-square approximation",statistic=h,p_value=p))
    for (rpm,g) in zip((30,40,60),groups)
        push!(summaries,(endpoint,rpm,volume_ml=100,n=length(g),mean=mean(g),sd=std(g),median=median(g)))
    end
end
# Vessel question is a two-group comparison, not the three-RPM omnibus test.
groups=[[r.value for r in units if r.endpoint=="aggregate_day5"&&r.volume_ml==v&&r.rpm==30] for v in (100,500)]
h,p=kw(groups);push!(tests,(endpoint="aggregate_day5",comparison="100 vs 500 mL at 30 RPM",method="Kruskal-Wallis two-group rank comparison, tie corrected",statistic=h,p_value=p))
CSV.write(joinpath(out,"tests.csv"),tests);CSV.write(joinpath(out,"summary.csv"),summaries)
write(joinpath(out,"README.txt"),"Independent unit: complete culture. Fold expansion is day5/day1 within the same culture, not a ratio of pooled means. Complete paired observations required. Missing cultures are excluded and listed by comparison through culture_values.csv. Original condition numbering. Three-RPM comparisons restricted to 100 mL, 80000 cells/mL. 60 RPM handling differs. Tests are observational comparisons, not isolated causal effects. Omnibus p-values are not pairwise p-values.\n")
write(joinpath(out,"COMPLETE"),"Culture-level statistics computed\n")
# Pairwise rank comparisons accompany, but do not replace, the three-rate omnibus test.
pairs=NamedTuple[]
for endpoint in ("aggregate_day5","cells_fold_day1_day5")
    localrows=NamedTuple[]
    for (a,b) in ((30,40),(30,60),(40,60))
        x=[r.value for r in units if r.endpoint==endpoint&&r.volume_ml==100&&r.rpm==a]
        y=[r.value for r in units if r.endpoint==endpoint&&r.volume_ml==100&&r.rpm==b]
        rr,ties=ranks(vcat(x,y));nx=length(x);ny=length(y);n=nx+ny
        u=sum(rr[1:nx])-nx*(nx+1)/2
        variance=nx*ny/12*((n+1)-sum(t^3-t for t in ties)/(n*(n-1)))
        z=max(0.,abs(u-nx*ny/2)-.5)/sqrt(variance)
        push!(localrows,(endpoint,rpm_a=a,rpm_b=b,n_a=nx,n_b=ny,U=u,p_value=2ccdf(Normal(),z)))
    end
    order=sortperm([r.p_value for r in localrows]);adjusted=zeros(3);running=0.
    for (j,i) in enumerate(order);running=max(running,min(1.,(4-j)*localrows[i].p_value));adjusted[i]=running;end
    append!(pairs,[merge(r,(p_holm=adjusted[i],method="Mann-Whitney normal approximation, tie and continuity corrections; Holm within three contrasts")) for (i,r) in enumerate(localrows)])
end
CSV.write(joinpath(out,"pairwise_rank_tests.csv"),pairs)
inputs=collect(CSV.File(joinpath(root,"data/processed/exposure_summary_inputs.csv")))
correlations=NamedTuple[]
for state in (:shear_pa,:edr_m2_s3)
    local pairs=[(only(r for r in inputs if r.scenario=="eul_va"&&r.volume_ml==v&&r.rpm==rpm),only(r for r in inputs if r.scenario=="lag_p50"&&r.size_um==200&&r.volume_ml==v&&r.rpm==rpm)) for (v,rpm) in ((100,30),(100,40),(100,60),(500,30))]
    x=[getproperty(a,state) for (a,b) in pairs];y=[getproperty(b,state) for (a,b) in pairs]
    r=cor(x,y);t=r*sqrt(2/(1-r*r))
    push!(correlations,(state=string(state),n=4,pearson_r=r,p_value=2ccdf(TDist(2),abs(t)),spearman_r=cor(first(ranks(x)),first(ranks(y)))))
end
CSV.write(joinpath(out,"exposure_correlations.csv"),correlations)
