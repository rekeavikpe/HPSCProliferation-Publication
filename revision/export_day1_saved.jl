include("October20.jl")
length(ARGS)==4||error("Usage: export_day1_saved.jl FIT SPLIT COMPLETED OUTPUT")
fit,split,source,out=ARGS
prep,external=source_contexts(fit,split);mkpath(out)
curves=NamedTuple[];observed=NamedTuple[];metrics=NamedTuple[];folds=NamedTuple[];terms=NamedTuple[]
for s in EIGHT
    m=s=="baseline" ? 0 : 6
    p=read_parameters(joinpath(source,"day1",s,"selected_parameters.csv"),m)
    for (label,base) in (("training60",prep.train),("selection40",prep.held),("500mL_day1_initialized",external))
        ctx=dayone(route(base,s));yp=prediction(ctx,m,p)
        append!(metrics,report_metrics(ctx,m,p;tag=s,cohort=label))
        for (i,r) in enumerate(ctx.observations)
            push!(observed,(scenario=s,cohort=label,condition=r.condition,day=r.day,state=r.state,mean=r.mean,sd=r.sd,predicted=yp[i]))
        end
        for (c,q) in ctx.conditions
            times=sort(unique(vcat(collect(24.:.25:maximum(q.times)),q.times)))
            pp=merge(p,(K_d=p.K_d*exp(-2),))
            sim=simulate(m,pp,q.u0,times.-24,[(time_h=f.time_h-24,fraction=f.fraction) for f in q.feeds];fresh_glucose=q.glucose,exposure=q.lookup,reltol=1e-12,abstol=1e-13)
            for (j,t) in enumerate(times)
                for (k,state) in enumerate(STATES)
                    push!(curves,(scenario=s,cohort=label,condition=c,day=t/24,state,side="pre",value=sim.predictions[k,j]))
                end
                for e in sim.events
                    isapprox(e.time_h+24,t;atol=1e-9)||continue
                    for (k,state) in enumerate(STATES)
                        push!(curves,(scenario=s,cohort=label,condition=c,day=t/24,state,side="post",value=e.after[k]))
                    end
                end
            end
            for (i,r) in enumerate(ctx.observations)
                r.condition==c||continue
                j=findfirst(==(24r.day),times);k=findfirst(==(r.state),STATES)
                isapprox(sim.predictions[k,j],yp[i];rtol=1e-7,atol=1e-7)||error("Dense curve differs from scored prediction")
                if r.state=="cells"&&(r.day==5||24r.day==maximum(q.times))
                    push!(folds,(scenario=s,cohort=label,condition=c,day=r.day,observed=r.mean/q.u0[1],predicted=yp[i]/q.u0[1]))
                end
            end
        end
    end
end
CSV.write(joinpath(out,"curves.csv"),curves);CSV.write(joinpath(out,"observations.csv"),observed)
CSV.write(joinpath(out,"errors.csv"),metrics);CSV.write(joinpath(out,"folds.csv"),folds)
write(joinpath(out,"README.txt"),"Saved day-1 fits only. No optimization. Measured initial conditions are partition-specific. Day1 excluded from errors. Dense curves checked against scored predictions. 500mL curves use measured 500mL day1 initial states.\n")
