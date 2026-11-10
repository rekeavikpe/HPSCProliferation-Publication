using CSV, TOML, SHA, Random, Statistics, ProfileLikelihood, Optimization, NLopt, SciMLBase
include("PackageWorkflow.jl")
using .PackageWorkflow
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)>=2||error("Usage: historical_lactate.jl FIT OUT [STARTS=8]")
fitdir,out=abspath.(ARGS[1:2]);n=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 8
n>=2||error("At least two starts required")
ispath(out)&&error("Use a fresh output directory")
base=checked_context(fitdir);model=legacy_model(5)
p=read_parameters(joinpath(fitdir,"model5/selected_parameters.csv"),model)
base_nll=gaussian_nll(base.y,prediction(base,model,p),base.sigma)
isapprox(base_nll,only(CSV.File(joinpath(fitdir,"model5/summary.csv"))).nll;atol=1e-6,rtol=0)||error("Source fit differs")
cfg=deepcopy(base.config);cfg["lactate_inhibition"]=true;cfg["fixed"]["K_Lac"]=30.40495
ctx=context(cfg=cfg);names=collect(free_names(ctx,model));length(names)==10||error("Expected ten free parameters")
mkpath(out);open(io->TOML.print(io,cfg),joinpath(out,"calibration.toml"),"w")
open(io->TOML.print(io,Dict("source_fit"=>fitdir,"fixed_K_Lac"=>30.40495,"starts"=>n,"seed"=>20261004,"optimizer"=>"ProfileLikelihood.jl MLE, NLopt.LN_BOBYQA","interpretation"=>"Historical combined-data value tested as a sensitivity assumption; not independent validation")),joinpath(out,"settings.toml"),"w")
files=[@__FILE__,joinpath(@__DIR__,"PackageWorkflow.jl"),joinpath(@__DIR__,"Manifest.toml")]
append!(files,[joinpath(ROOT,f) for f in ("src/Workflow.jl","src/Calibration.jl","src/Proliferation.jl","data/processed/measurements.csv","data/processed/series.csv","data/processed/feeding.csv","data/processed/lagrangian_exposure.csv")])
append!(files,[joinpath(fitdir,f) for f in ("calibration.toml","observations.csv","model5/selected_parameters.csv","model5/summary.csv")])
CSV.write(joinpath(out,"source_hashes.csv"),[(path=abspath(f),sha256=bytes2hex(sha256(read(f)))) for f in files])
function loglik(z,_)
    try
        -gaussian_nll(ctx.y,prediction(ctx,model,decode(ctx,model,z)),ctx.sigma)
    catch e
        e isa InterruptException&&rethrow()
        -1e100
    end
end
rng=MersenneTwister(20261004);z=encode(ctx,model,p)
starts=[i==1 ? copy(z) : clamp.(z .+ (i<=4 ? .04 : .15).*randn(rng,length(z)),.001,.999) for i in 1:n]
sols=Vector{Any}(undef,n)
Threads.@threads :greedy for i in 1:n
    prob=make_problem(loglik,starts[i],names)
    sol=mle(prob,NLopt.LN_BOBYQA();maxiters=20000,xtol_rel=1e-10,ftol_rel=1e-12)
    isapprox(loglik(sol.mle,nothing),sol.maximum;atol=1e-6,rtol=0)||error("MLE mismatch")
    sols[i]=sol;println("Start $i NLL=$(-sol.maximum) status=$(sol.retcode)");flush(stdout)
    save_parameters(joinpath(out,"start$(i)_parameters.csv"),decode(ctx,model,sol.mle))
end
best=minimum(-s.maximum for s in sols)
valid=findall(s->SciMLBase.successful_retcode(s.retcode)&&-s.maximum<=best+1e-6,sols)
isempty(valid)&&error("No converged start near best objective; inspect individual starts")
j=valid[argmin([-sols[i].maximum for i in valid])];sol=sols[j];q=decode(ctx,model,sol.mle)
CSV.write(joinpath(out,"starts.csv"),[(start=i,nll=-s.maximum,retcode=string(s.retcode),stopping_met=SciMLBase.successful_retcode(s.retcode)) for (i,s) in enumerate(sols)])
save_parameters(joinpath(out,"mle_parameters.csv"),q)
CSV.write(joinpath(out,"fit_comparison.csv"),[(formulation="100mL all-free reference",nll=base_nll,free_parameters=11,K_Lac=p.K_Lac),(formulation="historical K_Lac fixed; all other parameters free",nll=-sol.maximum,free_parameters=10,K_Lac=q.K_Lac)])
CSV.write(joinpath(out,"likelihood_penalty.csv"),[(delta_nll=-sol.maximum-base_nll,twice_delta_nll=2*(-sol.maximum-base_nll),converged_near_best=length(valid),starts=n)])
metrics=NamedTuple[];obs=NamedTuple[];curves=NamedTuple[]
for (cohort,conditions) in (("100mL calibration",(1,2,3,5,6)),("500mL retrospective",(4,)))
    a=context(cfg=base.config,conditions=conditions);b=context(cfg=cfg,conditions=conditions)
    y0=prediction(a,model,p);y1=prediction(b,model,q)
    for state in STATES
        ix=findall(r->r.state==state,a.observations);nz=filter(i->a.y[i]!=0,ix)
        push!(metrics,(cohort,state,reference_mrpe=100mean(abs.((y0[nz].-a.y[nz])./a.y[nz])),historical_mrpe=100mean(abs.((y1[nz].-a.y[nz])./a.y[nz])),reference_standardized_rmse=sqrt(mean(((y0[ix].-a.y[ix])./a.sigma[ix]).^2)),historical_standardized_rmse=sqrt(mean(((y1[ix].-a.y[ix])./a.sigma[ix]).^2)),delta_nll=sum(((y1[ix].-a.y[ix])./a.sigma[ix]).^2-((y0[ix].-a.y[ix])./a.sigma[ix]).^2)/2))
    end
    append!(obs,[(;cohort,r...,reference=y0[i],historical=y1[i]) for (i,r) in enumerate(a.observations)])
    for c in conditions
        setup=a.conditions[c];ts=sort(unique(vcat(collect(0.:3.:last(setup.times)),setup.times,[f.time_h for f in setup.feeds])))
        sims=[simulate(model,pp,setup.u0,ts,setup.feeds;fresh_glucose=setup.glucose,exposure=setup.lookup,lactate_inhibition=true,reltol=1e-12,abstol=1e-13) for pp in (p,q)]
        for (i,t) in enumerate(ts),(k,state) in enumerate(STATES)
            push!(curves,(cohort,condition=c,state,time_h=t,phase=0,reference=sims[1].predictions[k,i],historical=sims[2].predictions[k,i]))
        end
        for (i,e) in enumerate(sims[1].events),(k,state) in enumerate(STATES)
            push!(curves,(cohort,condition=c,state,time_h=e.time_h,phase=1,reference=e.after[k],historical=sims[2].events[i].after[k]))
        end
    end
end
sort!(curves;by=r->(r.condition,r.state,r.time_h,r.phase))
CSV.write(joinpath(out,"state_errors.csv"),metrics);CSV.write(joinpath(out,"observations.csv"),obs);CSV.write(joinpath(out,"trajectories.csv"),curves)
CSV.write(joinpath(out,"boundaries.csv"),[(parameter=string(name),value=getproperty(q,name),normalized_log_value=v,at_bound=v<1e-4||v>1-1e-4) for (name,v) in zip(names,encode(ctx,model,q))])
write(joinpath(out,"COMPLETE"),"Historical-value sensitivity complete. No 500-mL measurements used for optimization.\n")
println("Complete: $out")
