using CSV,TOML,Serialization,ProfileLikelihood,SciMLBase,SHA,Statistics
include("PackageWorkflow.jl");using .PackageWorkflow
include("../src/Workflow.jl");using .Workflow,.Workflow.Calibration,.Workflow.Calibration.Proliferation
include("../revision/CellContexts.jl")
length(ARGS)==2||error("Usage: run_october_cells.jl CELL_FIT OUTPUT")
fit,out=ARGS;cfg=TOML.parsefile(joinpath(fit,"calibration.toml"));line=strip(read(joinpath(fit,"cell_line.txt"),String));ctx=cell_context(line,cfg)
names=collect(free_names(ctx,6));p=read_parameters(joinpath(fit,"selected_parameters.csv"),6)
resolution=parse(Int,get(ENV,"HPSC_PROFILE_RESOLUTION","100"));mkpath(out)
fingerprint=bytes2hex(sha256(vcat(read(@__FILE__),read(joinpath(fit,"calibration.toml")),read(joinpath(fit,"selected_parameters.csv")),read(joinpath(@__DIR__,"PackageWorkflow.jl")),codeunits(string(resolution)*get(ENV,"HPSC_PROFILE_MAXITERS","10000")))))
lock=joinpath(out,"fingerprint.txt")
isfile(lock)&&strip(read(lock,String))!=fingerprint&&error("Changed profile source or resolution")
write(lock,fingerprint)
function ll(z,_)
    all(x->isfinite(x)&&0<=x<=1,z)||return -1e100
    try;-gaussian_nll(ctx.y,prediction(ctx,6,decode(ctx,6,z)),ctx.sigma)
    catch e;e isa InterruptException&&rethrow();-1e100;end
end
prob=make_problem(ll,encode(ctx,6,p),names)
solfile=joinpath(out,"mle.bin")
sol=if isfile(solfile);deserialize(solfile);else;s=fit_mle(prob);serialize(solfile,s);s;end
sol.maximum>=ll(encode(ctx,6,p),nothing)-1e-5||error("Profile MLE worse than fit")
save_parameters(joinpath(out,"mle_parameters.csv"),decode(ctx,6,sol.mle))
phys(i,z)=exp(log(cfg["bounds"][string(names[i])][1])+z*log(cfg["bounds"][string(names[i])][2]/cfg["bounds"][string(names[i])][1]))
Threads.@threads :greedy for i in eachindex(names)
    dir=joinpath(out,string(names[i]));isfile(joinpath(dir,"COMPLETE"))&&continue;mkpath(dir)
    println("Profiling $line $(names[i])");flush(stdout)
    other=setdiff(eachindex(names),[i])
    if isempty(other)
        grid=sort(unique(vcat(collect(range(0.,1.;length=2resolution+1)),[sol.mle[i]])))
        rows=[(parameter=string(names[i]),normalized=z,value=phys(i,z),relative_loglik=ll([z],nothing)-sol.maximum,verified_optimizer=true) for z in grid]
        CSV.write(joinpath(dir,"profile.csv"),rows)
        threshold=ProfileLikelihood.get_chisq_threshold(.95)
        left=any(r->r.normalized<sol.mle[i]&&r.relative_loglik<=threshold,rows)
        right=any(r->r.normalized>sol.mle[i]&&r.relative_loglik<=threshold,rows)
        inside=filter(r->r.relative_loglik>=threshold,rows)
        status=maximum(r.relative_loglik for r in rows)>1e-4 ? "numerical_review" : left&&right ? "bounded" : "not_bounded"
        CSV.write(joinpath(dir,"audit.csv"),[(parameter=string(names[i]),status,lower=minimum(r.value for r in inside),upper=maximum(r.value for r in inside),mle=phys(i,sol.mle[i]))])
        write(joinpath(dir,"COMPLETE"),status)
        continue
    end
    seeds=_ -> [copy(sol.mle[other]),fill(.25,length(other)),fill(.75,length(other))]
    prof,records=profile_parameter(sol.problem,sol,i;resolution,maxiters=parse(Int,get(ENV,"HPSC_PROFILE_MAXITERS","10000")),parallel=false,full_range=true,seed_provider=seeds)
    rows=NamedTuple[]
    for j in eachindex(prof.parameter_values[i])
        z=copy(sol.mle);z[i]=prof.parameter_values[i][j];z[other]=prof.other_mles[i][j]
        actual=ll(z,nothing)-sol.maximum
        isapprox(actual,prof.profile_values[i][j];atol=1e-5,rtol=0)||error("Profile likelihood verification failed")
        matches=filter(r->length(r.u)==length(other)&&maximum(abs.(r.u.-z[other]);init=0.)<1e-10&&abs(r.objective+actual)<1e-5,records)
        at_mle=maximum(abs.(z.-sol.mle))<1e-9
        ok=at_mle ? SciMLBase.successful_retcode(sol.retcode) : !isempty(matches)&&last(matches).stopping_met
        push!(rows,(parameter=string(names[i]),normalized=z[i],value=phys(i,z[i]),relative_loglik=actual,verified_optimizer=ok))
    end
    CSV.write(joinpath(dir,"profile.csv"),rows)
    threshold=ProfileLikelihood.get_chisq_threshold(.95)
    numerical=any(r->r.relative_loglik>1e-4||(!r.verified_optimizer&&r.relative_loglik>threshold-1),rows)
    left=any(r->r.normalized<sol.mle[i]&&r.relative_loglik<=threshold,rows)
    right=any(r->r.normalized>sol.mle[i]&&r.relative_loglik<=threshold,rows)
    status=numerical ? "numerical_review" : left&&right ? "bounded" : "not_bounded"
    ci=prof.confidence_intervals[i]
    CSV.write(joinpath(dir,"audit.csv"),[(parameter=string(names[i]),status,lower=phys(i,ci.lower),upper=phys(i,ci.upper),mle=phys(i,sol.mle[i]))])
    write(joinpath(dir,"COMPLETE"),status)
end
CSV.write(joinpath(out,"audit.csv"),vcat([collect(CSV.File(joinpath(out,string(n),"audit.csv"))) for n in names]...))
write(joinpath(out,"COMPLETE"),"Profiles computed. Audit distinguishes numerical failure from unbounded profiles.\n")
