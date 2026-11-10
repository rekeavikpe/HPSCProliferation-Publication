using CSV, DataFrames, TOML, SHA, Statistics, ProfileLikelihood, Optimization, NLopt, SciMLBase
include("../profiles/PackageWorkflow.jl")
using .PackageWorkflow
const ROOT=normpath(joinpath(@__DIR__,".."))
const SOURCE=read(joinpath(@__DIR__,"original_fixed3.jl"),String)
const FIXED_SOURCE=length(ARGS)>=3 ? abspath(ARGS[3]) : joinpath(@__DIR__,"inputs/new_100mL_selected_parameters.csv")
const NEW=Dict(String(r.parameter)=>Float64(r.value) for r in CSV.File(FIXED_SOURCE))
const CASES=[
    (name="A_original",parent="",change="Original six-condition script, including focused K_G profile range",hundred=false,newfixed=false,release=false,kg=40.,hydro=false,clamp=false,focused=true),
    (name="A_fullrange",parent="A_original",change="Profile grid uses complete original bounds",hundred=false,newfixed=false,release=false,kg=40.,hydro=false,clamp=false,focused=false),
    (name="B_100mL",parent="A_fullrange",change="Remove 500-mL condition only",hundred=true,newfixed=false,release=false,kg=40.,hydro=false,clamp=false,focused=false),
    (name="C_fixed_values",parent="B_100mL",change="Replace three fixed values with new 100-mL MLEs",hundred=true,newfixed=true,release=false,kg=40.,hydro=false,clamp=false,focused=false),
    (name="D_release_KAgg",parent="B_100mL",change="Release K_Agg; retain original fixed K_Lac and Y_GX",hundred=true,newfixed=false,release=true,kg=40.,hydro=false,clamp=false,focused=false),
    (name="E_KG400",parent="B_100mL",change="K_G upper bound 40 to 400",hundred=true,newfixed=false,release=false,kg=400.,hydro=false,clamp=false,focused=false),
    (name="F_hydro_bounds",parent="B_100mL",change="K_tau and K_eps upper bounds to 5",hundred=true,newfixed=false,release=false,kg=40.,hydro=true,clamp=false,focused=false),
    (name="G_exposure_clamp",parent="B_100mL",change="Clamp exposures below 200 um instead of extrapolating",hundred=true,newfixed=false,release=false,kg=40.,hydro=false,clamp=true,focused=false)]
for (name,change) in [("H_KLac_value","Replace only fixed K_Lac with the new 100-mL MLE"),("I_YGX_value","Replace only fixed Y_GX with the new 100-mL MLE"),("J_KAgg_value","Replace only fixed K_Agg with the new 100-mL MLE")]
    push!(CASES,merge(CASES[3],(;name,change,parent="B_100mL",newfixed=true)))
end

function build_case(case)
    firstpos=first(findfirst("const FIXED_K_Lac",SOURCE));lastpos=last(findfirst("const N_PARAMS = length(START)",SOURCE))
    text=SOURCE[firstpos:lastpos]
    case.hundred&&(text=replace(text,"for cond_idx in 1:NUM_CONDITIONS"=>"for cond_idx in (1,2,3,5,6)"))
    if case.newfixed
        replacements=Dict("K_Lac"=>("const FIXED_K_Lac = 30.40495","const FIXED_K_Lac = $(NEW["K_Lac"])"),"Y_GX"=>("const FIXED_Y_GX  = 3.620070e-06","const FIXED_Y_GX = $(NEW["Y_GX"])"),"K_Agg"=>("const FIXED_K_Agg = 6.174043e7","const FIXED_K_Agg = $(NEW["K_Agg"])"))
        subset=get(Dict("H_KLac_value"=>["K_Lac"],"I_YGX_value"=>["Y_GX"],"J_KAgg_value"=>["K_Agg"]),case.name,["K_Lac","Y_GX","K_Agg"])
        for name in subset
            before,after=replacements[name]
            occursin(before,text)||error("Original fixed-parameter definition changed")
            text=replace(text,before=>after)
        end
    end
    if case.release
        text=replace(text,"K_d, K_G, K_oG, Y_GLac, alpha, mu_max, K_tau, K_eps = p"=>"K_Agg, K_d, K_G, K_oG, Y_GLac, alpha, mu_max, K_tau, K_eps = p","        K_Agg = FIXED_K_Agg"=>"", "const START = ["=>"const START = [6.174043e7,", "const LB = ["=>"const LB = [1e5,", "const UB = ["=>"const UB = [1e10,", "const PARAM_NAMES = ["=>"const PARAM_NAMES = [\"K_Agg\",")
    end
    text=replace(text,"40.0,    # K_G"=>"$(case.kg),    # K_G")
    case.hydro&&(text=replace(text,"0.1,     # K_tau"=>"5.0,     # K_tau","0.05     # K_eps"=>"5.0     # K_eps"))
    case.clamp&&(text=replace(text,"return max(ys[1] + slope * (x - xs[1]), 1e-8)"=>"return ys[1]"))
    mod=Module(Symbol(case.name))
    Base.include_string(mod,"using CSV, DataFrames, Statistics, OrdinaryDiffEqVerner, DiffEqCallbacks, SciMLBase, XLSX, LaTeXStrings, Printf\ndatadir(args...)=joinpath($(repr(joinpath(@__DIR__,"inputs"))),last(args))\n"*text,case.name*".jl")
    Base.invokelatest(Base.invokelatest(getproperty,mod,:cache_experimental_data),CSV.read(joinpath(@__DIR__,"inputs/Combined100mLand500mLDataAllRPM.csv"),DataFrame)) |> x->(Base.invokelatest(getproperty,mod,:EXP_DATA_CACHE)[]=x)
    mod,text
end
function run_case(case,out)
    dest=joinpath(out,case.name);isfile(joinpath(dest,"COMPLETE"))&&return
    mkpath(dest);mod,text=build_case(case);write(joinpath(dest,"effective_original_definitions.jl"),text)
    names=copy(Base.invokelatest(getproperty,mod,:PARAM_NAMES));lo=copy(Base.invokelatest(getproperty,mod,:LB));hi=copy(Base.invokelatest(getproperty,mod,:UB));start=copy(Base.invokelatest(getproperty,mod,:START))
    ll=let fn=Base.invokelatest(getproperty,mod,:loglik)
        (p,data)->Base.invokelatest(fn,p,data)
    end
    settings=Dict("engine"=>"ProfileLikelihood.jl","model"=>5,"parameters"=>["K_G","mu_max"],"case"=>case.name,"parent"=>case.parent,"change"=>case.change,"coordinates"=>"original physical parameters","source_fit"=>"original fixed-three script reproduction","replacement_fixed_values_source"=>FIXED_SOURCE,"fixed_parameters"=>Dict("K_Lac"=>Base.invokelatest(getproperty,mod,:FIXED_K_Lac),"Y_GX"=>Base.invokelatest(getproperty,mod,:FIXED_Y_GX)))
    !case.release&&(settings["fixed_parameters"]["K_Agg"]=Base.invokelatest(getproperty,mod,:FIXED_K_Agg))
    open(io->TOML.print(io,settings),joinpath(dest,"settings.toml"),"w")
    candidates=Any[];starts=NamedTuple[]
    for f in (1.,.5,1.5,2.)
        seed=clamp.(start.*f,lo,hi)
        prob=LikelihoodProblem(ll,seed;data=nothing,syms=names,f_kwargs=(adtype=Optimization.AutoFiniteDiff(),),prob_kwargs=(lb=lo,ub=hi))
        sol=mle(prob,NLopt.LN_BOBYQA();maxiters=20000,xtol_rel=1e-10,ftol_rel=1e-12)
        push!(candidates,sol);push!(starts,(factor=f,nll=-ll(sol.mle,nothing),retcode=string(sol.retcode)))
        println(case.name," MLE start ",f," NLL ",last(starts).nll);flush(stdout)
    end
    best=argmin([r.nll for r in starts]);sol=candidates[best];prob=sol.problem
    CSV.write(joinpath(dest,"mle_starts.csv"),starts)
    CSV.write(joinpath(dest,"mle_parameters.csv"),[(parameter=n,value=v) for (n,v) in zip(names,sol.mle)])
    CSV.write(joinpath(dest,"mle_summary.csv"),[(model=5,nll=-sol.maximum,retcode=string(sol.retcode),stopping_met=SciMLBase.successful_retcode(sol.retcode))])
    CSV.write(joinpath(dest,"bounds.csv"),[(parameter=n,lower=l,upper=h) for (n,l,h) in zip(names,lo,hi)])
    ranges=construct_profile_ranges(sol,lo,hi,fill(100,length(names)))
    gi=findfirst(==("K_G"),names)
    if case.focused
        ranges[gi]=(LinRange(sol.mle[gi],2.,50),LinRange(sol.mle[gi],30.,50))
    end
    summary=NamedTuple[]
    for i in [gi,findfirst(==("mu_max"),names)]
        sub=joinpath(dest,names[i]);mkpath(sub);alg=AuditedBOBYQA()
        prof=ProfileLikelihood.profile(prob,sol,[i];param_ranges=ranges,parallel=true,conf_level=.95,alg,next_initial_estimate_method=:prev,progress=false)
        rows=NamedTuple[]
        for j in eachindex(prof.parameter_values[i])
            q=copy(sol.mle);q[i]=prof.parameter_values[i][j];others=prof.other_mles[i][j];q[setdiff(eachindex(q),[i])].=others
            actual=ll(q,nothing)-sol.maximum
            isapprox(actual,prof.profile_values[i][j];atol=1e-5,rtol=0)||error("Profile verification failed")
            matches=filter(r->r.u==others&&isapprox(r.objective,-actual;atol=1e-6,rtol=0),alg.records)
            stop=isempty(matches) ? missing : last(matches).stopping_met
            push!(rows,(parameter=names[i],parameter_value=q[i],profile_loglik=actual,stopping_met=stop))
        end
        CSV.write(joinpath(sub,"profile.csv"),rows)
        ci=prof.confidence_intervals[i]
        row=(parameter=names[i],mle=sol.mle[i],package_lower=ci.lower,package_upper=ci.upper,lower_hits_bound=isapprox(ci.lower,lo[i];rtol=1e-6),upper_hits_bound=isapprox(ci.upper,hi[i];rtol=1e-6),unconfirmed_points=count(r->r.stopping_met!==true,rows),better_minimum=maximum(r.profile_loglik for r in rows)>1e-5)
        CSV.write(joinpath(sub,"interval.csv"),[row]);push!(summary,row)
    end
    CSV.write(joinpath(dest,"profile_confidence_intervals.csv"),summary);write(joinpath(dest,"COMPLETE"),"Diagnostic original-script comparison, not validation\n")
end
out=length(ARGS)>=1 ? abspath(ARGS[1]) : joinpath(ROOT,"results/controlled_original")
mkpath(out)
CSV.write(joinpath(out,"design.csv"),CASES)
hashes=[(path=relpath(p,ROOT),sha256=bytes2hex(sha256(read(p)))) for p in [@__FILE__,joinpath(@__DIR__,"original_fixed3.jl"),joinpath(@__DIR__,"Manifest.toml"),joinpath(@__DIR__,"../profiles/PackageWorkflow.jl"),joinpath(@__DIR__,"inputs/Combined100mLand500mLDataAllRPM.csv"),joinpath(@__DIR__,"inputs/Particles_REBUILT_ALLNONAN_fromTXT.xlsx"),FIXED_SOURCE]]
hashfile=joinpath(out,"source_hashes.csv")
if isfile(hashfile)
    previous=Dict(String(r.path)=>String(r.sha256) for r in CSV.File(hashfile))
    previous==Dict(r.path=>r.sha256 for r in hashes)||error("Code or inputs changed. Choose a fresh output directory to preserve the previous comparison.")
end
CSV.write(hashfile,hashes)
selected=length(ARGS)>=2 ? filter(c->c.name==ARGS[2],CASES) : CASES
isempty(selected)&&error("Unknown comparison case")
Threads.@threads :greedy for case in selected
    run_case(case,out)
end
println("Completed controlled cases")
