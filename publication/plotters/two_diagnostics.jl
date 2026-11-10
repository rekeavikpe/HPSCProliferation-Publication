using Printf, CSV, TOML, ProfileLikelihood, Optimization, CairoMakie, LaTeXStrings, SciMLBase
# Makie 0.21 accesses a Julia internal removed in 1.12; use the public equivalent.
Makie.MakieCore.func2string(func::F) where {F<:Function} = string(nameof(func))
function render_profiles(input,out)
    input,out=abspath.((input,out));meta=TOML.parsefile(joinpath(input,"settings.toml"))
    meta["engine"]=="ProfileLikelihood.jl"||error("Expected ProfileLikelihood.jl results")
    isfile(joinpath(input,"COMPLETE"))||error("Profile suite is incomplete")
    original_order=["K_Agg","K_d","K_G","K_oG","K_Lac","Y_GLac","Y_GX","alpha","mu_max","K_tau","K_eps"]
    names=filter(n->n in meta["parameters"],original_order)
    length(names)==length(meta["parameters"])||error("Unsupported parameter labels")
    labelmap=Dict(zip(original_order,[L"K_{Agg}",L"K_d",L"K_G",L"K_{oG}",L"K_{Lac}",L"Y_{GLac}",L"Y_{GX}",L"\alpha",L"\mu_{\max}",L"K_{\tau}",L"K_{\varepsilon}"]))
    for (name,label) in get(meta,"parameter_latex",Dict())
        haskey(labelmap,name)&&(labelmap[name]=latexstring(label))
    end
    saved=Dict(String(r.parameter)=>Float64(r.value) for r in CSV.File(joinpath(input,"mle_parameters.csv")))
    mle_values=[saved[n] for n in names]
    # Restore the package plotting object from its exported results; no optimization.
    prob=LikelihoodProblem((theta,data)->error("Plot-only result adapter cannot fit data"),mle_values;syms=Symbol.(names))
    summary=only(CSV.File(joinpath(input,"mle_summary.csv")))
    code=getproperty(SciMLBase.ReturnCode,Symbol(summary.retcode))
    sol=ProfileLikelihood.LikelihoodSolution{length(names),typeof(mle_values),typeof(prob),Float64,typeof(code),Nothing}(mle_values,prob,nothing,-Float64(summary.nll),code)
    theta=Dict{Int,Vector{Float64}}();values=Dict{Int,Vector{Float64}}();splines=Dict{Int,Any}()
    intervals=Dict{Int,ProfileLikelihood.ConfidenceInterval{Float64,Float64}}();others=Dict{Int,Vector{Vector{Float64}}}()
    checks=NamedTuple[]
    for (i,name) in enumerate(names)
        rows=sort(collect(CSV.File(joinpath(input,name,"profile.csv")));by=r->r.parameter_value)
        theta[i]=Float64[r.parameter_value for r in rows];values[i]=Float64[r.profile_loglik for r in rows]
        ci=only(CSV.File(joinpath(input,name,"interval.csv")))
        lower_crossed=any(x<mle_values[i] && y<=-1.920729410347062 for (x,y) in zip(theta[i],values[i]))
        upper_crossed=any(x>mle_values[i] && y<=-1.920729410347062 for (x,y) in zip(theta[i],values[i]))
        intervals[i]=ProfileLikelihood.ConfidenceInterval(lower_crossed ? Float64(ci.package_lower) : minimum(theta[i]),upper_crossed ? Float64(ci.package_upper) : maximum(theta[i]),.95)
        ProfileLikelihood.spline_profile!(splines,i,theta[i],values[i],ProfileLikelihood.FritschCarlsonMonotonicInterpolation,ProfileLikelihood.Line)
        others[i]=Vector{Float64}[]
        push!(checks,(parameter=name,lower_crossed,upper_crossed,points=length(rows),lower_hits_bound=ci.lower_hits_bound,upper_hits_bound=ci.upper_hits_bound,unconfirmed_points=ci.unconfirmed_points,better_minimum=ci.better_minimum,lower_endpoint_profile_loglik=splines[i](Float64(ci.package_lower)),upper_endpoint_profile_loglik=splines[i](Float64(ci.package_upper))))
    end
    prof=ProfileLikelihood.ProfileLikelihoodSolution(theta,values,prob,sol,splines,intervals,others)
    CairoMakie.activate!()
    fig=ProfileLikelihood.plot_profiles(prof;
        latex_names=[labelmap[n] for n in names],show_mles=true,shade_ci=true,
        nrow=1,ncol=2,
        fig_kwargs=(figure_padding=(15,20,15,15),),
        axis_kwargs=(width=420,height=270,titlesize=23,xlabelsize=26,ylabelsize=26,xticklabelsize=21,yticklabelsize=21))
    axes=filter(ax->ax isa Axis,fig.content)
    for (i,ax) in enumerate(axes)
        font=get(ENV,"HPSC_FONT","Arial")
        ax.xlabelfont=font; ax.ylabelfont=font; ax.titlefont=font
        ax.xticklabelfont=font; ax.yticklabelfont=font
        parameter_labels=Dict("K_Lac"=>rich("K",subscript("Lac")), "Y_GX"=>rich("Y",subscript("GX")), "K_Agg"=>rich("K",subscript("Agg")), "K_d"=>rich("K",subscript("d")),
            "K_G"=>rich("K",subscript("G")), "K_oG"=>rich("K",subscript("oG")),
            "Y_GLac"=>rich("Y",subscript("GLac")), "alpha"=>"α", "mu_max"=>rich("µ",subscript("max")),
            "K_tau"=>rich("K",subscript("τ")), "K_eps"=>rich("K",subscript("ε")))
        ax.ylabel=rich("ℓ",subsup("p","*"),"(",parameter_labels[names[i]],") − ℓ",superscript("*");font=font)
        labels=Dict("K_Lac"=>rich("K",subscript("Lac")," (mM)"), "Y_GX"=>rich("Y",subscript("GX")," (µmol/cell)"), "K_Agg"=>rich("K",subscript("Agg")," (µm",superscript("3"),")"),
          "K_d"=>rich("K",subscript("d")," (h",superscript("−1"),")"),
          "K_G"=>rich("K",subscript("G")," (mM)"),
          "K_oG"=>rich("K",subscript("oG")," (mL/cell)"),
          "Y_GLac"=>rich("Y",subscript("GLac")," (mmol/mmol)"),
          "alpha"=>"α (µm)","mu_max"=>rich("µ",subscript("max")," (h",superscript("−1"),")"),
          "K_tau"=>rich("K",subscript("τ")," (Pa)"),
          "K_eps"=>rich("K",subscript("ε")," (m",superscript("2")," s",superscript("−3"),")"))
        ax.xlabel=labels[names[i]]
        ci=only(CSV.File(joinpath(input,names[i],"interval.csv")))
        function number(v)
            if names[i] in ("K_Agg","K_oG")
                exponent=floor(Int,log10(v))
                return rich(@sprintf("%.2f",v/10.0^exponent)," × 10",superscript(string(exponent));font=font)
            end
            names[i]=="K_G" ? @sprintf("%.2f",v) : @sprintf("%.3g",v)
        end
        ax.title=rich("($(Char(Int('a')+i-1))) 95% CI (",number(ci.package_lower),", ",number(ci.package_upper),")";font=font)
        ax.xticks=Makie.WilkinsonTicks(4)
        if names[i]=="K_Agg"
            ax.xticks=([1e8,2e8],[rich("1 × 10",superscript("8");font=font),rich("2 × 10",superscript("8");font=font)])
        elseif names[i]=="K_oG"
            ax.xticks=([2.8e-6,3.2e-6],[rich("2.8 × 10",superscript("−6");font=font),rich("3.2 × 10",superscript("−6");font=font)])
            ax.title=rich("(d) 95% CI (2.71, 3.27) × 10",superscript("−6");font=font)
        end
        if names[i]=="Y_GX"
            ax.xscale=log10
            ax.xticks=([1e-10,1e-8,1e-6,1e-5],[rich("10",superscript("−10");font=font),rich("10",superscript("−8");font=font),rich("10",superscript("−6");font=font),rich("10",superscript("−5");font=font)])
            xlims!(ax,8e-11,1.5e-5)
        elseif names[i]=="K_Lac"
            ax.xticks=[50.,75.,100.]
            xlims!(ax,47.,103.)
        end
        ylims!(ax,-2.7,.15)
        c=checks[i]
        if !c.lower_crossed || !c.upper_crossed
            side=!c.lower_crossed && !c.upper_crossed ? "both limits" : !c.lower_crossed ? "lower limit" : "upper limit"
            ax.title="($(Char(Int('a')+i-1))): $side unresolved"
        end
    end
    CairoMakie.resize_to_layout!(fig)
    mkpath(out)
    save(joinpath(out,"FigureS6_fixed_parameter_profiles.pdf"),fig;pt_per_unit=1)
    save(joinpath(out,"FigureS6_fixed_parameter_profiles.png"),fig;px_per_unit=3)

    CSV.write(joinpath(out,"profile_qualifications.csv"),checks)
    write(joinpath(out,"caption.txt"),"Generated by ProfileLikelihood.plot_profiles using the original manuscript script's styling (3-by-4 for eleven parameters; 3-by-3 for nine; 2-by-4 for eight), mathematical labels, default package splines, confidence-region shading, red reference lines and physical parameter axes. Existing package profile values and interval endpoints are reused without refitting. Where a computed profile does not cross the likelihood-ratio threshold, shading extends only to the tested endpoint and the title explicitly marks the corresponding limit unresolved. Raw package endpoints remain in the input CSVs. Unconfirmed nuisance optimizations require review. See profile_qualifications.csv and original fit diagnostics. Plot smoothing does not establish optimizer convergence.\n")
    if !isempty(get(meta,"fixed_parameters",Dict()))
        open(joinpath(out,"caption.txt"),"a") do io
            biological_fixed=Dict(n=>v for (n,v) in meta["fixed_parameters"] if !(n=="K_Lac"&&!get(meta,"lactate_inhibition",true)))
            !get(meta,"lactate_inhibition",true)&&println(io,"Lactate inhibition is exactly inactive; its unused K_Lac placeholder is not an estimated or fixed biological constant.")
            println(io,"Conditional profiles with fixed parameters: ",biological_fixed,". Intervals condition on these choices; consult the input settings and analysis provenance for the dataset used to estimate the fixed values.")
        end
    end
    println("Native ProfileLikelihood.plot_profiles figure saved: $out")
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2||error("Usage: --project=profile_figures profile_figures/render.jl PACKAGE_RESULTS OUTPUT")
    render_profiles(ARGS...)
end
