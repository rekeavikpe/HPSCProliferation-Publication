using CairoMakie, CSV, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==3||error("Usage: render_targeted.jl PROFILE_DIRECTORY LACTATE_LIMIT_DIRECTORY OUTPUT_DIRECTORY")
function main()
    profiledir=ARGS[1]=="-" ? nothing : abspath(ARGS[1])
    limitdir,out=abspath.(ARGS[2:3])
    CairoMakie.activate!();manuscript_theme!()
    if profiledir!==nothing
    profilemeta=TOML.parsefile(joinpath(profiledir,"run.toml"))
    limitmeta=TOML.parsefile(joinpath(limitdir,"run.toml"))
    profilemeta["source_fit"]==limitmeta["source_fit"]||error("Profile and limit used different source fits")
    profilemeta["parameter"]=="K_Lac"&&profilemeta["model"]==5||error("This renderer expects Model 6 K_Lac")
    rows=sort(collect(CSV.File(joinpath(profiledir,"profile.csv")));by=r->r.fixed_value)
    comparison=collect(CSV.File(joinpath(limitdir,"comparison.csv")))
    original=only(filter(r->r.formulation=="bounded lactate inhibition",comparison))
    inactive=only(filter(r->r.formulation=="exact inactive inhibition",comparison))
    limit=2(inactive.nll-original.nll)
    CairoMakie.activate!();manuscript_theme!()
    fig=Figure(size=(540,300))
    for col in 1:2
        ax=Axis(fig[1,col],xlabel="Lactate inhibition constant (mM)",ylabel="2 × NLL difference from saved fit",xscale=log10,title=col==1 ? "Full profile" : "Near-optimum view")
        lines!(ax,[r.fixed_value for r in rows],[r.twice_delta_nll for r in rows];color="#0072B2")
        for r in rows
            scatter!(ax,[r.fixed_value],[r.twice_delta_nll];color="#0072B2",marker=r.stopping_met ? :circle : :xcross,markersize=5)
        end
        hlines!(ax,[0.];color=:gray70,linestyle=:dot)
        hlines!(ax,[limit];color="#D55E00",linestyle=:dash)
        col==2&&ylims!(ax,min(-.15,limit-.1),10)
    end
    Legend(fig[2,1:2],[MarkerElement(color="#0072B2",marker=:circle),MarkerElement(color="#0072B2",marker=:xcross),LineElement(color="#D55E00",linestyle=:dash)],["Convergence reported","Convergence not confirmed","Exact inactive limit"];orientation=:horizontal,tellwidth=false,labelsize=6.5)
    export_figure(out,"lactate_profile",fig)
    write(joinpath(out,"lactate_profile_caption.txt"),"Model 6 multistart lactate profile from the saved 100-mL calibration. Other parameters are reoptimized at each fixed value using three initial starts and neighboring-point starts. Crosses flag the selected nuisance fit when convergence was not confirmed. The near-optimum panel clips only the displayed vertical range; full-range values remain in the companion panel. The dashed line is the separately reoptimized exact inactive-inhibition limit. No confidence-interval threshold is asserted for this bounded, working-likelihood analysis.\nAlt text: A lactate parameter profile is shown over its complete range and near the best likelihood, with the exact inactive limit marked for comparison.\n")
    end
    errors=collect(CSV.File(joinpath(limitdir,"prediction_comparison.csv")))
    states=("cells","glucose","lactate","aggregate")
    fig=Figure(size=(510,270));ax=Axis(fig[1,1],xlabel="State variable",ylabel="Change in MRPE (percentage points)",xticks=(1:4,["Cells","Glucose","Lactate","Diameter"]))
    hlines!(ax,[0.];color=:gray60)
    for (cohort,color,offset) in (("100mL calibration","#0072B2",-.08),("500mL retrospective","#CC79A7",.08))
        rr=[only(filter(r->r.cohort==cohort&&r.state==state,errors)) for state in states]
        scatter!(ax,collect(1:4).+offset,[r.inactive_mrpe-r.original_mrpe for r in rr];color,markersize=6,label=cohort)
    end
    Legend(fig[2,1],ax;orientation=:horizontal,tellwidth=false)
    export_figure(out,"lactate_limit_errors",fig)
    write(joinpath(out,"lactate_limit_errors_caption.txt"),"Change in MRPE after reoptimizing the exact inactive-inhibition model, relative to the original bounded-inhibition Model 6. Negative values mean lower error. Values are percentage-point differences, not relative percent changes. The 500-mL evaluation is retrospective and did not enter reoptimization.\nAlt text: Small changes in prediction error are displayed for four state variables in the training and scale-transfer cohorts.\n")
    println("Saved targeted-check figures to $out")
end
main()
