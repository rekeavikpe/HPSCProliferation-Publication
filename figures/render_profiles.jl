using CairoMakie, CSV, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==2||error("Usage: render_profiles.jl PROFILE_SUITE OUTPUT_DIRECTORY")
suite,out=abspath.(ARGS); meta=TOML.parsefile(joinpath(suite,"suite.toml"))
parameters=collect(CSV.File(joinpath(meta["source_fit"],"model$(meta["model"])/selected_parameters.csv")))
manuscript_theme!(); CairoMakie.activate!()
fig=Figure(size=(850,650))
for (i,p) in enumerate(parameters)
    rows=sort(collect(CSV.File(joinpath(suite,p.parameter,"profile.csv")));by=r->r.fixed_value)
    ax=Axis(fig[div(i-1,4)+1,mod(i-1,4)+1],title=p.parameter,xlabel="Parameter value",ylabel="Normalized log likelihood",xscale=log10)
    lines!(ax,[r.fixed_value for r in rows],[-r.twice_delta_nll/2 for r in rows];color="#0072B2")
    for good in (true,false)
        selected=filter(r->r.stopping_met==good,rows)
        isempty(selected)&&continue
        scatter!(ax,[r.fixed_value for r in selected],[-r.twice_delta_nll/2 for r in selected];marker=good ? :circle : :xcross,color=good ? "#0072B2" : "#D55E00",markersize=4)
    end
    hlines!(ax,[-3.841458820694124/2];color=:black,linestyle=:dash)
    vlines!(ax,[p.value];color=:gray50,linestyle=:dot)
    ylims!(ax,-5,max(.2,maximum(-r.twice_delta_nll/2 for r in rows)+.1))
end
Legend(fig[4,1:4],[LineElement(color=:black,linestyle=:dash),LineElement(color=:gray50,linestyle=:dot),MarkerElement(color="#D55E00",marker=:xcross)],["Nominal 95% likelihood-ratio reference","Saved calibration estimate","Stopping criterion not met"];orientation=:horizontal,tellwidth=false)
export_figure(out,"model$(meta["model"])_profiles",fig)
write(joinpath(out,"profile_caption.txt"),"One-parameter profiles, with all remaining parameters reoptimized using multiple starts. Logarithmic parameter axes; normalized log likelihood follows the original manuscript convention. The lower display limit is -5; full numerical curves remain in the source CSV files. The dashed line is a nominal asymptotic reference, not a guarantee of 95% coverage under bounds and correlated observations. Dotted lines indicate saved calibration estimates. Crosses mark unsuccessful stopping checks; search bounds must not be interpreted as confidence limits.\n")
