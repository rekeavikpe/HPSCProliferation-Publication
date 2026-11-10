using CSV, CairoMakie
include("ManuscriptStyle.jl")
using .ManuscriptStyle

root=only(ARGS)
# Use the requested font explicitly. Do not silently substitute a server font.
for (variable,macname,family) in (("HPSC_FONT","Arial.ttf","Arial"),
                                 ("HPSC_BOLD_FONT","Arial Bold.ttf","Arial:style=Bold"))
    if !haskey(ENV,variable)
        macpath=joinpath("/System/Library/Fonts/Supplemental",macname)
        if isfile(macpath)
            ENV[variable]=macpath
        elseif Sys.which("fc-match")!==nothing
            format="%{family}\n%{file}"
            found=split(strip(read(`fc-match -f $format $family`,String)),'\n')
            length(found)>=2 && occursin("arial",lowercase(found[1])) && isfile(found[2]) &&
                (ENV[variable]=found[2])
        end
    end
    haskey(ENV,variable) && isfile(ENV[variable]) ||
        error("Arial is unavailable. Render the saved results locally or set $variable to the Arial font file. Numerical results are already saved.")
end
manuscript_theme!()
mkpath(joinpath(root,"figures"))
curves=collect(CSV.File(joinpath(root,"curves.csv")))
points=collect(CSV.File(joinpath(root,"predictions.csv")))
states=["cells","glucose","lactate","aggregate"]
labels=["Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)"]
formulations=["Original D0 initialization","Aggregate initialized at D1"]
colors=[:gray45,:steelblue]
styles=[:dash,:solid]
for cohort in unique(String(r.cohort) for r in points)
    cs=sort(unique(r.condition for r in points if r.cohort==cohort))
    f=Figure(size=(1600,310length(cs)+120),fontsize=26)
    for (j,state) in enumerate(states)
        scale=j==1 ? 1e6 : 1.
        obs=filter(r->r.cohort==cohort && r.state==state && r.formulation==formulations[1],points)
        lines=filter(r->r.cohort==cohort && r.state==state,curves)
        ymax=max(maximum((r.observed+r.sd)/scale for r in obs),maximum(r.value/scale for r in lines))
        for (i,c) in enumerate(cs)
            ax=Axis(f[i,j],xlabel=i==length(cs) ? "Time (days)" : "",ylabel=labels[j],
                    title="($(Char('a'+(i-1)*4+j-1))) Condition $c",titlesize=24,
                    xlabelsize=25,ylabelsize=25,xticklabelsize=21,yticklabelsize=21)
            for k in 1:2
                rr=filter(r->r.cohort==cohort && r.condition==c && r.state==state && r.formulation==formulations[k],curves)
                lines!(ax,[r.day for r in rr],[r.value/scale for r in rr];color=colors[k],linestyle=styles[k],linewidth=3)
            end
            oo=filter(r->r.condition==c,obs)
            errorbars!(ax,[r.day for r in oo],[r.observed/scale for r in oo],[r.sd/scale for r in oo];
                       color=:black,linewidth=2.3,whiskerwidth=11)
            scatter!(ax,[r.day for r in oo],[r.observed/scale for r in oo];color=:black,markersize=10)
            lastday=maximum(r.day for r in oo)
            xlims!(ax,-.15,lastday+.3)
            ylims!(ax,0,ymax*1.12)
        end
    end
    elements=[LineElement(color=colors[k],linestyle=styles[k],linewidth=3) for k in 1:2]
    Legend(f[length(cs)+1,1:4],vcat(elements,[MarkerElement(color=:black,marker=:circle,markersize=10)]),
           vcat(formulations,["Experimental mean ± SD"]);orientation=:horizontal,labelsize=24)
    stem=joinpath(root,"figures","aggregate_day1_$cohort")
    save(stem*".png",f;px_per_unit=2)
    save(stem*".pdf",f)
end
write(joinpath(root,"FIGURES_COMPLETE"),"Aggregate day 1 comparison figures saved.\n")
println("Comparison figures saved: ",joinpath(root,"figures"))
