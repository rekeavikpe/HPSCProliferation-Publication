using CSV,CairoMakie
include("ManuscriptStyle.jl");using .ManuscriptStyle
root=only(ARGS);manuscript_theme!();mkpath(joinpath(root,"figures"))
rr=collect(CSV.File(joinpath(root,"curves.csv")));oo=collect(CSV.File(joinpath(root,"observations.csv")))
states=["cells","glucose","lactate","aggregate"]
ylabels=["Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)"]
scenarios=["baseline","shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both","full"]
names=["Model 0","Shear: growth","Shear: aggregate","Shear: both","EDR: growth","EDR: aggregate","EDR: both","Model 6"]
for cohort in unique(String(r.cohort) for r in rr)
    cs=sort(unique(r.condition for r in rr if r.cohort==cohort))
    for group in ("key","all")
        ids=group=="key" ? [1,6,8] : collect(1:8)
        f=Figure(size=(1600,350length(cs)+100),fontsize=18)
        colors=[:gray45,:darkorange,:goldenrod,:sienna,:purple,:steelblue,:teal,:crimson]
        styles=[:dash,:solid,:dash,:dot,:dashdot,:solid,:dot,:solid]
        for (i,c) in enumerate(cs),(j,state) in enumerate(states)
            scale=j==1 ? 1e6 : 1.
            a=Axis(f[i,j],xlabel="Time (days)",ylabel=ylabels[j],title="($(Char('a'+(i-1)*4+j-1))) Condition $c",titlesize=18)
            for k in ids
                data=filter(r->r.cohort==cohort&&r.condition==c&&r.state==state&&r.scenario==scenarios[k],rr)
                lines!(a,[r.day for r in data],[r.value/scale for r in data];color=colors[k],linestyle=styles[k],linewidth=2.5)
            end
            obs=filter(r->r.cohort==cohort&&r.condition==c&&r.state==state&&r.scenario=="full",oo)
            errorbars!(a,[r.day for r in obs],[r.mean/scale for r in obs],[r.sd/scale for r in obs];color=:black,linewidth=2,whiskerwidth=9)
            scatter!(a,[r.day for r in obs],[r.mean/scale for r in obs];color=:black,markersize=8)
        end
        Legend(f[length(cs)+1,1:4],[LineElement(color=colors[k],linestyle=styles[k],linewidth=2.5) for k in ids],names[ids];orientation=:horizontal,nbanks=group=="all" ? 2 : 1,labelsize=17)
        stem=joinpath(root,"figures","day1_$(cohort)_$group")
        save(stem*".png",f;px_per_unit=2);save(stem*".pdf",f)
    end
end
