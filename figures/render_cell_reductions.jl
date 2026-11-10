using CSV,CairoMakie
include("ManuscriptStyle.jl");using .ManuscriptStyle
root=only(ARGS);manuscript_theme!()
for (dir,_,files) in walkdir(root)
    "REDUCTION_COMPLETE" in files||continue
    obs=collect(CSV.File(joinpath(dir,"observations.csv")));rs=collect(CSV.File(joinpath(dir,"curves.csv")))
    f=Figure(size=(1100,750),fontsize=17)
    for (i,state) in enumerate(("cells","glucose","lactate","aggregate"))
        scale=i==1 ? 1e6 : 1.
        label=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")[i]
        ax=Axis(f[(i-1)÷2+1,(i-1)%2+1],xlabel="Time (days)",ylabel=label,title="($(Char('a'+i-1)))")
        for (strategy,color,style) in (("Unchanged ESI-017",:gray45,:dash),("Recalibrated",:steelblue,:solid))
            a=filter(r->r.state==state&&r.strategy==strategy,rs)
            lines!(ax,[r.day for r in a],[r.value/scale for r in a];color,linestyle=style,linewidth=2.5)
        end
        a=filter(r->r.state==state,obs)
        errorbars!(ax,[r.day for r in a],[r.mean/scale for r in a],[r.sd/scale for r in a];color=:black,linewidth=2,whiskerwidth=9)
        scatter!(ax,[r.day for r in a],[r.mean/scale for r in a];color=:black,markersize=8)
    end
    Legend(f[3,1:2],[LineElement(color=:gray45,linestyle=:dash),LineElement(color=:steelblue)], ["Unchanged ESI-017","Reduced recalibration"];orientation=:horizontal)
    save(joinpath(dir,"fit.png"),f;px_per_unit=2);save(joinpath(dir,"fit.pdf"),f)
    audit=collect(CSV.File(joinpath(dir,"profiles/audit.csv")));f=Figure(size=(1100,700),fontsize=17)
    for (i,a) in enumerate(audit)
        values=sort(collect(CSV.File(joinpath(dir,"profiles",a.parameter,"profile.csv")));by=r->r.value)
        ax=Axis(f[(i-1)÷3+1,(i-1)%3+1],xlabel=String(a.parameter),ylabel="Relative log-likelihood",xscale=log10,title="$(a.parameter): $(a.status)",titlesize=13)
        lines!(ax,[r.value for r in values],[r.relative_loglik for r in values];linewidth=2)
        hlines!(ax,[-1.920729410347062];color=:gray45,linestyle=:dash);ylims!(ax,-10,1)
    end
    save(joinpath(dir,"profiles.png"),f;px_per_unit=2);save(joinpath(dir,"profiles.pdf"),f)
end
