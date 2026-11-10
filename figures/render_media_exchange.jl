using CSV, CairoMakie, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1 || error("Usage: render_media_exchange.jl EVALUATION_DIRECTORY")
dir=abspath(ARGS[1]);isfile(joinpath(dir,"COMPLETE")) || error("Incomplete prediction run")
points=collect(CSV.File(joinpath(dir,"predictions.csv")))
curves=collect(CSV.File(joinpath(dir,"trajectories.csv")))
meta=TOML.parsefile(joinpath(dir,"analysis.toml"))
manuscript_theme!()
fig=Figure(size=(760,580))
states=("cells","glucose","lactate","aggregate")
labels=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
colors=Dict(20=>:blue,50=>:red);markers=Dict(20=>:circle,50=>:rect)
for (i,state) in enumerate(states)
    ax=Axis(fig[(i-1)÷2+1,(i-1)%2+1],xlabel="Time (days)",ylabel=labels[i],xlabelsize=11,ylabelsize=11,xticklabelsize=10,yticklabelsize=10,xticks=0:6,title="($(('a':'d')[i]))",titlealign=:left,titlesize=12)
    factor=state=="cells" ? 1e6 : 1.
    for protocol in (20,50)
        line=filter(r->r.model==5 && r.protocol_percent==protocol && r.state==state,curves)
        obs=filter(r->r.model==5 && r.protocol_percent==protocol && r.state==state,points)
        lines!(ax,[r.day for r in line],[r.value/factor for r in line];color=colors[protocol],linewidth=1.3)
        available=filter(r->!ismissing(r.sd),obs)
        errorbars!(ax,[r.day for r in available],[r.observed/factor for r in available],[r.sd/factor for r in available];color=(colors[protocol],.9),linewidth=.95,whiskerwidth=5.5)
        scatter!(ax,[r.day for r in obs],[r.observed/factor for r in obs];color=colors[protocol],marker=markers[protocol],markersize=5.2)
    end
    xlims!(ax,0,6.15)
end
entries=[LineElement(color=colors[20],linewidth=1.3),MarkerElement(color=colors[20],marker=markers[20],markersize=6),LineElement(color=colors[50],linewidth=1.3),MarkerElement(color=colors[50],marker=markers[50],markersize=6)]
Legend(fig[3,1:2],entries,["20% — Model prediction","20% — Experimental data","50% — Model prediction","50% — Experimental data"];orientation=:horizontal,nbanks=2,labelsize=11)
rowgap!(fig.layout,16);colgap!(fig.layout,24)
export_figure(joinpath(dir,"figures"),"media_exchange",fig)
source=basename(meta["source_fit"])
write(joinpath(dir,"index.html"),"""<!doctype html><meta charset="utf-8"><title>Media-exchange evaluation</title><style>body{font:17px system-ui;max-width:1100px;margin:36px auto;line-height:1.5}img{width:100%}</style><h1>20% and 50% media exchange</h1><p>Frozen Model 6 from $source. No kinetic parameters were fitted to these experiments. Both protocols exchange 50% on day 2; the assigned fraction applies on days 3–5. Observations are pre-exchange.</p><img src="figures/media_exchange.svg"><p>Initial glucose uses measured day-zero protocol means. Fresh-medium glucose is assumed equal to those means. Error bars show experimental SD, not prediction uncertainty.</p><p><a href="state_errors.csv">Errors by state and protocol, including baseline</a> · <a href="figures/media_exchange.pdf">Vector PDF</a> · <a href="analysis.toml">Analysis assumptions</a></p>""")
println("Saved six figure formats and review page in $dir")
