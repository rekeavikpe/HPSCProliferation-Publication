using CSV, CairoMakie
Makie.MakieCore.func2string(func::F) where {F<:Function}=string(nameof(func))
include("../figures/ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1||error("Usage: historical_lactate.jl RESULTS")
root=abspath(ARGS[1]);isfile(joinpath(root,"COMPLETE"))||error("Fit incomplete")
rows=collect(CSV.File(joinpath(root,"trajectories.csv")));obs=collect(CSV.File(joinpath(root,"observations.csv")))
manuscript_theme!();labels=("Cells (cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate diameter (μm)")
for c in sort(unique(r.condition for r in rows))
    fig=Figure(size=(540,400));condition=filter(r->r.condition==c,rows)
    Label(fig[0,1:2],"Condition $c · $(first(condition).cohort)";fontsize=11)
    for (i,state) in enumerate(("cells","glucose","lactate","aggregate"))
        r=filter(r->r.state==state,condition);ax=Axis(fig[(i-1)÷2+1,(i-1)%2+1];xlabel="Time (days)",ylabel=labels[i])
        t=[v.time_h/24 for v in r]
        lines!(ax,t,[v.historical for v in r];color=c==4 ? SCALE_COLOR : COLORS[c],label="Fixed K_Lac = 30.405")
        lines!(ax,t,[v.reference for v in r];color=:black,linestyle=:dash,label="All-free reference")
        measured=filter(v->v.condition==c&&v.state==state,obs)
        errorbars!(ax,[v.day for v in measured],[v.mean for v in measured],[v.sd for v in measured];color=c==4 ? SCALE_COLOR : COLORS[c],whiskerwidth=4)
        scatter!(ax,[v.day for v in measured],[v.mean for v in measured];color=c==4 ? SCALE_COLOR : COLORS[c],marker=MARKERS[c],label="Experimental data")
    end
    Legend(fig[3,1:2],first(filter(x->x isa Axis,fig.content));orientation=:horizontal)
    Label(fig[4,1:2],"Historical-value sensitivity · ten other parameters reoptimized on 100 mL only";fontsize=7)
    export_figure(joinpath(root,"figures"),"condition$c",fig)
end
println("Saved historical-lactate comparison plots")
