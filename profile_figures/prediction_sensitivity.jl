using CSV, CairoMakie
Makie.MakieCore.func2string(func::F) where {F<:Function}=string(nameof(func))
include("../figures/ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1||error("Usage: prediction_sensitivity.jl INACTIVE_PROFILE_DIRECTORY")
root=abspath(ARGS[1]);rows=collect(CSV.File(joinpath(root,"prediction_sensitivity.csv")))
obs=collect(CSV.File(joinpath(root,"sensitivity_observations.csv")))
manuscript_theme!()
labels=("Cells (cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate diameter (μm)")
for c in sort(unique(r.condition for r in rows))
    fig=Figure(size=(540,410));condition=filter(r->r.condition==c,rows)
    Label(fig[0,1:2],"Condition $c · $(first(condition).cohort)";fontsize=11)
    for (i,state) in enumerate(("cells","glucose","lactate","aggregate"))
        r=filter(r->r.state==state,condition);ax=Axis(fig[(i-1)÷2+1,(i-1)%2+1];xlabel="Time (days)",ylabel=labels[i])
        t=[v.time_h/24 for v in r]
        band!(ax,t,[v.profile_min for v in r],[v.profile_max for v in r];color=(c==4 ? SCALE_COLOR : COLORS[c],.2),label="Profile-path envelope")
        lines!(ax,t,[v.inactive_mle for v in r];color=c==4 ? SCALE_COLOR : COLORS[c],label="Inactive inhibition")
        lines!(ax,t,[v.active_mle for v in r];color=:black,linestyle=:dash,label="Fitted inhibition")
        measured=filter(v->v.condition==c&&v.state==state,obs)
        errorbars!(ax,[v.day for v in measured],[v.mean for v in measured],[v.sd for v in measured];color=c==4 ? SCALE_COLOR : COLORS[c],whiskerwidth=4)
        scatter!(ax,[v.day for v in measured],[v.mean for v in measured];color=c==4 ? SCALE_COLOR : COLORS[c],marker=MARKERS[c],label="Experimental data")
    end
    Legend(fig[3,1:2],first(filter(x->x isa Axis,fig.content));orientation=:horizontal,nbanks=2)
    Label(fig[4,1:2],"Conditional sensitivity envelope; not a 95% prediction interval";fontsize=7)
    export_figure(joinpath(root,"prediction_figures"),"condition$c",fig)
end
println("Saved prediction sensitivity figures")
