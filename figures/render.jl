using CairoMakie, CSV, Statistics, Random, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
function main()
CairoMakie.activate!();manuscript_theme!()
root=normpath(joinpath(@__DIR__,".."))
source=isempty(ARGS) ? joinpath(root,"reference/figure_data") : abspath(ARGS[1])
out=length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(@__DIR__,"output")
states=("cells","glucose","lactate","aggregate")
labels=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
points=collect(CSV.File(joinpath(source,"training_points.csv")))
curves=collect(CSV.File(joinpath(source,"training_curves.csv")))
source_meta=TOML.parsefile(joinpath(source,"source.toml"))
selected=source_meta["selected_model"]
source_run=get(source_meta,"source_run",basename(dirname(source)))
conditions=(1,2,3,5,6)
condition_labels=["Condition 1: 2×10⁴ cells/mL, 40 RPM","Condition 2: 8×10⁴ cells/mL, 40 RPM","Condition 3: 15×10⁴ cells/mL, 40 RPM","Condition 4: 8×10⁴ cells/mL, 30 RPM","Condition 5: 8×10⁴ cells/mL, 60 RPM"]
scale(state)=state=="cells" ? 1e6 : 1.
limits=Dict(state=>max((2.,18.,20.,360.)[k],1.08maximum(vcat([r.value/scale(state) for r in curves if r.state==state],[(r.mean+r.sd)/scale(state) for r in points if r.state==state]))) for (k,state) in enumerate(states))
for id in unique([0,selected])
    fig=Figure(size=(900,680))
    for (j,c) in enumerate(conditions),(i,state) in enumerate(states)
        ax=Axis(fig[i,j],xlabel=i==4 ? "Time (days)" : "",ylabel=j==1 ? labels[i] : "",xticks=0:6,xlabelsize=10,ylabelsize=10,xticklabelsize=8.5,yticklabelsize=8.5,title="($(Char(96+(i-1)*5+j)))",titlealign=:left,titlesize=10)
        rr=filter(r->r.model==id&&r.condition==c&&r.state==state,curves)
        pp=filter(r->r.condition==c&&r.state==state,points)
        color=COLORS[c]
        lines!(ax,[r.day for r in rr],[r.value/scale(state) for r in rr];color,linewidth=1.0)
        errorbars!(ax,[r.day for r in pp],[r.mean/scale(state) for r in pp],[r.sd/scale(state) for r in pp];color=(color,.9),whiskerwidth=5.5,linewidth=.95)
        scatter!(ax,[r.day for r in pp],[r.mean/scale(state) for r in pp];color,marker=MARKERS[c],markersize=5.2)
        xlims!(ax,-.15,maximum(r.day for r in rr)+.15);ylims!(ax,0,limits[state])
    end
    Legend(fig[5,1:5],[[LineElement(color=COLORS[c]),MarkerElement(color=COLORS[c],marker=MARKERS[c],markersize=5.2)] for c in conditions],condition_labels;orientation=:horizontal,nbanks=2,tellwidth=false,labelsize=11)
    colgap!(fig.layout,18);rowgap!(fig.layout,17)
    stem=id==0 ? "training_without_hydrodynamics" : "training_with_hydrodynamics"
    export_figure(out,stem,fig)
    write(joinpath(out,stem*"_caption.txt"),"100-mL calibration, Model $(publication_model(id)). Rows show viable cells, glucose, lactate and mean aggregate size. Columns are manuscript conditions 1–5 (legacy 1,2,3,5,6). Points are daily condition means; bars are between-bioreactor SD. Curves use saved fitted parameters. Metabolite jumps show exchange after sampling. Axis limits are shared between baseline and hydrodynamic figures. Day-zero states initialize the model; fitted observations start on day 1. Source calibration: $source_run.\nAlt text: Twenty panels show four measured culture outcomes and their fitted trajectories for five 100-mL conditions.\n")
end
if get(ENV,"HPSC_FIT_FIGURES_ONLY","0")!="1"
allstarts=collect(CSV.File(joinpath(source,"optimizer_bic.csv")))
ninitial=TOML.parsefile(joinpath(source,"source.toml"))["initial_starts_per_model"]
starts=filter(r->parse(Int,replace(r.run,"run"=>""))<=ninitial,allstarts)
refinements=filter(r->parse(Int,replace(r.run,"run"=>""))>ninitial,allstarts)
CSV.write(joinpath(out,"optimizer_refinements.csv"),refinements)
for view in (:full,:zoom)
fig=Figure(size=(720,440));ax=Axis(fig[1,1],xlabel="BIC value",ylabel="Models",yticks=(0:6,["Model $i" for i in 0:6]),ygridvisible=false,yreversed=true)
rng=MersenneTwister(20261003)
for id in publication_order()
    yid=publication_model(id)
    rows=filter(r->r.model==id,starts);v=[r.bic for r in rows];q1,q2,q3=quantile(v,[.25,.5,.75]);iqr=q3-q1
    lo=min(q1,minimum(filter(x->x>=q1-1.5iqr,v)));hi=max(q3,maximum(filter(x->x<=q3+1.5iqr,v)));color=MODEL_COLORS[id+1]
    poly!(ax,Point2f[(q1,yid-.25),(q3,yid-.25),(q3,yid+.25),(q1,yid+.25)];color=(color,.45),strokecolor=color,strokewidth=.8)
    lines!(ax,[lo,q1],[yid,yid];color);lines!(ax,[q3,hi],[yid,yid];color)
    for (x,height) in ((lo,.12),(hi,.12),(q2,.25))
        lines!(ax,[x,x],[yid-height,yid+height];color=x==q2 ? :black : color,linewidth=1)
    end
    for r in rows
        scatter!(ax,[r.bic],[yid+.13(rand(rng)-.5)];color=(color,.8),marker=r.stopping_met ? :circle : :xcross,markersize=3)
    end
end
ylims!(ax,-.6,6.6);ax.yreversed=true
Legend(fig[2,1],[MarkerElement(color=:gray,marker=:circle),MarkerElement(color=:gray,marker=:xcross)],["Stopping criterion met","Stopping criterion not met"];orientation=:horizontal,tellwidth=false)
if view==:zoom
    values=[r.bic for r in starts]
    low=minimum(values)
    upper=max(quantile(values,.90),maximum(minimum(r.bic for r in starts if r.model==id) for id in 0:6))
    upper=haskey(ENV,"HPSC_BIC_ZOOM_MAX") ? parse(Float64,ENV["HPSC_BIC_ZOOM_MAX"]) : upper
    upper>low || (upper=low+1)
    margin=.03*(upper-low)
    xlims!(ax,low-margin,upper+margin)
    clipped=[(model=id,above_zoom=count(r->r.model==id&&r.bic>upper+margin,starts),upper_limit=upper+margin) for id in 0:6]
    CSV.write(joinpath(out,"training_selection_zoom_counts.csv"),clipped)
    Label(fig[3,1],"Zoomed view: $(sum(r.above_zoom for r in clipped)) off-scale points; retained in full-range companion",fontsize=7,tellwidth=false)
    export_figure(out,"training_selection_zoom",fig)
else
    export_figure(out,"training_selection",fig)
end
end
write(joinpath(out,"training_selection_caption.txt"),"Training BIC across $ninitial initial calibration runs per model. Dependent refinements are excluded from these boxes and saved separately in optimizer_refinements.csv; the best-achieved ranking includes refinements. Each point uses the lower achieved NLL of its least-squares and polishing stages; crosses flag an unmet stopping criterion for that stage. The companion training_selection_zoom limits only the displayed x range, retaining boxes computed from all initial runs; its displayed note and CSV count off-scale points. Automatic zoom uses the pooled 90th percentile while retaining each model’s best initial fit, plus padding; HPSC_BIC_ZOOM_MAX overrides the upper focus bound. Boxes show actual quartiles, medians and 1.5-IQR whiskers; they may collapse when starts coincide. These are optimizer outcomes, not biological replicates or uncertainty intervals. Selection uses the best achieved BIC, not the box median. Duplicate historical Model 2 is removed.\nAlt text: Horizontal colored boxplots and jittered points compare training BIC across seven models; crosses mark stages that did not meet their stopping criterion.\n")
end
trajectories=collect(CSV.File(joinpath(source,"test_trajectories.csv")));pp=collect(CSV.File(joinpath(source,"test_predictions.csv")))
fig=Figure(size=(510,390));color=SCALE_COLOR
for (k,state) in enumerate(states)
    row,col=divrem(k-1,2)
    ax=Axis(fig[row+1,col+1],xlabel="Time (days)",ylabel=replace(labels[k],"\n"=>" "),xticks=0:5,title="($(Char(96+k)))",titlealign=:left)
    for id in unique([0,selected])
        rr=filter(r->r.publication_model==id&&r.state==state,trajectories)
        lines!(ax,[r.day for r in rr],[r.value/scale(state) for r in rr];color=id==0 ? :gray45 : color,linestyle=id==0 ? :dash : :solid)
    end
    data=filter(r->r.publication_model==selected&&r.state==state,pp)
    errorbars!(ax,[r.day for r in data],[r.observed/scale(state) for r in data],[r.sd/scale(state) for r in data];color=(color,.65),whiskerwidth=4,linewidth=.8)
    scatter!(ax,[r.day for r in data],[r.observed/scale(state) for r in data];color,marker=MARKERS[4])
    xlims!(ax,-.15,5.15)
end
Legend(fig[3,1:2],[LineElement(color=:gray45,linestyle=:dash),LineElement(color=color),MarkerElement(color=color,marker=MARKERS[4])],["Model 0","Model $(publication_model(selected))","Experimental data"];orientation=:horizontal,tellwidth=false)
export_figure(out,"scale_transfer",fig)
write(joinpath(out,"scale_transfer_caption.txt"),"Retrospective 500-mL predictions after 100-mL calibration. Solid dark teal: Model $(publication_model(selected)); dashed gray: Model 0; hexagons and bars: means and between-bioreactor SD. Dark teal distinguishes the 500-mL condition; the hexagon marker is retained. No refitting to 500-mL outcomes. Samples precede media exchange. Source calibration: $source_run.\nAlt text: Four panels compare baseline and hydrodynamic predictions with measured cell density, glucose, lactate and aggregate size over five days.\n")
if get(ENV,"HPSC_FIT_FIGURES_ONLY","0")!="1"
ranking=sort(collect(CSV.File(joinpath(source,"training_ranking.csv")));by=r->r.model)
fig=Figure(size=(510,280));ax=Axis(fig[1,1],xlabel="ΔBIC relative to best achieved fit",ylabel="Models",yticks=(0:6,["Model $i" for i in 0:6]),ygridvisible=false,yreversed=true)
best=minimum(r.bic for r in ranking)
for r in ranking
    delta=r.bic-best;color=MODEL_COLORS[Int(r.model)+1]
    lines!(ax,[0,delta],[publication_model(r.model),publication_model(r.model)];color,linewidth=1)
    scatter!(ax,[delta],[publication_model(r.model)];color,marker=r.stopping_met ? :circle : :xcross,markersize=6)
    text!(ax,delta+2,publication_model(r.model);text=string(round(r.bic;digits=2)),align=(:left,:center),fontsize=7)
end
xlims!(ax,-3,maximum(r.bic-best for r in ranking)+30);ylims!(ax,-.6,6.6);ax.yreversed=true
export_figure(out,"best_training_bic",fig)
write(joinpath(out,"best_training_bic_caption.txt"),"Best achieved BIC for each model; labels give absolute BIC and horizontal positions show ΔBIC. Crosses flag selected fits that did not meet the stopping criterion. This complements the full optimizer-start distribution; it does not discard or conceal poorly fitting starts.\nAlt text: Seven points compare best achieved model scores; Model $(publication_model(selected)) has the lowest BIC.\n")
end
println("Manuscript-style figures saved to $out")
end
main()
