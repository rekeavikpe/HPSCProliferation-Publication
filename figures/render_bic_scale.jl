using CairoMakie, CSV, Statistics, Random, TOML, FileIO, ImageIO
include("ManuscriptStyle.jl")
using .ManuscriptStyle

function main()
    2<=length(ARGS)<=3 || error("Usage: render_bic_scale.jl ORIGINAL_FIGURE_DATA OUTPUT [CFD_IMAGE]")
    source,out=abspath.(ARGS[1:2]);manuscript_theme!()
    with_cfd=length(ARGS)==3
    offset=with_cfd ? 2 : 0
    stem=with_cfd ? "cfd_bic_and_500mL" : "bic_and_500mL"
    meta=TOML.parsefile(joinpath(source,"source.toml"));selected=Int(meta["selected_model"])
    scope=get(meta,"bic_scope","Current all-100-mL analysis, not a 60/40 split.")
    prediction_source=get(meta,"prediction_source_run",meta["source_run"])
    runs=filter(r->parse(Int,replace(r.run,"run"=>""))<=meta["initial_starts_per_model"],collect(CSV.File(joinpath(source,"optimizer_bic.csv"))))
    curves=collect(CSV.File(joinpath(source,"test_trajectories.csv")))
    points=collect(CSV.File(joinpath(source,"test_predictions.csv")))
    fig=Figure(size=with_cfd ? (1950,720) : (1100,550))
    panel=fig.layout
    label_size=with_cfd ? 17 : 14
    tick_size=with_cfd ? 14 : 12
    panel_size=with_cfd ? 18 : 15
    if with_cfd
        prediction_grid=GridLayout(fig[1,3];tellwidth=false,tellheight=false)
        for (column,heading) in enumerate(("Hydrodynamic fields","Model selection","500 mL prediction"))
            Label(fig[0,column],heading;fontsize=22,font=:bold,tellwidth=false)
            colsize!(fig.layout,column,Relative((1/3, .26, 1-1/3-.26)[column]))
        end
        original=FileIO.load(abspath(ARGS[3]))
        image_axis=Axis(fig[1,1],aspect=DataAspect())
        image!(image_axis,rotr90(original))
        hidedecorations!(image_axis);hidespines!(image_axis)
        bic_slot=fig[1,2]
        legend_slot=fig[2,3]
        colgap!(fig.layout,50);rowgap!(fig.layout,24)
        colgap!(prediction_grid,24);rowgap!(prediction_grid,36)
    else
        Label(panel[0,1],"Model selection";fontsize=17,font=:bold)
        Label(panel[0,2:3],"500 mL prediction";fontsize=17,font=:bold)
        bic_slot=panel[1:2,1]
        legend_slot=panel[3,2:3]
    end
    ax=Axis(bic_slot,xlabel="Model",ylabel="BIC",xticks=(0:6,string.(0:6)),xgridvisible=false,xlabelsize=label_size,ylabelsize=label_size,xticklabelsize=tick_size,yticklabelsize=tick_size,title="($(Char('a'+offset)))",titlealign=:left,titlesize=panel_size)
    rng=MersenneTwister(20261005)
    stats=NamedTuple[]
    for pub in 0:6
        id=storage_model(pub);rr=filter(r->r.model==id,runs);v=Float64[r.bic for r in rr]
        q1,med,q3=quantile(v,[.25,.5,.75]);iqr=q3-q1
        lo=minimum(filter(x->x>=q1-1.5iqr,v));hi=maximum(filter(x->x<=q3+1.5iqr,v));color=MODEL_COLORS[id+1]
        poly!(ax,Point2f[(pub-.23,q1),(pub-.23,q3),(pub+.23,q3),(pub+.23,q1)];color=(color,.45),strokecolor=color,strokewidth=1)
        lines!(ax,[pub,pub],[lo,hi];color,linewidth=1)
        for (x,h) in ((lo,.1),(hi,.1),(med,.23))
            lines!(ax,[pub-h,pub+h],[x,x];color=x==med ? :black : color,linewidth=1.2)
        end
        for r in rr
            scatter!(ax,[pub+.22(rand(rng)-.5)],[r.bic];color=(color,.65),marker=:circle,markersize=4)
        end
        push!(stats,(publication_model=pub,storage_model=id,n=length(v),minimum=minimum(v),q1=q1,median=med,q3=q3,maximum=maximum(v)))
    end
    xlims!(ax,-.5,6.5)
    biclo,bichi=extrema(Float64[r.bic for r in runs])
    padding=max(5.,.12*(bichi-biclo))
    target=(bichi-biclo+2padding)/8
    magnitude=10.0^floor(log10(target))
    step=first(filter(x->x>=target,[1.,2.,5.,10.].*magnitude))
    lower=floor((biclo-padding)/step)*step
    upper=ceil((bichi+padding)/step)*step
    ylims!(ax,lower,upper);ax.yticks=collect(lower:step:upper)
    state_labels=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
    for (k,state) in enumerate(("cells","glucose","lactate","aggregate"))
        row,col=divrem(k-1,2);d=state=="cells" ? 1e6 : 1.
        slot=with_cfd ? prediction_grid[row+1,col+1] : panel[row+1,col+2]
        a=Axis(slot,xlabel="Time (days)",ylabel=state_labels[k],xticks=0:5,xlabelsize=label_size,ylabelsize=label_size,xticklabelsize=tick_size,yticklabelsize=tick_size,title="($(Char('a'+k+offset)))",titlealign=:left,titlesize=panel_size)
        for id in (0,selected)
            rr=filter(r->r.publication_model==id&&r.state==state,curves)
            lines!(a,[r.day for r in rr],[r.value/d for r in rr];color=id==0 ? :gray45 : SCALE_COLOR,linestyle=id==0 ? :dash : :solid,linewidth=1.6)
        end
        pp=filter(r->r.publication_model==selected&&r.state==state,points)
        errorbars!(a,[r.day for r in pp],[r.observed/d for r in pp],[r.sd/d for r in pp];color=SCALE_COLOR,whiskerwidth=8,linewidth=1.6)
        scatter!(a,[r.day for r in pp],[r.observed/d for r in pp];color=SCALE_COLOR,marker=MARKERS[4],markersize=7)
        xlims!(a,-.15,5.15)
    end
    Legend(legend_slot,[LineElement(color=:gray45,linestyle=:dash),LineElement(color=SCALE_COLOR),MarkerElement(color=SCALE_COLOR,marker=MARKERS[4])],["Model 0",model_label(selected),"Experimental data ± SD"];orientation=:horizontal,labelsize=with_cfd ? 16 : 12,tellwidth=false)
    if with_cfd
        colsize!(prediction_grid,1,Relative(.5))
        colsize!(prediction_grid,2,Relative(.5))
    else
        colsize!(panel,1,Relative(.38));colgap!(panel,30);rowgap!(panel,24)
    end
    export_figure(out,stem,fig)
    CSV.write(joinpath(out,"bic_box_statistics.csv"),stats)
    ranking=collect(CSV.File(joinpath(source,"training_ranking.csv")))
    convergence=map(0:6) do pub
        id=storage_model(pub);rr=filter(r->r.model==id,runs)
        selected_fit=only(filter(r->r.model==id,ranking))
        (publication_model=pub,storage_model=id,total_initial_runs=length(rr),
         stopping_met=count(r->r.stopping_met,rr),stopping_not_met=count(r->!r.stopping_met,rr),
         selected_fit_stopping_met=selected_fit.stopping_met,selected_fit_bic=selected_fit.bic)
    end
    CSV.write(joinpath(out,"supplementary_optimizer_convergence.csv"),convergence)
    body=join(["<tr><td>$(r.publication_model)</td><td>$(r.total_initial_runs)</td><td>$(r.stopping_met)</td><td>$(r.stopping_not_met)</td><td>$(r.selected_fit_stopping_met ? "Yes" : "No")</td></tr>" for r in convergence],"\n")
    write(joinpath(out,"supplementary_optimizer_convergence.html"),"""<!doctype html><meta charset="utf-8"><title>Optimizer convergence</title><style>body{font:16px Arial;max-width:900px;margin:50px auto;line-height:1.5}table{border-collapse:collapse;width:100%}th,td{padding:10px;border:1px solid #ccc;text-align:left}th{background:#edf2f5}</style><h1>Supplementary table: optimizer convergence</h1><table><tr><th>Model</th><th>Initial runs</th><th>Stopping criterion met</th><th>Not met</th><th>Selected fit: criterion met</th></tr>$body</table><p>Counts describe the saved stage chosen for each initial optimizer run. They are numerical diagnostics, not biological replicates. A finite score with an unmet stopping criterion is retained in the main figure using the same marker as other runs. The selected-fit status comes from the saved final ranking, which may include refinements. Meeting the criterion does not establish a global optimum or parameter identifiability.</p><p>Model 0 is included in both this table and the combined figure, which displays Models 0–6. $scope</p>""")
    write(joinpath(out,stem*"_caption.txt"),"$(with_cfd ? "(a,b) Original 100-mL CFD shear-stress and EDR panels at 30 and 60 rpm, reproduced unchanged from the supplied Figure1.jpg. CFD image remains raster; numerical panels are rendered from data. " : "")$scope ($(Char('a'+offset))) Vertical BIC boxplots for baseline Model 0 and six distinct hydrodynamic models, $(meta["initial_starts_per_model"]) optimizer runs per model; boxes use actual quartiles and observed 1.5-IQR whiskers without artificial width expansion. All finite saved starts are shown as circles regardless of stopping status; the supplementary optimizer-convergence table reports the breakdown. All seven models, including baseline Model 0, are displayed. ($(Char('b'+offset))–$(Char('e'+offset))) Retrospective 500-mL predictions for viable cells, glucose, lactate and aggregate size, respectively, using the all-100-mL calibrated parameters; baseline and selected hydrodynamic model are shown with experimental means and between-vessel SD. No 500-mL recalibration. BIC source: $(meta["source_run"]). Prediction source: $prediction_source.\n")
end
main()
