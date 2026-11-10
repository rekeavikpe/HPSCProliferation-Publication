using CSV, CairoMakie, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
1<=length(ARGS)<=2 || error("Usage: render_cell_lines.jl FROZEN_DIRECTORY [RECALIBRATION_DIRECTORY]")
dir=abspath(ARGS[1]); comparison=length(ARGS)==2; target=comparison ? abspath(ARGS[2]) : dir;
comparison && !isfile(joinpath(target,"COMPLETE")) && error("Incomplete recalibration");isfile(joinpath(dir,"COMPLETE")) || error("Incomplete cell-line evaluation")
points=collect(CSV.File(joinpath(dir,"predictions.csv")))
curves=collect(CSV.File(joinpath(dir,"trajectories.csv")))
meta=TOML.parsefile(joinpath(dir,"analysis.toml"))
meta["parameters_reestimated"]==false || error("This renderer is for frozen predictions")
adapted=comparison ? collect(CSV.File(joinpath(target,"trajectories.csv"))) : []
manuscript_theme!()
states=("cells","glucose","lactate","aggregate")
labels=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
colors=(RGBf(.165,.290,.514),RGBf(.486,.196,.165))
# Shared physical scales include every curve and experimental error bar.
shared_limits=Dict{String,Tuple{Float64,Float64,Float64}}()
for (state,step) in zip(states,(.5,5.,5.,50.))
    factor=state=="cells" ? 1e6 : 1.
    values=Float64[r.value/factor for r in curves if r.model==5 && r.state==state]
    append!(values,[r.value/factor for r in adapted if r.state==state])
    for r in points
        r.model==5 && r.state==state || continue
        sd=ismissing(r.sd) ? 0. : r.sd
        append!(values,[(r.observed-sd)/factor,(r.observed+sd)/factor])
    end
    low=min(0.,floor(minimum(values)/step)*step)
    high=ceil((maximum(values)+.04*(maximum(values)-low))/step)*step
    shared_limits[state]=(low,high,step)
end

function render_layout(orientation)
    portrait=orientation==:portrait
    fig=Figure(size=portrait ? (600,900) : (1080,650))
    for (col,line) in enumerate(("585A1","SCTi003-A"))
    header=portrait ? fig[1,col] : fig[2col-1,1:4]
    Box(header;color=colors[col],strokevisible=false)
    Label(header,line*(col==1 ? "  |  50 rpm" : "  |  40 rpm");fontsize=13,color=:white,font=:bold,padding=(8,8,7,7),tellwidth=false)
    for (row,state) in enumerate(states)
        slot=portrait ? fig[row+1,col] : fig[2col,row]
        ax=Axis(slot,xlabel=(!portrait || row==4) ? "Time (days)" : "",ylabel=labels[row],xlabelsize=11,ylabelsize=10,xticklabelsize=9,yticklabelsize=9,xticks=0:5,title="($(Char(96+(portrait ? 2*(row-1)+col : 4*(col-1)+row))))",titlealign=:left,titlesize=12)
        trajectory=filter(r->r.model==5 && r.cell_line==line && r.state==state,curves)
        obs=filter(r->r.model==5 && r.cell_line==line && r.state==state,points)
        isempty(trajectory) && error("Missing trajectory for $line/$state")
        isempty(obs) && error("Missing observations for $line/$state")
        factor=state=="cells" ? 1e6 : 1.
        lines!(ax,[r.day for r in trajectory],[r.value/factor for r in trajectory];color=:gray40,linestyle=:dash,linewidth=1.3)
        if comparison
            adapted_line=filter(r->r.cell_line==line && r.state==state,adapted)
            isempty(adapted_line) && error("Missing recalibrated trajectory")
            lines!(ax,[r.day for r in adapted_line],[r.value/factor for r in adapted_line];color=colors[col],linewidth=1.3)
        end
        available=filter(r->!ismissing(r.sd),obs)
        errorbars!(ax,[r.day for r in available],[r.observed/factor for r in available],[r.sd/factor for r in available];color=colors[col],linewidth=.95,whiskerwidth=5.5)
        scatter!(ax,[r.day for r in obs],[r.observed/factor for r in obs];color=colors[col],marker=:circle,markersize=5.2)
        xlims!(ax,0,maximum(r.day for r in trajectory)+.12)
        low,high,step=shared_limits[state]
        ylims!(ax,low,high)
        ax.yticks=collect(low:step:high)
    end
end
legend_groups=Vector{Any}[]
label_groups=Vector{String}[]
for col in 1:2
    entries=Any[LineElement(color=:gray40,linestyle=:dash,linewidth=1.3)]
    legend_labels=["Prediction without recalibration"]
    if comparison
        push!(entries,LineElement(color=colors[col],linewidth=1.3))
        push!(legend_labels,"Recalibrated fit")
    end
    push!(entries,MarkerElement(color=colors[col],marker=:circle,markersize=5.2))
    push!(legend_labels,"Experimental data (±SD)")
    push!(legend_groups,entries)
    push!(label_groups,legend_labels)
end
legend_slot=portrait ? fig[6,1:2] : fig[5,1:4]
Legend(legend_slot,legend_groups,label_groups,["585A1","SCTi003-A"];
    orientation=:horizontal,nbanks=comparison ? 3 : 2,titleposition=:top,labelsize=10,titlesize=11,
    framevisible=true,framecolor=:black,framewidth=.7,
    patchsize=(28,12),rowgap=5,colgap=20,groupgap=24,
    halign=:center,padding=(12,12,8,8),tellwidth=false)
rowgap!(fig.layout,12);colgap!(fig.layout,20)
stem=comparison ? "cell_lines_comparison" : "cell_lines_frozen"
export_figure(joinpath(target,"figures"),stem*"_"*string(orientation),fig)
end
render_layout(:portrait)
render_layout(:landscape)
stem=comparison ? "cell_lines_comparison" : "cell_lines_frozen"
# Preserve the existing portrait filenames used by the server workflow.
for extension in ("png","pdf","svg","jpeg","tiff","eps")
    cp(joinpath(target,"figures",stem*"_portrait."*extension),joinpath(target,"figures",stem*"."*extension);force=true)
end
source=basename(meta["source_fit"])
description=comparison ? "Gray dashed curves predict using the ESI-017 parameter estimates without recalibration. Colored solid curves refit seven biological parameters to each cell line; these are in-sample fits, not validation. Consult summary.csv for convergence and boundaries." : "Predictions use the ESI-017 parameters without refitting."
write(joinpath(target,"index.html"),"""<!doctype html><meta charset="utf-8"><title>iPSC predictions and fits</title><style>body{font:17px system-ui;max-width:800px;margin:36px auto;line-height:1.5}img{width:100%}</style><h1>Model 6 transfer to two iPSC lines</h1><p>Calibration source: $source. $description The no-dextran-sulfate cohorts match the original analysis; the 585A1 Methods description needs reconciliation.</p><h2>Portrait</h2><img src="figures/$(stem)_portrait.svg"><p><a href="figures/$(stem)_portrait.pdf">Portrait PDF</a></p><h2>Landscape</h2><img src="figures/$(stem)_landscape.svg"><p><a href="figures/$(stem)_landscape.pdf">Landscape PDF</a></p><p>585A1 exposure at 50 rpm is interpolated between 40 and 60 rpm CFD. Initial glucose is measured; fresh-medium glucose is assumed equal to its initial mean. SD is experimental dispersion, not a prediction interval.</p><p><a href="figures/$stem.pdf">Vector PDF</a> · <a href="predictions.csv">Predictions and measurements</a> · <a href="analysis.toml">Analysis assumptions</a></p>""")
println("Saved portrait and landscape versions of $stem in six formats to $target/figures")
