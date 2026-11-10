using CSV, CairoMakie, TOML
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1 || error("Usage: render_exposure_sensitivity.jl RESULTS")
dir=abspath(only(ARGS));isfile(joinpath(dir,"COMPLETE")) || error("Incomplete exposure comparison")
rank=collect(CSV.File(joinpath(dir,"ranking.csv")));runs=collect(CSV.File(joinpath(dir,"runs.csv")));errors=collect(CSV.File(joinpath(dir,"state_errors.csv")))
settings=TOML.parsefile(joinpath(dir,"analysis.toml"))
labels=Dict("lag_p50"=>"Lagrangian P50","lag_mean"=>"Lagrangian mean","lag_p25"=>"Lagrangian P25","lag_p75"=>"Lagrangian P75","lag_p90"=>"Lagrangian P90","lag_max"=>"Lagrangian maximum","eul_va"=>"Eulerian volume average","eul_p50"=>"Eulerian P50","eul_p75"=>"Eulerian P75","eul_p90"=>"Eulerian P90","eul_max"=>"Eulerian maximum")
manuscript_theme!()
fig=Figure(size=(760,510))
ax=Axis(fig[1,1],xlabel="ΔBIC relative to the best exposure summary",yticks=(1:length(rank),[labels[r.scenario] for r in rank]),xlabelsize=12,yticklabelsize=11,xticklabelsize=10)
for (i,r) in enumerate(rank)
    color=startswith(r.scenario,"lag_") ? "#245B83" : "#AA4E35"
    group=filter(x->x.scenario==r.scenario,runs)
    x=[g.bic-first(rank).bic for g in group]
    if length(group)>1
        boxplot!(ax,fill(i,length(group)),x;orientation=:horizontal,color=(color,.3),show_outliers=false,width=.5)
    end
    for g in group
        scatter!(ax,[g.bic-first(rank).bic],[i];color=(color,.6),marker=g.stopping_met ? :circle : :xcross,markersize=3.5)
    end
    scatter!(ax,[r.delta_bic],[i];color=color,marker=:diamond,markersize=8)
end
ylims!(ax,.5,length(rank)+.5);ax.yreversed=true
Legend(fig[2,1],[MarkerElement(color=:black,marker=:diamond,markersize=8),MarkerElement(color=:gray,marker=:circle,markersize=5),MarkerElement(color=:gray,marker=:xcross,markersize=5)],["Selected fit","Stopping criterion met","Stopping criterion not met"];orientation=:horizontal,labelsize=10,tellwidth=false)
export_figure(joinpath(dir,"figures"),"exposure_bic",fig)

actual=Figure(size=(900,380))
bicaxis=Axis(actual[1,1],xlabel="BIC",yticks=(1:length(rank),[labels[r.scenario] for r in rank]),xlabelsize=15,yticklabelsize=13,xticklabelsize=12)
for (i,r) in enumerate(rank)
    color=startswith(r.scenario,"lag_") ? "#245B83" : "#AA4E35"
    group=filter(x->x.scenario==r.scenario,runs)
    values=Float64[g.bic for g in group]
    isempty(values)&&error("No saved BIC values for $(r.scenario)")
    all(isfinite,values)||error("Nonfinite saved BIC values")
    if length(values)>1
        boxplot!(bicaxis,fill(i,length(values)),values;orientation=:horizontal,color=(color,.3),show_outliers=false,width=.5)
    end
    scatter!(bicaxis,values,fill(i,length(values));color=(color,.6),marker=:circle,markersize=2.5)
    scatter!(bicaxis,[r.bic],[i];color=color,marker=:diamond,markersize=6)
end
ylims!(bicaxis,.5,length(rank)+.5);bicaxis.yreversed=true
export_figure(joinpath(dir,"figures"),"exposure_bic_actual",actual)
write(joinpath(dir,"figures","exposure_bic_actual_caption.txt"),"Actual training BIC across saved optimization runs for each exposure summary. Diamonds mark selected fits; points show all saved runs. Rows retain the original ranking. The companion delta-BIC figure is preserved. Convergence flags remain in ranking.csv and runs.csv; marker shape in this figure does not encode convergence. No refitting or data filtering is performed.\n")

f=Figure(size=(700,500));a=Axis(f[1,1],xticks=(1:4,["Cells","Glucose","Lactate","Aggregate size"]),yticks=(1:length(rank),[labels[r.scenario] for r in rank]),xticklabelsize=11,yticklabelsize=10)
states=("cells","glucose","lactate","aggregate")
z=[only(filter(e->e.scenario==r.scenario && e.state==s,errors)).mrpe_percent for s in states,r in rank]
heatmap!(a,1:4,1:length(rank),z;colormap=:Blues)
for i in 1:4,j in 1:length(rank)
    text!(a,i,j;text=string(round(z[i,j];digits=1)),align=(:center,:center),fontsize=11,color=z[i,j]>.6maximum(z) ? :white : :black)
end
a.yreversed=true
Label(f[2,1],"Training MRPE (%) against condition–day means";fontsize=11,tellwidth=false)
export_figure(joinpath(dir,"figures"),"exposure_state_errors",f)
rows=join(["<tr><td>$(labels[r.scenario])</td><td>$(round(r.bic;digits=3))</td><td>$(round(r.delta_bic;digits=3))</td><td>$(r.stopping_met)</td><td>$(r.converged_starts)/$(r.total_starts)</td><td>$(r.boundaries)</td></tr>" for r in rank])
write(joinpath(dir,"index.html"),"""<!doctype html><meta charset="utf-8"><title>Hydrodynamic exposure sensitivity</title><style>body{font:16px system-ui;max-width:1100px;margin:36px auto;line-height:1.5}img{max-width:100%}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}</style><h1>Hydrodynamic exposure-summary sensitivity</h1><p>Model 6; five 100-mL conditions; 116 condition–day–state means; 11 estimated parameters per summary. $(settings["starts_per_summary"]) independent hybrid searches per summary. Training-only comparison; no 500-mL outcome selection. This is conditional on Model 6, not a new comparison of all model structures.</p><img src="figures/exposure_bic.svg"><h2>Actual BIC</h2><img src="figures/exposure_bic_actual.svg"><table><tr><th>Summary</th><th>BIC</th><th>ΔBIC</th><th>Selected stopping met</th><th>Converged starts</th><th>Boundaries</th></tr>$rows</table><img src="figures/exposure_state_errors.svg"><p>Boxplots describe optimizer variability, not biological uncertainty. All finite run scores are retained. Convergence does not establish identifiability. BIC uses the existing working diagonal likelihood; repeated-day correlation remains unmodelled.</p><p>Lagrangian summaries use predicted aggregate size with linear interpolation and endpoint clamping. Eulerian values are constant for a vessel/RPM. Shear and EDR summaries are marginal; maxima need not occur simultaneously. Source weighting definitions are inherited from the supplied CFD workbooks.</p><p><a href="ranking.csv">Ranking</a> · <a href="state_errors.csv">State errors and experimental CV</a> · <a href="runs.csv">All starts</a> · <a href="exposure_inputs.csv">CFD inputs and workbook cells</a></p>""")
