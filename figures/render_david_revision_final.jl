using CSV, CairoMakie, Statistics
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1||error("Usage: render_david_revision.jl RESULTS")
root=abspath(only(ARGS));out=joinpath(root,"reporting_corrected","figures");mkpath(out)
manuscript_theme!()
const STATES4=("cells","glucose","lactate","aggregate")
const LABELS=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
const DISPLAY=("Cells","Glucose","Lactate","Aggregate size","Overall")
function error_heatmap(stage,stem)
 rows=collect(CSV.File(joinpath(root,stage,"state_errors.csv")));rank=collect(CSV.File(joinpath(root,stage,"ranking.csv")))
 ids=[string(r.scenario) for r in rank];names=stage=="selection" ? ["Model $(r.publication_model)" for r in rank] : replace.(ids,"lag_"=>"Lagrangian ","eul_"=>"Eulerian ","p"=>"P","va"=>"volume average","max"=>"maximum")
 fig=Figure(size=(850,stage=="selection" ? 350 : 480),fontsize=12)
 for (col,cohort,title) in ((1,"calibration60","(a) Calibration: 60% of replicates"),(2,"selection40","(b) Selection: withheld 40%"))
  vals=[only(filter(r->string(r.scenario)==id&&r.cohort==cohort&&r.state==s,rows)).mrpe_percent for s in (STATES4...,"all"),id in ids]
  ax=Axis(fig[1,col],title=title,xticks=(1:5,collect(DISPLAY)),yticks=(1:length(ids),names),yreversed=true,xticklabelrotation=pi/4,xticklabelsize=11,yticklabelsize=11,titlesize=13)
  hm=heatmap!(ax,1:5,1:length(ids),vals;colormap=:Blues,colorrange=(0,max(1.,maximum(r.mrpe_percent for r in rows))))
  thresholds=[sort(collect(vals[i,:]))[min(3,size(vals,2))] for i in axes(vals,1)]
  for i in axes(vals,1),j in axes(vals,2);text!(ax,i,j;text=string(round(vals[i,j];digits=1)),align=(:center,:center),fontsize=11,font=vals[i,j]<=thresholds[i] ? :bold : :regular,color=vals[i,j]>.65maximum(r.mrpe_percent for r in rows) ? :white : :black);end
  col==2&&Colorbar(fig[1,3],hm;label="MRPE (%)",labelsize=12,ticklabelsize=10)
 end
 export_figure(out,stem,fig)
end
error_heatmap("selection","heldout_model_errors")
error_heatmap("exposure_selection","heldout_exposure_errors")
rows=collect(CSV.File(joinpath(root,"cell_lines/curves.csv")));obs=collect(CSV.File(joinpath(root,"cell_lines/observations.csv")))
# Two complementary figures avoid seven indistinguishable curves in one panel.
for (stem,strategies,labels) in (("cell_lines_biological_vs_metabolic",("transfer_model6","biological7","metabolic5"),("Model 6 transfer","Seven parameters refitted","Five metabolic parameters refitted")),("cell_lines_response_comparison",("transfer_model0","baseline_refit","transfer_model6","hydrodynamic2","all11"),("Model 0 transfer","Model 0 refitted","Model 6 transfer","Hydrodynamic parameters refitted","All Model 6 parameters refitted")))
 fig=Figure(size=(1100,560),fontsize=12);styles=(:solid,:dash,:dot,:dashdot,:dashdotdot);legend=Any[]
 for (i,line) in enumerate(("585A1","SCTi003-A")),(j,state) in enumerate(STATES4)
  color=i==1 ? RGBf(.165,.290,.514) : RGBf(.486,.196,.165);scale=state=="cells" ? 1e6 : 1.
  ax=Axis(fig[i,j],title="($(Char('a'+(i-1)*4+j-1))) $line",xlabel="Time (days)",ylabel=LABELS[j],xticks=0:maximum(r.day for r in obs if r.cell_line==line),titlesize=13,xlabelsize=12,ylabelsize=12,xticklabelsize=10,yticklabelsize=10)
  commonobs=filter(r->r.state==state,obs)
  commoncurves=filter(r->r.state==state&&r.strategy in strategies,rows)
  ymax=max(maximum((r.mean+r.sd)/scale for r in commonobs),maximum(r.value/scale for r in commoncurves));ylims!(ax,0,1.08ymax)
  for (k,strategy) in enumerate(strategies)
   rs=filter(r->r.cell_line==line&&r.state==state&&r.strategy==strategy,rows)
   l=lines!(ax,[r.day for r in rs],[r.value/scale for r in rs];color,linestyle=styles[k],linewidth=1.6)
   i==1&&j==1&&push!(legend,LineElement(color=:black,linestyle=styles[k],linewidth=1.6))
  end
  rs=filter(r->r.cell_line==line&&r.state==state,obs)
  errorbars!(ax,[r.day for r in rs],[r.mean/scale for r in rs],[r.sd/scale for r in rs];color,linewidth=1.3,whiskerwidth=6)
  scatter!(ax,[r.day for r in rs],[r.mean/scale for r in rs];color,markersize=5)
 end
 push!(legend,MarkerElement(color=:black,marker=:circle,markersize=6))
 Legend(fig[3,1:4],legend,[labels...,"Experimental data ± SD"];orientation=:horizontal,nbanks=2,labelsize=11)
 export_figure(out,stem,fig)
end
write(joinpath(root,"reporting_corrected","figures_report.html"),"""<!doctype html><meta charset="utf-8"><title>David revision analyses</title><style>body{font:16px Arial;max-width:1400px;margin:35px auto;line-height:1.5}img{width:100%}a{color:#005a70}</style><h1>Revision analyses</h1><p>The tables contain full precision. Figures show MRPE to one decimal place. Bold values identify the three lowest errors in each state column within each panel, including exact ties. Model and exposure comparisons use errors on the withheld 40%. Cell-line refits are assessed on their calibration data.</p><h2>Model selection</h2><img src="figures/heldout_model_errors.svg"><p><a href="../selection/ranking.csv">Calibration BIC and held-out likelihood</a> · <a href="../selection/state_errors.csv">Errors and experimental noise</a></p><h2>Exposure summaries</h2><img src="figures/heldout_exposure_errors.svg"><p><a href="../exposure_selection/ranking.csv">Ranking</a> · <a href="../exposure_selection/state_errors.csv">Errors</a></p><h2>Cell lines</h2><img src="figures/cell_lines_biological_vs_metabolic.svg"><img src="figures/cell_lines_response_comparison.svg"><p><a href="../cell_lines/state_errors.csv">MRPE and experimental noise</a> · <a href="../cell_lines/parameter_policy.csv">Estimated and fixed parameters</a> · <a href="../cell_lines/interpretation.txt">Analysis definitions</a></p><h2>Cross-validation and exclusion sensitivity</h2><p><a href="../pooled_validation/replicate_error_summary.csv">Mean, SD and variance of replicate errors</a> · <a href="../pooled_validation/fold_error_summary.csv">Fold-level summaries</a></p><h2>Shear and EDR</h2><p><a href="../separate_ablations/index.html">Eight formulations, with separate response parameters for both-term variants</a>. All-active and fixed-three profiles are under separate_ablations/profiles. Both profile sets are saved. Use the corrected interval table for reporting.</p><h2>Experimental comparisons</h2><p><a href="../rpm_statistics/welch_holm.csv">Welch tests, Holm-adjusted p-values and protocol checks</a> · <a href="../rpm_statistics/culture_endpoints.csv">Culture-level observations</a></p><p>Oxygen CFD is pending from Pouyan and is not included in these calculations.</p>""")
println("Report and publication figures saved: $root/reporting_corrected/figures_report.html")
