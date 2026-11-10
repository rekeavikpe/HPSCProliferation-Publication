using CSV, CairoMakie, Statistics, TOML, Random
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==2 || error("Usage: render_training60.jl OUTPUT_ROOT main|fits|cells")
root=abspath(ARGS[1]);out=joinpath(root,"figures");source=joinpath(root,"reference_fit/figure_data");mkpath(out);manuscript_theme!()
function save_response_figure(out,stem,fig)
 if get(ENV,"HPSC_PREVIEW_ONLY","0")=="1"
  save(joinpath(out,stem*".png"),fig;px_per_unit=1)
 else
  export_figure(out,stem,fig)
 end
end
if ARGS[2] in ("main","fits")
 if ARGS[2]=="main"
 fit=joinpath(root,"reference_fit")
 rows=sort(collect(CSV.File(joinpath(fit,"selection_scores.csv")));by=r->r.model)
 errors=collect(CSV.File(joinpath(fit,"selection_state_errors.csv")))
 runs=collect(CSV.File(joinpath(fit,"selection_optimizer_bic.csv")))
 f=Figure(size=(1150,520));states5=("cells","glucose","lactate","aggregate","all")
 ax=Axis(f[1,1],title="(a) Model selection",xlabel="Model",ylabel="BIC (40% selection data)",xticks=0:6,titlesize=23,xlabelsize=20,ylabelsize=20,xticklabelsize=17,yticklabelsize=17)
 rng=MersenneTwister(20261016)
 for r in rows
  values=[v.bic for v in runs if v.model==r.model]
  boxplot!(ax,fill(r.model,length(values)),values;show_outliers=false,color=(:steelblue,.45))
  scatter!(ax,r.model .+ .18 .*(rand(rng,length(values)) .- .5),values;markersize=3,color=(:steelblue,.5))
  scatter!(ax,[r.model],[r.bic];marker=:diamond,markersize=10,color=:black)
 end
 ax2=Axis(f[1,2],title="(b) Selection MRPE",xticks=(1:5,["Cells","Glucose","Lactate","Aggregate size","Overall"]),yticks=(1:7,["Model $i" for i in 0:6]),yreversed=true,xticklabelrotation=pi/4,titlesize=23,xticklabelsize=16,yticklabelsize=17)
 z=[only(filter(r->r.model==id&&r.state==st,errors)).mrpe_percent for st in states5,id in 0:6]
 h=heatmap!(ax2,1:5,1:7,z;colormap=:Blues,colorrange=(0,maximum(z)))
 for i in 1:5,j in 1:7
  text!(ax2,i,j;text=string(round(z[i,j];digits=1)),align=(:center,:center),fontsize=17,color=z[i,j]>.65maximum(z) ? :white : :black)
 end
 Colorbar(f[1,3],h;label="MRPE (%)",labelsize=18,ticklabelsize=16)
 colgap!(f.layout,35);export_figure(out,"model_selection_40",f)
 write(joinpath(out,"model_selection_40_caption.txt"),"BIC and MRPE evaluated on the 40% selection cultures using parameters fitted to the 60% training cultures. Boxes show repeated training optimization runs. Diamonds show the run selected by training likelihood for each model.\n")

 er=collect(CSV.File(joinpath(root,"exposure_transfer/selection_scores.csv")))
 em=collect(CSV.File(joinpath(root,"exposure_transfer/state_errors.csv")))
 names=replace.([String(r.scenario) for r in er],"lag_"=>"Lagrangian ","eul_"=>"Eulerian ","p"=>"P","va"=>"volume average","max"=>"maximum")
 ef=Figure(size=(1250,650))
 ea=Axis(ef[1,1],title="(a) Exposure selection",xlabel="BIC (40% selection data)",yticks=(1:length(er),names),yreversed=true,titlesize=21,xlabelsize=18,xticklabelsize=15,yticklabelsize=15)
 scatter!(ea,[r.bic for r in er],1:length(er);marker=:diamond,markersize=10,color="#005A70")
 for (col,cohort,title) in ((2,"selection40","(b) Selection MRPE"),(3,"500 mL prediction","(c) 500 mL prediction"))
  ex=Axis(ef[1,col],title=title,xticks=(1:5,["Cells","Glucose","Lactate","Aggregate size","Overall"]),yticks=(1:length(er),names),yticklabelsvisible=false,yreversed=true,xticklabelrotation=pi/4,titlesize=20,xticklabelsize=14)
  zz=[only(filter(v->v.scenario==r.scenario&&v.cohort==cohort&&v.state==st,em)).mrpe_percent for st in states5,r in er]
  hh=heatmap!(ex,1:5,1:length(er),zz;colormap=:Blues,colorrange=(0,maximum(v.mrpe_percent for v in em)))
  for i in 1:5,j in 1:length(er)
   text!(ex,i,j;text=string(round(zz[i,j];digits=1)),align=(:center,:center),fontsize=14,color=zz[i,j]>.65maximum(v.mrpe_percent for v in em) ? :white : :black)
  end
  col==3 && Colorbar(ef[1,4],hh;label="MRPE (%)",labelsize=18,ticklabelsize=15)
 end
 colgap!(ef.layout,25);export_figure(out,"exposure_selection_and_transfer",ef)
 end
states=("cells","glucose","lactate","aggregate")
statelabels=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
points=collect(CSV.File(joinpath(source,"training_points.csv")));curves=collect(CSV.File(joinpath(source,"training_curves.csv")))
testcurves=collect(CSV.File(joinpath(source,"test_trajectories.csv")));testpoints=collect(CSV.File(joinpath(source,"test_predictions.csv")))
@assert TOML.parsefile(joinpath(source,"source.toml"))["selected_model"]==5
conditions=(1,2,3,5,6);divide(s)=s=="cells" ? 1e6 : 1.
for nrows in (5,6)
 fig=Figure(size=(1100,nrows==6 ? 1220 : 1030),fontsize=13)
 for i in 1:nrows,j in 1:4
  state=states[j];d=divide(state);validation=i==6
  rr=validation ? filter(r->r.publication_model==5&&r.state==state,testcurves) : filter(r->r.model==5&&r.condition==conditions[i]&&r.state==state,curves)
  pp=validation ? filter(r->r.publication_model==5&&r.state==state,testpoints) : filter(r->r.condition==conditions[i]&&r.state==state,points)
  color=validation ? SCALE_COLOR : COLORS[conditions[i]]
  local ax=Axis(fig[i,j+1],title="($(Char(96+(i-1)*4+j)))",titlealign=:left,titlesize=14,xlabel=i==nrows ? "Time (days)" : "",ylabel=statelabels[j],xlabelsize=13,ylabelsize=12,xticklabelsize=11,yticklabelsize=11,xticks=0:6)
  lines!(ax,[r.day for r in rr],[r.value/d for r in rr];color,linewidth=1.8)
  yy=[(validation ? r.observed : r.mean)/d for r in pp]
  errorbars!(ax,[r.day for r in pp],yy,[r.sd/d for r in pp];color,linewidth=1.4,whiskerwidth=7)
  scatter!(ax,[r.day for r in pp],yy;color,markersize=6)
  xlims!(ax,-.15,maximum(r.day for r in rr)+.15)
  ylims!(ax,0,(2.2,18.,20.,380.)[j])
 end
 for i in 1:nrows
  Label(fig[i,1],i==6 ? "500 mL\nValidation" : "100 mL\nCondition $i";fontsize=14,font=:bold,rotation=0,tellheight=false)
 end
 Legend(fig[nrows+1,2:5],[LineElement(color=:black,linewidth=1.8),MarkerElement(color=:black,marker=:circle,markersize=6)],["Model 6","Experimental data ± SD"];orientation=:horizontal,tellwidth=false,labelsize=14)
 colgap!(fig.layout,6);rowgap!(fig.layout,16)
 save_response_figure(out,nrows==6 ? "model6_six_by_four" : "model6_five_by_four",fig)
end

elseif ARGS[2]=="cells"
const STATES4=("cells","glucose","lactate","aggregate")
const LABELS=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
rows=collect(CSV.File(joinpath(root,"cell_lines/curves.csv")));obs=collect(CSV.File(joinpath(root,"cell_lines/observations.csv")))
# Two complementary figures avoid seven indistinguishable curves in one panel.
for (stem,strategies,labels) in (("cell_lines_biological_vs_metabolic",("transfer_model6","biological7","metabolic5"),("Model 6 transfer","Seven parameters refitted","Five metabolic parameters refitted")),("cell_lines_response_comparison",("transfer_model0","baseline_refit","transfer_model6","hydrodynamic2","all11"),("Model 0 transfer","Model 0 refitted","Model 6 transfer","Hydrodynamic parameters refitted","All Model 6 parameters refitted")))
 fig=Figure(size=(1100,560),fontsize=12);styles=(:solid,:dash,:dot,:dashdot,:dashdotdot);legend=Any[]
 for (i,line) in enumerate(("585A1","SCTi003-A")),(j,state) in enumerate(STATES4)
  color=i==1 ? RGBf(.165,.290,.514) : RGBf(.486,.196,.165);scale=state=="cells" ? 1e6 : 1.
  local ax=Axis(fig[i,j],title="($(Char('a'+(i-1)*4+j-1))) $line",xlabel="Time (days)",ylabel=LABELS[j],xticks=0:maximum(r.day for r in obs if r.cell_line==line),titlesize=13,xlabelsize=12,ylabelsize=12,xticklabelsize=10,yticklabelsize=10)
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

else
 error("Unknown plotting mode")
end
