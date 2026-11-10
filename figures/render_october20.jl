using CSV,CairoMakie,Statistics
include("ManuscriptStyle.jl");using .ManuscriptStyle
root=only(ARGS);out=get(ENV,"HPSC_OCTOBER_FIGURE_OUT",joinpath(root,"figures"));mkpath(out);manuscript_theme!()
function writefig(stem,fig)
    save(joinpath(out,stem*".png"),fig;px_per_unit=3)
    save(joinpath(out,stem*".pdf"),fig)
end
candidate=collect(CSV.File(joinpath(root,"reuse/candidate_selection_optimizer_bic.csv")))
models=[0,1,3,4,5,6,7]
f=Figure(size=(1050,570),fontsize=18);ax=Axis(f[1,1],ylabel="BIC",xlabel="Model",title="Selection data (40%)",xticks=(1:7,string.(models)),xlabelsize=18,ylabelsize=18,xticklabelsize=16,yticklabelsize=16)
for (i,m) in enumerate(models)
    vals=Float64[r.bic for r in candidate if r.model==m]
    boxplot!(ax,fill(i,length(vals)),vals;color=(:steelblue,.55),show_outliers=false)
    scatter!(ax,[i+.11sin(j*2.4) for j in eachindex(vals)],vals;color=(:steelblue,.55),markersize=4)
end
writefig("candidate_selection_bic_original_numbers",f)
scenarios=["baseline","shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both","full"]
labels=["Model 0","Shear: cells","Shear: aggregate","Shear: both","EDR: cells","EDR: aggregate","EDR: both","Model 6"]
rows=collect(CSV.File(joinpath(root,"reuse/selection_optimizer_bic.csv")))
f=Figure(size=(1100,600),fontsize=18);ax=Axis(f[1,1],ylabel="BIC",xlabel="Formulation",title="Selection data (40%)",xticks=(1:8,labels),xticklabelrotation=pi/5)
for (i,s) in enumerate(scenarios)
    vals=Float64[r.bic for r in rows if r.scenario==s];isempty(vals)&&error("Missing optimizer distribution: $s")
    boxplot!(ax,fill(i,length(vals)),vals;color=(:steelblue,.55),show_outliers=false)
    scatter!(ax,[i+.11sin(j*2.4) for j in eachindex(vals)],vals;color=(:steelblue,.55),markersize=4)
end
writefig("ablation_selection_bic",f)
# Each density defines its own Da field. Overlay its level curves on shared axes.
oxy=sort(collect(CSV.File(joinpath(root,"oxygen/condition_summary.csv")));by=r->r.condition)
f=Figure(size=(1100,700),fontsize=17);ax=Axis(f[1,1],xlabel="kLa (h⁻¹)",ylabel="Oxygen uptake (amol cell⁻¹ s⁻¹)",title="Bulk oxygen demand / transfer capacity",limits=((.44,9.),(8.85,35.4)),xticks=[.5,1.,2.,3.,4.,5.,6.,7.,8.,9.])
colors=Makie.wong_colors();els=Any[];names=String[]
for (i,r) in enumerate(oxy)
    xx=collect(range(.44,9.;length=300));peak=r.peak_cells_ml;color=colors[mod1(i,length(colors))]
    for (level,style) in ((.1,:dot),(.5,:dash),(1.,:solid))
        yy=level .* xx .* .00021 ./ (peak*1000*3600*1e-18)
        lines!(ax,xx,yy;color,linestyle=style,linewidth=1.8)
    end
    push!(els,LineElement(color=color,linewidth=2));push!(names,"Condition $(r.condition) (Da $(round(r.Da;digits=2)))")
end
for kla in unique(r.kla_h for r in oxy)
    ids=[string(r.condition) for r in oxy if r.kla_h==kla]
    scatter!(ax,[kla],[17.7];color=:black,markersize=10)
    text!(ax,kla,ids==["5"]||ids==["6"] ? 16.4 : 18.5;text=join(ids,", "),fontsize=13,align=(:center,:bottom))
end
for (level,style) in ((.1,:dot),(.5,:dash),(1.,:solid))
    push!(els,LineElement(color=:black,linestyle=style,linewidth=2));push!(names,"Contour: Da = $level")
end
Legend(f[2,1],els,names;nbanks=3,framevisible=true,labelsize=13,tellheight=true)
rowsize!(f.layout,1,Relative(.80))
writefig("oxygen_contours",f)
# CFD exposure comparison uses the saved summaries, without fitting eleven models.
inputs=collect(CSV.File(joinpath(root,"reuse/exposure_distributions.csv")))
sweep=collect(CSV.File(joinpath(root,"reuse/exposure_rpm_sweep.csv")))
f=Figure(size=(1100,700),fontsize=16)
for (j,state,label) in ((1,:shear_pa,"Shear stress (Pa)"),(2,:edr_m2_s3,"EDR (m²/s³)"))
    local ax=Axis(f[1,j],xlabel="Eulerian volume average",ylabel="Lagrangian P50 (200 µm)",title="($(Char('a'+j-1))) $label")
    for (v,rpm) in ((100,30),(100,40),(100,60),(500,30))
        a=only(r for r in inputs if r.scenario=="eul_va"&&r.volume_ml==v&&r.rpm==rpm)
        b=only(r for r in inputs if r.scenario=="lag_p50"&&r.size_um==200&&r.volume_ml==v&&r.rpm==rpm)
        scatter!(ax,[getproperty(a,state)],[getproperty(b,state)];markersize=12)
        text!(ax,getproperty(a,state),getproperty(b,state);text=" $v mL, $rpm RPM",fontsize=10,align=(:left,:bottom))
    end
    ax2=Axis(f[2,j],xlabel="Agitation (RPM)",ylabel=label,title="($(Char('c'+j-1))) 100 mL")
    for (scenario,name,color) in (("eul_va","Eulerian volume average",:darkorange),("lag_p50","Lagrangian P50 (200 µm)",:steelblue))
        rr=sort(filter(r->r.scenario==scenario,sweep);by=r->r.rpm)
        lines!(ax2,[r.rpm for r in rr],[getproperty(r,state) for r in rr];color,label=name,linewidth=2)
    end
    axislegend(ax2;position=:lt,labelsize=11)
end
writefig("exposure_input_comparison",f)
# Requested plots are produced from whichever reduced cell model passed its audit.
const STATES4=("cells","glucose","lactate","aggregate")
const YLABELS=("Viable cells (×10⁶ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate size (µm)")
# final_rounds.csv is written by the cell driver, including numerical-review status.
finals=collect(CSV.File(joinpath(root,"cells/final_rounds.csv")))
allcurves=Dict(String(r.line)=>collect(CSV.File(joinpath(r.final_fit,"curves.csv"))) for r in finals)
allobs=Dict(String(r.line)=>collect(CSV.File(joinpath(r.final_fit,"observations.csv"))) for r in finals)
f=Figure(size=(1500,720),fontsize=17)
for (i,line) in enumerate(("585A1","SCTi003-A")),(j,state) in enumerate(STATES4)
    scale=j==1 ? 1e6 : 1.;rs=allcurves[line];obs=filter(r->r.state==state,allobs[line])
    ymax=maximum(vcat([r.value/scale for rows in values(allcurves) for r in rows if r.state==state],[(r.mean+r.sd)/scale for rows in values(allobs) for r in rows if r.state==state]))
    local ax=Axis(f[i,j],title="($(Char('a'+(i-1)*4+j-1))) $line",xlabel="Time (days)",ylabel=YLABELS[j],xlabelsize=16,ylabelsize=16,xticklabelsize=14,yticklabelsize=14,titlesize=18);ylims!(ax,0,1.06ymax)
    for (strategy,color,style) in (("Unchanged ESI-017",:gray45,:dash),("Recalibrated",i==1 ? :steelblue : :darkorange,:solid))
        rr=filter(r->r.state==state&&r.strategy==strategy,rs)
        lines!(ax,[r.day for r in rr],[r.value/scale for r in rr];color,linestyle=style,linewidth=2.4)
    end
    color=i==1 ? :steelblue : :darkorange
    errorbars!(ax,[r.day for r in obs],[r.mean/scale for r in obs],[r.sd/scale for r in obs];color,linewidth=2,whiskerwidth=9)
    scatter!(ax,[r.day for r in obs],[r.mean/scale for r in obs];color,markersize=8)
end
Legend(f[3,1:4],[LineElement(color=:gray45,linestyle=:dash),LineElement(color=:steelblue),LineElement(color=:darkorange),MarkerElement(color=:black,marker=:circle)], ["ESI-017 parameter transfer","585A1 recalibration","SCTi003-A recalibration","Experimental data ± SD"];orientation=:horizontal,framevisible=true,labelsize=16)
writefig("cell_lines_reduced_landscape",f)
for row in finals
    pdir=joinpath(row.final_fit,"profiles");isfile(joinpath(pdir,"audit.csv"))||continue
    audit=collect(CSV.File(joinpath(pdir,"audit.csv")));fig=Figure(size=(1100,650),fontsize=16)
    for (i,a) in enumerate(audit)
        vals=sort(collect(CSV.File(joinpath(pdir,a.parameter,"profile.csv")));by=r->r.value)
        axis=Axis(fig[(i-1)÷3+1,(i-1)%3+1],xlabel=String(a.parameter),ylabel="Relative log-likelihood",xscale=log10,title="$(a.parameter): $(a.status)",titlesize=12)
        lines!(axis,[r.value for r in vals],[r.relative_loglik for r in vals];linewidth=2)
        hlines!(axis,[-1.920729410347062];color=:gray40,linestyle=:dash);ylims!(axis,-10,1)
    end
    writefig("profiles_"*String(row.line),fig)
end
