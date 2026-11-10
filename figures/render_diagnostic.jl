using CairoMakie, CSV, Ghostscript_jll
include("../src/PublicationModels.jl")
using .PublicationModels
CairoMakie.activate!()
root=normpath(joinpath(@__DIR__,".."))
source=isempty(ARGS) ? joinpath(root,"reference/scale_transfer") : abspath(ARGS[1])
out=joinpath(@__DIR__,"output");mkpath(out)
set_theme!(Theme(font="TeX Gyre Heros Makie",fontsize=9,Axis=(xlabelsize=10,ylabelsize=10,titlesize=10,titlefont=:bold,xticklabelsize=9,yticklabelsize=9,spinewidth=.7,xgridvisible=false,ygridvisible=false),Lines=(linewidth=1.5,),Scatter=(markersize=4,)))
function export_figure(stem,fig;dpi=600)
    # Render off the synced Desktop, then replace complete output files.
    mktempdir() do staging
        pdf=joinpath(staging,stem*".pdf")
        save(pdf,fig;pt_per_unit=1)
        save(joinpath(staging,stem*".svg"),fig;pt_per_unit=1)
        gs=Ghostscript_jll.gs()
        for (extension,device) in (("png","png16m"),("jpeg","jpeg"),("tiff","tiff24nc"),("eps","eps2write"))
            dest=joinpath(staging,stem*"."*extension)
            run(`$gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=$device -r$dpi -dJPEGQ=98 -sOutputFile=$dest $pdf`)
        end
        for file in readdir(staging;join=true)
            mv(file,joinpath(out,basename(file));force=true)
        end
    end
end

states=("cells","glucose","lactate","aggregate")
labels=("Viable cells (10⁵ cells/mL)","Glucose (mM)","Lactate (mM)","Aggregate diameter (µm)")
trajectories=collect(CSV.File(joinpath(source,"test_trajectories.csv")))
points=collect(CSV.File(joinpath(source,"test_predictions.csv")))
candidates=sort(unique(r.publication_model for r in points if r.publication_model!=0))
length(candidates)<=1 || error("Expected at most one selected hydrodynamic model plus Model 0")
selected=isempty(candidates) ? 0 : only(candidates)
models=selected==0 ? ((0,"#BD5E2F",:dash),) : ((0,"#BD5E2F",:dash),(selected,"#176B87",:solid))
fig=Figure(size=(510,380))
for (k,state) in enumerate(states)
    row,col=divrem(k-1,2); scale=state=="cells" ? 1e5 : 1.
    ax=Axis(fig[row+1,col+1],xlabel="Culture day",ylabel=labels[k],xticks=0:5,title=string(Char('A'+k-1)),titlealign=:left)
    for (id,color,style) in models
        rr=filter(r->r.publication_model==id&&r.state==state,trajectories)
        lines!(ax,[r.day for r in rr],[r.value/scale for r in rr];color,linestyle=style)
    end
    pp=filter(r->r.publication_model==selected&&r.state==state,points)
    errorbars!(ax,[r.day for r in pp],[r.observed/scale for r in pp],[r.sd/scale for r in pp];color=:black,whiskerwidth=4,linewidth=.8)
    scatter!(ax,[r.day for r in pp],[r.observed/scale for r in pp];color=:black)
    xlims!(ax,-.1,5.15)
end
elements=Any[LineElement(color=color,linestyle=style,linewidth=1.5) for (id,color,style) in models]
labels_legend=[id==0 ? "No-hydrodynamics reference" : "Hydrodynamic model $(publication_model(id))" for (id,_,_) in models]
push!(elements,MarkerElement(color=:black,marker=:circle,markersize=5));push!(labels_legend,"Measured mean ± SD")
Legend(fig[3,1:2],elements,labels_legend,orientation=:horizontal,framevisible=false,labelsize=8)
export_figure("scale_transfer",fig)
write(joinpath(out,"scale_transfer_caption.txt"),"Retrospective prediction of the 500-mL condition after calibration on 100-mL data. Points show daily condition means; error bars show between-bioreactor SD, not prediction intervals. Lines use saved training-fit parameters without refitting to 500-mL outcomes. Vertical metabolite changes represent medium exchange; observations are pre-exchange. Panel labels: A, viable cells; B, glucose; C, lactate; D, aggregate diameter. These are reproducibility-run results, not the completed R2 analysis.\nAlt text: Four panels compare measured 500-mL culture outcomes with the displayed model predictions across five days. Glucose rises and lactate falls at media exchanges.\n")
println("Saved PDF, EPS, SVG, 600-dpi PNG/TIFF and JPEG convenience copy to $out")
rankfile=isfile(joinpath(source,"training_ranking.csv")) ? joinpath(source,"training_ranking.csv") : joinpath(dirname(source),"ranking.csv")
if isfile(rankfile)
    ranked=collect(CSV.File(rankfile))
    modelid(r)=publication_model(hasproperty(r,:model) ? Int(r.model) : Int(r.publication_model))
    sort!(ranked;by=modelid)
    best=minimum(r.bic for r in ranked)
    fig2=Figure(size=(510,230))
    ax=Axis(fig2[1,1],xlabel="Model",ylabel="ΔBIC",xticks=0:6)
    colors=[r.bic==best ? "#176B87" : "#A6B7C2" for r in ranked]
    barplot!(ax,modelid.(ranked),[r.bic-best for r in ranked];color=colors)
    for r in ranked
        !r.stopping_met && scatter!(ax,[modelid(r)],[r.bic-best+3];marker=:xcross,color="#A54532",markersize=7)
    end
    scatter!(ax,[modelid(only(filter(r->r.bic==best,ranked)))],[0.];color="#176B87",markersize=5)
    Legend(fig2[2,1],[MarkerElement(color="#A54532",marker=:xcross,markersize=7)],["Selected start did not meet stopping criterion"],framevisible=false,labelsize=8,tellwidth=false,tellheight=true)
    export_figure("training_selection",fig2)
    write(joinpath(out,"training_selection_caption.txt"),"Model comparison using the 100-mL working Gaussian likelihood. Bars show BIC relative to the lowest achieved BIC. Crosses mark models whose selected optimizer start did not meet its stopping criterion. These are model scores, not biological replicates or uncertainty intervals.\nAlt text: Bars compare seven model scores relative to the best achieved training BIC. Crosses identify unconverged selected fits.\n")
end
