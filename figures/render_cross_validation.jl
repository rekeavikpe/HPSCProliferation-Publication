using CairoMakie, CSV, Statistics
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1||error("Usage: render_cross_validation.jl RESULT_DIRECTORY")
root=abspath(ARGS[1]);isfile(joinpath(root,"COMPLETE"))||error("Results incomplete")
rows=collect(CSV.File(joinpath(root,"state_metrics.csv")))
rows=filter(r->r.cohort=="prediction",rows)
out=joinpath(root,"figures");manuscript_theme!()
styles=Dict(0=>(:gray,:circle),6=>(SCALE_COLOR,:rect))
states=("cells","glucose","lactate","aggregate","all")
labels=("Cells","Glucose","Lactate","Aggregate size","Overall")
for exp in ("loco","partial500")
    rr=filter(r->r.experiment==exp,rows);isempty(rr)&&continue
    folds=sort(unique(String(r.fold) for r in rr))
    f=Figure(size=(1080,630),fontsize=12)
    for (i,state) in enumerate(states)
        a=Axis(f[(i-1)÷3+1,(i-1)%3+1],title="($(Char(96+i))) $(labels[i])",titlealign=:left,titlesize=15,
            xlabel=exp=="loco" ? "Omitted 100 mL condition" : "500 mL culture split",ylabel="Prediction MRPE (%)",
            xlabelsize=13,ylabelsize=13,xticklabelsize=11,yticklabelsize=11,xticks=(1:length(folds),string.(1:length(folds))))
        for model in (0,6)
            variants=exp=="partial500" ? ("partial500_fit","100mL_only_reference") : ("fold_fit",)
            for variant in variants
                vals=[only(filter(r->r.fold==fold&&r.state==state&&r.publication_model==model&&r.variant==variant,rr)).mrpe_percent for fold in folds]
                color,marker=styles[model]
                scatterlines!(a,1:length(folds),vals;color,marker,markersize=7,linewidth=1.8,linestyle=variant=="100mL_only_reference" ? :dash : :solid,
                    label=exp=="partial500" ? "Model $model — $(variant=="partial500_fit" ? "with 4 larger-vessel cultures" : "100 mL only")" : "Model $model")
            end
        end
        noise=[first(filter(r->r.fold==fold&&r.state==state,rr)).experimental_noise_percent for fold in folds]
        scatterlines!(a,1:length(folds),noise;color=:black,marker=:diamond,markersize=5,linestyle=:dot,linewidth=1,label="Experimental noise")
    end
    Legend(f[3,1:3],f.content[1];orientation=:horizontal,nbanks=exp=="partial500" ? 2 : 1,labelsize=12)
    rowgap!(f.layout,25);colgap!(f.layout,30)
    export_figure(out,exp=="loco" ? "leave_condition_out" : "partial500_comparison",f)
end
write(joinpath(out,"captions.txt"),"Prediction MRPE by state and overall. Each point is a complete held-out condition or a held-out set of five 500 mL cultures. Partial-500 fits use all 100 mL cultures plus four 500 mL cultures. Dashed reference curves use a separate fit to 100 mL data only and are evaluated on exactly the same held-out cultures. Experimental noise is the mean sample SD/mean across eligible held-out measurements. Splits overlap and are not independent experiments. Lines connect discrete folds for readability. These are retrospective robustness checks of prespecified Models 0 and 6, not new independent validation.\n")
println("Saved cross-validation figures: $out")
