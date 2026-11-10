using CSV, CairoMakie
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1||error("Usage: render_exposure_transfer.jl EVALUATION")
dir=abspath(only(ARGS));isfile(joinpath(dir,"COMPLETE"))||error("Evaluation incomplete")
ranking=collect(CSV.File(joinpath(dir,"training_ranking.csv")))
metrics=collect(CSV.File(joinpath(dir,"state_errors.csv")))
labels=Dict("lag_p50"=>"Lagrangian P50","lag_mean"=>"Lagrangian mean","lag_p25"=>"Lagrangian P25","lag_p75"=>"Lagrangian P75","lag_p90"=>"Lagrangian P90","lag_max"=>"Lagrangian maximum","eul_va"=>"Eulerian volume average","eul_p50"=>"Eulerian P50","eul_p75"=>"Eulerian P75","eul_p90"=>"Eulerian P90","eul_max"=>"Eulerian maximum")
manuscript_theme!()
states=("cells","glucose","lactate","aggregate","all")
for (metric,stem,colorlabel) in ((:mrpe_percent,"exposure_transfer_mrpe","MRPE (%)"),(:experimental_noise_percent,"exposure_transfer_noise","Experimental Noise (%)"))
    matrices=[[Float64(getproperty(only(filter(r->r.scenario==rank.scenario&&r.cohort==cohort&&r.state==state,metrics)),metric)) for state in states,rank in ranking] for cohort in ("100 mL calibration","500 mL prediction")]
    upper=maximum(maximum(z) for z in matrices)
    fig=Figure(size=(1450,650))
    for (j,cohort) in enumerate(("100 mL calibration","500 mL prediction"))
        ax=Axis(fig[1,j],title="($(j==1 ? "a" : "b")) $cohort",titlesize=21,
            xticks=(1:5,["Cells","Glucose","Lactate","Aggregate size","Overall"]),xticklabelrotation=π/5,xticklabelsize=16,
            yticks=(1:length(ranking),[labels[r.scenario] for r in ranking]),yticklabelsize=15,yreversed=true)
        h=heatmap!(ax,1:5,1:length(ranking),matrices[j];colorrange=(0.,upper),colormap=:Blues)
        for x in 1:5,y in 1:length(ranking)
            v=matrices[j][x,y];text!(ax,x,y;text=string(round(v;digits=1)),align=(:center,:center),fontsize=17,color=v>.6upper ? :white : :black)
        end
        j==2&&Colorbar(fig[1,3],h;label=colorlabel,labelsize=18,ticklabelsize=15)
    end
    colgap!(fig.layout,28)
    export_figure(joinpath(dir,"figures"),stem,fig)
end
write(joinpath(dir,"index.html"),"""<!doctype html><meta charset="utf-8"><title>Exposure representation: scale transfer</title><style>body{font:16px Arial;max-width:1500px;margin:30px auto;line-height:1.5}img{width:100%}</style><h1>Exposure representation: saved-fit evaluation at 500 mL</h1><p>All eleven fits retain their 100 mL parameter estimates and use the matching 500 mL CFD summary. Rows retain the original training-BIC order. This is a retrospective prediction evaluation, with no optimization or 500 mL model selection.</p><img src="figures/exposure_transfer_mrpe.svg"><img src="figures/exposure_transfer_noise.svg"><p><a href="state_errors.csv">Per-state and overall metrics</a> · <a href="verification.csv">Training likelihood verification and convergence flags</a> · <a href="predictions.csv">Predictions</a> · <a href="metric_definitions.txt">Metric definitions</a></p><p>Experimental Noise is the mean sample SD relative to the experimental mean, expressed as a percentage. Overall pools eligible condition–day–state entries. Noise is identical across exposure representations because the experimental data are unchanged. Missing upper or lower identifiability limits and unsuccessful optimizer stopping criteria are not resolved by this evaluation.</p>""")
