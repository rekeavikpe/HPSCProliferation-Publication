using CairoMakie, CSV, TOML, Statistics, Random
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==1 || error("Usage: render_split_selection.jl SPLIT_RESULT_DIRECTORY")
source=abspath(ARGS[1]);isfile(joinpath(source,"COMPLETE")) || error("Split comparison incomplete")
out=joinpath(source,"figures");manuscript_theme!()
protocol=TOML.parsefile(joinpath(source,"protocol.toml"))
rank=collect(CSV.File(joinpath(source,"selection_ranking.csv")))
ids=sort([Int(r.publication_model) for r in rank])
f=Figure(size=(850,420))
a=Axis(f[1,1],xlabel="Model",ylabel="Calibration BIC",xticks=ids,title="(a) Calibration partition",titlealign=:left,titlesize=12)
b=Axis(f[1,2],xlabel="Model",ylabel="ΔNLL on selection partition",xticks=ids,title="(b) Held-out predictive comparison",titlealign=:left,titlesize=12)
rng=MersenneTwister(20261003);rows=NamedTuple[];counts=NamedTuple[]
for pub in ids
    r=only(filter(r->r.publication_model==pub,rank));storage=Int(r.storage_model)
    initial=NamedTuple[]
    for run in 1:protocol["starts"]
        stages=collect(CSV.File(joinpath(source,"model$storage","run$run","stages.csv")))
        best=argmin(s->s.nll,stages)
        push!(initial,(publication_model=pub,run=run,bic=2best.nll+r.k*log(r.training_n),stopping_met=best.stopping_met))
    end
    append!(rows,initial);color=MODEL_COLORS[storage+1];v=[s.bic for s in initial]
    boxplot!(a,fill(pub,length(v)),v;show_outliers=false,color=(color,.4),width=.5)
    scatter!(a,[pub+.2(rand(rng)-.5) for _ in v],v;color,markersize=4)
    scatter!(b,[pub],[r.selection_nll-minimum(s.selection_nll for s in rank)];color,markersize=8)
    push!(counts,(publication_model=pub,total=length(initial),stopping_met=count(s->s.stopping_met,initial),stopping_not_met=count(s->!s.stopping_met,initial),selected_training_fit_stopping_met=r.training_stopping_met))
end
colgap!(f.layout,35)
export_figure(out,"split_model_comparison",f)
CSV.write(joinpath(out,"optimizer_bic.csv"),rows)
CSV.write(joinpath(out,"supplementary_convergence.csv"),counts)
write(joinpath(out,"split_model_comparison_caption.txt"),"(a) Calibration BIC across $(protocol["starts"]) independent optimizer starts per model on the calibration cultures. Boxplots show actual quartiles; points are optimizer outcomes, not biological replicates. All runs share one locked culture split. (b) Predictive NLL differences on held-out selection cultures using the best training-fit parameters for each model, without selection refitting. Smaller is better. Selection scoring uses training-derived noise weights. The two panels represent different criteria; panel (b) is not held-out BIC. This is a retrospective internal comparison because the complete dataset was examined during prior development.\n")
println("Split figures saved to $out")
