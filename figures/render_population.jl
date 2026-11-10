using CairoMakie, CSV, Statistics, Random
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)>=3 || error("Usage: render_population.jl REFERENCE_BIC_CSV OUT POPULATION_RESULT_DIR...")
source,out=abspath.(ARGS[1:2])
old=filter(r->parse(Int,replace(r.run,"run"=>""))<=100,collect(CSV.File(source)))
new=vcat([collect(CSV.File(joinpath(dir,"runs.csv"))) for dir in ARGS[3:end]]...)
ids=sort(unique(Int[r.model for r in new]))
CairoMakie.activate!();manuscript_theme!()
ncols=min(4,length(ids));nrows=cld(length(ids),ncols)
fig=Figure(size=(225*ncols,230*nrows+60));rng=MersenneTwister(20261004)
table=NamedTuple[]
for (i,id) in enumerate(ids)
    before=[r.bic for r in old if r.model==id];after=[r.bic for r in new if r.model==id]
    best=min(minimum(before),minimum(after))
    row,col=divrem(i-1,ncols)
    ax=Axis(fig[row+1,col+1],title="Model $(publication_model(id))",xticks=([1,2],["Local starts\n(n=$(length(before)))","Population + local\n(n=$(length(after)))"]),
        ylabel=col==0 ? "ΔBIC from best achieved fit" : "",yscale=log10,
        yticks=([1.,11.,101.,1001.,10001.,100001.,1000001.],["0","10","100","1,000","10,000","100,000","1,000,000"]))
    for (j,values) in enumerate((before,after))
        y=1 .+ max.(0.,values.-best)
        scatter!(ax,j .+ .18 .* (rand(rng,length(y)).-.5),y;color=j==1 ? (:gray40,.6) : (SCALE_COLOR,.9),markersize=4)
        lines!(ax,[j-.2,j+.2],fill(median(y),2);color=:black,linewidth=1)
    end
    xlims!(ax,.6,2.4);ylims!(ax,.8,max(12.,1.3*(maximum(before)-best+1),1.3*(maximum(after)-best+1)))
    push!(table,(model=id,old_n=length(before),new_n=length(after),old_best=minimum(before),new_best=minimum(after),new_worst=maximum(after),new_range=maximum(after)-minimum(after)))
end
Label(fig[nrows+1,1:ncols],"All initial outcomes retained. Axis uses log₁₀(1 + ΔBIC); black bars are medians.\nFive-repeat diagnostic, not final BIC boxplots.",fontsize=9)
colgap!(fig.layout,22);rowgap!(fig.layout,22)
export_figure(out,"population_comparison",fig)
CSV.write(joinpath(out,"comparison.csv"),table)
