using CSV
include("../src/ReplicateSplits.jl")
using .ReplicateSplits
root=normpath(joinpath(@__DIR__,".."))
series=collect(CSV.File(joinpath(root,"data/processed/series.csv")))
out=isempty(ARGS) ? joinpath(root,"results/splits") : abspath(ARGS[1])
mkpath(out)
for repetition in 1:10
    split=random_split(series;seed=20261003+repetition-1)
    CSV.write(joinpath(out,"split$(repetition).csv"),split.assignments)
end
println("Saved 10 predetermined random 60/40 whole-culture splits. No outcomes were used to choose seeds.")
