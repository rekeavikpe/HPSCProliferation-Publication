module ReplicateSplits
using CSV, Random
export random_split
"Randomly partition whole bioreactors within each operating condition."
function random_split(series;fraction=.6,seed=20261003)
    0<fraction<1 || throw(ArgumentError("Fraction must lie between zero and one"))
    rng=MersenneTwister(seed)
    assignments=NamedTuple[]
    training=String[];selection=String[]
    for c in (1,2,3,5,6)
        ids=sort(String[r.series_id for r in series if r.condition==c && r.volume_ml==100])
        length(ids)>=2 || error("Too few cultures in condition $c")
        order=shuffle(rng,ids)
        ntrain=clamp(round(Int,fraction*length(ids)),1,length(ids)-1)
        for (i,id) in enumerate(order)
            role=i<=ntrain ? "calibration" : "selection"
            push!(i<=ntrain ? training : selection,id)
            push!(assignments,(legacy_condition=c,series_id=id,role=role,seed=seed))
        end
    end
    (;training,selection,assignments)
end
end
