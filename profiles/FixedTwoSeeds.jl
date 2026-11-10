module FixedTwoSeeds
using CSV
export load_seed_bank, seed_vectors, seed_provider

# A separate function call owns each bank. In particular, neither the profile
# rows nor the fixed coordinate can be a Core.Box shared by threaded jobs.
struct SeedBank
    coordinates::Vector{Float64}
    vectors::Vector{Vector{Float64}}
    reference::Vector{Float64}
end

function load_seed_bank(directory, names, index, bounds, mle)
    rows=collect(CSV.File(joinpath(directory,string(names[index]),"profile.csv")))
    saved=collect(CSV.File(joinpath(directory,string(names[index]),"nuisance_parameters.csv")))
    others=setdiff(eachindex(names),[index])
    isempty(rows) && error("Empty seed profile: $(names[index])")
    coordinates=Float64.(getproperty.(rows,:normalized_parameter))
    all(isfinite,coordinates) || error("Nonfinite seed coordinates")
    vectors=Vector{Float64}[]
    for point in eachindex(rows)
        entries=filter(v->v.point==point,saved)
        length(unique(v.parameter for v in entries))==length(entries) || error("Duplicate seed parameter")
        values=Dict(String(v.parameter)=>Float64(v.value) for v in entries)
        all(k->haskey(values,string(names[k])) && isfinite(values[string(names[k])]) && values[string(names[k])]>0,others) || error("Incomplete seed vector for $(names[index]) point $point")
        push!(vectors,[clamp(log(values[string(names[k])]/bounds[string(names[k])][1])/log(bounds[string(names[k])][2]/bounds[string(names[k])][1]),0.,1.) for k in others])
    end
    SeedBank(coordinates,vectors,copy(mle[others]))
end

function seed_vectors(bank::SeedBank, coordinate)
    isfinite(coordinate) || error("Nonfinite fixed coordinate")
    nearest=sortperm(abs.(bank.coordinates.-coordinate))[1:min(3,length(bank.coordinates))]
    # Optimizers mutate their starting vectors. Always return independent copies.
    [copy(bank.reference),[copy(bank.vectors[j]) for j in nearest]...]
end

function seed_provider(bank::SeedBank)
    function provide(fixed_problem)
        hasproperty(fixed_problem.f.f,Symbol("θₙ")) || error("ProfileLikelihood fixed-coordinate interface changed")
        coordinate=Float64(getproperty(fixed_problem.f.f,Symbol("θₙ")))
        seed_vectors(bank,coordinate)
    end
end
end
