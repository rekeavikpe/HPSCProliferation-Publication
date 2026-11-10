include("../profile_figures/render.jl")
length(ARGS)==1||error("Usage: --project=profile_figures controlled/render.jl RESULTS")
for case in CSV.File(joinpath(ARGS[1],"design.csv"))
    source=joinpath(ARGS[1],String(case.name))
    isfile(joinpath(source,"COMPLETE"))||error("Incomplete comparison: $(case.name)")
    render_profiles(source,joinpath(source,"figures"))
end
