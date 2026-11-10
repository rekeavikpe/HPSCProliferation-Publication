# Keep the existing ProfileLikelihood.jl line, shading and reference-line style.
include("../profile_figures/render.jl")
length(ARGS)==2 || error("Usage: render_fixed_two.jl PROFILE_RESULTS OUTPUT")
CairoMakie.set_theme!(fonts=(regular="Arial",bold="Arial Bold",italic="Arial Italic",bold_italic="Arial Bold Italic"))
render_profiles(ARGS...)
