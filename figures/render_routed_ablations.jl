using CSV, CairoMakie
include("ManuscriptStyle.jl")
using .ManuscriptStyle
include("AblationPlots.jl")
length(ARGS)==1||error("Usage: render_routed_ablations.jl RESULTS")
dir=abspath(only(ARGS));isfile(joinpath(dir,"COMPLETE"))||error("Calibration suite incomplete")
manuscript_theme!();stem=ablation_plot(dir;extended=true)
shared=get(TOML.parsefile(joinpath(dir,"manifest.toml"))["settings"],"shared_response_parameter",false)
scope_note=isfile(joinpath(dir,"selection_optimizer_bic.csv")) ? "Parameters fitted to the 60% training cultures. BIC evaluated on the 40% selection cultures." : "Exploratory all-100-mL comparison."
response_note=shared ? "Single-exposure effects on both equations share one response parameter." : "Single-exposure effects on both equations have separate response parameters."
write(joinpath(dir,"index.html"),"""<!doctype html><meta charset="utf-8"><title>Eight exposure assignments</title><style>body{font:16px Arial;max-width:1500px;margin:30px auto;line-height:1.5}img{width:100%}</style><h1>Eight exposure assignments</h1><p>$scope_note Assignments indicate direct inputs into the growth-inhibition and aggregate-loss terms; indirect coupling remains. The baseline and primary Model 6 are preserved. $response_note No 500 mL refitting.</p><img src="figures/$stem.svg"><p><a href="ranking.csv">BIC and parameter counts</a> · <a href="state_errors.csv">MRPE and Experimental Noise</a> · <a href="figures/$(stem)_optimizer_diagnostics.csv">Optimizer diagnostics</a></p><p>Profile results are saved separately under profiles/all-active and profiles/conditional for each formulation.</p>""")
