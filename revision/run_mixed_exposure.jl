include("MixedExposure.jl")
length(ARGS)==3 || error("Usage: run_mixed_exposure.jl REFERENCE_FIT SPLIT_FIT OUTPUT")
MixedExposure.run_mixed(ARGS...)
