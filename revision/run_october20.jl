include("October20.jl")
length(ARGS)==4||error("Usage: run_october20.jl reuse|p75|day1 FIT SPLIT OUTPUT")
run_october(ARGS...)
