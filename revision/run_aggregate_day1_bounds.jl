include("AggregateDayOneBounds.jl")
length(ARGS) in (4,5) || error("Usage: run_aggregate_day1_bounds.jl REFERENCE_FIT SPLIT_FIT COMPLETED_D1_RESULTS OUTPUT [--smoke]")
length(ARGS)==5 && ARGS[5]!="--smoke" && error("Unknown option")
AggregateDayOneBounds.run_checks(ARGS[1:4]...;smoke="--smoke" in ARGS)
