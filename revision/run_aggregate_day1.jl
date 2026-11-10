include("AggregateDayOne.jl")
length(ARGS)==3 || error("Usage: run_aggregate_day1.jl REFERENCE_FIT SPLIT_FIT OUTPUT")
AggregateDayOne.run_aggregate_dayone(ARGS...)
