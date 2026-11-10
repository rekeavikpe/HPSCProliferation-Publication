using CSV, TOML
include("../src/ProfileSummary.jl")
using .ProfileSummary
length(ARGS)==1||error("Usage: summarize_profiles.jl PROFILE_SUITE")
out=abspath(ARGS[1]); meta=TOML.parsefile(joinpath(out,"suite.toml")); id=meta["model"]
parameters=collect(CSV.File(joinpath(meta["source_fit"],"model$id/selected_parameters.csv")))
rows=NamedTuple[]
for p in parameters
    dir=joinpath(out,p.parameter)
    isfile(joinpath(dir,"COMPLETE"))||error("Incomplete profile: $(p.parameter)")
    profile=collect(CSV.File(joinpath(dir,"profile.csv")))
    s=summarize_profile(profile,p.value)
    push!(rows,(parameter=p.parameter,mle=p.value,lower_approx=s.better_minimum ? missing : s.lower.estimate,upper_approx=s.better_minimum ? missing : s.upper.estimate,lower_status=s.lower.status,upper_status=s.upper.status,lower_bracket_low=s.lower.bracket_low,lower_bracket_high=s.lower.bracket_high,upper_bracket_low=s.upper.bracket_low,upper_bracket_high=s.upper.bracket_high,failed_points=s.failed_points,better_minimum=s.better_minimum,minimum_twice_delta=s.minimum_twice_delta))
end
CSV.write(joinpath(out,"interval_diagnostics.csv"),rows)
write(joinpath(out,"interpretation.txt"),"Crossings use the nominal one-parameter 95% likelihood-ratio reference 3.841459, with log-linear interpolation inside the reported brackets. These are approximate working-likelihood intervals, not calibrated coverage guarantees. A bound is not a confidence limit. Missing crossings do not prove unboundedness outside the searched range. Boundary estimates, correlated observations and unsuccessful nuisance fits qualify interpretation. If better_minimum is true, refit the full model before interpreting intervals.\n")
