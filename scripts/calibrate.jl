using TOML
include("../src/Workflow.jl")
using .Workflow
length(ARGS)<=3 || error("Usage: calibrate.jl [OUTPUT_DIRECTORY] [STARTS_PER_MODEL] [CONFIG_FILE]")
out=isempty(ARGS) ? joinpath(ROOT,"results/calibration") : abspath(ARGS[1])
starts=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 6
starts>=1 || error("At least one start per model is required")
cfg=length(ARGS)>=3 ? TOML.parsefile(abspath(ARGS[3])) : config()
run_calibration(out;ctx=context(cfg=cfg),nstarts=starts)
