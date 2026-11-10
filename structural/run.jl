using StructuralIdentifiability, CSV, Random, Logging, TOML
include("models.jl")
using .StructuralModels
length(ARGS)>=1 || error("Usage: structural/run.jl MODEL [pooled|single40] [SEGMENT 0–5] [local|global]")
id=parse(Int,ARGS[1]);design=length(ARGS)>=2 ? Symbol(ARGS[2]) : :pooled;segment=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 0
mode=length(ARGS)>=4 ? Symbol(ARGS[4]) : :local
mode in (:local,:global)||error("Mode must be local or global")
base=length(ARGS)>=5 ? abspath(ARGS[5]) : joinpath(ROOT,"results/structural")
out=joinpath(base,"$(design)_segment$(segment)_$(mode)/model$id")
ispath(out)&&error("Output already exists; archive or use existing result")
mkpath(out)
write(joinpath(out,"model.jl"),model_source(id;design,segment)*"\n")
Random.seed!(20261004+100id+10segment)
ode=build_model(id;design,segment)
t=time()
result=if mode==:local
    assess_local_identifiability(ode;funcs_to_check=ode.parameters,prob_threshold=.999,loglevel=Logging.Info)
else
    assess_identifiability(ode;funcs_to_check=ode.parameters,prob_threshold=.999,loglevel=Logging.Info)
end
rows=[(model=id,design=string(design),segment=segment,parameter=string(k),classification=string(v),mode=string(mode)) for (k,v) in result]
CSV.write(joinpath(out,"parameters.csv"),rows)
open(joinpath(out,"run.toml"),"w") do io
    TOML.print(io,Dict("model"=>id,"design"=>string(design),"segment"=>segment,"mode"=>string(mode),"probability_threshold"=>.999,"seed"=>20261004+100id+10segment,"elapsed_seconds"=>time()-t,"julia"=>string(VERSION),"initial_conditions"=>"generic; not specialized to the nominal inoculum", "clock"=>"E is observed/known, E prime = -E/12; nonzero generic initial value", "scope"=>"Rational exposure segment, ideal continuous observation of all four states; no claim of practical identifiability or full hybrid global identifiability"))
end
write(joinpath(out,"COMPLETE"),"complete\n")
println(rows)
