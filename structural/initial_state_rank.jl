using StructuralIdentifiability, CSV, Random, TOML, LinearAlgebra
include("models.jl")
using .StructuralModels
const SI=StructuralIdentifiability
const NN=SI.Nemo
id=isempty(ARGS) ? 0 : parse(Int,ARGS[1])
base=length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(ROOT,"results/structural")
out=joinpath(base,"nominal_initial_state/model$id")
ispath(out)&&error("Output exists: $out")
mkpath(out)
prime=2147483647;precision=16;seed=20261004+id
ode=SI.reduce_ode_mod_p(build_model(id;design=:pooled,segment=0),prime)
field=NN.base_ring(ode.poly_ring);P=eltype(ode.x_vars);T=typeof(field(0))
rows=collect(CSV.File(joinpath(ROOT,"data/processed/measurements.csv");types=String))
initial=Dict("E"=>1//big(1))
for (j,c) in enumerate((1,2,3,5,6))
    glucose=[StructuralModels.decimal_rational(r.value) for r in rows if r.condition==string(c)&&r.day=="0"&&r.state=="glucose"&&r.status=="measured"]
    isempty(glucose)&&error("Missing glucose")
    initial["X$j"]=(20000,80000,150000,80000,80000)[j]//big(1)
    initial["G$j"]=sum(glucose)/length(glucose)
    initial["L$j"]=0//big(1);initial["A$j"]=15//big(1)
end
tofield(q)=field(numerator(q))/field(denominator(q))
state_name(x)=replace(string(x),"(t)"=>"")
ic=Dict{P,T}(x=>tofield(initial[state_name(x)]) for x in ode.x_vars)
rng=MersenneTwister(seed);params=Dict{P,T}(p=>field(rand(rng,1:prime-1)) for p in ode.parameters)
derivatives=SI.differentiate_output(ode,params,ic,Dict{P,Vector{T}}(),precision)
cols=[(y,j) for y in ode.y_vars for j in 0:precision-1]
matrix=NN.zero_matrix(field,length(ode.parameters),length(cols))
for (i,p) in enumerate(ode.parameters),(k,(y,j)) in enumerate(cols)
    matrix[i,k]=NN.coeff(derivatives[y][p],j)
end
r=rank(matrix);k=length(ode.parameters)
CSV.write(joinpath(out,"parameter_witness.csv"),[(parameter=string(p),value=string(params[p])) for p in ode.parameters])
CSV.write(joinpath(out,"initial_state.csv"),[(state=string(x),rational_value=string(initial[state_name(x)])) for x in ode.x_vars])
CSV.write(joinpath(out,"sensitivity_matrix.csv"),[(parameter=string(p),output=string(y),order=j,value=string(matrix[i,c])) for (i,p) in enumerate(ode.parameters) for (c,(y,j)) in enumerate(cols)])
open(joinpath(out,"result.toml"),"w") do io
    TOML.print(io,Dict("model"=>id,"rank"=>r,"parameters"=>k,"prime"=>prime,"precision"=>precision,"seed"=>seed,"full_rank"=>r==k,"scope"=>"Exact finite-field parameter sensitivity rank, fixed nominal initial states, pooled 100-mL lower exposure segment. Full rank certifies generic local identifiability at these known initial states; rank deficiency alone is inconclusive. Not practical or global identifiability."))
end
r==k || error("Rank $r/$k: inconclusive; investigate before claiming identifiability")
write(joinpath(out,"COMPLETE"),"complete\n")
println("Model $id: exact rank $r/$k at nominal initial states")
