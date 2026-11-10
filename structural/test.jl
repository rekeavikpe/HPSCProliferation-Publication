using StructuralIdentifiability, CSV, Test, Random
include("models.jl")
using .StructuralModels
const SI=StructuralIdentifiability
const NN=SI.Nemo
rows=collect(CSV.File(joinpath(ROOT,"reference/structural_fixture/rhs.csv");types=String))
parameters=collect(CSV.File(joinpath(ROOT,"reference/structural_fixture/parameters.csv");types=String))
qq(s)=NN.QQ(StructuralModels.decimal_rational(s))
@testset "Symbolic and numerical equation agreement" begin
    for id in 0:6,segment in 0:5
        ode=build_model(id;design=:pooled,segment)
        values=Dict(string(r.parameter)=>qq(r.value) for r in parameters if r.model==string(id))
        data=filter(r->r.model==string(id)&&r.segment==string(segment),rows)
        values["E(t)"]=qq(first(data).E)
        for r in data,s in ("X","G","L","A")
            values["$s$(r.experiment)(t)"]=qq(getproperty(r,Symbol(s)))
        end
        assignment=Dict(v=>values[string(v)] for v in vcat(ode.x_vars,ode.parameters))
        for r in data,s in ("X","G","L","A")
            x=only(filter(v->string(v)=="$s$(r.experiment)(t)",ode.x_vars))
            actual=Float64(SI.eval_at_dict(ode.x_equations[x],assignment))
            @test isapprox(actual,parse(Float64,getproperty(r,Symbol("d"*s)));rtol=2e-12,atol=1e-10)
        end
    end
end
@testset "Identifiability positive and negative controls" begin
    Random.seed!(20261004)
    identifiable=@ODEmodel(x'(t)=a*x(t),y(t)=x(t))
    confounded=@ODEmodel(x'(t)=a*b*x(t),y(t)=x(t))
    @test all(values(assess_local_identifiability(identifiable;funcs_to_check=identifiable.parameters,prob_threshold=.999)))
    @test !any(values(assess_local_identifiability(confounded;funcs_to_check=confounded.parameters,prob_threshold=.999)))
end
