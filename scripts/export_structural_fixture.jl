using CSV
include("../src/Workflow.jl")
using .Workflow, .Workflow.Calibration.Proliferation
out=joinpath(ROOT,"reference/structural_fixture");mkpath(out)
ctx=context();rows=NamedTuple[];parameter_rows=NamedTuple[]
for id in 0:6
    model=legacy_model(id)
    p=read_parameters(joinpath(ROOT,"reference/reproduction/model$id/selected_parameters.csv"),model)
    append!(parameter_rows,[(model=id,parameter=string(n),value=getproperty(p,n)) for n in keys(p)])
    for segment in 0:5,(j,c) in enumerate((1,2,3,5,6))
        size=(100.,300.,500.,700.,900.,1100.)[segment+1]
        u=[80000.,10.,3.,size];tau,edr=ctx.conditions[c].lookup(size)
        du=zeros(4);model_rhs!(du,u,p,12.,model,tau,edr)
        push!(rows,(model=id,segment=segment,experiment=j,X=u[1],G=u[2],L=u[3],A=u[4],E=exp(-1.),dX=du[1],dG=du[2],dL=du[3],dA=du[4]))
    end
end
CSV.write(joinpath(out,"rhs.csv"),rows)
CSV.write(joinpath(out,"parameters.csv"),parameter_rows)
println("Saved 210 equation fixtures from the numerical fitting kernel")
