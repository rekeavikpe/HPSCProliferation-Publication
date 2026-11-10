using CSV, TOML
include("../src/TargetedChecks.jl")
using .TargetedChecks, .TargetedChecks.Workflow, .TargetedChecks.Workflow.Calibration, .TargetedChecks.Workflow.Calibration.Proliferation
length(ARGS)==2||error("Usage: diagnose_precision.jl FIT_DIRECTORY OUTPUT_DIRECTORY")
fitdir=abspath(ARGS[1]);out=abspath(ARGS[2]);ctx=checked_context(fitdir)
prepare_output(out,fitdir;settings=Dict("analysis"=>"ODE precision and projected finite-difference gradients"))
scores=NamedTuple[];gradients=NamedTuple[]
for id in (1,2,5)
    m=legacy_model(id);p=read_parameters(joinpath(fitdir,"model$id/selected_parameters.csv"),m);z=encode(ctx,m,p)
    for tolerance in (1e-11,1e-12,1e-13,1e-14)
        values=prediction(ctx,m,p;reltol=tolerance,abstol=tolerance/10)
        push!(scores,(model=id,reltol=tolerance,nll=gaussian_nll(ctx.y,values,ctx.sigma)))
    end
    for tolerance in (1e-12,1e-13),step in (1e-4,1e-5,1e-6)
        residual=Workflow.residual_function(ctx,m;reltol=tolerance,abstol=tolerance/10)
        r=residual(z);g=zeros(length(z))
        for j in eachindex(z)
            a=copy(z);b=copy(z);a[j]=max(0,z[j]-step);b[j]=min(1,z[j]+step)
            g[j]=sum(r.*((residual(b).-residual(a))./(b[j]-a[j])))
        end
        projected=z.-clamp.(z.-g,0,1)
        for (j,name) in enumerate(free_names(ctx,m))
            push!(gradients,(model=id,reltol=tolerance,step=step,parameter=string(name),gradient=g[j],projected_gradient=projected[j],normalized_value=z[j]))
        end
    end
end
CSV.write(joinpath(out,"scores.csv"),scores)
CSV.write(joinpath(out,"gradients.csv"),gradients)
metadata=TOML.parsefile(joinpath(out,"run.toml"));metadata["status"]="complete"
open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,metadata);end
write(joinpath(out,"COMPLETE"),"complete; diagnosis only, no tolerance changes\n")
println(scores)
