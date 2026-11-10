module HybridRefinement
using CSV, TOML, Optim
using ..Workflow
using ..Workflow.Calibration
using ..Workflow.Calibration.Proliferation
export refine_original

function refine_original(ctx,model,z0,dir;iterations=3000)
    mkpath(dir)
    residual=Workflow.residual_function(ctx,model)
    # The manuscript's scaled_to_params clips to the unit hypercube.
    objective=z->all(isfinite,z) ? sum(abs2,residual(clamp.(z,0.,1.)))/2 : 1e100
    start=time()
    result=optimize(objective,copy(z0),NelderMead(),Optim.Options(iterations=iterations,g_abstol=1e-8))
    raw=Optim.minimizer(result);z=clamp.(raw,0.,1.)
    p=decode(ctx,model,z)
    nll=gaussian_nll(ctx.y,prediction(ctx,model,p),ctx.sigma)
    save_parameters(joinpath(dir,"polished_parameters.csv"),p)
    CSV.write(joinpath(dir,"coordinates.csv"),[(parameter=string(n),raw=x,feasible=y) for (n,x,y) in zip(free_names(ctx,model),raw,z)])
    stages=[(stage="polished",nll=nll,stopping_met=Optim.converged(result))]
    CSV.write(joinpath(dir,"stages.csv"),stages)
    open(joinpath(dir,"diagnostics.toml"),"w") do io
        TOML.print(io,Dict("method"=>"NelderMead with original clipped-coordinate objective",
            "iterations"=>Optim.iterations(result),"iteration_limit"=>iterations,
            "stopping_met"=>Optim.converged(result),"elapsed_seconds"=>time()-start,
            "clipped_coordinates"=>count(x->x<0||x>1,raw)))
    end
    write(joinpath(dir,"COMPLETE"),"saved actual optimizer status\n")
    stages
end
end
