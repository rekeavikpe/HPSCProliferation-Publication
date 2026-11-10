module BoundedLeastSquares
using LinearAlgebra
export bounded_lm

"Projected, damped Gauss-Newton for residuals with 0 <= x <= 1; diagnostic optimizer."
function bounded_lm(residual,x0;maxiter=250,h=3e-5,gtol=1e-4)
    x=clamp.(copy(x0),0.0,1.0);r=residual(x);f=sum(abs2,r)/2
    damping=1e-3;stalled=0;pg=Inf;calls=1
    for iter in 1:maxiter
        J=zeros(length(r),length(x))
        for j in eachindex(x)
            a=copy(x);b=copy(x);a[j]=max(0.0,x[j]-h);b[j]=min(1.0,x[j]+h)
            J[:,j]=(residual(b).-residual(a))./(b[j]-a[j]);calls+=2
        end
        g=J'*r;pg=maximum(abs.(x.-clamp.(x.-g,0.0,1.0)))
        pg<gtol && return x,(converged=true,iterations=iter,reason="projected_gradient",projected_gradient=pg,calls=calls)
        active=[(x[i]<=1e-10 && g[i]>0)||(x[i]>=1-1e-10 && g[i]<0) for i in eachindex(x)]
        free=findall(.!active);isempty(free) && return x,(converged=true,iterations=iter,reason="active_bounds",projected_gradient=pg,calls=calls)
        H=J[:,free]'*J[:,free];diagH=max.(diag(H),1e-6)
        accepted=false;gain=0.0
        for attempt in 1:12
            step=zeros(length(x));step[free]=-(H+damping*Diagonal(diagH))\g[free]
            maximum(abs.(step))>.25 && (step .*= .25/maximum(abs.(step)))
            trial=clamp.(x.+step,0.0,1.0)
            rr=residual(trial);calls+=1;ff=sum(abs2,rr)/2
            if isfinite(ff) && ff<f
                gain=f-ff;x=trial;r=rr;f=ff;damping=max(damping/3,1e-12);accepted=true;break
            end
            damping*=10
        end
        stalled=(accepted && gain>1e-8) ? 0 : stalled+1
        stalled>=5 && return x,(converged=false,iterations=iter,reason="objective_stalled",projected_gradient=pg,calls=calls)
    end
    x,(converged=false,iterations=maxiter,reason="iteration_limit",projected_gradient=pg,calls=calls)
end
end
