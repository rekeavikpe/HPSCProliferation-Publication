# Included by run_david_analyses.jl after importing RevisionCommon.
const BIOLOGICAL7=(:mu_max,:alpha,:K_G,:K_oG,:Y_GLac,:K_d,:K_Agg)
const METABOLIC5=(:K_G,:K_Lac,:K_oG,:Y_GX,:Y_GLac)
function cell_context(line,cfg)
    raw=filter(r->r.cell_line==line&&r.status=="measured",collect(CSV.File(joinpath(ROOT,"data/processed/cell_line_measurements.csv"))))
    rpm=only(unique(r.rpm for r in raw));days=Int(maximum(r.day for r in raw));initial(s)=mean(Float64[r.value for r in raw if r.day==0&&r.state==s])
    cfd=collect(CSV.File(joinpath(ROOT,"data/processed/lagrangian_exposure.csv")))
    grids=[sort(filter(r->r.volume_ml==100&&r.rpm==rpm0&&r.percentile==cfg["exposure_percentile"],cfd);by=r->r.size_um) for rpm0 in (40,60)]
    a,b=grids; @assert [r.size_um for r in a]==[r.size_um for r in b]
    w=(rpm-40)/20;lookup=linear_exposure([r.size_um for r in a],[(1-w)*x.shear_pa+w*y.shear_pa for (x,y) in zip(a,b)],[(1-w)*x.edr_m2_s3+w*y.edr_m2_s3 for (x,y) in zip(a,b)])
    obs=NamedTuple[]
    for day in 1:days,s in STATES
        vals=Float64[r.value for r in raw if r.day==day&&r.state==s];isempty(vals)&&continue
        push!(obs,(condition=1,day=Float64(day),state=s,mean=mean(vals),sd=length(vals)>1 ? std(vals) : 0.,n_vessels=length(vals)))
    end
    glucose=initial("glucose");feeds=[(time_h=24.0d,fraction=line=="SCTi003-A"&&d==2 ? .5 : .8) for d in 2:days-1]
    q=(u0=[initial("cells"),glucose,initial("lactate"),Float64(cfg["initial_aggregate_um"])],glucose,lookup,times=collect(24.:24.:24days),feeds)
    # Preserve original cell-line weighted SSR across every strategy. Not a comparison of BIC to ESI-017 fits.
    scales=Dict(zip(STATES,(5e5,5.,5.,100.)));weights=Dict(zip(STATES,(1.,1.,1.,1.5)))
    (;config=deepcopy(cfg),observations=obs,conditions=Dict(1=>q),y=Float64[r.mean for r in obs],sigma=[scales[r.state]/sqrt(weights[r.state]) for r in obs])
end
