using CSV, TOML, Statistics, SHA
include("Common.jl")
using .RevisionCommon, .RevisionCommon.Workflow, .RevisionCommon.Workflow.Calibration, .RevisionCommon.Workflow.Calibration.Proliferation
include("CellContexts.jl")
const EXPOSURES=("lag_p50","lag_mean","lag_p25","lag_p75","lag_p90","lag_max","eul_va","eul_p50","eul_p75","eul_p90","eul_max")
function exposure_context(base,scenario)
    inputs=collect(CSV.File(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))
    series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
    conditions=copy(base.conditions)
    for (c,q) in base.conditions
        info=first(filter(r->r.condition==c,series))
        rows=sort(filter(r->r.scenario==scenario&&r.volume_ml==info.volume_ml&&r.rpm==info.rpm,inputs);by=r->r.size_um)
        lookup=if startswith(scenario,"lag_")
            @assert [r.size_um for r in rows]==[200,400,600,800,1000]
            linear_exposure([r.size_um for r in rows],[r.shear_pa for r in rows],[r.edr_m2_s3 for r in rows])
        else
            r=only(rows);let tau=r.shear_pa,edr=r.edr_m2_s3;size->(tau,edr);end
        end
        conditions[c]=merge(q,(;lookup))
    end
    merge(base,(;conditions))
end
function verify_split(prep,splitfit)
    saved=collect(CSV.File(joinpath(splitfit,"split.csv")))
    Set(String(r.series_id) for r in saved if r.role=="calibration")==Set(prep.split.training)||error("Training split differs")
    Set(String(r.series_id) for r in saved if r.role=="selection")==Set(prep.split.selection)||error("Selection split differs")
    for id in 0:6
        m=legacy_model(id);p=read_parameters(joinpath(splitfit,"model$id/selected_parameters.csv"),m)
        actual=gaussian_nll(prep.train.y,prediction(prep.train,m,p),prep.train.sigma)
        recorded=only(CSV.File(joinpath(splitfit,"model$id/summary.csv"))).nll
        isapprox(actual,recorded;atol=1e-6,rtol=0)||error("Saved training likelihood differs for model$id")
    end
end
function main()
    length(ARGS)==4||error("Usage: run_david_analyses.jl selection|exposure|cells FINAL_FIT SPLIT_FIT OUTPUT")
    stage,fit,splitfit,out=ARGS;fit=abspath(fit);splitfit=abspath(splitfit);out=abspath(out)
    cfgfile=joinpath(fit,"calibration.toml");cfg=TOML.parsefile(cfgfile)
    isempty(cfg["fixed"])||error("Expected all-free reference fit")
    cfg["bounds"]["K_G"][2]==100&&cfg["bounds"]["K_Lac"][2]==100||error("Expected approved 100 mM bounds")
    starts=parse(Int,get(ENV,"HPSC_REVIEW_STARTS","100"));budget=parse(Int,get(ENV,"HPSC_REVIEW_EVALUATIONS","20000"));iterations=parse(Int,get(ENV,"HPSC_REVIEW_ITERATIONS","3000"))
    starts>0&&budget>=100&&iterations>0||error("Invalid budgets")
    files=[cfgfile,joinpath(fit,"model0/selected_parameters.csv"),joinpath(fit,"model5/selected_parameters.csv")]
    stage!="cells"&&append!(files,[joinpath(splitfit,"split.csv"),[joinpath(splitfit,"model$id/selected_parameters.csv") for id in 0:6]...])
    lock_inputs(out,Dict("stage"=>stage,"starts"=>starts,"budget"=>budget,"iterations"=>iterations,"julia"=>string(VERSION),"cell_policy"=>"biological7 and metabolic5, plus baseline, hydrodynamic-only and all-active"),files)
    rows=NamedTuple[];metrics=NamedTuple[];curves=NamedTuple[];observed=NamedTuple[];policies=NamedTuple[]
    if stage in ("selection","exposure")
        prep=SplitSelection.prepare(cfgfile);verify_split(prep,splitfit)
        CSV.write(joinpath(out,"split.csv"),prep.split.assignments)
        scenarios=stage=="selection" ? string.(0:6) : collect(EXPOSURES)
        for scenario in scenarios
            id=stage=="selection" ? parse(Int,scenario) : 5;m=legacy_model(id)
            train=stage=="selection" ? prep.train : exposure_context(prep.train,scenario)
            held=stage=="selection" ? prep.held : exposure_context(prep.held,scenario)
            p,summary=if stage=="selection"
                p=read_parameters(joinpath(splitfit,"model$id/selected_parameters.csv"),m)
                nll=gaussian_nll(train.y,prediction(train,m,p),train.sigma)
                p,(nll=nll,k=length(free_names(train,m)),bic=bic(nll,length(free_names(train,m)),length(train.y)))
            else
                # No all-data estimate is used as an anchor for held-out selection.
                fit_context(train,m,joinpath(out,scenario);starts,budget,iterations)
            end
            push!(rows,(scenario,publication_model=id==5 ? 6 : id==6 ? 5 : id,k=summary.k,training_nll=summary.nll,training_bic=summary.bic,selection_nll=gaussian_nll(held.y,prediction(held,m,p),held.sigma)))
            for (cohort,ctx) in (("calibration60",train),("selection40",held))
                append!(metrics,[merge((;scenario,cohort),r) for r in errors(ctx,prediction(ctx,m,p))])
            end
            CSV.write(joinpath(out,"ranking.csv"),sort(rows;by=r->r.selection_nll));CSV.write(joinpath(out,"state_errors.csv"),metrics)
        end
    elseif stage=="cells"
        for line in ("585A1","SCTi003-A")
            base=cell_context(line,cfg)
            append!(observed,[merge((;cell_line=line),r) for r in base.observations])
            for strategy in ("transfer_model0","transfer_model6","baseline_refit","biological7","metabolic5","hydrodynamic2","all11")
                m=strategy in ("transfer_model0","baseline_refit") ? 0 : 6
                parent=read_parameters(joinpath(fit,"model$(m==0 ? 0 : 5)/selected_parameters.csv"),m)
                free=strategy=="biological7" ? BIOLOGICAL7 : strategy=="metabolic5" ? METABOLIC5 : strategy=="hydrodynamic2" ? (:K_tau,:K_eps) : Tuple(parameter_names(m))
                transfer=startswith(strategy,"transfer_")
                c=deepcopy(cfg);c["fixed"]=Dict(string(n)=>getproperty(parent,n) for n in keys(parent) if transfer||!(n in free))
                ctx=merge(base,(config=c,));dir=joinpath(out,line,strategy);mkpath(dir)
                append!(policies,[(cell_line=line,strategy,parameter=string(n),estimated=!haskey(c["fixed"],string(n)),source_value=getproperty(parent,n)) for n in keys(parent)])
                p=if transfer;parent;else;first(fit_context(ctx,m,dir;starts,budget,iterations,anchor=parent));end
                save_parameters(joinpath(dir,"selected_parameters.csv"),p)
                yp=prediction(ctx,m,p)
                push!(rows,(cell_line=line,strategy,estimated_parameters=length(free_names(ctx,m)),weighted_ssr=sum(abs2,(yp.-ctx.y)./ctx.sigma),assessment=transfer ? "transfer without refitting" : "cell-line calibration"))
                append!(metrics,[merge((;cell_line=line,strategy),r) for r in errors(ctx,yp)])
                append!(curves,[merge((;cell_line=line,strategy),r) for r in export_curves(ctx,m,p)])
                CSV.write(joinpath(out,"fits.csv"),rows);CSV.write(joinpath(out,"state_errors.csv"),metrics);CSV.write(joinpath(out,"curves.csv"),curves)
            end
        end
        CSV.write(joinpath(out,"observations.csv"),observed);CSV.write(joinpath(out,"parameter_policy.csv"),policies)
        write(joinpath(out,"interpretation.txt"),"All refits use the same original state-scaled weighted SSR and final ESI-017 parameter bounds. Seven-parameter refitting therefore preserves the parameter set, not the old narrower bounds. Refitted errors are in-sample. Transfer errors use unchanged ESI-017 parameters. The author confirmed dextran sulfate was used experimentally for 585A1, superseding the raw CSV absent label. No explicit dextran-sulfate term is modeled. Overall MRPE pools eligible condition-day-state observations. Experimental Noise is mean sample SD / mean × 100.\n")
    else;error("Unknown stage $stage");end
    write(joinpath(out,"COMPLETE"),"Completed $(stage)\n")
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
