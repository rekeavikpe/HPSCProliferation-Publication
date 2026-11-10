using TOML, Statistics, Random
function ablation_plot(dir;extended=false)
    ranking=collect(CSV.File(joinpath(dir,"ranking.csv")))
    errors=collect(CSV.File(joinpath(dir,"state_errors.csv")))
    order=extended ? ("baseline","shear_cell","shear_aggregate","shear_both","edr_cell","edr_aggregate","edr_both","full") : ("baseline","shear_only","edr_only","full")
    labels=extended ? ["Model 0","Shear: cell","Shear: aggregate","Shear: both","EDR: cell","EDR: aggregate","EDR: both","Model 6"] : ["Model 0","Shear only","EDR only","Model 6"]
    if extended
        runs=collect(CSV.File(joinpath(dir,"optimizer_bic.csv")))
    else
        runs=[(scenario=String(r.scenario),run=Int(r.run),bic=Float64(r.bic),stopping_met=Bool(r.stopping_met)) for r in CSV.File(joinpath(dir,"runs.csv"))]
        fit=get(ENV,"HPSC_ABLATION_SOURCE_FIT",TOML.parsefile(joinpath(dir,"manifest.toml"))["settings"]["source_fit"])
        source=joinpath(fit,"figure_data")
        isfile(joinpath(source,"optimizer_bic.csv"))||error("Missing original optimizer distributions; export source figure_data or set HPSC_ABLATION_SOURCE_FIT")
        meta=TOML.parsefile(joinpath(source,"source.toml"))
        for r in CSV.File(joinpath(source,"optimizer_bic.csv"))
            r.model in (0,5)||continue
            num=parse(Int,replace(r.run,"run"=>""));num<=meta["initial_starts_per_model"]||continue
            push!(runs,(scenario=r.model==0 ? "baseline" : "full",run=num,bic=Float64(r.bic),stopping_met=Bool(r.stopping_met)))
        end
    end
    selection=isfile(joinpath(dir,"selection_optimizer_bic.csv"))
    if selection
        ranking=collect(CSV.File(joinpath(dir,"selection_ranking.csv")))
        runs=collect(CSV.File(joinpath(dir,"selection_optimizer_bic.csv")))
    end
    n=length(order);colors=extended ? ["#777777","#27689F","#6EA2CB","#004B73","#B65D39","#D79D65","#804126","#00866A"] : ["#777777","#4E79A7","#AA4E35","#006B7F"]
    f=Figure(size=extended ? (1500,650) : (1200,560))
    ax=Axis(f[1,1],title=selection ? "(a) Model selection (40%)" : "(a) 100 mL calibration",xlabel="Model formulation",ylabel="BIC",xticks=(1:n,labels),xticklabelrotation=π/5,titlesize=21,xlabelsize=18,ylabelsize=18,xticklabelsize=15,yticklabelsize=15,xgridvisible=false)
    rng=MersenneTwister(20261006);diagnostics=NamedTuple[]
    for (i,s) in enumerate(order)
        rr=filter(r->r.scenario==s,runs);isempty(rr)&&error("No optimizer values for $s")
        v=Float64[r.bic for r in rr];all(isfinite,v)||error("Nonfinite BIC")
        boxplot!(ax,fill(i,length(v)),v;show_outliers=false,width=.55,color=(colors[i],.5),strokecolor=colors[i],whiskercolor=colors[i],mediancolor=:black)
        scatter!(ax,i .+ .22 .* (rand(rng,length(v)).-.5),v;color=(colors[i],.65),markersize=4)
        selected=only(filter(r->r.scenario==s,ranking))
        push!(diagnostics,(scenario=s,runs=length(rr),stopping_met=count(r->r.stopping_met,rr),selected_stopping_met=selected.stopping_met,selected_bic=selected.bic))
    end
    xlims!(ax,.4,n+.6)
    lo,hi=extrema(Float64[r.bic for r in runs]);pad=max(5.,.08*(hi-lo));ylims!(ax,lo-pad,hi+pad)
    ax2=Axis(f[1,2],title="(b) 500 mL prediction",xticks=(1:4,["Cells","Glucose","Lactate","Aggregate size"]),yticks=(1:n,labels),xticklabelrotation=π/5,titlesize=21,xticklabelsize=15,yticklabelsize=15)
    states=("cells","glucose","lactate","aggregate")
    z=[only(filter(r->r.scenario==s&&r.state==state&&r.cohort=="500mL retrospective prediction",errors)).mrpe_percent for state in states,s in order]
    h=heatmap!(ax2,1:4,1:n,z;colormap=:Blues)
    for i in 1:4,j in 1:n
        text!(ax2,i,j;text=string(round(z[i,j];digits=1)),align=(:center,:center),fontsize=18,color=z[i,j]>.6maximum(z) ? :white : :black)
    end
    ax2.yreversed=true;Colorbar(f[1,3],h,label="MRPE (%)",labelsize=18,ticklabelsize=15)
    colgap!(f.layout,32)
    stem=extended ? "ablation_eight_comparison" : "ablation_comparison"
    export_figure(joinpath(dir,"figures"),stem,f)
    CSV.write(joinpath(dir,"figures",stem*"_optimizer_diagnostics.csv"),diagnostics)
    CSV.write(joinpath(dir,"figures",stem*"_bic_values.csv"),runs)
    write(joinpath(dir,"figures",stem*"_caption.txt"),"Vertical boxplots show actual BIC across saved optimization runs, not biological replicates or confidence intervals. All finite values are retained without clipping; supplementary diagnostics report convergence. Baseline and full Model 6 distributions are reused from the saved calibration. $(selection ? "Parameters fitted to 60% of cultures. BIC evaluated on the held-out 40%." : "All-100-mL sensitivity.") 500-mL predictions are made without refitting. Exposure assignments indicate direct equation inputs; indirect effects through growth and coupled states remain.\n")
    stem
end
