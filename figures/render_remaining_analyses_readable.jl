using CSV, CairoMakie
include("ManuscriptStyle.jl")
using .ManuscriptStyle
length(ARGS)==2 || error("Usage: render_remaining_analyses.jl oxygen|ablation RESULTS")
kind,dir=ARGS;isfile(joinpath(dir,"COMPLETE")) || error("Analysis incomplete")
manuscript_theme!()
table_html=""
if kind=="oxygen"
    rows=collect(CSV.File(joinpath(dir,"sensitivity.csv")))
    summary=collect(CSV.File(joinpath(dir,"condition_summary.csv")))
    sort!(summary;by=s->(s.volume_ml==500,s.source_condition))
    f=Figure(size=(1000,650))
    for (i,s) in enumerate(summary)
        oxygen_ax=Axis(f[(i-1)÷3+1,(i-1)%3+1],title="($(Char(96+i))) $(replace(s.label,"100mL"=>"100 mL","500mL"=>"500 mL"))",
            xlabel="kₗa (h⁻¹)",ylabel="Bulk oxygen demand / capacity",xscale=log10,yscale=log10,
            xticks=([0.1,0.2,0.5,1.,2.,5.,10.],["0.1","0.2","0.5","1","2","5","10"]),
            titlesize=13,xlabelsize=12,ylabelsize=12,xticklabelsize=10,yticklabelsize=10)
        for (mult,col) in zip((.5,1.,2.),("#4E79A7","#006B7F","#AA4E35"))
            q=sort(filter(x->x.source_condition==s.source_condition && x.q_multiplier==mult,rows);by=x->x.kla_h)
            lines!(oxygen_ax,[x.kla_h for x in q],[x.Da for x in q];color=col,linewidth=2)
        end
        hlines!(oxygen_ax,[1.];color=:black,linestyle=:dash)
        ylims!(oxygen_ax,.001,100);xlims!(oxygen_ax,.1,10)
    end
    Legend(f[3,1:3],[LineElement(color=c,linewidth=2) for c in ("#4E79A7","#006B7F","#AA4E35")],
        ["0.5 × reference uptake","Reference uptake: 17.7 amol cell⁻¹ s⁻¹","2 × reference uptake"];orientation=:horizontal,labelsize=11,tellwidth=false)
    export_figure(joinpath(dir,"figures"),"oxygen_sensitivity",f)
    title="Conditional bulk oxygen sensitivity"
    note="The dashed line is demand/capacity = 1. Uptake multipliers and kLa ranges are exploratory assumptions. This does not establish actual dissolved oxygen or exclude intra-aggregate hypoxia. Consult assumptions.toml and correlation_audit.csv before citing a numerical kLa."
    links="<a href='condition_summary.csv'>Six-condition table</a> · <a href='correlation_audit.csv'>Correlation audit</a> · <a href='assumptions.toml'>Assumptions and source URLs</a>"
    table_rows=join(["<tr><td>$(Char(96+i))</td><td>$(replace(s.label,"100mL"=>"100 mL","500mL"=>"500 mL"))</td><td>$(s.rpm)</td><td>$(round(s.conditional_kla_h;digits=3))</td><td>$(round(s.conditional_Da;digits=3))</td><td>$(round(s.required_kla_h_at_q_reference;digits=3))</td><td>$(s.kla_status)</td></tr>" for (i,s) in enumerate(summary)],"\n")
    table_html="""<h2>Six-condition oxygen assessment</h2><table><thead><tr><th>Panel</th><th>Condition</th><th>RPM</th><th>Conditional k<sub>L</sub>a (h⁻¹)</th><th>Demand / capacity</th><th>Required k<sub>L</sub>a at ratio 1 (h⁻¹)</th><th>Basis and assumptions</th></tr></thead><tbody>$table_rows</tbody></table><p>Ratios use the peak condition-day mean cell density, reference uptake 17.7 amol cell⁻¹ s⁻¹, and assumed saturation 0.21 mmol/L. The 100 mL transfer coefficients use Neto’s 41.5-mm wheel, source contact area and temperature-dependent water properties at 37 °C; the 500 mL value is an assumed input, not a measurement. Uptake and saturation were not measured in these cultures. Full inputs and individual-vessel peak checks are retained in the linked CSV and assumptions file.</p>"""
    if isfile(joinpath(dir,"temperature_comparison.csv"))
        comparisons=collect(CSV.File(joinpath(dir,"temperature_comparison.csv")))
        comparison_rows=join(["<tr><td>$(r.temperature_C)</td><td>$(r.volume_ml)</td><td>$(round(r.kla_h;digits=3))</td><td>$(r.label)</td></tr>" for r in comparisons],"\n")
        table_html*="""<h2>Temperature and volume check at 60 rpm</h2><table><tr><th>Temperature (°C)</th><th>Volume (mL)</th><th>Calculated k<sub>L</sub>a (h⁻¹)</th><th>Interpretation</th></tr>$comparison_rows</table><p>All entries are correlation calculations, not digitized measurements. The 37 °C calculation extends the experimental temperature range using Neto’s water-property functions. Source files and their pinned repository version are recorded in sources.csv and data/reference/neto2022/provenance.json.</p>"""
    end
    image="oxygen_sensitivity"
elseif kind=="ablation"
    include("AblationPlots.jl")
    ablation_plot(dir)
    title="Shear and EDR ablation sensitivity"
    note="This all-data sensitivity retains the primary Model 6 formulation. Each reduced variant has 10 active fitted parameters; baseline has 9 and full Model 6 has 11. Baseline/full fits are reused; both reduced variants are independently refitted. 500-mL outcomes are evaluated only after training selection. Convergence and boundaries remain in ranking.csv; no independent causal effect is established."
    links="<a href='ranking.csv'>Training BIC and convergence</a> · <a href='state_errors.csv'>Statewise and overall MRPE / Experimental Noise</a> · <a href='runs.csv'>All new starts</a>"
    image="ablation_comparison"
else
    error("Unknown analysis")
end
write(joinpath(dir,"index.html"),"""<!doctype html><meta charset="utf-8"><title>$title</title><style>body{font:16px Arial;max-width:1150px;margin:40px auto;line-height:1.5}img{width:100%}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:10px;border:1px solid #ddd;text-align:left}th{background:#edf2f5}</style><h1>$title</h1><p>$note</p><img src="figures/$image.svg"><p>$links</p>$table_html""")
