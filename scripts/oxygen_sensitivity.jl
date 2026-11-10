using CSV, TOML, Statistics, SHA, Test
const ROOT=normpath(joinpath(@__DIR__,".."))
length(ARGS)==1 || error("Usage: oxygen_sensitivity.jl OUTPUT")
out=abspath(only(ARGS));mkpath(out)
configfile=joinpath(ROOT,"config/oxygen_sensitivity.toml");cfg=TOML.parsefile(configfile)
series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
measurements=collect(CSV.File(joinpath(ROOT,"data/processed/measurements.csv")))
q0=cfg["q_reference_amol_cell_s"]*1e-18;cs=cfg["cstar_mol_L"]
# X is in cells/mL, q in mol/cell/s and kLa in h^-1.
capacity(q,x,kla)=q*x*1000/(kla/3600*cs)
required_kla(q,x)=q*x*1000/cs*3600
reference_dir=joinpath(ROOT,"data/reference/neto2022")
diffusion=collect(CSV.File(joinpath(reference_dir,"2PDabPlot.txt");delim='\t'))
temperatures=Float64[r.T for r in diffusion]
diffusivities=Float64[r.Dab*1e-4 for r in diffusion]
slope=sum((temperatures.-mean(temperatures)).*(diffusivities.-mean(diffusivities)))/sum(abs2,temperatures.-mean(temperatures))
intercept=mean(diffusivities)-slope*mean(temperatures)

# Transcribed from Neto's Constants.m, calcRhoH2O.m and calcNuH2O.m.
function water_properties(t,pressure)
    rho0=evalpoly(t,(999.84847,6.337563e-2,-8.523829e-3,6.943248e-5,-3.821216e-7))
    kt=evalpoly(t,(50.83101e-8,-3.68293e-9,7.263725e-11,-6.597702e-13,2.87767e-15))
    rho=rho0*(1+kt*(pressure*100-101.325))
    dt=20-t
    nu=1.002*10^(dt/(96+t)*evalpoly(dt,(1.2364,-1.37e-3,5.7e-6,2.55e-8)))/rho*1e-3
    dab=intercept+slope*t
    (;rho,nu,mu=rho*nu,dab)
end
function neto(rpm,cfg;temperature=cfg["temperature_C"],volume=1e-4,area=cfg["area_100_m2"])
    d=cfg["wheel_diameter_assumed_m"];a=area/volume
    props=water_properties(temperature,cfg["pressure_bar"])
    re=(rpm/60)*d^2/props.nu
    sc=props.nu/props.dab
    g=d*a
    sh=1.68*re^0.551*sc^(1/3)*g^1.18
    kla=sh*props.dab/d*a*3600
    (;kla,re,sc,g,props)
end
@testset "Oxygen dimensional and threshold checks" begin
    @test required_kla(q0,1e6)≈0.30342857142857144
    @test capacity(q0,1e6,required_kla(q0,1e6))≈1
    @test capacity(2q0,1e6,1.)≈2capacity(q0,1e6,1.)
    @test capacity(q0,1e6,2.)≈capacity(q0,1e6,1.)/2
    @test cfg["area_100_m2"]/1e-4≈16.37976
    @test neto(40,cfg).g≈.67976004
    @test length(diffusion)==21
    @test water_properties(20.,1.01325).mu≈.001002
    @test water_properties(37.,1.01325).dab≈3.428867894794901e-9
    @test neto(60,cfg;temperature=37.).kla>neto(60,cfg;temperature=24.2).kla
end
peaks=NamedTuple[];grid=NamedTuple[];audit=NamedTuple[]
for c in 1:6
    info=first(filter(r->r.condition==c,series))
    rows=filter(r->r.condition==c && r.state=="cells" && r.status=="measured",measurements)
    daily=[(day=d,value=mean(Float64[r.value for r in rows if r.day==d])) for d in sort(unique(r.day for r in rows))]
    best=argmax(r->r.value,daily);peak=best.value;individual=maximum(r.value for r in rows)
    z=neto(info.rpm,cfg)
    nominal=info.volume_ml==100 ? z.kla : cfg["kla_500_assumed_h"]
    role=info.volume_ml==500 ? "500mL evaluation" : "100mL condition $(c<4 ? c : c-1)"
    push!(peaks,(source_condition=c,label=role,rpm=info.rpm,volume_ml=info.volume_ml,
        peak_day=best.day,peak_condition_mean_cells_ml=peak,max_individual_cells_ml=individual,
        required_kla_h_at_q_reference=required_kla(q0,peak),required_kla_h_individual=required_kla(q0,individual),
        conditional_kla_h=nominal,conditional_Da=capacity(q0,peak,nominal),
        kla_status=info.volume_ml==100 ? "Neto correlation, verified wheel and area, original water-property functions at 37 C" : "Unvalidated inherited 1.8/h assumption"))
    for mult in cfg["q_multipliers"],kla in exp.(range(log(cfg["kla_grid_min_h"]),log(cfg["kla_grid_max_h"]);length=cfg["kla_grid_points"]))
        push!(grid,(source_condition=c,q_multiplier=mult,q_amol_cell_s=mult*cfg["q_reference_amol_cell_s"],kla_h=kla,Da=capacity(q0*mult,peak,kla)))
    end
    if info.volume_ml==100
        d=cfg["wheel_diameter_assumed_m"]
        push!(audit,(source_condition=c,Re=z.re,Sc=z.sc,G=z.g,in_original_Re_range=624<=z.re<=3172,
            in_original_Sc_range=366<=z.sc<=395,in_original_G_range=.651<=z.g<=1.08,
            temperature_C=cfg["temperature_C"],rho_kg_m3=z.props.rho,nu_m2_s=z.props.nu,mu_Pa_s=z.props.mu,diffusivity_m2_s=z.props.dab,conditional_kla_h=z.kla))
    end
end
sort!(peaks;by=r->(r.volume_ml==500,r.source_condition))
checks=[(label=label,temperature_C=t,volume_ml=v*1e6,rpm=60,kla_h=neto(60,cfg;temperature=t,volume=v,area=a).kla) for (label,t,v,a) in (("Neto 60-rpm experimental settings: correlation estimate, not measured value",24.2,1.05e-4,.001647072),("100-mL room-temperature comparison",24.2,1e-4,cfg["area_100_m2"]),("100-mL culture-temperature estimate",37.,1e-4,cfg["area_100_m2"]))]
CSV.write(joinpath(out,"temperature_comparison.csv"),checks)
CSV.write(joinpath(out,"condition_summary.csv"),peaks)
CSV.write(joinpath(out,"sensitivity.csv"),grid);CSV.write(joinpath(out,"correlation_audit.csv"),audit)
cp(configfile,joinpath(out,"assumptions.toml");force=true)
files=vcat([configfile,@__FILE__,joinpath(ROOT,"data/processed/series.csv"),joinpath(ROOT,"data/processed/measurements.csv")],readdir(reference_dir;join=true))
CSV.write(joinpath(out,"sources.csv"),[(path=relpath(f,ROOT),sha256=bytes2hex(sha256(read(f)))) for f in files])
write(joinpath(out,"interpretation.txt"),"Primary output is required kLa at Da=1, conditional on oxygen uptake and saturation. Da=OUR/(kLa*Cstar) describes bulk transfer capacity at maximal driving force, not actual dissolved oxygen or aggregate-core oxygen. Peak density is the maximum condition-day mean; the maximum individual value is supplied separately. The 100-mL correlation uses Neto's verified 41.5-mm wheel, exact source contact area and original temperature-dependent water-property functions. The 37 C Schmidt value is outside the experimental fitting range; the paper supports temperature extension numerically. Estimates remain conditional on water-to-medium and headspace transferability. The 500-mL kLa is assumed. Sensitivity ranges are exploratory, not empirical uncertainty bounds. No hypoxia exclusion or oxygen-replete classification is made. Neto's measurements used uncapped water vessels; applicability to culture medium/headspace conditions is limited.\n")
write(joinpath(out,"COMPLETE"),"Conditional oxygen sensitivity calculated\n")
println("Oxygen tables saved: ",out)
