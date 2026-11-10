module StructuralModels
using StructuralIdentifiability, CSV
export build_model, model_source, exposure_coefficients, ROOT
const ROOT=normpath(joinpath(@__DIR__,".."))
const LEGACY=(0,1,3,4,5,6,7)
function decimal_rational(text)
    s=lowercase(string(text)); parts=split(s,'e'); exponent=length(parts)==2 ? parse(Int,parts[2]) : 0
    mantissa=parts[1]; dot=findfirst(==('.'),mantissa)
    decimals=dot===nothing ? 0 : length(mantissa)-dot
    numerator=parse(BigInt,replace(mantissa,"."=>""));power=exponent-decimals
    power>=0 ? numerator*big(10)^power//big(1) : numerator//big(10)^(-power)
end
function exposure_coefficients(rpm,segment)
    rows=sort(filter(r->r.volume_ml=="100"&&r.rpm==string(rpm)&&r.percentile=="50",collect(CSV.File(joinpath(ROOT,"data/processed/lagrangian_exposure.csv");types=String)));by=r->parse(Float64,r.size_um))
    isempty(rows) && error("No CFD data for $rpm RPM")
    sizes=decimal_rational.([r.size_um for r in rows]); tau=decimal_rational.([r.shear_pa for r in rows]); edr=decimal_rational.([r.edr_m2_s3 for r in rows])
    if segment==0
        return (tau_intercept=tau[1],tau_slope=0//big(1),edr_intercept=edr[1],edr_slope=0//big(1))
    elseif segment==5
        return (tau_intercept=tau[end],tau_slope=0//big(1),edr_intercept=edr[end],edr_slope=0//big(1))
    end
    segment in 1:4 || error("Segments: 0 lower clamp; 1–4 interpolation; 5 upper clamp")
    i=segment;st=(tau[i+1]-tau[i])/(sizes[i+1]-sizes[i]);se=(edr[i+1]-edr[i])/(sizes[i+1]-sizes[i])
    (tau_intercept=tau[i]-st*sizes[i],tau_slope=st,edr_intercept=edr[i]-se*sizes[i],edr_slope=se)
end
ratstr(x)="($(numerator(x))//$(denominator(x)))"
function model_source(public_id;design=:pooled,segment=0)
    public_id in 0:6 || error("Public model ID must be 0–6")
    model=LEGACY[public_id+1]
    rpms=design==:pooled ? [40,40,40,30,60] : design==:single40 ? [40] : error("Unknown design")
    equations=String["E'(t)=-E(t)//12"]
    for (c,rpm) in enumerate(rpms)
        X,G,L,A=["$(s)$(c)(t)" for s in ("X","G","L","A")]
        q=exposure_coefficients(rpm,segment)
        tau="($(ratstr(q.tau_intercept))+$(ratstr(q.tau_slope))*$A)"
        edr="($(ratstr(q.edr_intercept))+$(ratstr(q.edr_slope))*$A)"
        base="(mu_max*$G/($G+K_G)*K_Lac/($L+K_Lac)*K_Agg/($A^3+K_Agg))"
        shear=if model==0
            "1"
        elseif model==3
            "(1/(1+$tau/K_tau+($tau/K_tau)^2/2))"
        elseif model==5
            "(K_tau^2/($tau^2+K_tau^2))"
        elseif model==7
            effective="(K_tau*(1+gamma*$edr/($edr+K_eps)))"
            "($effective/($tau+$effective))"
        else
            "(K_tau/($tau+K_tau))"
        end
        growth="($base*$shear)";metabolic=model==4 ? base : growth
        lactate="(K_oG*$metabolic*$G)"
        aggregate=model==0 ? "alpha*$base" : model==6 ? "$growth*(alpha-$edr/K_eps*$A)" : "alpha*$metabolic-delta_max*$edr/($edr+K_eps)*$A"
        append!(equations,["X$(c)'(t)=($growth-K_d*E(t))*$X","G$(c)'(t)=-(Y_GX*$metabolic+Y_GLac*$lactate)*$X","L$(c)'(t)=$lactate*$X","A$(c)'(t)=$aggregate"])
    end
    push!(equations,"clock(t)=E(t)")
    for c in eachindex(rpms),s in ("X","G","L","A")
        push!(equations,"obs$(s)$(c)(t)=$(s)$(c)(t)")
    end
    "@ODEmodel(\n"*join(equations,",\n")*"\n)"
end
function build_model(public_id;design=:pooled,segment=0)
    source=model_source(public_id;design,segment)
    Core.eval(@__MODULE__,Meta.parse(source))
end
end
