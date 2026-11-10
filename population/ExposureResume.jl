# Only these audited earlier drivers are compatible with this resume implementation.
const COMPATIBLE_EXPOSURE_DRIVERS=Set([
    "abe8f03b47a58f7102e1786e9585ccbef9ea54e4f48ad772664dd8ccfc97bfd6",
    "d02711ffdcc7df7130fac21bf003e35d69b5eb46cea3a03bd03d2fd20ae7f54e",
])
function verify_resume_inputs(out,settings,files)
    previous=TOML.parsefile(joinpath(out,"analysis.toml"))
    for key in ("model","legacy_model","n","k","scenarios","starts_per_summary","population_size","max_evaluations","local_iterations","seed_policy","method","inputs")
        previous[key]==settings[key] || error("Resume settings differ: $key")
    end
    previous["julia"]==string(VERSION) || error("Resume with original Julia version $(previous["julia"])")
    saved=collect(CSV.File(joinpath(out,"source_hashes.csv")))
    for r in saved
        path=normpath(joinpath(ROOT,r.path));isfile(path) || error("Resume input missing: $(r.path)")
        current=bytes2hex(sha256(read(path)))
        compatible=path==joinpath(@__DIR__,"exposure_sensitivity.jl") && r.sha256 in COMPATIBLE_EXPOSURE_DRIVERS
        (current==r.sha256 || compatible) || error("Resume input changed: $(r.path). Existing fits preserved.")
    end
    for f in files
        any(r->normpath(joinpath(ROOT,r.path))==normpath(f),saved) || error("Resume source not present in original manifest: $f")
    end
    bytes2hex(sha256(read(joinpath(out,"calibration.toml"))))==bytes2hex(sha256(read(joinpath(settings["source_fit"],"calibration.toml")))) || error("Saved calibration differs")
    bytes2hex(sha256(read(joinpath(out,"exposure_inputs.csv"))))==bytes2hex(sha256(read(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))) || error("Saved exposure inputs differ")
end
function verify_saved_fit(dir,job,ctx)
    try
        isfile(joinpath(dir,"FAILED")) && return (false,"previous failure; preserve and retry")
        r=only(CSV.File(joinpath(dir,"summary.csv")))
        r.scenario==job.scenario && r.run==job.run || error("Search identity mismatch")
        r.stopping_met isa Bool || error("Missing optimizer status")
        p=read_parameters(joinpath(dir,"selected_parameters.csv"),6)
        z=encode(ctx,6,p)
        all(x->isfinite(x)&&-1e-8<=x<=1+1e-8,z) || error("Parameters outside original bounds")
        nll=gaussian_nll(ctx.y,prediction(ctx,6,p),ctx.sigma)
        isfinite(nll) && isapprox(nll,r.nll;atol=1e-6,rtol=0) || error("Saved parameters do not reproduce likelihood")
        isapprox(bic(nll,11,116),r.bic;atol=2e-6,rtol=0) || error("Saved BIC mismatch")
        stages=collect(CSV.File(joinpath(dir,"stages.csv")))
        any(s->s.stage==r.stage && s.stopping_met==r.stopping_met && isapprox(s.nll,r.nll;atol=1e-6,rtol=0),stages) || error("Saved stage mismatch")
        (true,"verified saved likelihood, BIC, bounds and optimizer status")
    catch e
        e isa InterruptException && rethrow()
        (false,sprint(showerror,e))
    end
end
