#= Recompute K_Lac and Y_GX profiles with all ten nuisance parameters free.

Independent grid points own their objective, starts and optimizer state. Only
the new output directory is written; the primary fit and nine-profile run are
read-only. Saved ProfileLikelihood grids are reused and extended if necessary.
=#
include("run_fixed_two.jl")
using Dates, Statistics

const DIAGNOSTIC_NAMES = ("K_Lac", "Y_GX")
const DIAGNOSTIC_THRESHOLD = ProfileLikelihood.get_chisq_threshold(.95)

function diagnostic_setup(fitdir, seed_dir)
    validated = fixed_two_context(fitdir)
    ctx, model = validated.raw, validated.model
    names = collect(free_names(ctx, model))
    isempty(ctx.config["fixed"]) && length(names) == 11 || error("Diagnostics require all eleven parameters free")
    seedcfg = TOML.parsefile(joinpath(seed_dir, "calibration.toml"))
    seedctx = context(cfg=seedcfg)
    seedmeta = TOML.parsefile(joinpath(seed_dir, "settings.toml"))
    seedmeta["model"] == 5 && seedmeta["profile_scope"] == "all-active" || error("Expected original all-active Model 6 seeds")
    isempty(seedcfg["fixed"]) || error("Seeds must not fix a parameter")
    seedcfg["bounds"] == ctx.config["bounds"] || error("Seed bounds differ")
    seedctx.observations == ctx.observations && seedctx.y == ctx.y && seedctx.sigma == ctx.sigma || error("Seed observations or weights differ")
    # Reuse the already converged all-active profiling reference. The original
    # calibration selected the same optimum but its final polish hit its budget;
    # do not invent a successful retcode for that earlier optimizer invocation.
    reference_summary = only(CSV.File(joinpath(seed_dir,"mle_summary.csv")))
    reference_summary.model==5 && reference_summary.free_parameters==11 && reference_summary.stopping_met===true && reference_summary.smoke_test===false || error("Expected a converged all-active profile reference")
    p = read_parameters(joinpath(seed_dir,"mle_parameters.csv"),model)
    z = clamp.(encode(ctx, model, p), 0., 1.)
    nll = gaussian_nll(ctx.y, prediction(ctx, model, p), ctx.sigma)
    isapprox(nll,reference_summary.nll;atol=1e-6,rtol=0) || error("Saved all-active reference likelihood differs")
    isapprox(nll,validated.saved_nll;atol=1e-5,rtol=0) || error("Diagnostic reference differs materially from the primary optimum")
    (; ctx, model, p, names, z, nll, saved_nll=validated.saved_nll, reference_retcode=String(reference_summary.retcode))
end

function diagnostic_loglik(ctx, model, z)
    all(v -> isfinite(v) && 0 <= v <= 1, z) || return -1e100
    try
        -gaussian_nll(ctx.y, prediction(ctx, model, decode(ctx, model, z)), ctx.sigma)
    catch e
        e isa InterruptException && rethrow()
        -1e100
    end
end

# Coordinates separated only by rounding at a bound/reference do not provide
# distinct profile information. Do not send those duplicate points to optimizers.
function diagnostic_grid(coordinates, mle)
    canonical(x) = abs(x) < 1e-9 ? 0. : abs(x-1) < 1e-9 ? 1. : abs(x-mle) < 1e-9 ? mle : x
    sort(unique(canonical.(vcat(Float64.(coordinates), [Float64(mle)]))))
end

point_id(z) = "point_" * string(reinterpret(UInt64, Float64(z)); base=16, pad=16)
point_directory(out, name, z) = joinpath(out, name, "points", point_id(z))
point_complete(out, name, z) = isfile(joinpath(point_directory(out, name, z), "COMPLETE"))

function interpolated_seed(bank, coordinate)
    left = findlast(v -> v < coordinate-1e-12, bank.coordinates)
    right = findfirst(v -> v > coordinate+1e-12, bank.coordinates)
    (left === nothing || right === nothing) && return nothing
    weight = (coordinate-bank.coordinates[left])/(bank.coordinates[right]-bank.coordinates[left])
    (1-weight).*bank.vectors[left] .+ weight.*bank.vectors[right]
end

function diagnostic_starts(bank, coordinate, extra=Vector{Float64}[])
    nearest = argmin(abs.(bank.coordinates.-coordinate))
    candidates = [copy(bank.vectors[nearest]), seed_vectors(bank, coordinate)..., deepcopy(extra)...]
    interpolated = interpolated_seed(bank, coordinate)
    interpolated === nothing || push!(candidates, interpolated)
    # A saved nuisance optimum can lie almost exactly on a bound. Retain those
    # starts and also start slightly inside the box so BOBYQA can explore with
    # a useful initial radius. The optimization bounds themselves do not change.
    append!(candidates, [clamp.(s, 1e-3, 1-1e-3) for s in deepcopy(candidates)])
    starts = Vector{Float64}[]
    for start in candidates
        all(isfinite, start) || error("Nonfinite starting vector")
        start = clamp.(start, 0., 1.)
        any(s -> maximum(abs.(s.-start)) < 1e-12, starts) || push!(starts, start)
    end
    starts
end

function optimize_diagnostic_point(loglik, reference, ref_nll, index, coordinate, starts; maxiters=10000)
    others = setdiff(eachindex(reference), [index])
    function objective(u, _)
        z = copy(reference)
        z[index] = coordinate
        z[others] .= u
        -loglik(z)-ref_nll
    end
    f = PackageWorkflow.Optimization.OptimizationFunction(objective, PackageWorkflow.Optimization.AutoFiniteDiff())
    prob = PackageWorkflow.Optimization.OptimizationProblem(f, copy(first(starts)); lb=zeros(length(others)), ub=ones(length(others)))
    # This closure and optimizer belong to this invocation, not an outer loop.
    alg = AuditedBOBYQA(NamedTuple[], ReentrantLock(), _ -> deepcopy(starts[2:end]))
    sol = SciMLBase.solve(prob, alg; maxiters, abstol=1e-8)
    z = copy(reference); z[index] = coordinate; z[others] .= sol.u
    actual = loglik(z)+ref_nll
    isfinite(actual) && isapprox(actual, -Float64(sol.objective); atol=1e-6, rtol=0) || error("Point objective replay differs")
    actual <= 1e-5 || error("Found a better reference likelihood; do not publish this scan without review")
    (; z, relative_loglik=actual, stopping_met=SciMLBase.successful_retcode(sol.retcode), retcode=string(sol.retcode), records=alg.records)
end

function run_diagnostic_point(prep, job, bank, out; maxiters=10000)
    name, coordinate = job.name, job.coordinate
    dest = point_directory(out, name, coordinate)
    point_complete(out, name, coordinate) && return false
    index = something(findfirst(==(Symbol(name)), prep.names))
    starts = diagnostic_starts(bank, coordinate, job.extra)
    loglik = z -> diagnostic_loglik(prep.ctx, prep.model, z)
    result = optimize_diagnostic_point(loglik, prep.z, prep.nll, index, coordinate, starts; maxiters)
    # Do not write COMPLETE until all numerical files are safely in place.
    mkpath(dest)
    p = decode(prep.ctx, prep.model, result.z)
    save_parameters(joinpath(dest, "parameters.csv"), p)
    row = (parameter=name, parameter_value=getproperty(p, Symbol(name)), normalized_parameter=coordinate,
           profile_loglik=result.relative_loglik, stopping_met=result.stopping_met, retcode=result.retcode,
           starts=length(starts), best_seed=only(result.records).best_seed)
    CSV.write(joinpath(dest, "point.csv"), [row])
    CSV.write(joinpath(dest, "optimizer_calls.csv"), [(objective=r.objective, stopping_met=r.stopping_met, retcode=r.retcode, starts=r.starts, best_seed=r.best_seed, start_objective_spread=r.start_objective_spread) for r in result.records])
    write(joinpath(dest, "COMPLETE"), "Point saved. Review stopping status and curve audit before inference.\n")
    true
end

function run_diagnostic_batch(prep, jobs, banks, out; maxiters=10000)
    pending = filter(j -> !point_complete(out, j.name, j.coordinate), jobs)
    isempty(pending) && return
    counter = Threads.Atomic{Int}(0)
    print_lock = ReentrantLock()
    errors = String[]
    println("Optimizing ", length(pending), " independent profile points with up to ", Threads.nthreads(), " threads"); flush(stdout)
    Threads.@threads :greedy for k in eachindex(pending)
        job = pending[k]
        try
            run_diagnostic_point(prep, job, banks[job.name], out; maxiters)
            finished = Threads.atomic_add!(counter, 1)+1
            lock(print_lock) do
                println("Completed profile point ", job.name, " ", finished, "/", length(pending), " in this batch"); flush(stdout)
            end
        catch e
            e isa InterruptException && rethrow()
            message = sprint(showerror, e, catch_backtrace())
            dest = point_directory(out, job.name, job.coordinate); mkpath(dest)
            write(joinpath(dest, "FAILED.txt"), message*"\n")
            lock(print_lock) do
                push!(errors, job.name*": "*message)
                println(stderr, "FAILED profile point ", job.name, ": ", sprint(showerror, e)); flush(stderr)
            end
        end
    end
    isempty(errors) || error("$(length(errors)) point(s) failed; completed point checkpoints are preserved. See FAILED.txt files.")
end

function completed_points(out, name)
    parent = joinpath(out, name, "points")
    isdir(parent) || return NamedTuple[]
    rows = [merge(NamedTuple(only(CSV.File(joinpath(parent, d, "point.csv")))), (directory=joinpath(parent, d),))
            for d in readdir(parent) if isfile(joinpath(parent, d, "COMPLETE"))]
    sort(rows; by=r -> r.normalized_parameter)
end

function extension_jobs(prep, out, name)
    rows = completed_points(out, name)
    index = something(findfirst(==(Symbol(name)), prep.names))
    mle = prep.z[index]
    result = NamedTuple[]
    for side in (:lower, :upper)
        branch = filter(r -> side == :lower ? r.normalized_parameter < mle : r.normalized_parameter > mle, rows)
        any(r -> r.profile_loglik <= DIAGNOSTIC_THRESHOLD, branch) && continue
        endpoint = side == :lower ? first(rows) : last(rows)
        bound = side == :lower ? 0. : 1.
        abs(endpoint.normalized_parameter-bound) < 1e-9 && continue
        zvals = [r.normalized_parameter for r in rows]
        gaps = filter(v -> v > 1e-8, side == :lower ? diff(zvals) : reverse(diff(zvals)))
        step = isempty(gaps) ? .01 : first(gaps)
        p = read_parameters(joinpath(endpoint.directory, "parameters.csv"), prep.model)
        warm = encode(prep.ctx, prep.model, p)[setdiff(eachindex(prep.names), [index])]
        for j in 1:4
            z = clamp(endpoint.normalized_parameter+(side == :lower ? -1 : 1)*j*step, 0., 1.)
            any(t -> t.coordinate == z, result) || push!(result, (name=name, coordinate=z, extra=[copy(warm)]))
            z == bound && break
        end
    end
    result
end

function diagnostic_limits(coordinates, values, mle)
    center = argmin(abs.(coordinates.-mle))
    splines = Dict{Int, Any}()
    ProfileLikelihood.spline_profile!(splines, 1, coordinates, values, ProfileLikelihood.FritschCarlsonMonotonicInterpolation, ProfileLikelihood.Line)
    function endpoint(order)
        for (a,b) in zip(order[1:end-1], order[2:end])
            if values[a] >= DIAGNOSTIC_THRESHOLD && values[b] <= DIAGNOSTIC_THRESHOLD
                lo,hi = minmax(coordinates[a], coordinates[b])
                flo = splines[1](lo)-DIAGNOSTIC_THRESHOLD
                for _ in 1:60
                    mid = (lo+hi)/2
                    if (splines[1](mid)-DIAGNOSTIC_THRESHOLD >= 0) == (flo >= 0); lo=mid; else; hi=mid; end
                end
                return ((lo+hi)/2, true)
            end
        end
        (coordinates[last(order)], false)
    end
    lower,lc = endpoint(collect(center:-1:1))
    upper,uc = endpoint(collect(center:length(coordinates)))
    (; lower, upper, lower_crossed=lc, upper_crossed=uc)
end

function diagnostic_curve_issues(values, stops)
    issues = String[]
    all(isfinite, values) || push!(issues, "nonfinite likelihood")
    maximum(values) <= 1e-5 || push!(issues, "better likelihood than reference")
    all(stops) || push!(issues, "nuisance optimization requires review")
    inside = values .>= DIAGNOSTIC_THRESHOLD
    count(i -> inside[i] && (i==1 || !inside[i-1]), eachindex(inside)) <= 1 || push!(issues, "disconnected confidence region")
    dips = count(i -> min(values[i-1],values[i+1])-values[i] > .2 && max(values[i-1],values[i+1]) > -3, 2:length(values)-1)
    dips == 0 || push!(issues, "$(dips) sharp dips near confidence region")
    issues
end

function export_diagnostic(prep, out, name)
    rows = completed_points(out, name)
    isempty(rows) && error("No completed points for $name")
    nameindex = something(findfirst(==(Symbol(name)), prep.names))
    profile = NamedTuple[]; nuisance = NamedTuple[]
    for (j,r) in enumerate(rows)
        p = read_parameters(joinpath(r.directory,"parameters.csv"),prep.model)
        actual = prep.nll-gaussian_nll(prep.ctx.y,prediction(prep.ctx,prep.model,p),prep.ctx.sigma)
        isapprox(actual,r.profile_loglik;atol=1e-5,rtol=0) || error("Saved point replay differs: $name point $j")
        push!(profile,(parameter=name,parameter_value=r.parameter_value,profile_loglik=actual,normalized_parameter=r.normalized_parameter,stopping_met=r.stopping_met,retcode=r.retcode))
        append!(nuisance,[(point=j,parameter=string(n),value=getproperty(p,n)) for n in prep.names])
    end
    coordinates=Float64[r.normalized_parameter for r in profile]; values=Float64[r.profile_loglik for r in profile]
    limits=diagnostic_limits(coordinates,values,prep.z[nameindex])
    lo,hi=prep.ctx.config["bounds"][name]
    physical(z)=exp(log(lo)+z*log(hi/lo))
    issues=diagnostic_curve_issues(values,[r.stopping_met===true for r in profile])
    qualifications=String[]
    limits.lower_crossed || push!(qualifications, "lower limit unresolved within tested bounds")
    limits.upper_crossed || push!(qualifications, "upper limit unresolved within tested bounds")
    interval=(parameter=name,mle=getproperty(prep.p,Symbol(name)),package_lower=physical(limits.lower),package_upper=physical(limits.upper),lower_hits_bound=!limits.lower_crossed,upper_hits_bound=!limits.upper_crossed,unconfirmed_points=count(r->r.stopping_met!==true,profile),better_minimum=maximum(values)>1e-5)
    audit=(parameter=name,points=length(profile),lower_crossed=limits.lower_crossed,upper_crossed=limits.upper_crossed,numerical_review_required=!isempty(issues),status=join(vcat(issues,qualifications)," | "),lower=limits.lower_crossed ? physical(limits.lower) : NaN,upper=limits.upper_crossed ? physical(limits.upper) : NaN)
    dest=joinpath(out,name)
    CSV.write(joinpath(dest,"profile.csv"),profile)
    CSV.write(joinpath(dest,"nuisance_parameters.csv"),nuisance)
    CSV.write(joinpath(dest,"interval.csv"),[interval])
    write(joinpath(dest,"COMPLETE"),"Numerical curve exported. Consult audit.csv before use.\n")
    println("FINISHED diagnostic ",name,": ", isempty(audit.status) ? "thresholds crossed on both sides" : audit.status);flush(stdout)
    (;interval,audit)
end

function fixed_parameter_diagnostics_main(args=ARGS)
    length(args)==3 || error("Usage: run_fixed_parameter_diagnostics.jl REFERENCE_FIT ALL_ACTIVE_SEEDS NEW_OUTPUT")
    fit,seed,out=abspath.(args)
    safe_output(out,[fit,seed,joinpath(ROOT,"results/training60_final/profile_fixed_two_intervals_recovered"),joinpath(ROOT,"results/training60_final/profile_fixed_two")])
    prep=diagnostic_setup(fit,seed)
    maxiters=parse(Int,get(ENV,"HPSC_DIAGNOSTIC_MAXITERS","10000"))
    maxiters>0 || error("Invalid optimizer budget")
    smoke=get(ENV,"HPSC_DIAGNOSTIC_SMOKE","0")=="1"
    settings=Dict("engine"=>"ProfileLikelihood.jl","model"=>5,"publication_model"=>6,"parameters"=>collect(DIAGNOSTIC_NAMES),"free_parameters"=>11,"fixed_parameters"=>Dict{String,Float64}(),"training_cultures"=>19,"selection_cultures_used"=>0,"source_fit"=>fit,"seed_directory"=>seed,"maxiters"=>maxiters,"smoke_test"=>smoke,"implementation"=>"Independent grid-point nuisance optimization; package likelihood-ratio threshold and monotone interpolation","grid_policy"=>"Saved grid, merging coordinates within 1e-9 of reference or bounds; extend by original edge spacing until threshold or bound","primary_fit_reestimated"=>false)
    files=[joinpath(ROOT,f) for f in ("profiles/run_fixed_parameter_diagnostics.jl","profiles/run_fixed_two.jl","profiles/PackageWorkflow.jl","profiles/FixedTwoSeeds.jl","profiles/Project.toml","profiles/Manifest.toml","src/TargetedChecks.jl","src/Workflow.jl","src/Calibration.jl","src/Proliferation.jl","data/processed/measurements.csv","data/processed/series.csv","data/processed/feeding.csv","data/processed/lagrangian_exposure.csv")]
    append!(files,[joinpath(fit,f) for f in ("calibration.toml","split.csv","observations.csv","model5/selected_parameters.csv","model5/summary.csv")])
    append!(files,[joinpath(seed,f) for f in ("calibration.toml","settings.toml","mle_parameters.csv","mle_summary.csv")])
    for name in DIAGNOSTIC_NAMES, file in ("profile.csv","nuisance_parameters.csv");push!(files,joinpath(seed,name,file));end
    hashes=[(path=f,sha256=bytes2hex(sha256(read(f)))) for f in files]
    if isdir(out)
        TOML.parsefile(joinpath(out,"settings.toml"))==settings || error("Resume settings differ")
        old=collect(CSV.File(joinpath(out,"source_hashes.csv")))
        length(old)==length(hashes) && all(a.path==b.path && a.sha256==b.sha256 for (a,b) in zip(old,hashes)) || error("Resume sources differ")
        if isfile(joinpath(out,"COMPLETE"));println("Already complete: $out");return;end
    else
        mkpath(out)
        open(io->TOML.print(io,settings),joinpath(out,"settings.toml"),"w")
        CSV.write(joinpath(out,"source_hashes.csv"),hashes)
        open(io->TOML.print(io,prep.ctx.config),joinpath(out,"calibration.toml"),"w")
        save_parameters(joinpath(out,"mle_parameters.csv"),prep.p)
        CSV.write(joinpath(out,"mle_summary.csv"),[(model=5,saved_nll=prep.saved_nll,initial_nll=prep.nll,nll=prep.nll,free_parameters=11,retcode=prep.reference_retcode,stopping_met=true,smoke_test=smoke)])
        write(joinpath(out,"fixed_parameters.csv"),"parameter,value,source\n")
    end
    banks=Dict{String,FixedTwoSeeds.SeedBank}()
    jobs=NamedTuple[]
    for name in DIAGNOSTIC_NAMES
        index=something(findfirst(==(Symbol(name)),prep.names))
        bank=load_seed_bank(seed,prep.names,index,prep.ctx.config["bounds"],prep.z)
        banks[name]=bank
        grid=diagnostic_grid(bank.coordinates,prep.z[index])
        if smoke
            # Include a known formerly underoptimized Y_GX point in a real-data test.
            grid=name=="Y_GX" ? [bank.coordinates[108]] : [first(grid)]
        end
        append!(jobs,[(name=name,coordinate=z,extra=Vector{Float64}[]) for z in grid])
    end
    CSV.write(joinpath(out,"initial_grid.csv"),[(parameter=j.name,normalized_parameter=j.coordinate) for j in jobs])
    open(io->println(io,Dates.now()," START independent-point diagnostics"),joinpath(out,"progress.log"),"a")
    run_diagnostic_batch(prep,jobs,banks,out;maxiters)
    if smoke
        write(joinpath(out,"SMOKE_COMPLETE"),"Real-data test points only; not publication curves.\n")
        println("Real-data smoke points finished: $out");return
    end
    for batch in 1:250
        additional=vcat([extension_jobs(prep,out,name) for name in DIAGNOSTIC_NAMES]...)
        isempty(additional) && break
        batch<250 || error("Profile extension did not resolve; inspect output")
        println("Extending profile grid: ",length(additional)," points");flush(stdout)
        run_diagnostic_batch(prep,additional,banks,out;maxiters)
    end
    exports=[export_diagnostic(prep,out,name) for name in DIAGNOSTIC_NAMES]
    CSV.write(joinpath(out,"profile_confidence_intervals.csv"),[e.interval for e in exports])
    CSV.write(joinpath(out,"audit.csv"),[e.audit for e in exports])
    write(joinpath(out,"COMPLETE"),"Numerical diagnostics complete; inspect audit.csv. Unresolved confidence limits are reported, not forced to cross.\n")
    if any(e->e.audit.numerical_review_required,exports)
        write(joinpath(out,"NUMERICAL_REVIEW_REQUIRED"),"Inspect audit.csv before using these curves.\n")
    end
    open(io->println(io,Dates.now()," FINISHED diagnostics and audit"),joinpath(out,"progress.log"),"a")
    println("Saved two diagnostic profiles: $out");flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
    fixed_parameter_diagnostics_main()
end
