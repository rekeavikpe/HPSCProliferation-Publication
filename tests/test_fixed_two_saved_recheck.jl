using Test
include("../profiles/run_fixed_two.jl")

@testset "Recovered profiles keep valid values and flag improved nuisance fits" begin
    names=[:x,:y];bounds=Dict("x"=>[1.,exp(1.)],"y"=>[1.,exp(1.)])
    loglik(z,_)=-.5*(((z[1]-.5)/.05)^2+((z[2]-.4-.7*(z[1]-.5))/.1)^2)
    problem=make_problem(loglik,[.5,.4],names)
    sol=fit_mle(problem;maxiters=100)
    for offset in (0.,.02)
        mktempdir() do dir
            dest=joinpath(dir,"x");mkpath(dest)
            xs=collect(range(.35,.65;length=21))
            vectors=[[x,.4+.7*(x-.5)+offset] for x in xs]
            rows=[(parameter="x",parameter_value=exp(z[1]),normalized_parameter=z[1],profile_loglik=loglik(z,nothing)-sol.maximum,stopping_met=true,retcode="Success") for z in vectors]
            CSV.write(joinpath(dest,"profile.csv"),rows)
            CSV.write(joinpath(dest,"nuisance_parameters.csv"),[(point=j,parameter=string(names[k]),value=exp(z[k])) for (j,z) in enumerate(vectors) for k in eachindex(names)])
            write(joinpath(dest,"interval.csv"),"saved interval")
            write(joinpath(dest,"optimizer_calls.csv"),"saved calls")
            write(joinpath(dest,"COMPLETE"),"complete")
            write(joinpath(dest,"RECOVERED_NEEDS_CHECK"),"check")
            before=sha256(read(joinpath(dest,"profile.csv")))
            bank=load_seed_bank(dir,names,1,bounds,sol.mle)
            accepted=recheck_saved_profile(dest,names,1,bounds,sol,loglik,bank;maxiters=500)
            @test accepted==(offset==0.)
            @test sha256(read(joinpath(dest,"profile.csv")))==before
            @test isfile(joinpath(dest,"COMPLETE"))==(offset==0.)
            @test !isfile(joinpath(dest,"RECOVERED_NEEDS_CHECK"))
            @test isfile(joinpath(dest,"RECOMPUTE_REQUIRED"))==(offset!=0.)
            if offset!=0.
                @test sha256(read(joinpath(dest,"before_seed_recheck/profile.csv")))==before
            else
                @test isfile(joinpath(dest,"RECOVERY_CHECKED"))
            end
            checks=collect(CSV.File(joinpath(dest,"recovery_recheck.csv")))
            @test length(checks)==5
            @test all(r->r.stopping_met,checks)
        end
    end
end
