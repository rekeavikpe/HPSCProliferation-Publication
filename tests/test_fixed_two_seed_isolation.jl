using Test, CSV
include("../profiles/FixedTwoSeeds.jl")
using .FixedTwoSeeds

@testset "Parallel parameter and direction jobs keep separate seed banks" begin
    mktempdir() do dir
        names=[:a,:b,:c]
        bounds=Dict(string(n)=>[1.,exp(1.)] for n in names)
        reference=[.2,.3,.4]
        counts=[10,20,7]
        for i in eachindex(names)
            sub=joinpath(dir,string(names[i]));mkpath(sub)
            CSV.write(joinpath(sub,"profile.csv"),[(normalized_parameter=(j-1)/(counts[i]-1),) for j in 1:counts[i]])
            CSV.write(joinpath(sub,"nuisance_parameters.csv"),[(point=j,parameter=string(n),value=exp(.1i+.001j+.0001k)) for j in 1:counts[i] for (k,n) in enumerate(names)])
        end
        banks=Vector{FixedTwoSeeds.SeedBank}(undef,3)
        providers=Vector{Any}(undef,3)
        Threads.@threads for i in eachindex(names)
            banks[i]=load_seed_bank(dir,names,i,bounds,reference)
            providers[i]=seed_provider(banks[i])
        end
        @test all(p->all(n->!(getfield(p,n) isa Core.Box),fieldnames(typeof(p))),providers)
        checks=fill(false,900)
        Threads.@threads for task in eachindex(checks)
            i=mod1(task,3);coordinate=isodd(task) ? 0. : 1.
            # The package supplies its fixed coordinate through this interface.
            problem=(f=(f=(θₙ=coordinate,),),)
            yield()
            obtained=providers[i](problem)
            js=coordinate==0 ? [1,2,3] : [counts[i],counts[i]-1,counts[i]-2]
            others=setdiff(eachindex(names),[i])
            expected=[copy(reference[others]),[[.1i+.001j+.0001k for k in others] for j in js]...]
            checks[task]=all(isapprox.(reduce(vcat,obtained),reduce(vcat,expected);atol=1e-12))
            # A solve must not modify the saved reference or another job's seeds.
            obtained[1].=999.;obtained[end].=-999.
        end
        @test all(checks)
        @test reference==[.2,.3,.4]
        for i in eachindex(names)
            @test seed_vectors(banks[i],0.)[1]==reference[setdiff(eachindex(names),[i])]
            @test length(seed_vectors(banks[i],.5))==4
        end
        @test_throws ErrorException seed_vectors(banks[1],NaN)
        @test_throws ErrorException providers[1]((f=(f=(wrong=.5,),),))
    end
end

@testset "Regression reproduces the old captured-variable failure" begin
    function old_shared_capture()
        old=nothing
        fixed=nothing
        providers=Function[]
        for i in 1:2
            old=collect(1:(10i))
            vectors=Dict(j=>j for j in eachindex(old))
            push!(providers,function(coordinate)
                fixed=coordinate
                nearest=sortperm(abs.(old.-fixed))[1:3]
                [vectors[j] for j in nearest]
            end)
        end
        providers[1](20.)
    end
    # Even sequential invocation of the two escaped closures is enough to
    # reproduce key 20 not found. Thread scheduling made its timing variable.
    @test_throws KeyError old_shared_capture()
end
