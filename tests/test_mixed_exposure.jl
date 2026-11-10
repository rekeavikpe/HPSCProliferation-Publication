using Test, CSV, TOML, Statistics
include("../revision/MixedExposure.jl")
using .MixedExposure
using .MixedExposure.RevisionCommon, .RevisionCommon.Workflow, .RevisionCommon.Workflow.Calibration
using .RevisionCommon.Workflow.Calibration.Proliferation

@testset "Independent P50 shear and P75 EDR routing" begin
    rows=collect(CSV.File(joinpath(ROOT,"data/processed/exposure_summary_inputs.csv")))
    for (v,rpm) in ((100,30),(100,40),(100,60),(500,30))
        grid(p)=sort(filter(r->r.scenario=="lag_p$p"&&r.volume_ml==v&&r.rpm==rpm,rows);by=r->r.size_um)
        a,b=grid(50),grid(75);lookup=mixed_lookup(v,rpm)
        for i in eachindex(a)
            @test lookup(a[i].size_um)==(a[i].shear_pa,b[i].edr_m2_s3)
        end
        @test collect(lookup(300)) ≈ [(a[1].shear_pa+a[2].shear_pa)/2,(b[1].edr_m2_s3+b[2].edr_m2_s3)/2]
        @test lookup(15)==lookup(200)
        @test lookup(1200)==lookup(1000)
    end
    @test collect(mixed_lookup(100,50)(300)) ≈ (collect(mixed_lookup(100,40)(300)).+collect(mixed_lookup(100,60)(300)))./2
    @test_throws ErrorException mixed_lookup(100,100)
end

if !isempty(ARGS)
    length(ARGS)==2 || error("Supply REFERENCE_FIT SPLIT_FIT")
    prep,p=source_contexts(ARGS...)
    @testset "Preserved split, objective and evaluation data" begin
        @test length(prep.split.training)==19
        @test length(prep.split.selection)==12
        @test isempty(intersect(prep.split.training,prep.split.selection))
        for (_,base,settings) in evaluations(prep)
            same=replace_exposure(base;edr_percentile=50,settings)
            mixed=replace_exposure(base;settings)
            @test same.y==mixed.y==base.y
            @test same.sigma==mixed.sigma==base.sigma
            @test mixed.observations==base.observations
            @test mixed.config==base.config
            for (c,q) in base.conditions
                @test mixed.conditions[c].u0==q.u0
                @test mixed.conditions[c].feeds==q.feeds
                @test mixed.conditions[c].times==q.times
                for size in (15.,200.,300.,800.,1200.)
                    @test collect(same.conditions[c].lookup(size)) ≈ collect(q.lookup(size)) rtol=1e-12
                    @test mixed.conditions[c].lookup(size)[1] ≈ q.lookup(size)[1] rtol=1e-12
                end
            end
            @test prediction(same,6,p) ≈ prediction(base,6,p) rtol=1e-8
        end
        @test length(free_names(replace_exposure(prep.train),6))==11
        nll=gaussian_nll(prep.held.y,prediction(prep.held,6,p),prep.held.sigma)
        @test round(bic(nll,11,length(prep.held.y));digits=2)==1458.73
        for (cohort,ctx,_) in evaluations(prep)
            cohort in ("training60","selection40","validation500","media20","media50") || continue
            target=Dict("training60"=>14.2,"selection40"=>17.4,"validation500"=>16.0,"media20"=>17.9,"media50"=>15.3)[cohort]
            actual=only(r.mrpe_percent for r in errors(ctx,prediction(ctx,6,p)) if r.state=="all")
            @test round(actual;digits=1)==target
        end
    end
end
