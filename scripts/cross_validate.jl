using CSV, Statistics
include("../src/Workflow.jl")
include("../src/ReplicateSplits.jl")
using .Workflow, .ReplicateSplits
length(ARGS)==2 || error("Usage: cross_validate.jl SPLIT_NUMBER OUTPUT_DIRECTORY")
repetition=parse(Int,ARGS[1]); repetition>0 || error("Split number must be positive")
out=abspath(ARGS[2]); ispath(out) && error("Output already exists")
series=collect(CSV.File(joinpath(ROOT,"data/processed/series.csv")))
split=random_split(series;seed=20261003+repetition-1)
train=context(ids=split.training)
run_calibration(out;ctx=train,seed=20261003+1000repetition)
CSV.write(joinpath(out,"split.csv"),split.assignments)
# Selection outcomes enter only after all training fits are saved.
held=context(ids=split.selection,fallback_glucose=Dict(c=>q.glucose for (c,q) in train.conditions))
results=NamedTuple[]
for id in 0:6
    m=legacy_model(id)
    p=read_parameters(joinpath(out,"model$id/selected_parameters.csv"),m)
    pred=prediction(held,m,p)
    for state in STATES
        inds=findall(r->r.state==state,held.observations)
        isempty(inds) && continue
        y=held.y[inds];yp=pred[inds]
        nonzero=findall(!iszero,y)
        # Scale errors using calibration-set SDs/floors, never selection-set SDs.
        scales=[train.sigma[only(findall(r->r.condition==held.observations[i].condition && r.day==held.observations[i].day && r.state==state,train.observations))] for i in inds]
        push!(results,(model=id,state=state,n=length(y),mrpe_percent=100mean(abs.((yp[nonzero].-y[nonzero])./y[nonzero])),rmse=sqrt(mean((yp.-y).^2)),bias=mean(yp.-y),scaled_mse=mean(((yp.-y)./scales).^2)))
    end
end
CSV.write(joinpath(out,"selection_errors.csv"),results)
println("Selection errors saved; they are predictive scores, not held-out BIC.")
