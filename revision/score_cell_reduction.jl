include("October20.jl")
dir,dest=ARGS
cfg=TOML.parsefile(joinpath(dir,"calibration.toml"));line=strip(read(joinpath(dir,"cell_line.txt"),String));ctx=cell_context(line,cfg)
p=read_parameters(joinpath(dir,"reported_parameters.csv"),6)
CSV.write(dest,[score_row(ctx,6,p,line,"cell_line_fit")])
