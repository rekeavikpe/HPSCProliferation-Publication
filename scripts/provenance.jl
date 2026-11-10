using SHA, CSV
root=normpath(joinpath(@__DIR__,".."))
files=String[]
for directory in ("src","scripts","test","config","data/processed","reference","figures","structural","profiles","profile_figures","tools")
    for (folder,_,names) in walkdir(joinpath(root,directory))
        occursin(joinpath("figures","output"),folder) && continue
        append!(files,[joinpath(folder,n) for n in names if any(ext->endswith(n,ext),(".jl",".toml",".csv",".json",".py",".txt"))])
    end
end
append!(files,[joinpath(root,"Project.toml"),joinpath(root,"Manifest.toml")])
CSV.write(joinpath(root,"reference/file_hashes.csv"),[(path=relpath(f,root),sha256=bytes2hex(sha256(read(f)))) for f in sort(files) if f!=joinpath(root,"reference/file_hashes.csv")])
