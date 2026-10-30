# Shared file/provenance utilities only; no sampling side effects on include.
using TOML, SHA
const REPOSITORY = normpath(joinpath(@__DIR__,"../.."))
file_hash(path)=bytes2hex(sha256(read(path)))
function source_hashes(paths)
    Dict(path=>file_hash(joinpath(REPOSITORY,path)) for path in paths)
end
function generation_sources()
    paths=["Project.toml","experiments/scaling/run.jl","experiments/scaling/common.jl"]
    append!(paths,[relpath(joinpath(d,f),REPOSITORY) for (d,_,fs) in walkdir(joinpath(REPOSITORY,"src")) for f in fs if endswith(f,".jl")])
    source_hashes(paths)
end
save_toml(path,data)=open(io->TOML.print(io,data),path,"w")

function critical_cases(config)
    cases=NamedTuple[]
    for q in config["q_values"], L in config["sizes"], (chain,initial) in enumerate(config["rc_initials"])
        push!(cases,(;family="RandomCluster",algorithm="sweeny",L,q=Float64(q),chain,initial))
    end
    for L in config["spin_sizes"], (chain,initial) in enumerate(config["potts_initials"])
        push!(cases,(;family="Potts",algorithm="swendsen_wang",L,q=Float64(config["spin_q"]),chain,initial))
    end
    cases
end
