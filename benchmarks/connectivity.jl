# Paired timing repeats: both backends start from the same warmed configuration
# and random stream, and must finish with identical occupations and future RNG.
using ClassicalMonteCarlo, Random, Statistics, TOML, SHA, Dates

function sweep_batch!(rng,state,p,q,workspace,sweeps)
    accepted=0
    for _ in 1:sweeps
        accepted+=sweeny_sweep!(rng,state,p,q;workspace)
    end
    accepted
end

function benchmark(output)
    ispath(output) && error("choose a new benchmark output directory")
    mkpath(output)
    write(joinpath(output,"timings.csv"),"q,L,backend,repeat,sweeps,seconds_per_sweep\n")
    rows=NamedTuple[]
    for q in (2.,3.), L in (8,16,32,64)
        p=sqrt(q)/(1+sqrt(q))
        graph=square_bond_graph(L;periodic=true)
        equilibrium=BondConfiguration(graph)
        rng=Xoshiro(20261031+100L+Int(q))
        sweep_batch!(rng,equilibrium,p,q,SweenyWorkspace(graph),2048)
        sweeps=max(4,4096÷L^2)
        times=Dict("union_find_rebuild"=>Float64[],"interleaved_bfs"=>Float64[])
        allocations=Dict{String,Int}()
        for repeat in 1:5
            states=[BondConfiguration(graph,equilibrium.occupied) for _ in 1:2]
            streams=[copy(rng),copy(rng)]
            workspaces=(BondWorkspace(graph),SweenyWorkspace(graph))
            names=("union_find_rebuild","interleaved_bfs")
            for i in (isodd(repeat) ? (1,2) : (2,1))
                # Warm the specialization on separate occupations and RNG.
                sweep_batch!(copy(rng),BondConfiguration(graph,equilibrium.occupied),p,q,workspaces[i],2)
                GC.gc()
                elapsed=@elapsed sweep_batch!(streams[i],states[i],p,q,workspaces[i],sweeps)
                push!(times[names[i]],elapsed/sweeps)
                allocations[names[i]]=@allocated sweep_batch!(copy(rng),states[i],p,q,workspaces[i],0)
            end
            states[1].occupied==states[2].occupied || error("backend trajectory mismatch")
            rand(streams[1],UInt64,16)==rand(streams[2],UInt64,16) || error("backend RNG mismatch")
        end
        # Dedicated allocation check excludes constructing/copying the RNG.
        for (name,w) in (("union_find_rebuild",BondWorkspace(graph)),("interleaved_bfs",SweenyWorkspace(graph)))
            r=copy(rng); s=BondConfiguration(graph,equilibrium.occupied)
            sweep_batch!(r,s,p,q,w,1)
            allocations[name]=@allocated sweep_batch!(r,s,p,q,w,1)
        end
        ref=median(times["union_find_rebuild"]); fast=median(times["interleaved_bfs"])
        push!(rows,(;q,L,sweeps,repeats=5,reference_seconds=ref,ibfs_seconds=fast,
                    speedup=ref/fast,reference_bytes=allocations["union_find_rebuild"],ibfs_bytes=allocations["interleaved_bfs"]))
        println(last(rows)); flush(stdout)
        for name in keys(times)
            open(joinpath(output,"timings.csv"),"a") do io
                for (repeat,seconds) in enumerate(times[name])
                    println(io,join((q,L,name,repeat,sweeps,seconds),','))
                end
            end
        end
    end
    open(joinpath(output,"summary.csv"),"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row->println(io,join(values(row),',')),rows)
    end
    root=dirname(@__DIR__)
    paths=vcat(["Project.toml","benchmarks/connectivity.jl"],
        [relpath(joinpath(d,f),root) for (d,_,fs) in walkdir(joinpath(root,"src")) for f in fs if endswith(f,".jl")])
    metadata=Dict("julia_version"=>string(VERSION),"cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,
        "timestamp_utc"=>string(now(UTC)),"scope"=>"update-only; lazy skip in both; no measurement, compile or setup; medians of five alternating-order repeats",
        "paired_trajectories_and_rng_equal"=>true,"source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(root,p)))) for p in paths))
    open(io->TOML.print(io,metadata),joinpath(output,"metadata.toml"),"w")
end

if abspath(PROGRAM_FILE)==@__FILE__
    benchmark(isempty(ARGS) ? joinpath(@__DIR__,"results","connectivity") : abspath(ARGS[1]))
end
