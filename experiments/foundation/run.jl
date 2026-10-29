# Generate raw evidence only. Analysis and plotting are separate, repeatable jobs.
using ClassicalMonteCarlo, TOML, SHA, Dates

function source_hashes()
    root = normpath(joinpath(@__DIR__, "../.."))
    paths = ["Project.toml",
             joinpath("experiments","foundation","run.jl"),
             joinpath("experiments","foundation","exact.jl")]
    for (directory, _, files) in walkdir(joinpath(root,"src"))
        append!(paths, [relpath(joinpath(directory,f),root)
                        for f in files if endswith(f,".jl")])
    end
    Dict(p=>bytes2hex(sha256(read(joinpath(root,p)))) for p in paths)
end

function run_experiments(config_path, output)
    cfg = TOML.parsefile(config_path)
    ispath(output) && (!isdir(output) || !isempty(readdir(output))) &&
        error("output must be a new or empty directory: $output")
    all(b -> b > 0 && cfg["samples"] % b == 0 && cfg["samples"]÷b >= 2,
        cfg["block_sizes"]) || error("all block sizes must divide samples into at least two blocks")
    length(cfg["initial_states"]) >= 2 || error("multiple independent chains required")
    mkpath(output)
    cp(config_path,joinpath(output,"config.toml"))
    hashes = source_hashes()
    cases = NamedTuple[]
    for algorithm in cfg["algorithms"], (chain,initial) in enumerate(cfg["initial_states"])
        push!(cases,(experiment="exact",L=cfg["exact_L"],T=cfg["exact_T"],algorithm,chain,initial))
    end
    for L in cfg["sizes"], algorithm in cfg["algorithms"], (chain,initial) in enumerate(cfg["initial_states"])
        push!(cases,(experiment="critical",L,T=2/log1p(sqrt(2)),algorithm,chain,initial))
    end
    # Compile all selected kernel/run paths on throwaway states before timing.
    for algorithm in cfg["algorithms"]
        warm = start_run(RunSpec(L=3,beta=1/3,algorithm=Symbol(algorithm),burn=4,samples=8,seed=77))
        advance!(warm,4)
        run!(warm)
    end
    manifest = Dict{String,Any}("schema_version"=>1,"started_utc"=>string(now(UTC)),
        "julia_version"=>string(VERSION),"machine"=>Sys.MACHINE,"cpu_name"=>Sys.CPU_NAME,
        "threads"=>Threads.nthreads(),"logical_cpu_threads"=>Sys.CPU_THREADS,
        "source_sha256"=>hashes,"runs"=>Any[],
        "timing"=>"production kernel updates plus E/M recording; excludes compilation, warmup, storage and analysis; explicit GC before each measurement",
        "clock"=>"one metropolis_sweep, one wolff_cluster, or one swendsen_wang_sweep; do not equate these costs")
    first_state = nothing
    for (id,c) in enumerate(cases)
        runid = "run_"*lpad(id,3,'0')
        spec = RunSpec(L=c.L,beta=1/c.T,algorithm=Symbol(c.algorithm),initial=Symbol(c.initial),
            seed=cfg["seed"]+1009 * id,burn=cfg["burn"],samples=cfg["samples"],stride=cfg["stride"])
        state = start_run(spec)
        advance!(state,spec.burn)
        # Measurement functions are compiled on the throwaway states above.
        GC.gc()
        seconds = @elapsed run!(state)
        export_run(joinpath(output,"runs",runid),state)
        record = Dict{String,Any}(string(k)=>v for (k,v) in pairs(c))
        merge!(record,Dict("run_id"=>runid,"seed"=>spec.seed,"samples"=>spec.samples,
            "burn"=>spec.burn,"stride"=>spec.stride,"production_seconds"=>seconds,
            "update_unit"=>string(update_unit(state.kernel)),"work_unit"=>string(work_unit(state.kernel))))
        push!(manifest["runs"],record)
        id == 1 && (first_state = state)
        println("[$id/$(length(cases))] $runid $(c.experiment) L=$(c.L) $(c.algorithm): $(round(seconds;digits=3)) s")
    end
    # A visible integration audit: restore a partially completed production run
    # and compare all recorded values, counters and final configuration.
    original = first_state
    partial = start_run(original.spec)
    advance!(partial,original.spec.burn+(original.spec.samples÷2)*original.spec.stride)
    checkpoint = joinpath(output,"restart.checkpoint")
    save_checkpoint(checkpoint,partial)
    resumed = load_checkpoint(checkpoint;expected_spec=original.spec)
    run!(resumed)
    matched = resumed.energies == original.energies && resumed.magnetizations == original.magnetizations &&
        resumed.measurement_updates == original.measurement_updates && resumed.kernel_work == original.kernel_work &&
        resumed.model.spins == original.model.spins && resumed.updates == original.updates
    matched || error("continuous and resumed runs disagree")
    manifest["restart_audit"] = Dict("run_id"=>"run_001","exactly_equal"=>matched,
        "checkpoint_update"=>partial.updates,"final_update"=>resumed.updates,
        "scope"=>"same Julia/code/config environment; arrays, counters, final spins")
    hashes == source_hashes() || error("source changed during generation; rerun from stable code")
    manifest["finished_utc"] = string(now(UTC))
    manifest["measurements"] = sum(r["samples"] for r in manifest["runs"])
    open(io->TOML.print(io,manifest),joinpath(output,"manifest.toml"),"w")
    println("Saved $(manifest["measurements"]) measurements to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    config = length(ARGS)>=1 ? abspath(ARGS[1]) : joinpath(@__DIR__,"configs","demo.toml")
    output = length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(@__DIR__,"results","demo")
    run_experiments(config,output)
end
