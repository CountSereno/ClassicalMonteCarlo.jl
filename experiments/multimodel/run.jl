# Sampling writes evidence once; analysis and plotting never restart a sampler.
using ClassicalMonteCarlo, TOML, SHA, Dates, Random

function generation_hashes()
    repository = normpath(joinpath(@__DIR__, "../.."))
    paths = ["Project.toml", "experiments/multimodel/run.jl"]
    for (directory, _, files) in walkdir(joinpath(repository, "src"))
        append!(paths, [relpath(joinpath(directory, file), repository)
                        for file in files if endswith(file, ".jl")])
    end
    Dict(path => bytes2hex(sha256(read(joinpath(repository, path)))) for path in paths)
end

function experiment_cases(config)
    cases = NamedTuple[]
    q = config["potts_q"]
    beta = config["exact_beta"]
    p = -expm1(-beta)
    for algorithm in config["potts_algorithms"], (chain, initial) in enumerate(config["potts_initials"])
        push!(cases, (; family="Potts", experiment="exact", L=config["exact_L"],
                      q=Float64(q), beta, p, algorithm, chain, initial))
    end
    for rc_q in config["rc_q_values"]
        algorithms = rc_q == 1 ? ["sweeny", "independent"] : ["sweeny"]
        for algorithm in algorithms, (chain, initial) in enumerate(config["rc_initials"])
            push!(cases, (; family="RandomCluster", experiment="exact", L=config["exact_L"],
                          q=Float64(rc_q), beta, p, algorithm, chain, initial))
        end
    end
    beta = log1p(sqrt(q))
    p = -expm1(-beta)
    for L in config["critical_sizes"], (family, algorithm, initials) in
        (("Potts", "swendsen_wang", config["potts_initials"]),
         ("RandomCluster", "sweeny", config["rc_initials"]))
        for (chain, initial) in enumerate(initials)
            push!(cases, (; family, experiment="critical", L, q=Float64(q), beta, p,
                          algorithm, chain, initial))
        end
    end
    cases
end

function specification(case, config, seed; burn=config["burn"], samples=config["samples"])
    common = (; L=case.L, algorithm=Symbol(case.algorithm), initial=Symbol(case.initial),
               seed, burn, samples, stride=config["stride"])
    if case.family == "Potts"
        PottsRunSpec(; common..., q=Int(case.q), beta=case.beta, J=1)
    else
        RandomClusterRunSpec(; common..., p=case.p, q=case.q, periodic=true)
    end
end

same_configuration(a::Potts, b::Potts) = a.spins == b.spins
same_configuration(a::BondConfiguration, b::BondConfiguration) =
    a.occupied == b.occupied && a.graph.nvertices == b.graph.nvertices && a.graph.edges == b.graph.edges

function audit_restart(original, directory, run_id)
    partial = start_run(original.spec)
    midpoint = original.spec.burn + (original.spec.samples ÷ 2) * original.spec.stride
    advance!(partial, midpoint)
    checkpoint = joinpath(directory, "checkpoints", run_id * ".checkpoint")
    save_checkpoint(checkpoint, partial)
    restored = load_checkpoint(checkpoint; expected_spec=original.spec)
    run!(restored)
    # The generic observations contract avoids assumptions about spin/bond fields.
    equal = observations(restored) == observations(original) &&
            same_configuration(restored.model, original.model) &&
            restored.updates == original.updates &&
            restored.measurement_updates == original.measurement_updates &&
            restored.kernel_work == original.kernel_work &&
            restored.pending_work == original.pending_work &&
            rand(copy(restored.rng), UInt64, 16) == rand(copy(original.rng), UInt64, 16)
    equal || error("checkpoint changed the trajectory or future random stream")
    Dict("run_id" => run_id, "exactly_equal" => equal,
         "checkpoint_update" => midpoint, "final_update" => restored.updates,
         "scope" => "paired observables, final physical configuration, clock, work counters and future RNG output; same environment")
end

function run_experiments(config_path, output)
    config = TOML.parsefile(config_path)
    ispath(output) && (!isdir(output) || !isempty(readdir(output))) &&
        error("output must be a new or empty directory: $output")
    all(b -> b > 0 && config["samples"] % b == 0 && config["samples"] ÷ b >= 2,
        config["block_sizes"]) || error("block sizes must divide samples into at least two blocks")
    length(config["potts_initials"]) == length(config["rc_initials"]) >= 2 ||
        error("use equal numbers of independent chains in both representations")
    hashes = generation_hashes()
    cases = experiment_cases(config)
    # Compile all kernel/measurement paths on throwaway states, outside timing.
    for index in unique(i -> (cases[i].family, cases[i].algorithm), eachindex(cases))
        warm = start_run(specification(cases[index], config, 17; burn=4, samples=8))
        advance!(warm, warm.spec.burn)
        run!(warm)
    end
    mkpath(output)
    cp(config_path, joinpath(output, "config.toml"))
    manifest = Dict{String,Any}(
        "schema_version" => 1, "started_utc" => string(now(UTC)),
        "julia_version" => string(VERSION), "machine" => Sys.MACHINE,
        "cpu_name" => Sys.CPU_NAME, "threads" => Threads.nthreads(),
        "logical_cpu_threads" => Sys.CPU_THREADS, "source_sha256" => hashes,
        "timing" => "production updates plus all recorded observables; excludes compile, burn, I/O and analysis; explicit GC before each timer",
        "boundary" => "periodic square", "J" => 1, "kB" => 1,
        "runs" => Any[], "restart_audits" => Any[])
    audited = Set{Tuple{String,String}}()
    for (index, case) in enumerate(cases)
        run_id = "run_" * lpad(index, 3, '0')
        seed = config["seed"] + 1009 * index
        state = start_run(specification(case, config, seed))
        advance!(state, state.spec.burn)
        GC.gc()
        seconds = @elapsed run!(state)
        export_run(joinpath(output, "runs", run_id), state)
        row = Dict{String,Any}(string(key) => value for (key, value) in pairs(case))
        merge!(row, Dict("run_id" => run_id, "seed" => seed, "samples" => config["samples"],
            "burn" => config["burn"], "stride" => config["stride"],
            "nvertices" => case.L^2, "nedges" => 2 * case.L^2,
            "production_seconds" => seconds,
            "update_unit" => update_unit(state.kernel), "work_unit" => work_unit(state.kernel)))
        push!(manifest["runs"], row)
        key = (case.family, case.algorithm)
        if !(key in audited)
            push!(manifest["restart_audits"], audit_restart(state, output, run_id))
            push!(audited, key)
        end
        println("[$index/$(length(cases))] $run_id $(case.family) $(case.experiment) L=$(case.L) q=$(case.q) $(case.algorithm): $(round(seconds; digits=3)) s")
    end
    hashes == generation_hashes() || error("source changed during production; do not use this dataset")
    manifest["finished_utc"] = string(now(UTC))
    manifest["measurements"] = sum(row["samples"] for row in manifest["runs"])
    open(io -> TOML.print(io, manifest), joinpath(output, "manifest.toml"), "w")
    println("Saved $(manifest["measurements"]) measurements to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    config = isempty(ARGS) ? joinpath(@__DIR__, "configs", "demo.toml") : abspath(ARGS[1])
    output = length(ARGS) < 2 ? joinpath(@__DIR__, "results", "demo") : abspath(ARGS[2])
    run_experiments(config, output)
end
