#=
Reproducible Ising algorithm comparison and finite-size diagnostic scan.

Run from the repository root:
    julia --project=. experiments/run_ising.jl experiments/configs/demo.toml
Only Julia standard libraries and the parent package are required.
=#
using ClassicalMonteCarlo, Statistics, TOML, SHA, Dates, Random

const TC = 2 / log(1 + sqrt(2))
const FIELDS = (:energy_density, :abs_magnetization_density, :heat_capacity,
                :susceptibility, :susceptibility_abs, :binder)

function update!(model, beta, algorithm, rng)
    algorithm == "metropolis" && return metropolis!(model, beta; rng=rng)
    algorithm == "wolff" && return wolff!(model, beta; rng=rng)
    error("unknown algorithm: $algorithm")
end

function diagnostic(x, maxlag)
    if var(x) == 0
        return (tau=NaN, window=0, stderr_mean=NaN)
    end
    integrated_autocorrelation_time(x; maxlag=maxlag)
end

function one_chain(L, T, algorithm, initial, seed, cfg)
    rng = new_rng(seed)
    model = Ising(L; init=initial, rng=rng)
    beta = 1 / T
    for _ in 1:cfg["burn"]
        update!(model, beta, algorithm, rng)
    end
    n = cfg["samples"]
    E, M, work = zeros(n), zeros(n), zeros(Int, n)
    # Compile the selected update and both measurement methods before timing.
    # These read-only calls do not change the RNG stream or state.
    energy(model)
    magnetization(model)
    elapsed = @elapsed for i in 1:n
        work[i] = update!(model, beta, algorithm, rng)
        E[i] = energy(model)
        M[i] = magnetization(model)
    end
    obs = ising_observables(E, M; beta=beta, nsites=L^2, blocksize=cfg["blocksize"])
    dE = diagnostic(E, cfg["maxlag"])
    dM = diagnostic(abs.(M), cfg["maxlag"])
    # A second, longer block length helps reveal unstable error estimates.
    larger_block = 4 * cfg["blocksize"]
    obs4 = ising_observables(E, M; beta=beta, nsites=L^2, blocksize=larger_block)
    (E=E, M=M, work=work, elapsed=elapsed, obs=obs, obs4=obs4, dE=dE, dM=dM)
end

function main()
    config_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "configs", "demo.toml")
    cfg = TOML.parsefile(config_path)
    out = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(@__DIR__, "results", cfg["name"])
    ispath(out) && !isempty(readdir(out)) && error("output directory is not empty; choose a new output path")
    cfg["samples"] % (4 * cfg["blocksize"]) == 0 || error("samples must be divisible by 4 * blocksize")
    cfg["samples"] >= 8 * cfg["blocksize"] || error("at least two large blocks required")
    1 <= cfg["maxlag"] < cfg["samples"] || error("maxlag must lie in [1, samples-1]")
    cfg["chains"] >= 2 || error("at least two independent chains are required")
    mkpath(joinpath(out, "traces"))
    cp(config_path, joinpath(out, "config.toml"))
    src = joinpath(@__DIR__, "..", "src")
    hashes = Dict(basename(p) => bytes2hex(sha256(read(p))) for p in sort(readdir(src; join=true)) if endswith(p, ".jl"))
    hashes["experiments/run_ising.jl"] = bytes2hex(sha256(read(@__FILE__)))
    metadata = Dict("config"=>cfg, "julia_version"=>string(VERSION),
        "started_utc"=>string(now(UTC)), "platform"=>string(Sys.MACHINE),
        "threads"=>Threads.nthreads(), "source_sha256"=>hashes,
        "units"=>"J = kB = 1; total E and signed total M in traces; periodic square lattice",
        "measurement_interval"=>"one Metropolis sweep or one Wolff cluster, depending on algorithm",
        "runtime"=>"production updates plus E/M measurement; excludes compilation, burn-in, I/O, analysis",
        "work"=>"Metropolis accepted flips; Wolff flipped cluster size; neither equals CPU work",
        "tau_method"=>"positive-window rho sum, stops at first nonpositive lag; cannot capture negative tails",
        "initial_states"=>"alternating random/up across independent chains")
    cases = NamedTuple[]
    for L in cfg["sizes"], algorithm in ("metropolis", "wolff"), chain in 1:cfg["chains"]
        push!(cases, (experiment="benchmark", L=L, T=TC, algorithm=algorithm, chain=chain))
    end
    for L in cfg["sizes"], T in cfg["temperatures"], chain in 1:cfg["chains"]
        push!(cases, (experiment="scan", L=L, T=T, algorithm="wolff", chain=chain))
    end
    basefields = ["run_id", "experiment", "L", "T", "algorithm", "chain", "initial", "seed",
        "samples", "burn", "blocksize", "production_seconds", "mean_flipped_sites", "tau_E", "tau_absM",
        "window_E", "window_absM", "ESS_E", "ESS_absM", "ESS_E_per_second", "ESS_absM_per_second",
        "tau_E_spin_sweeps", "tau_absM_spin_sweeps", "window_truncated", "block_warning"]
    statfields = reduce(vcat, [[string(f), string(f)*"_stderr", string(f)*"_stderr_large_block"] for f in FIELDS])
    open(joinpath(out, "summary.csv"), "w") do summary
        println(summary, join(vcat(basefields, statfields), ','))
        for (i, c) in enumerate(cases)
            seed = cfg["seed"] + i * 1009
            initial = isodd(c.chain) ? :random : :up
            runid = "run_" * lpad(i, 3, '0')
            result = one_chain(c.L, c.T, c.algorithm, initial, seed, cfg)
            E, M, work = result.E, result.M, result.work
            n = length(E)
            dE, dM = result.dE, result.dM
            # Wolff work-normalized time uses average flipped spins / N;
            # it is an approximate convention, not an elapsed CPU-time measure.
            scale = c.algorithm == "wolff" ? mean(work) / c.L^2 : 1.0
            essE, essM = n/(2 * dE.tau), n/(2 * dM.tau)
            truncated = dE.window == cfg["maxlag"] || dM.window == cfg["maxlag"]
            block_warning = cfg["blocksize"] < 10 * max(dE.tau, dM.tau)
            row = Any[runid, c.experiment, c.L, c.T, c.algorithm, c.chain, string(initial), seed,
                n, cfg["burn"], cfg["blocksize"], result.elapsed, mean(work), dE.tau, dM.tau,
                dE.window, dM.window, essE, essM, essE/result.elapsed, essM/result.elapsed,
                dE.tau*scale, dM.tau*scale, truncated, block_warning]
            for f in FIELDS
                s = getproperty(result.obs, f)
                s4 = getproperty(result.obs4, f)
                append!(row, [s.raw, s.stderr, s4.stderr])
            end
            println(summary, join(row, ','))
            flush(summary)
            open(joinpath(out, "traces", runid*".csv"), "w") do trace
                println(trace, "measurement,total_energy,total_magnetization,flipped_sites")
                for j in eachindex(E)
                    println(trace, j, ',', E[j], ',', M[j], ',', work[j])
                end
            end
            println("[$i/$(length(cases))] $runid $(c.experiment) L=$(c.L) T=$(round(c.T; digits=4)) $(c.algorithm): tauE=$(round(dE.tau; digits=2)), tau|M|=$(round(dM.tau; digits=2))")
        end
    end
    metadata["finished_utc"] = string(now(UTC))
    metadata["runs"] = length(cases)
    open(io -> TOML.print(io, metadata), joinpath(out, "metadata.toml"), "w")
    println("Saved to $out")
end

main()
