# Follow up a finite-window warning without generating or changing any samples.
using ClassicalMonteCarlo, Statistics, TOML, SHA, Printf
include("tables.jl")
using .FoundationTables

function sensitivity_source_hashes()
    repository = normpath(joinpath(@__DIR__, "../.."))
    sources = ["Project.toml", "src/analysis.jl",
               "experiments/foundation/window_sensitivity.jl",
               "experiments/foundation/tables.jl"]
    Dict(path => bytes2hex(sha256(read(joinpath(repository, path)))) for path in sources)
end

function window_sensitivity(root)
    hashes = sensitivity_source_hashes()
    config = TOML.parsefile(joinpath(root, "config.toml"))
    manifest = TOML.parsefile(joinpath(root, "manifest.toml"))
    analysis = joinpath(root, "analysis")
    primary_path = joinpath(analysis, "chains.csv")
    primary = Dict(row["run_id"] => row for row in read_table(primary_path))
    L = maximum(config["sizes"])
    base_lag = config["maxlag"]
    limits = [base_lag, 2base_lag, 4base_lag]
    runs = filter(row -> row["experiment"] == "critical" && row["L"] == L,
                  manifest["runs"])
    rows = Dict{String,Any}[]
    input_hashes = Dict("analysis/chains.csv" => bytes2hex(sha256(read(primary_path))))
    input_hashes["config.toml"] = bytes2hex(sha256(read(joinpath(root, "config.toml"))))
    input_hashes["manifest.toml"] = bytes2hex(sha256(read(joinpath(root, "manifest.toml"))))

    for run in runs
        relative_path = joinpath("runs", run["run_id"], "trace.csv")
        trace_path = joinpath(root, relative_path)
        metadata = TOML.parsefile(joinpath(dirname(trace_path), "metadata.toml"))
        trace_hash = bytes2hex(sha256(read(trace_path)))
        trace_hash == metadata["trace_sha256"] || error("trace checksum mismatch: $relative_path")
        input_hashes[relative_path] = trace_hash
        trace = read_trace(trace_path)
        maximum(limits) < length(trace.E) || error("sensitivity lag must be shorter than trace")
        for (observable, values, baseline_key) in (("energy", trace.E, "tau_energy_records"),
                                                  ("absM", abs.(trace.M), "tau_absM_records"))
            baseline_tau = parse(Float64, primary[run["run_id"]][baseline_key])
            for limit in limits
                result = integrated_autocorrelation_time(values; maxlag=limit)
                # Re-analysis must reproduce the already published baseline exactly.
                limit == base_lag && result.tau != baseline_tau &&
                    error("baseline analysis disagrees for $(run["run_id"]): $observable")
                push!(rows, Dict("run_id" => run["run_id"], "L" => L,
                    "algorithm" => run["algorithm"], "chain" => run["chain"],
                    "observable" => observable, "maxlag" => limit,
                    "tau_records" => result.tau, "tau_updates" => result.tau * run["stride"],
                    "window" => result.window, "truncated" => result.window == limit,
                    "tau_over_baseline" => result.tau / baseline_tau,
                    "ESS_over_baseline" => baseline_tau / result.tau,
                    "ESS_per_second" => length(values) / (2 * result.tau * run["production_seconds"])))
            end
        end
    end
    hashes == sensitivity_source_hashes() || error("analysis source changed during computation")
    write_table(joinpath(analysis, "window_sensitivity.csv"), rows,
        ["run_id", "L", "algorithm", "chain", "observable", "maxlag", "tau_records",
         "tau_updates", "window", "truncated", "tau_over_baseline", "ESS_over_baseline", "ESS_per_second"])
    info = Dict("julia_version" => string(VERSION), "L" => L,
        "runs" => length(runs), "maxlags" => limits, "new_samples" => 0,
        "source_sha256" => hashes, "input_sha256" => input_hashes,
        "rule" => "same first-nonpositive ACF estimator as primary analysis; only maximum lag changes",
        "baseline_exactly_reproduced" => true,
        "truncations_by_limit" => Dict(string(limit) => count(row -> row["maxlag"] == limit && row["truncated"], rows)
                                       for limit in limits))
    open(io -> TOML.print(io, info), joinpath(analysis, "window_sensitivity.toml"), "w")
    open(joinpath(root, "WINDOW_SENSITIVITY.md"), "w") do io
        println(io, "# Offline ACF-window sensitivity at L=$L\n")
        println(io, "Reanalysed $(length(runs)) saved critical-temperature chains with maximum lags $(join(limits, ", ")). No samples were generated and the primary analysis remains at lag $base_lag. All baseline tau values were reproduced exactly; input and analysis-source hashes are recorded in analysis/window_sensitivity.toml.\n")
        println(io, "| Algorithm | Observable | Mean tau at $(limits[1]) | At $(limits[2]) | At $(limits[3]) | Truncated at $(limits[3]) |")
        println(io, "|---|---|---:|---:|---:|---:|")
        for algorithm in config["algorithms"], observable in ("energy", "absM")
            group = filter(row -> row["algorithm"] == algorithm && row["observable"] == observable, rows)
            means = [mean(row["tau_updates"] for row in group if row["maxlag"] == limit) for limit in limits]
            ntruncated = count(row -> row["maxlag"] == last(limits) && row["truncated"], group)
            println(io, @sprintf("| %s | %s | %.3f | %.3f | %.3f | %d |", algorithm, observable, means..., ntruncated))
        end
        println(io, "\nTau uses native algorithm updates, so different algorithms have different clocks. For a fixed chain the inferred ESS varies inversely with tau; timing data are reused without retiming.\n")
        changed = filter(row -> row["maxlag"] == last(limits) && row["tau_over_baseline"] > 1 + 1e-12, rows)
        if !isempty(changed)
            println(io, "| Changed chain | Observable | Tau ratio $(last(limits)) / $base_lag | ESS ratio | Last included lag |")
            println(io, "|---|---|---:|---:|---:|")
            for row in changed
                println(io, @sprintf("| %s (%s) | %s | %.4f | %.4f | %d |", row["run_id"], row["algorithm"], row["observable"], row["tau_over_baseline"], row["ESS_over_baseline"], row["window"]))
            end
        end
        println(io, "\nThis follow-up tests one window choice on the same finite traces. Reaching a first zero before the new limit does not establish stationarity, remove long-tail noise, or validate signed-M susceptibility errors. It does not replace longer independent chains or a separate equilibration study.")
    end
    println("Window-sensitivity flags by maximum lag: ", info["truncations_by_limit"])
    println("Saved sensitivity table, source/input hashes and WINDOW_SENSITIVITY.md; no new samples.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    window_sensitivity(isempty(ARGS) ? joinpath(@__DIR__, "results", "demo") : abspath(first(ARGS)))
end
