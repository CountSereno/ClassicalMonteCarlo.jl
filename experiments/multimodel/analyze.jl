# Offline inference from paired records; no sampler is invoked here.
using ClassicalMonteCarlo, Statistics, TOML, SHA, Printf
include("tables.jl")
include("exact.jl")
using .MultiModelTables, .ExactModels

const POTTS_FIELDS = (:energy_density, :order_parameter, :heat_capacity,
                      :order_parameter_fluctuation, :moment_ratio, :simplex_second_moment)
const RC_FIELDS = (:bond_density, :cluster_density, :largest_cluster_fraction,
                   :mean_cluster_size, :bond_variance_per_edge)

function analysis_hashes()
    repository = normpath(joinpath(@__DIR__, "../.."))
    paths = ["Project.toml", "src/analysis.jl", "src/model_observables.jl",
             "experiments/multimodel/analyze.jl", "experiments/multimodel/exact.jl",
             "experiments/multimodel/tables.jl"]
    Dict(path => bytes2hex(sha256(read(joinpath(repository, path)))) for path in paths)
end

function chain_statistics(row, trace, blocksize)
    if row["family"] == "Potts"
        potts_observables(trace["total_energy"], trace["order_parameter"];
            beta=row["beta"], nsites=row["nvertices"], blocksize,
            simplex_squared=trace["simplex_squared"])
    else
        random_cluster_observables(trace["occupied_bonds"], trace["cluster_count"],
            trace["largest_cluster"], trace["cluster_size_squared_sum"];
            nedges=row["nedges"], nvertices=row["nvertices"], blocksize)
    end
end

function correlation_diagnostic(values, maxlag)
    if var(values) == 0
        return (; tau=NaN, window=0, stderr_mean=NaN)
    end
    integrated_autocorrelation_time(values; maxlag)
end

function combine_chains(rows)
    n = length(rows)
    within = sqrt(sum(row["stderr_max"]^2 for row in rows)) / n
    between = std([row["raw"] for row in rows]) / sqrt(n)
    (; mean=mean(row["raw"] for row in rows), diagnostic_se=max(within, between),
       within_se=within, between_se=between, chains=n)
end

function analyze(root)
    hashes = analysis_hashes()
    config_path = joinpath(root, "config.toml")
    manifest_path = joinpath(root, "manifest.toml")
    config, manifest = TOML.parsefile(config_path), TOML.parsefile(manifest_path)
    input_hashes = Dict("config.toml" => bytes2hex(sha256(read(config_path))),
                        "manifest.toml" => bytes2hex(sha256(read(manifest_path))))
    output = joinpath(root, "analysis"); mkpath(output)
    L, q = config["exact_L"], config["potts_q"]
    beta = config["exact_beta"]
    p = -expm1(-beta)
    potts_exact = enumerate_potts(L, q, beta)
    edges = square_edges(L; periodic=true)
    rc_exact = Dict(string(Float64(rcq)) => enumerate_random_cluster(L^2, edges, p, rcq)
                    for rcq in config["rc_q_values"])
    integer_rc = rc_exact[string(Float64(q))]
    exact_energy_difference = potts_exact["energy_density"] + 2 * integer_rc["bond_density"] / p
    exact_simplex_difference = potts_exact["simplex_second_moment"] - integer_rc["mean_cluster_size"] / L^2
    abs(exact_energy_difference) < 1e-10 || error("independent exact FK energy relation failed")
    abs(exact_simplex_difference) < 1e-10 || error("independent exact FK simplex relation failed")
    open(io -> TOML.print(io, Dict("potts" => potts_exact, "random_cluster" => rc_exact)),
         joinpath(output, "exact.toml"), "w")

    estimates = Dict{String,Any}[]
    blocks = Dict{String,Any}[]
    diagnostics = Dict{String,Any}[]
    acfs = Dict{String,Any}[]
    for row in manifest["runs"]
        relative = joinpath("runs", row["run_id"], "trace.csv")
        path = joinpath(root, relative)
        metadata_path = joinpath(dirname(path), "metadata.toml")
        metadata = TOML.parsefile(metadata_path)
        digest = bytes2hex(sha256(read(path)))
        digest == metadata["trace_sha256"] || error("trace checksum mismatch")
        input_hashes[relative] = digest
        input_hashes[relpath(metadata_path, root)] = bytes2hex(sha256(read(metadata_path)))
        trace = read_trace(path)
        length(trace["measurement"]) == row["samples"] || error("sample count mismatch")
        fields = row["family"] == "Potts" ? POTTS_FIELDS : RC_FIELDS
        raw = Dict{String,Float64}()
        errors = Dict(string(field) => 0.0 for field in fields)
        identity = Dict(key => row[key] for key in ("run_id", "family", "experiment", "L", "q", "algorithm", "chain"))
        for blocksize in config["block_sizes"]
            result = chain_statistics(row, trace, blocksize)
            for field in fields
                name = string(field)
                stat = getproperty(result, field)
                raw[name] = stat.raw
                errors[name] = max(errors[name], stat.stderr)
                push!(blocks, merge(copy(identity), Dict("observable" => name,
                    "blocksize" => blocksize, "nblocks" => row["samples"] ÷ blocksize,
                    "raw" => stat.raw, "stderr" => stat.stderr)))
            end
        end
        for field in fields
            name = string(field)
            push!(estimates, merge(copy(identity), Dict("observable" => name,
                  "raw" => raw[name], "stderr_max" => errors[name])))
        end
        primary, secondary = row["family"] == "Potts" ?
            ("total_energy", "simplex_squared") : ("occupied_bonds", "cluster_size_squared_sum")
        first = correlation_diagnostic(trace[primary], config["maxlag"])
        second = correlation_diagnostic(trace[secondary], config["maxlag"])
        push!(diagnostics, merge(copy(row), Dict(
            "primary_observable" => primary, "secondary_observable" => secondary,
            "tau_primary_records" => first.tau, "tau_secondary_records" => second.tau,
            "tau_primary_updates" => first.tau * row["stride"],
            "tau_secondary_updates" => second.tau * row["stride"],
            "primary_window" => first.window, "secondary_window" => second.window,
            "primary_ESS_per_second" => row["samples"] / (2 * first.tau * row["production_seconds"]),
            "window_truncated" => first.window == config["maxlag"] || second.window == config["maxlag"],
            "small_block_warning" => minimum(config["block_sizes"]) < 10 * max(first.tau, second.tau),
            "large_block_warning" => maximum(config["block_sizes"]) < 10 * max(first.tau, second.tau),
            "mean_kernel_work" => mean(trace["kernel_work"]))))
        retain_acf = (row["experiment"] == "critical" && row["L"] == maximum(config["critical_sizes"])) ||
                     (row["family"] == "RandomCluster" && row["experiment"] == "exact" && row["q"] == 1)
        if retain_acf
            for observable in (primary, secondary)
                rho = autocorrelation(trace[observable]; maxlag=config["maxlag"])
                for lag in 0:config["maxlag"]
                    push!(acfs, merge(copy(identity), Dict("observable" => observable,
                          "lag" => lag, "acf" => rho[lag + 1])))
                end
            end
        end
    end
    identity_fields = ["run_id", "family", "experiment", "L", "q", "algorithm", "chain"]
    write_table(joinpath(output, "estimates.csv"), estimates,
                vcat(identity_fields, ["observable", "raw", "stderr_max"]))
    write_table(joinpath(output, "block_errors.csv"), blocks,
                vcat(identity_fields, ["observable", "blocksize", "nblocks", "raw", "stderr"]))
    write_table(joinpath(output, "acf.csv"), acfs,
                vcat(identity_fields, ["observable", "lag", "acf"]))
    write_table(joinpath(output, "chains.csv"), diagnostics,
        vcat(identity_fields, ["initial", "seed", "burn", "samples", "stride", "production_seconds",
            "update_unit", "work_unit", "mean_kernel_work", "primary_observable", "secondary_observable",
            "tau_primary_records", "tau_secondary_records", "tau_primary_updates", "tau_secondary_updates",
            "primary_window", "secondary_window", "primary_ESS_per_second", "window_truncated",
            "small_block_warning", "large_block_warning"]))

    group_fields = ["family", "experiment", "L", "q", "algorithm", "observable"]
    keys = unique([Tuple(row[field] for field in group_fields) for row in estimates])
    aggregate = Dict{String,Any}[]
    for key in keys
        rows = filter(row -> Tuple(row[field] for field in group_fields) == key, estimates)
        identity = Dict(field => value for (field, value) in zip(group_fields, key))
        reference = identity["family"] == "Potts" ? potts_exact : rc_exact[string(identity["q"])]
        exact_value = identity["experiment"] == "exact" ? reference[identity["observable"]] : NaN
        statistic = combine_chains(rows)
        push!(aggregate, merge(identity, Dict(string(k) => v for (k, v) in pairs(statistic)),
            Dict("exact" => exact_value, "difference_over_se" => (statistic.mean - exact_value) / statistic.diagnostic_se)))
    end
    write_table(joinpath(output, "aggregate.csv"), aggregate,
        vcat(group_fields, ["mean", "diagnostic_se", "within_se", "between_se", "chains", "exact", "difference_over_se"]))
    fk = Dict{String,Any}[]
    for size in config["critical_sizes"]
        select(family, observable) = only(filter(row -> row["family"] == family &&
            row["experiment"] == "critical" && row["L"] == size && row["observable"] == observable, aggregate))
        critical_p = sqrt(q) / (1 + sqrt(q))
        for (name, spin_name, bond_name, scale) in
            (("energy_density", "energy_density", "bond_density", -2 / critical_p),
             ("simplex_second_moment", "simplex_second_moment", "mean_cluster_size", 1 / size^2))
            spin, bond = select("Potts", spin_name), select("RandomCluster", bond_name)
            spin_value, bond_value = spin["mean"], scale * bond["mean"]
            spin_se, bond_se = spin["diagnostic_se"], abs(scale) * bond["diagnostic_se"]
            error = hypot(spin_se, bond_se)
            push!(fk, Dict("L" => size, "q" => q, "observable" => name,
                "spin_mean" => spin_value, "spin_se" => spin_se,
                "bond_mean" => bond_value, "bond_se" => bond_se,
                "difference" => spin_value - bond_value,
                "combined_se" => error, "difference_over_se" => (spin_value - bond_value) / error))
        end
    end
    write_table(joinpath(output, "fk_comparison.csv"), fk,
        ["L", "q", "observable", "spin_mean", "spin_se", "bond_mean", "bond_se", "difference", "combined_se", "difference_over_se"])
    report = Dict("runs" => length(manifest["runs"]), "measurements" => manifest["measurements"],
        "julia_version" => string(VERSION), "source_sha256" => hashes, "input_sha256" => input_hashes,
        "production_seconds" => sum(row["production_seconds"] for row in manifest["runs"]),
        "max_exact_discrepancy_over_se" => maximum(abs(row["difference_over_se"]) for row in aggregate if row["experiment"] == "exact"),
        "max_fk_discrepancy_over_se" => maximum(abs(row["difference_over_se"]) for row in fk),
        "exact_fk_energy_difference" => exact_energy_difference,
        "exact_fk_simplex_difference" => exact_simplex_difference,
        "small_block_warnings" => count(row["small_block_warning"] for row in diagnostics),
        "large_block_warnings" => count(row["large_block_warning"] for row in diagnostics),
        "window_truncations" => count(row["window_truncated"] for row in diagnostics),
        "error_rule" => "max(between-chain SEM, quadrature-combined largest paired-block SE/nchains); diagnostic, not calibrated coverage",
        "tau_rule" => "first-nonpositive ACF truncation; retained record units; no convergence guarantee")
    hashes == analysis_hashes() || error("analysis source changed while calculating")
    open(io -> TOML.print(io, report), joinpath(output, "report.toml"), "w")
    write_report(root, config, manifest, report, aggregate, fk, diagnostics)
    println("Analysed $(report["runs"]) chains; exact max discrepancy=$(report["max_exact_discrepancy_over_se"]), FK=$(report["max_fk_discrepancy_over_se"]) diagnostic SE.")
end

function write_report(root, config, manifest, report, aggregate, fk, diagnostics)
    open(joinpath(root, "RESULTS.md"), "w") do io
        println(io, "# Executed Potts and random-cluster study\n")
        println(io, "$(report["runs"]) independent-seed chains, $(report["measurements"]) recorded measurements; Julia $(VERSION). Each chain has $(config["burn"]) warmup updates and $(config["samples"]) production records, stride $(config["stride"]). Raw data, model-specific units, source hashes and restart evidence are saved.\n")
        println(io, "## Exact finite-system validation\n")
        println(io, "Periodic L=$(config["exact_L"]) Potts q=$(config["potts_q"]), beta=$(config["exact_beta"]), is enumerated independently over all q^(L²) spin assignments. The same graph's random-cluster measure is enumerated over all 2^(2L²) bond subsets at q=$(join(config["rc_q_values"], ", ")). The q=1 case also has independent Bernoulli redraws.\n")
        println(io, "| Family | q | Algorithm | Largest absolute discrepancy / diagnostic SE |\n|---|---:|---|---:|")
        conditions = unique((row["family"], row["q"], row["algorithm"]) for row in aggregate if row["experiment"] == "exact")
        for (family, q, algorithm) in conditions
            selected = filter(row -> row["family"] == family && row["q"] == q && row["algorithm"] == algorithm && row["experiment"] == "exact", aggregate)
            maximum_discrepancy = maximum(abs(row["difference_over_se"]) for row in selected)
            println(io, @sprintf("| %s | %.1f | %s | %.3f |", family, q, algorithm, maximum_discrepancy))
        end
        println(io, "\nThese many correlated comparisons are validation diagnostics, not a simultaneous confidence-coverage test. Full estimates and references are in analysis/aggregate.csv. Potts max-colour fluctuations are not relabelled as a signed magnetic susceptibility.\n")
        println(io, "## Cross-representation checks at the known q=$(config["potts_q"]) critical point\n")
        println(io, "Spin SW and bond Sweeny run independently on identical periodic graphs at beta=log(1+sqrt(q)) and p=sqrt(q)/(1+sqrt(q)). The checked identities are <E>/N=-<b>/(pN) and <simplex_squared>=<sum(cluster_size²)>/N².\n")
        println(io, "| L | Quantity | Spin estimate ± SE | Bond-derived estimate ± SE | Difference / combined SE |\n|---:|---|---:|---:|---:|")
        for row in fk
            println(io, @sprintf("| %d | %s | %.6f ± %.6f | %.6f ± %.6f | %.3f |", row["L"], row["observable"], row["spin_mean"], row["spin_se"], row["bond_mean"], row["bond_se"], row["difference_over_se"]))
        end
        println(io, "\nAt L=3 the two independent enumerators satisfy both identities to absolute error below 1e-10. The spin mapping requires integer q; q=1.5 is sampled and validated only in the random-cluster representation.\n")
        println(io, "## Correlation, cost and limitations\n")
        println(io, "| L | Family / algorithm | Mean primary tau (native updates) | Mean primary ESS / second |\n|---:|---|---:|---:|")
        for size in config["critical_sizes"], family in ("Potts", "RandomCluster")
            rows = filter(row -> row["experiment"] == "critical" && row["L"] == size && row["family"] == family, diagnostics)
            println(io, @sprintf("| %d | %s / %s | %.3f | %.1f |", size, family, first(rows)["algorithm"], mean(row["tau_primary_updates"] for row in rows), mean(row["primary_ESS_per_second"] for row in rows)))
        end
        println(io, "\nPrimary traces are Potts energy and RC occupied-bond count. They estimate related equilibrium means but have different transition processes and conditional noise; these ESS values are observable-specific. One spin sweep, one Wolff cluster, one Sweeny m-attempt sweep and one independent full redraw are different clocks. Production timings include measurement and are short, sequential, single-machine implementation timings. They exclude compilation, warmup, storage and analysis.\n")
        println(io, "The reference Sweeny kernel rebuilds connectivity for edge proposals; scratch reuse does not turn its O(m²)-per-sweep method into a dynamic-connectivity algorithm. No critical point, dynamic exponent, or phase-transition order is inferred from these sizes.\n")
        println(io, "Sweeny observations are one full sweep apart. This study cannot resolve sub-sweep relaxation or test published critical speeding-up; tau near the positive-window estimator's 0.5-record floor does not rule out faster modes.\n")
        println(io, "There are $(report["small_block_warnings"]) small-block warnings, $(report["large_block_warnings"]) largest-block warnings and $(report["window_truncations"]) ACF-window flags. These compare selected linear-observable correlations with block/window choices; they neither establish equilibrium nor certify errors of nonlinear quantities. Block sizes $(join(config["block_sizes"], ", ")) and all per-chain diagnostics remain available.\n")
        println(io, "$(length(manifest["restart_audits"])) model/algorithm restart audits reproduce paired observations, final physical configurations, update/work counters and future random draws exactly in this environment.\n")
        println(io, "![Exact benchmarks](figures/exact_validation.png)\n\n![FK identities](figures/fk_comparison.png)\n\n![Correlation diagnostics](figures/correlation_diagnostics.png)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze(isempty(ARGS) ? joinpath(@__DIR__, "results", "demo") : abspath(first(ARGS)))
end
