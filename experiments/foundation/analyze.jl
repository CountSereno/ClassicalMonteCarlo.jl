# Offline analysis: never starts a sampler or overwrites raw data.
using ClassicalMonteCarlo, Statistics, TOML, SHA, Printf
include("exact.jl")
include("tables.jl")
using .ExactIsing, .FoundationTables

const OBSERVABLES = (:energy_density,:abs_magnetization_density,:heat_capacity,
                    :susceptibility,:susceptibility_abs,:binder)

function analysis_source_hashes()
    repository=normpath(joinpath(@__DIR__,"../.."))
    paths=vcat([joinpath("experiments","foundation",f) for f in ("analyze.jl","exact.jl","tables.jl")],
               ["src/analysis.jl","Project.toml"])
    Dict(path=>bytes2hex(sha256(read(joinpath(repository,path)))) for path in paths)
end

function diagnostic(x,maxlag)
    if var(x) == 0
        return (tau=NaN,window=0,stderr_mean=NaN)
    end
    integrated_autocorrelation_time(x;maxlag)
end

function summarize_group(rows,key)
    n = length(rows)
    values = [r[key] for r in rows]
    within = sqrt(sum(r[key*"_stderr_max"]^2 for r in rows))/n
    between = std(values)/sqrt(n)
    (value=mean(values),se=max(within,between),within_se=within,between_se=between)
end

function analyze(output)
    hashes=analysis_source_hashes()
    manifest = TOML.parsefile(joinpath(output,"manifest.toml"))
    cfg = TOML.parsefile(joinpath(output,"config.toml"))
    exact = enumerate_ising(cfg["exact_L"],1/cfg["exact_T"])
    analysis = joinpath(output,"analysis"); mkpath(analysis)
    open(io->TOML.print(io,exact),joinpath(analysis,"exact.toml"),"w")
    summaries = Dict{String,Any}[]
    blocks = Dict{String,Any}[]
    correlations = Dict{String,Any}[]
    for r in manifest["runs"]
        folder = joinpath(output,"runs",r["run_id"])
        metadata = TOML.parsefile(joinpath(folder,"metadata.toml"))
        tracepath = joinpath(folder,"trace.csv")
        bytes2hex(sha256(read(tracepath))) == metadata["trace_sha256"] || error("trace hash mismatch")
        trace = read_trace(tracepath)
        length(trace.E) == r["samples"] || error("sample count mismatch")
        trace.updates == collect((r["burn"]+r["stride"]):r["stride"]:(r["burn"]+r["samples"]*r["stride"])) ||
            error("measurement update clock mismatch")
        N = r["L"]^2
        all(abs.(trace.E) .<= 2N) && all(abs.(trace.M) .<= N) || error("unphysical spin observable")
        all(mod.(trace.M .- N,2).==0) || error("magnetization parity violation")
        row = copy(r)
        for field in OBSERVABLES
            row[string(field)*"_stderr_max"] = 0.0
        end
        for b in cfg["block_sizes"]
            obs = ising_observables(trace.E,trace.M;beta=1/r["T"],nsites=N,blocksize=b)
            for field in OBSERVABLES
                f = string(field); stat = getproperty(obs,field)
                row[f] = stat.raw
                row[f*"_stderr_max"] = max(row[f*"_stderr_max"],stat.stderr)
                push!(blocks,Dict("run_id"=>r["run_id"],"experiment"=>r["experiment"],"L"=>r["L"],
                    "algorithm"=>r["algorithm"],"chain"=>r["chain"],"observable"=>f,
                    "blocksize"=>b,"nblocks"=>length(trace.E)÷b,"raw"=>stat.raw,"stderr"=>stat.stderr))
            end
        end
        dE = diagnostic(trace.E,cfg["maxlag"])
        dM = diagnostic(abs.(trace.M),cfg["maxlag"])
        # tau in units of stored records. Convert to named algorithm updates,
        # then, for Wolff only, a rough changed-spin/N work convention.
        row["tau_energy_records"] = dE.tau
        row["tau_absM_records"] = dM.tau
        row["tau_energy_updates"] = dE.tau*r["stride"]
        row["tau_absM_updates"] = dM.tau*r["stride"]
        row["tau_energy_spin_sweeps"] = dE.tau * (r["algorithm"]=="wolff" ? mean(trace.work)/N : r["stride"])
        row["ESS_energy_per_second"] = length(trace.E)/(2 * dE.tau*r["production_seconds"])
        row["ESS_absM_per_second"] = length(trace.E)/(2 * dM.tau*r["production_seconds"])
        row["window_energy"] = dE.window
        row["window_absM"] = dM.window
        row["window_truncated"] = dE.window==cfg["maxlag"] || dM.window==cfg["maxlag"]
        row["short_block_warning"] = minimum(cfg["block_sizes"]) < 10 * max(dE.tau,dM.tau)
        row["large_block_warning"] = maximum(cfg["block_sizes"]) < 10 * max(dE.tau,dM.tau)
        row["signed_m_mean"] = mean(trace.M)/N
        row["positive_m_fraction"] = mean(trace.M .> 0)
        push!(summaries,row)
        if r["experiment"]=="critical" && r["L"]==maximum(cfg["sizes"])
            for (name,x) in (("energy",trace.E),("absM",abs.(trace.M)),("signedM",trace.M))
                rho=autocorrelation(x;maxlag=cfg["maxlag"])
                for lag in 0:cfg["maxlag"]
                    push!(correlations,Dict("run_id"=>r["run_id"],"algorithm"=>r["algorithm"],
                        "chain"=>r["chain"],"observable"=>name,"lag"=>lag,"acf"=>rho[lag+1]))
                end
            end
        end
    end
    summaryfields = ["run_id","experiment","L","T","algorithm","chain","initial","seed","samples","burn","stride",
        "production_seconds","update_unit","work_unit","signed_m_mean","positive_m_fraction",
        "tau_energy_records","tau_absM_records","tau_energy_updates","tau_absM_updates","tau_energy_spin_sweeps",
        "ESS_energy_per_second","ESS_absM_per_second","window_energy","window_absM","window_truncated",
        "short_block_warning","large_block_warning"]
    for field in OBSERVABLES
        append!(summaryfields,[string(field),string(field)*"_stderr_max"])
    end
    write_table(joinpath(analysis,"chains.csv"),summaries,summaryfields)
    write_table(joinpath(analysis,"block_errors.csv"),blocks,["run_id","experiment","L","algorithm","chain","observable","blocksize","nblocks","raw","stderr"])
    write_table(joinpath(analysis,"acf.csv"),correlations,["run_id","algorithm","chain","observable","lag","acf"])
    aggregate = Dict{String,Any}[]
    for ex in ("exact","critical"),L in sort(unique([r["L"] for r in summaries if r["experiment"]==ex])),algorithm in cfg["algorithms"]
        selected = filter(r->r["experiment"]==ex && r["L"]==L && r["algorithm"]==algorithm,summaries)
        for field in OBSERVABLES
            key=string(field); stat=summarize_group(selected,key)
            exact_value=ex=="exact" ? exact[key] : NaN
            push!(aggregate,Dict("experiment"=>ex,"L"=>L,"T"=>selected[1]["T"],"algorithm"=>algorithm,
                "observable"=>key,"mean"=>stat.value,"diagnostic_se"=>stat.se,"within_se"=>stat.within_se,
                "between_se"=>stat.between_se,"exact"=>exact_value,
                "difference_over_se"=>(stat.value-exact_value)/stat.se))
        end
    end
    write_table(joinpath(analysis,"aggregate.csv"),aggregate,["experiment","L","T","algorithm","observable","mean","diagnostic_se","within_se","between_se","exact","difference_over_se"])
    report=Dict("runs"=>length(summaries),"measurements"=>manifest["measurements"],
        "production_seconds"=>sum(r["production_seconds"] for r in summaries),
        "short_block_warnings"=>sum(r["short_block_warning"] for r in summaries),
        "large_block_warnings"=>sum(r["large_block_warning"] for r in summaries),
        "acf_truncations"=>sum(r["window_truncated"] for r in summaries),
        "max_exact_discrepancy_over_diagnostic_se"=>maximum(abs(r["difference_over_se"]) for r in aggregate if r["experiment"]=="exact"),
        "error_rule"=>"max(between-chain SEM, quadrature of largest within-chain block-jackknife SE divided by nchains); diagnostic, not calibrated CI",
        "tau_rule"=>"first nonpositive ACF truncation; positive tail only; not equilibration proof",
        "julia_version"=>string(VERSION),
        "source_sha256"=>hashes)
    hashes==analysis_source_hashes() || error("analysis source changed while computing; rerun analysis")
    open(io->TOML.print(io,report),joinpath(analysis,"report.toml"),"w")
    open(joinpath(output,"RESULTS.md"),"w") do io
        println(io,"# Executed foundation study\n")
        println(io,"$(report["runs"]) independent-seed chains, $(report["measurements"]) measurements; Julia $(VERSION).\n")
        println(io,"Each chain uses $(cfg["burn"]) warmup updates, $(cfg["samples"]) retained measurements and stride $(cfg["stride"]). Block sizes: $(join(cfg["block_sizes"],", ")); maximum ACF lag: $(cfg["maxlag"]). Full traces, configuration, source hashes and exact-reference values are saved alongside this report.\n")
        println(io,"## Independent finite-system reference\n")
        println(io,"Exact enumeration covers all $(exact["states"]) states of the periodic L=$(cfg["exact_L"]) Ising model at T=$(cfg["exact_T"]).\n")
        println(io,"| Algorithm | Energy / N | Absolute magnetization / N | Binder U4 |\n|---|---:|---:|---:|")
        for algorithm in cfg["algorithms"]
            vals=String[]
            for observable in ("energy_density","abs_magnetization_density","binder")
                row=only(filter(r->r["experiment"]=="exact"&&r["algorithm"]==algorithm&&r["observable"]==observable,aggregate))
                push!(vals,@sprintf("%.6f ± %.6f",row["mean"],row["diagnostic_se"]))
            end
            println(io,"| $algorithm | ",join(vals," | ")," |")
        end
        println(io,@sprintf("| Exact | %.6f | %.6f | %.6f |",exact["energy_density"],exact["abs_magnetization_density"],exact["binder"]))
        println(io,"\nAcross the six observables and three algorithms, the largest discrepancy is ",@sprintf("%.2f",report["max_exact_discrepancy_over_diagnostic_se"])," diagnostic standard errors. This is a validation comparison, not a simultaneous coverage guarantee.\n")
        println(io,"## Critical algorithm comparison\n")
        println(io,"| L | Algorithm | Mean energy tau (native updates) | Mean energy ESS / second |\n|---:|---|---:|---:|")
        for L in cfg["sizes"],algorithm in cfg["algorithms"]
            selected=filter(r->r["experiment"]=="critical"&&r["L"]==L&&r["algorithm"]==algorithm,summaries)
            println(io,@sprintf("| %d | %s | %.3f | %.1f |",L,algorithm,mean(r["tau_energy_updates"] for r in selected),mean(r["ESS_energy_per_second"] for r in selected)))
        end
        println(io,"\nNative updates are Metropolis sweeps, Wolff cluster updates, or SW sweeps. They do not have equal cost. Production timing includes measurements, excludes compile/warmup/storage, and is a limited single-machine implementation benchmark.\n")
        println(io,"There are $(report["short_block_warnings"]) short-block warnings, $(report["large_block_warnings"]) large-block warnings and $(report["acf_truncations"]) positive-window truncation flags. No flag is an equilibration certificate. Only $(length(cfg["initial_states"])) chains per condition were used; error bars are diagnostic. No critical temperature or exponent was fitted.\n")
        println(io,"Block and window flags use energy and absolute M only. They do not certify the errors of signed-M susceptibility or nonlinear observables: slow magnetization-sector transitions can remain hidden. Signed-M means, sign fractions and its L=$(maximum(cfg["sizes"])) ACF are saved for further inspection.\n")
        println(io,"The saved partial run resumed to exactly the same measurements, counters and final spins as uninterrupted execution in this environment. See manifest.toml.\n")
        println(io,"![Exact benchmark](figures/exact_validation.png)\n\n![Algorithm efficiency](figures/algorithm_efficiency.png)\n\n![Correlation and block diagnostics](figures/correlation_diagnostics.png)")
    end
    println(report)
end

if abspath(PROGRAM_FILE)==@__FILE__
    analyze(length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(@__DIR__,"results","demo"))
end
