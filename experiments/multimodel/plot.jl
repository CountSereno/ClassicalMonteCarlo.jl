# Pure Julia figures from saved analysis tables. No model updates occur here.
using CairoMakie, Statistics, TOML, SHA
include("tables.jl")
using .MultiModelTables

const COLOURS = Dict("metropolis" => "#2864a0", "wolff" => "#d07b35",
                     "swendsen_wang" => "#248b81", "sweeny" => "#7855a0",
                     "independent" => "#c75b42")
const LABELS = Dict("metropolis" => "Metropolis", "wolff" => "Wolff",
                    "swendsen_wang" => "Swendsen–Wang", "sweeny" => "Sweeny",
                    "independent" => "Independent redraw")
number(row, key) = parse(Float64, row[key])
select(rows; conditions...) = filter(row -> all(row[string(k)] == string(v) for (k, v) in conditions), rows)

function save_figure(directory, name, figure)
    save(joinpath(directory, name * ".png"), figure; px_per_unit=1.8)
    save(joinpath(directory, name * ".svg"), figure)
    save(joinpath(directory, name * ".pdf"), figure)
end

function points_with_errors!(axis, x, y, error; colour, label=nothing, marker=:circle)
    errorbars!(axis, x, y, error; color=colour, whiskerwidth=10, linewidth=1.5)
    scatter!(axis, x, y; color=colour, marker, markersize=10, label)
end

function plot_all(root)
    analysis = joinpath(root, "analysis")
    aggregate = read_table(joinpath(analysis, "aggregate.csv"))
    fk = read_table(joinpath(analysis, "fk_comparison.csv"))
    chains = read_table(joinpath(analysis, "chains.csv"))
    blocks = read_table(joinpath(analysis, "block_errors.csv"))
    acf = read_table(joinpath(analysis, "acf.csv"))
    config = TOML.parsefile(joinpath(root, "config.toml"))
    output = joinpath(root, "figures"); mkpath(output)
    set_theme!(Theme(fontsize=14, Axis=(xgridvisible=false,
        ygridcolor=(:black, .07), topspinevisible=false, rightspinevisible=false),
        Legend=(framevisible=false,)))

    figure = Figure(size=(1240, 760))
    Label(figure[0, 1:3], "Potts and random-cluster models: independent exact benchmarks", fontsize=23, font=:bold)
    algorithms = config["potts_algorithms"]
    for (column, (observable, title, ylabel)) in enumerate([
        ("energy_density", "Potts q=$(config["potts_q"]): energy", "⟨E⟩ / N"),
        ("order_parameter", "Potts: max-colour order", "⟨m_max⟩"),
        ("simplex_second_moment", "Potts: simplex second moment", "⟨m_simplex²⟩")])
        axis = Axis(figure[1, column]; title, ylabel,
            xticks=(1:length(algorithms), ["Metro.", "Wolff", "SW"]))
        rows = select(aggregate; family="Potts", experiment="exact", observable)
        hlines!(axis, [number(first(rows), "exact")]; color=:gray30, linestyle=:dash)
        for (x, algorithm) in enumerate(algorithms)
            row = only(select(rows; algorithm))
            points_with_errors!(axis, [x], [number(row, "mean")], [number(row, "diagnostic_se")]; colour=COLOURS[algorithm])
        end
        xlims!(axis, .5, length(algorithms) + .5)
    end
    axes = Axis[]
    for (column, (observable, title, ylabel)) in enumerate([
        ("bond_density", "Random cluster: occupied bonds", "⟨b⟩ / number of edges"),
        ("cluster_density", "Random cluster: component density", "⟨k⟩ / N"),
        ("mean_cluster_size", "Random cluster: size-biased mean", "⟨Σ cluster_size²⟩ / N")])
        axis = Axis(figure[2, column]; title, ylabel, xlabel="Cluster weight q", xticks=config["rc_q_values"])
        push!(axes, axis)
        rows = select(aggregate; family="RandomCluster", experiment="exact", observable, algorithm="sweeny")
        sort!(rows; by=row -> number(row, "q"))
        x = [number(row, "q") for row in rows]
        lines!(axis, x, [number(row, "exact") for row in rows]; color=:gray30,
               linestyle=:dash, label="Exact enumeration")
        points_with_errors!(axis, x, [number(row, "mean") for row in rows],
            [number(row, "diagnostic_se") for row in rows]; colour=COLOURS["sweeny"], label="Sweeny")
        redraw = only(select(aggregate; family="RandomCluster", experiment="exact", observable, algorithm="independent"))
        points_with_errors!(axis, [1.03], [number(redraw, "mean")], [number(redraw, "diagnostic_se")];
            colour=COLOURS["independent"], label="q=1 independent redraw", marker=:diamond)
    end
    Legend(figure[3, 1:3], axes[1]; orientation=:horizontal)
    Label(figure[4, 1:3], "Periodic L=$(config["exact_L"]); β=$(config["exact_beta"]), p=1−exp(−β) · four chains per condition · bars: diagnostic ±1 SE · redraw marker offset for visibility", fontsize=12, color=:gray35)
    save_figure(output, "exact_validation", figure)

    figure = Figure(size=(1280, 470))
    Label(figure[0, 1:3], "q=$(config["potts_q"]) at the known critical point: two independent representations", fontsize=23, font=:bold)
    quantity_labels = Dict("energy_density" => ("Energy identity", "Energy per site"),
                           "simplex_second_moment" => ("Simplex / FK-cluster identity", "Second moment"))
    comparisons = Axis[]
    for (column, observable) in enumerate(("energy_density", "simplex_second_moment"))
        title, ylabel = quantity_labels[observable]
        axis = Axis(figure[1, column]; title, ylabel, xlabel="Linear size L", xticks=config["critical_sizes"])
        push!(comparisons, axis)
        rows = select(fk; observable)
        sort!(rows; by=row -> number(row, "L"))
        x = [number(row, "L") for row in rows]
        for (prefix, offset, colour, label, marker) in
            (("spin", -.08, COLOURS["swendsen_wang"], "Potts / SW", :circle),
             ("bond", .08, COLOURS["sweeny"], "FK bonds / Sweeny", :diamond))
            y = [number(row, prefix * "_mean") for row in rows]
            lines!(axis, x .+ offset, y; color=colour, linewidth=1.5)
            points_with_errors!(axis, x .+ offset, y, [number(row, prefix * "_se") for row in rows]; colour, label, marker)
        end
    end
    axis = Axis(figure[1, 3]; title="Difference relative to uncertainty", xlabel="Linear size L",
                ylabel="(spin − bond) / combined SE", xticks=config["critical_sizes"])
    hlines!(axis, [0.]; color=:gray35, linestyle=:dash)
    for (observable, colour, marker) in (("energy_density", "#2864a0", :circle),
                                        ("simplex_second_moment", "#d07b35", :diamond))
        rows = select(fk; observable)
        sort!(rows; by=row -> number(row, "L"))
        scatterlines!(axis, [number(row, "L") for row in rows], [number(row, "difference_over_se") for row in rows];
            color=colour, marker, markersize=9, label=observable == "energy_density" ? "Energy" : "Simplex²")
    end
    Legend(figure[2, 1:2], comparisons[1]; orientation=:horizontal)
    Legend(figure[2, 3], axis; orientation=:horizontal)
    Label(figure[3, 1:3], "⟨E⟩/N = −⟨b⟩/(pN); ⟨m_simplex²⟩ = ⟨Σ cluster_size²⟩/N² · bars: diagnostic ±1 SE · no fitted critical point or exponent", fontsize=12, color=:gray35)
    save_figure(output, "fk_comparison", figure)

    figure = Figure(size=(1190, 800))
    maximum_size = maximum(config["critical_sizes"])
    Label(figure[0, 1:2], "Correlation, block errors and measured implementation cost", fontsize=23, font=:bold)
    free_axis = Axis(figure[1, 1]; title="q=1: independent bonds and local updates", xlabel="Lag in native recorded updates", ylabel="Bond-count ACF")
    critical_axis = Axis(figure[1, 2]; title="q=$(config["potts_q"]) critical, L=$maximum_size", xlabel="Lag in native recorded updates", ylabel="Energy / bond-count ACF")
    error_axis = Axis(figure[2, 1]; title="L=$maximum_size: energy-density block errors", xlabel="Block size (records)",
                      ylabel="Mean within-chain jackknife SE", xscale=log2, yscale=log10, xticks=config["block_sizes"])
    cost_axis = Axis(figure[2, 2]; title="Production updates plus observations", xlabel="Linear size L",
                     ylabel="Microseconds / retained record", yscale=log10, xticks=config["critical_sizes"])
    function draw_acf!(axis, rows, algorithm, name)
        lags = sort(unique(number(row, "lag") for row in rows))
        for chain in 1:length(config["rc_initials"])
            values = select(rows; chain)
            lines!(axis, [number(row, "lag") for row in values], [number(row, "acf") for row in values];
                   color=(COLOURS[algorithm], .18), linewidth=.8)
        end
        means = [mean(number(row, "acf") for row in rows if number(row, "lag") == lag) for lag in lags]
        lines!(axis, lags, means; color=COLOURS[algorithm], linewidth=2, label=name)
    end
    for algorithm in ("sweeny", "independent")
        rows = select(acf; family="RandomCluster", experiment="exact", q=1.0,
                      algorithm, observable="occupied_bonds")
        draw_acf!(free_axis, rows, algorithm, LABELS[algorithm])
    end
    for (family, algorithm, observable) in (("Potts", "swendsen_wang", "total_energy"),
                                           ("RandomCluster", "sweeny", "occupied_bonds"))
        rows = select(acf; family, experiment="critical", L=maximum_size, algorithm, observable)
        draw_acf!(critical_axis, rows, algorithm, family == "Potts" ? "Potts energy / SW" : "Bond count / Sweeny")
        field = family == "Potts" ? "energy_density" : "bond_density"
        scale = family == "Potts" ? 1. : 2 / (sqrt(config["potts_q"]) / (1 + sqrt(config["potts_q"])))
        rows = select(blocks; family, experiment="critical", L=maximum_size, algorithm, observable=field)
        for chain in 1:length(config["potts_initials"])
            values = select(rows; chain)
            errors = [number(only(select(values; blocksize)), "stderr") * scale for blocksize in config["block_sizes"]]
            lines!(error_axis, config["block_sizes"], errors; color=(COLOURS[algorithm], .2), linewidth=.8)
        end
        errors = [mean(number(row, "stderr") for row in select(rows; blocksize)) * scale for blocksize in config["block_sizes"]]
        scatterlines!(error_axis, config["block_sizes"], errors; color=COLOURS[algorithm], markersize=9)
        rows = select(chains; family, experiment="critical", algorithm)
        x = config["critical_sizes"]
        y = [mean(number(row, "production_seconds") / number(row, "samples") * 1e6 for row in select(rows; L)) for L in x]
        scatterlines!(cost_axis, x, y; color=COLOURS[algorithm], markersize=9, label=LABELS[algorithm])
    end
    for axis in (free_axis, critical_axis)
        hlines!(axis, [0.]; color=:gray50, linewidth=.8)
        ylims!(axis, -.1, 1.02)
    end
    xlims!(free_axis, 0, min(30, config["maxlag"]))
    xlims!(critical_axis, 0, min(120, config["maxlag"]))
    axislegend(free_axis; position=:rt)
    axislegend(critical_axis; position=:rt)
    Legend(figure[3, 1:2], cost_axis; orientation=:horizontal)
    Label(figure[4, 1:2], "Faint curves: individual chains · sweep / redraw clocks and per-record work differ · Sweeny uses reference connectivity rebuilding", fontsize=12, color=:gray35)
    Label(figure[5, 1:2], "Short single-machine timings; no dynamic-exponent fit · block/window diagnostics do not establish equilibrium", fontsize=12, color=:gray35)
    save_figure(output, "correlation_diagnostics", figure)

    plotting = normpath(joinpath(@__DIR__, "../plotting"))
    input_names = ["aggregate.csv", "fk_comparison.csv", "chains.csv", "block_errors.csv", "acf.csv"]
    metadata = Dict("julia_version" => string(VERSION), "cairomakie_version" => string(Base.pkgversion(CairoMakie)),
        "plot_source_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "table_reader_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "tables.jl")))),
        "config_sha256" => bytes2hex(sha256(read(joinpath(root, "config.toml")))),
        "input_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(analysis, name)))) for name in input_names),
        "plot_project_sha256" => bytes2hex(sha256(read(joinpath(plotting, "Project.toml")))),
        "plot_manifest_sha256" => bytes2hex(sha256(read(joinpath(plotting, "Manifest.toml")))))
    open(io -> TOML.print(io, metadata), joinpath(output, "metadata.toml"), "w")
    println("Saved three PNG/SVG/PDF figures to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_all(isempty(ARGS) ? joinpath(@__DIR__, "results", "demo") : abspath(first(ARGS)))
end
