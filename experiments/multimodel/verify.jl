# Independent checks of saved observations. Does not generate new samples.
using TOML, SHA, Test
include("tables.jl")
include("exact.jl")
using .MultiModelTables, .ExactModels

function verify(root; require_current_source=false)
    manifest = TOML.parsefile(joinpath(root, "manifest.toml"))
    config = TOML.parsefile(joinpath(root, "config.toml"))
    records = manifest["runs"]
    hashes = copy(manifest["source_sha256"])
    @test length(unique(row["seed"] for row in records)) == length(records)
    measurements = 0
    @testset "model-specific traces and provenance" begin
        for row in records
            folder = joinpath(root, "runs", row["run_id"])
            meta = TOML.parsefile(joinpath(folder, "metadata.toml"))
            trace_path = joinpath(folder, "trace.csv")
            @test bytes2hex(sha256(read(trace_path))) == meta["trace_sha256"]
            @test meta["schema_version"] == 2
            @test meta["complete"]
            @test meta["model"] == row["family"]
            @test meta["config"]["model"] == row["family"]
            @test meta["boundary"] == "periodic square"
            @test meta["trace_columns"] == String.(split(open(readline, trace_path), ','))
            @test meta["samples_retained"] == row["samples"]
            @test meta["updates_completed"] == row["burn"] + row["samples"] * row["stride"]
            @test meta["update_unit"] == row["update_unit"]
            @test meta["work_unit"] == row["work_unit"]
            @test meta["julia_version"] == manifest["julia_version"]
            for key in ("L", "q", "seed", "burn", "samples", "stride", "algorithm", "initial")
                @test meta["config"][key] == row[key]
            end
            if row["family"] == "Potts"
                @test meta["config"]["beta"] == row["beta"]
                @test meta["config"]["J"] == 1
            else
                @test meta["config"]["p"] == row["p"]
                @test meta["config"]["periodic"]
            end
            for (path, hash) in meta["source_sha256"]
                @test haskey(hashes, path) && hashes[path] == hash
            end
            trace = read_trace(trace_path)
            n = row["samples"]
            N, m = row["nvertices"], row["nedges"]
            @test all(values -> length(values) == n, values(trace))
            @test trace["update"] == collect(row["burn"] + row["stride"]:row["stride"]:row["burn"] + n * row["stride"])
            @test all(isinteger, trace["kernel_work"])
            if row["family"] == "Potts"
                @test Set(keys(trace)) == Set(["measurement", "update", "total_energy", "order_parameter", "simplex_squared", "kernel_work"])
                @test all(e -> -m <= e <= 0 && isinteger(e), trace["total_energy"])
                @test all(x -> 0 <= x <= 1, trace["order_parameter"])
                @test all(x -> -1e-12 <= x <= 1 + 1e-12, trace["simplex_squared"])
                max_counts = (trace["order_parameter"] .* (row["q"] - 1) .+ 1) .* N ./ row["q"]
                @test all(x -> isapprox(x, round(x); atol=1e-10), max_counts)
                @test all(x -> 0 <= x <= N * row["stride"], trace["kernel_work"])
            else
                @test Set(keys(trace)) == Set(["measurement", "update", "occupied_bonds", "cluster_count", "largest_cluster", "cluster_size_squared_sum", "kernel_work"])
                bonds = trace["occupied_bonds"]
                clusters = trace["cluster_count"]
                largest = trace["largest_cluster"]
                second = trace["cluster_size_squared_sum"]
                @test all(isinteger, bonds) && all(b -> 0 <= b <= m, bonds)
                @test all(isinteger, clusters) && all(k -> 1 <= k <= N, clusters)
                @test all(isinteger, largest) && all(c -> 1 <= c <= N, largest)
                @test all(isinteger, second) && all(s -> N <= s <= N^2, second)
                @test all(clusters .>= N .- bonds)
                @test all(N .<= clusters .* largest)
                @test all(second .>= largest.^2 .+ N .- largest)
                @test all(second .<= N .* largest)
                @test all(mod.(second .- N, 2) .== 0)
                @test all(x -> 0 <= x <= m * row["stride"], trace["kernel_work"])
                if row["algorithm"] == "independent"
                    @test all(==(m * row["stride"]), trace["kernel_work"])
                end
            end
            measurements += n
        end
        @test measurements == manifest["measurements"]
        @test all(audit["exactly_equal"] for audit in manifest["restart_audits"])
    end
    @testset "independent exact reference limits" begin
        L, q = config["exact_L"], config["potts_q"]
        free = enumerate_potts(L, q, 0.)
        @test free["energy_density"] ≈ -2 / q atol=1e-12
        @test free["simplex_second_moment"] ≈ 1 / L^2 atol=1e-12
        @test free["heat_capacity"] == 0
        frozen = enumerate_potts(L, q, 100.)
        @test frozen["energy_density"] ≈ -2 atol=1e-12
        @test frozen["order_parameter"] ≈ 1 atol=1e-12
        @test frozen["simplex_second_moment"] ≈ 1 atol=1e-12
        edges = square_edges(2; periodic=false)
        percolation = enumerate_random_cluster(4, edges, .37, 1.)
        @test percolation["bond_density"] ≈ .37 atol=1e-12
        @test percolation["bond_variance_per_edge"] ≈ .37 * .63 atol=1e-12
    end
    repository = normpath(joinpath(@__DIR__, "../.."))
    mismatches = [path for (path, hash) in hashes if !isfile(joinpath(repository, path)) || bytes2hex(sha256(read(joinpath(repository, path)))) != hash]
    if require_current_source
        isempty(mismatches) || error("current source differs: $(join(mismatches, ", "))")
    elseif !isempty(mismatches)
        @warn "Historical data passed; current source differs from the generation source" mismatches
    end
    println("Verified $(length(records)) unique-seed chains, $measurements records, exact reference limits and saved restart audits.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    paths = filter(!=("--require-current-source"), ARGS)
    root = isempty(paths) ? joinpath(@__DIR__, "results", "demo") : abspath(first(paths))
    verify(root; require_current_source="--require-current-source" in ARGS)
end
