# Read-only integrity checks and independent reference limits. No resampling.
using TOML, SHA, Test
include("exact.jl")
include("tables.jl")
using .ExactIsing,.FoundationTables

function verify(root;require_current_source=false)
    manifest=TOML.parsefile(joinpath(root,"manifest.toml"))
    repository=normpath(joinpath(@__DIR__,"../.."))
    recorded_hashes=copy(manifest["source_sha256"])
    seeds=[r["seed"] for r in manifest["runs"]]
    @test length(seeds)==length(unique(seeds))
    count=0
    @testset "portable traces and physical limits" begin
        for row in manifest["runs"]
            folder=joinpath(root,"runs",row["run_id"])
            meta=TOML.parsefile(joinpath(folder,"metadata.toml"))
            @test meta["complete"]
            @test meta["samples_retained"]==row["samples"]
            @test meta["updates_completed"]==row["burn"]+row["samples"]*row["stride"]
            @test meta["update_unit"]==row["update_unit"]
            @test meta["work_unit"]==row["work_unit"]
            @test meta["julia_version"]==manifest["julia_version"]
            for key in ("L","algorithm","initial","seed","burn","samples","stride")
                @test meta["config"][key]==row[key]
            end
            @test meta["config"]["beta"]==1/row["T"]
            @test meta["config"]["J"]==1.0
            for (source,hash) in meta["source_sha256"]
                @test get(recorded_hashes,source,nothing)==hash
                recorded_hashes[source]=hash
            end
            path=joinpath(folder,"trace.csv")
            @test bytes2hex(sha256(read(path)))==meta["trace_sha256"]
            trace=read_trace(path);N=row["L"]^2
            @test length(trace.E)==row["samples"]
            @test all(isfinite,trace.E)
            @test all(abs.(trace.E).<=2N)
            @test all(mod.(trace.E.+2N,4).==0)
            @test all(abs.(trace.M).<=N)
            @test all(mod.(trace.M.-N,2).==0)
            @test all(x -> 0 <= x <= N*row["stride"], trace.work)
            @test trace.updates==collect((row["burn"]+row["stride"]):row["stride"]:(row["burn"]+row["samples"]*row["stride"]))
            count+=length(trace.E)
        end
        @test count==manifest["measurements"]
        @test manifest["restart_audit"]["exactly_equal"]
    end
    @testset "independent enumerator limits" begin
        free=enumerate_ising(3,0)
        @test free["energy_density"]≈0 atol=1e-14
        @test free["heat_capacity"]==0
        @test free["binder"]≈2/27
        @test sum(values(free["density_of_states"]))==512
        frozen=enumerate_ising(3,100)
        @test frozen["energy_density"]≈-2
        @test frozen["abs_magnetization_density"]≈1
        @test frozen["binder"]≈2/3
    end
    mismatches=[path for (path,hash) in recorded_hashes if !isfile(joinpath(repository,path)) || bytes2hex(sha256(read(joinpath(repository,path))))!=hash]
    if require_current_source
        isempty(mismatches) || error("generation source differs from recorded revision: $(join(mismatches,", "))")
    elseif !isempty(mismatches)
        @warn "Raw data checks passed; current generator source differs from recorded version" mismatches
    end
    println("Verified $(length(seeds)) unique-seed chains, $count records, checksums, clocks, physical limits, and exact enumeration limits.")
end

if abspath(PROGRAM_FILE)==@__FILE__
    paths=filter(!=("--require-current-source"),ARGS)
    verify(isempty(paths) ? joinpath(@__DIR__,"results","demo") : abspath(first(paths));require_current_source="--require-current-source" in ARGS)
end
