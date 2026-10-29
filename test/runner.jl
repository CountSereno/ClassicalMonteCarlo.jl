using Test, Random, TOML, SHA

function assert_same_run(a, b)
    @test a.updates == b.updates
    @test a.pending_work == b.pending_work
    @test a.model.spins == b.model.spins
    @test a.measurement_updates == b.measurement_updates
    @test a.energies == b.energies
    @test a.magnetizations == b.magnetizations
    @test a.kernel_work == b.kernel_work
    # Compare future draws without changing either running stream.
    @test rand(copy(a.rng), UInt64, 16) == rand(copy(b.rng), UInt64, 16)
end

@testset "validated experiment schedules" begin
    for overrides in ((L=2,), (beta=NaN,), (beta=-1,), (beta=Inf,),
                      (J=Inf,), (samples=0,), (burn=-1,), (stride=0,),
                      (seed=-1,), (algorithm=:unknown,), (initial=:unknown,),
                      (algorithm=:wolff, J=-1), (algorithm=:swendsen_wang, J=-1))
        options = merge((L=3, beta=.4), overrides)
        @test_throws ArgumentError RunSpec(; options...)
    end
    @test_throws OverflowError RunSpec(L=typemax(Int), beta=.4)
    @test_throws OverflowError RunSpec(L=3, beta=.4, samples=typemax(Int), stride=2)
    @test RunSpec(L=3, beta=0, J=-1).J == -1
    r = start_run(RunSpec(L=3, beta=.4, burn=2, samples=3, stride=4))
    @test advance!(r, 0) === r
    @test_throws ArgumentError advance!(r, -1)
    advance!(r, 5)
    @test isempty(r.energies)
    advance!(r, 1)
    @test r.measurement_updates == [6]
    run!(r)
    @test r.measurement_updates == [6, 10, 14]
    @test iscomplete(r)
    done = deepcopy(r)
    advance!(r, typemax(Int))
    run!(r)
    assert_same_run(r, done)
end

@testset "pause and checkpoint preserve the Markov clock" begin
    mktempdir() do dir
        for algorithm in (:metropolis, :wolff, :swendsen_wang)
            spec = RunSpec(L=4, beta=.44, algorithm=algorithm,
                           seed=737, burn=7, samples=18, stride=3)
            continuous = run!(start_run(spec))
            @test length(continuous.energies) == spec.samples
            @test continuous.measurement_updates == collect(10:3:61)
            # Pause inside warm-up, at its boundary, inside a measurement
            # interval, immediately after a measurement, and after completion.
            for pause_at in (0, 3, 7, 8, 10, 60, 61)
                partial = advance!(start_run(spec), pause_at)
                checkpoint = joinpath(dir, string(algorithm) * ".checkpoint")
                save_checkpoint(checkpoint, partial)
                restored = load_checkpoint(checkpoint; expected_spec=spec)
                assert_same_run(restored, partial)
                run!(restored)
                assert_same_run(restored, continuous)
            end
            # A different execution order cannot perturb another explicit RNG.
            a, b = start_run(spec), start_run(spec)
            other = start_run(RunSpec(L=3, beta=.1, seed=888, burn=2, samples=20))
            advance!(a, 20); run!(other); advance!(a, 17); run!(a)
            run!(b)
            assert_same_run(a, b)
        end
    end
end

# Rewrite only the safe TOML envelope to test rejection paths before native
# deserialization. Payload bytes remain untouched unless explicitly requested.
function rewrite_checkpoint_header(source, destination, change!)
    magic, header, payload = open(source, "r") do io
        magic = readline(io)
        n = parse(Int, readline(io))
        header = TOML.parse(String(read(io, n)))
        magic, header, read(io)
    end
    change!(header)
    text = sprint(io -> TOML.print(io, header))
    open(destination, "w") do io
        println(io, magic)
        println(io, ncodeunits(text))
        write(io, text)
        write(io, payload)
    end
end

@testset "checkpoint compatibility and corruption checks" begin
    mktempdir() do dir
        spec = RunSpec(L=3, beta=.5, burn=3, samples=5, stride=2, seed=101)
        state = advance!(start_run(spec), 8)
        original = joinpath(dir, "original.checkpoint")
        altered = joinpath(dir, "altered.checkpoint")
        save_checkpoint(original, state)
        @test_throws ArgumentError load_checkpoint(original;
            expected_spec=RunSpec(L=3, beta=.6, burn=3, samples=5, stride=2, seed=101))
        for change! in (h -> h["schema_version"] = 99,
                        h -> h["julia_version"] = "0.0.0",
                        h -> h["platform"] = "different-platform",
                        h -> h["rng_type"] = "UnknownRNG",
                        h -> h["source_sha256"]["src/runner.jl"] = "changed",
                        h -> h["config"]["beta"] = 9.0,
                        h -> h["payload_sha256"] = "corrupt",
                        h -> h["updates_completed"] = 100)
            rewrite_checkpoint_header(original, altered, change!)
            @test_throws ArgumentError load_checkpoint(altered)
        end
        bytes = read(original)
        bytes[end] ⊻= 0x01
        write(altered, bytes)
        @test_throws ArgumentError load_checkpoint(altered)
        write(altered, "not a checkpoint\n")
        @test_throws ArgumentError load_checkpoint(altered)
        write(altered, "ClassicalMonteCarlo checkpoint 1\n999999999\n")
        @test_throws ArgumentError load_checkpoint(altered)
        broken = deepcopy(state)
        push!(broken.energies, 0.0)
        @test_throws ArgumentError save_checkpoint(altered, broken)
        changed = deepcopy(state)
        changed.source_sha256["src/runner.jl"] = "not-the-loaded-source"
        @test_throws ArgumentError save_checkpoint(altered, changed)
        @test_throws ArgumentError export_run(joinpath(dir, "changed"), changed)
        @test !ispath(joinpath(dir, "changed"))
        @test_throws ArgumentError run_metadata(changed)
        for algorithm in (:wolff, :swendsen_wang)
            wrong_workspace = start_run(RunSpec(L=3, beta=.4, algorithm=algorithm))
            if algorithm == :wolff
                resize!(wrong_workspace.kernel.workspace.queue, 1)
            else
                resize!(wrong_workspace.kernel.workspace.parent, 1)
            end
            @test_throws ArgumentError save_checkpoint(altered, wrong_workspace)
        end
    end
end

@testset "portable data export and explicit work units" begin
    mktempdir() do dir
        for algorithm in (:metropolis, :wolff, :swendsen_wang)
            spec = RunSpec(L=3, beta=.3, algorithm=algorithm, burn=2,
                           samples=4, stride=2, seed=91)
            state = run!(start_run(spec))
            output = joinpath(dir, string(algorithm))
            export_run(output, state)
            @test_throws ArgumentError export_run(output, state)
            metadata = TOML.parsefile(joinpath(output, "metadata.toml"))
            @test metadata["complete"]
            @test metadata["updates_completed"] == 10
            @test metadata["samples_retained"] == 4
            @test metadata["update_unit"] == update_unit(state.kernel)
            @test metadata["work_unit"] == work_unit(state.kernel)
            @test metadata["trace_sha256"] == bytes2hex(sha256(read(joinpath(output, "trace.csv"))))
            lines = readlines(joinpath(output, "trace.csv"))
            @test lines[1] == "measurement,update,total_energy,total_magnetization,kernel_work"
            @test length(lines) == 5
            for (i, line) in enumerate(lines[2:end])
                fields = split(line, ',')
                @test parse(Int, fields[2]) == state.measurement_updates[i]
                @test parse(Float64, fields[3]) == state.energies[i]
                @test parse(Int, fields[4]) == state.magnetizations[i]
                @test parse(Int, fields[5]) == state.kernel_work[i]
            end
        end
        partial = advance!(start_run(RunSpec(L=3, beta=.2)), 1)
        export_run(joinpath(dir, "partial"), partial)
        @test !TOML.parsefile(joinpath(dir, "partial", "metadata.toml"))["complete"]
    end
end

@testset "source drift is rejected before labeling a run" begin
    # Exercise real file drift in an isolated package copy, never by editing the
    # user's source tree or files used by concurrently running experiments.
    mktempdir() do dir
        root = dirname(@__DIR__)
        cp(joinpath(root, "src"), joinpath(dir, "src"))
        cp(joinpath(root, "Project.toml"), joinpath(dir, "Project.toml"))
        sandbox = Module(:RunnerSourceDriftFixture)
        Base.include(sandbox, joinpath(dir, "src", "ClassicalMonteCarlo.jl"))
        # Enter the new definition world before looking up freshly included
        # bindings as well as calling their methods (Julia 1.12+).
        Base.invokelatest() do
            fixture = getfield(sandbox, :ClassicalMonteCarlo)
            spec = fixture.RunSpec(L=3, beta=.4, burn=1, samples=2)
            state = fixture.start_run(spec)
            open(joinpath(dir, "src", "spins.jl"), "a") do io
                println(io, "\n# Simulated edit after this package was loaded.")
            end
            @test_throws ArgumentError fixture.start_run(spec)
            @test_throws ArgumentError fixture.run_metadata(state)
            @test_throws ArgumentError fixture.export_run(joinpath(dir, "output"), state)
            @test_throws ArgumentError fixture.save_checkpoint(joinpath(dir, "state.cmc"), state)
            @test !ispath(joinpath(dir, "output"))
            @test !ispath(joinpath(dir, "state.cmc"))
        end
    end
end

@testset "source fingerprints invalidate compiled snapshots" begin
    # A fresh process must compile a new snapshot after Project-only edits or
    # edits to source files not included by the package. In-process edits still
    # fail closed, as tested above. Use an isolated writable compilation cache.
    mktempdir() do dir
        root = dirname(@__DIR__)
        cp(joinpath(root, "src"), joinpath(dir, "src"))
        cp(joinpath(root, "Project.toml"), joinpath(dir, "Project.toml"))
        script = "using ClassicalMonteCarlo; " *
                 "run_metadata(start_run(RunSpec(L=3, beta=0.4, burn=0, samples=1))); " *
                 "print(\"validated\")"
        command = `$(Base.julia_cmd()) --startup-file=no --project=$dir -e $script`
        separator = Sys.iswindows() ? ';' : ':'
        depot = join([joinpath(dir, "depot"); DEPOT_PATH], separator)
        isolated = addenv(command, "JULIA_LOAD_PATH" => "@$(separator)@stdlib",
                         "JULIA_DEPOT_PATH" => depot)
        @test read(isolated, String) == "validated"
        open(joinpath(dir, "Project.toml"), "a") do io
            println(io, "\n# Project-only provenance edit.")
        end
        @test read(isolated, String) == "validated"
        unused = joinpath(dir, "src", "audit_unused.jl")
        write(unused, "# Un-included provenance input.\n")
        @test read(isolated, String) == "validated"
        open(unused, "a") do io
            println(io, "# Edit existing un-included source.")
        end
        @test read(isolated, String) == "validated"
    end
end
