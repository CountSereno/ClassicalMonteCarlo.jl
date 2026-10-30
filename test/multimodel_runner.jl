using Test, Random, TOML, SHA

function assert_same_multimodel_run(a, b)
    @test typeof(a.spec) === typeof(b.spec)
    @test a.updates == b.updates
    @test a.pending_work == b.pending_work
    @test a.measurement_updates == b.measurement_updates
    @test a.kernel_work == b.kernel_work
    @test keys(observations(a)) == keys(observations(b))
    for name in keys(observations(a))
        @test getproperty(observations(a), name) == getproperty(observations(b), name)
    end
    if a.model isa BondConfiguration
        @test a.model.occupied == b.model.occupied
    else
        @test a.model.spins == b.model.spins
    end
    @test rand(copy(a.rng), UInt64, 16) == rand(copy(b.rng), UInt64, 16)
end

@testset "separate validated model specifications" begin
    for overrides in ((L=2,), (q=1,), (beta=-1,), (beta=Inf,), (J=Inf,),
                      (initial=:up,), (algorithm=:sweeny,), (seed=-1,),
                      (burn=-1,), (samples=0,), (stride=0,),
                      (algorithm=:wolff, J=-1))
        @test_throws ArgumentError PottsRunSpec(; merge((L=3, q=3, beta=.6), overrides)...)
    end
    @test_throws TypeError PottsRunSpec(L=3, q=2.5, beta=.6)
    @test_throws OverflowError PottsRunSpec(L=typemax(Int), q=3, beta=.6)
    @test PottsRunSpec(L=3, q=2, beta=0, J=-1).J == -1
    for overrides in ((L=2,), (p=-.1,), (p=NaN,), (p=0,), (p=1,),
                      (q=0,), (q=Inf,), (initial=:ordered,), (algorithm=:wolff,),
                      (algorithm=:independent, q=2), (seed=-1,),
                      (burn=-1,), (samples=0,), (stride=0,))
        @test_throws ArgumentError RandomClusterRunSpec(; merge((L=3, p=.4), overrides)...)
    end
    @test_throws ArgumentError RandomClusterRunSpec(L=1, p=.5, periodic=false)
    @test_throws OverflowError RandomClusterRunSpec(L=typemax(Int), p=.5)
    @test_throws OverflowError RunSpec(L=3, beta=.4, samples=1, stride=typemax(Int)÷2)
    @test RandomClusterRunSpec(L=3, p=.4, q=.7).q == .7
    for p in (0, 1)
        @test RandomClusterRunSpec(L=1, periodic=false, p=p, algorithm=:independent).p == p
    end
end

@testset "paired observations preserve physical meanings" begin
    old = run!(start_run(RunSpec(L=3, beta=.4, burn=0, samples=2)))
    @test old.energies === observations(old).total_energy
    @test old.magnetizations === observations(old).total_magnetization
    @test :energies in propertynames(old)
    potts = start_run(PottsRunSpec(L=3, q=3, beta=.5, burn=0, samples=1))
    @test keys(observations(potts)) == (:total_energy, :order_parameter, :simplex_squared)
    @test !hasproperty(potts, :magnetizations)
    # Ordered and balanced populations have exact known endpoint moments;
    # colour permutations cannot change either invariant order observable.
    for (spins, order, simplex) in ((ones(Int, 9), 1.0, 1.0),
                                   (repeat(1:3, 3), 0.0, 0.0),
                                   ([1,1,1,1,1,1,2,2,2], .5, 1/3))
        potts.model.spins .= spins
        before = copy(potts.rng)
        ClassicalMonteCarlo._record_observations!(observations(potts), potts.model, potts.spec, potts.kernel)
        @test last(observations(potts).order_parameter) ≈ order atol=1e-14
        @test last(observations(potts).simplex_squared) ≈ simplex atol=1e-14
        @test rand(copy(potts.rng), UInt64, 8) == rand(before, UInt64, 8)
    end
    binary = start_run(PottsRunSpec(L=3, q=2, beta=.5, burn=0, samples=1))
    binary.model.spins .= [ones(Int, 6); fill(2, 3)]
    ClassicalMonteCarlo._record_observations!(observations(binary), binary.model, binary.spec, binary.kernel)
    @test only(observations(binary).simplex_squared) ≈ (3/9)^2
    @test only(observations(binary).order_parameter) ≈ 3/9

    for p in (0, 1)
        state = run!(start_run(RandomClusterRunSpec(L=2, periodic=false, p=p,
            algorithm=:independent, initial=:random, burn=0, samples=3)))
        obs = observations(state)
        @test obs.occupied_bonds == fill(p == 0 ? 0 : 4, 3)
        @test obs.cluster_count == fill(p == 0 ? 4 : 1, 3)
        @test obs.largest_cluster == fill(p == 0 ? 1 : 4, 3)
        @test obs.cluster_size_squared_sum == fill(p == 0 ? 4 : 16, 3)
        @test state.kernel_work == fill(4, 3)
        @test update_unit(state.kernel) == "independent_bond_redraw"
        @test work_unit(state.kernel) == "resampled_bonds"
    end
    singleton = run!(start_run(RandomClusterRunSpec(L=1, periodic=false, p=.5,
        algorithm=:independent, burn=0, samples=2)))
    @test observations(singleton).cluster_count == [1, 1]
    @test singleton.kernel_work == [0, 0]
end

@testset "all model families share exact pause and restart" begin
    specs = AbstractRunSpec[
        (PottsRunSpec(L=3, q=3, beta=.8, algorithm=a, seed=719,
                      burn=4, samples=12, stride=3) for a in (:metropolis, :wolff, :swendsen_wang))...,
        RandomClusterRunSpec(L=2, periodic=false, p=.4, q=1.5, seed=719,
                             burn=4, samples=12, stride=3),
        RandomClusterRunSpec(L=3, p=.6, q=3, initial=:full, seed=719,
                             burn=4, samples=12, stride=3),
        RandomClusterRunSpec(L=3, p=.5, algorithm=:independent, seed=719,
                             burn=4, samples=12, stride=3),
    ]
    mktempdir() do dir
        for (index, spec) in enumerate(specs)
            continuous = run!(start_run(spec))
            @test continuous.measurement_updates == collect(7:3:40)
            for boundary in (0, 2, 4, 5, 7, 39, 40)
                partial = advance!(start_run(spec), boundary)
                path = joinpath(dir, "checkpoint-$index")
                save_checkpoint(path, partial)
                restored = load_checkpoint(path; expected_spec=spec)
                assert_same_multimodel_run(restored, partial)
                run!(restored)
                assert_same_multimodel_run(restored, continuous)
                @test run!(restored) === restored
            end
        end
    end
end

@testset "model-specific export and checkpoint guards" begin
    mktempdir() do dir
        specs = (PottsRunSpec(L=3, q=3, beta=.7, burn=1, samples=4, stride=2),
                 RandomClusterRunSpec(L=2, periodic=false, p=.4, q=1.5,
                                      burn=1, samples=4, stride=2))
        for (index, spec) in enumerate(specs)
            state = run!(start_run(spec))
            output = joinpath(dir, "run-$index")
            export_run(output, state)
            metadata = TOML.parsefile(joinpath(output, "metadata.toml"))
            @test metadata["schema_version"] == ClassicalMonteCarlo.RUN_SCHEMA_VERSION
            @test metadata["model"] == metadata["config"]["model"]
            @test metadata["samples_retained"] == 4
            @test metadata["complete"]
            @test metadata["trace_columns"][3:end-1] == string.(collect(keys(observations(state))))
            @test metadata["trace_sha256"] == bytes2hex(sha256(read(joinpath(output, "trace.csv"))))
            lines = readlines(joinpath(output, "trace.csv"))
            @test split(lines[1], ',') == metadata["trace_columns"]
            for (i, line) in enumerate(lines[2:end])
                fields = split(line, ',')
                @test parse(Int, fields[1]) == i
                @test parse(Int, fields[2]) == state.measurement_updates[i]
                for (j, column) in enumerate(values(observations(state)))
                    @test parse(eltype(column), fields[j+2]) == column[i]
                end
            end
            checkpoint = joinpath(dir, "state-$index")
            save_checkpoint(checkpoint, state)
            @test_throws ArgumentError load_checkpoint(checkpoint; expected_spec=RunSpec(L=3, beta=.7))
            corrupted = deepcopy(state)
            if spec isa PottsRunSpec
                observations(corrupted).order_parameter[1] = -1
            else
                observations(corrupted).cluster_count[1] = 0
            end
            @test_throws ArgumentError save_checkpoint(joinpath(dir, "bad"), corrupted)
            @test_throws ArgumentError export_run(joinpath(dir, "bad-export"), corrupted)
            @test !ispath(joinpath(dir, "bad-export"))
            corrupted = deepcopy(state)
            corrupted.pending_work = 1
            @test_throws ArgumentError save_checkpoint(joinpath(dir, "bad"), corrupted)
        end
        graphstate = start_run(RandomClusterRunSpec(L=2, periodic=false, p=.5))
        graphstate.model.graph.edges[1] = (1, 4)
        @test_throws ArgumentError save_checkpoint(joinpath(dir, "bad-graph"), graphstate)
        pottsstate = start_run(PottsRunSpec(L=3, q=3, beta=.7))
        pottsstate.model.spins[1] = 4
        @test_throws ArgumentError save_checkpoint(joinpath(dir, "bad-spin"), pottsstate)
        impossible = run!(start_run(RandomClusterRunSpec(L=3, p=.5, q=1.5,
            burn=0, samples=1)))
        observations(impossible).occupied_bonds[1] = 8
        observations(impossible).cluster_count[1] = 2
        observations(impossible).largest_cluster[1] = 1
        observations(impossible).cluster_size_squared_sum[1] = 9
        @test_throws ArgumentError save_checkpoint(joinpath(dir, "bad-partition"), impossible)
        redraw = run!(start_run(RandomClusterRunSpec(L=2, periodic=false, p=.5,
            algorithm=:independent, burn=0, samples=2, stride=3)))
        redraw.kernel_work[1] -= 1
        @test_throws ArgumentError export_run(joinpath(dir, "bad-redraw"), redraw)
        old = joinpath(dir, "schema-one")
        write(old, "ClassicalMonteCarlo checkpoint 1\n1\nx")
        @test_throws ArgumentError load_checkpoint(old)
    end
end
