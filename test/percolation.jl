using Test, Random

@testset "bond percolation and random cluster" begin
    open_graph = square_bond_graph(3)
    torus = square_bond_graph(3; periodic=true)
    @test (open_graph.nvertices, length(open_graph.edges)) == (9, 12)
    @test (torus.nvertices, length(torus.edges)) == (9, 18)
    @test_throws ArgumentError square_bond_graph(2; periodic=true)
    @test_throws ArgumentError BondGraph(2, [(1, 2), (2, 1)])
    @test_throws ArgumentError BondGraph(2, [(1, 1)])
    @test_throws ArgumentError BondGraph(2, [(1, 3)])

    state = BondConfiguration(open_graph)
    @test cluster_sizes(state) == ones(Int, 9)
    sample_bond_percolation!(MersenneTwister(1), state, 1)
    @test all(state.occupied)
    @test cluster_sizes(state) == [9]
    sample_bond_percolation!(MersenneTwister(1), state, 0)
    @test !any(state.occupied)
    @test_throws ArgumentError sample_bond_percolation!(state, NaN)
    @test_throws ArgumentError sweeny_update!(state, 0, 2)
    @test_throws ArgumentError sweeny_update!(state, 0.5, 0)

    # Verify exact detailed balance for every configuration and every proposed
    # edge on a three-edge cycle, including occupied bridges and cycle edges.
    triangle = BondGraph(3, [(1, 2), (2, 3), (1, 3)])
    p, q = 0.37, 1.5
    for mask in 0:7, i in 1:3
        bits = BitVector([!iszero(mask & (1 << (j - 1))) for j in 1:3])
        s = BondConfiguration(triangle, bits)
        t = BondConfiguration(triangle, bits)
        t.occupied[i] = !t.occupied[i]
        weight(x) = p^count(x.occupied) * (1 - p)^(3 - count(x.occupied)) * q^nclusters(x)
        logr = ClassicalMonteCarlo._sweeny_logratio(s, i, p, q)
        @test exp(logr) ≈ weight(t) / weight(s)
        reverse = ClassicalMonteCarlo._sweeny_logratio(t, i, p, q)
        @test weight(s) * min(1, exp(logr)) ≈ weight(t) * min(1, exp(reverse))
    end

    # Independent q=1 samples have E[b]=m*p, separately from the Markov chain.
    rng = MersenneTwister(320)
    state = BondConfiguration(triangle)
    total = 0
    for _ in 1:20_000
        sample_bond_percolation!(rng, state, p)
        total += count(state.occupied)
    end
    @test abs(total / 20_000 - 3p) < 0.03

    # Compare a fixed-seed run with the exact ensemble probability on a small
    # graph. The deliberately wide tolerance allows autocorrelation.
    state = BondConfiguration(triangle)
    rng = MersenneTwister(107)
    for _ in 1:2_000
        sweeny_sweep!(rng, state, p, q)
    end
    occupied_sum = 0
    for _ in 1:30_000
        sweeny_sweep!(rng, state, p, q)
        occupied_sum += count(state.occupied)
    end
    weights = Float64[]
    ns = Int[]
    for mask in 0:7
        bits = BitVector([!iszero(mask & (1 << (j - 1))) for j in 1:3])
        s = BondConfiguration(triangle, bits)
        b = count(bits)
        push!(ns, b)
        push!(weights, p^b * (1-p)^(3-b) * q^nclusters(s))
    end
    @test abs(occupied_sum / 30_000 - sum(ns .* weights) / sum(weights)) < 0.05
end
