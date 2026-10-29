using Test, Random
using ClassicalMonteCarlo

@testset "periodic spin models" begin
    rng = MersenneTwister(314159)
    @test_throws ArgumentError Ising(2)
    @test_throws ArgumentError Potts(3, 1)
    @test_throws ArgumentError Ising(3; J=Inf)
    @test_throws ArgumentError Ising(3; init=:bad)
    @test_throws ArgumentError metropolis!(Ising(3), -1)
    @test_throws ArgumentError wolff!(Ising(3; J=-1), 1)
    @test_throws ArgumentError swendsen_wang!(Potts(3, 3; J=-1), 1)

    i = Ising(3, 4; init=:up)
    p = Potts(3, 4, 3; init=:ordered)
    @test energy(i) == -24
    @test magnetization(i) == 12
    @test energy(p) == -24
    @test magnetization(p) == 1
    @test metropolis!(i, 0; rng, sweeps=0) == 0
    @test 0 <= metropolis!(i, 0; rng) <= 12
    @test wolff!(i, 0; rng) in (0, 1)
    @test swendsen_wang!(p, 0; rng) == 12
    @test all(s -> 1 <= s <= 3, p.spins)

    # The local energy increment must agree with recomputation for every site
    # and every possible replacement, including periodic boundary sites.
    for m in (Ising(3, 4; rng), Potts(3, 4, 3; rng))
        for idx in eachindex(m.spins)
            old = m.spins[idx]
            choices = m isa Ising ? (-old,) : (1, 2, 3)
            for new in choices
                before = energy(m)
                dE = ClassicalMonteCarlo._delta_energy(m, idx, new)
                m.spins[idx] = new
                @test energy(m) - before ≈ dE
                m.spins[idx] = old
            end
        end
    end
    for construct in (rng -> Ising(4; rng), rng -> Potts(4, 3; rng))
        a, b = construct(MersenneTwister(17)), construct(MersenneTwister(17))
        ra, rb = MersenneTwister(19), MersenneTwister(19)
        for update! in (metropolis!, wolff!, swendsen_wang!)
            @test update!(a, 0.4; rng=ra) == update!(b, 0.4; rng=rb)
            @test a.spins == b.spins
        end
    end
end

@testset "cluster workspaces preserve updates and reset scratch" begin
    for make in (rng -> Ising(3, 4; rng), rng -> Potts(3, 4, 3; rng))
        for (update!, make_workspace) in ((wolff!, WolffWorkspace),
                                          (swendsen_wang!, SwendsenWangWorkspace))
            cached = make(MersenneTwister(913))
            fresh = deepcopy(cached)
            rc, rf = MersenneTwister(27), MersenneTwister(27)
            workspace = make_workspace(cached)
            # Deliberately dirty scratch before reuse. It must not become part
            # of the physical Markov state or affect the random-number order.
            if workspace isa WolffWorkspace
                fill!(workspace.visited, true)
                fill!(workspace.queue, 0)
            else
                fill!(workspace.parent, -1)
                fill!(workspace.sizes, -1)
                fill!(workspace.colours, 0)
            end
            for beta in (0.0, 0.3, 0.8, 0.0, 0.5), _ in 1:20
                @test update!(cached, beta; rng=rc, workspace) == update!(fresh, beta; rng=rf)
                @test cached.spins == fresh.spins
            end
            @test rand(rc) == rand(rf)
        end
    end
    @test_throws ArgumentError wolff!(Ising(4), 0.4; workspace=WolffWorkspace(Ising(3)))
    @test_throws ArgumentError swendsen_wang!(Ising(4), 0.4; workspace=SwendsenWangWorkspace(Ising(3)))
    @test_throws ArgumentError swendsen_wang!(Potts(3, 3), 0.4; workspace=SwendsenWangWorkspace(Ising(3)))
    huge = big"1e10000"
    @test_throws ArgumentError Ising(3; J=huge)
    @test_throws ArgumentError metropolis!(Ising(3), huge)
    @test_throws ArgumentError wolff!(Ising(3), huge)
end

# Exhaustive 3×3 equilibrium means test the kernels against the distribution,
# not against another implementation of the updates. Fixed seeds and a generous
# tolerance reduce false failures while still catching wrong bond probabilities.
function exact_mean_energy(m, beta)
    n = length(m.spins)
    q = m isa Ising ? 2 : m.q
    z = 0.0
    ez = 0.0
    for code in 0:q^n-1
        x = code
        for i in 1:n
            digit = x % q
            m.spins[i] = m isa Ising ? 2digit - 1 : digit + 1
            x ÷= q
        end
        e = energy(m)
        w = exp(-beta * e)
        z += w
        ez += e * w
    end
    ez / z
end

@testset "spin algorithms sample exact small-system energy" begin
    beta = 0.3
    for make in (() -> Ising(3; init=:up), () -> Potts(3, 3; init=:ordered))
        expected = exact_mean_energy(make(), beta)
        for (seed, update!) in enumerate((metropolis!, wolff!, swendsen_wang!))
            m = make()
            rng = MersenneTwister(100 + seed)
            for _ in 1:2_000
                update!(m, beta; rng)
            end
            total = 0.0
            for _ in 1:30_000
                update!(m, beta; rng)
                total += energy(m)
            end
            @test abs(total / 30_000 - expected) < 0.22
        end
    end
end

@testset "zero-coupling spin chains are not parity restricted" begin
    # A non-lazy sweep on an even lattice flips exactly N=16 spins and
    # cannot change the parity of the number of up/colour-1 spins.
    for (make, beta) in ((() -> Ising(4; init=:up), 0.0),
                         (() -> Ising(4; J=-1, init=:up), 0.0),
                         (() -> Ising(4; J=0, init=:up), 0.4),
                         (() -> Potts(4, 2; init=:ordered), 0.0),
                         (() -> Potts(4, 2; J=0, init=:ordered), 0.4))
        for update! in (metropolis!, wolff!)
            update! === wolff! && make().J < 0 && continue
            a, b = make(), make()
            ra, rb = MersenneTwister(901), MersenneTwister(901)
            parities = Set{Bool}()
            outcomes = Set{Int}()
            for _ in 1:200
                outcome = update!(a, beta; rng=ra)
                @test outcome == update!(b, beta; rng=rb)
                @test a.spins == b.spins
                push!(outcomes, outcome)
                push!(parities, isodd(count(==(1), a.spins)))
            end
            @test length(parities) == 2
            update! === wolff! && @test outcomes == Set((0, 1))
        end
    end
end
