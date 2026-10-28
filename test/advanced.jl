@testset "RNG and advanced sampling" begin
    @test rand(new_rng(42), UInt64, 8) == rand(new_rng(42), UInt64, 8)
    @test rand(independent_rngs(9, 2)[1], UInt64) != rand(independent_rngs(9, 2)[2], UInt64)
    @test_throws ArgumentError new_rng(-1)

    rng = new_rng(42)
    w = IsingWorm(3, 0.3)
    @test_throws ArgumentError IsingWorm(2, 0.3)
    @test_throws ArgumentError IsingWorm(3, 0.0)
    # Bond parity must equal the current defect set at every step.
    for _ in 1:3000
        worm_step!(rng, w)
        parity = falses(9)
        for site in 1:9
            row, col = divrem(site - 1, 3)
            for (edge, neighbor) in ((2site-1, row*3 + mod(col+1,3)+1),
                                     (2site, mod(row+1,3)*3 + col+1))
                if w.bonds[edge]
                    parity[site] = !parity[site]
                    parity[neighbor] = !parity[neighbor]
                end
            end
        end
        expected = falses(9)
        if w.head != w.tail
            expected[w.head] = expected[w.tail] = true
        end
        @test parity == expected
    end
    # Exactly enumerate 3×3 spins: dimensionless susceptibility <M²>/N at zero field (physical χ includes β).
    beta = 0.3
    Z = 0.0; M2 = 0.0
    for mask in 0:(2^9-1)
        spins = [isodd(mask >> (i-1)) ? 1 : -1 for i in 1:9]
        bonds = sum(spins[i] * (spins[div(i-1,3)*3+mod(i,3)+1] +
                    spins[mod(div(i-1,3)+1,3)*3+mod(i-1,3)+1]) for i in 1:9)
        weight = exp(beta*bonds)
        Z += weight
        M2 += weight * sum(spins)^2 / 9
    end
    w = IsingWorm(3, beta)
    rng = new_rng(739)
    closed = 0
    for step in 1:500_000
        worm_step!(rng, w)
        step > 10_000 && (closed += worm_closed(w))
    end
    @test 1 / (closed / 490_000) ≈ M2 / Z rtol=0.07

    # One state in bin 1, two in bin 2; uniform independent proposal is symmetric.
    result = wang_landau(new_rng(25), 1; nbins=2, bin=s -> s == 1 ? 1 : 2,
                         propose=(r,s) -> rand(r, 1:3), check_interval=500,
                         logf_final=1e-5, max_steps=1_000_000)
    @test result.converged
    @test result.log_density[2] - result.log_density[1] ≈ log(2) atol=0.15
    blocked = wang_landau(new_rng(1), 1; nbins=2, bin=identity,
                          propose=(r,s) -> s, check_interval=10, max_steps=100)
    @test !blocked.converged && blocked.steps == 100
    @test_throws ArgumentError wang_landau(new_rng(1), 1; nbins=0, bin=identity, propose=(r,s)->s)
end
