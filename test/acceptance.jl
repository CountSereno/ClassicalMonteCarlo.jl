using Test, Random
using ClassicalMonteCarlo

@testset "portable Metropolis-Hastings filter" begin
    @test metropolis_accept(0.0, 0.0)
    @test metropolis_accept(0.0, prevfloat(1.0))
    @test metropolis_accept(Inf, 0.9)
    @test !metropolis_accept(-Inf, 0.0)
    @test metropolis_accept(-1000.0, 0.0)
    @test !metropolis_accept(-1000.0, nextfloat(0.0))
    @test metropolis_accept(log(0.25), 0.2)
    @test !metropolis_accept(log(0.25), 0.25)
    @test !metropolis_accept(log(0.25), 0.3)
    for u in (-0.1, 1.0, Inf, NaN)
        @test_throws ArgumentError metropolis_accept(Inf, u)
    end
    @test_throws ArgumentError metropolis_accept(NaN, 0.5)

    @test metropolis_logratio(2, 3) == -6
    @test metropolis_logratio(0, typemax(Int)) == 0
    @test metropolis_logratio(2, typemax(Int)) < 0 # no integer wraparound
    @test metropolis_logratio(2, -3; log_proposal_ratio=-1) == 5
    @test metropolis_logratio(1e300, 1e300) == -Inf
    @test metropolis_logratio(1e300, -1e300) == Inf
    @test metropolis_logratio(1, 1; log_proposal_ratio=-Inf) == -Inf
    for beta in (-1, Inf, NaN)
        @test_throws ArgumentError metropolis_logratio(beta, 1)
    end
    for de in (Inf, -Inf, NaN)
        @test_throws ArgumentError metropolis_logratio(1, de)
    end
    @test_throws ArgumentError metropolis_logratio(1, 1; log_proposal_ratio=NaN)
    @test_throws ArgumentError metropolis_logratio(1e300, 1e300; log_proposal_ratio=Inf)

    # The random adapter's draw count is explicit and independently testable.
    for logratio in (0.0, 1.0, Inf, -Inf)
        actual, reference = MersenneTwister(314), MersenneTwister(314)
        @test metropolis_accept(actual, logratio) == (logratio >= 0)
        @test rand(actual) == rand(reference)
    end
    actual, reference = MersenneTwister(12), MersenneTwister(12)
    @test metropolis_accept(actual, -0.7) == metropolis_accept(-0.7, rand(reference))
    @test rand(actual) == rand(reference)

    # Numerically integrate the pure decisions over a fixed uniform grid. The
    # resulting asymmetric-proposal transition matrix must preserve the chosen
    # distribution. A sign error in the Hastings term fails detailed balance.
    target = [0.2, 0.3, 0.5]
    proposal = [0.0 0.8 0.2; 0.1 0.0 0.9; 0.6 0.4 0.0]
    energies = -log.(target)
    transition = zeros(3, 3)
    grid_length = 10_000
    for i in 1:3, j in 1:3
        i == j && continue
        correction = log(proposal[j, i]) - log(proposal[i, j])
        ratio = metropolis_logratio(1, energies[j]-energies[i]; log_proposal_ratio=correction)
        accepted = count(k -> metropolis_accept(ratio, (k-0.5)/grid_length), 1:grid_length)
        transition[i, j] = proposal[i, j] * accepted/grid_length
    end
    for i in 1:3
        transition[i, i] = 1-sum(transition[i, :])
    end
    flux = target .* transition
    @test maximum(abs.(flux-transpose(flux))) < 5e-5
    @test maximum(abs.(transpose(transition)*target-target)) < 5e-5
end

@testset "acceptance source is standalone" begin
    # Include in an otherwise empty namespace: no model or package globals may
    # be required by this small, portable component.
    host = Module(:StandaloneAcceptanceTest)
    Base.include(host, joinpath(@__DIR__, "..", "src", "acceptance.jl"))
    filter = getproperty(getproperty(host, :Acceptance), :metropolis_accept)
    @test Base.invokelatest(filter, -1.0, 0.1)
    @test !Base.invokelatest(filter, -1.0, 0.9)
end
