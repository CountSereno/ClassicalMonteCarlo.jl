using Test, Statistics

@testset "analysis" begin
    @test block_means(collect(1:6), 2) == [1.5, 3.5, 5.5]
    @test_throws ArgumentError block_means([1, 2, 3], 2)
    @test_throws ArgumentError block_means([1.0, NaN], 1)
    @test_throws ArgumentError block_means([1, 2], 0)

    v = [1.0, 2.0, 3.0, 4.0]
    jk = jackknife(mean, v)
    @test jk.estimate ≈ mean(v)
    @test jk.bias ≈ 0 atol=1e-12
    @test jk.stderr ≈ std(v)/sqrt(length(v))
    @test jackknife(mean, v; blocksize=2).stderr ≈ 1.0
    @test_throws ArgumentError jackknife(mean, v; blocksize=4)
    @test_throws ArgumentError jackknife(mean, v; blocksize=3)
    @test_throws ArgumentError jackknife(x -> Inf, v)

    @test autocorrelation([1.0, 2, 3]; maxlag=2) ≈ [1.0, 0.0, -0.5]
    tau = integrated_autocorrelation_time([1.0, 2, 3]; maxlag=2)
    @test tau.tau == 0.5
    @test tau.window == 0
    @test tau.stderr_mean ≈ sqrt(1/3)
    @test_throws ArgumentError autocorrelation(ones(5))
    @test_throws ArgumentError autocorrelation(v; maxlag=4)
    @test_throws ArgumentError integrated_autocorrelation_time(v; maxlag=0)

    L = [4.0, 8, 16, 32]
    y = 3 .+ 2 ./ L
    fit = finite_size_fit(L, y; correction_exponent=1)
    @test fit.intercept ≈ 3
    @test fit.slope ≈ 2
    @test fit.chi2 ≈ 0 atol=1e-28
    @test fit.dof == 2
    weighted = finite_size_fit(L, y; correction_exponent=1, sigma=fill(0.1, 4))
    @test weighted.intercept ≈ 3
    @test weighted.slope ≈ 2
    @test weighted.weighted && weighted.intercept_stderr > 0
    @test_throws ArgumentError finite_size_fit(L, y; correction_exponent=0)
    @test_throws ArgumentError finite_size_fit([1, 1, 1], [1, 2, 3]; correction_exponent=1)
    @test_throws ArgumentError finite_size_fit(L, y; correction_exponent=1, sigma=[0.1, 0.1, 0, 0.1])
end

@testset "paired Ising observables and block uncertainty" begin
    E = [-8.0, -4, -4, 0, -8, 0, -4, 0]
    M = [4.0, 2, -2, 0, -4, 0, 2, -2]
    beta = 0.4
    N = 4
    result = ising_observables(E, M; beta, nsites=N, blocksize=2)
    @test (result.n_samples, result.n_blocks, result.blocksize) == (8, 4, 2)
    @test result.energy_density.raw ≈ mean(E)/N
    @test result.abs_magnetization_density.raw ≈ mean(abs.(M))/N
    @test result.heat_capacity.raw ≈ beta^2*var(E; corrected=false)/N
    @test result.susceptibility.raw ≈ beta*var(M; corrected=false)/N
    @test result.susceptibility_abs.raw ≈ beta*var(abs.(M); corrected=false)/N
    @test result.binder.raw ≈ 1-mean(M.^4)/(3mean(M.^2)^2)
    @test result.energy_density.stderr ≈ std(block_means(E/N, 2))/2

    # A separate, direct delete-block calculation checks that correlated
    # moment estimates are resampled together, including signed/absolute M.
    definitions = (
        (e,m) -> mean(e)/N,
        (e,m) -> mean(abs.(m))/N,
        (e,m) -> beta^2*var(e; corrected=false)/N,
        (e,m) -> beta*var(m; corrected=false)/N,
        (e,m) -> beta*var(abs.(m); corrected=false)/N,
        (e,m) -> 1-mean(m.^4)/(3mean(m.^2)^2),
    )
    fields = (:energy_density, :abs_magnetization_density, :heat_capacity,
              :susceptibility, :susceptibility_abs, :binder)
    for (name, f) in zip(fields, definitions)
        replicas = [f(E[setdiff(1:8, (2k-1):2k)], M[setdiff(1:8, (2k-1):2k)]) for k in 1:4]
        value = getproperty(result, name)
        expected_bias = 3*(mean(replicas)-f(E,M))
        @test value.bias ≈ expected_bias atol=1e-13
        @test value.estimate ≈ f(E,M)-expected_bias atol=1e-13
        @test value.stderr ≈ sqrt(3/4*sum((replicas .- mean(replicas)).^2)) atol=1e-13
    end

    # Thermodynamics must be unchanged by an energy reference shift or a
    # global spin inversion. The large shift catches unstable raw E² methods.
    shifted = ising_observables(E .+ 1e12, -M; beta, nsites=N, blocksize=2)
    @test shifted.heat_capacity.raw ≈ result.heat_capacity.raw
    @test shifted.heat_capacity.stderr ≈ result.heat_capacity.stderr
    for name in (:abs_magnetization_density, :susceptibility, :susceptibility_abs, :binder)
        @test getproperty(shifted, name) == getproperty(result, name)
    end

    # Exhaust all 2^9 equally probable configurations at infinite temperature.
    # Independent spins give <M²>=N and U4=2/(3N), including finite-size terms.
    model = Ising(3)
    all_E = Float64[]
    all_M = Float64[]
    for mask in 0:511
        for i in 1:9
            model.spins[i] = isodd(mask >> (i-1)) ? 1 : -1
        end
        push!(all_E, energy(model))
        push!(all_M, magnetization(model))
    end
    exact = ising_observables(all_E, all_M; beta=0, nsites=9, blocksize=16)
    @test exact.energy_density.raw ≈ 0 atol=1e-14
    @test exact.binder.raw ≈ 2/27 atol=1e-13
    @test exact.heat_capacity.raw == exact.susceptibility.raw == 0

    @test_throws ArgumentError ising_observables(E, M[1:7]; beta, nsites=N)
    @test_throws ArgumentError ising_observables(E, M; beta, nsites=0)
    @test_throws ArgumentError ising_observables(E, M; beta=-1, nsites=N)
    @test_throws ArgumentError ising_observables(E, M; beta, nsites=N, blocksize=3)
    @test_throws ArgumentError ising_observables(E, M; beta, nsites=N, blocksize=8)
    @test_throws ArgumentError ising_observables(E, 2M; beta, nsites=N)
    @test_throws ArgumentError ising_observables(E, zeros(8); beta, nsites=N)
    @test_throws ArgumentError ising_observables([0., 1.], [0., 1.]; beta, nsites=N)
end
