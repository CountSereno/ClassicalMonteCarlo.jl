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
