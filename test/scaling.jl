@testset "power laws, independent weighted reference and effective covariance" begin
    sizes=[8.,12.,16.,24.,32.,48.,64.]
    values=2.3 .* sizes.^1.75
    fit=power_law_fit(sizes,values;sigma=.01values)
    @test fit.exponent ≈ 1.75
    @test fit.amplitude ≈ 2.3
    @test fit.chi2 < 1e-20
    @test fit.dof == 5
    design=hcat(ones(7),log.(sizes))
    expected=inv(transpose(design)*design)*.01^2
    @test fit.covariance ≈ expected
    noisy=values .* exp.([.01,-.02,.015,-.013,.005,.01,-.008])
    @test power_law_fit(sizes,noisy;sigma=.01noisy).exponent ≈
        ((transpose(design)*design) \ (transpose(design)*log.(noisy)))[2]
    eff=effective_exponents(sizes,values;sigma=.01values)
    @test all(isapprox.(eff.exponents,1.75))
    @test eff.covariance[1,2] ≈ -.01^2/(log(12/8)*log(16/12))
    @test eff.covariance[1,3] == 0
    @test all(eff.stderr .> 0)
    @test power_law_fit(sizes,values).exponent ≈ 1.75
    for (L,y) in (([1,2],[2,3]),([1,1,3],[2,3,4]),([1,2,3],[2,0,4]),
                  ([1,2,Inf],[1,2,3]),([1,2,3],[1,2]),([3,2,1],[1,2,3]))
        @test_throws ArgumentError power_law_fit(L,y)
    end
    for sigma in ([0.,1.,2.],[-1.,1.,1.],[1.,1.], [1.,1.,Inf])
        @test_throws ArgumentError power_law_fit([1,2,3],[1,2,3];sigma)
    end
end
