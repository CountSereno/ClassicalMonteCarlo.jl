using ClassicalMonteCarlo, Random
# Toy density of states: one microscopic state in bin 1, two in bin 2.
# The exact log density difference is log(2). Proposals must be symmetric.
rng = new_rng(25)
result = wang_landau(rng, 1; nbins=2, bin=s -> s == 1 ? 1 : 2,
                    propose=(rng, state) -> rand(rng, 1:3),
                    check_interval=500, logf_final=1e-5, max_steps=1_000_000)
println("Completed schedule: ", result.converged, "; steps: ", result.steps)
println("log g(2) - log g(1): ", result.log_density[2] - result.log_density[1],
        "; exact: ", log(2))
