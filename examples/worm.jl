using ClassicalMonteCarlo

function main()
    rng = new_rng(13)
    beta = 0.3
    worm = IsingWorm(8, beta)
    for _ in 1:20_000
        worm_step!(rng, worm)
    end
    closed = 0
    nsteps = 500_000
    for _ in 1:nsteps
        worm_step!(rng, worm)
        closed += worm_closed(worm)
    end
    println("<M²>/N = ", nsteps / closed, " (dimensionless; physical susceptibility is beta times this)")
    println("Estimate uncertainty using independent runs or blocking of the closure indicator.")
end
main()
