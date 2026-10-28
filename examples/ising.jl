using ClassicalMonteCarlo, Statistics

# A small demonstration, not a convergence or precision study.
rng = new_rng(20261028)
model = Ising(16; init=:random, rng=rng)
beta = log(1 + sqrt(2)) / 2
for _ in 1:1_000
    wolff!(model, beta; rng=rng)
end
energies = Float64[]
magnetizations = Float64[]
for _ in 1:10_000
    wolff!(model, beta; rng=rng)
    push!(energies, energy(model) / length(model.spins))
    push!(magnetizations, abs(magnetization(model)) / length(model.spins))
end
stats = integrated_autocorrelation_time(energies; maxlag=500)
println("Julia version: ", VERSION, "; seed: 20261028; L: 16; beta: ", beta)
println("<E>/N = ", mean(energies), " ± ", stats.stderr_mean)
println("<|M|>/N = ", mean(magnetizations))
println("tau_int(E) = ", stats.tau, " Wolff cluster steps; window = ", stats.window)
