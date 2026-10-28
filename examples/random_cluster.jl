using ClassicalMonteCarlo, Statistics
rng = new_rng(42)
state = BondConfiguration(square_bond_graph(6; periodic=true))
p, q = 0.5, 1.5
for _ in 1:300
    sweeny_sweep!(rng, state, p, q)
end
sizes = Float64[]
for _ in 1:1000
    sweeny_sweep!(rng, state, p, q)
    push!(sizes, first(cluster_sizes(state)))
end
println("Mean largest cluster size: ", mean(sizes))
println("Reference implementation: connectivity is rebuilt on each proposal.")
