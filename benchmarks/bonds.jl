# Allocation probe, not a claim of wall-clock or asymptotic speedup.
using ClassicalMonteCarlo, Random

function repeat_bond_updates!(update!::F, repetitions) where F
    for _ in 1:repetitions
        update!()
    end
end

function main()
    L, repetitions = 8, 1000
    p, q = sqrt(3)/(1+sqrt(3)), 3.0
    graph = square_bond_graph(L;periodic=true)
    fresh, reused = BondConfiguration(graph), BondConfiguration(graph)
    rng_fresh, rng_reused = new_rng(813), new_rng(813)
    scratch = BondWorkspace(reused)
    update_fresh! = () -> sweeny_sweep!(rng_fresh,fresh,p,q;workspace=BondWorkspace(fresh))
    update_reused! = () -> sweeny_sweep!(rng_reused,reused,p,q;workspace=scratch)
    repeat_bond_updates!(update_fresh!,100)
    repeat_bond_updates!(update_reused!,100)
    @assert fresh.occupied == reused.occupied
    GC.gc()
    fresh_bytes = @allocated repeat_bond_updates!(update_fresh!,repetitions)
    reused_bytes = @allocated repeat_bond_updates!(update_reused!,repetitions)
    @assert fresh.occupied == reused.occupied
    @assert rand(rng_fresh,UInt64,16) == rand(rng_reused,UInt64,16)
    # Measurement scratch is reusable independently of the update kernel.
    measure! = () -> cluster_statistics(reused;workspace=scratch)
    repeat_bond_updates!(measure!,100)
    measured_bytes = @allocated repeat_bond_updates!(measure!,repetitions)
    println("julia_version,L,q,updates,fresh_sweep_bytes,reused_sweep_bytes,reused_measurement_bytes")
    println(join((VERSION,L,q,repetitions,fresh_bytes,reused_bytes,measured_bytes),','))
end

main()
