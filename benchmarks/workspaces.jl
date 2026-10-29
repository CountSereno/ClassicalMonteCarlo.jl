# Reproduce temporary allocation counts without adding a benchmarking dependency.
# Usage: julia --project=. benchmarks/workspaces.jl > allocations.csv
using ClassicalMonteCarlo, Random

function repeat_updates!(f::F, n) where F
    for _ in 1:n
        f()
    end
    nothing
end

function compare_storage(model, algorithm, updates)
    make_workspace = algorithm === :wolff ? WolffWorkspace : SwendsenWangWorkspace
    update! = algorithm === :wolff ? wolff! : swendsen_wang!
    beta = log(1 + sqrt(2)) / 2
    a, b = deepcopy(model), deepcopy(model)
    ra, rb = new_rng(4161), new_rng(4161)
    workspace = make_workspace(b)
    fresh = () -> update!(a, beta; rng=ra)
    reused = () -> update!(b, beta; rng=rb, workspace=workspace)

    # Compile both paths before measuring. Identical calls preserve identical
    # random streams, which also checks that scratch reuse changes no physics.
    repeat_updates!(fresh, 100)
    repeat_updates!(reused, 100)
    a.spins == b.spins || error("workspace changed warm-up trajectory")
    GC.gc()
    fresh_bytes = @allocated repeat_updates!(fresh, updates)
    reused_bytes = @allocated repeat_updates!(reused, updates)
    a.spins == b.spins || error("workspace changed production trajectory")
    rand(ra, UInt64) == rand(rb, UInt64) || error("workspace changed RNG consumption")
    (fresh_bytes, reused_bytes)
end

function main()
    println("julia_version,model,L,algorithm,updates,fresh_bytes,reused_bytes")
    for (name, model) in (("Ising", Ising(16; rng=new_rng(17))),
                           ("Potts_q3", Potts(16,3; rng=new_rng(17)))),
        algorithm in (:wolff, :swendsen_wang)
        fresh, reused = compare_storage(model, algorithm, 1000)
        println(join((VERSION,name,16,algorithm,1000,fresh,reused),','))
    end
end

main()
