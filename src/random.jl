export new_rng, independent_rngs

"""
    new_rng(seed::Integer)

Create an independently usable Julia `Xoshiro` pseudorandom generator. Give each
simulation its own generator and record the Julia version as well as the seed;
Julia does not promise identical streams across versions. Not cryptographic.
"""
function new_rng(seed::Integer)
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    return Random.Xoshiro(seed)
end

"""
    independent_rngs(seed, n)

Create `n` deterministically seeded generators. Seeds are drawn from a master
Xoshiro stream; this is convenient stream separation, not a proof of disjoint
substreams. Never share a mutable RNG across concurrently running simulations.
"""
function independent_rngs(seed::Integer, n::Integer)
    n >= 0 || throw(ArgumentError("n must be nonnegative"))
    master = new_rng(seed)
    return [Random.Xoshiro(rand(master, UInt64)) for _ in 1:n]
end
