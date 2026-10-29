# Specifications contain physical parameters and a shared update schedule.
# A model gets its own type: no irrelevant beta/field/q placeholders are stored.
abstract type AbstractRunSpec end
abstract type AbstractSpinRunSpec <: AbstractRunSpec end

function _run_schedule(L, seed, burn, samples, stride; minimum_L=3)
    L >= minimum_L || throw(ArgumentError("L must be at least $minimum_L"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    burn >= 0 || throw(ArgumentError("burn must be nonnegative"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    stride > 0 || throw(ArgumentError("stride must be positive"))
    l, s, w, n, d = Int(L), Int(seed), Int(burn), Int(samples), Int(stride)
    vertices = Base.checked_mul(l, l)
    Base.checked_mul(vertices, d)
    Base.checked_add(w, Base.checked_mul(n, d))
    l, s, w, n, d
end

function _spin_run_parameters(beta, J, algorithm)
    b, j = Float64(beta), Float64(J)
    isfinite(b) && b >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    isfinite(j) || throw(ArgumentError("J must be finite"))
    algorithm in (:metropolis, :wolff, :swendsen_wang) || throw(ArgumentError("unknown spin algorithm"))
    algorithm != :metropolis && j < 0 && throw(ArgumentError("cluster kernels require J >= 0"))
    b, j
end

"""
    RunSpec(; L, beta, J=1, algorithm=:metropolis, initial=:random,
              seed=0, burn=1000, samples=1000, stride=1)

Periodic square, zero-field Ising specification. This original constructor is
kept for compatibility. One update is a Metropolis sweep, a Wolff cluster, or
an SW sweep. Initial states are `:random`, `:up` or `:down`.

All specifications measure first at `burn + stride`, then every `stride`
updates, retaining exactly `samples` paired observations. Warm-up is not stored.
"""
struct RunSpec <: AbstractSpinRunSpec
    L::Int
    beta::Float64
    J::Float64
    algorithm::Symbol
    initial::Symbol
    seed::Int
    burn::Int
    samples::Int
    stride::Int

    function RunSpec(; L::Integer, beta::Real, J::Real=1,
                     algorithm::Symbol=:metropolis, initial::Symbol=:random,
                     seed::Integer=0, burn::Integer=1000,
                     samples::Integer=1000, stride::Integer=1)
        l, s, w, n, d = _run_schedule(L, seed, burn, samples, stride)
        b, j = _spin_run_parameters(beta, J, algorithm)
        initial in (:random, :up, :down) || throw(ArgumentError("unknown Ising initial state"))
        new(l, b, j, algorithm, initial, s, w, n, d)
    end
end

"""
    PottsRunSpec(; L, q, beta, J=1, algorithm=:metropolis, initial=:random,
                   seed=0, burn=1000, samples=1000, stride=1)

Periodic square, zero-field integer-q Potts specification, with `q >= 2`.
Initial states are `:random` or `:ordered` (all colour 1). The Hamiltonian is
`-J * sum(delta(s_i,s_j))`; q=2 does not use the Ising coupling convention.
The stored order parameter is `(q*max_a(N_a)/N - 1)/(q-1)`, not signed total M.
`simplex_squared = (q*sum_a((N_a/N)^2)-1)/(q-1)` records the squared simplex
magnetization for an independent Fortuin–Kasteleyn moment comparison.
Update and sampling units are the same named spin-kernel units as `RunSpec`.
"""
struct PottsRunSpec <: AbstractSpinRunSpec
    L::Int
    q::Int
    beta::Float64
    J::Float64
    algorithm::Symbol
    initial::Symbol
    seed::Int
    burn::Int
    samples::Int
    stride::Int

    function PottsRunSpec(; L::Integer, q::Integer, beta::Real, J::Real=1,
                          algorithm::Symbol=:metropolis, initial::Symbol=:random,
                          seed::Integer=0, burn::Integer=1000,
                          samples::Integer=1000, stride::Integer=1)
        l, s, w, n, d = _run_schedule(L, seed, burn, samples, stride)
        b, j = _spin_run_parameters(beta, J, algorithm)
        q >= 2 || throw(ArgumentError("Potts q must be at least 2"))
        initial in (:random, :ordered) || throw(ArgumentError("unknown Potts initial state"))
        new(l, Int(q), b, j, algorithm, initial, s, w, n, d)
    end
end

"""
    RandomClusterRunSpec(; L, p, q=1, periodic=true, algorithm=:sweeny,
                           initial=:empty, seed=0, burn=1000,
                           samples=1000, stride=1)

Random-cluster measure `p^b*(1-p)^(m-b)*q^k` on an L×L square graph. `k`
includes isolated vertices. Periodic boundaries require L≥3; open ones L≥1.
Sweeny uses `0 < p < 1`, `q > 0` (including nonintegers), and a graph with
at least one edge. One Sweeny sweep attempts m lazy edge updates with replacement.

`:independent` redraws all m bonds and is valid only for `q == 1`; p=0 and p=1
are allowed. Its update unit is a whole-graph Bernoulli redraw, not a Sweeny sweep.
Initial states are `:empty`, `:full`, or `:random` (independent Bernoulli(p)
bonds, which is only an equilibrium draw at q=1).
"""
struct RandomClusterRunSpec <: AbstractRunSpec
    L::Int
    p::Float64
    q::Float64
    periodic::Bool
    algorithm::Symbol
    initial::Symbol
    seed::Int
    burn::Int
    samples::Int
    stride::Int

    function RandomClusterRunSpec(; L::Integer, p::Real, q::Real=1,
                                  periodic::Bool=true, algorithm::Symbol=:sweeny,
                                  initial::Symbol=:empty, seed::Integer=0,
                                  burn::Integer=1000, samples::Integer=1000,
                                  stride::Integer=1)
        l, s, w, n, d = _run_schedule(L, seed, burn, samples, stride;
                                      minimum_L=periodic ? 3 : 1)
        pp, qq = Float64(p), Float64(q)
        isfinite(pp) && 0 <= pp <= 1 || throw(ArgumentError("p must be finite and in [0, 1]"))
        isfinite(qq) && qq > 0 || throw(ArgumentError("q must be finite and positive"))
        algorithm in (:sweeny, :independent) || throw(ArgumentError("unknown random-cluster algorithm"))
        if algorithm == :sweeny
            0 < pp < 1 || throw(ArgumentError("Sweeny requires 0 < p < 1"))
            l >= 2 || throw(ArgumentError("Sweeny requires at least one edge"))
        else
            qq == 1 || throw(ArgumentError("independent redraws require q == 1"))
        end
        initial in (:empty, :full, :random) || throw(ArgumentError("unknown random-cluster initial state"))
        # These integer observables and work counters must be representable.
        vertices = Base.checked_mul(l, l)
        Base.checked_mul(vertices, vertices)
        edges = periodic ? Base.checked_mul(2, vertices) : Base.checked_mul(2l, l-1)
        Base.checked_mul(edges, d)
        new(l, pp, qq, periodic, algorithm, initial, s, w, n, d)
    end
end

_model_name(::RunSpec) = "Ising"
_model_name(::PottsRunSpec) = "Potts"
_model_name(::RandomClusterRunSpec) = "RandomCluster"

function _spec_record(spec::AbstractRunSpec)
    record = Dict{String,Any}("model" => _model_name(spec))
    for name in fieldnames(typeof(spec))
        value = getfield(spec, name)
        record[string(name)] = value isa Symbol ? string(value) : value
    end
    record
end
