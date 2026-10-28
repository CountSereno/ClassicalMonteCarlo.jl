export Ising, Potts, energy, magnetization, metropolis!, wolff!, swendsen_wang!

"""
    Ising(Lx, Ly=Lx; J=1, init=:random, rng=Random.default_rng())

A periodic square-lattice Ising configuration. The Hamiltonian is
`-J * sum(s[i]*s[j])` over bonds in the positive x and y directions.
Both dimensions must be at least 3. `init` is `:random`, `:up`, or `:down`.
"""
struct Ising
    Lx::Int
    Ly::Int
    J::Float64
    spins::Vector{Int8}
end

"""
    Potts(Lx, Ly, q; J=1, init=:random, rng=Random.default_rng())
    Potts(L, q; kwargs...)

A periodic square-lattice q-state Potts configuration with Hamiltonian
`-J * sum(spins[i] == spins[j])` over positive x and y bonds. Both
dimensions must be at least 3, and q must be at least 2. `init` is
`:random` or `:ordered`.
"""
struct Potts
    Lx::Int
    Ly::Int
    q::Int
    J::Float64
    spins::Vector{Int}
end

function _spin_parameters(Lx::Integer, Ly::Integer, J::Real)
    Lx >= 3 && Ly >= 3 || throw(ArgumentError("both lattice dimensions must be at least 3"))
    isfinite(J) || throw(ArgumentError("J must be finite"))
    lx, ly = Int(Lx), Int(Ly)
    Base.checked_mul(lx, ly)
    lx, ly, Float64(J)
end

function Ising(Lx::Integer, Ly::Integer=Lx; J::Real=1, init::Symbol=:random,
               rng::AbstractRNG=Random.default_rng())
    lx, ly, j = _spin_parameters(Lx, Ly, J)
    init in (:random, :up, :down) || throw(ArgumentError("init must be :random, :up, or :down"))
    spins = init === :random ? Int8.(2 .* rand(rng, Bool, lx * ly) .- 1) :
            fill(init === :up ? Int8(1) : Int8(-1), lx * ly)
    Ising(lx, ly, j, spins)
end

function Potts(Lx::Integer, Ly::Integer, q::Integer; J::Real=1, init::Symbol=:random,
               rng::AbstractRNG=Random.default_rng())
    lx, ly, j = _spin_parameters(Lx, Ly, J)
    q >= 2 || throw(ArgumentError("q must be at least 2"))
    qq = Int(q)
    init in (:random, :ordered) || throw(ArgumentError("init must be :random or :ordered"))
    spins = init === :random ? rand(rng, 1:qq, lx * ly) : ones(Int, lx * ly)
    Potts(lx, ly, qq, j, spins)
end
Potts(L::Integer, q::Integer; kwargs...) = Potts(L, L, q; kwargs...)

# Periodic neighbours, with x as the fast index. Dimensions >= 3 avoid self-loops
# and ambiguous double bonds in the single-site local energy.
@inline function _neighbours(m::Union{Ising,Potts}, i::Int)
    x = mod1(i, m.Lx)
    y = (i - 1) ÷ m.Lx + 1
    left = x == 1 ? i + m.Lx - 1 : i - 1
    right = x == m.Lx ? i - m.Lx + 1 : i + 1
    down = y == 1 ? i + (m.Ly - 1) * m.Lx : i - m.Lx
    up = y == m.Ly ? i - (m.Ly - 1) * m.Lx : i + m.Lx
    (left, right, down, up)
end

"""Total Hamiltonian, counting each positive-direction bond exactly once."""
function energy(m::Union{Ising,Potts})
    E = 0.0
    @inbounds for i in eachindex(m.spins)
        _, right, _, up = _neighbours(m, i)
        E += _bondenergy(m, m.spins[i], m.spins[right])
        E += _bondenergy(m, m.spins[i], m.spins[up])
    end
    E
end
@inline _bondenergy(m::Ising, a, b) = -m.J * a * b
@inline _bondenergy(m::Potts, a, b) = -m.J * (a == b)

"""Total signed Ising magnetization, or normalized Potts order parameter in [0,1]."""
magnetization(m::Ising) = sum(m.spins)
function magnetization(m::Potts)
    counts = zeros(Int, m.q)
    for s in m.spins
        counts[s] += 1
    end
    (m.q * maximum(counts) / length(m.spins) - 1) / (m.q - 1)
end

@inline function _delta_energy(m::Union{Ising,Potts}, i::Int, new)
    old = m.spins[i]
    delta = 0.0
    @inbounds for j in _neighbours(m, i)
        delta += _bondenergy(m, new, m.spins[j]) - _bondenergy(m, old, m.spins[j])
    end
    delta
end
@inline _proposal(m::Ising, i, rng) = -m.spins[i]
@inline function _proposal(m::Potts, i, rng)
    v = rand(rng, 1:(m.q - 1))
    v >= m.spins[i] ? v + 1 : v
end
function _check_beta(beta::Real)
    isfinite(beta) && beta >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    Float64(beta)
end

"""
    metropolis!(model, beta; rng=Random.default_rng(), sweeps=1)

Run `sweeps * N` proposals, selecting sites independently with replacement.
Returns the number of accepted proposals. Supports positive and negative J.
For the zero-coupling Ising and q=2 Potts cases, each proposal is skipped
with probability 1/2 to avoid parity-restricted or periodic sampled chains.
"""
function metropolis!(m::Union{Ising,Potts}, beta::Real;
                     rng::AbstractRNG=Random.default_rng(), sweeps::Integer=1)
    b = _check_beta(beta)
    sweeps >= 0 || throw(ArgumentError("sweeps must be nonnegative"))
    n = Base.checked_mul(Int(sweeps), length(m.spins))
    accepted = 0
    lazy = b * m.J == 0 && (m isa Ising || m.q == 2)
    for _ in 1:n
        lazy && rand(rng, Bool) && continue
        i = rand(rng, eachindex(m.spins))
        new = _proposal(m, i, rng)
        dE = _delta_energy(m, i, new)
        if dE <= 0 || rand(rng) < exp(-b * dE)
            m.spins[i] = new
            accepted += 1
        end
    end
    accepted
end

@inline _cluster_probability(m::Ising, b) = -expm1(-2 * b * m.J)
@inline _cluster_probability(m::Potts, b) = -expm1(-b * m.J)
@inline _other_colour(m::Ising, old, rng) = -old
@inline function _other_colour(m::Potts, old, rng)
    v = rand(rng, 1:(m.q - 1))
    v >= old ? v + 1 : v
end
@inline _random_colour(m::Ising, rng) = rand(rng, Bool) ? Int8(1) : Int8(-1)
@inline _random_colour(m::Potts, rng) = rand(rng, 1:m.q)
function _cluster_beta(m::Union{Ising,Potts}, beta::Real)
    b = _check_beta(beta)
    m.J >= 0 || throw(ArgumentError("cluster algorithms require J >= 0"))
    b
end

"""
    wolff!(model, beta; rng=Random.default_rng())

Flip one Wolff cluster and return its size. Cluster updates require J >= 0.
One call is one cluster, not one lattice sweep. At zero coupling for Ising
and q=2 Potts, a probability-1/2 no-op is used to avoid a periodic chain;
the function returns 0 on that no-op.
"""
function wolff!(m::Union{Ising,Potts}, beta::Real; rng::AbstractRNG=Random.default_rng())
    b = _cluster_beta(m, beta)
    p = _cluster_probability(m, b)
    if b * m.J == 0 && (m isa Ising || m.q == 2) && rand(rng, Bool)
        return 0
    end
    seed = rand(rng, eachindex(m.spins))
    old = m.spins[seed]
    new = _other_colour(m, old, rng)
    visited = falses(length(m.spins))
    visited[seed] = true
    cluster = Int[seed]
    head = 1
    while head <= length(cluster)
        i = cluster[head]
        for j in _neighbours(m, i)
            if !visited[j] && m.spins[j] == old && rand(rng) < p
                visited[j] = true
                push!(cluster, j)
            end
        end
        head += 1
    end
    for i in cluster
        m.spins[i] = new
    end
    length(cluster)
end

@inline function _findroot!(parent, i)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    i
end
@inline function _union!(parent, sizes, i, j)
    a, b = _findroot!(parent, i), _findroot!(parent, j)
    if a != b
        if sizes[a] < sizes[b]
            a, b = b, a
        end
        parent[b] = a
        sizes[a] += sizes[b]
    end
    nothing
end

"""
    swendsen_wang!(model, beta; rng=Random.default_rng())

One Swendsen–Wang bond-and-recolour sweep; returns the number of clusters.
Cluster updates require J >= 0.
"""
function swendsen_wang!(m::Union{Ising,Potts}, beta::Real; rng::AbstractRNG=Random.default_rng())
    b = _cluster_beta(m, beta)
    p = _cluster_probability(m, b)
    n = length(m.spins)
    parent = collect(1:n)
    sizes = ones(Int, n)
    for i in 1:n
        _, right, _, up = _neighbours(m, i)
        for j in (right, up)
            if m.spins[i] == m.spins[j] && rand(rng) < p
                _union!(parent, sizes, i, j)
            end
        end
    end
    colours = similar(m.spins)
    nclusters = 0
    for i in 1:n
        if _findroot!(parent, i) == i
            colours[i] = _random_colour(m, rng)
            nclusters += 1
        end
    end
    for i in 1:n
        m.spins[i] = colours[_findroot!(parent, i)]
    end
    nclusters
end
