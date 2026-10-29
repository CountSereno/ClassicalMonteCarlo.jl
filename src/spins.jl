export Ising, Potts, energy, magnetization, metropolis!, wolff!, swendsen_wang!,
       WolffWorkspace, SwendsenWangWorkspace

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
    j = Float64(J)
    isfinite(j) || throw(ArgumentError("J must have a finite Float64 representation"))
    lx, ly, j
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
    b = Float64(beta)
    isfinite(b) || throw(ArgumentError("beta must have a finite Float64 representation"))
    b
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
        if metropolis_accept(rng, metropolis_logratio(b, dE))
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
    WolffWorkspace(model)

Reusable storage for one Wolff chain: a visited mask and a bounded FIFO queue.
The storage holds no physical state and is reset before cluster construction. It may be reused
for any Ising/Potts lattice with the same number of sites, but must not be shared
between concurrent updates. Passing it explicitly avoids per-cluster arrays.
"""
struct WolffWorkspace
    visited::BitVector
    queue::Vector{Int}
end
function WolffWorkspace(m::Union{Ising,Potts})
    n = length(m.spins)
    WolffWorkspace(falses(n), Vector{Int}(undef, n))
end

function _check_workspace(workspace::WolffWorkspace, m::Union{Ising,Potts})
    n = length(m.spins)
    length(workspace.visited) == n && length(workspace.queue) == n ||
        throw(ArgumentError("Wolff workspace must have one visited flag and queue slot per site"))
    nothing
end

"""
    wolff!(model, beta; rng=Random.default_rng(), workspace=WolffWorkspace(model))

Flip one Wolff cluster and return its size. Cluster updates require J >= 0.
One call is one cluster, not one lattice sweep. At zero coupling for Ising
and q=2 Potts, a probability-1/2 no-op is used to avoid a periodic chain;
the function returns 0 on that no-op. Reuse a workspace across calls to avoid
allocating the visited mask and cluster queue for every update.
"""
function wolff!(m::Union{Ising,Potts}, beta::Real;
                rng::AbstractRNG=Random.default_rng(),
                workspace::WolffWorkspace=WolffWorkspace(m))
    b = _cluster_beta(m, beta)
    _check_workspace(workspace, m)
    p = _cluster_probability(m, b)
    if b * m.J == 0 && (m isa Ising || m.q == 2) && rand(rng, Bool)
        return 0
    end
    seed = rand(rng, eachindex(m.spins))
    old = m.spins[seed]
    new = _other_colour(m, old, rng)
    visited, queue = workspace.visited, workspace.queue
    fill!(visited, false)
    visited[seed] = true
    queue[1] = seed
    cluster_size = 1
    head = 1
    while head <= cluster_size
        i = queue[head]
        for j in _neighbours(m, i)
            if !visited[j] && m.spins[j] == old && rand(rng) < p
                visited[j] = true
                cluster_size += 1
                queue[cluster_size] = j
            end
        end
        head += 1
    end
    for k in 1:cluster_size
        m.spins[queue[k]] = new
    end
    cluster_size
end

"""
    SwendsenWangWorkspace(model)

Reusable union-find and recolouring arrays for one SW chain. Each update resets
the disjoint sets and assigns new root colours. The number of sites and spin
element type must match the model. Do not share mutable workspace arrays across
concurrent chains.
"""
struct SwendsenWangWorkspace{T<:Integer}
    parent::Vector{Int}
    sizes::Vector{Int}
    colours::Vector{T}
end
function SwendsenWangWorkspace(m::Union{Ising,Potts})
    n = length(m.spins)
    SwendsenWangWorkspace(collect(1:n), ones(Int, n), similar(m.spins))
end

function _check_workspace(workspace::SwendsenWangWorkspace, m::Union{Ising,Potts})
    n = length(m.spins)
    length(workspace.parent) == n && length(workspace.sizes) == n &&
        length(workspace.colours) == n || throw(ArgumentError("SW workspace arrays must match the number of sites"))
    eltype(workspace.colours) === eltype(m.spins) ||
        throw(ArgumentError("SW workspace colour type must match the model spin type"))
    nothing
end

"""
    swendsen_wang!(model, beta; rng=Random.default_rng(), workspace=SwendsenWangWorkspace(model))

One Swendsen–Wang bond-and-recolour sweep; returns the number of clusters.
Cluster updates require J >= 0. Pass a persistent workspace to reuse all
union-find and colour arrays without changing the sampling or RNG order.
"""
function swendsen_wang!(m::Union{Ising,Potts}, beta::Real;
                       rng::AbstractRNG=Random.default_rng(),
                       workspace::SwendsenWangWorkspace=SwendsenWangWorkspace(m))
    b = _cluster_beta(m, beta)
    _check_workspace(workspace, m)
    p = _cluster_probability(m, b)
    n = length(m.spins)
    parent, sizes, colours = workspace.parent, workspace.sizes, workspace.colours
    _reset_sets!(parent, sizes)
    for i in 1:n
        _, right, _, up = _neighbours(m, i)
        for j in (right, up)
            if m.spins[i] == m.spins[j] && rand(rng) < p
                _union!(parent, sizes, i, j)
            end
        end
    end
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
