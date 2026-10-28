export BondGraph, square_bond_graph, BondConfiguration, sample_bond_percolation!,
       cluster_sizes, nclusters, sweeny_update!, sweeny_sweep!

"""
    BondGraph(nvertices, edges)

Undirected finite simple graph with vertices numbered `1:nvertices`. Every edge
is an unordered pair of distinct vertices; duplicate edges are rejected. This
representation intentionally excludes multigraphs and self-loops.
"""
struct BondGraph
    nvertices::Int
    edges::Vector{Tuple{Int,Int}}
    function BondGraph(nvertices::Integer, edges)
        nvertices > 0 || throw(ArgumentError("nvertices must be positive"))
        pairs = Tuple{Int,Int}[]
        seen = Set{Tuple{Int,Int}}()
        for e in edges
            length(e) == 2 || throw(ArgumentError("each edge must have two endpoints"))
            u, v = Int(e[1]), Int(e[2])
            1 <= u <= nvertices && 1 <= v <= nvertices || throw(ArgumentError("edge endpoint out of range"))
            u != v || throw(ArgumentError("self-loops are not supported"))
            edge = minmax(u, v)
            edge ∉ seen || throw(ArgumentError("duplicate edge"))
            push!(seen, edge)
            push!(pairs, edge)
        end
        new(Int(nvertices), pairs)
    end
end

"""
    square_bond_graph(L; periodic=false)

Two-dimensional `L × L` nearest-neighbor graph. Open boundaries allow `L ≥ 1`;
periodic boundaries require `L ≥ 3`, so distinct edges remain unambiguous.
"""
function square_bond_graph(L::Integer; periodic::Bool=false)
    L >= (periodic ? 3 : 1) || throw(ArgumentError("L is too small for the chosen boundary"))
    id(x, y) = x + (y - 1) * L
    edges = Tuple{Int,Int}[]
    for y in 1:L, x in 1:L
        if x < L
            push!(edges, (id(x, y), id(x + 1, y)))
        elseif periodic
            push!(edges, (id(x, y), id(1, y)))
        end
        if y < L
            push!(edges, (id(x, y), id(x, y + 1)))
        elseif periodic
            push!(edges, (id(x, y), id(x, 1)))
        end
    end
    BondGraph(L * L, edges)
end

"""An edge-occupation configuration on a `BondGraph`, initially empty."""
struct BondConfiguration
    graph::BondGraph
    occupied::BitVector
    function BondConfiguration(graph::BondGraph, occupied::AbstractVector{Bool}=falses(length(graph.edges)))
        length(occupied) == length(graph.edges) || throw(ArgumentError("one occupation value is required per edge"))
        new(graph, BitVector(occupied))
    end
end

"""
    sample_bond_percolation!(rng, state, p)

Replace all edge occupations with independent Bernoulli(`p`) draws. Site
percolation and spanning/wrapping observables are not implemented here.
"""
function sample_bond_percolation!(rng::AbstractRNG, state::BondConfiguration, p::Real)
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must be finite and in [0, 1]"))
    for i in eachindex(state.occupied)
        state.occupied[i] = rand(rng) < p
    end
    state
end
sample_bond_percolation!(state::BondConfiguration, p::Real) = sample_bond_percolation!(Random.default_rng(), state, p)

# Straightforward union-find rebuilt after each update: correct for deletions,
# though not intended for large-scale performance.
function _bond_components(state::BondConfiguration; skip::Int=0)
    parent = collect(1:state.graph.nvertices)
    sizes = ones(Int, state.graph.nvertices)
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        i
    end
    for (i, (u, v)) in enumerate(state.graph.edges)
        if i != skip && state.occupied[i]
            a, b = root(u), root(v)
            if a != b
                if sizes[a] < sizes[b]
                    a, b = b, a
                end
                parent[b] = a
                sizes[a] += sizes[b]
            end
        end
    end
    [root(i) for i in eachindex(parent)]
end

"""Sizes of all connected components, including isolated vertices, descending."""
function cluster_sizes(state::BondConfiguration)
    labels = _bond_components(state)
    counts = Dict{Int,Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    sort!(collect(values(counts)); rev=true)
end
"""Number of connected components, including isolated vertices."""
nclusters(state::BondConfiguration) = length(cluster_sizes(state))

# log of the proposed/current random-cluster weight. Skipping the chosen edge
# distinguishes bridges from edges whose endpoints remain connected.
function _sweeny_logratio(state::BondConfiguration, edge_index::Integer, p::Real, q::Real)
    isfinite(p) && 0 < p < 1 || throw(ArgumentError("Sweeny requires finite 0 < p < 1"))
    isfinite(q) && q > 0 || throw(ArgumentError("q must be finite and positive"))
    checkbounds(state.occupied, edge_index)
    labels = _bond_components(state; skip=Int(edge_index))
    u, v = state.graph.edges[edge_index]
    disconnected = labels[u] != labels[v]
    log_add_ratio = log(p) - log1p(-p) - (disconnected ? log(q) : 0)
    state.occupied[edge_index] ? -log_add_ratio : log_add_ratio
end

"""
    sweeny_update!(rng, state, p, q; edge=nothing) -> Bool

One lazy symmetric single-edge Metropolis proposal for the random-cluster measure
`p^b (1-p)^(m-b) q^k`, with `b` occupied bonds, `m` total bonds, and `k`
components including isolated vertices. `q > 0` may be noninteger. A uniformly
random edge is selected unless `edge` is supplied. The proposal is skipped
with probability 1/2, ensuring an aperiodic chain even at `p=1/2, q=1`.
Returns whether it flipped.
This reference implementation rebuilds connectivity for each proposal and is
therefore slow on large graphs. Requires `0 < p < 1`.
"""
function sweeny_update!(rng::AbstractRNG, state::BondConfiguration, p::Real, q::Real; edge::Union{Nothing,Integer}=nothing)
    isempty(state.occupied) && throw(ArgumentError("Sweeny update requires at least one edge"))
    i = isnothing(edge) ? rand(rng, eachindex(state.occupied)) : edge
    logratio = _sweeny_logratio(state, i, p, q)
    rand(rng) < 0.5 && return false
    if logratio >= 0 || log(rand(rng)) < logratio
        state.occupied[i] = !state.occupied[i]
        return true
    end
    false
end
sweeny_update!(state::BondConfiguration, p::Real, q::Real; kwargs...) = sweeny_update!(Random.default_rng(), state, p, q; kwargs...)

"""Attempt `m` uniformly random single-edge updates with replacement; return the number accepted."""
function sweeny_sweep!(rng::AbstractRNG, state::BondConfiguration, p::Real, q::Real)
    isempty(state.occupied) && throw(ArgumentError("Sweeny sweep requires at least one edge"))
    accepted = 0
    for _ in eachindex(state.occupied)
        accepted += sweeny_update!(rng, state, p, q)
    end
    accepted
end
sweeny_sweep!(state::BondConfiguration, p::Real, q::Real) = sweeny_sweep!(Random.default_rng(), state, p, q)
