export BondGraph, square_bond_graph, BondConfiguration, sample_bond_percolation!,
       cluster_sizes, nclusters, sweeny_update!, sweeny_sweep!,
       BondWorkspace, cluster_statistics, sweeny_logratio

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

"""
    BondWorkspace(graph_or_state)

Reusable union-find scratch arrays for Sweeny updates and bond measurements.
They contain neither physical occupations nor an RNG. Connectivity is rebuilt
on each call, so deleting a bridge is handled correctly. This avoids temporary
array allocation but is not a dynamic-connectivity algorithm. Each concurrent
chain must own its workspace; only the number of vertices must match.
"""
struct BondWorkspace
    parent::Vector{Int}
    sizes::Vector{Int}
end
BondWorkspace(graph::BondGraph) = BondWorkspace(collect(1:graph.nvertices), ones(Int, graph.nvertices))
BondWorkspace(state::BondConfiguration) = BondWorkspace(state.graph)

function _check_workspace(workspace::BondWorkspace, state::BondConfiguration)
    n = state.graph.nvertices
    length(workspace.parent) == length(workspace.sizes) == n ||
        throw(ArgumentError("bond workspace arrays must match the number of vertices"))
    workspace.parent !== workspace.sizes || throw(ArgumentError("bond workspace arrays must not alias"))
    nothing
end

# Return root labels in workspace.parent. The returned vector is borrowed:
# another call using the same workspace replaces its contents.
function _bond_components(state::BondConfiguration; skip::Int=0,
                          workspace::BondWorkspace=BondWorkspace(state))
    _check_workspace(workspace, state)
    0 <= skip <= length(state.occupied) || throw(ArgumentError("invalid skipped edge"))
    parent, sizes = workspace.parent, workspace.sizes
    _reset_sets!(parent, sizes)
    for (i, (u, v)) in enumerate(state.graph.edges)
        if i != skip && state.occupied[i]
            _union!(parent, sizes, u, v)
        end
    end
    for i in eachindex(parent)
        parent[i] = _findroot!(parent, i)
    end
    parent
end

"""Sizes of all connected components, including isolated vertices, descending."""
function cluster_sizes(state::BondConfiguration; workspace::BondWorkspace=BondWorkspace(state))
    labels = _bond_components(state; workspace=workspace)
    sort!([workspace.sizes[i] for i in eachindex(labels) if labels[i] == i]; rev=true)
end
"""Number of connected components, including isolated vertices."""
nclusters(state::BondConfiguration; workspace::BondWorkspace=BondWorkspace(state)) =
    cluster_statistics(state; workspace=workspace).cluster_count

"""
    cluster_statistics(state; workspace=BondWorkspace(state))

Measure `occupied_bonds`, `cluster_count`, `largest_cluster`, and
`cluster_size_squared_sum = sum(size^2)` over **all** FK components, including
isolated vertices and the largest cluster. Counts are integers, not normalized
densities. The physical occupations and random streams are untouched. Only
scratch arrays are overwritten. Integer overflow is rejected explicitly.
"""
function cluster_statistics(state::BondConfiguration; workspace::BondWorkspace=BondWorkspace(state))
    labels = _bond_components(state; workspace=workspace)
    components, largest, squared_sum = 0, 0, 0
    for i in eachindex(labels)
        if labels[i] == i
            size = workspace.sizes[i]
            components += 1
            largest = max(largest, size)
            squared_sum = Base.checked_add(squared_sum, Base.checked_mul(size, size))
        end
    end
    (occupied_bonds=count(state.occupied), cluster_count=components,
     largest_cluster=largest, cluster_size_squared_sum=squared_sum)
end

"""
    sweeny_logratio(occupied, endpoints_disconnected, p, q)

Pure log weight ratio for toggling one edge in `p^b (1-p)^(m-b) q^k`.
`endpoints_disconnected` describes connectivity **with that edge removed**,
even for an insertion. Adding a disconnected edge merges two components:
the log ratio is `log(p/(1-p)) - log(q)`; otherwise it is `log(p/(1-p))`.
Deletion reverses that sign. Requires finite `0<p<1`, `q>0` (noninteger allowed).

This function has no graph, RNG or mutation dependency. The caller determines
connectivity, calls `metropolis_accept` and applies the accepted toggle.
"""
function sweeny_logratio(occupied::Bool, endpoints_disconnected::Bool, p::Real, q::Real)
    isfinite(p) && 0 < p < 1 || throw(ArgumentError("Sweeny requires finite 0 < p < 1"))
    isfinite(q) && q > 0 || throw(ArgumentError("q must be finite and positive"))
    log_add_ratio = log(p) - log1p(-p) - (endpoints_disconnected ? log(q) : 0)
    occupied ? -log_add_ratio : log_add_ratio
end

# log of the proposed/current random-cluster weight. Skipping the chosen edge
# distinguishes bridges from edges whose endpoints remain connected.
function _sweeny_logratio(state::BondConfiguration, edge_index::Integer, p::Real, q::Real;
                         workspace::BondWorkspace=BondWorkspace(state))
    isfinite(p) && 0 < p < 1 || throw(ArgumentError("Sweeny requires finite 0 < p < 1"))
    isfinite(q) && q > 0 || throw(ArgumentError("q must be finite and positive"))
    checkbounds(state.occupied, edge_index)
    labels = _bond_components(state; skip=Int(edge_index), workspace=workspace)
    u, v = state.graph.edges[edge_index]
    disconnected = labels[u] != labels[v]
    sweeny_logratio(state.occupied[edge_index], disconnected, p, q)
end

"""
    sweeny_update!(rng, state, p, q; edge=nothing, workspace=BondWorkspace(state)) -> Bool

One lazy symmetric single-edge Metropolis proposal for the random-cluster measure
`p^b (1-p)^(m-b) q^k`, with `b` occupied bonds, `m` total bonds, and `k`
components including isolated vertices. `q > 0` may be noninteger. A uniformly
random edge is selected unless `edge` is supplied. The proposal is skipped
with probability 1/2, ensuring an aperiodic chain even at `p=1/2, q=1`.
Returns whether it flipped.
This reference implementation rebuilds connectivity for each proposal and is
therefore slow on large graphs. Requires `0 < p < 1`.
"""
function sweeny_update!(rng::AbstractRNG, state::BondConfiguration, p::Real, q::Real;
                        edge::Union{Nothing,Integer}=nothing,
                        workspace::BondWorkspace=BondWorkspace(state))
    isempty(state.occupied) && throw(ArgumentError("Sweeny update requires at least one edge"))
    i = isnothing(edge) ? rand(rng, eachindex(state.occupied)) : edge
    logratio = _sweeny_logratio(state, i, p, q; workspace=workspace)
    rand(rng) < 0.5 && return false
    if metropolis_accept(rng, logratio)
        state.occupied[i] = !state.occupied[i]
        return true
    end
    false
end
sweeny_update!(state::BondConfiguration, p::Real, q::Real; kwargs...) = sweeny_update!(Random.default_rng(), state, p, q; kwargs...)

"""
Attempt `m` uniformly random lazy single-edge updates with replacement; return
the number accepted. Reuse `workspace` across sweeps to retain scratch storage.
A full connectivity rebuild per proposal remains: roughly O(m*(n+m)) work per
sweep, apart from near-constant amortized union-find factors.
"""
function sweeny_sweep!(rng::AbstractRNG, state::BondConfiguration, p::Real, q::Real;
                       workspace::BondWorkspace=BondWorkspace(state))
    isempty(state.occupied) && throw(ArgumentError("Sweeny sweep requires at least one edge"))
    accepted = 0
    for _ in eachindex(state.occupied)
        accepted += sweeny_update!(rng, state, p, q; workspace=workspace)
    end
    accepted
end
sweeny_sweep!(state::BondConfiguration, p::Real, q::Real; kwargs...) =
    sweeny_sweep!(Random.default_rng(), state, p, q; kwargs...)
