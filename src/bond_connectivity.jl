export SweenyWorkspace, endpoints_connected

"""
    SweenyWorkspace(graph_or_state)

Reusable interleaved breadth-first-search (IBFS) scratch space. Static adjacency
uses contiguous edge-index and neighbour arrays; occupations remain owned by
the configuration. Two queues explore opposite endpoints, stopping when they
meet or either component is exhausted. Signed visit stamps avoid clearing an
N-element array for each query. No connectivity survives between queries.

The graph must be the same object and its topology must not be mutated while
this workspace is in use. Occupations may change freely. Each concurrent chain
needs its own workspace. `measurements` holds independent union-find scratch
for full cluster statistics. This is not a dynamic spanning-forest algorithm;
a query still has O(N+m) worst-case cost and storage is O(N+m).
"""
mutable struct SweenyWorkspace
    graph::BondGraph
    offsets::Vector{Int}
    neighbours::Vector{Int}
    edge_indices::Vector{Int}
    marks::Vector{Int}
    left_queue::Vector{Int}
    right_queue::Vector{Int}
    stamp::Int
    measurements::BondWorkspace
end

function SweenyWorkspace(graph::BondGraph)
    n, m = graph.nvertices, length(graph.edges)
    degrees = zeros(Int, n)
    for (u, v) in graph.edges
        degrees[u] += 1
        degrees[v] += 1
    end
    offsets = vcat(1, 1 .+ cumsum(degrees))
    neighbours, edge_indices = Vector{Int}(undef, 2m), Vector{Int}(undef, 2m)
    cursor = copy(offsets)
    for (edge, (u, v)) in enumerate(graph.edges)
        for (a, b) in ((u, v), (v, u))
            neighbours[cursor[a]] = b
            edge_indices[cursor[a]] = edge
            cursor[a] += 1
        end
    end
    SweenyWorkspace(graph, offsets, neighbours, edge_indices, zeros(Int, n),
                    Vector{Int}(undef, n), Vector{Int}(undef, n), 0, BondWorkspace(graph))
end
SweenyWorkspace(state::BondConfiguration) = SweenyWorkspace(state.graph)

function _check_workspace(w::SweenyWorkspace, state::BondConfiguration)
    w.graph === state.graph || throw(ArgumentError("Sweeny workspace belongs to a different graph"))
    n, m = state.graph.nvertices, length(state.graph.edges)
    length(w.offsets) == n+1 && length(w.neighbours) == length(w.edge_indices) == 2m &&
        length(w.marks) == length(w.left_queue) == length(w.right_queue) == n ||
        throw(ArgumentError("invalid Sweeny scratch dimensions"))
    0 <= w.stamp <= typemax(Int) || throw(ArgumentError("invalid visit stamp"))
    nothing
end

# Persistence boundaries validate topology and scratch aliasing once, outside
# the proposal loop. Scratch contents are otherwise immaterial to the physics.
function _validate_connectivity(w::SweenyWorkspace, state::BondConfiguration)
    _check_workspace(w, state)
    expected = SweenyWorkspace(state.graph)
    w.offsets == expected.offsets && w.neighbours == expected.neighbours &&
        w.edge_indices == expected.edge_indices || throw(ArgumentError("adjacency/topology mismatch"))
    arrays = (w.offsets, w.neighbours, w.edge_indices, w.marks, w.left_queue,
              w.right_queue, w.measurements.parent, w.measurements.sizes)
    all(arrays[i] !== arrays[j] for i in eachindex(arrays) for j in 1:i-1) ||
        throw(ArgumentError("connectivity scratch arrays must not alias"))
    _check_workspace(w.measurements, state)
    # Valid marks were written only in earlier generations. Otherwise a future
    # query could confuse corrupted checkpoint scratch with a visited vertex.
    all(x -> -w.stamp <= x <= w.stamp, w.marks) || throw(ArgumentError("invalid visit marks"))
    nothing
end
_validate_connectivity(w::BondWorkspace, state::BondConfiguration) = _check_workspace(w, state)
_measurement_workspace(w::BondWorkspace) = w
_measurement_workspace(w::SweenyWorkspace) = w.measurements

# Expand one vertex, then let the other queue take a turn. The signed stamp
# identifies the owning search without an extra ownership array.
@inline function _expand_front!(w, state, skip, vertex, queue, tail, stamp)
    for index in w.offsets[vertex]:(w.offsets[vertex+1]-1)
        edge = w.edge_indices[index]
        (edge == skip || !state.occupied[edge]) && continue
        neighbour = w.neighbours[index]
        mark = w.marks[neighbour]
        mark == -stamp && return tail, true
        if mark != stamp
            w.marks[neighbour] = stamp
            tail += 1
            queue[tail] = neighbour
        end
    end
    tail, false
end

function _endpoints_connected!(state, edge, w::SweenyWorkspace)
    if w.stamp == typemax(Int)
        fill!(w.marks, 0)
        w.stamp = 0
    end
    w.stamp += 1
    stamp = w.stamp
    u, v = state.graph.edges[edge]
    w.left_queue[1], w.right_queue[1] = u, v
    w.marks[u], w.marks[v] = stamp, -stamp
    left_head = right_head = left_tail = right_tail = 1
    while left_head <= left_tail && right_head <= right_tail
        left_tail, connected = _expand_front!(w, state, edge,
            w.left_queue[left_head], w.left_queue, left_tail, stamp)
        connected && return true
        left_head += 1
        left_head > left_tail && return false
        right_tail, connected = _expand_front!(w, state, edge,
            w.right_queue[right_head], w.right_queue, right_tail, -stamp)
        connected && return true
        right_head += 1
    end
    false
end

"""
    endpoints_connected(state, edge; workspace=SweenyWorkspace(state))

Are the endpoints connected when `edge` is excluded, regardless of its current
occupation? Handles bridge deletions and cycle edges without changing the
configuration or consuming randomness. A `BondWorkspace` selects the full
union-find rebuild reference. Reuse a workspace for repeated queries.
"""
function endpoints_connected(state::BondConfiguration, edge::Integer;
                             workspace::Union{BondWorkspace,SweenyWorkspace}=SweenyWorkspace(state))
    checkbounds(state.occupied, edge)
    _check_workspace(workspace, state)
    _endpoints_connected!(state, Int(edge), workspace)
end
function _endpoints_connected!(state, edge, workspace::BondWorkspace)
    labels = _bond_components(state; skip=edge, workspace)
    u, v = state.graph.edges[edge]
    labels[u] == labels[v]
end
