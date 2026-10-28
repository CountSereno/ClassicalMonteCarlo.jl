export IsingWorm, worm_step!, worm_closed, worm_bond_count

"""
    IsingWorm(L, beta; J=1.0, origin=1)

High-temperature graph worm on an L×L periodic square lattice, with L ≥ 3,
zero field and ferromagnetic J > 0. Edges carry occupation 0 or 1 with
weight tanh(beta*J)^number_of_occupied_edges. `head` and `tail` mark odd-degree
vertices; both equal `origin` initially. A step proposes one head move and
bond toggle. `beta` must be positive and finite, with tanh(beta*J) < 1.
This is a graph ensemble, not a spin configuration. The observable and time
conventions must be chosen for that ensemble.
"""
mutable struct IsingWorm
    L::Int
    x::Float64
    bonds::BitVector
    head::Int
    tail::Int
end

function IsingWorm(L::Integer, beta::Real; J::Real=1.0, origin::Integer=1)
    L >= 3 || throw(ArgumentError("L must be at least 3"))
    isfinite(beta) && beta > 0 && isfinite(J) && J > 0 ||
        throw(ArgumentError("beta and J must be finite and positive"))
    1 <= origin <= L*L || throw(ArgumentError("origin out of bounds"))
    x = tanh(Float64(beta)*Float64(J))
    0 < x < 1 || throw(ArgumentError("tanh(beta*J) must be representable strictly between 0 and 1"))
    return IsingWorm(Int(L), x, falses(2*L*L), Int(origin), Int(origin))
end

worm_closed(w::IsingWorm) = w.head == w.tail
worm_bond_count(w::IsingWorm) = count(w.bonds)

"""Propose one local worm move; return whether the move was accepted."""
function worm_step!(rng::AbstractRNG, w::IsingWorm)
    L = w.L
    site = w.head
    row, col = divrem(site - 1, L)
    direction = rand(rng, 1:4)
    if direction == 1 # positive column, forward horizontal edge
        next = row*L + mod(col+1, L) + 1
        edge = 2*site - 1
    elseif direction == 2 # negative column
        next = row*L + mod(col-1, L) + 1
        edge = 2*next - 1
    elseif direction == 3 # positive row, forward vertical edge
        next = mod(row+1, L)*L + col + 1
        edge = 2*site
    else
        next = mod(row-1, L)*L + col + 1
        edge = 2*next
    end
    # Adding a bond costs a factor x; removing one has acceptance 1.
    if w.bonds[edge] || rand(rng) < w.x
        w.bonds[edge] = !w.bonds[edge]
        w.head = next
        return true
    end
    return false
end
