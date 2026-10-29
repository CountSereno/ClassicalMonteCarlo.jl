# Small, model-independent disjoint-set primitives. Callers own the arrays;
# no lattice, RNG, physical state, allocation or file I/O is hidden here.
# These routines are internal: callers validate dimensions before hot loops.
function _reset_sets!(parent, sizes)
    for i in eachindex(parent)
        parent[i] = i
        sizes[i] = 1
    end
    nothing
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
