export potts_observables, random_cluster_observables

# These estimators only depend on the statistical helpers in analysis.jl. They
# do not construct models, consume random numbers or know about the runner.
# A row is one simultaneous observation; all moments in a time block must be
# removed together when propagating errors through a nonlinear statistic.
function _paired_moment_jackknife(names::NTuple{N,Symbol}, nsamples::Int,
                                  blocksize::Integer, nmoments::Int,
                                  moment_at, evaluate) where N
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    nsamples % blocksize == 0 || throw(ArgumentError("blocksize must divide the number of observations"))
    nblocks = nsamples ÷ blocksize
    nblocks >= 2 || throw(ArgumentError("at least two blocks are required"))
    moments = zeros(nblocks, nmoments)
    for block in 1:nblocks
        for i in ((block-1)*blocksize+1):(block*blocksize)
            values = moment_at(i)
            for j in 1:nmoments
                moments[block, j] += values[j]
            end
        end
    end
    moments ./= blocksize
    all(isfinite, moments) || throw(ArgumentError("moments overflowed; rescale the input units"))
    total = vec(sum(moments; dims=1))
    raw = evaluate(total/nblocks)
    all(isfinite, raw) || throw(ArgumentError("observable calculation is nonfinite"))
    replicas = Matrix{Float64}(undef, nblocks, N)
    for block in 1:nblocks
        values = evaluate((total - @view(moments[block, :]))/(nblocks-1))
        all(isfinite, values) || throw(ArgumentError("delete-block observable calculation is nonfinite"))
        replicas[block, :] .= values
    end
    statistics = ntuple(j -> _jackknife_result(raw[j], @view(replicas[:, j])), N)
    merge(NamedTuple{names}(statistics),
          (n_samples=nsamples, n_blocks=nblocks, blocksize=Int(blocksize)))
end

function _observable_count(value::Integer, name)
    value > 0 || throw(ArgumentError("$name must be positive"))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("$name must have a finite Float64 representation"))
    result
end

function _paired_observable_series(series...)
    values = map(x -> _finite_series(x, 2), series)
    all(x -> all(isfinite, x), values) || throw(ArgumentError("observations must have finite Float64 representations"))
    all(x -> length(x) == length(first(values)), values) || throw(ArgumentError("paired observations must have equal lengths"))
    values
end

"""
    potts_observables(E, order; beta, nsites, blocksize=1,
                      simplex_squared=nothing)

Statistics for paired stationary Potts observations. `E` is **total** energy;
`order` is the normalized largest-colour population
`(q*maximum(N_a)/N - 1)/(q-1)` in `[0,1]`, not a signed total magnetization.
Returns `energy_density`, `order_parameter`, `heat_capacity`,
`order_parameter_fluctuation`, and `moment_ratio`. Each field contains
`(estimate, stderr, bias, raw)`, followed by sampling/block metadata.

Raw definitions are `<E>/N`, `<order>`, `beta^2*Var(E)/N`,
`beta*N*Var(order)`, and `<order^4>/<order^2>^2`. The fluctuation is an
order-parameter proxy, **not** a magnetic linear-response susceptibility.
The moment ratio is not the scalar Ising Binder cumulant; its critical value
depends on the observable definition and geometry.

Optionally pass a paired `simplex_squared` series defined by
`(q*sum(N_a^2)/N^2 - 1)/(q-1)`. Its mean is returned as
`simplex_second_moment`. This is a different observable from `order^2` for
`q>2`. In the zero-field ferromagnetic FK representation its expectation is
`<sum(cluster_size^2)>/N^2` at `p=1-exp(-beta*J)`.

Every delete-block replica removes all paired moments together. Energy and
order fluctuations are centered before forming squared moments; no resampling
or equilibration decisions occur here. A full or delete-block sample with
`<order^2> == 0` has undefined moment ratio and raises `ArgumentError`.
Constant nonzero order and zero fluctuation are supported. Blocks must be
long compared with relevant correlations; compare block lengths and chains.
"""
function potts_observables(E::AbstractVector{<:Real}, order::AbstractVector{<:Real};
                           beta::Real, nsites::Integer, blocksize::Integer=1,
                           simplex_squared::Union{Nothing,AbstractVector{<:Real}}=nothing)
    energies, orders = _paired_observable_series(E, order)
    n = _observable_count(nsites, "nsites")
    b = Float64(beta)
    isfinite(b) && b >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    all(x -> 0 <= x <= 1, orders) || throw(ArgumentError("normalized Potts order must lie in [0,1]"))
    simplex = if isnothing(simplex_squared)
        nothing
    else
        values, = _paired_observable_series(simplex_squared)
        length(values) == length(energies) || throw(ArgumentError("simplex_squared must match the paired trace length"))
        all(x -> 0 <= x <= 1, values) || throw(ArgumentError("simplex_squared must lie in [0,1]"))
        values
    end
    energy_offset, order_offset = first(energies), first(orders)
    function moment_at(i)
        e = (energies[i] - energy_offset)/n
        m = orders[i]
        centered_order = m - order_offset
        base = (e, e^2, m, centered_order, centered_order^2, m^2, m^4)
        isnothing(simplex) ? base : (base..., simplex[i])
    end
    function evaluate(v)
        v[6] > 0 || throw(ArgumentError("Potts moment ratio is undefined for a full or delete-block sample with zero mean order^2"))
        base = (energy_offset/n + v[1], v[3],
                b^2*n*(v[2]-v[1]^2), b*n*(v[5]-v[4]^2), v[7]/v[6]^2)
        isnothing(simplex) ? base : (base..., v[8])
    end
    names = (:energy_density, :order_parameter, :heat_capacity,
             :order_parameter_fluctuation, :moment_ratio)
    if !isnothing(simplex)
        names = (names..., :simplex_second_moment)
    end
    _paired_moment_jackknife(names, length(energies), blocksize,
                            isnothing(simplex) ? 7 : 8, moment_at, evaluate)
end

"""
    random_cluster_observables(bonds, clusters, largest, size_squared_sum;
                               nedges, nvertices, blocksize=1)

Paired geometric statistics of a random-cluster graph. Each row contains the
number of occupied bonds `b`, the number of components `k` including isolated
vertices, the largest component size, and `sum(component_size^2)` including
the largest component. Inputs are unnormalised integer-valued counts.

Returns `bond_density=<b>/nedges`, `cluster_density=<k>/nvertices`,
`largest_cluster_fraction=<largest>/nvertices`,
`mean_cluster_size=<size_squared_sum>/nvertices`, and
`bond_variance_per_edge=Var(b)/nedges`, each with joint-block jackknife
`(estimate, stderr, bias, raw)` and sample/block metadata. `mean_cluster_size`
is the expected size of the cluster containing a uniformly selected vertex;
it includes the largest cluster and is not a finite-cluster-only susceptibility.
The bond variance is a geometric fluctuation, not a heat capacity.

Positive `nedges` is required because bond density is undefined for an edgeless
graph. Constant endpoint samples are allowed and have zero estimated variance.
Bounds and simple consistency conditions are checked; these checks do not
reconstruct the graph or establish equilibration. No `p`, `q`, spin energy or
temperature is silently inferred from geometric observations.
"""
function random_cluster_observables(bonds::AbstractVector{<:Real},
                                    clusters::AbstractVector{<:Real},
                                    largest::AbstractVector{<:Real},
                                    size_squared_sum::AbstractVector{<:Real};
                                    nedges::Integer, nvertices::Integer,
                                    blocksize::Integer=1)
    b, k, c, s2 = _paired_observable_series(bonds, clusters, largest, size_squared_sum)
    m = _observable_count(nedges, "nedges")
    n = _observable_count(nvertices, "nvertices")
    isfinite(n^2) || throw(ArgumentError("nvertices^2 must have a finite Float64 representation"))
    all(x -> isinteger(x) && 0 <= x <= m, b) || throw(ArgumentError("occupied bond counts must be integers in [0,nedges]"))
    all(x -> isinteger(x) && 1 <= x <= n, k) || throw(ArgumentError("cluster counts must be integers in [1,nvertices]"))
    all(x -> isinteger(x) && 1 <= x <= n, c) || throw(ArgumentError("largest cluster sizes must be integers in [1,nvertices]"))
    all(x -> isinteger(x) && n <= x <= n^2, s2) || throw(ArgumentError("squared size sums must be integers in [nvertices,nvertices^2]"))
    for i in eachindex(b)
        b[i] >= n-k[i] || throw(ArgumentError("too few occupied bonds for the component count"))
        c[i]+k[i]-1 <= n <= k[i]*c[i] || throw(ArgumentError("inconsistent cluster count and largest size"))
        c[i]^2+n-c[i] <= s2[i] <= n*c[i] || throw(ArgumentError("inconsistent largest size and squared size sum"))
        iszero(mod(s2[i]-n, 2)) || throw(ArgumentError("squared cluster sizes must have the parity of nvertices"))
    end
    offset = first(b)
    moment_at(i) = ((b[i]-offset)/m, ((b[i]-offset)/m)^2,
                    k[i]/n, c[i]/n, s2[i]/n)
    evaluate(v) = (offset/m+v[1], v[3], v[4], v[5], m*(v[2]-v[1]^2))
    names = (:bond_density, :cluster_density, :largest_cluster_fraction,
             :mean_cluster_size, :bond_variance_per_edge)
    _paired_moment_jackknife(names, length(b), blocksize, 5, moment_at, evaluate)
end
