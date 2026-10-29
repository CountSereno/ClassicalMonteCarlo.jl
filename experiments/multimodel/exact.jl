"""
Independent small-system reference distributions for Potts and random-cluster
experiments. No production model, energy, connectivity or estimator code is
called: spin energies use explicit coordinates and bond connectivity uses BFS.
Only standard libraries are needed; exhaustive state counts are capped.
"""
module ExactModels

export enumerate_potts, enumerate_random_cluster, square_edges

"""Independent right/up square-lattice edge list; periodic L >= 3, open L >= 1."""
function square_edges(L::Integer; periodic::Bool=true)
    L >= (periodic ? 3 : 1) || throw(ArgumentError("L is too small for this boundary"))
    edges = Tuple{Int,Int}[]
    index(x,y) = x + L*(y-1)
    for y in 1:L, x in 1:L
        if x < L
            push!(edges, (index(x,y), index(x+1,y)))
        elseif periodic
            push!(edges, (index(x,y), index(1,y)))
        end
        if y < L
            push!(edges, (index(x,y), index(x,y+1)))
        elseif periodic
            push!(edges, (index(x,y), index(x,1)))
        end
    end
    edges
end

function state_count(base, exponent; limit=2_000_000)
    states = 1
    for _ in 1:exponent
        states <= limit÷base || throw(ArgumentError("exact enumeration exceeds $limit states"))
        states *= base
    end
    states
end

function normalized_weights(logweights)
    shift = maximum(logweights)
    isfinite(shift) || throw(ArgumentError("no finite-weight state in exact ensemble"))
    weights = exp.(logweights .- shift)
    partition = sum(weights)
    weights ./= partition
    weights, shift + log(partition)
end

"""
    enumerate_potts(L, q, beta; J=1)

Exhaust a periodic square, zero-field integer-q Potts model. Returns population
observables in the same normalizations as `potts_observables`, plus configuration
and state-count metadata. `simplex_second_moment` uses all colour populations,
whereas `moment_ratio` uses the normalized largest-colour population. At most
two million configurations are permitted (e.g. L=3, q=3 has 19,683 states).
"""
function enumerate_potts(L::Integer, q::Integer, beta::Real; J::Real=1)
    L >= 3 || throw(ArgumentError("periodic Potts reference requires L >= 3"))
    q >= 2 || throw(ArgumentError("q must be at least two"))
    b, coupling = Float64(beta), Float64(J)
    isfinite(b) && b >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    isfinite(coupling) || throw(ArgumentError("J must be finite"))
    N = Base.checked_mul(Int(L), Int(L))
    count = state_count(q, N)
    energies = Vector{Float64}(undef, count)
    orders = similar(energies)
    simplex = similar(energies)
    colours = zeros(Int, N)
    populations = zeros(Int, q)
    for code in 0:(count-1)
        value = code
        fill!(populations, 0)
        for site in 1:N
            colour = value % q + 1
            value ÷= q
            colours[site] = colour
            populations[colour] += 1
        end
        equal_pairs = 0
        for y in 1:L, x in 1:L
            site = x + L*(y-1)
            right = mod(x,L)+1 + L*(y-1)
            up = x + L*mod(y,L)
            equal_pairs += colours[site] == colours[right]
            equal_pairs += colours[site] == colours[up]
        end
        energies[code+1] = -coupling*equal_pairs
        orders[code+1] = (q*maximum(populations)/N-1)/(q-1)
        simplex[code+1] = (q*sum(abs2, populations)/N^2-1)/(q-1)
    end
    all(isfinite, energies) || throw(ArgumentError("energy overflow in exact enumeration"))
    # Subtract the minimum before multiplying by beta to avoid needless
    # Boltzmann overflow at low temperature. The physical log Z is restored.
    ground = minimum(energies)
    weights, shifted_logZ = normalized_weights(-b .* (energies .- ground))
    average(values) = sum(weights .* values)
    energy_mean, order_mean = average(energies), average(orders)
    order2 = average(orders.^2)
    order2 > 0 || throw(ArgumentError("Potts order moment ratio is undefined"))
    Dict{String,Any}(
        "model"=>"Potts", "L"=>Int(L), "q"=>Int(q), "beta"=>b, "J"=>coupling,
        "nvertices"=>N, "nedges"=>2N, "states"=>count,
        "log_partition"=>shifted_logZ-b*ground,
        "energy_density"=>energy_mean/N,
        "order_parameter"=>order_mean,
        "heat_capacity"=>b^2*average((energies .- energy_mean).^2)/N,
        "order_parameter_fluctuation"=>b*N*average((orders .- order_mean).^2),
        "moment_ratio"=>average(orders.^4)/order2^2,
        "simplex_second_moment"=>average(simplex))
end

function checked_edges(nvertices, edges)
    nvertices > 0 || throw(ArgumentError("nvertices must be positive"))
    pairs = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for edge in edges
        length(edge) == 2 || throw(ArgumentError("each edge needs two endpoints"))
        u, v = Int(edge[1]), Int(edge[2])
        1 <= u <= nvertices && 1 <= v <= nvertices || throw(ArgumentError("endpoint outside graph"))
        u != v || throw(ArgumentError("self loops are unsupported"))
        pair = minmax(u,v)
        pair ∉ seen || throw(ArgumentError("parallel edges are unsupported"))
        push!(seen, pair)
        push!(pairs, pair)
    end
    pairs
end

"""
    enumerate_random_cluster(nvertices, edges, p, q)

Enumerate every bond subset of a finite simple graph with probabilities
proportional to `p^b*(1-p)^(m-b)*q^k`. Connectivity is evaluated by independent
breadth-first searches. Real `q>0` and `0<=p<=1` are supported, including exact
endpoint distributions. The graph must have at least one edge; at most two
million subsets are permitted. No spin-energy or temperature interpretation is
assigned to its geometric moments.
"""
function enumerate_random_cluster(nvertices::Integer, edges, p::Real, q::Real)
    probability, cluster_weight = Float64(p), Float64(q)
    isfinite(probability) && 0 <= probability <= 1 || throw(ArgumentError("p must be finite and in [0,1]"))
    isfinite(cluster_weight) && cluster_weight > 0 || throw(ArgumentError("q must be finite and positive"))
    pairs = checked_edges(nvertices, edges)
    N, m = Int(nvertices), length(pairs)
    m > 0 || throw(ArgumentError("bond-density reference requires at least one edge"))
    count = state_count(2, m)
    adjacency = [Tuple{Int,Int}[] for _ in 1:N]
    for (index, (u,v)) in enumerate(pairs)
        push!(adjacency[u], (v,index))
        push!(adjacency[v], (u,index))
    end
    bonds = Vector{Int}(undef, count)
    clusters = similar(bonds)
    largest = similar(bonds)
    squared_sum = similar(bonds)
    logweights = Vector{Float64}(undef, count)
    visited, queue = falses(N), zeros(Int,N)
    logq = log(cluster_weight)
    for subset in 0:(count-1)
        b = count_ones(subset)
        fill!(visited, false)
        k, maximum_size, s2 = 0, 0, 0
        for start in 1:N
            visited[start] && continue
            k += 1
            visited[start] = true
            queue[1] = start
            first, last = 1, 1
            while first <= last
                vertex = queue[first]
                first += 1
                for (neighbour, edge) in adjacency[vertex]
                    if !visited[neighbour] && !iszero(subset & (1 << (edge-1)))
                        visited[neighbour] = true
                        last += 1
                        queue[last] = neighbour
                    end
                end
            end
            maximum_size = max(maximum_size, last)
            s2 += last^2
        end
        i = subset+1
        bonds[i], clusters[i], largest[i], squared_sum[i] = b, k, maximum_size, s2
        logweights[i] = if probability == 0
            b == 0 ? k*logq : -Inf
        elseif probability == 1
            b == m ? k*logq : -Inf
        else
            b*log(probability) + (m-b)*log1p(-probability) + k*logq
        end
    end
    weights, logZ = normalized_weights(logweights)
    average(values) = sum(weights .* values)
    mean_bonds, mean_clusters = average(bonds), average(clusters)
    Dict{String,Any}(
        "model"=>"RandomCluster", "p"=>probability, "q"=>cluster_weight,
        "nvertices"=>N, "nedges"=>m, "states"=>count, "log_partition"=>logZ,
        "mean_bonds"=>mean_bonds, "mean_clusters"=>mean_clusters,
        "bond_density"=>mean_bonds/m,
        "cluster_density"=>mean_clusters/N,
        "largest_cluster_fraction"=>average(largest)/N,
        "mean_cluster_size"=>average(squared_sum)/N,
        "bond_variance_per_edge"=>average((bonds .- mean_bonds).^2)/m)
end

end # module ExactModels
