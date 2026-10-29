using Test, Random

# Independent graph traversal: never calls production union-find or measurements.
function traversal_component_sizes(graph, occupied)
    adjacency = [Int[] for _ in 1:graph.nvertices]
    for (i, (u, v)) in enumerate(graph.edges)
        if occupied[i]
            push!(adjacency[u], v)
            push!(adjacency[v], u)
        end
    end
    visited = falses(graph.nvertices)
    sizes = Int[]
    for start in 1:graph.nvertices
        visited[start] && continue
        visited[start] = true
        stack, size = [start], 0
        while !isempty(stack)
            vertex = pop!(stack)
            size += 1
            for neighbour in adjacency[vertex]
                if !visited[neighbour]
                    visited[neighbour] = true
                    push!(stack, neighbour)
                end
            end
        end
        push!(sizes, size)
    end
    sort!(sizes; rev=true)
end

@testset "bond scratch storage and independent connectivity reference" begin
    graphs = [BondGraph(4, Tuple{Int,Int}[]),
              BondGraph(5, [(1,2), (2,3), (1,3)]),
              square_bond_graph(2)]
    for graph in graphs
        workspace = BondWorkspace(graph)
        m = length(graph.edges)
        for mask in 0:(1 << m)-1
            occupied = BitVector(!iszero(mask & (1 << (i-1))) for i in 1:m)
            state = BondConfiguration(graph, occupied)
            reference = traversal_component_sizes(graph, occupied)
            # Stale scratch contents must not enter a later connectivity query.
            fill!(workspace.parent, -9)
            fill!(workspace.sizes, -9)
            stats = cluster_statistics(state; workspace=workspace)
            @test stats == (occupied_bonds=count(occupied),
                cluster_count=length(reference), largest_cluster=maximum(reference),
                cluster_size_squared_sum=sum(abs2, reference))
            @test cluster_sizes(state; workspace=workspace) == reference
            @test nclusters(state; workspace=workspace) == length(reference)
            @test state.occupied == occupied
            for i in 1:m
                proposed = copy(occupied)
                proposed[i] = !proposed[i]
                delta_k = length(traversal_component_sizes(graph, proposed)) - length(reference)
                delta_b = occupied[i] ? -1 : 1
                p, q = .37, 1.5
                exact_logratio = delta_b*(log(p)-log1p(-p)) + delta_k*log(q)
                @test ClassicalMonteCarlo._sweeny_logratio(state,i,p,q;workspace=workspace) ≈ exact_logratio
                # A bridge deletion or component-merging insertion has |Δk|=1.
                @test sweeny_logratio(occupied[i],delta_k!=0,p,q) ≈ exact_logratio
            end
        end
    end
    state = BondConfiguration(square_bond_graph(3))
    @test_throws ArgumentError cluster_statistics(state;workspace=BondWorkspace(BondGraph(2, [(1,2)])))
    same = ones(Int, 9)
    @test_throws ArgumentError cluster_statistics(state;workspace=BondWorkspace(same,same))
    for p in (0,1,NaN,Inf), occupied in (false,true)
        @test_throws ArgumentError sweeny_logratio(occupied,false,p,2)
    end
    for q in (0,-1,NaN,Inf)
        @test_throws ArgumentError sweeny_logratio(false,true,.5,q)
    end
    @test sweeny_logratio(false,false,.5,3) == 0
    @test sweeny_logratio(true,true,.5,3) ≈ log(3)
    @test sweeny_logratio(false,true,.5,3) ≈ -log(3)
end

@testset "cached Sweeny updates preserve trajectory and future RNG" begin
    for (p,q) in ((.37,1.5),(.5,1.0),(.63,3.0))
        graph = square_bond_graph(4;periodic=true)
        fresh, cached = BondConfiguration(graph), BondConfiguration(graph)
        a, b = Xoshiro(417), Xoshiro(417)
        workspace = BondWorkspace(cached)
        for _ in 1:64
            accepted_a = sweeny_sweep!(a,fresh,p,q)
            accepted_b = sweeny_sweep!(b,cached,p,q;workspace=workspace)
            @test accepted_a == accepted_b
            @test fresh.occupied == cached.occupied
            # Interspersing measurements changes scratch only.
            @test cluster_statistics(cached;workspace=workspace) == cluster_statistics(fresh)
        end
        @test rand(a,UInt64,16) == rand(b,UInt64,16)
    end
end
