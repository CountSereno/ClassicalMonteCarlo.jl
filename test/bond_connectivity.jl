# Independent single-source traversal, separate from production queues and UF.
function reference_connected(state, skip)
    u, target = state.graph.edges[skip]
    seen, stack = Set([u]), [u]
    while !isempty(stack)
        vertex = pop!(stack)
        for (edge, (a,b)) in enumerate(state.graph.edges)
            (edge == skip || !state.occupied[edge]) && continue
            neighbour = a == vertex ? b : b == vertex ? a : 0
            neighbour == target && return true
            if neighbour != 0 && neighbour ∉ seen
                push!(seen, neighbour)
                push!(stack, neighbour)
            end
        end
    end
    false
end

@testset "interleaved connectivity: exhaustive graphs and visit stamps" begin
    graphs = [BondGraph(5, [(i,j) for i in 1:5 for j in i+1:5]),
              BondGraph(7, [(1,2),(2,3),(1,3),(3,4),(5,6)]), square_bond_graph(2)]
    for graph in graphs
        state = BondConfiguration(graph)
        workspace = SweenyWorkspace(state)
        m = length(state.occupied)
        for mask in 0:2^m-1
            state.occupied .= [!iszero(mask & (1 << (i-1))) for i in 1:m]
            before = copy(state.occupied)
            for edge in 1:m
                @test endpoints_connected(state,edge;workspace) == reference_connected(state,edge)
            end
            @test state.occupied == before
        end
        workspace.stamp = typemax(Int)
        fill!(workspace.marks, typemax(Int))
        @test endpoints_connected(state,1;workspace) == reference_connected(state,1)
        @test workspace.stamp == 1
        @test_throws BoundsError endpoints_connected(state,0;workspace)
        @test_throws ArgumentError endpoints_connected(BondConfiguration(deepcopy(graph)),1;workspace)
    end
end

@testset "IBFS and original Sweeny decisions are pathwise identical" begin
    for L in (3,8,16), (p,q) in ((.2,.3),(.5,1.0),(.6,1.5),(.72,3.0)), initial in (false,true)
        graph = square_bond_graph(L;periodic=true)
        a,b = BondConfiguration(graph),BondConfiguration(graph)
        fill!(a.occupied,initial); fill!(b.occupied,initial)
        ra,rb = Xoshiro(34),Xoshiro(34)
        wa,wb = BondWorkspace(a),SweenyWorkspace(b)
        for _ in 1:500
            # Preserve the old ordering explicitly: graph query precedes lazy draw.
            edge = rand(ra,eachindex(a.occupied))
            ratio = ClassicalMonteCarlo._sweeny_logratio(a,edge,p,q;workspace=wa)
            accepted = !(rand(ra)<.5) && metropolis_accept(ra,ratio)
            accepted && (a.occupied[edge] = !a.occupied[edge])
            @test sweeny_update!(rb,b,p,q;workspace=wb) == accepted
            @test a.occupied == b.occupied
        end
        @test rand(ra,UInt64,16) == rand(rb,UInt64,16)
        @test cluster_statistics(a) == cluster_statistics(b;workspace=wb.measurements)
    end
    state=BondConfiguration(square_bond_graph(3)); w=SweenyWorkspace(state)
    for p in (0,1,NaN,Inf), q in (.5,1.,3.)
        @test_throws ArgumentError sweeny_update!(Xoshiro(1),state,p,q;workspace=w)
    end
    for q in (0,-1,NaN,Inf)
        @test_throws ArgumentError sweeny_update!(Xoshiro(1),state,.5,q;workspace=w)
    end
    # q=1 bypasses connectivity completely (even for non-lazy moves).
    sweeny_sweep!(Xoshiro(4),state,.5,1.;workspace=w)
    @test w.stamp == 0
end

@testset "IBFS checkpoints validate topology and resume scratch safely" begin
    mktempdir() do dir
        spec=RandomClusterRunSpec(L=8,p=.6,q=2.3,seed=53,burn=5,samples=30,stride=2)
        direct=run!(start_run(spec)); partial=start_run(spec)
        advance!(partial,18)
        path=joinpath(dir,"state.checkpoint")
        save_checkpoint(path,partial); resumed=run!(load_checkpoint(path))
        @test observations(resumed)==observations(direct)
        @test resumed.model.occupied==direct.model.occupied
        @test rand(resumed.rng,UInt64,16)==rand(direct.rng,UInt64,16)
        @test run_metadata(resumed)["connectivity_backend"]=="interleaved_bfs"
        bad=start_run(spec); bad.kernel.workspace.edge_indices[1]=0
        @test_throws ArgumentError save_checkpoint(path,bad)
        bad=start_run(spec); bad.kernel.workspace.marks[1]=99
        @test_throws ArgumentError save_checkpoint(path,bad)
        bad=start_run(spec); bad.kernel.workspace.left_queue=bad.kernel.workspace.marks
        @test_throws ArgumentError save_checkpoint(path,bad)
        old=joinpath(dir,"old")
        write(old,"ClassicalMonteCarlo checkpoint 2\n")
        @test_throws ArgumentError load_checkpoint(old)
    end
end
