using Test, Statistics
include(joinpath(@__DIR__, "../experiments/multimodel/exact.jl"))
include(joinpath(@__DIR__, "../experiments/foundation/exact.jl"))
using .ExactModels

# Direct deletion of the original paired rows is independent of the optimized
# block-moment implementation. Test nonlinear errors, not only raw averages.
function check_paired_jackknife(statistic, values, result; blocksize)
    count = length(first(values))
    nblocks = count÷blocksize
    raw = statistic(values...)
    replicas = [begin
        keep = setdiff(1:count, ((block-1)*blocksize+1):(block*blocksize))
        statistic((series[keep] for series in values)...)
    end for block in 1:nblocks]
    bias = (nblocks-1)*(mean(replicas)-raw)
    stderr = sqrt((nblocks-1)/nblocks*sum((replicas .- mean(replicas)).^2))
    @test result.raw ≈ raw atol=1e-12
    @test result.bias ≈ bias atol=1e-12
    @test result.estimate ≈ raw-bias atol=1e-12
    @test result.stderr ≈ stderr atol=1e-12
end

@testset "paired Potts moments retain model definitions" begin
    E = [-18., -12, -10, -8, -18, -12, -6, -10]
    order = [1., .5, .25, .25, 1, .5, .25, .5]
    simplex = [1., .25, .0625, .25, 1, .25, .0625, .5]
    beta, N = .7, 9
    values = (E, order, simplex)
    result = potts_observables(E, order; beta, nsites=N, blocksize=2, simplex_squared=simplex)
    @test (result.n_samples, result.n_blocks, result.blocksize) == (8,4,2)
    definitions = (
        :energy_density => (e,m,s) -> mean(e)/N,
        :order_parameter => (e,m,s) -> mean(m),
        :heat_capacity => (e,m,s) -> beta^2*var(e; corrected=false)/N,
        :order_parameter_fluctuation => (e,m,s) -> beta*N*var(m; corrected=false),
        :moment_ratio => (e,m,s) -> mean(m.^4)/mean(m.^2)^2,
        :simplex_second_moment => (e,m,s) -> mean(s),
    )
    for (name, definition) in definitions
        check_paired_jackknife(definition, values, getproperty(result,name); blocksize=2)
    end
    shifted = potts_observables(E .+ 1e12, order; beta, nsites=N, blocksize=2)
    @test shifted.heat_capacity.raw ≈ result.heat_capacity.raw
    @test shifted.heat_capacity.stderr ≈ result.heat_capacity.stderr
    @test shifted.order_parameter_fluctuation == result.order_parameter_fluctuation
    @test !hasproperty(shifted, :simplex_second_moment)
    @test result.simplex_second_moment.raw != mean(order.^2)

    frozen = potts_observables(fill(-18.,8), ones(8); beta, nsites=N, blocksize=2, simplex_squared=ones(8))
    @test frozen.energy_density.raw == -2
    @test frozen.order_parameter.raw == frozen.moment_ratio.raw == 1
    @test frozen.heat_capacity.raw == frozen.heat_capacity.stderr == 0
    @test frozen.order_parameter_fluctuation.raw == frozen.order_parameter_fluctuation.stderr == 0
    @test frozen.simplex_second_moment.raw == 1
    @test_throws ArgumentError potts_observables(E, zeros(8); beta, nsites=N)
    @test_throws ArgumentError potts_observables([0.,0], [0.,1]; beta, nsites=N)
    @test_throws ArgumentError potts_observables(E, order[1:7]; beta, nsites=N)
    @test_throws ArgumentError potts_observables(E, order; beta=-1, nsites=N)
    @test_throws ArgumentError potts_observables(E, order; beta=Inf, nsites=N)
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=0)
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=N, blocksize=0)
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=N, blocksize=3)
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=N, blocksize=8)
    @test_throws ArgumentError potts_observables(E, 2order; beta, nsites=N)
    @test_throws ArgumentError potts_observables([Inf; E[2:end]], order; beta, nsites=N)
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=N, simplex_squared=simplex[1:7])
    @test_throws ArgumentError potts_observables(E, order; beta, nsites=N, simplex_squared=2simplex)
end

@testset "paired random-cluster geometric moments" begin
    # Empty, one bond, adjacent pair, disjoint pair, path, cycle, edge, empty.
    b = [0,1,2,2,3,4,1,0]
    k = [4,3,2,2,1,1,3,4]
    c = [1,2,3,2,4,4,2,1]
    s2 = [4,6,10,8,16,16,6,4]
    N, edges = 4,4
    values = (b,k,c,s2)
    result = random_cluster_observables(values...; nedges=edges,nvertices=N,blocksize=2)
    definitions = (
        :bond_density => (b,k,c,s2) -> mean(b)/edges,
        :cluster_density => (b,k,c,s2) -> mean(k)/N,
        :largest_cluster_fraction => (b,k,c,s2) -> mean(c)/N,
        :mean_cluster_size => (b,k,c,s2) -> mean(s2)/N,
        :bond_variance_per_edge => (b,k,c,s2) -> var(b; corrected=false)/edges,
    )
    for (name, definition) in definitions
        check_paired_jackknife(definition, values, getproperty(result,name); blocksize=2)
    end
    @test (result.n_samples,result.n_blocks,result.blocksize) == (8,4,2)
    empty = random_cluster_observables(zeros(4),fill(N,4),ones(4),fill(N,4);nedges=edges,nvertices=N)
    @test empty.bond_density.raw == empty.bond_variance_per_edge.raw == 0
    @test empty.mean_cluster_size.raw == empty.cluster_density.raw == 1
    full = random_cluster_observables(fill(edges,4),ones(4),fill(N,4),fill(N^2,4);nedges=edges,nvertices=N)
    @test full.bond_density.raw == full.largest_cluster_fraction.raw == 1
    @test full.mean_cluster_size.raw == N
    @test full.bond_variance_per_edge.raw == full.bond_variance_per_edge.stderr == 0
    @test_throws ArgumentError random_cluster_observables(b,k,c,s2;nedges=0,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b,k,c,s2;nedges=edges,nvertices=0)
    @test_throws ArgumentError random_cluster_observables(b[1:7],k,c,s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b .+ .5,k,c,s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(2b,k,c,s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b,zeros(8),c,s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b,k,zeros(8),s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b,k,c,2s2;nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(zeros(8),ones(8),fill(N,8),fill(N^2,8);nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(fill(2,4),fill(2,4),fill(2,4),fill(7,4);nedges=edges,nvertices=N)
    @test_throws ArgumentError random_cluster_observables(b,k,c,s2;nedges=edges,nvertices=N,blocksize=8)
end

@testset "independent Potts and random-cluster exact references" begin
    edges = square_edges(2; periodic=false)
    @test length(edges) == 4
    @test length(square_edges(3)) == 18
    @test isempty(square_edges(1;periodic=false))
    @test_throws ArgumentError square_edges(2)
    for q in (2,3)
        free = enumerate_potts(3,q,0)
        frozen = enumerate_potts(3,q,100)
        @test free["states"] == q^9
        @test free["energy_density"] ≈ -2/q atol=1e-13
        @test free["simplex_second_moment"] ≈ 1/9 atol=1e-13
        @test free["log_partition"] ≈ 9log(q)
        @test free["heat_capacity"] == free["order_parameter_fluctuation"] == 0
        @test frozen["energy_density"] ≈ -2
        @test frozen["order_parameter"] ≈ 1
        @test frozen["moment_ratio"] ≈ 1
        @test frozen["simplex_second_moment"] ≈ 1
    end
    # The two-state Potts Hamiltonian with J=1 is H_Ising/2-N on this
    # periodic square graph. Its temperature convention must not be reused
    # unchanged as an Ising beta; this also distinguishes fluctuation units.
    beta_two = .7
    two_colour = enumerate_potts(3,2,beta_two)
    ising = ExactIsing.enumerate_ising(3,beta_two/2)
    @test two_colour["energy_density"] ≈ ising["energy_density"]/2-1
    @test two_colour["heat_capacity"] ≈ ising["heat_capacity"]
    @test two_colour["order_parameter"] ≈ ising["abs_magnetization_density"]
    @test two_colour["order_parameter_fluctuation"] ≈ 2ising["susceptibility_abs"]
    @test two_colour["moment_ratio"] ≈ 3*(1-ising["binder"])
    @test two_colour["simplex_second_moment"] ≈ ising["susceptibility"]/(9beta_two/2)
    for p in (0.,.13,.5,1.)
        independent = enumerate_random_cluster(4,edges,p,1.)
        @test independent["states"] == 16
        @test independent["bond_density"] ≈ p atol=1e-14
        @test independent["bond_variance_per_edge"] ≈ p*(1-p) atol=1e-14
        @test independent["log_partition"] ≈ 0 atol=1e-14
    end
    for q in (.3,1.,1.5,3.)
        empty = enumerate_random_cluster(4,edges,0,q)
        full = enumerate_random_cluster(4,edges,1,q)
        @test empty["cluster_density"] == 1
        @test empty["mean_cluster_size"] == 1
        @test empty["largest_cluster_fraction"] == .25
        @test full["cluster_density"] == .25
        @test full["mean_cluster_size"] == 4
        @test full["largest_cluster_fraction"] == 1
        @test empty["log_partition"] ≈ 4log(q)
        @test full["log_partition"] ≈ log(q)
    end
    # On a tree k=N-b, so even noninteger q gives independent occupied bonds
    # with probability p/(p+q*(1-p)). This is an exact noninteger-q check.
    p,q = .37,1.5
    effective_p = p/(p+q*(1-p))
    tree = enumerate_random_cluster(4,[(1,2),(2,3),(3,4)],p,q)
    @test tree["bond_density"] ≈ effective_p
    @test tree["bond_variance_per_edge"] ≈ effective_p*(1-effective_p)
    @test tree["mean_clusters"] ≈ 4-3effective_p
    # The FK identity independently connects colour enumeration with BFS bond
    # enumeration. q=3 uses 19,683 spins versus 262,144 distinct bond subsets.
    beta = .7
    potts = enumerate_potts(3,3,beta)
    bond_probability = -expm1(-beta)
    fk = enumerate_random_cluster(9,square_edges(3),bond_probability,3.)
    @test potts["energy_density"] ≈ -fk["mean_bonds"]/(bond_probability*9) rtol=1e-12
    @test potts["simplex_second_moment"] ≈ fk["mean_cluster_size"]/9 rtol=1e-12
    @test potts["log_partition"] ≈ fk["log_partition"] + 18beta rtol=1e-12
    @test_throws ArgumentError enumerate_potts(2,3,.5)
    @test_throws ArgumentError enumerate_potts(3,1,.5)
    @test_throws ArgumentError enumerate_potts(3,3,-.5)
    @test_throws ArgumentError enumerate_potts(8,3,.5)
    @test_throws ArgumentError enumerate_random_cluster(4,edges,.5,0)
    @test_throws ArgumentError enumerate_random_cluster(4,edges,-.1,1)
    @test_throws ArgumentError enumerate_random_cluster(4,[(1,2),(2,1)],.5,1)
    @test_throws ArgumentError enumerate_random_cluster(4,[(1,1)],.5,1)
    @test_throws ArgumentError enumerate_random_cluster(4,[(1,5)],.5,1)
    @test_throws ArgumentError enumerate_random_cluster(1,Tuple{Int,Int}[],.5,1)
end
