# Model-specific adapters. The execution loop knows only how to step a kernel
# and append one paired observation; statistical inference stays outside it.
abstract type AbstractUpdateKernel end
struct MetropolisKernel <: AbstractUpdateKernel end
struct WolffKernel{W} <: AbstractUpdateKernel
    workspace::W
end
struct SwendsenWangKernel{W} <: AbstractUpdateKernel
    workspace::W
end
struct SweenyKernel{W} <: AbstractUpdateKernel
    workspace::W
end
struct IndependentBondKernel{W} <: AbstractUpdateKernel
    workspace::W
end

WolffKernel(model::Union{Ising,Potts}) = WolffKernel(WolffWorkspace(model))
SwendsenWangKernel(model::Union{Ising,Potts}) = SwendsenWangKernel(SwendsenWangWorkspace(model))
SweenyKernel(model::BondConfiguration) = SweenyKernel(BondWorkspace(model))
IndependentBondKernel(model::BondConfiguration) = IndependentBondKernel(BondWorkspace(model))

update_unit(::MetropolisKernel) = "metropolis_sweep"
update_unit(::WolffKernel) = "wolff_cluster"
update_unit(::SwendsenWangKernel) = "swendsen_wang_sweep"
update_unit(::SweenyKernel) = "sweeny_sweep"
update_unit(::IndependentBondKernel) = "independent_bond_redraw"
work_unit(::MetropolisKernel) = "accepted_flips"
work_unit(::WolffKernel) = "flipped_sites"
work_unit(::SwendsenWangKernel) = "formed_clusters"
work_unit(::SweenyKernel) = "accepted_bond_flips"
work_unit(::IndependentBondKernel) = "resampled_bonds"
_algorithm(::MetropolisKernel) = :metropolis
_algorithm(::WolffKernel) = :wolff
_algorithm(::SwendsenWangKernel) = :swendsen_wang
_algorithm(::SweenyKernel) = :sweeny
_algorithm(::IndependentBondKernel) = :independent

_initialize_model(rng, spec::RunSpec) = Ising(spec.L; J=spec.J, init=spec.initial, rng=rng)
_initialize_model(rng, spec::PottsRunSpec) = Potts(spec.L, spec.q; J=spec.J, init=spec.initial, rng=rng)
function _initialize_model(rng, spec::RandomClusterRunSpec)
    model = BondConfiguration(square_bond_graph(spec.L; periodic=spec.periodic))
    spec.initial == :full && fill!(model.occupied, true)
    spec.initial == :random && sample_bond_percolation!(rng, model, spec.p)
    model
end

function _make_kernel(model, spec::AbstractSpinRunSpec)
    spec.algorithm == :metropolis && return MetropolisKernel()
    spec.algorithm == :wolff && return WolffKernel(model)
    spec.algorithm == :swendsen_wang && return SwendsenWangKernel(model)
    throw(ArgumentError("unknown spin kernel"))
end
_make_kernel(model, spec::RandomClusterRunSpec) = spec.algorithm == :sweeny ? SweenyKernel(model) : IndependentBondKernel(model)

_runner_step!(rng, model, spec::AbstractSpinRunSpec, ::MetropolisKernel) =
    metropolis!(model, spec.beta; rng=rng)
_runner_step!(rng, model, spec::AbstractSpinRunSpec, k::WolffKernel) =
    wolff!(model, spec.beta; rng=rng, workspace=k.workspace)
_runner_step!(rng, model, spec::AbstractSpinRunSpec, k::SwendsenWangKernel) =
    swendsen_wang!(model, spec.beta; rng=rng, workspace=k.workspace)
_runner_step!(rng, model, spec::RandomClusterRunSpec, kernel::SweenyKernel) =
    sweeny_sweep!(rng, model, spec.p, spec.q; workspace=kernel.workspace)
function _runner_step!(rng, model, spec::RandomClusterRunSpec, ::IndependentBondKernel)
    sample_bond_percolation!(rng, model, spec.p)
    length(model.occupied)
end

_empty_observations(::RunSpec) = (total_energy=Float64[], total_magnetization=Int[])
_empty_observations(::PottsRunSpec) =
    (total_energy=Float64[], order_parameter=Float64[], simplex_squared=Float64[])
_empty_observations(::RandomClusterRunSpec) =
    (occupied_bonds=Int[], cluster_count=Int[], largest_cluster=Int[], cluster_size_squared_sum=Int[])

function _record_observations!(trace, model::Ising, ::RunSpec, kernel)
    push!(trace.total_energy, energy(model))
    push!(trace.total_magnetization, magnetization(model))
end
function _record_observations!(trace, model::Potts, ::PottsRunSpec, kernel)
    counts = zeros(Int, model.q)
    for colour in model.spins
        counts[colour] += 1
    end
    n, q = length(model.spins), model.q
    order = (q * (maximum(counts)/n) - 1)/(q-1)
    simplex_squared = (q * sum(count -> (count/n)^2, counts) - 1)/(q-1)
    push!(trace.total_energy, energy(model))
    # Clamp only floating round-off at the exact physical endpoints.
    push!(trace.order_parameter, clamp(order, 0.0, 1.0))
    push!(trace.simplex_squared, clamp(simplex_squared, 0.0, 1.0))
end
function _record_observations!(trace, model::BondConfiguration, ::RandomClusterRunSpec, kernel)
    # One component calculation supplies all cluster observables. It consumes
    # no random values and never recolours a random-cluster configuration.
    result = cluster_statistics(model; workspace=kernel.workspace)
    for name in keys(trace)
        push!(getproperty(trace, name), getproperty(result, name))
    end
end

_boundary(::AbstractSpinRunSpec) = "periodic square"
_boundary(spec::RandomClusterRunSpec) = spec.periodic ? "periodic square" : "open square"
_units(::RunSpec) = "kB=1; total E; signed total M"
_units(::PottsRunSpec) = "kB=1; total E; normalized maximum-colour order parameter; squared simplex magnetization"
_units(::RandomClusterRunSpec) = "bond/component/site counts; no spin energy or signed magnetization"
_observable_units(::RunSpec) = Dict("total_energy" => "total Hamiltonian", "total_magnetization" => "signed total spin")
_observable_units(::PottsRunSpec) = Dict("total_energy" => "total Hamiltonian", "order_parameter" => "(q*max_colour_count/N-1)/(q-1)", "simplex_squared" => "(q*sum_colour(count/N)^2-1)/(q-1)")
_observable_units(::RandomClusterRunSpec) = Dict("occupied_bonds" => "occupied edge count", "cluster_count" => "components including isolated vertices", "largest_cluster" => "sites in largest component", "cluster_size_squared_sum" => "sum of squared component site counts")

function _validate_model(model::Ising, spec::RunSpec)
    model.Lx == model.Ly == spec.L && model.J == spec.J || throw(ArgumentError("Ising model/configuration mismatch"))
    length(model.spins) == spec.L^2 && all(s -> s == -1 || s == 1, model.spins) || throw(ArgumentError("invalid Ising state"))
end
function _validate_model(model::Potts, spec::PottsRunSpec)
    model.Lx == model.Ly == spec.L && model.J == spec.J && model.q == spec.q || throw(ArgumentError("Potts model/configuration mismatch"))
    length(model.spins) == spec.L^2 && all(s -> 1 <= s <= spec.q, model.spins) || throw(ArgumentError("invalid Potts state"))
end
function _validate_model(model::BondConfiguration, spec::RandomClusterRunSpec)
    expected = square_bond_graph(spec.L; periodic=spec.periodic)
    model.graph.nvertices == expected.nvertices && model.graph.edges == expected.edges || throw(ArgumentError("graph/configuration mismatch"))
    length(model.occupied) == length(expected.edges) || throw(ArgumentError("invalid bond state"))
end
_validate_model(model, spec) = throw(ArgumentError("model type and specification disagree"))

function _validate_observations(trace, spec::RunSpec)
    nsites = spec.L^2
    all(isfinite, trace.total_energy) || throw(ArgumentError("nonfinite energy trace"))
    bound = 2 * nsites * abs(spec.J)
    tolerance = isfinite(bound) ? 8 * eps(Float64) * nsites * bound : 0.0
    all(e -> abs(e) <= bound + tolerance, trace.total_energy) || throw(ArgumentError("Ising energy out of bounds"))
    all(m -> -nsites <= m <= nsites && iseven(m + nsites), trace.total_magnetization) || throw(ArgumentError("invalid signed magnetization trace"))
end
function _validate_observations(trace, spec::PottsRunSpec)
    all(isfinite, trace.total_energy) || throw(ArgumentError("nonfinite energy trace"))
    lo, hi = minmax(-2 * spec.L^2 * spec.J, 0.0)
    bound = max(abs(lo), abs(hi))
    tolerance = isfinite(bound) ? 8 * eps(Float64) * spec.L^2 * bound : 0.0
    all(e -> lo-tolerance <= e <= hi+tolerance, trace.total_energy) || throw(ArgumentError("Potts energy out of bounds"))
    all(m -> isfinite(m) && 0 <= m <= 1, trace.order_parameter) || throw(ArgumentError("invalid Potts order parameter trace"))
    all(m -> isfinite(m) && 0 <= m <= 1, trace.simplex_squared) || throw(ArgumentError("invalid Potts simplex moment trace"))
end
function _validate_observations(trace, spec::RandomClusterRunSpec)
    nsites = spec.L^2
    nedges = spec.periodic ? 2 * nsites : 2 * spec.L * (spec.L-1)
    for (b, k, largest, s2) in zip(values(trace)...)
        0 <= b <= nedges || throw(ArgumentError("invalid occupied-bond count"))
        1 <= k <= nsites || throw(ArgumentError("invalid cluster count"))
        1 <= largest <= nsites-k+1 || throw(ArgumentError("invalid largest-cluster size"))
        nsites <= s2 <= nsites^2 && iseven(s2-nsites) || throw(ArgumentError("invalid cluster-size second moment"))
        nsites <= k*largest || throw(ArgumentError("largest component cannot cover all sites with the given cluster count"))
        largest^2 + (nsites-largest) <= s2 <= largest*nsites || throw(ArgumentError("inconsistent cluster-size observables"))
        k >= nsites-b || throw(ArgumentError("cluster count and occupied bonds are inconsistent"))
        b > 0 && k == nsites && throw(ArgumentError("occupied bonds must join some vertices"))
    end
end

_maximum_work(model::Union{Ising,Potts}) = length(model.spins)
_maximum_work(model::BondConfiguration) = length(model.occupied)

function _validate_kernel(kernel, model, spec)
    _algorithm(kernel) == spec.algorithm || throw(ArgumentError("kernel/configuration mismatch"))
    if spec isa AbstractSpinRunSpec
        kernel isa Union{MetropolisKernel,WolffKernel,SwendsenWangKernel} || throw(ArgumentError("invalid spin kernel"))
        kernel isa Union{WolffKernel,SwendsenWangKernel} && _check_workspace(kernel.workspace, model)
    else
        kernel isa Union{SweenyKernel,IndependentBondKernel} || throw(ArgumentError("invalid random-cluster kernel"))
        _check_workspace(kernel.workspace, model)
    end
end
