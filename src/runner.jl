# A deliberately small execution layer: kernels update physics; the runner owns
# the clock, sampling schedule, random stream and recoverable experiment state.
using Serialization, SHA, TOML

export RunSpec, RunState, MetropolisKernel, WolffKernel, SwendsenWangKernel,
       start_run, advance!, run!, iscomplete, update_unit, work_unit,
       save_checkpoint, load_checkpoint, export_run, run_metadata

const RUN_SCHEMA_VERSION = 1
const _CHECKPOINT_MAGIC = "ClassicalMonteCarlo checkpoint 1"

"""
    RunSpec(; L, beta, J=1, algorithm=:metropolis, initial=:random,
              seed=0, burn=1000, samples=1000, stride=1)

Validated schedule for a periodic square, zero-field Ising experiment. `burn`
and `stride` count kernel updates: one Metropolis sweep, one Wolff cluster, or
one Swendsen–Wang sweep. The first measurement follows `burn + stride` updates;
initial and warm-up states are not measured. `samples` is the retained count.

The runner supports this one model family intentionally. New model families
should add their own validated specification rather than overload these units.
"""
struct RunSpec
    L::Int
    beta::Float64
    J::Float64
    algorithm::Symbol
    initial::Symbol
    seed::Int
    burn::Int
    samples::Int
    stride::Int

    function RunSpec(; L::Integer, beta::Real, J::Real=1,
                     algorithm::Symbol=:metropolis, initial::Symbol=:random,
                     seed::Integer=0, burn::Integer=1000,
                     samples::Integer=1000, stride::Integer=1)
        L >= 3 || throw(ArgumentError("L must be at least 3"))
        b, j = Float64(beta), Float64(J)
        isfinite(b) && b >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
        isfinite(j) || throw(ArgumentError("J must be finite"))
        algorithm in (:metropolis, :wolff, :swendsen_wang) || throw(ArgumentError("unknown algorithm"))
        algorithm != :metropolis && j < 0 && throw(ArgumentError("cluster kernels require J >= 0"))
        initial in (:random, :up, :down) || throw(ArgumentError("unknown initial state"))
        seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
        burn >= 0 || throw(ArgumentError("burn must be nonnegative"))
        samples > 0 || throw(ArgumentError("samples must be positive"))
        stride > 0 || throw(ArgumentError("stride must be positive"))
        l, s, w, n, d = Int(L), Int(seed), Int(burn), Int(samples), Int(stride)
        Base.checked_mul(l, l)
        Base.checked_add(w, Base.checked_mul(n, d))
        new(l, b, j, algorithm, initial, s, w, n, d)
    end
end

abstract type AbstractIsingKernel end
struct MetropolisKernel <: AbstractIsingKernel end
struct WolffKernel{W} <: AbstractIsingKernel
    workspace::W
end
struct SwendsenWangKernel{W} <: AbstractIsingKernel
    workspace::W
end
WolffKernel(model::Ising) = WolffKernel(WolffWorkspace(model))
SwendsenWangKernel(model::Ising) = SwendsenWangKernel(SwendsenWangWorkspace(model))

update_unit(::MetropolisKernel) = "metropolis_sweep"
update_unit(::WolffKernel) = "wolff_cluster"
update_unit(::SwendsenWangKernel) = "swendsen_wang_sweep"
work_unit(::MetropolisKernel) = "accepted_flips"
work_unit(::WolffKernel) = "flipped_sites"
work_unit(::SwendsenWangKernel) = "formed_clusters"
_algorithm(::MetropolisKernel) = :metropolis
_algorithm(::WolffKernel) = :wolff
_algorithm(::SwendsenWangKernel) = :swendsen_wang

_runner_step!(rng, model, beta, ::MetropolisKernel) = metropolis!(model, beta; rng=rng)
_runner_step!(rng, model, beta, k::WolffKernel) = wolff!(model, beta; rng=rng, workspace=k.workspace)
_runner_step!(rng, model, beta, k::SwendsenWangKernel) = swendsen_wang!(model, beta; rng=rng, workspace=k.workspace)

function _make_kernel(model, algorithm)
    algorithm === :metropolis && return MetropolisKernel()
    algorithm === :wolff && return WolffKernel(model)
    algorithm === :swendsen_wang && return SwendsenWangKernel(model)
    throw(ArgumentError("unknown kernel"))
end

"""
Recoverable state of one independent chain. Paired `energies` and
`magnetizations` contain total E and signed total M. `measurement_updates`
contains the total update index including warm-up. `kernel_work` is the sum of
the kernel's returned diagnostic over the measurement interval, whose meaning
is given by `work_unit(state.kernel)`; it is not a generic CPU-work measure.

Use `start_run` rather than constructing this type manually. Vectors are kept
in memory in this initial implementation, so storage is O(retained samples).
"""
mutable struct RunState{K<:AbstractIsingKernel,R<:AbstractRNG}
    spec::RunSpec
    model::Ising
    kernel::K
    rng::R
    updates::Int
    pending_work::Int
    measurement_updates::Vector{Int}
    energies::Vector{Float64}
    magnetizations::Vector{Int}
    kernel_work::Vector{Int}
    source_sha256::Dict{String,String}
end

function start_run(spec::RunSpec)
    _assert_runner_sources(_LOADED_RUNNER_SOURCES)
    rng = new_rng(spec.seed)
    model = Ising(spec.L; J=spec.J, init=spec.initial, rng=rng)
    kernel = _make_kernel(model, spec.algorithm)
    RunState(spec, model, kernel, rng, 0, 0, Int[], Float64[], Int[], Int[], copy(_LOADED_RUNNER_SOURCES))
end

_total_updates(spec::RunSpec) = spec.burn + spec.samples * spec.stride
iscomplete(state::RunState) = state.updates == _total_updates(state.spec)

"""
    advance!(state, nupdates)

Advance by at most `nupdates` kernel calls, stopping at the scheduled end.
Pausing inside warm-up or between measurements is supported. RNG state and
the unfinished interval's work counter remain in `state`. Returns `state`.
"""
function advance!(state::RunState, nupdates::Integer)
    nupdates >= 0 || throw(ArgumentError("nupdates must be nonnegative"))
    steps = Int(min(nupdates, _total_updates(state.spec) - state.updates))
    spec = state.spec
    for _ in 1:steps
        work = _runner_step!(state.rng, state.model, spec.beta, state.kernel)
        state.updates += 1
        state.updates <= spec.burn && continue
        state.pending_work += work
        if (state.updates - spec.burn) % spec.stride == 0
            push!(state.measurement_updates, state.updates)
            push!(state.energies, energy(state.model))
            push!(state.magnetizations, magnetization(state.model))
            push!(state.kernel_work, state.pending_work)
            state.pending_work = 0
        end
    end
    state
end

"""Complete the remaining warm-up and production schedule; safe to call twice."""
run!(state::RunState) = advance!(state, _total_updates(state.spec) - state.updates)

function _spec_record(spec::RunSpec)
    Dict{String,Any}("L" => spec.L, "beta" => spec.beta, "J" => spec.J,
        "algorithm" => string(spec.algorithm), "initial" => string(spec.initial),
        "seed" => spec.seed, "burn" => spec.burn, "samples" => spec.samples,
        "stride" => spec.stride)
end

function _runner_source_hashes(; register_dependencies=false)
    root = dirname(@__DIR__)
    paths = String[joinpath(root, "Project.toml")]
    for (directory, _, files) in walkdir(@__DIR__)
        # Directory dependencies detect added/deleted source files, including
        # files not included by the module. Avoid 1.11-only keyword arguments.
        register_dependencies && Base.include_dependency(directory)
        append!(paths, [joinpath(directory, f) for f in files if endswith(f, ".jl")])
    end
    if register_dependencies
        foreach(Base.include_dependency, paths)
    end
    Dict(relpath(p, root) => bytes2hex(sha256(read(p))) for p in sort(paths))
end

# Captured while the package code is loaded/precompiled, not freshly attached to
# an already-running chain at export time. Refuse disk drift rather than label
# old in-memory kernels with hashes of newly edited source files.
const _LOADED_RUNNER_SOURCES = _runner_source_hashes(; register_dependencies=true)

function _assert_runner_sources(snapshot)
    snapshot == _LOADED_RUNNER_SOURCES || throw(ArgumentError("run source snapshot differs from loaded code"))
    _runner_source_hashes() == _LOADED_RUNNER_SOURCES || throw(ArgumentError("source files changed after package load; restart Julia before running/exporting"))
    nothing
end

"""Machine-readable provenance; no timing values are inferred by the runner."""
function run_metadata(state::RunState)
    _assert_runner_sources(state.source_sha256)
    Dict{String,Any}("schema_version" => RUN_SCHEMA_VERSION,
        "julia_version" => string(VERSION), "platform" => string(Sys.MACHINE),
        "rng_type" => string(typeof(state.rng)), "model" => "Ising",
        "boundary" => "periodic square", "units" => "kB=1; total E; signed total M",
        "update_unit" => update_unit(state.kernel), "work_unit" => work_unit(state.kernel),
        "work_aggregation" => "sum over each measurement interval; not elapsed CPU work",
        "updates_completed" => state.updates, "samples_retained" => length(state.energies),
        "complete" => iscomplete(state), "config" => _spec_record(state.spec),
        "source_sha256" => copy(state.source_sha256))
end

function _validate_run_state(state::RunState)
    spec = state.spec
    0 <= state.updates <= _total_updates(spec) || throw(ArgumentError("invalid update counter"))
    n = max(0, state.updates - spec.burn) ÷ spec.stride
    all(v -> length(v) == n, (state.measurement_updates, state.energies,
                              state.magnetizations, state.kernel_work)) || throw(ArgumentError("inconsistent retained traces"))
    state.measurement_updates == [spec.burn + i*spec.stride for i in 1:n] || throw(ArgumentError("invalid measurement clock"))
    state.model.Lx == state.model.Ly == spec.L && state.model.J == spec.J || throw(ArgumentError("model/configuration mismatch"))
    length(state.model.spins) == spec.L^2 && all(s -> s == -1 || s == 1, state.model.spins) || throw(ArgumentError("invalid Ising state"))
    _algorithm(state.kernel) == spec.algorithm || throw(ArgumentError("kernel/configuration mismatch"))
    if state.kernel isa Union{WolffKernel,SwendsenWangKernel}
        _check_workspace(state.kernel.workspace, state.model)
    end
    state.pending_work >= 0 && all(>=(0), state.kernel_work) || throw(ArgumentError("invalid work counter"))
    all(isfinite, state.energies) || throw(ArgumentError("nonfinite energy trace"))
    all(m -> abs(m) <= spec.L^2 && iseven(m + spec.L^2), state.magnetizations) || throw(ArgumentError("invalid magnetization trace"))
    if state.updates <= spec.burn || (state.updates - spec.burn) % spec.stride == 0
        state.pending_work == 0 || throw(ArgumentError("unexpected unfinished interval"))
    end
    state
end

"""
    save_checkpoint(path, state)

Atomically replace a local checkpoint in the destination directory. The header
records schema, Julia version, source fingerprints, RNG type and full run
specification; the serialized state is protected by SHA-256 against corruption.
This format is for trusted, same-environment restart, not portable data exchange.
"""
function save_checkpoint(path::AbstractString, state::RunState)
    _validate_run_state(state)
    _assert_runner_sources(state.source_sha256)
    buffer = IOBuffer()
    serialize(buffer, state)
    payload = take!(buffer)
    header = run_metadata(state)
    header["payload_sha256"] = bytes2hex(sha256(payload))
    header_text = sprint(io -> TOML.print(io, header))
    parent = dirname(abspath(path))
    mkpath(parent)
    temporary, io = mktemp(parent)
    try
        println(io, _CHECKPOINT_MAGIC)
        println(io, ncodeunits(header_text))
        write(io, header_text)
        write(io, payload)
        close(io)
        # rename replaces the same-directory destination atomically on the
        # supported local filesystems. No delete-then-move gap is introduced.
        Base.Filesystem.rename(temporary, abspath(path))
    catch
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
        rethrow()
    end
    path
end

"""
    load_checkpoint(path; expected_spec=nothing)

Restore a **trusted local** checkpoint. Never load untrusted files: Julia
Serialization is not a safe untrusted-data format. Before deserializing, reject
an incompatible schema, Julia version, source tree, RNG type, expected run
specification, or payload checksum. No implicit version migration is performed.
"""
function load_checkpoint(path::AbstractString; expected_spec::Union{Nothing,RunSpec}=nothing)
    _assert_runner_sources(_LOADED_RUNNER_SOURCES)
    header, payload = open(path, "r") do io
        readline(io) == _CHECKPOINT_MAGIC || throw(ArgumentError("unrecognized checkpoint format"))
        header_length = tryparse(Int, readline(io))
        header_length !== nothing && 0 < header_length <= 1_048_576 || throw(ArgumentError("invalid checkpoint header length"))
        raw_header = read(io, header_length)
        length(raw_header) == header_length || throw(ArgumentError("truncated checkpoint header"))
        TOML.parse(String(raw_header)), read(io)
    end
    get(header, "schema_version", nothing) == RUN_SCHEMA_VERSION || throw(ArgumentError("incompatible checkpoint schema"))
    get(header, "julia_version", nothing) == string(VERSION) || throw(ArgumentError("incompatible Julia version"))
    get(header, "platform", nothing) == string(Sys.MACHINE) || throw(ArgumentError("incompatible checkpoint platform"))
    get(header, "source_sha256", nothing) == _runner_source_hashes() || throw(ArgumentError("checkpoint source/configuration code differs"))
    get(header, "rng_type", nothing) == string(typeof(new_rng(0))) || throw(ArgumentError("incompatible RNG type"))
    expected_spec !== nothing && get(header, "config", nothing) != _spec_record(expected_spec) && throw(ArgumentError("checkpoint run specification differs"))
    get(header, "payload_sha256", nothing) == bytes2hex(sha256(payload)) || throw(ArgumentError("checkpoint payload checksum mismatch"))
    state = deserialize(IOBuffer(payload))
    state isa RunState || throw(ArgumentError("checkpoint does not contain a RunState"))
    _spec_record(state.spec) == get(header, "config", nothing) || throw(ArgumentError("checkpoint header/state mismatch"))
    string(typeof(state.rng)) == header["rng_type"] || throw(ArgumentError("checkpoint RNG/header mismatch"))
    state.source_sha256 == header["source_sha256"] || throw(ArgumentError("checkpoint source/header mismatch"))
    for (key, actual) in (("updates_completed", state.updates),
                          ("samples_retained", length(state.energies)),
                          ("complete", iscomplete(state)),
                          ("update_unit", update_unit(state.kernel)),
                          ("work_unit", work_unit(state.kernel)))
        get(header, key, nothing) == actual || throw(ArgumentError("checkpoint $key/header mismatch"))
    end
    _validate_run_state(state)
end

"""
    export_run(directory, state)

Write portable `trace.csv` and `metadata.toml` into a new or empty directory.
Partial runs are allowed and explicitly marked. Existing files are never
overwritten. This data export is independent of the restart checkpoint.
"""
function export_run(directory::AbstractString, state::RunState)
    _validate_run_state(state)
    _assert_runner_sources(state.source_sha256)
    if ispath(directory)
        isdir(directory) && isempty(readdir(directory)) || throw(ArgumentError("output must be a new or empty directory"))
    else
        mkpath(directory)
    end
    trace_path = joinpath(directory, "trace.csv")
    open(trace_path, "w") do io
        println(io, "measurement,update,total_energy,total_magnetization,kernel_work")
        for i in eachindex(state.energies)
            println(io, join((i, state.measurement_updates[i], state.energies[i],
                             state.magnetizations[i], state.kernel_work[i]), ','))
        end
    end
    metadata = run_metadata(state)
    metadata["trace_sha256"] = bytes2hex(sha256(read(trace_path)))
    open(io -> TOML.print(io, metadata), joinpath(directory, "metadata.toml"), "w")
    directory
end
