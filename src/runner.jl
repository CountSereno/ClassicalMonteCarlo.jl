# Shared execution and persistence. Model definitions, kernel updates and
# statistical inference remain independent of this scheduling layer.
using Serialization, SHA, TOML

export AbstractRunSpec, RunSpec, PottsRunSpec, RandomClusterRunSpec, RunState,
       MetropolisKernel, WolffKernel, SwendsenWangKernel, SweenyKernel,
       IndependentBondKernel, observations, start_run, advance!, run!, iscomplete,
       update_unit, work_unit, save_checkpoint, load_checkpoint, export_run,
       run_metadata

# Schema 3 records the Sweeny connectivity backend and serializes IBFS scratch.
# Older native checkpoints require their original source environment.
const RUN_SCHEMA_VERSION = 3
const _CHECKPOINT_MAGIC = "ClassicalMonteCarlo checkpoint 3"

include("runner/specifications.jl")
include("runner/adapters.jl")

"""
One recoverable independent chain with concrete specification, model, kernel,
RNG and observation-table types. `observations(state)` returns paired vectors
whose names and physical meanings depend on the model. `measurement_updates`
is the total update count including warm-up; `kernel_work` sums the kernel's
named diagnostic over the corresponding interval, not elapsed CPU work.

Construct with `start_run`. Storage is O(retained samples). For old Ising code,
`state.energies` and `state.magnetizations` remain aliases of the same vectors;
other models use the explicit `observations` interface.
"""
mutable struct RunState{S<:AbstractRunSpec,M,K<:AbstractUpdateKernel,R<:AbstractRNG,T<:NamedTuple}
    spec::S
    model::M
    kernel::K
    rng::R
    updates::Int
    pending_work::Int
    measurement_updates::Vector{Int}
    trace::T
    kernel_work::Vector{Int}
    source_sha256::Dict{String,String}
end

"""
    observations(state)

The live, paired observation vectors as a typed NamedTuple. This accessor does
not copy, sample, or consume RNG state. Names are the physical CSV column names:
Ising `(total_energy, total_magnetization)`, Potts `(total_energy,
order_parameter, simplex_squared)`, and random cluster `(occupied_bonds,
cluster_count, largest_cluster, cluster_size_squared_sum)`. Treat retained
observations as read-only evidence during a run; copy them for destructive work.
"""
observations(state::RunState) = getfield(state, :trace)

# Narrow compatibility adapter: graph counts and Potts order parameters never
# masquerade as signed Ising magnetizations. Scheduling uses observations().
function Base.getproperty(state::RunState{RunSpec}, name::Symbol)
    name === :energies && return getfield(state, :trace).total_energy
    name === :magnetizations && return getfield(state, :trace).total_magnetization
    getfield(state, name)
end
Base.propertynames(state::RunState{RunSpec}, private::Bool=false) =
    (fieldnames(typeof(state))..., :energies, :magnetizations)

function start_run(spec::AbstractRunSpec)
    _assert_runner_sources(_LOADED_RUNNER_SOURCES)
    rng = new_rng(spec.seed)
    model = _initialize_model(rng, spec)
    kernel = _make_kernel(model, spec)
    RunState(spec, model, kernel, rng, 0, 0, Int[], _empty_observations(spec),
             Int[], copy(_LOADED_RUNNER_SOURCES))
end

_total_updates(spec::AbstractRunSpec) = spec.burn + spec.samples * spec.stride
iscomplete(state::RunState) = state.updates == _total_updates(state.spec)

"""
    advance!(state, nupdates)

Advance by at most `nupdates` kernel calls and stop at the scheduled end. One
call is a named sweep/cluster/redraw, defined by `update_unit(state.kernel)`.
Pauses inside warm-up and inside a measurement interval preserve both the RNG
and unfinished work count. The first measurement is at `burn + stride`.
"""
function advance!(state::RunState, nupdates::Integer)
    nupdates >= 0 || throw(ArgumentError("nupdates must be nonnegative"))
    steps = Int(min(nupdates, _total_updates(state.spec) - state.updates))
    spec = state.spec
    for _ in 1:steps
        work = _runner_step!(state.rng, state.model, spec, state.kernel)
        state.updates += 1
        state.updates <= spec.burn && continue
        state.pending_work += work
        if (state.updates - spec.burn) % spec.stride == 0
            _record_observations!(observations(state), state.model, spec, state.kernel)
            push!(state.measurement_updates, state.updates)
            push!(state.kernel_work, state.pending_work)
            state.pending_work = 0
        end
    end
    state
end

"""Complete the remaining warm-up and production schedule; safe to call twice."""
run!(state::RunState) = advance!(state, _total_updates(state.spec) - state.updates)

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

_trace_columns(state::RunState) =
    ["measurement", "update", string.(collect(keys(observations(state))))..., "kernel_work"]

"""Machine-readable physical units, schedule, data schema and source provenance."""
function run_metadata(state::RunState)
    _assert_runner_sources(state.source_sha256)
    Dict{String,Any}("schema_version" => RUN_SCHEMA_VERSION,
        "julia_version" => string(VERSION), "platform" => string(Sys.MACHINE),
        "rng_type" => string(typeof(state.rng)), "model" => _model_name(state.spec),
        "boundary" => _boundary(state.spec), "units" => _units(state.spec),
        "trace_columns" => _trace_columns(state),
        "observable_units" => _observable_units(state.spec),
        "update_unit" => update_unit(state.kernel), "work_unit" => work_unit(state.kernel),
        "connectivity_backend" => _connectivity_backend(state.kernel),
        "work_aggregation" => "sum over each measurement interval; not elapsed CPU work",
        "updates_completed" => state.updates,
        "samples_retained" => length(state.measurement_updates),
        "complete" => iscomplete(state), "config" => _spec_record(state.spec),
        "source_sha256" => copy(state.source_sha256))
end

function _validate_run_state(state::RunState)
    spec, trace = state.spec, observations(state)
    0 <= state.updates <= _total_updates(spec) || throw(ArgumentError("invalid update counter"))
    n = max(0, state.updates - spec.burn) ÷ spec.stride
    typeof(trace) == typeof(_empty_observations(spec)) || throw(ArgumentError("observation table/model mismatch"))
    all(v -> length(v) == n, values(trace)) && length(state.measurement_updates) == n &&
        length(state.kernel_work) == n || throw(ArgumentError("inconsistent retained traces"))
    state.measurement_updates == [spec.burn + i*spec.stride for i in 1:n] || throw(ArgumentError("invalid measurement clock"))
    _validate_model(state.model, spec)
    _validate_kernel(state.kernel, state.model, spec)
    _validate_observations(trace, spec)
    maximum_work = _maximum_work(state.model)
    all(w -> 0 <= w <= maximum_work*spec.stride, state.kernel_work) || throw(ArgumentError("invalid retained work counter"))
    interval_updates = state.updates <= spec.burn ? 0 : (state.updates-spec.burn) % spec.stride
    0 <= state.pending_work <= maximum_work*interval_updates || throw(ArgumentError("invalid unfinished interval work"))
    if state.kernel isa IndependentBondKernel
        all(==(maximum_work*spec.stride), state.kernel_work) || throw(ArgumentError("independent redraw work must equal the number of resampled edges"))
        state.pending_work == maximum_work*interval_updates || throw(ArgumentError("incomplete independent redraw work counter"))
    end
    state
end

"""
    save_checkpoint(path, state)

Atomically replace a local checkpoint via a same-directory temporary file and
rename. The TOML header records schema, Julia/platform, model/configuration,
physical column units, source fingerprints and RNG type. A payload SHA-256
check detects corruption before native deserialization.

This format is only for trusted, matching-environment restart, not portable
exchange or untrusted input. Schema 1 checkpoints require the old source tree.
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

Restore a trusted local checkpoint in the matching Julia/code/platform
context. The header, expected model specification and payload checksum are
checked before deserialization; the restored model, observation schema,
workspace and counters are then validated. No implicit migration is performed.
Never load untrusted Julia Serialization data.
"""
function load_checkpoint(path::AbstractString; expected_spec::Union{Nothing,AbstractRunSpec}=nothing)
    _assert_runner_sources(_LOADED_RUNNER_SOURCES)
    header, payload = open(path, "r") do io
        readline(io) == _CHECKPOINT_MAGIC || throw(ArgumentError("unrecognized or incompatible checkpoint format"))
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
    if expected_spec !== nothing
        get(header, "model", nothing) == _model_name(expected_spec) &&
            get(header, "config", nothing) == _spec_record(expected_spec) || throw(ArgumentError("checkpoint run specification differs"))
    end
    get(header, "payload_sha256", nothing) == bytes2hex(sha256(payload)) || throw(ArgumentError("checkpoint payload checksum mismatch"))
    state = deserialize(IOBuffer(payload))
    state isa RunState || throw(ArgumentError("checkpoint does not contain a RunState"))
    _validate_run_state(state)
    actual = run_metadata(state)
    for key in keys(actual)
        get(header, key, nothing) == actual[key] || throw(ArgumentError("checkpoint $key/header mismatch"))
    end
    state
end

"""
    export_run(directory, state)

Write portable `trace.csv` and `metadata.toml` to a new or empty directory.
Columns are `measurement,update`, the named model observations, then
`kernel_work`. The Ising column names/order remain compatible with schema 1;
Potts and random-cluster exports explicitly use their own physical quantities.
Partial runs are allowed and marked. Existing files are never overwritten.
"""
function export_run(directory::AbstractString, state::RunState)
    _validate_run_state(state)
    metadata = run_metadata(state)
    if ispath(directory)
        isdir(directory) && isempty(readdir(directory)) || throw(ArgumentError("output must be a new or empty directory"))
    else
        mkpath(directory)
    end
    trace_path = joinpath(directory, "trace.csv")
    columns = values(observations(state))
    open(trace_path, "w") do io
        println(io, join(_trace_columns(state), ','))
        for i in eachindex(state.measurement_updates)
            values_at_i = map(column -> column[i], columns)
            println(io, join((i, state.measurement_updates[i], values_at_i..., state.kernel_work[i]), ','))
        end
    end
    metadata["trace_sha256"] = bytes2hex(sha256(read(trace_path)))
    open(io -> TOML.print(io, metadata), joinpath(directory, "metadata.toml"), "w")
    directory
end
