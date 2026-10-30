# A shared runner with model-specific observations

The execution layer separates five responsibilities:

1. Each specification validates one model's physical parameters and a shared schedule.
2. Kernel adapters call the existing update functions and own scratch workspaces.
3. One typed `RunState` owns the model, RNG, clock and paired observation vectors.
4. Export/checkpoint code handles evidence and restart without knowing statistical estimators.
5. Offline analysis interprets the named physical quantities and estimates uncertainty.

`src/runner/specifications.jl` defines configurations;
`src/runner/adapters.jl` defines initialization, stepping, measurements and model
validation; `src/runner.jl` implements the common scheduling/persistence path.
The core uses only Julia standard libraries. Kernels remain usable without this
runner, and the standalone acceptance rule does not depend on it.

## Specify a model, then use the same execution API

```julia
using ClassicalMonteCarlo

ising = RunSpec(L=16, beta=log(1 + sqrt(2))/2,
                algorithm=:wolff, seed=20261030,
                burn=2000, samples=8000, stride=1)

potts = PottsRunSpec(L=8, q=3, beta=log(1 + sqrt(3)),
                     algorithm=:swendsen_wang, initial=:ordered,
                     seed=20261031, burn=2000, samples=8000, stride=1)

rc = RandomClusterRunSpec(L=8, p=sqrt(3)/(1 + sqrt(3)), q=3,
                          periodic=true, algorithm=:sweeny, initial=:empty,
                          seed=20261032, burn=2000, samples=8000, stride=1)

state = start_run(potts)
advance!(state, 2137)
save_checkpoint("checkpoints/potts.cmc", state)
resumed = load_checkpoint("checkpoints/potts.cmc"; expected_spec=potts)
run!(resumed)
export_run("results/potts", resumed)
E = observations(resumed).total_energy
m = observations(resumed).order_parameter
```

These example lengths are not a convergence prescription. Sweeny now uses a
reusable interleaved BFS workspace; metadata records `interleaved_bfs`. Full
cluster measurements have separate union-find scratch. The graph topology is
fixed while occupations change. See the [connectivity/scaling study](../experiments/scaling/README.md).

`RunSpec` keeps the original Ising constructor. `PottsRunSpec` is periodic square,
integer q≥2, zero field; `:random` and `:ordered` initial states are supported.
The Hamiltonian is `-J Σ δ(s_i,s_j)`; q=2 uses a different coupling convention
from `-J Σ s_i s_j` Ising. Metropolis accepts either sign of J; both spin cluster
algorithms require J≥0.

`RandomClusterRunSpec` supports open or periodic square graphs. Sweeny requires
0<p<1, finite q>0 (including nonintegers) and at least one edge. The independent
Bernoulli kernel is allowed only at q=1, including p=0 and p=1. For example:

```julia
percolation = RandomClusterRunSpec(L=8, p=.5, q=1, algorithm=:independent,
                                   burn=0, samples=8000, seed=17)
```

Graph initial states are `:empty`, `:full` or `:random`. Random initialization
uses independent Bernoulli(p) bonds and is an equilibrium draw only for q=1.
Open boundaries allow L≥1; periodic boundaries require L≥3, with no self-loops
or multiple edges. The graph's ordered edge list is part of restart validation.
The runner currently constructs square graphs; arbitrary `BondGraph` objects
remain available through the lower-level kernel API.

## Explicit clocks and paired data

The first measurement follows `burn + stride`, then `burn + 2stride`, and so on.
`advance!` can stop inside warm-up or partway through a measurement interval and
stops at the scheduled end. `run!` on a completed state is a no-op. Each chain
owns its RNG: interleaving other chains does not consume this chain's stream.

| Kernel | One update | `kernel_work`, summed over the measurement interval |
|---|---|---|
| Metropolis | N random site proposals with replacement | accepted flips |
| Wolff | one single-cluster update | flipped sites |
| Swendsen–Wang | one full bond-and-recolour sweep | formed clusters |
| Sweeny | m lazy edge proposals with replacement | accepted bond flips |
| Independent bonds | a Bernoulli redraw of all m bonds | resampled bonds (m) |

These counts are distinct diagnostics, not interchangeable CPU-work units. In
particular, an SW cluster count is not a flipped-spin count, and a full independent
redraw is not a Sweeny sweep. Kernel parameters and return contracts are unchanged
by scheduling. Measurements consume no RNG values.

`observations(state)` returns a typed NamedTuple of live vectors. All columns
have the same measurement indices and must remain paired in nonlinear analysis.
This accessor does not allocate copies; copy columns before destructive analysis.

| Model | Columns after `measurement,update` and before `kernel_work` |
|---|---|
| Ising | `total_energy,total_magnetization` |
| Potts | `total_energy,order_parameter,simplex_squared` |
| Random cluster | `occupied_bonds,cluster_count,largest_cluster,cluster_size_squared_sum` |

For Potts, `order_parameter = (q max_a N_a/N - 1)/(q-1)` is a maximum-colour
quantity, not signed total magnetization. `simplex_squared =
(q Σ_a (N_a/N)^2 - 1)/(q-1)` is the squared simplex magnetization. Both are
computed from the same colour counts; floating round-off at exact endpoints is
clamped to [0,1]. For the zero-field ferromagnetic model with integer q and
`p = 1-exp(-beta*J)`, its equilibrium mean can be compared with `⟨Σ_c |c|²⟩/N²`
in the corresponding Fortuin–Kasteleyn ensemble. Neither a
maximum-colour variance nor an Ising Binder formula should be silently relabelled
as a general Potts susceptibility or Binder definition.

Random-cluster components include isolated vertices. `cluster_size_squared_sum`
is the unnormalized sum Σ_c |c|², not a spin magnetization or an energy. Physical
mappings to Potts require the matching p, q, graph and boundary conditions.

Old Ising code may continue using `state.energies` and `state.magnetizations`;
these are compatibility aliases of the same Ising vectors. Other models expose
the explicit observation names. Ising CSV column names and order are unchanged.

## Checkpoint contract

The checkpoint stores the physical model, RNG, kernel workspace, update count,
unfinished measurement interval and all retained observations. Its TOML envelope
records schema, Julia version, platform, model/specification, units, observation
columns, RNG type and SHA-256 hashes of source files and package configuration.

Source fingerprints are captured when the package is loaded/precompiled and
copied into each chain. Start/save/export reject source changes since that load.
The compilation cache tracks the Project file and source directory contents as
well as included files, so a normal restart refreshes the snapshot after edits.
Restore checks schema, Julia/platform, source, expected specification and payload
checksum before native deserialization; it then validates model/graph, workspace,
observation schema, clocks, work counters and all recorded metadata.

**Only load trusted local checkpoints.** Julia Serialization is not a safe format
for untrusted input. A checksum detects corruption; it does not authenticate a
file. Checkpoints are written through a temporary file and same-directory rename.
Atomic replacement on supported local filesystems is not a guarantee against all
power-loss or filesystem failures.

The shared model-specific state uses checkpoint **schema 3**. Schema 1 and 2 native
checkpoints require the original code environment and are rejected here before
deserialization. There is no implicit migration. Old portable CSV/TOML datasets
remain readable independently of native restart, and their recorded source
hashes remain historical evidence.

Tests compare uninterrupted runs with pause/save/load/continue at warm-up,
stride and completion boundaries for Ising, Potts and random cluster, including
noninteger q and independent percolation. They compare complete observations,
final configurations and future RNG draws, not only final averages.

## Portable export and current limits

`export_run` writes `trace.csv` and `metadata.toml` into a new or empty directory.
Metadata includes `trace_columns`, `observable_units`, the model-discriminated
configuration, clock/work meanings, completion status, provenance and a trace
checksum. Partial exports are marked `complete=false`; export never overwrites
existing evidence. Native checkpoints are separate from these portable tables.

The runner retains O(samples) vectors in memory. It does not yet stream chunks,
schedule processes, choose equilibration lengths, choose block sizes or infer
critical exponents. Wall-clock timing is external to state, so it cannot affect
restart equality. Benchmarks must separate compilation, warm-up, production,
measurements and I/O.

A new model should add a validated specification plus initialization, kernel,
measurement and validation adapters. The common scheduling and persistence code
should not acquire model-specific branches or unrelated optional fields.
