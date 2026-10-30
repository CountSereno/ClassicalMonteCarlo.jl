# Core design: portable kernels and a shared model runner

The package connects portable weight/acceptance rules, model updates, an explicit
simulation clock, recoverable runs, offline analysis and Julia figures. Ising,
Potts and random-cluster models share execution and persistence while keeping
their physical parameters and observations distinct. Existing model/update APIs
and Ising observation aliases remain usable. The larger architecture in
[the research plan](research-plan.md) remains a development direction.

## Separate numerical decisions from physical state

`src/acceptance.jl` is a self-contained module that depends only on `Random`.
It can be copied or included in another Julia project without the lattice,
models, runner or plotting dependencies.

```julia
include("acceptance.jl")
using .Acceptance

# Energy and proposal generation belong to the model/caller.
logratio = metropolis_logratio(beta, proposed_energy - current_energy;
    log_proposal_ratio=log_q_reverse - log_q_forward)

# Pure form: useful in tests and ports; no random state is touched.
accepted = metropolis_accept(logratio, uniform_variate)

# Adapter form: explicitly receives the chain's RNG.
accepted = metropolis_accept(rng, logratio)
```

The sign convention is `log(target_new/target_old) + log(proposal_reverse/proposal_forward)`.
The pure function accepts a uniform value in `[0,1)` and compares in log space;
it rejects NaN inputs, handles infinite ratios explicitly, and never mutates a
configuration. The RNG adapter consumes no value for deterministic outcomes and
one for a finite negative log ratio. Tests cover this contract and an asymmetric
proposal's exact detailed balance. [Runnable standalone example](../examples/acceptance.jl).

Ising/Potts Metropolis, Sweeny, Worm and Wang–Landau reuse this decision. They
still compute their own weight ratios and implement their own state changes.
Wang–Landau deliberately draws one value each proposal, including deterministic
acceptance, matching that algorithm's existing RNG-consumption convention.
Using the same acceptance filter does not make their ensembles or clocks the same.

For Sweeny, `sweeny_logratio(occupied, endpoints_disconnected, p, q)` is a pure
weight calculation. Connectivity is determined with the selected edge removed,
so it distinguishes deleting a bridge from deleting an edge on a cycle. The
reference update first queries connectivity, then computes the ratio, applies
the lazy proposal and acceptance decision, and finally toggles the occupation.
This keeps graph traversal, probability and state mutation independently
reviewable. The small disjoint-set primitives in `src/disjoint_sets.jl` contain
no model or RNG dependency and are shared by SW and bond connectivity.

## Reuse scratch memory without hiding model updates

`WolffWorkspace(model)` contains a visited mask and a bounded cluster queue.
`SwendsenWangWorkspace(model)` contains union-find and recolouring arrays.
Neither owns the physical spins or an RNG. Updates reset the necessary fields
on each call and reject incompatible dimensions/types. Each concurrent chain
must own a workspace.

```julia
workspace = WolffWorkspace(model)
for _ in 1:1000
    wolff!(model, beta; rng=rng, workspace=workspace)
end
```

Omitting the keyword retains the simple existing API by constructing temporary
storage. Tests compare the complete trajectories and RNG streams of both paths,
including repeated calls and zero coupling. `benchmarks/workspaces.jl` measures
allocations after compilation; it does not infer a wall-clock speedup from fewer
allocated bytes.

`BondWorkspace(graph_or_state)` similarly owns reusable connectivity arrays.
`cluster_statistics` obtains all four geometric observations from one component
calculation and consumes no RNG values. Measurements may reuse the kernel's
scratch arrays because each later query resets them. This is the full-rebuild
reference. The default `SweenyWorkspace` instead uses static contiguous adjacency,
two BFS queues and signed visit stamps. Searches alternate endpoints and stop on
contact or exhaustion of one component. Only graph topology is cached, never
connectivity across bond changes. Graph topology must remain fixed.

`endpoints_connected`, the pure `sweeny_logratio`, and the acceptance filter have
separate responsibilities. Lazy skips and q=1 avoid queries, with exactly the
same RNG stream and trajectory as the original algorithm. Cluster measurement
uses a separate `BondWorkspace` inside the Sweeny workspace. Both algorithms
retain O(m*(n+m)) worst-case sweep cost; IBFS is not a dynamic forest.
The [paired timing study](../experiments/scaling/README.md) measures actual
speedups independently of allocation and checks identical final states/RNG.

The separate [bond allocation probe](../benchmarks/bonds.jl) reports allocation
for updates and measurements independently. Its measured [L=8 result](../benchmarks/results/bond_allocations.csv)
from commit 44ddbc7 on Julia 1.12.6 uses 1,000 reference sweeps: fresh per-sweep scratch allocated 1,152,000 B,
while reused sweeps and reused geometric measurements each allocated 0 measured
heap bytes. Production traces and checkpoint storage still allocate memory.

The recorded [allocation probe](../benchmarks/results/workspace_allocations.csv)
on Julia 1.12.6 used L=16 and 1,000 updates after warm-up. Reused workspaces
allocated 0 measured heap bytes for both Ising and q=3 Potts, with both cluster
algorithms. Fresh storage allocated 2.24 MB for Wolff, 4.544 MB for Ising SW and
6.336 MB for Potts SW. The probe checks matching spins and future RNG output;
these figures concern scratch allocation only, not total application memory.

```bash
julia --startup-file=no --project=. benchmarks/workspaces.jl
```

## Keep execution and inference separate

The [runner](runner.md) uses one configuration type per model: `RunSpec` for
Ising, `PottsRunSpec` and `RandomClusterRunSpec`. The shared `RunState` has concrete
types for specification, model, kernel, RNG and observation table. Its typed
NamedTuple contains only that model's vectors; it has no irrelevant optional
beta, energy or magnetization placeholders. `observations(state)` exposes these
paired vectors; the old Ising `energies`/`magnetizations` properties are narrow
compatibility aliases.

Model adapters construct configurations, invoke kernels, append measurements
and validate physical bounds. The shared loop only manages the update count,
warm-up, measurement interval and unfinished work count. `advance!` pauses at
any update boundary; `run!` completes the schedule. Kernels do not write files
or calculate error bars. `src/model_observables.jl` operates on saved numeric
vectors without any model construction or random sampling.

One update is a Metropolis sweep, Wolff cluster, SW sweep, lazy Sweeny sweep or
whole-graph independent redraw. Its work counter has a separate named meaning:
accepted flips, flipped sites, formed clusters, accepted bond flips or resampled
bonds. These are never treated as a common CPU-work unit.

Portable CSV/TOML exports support independent analysis. A checksum-guarded local
checkpoint preserves the in-memory state for trusted, same-environment restart.
Restoration rejects incompatible versions/configuration/code; it is not a
general cross-version interchange format. Schema 3 adds IBFS scratch/backend metadata and rejects schema 1/2
native checkpoints; their original source environment is needed for restart.
Portable historical CSV/TOML datasets remain readable. The runner keeps traces in
memory, O(samples), and does not yet implement streaming, MPI, arbitrary models,
or a generic observable plugin system.

## Dependencies and review boundaries

The core remains standard-library-only. CairoMakie and its resolved environment
live under `experiments/plotting`, separate from the simulation package. The
[foundation experiments](../experiments/foundation/README.md) produce raw traces,
analysis tables and figures in separate steps. Prior experiment artifacts remain
historical evidence for their recorded source hashes; a changed source tree is
not expected to reproduce an old checkpoint without the original environment.
The [multi-model experiments](../experiments/multimodel/README.md) add independent
enumeration and finite-volume FK energy and simplex-moment identities. They do
not equate maximum-colour order with the largest FK cluster or relabel a
geometric variance as a spin heat capacity.

When extending a kernel, first add its mathematical ratio, then a readable
reference transition, exact/invariant tests, and only then an optimized path.
When adding a model to the runner, define its state and estimator units explicitly
rather than adding unrelated optional fields to the existing Ising specification.

## Finite-size inference

`src/scaling.jl` accepts numeric arrays independently of models and runners.
Positive-data power fits use centered QR; adjacent effective exponents retain
their shared-data covariance. The [scaling experiment](../experiments/scaling/README.md)
keeps fit cutoffs, fixed correction ansatzes, block/bootstrap uncertainty and
source/input hashes in analysis artifacts. A fit standard error does not include
corrections to scaling or failure to equilibrate.
