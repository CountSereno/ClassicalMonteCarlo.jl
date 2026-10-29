# Core design and first implementation increment

The first implementation increment follows one complete path: a portable
acceptance rule, existing model updates, an explicit simulation clock, recoverable
Ising runs, offline analysis and Julia figures. It preserves the current model
constructors and update APIs. The larger architecture in [the research plan](research-plan.md)
is a development direction; it is not all implemented in this increment.

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

The [runner](runner.md) intentionally implements one well-defined Ising workflow.
`RunSpec` validates the physical settings and schedule; `RunState` owns the model,
RNG, kernel workspace, update count, unfinished measurement interval and paired
time series. `advance!` pauses at any update boundary; `run!` completes the schedule.
Kernels do not write files or calculate error bars.

One runner update is a Metropolis sweep, a Wolff cluster, or an SW sweep.
The meaning of the returned work counter is recorded separately: accepted flips,
flipped sites, or formed clusters. It is never treated as a common CPU-work unit.

Portable CSV/TOML exports support independent analysis. A checksum-guarded local
checkpoint preserves the in-memory state for trusted, same-environment restart.
Restoration rejects incompatible versions/configuration/code; it is not a
general cross-version interchange format. This initial runner keeps traces in
memory, O(samples), and does not yet implement streaming, MPI, arbitrary models,
or a generic observable plugin system.

## Dependencies and review boundaries

The core remains standard-library-only. CairoMakie and its resolved environment
live under `experiments/plotting`, separate from the simulation package. The
[foundation experiments](../experiments/foundation/README.md) produce raw traces,
analysis tables and figures in separate steps. Prior experiment artifacts remain
historical evidence for their recorded source hashes; a changed source tree is
not expected to reproduce an old checkpoint without the original environment.

When extending a kernel, first add its mathematical ratio, then a readable
reference transition, exact/invariant tests, and only then an optimized path.
When adding a model to the runner, define its state and estimator units explicitly
rather than adding unrelated optional fields to the existing Ising specification.
