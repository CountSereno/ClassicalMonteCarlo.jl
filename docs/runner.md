# A readable and restartable Ising experiment

The execution layer separates four responsibilities:

1. `RunSpec` validates the physical parameters and sampling schedule.
2. A kernel calls the model's existing update function and owns reusable scratch arrays.
3. `RunState` owns one explicit RNG, the exact update clock and paired observations.
4. Export and analysis operate on saved data; they do not silently sample again.

This first implementation intentionally supports periodic square, zero-field
Ising. It is a complete small experiment path, not a universal simulation
framework. The existing model and update APIs remain usable without the runner.

## Run, pause and continue

```julia
using ClassicalMonteCarlo

spec = RunSpec(L=16, beta=log(1 + sqrt(2))/2,
               algorithm=:wolff, initial=:random, seed=20261030,
               burn=2000, samples=8000, stride=1)
state = start_run(spec)
advance!(state, 2000) # finish the scheduled warm-up
production_seconds = @elapsed run!(state)
export_run("results/my-run", state)
```

`algorithm` is `:metropolis`, `:wolff` or `:swendsen_wang`. One update is
respectively one N-proposal sweep, one single-cluster update or one full
bond-and-recolour sweep. `stride` counts these updates. A measurement follows
`burn + stride`, then `burn + 2stride`, and so on. `advance!` can stop during
warm-up or partway through a measurement interval. It stops at the scheduled
end; `run!` on a completed state is a no-op.

Each chain has its own RNG; interleaving another chain cannot consume its random
numbers. For multi-chain experiments, assign and record each chain's seed before
dispatch. The runner currently does not implement a parallel scheduler or
automatic seed allocation.

## Checkpoint contract

```julia
partial = start_run(spec)
advance!(partial, 2137)
save_checkpoint("checkpoints/my-run.cmc", partial)
resumed = load_checkpoint("checkpoints/my-run.cmc"; expected_spec=spec)
run!(resumed)
```

The checkpoint stores the model, random stream, kernel workspace, update count,
unfinished measurement interval and all retained traces. An ASCII/TOML envelope
records schema, Julia version, platform, RNG type, configuration and SHA-256
hashes of the source tree and package configuration. Source fingerprints are
captured when the package is loaded/precompiled and copied into each chain;
starting, saving or exporting rejects source files edited since that load.
Restart Julia after a source change rather than attaching new hashes to old
in-memory kernels. Every fingerprint input and source directory is registered
as a precompilation dependency, so the fresh process rebuilds its snapshot after
Project-only edits or adding a source file. The payload checksum is checked before
native deserialization. Incompatible checkpoints are rejected; there is no
implicit migration. `expected_spec` additionally protects against accidentally
resuming the wrong experiment.

**Only load trusted local checkpoints.** Native Julia Serialization is not a safe
format for untrusted input. A checksum detects corruption; it does not establish
authenticity. Restart is promised only in the matching software environment,
not across Julia versions. A checkpoint is written through a temporary file in
the same directory followed by rename; an existing checkpoint is atomically
replaced on supported local filesystems. This is not a guarantee against every
power-loss or filesystem failure.

Tests compare continuous execution with pause/save/load/continue at warm-up,
stride and completion boundaries for all three kernels. They compare both the
full trajectories and future RNG output, not only final means.

## Portable data contract

`export_run` accepts a new or empty directory and writes:

- `trace.csv`: `measurement,update,total_energy,total_magnetization,kernel_work`.
- `metadata.toml`: configuration, physical/update/work units, completion status,
  RNG/Julia/platform identity, source hashes and the trace checksum.

E is total Hamiltonian and M is signed total magnetization. `update` includes
warm-up. `kernel_work` sums the update return values over each stride:

| Kernel | Work field meaning |
|---|---|
| Metropolis | accepted single-spin flips |
| Wolff | flipped sites in all clusters in the interval |
| Swendsen–Wang | number of formed clusters summed over sweeps |

These are distinct diagnostics, not interchangeable CPU-work units. In
particular, the SW cluster count is not the number of flipped sites. Wall-clock
timing is measured externally so it cannot affect simulation state or restart
equality. Separate compilation, warm-up, production, measurement and I/O when
reporting performance.

Partial exports are allowed and marked `complete=false`. They are useful for
inspection, but do not satisfy a completed experiment's sampling schedule.
Export never overwrites existing evidence. Portable CSV/TOML data is separate
from the version-sensitive checkpoint.

## Deliberate limits

The runner keeps O(samples) traces in memory. It does not yet stream chunked data,
compute online accumulators, schedule multiple processes, diagnose equilibration,
choose block lengths or infer critical exponents. Statistics belong in an
explicit analysis step, using `ising_observables` on the retained paired E/M
records and checking multiple block lengths and independent chains.

To extend it, first add a kernel adapter with an explicit update and diagnostic
unit. A new representation should define its own validated state and measurement
contract rather than relabel graph observables as spin energies or magnetization.
