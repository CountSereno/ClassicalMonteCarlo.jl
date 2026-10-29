# A reproducible first development cycle

This study exercises the new portable Metropolis filter, reusable cluster
workspaces, and restartable Julia runner in one inspectable workflow:
**generate → validate → analyse → plot**. Generation never silently reruns
when analysing or plotting a saved dataset.

The numerical questions are deliberately narrow:

1. Do three update algorithms reproduce an independently enumerated 3×3 Ising
   distribution, including nonlinear thermodynamic observables?
2. At the known critical temperature, how do energy and magnetization
   correlations, block errors and measured throughput differ across algorithms?
3. Does checkpoint/resume reproduce uninterrupted execution in the same environment?

These are validation and finite-size benchmarking questions, not claims of new
critical exponents. See [the actual results](results/demo/RESULTS.md).

## Run from the repository root

Generation, verification and analysis use the package and Julia standard libraries only:

```bash
julia --startup-file=no --project=. experiments/foundation/run.jl \
  experiments/foundation/configs/demo.toml experiments/foundation/results/my-run
julia --startup-file=no --project=. experiments/foundation/verify.jl \
  experiments/foundation/results/my-run --require-current-source
julia --startup-file=no --project=. experiments/foundation/analyze.jl \
  experiments/foundation/results/my-run
julia --startup-file=no --project=. experiments/foundation/window_sensitivity.jl \
  experiments/foundation/results/my-run
```

Generation requires a new or empty output directory. Analysis only writes its
own derived tables/report. Every individual trace has a checksum and the runner's
portable metadata. A restart checkpoint is retained for the first chain; it is
only suitable for the recorded code/Julia environment. Never load an untrusted
Julia serialization file.

Plotting has its own locked environment, keeping GUI/graphics dependencies out of
the simulation core:

```bash
julia --project=experiments/plotting -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=experiments/plotting \
  experiments/foundation/plot.jl experiments/foundation/results/my-run
```

For this workspace, the temporary working depot is
`JULIA_DEPOT_PATH=/tmp/cmc-foundation-depot`.
On another machine, use the normal writable Julia depot instead. Existing saved
data can be verified, analysed and plotted without running the generator.
Omit `--require-current-source` when inspecting historical data after changing
code: checksums and physical checks still run, with any source drift reported.

## Model and sampling conventions

- Ferromagnetic Ising on a periodic square lattice, J=kB=1, zero field.
- Exact benchmark: L=3, T=3; all 512 configurations enumerated by a separate
  coordinate-based implementation that does not call the package's energy or
  neighbour routines.
- Critical benchmark: L=8,16,24, with the known Tc=2/log(1+sqrt(2)) inserted.
- Metropolis, Wolff, and Swendsen–Wang, four independent seeds per condition;
  starts random/up/down/random. Each chain has 8,192 warmup updates and 32,768
  production records, one recorded algorithm update apart. Total: 48 chains
  and 1,572,864 measurements.
- One Metropolis update is N random site proposals; one Wolff update flips one
  cluster; one SW update constructs/recolors all FK clusters. These clocks must
  not be compared as if each step had equal computational cost.
- E and M in the traces are total energy and signed total magnetization. SW
  `kernel_work` counts formed clusters, not flipped spins. Wolff spin-sweep
  conversion uses average changed-spin count/N and is only a work convention.

## Analysis and limitations

All analysis is written in Julia. `ising_observables` deletes paired E/M time
blocks together, preserving moment covariance. Block sizes 128, 512 and 2,048
are retained; each group uses the larger of the independent-chain SEM and the
quadrature-combined largest within-chain block error. This is a descriptive
error rule, not a calibrated simultaneous confidence interval. Four chains are
insufficient for precision coverage claims.

The per-site definitions are `cᵥ = β² Var(E)/N`, `χ = β Var(M)/N`,
`χ_abs = β(⟨M²⟩ − ⟨|M|⟩²)/N`, and `U₄ = 1 − ⟨M⁴⟩/(3⟨M²⟩²)`.
The absolute-magnetization fluctuation proxy `χ_abs` is distinct from the
zero-field susceptibility `χ`. Each chain uses its measured signed mean M in
the latter; failure to explore both magnetization signs can bias that estimate.

The positive-window autocorrelation estimator stops at the first nonpositive
ACF. Full energy, absolute-M and signed-M ACFs are retained at L=24. The positive
window can miss negative tails and slow modes; the code does not apply its ESS
formula to signed M, Binder or other nonlinear observables. Block and window
warnings remain in the report. `ESS = n/(2τ)` uses `τ = 1/2 + Σρ(t)` over the
retained positive window. The block warnings compare block length to `10τ` for
energy and absolute M only; this is a heuristic, not a bound for every observable.
Warmup sensitivity and larger sizes remain future
work; neither a small tau nor agreement among a few chains proves equilibration.

The executed run produced two chains with ACF-window flags at lag 512. A
[separate offline follow-up](results/demo/WINDOW_SENSITIVITY.md) recalculates
the L=24 diagnostics at maximum lags 512, 1,024 and 2,048 from the same traces.
It preserves the primary analysis and records its own source/input hashes.

Timing includes production updates plus E/M measurement, excludes JIT compilation,
warmup, storage and analysis. All three algorithm paths are compiled on separate
throwaway states first; an explicit garbage collection precedes each timing.
Each chain is timed once, sequentially on one machine;
throughput is a preliminary implementation measurement, not a universal hardware
or algorithm ranking. Figure lines only guide the eye; no Tc or exponent is fitted.

## File responsibilities

| File | Responsibility |
|---|---|
| `run.jl` | Read configuration, invoke the runner, persist evidence, audit restart |
| `exact.jl` | Independent finite-system enumeration |
| `tables.jl` | Small explicit reader for this unquoted numeric CSV schema |
| `verify.jl` | Read-only checksums, clocks, physical limits and exact limits |
| `analyze.jl` | Joint statistics, ACF, chain aggregation and result report |
| `window_sensitivity.jl` | Saved-trace follow-up at three ACF maximum lags |
| `plot.jl` | CairoMakie PNG/PDF/SVG from saved analysis tables |
| `configs/demo.toml` | The executed, fixed design |

Source hashes identify the actual implementation used. Reusing a seed after an
algorithm/refactor change does not guarantee the same trajectory, and Julia does
not promise identical RNG streams across versions. Preserve the raw records as
well as the seed and environment.
