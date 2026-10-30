# Sweeny connectivity and finite-size scaling

This stage connects an independently checked optimization to a reproducible
critical-scaling study. The [executed report](results/critical/RESULTS.md)
contains the numerical results; [configuration](configs/critical.toml) records
the full design. The core uses only Julia standard libraries.

## Connectivity design

`SweenyWorkspace` builds a static contiguous adjacency table once. Each query
excludes the candidate edge and alternates BFS expansions from its endpoints.
It stops on contact or when either queue empties. Signed generation stamps
avoid a full visit-array clear; wraparound explicitly clears the marks.
There is no persistent connectivity cache to invalidate after a deletion.
The graph topology must stay fixed; occupations can change freely. A workspace
is owned by one chain, and its separate `measurements` workspace computes full
cluster statistics only when requested.

`endpoints_connected` owns graph traversal, `sweeny_logratio` owns the pure
weight ratio, and `metropolis_accept` owns the acceptance decision. Passing a
`BondWorkspace` selects the full-rebuild reference; passing `SweenyWorkspace`
selects IBFS. Lazy proposals and q=1 skip unnecessary searches. RNG consumption
matches the previous algorithm exactly, including deterministic acceptance.
The runner now defaults to IBFS and records its backend in metadata.

The method follows the interleaved-search approach described in
[Elçi and Weigel (2013), §II.2.1](https://arxiv.org/html/1307.6647v4#S2.SS2.SSS1).
It is independently implemented here. It is **not** the fully dynamic forest
algorithm discussed later in that paper: worst-case query work remains O(N+m).

Validation includes exhaustive bond subsets on a complete five-vertex graph,
disconnected graphs, bridges/cycles, visit-stamp rollover, corrupted-checkpoint
guards, and stepwise comparison with the original implementation. Both accepted
updates and future random draws must match, not just distributional averages.

## Paired performance benchmark

Run from the repository root:

```bash
julia --startup-file=no --project=. benchmarks/connectivity.jl /tmp/cmc-benchmark
```

The [recorded benchmark](../../benchmarks/results/connectivity/summary.csv)
uses five alternating-order repeats per size/weight, after compilation and
warmup. Both backends start each repeat from the same configuration and RNG;
final occupations and future RNG output are checked. Both skip lazy queries.
At L=64 the measured median update-only speedups are about 31× (q=2) and 47×
(q=3); both allocate zero heap bytes per reused-workspace sweep. These are
single-machine short timings, not asymptotic complexity estimates. Full study
timings additionally include cluster measurements. Raw repeat times, hardware,
Julia version and source hashes are saved alongside the summary.

## Experiment and observable definitions

- Periodic square lattices, L=8,12,16,24,32,48,64, q=2 and 3: four Sweeny
  chains per point, from empty/full/random/random occupations.
- Independent q=3 Potts SW chains at L=16,32,64: four random/ordered starts.
- Each chain discards 8192 sweeps, then retains 16384 measurements at stride 1.
  There are 68 chains and 1,114,112 records. Seeds differ at every parameter
  point and chain; there is no common-random-number coupling across sizes.
- `C1 = largest_cluster`, including a possible wrapping component, is a site
  count. `S2_over_N = sum_C |C|² / N` includes all components and isolated sites.
  It is not the finite-cluster-only percolation susceptibility.
- All parameters are set at the known critical point
  `p_c = sqrt(q)/(1+sqrt(q))`, `beta_c*J = log(1+sqrt(q))`, J=1.
  This identity comes from the [square-lattice critical-point result](https://arxiv.org/abs/1006.5073);
  this experiment does not estimate p_c or nu.

At criticality, the leading powers used as reference are
`C1 ∝ L^d_f` and `S2_over_N ∝ L^(2d_f-2)`:

| q | d_f | 2d_f-2 |
|---:|---:|---:|
| 2 | 15/8 | 7/4 |
| 3 | 28/15 | 26/15 |

Definitions, the leading dimensions and multiple correction mechanisms are
discussed by [Xu, Salas and Deng (2025)](https://doi.org/10.3390/e27040418).
Their numerical study uses a loop model with different correction amplitudes;
we do not transfer those amplitudes or a fitted correction exponent to this
square-lattice bond model.

## Statistical analysis

Each retained trajectory is analyzed separately. Consecutive blocks of 128,
512 and 2048 sweeps give block errors; 512 is the preset primary choice.
Group means use four equally weighted chains. The displayed diagnostic error
is the larger of the independently combined primary-block errors and the SE
of the four chain means. This is a conservative diagnostic rule, not a proven
confidence-coverage guarantee. The analysis additionally repeats the L>=16 fits
with each of the three block lengths, retaining all individual large-block flags.

`power_law_fit` uses centered QR in log space, with independent relative errors
sigma/y. Leading-power fits use every cutoff leaving at least three sizes.
400 within-chain block-bootstrap replicas provide conditional interval checks
with fixed original fit weights. They cannot capture finite-size model bias.
`effective_exponents` reports adjacent log ratios with their full covariance:
neighbouring ratios share data and are anticorrelated. Separate fixed-reference
power fits test `y/L^x = a + b*L^-omega` at omega=0.8,1.2,4/3,2.0. These are nuisance
ansatz sensitivity checks; the data do not determine a correction exponent.

Per-chain diagnostics include ACF windows, block/tau, block-error stability,
and early/late production differences. They do not prove equilibrium. The burn
trajectory is not saved; the analysis also repeats L>=16 fits after discarding
one-quarter or one-half of production. These overlapping fits are correlated
and provide a sensitivity check, not a substitute for longer preset burns.
The q=3 SW comparison tests `<E>/N = -<b>/(pN)` and
`N*<simplex_squared> = <S2>/N`. These are independent equilibrium-mean checks,
not identifications of individual spin and bond configurations.

## Reproduce without overwriting evidence

```bash
# Sampling and strict read-only checks. Choose a fresh output directory.
julia --startup-file=no --project=. experiments/scaling/run.jl \
  experiments/scaling/configs/critical.toml experiments/scaling/results/my-run
julia --startup-file=no --project=. experiments/scaling/verify.jl \
  experiments/scaling/results/my-run --require-current-source

# Analyze saved observations; never runs a sampler.
julia --startup-file=no --project=. experiments/scaling/analyze.jl \
  experiments/scaling/results/my-run

# A separate pinned environment keeps CairoMakie out of the simulation package.
julia --startup-file=no --project=experiments/plotting -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=experiments/plotting experiments/scaling/plot.jl \
  experiments/scaling/results/my-run
julia --startup-file=no --project=. experiments/scaling/verify_analysis.jl \
  experiments/scaling/results/my-run --require-current-source
```

For the included results use `results/critical` in these commands and omit the
sampling step. Native checkpoints use schema 3 and require their recorded Julia
and source environment. Historical CSV/TOML data from earlier stages remain
unchanged and readable; source drift must not be presented as an exact replay.

| File | Responsibility |
|---|---|
| `common.jl` | Case design and file/source provenance helpers |
| `run.jl` | Sampling, production timing, export and exact restart audits |
| `verify.jl` | Shared trace checks plus design/backend/checkpoint provenance |
| `analyze.jl` | Blocks, fits, bootstrap, cross-representation checks and report |
| `plot.jl` | Four PNG/SVG/PDF figure sets from saved tables only |
| `verify_analysis.jl` | Artifact hashes and independent closed-form fit checks |

Next steps are larger sizes/longer chains with multiple preset burn lengths,
explicit sub-sweep sampling when studying dynamics, and a separately verified
dynamic-connectivity backend if traversal cost warrants the complexity.
