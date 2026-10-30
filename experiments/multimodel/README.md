# Potts and random-cluster experiments with one runner

This executed study exercises model-specific specifications and observations
through one shared scheduling, restart and export interface. Sampling, checking,
analysis and plotting are separate commands. [Read the numerical report](results/demo/RESULTS.md).

The three questions are:

1. Do Potts Metropolis, Wolff and Swendsen–Wang reproduce independently
   enumerated finite-volume observables?
2. Does Sweeny sample integer and noninteger random-cluster weights correctly,
   and does its q=1 limit agree with independent Bernoulli redraws?
3. On the same lattice, do independently sampled spin and FK-bond ensembles
   satisfy their energy and simplex-moment identities?

## Reproduce from the repository root

The first three steps require only Julia standard libraries and this package:

```bash
julia --startup-file=no --project=. experiments/multimodel/run.jl \
  experiments/multimodel/configs/demo.toml experiments/multimodel/results/my-run
julia --startup-file=no --project=. experiments/multimodel/verify.jl \
  experiments/multimodel/results/my-run --require-current-source
julia --startup-file=no --project=. experiments/multimodel/analyze.jl \
  experiments/multimodel/results/my-run
```

Generation accepts only a new or empty output directory. Analysis reads the
saved traces without resampling. Every trace has a checksum, full run
configuration, named update/work units and hashes of the loaded source. The
study manifest records unique seeds, initials, timing boundaries and restart
audits. Checkpoints are trusted, same-environment restart files; CSV/TOML are the
portable analysis records. Never load a checkpoint from an untrusted source.

Julia/CairoMakie plotting has a separate resolved environment:

```bash
julia --project=experiments/plotting -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=experiments/plotting \
  experiments/multimodel/plot.jl experiments/multimodel/results/my-run
```

These commands also inspect or plot the supplied `results/demo` without
rerunning Monte Carlo. Omit `--require-current-source` only when inspecting a
historical dataset after editing the code; the verifier then reports source
drift while retaining data-integrity and physical checks. The source hashes in
the data are never silently rewritten.

## Executed design and model conventions

All lattices are periodic squares; J=kB=1 for spins. Each condition has four
independent seeds. Potts starts alternate random/ordered; bond starts are
empty/full/random/random. The configuration records 4,096 warmup updates,
16,384 production records and stride one, for **44 chains and 720,896 records**.

- Potts exact benchmark: L=3, q=3, beta=0.7, three update algorithms; enumerate
  all 19,683 spin assignments independently using coordinate-based bonds.
- Random-cluster exact benchmark: the same L=3 graph, p=1-exp(-0.7), q=1,1.5,3;
  enumerate all 262,144 occupied-edge subsets using independent BFS connectivity.
  Sweeny is compared with these references; q=1 also uses independent full redraws.
- Cross-representation benchmark: q=3, L=4,8, beta=log(1+sqrt(q)),
  p=sqrt(q)/(1+sqrt(q)); compare Potts SW with random-cluster Sweeny.

The square-lattice critical values are supplied from the known theoretical
result, not estimated from this small study. See
[Beffara and Duminil-Copin](https://arxiv.org/abs/1006.5073).

One Metropolis update is N randomly selected spin proposals; one Wolff update
is one cluster; one SW update constructs and recolours all FK clusters. One
Sweeny update in this runner is m randomly selected lazy edge proposals, where
m is the number of graph edges. An independent update redraws every edge.
The returned diagnostic counts have separate units in each run's metadata;
they are not a universal amount of CPU work.

## Observables and exact FK relations

Potts records total E, the normalized maximum-colour order
`m_max=(q*max_a(n_a)/N-1)/(q-1)`, and
`simplex_squared=(q*sum_a(n_a²)/N²-1)/(q-1)`.
The largest-colour population and the largest FK component are distinct finite
system observables; this study does not identify their means.

Potts analysis reports E/N, m_max, heat capacity `beta² Var(E)/N`, the
explicitly named order fluctuation proxy `beta*N*Var(m_max)`, moment ratio
`<m_max⁴>/<m_max²>²`, and `<simplex_squared>`. The proxy is not a signed linear
response susceptibility, and the moment ratio is not labelled an Ising Binder
cumulant. The runner records FK-compatible simplex information separately.

Random-cluster records b occupied edges, k connected components, largest
component size and `S2=sum_c(size_c²)`, including isolated vertices. Analysis
reports `<b>/m`, `<k>/N`, `<largest>/N`, `<S2>/N` and `Var(b)/m`.
`<S2>/N` is the mean size of the cluster containing a uniformly selected vertex;
the largest component is included, so it is not a connected fluctuation with
that component removed.

For integer q, assigning an independent uniform Potts colour to each FK
component gives two finite-graph identities:

```
E[b | spins] = p * number of like-colour edges
<E>/N = -J <b>/(pN)
E[simplex_squared | bonds] = S2/N²
<simplex_squared> = <S2>/N²
```

The same graph, boundary, q and p=1-exp(-beta*J) are essential. The q=1.5
ensembles have a valid bond measure and exact reference, but no integer-colour
spin model is assigned to them. For the random-cluster representation and
connectivity costs see [Elçi and Weigel](https://arxiv.org/html/1307.6647v4).
The code used for this recorded dataset rebuilds connectivity with reusable scratch memory. It does not
implement the dynamic-connectivity algorithm studied in that work.

## Error bars, correlations and costs

All moments for one observation are removed together in block jackknife
estimates. Block sizes 128,512,2048 are retained. The grouped diagnostic SE is
the larger of the between-chain SEM and the quadrature-combined largest
within-chain block SE. With four chains this rule does not supply calibrated
confidence coverage, particularly for nonlinear moment ratios and variances.

ACFs are retained for both critical representations at L=8 and for both q=1
bond algorithms. Primary diagnostics use Potts energy or RC bond count;
secondary diagnostics use simplex_squared or S2. The positive-window estimator
stops at the first nonpositive ACF or lag 512. A flag is retained when the
window reaches its limit or a chosen block is shorter than 10 estimated tau.
No flag, agreement between a few chains, or flat error curve proves equilibrium.
Slow tails, initial-state dependence and larger systems need further work.

Reported timing includes production transitions and all recorded measurements,
but excludes compilation, warmup, disk output and analysis. Selected kernel
paths are compiled on throwaway states first; each sequential chain is timed
once after explicit garbage collection. Times are short and depend on machine
state and implementation. Energy and bond count have different conditional
noise and transition processes, so their ESS/second is observable-specific.
Lines guide the eye; no dynamic exponent or transition order is fitted.

The Sweeny traces are measured once per full sweep. They cannot resolve
sub-sweep relaxation or the critical speeding-up discussed by
[Deng et al.](https://arxiv.org/abs/cond-mat/0701113). A positive-window estimate
near its 0.5-record floor does not exclude faster modes; studying them would
require explicitly recording edge-attempt time and a suitable estimator.

## Files and responsibilities

| File | Responsibility |
|---|---|
| `run.jl` | Generate explicit cases, call runner, persist traces and audit restart |
| `exact.jl` | Independent spin enumeration and BFS bond-subset enumeration |
| `verify.jl` | Read-only schema, checksum, provenance and physical constraints |
| `tables.jl` | Explicit reader/writer for this numeric CSV schema |
| `analyze.jl` | Paired moments, block errors, ACFs and FK comparisons |
| `plot.jl` | PNG/SVG/PDF figures from analysis tables only |
| `configs/demo.toml` | Fixed design, seeds and statistical settings |

The core package contains the reusable estimators and runner. This directory
contains the study-specific design, exact references and presentation, keeping
physics state updates independent of tables, plotting and filesystem layout.

## Historical source environment

The included dataset was generated at [commit 44ddbc7](https://github.com/CountSereno/ClassicalMonteCarlo.jl/tree/44ddbc750e7fa32878ee9a5b42ddabf921326674).
Its schema-2 native checkpoints require that environment. Current code uses
IBFS and schema 3; historical CSV/TOML remains readable. Omit
`--require-current-source` when verifying historical data against a newer tree.
The [next-stage study](../scaling/README.md) records the new backend and larger sizes.
