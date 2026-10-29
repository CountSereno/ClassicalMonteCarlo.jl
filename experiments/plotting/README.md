# Julia figure environment

This environment is separate from the standard-library-only simulation package.
The checked-in Manifest records the resolved graphics stack used for the
foundation figures: Julia 1.12.6, CairoMakie 0.15.16 and Makie 0.24.16.

From the repository root:

```bash
julia --project=experiments/plotting -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=experiments/plotting \
  experiments/foundation/plot.jl experiments/foundation/results/demo
```

The script reads saved analysis tables and writes PNG, SVG and PDF figures; it
never launches a simulation. Figure metadata records the input-table hashes,
plot-script hash and Manifest hash. Keep the resolved Manifest when reproducing
the saved figures; updating dependencies is a separate, deliberate operation.

Use Julia 1.12.6 to reproduce this graphics environment. The core package's
Julia 1.10 compatibility target does not imply that this resolved plotting
environment was tested on Julia 1.10.

In a sandbox with a read-only default font cache, set `XDG_CACHE_HOME` to a
writable temporary directory before starting Julia. The local run used
`XDG_CACHE_HOME=/tmp/cmc-foundation-cache`; this changes cache placement, not the
simulation or figure inputs.
