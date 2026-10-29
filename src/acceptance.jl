"""
Model-independent Metropolis–Hastings acceptance rules.

This file can also be used outside ClassicalMonteCarlo:
`include("acceptance.jl"); using .Acceptance`.
It depends only on Julia's Random standard library. The pure scalar method is
the reference contract; RNG management, proposals and state mutation belong to
the caller.
"""
module Acceptance

using Random: AbstractRNG, rand

export metropolis_accept, metropolis_logratio

@inline function check_logratio(logratio::Real)
    isnan(logratio) && throw(ArgumentError("log acceptance ratio must not be NaN"))
    logratio
end

"""
    metropolis_accept(logratio, u) -> Bool

Pure Metropolis–Hastings decision for a supplied uniform variate `0 <= u < 1`.
`logratio` means

    log π(proposed) - log π(current)
      + log Q(current | proposed) - log Q(proposed | current).

The acceptance probability is `min(1, exp(logratio))`. Comparison in log space
avoids overflow/underflow from exponentiating extreme ratios. The comparison is
strict: equality with the acceptance threshold rejects. A finite negative ratio
accepts `u=0`, while `logratio=-Inf` always rejects (zero acceptance probability).
`+Inf` always accepts. NaN ratios and variates outside `[0,1)` are errors.

No random numbers are generated and no input is mutated. Keeping this operation
separate makes acceptance tests and ports to another model or language simple.
"""
@inline function metropolis_accept(logratio::Real, u::Real)
    check_logratio(logratio)
    isfinite(u) && 0 <= u < 1 || throw(ArgumentError("uniform variate must lie in [0, 1)"))
    logratio >= 0 && return true
    logratio == -Inf && return false
    log(u) < logratio
end

"""
    metropolis_accept(rng::AbstractRNG, logratio) -> Bool

RNG adapter to the pure acceptance rule. Deterministic outcomes (`logratio >= 0`
or `-Inf`) consume no random value; a finite negative ratio consumes one uniform
value. This consumption convention is part of the interface. Pass a separate
RNG to each independent chain.
"""
@inline function metropolis_accept(rng::AbstractRNG, logratio::Real)
    check_logratio(logratio)
    logratio >= 0 && return true
    logratio == -Inf && return false
    metropolis_accept(logratio, rand(rng))
end

"""
    metropolis_logratio(beta, delta_energy; log_proposal_ratio=0)

Build the log acceptance ratio for a Boltzmann target:
`-beta * delta_energy + log_proposal_ratio`, with `delta_energy = E_new - E_old`.
`log_proposal_ratio` is **reverse minus forward** log proposal probability; it
is zero only for a symmetric proposal. `beta` must be finite and nonnegative,
and `delta_energy` finite. The proposal term may be ±Inf but not NaN.

Integer inputs are converted to floating arithmetic before multiplication,
avoiding integer overflow. Overflow of finite floating inputs to ±Inf represents
a deterministic decision; an indeterminate sum of opposite infinities is
rejected. This helper encodes the target ratio, not proposal generation or
Boltzmann equilibration.
"""
@inline function metropolis_logratio(beta::Real, delta_energy::Real;
                                    log_proposal_ratio::Real=0)
    isfinite(beta) && beta >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    isfinite(delta_energy) || throw(ArgumentError("energy difference must be finite"))
    check_logratio(log_proposal_ratio)
    ratio = -float(beta) * float(delta_energy) + float(log_proposal_ratio)
    check_logratio(ratio)
end

end # module Acceptance
