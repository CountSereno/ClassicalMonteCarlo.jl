# A model-independent Metropolis–Hastings filter.
# This also works in a separate project with only Julia's Random standard library.
include(joinpath(@__DIR__, "..", "src", "acceptance.jl"))
using .Acceptance
using Random

# Physical weights, proposal probabilities, acceptance and state mutation are
# deliberately separate operations. Delta E means proposed minus current energy.
beta, delta_energy = 0.5, 2.0
logratio = metropolis_logratio(beta, delta_energy)
@assert metropolis_accept(logratio, 0.2)
@assert !metropolis_accept(logratio, 0.8)

# For a nonsymmetric proposal, include reverse / forward proposal probability.
logratio_asymmetric = metropolis_logratio(beta, delta_energy;
    log_proposal_ratio=log(0.4) - log(0.2))
@assert isapprox(logratio_asymmetric, -1 + log(2))

rng = Xoshiro(20261030)
println((symmetric_logratio=logratio,
         asymmetric_logratio=logratio_asymmetric,
         accepted=metropolis_accept(rng, logratio_asymmetric)))
