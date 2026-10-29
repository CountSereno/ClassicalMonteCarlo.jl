"""Independent exhaustive Ising reference, separate from production kernels.

Only use on small periodic lattices (3 <= L <= 4). The right/up bonds are
written explicitly; this does not call energy, magnetization or the sampler's
neighbour helper. J=kB=1 and the field is zero.
"""
module ExactIsing
using Statistics
export enumerate_ising

function enumerate_ising(L::Integer, beta::Real)
    3 <= L <= 4 || throw(ArgumentError("exact enumeration supports L=3 or 4"))
    isfinite(beta) && beta >= 0 || throw(ArgumentError("finite beta >= 0 required"))
    N = L^2
    energies = Float64[]
    magnetic = Float64[]
    spectrum = Dict{Int,Int}()
    for code in 0:(1 << N)-1
        spin = [iszero(code & (1 << (site-1))) ? -1 : 1 for site in 1:N]
        e = 0
        for y in 1:L, x in 1:L
            i = x + L*(y-1)
            right = mod(x,L)+1 + L*(y-1)
            up = x + L*mod(y,L)
            e -= spin[i]*(spin[right]+spin[up])
        end
        spectrum[e] = get(spectrum,e,0)+1
        push!(energies,e)
        push!(magnetic,sum(spin))
    end
    weights = exp.(-beta .* (energies .- minimum(energies)))
    weights ./= sum(weights)
    moment(v) = sum(weights .* v)
    E = moment(energies)
    M = moment(magnetic)
    M2 = moment(magnetic.^2)
    result = Dict(
        "L"=>L,"beta"=>Float64(beta),"states"=>length(energies),
        "energy_density"=>E/N,
        "abs_magnetization_density"=>moment(abs.(magnetic))/N,
        "heat_capacity"=>beta^2*(moment(energies.^2)-E^2)/N,
        "susceptibility"=>beta*(M2-M^2)/N,
        "susceptibility_abs"=>beta*(M2-moment(abs.(magnetic))^2)/N,
        "binder"=>1-moment(magnetic.^4)/(3M2^2),
        "density_of_states"=>Dict(string(k)=>v for (k,v) in spectrum))
    @assert sum(values(spectrum)) == 2^N
    @assert spectrum[-2N] == 2
    @assert abs(M) < 1e-12
    result
end
end
