export block_means, jackknife, autocorrelation, integrated_autocorrelation_time,
       finite_size_fit, ising_observables

# These routines deliberately keep the statistical assumptions explicit. They
# operate on stationary post-equilibration time series, not raw warmup data.
function _finite_series(x::AbstractVector{<:Real}, min_length::Int=1)
    length(x) >= min_length || throw(ArgumentError("at least $min_length observations are required"))
    all(isfinite, x) || throw(ArgumentError("observations must be finite"))
    return Float64.(x)
end

function _jackknife_result(raw::Real, replicas::AbstractVector{<:Real})
    nb = length(replicas)
    avg = mean(replicas)
    bias = (nb-1)*(avg-raw)
    stderr = sqrt((nb-1)/nb * sum((z-avg)^2 for z in replicas))
    return (estimate=raw-bias, stderr=stderr, bias=bias, raw=Float64(raw))
end

"""
    block_means(x, blocksize)

Mean of each consecutive, non-overlapping block. `length(x)` must be divisible
by `blocksize`, so no observations are silently discarded. For error estimates,
blocks must be long compared with the correlation time and approximately
independent; checking this is the caller's responsibility.
"""
function block_means(x::AbstractVector{<:Real}, blocksize::Integer)
    v = _finite_series(x)
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    length(v) % blocksize == 0 || throw(ArgumentError("blocksize must divide the number of observations"))
    return [mean(@view v[i:i+blocksize-1]) for i in 1:blocksize:length(v)]
end

"""
    jackknife(f, x; blocksize=1)

Delete-one-block jackknife for a scalar statistic `f`, which is called on the
full vector and on vectors with each contiguous block removed. Returns
`(estimate, stderr, bias, raw)`; `estimate = raw - bias`. The length must be
divisible by `blocksize` and there must be at least two blocks. `f` must return
a finite real scalar. Blocks should be approximately independent and `f`
should be sufficiently smooth for the jackknife approximation to be reliable.
"""
function jackknife(f, x::AbstractVector{<:Real}; blocksize::Integer=1)
    v = _finite_series(x, 2)
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    length(v) % blocksize == 0 || throw(ArgumentError("blocksize must divide the number of observations"))
    nb = length(v) ÷ blocksize
    nb >= 2 || throw(ArgumentError("at least two blocks are required"))
    scalar(y) = y isa Real && isfinite(y) ? Float64(y) : throw(ArgumentError("statistic must return a finite real scalar"))
    raw = scalar(f(v))
    replicas = Vector{Float64}(undef, nb)
    for k in 1:nb
        start = (k-1)*blocksize + 1
        stop = k*blocksize
        replicas[k] = scalar(f(vcat(@view(v[1:start-1]), @view(v[stop+1:end]))))
    end
    return _jackknife_result(raw, replicas)
end

"""
    ising_observables(E, M; beta, nsites, blocksize=1)

Estimate zero-field Ising thermodynamic observables from paired stationary
time series of **total** energy `E` and **signed total** magnetization `M`.
The six result fields are `energy_density`, `abs_magnetization_density`,
`heat_capacity`, `susceptibility`, `susceptibility_abs`, and `binder`. Each is
`(estimate, stderr, bias, raw)` with a delete-contiguous-block jackknife;
`estimate = raw - bias`. Metadata `n_samples`, `n_blocks`, and `blocksize`
record the data actually used. No observations are silently discarded.

With `N=nsites` and k_B=1, the raw definitions are
`<E>/N`, `<|M|>/N`, `beta^2*(<E^2>-<E>^2)/N`,
`beta*(<M^2>-<M>^2)/N`, `beta*(<M^2>-<|M|>^2)/N`, and
`U4=1-<M^4>/(3<M^2>^2)`, respectively. `susceptibility_abs` is an
absolute-magnetization fluctuation proxy, **not** the finite-system linear
response susceptibility. The signed definition requires sampling both
magnetization sectors at zero field.

All moments from each time block are deleted together, preserving their
covariance in the nonlinear errors. Energy is centered before moments are
formed to avoid subtracting nearly equal large squared energies. Cost is
linear in sample count plus block count. Blocks must be approximately
independent: compare multiple block lengths and independent equilibrated
chains; this function does not diagnose equilibration or choose block sizes.
At least two equal blocks are required. `binder` must be defined both for
the full sample and every delete-block sample (`<M^2> > 0`), or an
`ArgumentError` is raised. Bias-corrected finite-sample estimates can lie
outside the physical range; inspect `raw` and the uncertainty as well.
"""
function ising_observables(E::AbstractVector{<:Real}, M::AbstractVector{<:Real};
                           beta::Real, nsites::Integer, blocksize::Integer=1)
    energies = _finite_series(E, 2)
    magnetic = _finite_series(M, 2)
    length(energies) == length(magnetic) || throw(ArgumentError("E and M must have equal lengths"))
    nsites > 0 || throw(ArgumentError("nsites must be positive"))
    n = Float64(nsites)
    isfinite(n) || throw(ArgumentError("nsites must have a finite Float64 representation"))
    b = _check_beta(beta)
    isfinite(b) || throw(ArgumentError("beta must have a finite Float64 representation"))
    all(x -> abs(x) <= n, magnetic) || throw(ArgumentError("total magnetization must lie in [-nsites, nsites]"))
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    length(energies) % blocksize == 0 || throw(ArgumentError("blocksize must divide the number of observations"))
    nb = length(energies) ÷ blocksize
    nb >= 2 || throw(ArgumentError("at least two blocks are required"))

    offset = mean(energies)
    moments = zeros(nb, 6)
    for k in 1:nb
        for i in ((k-1)*blocksize+1):(k*blocksize)
            e = (energies[i] - offset)/n
            m = magnetic[i]/n
            moments[k, 1] += e
            moments[k, 2] += e^2
            moments[k, 3] += m
            moments[k, 4] += abs(m)
            moments[k, 5] += m^2
            moments[k, 6] += m^4
        end
    end
    moments ./= blocksize
    all(isfinite, moments) || throw(ArgumentError("moments overflowed; rescale the energy units"))
    total = vec(sum(moments; dims=1))
    function evaluate(v)
        v[5] > 0 || throw(ArgumentError("Binder cumulant is undefined when a full or delete-block sample has zero mean M^2"))
        values = (offset/n + v[1], v[4],
                  b^2*n*(v[2]-v[1]^2),
                  b*n*(v[5]-v[3]^2),
                  b*n*(v[5]-v[4]^2),
                  1-v[6]/(3v[5]^2))
        all(isfinite, values) || throw(ArgumentError("observable calculation overflowed; rescale the energy units"))
        values
    end
    raw = evaluate(total/nb)
    replicas = Matrix{Float64}(undef, nb, 6)
    for k in 1:nb
        replicas[k, :] .= evaluate((total - @view(moments[k, :]))/(nb-1))
    end
    results = ntuple(j -> _jackknife_result(raw[j], @view(replicas[:, j])), 6)
    return (energy_density=results[1], abs_magnetization_density=results[2],
            heat_capacity=results[3], susceptibility=results[4],
            susceptibility_abs=results[5], binder=results[6],
            n_samples=length(energies), n_blocks=nb, blocksize=Int(blocksize))
end

"""
    autocorrelation(x; maxlag=length(x)-1)

Return normalized autocorrelation estimates at lags 0 through `maxlag`, using
one common denominator (the sum of squared deviations). Constant series have
undefined normalized autocorrelation and are rejected. Results assume a
stationary time series; noise grows at large lags.
"""
function autocorrelation(x::AbstractVector{<:Real}; maxlag::Integer=length(x)-1)
    v = _finite_series(x, 2)
    0 <= maxlag < length(v) || throw(ArgumentError("maxlag must lie between 0 and length(x)-1"))
    centered = v .- mean(v)
    denominator = sum(abs2, centered)
    denominator > 0 && isfinite(denominator) || throw(ArgumentError("series must have finite nonzero variance"))
    return [sum(centered[i]*centered[i+lag] for i in 1:length(v)-lag)/denominator for lag in 0:maxlag]
end

"""
    integrated_autocorrelation_time(x; maxlag=min(length(x)-1, 1000))

Estimate `tau_int = 1/2 + sum(rho(t))` by summing positive-lag estimates up to
(but excluding) the first nonpositive value or `maxlag`, whichever comes first.
Returns `(tau, window, stderr_mean)`, where `stderr_mean` uses
`var(x) * 2tau / length(x)`. This simple positive-window estimate is noisy for
short chains and is not a convergence test. Strongly anticorrelated chains need
a different windowing method; the truncation cannot capture negative tails.
"""
function integrated_autocorrelation_time(x::AbstractVector{<:Real}; maxlag::Integer=min(length(x)-1, 1000))
    v = _finite_series(x, 2)
    maxlag >= 1 || throw(ArgumentError("maxlag must be at least 1"))
    rho = autocorrelation(v; maxlag=maxlag)
    tau = 0.5
    window = 0
    for lag in 1:maxlag
        rho[lag+1] > 0 || break
        tau += rho[lag+1]
        window = lag
    end
    return (tau=tau, window=window, stderr_mean=sqrt(var(v)*2*tau/length(v)))
end

"""
    finite_size_fit(L, y; correction_exponent, sigma=nothing)

Fit the fixed-exponent ansatz `y(L) = intercept + slope * L^(-correction_exponent)`
by least squares. This is a basic extrapolation tool, not a general nonlinear
finite-size-scaling fit: choose the exponent and size range based on physics,
and inspect residuals/corrections. With `sigma`, independent known one-sigma
uncertainties are used as weights and parameter errors use the unscaled
inverse normal matrix. Without `sigma`, residual variance is estimated with
`n-2` degrees of freedom. At least three finite observations and distinct
transformed sizes are required. Returns coefficients, standard errors, residuals,
and `chi2` (weighted only when `sigma` is supplied).
"""
function finite_size_fit(L::AbstractVector{<:Real}, y::AbstractVector{<:Real}; correction_exponent::Real, sigma=nothing)
    lengths = _finite_series(L, 3)
    vals = _finite_series(y, 3)
    length(lengths) == length(vals) || throw(ArgumentError("L and y must have equal lengths"))
    all(>(0), lengths) || throw(ArgumentError("sizes must be positive"))
    isfinite(correction_exponent) && correction_exponent > 0 || throw(ArgumentError("correction_exponent must be finite and positive"))
    x = lengths .^ (-correction_exponent)
    all(isfinite, x) || throw(ArgumentError("transformed sizes must be finite"))
    weights = if sigma === nothing
        ones(length(vals))
    else
        sigma isa AbstractVector{<:Real} || throw(ArgumentError("sigma must be a vector of real uncertainties"))
        s = _finite_series(sigma, 3)
        length(s) == length(vals) || throw(ArgumentError("sigma and y must have equal lengths"))
        all(>(0), s) || throw(ArgumentError("uncertainties must be positive"))
        w = 1.0 ./ (s .^ 2)
        all(isfinite, w) && all(>(0), w) || throw(ArgumentError("uncertainties yield invalid weights"))
        w
    end
    sw = sum(weights)
    xbar = sum(weights .* x)/sw
    ybar = sum(weights .* vals)/sw
    sxx = sum(weights .* (x .- xbar).^2)
    isfinite(sxx) && sxx > 0 || throw(ArgumentError("transformed sizes must not all be equal"))
    slope = sum(weights .* (x .- xbar) .* (vals .- ybar))/sxx
    intercept = ybar - slope*xbar
    residuals = vals .- (intercept .+ slope .* x)
    chi2 = sum(weights .* residuals.^2)
    scale = sigma === nothing ? chi2/(length(vals)-2) : 1.0
    return (intercept=intercept, slope=slope,
            intercept_stderr=sqrt(scale*(1/sw+xbar^2/sxx)),
            slope_stderr=sqrt(scale/sxx), residuals=residuals,
            chi2=chi2, dof=length(vals)-2, weighted=sigma !== nothing)
end
