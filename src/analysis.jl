export block_means, jackknife, autocorrelation, integrated_autocorrelation_time,
       finite_size_fit

# These routines deliberately keep the statistical assumptions explicit. They
# operate on stationary post-equilibration time series, not raw warmup data.
function _finite_series(x::AbstractVector{<:Real}, min_length::Int=1)
    length(x) >= min_length || throw(ArgumentError("at least $min_length observations are required"))
    all(isfinite, x) || throw(ArgumentError("observations must be finite"))
    return Float64.(x)
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
    avg = mean(replicas)
    bias = (nb-1)*(avg-raw)
    stderr = sqrt((nb-1)/nb * sum((z-avg)^2 for z in replicas))
    return (estimate=raw-bias, stderr=stderr, bias=bias, raw=raw)
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
