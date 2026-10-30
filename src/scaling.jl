export power_law_fit, effective_exponents

"""
    power_law_fit(L, y; sigma=nothing)

Fit positive data to `y=A*L^exponent` by linear least squares in log space.
With independent known standard errors `sigma`, the delta-method log errors
are `sigma/y`; these must be small for the approximation to be useful. QR is
used on a centered design matrix. The returned covariance is for
`(log_amplitude, exponent)`; it is not rescaled by reduced chi-square when
errors are supplied. Unweighted errors use residual variance with n-2 degrees
of freedom. `residuals` and `chi2` are in log space.

This is a leading-power diagnostic, not an automatic asymptotic exponent
determination. Inspect size-cutoff stability, residuals and correction terms.
Correlated sizes require a covariance-aware fit outside this interface.
"""
function power_law_fit(L::AbstractVector{<:Real}, y::AbstractVector{<:Real}; sigma=nothing)
    lengths, values = _positive_scaling_data(L,y,3)
    x,z = log.(lengths),log.(values)
    errors = _scaling_errors(sigma,values)
    x0=mean(x)
    design=hcat(ones(length(x)),x .- x0)
    weighted_design=design ./ errors
    factor=qr(weighted_design)
    coefficients=factor \ (z ./ errors)
    r=Matrix(factor.R)
    all(isfinite,coefficients) && minimum(abs.(diag(r))) > eps(Float64)*maximum(abs.(r)) ||
        throw(ArgumentError("sizes give a numerically singular design"))
    residuals=z-design*coefficients
    chi2=sum(abs2,residuals ./ errors)
    dof=length(x)-2
    scale=sigma === nothing ? chi2/dof : 1.0
    inverse_r=r \ Matrix{Float64}(I,2,2)
    transform=[1.0 -x0; 0.0 1.0]
    covariance=scale .* (transform*inverse_r*transpose(inverse_r)*transpose(transform))
    log_amplitude=coefficients[1]-x0*coefficients[2]
    (amplitude=exp(log_amplitude),log_amplitude,exponent=coefficients[2],
     exponent_stderr=sqrt(covariance[2,2]),log_amplitude_stderr=sqrt(covariance[1,1]),
     covariance,residuals,chi2,dof,weighted=sigma !== nothing)
end

function _positive_scaling_data(L,y,minimum_length)
    lengths,values=_finite_series(L,minimum_length),_finite_series(y,minimum_length)
    length(lengths)==length(values) || throw(ArgumentError("L and y must have equal lengths"))
    all(>(0),lengths) && all(>(0),values) || throw(ArgumentError("sizes and observations must be positive"))
    all(diff(lengths).>0) || throw(ArgumentError("sizes must be strictly increasing"))
    lengths,values
end
function _scaling_errors(sigma,values)
    sigma === nothing && return ones(length(values))
    sigma isa AbstractVector{<:Real} || throw(ArgumentError("sigma must be a real vector"))
    errors=_finite_series(sigma,length(values))
    length(errors)==length(values) && all(>(0),errors) || throw(ArgumentError("sigma must be positive and match y"))
    errors ./= values
    all(x -> isfinite(x) && x>0,errors) || throw(ArgumentError("invalid relative errors"))
    errors
end

"""
    effective_exponents(L, y; sigma)

Adjacent-size slopes `log(y[i+1]/y[i])/log(L[i+1]/L[i])` and their full
delta-method covariance, for independent measurements at each size. Adjacent
slopes share one datum and are anticorrelated; do not fit them as independent.
Requires positive means/errors and strictly increasing sizes.
"""
function effective_exponents(L::AbstractVector{<:Real},y::AbstractVector{<:Real};sigma)
    lengths,values=_positive_scaling_data(L,y,2)
    errors=_scaling_errors(sigma,values)
    sigma === nothing && throw(ArgumentError("effective exponents require uncertainties"))
    n=length(values)-1
    differences=diff(log.(lengths))
    jacobian=zeros(n,n+1)
    for i in 1:n
        jacobian[i,i]=-1/differences[i]
        jacobian[i,i+1]=1/differences[i]
    end
    covariance=jacobian*Diagonal(errors.^2)*transpose(jacobian)
    (sizes=sqrt.(lengths[1:end-1].*lengths[2:end]),
     exponents=diff(log.(values))./differences,
     stderr=sqrt.(diag(covariance)),covariance)
end
