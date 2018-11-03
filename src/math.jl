export binom2, 
       xlogx, xlogy,
       xexpx, xexpy,
       log1pexp, log1mexp,
       midpoint


"fast binomial(n,2)"
@inline function binom2(n::Int)
	@boundscheck @assert n ≥ 0
	((n - 1) * n) >> 1
end


#= Many of the following functions are based on
StatsFuns.jl, with minor modifications. =#

"x * log(x), giving zero for x = 0"
xlogx(x::Real) = iszero(x) ? float(x) : x * log(x)

"x*log(y), giving zero if x is zero and y non-negative"
xlogy(x::T, y::T) where {T<:Real} = iszero(x) && y ≥ 0 ? float(x) : x * log(y)
xlogy(x::Real, y::Real) = xlogy(promote(x, y)...)

"x * exp(x), giving zero for x = -Inf"
xexpx(x::Real) = isfinite(x) ? x * exp(x) : exp(x)

"x * exp(y), giving zero for y = -Inf"
xexpy(x::T, y::T) where {T<:Real} = isnan(x) || isfinite(y) ? x * exp(y) : exp(y)
xexpy(x::Real, y::Real) = xexpy(promote(x, y)...)

"log(1+exp(x))"
log1pexp(x::Real) = x ≤ -37. ? exp(x) : x ≤ 18. ? log1p(exp(x)) : x ≤ 33.3 ? x + exp(-x) : float(x)

"log(1-exp(x))"
log1mexp(x::Real) = x < logofhalf ? log1p(-exp(x)) : log(-expm1(x))


#= this is called loghalf in StatsFuns, and I was always getting
name conflicts, so I renamed it here to logofhalf =#
Base.@irrational logofhalf -0.6931471805599453094 log(big(0.5))


"""midpoint of the interval [a,b]"""
function midpoint(a::Float64, b::Float64)
    #= Based on F. Goualard. 2014, Table VII
    DOI: 10.1145/2493882 =#

    if !(a ≤ b)
        return NaN
    elseif a == -Inf
        if b == Inf
            return 0.
        else
            return nextfloat(-Inf)
        end
    elseif b == Inf
        return prevfloat(Inf)
    end
    
    mid = 0.5*(a + b)
    if isfinite(mid)
        return mid
    else
        return 0.5 * a + 0.5 * b
    end
end
