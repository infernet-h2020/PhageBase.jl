#= Many of the following functions are based on StatsFuns.jl =#

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
