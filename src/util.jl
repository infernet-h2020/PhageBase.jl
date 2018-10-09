export @checkposint, @checknonnegint
export fermi_dirac_prob, fermi_dirac_1mp, fermi_dirac_logp, fermi_dirac_l1mp
export xlogx, xlogy, xexpx, log1pexp, log1mexp


"throws an error if an argument is not a positive integer"
macro checkposint(X::Union{Expr,Symbol}...)
	ex = :(nothing)
	for x in X
		thisexpr = :($x isa $Int && $x > 0 || throw(ArgumentError(string($(string(x)), " must be a positive integer; got ", $x))))
		ex = :($ex; $thisexpr)
	end
	return esc(:($ex; $nothing))
end

"throws an error if an argument is not a non-negative integer"
macro checknonnegint(X::Union{Expr,Symbol}...)
	ex = :(nothing)
	for x in X
		thisexpr = :($x isa $Int && $x ≥ 0 || throw(ArgumentError(string($(string(x)), " must be a non-negative integer; got ", $x))))
		ex = :($ex; $thisexpr)
	end
	return esc(:($ex; $nothing))
end


"fast binomial(n,2)"
@inline function binom2(n::Int)
	@boundscheck @assert n ≥ 0
	((n - 1) * n) >> 1
end


"Fermi-Dirac binding probability"
fermi_dirac_prob(φ::Real) = 1 / (1 + exp(-φ))
fermi_dirac_logp(φ::Real) = -log1pexp(-φ) # log(p)

fermi_dirac_1mp(φ::Real) = fermi_dirac_prob(-φ) # 1 - p
fermi_dirac_l1mp(φ::Real) = fermi_dirac_logp(-φ) # log(1-p)


for f in (:fermi_dirac_prob, 
		  :fermi_dirac_logp, 
		  :fermi_dirac_1mp, 
		  :fermi_dirac_l1mp)
	
	@eval begin
		($f)(μ::Real, E::Real) = ($f)(μ - E)
	end

end


#= Many of the following functions are based on StatsFuns.jl =#

"x * log(x), giving zero for x = 0"
xlogx(x::Real) = iszero(x) ? float(x) : x * log(x)

"x*log(y), giving zero if x is zero and y non-negative"
xlogy(x::T, y::T) where {T<:Real} = iszero(x) && y ≥ 0 ? float(x) : x * log(y)
xlogy(x::Real, y::Real) = xlogy(promote(x, y)...)

"x * exp(x), giving zero for x = -Inf"
xexpx(x::Real) = isfinite(x) ? x * exp(x) : exp(x)

"log(1+exp(x))"
log1pexp(x::Real) = x ≤ -37. ? exp(x) : x ≤ 18. ? log1p(exp(x)) : x ≤ 33.3 ? x + exp(-x) : float(x)

Base.@irrational loghalf -0.6931471805599453094 log(big(0.5))

"log(1-exp(x))"
log1mexp(x::Real) = x < loghalf ? log1p(-exp(x)) : log(-expm1(x))
