export @checkposint, @checknonnegint
export fermi_dirac_prob, fermi_dirac_1mp, fermi_dirac_logp, fermi_dirac_l1mp
export log1pexp, xlogx

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


"accurate log(1 + exp(x))"
@inline function log1pexp(z::Real)
	x = float(z)
	x ≤ -37 ? exp(x) : x ≤ 18 ? log1p(exp(x)) : x ≤ 33.3 ? x + exp(-x) : x;
end


"x * log(x), giving zero if x == 0"
@inline xlogx(x::Real) = x == 0 ? zero(x) : x * log(x)


"Fermi-Dirac binding probability"
fermi_dirac_prob(φ::Real) = 1 / (1 + exp(-φ))
fermi_dirac_logp(φ::Real) = -log1pexp(-φ) # log(p)

fermi_dirac_1mp(φ::Real) = fermi_dirac_prob(-φ) # 1 - p
fermi_dirac_l1mp(φ::Real) = fermi_dirac_logp(-φ) # log(1-p)

# fermi_dirac_prob(μ::Real, E::Real) = fermi_dirac_prob(μ - E)
# fermi_dirac_logp(μ::Real, E::Real) = fermi_dirac_logp(μ - E)
# fermi_dirac_1mp(μ::Real, E::Real) = fermi_dirac_1mp(μ - E)
# fermi_dirac_l1mp(μ::Real, E::Real) = fermi_dirac_l1mp(μ - E)


for f in (:fermi_dirac_prob, 
		  :fermi_dirac_logp, 
		  :fermi_dirac_1mp, 
		  :fermi_dirac_l1mp)
	
	@eval begin
		($f)(μ::Real, E::Real) = ($f)(μ - E)
	end

end