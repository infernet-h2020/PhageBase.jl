export @checkposint, @checknonnegint
export fermi_dirac_prob, fermi_dirac_1mp, fermi_dirac_logp, fermi_dirac_l1mp
export xlogx, xlogy, xexpx, xexpy, log1pexp, log1mexp


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


"""matrix of diversities. The diversity at round t
of replicate v is defined as the number of sequences
s for which N[s,v,t] > 0."""
function diversities(N::AbstractArray{<:Real, 3})
	sum(x -> x > 0, N; dims=1)
end