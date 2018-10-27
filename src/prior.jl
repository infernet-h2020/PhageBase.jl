export FieldsPrior, GaussPrior, 
       log_prior, log_prior_grad!, log_prior_hess!


abstract type FieldsPrior{A,L} end


struct GaussPrior{A,L} <: FieldsPrior{A,L}
    η::Vector{Float64}
	function GaussPrior{A,L}(η::AbstractVector{Float64}) where {A,L}
		@checkposint A L
		@assert length(η) == fieldslen(A,L)
		for x in η @assert 0 ≤ x < Inf end
        new(η)
    end
end


function GaussPrior{A,L}(ηh::Real, 
                         ηJ::Real) where {A,L}
	@checkposint A L
	@assert 0 ≤ ηh < Inf
	@assert 0 ≤ ηJ < Inf

	η = [fill(ηh, A*L);
		 fill(ηJ, binom2(L)*A^2)]
	GaussPrior{A,L}(float(η))
end


GaussPrior{A,L}(η::Real) where {A,L} = GaussPrior{A,L}(η, η)


function GaussPrior{A,L}() where {A,L}
	@checkposint A L
	η = zeros(fieldslen(A,L))
	GaussPrior{A,L}(η)
end


GaussPrior(::Fields{A,L}, η::Union{Real,AbstractVector{<:Real}}) where {A,L} = GaussPrior{A,L}(η)
GaussPrior(::Fields{A,L}, ηh::Real, ηJ::Real) where {A,L} = GaussPrior{A,L}(ηh, ηJ)


"fields log prior"
function log_prior(fields::Fields{A,L,U},
				   prior::GaussPrior{A,L}
				   ) where {A,L,U}
	p = zero(U)
	for f = 1 : fieldslen(A,L)
		p -= prior.η[f] * fields[f]^2
	end
	p/2
end


"adds to G the gradient of the fields log prior"
function log_prior_grad!(G::AbstractVector{Float64},
                         fields::Fields{A,L,Float64},
                         prior::GaussPrior{A,L}
                         ) where {A,L,V,T}
    @assert length(G) == length(fields.x)
    for f = 1 : length(prior.η)    # TODO: @inbounds
        G[f] -= prior.η[f] * fields[f]
    end
    nothing
end


"adds to H the Hessian of the fields log prior"
function log_prior_hess!(H::AbstractMatrix{Float64},
                         prior::GaussPrior{A,L}
                         ) where {A,L}
    flen = fieldslen(A,L)
    @assert size(H,1) == size(H,2) ≥ flen
    for f = 1 : flen    # TODO: @inbounds
        H[f,f] -= prior.η[f]
    end
    nothing
end