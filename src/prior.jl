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


function GaussPrior{A,L}(ηh::Real, ηJ::Real) where {A,L}
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
function log_prior(fields::Union{Fields{A,L,U},FieldsChem{A,L,U}},
				   prior::GaussPrior{A,L}
				   ) where {A,L,U}
	p = zero(U)
	@inbounds for f = 1 : length(prior.η)
		p -= prior.η[f] * fields[f]^2
	end
	p/2
end


"adds to G the gradient of the fields log prior"
function log_prior_grad!(G::AbstractVector{Float64},
                         fields::Union{Fields{A,L},FieldsChem{A,L}},
                         prior::GaussPrior{A,L}
                         ) where {A,L,V,T}
    @boundscheck @assert length(G) == length(fields.x)
    @inbounds for f = 1 : length(prior.η)
        G[f] -= prior.η[f] * fields[f]
    end
    nothing
end


"adds to H the Hessian of the fields log prior"
function log_prior_hess!(H::AbstractMatrix{Float64},
                         prior::GaussPrior{A,L}
                         ) where {A,L}
    flen = fieldslen(A,L)
    @boundscheck @assert size(H,1) == size(H,2) ≥ flen
    @inbounds for f = 1 : flen
        H[f,f] -= prior.η[f]
    end
    nothing
end