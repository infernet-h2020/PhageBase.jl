export AbstractPrior, GaussPrior, L1Prior,
       log_prior, log_prior_grad!, log_prior_hess!


abstract type AbstractPrior{A,L} end


# **********
# Gauss Prior
# **********

struct GaussPrior{A,L} <: AbstractPrior{A,L}
    η::Vector{Float64}  # regularization weight (= 1 / prior variance)
    ξ::Vector{Float64}  # center of the Guassian prior
    function GaussPrior{A,L}(η::AbstractVector{Float64},
                             ξ::AbstractVector{Float64}) where {A,L}
		@checkposint A L
		@assert length(η) == length(ξ) == fieldslen(A,L)
        for x in η @assert 0 ≤ x < Inf end
        for x in ξ @assert -Inf < x < Inf end
        new(η, ξ)
    end
end

function GaussPrior{A,L}(η::AbstractVector{Float64}) where {A,L}
    ξ = zeros(fieldslen(A,L))
    GaussPrior{A,L}(η, ξ)
end

function GaussPrior{A,L}(ηh::Real, ηJ::Real) where {A,L}
	@checkposint A L
	@assert 0 ≤ ηh < Inf
	@assert 0 ≤ ηJ < Inf

	η = [fill(Float64(ηh), A*L);
         fill(Float64(ηJ), binom2(L)*A^2)]
	GaussPrior{A,L}(η)
end


GaussPrior{A,L}(η::Real) where {A,L} = GaussPrior{A,L}(η, η)


function GaussPrior{A,L}() where {A,L}
	@checkposint A L
    η = zeros(fieldslen(A,L))
    ξ = zeros(fieldslen(A,L))
	GaussPrior{A,L}(η, ξ)
end


GaussPrior(::Fields{A,L}, η::Union{Real,AbstractVector{<:Real}}) where {A,L} = GaussPrior{A,L}(η)
GaussPrior(::Fields{A,L}, ηh::Real, ηJ::Real) where {A,L} = GaussPrior{A,L}(ηh, ηJ)
GaussPrior(::Fields{A,L},
		   η::AbstractVector{<:Real},
		   ξ::AbstractVector{<:Real}) where {A,L} = GaussPrior{A,L}(η, ξ)


"fields log prior"
function log_prior(fields::FieldsAny{A,L,U},
				   prior::GaussPrior{A,L}) where {A,L,U}
	p = zero(U)
	@inbounds for f = 1 : length(prior.η)
		p -= prior.η[f] * (fields.x[f] - prior.ξ[f])^2
	end
	p/2
end


"adds to G the gradient of the fields log prior"
function log_prior_grad!(G::AbstractVector{Float64},
                         fields::FieldsAny{A,L,Float64},
                         prior::GaussPrior{A,L}) where {A,L,V,T}
    @boundscheck @assert length(G) == length(fields.x)
    @inbounds for f = 1 : length(prior.η)
        G[f] -= prior.η[f] * (fields.x[f] - prior.ξ[f])
    end
    nothing
end


"adds to H the Hessian of the fields log prior"
function log_prior_hess!(H::AbstractMatrix{Float64},
                         prior::GaussPrior{A,L}) where {A,L}
    flen = fieldslen(A,L)
    @boundscheck @assert size(H,1) == size(H,2) ≥ flen
    @inbounds for f = 1 : flen
        H[f,f] -= prior.η[f]
    end
    nothing
end


# **********
# L1 Prior
# **********

struct L1Prior{A,L} <: AbstractPrior{A,L}
    η::Vector{Float64}  # regularization weight
    ξ::Vector{Float64}  # center
    function L1Prior{A,L}(η::AbstractVector{Float64},
                          ξ::AbstractVector{Float64}) where {A,L}
		@checkposint A L
		@assert length(η) == length(ξ) == fieldslen(A,L)
        for x in η @assert 0 ≤ x < Inf end
        for x in ξ @assert -Inf < x < Inf end
        new(η, ξ)
    end
end

function L1Prior{A,L}(η::AbstractVector{Float64}) where {A,L}
    ξ = zeros(fieldslen(A,L))
    L1Prior{A,L}(η, ξ)
end

function L1Prior{A,L}(ηh::Real, ηJ::Real) where {A,L}
	@checkposint A L
	@assert 0 ≤ ηh < Inf
	@assert 0 ≤ ηJ < Inf

	η = [fill(Float64(ηh), A*L);
         fill(Float64(ηJ), binom2(L)*A^2)]
    L1Prior{A,L}(η)
end

L1Prior{A,L}(η::Real) where {A,L} = L1Prior{A,L}(η, η)


"fields log prior"
function log_prior(fields::Union{Fields{A,L,U}, FieldsChem{A,L,U}},
				   prior::L1Prior{A,L}) where {A,L,U}
	p = zero(U)
	@inbounds for f = 1 : length(prior.η)
		p -= prior.η[f] * abs(fields.x[f] - prior.ξ[f])
	end
	p/2
end
