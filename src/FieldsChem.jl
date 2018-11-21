export FieldsChem, bindφ, getμ, setμ!


"wrapper around fields vector that includes chemical potential (h,J,μ)"
struct FieldsChem{A,L,U} <: AbstractFields{A,L,U}
    x::Vector{U}
    function FieldsChem{A,L,U}(x::AbstractVector{U}) where {A, L, U<:Real}
        @boundscheck @assert length(x) == fieldslen(A,L,) + 1
        new(x)
    end
end

FieldsChem{A,L}(x::Vector{U}) where {A,L,U<:Real} = FieldsChem{A,L,U}(x)
FieldsChem{A,L,U}() where {A,L,U<:Real} = FieldsChem{A,L}(zeros(U, fieldslen(A,L) + 1))
FieldsChem{A,L}() where {A,L} = FieldsChem{A,L,Float64}()

FieldsChem(f::Fields{A,L,U}) where {A,L,U} = FieldsChem{A,L,U}(f.x)
Fields(f::FieldsChem{A,L,U}) where {A,L,U} = Fields{A,L,U}(f.x)

FieldsAny = Union{Fields{A,L}, FieldsChem{A,L}} where {A,L}


"get the chemical potential (μ)"
getμ(fields::FieldsChem) = @inbounds fields.x[end]
"set the chemical potential (μ)"
setμ!(fields::FieldsChem, μ::Real) = @inbounds fields.x[end] = μ

energy(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = energy(Fields(fields), s)

"φ = μ - E"
bindφ(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = getμ(fields) - energy(fields, s)


for f in (:fermi_dirac_prob, 
		  :fermi_dirac_logp, 
		  :fermi_dirac_1mp, 
          :fermi_dirac_l1mp)
    
	@eval begin
        ($f)(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = ($f)(bindφ(fields, s))
	end
end
