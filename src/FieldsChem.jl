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

FieldsChem(f::Fields{A,L,U}, μ::Real) where {A,L,U} = FieldsChem{A,L,U}([f.x; μ])
Fields(f::FieldsChem{A,L,U}) where {A,L,U} = Fields{A,L,U}(f.x[1:end-1]) # this copies!

FieldsAny = Union{Fields{A,L,U}, FieldsChem{A,L,U}} where {A,L,U}


"get the chemical potential (μ)"
getμ(fields::FieldsChem) = @inbounds fields.x[end]
"set the chemical potential (μ)"
setμ!(fields::FieldsChem, μ::Real) = @inbounds fields.x[end] = μ


"energy of sequence s using 'fields' (h,J)"
function energy(fields::FieldsAny{A,L,U}, s::Sequence{A,L}) where {A,L,U}
    E = zero(U)

    offset = 0
	@inbounds for i=1:L
		E -= fields.x[s[i] + offset]
		offset += A
    end

    @assert offset == A*L

    A2 = A^2;

	@inbounds for j=1:L, i=1:j-1
		E -= fields.x[s[i] + (s[j]-1)*A + offset]
		offset += A2
	end

    E
end


"fast computation of energy of sequence s using 'fields' (h,J)"
function energy(fields::FieldsAny{A,L,U}, s::FastSeq{A,L}) where {A,L,U}
	E = zero(U)
    @inbounds for f in s.fieldidx
		E -= fields.x[f]
    end
	E
end


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
