export FieldsChem, bindφ


"wrapper around fields vector that includes chemical potential (h,J,μ)"
struct FieldsChem{A,L,U} <: AbstractFields{A,L,U}
    x::Vector{U}
    function FieldsChem{A,L,U}(x::AbstractVector{U}) where {A, L, U<:Real}
        @boundscheck @assert length(x) ≥ fieldslen(A,L,) + 1
        new(x)
    end
end

FieldsChem{A,L}(x::Vector{U}) where {A,L,U<:Real} = FieldsChem{A,L,U}(x)
FieldsChem{A,L,U}() where {A,L,U<:Real} = FieldsChem{A,L}(zeros(U, fieldslen(A,L) + 1))
FieldsChem{A,L}() where {A,L} = FieldsChem{A,L,Float64}()

FieldsChem(f::Fields{A,L,U}) where {A,L,U} = FieldsChem{A,L,U}(f.x)
Fields(f::FieldsChem{A,L,U}) where {A,L,U} = Fields{A,L,U}(f.x)

FieldsAny = Union{Fields{A,L}, FieldsChem{A,L}} where {A,L}


Base.propertynames(::FieldsChem) = (:x, :μ)


#= so that I can dispatch on the property symbol =#
Base.getproperty(fields::FieldsChem, sym::Symbol) = getproperty(fields, Val(sym))
Base.setproperty!(fields::FieldsChem, sym::Symbol, v) = setproperty!(fields, Val(sym), v)


"get the chemical potential (μ)"
Base.getproperty(fields::FieldsChem, ::Val{:μ}) = fields[end] #TODO: @inbounds
"set the chemical potential (μ)"
Base.setproperty!(fields::FieldsChem, ::Val{:μ}, μ::Real) = fields[end] = μ #TODO: @inbounds

"get the fields vector"
Base.getproperty(fields::FieldsChem, ::Val{:x}) = getfield(fields, :x) #TODO: @inbounds


energy(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = energy(Fields(fields), s)

"φ = μ - E"
bindφ(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = fields.μ - energy(fields, s)


for f in (:fermi_dirac_prob, 
		  :fermi_dirac_logp, 
		  :fermi_dirac_1mp, 
          :fermi_dirac_l1mp)
    
	@eval begin
        ($f)(fields::FieldsChem{A,L}, s::SeqAny{A,L}) where {A,L} = ($f)(bindφ(fields, s))
	end
end
