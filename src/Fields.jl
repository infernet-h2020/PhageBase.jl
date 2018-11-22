export AbstractFields, Fields, field_index, fieldslen


abstract type AbstractFields{A,L,U<:Real} end


"wrapper around fields vector (h,J)"
struct Fields{A,L,U} <: AbstractFields{A,L,U}
    x::Vector{U}
    function Fields{A,L,U}(x::AbstractVector{U}) where {A, L, U<:Real}
        @boundscheck @assert length(x) == fieldslen(A,L)
        new(x)
    end
end

Fields{A,L}(x::Vector{U}) where {A, L, U<:Real} = Fields{A,L,U}(x)
Fields{A,L,U}() where {A, L, U<:Real} = Fields{A,L}(zeros(U, fieldslen(A,L)))
Fields{A,L}() where {A,L} = Fields{A,L,Float64}()


"get h[a,i]"
function Base.getindex(fields::AbstractFields, a::Int, i::Int)
    idx = field_index(fields, a, i)
    @inbounds fields.x[idx]
end

"get J[a,b,i,j]"
function Base.getindex(fields::AbstractFields, 
                       a::Int, b::Int, i::Int, j::Int)
    idx = field_index(fields, a, b, i, j)
    @inbounds fields.x[idx]
end

"set h[a,i]"
function Base.setindex!(fields::AbstractFields, value::Real, 
                        a::Int, i::Int)
    idx = field_index(fields, a, i)
    @inbounds fields.x[idx] = value
end

"set J[a,b,i,j]"
function Base.setindex!(fields::AbstractFields, value::Real,
                        a::Int, b::Int, i::Int, j::Int)
    idx = field_index(fields, a, b, i, j)
    @inbounds fields.x[idx] = value
end


Base.length(fields::AbstractFields) = length(fields.x)
Base.iterate(fields::AbstractFields, state = 1) = iterate(fields.x, state)
Base.eltype(::Type{AbstractFields{A,L,U}}) where {A,L,U<:Real} = @isdefined(U) ? U : Real

Base.firstindex(fields::AbstractFields) = firstindex(fields.x)
Base.lastindex(fields::AbstractFields)  = lastindex(fields.x)


"length of fields vector (h,J)"
@inline function fieldslen(A::Int, L::Int)
    @assert A ≥ 1 && L ≥ 1
    (L + binom2(L)A)A
end


"index of h[a,i] in fields vector"
function field_index(::AbstractFields{A,L}, a::Int, i::Int) where {A,L}
    @boundscheck @assert 1 ≤ a ≤ A
    @boundscheck @assert 1 ≤ i ≤ L
    a + (i-1)A
end

"index of J[a,b,i,j] in fields vector"
function field_index(::AbstractFields{A,L}, a::Int, b::Int, i::Int, j::Int) where {A,L}
    @boundscheck @assert 1 ≤ a ≤ A && 1 ≤ b ≤ A
    @boundscheck @assert 1 ≤ i < j ≤ L
    a + (b-1 + L + (i-1 + binom2(j-1))A)A
end

