export Fields, field_index, fieldslen

struct Fields{A, L, U<:Real}
    x::Vector{U}
    function Fields{A,L,U}(x::AbstractVector{U}) where {A, L, U<:Real}
        @boundscheck @assert length(x) ≥ fieldslen(A,L,)
        new(x)
    end
end

Fields{A,L}(x::Vector{U}) where {A, L, U<:Real} = Fields{A,L,U}(x)
Fields{A,L,U}() where {A, L, U<:Real} = Fields{A,L}(zeros(U, fieldslen(A,L)))
Fields{A,L}() where {A,L} = Fields{A,L,Float64}()


"gets field at linear index idx"
Base.getindex(fields::Fields, idx::Int) = fields.x[idx]

"sets the field at linear index idx"
Base.setindex!(fields::Fields, value::Real, 
               idx::Int) = fields.x[idx] = value

"get h[a,i]"
function Base.getindex(fields::Fields, a::Int, i::Int)
    idx = field_index(fields, a, i)
    fields[idx]
end

"get J[a,b,i,j]"
function Base.getindex(fields::Fields, 
                       a::Int, b::Int, i::Int, j::Int)
    idx = field_index(fields, a, b, i, j)
    fields[idx]
end

"set h[a,i]"
function Base.setindex!(fields::Fields, value::Real, 
                        a::Int, i::Int)
    idx = field_index(fields, a, i)
    fields[idx] = value
end

"set J[a,b,i,j]"
function Base.setindex!(fields::Fields, value::Real,
                        a::Int, b::Int, i::Int, j::Int)
    idx = field_index(fields, a, b, i, j)
    fields[idx] = value
end


Base.length(fields::Fields) = length(fields.x)
Base.iterate(fields::Fields, state = 1) = iterate(fields.x, state)
Base.eltype(::Type{Fields{A,L,U}}) where {A,L,U<:Real} = @isdefined(U) ? U : Real

Base.firstindex(fields::Fields) = firstindex(fields.x)
Base.lastindex(fields::Fields)  = lastindex(fields.x)


"length of fields vector (h,J)"
@inline function fieldslen(A::Int, L::Int)
    @assert A ≥ 1 && L ≥ 1
    (L + binom2(L)A)A
end


"index of h[a,i] in fields vector"
function field_index(::Fields{A,L}, a::Int, i::Int) where {A,L}
    field_index(A, L, a, i)
end

"index of h[a,i] in fields vector"
function field_index(::Fields{A,L}, a::Int, b::Int, i::Int, j::Int) where {A,L}
    field_index(A, L, a, b, i, j)
end

"index of J[a,b,i,j] in fields vector"
function field_index(A::Int, L::Int, a::Int, i::Int)
    @boundscheck @assert 1 ≤ a ≤ A
    @boundscheck @assert 1 ≤ i ≤ L
    a + (i-1)A
end

"index of J[a,b,i,j] in fields vector"
function field_index(A::Int, L::Int, a::Int, b::Int, i::Int, j::Int)
    @boundscheck @assert 1 ≤ a ≤ A && 1 ≤ b ≤ A
    @boundscheck @assert 1 ≤ i < j ≤ L
    a + (b-1 + L + (i-1 + binom2(j-1))A)A
end
