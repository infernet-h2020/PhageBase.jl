export Fields, get_field, field_index, set_field!

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

Base.getindex(fields::Fields, inds...) = fields.x[inds...]
Base.setindex!(fields::Fields, v, inds...) = setindex!(fields.x, v, inds...)

"fast binomial(n,2)"
@inline function binom2(n::Int)
	@boundscheck @assert n ≥ 0
	((n - 1) * n) >> 1
end

"length of fields vector (h,J)"
@inline function fieldslen(A::Int, L::Int)
    @boundscheck @assert A ≥ 1 && L ≥ 1
    @inbounds (L + binom2(L)A)A
end

"index of h[a,i] in fields vector"
function field_index(::Fields{A,L}, a::Int, i::Int) where {A,L}
    @boundscheck @assert 1 ≤ a ≤ A && 1 ≤ i ≤ L
    a + (i-1)A
end

"index of J[a,b,i,j] in fields vector"
function field_index(::Fields{A,L}, a::Int, b::Int, i::Int, j::Int) where {A,L}
    @boundscheck @assert 1 ≤ a ≤ A && 1 ≤ b ≤ A && 1 ≤ i < j ≤ L
    a + (b-1 + L + (i-1 + binom2(j-1))A)A
end

"get h[a,i]"
function get_field(fields::Fields, a::Int, i::Int)
    idx = field_index(fields, a, i)
    @inbounds fields[idx]
end

"get J[a,b,i,j]"
function get_field(fields::Fields, a::Int, b::Int, i::Int, j::Int)
    idx = field_index(fields, a, b, i, j)
    @inbounds fields[idx]    
end

"set h[a,i]"
function set_field!(fields::Fields, 
                    a::Int, i::Int, 
                    val::Real)
    idx = field_index(fields,a,i)
    @inbounds fields.x[idx] = val
    nothing
end

"set J[a,b,i,j]"
function set_field!(fields::Fields, 
                    a::Int, b::Int, i::Int, j::Int,
                    val::Real)
    idx = field_index(fields,a,b,i,j)
    @inbounds fields.x[idx] = val
    nothing
end
