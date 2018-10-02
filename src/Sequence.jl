export Sequence, energy

import Random

"sequence type, with all letters guaranteed to be in the range 1 ≤ a ≤ A"
struct Sequence{A,L}
    s::NTuple{L,Int}
    function Sequence{A,L}(s::NTuple{L,Integer}) where {A,L}
		@checkposint A
        for i = 1:L
			if !(1 ≤ s[i] ≤ A)
				throw(ArgumentError("sequence letters must be 1 ≤ s[i] ≤ A; got s[$i]=$(s[i]), A=$A"))
			end
		end
		
        new(s)
    end
end

Sequence{A}(s::NTuple{L,Integer}) where {A,L} = Sequence{A,L}(s)
Sequence{A,L}(s::Integer...) where {A,L} = Sequence{A,L}(s)
Sequence{A}(s::Integer...) where {A} = Sequence{A}(s)
Sequence{A,L}(s::AbstractVector{<:Integer}) where {A,L} = Sequence{A,L}(Tuple(s))

Base.convert(::Type{Sequence{A,L}}, s::NTuple{L,Int}) where {A,L} = Sequence{A,L}(s)
Base.convert(::Type{NTuple{L,Int}}, s::Sequence{A,L}) where {A,L} = s.s

Base.:(==)(s::Sequence, r::Sequence) = s.s == r.s
Base.hash(s::Sequence) = hash(s.s)

Base.getindex(s::Sequence, i) = s.s[i]
Base.length(s::Sequence) = length(s.s)
Base.eltype(::Type{<:Sequence}) = eltype(s.s)
Base.iterate(s::Sequence) = iterate(s.s)
Base.iterate(s::Sequence, state) = iterate(s.s, state)

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Sequence{A,L}}) where {A,L}
	Sequence{A,L}(ntuple(i -> rand(1:A), Val(L)))
end

"energy of sequence s using 'fields' (h,J)"
function energy(fields::Fields{A,L,U}, s::Sequence{A,L}) where {A,L,U}   
    E = zero(U)
    
    offset = 0
	@inbounds for i=1:L
		E -= fields[s[i] + offset]
		offset += A
    end
    
    @assert offset == A*L
    
    A2 = A^2;

	@inbounds for j=1:L, i=1:j-1
		E -= fields[s[i] + (s[j]-1)*A + offset]
		offset += A2
	end

    E
end
