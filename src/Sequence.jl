export Sequence, energy, hamming

import Random

"sequence type, with all letters guaranteed to be in the range 1 ≤ a ≤ A"
struct Sequence{A,L}
    s::NTuple{L,Int}
    function Sequence{A,L}(s::NTuple{L,Integer}) where {A,L}
		@checkposint A
		for i = 1:L
			@assert 1 ≤ s[i] ≤ A
		end
		
        new(s)
    end
end

Sequence{A,L}(s::Integer...) where {A,L} = Sequence{A,L}(s)
Sequence{A,L}(s::AbstractVector{<:Integer}) where {A,L} = Sequence{A,L}(Tuple(s))
Sequence{A}(s::Integer...) where {A} = Sequence{A}(s)
Sequence{A}(s::NTuple{L,Integer}) where {A,L} = Sequence{A,L}(s)
Sequence{A}(s::AbstractVector{<:Integer}) where {A} = Sequence{A}(Tuple(s))

Base.convert(::Type{Sequence{A,L}}, s::NTuple{L,Int}) where {A,L} = Sequence{A,L}(s)
Base.convert(::Type{NTuple{L,Int}}, s::Sequence{A,L}) where {A,L} = s.s
Base.convert(::Type{Sequence{A,L}}, s::Sequence{B,L}) where {A,B,L} = Sequence{A,L}(s.s)

Base.:(==)(s::Sequence, r::Sequence) = s.s == r.s
Base.hash(s::Sequence) = hash(s.s)

Base.getindex(s::Sequence, i) = s.s[i]
Base.length(s::Sequence) = length(s.s)
Base.eltype(s::Type{<:Sequence}) = Int
Base.iterate(s::Sequence) = iterate(s.s)
Base.iterate(s::Sequence, state) = iterate(s.s, state)

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Sequence{A,L}}) where {A,L}
	Sequence{A,L}(ntuple(i -> rand(1:A), Val(L)))
end


"""Hamming distance between two sequences. 
Same as number of positions where they differ."""
function hamming(s1::Sequence{A,L},
				 s2::Sequence{B,L}) where {A,B,L}
	d = 0
	for i = 1:L
		if s1[i] ≠ s2[i]
			d += 1
		end
	end
	return d
end