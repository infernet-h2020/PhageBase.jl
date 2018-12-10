export Sequence, energy, hamming,
	   subseq, seqinsert

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
Base.convert(::Type{Sequence{A}}, s::Sequence{B,L}) where {A,B,L} = Sequence{A,L}(s.s)

Base.collect(s::Sequence{A,L}) = collect(s.s)

Base.:(==)(s::Sequence, r::Sequence) = s.s == r.s
Base.hash(s::Sequence) = hash(s.s)


function Base.isless(s1::Sequence{A,L},
					 s2::Sequence{B,L}) where {A,B,L}
	@inbounds for i = 1:L
		if s1[i] < s2[i]
			return true
		elseif s1[i] > s2[i]
			return false
		end
	end
	return false
end


Base.getindex(s::Sequence, i) = s.s[i]
Base.length(s::Sequence) = length(s.s)
Base.eltype(s::Type{<:Sequence}) = Int
Base.iterate(s::Sequence) = iterate(s.s)
Base.iterate(s::Sequence, state) = iterate(s.s, state)


function Random.rand(rng::Random.AbstractRNG,
					 ::Random.SamplerType{Sequence{A,L}}
					 ) where {A,L}
	Sequence{A,L}(ntuple(i -> rand(1:A), Val(L)))
end


"""
	hamming(s1, s2)

Hamming distance between two sequences.
"""
function hamming(s1::Sequence{A,L},
				 s2::Sequence{B,L}) where {A,B,L}
	d = 0
	@inbounds for i = 1:L
		if s1[i] ≠ s2[i]
			d += 1
		end
	end
	return d
end


"""
	subseq(s, i1, i2)

subsequence of s, from position i1 to i2
"""
function subseq(s::Sequence{A,L}, i1::Int, i2::Int) where {A,L}
	@assert 0 < i1 < i2 ≤ L
	Sequence{A}(s.s[i1:i2])
end


"""
	seqinsert(seq0, i0, seq1)

positions i0, i0+1, ..., i0+L1-1 in seq0 by seq1
"""
function seqinsert(seq0::Sequence{A,L}, i0::Int,
				   seq1::Sequence{A,L1}) where {A,L,L1}
	@assert 1 ≤ i0 ≤ L
	@assert i0 + L1 - 1 ≤ L
	snew = tuplejoin(seq0.s[1:i0-1], seq1.s, seq0.s[i0+L1:L])
	Sequence{A,L}(snew)
end
