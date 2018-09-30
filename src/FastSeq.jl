export FastSeq

"representation of sequence that enables fast energy computation"
struct FastSeq{A,L,FIdx}
    sequence::Sequence{A,L}
    fieldidx::NTuple{FIdx,Int}  # indices into fields vector that contribute to the energy of this sequence

    function FastSeq{A,L,FIdx}(s::Sequence{A,L}) where {A,L,FIdx}
        @checkposint A L FIdx
        if FIdx ≠ fidxlen(L)
            throw(ArgumentError("FIdx=$FIdx inconsistent with L=$L; expected FIdx $(fidxlen(L::Int))"))
        end

        offset = 0; A2 = A^2;
        fieldidx = Vector{Int}(undef, FIdx);
        @inbounds for i = 1:L
            fieldidx[i] = s[i] + offset
            offset += A
        end
        pos = L
        @inbounds for j = 1:L
            local_offset = (s[j]-1)*A
            for i=1:j-1
                fieldidx[pos+=1] = s[i] + local_offset + offset
                offset += A2
            end
        end
        
        new{A,L,FIdx}(s, (fieldidx...,))
    end
end

FastSeq(s::Sequence{A,L}) where {A,L} = FastSeq{A,L,fidxlen(L)}(s)

Base.convert(::Type{NTuple{L,Int}}, s::FastSeq{A,L}) where {A,L} = s.sequence.s
Base.convert(::Type{Sequence{A,L}}, s::FastSeq{A,L}) where {A,L} = s.sequence

"fast computation of energy of sequence s using 'fields' (h,J)"
function energy(fields::Fields{A,L,U}, s::FastSeq{A,L}) where {A,L,U}
	E = zero(U)
    
    @inbounds @simd for f in s.fieldidx
		E -= fields[f]
    end
    
	E
end

"length of field indices vector"
@inline function fidxlen(L::Int)
	@boundscheck @checkposint L
	binom2(L+1)
end