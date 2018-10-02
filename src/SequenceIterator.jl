export SequenceIterator

"an iterator over all sequences of length L with letters 1,...,A"
struct SequenceIterator{A,L}
	function SequenceIterator{A,L}() where {A,L}
		@checkposint A L
		new()
	end
end

function Base.iterate(::SequenceIterator{A,L}) where {A,L}
	state = ones(Int,L)
	return (Sequence{A,L}(state), state)
end

function Base.iterate(::SequenceIterator{A,L}, state::Vector{Int}) where {A,L}
	for i = L:-1:1
		if state[i] < A
			state[i] += 1
			for j = i+1:L
				state[j] = 1
			end
			return (Sequence{A,L}(state), state)
		end
	end
	return nothing
end

Base.eltype(::Type{SequenceIterator{A,L}}) where {A,L} = Sequence{A,L}
Base.length(::SequenceIterator{A,L}) where {A,L} = A^L
