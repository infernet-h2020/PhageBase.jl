export SequencesIterator

"an iterator over all sequences of length L with letters 1,...,A"
struct SequencesIterator{A,L}
	function SequencesIterator{A,L}() where {A,L}
		@checkposint A L
		new()
	end
end

function Base.iterate(::SequencesIterator{A,L}) where {A,L}
	state = ones(Int,L)
	(Sequence{A,L}(state), state)
end

function Base.iterate(::SequencesIterator{A,L}, state::Vector{Int}) where {A,L}
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

Base.eltype(::Type{SequencesIterator{A,L}}) where {A,L} = Sequence{A,L}
Base.length(::SequencesIterator{A,L}) where {A,L} = A^L
