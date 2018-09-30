export Dataset, selectivities, randtrain, seqcounts

import Random

"phage display experiment dataset"
struct Dataset{A, L, V, T, C<:Real}
	sequences::Vector{Sequence{A,L}}
	N::Array{C,3} # S×V×T tensor of counts, where S is sequences, V is replicates, and T is rounds
	
	#= In simulations we need N to be integer, but in the inference sometimes it is better that it is real.
	That's why I allow C<:Real here =#

	function Dataset{A,L,V,T,C}(sequences::AbstractVector{Sequence{A,L}},
								N::AbstractArray{C,3}) where {A, L, V, T, C<:Real}
		S = length(sequences)
		@assert size(N) == (S,V,T)
		@assert all(n -> 0 ≤ n < Inf, N)
		new(copy(sequences), copy(N))
	end
end

function Dataset(sequences::AbstractVector{Sequence{A,L}}, N::AbstractArray{C,3}) where {A, L, C<:Real}
	S,V,T = size(N)
	Dataset{A,L,V,T,C}(sequences, N)
end

"average selectivities of all sequences in an experiment"
function selectivities(d::Dataset)
	S,V,T = size(d.N)
	θ = zeros(S)
	for t=1:T-1, v=1:V
		normalization = (T-1)V/S * sum(d.N[s,v,t+1]/d.N[s,v,t] for s=1:S if d.N[s,v,t] > 0)
		for s=1:S
			if d.N[s,v,t] > 0
				θ[s] += d.N[s,v,t+1] / d.N[s,v,t] / normalization
			end
		end
	end
	nothing
end

"extracts the counts of a sequence from a dataset"
function seqcounts(data::Dataset{A,L,V,T,C}, 
				   sequence::Sequence{A,L})::Array{C,3} where {A,L,V,T,C}
	S = length(data.sequences)
	N = zeros(C,1,V,T)
	while true
		s = findnext(data.sequences, sequence)
		if s == 0
			return N
		else
			for t=1:T, v=1:V
				N[1,v,t] += data.N[s,v,t]
			end
		end
	end
end

"splits data into training and test datasets"
function randtrain(d::Dataset, trainsize::Integer)
	S = length(d.sequences)
	if !(0 ≤ trainsize ≤ S)
		throw(ArgumentError("trainsize ($trainsize) negative or larger than number of sequences ($S)"))
	end

	P = Random.randperm(S)

	trainidx = P[1 : trainsize]
	testsidx = P[trainsize + 1 : end]

	@assert length(trainidx) == trainsize

	traindata = Dataset(d.sequences[trainidx], d.N[trainidx,:,:])
	testsdata = Dataset(d.sequences[testsidx], d.N[testsidx,:,:])

    return traindata, testsdata
end

"splits data into training, test and validation datasets"
function randtrain(d::Dataset, trainsize::Integer, testssize::Integer)
	S = length(d.sequences)

	if trainsize < 0 || testssize < 0 || trainsize + testssize > S
		throw(ArgumentError("cannot satisfy trainsize ($trainsize), testssize ($testsize) with $S sequences"))
	end

    P = Random.randperm(S)

    trainidx = P[1 : trainsize]
	testsidx = P[trainsize + 1 : trainsize + testssize]
	valididx = P[trainsize + testssize + 1 : end]

	@assert length(trainidx) == trainsize
	@assert length(testsidx) == testssize

	traindata = Dataset(d.sequences[trainidx], d.N[trainidx,:,:])
	testsdata = Dataset(d.sequences[testsidx], d.N[testsidx,:,:])
	validdata = Dataset(d.sequences[valididx], d.N[valididx,:,:])

    return traindata, testsdata, validdata
end
