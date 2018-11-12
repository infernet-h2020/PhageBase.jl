export Dataset, number_of_sequences
export selectivities, mean_selectivities,
	   enrichments_freq, enrichments_fold
export randtrain, seqcounts, normalize_counts!, scale_counts!,
	   subdata, randsplit


using Random, Statistics


"phage display experiment dataset"
struct Dataset{A, L, V, T, C<:Real}
	# all sequences in all rounds and replicates
	sequences::Vector{Sequence{A,L}}
	
	#= S×V×T tensor of counts, where S is sequences, 
	V is replicates, and T is rounds. A value 
	N[s,v,t] == NaN indicates that this count was not
	measured. TODO: use missing?? =#
	N::Array{C,3}
	
	#= In simulations we need N to be integer, but in 
	the inference it might be better that it is real.
	That's why I allow C<:Real here =#

	function Dataset{A,L,V,T,C}(sequences::AbstractVector{Sequence{A,L}},
								N::AbstractArray{C,3}) where {A,L,V,T,C<:Real}
		@checkposint A L V T
		S = length(sequences)
		if size(N) ≠ (S,V,T)
			throw(ArgumentError(string("bad N shape; got ", size(N), ", expected ", (S,V,T))))
		elseif !all(n-> 0 ≤ n < Inf, N)
			throw(ArgumentError("counts must be non-negative and finite"))
		end
		new(sequences, N)
	end
end

function Dataset(sequences::AbstractVector{Sequence{A,L}}, 
				 N::AbstractArray{C,3}) where {A,L,C<:Real}
	S,V,T = size(N)
	Dataset{A,L,V,T,C}(sequences, N)
end


number_of_sequences(d::Dataset) = length(d.sequences)
diversities(d::Dataset) = diversities(d.N)


"""
	enrychments_freq(data)

frequency enrichment ratios of each sequence
in every round and replicates
"""
function enrichments_freq(d::Dataset)
	S,V,T = size(d.N)
	enrich = enrichments_fold(d)
	# normalize by total counts in each round to get frequency ratios
	enrich .* sum(d.N[:,:,1:T-1]; dims=1) / sum(d.N[:,:,2:T]; dims=1)
end


"""
	enrichments_fold(data)

fold enrichment ratios of each sequence
in every round and replicates
"""
function enrichments_fold(d::Dataset)
	S,V,T = size(d.N)
	enrich = zeros(S,V,T-1)
	for t = 1:T-1, v = 1:V, s = 1:S
		if iszero(d.N[s,v,t]) && iszero(d.N[s,v,t+1])
			enrich[s,v,t] = 0.
		else
			enrich[s,v,t] = d.N[s,v,t+1] / d.N[s,v,t]
		end
	end
	return enrich
end


"""
	selectivities(data)

selectivities as defined in Eq. (1) of Boyer et al 2016 PNAS,
except that we divide by library size to obtain a quantity that is
approximately independent of library size
"""
function selectivities(d::Dataset)
	enrich = enrichments_fold(d)
	θ = enrich ./ sum(enrich; dims=1)  # Eq. (1) of Boyer et al 2016 PNAS
	θ .* length(d.sequences) # divide by library size
end


"""
	mean_selectivities(data)

mean selectivities of each sequence over all rounds and replicates
"""
mean_selectivities(d::Dataset) = vec(mean(selectivities(d); dims=(2,3)))


"extracts the counts of a sequence from a dataset"
function seqcounts(data::Dataset{A,L,V,T,C}, 
				   sequence::Sequence{A,L}
				   )::Array{C,3} where {A,L,V,T,C}
	S = length(data.sequences)
	N = zeros(C,1,V,T)
	while true
		s = findnext(data.sequences, sequence)
		if s == 0
			return N
		else
			for t = 1:T, v = 1:V
				N[1,v,t] += data.N[s,v,t]
			end
		end
	end
end


"""randomly split a list of length 'len' into sublists of lengths
'sublens'. Returns lists of indices into the list to effect the splitting."""
function randsplit(len::Int, sublens::Int...)
	@assert sum(sublens) == len
	for sublen in sublens @assert sublen ≥ 0 end
	P = randperm(len)
	offset = 0
	[P[offset + 1 : (offset += sublen)] for (i, sublen) in enumerate(sublens)]
end


"dataset corresponding to sequence indices idx"
subdata(d::Dataset, idx::AbstractVector{Int}) = Dataset(d.sequences[idx], d.N[idx,:,:])


"splits data into training and test datasets"
function randtrain(d::Dataset,
				   trainsize::Integer)
	S = length(d.sequences)
	@assert 0 ≤ trainsize ≤ S

	trainidx, testsidx = randsplit(S, trainsize, S - trainsize)

	traindata = subdata(d, trainidx)
	testsdata = subdata(d, testsidx)

	traindata, testsdata
end


"splits data into training, test and validation datasets"
function randtrain(d::Dataset, 
				   trainsize::Integer,
				   testssize::Integer)
	S = length(d.sequences)
	@assert trainsize ≥ 0 && testssize ≥ 0 && trainsize + testssize ≤ S

	trainidx, testsidx, valididx = randsplit(S, trainsize, testssize, 
											 S - trainsize - testssize)

	traindata = subdata(d, trainidx)
	testsdata = subdata(d, testsidx)
	validdata = subdata(d, valididx)

	traindata, testsdata, validdata
end


"Normalize counts so that Ntot at every round equals N"
function normalize_counts!(data::Dataset{A,L,V,T,Float64},
						   N::Real) where {A,L,V,T}
	@assert 0 < N < Inf
	for t = 1:T, v = 1:V
		N0 = sum(data.N[:,v,t])
		data.N[:,v,t] .*= N / N0
	end
	nothing
end


"Scale all counts by a constant factor"
function scale_counts!(data::Dataset, C::Real) where {A,L,V,T}
	@assert 0 < C < Inf
	for idx in eachindex(data.N)
		data.N[idx] *= C
	end
end