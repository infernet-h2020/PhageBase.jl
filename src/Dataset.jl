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
	counts[s,v,t] == NaN indicates that this count was not
	measured. TODO: use missing?? =#
	counts::Array{C,3}
	
	#= In simulations we need counts to be integer, but in 
	the inference it might be better that it is real.
	That's why I allow C<:Real here =#

	function Dataset{A,L,V,T,C}(sequences::AbstractVector{Sequence{A,L}},
								counts::AbstractArray{C,3}) where {A,L,V,T,C<:Real}
		@checkposint A L V T
		S = length(sequences)
		if size(counts) ≠ (S,V,T)
			throw(ArgumentError(string("bad 'counts' shape; got ", size(counts), ", expected ", (S,V,T))))
		end
		for x in counts
			if !(0 ≤ x < Inf)
				throw(ArgumentError("counts must be non-negative and finite; got $x"))
			end
		end
		new(sequences, counts)
	end
end

function Dataset(sequences::AbstractVector{Sequence{A,L}}, 
				 counts::AbstractArray{C,3}) where {A,L,C<:Real}
	S,V,T = size(counts)
	Dataset{A,L,V,T,C}(sequences, counts)
end


number_of_sequences(d::Dataset) = length(d.sequences)
diversities(d::Dataset) = diversities(d.counts)


"""
	enrychments_freq(data)

frequency enrichment ratios of each sequence
in every round and replicates
"""
function enrichments_freq(d::Dataset)
	S,V,T = size(d.counts)
	freq = d.counts ./ sum(d.counts; dims=1) # sequence frequencies per round
	freq[:,:,2:T] ./ freq[:,:,1:T-1]	
	#enrich = enrichments_fold(d)
	# normalize by total counts in each round to get frequency ratios
	#enrich .* sum(d.counts[:,:,1:T-1]; dims=1) / sum(d.counts[:,:,2:T]; dims=1)
end


"""
	enrichments_fold(data)

fold enrichment ratios of each sequence
in every round and replicates
"""
function enrichments_fold(d::Dataset)
	S,V,T = size(d.counts)
	d.counts[:,:,2:T] ./ d.counts[:,:,1:T-1]
	# enrich = zeros(S,V,T-1)
	# for t = 1:T-1, v = 1:V, s = 1:S
	# 	if iszero(d.counts[s,v,t]) && iszero(d.counts[s,v,t+1])
	# 		enrich[s,v,t] = 0.
	# 	else
	# 		enrich[s,v,t] = d.counts[s,v,t+1] / d.counts[s,v,t]
	# 	end
	# end
	# return enrich
end


"""
	selectivities(data)

selectivities as defined in Eq. (1) of Boyer et al 2016 PNAS,
except that we multiply by library size to obtain a quantity that is
approximately independent of library size
"""
function selectivities(data::Dataset)
	S, V, T = size(data.counts)
	freq = data.counts ./ sum(data.counts; dims=1)
	sel = freq[:,:,2:T] ./ freq[:,:,1:T-1]
	enrich = enrichments_fold(data)

	#= Eq. (1) of Boyer et al 2016 PNAS.
	sequences not present at round t are given an θ = NaN at round t. =#
	θ = enrich ./ sum(x -> isnan(x) ? 0. : x, enrich; dims=1)
	
	#= By multiplying by library size, we get a quantity that can
	be compared across multiple libraries. For example, we can put
	train and test data in the same scatter plot. =#
	θ .* length(data.sequences)
end


"""
	mean_selectivities(data)

mean selectivities of each sequence over all
rounds and replicates, ignoring NaNs
"""
function mean_selectivities(data::Dataset)
	S, V, T = size(data.counts)
	θ = selectivities(data)
	θavg = zeros(S)
	for s = 1:S
		if all(isnan, θ[s,:,:])
			θavg[s] = NaN
		else
			θavg[s] = mean(x for x in θ[s,:,:] if !isnan(x))
		end
	end
	θavg
end


"""
	enrichscores(data, wt)

Functional scores as defined by Rubin et al 2017 Gen. Bio.,
where wt is the index of the wildtype sequence in data.sequences.
"""
function enrichscores(data::Dataset, wt::Int)
	@assert 1 ≤ wt ≤ length(data.sequences)
	@assert false # TODO: implement
end


"extracts the counts of a sequence from a dataset"
function seqcounts(data::Dataset{A,L,V,T,C}, 
				   sequence::Sequence{A,L}
				   )::Array{C,3} where {A,L,V,T,C}
	S = length(data.sequences)
	counts = zeros(C,1,V,T)
	while true
		s = findnext(data.sequences, sequence)
		if s == 0
			return counts
		else
			for t = 1:T, v = 1:V
				counts[1,v,t] += data.counts[s,v,t]
			end
		end
	end
end


"""
	randsplit(len, sublens...)

split a list of length 'len' into sublists of lengths
'sublens'. Returns lists of indices into the list to
effect the splitting.
"""
function randsplit(len::Int, sublens::Int...)
	@assert sum(sublens) == len
	for sublen in sublens
		@assert sublen ≥ 0
	end
	
	P = randperm(len)
	offset = 0
	[P[offset + 1 : (offset += sublen)]
		for (i, sublen) in enumerate(sublens)]
end


"""
	subdata(data, idx)

sub dataset corresponding to sequence indices idx
"""
subdata(data::Dataset, idx::AbstractVector{Int}) =
	Dataset(data.sequences[idx], data.counts[idx,:,:])


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


"""
	randtrain(data, trainsize, testssize)

splits data into train and tests datasets
"""
function randtrain(data::Dataset, 
				   trainsize::Integer,
				   testssize::Integer)
	S = length(data.sequences)
	@assert trainsize ≥ 0 && testssize ≥ 0 && trainsize + testssize ≤ S

	trainidx, testsidx, valididx = randsplit(S, trainsize, testssize, 
											 S - trainsize - testssize)

	traindata = subdata(data, trainidx)
	testsdata = subdata(data, testsidx)
	validdata = subdata(data, valididx)

	traindata, testsdata, validdata
end


"""
	normalize_counts!(data, N)

Normalize counts so that total counts at every round
and replicate equals N.
"""
function normalize_counts!(data::Dataset{A,L,V,T,Float64},
						   N::Real) where {A,L,V,T}
	@assert 0 < N < Inf
	for t = 1:T, v = 1:V
		N0 = sum(data.counts[:,v,t])
		data.counts[:,v,t] .*= N / N0
	end
	nothing
end


"Scale all counts by a constant factor"
function scale_counts!(data::Dataset, C::Real)
	@assert 0 < C < Inf
	for idx in eachindex(data.counts)
		data.counts[idx] *= C
	end
end