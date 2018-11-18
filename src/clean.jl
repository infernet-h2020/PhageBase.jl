export clean, clean!, rmzero


"""
    clean!(data)

Modifies data in-place, setting N[s,v,t] to zero if
N[s,v,t0] was zero for any 0 < t0 < t.
"""
function clean!(data::Dataset)
    S, V, T = size(data.counts)
    for s = 1:S, v = 1:V
        for t0 = 1:T-1
            if iszero(data.counts[s,v,t0])
                for t = t0 + 1 : T
                    data.counts[s,v,t] = data.counts[s,v,t0]
                end
                break
            end
        end
    end
end


"""
    clean(data)

Like clean!(data), but returns a copy of the data,
leaving the original unmodified.
"""
function clean(data::Dataset)
    d = deepcopy(data)
    clean!(d)
    d
end


"""
    rmzero(data)

Returns a copy of data where sequences with
zero counts in all rounds and replicates have
been removed.
"""
function rmzero(data::Dataset)
    S, V, T = size(data.counts)
    flag = fill(false, S)
    for s = 1:S
        for v = 1:V, t = 1:T
            if !iszero(data.counts[s,v,t])
                flag[s] = true
                break
            end
        end
    end
    
    Dataset(data.sequences[flag], data.counts[flag,:,:])
end

