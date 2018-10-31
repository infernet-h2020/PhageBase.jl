export zero_sum_gauge!,
       zero_sum_gauge,
       is_zero_sum_gauge


"""transforms fields to zero-sum gauge, returning the
change in μ. If length(fields) > fieldslen(A,L), this
assumes that fields[end] is the chemical potential and
modifies it too."""
function zero_sum_gauge!(fields::Fields{A,L}) where {A,L}
    #= This is super inefficient.
    TODO: probably this can be made faster? =#
    f0 = deepcopy(fields)
    for a = 1:A, b = 1:A, i = 1:L, j = i+1:L
        fields[a,b,i,j] -= (1/A)sum(f0[x,b,i,j] for x=1:A)
        fields[a,b,i,j] -= (1/A)sum(f0[a,y,i,j] for y=1:A)
        fields[a,b,i,j] += (1/A^2)sum(f0[x,y,i,j] for x=1:A, y=1:A)
    end
    for a = 1:A, i = 1:L
        fields[a,i] -= (1/A)sum(f0[a,i] for a=1:A)
        if i < L
            fields[a,i] += (1/A)sum(f0[a,x,i,j] for x=1:A, j=i+1:L)
            fields[a,i] -= (1/A^2)sum(f0[x,y,i,j] for x=1:A, y=1:A, j=i+1:L)
        end
        if i > 1
            fields[a,i] += (1/A)sum(f0[x,a,j,i] for x=1:A, j=1:i-1)
            fields[a,i] -= (1/A^2)sum(f0[x,y,j,i] for x=1:A, y=1:A, j=1:i-1)
        end
    end
    flen = fieldslen(A,L)
    Δμ = (1/A)sum(f0[a,i] for a=1:A, i=1:L)
    Δμ += (1/A^2)sum(f0[a,b,i,j] for a=1:A, b=1:A, i=1:L for j=i+1:L)
    if length(fields) > flen
        fields[flen + 1] += Δμ
    end
    Δμ
end


function zero_sum_gauge(fields::Fields{A,L}) where {A,L}
    f = deepcopy(fields)
    zero_sum_gauge!(f)
    f
end


"""return true if the fields are in zero-sum gauge, or false otherwise.
Uses atol for comparison with zero."""
function is_zero_sum_gauge(fields::Fields{A,L}; atol=1e-10) where {A,L}
    @assert atol ≥ 0
    for i = 1:L
        if abs(sum(fields[a,i] for a=1:A)) > atol
            return false
        end
    end
    for i = 1:L, j = i+1:L
        for a = 1:A
            if abs(sum(fields[a,b,i,j] for b=1:A)) > atol
                return false
            end
        end
        for b = 1:A
            if abs(sum(fields[a,b,i,j] for a=1:A)) > atol
                return false
            end
        end
    end
    return true
end