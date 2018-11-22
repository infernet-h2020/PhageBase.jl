export zero_sum_gauge!,
       zero_sum_gauge,
       is_zero_sum_gauge


"""
    zero_sum_gauge!(fields)

transforms fields to zero-sum gauge
"""
zero_sum_gauge!(fields::Fields) = zero_sum_gauge_fieldsonly!(fields)

"""
    zero_sum_gauge!(fields)

transforms fields to zero-sum gauge, also changing
the chemical potential
"""
function zero_sum_gauge!(fields::FieldsChem{A,L}) where {A,L}
    Δμ = (1/A)sum(fields[a,i] for a=1:A, i=1:L)
    Δμ += (1/A^2)sum(fields[a,b,i,j] for a=1:A, b=1:A, i=1:L for j=i+1:L)
    zero_sum_gauge_fieldsonly!(fields)
    setμ!(fields, getμ(fields) + Δμ)
    nothing
end


"""
    zero_sum_gauge_fieldsonly!(fields)

Transforms fields (h,J) to zero sum gauge.
"""
function zero_sum_gauge_fieldsonly!(fields::FieldsAny{A,L}) where {A,L}
    #= This is a dumb way to do this. It can probably
    be made more efficient, but I think I never use
    this in any performance critical code. =#

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
end


function zero_sum_gauge(fields::FieldsAny)
    f = deepcopy(fields)
    zero_sum_gauge!(f)
    f
end


"""return true if the fields are in zero-sum gauge, or false otherwise.
Uses atol for comparison with zero."""
function is_zero_sum_gauge(fields::FieldsAny{A,L}; atol=1e-10) where {A,L}
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
