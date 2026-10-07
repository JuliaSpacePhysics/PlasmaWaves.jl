"""
    spectral_matrix(Xf)

Compute the spectral matrix ``S`` defined by

```math
S_{ij}(f) = X_i(f) X_j^*(f),
```

where ``X_i(f)``=`Xf[f, i]` is the FFT of the ``i``-th component and ``*`` denotes complex conjugation.
"""
function spectral_matrix(Xf::AbstractMatrix{T}) where {T <: Complex}
    S = Array{T}(undef, size(Xf, 2), size(Xf, 2), size(Xf, 1))
    return spectral_matrix!(S, Xf)
end

function spectral_matrix!(S, Xf)
    @inbounds for f in axes(Xf, 1), j in axes(Xf, 2)
        xj = conj(Xf[f, j])
        for i in axes(Xf, 2)
            S[i, j, f] = Xf[f, i] * xj
        end
    end
    return S
end

"""
    spectral_matrix(X, dim = 1)

Compute the spectral matrix ``S(f)`` given the time series data `X` along dimension `dim`.

Returns a 3-D array of size ``n × n × (\\lfloor N/2 \\rfloor + 1)`` for `n` components of length `N`.
"""
function spectral_matrix(X::AbstractMatrix{<:Real}, dim = 1)
    Xf = rfft(X, dim)
    return dim == 1 ? spectral_matrix(Xf) : spectral_matrix(transpose(Xf))
end

# `S = Σₖ aa[k] Xf[g, :] * Xf[g, :]'` over the bins g centred on `f`.
function smoothed_spectral_matrix!(S, Xf, aa, f)
    h = length(aa) ÷ 2
    checkbounds(Xf, (f - h):(f + h), axes(Xf, 2))
    @inbounds for j in axes(Xf, 2)
        d = zero(real(eltype(S)))
        for k in eachindex(aa)
            d += aa[k] * abs2(Xf[f - h + k - 1, j])
        end
        S[j, j] = d
        for i in 1:(j - 1)
            acc = zero(eltype(S))
            for k in eachindex(aa)
                g = f - h + k - 1
                acc += aa[k] * (Xf[g, i] * conj(Xf[g, j]))
            end
            S[i, j] = acc
            S[j, i] = conj(acc)
        end
    end
    return S
end
