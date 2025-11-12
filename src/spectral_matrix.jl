"""
    spectral_matrix(Xf)

Compute the spectral matrix ``S`` defined by

```math
S_{ij}(f) = X_i(f) X_j^*(f),
```

where ``X_i(f)``=`Xf[f, i]` is the FFT of the i-th component and * denotes complex conjugation.
"""
function spectral_matrix(Xf::AbstractMatrix{<:Complex})
    return @tullio S[f, i, j] := Xf[f, i] * conj(Xf[f, j])
end

spectral_matrix!(S, Xf) = (@tullio S[i, j, f] = Xf[f, i] * conj(Xf[f, j]))

"""
    spectral_matrix(X, window)

Compute the spectral matrix ``S(f)`` given the time series data `X`.

Returns a 3-D array of size ``N_{freq}, n, n``, where ``N_{freq} = \\lfloor N/2 \\rfloor`` 
    and `n` is the dimensionality (number of components).

# Arguments
- `X`: Matrix where each column is a component of the multivariate time series, or a vector of vectors.
- `window`: A window function (optional). If not provided, a rectangular window (no windowing) is used.
"""
function spectral_matrix(X::AbstractMatrix{<:Real})
    nfft = size(X, 1)
    # Compute FFTs and normalize
    Xf = rfft(X, 1) ./ sqrt(nfft)
    return spectral_matrix(Xf)
end

"""
    smooth_spectral_matrix!(S_smooth, S, aa)

In-place version of `smooth_spectral_matrix` that writes results to a pre-allocated array.
"""
function smooth_spectral_matrix!(S_smooth, S, aa)
    Nfreq = size(S, 3)
    M = length(aa)
    halfM = div(M, 2)
    # For boundary frequencies, copy original S
    S_smooth[:, :, 1:halfM] .= @view S[:, :, 1:halfM]
    @inbounds for f in (halfM + 1):(Nfreq - halfM), j in axes(S, 2), i in axes(S, 1)
        S_smooth[i, j, f] = sum(1:M) do k
            aa[k] * S[i, j, f - halfM + k - 1]
        end
    end
    slices = (Nfreq - halfM + 1):Nfreq
    S_smooth[:, :, slices] .= @view S[:, :, slices]
    return S_smooth
end
