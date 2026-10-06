# polarisation analysis
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/twavpol.pro
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/wavpol.pro
# https://pyspedas.readthedocs.io/en/latest/_modules/pyspedas/analysis/twavpol.html
# https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/twavpol.py

"""
    degree_of_polarization(S)

Degree of polarization ``p^2`` of the ``n×n`` spectral matrix `S` [samsonCommentsDescriptionsPolarization1980](@cite):
1 for a pure state, 0 for ``S ∝ I``.

```math
\\begin{aligned}
p^2  &= 1-\\frac{(tr 𝐒)^2-(tr 𝐒^2)}{(tr 𝐒)^2-n^{-1}(tr 𝐒)^2} \\\\
    &= \\frac{n(tr 𝐒^2)-(tr 𝐒)^2}{(n-1)(tr 𝐒)^2}
\\end{aligned}
```
"""
function degree_of_polarization(S)
    n = size(S, 1)
    trS = zero(eltype(S))
    trS2 = zero(eltype(S))
    @inbounds for i in axes(S, 1)
        trS += S[i, i]
        for j in axes(S, 2)
            trS2 += S[i, j] * S[j, i]
        end
    end
    return real((n * trS2 - trS^2) / ((n - 1) * trS^2))
end


"""
    degree_of_polarization(S0, S1, S2, S3)

Degree of polarization ``p^2 = (S_1^2 + S_2^2 + S_3^2) / S_0^2`` from the [Stokes parameters](https://en.wikipedia.org/wiki/Stokes_parameters),
equal to [`degree_of_polarization(S)`](@ref) of the corresponding 2×2 spectral matrix.
"""
degree_of_polarization(S0, S1, S2, S3) = (S1^2 + S2^2 + S3^2) / S0^2

# Angle of the line along (x, y, z) from the z axis, in [0, π/2] since the sign of a wave vector is undetermined; NaN for the zero vector.
_angle_from_z(z, x, y) = iszero(x) && iszero(y) && iszero(z) ? oftype(z, NaN) : atan(hypot(x, y), abs(z))

# Handedness of rotation about z, never 0: at exactly perpendicular propagation it is undefined, but the axial ratio is not.
_handedness(x) = ifelse(x < 0, -one(x), one(x))

_smooth_t(nfft) = let xs = 0:(nfft - 1)
    @. 0.54 - 0.46 * cos(2π * (xs / nfft))
end
_hamming3() = (0.08, 1, 0.08)

"""
    wavpol(X, fs = 1; method = Means(), nfft = 256, noverlap = div(nfft, 2), smooth_t, smooth_f, dim = 1)

Polarization analysis of 3-component time series `X` (components along dimension 2, or 1 with `dim = 2`) sampled at `fs`.

For each FFT window, the spectral matrix ``S(f)`` of the windowed data is smoothed over frequency, scaled to unit trace and passed to `method`,
which turns it into wave parameters. Methods: [`Means`](@ref) (SPEDAS `wavpol`), [`Samson`](@ref) (principal eigenvector)
and [`Santolik`](@ref) (SVD); any callable mapping the 3×3 ``S`` to a `NamedTuple` of reals also works.

The data are assumed to be in a right-handed, field-aligned coordinate system with z along the ambient magnetic field.

# Keywords
- `noverlap`: Number of samples shared by consecutive windows
- `smooth_t`: Time-domain window of length `nfft` (default: Hamming)
- `smooth_f`: Frequency-domain smoothing weights of odd length (default: 3-point Hamming)

# Returns
A named tuple with `indices` (centre sample of each window), `freqs`, `power` (one-sided power spectral density, input units² / Hz),
and the fields of `method`, each a (window × frequency) matrix. A field name means the same quantity whichever method returns it:
- `degpol`: Degree of polarization ``p^2`` in [0, 1]; see [`degree_of_polarization`](@ref)
- `waveangle`: Angle between the wave normal and z, in [0, π/2]
- `ellipticity`: Minor-to-major axis ratio of the polarization ellipse in its own plane, in [-1, 1], negative for left-hand rotation about z
- `ellipticity_perp`: The same for the ellipse projected onto the (x, y) plane
- `planarity`: ``1 - \\sqrt{W_3 / W_1}`` from the singular values ``W_1 ≥ W_2 ≥ W_3`` of ``[Re S; Im S]``

# Notes
`smooth_f` is needed because otherwise the rank of the spectral matrix ``S̃(f)`` is 1, yielding a constant (fully polarized) result ``degpol(f) = 1``. Frequency smoothing introduces ensemble averaging, corresponding to different realizations, so ``S̃(f)`` gains fuller rank. The `length(smooth_f) ÷ 2` bins at each end, where the smoothing window does not fit, report only `power`; the other outputs are NaN there.
"""
wavpol(X, args...; kw...) = _transpose(_wavpol, X, args...; kw...)

function _wavpol(X::AbstractMatrix{T}, fs = 1; method = Means(), nfft = 256, noverlap = div(nfft, 2), smooth_t = _smooth_t(nfft), smooth_f = _hamming3()) where {T}
    n = 3
    @assert size(X, 2) == n
    N = size(X, 1)
    Nfreq = div(nfft, 2) + 1
    freqs = (fs / nfft) * (0:(Nfreq - 1))

    0 <= noverlap < nfft || throw(ArgumentError("need 0 ≤ noverlap < nfft, got noverlap = $noverlap, nfft = $nfft"))
    isodd(length(smooth_f)) || throw(ArgumentError("smooth_f needs odd length to centre on each bin, got $(length(smooth_f))"))
    step = nfft - noverlap
    nsteps = N < nfft ? 0 : fld(N - nfft, step) + 1
    indices = (1 + div(nfft, 2)) .+ step .* (0:(nsteps - 1))
    aa = map(T, smooth_f ./ sum(smooth_f))
    window = smooth_t ./ nfft # FFT normalization folded in
    h = length(aa) ÷ 2
    psd = 1 / (fs * sum(abs2, window))
    # One-sided PSD: DC and Nyquist have no negative-frequency twin to fold in.
    fold(f) = f == 1 || 2(f - 1) == nfft ? psd : 2psd

    power = zeros(T, nsteps, Nfreq)
    # The probe call only fixes the output names.
    outs = map(_ -> fill(T(NaN), nsteps, Nfreq), method(Matrix{Complex{T}}(I, n, n) / n))

    plan = plan_rfft(zeros(T, nfft, n), 1)

    Threads.@threads for j in 1:nsteps
        @no_escape begin
            Xw = @alloc(T, nfft, n)
            Xf = @alloc(Complex{T}, Nfreq, n)
            Sf = @alloc(Complex{T}, n, n)
            start = 1 + (j - 1) * step
            Xw .= view(X, start:(start + nfft - 1), :) .* window
            mul!(Xf, plan, Xw)
            for f in 1:Nfreq
                if h < f <= Nfreq - h
                    smoothed_spectral_matrix!(Sf, Xf, aa, f)
                    trS = real(tr(Sf))
                    power[j, f] = fold(f) * trS
                    # The methods raise S to up to the 8th power, which underflows Float32 for small fields.
                    Sf ./= trS
                    map((o, v) -> (o[j, f] = v), outs, method(Sf))
                else
                    power[j, f] = fold(f) * sum(i -> abs2(Xf[f, i]), 1:n)
                end
            end
        end
    end
    return (; indices, freqs, power, outs...)
end
