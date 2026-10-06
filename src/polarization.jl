# polarisation analysis
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/twavpol.pro
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/wavpol.pro
# https://pyspedas.readthedocs.io/en/latest/_modules/pyspedas/analysis/twavpol.html
# https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/twavpol.py

"""
    polarization(S)

Compute the degree of polarization (DOP) `p^2` from spectral matrix `S`.

```math
\\begin{aligned}
p^2  &= 1-\\frac{(tr 𝐒)^2-(tr 𝐒^2)}{(tr 𝐒)^2-n^{-1}(tr 𝐒)^2} \\\\
    &= \\frac{n(tr 𝐒^2)-(tr 𝐒)^2}{(n-1)(tr 𝐒)^2}
\\end{aligned}
```
"""
function polarization(S)
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
    wave_normal_angle(S)

Angle between the wave vector and z from the imaginary part of the spectral matrix `S` (Means 1972):
``𝐤 ∥ (Im S_{23}, -Im S_{13}, Im S_{12})``, folded into [0, π/2] since the sign of ``𝐤`` is undetermined.
NaN when ``Im S = 0``.
"""
function wave_normal_angle(S)
    A, B, C = imag(S[1, 2]), imag(S[1, 3]), imag(S[2, 3])
    return iszero(A) && iszero(B) && iszero(C) ? oftype(A, NaN) : atan(hypot(B, C), abs(A))
end

# https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/twavpol.py#L450
# Reduced spectral leakage: The FFT spectrum becomes smoother, peaks clearer.
_smooth_t(nfft) = let xs = 0:(nfft - 1)
    @. 0.54 - 0.46 * cos(2π * (xs / nfft))
end
_hamming3() = (0.08, 1, 0.08)

"""
    wavpol(X, fs=1; nfft=256, noverlap=div(nfft, 2), smooth_t=_smooth_t(nfft), smooth_f=_hamming3())

Perform polarization analysis of `n`-component time series data `X` (each column is a component) of sampling frequency `fs`.

For each FFT window (with specified overlap), the routine:
1. Applies a time-domain window function and computes the FFT to construct the spectral matrix ``S(f)``
2. Applies frequency smoothing using a window function
3. Computes wave parameters: power, degree of polarization, wave normal angle, ellipticity, and helicity

The analysis assumes the data are in a right-handed, field-aligned coordinate system 
(with Z along the ambient magnetic field).

# Keywords
- `nfft`: Number of points for FFT (default: 256)
- `noverlap`: Number of samples shared by consecutive windows (default: nfft÷2)
- `smooth_t`: Time domain window function (default: Hann window)
- `smooth_f`: Frequency domain smoothing window of odd length (default: 3-point Hamming window)

# Returns
A named tuple containing:
- `indices`: Time indices for each FFT window
- `freqs`: Frequency array
- `power`: One-sided power spectral density (input units² / Hz)
- `degpol`: Degree of polarization [0,1]
- `waveangle`: Wave normal angle [0,π/2]
- `ellipticity`: Wave ellipticity [-1,1], negative for left-hand polarized
- `helicity`: Wave helicity

# Notes
- `smooth_f` is needed because otherwise the rank of the spectral matrix ``S̃(f)`` is 1, yielding a constant (fully polarized) result ``degpol(f) = 1``. Frequency smoothing introduces ensemble averaging, corresponding to different realizations, so ``S̃(f)`` gains fuller rank. The `length(smooth_f) ÷ 2` bins at each end, where the smoothing window does not fit, report only `power`; the other outputs are NaN there.
-  The cross-spectral density matrix ``S(f)`` is the Fourier transform of ``R(τ) = <X(t) X(t+τ)^†>`` ([Wiener-Khinchin theorem](https://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem)).

See also: [`polarization`](@ref), [`wave_normal_angle`](@ref), [`wpol_helicity`](@ref)
"""
wavpol(X, args...; kw...) = _transpose(_wavpol, X, args...; kw...)

_wavpol(X, fs = 1; kw...) = _spectral_analysis(_wavpol_kernel, X, fs; kw...)

function _wavpol_kernel(S)
    helicity, ellipticity = wpol_helicity(S)
    return (; degpol = polarization(S), waveangle = wave_normal_angle(S), ellipticity, helicity)
end

# Windowed FFT → smoothed spectral matrix → `kernel(S)::NamedTuple` per (window, frequency).
function _spectral_analysis(kernel, X::AbstractMatrix{T}, fs; nfft = 256, noverlap = div(nfft, 2), smooth_t = _smooth_t(nfft), smooth_f = _hamming3()) where {T}
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
    outs = map(_ -> fill(T(NaN), nsteps, Nfreq), kernel(zeros(Complex{T}, n, n)))

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
                    power[j, f] = fold(f) * real(tr(Sf))
                    map((o, v) -> (o[j, f] = v), outs, kernel(Sf))
                else
                    power[j, f] = fold(f) * sum(i -> abs2(Xf[f, i]), 1:n)
                end
            end
        end
    end
    return (; indices, freqs, power, outs...)
end
