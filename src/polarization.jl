# polarisation analysis
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/twavpol.pro
# https://github.com/spedas/bleeding_edge/blob/master/general/science/wavpol/wavpol.pro
# https://pyspedas.readthedocs.io/en/latest/_modules/pyspedas/analysis/twavpol.html
# https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/twavpol.py

include("./svd.jl")

struct StokesVector{T}
    S0::T
    S1::T
    S2::T
    S3::T
end

"""
    polarization(S0, S1, S2, S3)
    polarization(S::StokesVector)

Compute the degree of polarization (p) from Stoke parameters or a Stokes vector.

# Reference
- [Wikipedia](https://en.wikipedia.org/wiki/Polarization_(waves))
- [Stokes parameters](https://en.wikipedia.org/wiki/Stokes_parameters)
"""
function polarization(S0, S1, S2, S3)
    return sqrt(S1^2 + S2^2 + S3^2) / S0
end

polarization(S::StokesVector) = polarization(S.S0, S.S1, S.S2, S.S3)

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
    trS2 = tr(S * S)
    trS = tr(S)
    return real((n * trS2 - trS^2) / ((n - 1) * trS^2))
end


"""
Wave normal angle is the angle between (wnx, wny) and the vertical |wnz|
Use the imaginary parts of off-diagonals.
Define:``A = Im(S₁₂), B = Im(S₁₃), C = Im(S₂₃)``
"""
function wave_normal_angle(S)
    A = imag(S[1, 2])
    B = imag(S[1, 3])
    C = imag(S[2, 3])
    aaa2 = sqrt(A^2 + B^2 + C^2)
    return if aaa2 != 0
        # Normalize contributions to get directional cosines.
        wnx = abs(C / aaa2)
        wny = -abs(B / aaa2)
        wnz = A / aaa2
        atan(sqrt(wnx^2 + wny^2), abs(wnz))
    else
        NaN
    end
end

# https://github.com/spedas/pyspedas/blob/master/pyspedas/analysis/twavpol.py#L450
# Reduced spectral leakage: The FFT spectrum becomes smoother, peaks clearer.
_smooth_t(nfft) = let xs = 0:(nfft - 1)
    @. 0.54 - 0.46 * cos(2π * (xs / nfft))
end
_hamming3() = SA[0.08, 1, 0.08]

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
- `noverlap`: Number of overlapping points between windows (default: nfft÷2)
- `smooth_t`: Time domain window function (default: Hann window)
- `smooth_f`: Frequency domain smoothing window (default: 3-point Hamming window)

# Returns
A named tuple containing:
- `indices`: Time indices for each FFT window
- `freqs`: Frequency array
- `power`: Power spectral density, normalized by frequency bin width and window function
- `degpol`: Degree of polarization [0,1]
- `waveangle`: Wave normal angle [0,π/2]
- `ellipticity`: Wave ellipticity [-1,1], negative for left-hand polarized
- `helicity`: Wave helicity

# Notes
- `smooth_f` is needed because otherwise the rank of the spectral matrix ``S̃(f)`` is 1, yielding a constant (fully polarized) result ``degpol(f) = 1``. Frequency smoothing introduces ensemble averaging, corresponding to different realizations, so ``S̃(f)`` gains fuller rank.
-  The cross-spectral density matrix ``S(f)`` is the Fourier transform of ``R(τ) = <X(t) X(t+τ)^†>`` ([Wiener-Khinchin theorem](https://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem)).

See also: [`polarization`](@ref), [`wave_normal_angle`](@ref), [`wpol_helicity`](@ref)
"""
function wavpol(X::AbstractMatrix{T}, fs = 1; nfft = 256, noverlap = div(nfft, 2), smooth_t = _smooth_t(nfft), smooth_f = _hamming3()) where {T}
    n = 3
    @assert size(X, 2) == n
    N = size(X, 1)
    Nfreq = div(nfft, 2) + 1
    freqs = (fs / nfft) * (0:(Nfreq - 1))

    # Define the number of FFT windows
    nsteps = floor(Int, (N - nfft) / noverlap) + 1
    indices = 1 .+ (0:(nsteps - 1)) * noverlap .+ div(nfft, 2)
    # normalize the smooth window for frequency smoothing
    smooth_f = smooth_f / sum(smooth_f)

    # Preallocate arrays for the results.
    power, degpol, waveangle, ellipticity, helicity =
        ntuple(_ -> zeros(T, nsteps, Nfreq), 5)

    plan = plan_rfft(zeros(T, nfft, n), 1)

    SfType = SMatrix{n, n, Complex{T}}
    Threads.@threads for j in 1:nsteps
        @no_escape begin
            Xw = @alloc(T, nfft, n)
            Xf = @alloc(Complex{T}, Nfreq, n)
            S = @alloc(Complex{T}, n, n, Nfreq)
            Sm = @alloc(Complex{T}, n, n, Nfreq)

            start_idx = 1 + (j - 1) * noverlap
            end_idx = start_idx + nfft - 1
            Xw .= view(X, start_idx:end_idx, :) .* smooth_t
            mul!(Xf, plan, Xw)
            Xf ./= sqrt(nfft) # Normalize
            spectral_matrix!(S, Xf)
            smooth_spectral_matrix!(Sm, S, smooth_f)
            # Compute the following polarization parameters from the spectral matrix ``S``:
            for f in 1:Nfreq
                Sf = SfType(view(Sm, :, :, f))
                power[j, f] = real(tr(Sf))
                degpol[j, f] = polarization(Sf)
                waveangle[j, f] = wave_normal_angle(Sf)
                helicity[j, f], ellipticity[j, f] = wpol_helicity(Sf, waveangle[j, f])
            end
        end
    end

    # Scaling power results to units with meaning
    binwidth = fs / nfft
    W = sum(smooth_t .^ 2) / nfft
    power_s = power * 2 / (binwidth * W)
    return (; indices, freqs, power = power_s, degpol, waveangle, ellipticity, helicity)
end

function twavpol end
function twavpol_svd end
