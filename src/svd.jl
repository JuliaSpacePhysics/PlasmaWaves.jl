svd_S(S) = begin
    A = @SMatrix [
        real(S[1, 1])   real(S[1, 2])   real(S[1, 3]);
        real(S[1, 2])   real(S[2, 2])   real(S[2, 3]);
        real(S[1, 3])   real(S[2, 3])   real(S[3, 3]);
        0               -imag(S[1, 2]) -imag(S[1, 3]);
        imag(S[1, 2])   0               -imag(S[2, 3]);
        imag(S[1, 3])   imag(S[2, 3])   0;
    ]
    return svd(A; full = false)
end

function svd_polarization(S::AbstractMatrix)
    _, Svals, V = svd_S(S)
    # singular values in Svals (vector length = 3)
    # take v = column 3 of V (associated with smallest singular value?) or as in your algorithm
    v = V[:, 3]
    v = v .* sign(v[3])
    # compute θ, φ
    theta = atan(sqrt(v[1]^2 + v[2]^2) / v[3])
    phi = atan(v[2], v[1])
    # compute planarity and ellipticity
    # ellipticity: ratio of the two axes of polarization ellipse * sign of polarization
    w1, w2, w3 = sort(Svals)
    planarity = 1.0 - sqrt(w1 / w3)
    ellipticity = (w2 / w3) * sign(imag(S[1, 2]))
    return (; theta, phi, planarity, ellipticity)
end

wavpol_svd(X, args...; kw...) = _transpose(_wavpol_svd, X, args...; kw...)

function _wavpol_svd(X::AbstractMatrix{T}, fs = 1; nfft = 256, noverlap = div(nfft, 2), smooth_t = _smooth_t(nfft), smooth_f = _hamming3()) where {T}
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
    power = zeros(T, nsteps, Nfreq)
    planarity = zeros(T, nsteps, Nfreq)
    waveangle = zeros(T, nsteps, Nfreq)
    ellipticity = zeros(T, nsteps, Nfreq)

    plan = plan_rfft(zeros(T, nfft, n), 1)

    SfType = SMatrix{n, n, Complex{T}}
    # tforeach(1:nsteps) do j
    Threads.@threads for j in 1:nsteps
        # tforeach(1:nsteps) do j
        @no_escape begin
            Xw = @alloc(T, nfft, n)
            Xf = @alloc(Complex{T}, Nfreq, n)
            S = @alloc(Complex{T}, n, n, Nfreq)
            Sm = @alloc(Complex{T}, n, n, Nfreq)

            start_idx = 1 + (j - 1) * noverlap
            end_idx = start_idx + nfft - 1
            Xw .= view(X, start_idx:end_idx, :) .* smooth_t
            mul!(Xf, plan, Xw)
            Xf ./= nfft # normalize
            spectral_matrix!(S, Xf)
            smooth_spectral_matrix!(Sm, S, smooth_f)
            for f in 1:Nfreq
                Sf = @views SfType(Sm[:, :, f])
                res = svd_polarization(Sf)
                power[j, f] = real(tr(Sf))
                planarity[j, f] = res.planarity
                waveangle[j, f] = res.theta
                ellipticity[j, f] = res.ellipticity
            end
        end
    end

    # Scaling power results to units with meaning
    binwidth = fs / nfft
    W = sum(smooth_t .^ 2) / nfft
    power_s = power * 2 / (binwidth * W)

    return (; indices, freqs, power = power_s, planarity, waveangle, ellipticity)
end
