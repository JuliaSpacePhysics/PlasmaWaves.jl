function _smallest_eigenvector(a11, a22, a33, a12, a13, a23, λ)
    b11, b22, b33 = a11 - λ, a22 - λ, a33 - λ
    # Three columns of adj(A − λI); each lies in the null space. Pick the largest for stability.
    c1 = (b22 * b33 - a23^2, a23 * a13 - a12 * b33, a12 * a23 - b22 * a13)
    c2 = (a12 * b33 - a13 * a23, a13^2 - b11 * b33, b11 * a23 - a12 * a13)
    c3 = (a12 * a23 - a13 * b22, a13 * a12 - b11 * a23, b11 * b22 - a12^2)
    v = argmax(c -> c[1]^2 + c[2]^2 + c[3]^2, (c1, c2, c3))
    invn = inv(sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    return v[1] * invn, v[2] * invn, v[3] * invn
end

# Upper triangle of B^T B, where B is the 6×3 real matrix [Re(S); skew(Im(S))] for 3×3 Hermitian S.
# This equals the Gram matrix whose eigenvalues give SVD singular values squared.
@inbounds function _gram_upper(S)
    r11, r22, r33 = real(S[1, 1]), real(S[2, 2]), real(S[3, 3])
    r12, i12 = reim(S[1, 2])
    r13, i13 = reim(S[1, 3])
    r23, i23 = reim(S[2, 3])

    a11 = r11^2 + r12^2 + r13^2 + i12^2 + i13^2
    a22 = r12^2 + r22^2 + r23^2 + i12^2 + i23^2
    a33 = r13^2 + r23^2 + r33^2 + i13^2 + i23^2
    a12 = r11 * r12 + r12 * r22 + r13 * r23 + i13 * i23
    a13 = r11 * r13 + r12 * r23 + r13 * r33 - i12 * i23
    a23 = r12 * r13 + r22 * r23 + r23 * r33 + i12 * i13
    return a11, a22, a33, a12, a13, a23
end

function svd_polarization(S::AbstractMatrix)
    a11, a22, a33, a12, a13, a23 = _gram_upper(S)

    q = (a11 + a22 + a33) / 3
    p1 = a12^2 + a13^2 + a23^2
    if p1 == 0
        λ1 = max(a11, a22, a33)
        λ3 = min(a11, a22, a33)
        λ2 = a11 + a22 + a33 - λ1 - λ3
        imin = argmin((a11, a22, a33))
        v1 = imin == 1 ? one(a11) : zero(a11)
        v2 = imin == 2 ? one(a22) : zero(a22)
        v3 = imin == 3 ? one(a33) : zero(a33)
    else
        p2 = (a11 - q)^2 + (a22 - q)^2 + (a33 - q)^2 + 2p1
        p = sqrt(p2 / 6)
        b11 = (a11 - q) / p
        b22 = (a22 - q) / p
        b33 = (a33 - q) / p
        b12 = a12 / p
        b13 = a13 / p
        b23 = a23 / p
        r = (b11 * b22 * b33 + 2 * b12 * b13 * b23 - b11 * b23^2 - b22 * b13^2 - b33 * b12^2) / 2
        ϕ = acos(clamp(r, -1, 1)) / 3
        λ1 = q + 2p * cos(ϕ)
        λ3 = q + 2p * cos(ϕ + 2π / 3)
        λ2 = 3q - λ1 - λ3
        v1, v2, v3 = _smallest_eigenvector(a11, a22, a33, a12, a13, a23, λ3)
    end

    s = ifelse(v3 < 0, -one(v3), one(v3))
    v1 *= s
    v2 *= s
    v3 *= s
    theta = atan(sqrt(v1^2 + v2^2), v3)
    phi = atan(v2, v1)
    planarity = 1.0 - (max(λ3, 0) / λ1)^(1 // 4)
    ellipticity = sqrt(max(λ2, 0) / λ1) * sign(imag(S[1, 2]))
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
    smooth_f = smooth_f ./ sum(smooth_f)

    # Preallocate arrays for the results.
    power = zeros(T, nsteps, Nfreq)
    planarity = zeros(T, nsteps, Nfreq)
    waveangle = zeros(T, nsteps, Nfreq)
    ellipticity = zeros(T, nsteps, Nfreq)

    plan = plan_rfft(zeros(T, nfft, n), 1)

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
            Xf ./= nfft # normalize
            spectral_matrix!(S, Xf)
            smooth_spectral_matrix!(Sm, S, smooth_f)
            for f in 1:Nfreq
                Sf = @views Sm[:, :, f]
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
