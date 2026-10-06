"""
    Means()

Polarization `method` for [`wavpol`](@ref) reading the polarization ellipse off the rows of the spectral matrix `S`, as SPEDAS `wavpol` does.
Returns `degpol`, `waveangle`, `ellipticity` and `ellipticity_perp`.

The wave normal follows [meansUseThreedimensionalCovariance1972](@citet): ``𝐤 ∥ (Im S_{23}, -Im S_{13}, Im S_{12})``.

Row `c` gives the field vector ``𝐮 ∝ S_{c,:}^*`` up to a phase. The polarization ellipse traced by
``Re(𝐮 e^{iφ})`` has semi-axes ``\\sqrt{(|𝐮|^2 ± |𝐮^T 𝐮|)/2}`` and area ``π |Re 𝐮 × Im 𝐮|``, so its axial ratio is

```math
\\frac{2 |Re 𝐮 × Im 𝐮|}{|𝐮|^2 + |𝐮^T 𝐮|}.
```

`ellipticity_perp` applies the same to the (x, y) components. Both are signed by ``Im S_{12}``.
Each ratio is formed from numerator and denominator summed over rows, so row `c` counts ∝ ``S_{cc}`` for a pure wave.

Differences from SPEDAS, whose `helicity` is `abs(ellipticity)` and whose `ellipticity` is `ellipticity_perp`:
- SPEDAS averages the rows equally; for a parallel wave over white noise the δB_z row then holds only noise,
  and helicity tends to 2/3 of its true value however high the SNR.
- SPEDAS drops a factor 2 in the rotation angle of the (x, y) projection, so its ellipticity differs for tilted ellipses.
"""
struct Means end

"""
    Samson()

Polarization `method` for [`wavpol`](@ref) reading the polarization ellipse off the principal eigenvector ``𝐯`` of the spectral matrix,
its pure-state part [samsonCommentsDescriptionsPolarization1980](@cite). Returns the same fields as [`Means`](@ref), with the wave normal ``𝐤 ∥ Re 𝐯 × Im 𝐯``.

Isotropic noise ``σ^2 I`` leaves ``𝐯`` unchanged, so in the ensemble limit `ellipticity` is unbiased at any SNR.
On finite samples it stays closer to the truth than [`Means`](@ref) below `degpol` ≈ 0.5; above that the two agree,
and `Samson` takes about 1.4× as long.
"""
struct Samson end

# For field vector 𝐮: Re 𝐮 × Im 𝐮 as (z, x, y) components (up to sign), and the denominators of the axial ratios
# of the ellipse in 3-D and of its (x, y) projection; see `Means`.
@inline function _ellipse(u1, u2, u3)
    x = (imag(u1 * conj(u2)), imag(u2 * conj(u3)), imag(u3 * conj(u1)))
    d = abs2(u1) + abs2(u2) + abs2(u3) + abs(u1^2 + u2^2 + u3^2)
    d⊥ = abs2(u1) + abs2(u2) + abs(u1^2 + u2^2)
    return x, d, d⊥
end

function (::Means)(S)
    T = real(eltype(S))
    n = d = n⊥ = d⊥ = zero(T)
    for c in 1:3
        (x12, x23, x31), dc, d⊥c = _ellipse(S[c, 1], S[c, 2], S[c, 3])
        n += 2 * sqrt(x12^2 + x23^2 + x31^2)
        d += dc
        n⊥ += 2 * abs(x12)
        d⊥ += d⊥c
    end
    s = _handedness(imag(S[1, 2]))
    return (;
        degpol = polarization(S), waveangle = _angle_from_z(imag(S[1, 2]), imag(S[1, 3]), imag(S[2, 3])),
        ellipticity = s * n / d, ellipticity_perp = s * n⊥ / d⊥,
    )
end

function (::Samson)(S)
    A = (real(S[1, 1]), real(S[2, 2]), real(S[3, 3]), S[1, 2], S[1, 3], S[2, 3])
    (x12, x23, x31), d, d⊥ = _ellipse(_eigvec3(A..., first(_eigvals3(A...)))...)
    # x12 = (Re 𝐯 × Im 𝐯)_z is positive for right-hand rotation about z
    s = _handedness(x12)
    return (;
        degpol = polarization(S), waveangle = _angle_from_z(x12, x23, x31),
        ellipticity = s * 2 * sqrt(x12^2 + x23^2 + x31^2) / d, ellipticity_perp = s * 2 * abs(x12) / d⊥,
    )
end

# Eigenvalues λ1 ≥ λ2 ≥ λ3 of the 3×3 Hermitian matrix with diagonal a11, a22, a33 and upper triangle a12, a13, a23,
# from the trigonometric solution of the characteristic cubic.
@inline function _eigvals3(a11, a22, a33, a12, a13, a23)
    p1 = abs2(a12) + abs2(a13) + abs2(a23)
    # Exact for diagonal A: the cubic leaves an eps-level λ3, which the 4th root in Santolik planarity amplifies to ~1e-4.
    if iszero(p1)
        λ1, λ3 = max(a11, a22, a33), min(a11, a22, a33)
        return λ1, a11 + a22 + a33 - λ1 - λ3, λ3
    end
    q = (a11 + a22 + a33) / 3
    p = sqrt(((a11 - q)^2 + (a22 - q)^2 + (a33 - q)^2 + 2p1) / 6)
    b11, b22, b33 = a11 - q, a22 - q, a33 - q
    # det((A - qI) / p)
    detB = (b11 * b22 * b33 + 2real(a12 * a23 * conj(a13)) - b11 * abs2(a23) - b22 * abs2(a13) - b33 * abs2(a12)) / p^3
    ϕ = acos(clamp(detB / 2, -1, 1)) / 3
    λ1 = q + 2p * cos(ϕ)
    λ3 = q + 2p * cos(ϕ + 2π / 3)
    return λ1, 3q - λ1 - λ3, λ3
end

# Eigenvector of the eigenvalue λ, unnormalized; zero when λ is degenerate.
# (A - λI)𝐯 = 0 gives r ⋅ 𝐯 = 0 (no conjugate) for every row r, so 𝐯 ∥ rᵢ × rⱼ; take the largest for stability.
function _eigvec3(a11, a22, a33, a12, a13, a23, λ)
    r1, r2, r3 = (a11 - λ, a12, a13), (conj(a12), a22 - λ, a23), (conj(a13), conj(a23), a33 - λ)
    return argmax(_norm2, (_cross(r1, r2), _cross(r1, r3), _cross(r2, r3)))
end

_cross(x, y) = (x[2] * y[3] - x[3] * y[2], x[3] * y[1] - x[1] * y[3], x[1] * y[2] - x[2] * y[1])
_norm2(v) = abs2(v[1]) + abs2(v[2]) + abs2(v[3])
