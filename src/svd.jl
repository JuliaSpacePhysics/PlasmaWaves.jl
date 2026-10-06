# Upper triangle of B^T B, where B is the 6×3 real matrix [Re(S); skew(Im(S))] for 3×3 Hermitian S.
# This equals the Gram matrix whose eigenvalues give SVD singular values squared.
# Squaring also squares the condition number, so it is formed in at least Float64:
# in Float32, planarity of an exactly planar wave came out as low as 0.9.
@inbounds function _gram_upper(S)
    R = promote_type(Float64, real(eltype(S)))
    r11, r22, r33 = R(real(S[1, 1])), R(real(S[2, 2])), R(real(S[3, 3]))
    r12, i12 = reim(Complex{R}(S[1, 2]))
    r13, i13 = reim(Complex{R}(S[1, 3]))
    r23, i23 = reim(Complex{R}(S[2, 3]))

    a11 = r11^2 + r12^2 + r13^2 + i12^2 + i13^2
    a22 = r12^2 + r22^2 + r23^2 + i12^2 + i23^2
    a33 = r13^2 + r23^2 + r33^2 + i13^2 + i23^2
    a12 = r11 * r12 + r12 * r22 + r13 * r23 + i13 * i23
    a13 = r11 * r13 + r12 * r23 + r13 * r33 - i12 * i23
    a23 = r12 * r13 + r22 * r23 + r23 * r33 + i12 * i13
    return a11, a22, a33, a12, a13, a23
end

"""
    Santolik()

Polarization `method` for [`wavpol`](@ref) from the singular value decomposition of the real 6×3 matrix ``[Re S; Im S]``
[santolikSingularValueDecomposition2003](@cite). Returns `degpol`, `planarity`, `waveangle` and `ellipticity`.

The wave normal is the right singular vector of the smallest singular value ``W_3``, and `ellipticity` is ``W_2 / W_1`` signed by ``Im S_{12}``.
`planarity` comes from the squared singular values ``W_i^2``, so near 1 it is resolved only to about ``ε^{1/4}`` (``10^{-4}`` in Float64).
"""
struct Santolik end

function (::Santolik)(S)
    A = _gram_upper(S)
    λ1, λ2, λ3 = _eigvals3(A...)
    k1, k2, k3 = _eigvec3(A..., λ3)
    return (;
        degpol = degree_of_polarization(S), planarity = 1 - sqrt(sqrt(max(λ3, 0) / λ1)),
        waveangle = _angle_from_z(k3, k1, k2), ellipticity = sqrt(max(λ2, 0) / λ1) * _handedness(imag(S[1, 2])),
    )
end
