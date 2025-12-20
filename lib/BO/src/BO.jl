module BO

using SpecialFunctions
using FFTW
using QuadGK
using LinearAlgebra

export fDrHH, HHSolverParams
export funAn, funBn, funCn, funZl, funIn
export gen_fv2d, expand_fv2d, hermite_coefficients
export Species, create_solver_params
export solve_dispersion_matrix, JPoleCoefficients, get_jpole_coefficients
export FluidSpecies, FluidSolverParams, create_fluid_params, solve_fluid_dispersion

include("constants.jl")
include("distributions.jl")

"""
    HHSolverParams{T}

Parameters for the Hermite-Hermite dispersion relation solver.

See also: [`fDrHH`](@ref)
"""
struct HHSolverParams{T}
    S::Int                      # Number of species
    c2::T                       # Speed of light squared
    wcs::Vector{T}              # Cyclotron frequencies
    wps2::Vector{T}             # Plasma frequencies squared
    rhocs::Vector{T}            # Cyclotron radii
    kx::T                       # Perpendicular wave vector
    kz::T                       # Parallel wave vector
    vtzs::Vector{T}             # Parallel thermal velocities
    vtps::Vector{T}             # Perpendicular thermal velocities
    vdsz::Vector{T}             # Parallel drift velocities
    ds::Vector{T}               # Ring beam drift parameters
    As::Vector{T}               # Normalization parameters
    Nss::Vector{Int}            # Maximum harmonic numbers
    aslm::Vector{Matrix{T}}     # Hermite expansion coefficients
    msmax::Vector{Int}          # Max perpendicular Hermite index
    lsmax::Vector{Int}          # Max parallel Hermite index
end

"""
    funIn(n)

Normalization integral I_n = ∫ v^n exp(-v²) dv / √π.

Returns Γ((n+1)/2)/√π for even n, 0 for odd n.
"""
function funIn(n::Int)
    return isodd(n) ? 0.0 : gamma(0.5 * (n + 1)) / sqrt(π)
end

"""
    funAn(n, a, d, m, p)

Perpendicular integral for J_n²: ∫ J_n²(ay) exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funAn(n::Int, a, d, m::Int, p::Int)
    m < 0 && return 0.0

    ym = 10.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    integrand(y) = besselj(n, a * y)^2 * exp(-(y - d)^2) * (y - d)^m * y^p
    result, _ = quadgk(integrand, ymin, ymax, rtol=1e-15, atol=1e-10)
    return result
end

"""
    funBn(n, a, d, m, p)

Perpendicular integral for J_n·J_n': ∫ J_n(ay) J_n'(ay) exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funBn(n::Int, a, d, m::Int, p::Int)
    m < 0 && return 0.0

    ym = 10.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    function integrand(y)
        jn = besselj(n, a * y)
        jn_deriv = 0.5 * (besselj(n - 1, a * y) - besselj(n + 1, a * y))
        return jn * jn_deriv * exp(-(y - d)^2) * (y - d)^m * y^p
    end

    result, _ = quadgk(integrand, ymin, ymax, rtol=1e-15, atol=1e-10)
    return result
end

"""
    funCn(n, a, d, m, p)

Perpendicular integral for (J_n')²: ∫ [J_n'(ay)]² exp(-(y-d)²) (y-d)^m y^p dy.
"""
function funCn(n::Int, a, d, m::Int, p::Int)
    m < 0 && return 0.0

    ym = 10.0
    ymin = max(0.0, d - ym)
    ymax = ym + d

    function integrand(y)
        jn_deriv = 0.5 * (besselj(n - 1, a * y) - besselj(n + 1, a * y))
        return 0.25 * jn_deriv^2 * exp(-(y - d)^2) * (y - d)^m * y^p
    end

    result, _ = quadgk(integrand, ymin, ymax, rtol=1e-15, atol=1e-10)
    return result
end

"""
    funZl(z, l=0, N=128)

Generalized Plasma Dispersion Function: Z_l(z) = (1/√π) ∫ v^l exp(-v²)/(v-z) dv.

Uses Weideman's FFT-based Hilbert transform method.

See also: Xie, Phys. Plasmas 20, 092125 (2013)
"""
function funZl(z::ComplexF64, l::Int=0, N::Int=128)
    F_func(v) = v^l * exp(-v^2) / sqrt(π)
    return calZ(z, F_func, N)
end

# Hilbert transform for arbitrary distribution function
function calZ(z::ComplexF64, F_func::Function, N::Int, del::Int=1)
    Z = hilbert_transform(z, F_func, N, del)
    if isnan(Z)
        Z = hilbert_transform(z + 1e-10, F_func, N, del)
    end
    return Z
end

# Weideman's FFT-based Hilbert transform
function hilbert_transform(z::ComplexF64, F_func::Function, N::Int, del::Int=1,
                           L=sqrt(N / sqrt(2)))
    if imag(z) == 0.0
        return hilbert_real_line(real(z), F_func, N, L)
    else
        return hilbert_half_plane(z, F_func, N, L, del)
    end
end

# Weideman (1995) method for real line
function hilbert_real_line(x, F_func::Function, N::Int, L)
    n = collect(-N:(N-1))
    v = L * tan.(π * (n .+ 0.5) / (2N))

    FF = F_func.(v)
    FF[isnan.(FF)] .= 0.0

    a = fft(fftshift(FF .* (L .- 1im * v)))
    a = exp.(-1im * n * π / (2N)) .* fftshift(a)
    a = reverse(1im * (sign.(n .+ 0.5) .* a)) / (2N)

    t = (L + 1im * x) / (L - 1im * x)
    h = evalpoly(t, a) / (t^N * (L - 1im * x))
    Z = h + 1im * F_func(x)

    return π * Z
end

# Weideman (1994) method for upper/lower half-plane
function hilbert_half_plane(z::ComplexF64, F_func::Function, N::Int, L, del::Int)
    M = 2N
    M2 = 2M

    k = collect((-M+1):(M-1))
    theta = k * π / M
    v = L * tan.(theta / 2)

    FF = F_func.(v)
    FF[isnan.(FF)] .= 0.0

    W = L^2 .+ v.^2
    FF = vcat([0.0], FF .* W)

    a = fft(fftshift(FF)) / M2
    a0 = a[1]
    a = reverse(a[2:(N+1)])

    z1 = imag(z) > 0 ? z : conj(z)
    t = (L + 1im * z1) / (L - 1im * z1)
    p = evalpoly(t, a)
    h = 1im * (2p / (L - 1im * z1)^2 + (a0 / L) / (L - 1im * z1))

    Z = imag(z) > 0 ? h : conj(h) + del * 2im * F_func(z)
    return π * Z
end

"""
    fDrHH(w, params)

Electromagnetic dispersion relation for hot magnetized plasma with arbitrary
velocity distributions using Hermite-Hermite expansion.

Returns det(D(ω,k)) where D is the dispersion tensor. Roots give plasma wave
dispersion relations.

See also: Xie, Phys. Plasmas 32, 060702 (2025)
"""
function fDrHH(w::ComplexF64, params::HHSolverParams)
    (; S, c2, wcs, wps2, rhocs, kx, kz, vtzs, vtps, vdsz, ds, As, Nss, aslm, msmax, lsmax) = params

    # Handle singularities
    kx_local = kx == 0.0 ? 1e-30 : kx
    kz_local = kz == 0.0 ? 1e-4 * kx_local : kz

    # Refractive index components
    nx = sqrt(c2) * kx_local / w
    nz = sqrt(c2) * kz_local / w

    # Perpendicular wavenumber parameter
    as = kx_local * rhocs * sqrt(2)

    # Initialize dielectric tensor
    Kxx, Kxy, Kxz = 1.0 + 0im, 0.0im, 0.0im
    Kyx, Kyy, Kyz = 0.0im, 1.0 + 0im, 0.0im
    Kzx, Kzy, Kzz = 0.0im, 0.0im, 1.0 + 0im

    for s in 1:S
        for n in -Nss[s]:Nss[s]
            # Perpendicular integrals
            Ans = zeros(msmax[s] + 4, 2)
            Bns = zeros(msmax[s] + 4, 2)
            Cns = zeros(msmax[s] + 4, 2)

            for m in 0:(msmax[s] + 2)
                idm = m + 2
                Ans[idm, 1] = 2.0 / As[s] * funAn(n, as[s], ds[s], m, 0)
                Ans[idm, 2] = 2.0 / As[s] * funAn(n, as[s], ds[s], m, 1)
                Bns[idm, 1] = 2.0 / As[s] * funBn(n, as[s], ds[s], m, 1)
                Bns[idm, 2] = 2.0 / As[s] * funBn(n, as[s], ds[s], m, 2)
                Cns[idm, 1] = 2.0 / As[s] * funCn(n, as[s], ds[s], m, 2)
                Cns[idm, 2] = 2.0 / As[s] * funCn(n, as[s], ds[s], m, 3)
            end

            # Parallel dispersion functions
            zetasn = (w - kz_local * vdsz[s] - n * wcs[s]) / (kz_local * vtzs[s])
            Zlns = zeros(ComplexF64, lsmax[s] + 5)
            Il = zeros(lsmax[s] + 5)

            for l in 0:(lsmax[s] + 3)
                Zlns[l+2] = funZl(zetasn, l)
                Il[l+2] = funIn(l)
            end

            # Susceptibility tensor elements
            Xsn_xx, Xsn_xy, Xsn_xz = 0.0im, 0.0im, 0.0im
            Xsn_yy, Xsn_yz, Xsn_zz = 0.0im, 0.0im, 0.0im

            vr = vtps[s] / vtzs[s]
            dr = vdsz[s] / vtzs[s]

            for l in 0:lsmax[s], m in 0:msmax[s]
                coeff = aslm[s][l+1, m+1]
                idx_m, idx_mm = m + 3, m + 1  # m+1+2 and m-1+2

                # Common terms
                nwz = n * wcs[s] / (kz_local * vtzs[s])
                nwp = n * wcs[s] / (kz_local * vtps[s])
                Zl = Zlns[l+2]
                dAm = 2Ans[idx_m, 1] - m * Ans[idx_mm, 1]
                dBm = 2Bns[idx_m, 1] - m * Bns[idx_mm, 1]
                dCm = 2Cns[idx_m, 1] - m * Cns[idx_mm, 1]
                dZl = 2Zlns[l+3] - l * Zlns[l+1]

                # XX
                term1 = (nwz * Zl - Il[l+2]) * dAm
                term2 = Ans[m+2, 2] * vr^2 * dZl
                Xsn_xx += coeff * (term1 + term2)

                # XY
                term1 = (nwz * Zl - Il[l+2]) * dBm
                term2 = Bns[m+2, 2] * vr^2 * dZl
                Xsn_xy += coeff * (term1 + term2)

                # YY
                term1 = (nwz * Zl - Il[l+2]) * dCm
                term2 = Cns[m+2, 2] * vr^2 * dZl
                Xsn_yy += coeff * (term1 + term2)

                # XZ
                term1 = nwp * (dr * Zl + Zlns[l+3]) * dAm
                term2 = Ans[m+2, 2] * vr * ((2Zlns[l+4] - l * Zlns[l+2]) + dr * dZl)
                Xsn_xz += coeff * (term1 + term2)

                # YZ
                term1 = nwp * (dr * Zl + Zlns[l+3]) * dBm
                term2 = Bns[m+2, 2] * vr * ((2Zlns[l+4] - l * Zlns[l+2]) + dr * dZl)
                Xsn_yz += coeff * (term1 + term2)

                # ZZ
                term1_inner = dr^2 * Zl + 2dr * Zlns[l+3] + Zlns[l+4]
                term1 = nwz * term1_inner * dAm

                term2_inner1 = 2dr * (Zlns[l+4] - Il[l+3]) + 2(Zlns[l+5] - Il[l+4]) -
                               l * dr * (Zlns[l+2] - Il[l+1]) - l * (Zlns[l+3] - Il[l+2])
                term2_inner2 = dr * (dr * dZl + 2(Zlns[l+4] - l * Zlns[l+2]))
                term2 = Ans[m+2, 2] * vr^2 * (term2_inner1 + term2_inner2)
                Xsn_zz += coeff * (term1 + term2)
            end

            # Scaling factors
            nwkp = n * wcs[s] / (kx_local * vtps[s])
            Xsn_xx *= nwkp^2
            Xsn_xy *= 1im * nwkp
            Xsn_xz *= nwkp
            Xsn_yz *= -1im
            Xsn_zz *= (vtzs[s] / vtps[s])^2

            # Symmetry
            Xsn_yx = -Xsn_xy
            Xsn_zx = Xsn_xz
            Xsn_zy = -Xsn_yz

            # Add to dielectric tensor
            wp2w2 = wps2[s] / w^2
            Kxx += wp2w2 * Xsn_xx
            Kxy += wp2w2 * Xsn_xy
            Kxz += wp2w2 * Xsn_xz
            Kyx += wp2w2 * Xsn_yx
            Kyy += wp2w2 * Xsn_yy
            Kyz += wp2w2 * Xsn_yz
            Kzx += wp2w2 * Xsn_zx
            Kzy += wp2w2 * Xsn_zy
            Kzz += wp2w2 * Xsn_zz
        end
    end

    # Dispersion tensor D = K - n⊗n
    Dxx = Kxx - nz^2
    Dxy = Kxy
    Dxz = Kxz + nx * nz
    Dyx = Kyx
    Dyy = Kyy - (nx^2 + nz^2)
    Dyz = Kyz
    Dzx = Kzx + nx * nz
    Dzy = Kzy
    Dzz = Kzz - nx^2

    # Determinant
    return Dxx * Dyy * Dzz + Dyx * Dzy * Dxz + Dzx * Dyz * Dxy -
           Dxz * Dyy * Dzx - Dyz * Dzy * Dxx - Dzz * Dyx * Dxy
end

include("matrix_solver.jl")
include("fluid_solver.jl")

end # module BO
