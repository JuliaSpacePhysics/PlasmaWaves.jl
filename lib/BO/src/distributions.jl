# Distribution function generation and Hermite expansion utilities
# Ported from MATLAB boarbitrary code by Hua-sheng XIE

"""
    Species

Species parameters for plasma dispersion relation solver.

# Fields
- `q`: Charge in units of electron charge (e.g., -1 for electron, 1 for proton)
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity (in units of c)
- `vdr`: Perpendicular ring beam drift velocity (in units of c)
- `aslm`: Hermite expansion coefficients (optional, computed if not provided)
"""
struct Species{T}
    q::T           # Charge (units of e)
    m::T           # Mass (units of m_p)
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    vdr::T         # Perpendicular drift velocity (c)
    aslm::Matrix{T}  # Hermite expansion coefficients
end

"""
    Species(q, m, n, Tz, Tp; vdz=0.0, vdr=0.0, distribution=:maxwellian, kwargs...)

Create a Species with automatic Hermite coefficient computation.

# Distributions
- `:maxwellian` or `:bi_maxwellian`: Drift bi-Maxwellian ring beam
- `:bi_kappa`: Bi-kappa distribution
- `:product_bi_kappa`: Product bi-kappa distribution
- `:kappa_maxwellian`: Kappa-Maxwellian distribution
- `:shell`: Shell distribution
- `:slowing_down`: Slowing down distribution

# Keyword Arguments
- `kappa`: Kappa index for kappa distributions (default: 5.5)
- `Nz`, `Nx`: Maximum Hermite indices (default: 16)

See also: [`expand_fv2d`](@ref), [`gen_fv2d`](@ref)
"""
function Species(q, m, n, Tz, Tp; vdz=0.0, vdr=0.0, distribution=:maxwellian, kwargs...)
    T = promote_type(typeof(q), typeof(m), typeof(n), typeof(Tz), typeof(Tp), typeof(vdz), typeof(vdr))

    if distribution == :maxwellian || distribution == :bi_maxwellian
        # For Maxwellian, aslm is simply 1.0
        aslm = ones(T, 1, 1)
    else
        # Generate distribution on grid and expand in Hermite basis
        aslm = hermite_coefficients(Tz, Tp, vdz, vdr, m; distribution, kwargs...)
    end

    Species{T}(T(q), T(m), T(n), T(Tz), T(Tp), T(vdz), T(vdr), aslm)
end

"""
    hermite_coefficients(Tz, Tp, vdz, vdr, mass; distribution=:maxwellian, kappa=5.5, Nz=16, Nx=16)

Compute Hermite expansion coefficients for a velocity distribution function.

Returns matrix `aslm[l+1, m+1]` of coefficients for the expansion
``f(v_z, v_\\perp) = \\sum_{l,m} a_{lm} g_l(v_z) g_m(v_\\perp)``
where ``g_n(v) = v^n \\exp(-v^2)``.
"""
function hermite_coefficients(Tz, Tp, vdz, vdr, mass;
                               distribution=:maxwellian, kappa=5.5, Nz=16, Nx=16)
    if distribution == :maxwellian || distribution == :bi_maxwellian
        return ones(1, 1)
    end

    # Generate distribution and expand
    fvdat = gen_fv2d(Tz, Tp, vdz, vdr, mass; distribution, kappa)
    aslm = expand_fv2d(fvdat; Nz, Nx)
    return aslm
end

"""
    gen_fv2d(Tz, Tp, vdz, vdr, mass; distribution=:bi_kappa, kappa=5.5)

Generate 2D velocity distribution function on a grid.

Returns a NamedTuple with fields: `vz`, `vx`, `fv`, `dvz`, `dvx`, `vtz`, `vtx`, `vdz`, `vdx`.

# Distributions
- `:bi_maxwellian`: Drift bi-Maxwellian ring beam
- `:bi_kappa`: Bi-kappa distribution
- `:product_bi_kappa`: Product bi-kappa distribution
- `:kappa_maxwellian`: Kappa-Maxwellian
- `:shell`: Shell distribution
- `:slowing_down`: Slowing down distribution
"""
function gen_fv2d(Tz, Tp, vdz_c, vdr_c, mass; distribution=:bi_kappa, kappa=5.5)
    ms = mass * M_P
    Tzs = Tz * Q_E / K_B  # eV -> K
    Tps = Tp * Q_E / K_B
    vdsz = vdz_c * C_LIGHT  # c -> m/s
    vdsr = vdr_c * C_LIGHT

    # Thermal velocities with kappa correction if needed
    if distribution == :bi_kappa
        vtzs = sqrt(2 * (1 - 1.5/kappa) * K_B * Tzs / ms)
        vtps = sqrt(2 * (1 - 1.5/kappa) * K_B * Tps / ms)
    elseif distribution == :product_bi_kappa
        vtzs = sqrt(2 * (1 - 0.5/kappa) * K_B * Tzs / ms)
        vtps = sqrt(2 * (1 - 1/kappa) * K_B * Tps / ms)
    else
        vtzs = sqrt(2 * K_B * Tzs / ms)
        vtps = sqrt(2 * K_B * Tps / ms)
    end

    vtx, vtz = vtps, vtzs
    vdx, vdz_val = vdsr, vdsz

    # Create velocity grid
    dvx = 0.05 * vtx
    dvz = 0.05 * vtz
    vz_range = range(-10*vtz + vdz_val, 10*vtz + vdz_val, step=dvz)
    vx_range = range(0, 10*vtx, step=dvx)

    vz = [z for z in vz_range, _ in vx_range]
    vx = [x for _ in vz_range, x in vx_range]

    # Generate distribution
    fv = if distribution == :bi_maxwellian || distribution == :maxwellian
        bs = vdsr / vtps
        As = exp(-bs^2) + sqrt(π) * bs * erfc(-bs)
        coef = 1 / (sqrt(π^3) * vtz * vtx^2 * As)
        @. coef * exp(-(vz - vdz_val)^2 / vtz^2 - (vx - vdx)^2 / vtx^2)
    elseif distribution == :bi_kappa
        coef = 1 / (sqrt(π^3 * kappa^3) * vtz * vtx^2) * exp(loggamma(kappa + 1) - loggamma(kappa - 0.5))
        @. coef * (1 + (vz - vdz_val)^2 / (kappa * vtz^2) + vx^2 / (kappa * vtx^2))^(-kappa - 1)
    elseif distribution == :product_bi_kappa
        kappaz, kappax = kappa, kappa
        coef = 1 / (sqrt(π^3 * kappaz) * vtz * vtx^2) * exp(loggamma(kappaz + 1) - loggamma(kappaz + 0.5))
        @. coef * (1 + (vz - vdz_val)^2 / (kappaz * vtz^2))^(-kappaz - 1) * (1 + vx^2 / (kappax * vtx^2))^(-kappax - 1)
    elseif distribution == :kappa_maxwellian
        coef = 1 / (sqrt(π^3 * kappa) * vtz * vtx^2) * exp(loggamma(kappa) - loggamma(kappa - 0.5))
        @. coef * (1 + (vz - vdz_val)^2 / (kappa * vtz^2))^(-kappa) * exp(-vx^2 / vtx^2)
    elseif distribution == :shell
        bs = vdx / vtx
        As = 2/sqrt(π) * bs * exp(-bs^2) + (2*bs^2 + 1) * erfc(-bs)
        coef = 1 / (sqrt(π^3) * vtx^3 * As)
        @. coef * exp(-(sqrt(vz^2 + vx^2) - vdx)^2 / vtx^2)
    elseif distribution == :slowing_down
        coef = 3 / (4π * log(1 + vdx^3 / vtx^3))
        r = @. sqrt(vz^2 + vx^2)
        @. coef * (r <= vdx) / (r^3 + vtx^3)
    else
        error("Unknown distribution: $distribution")
    end

    return (vz=vz, vx=vx, fv=fv, dvz=dvz, dvx=dvx, vtz=vtz, vtx=vtx, vdz=vdz_val, vdx=vdx)
end

"""
    expand_fv2d(fvdat; Nz=16, Nx=16)

Expand velocity distribution in Hermite basis.

Returns the coefficient matrix `alm[l+1, m+1]` for the power-exponential basis:
``f(v_z, v_\\perp) = c_0 \\sum_{l,m} a_{lm} ((v_z - d_z)/L_z)^l e^{-((v_z - d_z)/L_z)^2} \\cdot ((v_\\perp - d_\\perp)/L_\\perp)^m e^{-((v_\\perp - d_\\perp)/L_\\perp)^2}``
"""
function expand_fv2d(fvdat; Nz=16, Nx=16)
    (; vz, vx, fv, dvz, dvx, vtz, vtx, vdz, vdx) = fvdat

    # Scaling parameters
    Lz, Lx = vtz, vtx
    dz, dx = 0.0, 0.0

    # Hermite basis functions
    function frhol(z, l)
        arg = sqrt(2) * (z - dz) / Lz
        1 / sqrt(2^l * factorial(l) * sqrt(π)) * hermiteH(l, arg) * exp(-arg^2 / 2)
    end

    function fum(x, m)
        arg = sqrt(2) * (x - dx) / Lx
        1 / sqrt(2^m * factorial(m) * sqrt(π)) * hermiteH(m, arg) * exp(-arg^2 / 2)
    end

    # Compute raw coefficients a0lm
    a0lm = zeros(Nz + 1, Nx + 1)
    for jz in 0:Nz, jx in 0:Nx
        l, m = jz, jx
        # Numerical integration including vx < 0 region by symmetry
        sum_val = sum(fv[:, 1] .* frhol.(vz[:, 1], l) .* fum.(vx[:, 1], m))
        for ix in 2:size(fv, 2)
            sum_val += sum(fv[:, ix] .* frhol.(vz[:, ix], l) .* fum.(vx[:, ix], m))
            sum_val += sum(fv[:, ix] .* frhol.(vz[:, ix], l) .* fum.(-vx[:, ix], m))
        end
        a0lm[jz + 1, jx + 1] = sum_val * dvx * dvz * 2 / (Lx * Lz)
    end

    # Normalization coefficient
    bs = dx / Lx
    As = exp(-bs^2) + sqrt(π) * bs * erfc(-bs)
    cs0 = 1 / (sqrt(π^3) * Lz * Lx^2 * As)

    # Convert to power-exponential basis
    alm = a0lm_to_alm(a0lm / cs0)

    return alm
end

"""
    hermiteH(n, x)

Compute the physicist's Hermite polynomial H_n(x).
"""
function hermiteH(n::Int, x)
    if n == 0
        return one(x)
    elseif n == 1
        return 2x
    else
        # Recurrence: H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
        H_prev, H_curr = one(x), 2x
        for k in 1:(n-1)
            H_prev, H_curr = H_curr, 2x * H_curr - 2k * H_prev
        end
        return H_curr
    end
end

"""
    a0lm_to_alm(a0lm)

Convert Hermite function coefficients to power-exponential basis coefficients.
"""
function a0lm_to_alm(a0lm)
    lmax, mmax = size(a0lm)
    nmax = max(lmax, mmax)

    # Build Hermite polynomial coefficient matrix
    cHn0 = zeros(nmax, nmax)
    cHn0[1, 1] = 1.0  # H_0(x) = 1
    if nmax >= 2
        cHn0[2, 2] = 2.0  # H_1(x) = 2x
    end

    # H_n(x) = 2x H_{n-1}(x) - 2(n-1) H_{n-2}(x)
    for n in 2:(nmax-1)
        cHn0[n+1, 1] = -cHn0[n, 2]
        for k in 1:n
            cHn0[n+1, k+1] += 2 * cHn0[n, k] - 2 * (n - 1) * cHn0[n-1, k+1]
        end
    end

    # Normalize
    cHn = zeros(nmax, nmax)
    for n in 1:nmax
        for k in 1:n
            cHn[n, k] = cHn0[n, k] / sqrt(2^((n-1)-(k-1)) * factorial(n-1) * sqrt(π))
        end
    end

    # Transform coefficients
    alm = zeros(lmax, mmax)
    for jz in 1:lmax, jx in 1:mmax
        alm .+= a0lm[jz, jx] * (cHn[jz, 1:lmax] * cHn[jx, 1:mmax]')
    end

    return alm
end

"""
    create_solver_params(species::Vector{Species}, B0, kx, kz; N=2)

Create HHSolverParams from a list of Species.

# Arguments
- `species`: Vector of Species objects
- `B0`: Magnetic field strength (Tesla)
- `kx`: Perpendicular wave vector (m⁻¹)
- `kz`: Parallel wave vector (m⁻¹)
- `N`: Maximum harmonic number (default: 2)

See also: [`Species`](@ref), [`fDrHH`](@ref)
"""
function create_solver_params(species::Vector{<:Species}, B0, kx, kz; N=2)
    S = length(species)
    T = Float64

    c2 = C_LIGHT^2

    # Compute derived quantities for each species
    qs = [s.q * Q_E for s in species]
    ms = [s.m * M_P for s in species]
    ns = [s.n for s in species]
    Tzs = [s.Tz * Q_E / K_B for s in species]  # eV -> K
    Tps = [s.Tp * Q_E / K_B for s in species]

    vtzs = [sqrt(2 * K_B * Tzs[i] / ms[i]) for i in 1:S]
    vtps = [sqrt(2 * K_B * Tps[i] / ms[i]) for i in 1:S]

    wps = [sqrt(ns[i] * qs[i]^2 / (ms[i] * EPS_0)) for i in 1:S]
    wps2 = wps .^ 2

    wcs = [B0 * qs[i] / ms[i] for i in 1:S]
    rhocs = [sqrt(K_B * Tps[i] / ms[i]) / wcs[i] for i in 1:S]

    vdsz = [sp.vdz * C_LIGHT for sp in species]
    ds = [species[i].vdr * C_LIGHT / vtps[i] for i in 1:S]

    As = [exp(-ds[i]^2) + sqrt(π) * ds[i] * erfc(-ds[i]) for i in 1:S]

    Nss = fill(N, S)
    aslm = [s.aslm for s in species]
    msmax = [size(s.aslm, 2) - 1 for s in species]
    lsmax = [size(s.aslm, 1) - 1 for s in species]

    HHSolverParams{T}(
        S, T(c2), T.(wcs), T.(wps2), T.(rhocs),
        T(kx), T(kz), T.(vtzs), T.(vtps), T.(vdsz),
        T.(ds), T.(As), Nss, aslm, msmax, lsmax
    )
end
