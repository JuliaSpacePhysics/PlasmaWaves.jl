# Electromagnetic fluid dispersion relation solver
# Ported from MATLAB bo_fluid_matrix.m by Hua-sheng XIE
#
# References:
# - Xie, H. S. (2014). PDRF: A general dispersion relation solver for
#   magnetized multi-fluid plasma. Comput. Phys. Comm. 185, 670-675.

using SparseArrays

export FluidSpecies, FluidSolverParams, solve_fluid_dispersion

"""
    FluidSpecies{T}

Fluid species parameters for the multi-fluid dispersion relation solver.

# Fields
- `q`: Charge in units of electron charge
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity (in units of c)
- `gamma_z`: Parallel polytrope exponent (default: 1.0 for isothermal)
- `gamma_p`: Perpendicular polytrope exponent (default: 1.0 for isothermal)
"""
struct FluidSpecies{T}
    q::T           # Charge (units of e)
    m::T           # Mass (units of m_p)
    n::T           # Number density (m^-3)
    Tz::T          # Parallel temperature (eV)
    Tp::T          # Perpendicular temperature (eV)
    vdz::T         # Parallel drift velocity (c)
    gamma_z::T     # Parallel polytrope exponent
    gamma_p::T     # Perpendicular polytrope exponent
end

"""
    FluidSpecies(q, m, n, Tz, Tp; vdz=0.0, gamma_z=1.0, gamma_p=1.0)

Create a FluidSpecies with specified parameters.

# Arguments
- `q`: Charge in units of e (-1 for electron, 1 for proton)
- `m`: Mass in units of proton mass
- `n`: Number density (m⁻³)
- `Tz`: Parallel temperature (eV)
- `Tp`: Perpendicular temperature (eV)
- `vdz`: Parallel drift velocity in units of c (default: 0)
- `gamma_z`: Parallel polytrope exponent (default: 1.0, isothermal)
- `gamma_p`: Perpendicular polytrope exponent (default: 1.0, isothermal)

# Fluid Closure Models
- `gamma_z = gamma_p = 1.0`: Isothermal
- `gamma_z = 3.0, gamma_p = 2.0`: CGL (Chew-Goldberger-Low) model

See also: [`solve_fluid_dispersion`](@ref)
"""
function FluidSpecies(q, m, n, Tz, Tp; vdz=0.0, gamma_z=1.0, gamma_p=1.0)
    T = promote_type(typeof(q), typeof(m), typeof(n), typeof(Tz), typeof(Tp),
                     typeof(vdz), typeof(gamma_z), typeof(gamma_p))
    FluidSpecies{T}(T(q), T(m), T(n), T(Tz), T(Tp), T(vdz), T(gamma_z), T(gamma_p))
end

"""
    FluidSolverParams{T}

Parameters for the fluid dispersion relation solver.
"""
struct FluidSolverParams{T}
    S::Int              # Number of species
    c2::T               # Speed of light squared
    qs::Vector{T}       # Charges (C)
    ms::Vector{T}       # Masses (kg)
    ns0::Vector{T}      # Number densities (m^-3)
    wcs::Vector{T}      # Cyclotron frequencies (rad/s)
    vdsz::Vector{T}     # Parallel drift velocities (m/s)
    Psz::Vector{T}      # Parallel pressures (Pa)
    Psp::Vector{T}      # Perpendicular pressures (Pa)
    rhoms::Vector{T}    # Mass densities (kg/m^3)
    csz::Vector{T}      # Parallel sound speeds (m/s)
    csp::Vector{T}      # Perpendicular sound speeds (m/s)
    gamma_z::Vector{T}  # Parallel polytrope exponents
    gamma_p::Vector{T}  # Perpendicular polytrope exponents
    B0::T               # Magnetic field (T)
end

"""
    create_fluid_params(species::Vector{<:FluidSpecies}, B0)

Create FluidSolverParams from a list of FluidSpecies.

# Arguments
- `species`: Vector of FluidSpecies objects
- `B0`: Magnetic field strength (Tesla)
"""
function create_fluid_params(species::Vector{<:FluidSpecies}, B0)
    S = length(species)
    T = Float64

    c2 = C_LIGHT^2

    qs = [s.q * Q_E for s in species]
    ms = [s.m * M_P for s in species]
    ns0 = [s.n for s in species]

    # Temperatures in Kelvin
    Tzs = [s.Tz * Q_E / K_B for s in species]
    Tps = [s.Tp * Q_E / K_B for s in species]

    # Cyclotron frequencies
    wcs = [B0 * qs[i] / ms[i] for i in 1:S]

    # Drift velocities
    vdsz = [s.vdz * C_LIGHT for s in species]

    # Pressures
    Psz = [K_B * ns0[i] * Tzs[i] for i in 1:S]
    Psp = [K_B * ns0[i] * Tps[i] for i in 1:S]

    # Mass densities
    rhoms = [ms[i] * ns0[i] for i in 1:S]

    # Polytrope exponents
    gamma_z = [s.gamma_z for s in species]
    gamma_p = [s.gamma_p for s in species]

    # Sound speeds (adiabatic model)
    csz = [sqrt(gamma_z[i] * Psz[i] / rhoms[i]) for i in 1:S]
    csp = [sqrt(gamma_p[i] * Psp[i] / rhoms[i]) for i in 1:S]

    FluidSolverParams{T}(
        S, T(c2), T.(qs), T.(ms), T.(ns0), T.(wcs), T.(vdsz),
        T.(Psz), T.(Psp), T.(rhoms), T.(csz), T.(csp),
        T.(gamma_z), T.(gamma_p), T(B0)
    )
end

"""
    solve_fluid_dispersion(params::FluidSolverParams, kx, kz)

Solve the multi-fluid electromagnetic dispersion relation using matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

The fluid model includes:
- Continuity equation: ∂n/∂t + ∇·(nv) = 0
- Momentum equation: m(∂v/∂t + v·∇v) = q(E + v×B) - ∇P/n
- Maxwell's equations for E and B

# Arguments
- `params`: FluidSolverParams with plasma parameters
- `kx`: Perpendicular wave vector component (m⁻¹)
- `kz`: Parallel wave vector component (m⁻¹)

# Returns
- Vector of complex eigenfrequencies ω (rad/s)

See also: [`FluidSpecies`](@ref), [`create_fluid_params`](@ref)
"""
function solve_fluid_dispersion(params::FluidSolverParams, kx, kz)
    (; S, c2, qs, ms, ns0, wcs, vdsz, Psz, Psp, rhoms, csz, csp, gamma_z, gamma_p, B0) = params

    # Matrix dimensions: 4 variables per species (n, vx, vy, vz) + 6 fields (Ex, Ey, Ez, Bx, By, Bz)
    SJ = 4 * S
    NN = SJ + 6

    # k dot v_drift
    kvds = kz .* vdsz  # kx*vdsx + kz*vdsz, but vdsx = 0

    # Pressure/B0 terms
    if B0 != 0
        PspB0 = Psp ./ B0
        PszB0 = Psz ./ B0
    else
        PspB0 = zeros(S)
        PszB0 = zeros(S)
    end

    # Build sparse matrix
    I_idx = Int[]
    J_idx = Int[]
    V_val = ComplexF64[]

    function add_entry!(i, j, v)
        if v != 0
            push!(I_idx, i)
            push!(J_idx, j)
            push!(V_val, ComplexF64(v))
        end
    end

    for s in 1:S
        ind = (s - 1) * 4

        # dn ~ n & v (continuity equation)
        add_entry!(ind + 1, ind + 1, kvds[s])
        add_entry!(ind + 1, ind + 2, kx * ns0[s])
        add_entry!(ind + 1, ind + 4, kz * ns0[s])

        # dv ~ n & v (momentum equation)
        add_entry!(ind + 2, ind + 1, kx * csp[s]^2 / ns0[s])
        add_entry!(ind + 4, ind + 1, kz * csz[s]^2 / ns0[s])
        add_entry!(ind + 2, ind + 2, kvds[s])
        add_entry!(ind + 3, ind + 3, kvds[s])
        add_entry!(ind + 4, ind + 4, kvds[s])
        add_entry!(ind + 3, ind + 2, -1im * wcs[s])
        add_entry!(ind + 2, ind + 3, 1im * wcs[s])

        # dv ~ E (Lorentz force)
        add_entry!(ind + 2, SJ + 1, 1im * qs[s] / ms[s])
        add_entry!(ind + 3, SJ + 2, 1im * qs[s] / ms[s])
        add_entry!(ind + 4, SJ + 3, 1im * qs[s] / ms[s])

        # dv ~ B (magnetic force and pressure anisotropy)
        add_entry!(ind + 2, SJ + 4, kz * (PszB0[s] - PspB0[s]) / rhoms[s])
        add_entry!(ind + 2, SJ + 5, -1im * qs[s] / ms[s] * vdsz[s])
        add_entry!(ind + 4, SJ + 4, kx * (PszB0[s] - PspB0[s]) / rhoms[s])
        add_entry!(ind + 3, SJ + 4, 1im * qs[s] / ms[s] * vdsz[s])
        add_entry!(ind + 3, SJ + 5, kz * (PszB0[s] - PspB0[s]) / rhoms[s])

        # dE ~ n (charge density contribution)
        add_entry!(SJ + 3, ind + 1, -1im * qs[s] * vdsz[s] / EPS_0)

        # dE ~ v (current density contribution)
        add_entry!(SJ + 1, ind + 2, -1im * qs[s] * ns0[s] / EPS_0)
        add_entry!(SJ + 2, ind + 3, -1im * qs[s] * ns0[s] / EPS_0)
        add_entry!(SJ + 3, ind + 4, -1im * qs[s] * ns0[s] / EPS_0)
    end

    # E(B): Faraday's law contribution to E equation
    add_entry!(SJ + 1, SJ + 5, c2 * kz)
    add_entry!(SJ + 2, SJ + 4, -c2 * kz)
    add_entry!(SJ + 2, SJ + 6, c2 * kx)
    add_entry!(SJ + 3, SJ + 5, -c2 * kx)

    # B(E): Faraday's law
    add_entry!(SJ + 4, SJ + 2, -kz)
    add_entry!(SJ + 5, SJ + 1, kz)
    add_entry!(SJ + 5, SJ + 3, -kx)
    add_entry!(SJ + 6, SJ + 2, kx)

    M = sparse(I_idx, J_idx, V_val, NN, NN)

    # Solve eigenvalue problem
    eigenvalues = eigen(Matrix(M)).values

    return eigenvalues
end
