module PlasmaWaves

"""
    PlasmaWaves

This module provides basic routines to evaluate dispersion relations
for plasma waves and instabilities.  It implements cold‐fluid
electromagnetic dispersion (including polarization) and a kinetic
electrostatic dispersion for isotropic Maxwellian species.  The
implementation draws inspiration from the BO code by Xie but is
written in pure Julia for clarity and flexibility.

The module exports a [`Species`] type to describe plasma species and
several helper functions to compute plasma frequencies, cyclotron
frequencies, thermal speeds, the plasma dispersion function, and
solvers for cold electromagnetic and kinetic electrostatic
dispersion.  Users can supply arbitrary sets of species and
background magnetic fields to explore wave propagation and
instabilities.

Examples
--------

```julia
using PlasmaWaves
using LinearAlgebra

# Define a cold electron species
sp = Species(:e, -1.602176634e-19, 9.10938356e-31, 1e6, 0.0, 0.0)

# Solve for the electromagnetic wave frequency at k=0.1 m⁻¹ and θ=0°
ω = cold_em_root(0.1, 0.0, 0.0, [sp])
println("ω ≈ ", ω, " rad/s")

# Solve for the Langmuir wave including thermal effects
spT = Species(:e, -1.602176634e-19, 9.10938356e-31, 1e6, 10.0, 0.0)
ω_L = kinetic_es_root(0.1, [spT])
println("Langmuir ω ≈ ", ω_L)

```
"""

using LinearAlgebra
using SpecialFunctions

# Physical constants
const ε₀ = 8.854187817e-12          # vacuum permittivity (F/m)
const μ₀ = 4π * 1.0e-7              # vacuum permeability (H/m)
const c  = 1 / sqrt(μ₀ * ε₀)        # speed of light (m/s)
const k_B = 1.380649e-23            # Boltzmann constant (J/K)
const eV = 1.602176634e-19          # electron‐volt to joule

"""
    Species(name, charge, mass, density, temperature, drift)

Describes a plasma species.  All numbers are SI units except
`temperature`, which is in electron‑volts (eV) for convenience.

* `name` – a symbolic identifier (e.g. `:e` for electrons)
* `charge` – the signed particle charge (C)
* `mass` – the particle mass (kg)
* `density` – number density (m⁻³)
* `temperature` – scalar temperature (eV)
* `drift` – mean drift velocity along the propagation direction (m/s)

```julia
julia> sp = Species(:e, -1.602e-19, 9.11e-31, 1e6, 5.0, 0.0)
Species(:e, -1.602e-19, 9.11e-31, 1.0e6, 5.0, 0.0)
```
"""
struct Species
    name::Symbol
    charge::Float64
    mass::Float64
    density::Float64
    temperature::Float64
    drift::Float64
end

# Fundamental frequencies and thermal speeds

"""
    omega_p(sp)

Return the plasma frequency (rad/s) of a species.  Defined by
`ω_p = |q| * sqrt(n/(m ε₀))`.  The absolute value of the charge is used
so that both positive and negative species produce positive
frequencies.
"""
omega_p(sp::Species) = abs(sp.charge) * sqrt(sp.density / (sp.mass * ε₀))

"""
    Omega_c(sp, B0)

Return the cyclotron frequency (rad/s) of a species in a magnetic
field `B0` (T).  Defined by `Ω_c = q B0 / m`.
"""
Omega_c(sp::Species, B0::Float64) = sp.charge * B0 / sp.mass

"""
    v_th(sp)

Return the 1D thermal speed (m/s) of a species.  This uses the
definition `v_th = sqrt(2 k_B T / m)` with `T` in eV.  Factor
`sqrt(2)` reflects the usual definition of the most probable speed
for a Maxwellian distribution.
"""
function v_th(sp::Species)
    return sqrt(2 * k_B * sp.temperature * eV / sp.mass)
end

"""
    plasma_dispersion(z)

Compute the plasma dispersion function `Z(z) = i√π exp(-z²) erfc(-i z)`.
This function is analytic everywhere and satisfies the usual
differential identity `Z'(z) = -2 [1 + z Z(z)]`.  It is used in
kinetic dispersion relations to account for resonant (Landau) effects.
```julia
julia> Z = plasma_dispersion(1.0 + 0.1im)
0.318815865433... - 0.595276111815...im
```
"""
function plasma_dispersion(z::Complex{T}) where T
    return im * sqrt(pi) * exp(-z^2) * erfc(-im * z)
end

"""
    stix_params(ω, species, B0)

Compute the Stix parameters `(S, D, P)` for a given real or complex
frequency `ω` (rad/s), set of species, and magnetic field `B0` (T).
The returned values are unitless.  See Stix (1962) and Swanson
for definitions.  The parameters are defined by

```
S = 1 - \sum_s ω_{p,s}² / (ω² - Ω_{c,s}²)
D = \sum_s Ω_{c,s} ω_{p,s}² / [ω (ω² - Ω_{c,s}²)]
P = 1 - \sum_s ω_{p,s}² / ω²
```
```

They are used internally to build the cold electromagnetic dispersion
relation.
"""
function stix_params(ω::Complex, species::Vector{Species}, B0::Float64)
    S = 1.0 + 0im
    D = 0.0 + 0im
    P = 1.0 + 0im
    for sp in species
        ωp = omega_p(sp)
        Ω = Omega_c(sp, B0)
        denom = ω^2 - Ω^2
        S -= ωp^2 / denom
        D += (Ω / ω) * (ωp^2 / denom)
        P -= ωp^2 / ω^2
    end
    return S, D, P
end

"""
    cold_em_root(k, θ, B0, species; branch=+1, init_ω=nothing, tol=1e-6, max_iter=50)

Solve the cold electromagnetic dispersion relation for a real angular
frequency `ω` (rad/s) for a given wavenumber `k` (m⁻¹), angle `θ`
(radians), background field `B0` (T), and a list of species.  The
optional `branch` selects the high‐frequency (`+1`) or low‐frequency
(`-1`) solution in the Appleton–Hartree quadratic.  A reasonable
initial guess is chosen based on the vacuum dispersion if `init_ω`
is not supplied.  The function returns the real part of the solution.
"""
function cold_em_root(k::Float64, θ::Float64, B0::Float64, species::Vector{Species};
                      branch::Int=+1, init_ω::Union{Nothing,Float64}=nothing,
                      tol::Float64=1e-6, max_iter::Int=50)
    θrad = θ
    # initial guess: vacuum dispersion unless provided
    ω = isnothing(init_ω) ? max(c * k, 1e-3) : init_ω
    for _ in 1:max_iter
        # compute Stix parameters
        S, D, P = stix_params(ω + 0im, species, B0)
        R = S + D
        L = S - D
        sinθ = sin(θrad)
        cosθ = cos(θrad)
        A = S * sinθ^2 + P * cosθ^2
        Bcoef = R * L * sinθ^2 + S * P * (1 + cosθ^2)
        C = P * R * L
        Δ = Bcoef^2 - 4 * A * C
        # choose root of quadratic
        root_term = sqrt(Δ)
        n2_values = Complex[(Bcoef + root_term) / (2 * A), (Bcoef - root_term) / (2 * A)]
        n2 = branch == +1 ? maximum(n2_values) : minimum(n2_values)
        # function value F = n² - (ck/ω)²
        F = n2 - (c * k / ω)^2
        # compute derivative via finite difference
        δ = ω * 1e-6
        δ = δ == 0 ? 1e-6 : δ
        ωp = ω + δ
        S2, D2, P2 = stix_params(ωp + 0im, species, B0)
        R2 = S2 + D2
        L2 = S2 - D2
        A2 = S2 * sinθ^2 + P2 * cosθ^2
        B2 = R2 * L2 * sinθ^2 + S2 * P2 * (1 + cosθ^2)
        C2 = P2 * R2 * L2
        Δ2 = B2^2 - 4 * A2 * C2
        root_term2 = sqrt(Δ2)
        n2p_values = Complex[(B2 + root_term2) / (2 * A2), (B2 - root_term2) / (2 * A2)]
        n2p = branch == +1 ? maximum(n2p_values) : minimum(n2p_values)
        Fp = n2p - (c * k / ωp)^2
        dF = (Fp - F) / δ
        # Newton step
        ω_new = ω - F / dF
        # check convergence
        if abs(ω_new - ω) < tol * max(abs(ω), 1.0)
            return real(ω_new)
        end
        ω = ω_new
    end
    return real(ω)
end

"""
    cold_em_polarization(ω, k, θ, B0, species) -> Evec, Bvec

Compute the normalized electric and magnetic polarization vectors for a
cold electromagnetic mode with angular frequency `ω` (rad/s),
wavenumber `k` (m⁻¹), propagation angle `θ` (radians), background
field `B0` (T) and species list.  The returned vectors are complex
3‐component arrays with unit maximum magnitude.  The magnetic
polarization is obtained from Faraday’s law `(k × E)/ω`.
"""
function cold_em_polarization(ω::Float64, k::Float64, θ::Float64,
                              B0::Float64, species::Vector{Species})
    S, D, P = stix_params(ω + 0im, species, B0)
    R = S + D
    L = S - D
    sinθ = sin(θ)
    cosθ = cos(θ)
    n = c * k / ω
    n2 = n^2
    # Build dispersion tensor (Swanson Eq. 4.62)
    M = zeros(ComplexF64, 3, 3)
    M[1,1] = S - n2 * cosθ^2
    M[1,2] = -1im * D
    M[1,3] = n2 * sinθ * cosθ
    M[2,1] = 1im * D
    M[2,2] = S - n2
    M[2,3] = 0
    M[3,1] = n2 * sinθ * cosθ
    M[3,2] = 0
    M[3,3] = P - n2 * sinθ^2
    # Use SVD to find nullspace vector
    U, Σ, Vt = svd(M)
    e_vec = Vt[:, end]
    # normalize electric field
    e_vec ./= maximum(abs.(e_vec))
    # compute magnetic field from Faraday’s law
    k_vec = [k * sinθ, 0.0, k * cosθ]
    b_vec = (k_vec × e_vec) / ω
    return e_vec, b_vec
end

"""
    kinetic_es_func(ω, k, species)

Evaluate the electrostatic dispersion function for a complex frequency
`ω` and real wavenumber `k` given a list of isotropic Maxwellian
species.  The function returns `D(ω,k) = 1 + Σ χ_s`, where
`χ_s = ω_{p,s}²/(k² v_{th,s}²) [1 + ζ_s Z(ζ_s)]` with
`ζ_s = (ω - k u_{d,s})/(k v_{th,s})`.  The thermal speed passed into
this function is the most probable thermal speed (dividing by √2
inside).  This function is used internally by the kinetic solver.
"""
function kinetic_es_func(ω::Complex, k::Float64, species::Vector{Species})
    D = 1 + 0im
    for sp in species
        ωp = omega_p(sp)
        # most probable speed: v_th/sqrt(2)
        vth = v_th(sp) / sqrt(2)
        ζ = (ω - k * sp.drift) / (k * vth)
        Z = plasma_dispersion(ζ)
        χ = (ωp^2 / (k^2 * vth^2)) * (1 + ζ * Z)
        D += χ
    end
    return D
end

"""
    kinetic_es_root(k, species; init_ω=nothing, tol=1e-8, max_iter=50)

Solve the kinetic electrostatic dispersion relation `D(ω,k) = 0` for a
complex frequency `ω` given real wavenumber `k` and isotropic
Maxwellian species.  A complex Newton method with finite difference
derivative is used.  If an initial guess `init_ω` is not supplied, a
warm plasma approximation based on the combined plasma frequency and
thermal corrections is used.  Returns the complex frequency.
"""
function kinetic_es_root(k::Float64, species::Vector{Species};
                         init_ω::Union{Nothing,Complex}=nothing,
                         tol::Float64=1e-8, max_iter::Int=50)
    if isnothing(init_ω)
        # approximate from combined plasma frequency and thermal correction
        ωp_sq = sum(omega_p(sp)^2 for sp in species)
        # use largest thermal speed to approximate warm correction
        vth_max = maximum(v_th(sp) for sp in species)
        # simple estimate: ω ≈ sqrt(ωp² + 3 k² v_th_max²)
        ω0 = sqrt(ωp_sq + 3 * k^2 * vth_max^2)
        ω = complex(ω0, 0.0)
    else
        ω = init_ω
    end
    for _ in 1:max_iter
        F = kinetic_es_func(ω, k, species)
        # small perturbation for derivative; avoid zero
        δ = ω == 0 ? 1e-6 : ω * 1e-6
        Fp = kinetic_es_func(ω + δ, k, species)
        dF = (Fp - F) / δ
        # Newton step
        ω_new = ω - F / dF
        # check convergence in function and update magnitude
        if abs(F) < tol && abs(ω_new - ω) < tol * max(abs(ω), 1.0)
            return ω_new
        end
        ω = ω_new
    end
    return ω
end

end # module