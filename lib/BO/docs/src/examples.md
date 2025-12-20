# Examples

This page demonstrates how to use the BO.jl package to solve kinetic dispersion relations for plasmas with arbitrary velocity distributions.

## Basic Usage: Bi-Maxwellian Plasma

The simplest case is a bi-Maxwellian distribution (isotropic or anisotropic Maxwellian).

```julia
using BO
using NLsolve  # For root finding

# Physical parameters
B0 = 0.1  # Magnetic field (Tesla)

# Define plasma species
# Protons (ions)
proton = Species(
    1.0,      # charge (units of e)
    1.0,      # mass (units of m_p)
    5e19,     # density (m^-3)
    1986.734, # parallel temperature (eV)
    993.367,  # perpendicular temperature (eV)
)

# Electrons (Maxwellian)
electron = Species(
    -1.0,      # charge
    5.447e-4,  # mass (m_e/m_p)
    5e19,      # density
    496.683,   # parallel temperature (eV)
    496.683,   # perpendicular temperature (eV)
)

species = [proton, electron]

# Wave vector parameters
k = 0.2  # normalized wavenumber
theta = 45 * π / 180  # propagation angle

# Compute normalizations
wn = abs(B0 * BO.Q_E / BO.M_P)  # proton cyclotron frequency
vA = B0 / sqrt(BO.MU_0 * sum(s.m * BO.M_P * s.n for s in species))  # Alfvén speed
kn = wn / vA

kx = k * kn * sin(theta)
kz = k * kn * cos(theta)

# Create solver parameters
params = create_solver_params(species, B0, kx, kz; N=2)

# Solve dispersion relation
# Initial guess for frequency (normalized to wn)
w0 = (0.0 + 0.01im) * wn

# Root finding
function dispersion(w_vec)
    w = w_vec[1] + im * w_vec[2]
    f = fDrHH(complex(w), params)
    [real(f), imag(f)]
end

result = nlsolve(dispersion, [real(w0), imag(w0)], ftol=1e-14)
w_solution = result.zero[1] + im * result.zero[2]

println("ω/ωₙ = ", w_solution / wn)
```

## Kappa Distribution

Kappa distributions are common in space plasmas and exhibit power-law tails. BO.jl supports several kappa-type distributions.

```julia
using BO
using NLsolve

B0 = 0.1

# Ion with product bi-kappa distribution
ion_kappa = Species(
    1.0, 1.0, 5e19, 1986.734, 993.367;
    distribution=:product_bi_kappa,
    kappa=5.5,  # Kappa index
    Nz=16, Nx=16  # Hermite expansion order
)

# Maxwellian electrons
electron = Species(-1.0, 5.447e-4, 5e19, 496.683, 496.683)

species = [ion_kappa, electron]

# Setup wave vector (same as above)
wn = abs(B0 * BO.Q_E / BO.M_P)
vA = B0 / sqrt(BO.MU_0 * sum(s.m * BO.M_P * s.n for s in species))
kn = wn / vA
k, theta = 0.2, 45π/180
kx, kz = k * kn * sin(theta), k * kn * cos(theta)

params = create_solver_params(species, B0, kx, kz)
```

## Ring Beam Electron Instabilities (Umeda 2012)

This example reproduces the ring-beam electron instabilities from Umeda et al. (2012),
demonstrating unstable wave solutions under a ring beam electron distribution at θ = 40°.

The plasma consists of:
- Ring beam electrons (10% density): parallel drift `vdz = 0.1c`, perpendicular ring drift `vdr = 0.05c`
- Background electrons (90% density): Maxwellian

```julia
using BO
using NLsolve

# Physical parameters from Umeda 2012
B0 = 96.24e-9  # Magnetic field (Tesla)
T_eV = 51.0    # Temperature (eV) for both populations
m_e_mp = 5.447e-4  # Electron mass ratio m_e/m_p

# Ring beam electron population (10% density)
ring_beam = Species(
    -1.0,       # charge (electron)
    m_e_mp,     # mass (m_e/m_p)
    1e5,        # density (m^-3)
    T_eV,       # parallel temperature (eV)
    T_eV,       # perpendicular temperature (eV);
    vdz=0.1,    # parallel drift velocity (0.1c)
    vdr=0.05,   # perpendicular ring beam drift (0.05c)
)

# Background Maxwellian electrons (90% density)
background = Species(
    -1.0,
    m_e_mp,
    9e5,        # 90% of total density
    T_eV,
    T_eV,
)

species = [ring_beam, background]

# Compute normalization (electron cyclotron frequency and Debye length)
wce = abs(B0 * BO.Q_E / BO.M_E)  # Electron cyclotron frequency
ms = [s.m * BO.M_P for s in species]
ns = [s.n for s in species]
Ts = [s.Tz * BO.Q_E / BO.K_B for s in species]
lambdaDs = [sqrt(BO.EPS_0 * BO.K_B * Ts[i] / (ns[i] * BO.Q_E^2)) for i in 1:2]
lambdaD = sqrt(1.0 / sum(1.0 ./ lambdaDs.^2))  # Total Debye length
kn = 1 / lambdaD
wn = wce

# Propagation angle θ = 40°
theta = 40 * π / 180

# Scan over normalized wavenumber k*λD
k_range = 0.01:0.0025:0.3

function solve_dispersion(species, B0, k, theta, kn, wn; N=6)
    kx = k * kn * sin(theta)
    kz = k * kn * cos(theta)
    params = create_solver_params(species, B0, kx, kz; N)

    function dispersion(w_vec)
        w = w_vec[1] + im * w_vec[2]
        f = fDrHH(complex(w), params)
        [real(f), imag(f)]
    end

    return dispersion, params
end

# Find unstable modes
results = ComplexF64[]
w_guess = (0.0 + 0.01im) * wn

for k in k_range
    dispersion, params = solve_dispersion(species, B0, k, theta, kn, wn)

    try
        result = nlsolve(dispersion, [real(w_guess), imag(w_guess)], ftol=1e-14)
        w_solution = result.zero[1] + im * result.zero[2]
        push!(results, w_solution / wn)
        w_guess = w_solution
    catch
        push!(results, NaN + im * NaN)
    end
end

# Results: real part is ωᵣ/ωce, imaginary part is growth rate γ/ωce
# Positive imaginary part indicates instability
println("k*λD = ", collect(k_range)[1:5], " ...")
println("ωᵣ/ωce = ", real.(results[1:5]), " ...")
println("γ/ωce = ", imag.(results[1:5]), " ...")
```

The analytical drift bi-Maxwellian ring beam distribution is given by:

```math
f(v_\\parallel, v_\\perp) = \\frac{1}{\\pi^{3/2} v_{t\\parallel} v_{t\\perp}^2 A_s}
\\exp\\left(-\\frac{(v_\\parallel - v_{d\\parallel})^2}{v_{t\\parallel}^2} - \\frac{(v_\\perp - v_{d\\perp})^2}{v_{t\\perp}^2}\\right)
```

where ``A_s = e^{-d_s^2} + \\sqrt{\\pi} d_s \\mathrm{erfc}(-d_s)`` is the normalization factor
with ``d_s = v_{d\\perp}/v_{t\\perp}``.

### Ring Beam Distribution Visualization

```julia
using BO

# Generate the ring beam distribution on a velocity grid
T_eV = 51.0
vdz_c = 0.1   # parallel drift (c)
vdr_c = 0.05  # perpendicular drift (c)
m_e_mp = 5.447e-4

# Generate distribution using gen_fv2d
fvdat = gen_fv2d(T_eV, T_eV, vdz_c, vdr_c, m_e_mp; distribution=:bi_maxwellian)

# fvdat contains: vz, vx (velocity grids), fv (distribution values)
# Plot with contourf(fvdat.vz, fvdat.vx, fvdat.fv)
```

## Shell Distribution

Shell distributions, where particles are concentrated on a spherical shell in velocity space, are relevant for solar wind and magnetospheric physics.

```julia
using BO

# Shell distribution
shell_species = Species(
    1.0, 1.0, 1e19, 100.0, 100.0;
    vdr=0.01,  # shell radius parameter
    distribution=:shell,
    Nz=16, Nx=16
)
```

## Computing Dispersion Curves

To compute dispersion curves over a range of wavenumbers:

```julia
using BO
using NLsolve

function compute_dispersion_curve(species, B0, k_range, theta; N=2)
    wn = abs(B0 * BO.Q_E / BO.M_P)
    vA = B0 / sqrt(BO.MU_0 * sum(s.m * BO.M_P * s.n for s in species))
    kn = wn / vA

    results = ComplexF64[]
    w_guess = (0.0 + 0.01im) * wn

    for k in k_range
        kx = k * kn * sin(theta)
        kz = k * kn * cos(theta)
        params = create_solver_params(species, B0, kx, kz; N)

        function dispersion(w_vec)
            w = w_vec[1] + im * w_vec[2]
            f = fDrHH(complex(w), params)
            [real(f), imag(f)]
        end

        try
            result = nlsolve(dispersion, [real(w_guess), imag(w_guess)], ftol=1e-14)
            w_solution = result.zero[1] + im * result.zero[2]
            push!(results, w_solution / wn)
            w_guess = w_solution  # Use previous solution as next guess
        catch
            push!(results, NaN + im * NaN)
        end
    end

    return results
end

# Example usage
B0 = 0.1
proton = Species(1.0, 1.0, 5e19, 1986.734, 993.367)
electron = Species(-1.0, 5.447e-4, 5e19, 496.683, 496.683)

k_range = 0.1:0.025:0.4
theta = 45π/180

omega = compute_dispersion_curve([proton, electron], B0, k_range, theta)

# Plot results
# using Plots
# plot(k_range, real.(omega), label="ωᵣ/ωₙ", xlabel="k/kₙ", ylabel="ω/ωₙ")
# plot!(k_range, imag.(omega), label="ωᵢ/ωₙ")
```

## Matrix Eigenvalue Solver

The matrix eigenvalue method finds all wave modes simultaneously by transforming the dispersion relation
into a matrix eigenvalue problem using J-pole approximation for the plasma dispersion function.

This approach is more efficient when you need to find multiple modes at once, and doesn't require
initial guesses for the root finder.

```julia
using BO

# Umeda 2012 ring beam configuration
B0 = 96.24e-9  # Tesla

# Ring beam electrons (10% density)
ring_beam = Species(-1.0, 5.447e-4, 1e5, 51.0, 51.0; vdz=0.1, vdr=0.05)
# Background electrons (90% density)
background = Species(-1.0, 5.447e-4, 9e5, 51.0, 51.0)

species = [ring_beam, background]

# Compute normalization
wce = abs(B0 * BO.Q_E / BO.M_E)

# Calculate Debye length for normalization
Ts = [s.Tz * BO.Q_E / BO.K_B for s in species]
ns = [s.n for s in species]
lambdaDs = [sqrt(BO.EPS_0 * BO.K_B * Ts[i] / (ns[i] * BO.Q_E^2)) for i in eachindex(species)]
lambdaD = sqrt(1.0 / sum(1.0 ./ lambdaDs.^2))

# Wave vector: k*λD = 0.03, θ = 40°
k = 0.03 / lambdaD
theta = 40 * π / 180
kx = k * sin(theta)
kz = k * cos(theta)

# Create solver parameters
params = create_solver_params(species, B0, kx, kz; N=6)

# Solve using matrix eigenvalue method
# J=12 provides good accuracy (J-pole approximation order)
eigenvalues = solve_dispersion_matrix(params, kx, kz; J=12)

# Filter for unstable modes (positive growth rate)
unstable = filter(e -> isfinite(e) && imag(e) > 0.001*wce, eigenvalues)

# Display results
for e in sort(unstable, by=e->imag(e), rev=true)[1:min(5, end)]
    println("ω/ωce = ", real(e)/wce, " + ", imag(e)/wce, "i")
end
# Expected: ω/ωce ≈ 0.62 + 0.16i (unstable mode)
```

### J-pole Coefficients

The J-pole approximation expresses the plasma dispersion function as:

```math
Z(\zeta) \approx \sum_{j=1}^{J} \frac{b_j}{\zeta - c_j}
```

Available J values: 4, 6, 8, 10, 12, 16, 20, 24, 28, 32. Larger J provides higher accuracy
but increases matrix size and computation time.

```julia
using BO

# Get J-pole coefficients
jpole = get_jpole_coefficients(12)
println("J = ", jpole.J)
println("b coefficients: ", jpole.bzj[1:3], " ...")
println("c coefficients: ", jpole.czj[1:3], " ...")
```

### Dispersion Curve Scan with Matrix Method

```julia
using BO

# Setup species (Umeda 2012)
B0 = 96.24e-9
ring_beam = Species(-1.0, 5.447e-4, 1e5, 51.0, 51.0; vdz=0.1, vdr=0.05)
background = Species(-1.0, 5.447e-4, 9e5, 51.0, 51.0)
species = [ring_beam, background]

# Normalizations
wce = abs(B0 * BO.Q_E / BO.M_E)
Ts = [s.Tz * BO.Q_E / BO.K_B for s in species]
ns = [s.n for s in species]
lambdaDs = [sqrt(BO.EPS_0 * BO.K_B * Ts[i] / (ns[i] * BO.Q_E^2)) for i in eachindex(species)]
lambdaD = sqrt(1.0 / sum(1.0 ./ lambdaDs.^2))

# Scan k*λD from 0.01 to 0.3
theta = 40 * π / 180
k_range = 0.01:0.005:0.3

results = []
for k_norm in k_range
    k = k_norm / lambdaD
    kx = k * sin(theta)
    kz = k * cos(theta)

    params = create_solver_params(species, B0, kx, kz; N=6)
    eigenvalues = solve_dispersion_matrix(params, kx, kz; J=12)

    # Find most unstable mode
    unstable = filter(e -> isfinite(e) && imag(e) > 0, eigenvalues)
    if !isempty(unstable)
        w = sort(unstable, by=e->imag(e), rev=true)[1]
        push!(results, (k=k_norm, wr=real(w)/wce, wi=imag(w)/wce))
    end
end

# Plot: real(ω)/ωce vs k*λD, imag(ω)/ωce vs k*λD
# using Plots
# plot([r.k for r in results], [r.wr for r in results], label="ωᵣ/ωce")
# plot!([r.k for r in results], [r.wi for r in results], label="γ/ωce")
```

## Direct Parameter Construction

For advanced usage, you can construct `HHSolverParams` directly:

```julia
using BO

# Directly specify all parameters
params = HHSolverParams{Float64}(
    2,                    # S: number of species
    BO.C_LIGHT^2,         # c2: speed of light squared
    [1e8, -1e11],         # wcs: cyclotron frequencies
    [1e16, 1e20],         # wps2: plasma frequencies squared
    [1e-3, 1e-6],         # rhocs: cyclotron radii
    1e3,                  # kx: perpendicular wavenumber
    1e3,                  # kz: parallel wavenumber
    [1e6, 1e7],           # vtzs: parallel thermal velocities
    [1e6, 1e7],           # vtps: perpendicular thermal velocities
    [0.0, 0.0],           # vdsz: parallel drift velocities
    [0.0, 0.0],           # ds: ring beam parameters
    [1.0, 1.0],           # As: normalization parameters
    [2, 2],               # Nss: maximum harmonic numbers
    [ones(1,1), ones(1,1)], # aslm: Hermite coefficients
    [0, 0],               # msmax: max perpendicular Hermite index
    [0, 0],               # lsmax: max parallel Hermite index
)

# Evaluate dispersion relation
w = 1e8 + 1e6im
D = fDrHH(w, params)
```

## Fluid Dispersion Solver

For cases where kinetic effects are not important, BO.jl provides a multi-fluid electromagnetic
dispersion solver. This is faster than the kinetic solver and useful for MHD-scale phenomena.

```julia
using BO

# Define fluid species
B0 = 1e-4  # Magnetic field (Tesla)

# Protons (isothermal)
proton = FluidSpecies(
    1.0,      # charge (units of e)
    1.0,      # mass (units of m_p)
    1e18,     # density (m^-3)
    100.0,    # parallel temperature (eV)
    100.0,    # perpendicular temperature (eV)
)

# Electrons
electron = FluidSpecies(
    -1.0,
    5.447e-4,  # m_e/m_p
    1e18,
    100.0,
    100.0,
)

species = [proton, electron]
params = create_fluid_params(species, B0)

# Wave vector
kx = 1e-3  # m^-1
kz = 1e-2  # m^-1

# Solve fluid dispersion relation
eigenvalues = solve_fluid_dispersion(params, kx, kz)

# Normalize by ion cyclotron frequency
wci = abs(params.wcs[1])
println("Eigenvalues (ω/ωci):")
for e in sort(eigenvalues, by=e->abs(real(e)))
    println("  ", real(e)/wci, " + ", imag(e)/wci, "i")
end
```

### Fluid Closure Models

The fluid solver supports different polytrope closure models through the `gamma_z` and `gamma_p`
parameters:

| Model | gamma_z | gamma_p | Description |
|-------|---------|---------|-------------|
| Isothermal | 1.0 | 1.0 | Constant temperature (default) |
| Adiabatic | 5/3 | 5/3 | Adiabatic compression |
| CGL | 3.0 | 2.0 | Chew-Goldberger-Low double adiabatic |

```julia
using BO

# CGL closure for anisotropic plasma
proton_cgl = FluidSpecies(
    1.0, 1.0, 1e18, 100.0, 50.0;  # Tp < Tz (temperature anisotropy)
    gamma_z=3.0, gamma_p=2.0      # CGL closure
)
```

### Including Drift

Parallel drift velocities can be specified for beam-plasma systems:

```julia
using BO

# Drifting proton beam
beam = FluidSpecies(
    1.0, 1.0, 1e17, 100.0, 100.0;
    vdz=0.01  # 1% of c parallel drift
)

# Background plasma
background = FluidSpecies(1.0, 1.0, 9e17, 100.0, 100.0)
electron = FluidSpecies(-1.0, 5.447e-4, 1e18, 100.0, 100.0)

species = [beam, background, electron]
```

## Supported Distributions

| Distribution | Symbol | Description |
|--------------|--------|-------------|
| Bi-Maxwellian | `:maxwellian`, `:bi_maxwellian` | Standard Maxwellian with temperature anisotropy |
| Bi-kappa | `:bi_kappa` | Kappa distribution in both directions |
| Product bi-kappa | `:product_bi_kappa` | Product of parallel and perpendicular kappa |
| Kappa-Maxwellian | `:kappa_maxwellian` | Kappa parallel, Maxwellian perpendicular |
| Shell | `:shell` | Particles on spherical shell in velocity space |
| Slowing down | `:slowing_down` | Slowing down distribution for energetic particles |

## References

- Xie, H. S. (2025). Efficient Framework for Solving Plasma Waves with Arbitrary Distributions. [arXiv:2501.06477](https://arxiv.org/abs/2501.06477)
- Xie, H. S. (2025). Phys. Plasmas 32, 060702. [doi:10.1063/5.0275307](https://doi.org/10.1063/5.0275307)
- Xie, H. S. (2014). PDRF: A general dispersion relation solver for magnetized multi-fluid plasma. Comput. Phys. Comm. 185, 670-675. [doi:10.1016/j.cpc.2013.10.012](https://doi.org/10.1016/j.cpc.2013.10.012)
