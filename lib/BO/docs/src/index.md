```@meta
CurrentModule = BO
```

# BO.jl

Efficient framework for solving plasma wave kinetic dispersion relations with arbitrary velocity distributions.

## Overview

BO.jl is a Julia implementation of the BO-Arbitrary framework for solving electromagnetic kinetic dispersion relations in magnetized plasmas. It extends the standard Maxwellian-based solvers to support arbitrary velocity distribution functions through Hermite-Hermite (HH) expansion.

### Features

- Solve kinetic electromagnetic dispersion relations for hot magnetized plasmas
- Support for arbitrary velocity distributions via Hermite expansion:
  - Bi-Maxwellian (ring beam)
  - Bi-kappa distributions
  - Shell distributions
  - Slowing down distributions
  - Custom distributions via grid input
- Multiple solver methods:
  - Root finding: `fDrHH` for targeted kinetic mode search
  - Kinetic matrix eigenvalue: `solve_dispersion_matrix` for finding all kinetic modes
  - Fluid matrix eigenvalue: `solve_fluid_dispersion` for multi-fluid MHD modes
- Efficient FFT-based plasma dispersion function computation
- J-pole approximation for efficient matrix eigenvalue computation

## Quick Start

```julia
using BO
using NLsolve

# Define plasma species
proton = Species(1.0, 1.0, 5e19, 1000.0, 500.0)  # q, m, n, Tz, Tp
electron = Species(-1.0, 5.447e-4, 5e19, 500.0, 500.0)

# Setup parameters
B0 = 0.1  # Tesla
kx, kz = 1e3, 1e3  # m^-1
params = create_solver_params([proton, electron], B0, kx, kz)

# Method 1: Root finding for specific mode
w0 = 1e6 + 1e4im
result = nlsolve(w -> let f = fDrHH(w[1]+im*w[2], params); [real(f), imag(f)] end,
                 [real(w0), imag(w0)])

# Method 2: Matrix eigenvalue for all kinetic modes
eigenvalues = solve_dispersion_matrix(params, kx, kz; J=8)

# Method 3: Fluid solver for MHD modes
fluid_proton = FluidSpecies(1.0, 1.0, 5e19, 1000.0, 500.0)
fluid_electron = FluidSpecies(-1.0, 5.447e-4, 5e19, 500.0, 500.0)
fluid_params = create_fluid_params([fluid_proton, fluid_electron], B0)
fluid_eigenvalues = solve_fluid_dispersion(fluid_params, kx, kz)
```

See the [Examples](@ref) page for more detailed usage.

## References

- Xie, H. S. (2025). Efficient Framework for Solving Plasma Waves with Arbitrary Distributions. [arXiv:2501.06477](https://arxiv.org/abs/2501.06477)
- Xie, H. S. (2025). Phys. Plasmas 32, 060702. [doi:10.1063/5.0275307](https://doi.org/10.1063/5.0275307)
- Xie, H. S. (2014). PDRF: A general dispersion relation solver for magnetized multi-fluid plasma. Comput. Phys. Comm. 185, 670-675.

## API Reference

```@index
```

```@autodocs
Modules = [BO]
```
