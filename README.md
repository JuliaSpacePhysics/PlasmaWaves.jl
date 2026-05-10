# PlasmaWaves

Wave analysis for (space) plasmas.

## Quick start

```julia
using Pkg; Pkg.add("PlasmaWaves")
using PlasmaWaves

# X is an N×3 array of field-aligned magnetic fluctuations
res = wavpol(X, fs = 128.0; nfft = 256)

res.degpol    # degree of polarization across time-frequency bins
res.waveangle # wave normal angle estimates
```

For SVD-derived planarity metrics, call `wavpol_svd` or `twavpol_svd`.

## Features and Roadmap

- [x] Wave polarization analysis with degree of polarization, wave normal angle, helicity, ellipticity, and planarity metrics
- [ ] Wave propagation analysis
  - [x] SVD of the magnetic spectral matrix
  - [ ] Electromagnetic SVD

See [Wave polarization cross-validation with PySPEDAS](https://juliaspacephysics.github.io/SPEDAS.jl/dev/validation/pyspedas/) for comparison and benchmarking against PySPEDAS implementation.

## Elsewhere

- [PlasmaBO.jl](https://github.com/JuliaSpacePhysics/PlasmaBO.jl) for wave dispersion relation analysis.

## Status

⚠️ **Development Status**: This package is in active development. While functional, the functionality is not fully tested (it has been cross-validated with a Python implementation in `PySPEDAS`). Please test thoroughly for scientific work.

[![DOI](https://zenodo.org/badge/1094801941.svg)](https://doi.org/10.5281/zenodo.17657870)
[![version](https://juliahub.com/docs/General/PlasmaWaves/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlasmaWaves)

[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaWaves.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaWaves.jl)