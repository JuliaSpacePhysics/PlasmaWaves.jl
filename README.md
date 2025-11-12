# PlasmaWaves

[![Build Status](https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/PlasmaWaves.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/PlasmaWaves.jl)

Wave analysis for (space) plasmas, built in Julia.

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("PlasmaWaves")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/PlasmaWaves.jl/dev/)

## Features and Roadmap

- [x] Wave polarization analysis with degree of polarization, wave normal angle, helicity, ellipticity, and planarity metrics
- [ ] Wave propagation analysis
  - [x] SVD of the magnetic spectral matrix
  - [ ] Electromagnetic SVD
- [ ] Wave dispersion relation analysis

## Quick start

```julia
using PlasmaWaves

# X is an N×3 array of field-aligned magnetic fluctuations
res = wavpol(X, fs = 128.0; nfft = 256)

res.degpol    # degree of polarization across time-frequency bins
res.waveangle # wave normal angle estimates
```

For SVD-derived planarity metrics, call `wavpol_svd` or `twavpol_svd`.

⚠️ **Development Status**: This package is in active development. While functional, the functionality is not well tested and the API may undergo changes in future releases. Please test thoroughly before using in scientific work.
