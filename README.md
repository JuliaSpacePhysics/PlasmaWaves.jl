# PlasmaWaves

[![DOI](https://zenodo.org/badge/1094801941.svg)](https://doi.org/10.5281/zenodo.17657870)

Wave analysis for (space) plasmas.

## Quick start

```julia
using Pkg; Pkg.add("PlasmaWaves")
using PlasmaWaves

# X is an N×3 array of field-aligned magnetic fluctuations
res = wavpol(X, 128.0; nfft = 256)

res.degpol    # degree of polarization across time-frequency bins
res.waveangle # wave normal angle estimates
```

Pass `method = Samson()` (principal eigenvector) or `method = Santolik()` (SVD, with planarity) for other estimators.

## Features and

- [x] Wave polarization analysis with degree of polarization, wave normal angle, ellipticity, and planarity metrics
- [x] Wave propagation analysis using SVD of magnetic field

See [Wave polarization cross-validation with PySPEDAS](https://juliaspacephysics.github.io/PlasmaWaves.jl/dev/pyspedas/) for comparison and benchmarking against PySPEDAS implementation.

## Elsewhere

- [PlasmaBO.jl](https://github.com/JuliaSpacePhysics/PlasmaBO.jl) for wave dispersion relation analysis.
