# Changelog

## [0.2.0]

### Changed

- **Breaking**: the polarization method is the `method` keyword of `wavpol` and `twavpol`: `Means()` (default, the SPEDAS row method), `Samson()`, `Santolik()`, or any callable mapping the spectral matrix, scaled to unit trace, to a `NamedTuple`. `twavpol_svd(X)` and `wavpol_svd(X)` are deprecated for `method = Santolik()` ([#26]).
- **Breaking**: a result field means the same quantity for every method ([#26]):
  - `ellipticity` is the signed axial ratio of the polarization ellipse in its own plane, as `Santolik` already returned. It replaces `helicity` (unsigned) of the default method.
  - `ellipticity_perp` is the axial ratio of the ellipse projected onto the plane ⊥ z, previously `ellipticity` of the default method (PySPEDAS `elliptict`).
- **Breaking**: `wpol_helicity` is removed; `Means()(S)` returns the same quantities ([#26]).
- **Breaking**: `polarization` is renamed `degree_of_polarization` and no longer exported. The Stokes form `degree_of_polarization(S0, S1, S2, S3)` returns p², as the spectral-matrix form does; `polarization(S0, S1, S2, S3)` returned p ([#26]).
- **Breaking**: `smooth_f` must have odd length; an even window placed results half a bin off ([#25]).

### Added

- `Samson()`: ellipticity and wave normal from the principal eigenvector of the spectral matrix, less biased than `Means()` for weakly polarized waves ([#26]).
- `Santolik()` also returns `degpol` ([#26]).

### Fixed

- `noverlap` is the number of samples shared by consecutive windows, as documented; it was used as the step ([#25]).
- `power` is a one-sided PSD that no longer doubles the DC and Nyquist bins. Polarization parameters are NaN at the edge bins where the frequency smoothing does not fit, instead of values from an unsmoothed, rank-1 spectral matrix ([#25]).
- Ellipticity is no longer zero at exactly parallel propagation, nor NaN when a component has no power. Rows of the spectral matrix are weighted by power, so a noise-only component no longer biases ellipticity toward 2/3 of its value ([#25]).
- `Santolik()` on Float32 input forms its Gram matrix in Float64; planarity of a planar wave came out as low as 0.9 ([#25]).
- Ellipticity is no longer zero at exactly perpendicular propagation, nor 0 or NaN for Float32 input of small amplitude (about 1e-6), where powers of the spectral matrix underflowed ([#26]).

[#25]: https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/pull/25
[#26]: https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/pull/26
