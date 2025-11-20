module PlasmaWaves
using FFTW
using LinearAlgebra
using StaticArrays
using SpaceDataModel: times, cadence

using Tullio, Bumper
export spectral_matrix, wavpol, twavpol, twavpol_svd, wpol_helicity, polarization

include("spectral_matrix.jl")
include("polarization.jl")
include("helicty.jl")
include("Stokes.jl")

"""
    twavpol(X; fs = nothing, nfft = 256, noverlap = div(nfft, 2))

Polarization analysis of time series data `X` (each column is a component) of sampling frequency `fs`.

If `fs` is not provided, it will be inferred from the `times` dimension of `X`.

See [`wavpol`](@ref) for details.
"""
twavpol(X; kwargs...) = _twavpol(wavpol, X; kwargs...)

"""
    twavpol_svd(X; fs = nothing, nfft = 256, noverlap = div(nfft, 2))

Polarization analysis of time series data `X` (each column is a component) of sampling frequency `fs`, using singular value decomposition (SVD) method.
"""
twavpol_svd(x; kwargs...) = _twavpol(wavpol_svd, x; kwargs...)

# Internal function for dispatch
function _twavpol(f, x; fs = nothing, nfft = 256, noverlap = div(nfft, 2), kwargs...)
    t = times(x)
    fs = @something fs 1 / cadence(Float64, t)
    res = f(x, fs; nfft, noverlap, kwargs...)
    return (times = t[res.indices], res...)
end


end # module
