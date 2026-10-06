"""
    PlasmaWaves

Plasma wave analysis: polarization with [`wavpol`](@ref) and [`twavpol`](@ref).
"""
module PlasmaWaves
using FFTW
using LinearAlgebra
using SpaceDataModel: unwrap, times, cadence, SpaceDataModel

using Bumper
using PrecompileTools
export spectral_matrix, wavpol, twavpol
export Means, Samson, Santolik

include("utils.jl")
include("spectral_matrix.jl")
include("polarization.jl")
include("svd.jl")
include("methods.jl")

"""
    twavpol(X; fs, kw...)

[`wavpol`](@ref) on time series `X`, with `fs` inferred from the `times` dimension of `X` if not given.
Adds the window centre `times` to the result.
"""
twavpol(X; kwargs...) = _twavpol(X; kwargs...)

@deprecate twavpol_svd(X; kwargs...) twavpol(X; method = Santolik(), kwargs...)
@deprecate wavpol_svd(X, args...; kwargs...) wavpol(X, args...; method = Santolik(), kwargs...) false

function _twavpol(x; fs = nothing, dim = 1, kwargs...)
    t = unwrap(SpaceDataModel.dim(x, dim))
    fs = @something fs 1 / cadence(Float64, t)
    res = wavpol(x, fs; dim, kwargs...)
    return (; times = t[res.indices], res...)
end

include("workload.jl")

end # module
