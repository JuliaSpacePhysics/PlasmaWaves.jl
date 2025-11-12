module PlasmaWavesDimensionalDataExt

using DimensionalData
using PlasmaWaves
using PlasmaWaves: wavpol_svd
import PlasmaWaves: twavpol, twavpol_svd
using DimensionalData.Dimensions: @dim, Dimension, TimeDim

function timedim(x, query = nothing)
    query = something(query, TimeDim)
    qdim = dims(x, query)
    return isnothing(qdim) ? dims(x, 1) : qdim
end

times(x::AbstractDimArray, args...) = lookup(timedim(x, args...))

abstract type FrequencyDim{T} <: Dimension{T} end
@dim 𝑓 FrequencyDim "Frequency"

_ellipticity(res, dims) = DimArray(res.ellipticity, dims; name = "Ellipticity")
_power(res, dims) = DimArray(res.power, dims; name = "Power")
_waveangle(res, dims) = DimArray(res.waveangle, dims; name = "Wave normal angle")

"""
    twavpol(x; nfft = 256, noverlap = div(nfft, 2))

A convenience wrapper around [`wavpol`](@ref) that works with DimensionalData arrays.

It automatically extracts the time dimension and returns the results as a DimStack with properly labeled dimensions.
"""
function PlasmaWaves.twavpol(x::AbstractDimArray; fs = nothing, nfft = 256, noverlap = div(nfft, 2), kwargs...)
    t = times(x)
    fs = @something fs samplingrate(t)
    res = wavpol(parent(x), fs; nfft, noverlap, kwargs...)
    dims = (Ti(t[res.indices]), 𝑓(res.freqs))
    return DimStack(
        (
            power = _power(res, dims),
            degpol = DimArray(res.degpol, dims; name = "Degree of polarization"),
            waveangle = _waveangle(res, dims),
            ellipticity = _ellipticity(res, dims),
            helicity = DimArray(res.helicity, dims; name = "Helicity"),
        )
    )
end

function PlasmaWaves.twavpol_svd(x::AbstractDimArray; fs = nothing, nfft = 256, noverlap = div(nfft, 2), kwargs...)
    t = times(x)
    fs = @something fs samplingrate(t)
    res = wavpol_svd(parent(x), fs; nfft, noverlap, kwargs...)
    dims = (Ti(t[res.indices]), 𝑓(res.freqs))
    return DimStack(
        (
            power = _power(res, dims),
            planarity = DimArray(res.planarity, dims; name = "Planarity"),
            waveangle = _waveangle(res, dims),
            ellipticity = _ellipticity(res, dims),
        )
    )
end

end
