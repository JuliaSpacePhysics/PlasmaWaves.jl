module PlasmaWavesDimensionalDataExt

using DimensionalData
using PlasmaWaves
using PlasmaWaves: wavpol, wavpol_svd, _twavpol
import PlasmaWaves: twavpol, twavpol_svd
using DimensionalData.Dimensions: @dim, Dimension, TimeDim, hasdim, dimnum

abstract type FrequencyDim{T} <: Dimension{T} end
@dim 𝑓 FrequencyDim "Frequency"

# A no-error version of `dimnum`
_dimnum(x, dim) = hasdim(x, dim) ? dimnum(x, dim) : nothing

_ellipticity(res, dims) = DimArray(res.ellipticity, dims; name = "Ellipticity")
_power(res, dims) = DimArray(res.power, dims; name = "Power")
_waveangle(res, dims) = DimArray(res.waveangle, dims; name = "Wave normal angle")

function PlasmaWaves.twavpol(x::AbstractDimArray; dim = nothing, kwargs...)
    dim = @something dim _dimnum(x, TimeDim) _dimnum(x, Dim{:time}) 1
    res = _twavpol(wavpol, x; dim = dim, kwargs...)
    dims = dim == 1 ? (Ti(res.times), 𝑓(res.freqs)) : (𝑓(res.freqs), Ti(res.times))
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

function PlasmaWaves.twavpol_svd(x::AbstractDimArray; dim = nothing, kwargs...)
    dim = @something dim _dimnum(x, TimeDim) _dimnum(x, Dim{:time}) 1
    res = _twavpol(wavpol_svd, x; dim = dim, kwargs...)
    dims = dim == 1 ? (Ti(res.times), 𝑓(res.freqs)) : (𝑓(res.freqs), Ti(res.times))
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
