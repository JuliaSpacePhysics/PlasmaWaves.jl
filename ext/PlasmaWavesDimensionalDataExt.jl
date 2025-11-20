module PlasmaWavesDimensionalDataExt

using DimensionalData
using PlasmaWaves
using PlasmaWaves: wavpol, wavpol_svd, samplingrate, _twavpol
import PlasmaWaves: twavpol, twavpol_svd
using DimensionalData.Dimensions: @dim, Dimension

abstract type FrequencyDim{T} <: Dimension{T} end
@dim 𝑓 FrequencyDim "Frequency"

_ellipticity(res, dims) = DimArray(res.ellipticity, dims; name = "Ellipticity")
_power(res, dims) = DimArray(res.power, dims; name = "Power")
_waveangle(res, dims) = DimArray(res.waveangle, dims; name = "Wave normal angle")

function PlasmaWaves.twavpol(x::AbstractDimArray; kwargs...)
    res = _twavpol(wavpol, x; kwargs...)
    dims = (Ti(res.times), 𝑓(res.freqs))
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

function PlasmaWaves.twavpol_svd(x::AbstractDimArray; kwargs...)
    res = _twavpol(wavpol_svd, x; kwargs...)
    dims = (Ti(res.times), 𝑓(res.freqs))
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
