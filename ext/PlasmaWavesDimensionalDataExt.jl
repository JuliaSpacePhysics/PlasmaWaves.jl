module PlasmaWavesDimensionalDataExt

using DimensionalData
using PlasmaWaves
using PlasmaWaves: _twavpol
using DimensionalData.Dimensions: @dim, Dimension, TimeDim, hasdim, dimnum

abstract type FrequencyDim{T} <: Dimension{T} end
@dim 𝑓 FrequencyDim "Frequency"

# A no-error version of `dimnum`
_dimnum(x, dim) = hasdim(x, dim) ? dimnum(x, dim) : nothing

const LABELS = (;
    power = "Power", degpol = "Degree of polarization", waveangle = "Wave normal angle",
    ellipticity = "Ellipticity", ellipticity_perp = "Ellipticity ⊥ B", planarity = "Planarity",
)

function PlasmaWaves.twavpol(x::AbstractDimArray; dim = nothing, kwargs...)
    dim = @something dim _dimnum(x, TimeDim) _dimnum(x, Dim{:time}) 1
    res = _twavpol(x; dim, kwargs...)
    dims = dim == 1 ? (Ti(res.times), 𝑓(res.freqs)) : (𝑓(res.freqs), Ti(res.times))
    fields = Base.structdiff(res, NamedTuple{(:times, :indices, :freqs)})
    layers = map(keys(fields), values(fields)) do k, v
        DimArray(v, dims; name = string(get(LABELS, k, k)))
    end
    return DimStack(NamedTuple{keys(fields)}(layers))
end

end
