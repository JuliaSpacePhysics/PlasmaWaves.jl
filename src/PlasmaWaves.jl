module PlasmaWaves
using FFTW
using LinearAlgebra
using StaticArrays

using Tullio, Bumper
export spectral_matrix, wavpol, twavpol, twavpol_svd, wpol_helicity, polarization

include("spectral_matrix.jl")
include("polarization.jl")
include("helicty.jl")
include("Stokes.jl")

function twavpol end
function twavpol_svd end

end # module
