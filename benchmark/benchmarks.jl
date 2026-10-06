using BenchmarkTools
using PlasmaWaves
using Random

const SUITE = BenchmarkGroup()

const rng = MersenneTwister(1)
const Xf = randn(rng, ComplexF64, 129, 3)
const X = randn(rng, 4096, 3)

SUITE["spectral_matrix"] = BenchmarkGroup()
SUITE["spectral_matrix"]["complex"] = @benchmarkable spectral_matrix($Xf)
SUITE["spectral_matrix"]["real"] = @benchmarkable spectral_matrix($X)

SUITE["wavpol"] = BenchmarkGroup()
SUITE["wavpol"]["means"] = @benchmarkable wavpol($X; nfft = 256)
SUITE["wavpol"]["samson"] = @benchmarkable wavpol($X; nfft = 256, method = Samson())
SUITE["wavpol"]["santolik"] = @benchmarkable wavpol($X; nfft = 256, method = Santolik())
