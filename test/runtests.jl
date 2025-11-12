using Test
using PlasmaWaves

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(PlasmaWaves)
end
