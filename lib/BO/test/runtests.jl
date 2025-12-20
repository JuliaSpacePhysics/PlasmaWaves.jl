using BO
using Test
using Aqua

@testset "BO.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BO)
    end
    # Write your tests here.
end
