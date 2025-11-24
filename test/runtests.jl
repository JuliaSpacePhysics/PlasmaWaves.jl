using Test
using PlasmaWaves
using LinearAlgebra
using JET: @test_call

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(PlasmaWaves)
end

@testset "JET" begin
    @test_nowarn PlasmaWaves.workload()
    @test_call PlasmaWaves.workload()
end

@testset "svd of spectral matrix" begin
    a = 3.5
    b = 1.2
    Xf = zeros(ComplexF64, 1, 3)
    Xf[1, 1] = a
    Xf[1, 2] = 1im * b
    S = spectral_matrix(Xf)
    Sf = @view S[:, :, 1]

    @test Sf ≈ [a^2 -im * a * b 0; im * a * b b^2 0; 0 0 0]
    A = [
        a^2     0     0;
        0       b^2   0;
        0       0     0;
        0       a * b 0;
        -a * b  0     0;
        0       0     0
    ]
    U, S, V = PlasmaWaves.svd_S(Sf)
    @test U * Diagonal(S) * V' ≈ A
    @test S ≈ [a * sqrt(a^2 + b^2), b * sqrt(a^2 + b^2), 0]
    s = Base.sign(U[1, 1])
    @test s * U[1, 1] ≈ a / sqrt(a^2 + b^2) + 0im
    @test s * U[2, 2] ≈ b / sqrt(a^2 + b^2) + 0im
    @test s * U[4, 2] ≈ a / sqrt(a^2 + b^2) + 0im
end

@testset "Spectral matrix from time sequence" begin
    # Note: the sign of the off-diagonal elements is opposite to that in the reference due to the convention of DFT they use compared to FFTW
    # $\omega =\frac{2 π k}{N τ}$
    a = 2.8
    b = 0.9
    f = 5.0
    K = 4          # controls sampling cadence (tau = 1 / (4 K f))
    Q = 3          # number of full wave periods observed
    τ = 1 / (4 * K * f)
    M = 4 * K * Q
    t = collect(0:(M - 1)) .* τ

    B1 = @. a * cos(2π * f * t)
    B2 = @. b * sin(2π * f * t)
    B3 = zeros(length(t))
    X = hcat(B1, B2, B3)

    S = spectral_matrix(X)
    freq_bin = Q + 1            # bin aligned with frequency f : k = f * N * τ = Q
    scale = (M / 2)^2
    Sf = @view S[:, :, freq_bin]
    Sf ./= scale
    @test Sf ≈ [
        a^2             im * a * b  0;
        -im * a * b     b^2         0;
        0               0           0
    ]
    @test spectral_matrix(X) ≈ spectral_matrix(X', 2)
end

using Downloads, JLD2

function get_test_data(url)
    filename = splitpath(url)[end]
    localpath = joinpath(@__DIR__, filename)
    isfile(localpath) || Downloads.download(url, localpath)
    return localpath
end

@testset "DimensionalData Integration" begin
    using DimensionalData, UnixTimes
    thc_scf_fac_url = "https://github.com/JuliaSpacePhysics/PlasmaWaves.jl/releases/download/v0.1.1/thc_scf_fac.jld2"
    fpath = get_test_data(thc_scf_fac_url)
    @load fpath thc_scf_fac
    result = twavpol(thc_scf_fac)
    @test result.power == twavpol(thc_scf_fac').power'

    result2 = twavpol_svd(thc_scf_fac)
    @test result2.power == twavpol_svd(thc_scf_fac').power'
end
