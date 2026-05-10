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

@testset "svd_polarization" begin
    # RHCP wave in xy-plane propagating along z: Xf = [a, ib, 0]
    # → S[1,2] = -iab, wave normal = ẑ, planarity = 1, ellipticity = -b/a
    a, b = 3.5, 1.2
    Xf = zeros(ComplexF64, 1, 3)
    Xf[1, 1] = a
    Xf[1, 2] = 1im * b
    S3d = spectral_matrix(Xf)
    Sf = @view S3d[:, :, 1]
    @test Sf ≈ [a^2 -im*a*b 0; im*a*b b^2 0; 0 0 0]
    res = PlasmaWaves.svd_polarization(Sf)
    @test res.theta ≈ 0
    @test res.planarity ≈ 1.0
    @test res.ellipticity ≈ -b / a

    # Wave normal along x: Gram matrix has a11=λ₃=0, a12=a13=0
    S_x = ComplexF64[0 0 0; 0 4.0 2.0; 0 2.0 8.0]
    res_x = PlasmaWaves.svd_polarization(S_x)
    @test !any(isnan, (res_x.theta, res_x.phi, res_x.planarity, res_x.ellipticity))
    @test res_x.theta ≈ π / 2
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
