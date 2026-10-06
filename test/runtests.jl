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

    # Planar wave (S has no component along its normal n̂) in Float32
    planarity = map(1:1000) do _
        n̂ = normalize(randn(Float32, 3))
        P = I - n̂ * n̂'
        u, v = P * randn(ComplexF32, 3), P * randn(ComplexF32, 3)
        PlasmaWaves.svd_polarization(u * u' + v * v').planarity
    end
    @test minimum(planarity) > 0.99
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

@testset "elliptical wave" begin
    # Ellipse with axes 1 and b, rotating right-handed about k̂, which is tilted by θ from z
    b, nfft = 0.6, 256
    t = 0:2047
    ω = 2π * 20 / nfft # centred on bin 21
    for θ in (0.0, π / 6, π / 3)
        e1, e2 = [cos(θ), 0, -sin(θ)], [0, 1, 0]
        X = cos.(ω .* t) * e1' .+ b .* sin.(ω .* t) * e2'
        r = wavpol(X; nfft)
        @test all(≈(1), r.degpol[:, 21])
        @test all(x -> isapprox(x, θ; atol = 1.0e-8), r.waveangle[:, 21])
        @test all(≈(b), r.helicity[:, 21])
        # (x, y) projection has axes cos θ and b; which is major flips at θ = π/3
        @test all(≈(min(b / cos(θ), cos(θ) / b)), r.ellipticity[:, 21])
        r = PlasmaWaves.wavpol_svd(X; nfft)
        @test all(x -> isapprox(x, 1; atol = 1.0e-3), r.planarity[:, 21]) # Gram matrix: ~eps^¼
        @test all(x -> isapprox(x, θ; atol = 1.0e-6), r.waveangle[:, 21])
        @test all(≈(b), r.ellipticity[:, 21])
    end
    # Expected spectral matrix of a parallel wave over white noise: the z row holds only noise
    v = [1, -b * im, 0]
    h, e = PlasmaWaves.wpol_helicity(1.0e4 * v * v' + I)
    @test h ≈ b rtol = 1.0e-3
    @test e ≈ b rtol = 1.0e-3
end

@testset "wavpol windows" begin
    X = randn(1024, 3)
    @test wavpol(X; nfft = 256, noverlap = 192).indices == 129:64:897
    @test wavpol(X; nfft = 256, noverlap = 0).indices == 129:256:897
end

@testset "wavpol power obeys Parseval" begin
    nfft, fs = 64, 4.0
    X = randn(nfft, 3)
    w = PlasmaWaves._smooth_t(nfft)
    res = wavpol(X, fs; nfft, smooth_f = (1,))
    @test sum(res.power) * fs / nfft ≈ sum(abs2, w .* X) / sum(abs2, w)
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
