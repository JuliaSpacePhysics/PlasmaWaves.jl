using Test
using PlasmaWaves

"""
This test suite validates a few basic limits of the `PlasmaWaves`
package.  These tests compare numerical solutions of the dispersion
relations against simple analytic approximations in limiting cases.
They are meant to catch gross errors in the implementation rather
than provide exhaustive coverage.  As with all plasma problems,
appropriate tolerances are used because the Newton solver may
converge to slightly different values depending on the initial
guesses and the model assumptions.
"""

# Physical constants for convenience
const e_charge = 1.602176634e-19
const m_e = 9.10938356e-31
const m_p = 1.67262192369e-27

@testset "Cold electromagnetic plasma" begin
    # Single cold electron species, B=0.  The solution should reduce
    # to the vacuum dispersion corrected by the plasma frequency:
    # ω² = ω_p² + (c k)².  Use k small enough such that frequency is
    # dominated by the plasma term to avoid numerical cancellations.
    n = 1e6
    sp = Species(:e, -e_charge, m_e, n, 0.0, 0.0)
    k = 0.05  # m^-1
    ω_num = cold_em_root(k, 0.0, 0.0, [sp])
    ωp = omega_p(sp)
    ω_pred = sqrt(ωp^2 + (PlasmaWaves.c * k)^2)
    @test isapprox(ω_num, ω_pred; rtol=5e-3)
end

@testset "Langmuir wave with warm correction" begin
    # Warm Langmuir wave (electron only).  At small kλ_D the real part
    # of the frequency satisfies the Bohm–Gross dispersion relation
    # ω² ≈ ω_p² + 3 k² v_th².  Imaginary part is negative (Landau
    # damping).  Use a modest k so that approximation holds.
    n = 1e6
    Te = 10.0 # eV
    sp = Species(:e, -e_charge, m_e, n, Te, 0.0)
    # compute Debye length and choose k such that kλ_D ≈ 0.3
    vth = v_th(sp)
    ωp = omega_p(sp)
    λ_D = vth / (sqrt(2) * ωp)
    kλ = 0.3
    k = kλ / λ_D
    ω_num = kinetic_es_root(k, [sp])
    # predicted real frequency from Bohm–Gross
    ω_pred = sqrt(ωp^2 + 3 * k^2 * vth^2)
    @test isapprox(real(ω_num), ω_pred; rtol=1e-2)
    @test imag(ω_num) < 0.0
end

@testset "Ion acoustic wave" begin
    # Electrostatic wave in electron–ion plasma.  For moderate k, the
    # real part of the frequency should follow the ion acoustic speed:
    # ω ≈ k c_s with c_s² = (k_B/m_i)(T_e + 3 T_i).  Choose kλ_D small
    # so that ions are slow compared to electrons but warm corrections
    # still apply.  Note that tolerances are fairly loose because the
    # kinetic model includes Landau damping and other effects not
    # captured by the fluid expression.
    n = 1e6
    Te = 5.0  # eV
    Ti = 1.0  # eV
    spe = Species(:e, -e_charge, m_e, n, Te, 0.0)
    spi = Species(:p, e_charge,  m_p, n, Ti, 0.0)
    vth_e = v_th(spe)
    ωpe = omega_p(spe)
    λ_De = vth_e / (sqrt(2) * ωpe)
    kλ = 0.1
    k = kλ / λ_De
    ω_num = kinetic_es_root(k, [spe, spi])
    # Ion acoustic speed (adiabatic electrons, warm ions factor 3)
    c_s = sqrt((PlasmaWaves.k_B * (Te * PlasmaWaves.eV + 3 * Ti * PlasmaWaves.eV)) / m_p)
    ω_pred = k * c_s
    @test isapprox(real(ω_num), ω_pred; rtol=0.1)
    @test imag(ω_num) < 0.0
end