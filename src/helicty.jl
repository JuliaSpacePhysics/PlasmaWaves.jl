"""
Phase factor `exp (i φ)` satisfies the following equation

``\\exp (4 i φ) = \\exp (-2 i γ)``

where

``γ = \\arctan (2 Re(𝐮)^𝐓 Im(𝐮) /(Re(𝐮)^2-Im(𝐮)^2))``
"""
@inline function phase_factor(u1, u2, u3)
    upper = 2 * (real(u1) * imag(u1) + real(u2) * imag(u2) + real(u3) * imag(u3))
    lower = real(u1)^2 + real(u2)^2 + real(u3)^2 - imag(u1)^2 - imag(u2)^2 - imag(u3)^2
    γ = atan(upper, lower)
    return exp(-im * γ / 2)
end

"""
    wpol_helicity(S, waveangle)

Compute helicity and ellipticity for a single frequency.

# Arguments
- `S`: Spectral matrix for a single frequency, size (3,3)
- `waveangle`: Wave normal angle for this frequency

# Returns
- `helicity`: Average helicity across the three components
- `ellipticity`: Average ellipticity across the three components
"""
function wpol_helicity(S, waveangle)
    helicity = 0
    ellipticity = 0

    for comp in 1:3
        alph = sqrt(real(S[comp, comp]))
        if comp == 1
            l1 = complex(alph)
            l2 = (real(S[1, 2]) / alph) + im * (-imag(S[1, 2]) / alph)
            l3 = (real(S[1, 3]) / alph) + im * (-imag(S[1, 3]) / alph)
        elseif comp == 2
            l1 = (real(S[2, 1]) / alph) + im * (-imag(S[2, 1]) / alph)
            l2 = complex(alph)
            l3 = (real(S[2, 3]) / alph) + im * (-imag(S[2, 3]) / alph)
        else
            l1 = (real(S[3, 1]) / alph) + im * (-imag(S[3, 1]) / alph)
            l2 = (real(S[3, 2]) / alph) + im * (-imag(S[3, 2]) / alph)
            l3 = complex(alph)
        end

        p = phase_factor(l1, l2, l3)
        y1 = p * l1
        y2 = p * l2
        y3 = p * l3

        norm_real = sqrt(real(y1)^2 + real(y2)^2 + real(y3)^2)
        norm_imag = sqrt(imag(y1)^2 + imag(y2)^2 + imag(y3)^2)
        helicity += norm_imag / norm_real / 3

        # TODO: why there is no 2 in front of uppere for `PySPEDAS`
        uppere = 2 * (imag(y1) * real(y1) + imag(y2) * real(y2))
        lowere = -imag(y1)^2 + real(y1)^2 - imag(y2)^2 + real(y2)^2
        gammarot = atan(uppere, lowere)
        p_rot = exp(-1im * 0.5 * gammarot)
        y1_rot = p_rot * y1
        y2_rot = p_rot * y2

        num = sqrt(imag(y1_rot)^2 + imag(y2_rot)^2)
        den = sqrt(real(y1_rot)^2 + real(y2_rot)^2)
        ellip_val = (den != 0) ? num / den : NaN
        # Adjust sign using the off-diagonal of ematspec and the wave normal angle
        sign_factor = sign(imag(S[1, 2]) * sin(waveangle))
        ellipticity += ellip_val * sign_factor / 3
    end
    return helicity, ellipticity
end
