"""
    wpol_helicity(S)

Helicity and ellipticity of the 3×3 spectral matrix `S` at one frequency, each a power-weighted mean over the rows of `S`.

Row `c` gives the field vector ``𝐮 ∝ S_{c,:}^*`` up to a phase. The polarization ellipse traced by
``Re(𝐮 e^{iφ})`` has semi-axes ``\\sqrt{(|𝐮|^2 ± |𝐮^T 𝐮|)/2}`` and area ``π |Re 𝐮 × Im 𝐮|``, so its axial ratio is

```math
\\frac{2 |Re 𝐮 × Im 𝐮|}{|𝐮|^2 + |𝐮^T 𝐮|}.
```

Helicity uses all three components; ellipticity uses the (x, y) projection, signed by ``Im S_{12}``
(negative for left-hand rotation about z).

Each ratio is formed from numerator and denominator summed over rows, so row `c` counts ∝ ``S_{cc}`` for a pure wave.
PySPEDAS/IDL `wavpol` instead average the rows equally; for a parallel wave over white noise the δB_z row then holds only noise,
and helicity tends to 2/3 of its true value however high the SNR.
They also drop a factor 2 in the rotation angle of the (x, y) projection, so their ellipticity differs for tilted ellipses.
"""
function wpol_helicity(S)
    T = real(eltype(S))
    hn = hd = en = ed = zero(T)
    for c in 1:3
        s1, s2, s3 = S[c, 1], S[c, 2], S[c, 3]
        # Re 𝐮 × Im 𝐮 up to sign, components (z, x, y)
        x12, x23, x31 = imag(s1 * conj(s2)), imag(s2 * conj(s3)), imag(s3 * conj(s1))
        hn += 2 * sqrt(x12^2 + x23^2 + x31^2)
        hd += abs2(s1) + abs2(s2) + abs2(s3) + abs(s1^2 + s2^2 + s3^2)
        en += 2 * abs(x12)
        ed += abs2(s1) + abs2(s2) + abs(s1^2 + s2^2)
    end
    return hn / hd, sign(imag(S[1, 2])) * en / ed
end

@deprecate wpol_helicity(S, waveangle) wpol_helicity(S) false
