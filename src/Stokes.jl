"""
Stokes Parameters

Set of values that describe the polarization state of electromagnetic radiation

- [Wikipedia](https://en.wikipedia.org/wiki/Stokes_parameters)
"""
struct StokesParameters{T}
    S0::T
    S1::T
    S2::T
    S3::T
end


"""
    polarization(S0, S1, S2, S3)
    polarization(S::StokesParameters)

Compute the degree of polarization (p) from Stoke parameters or a Stokes vector.

# Reference
- [Wikipedia](https://en.wikipedia.org/wiki/Polarization_(waves))
- [Stokes parameters](https://en.wikipedia.org/wiki/Stokes_parameters)
"""
function polarization(S0, S1, S2, S3)
    return sqrt(S1^2 + S2^2 + S3^2) / S0
end

polarization(S::StokesParameters) = polarization(S.S0, S.S1, S.S2, S.S3)
