"""
    polarization(S0, S1, S2, S3)

Degree of polarization from the [Stokes parameters](https://en.wikipedia.org/wiki/Stokes_parameters).
"""
function polarization(S0, S1, S2, S3)
    return sqrt(S1^2 + S2^2 + S3^2) / S0
end
