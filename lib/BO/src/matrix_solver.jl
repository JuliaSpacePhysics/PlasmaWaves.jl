# Matrix eigenvalue solver for plasma dispersion relations
# Ported from MATLAB boarbitrary/em3dHH/boHH code by Hua-sheng XIE

using SparseArrays

export solve_dispersion_matrix, JPoleCoefficients, get_jpole_coefficients

"""
    JPoleCoefficients{T}

J-pole approximation coefficients for the plasma dispersion function.
Z(ζ) ≈ Σⱼ bⱼ/(ζ - cⱼ)

See: Xie & Xiao, Plasma Science and Technology 18(2), 97 (2016)
"""
struct JPoleCoefficients{T <: Complex}
    J::Int
    bzj::Vector{T}
    czj::Vector{T}
end

"""
    get_jpole_coefficients(J=8)

Get J-pole approximation coefficients for the plasma dispersion function.

Supported values: J ∈ {4, 6, 8, 10, 12, 16, 20, 24, 28, 32}
"""
function get_jpole_coefficients(J::Int = 8)
    T = ComplexF64

    if J == 8
        bzj = T[
            -0.017340112270401 - 0.046306439626294im,
            -0.739917811220052 + 0.839518284620274im,
            5.840632105105495 + 0.95360275132204im,
            -5.583374181615043 - 11.208550459628098im,
        ]
        czj = T[
            2.237687725134293 - 1.625941024120362im,
            1.465234091939142 - 1.789620299603315im,
            0.839253966367922 - 1.891995211531426im,
            0.273936218055381 - 1.941787037576095im,
        ]
    elseif J == 12
        bzj = T[
            -10.020983259474214017 - 14.728932929429874883im,
            -0.58878169153449514493 + 0.19067303610080007359im,
            -0.27475707659732384029 + 3.617920717493884482im,
            0.00045713742777499515344 + 0.00027155393843737098852im,
            0.017940627032508378515 - 0.036436053276701248142im,
            10.366124263145749629 - 2.5069048649816145967im,
        ]
        czj = T[
            0.22660012611958088507627 - 2.0716877594897791206264im,
            -1.70029215163003500750575 - 1.8822474221612724460388im,
            1.17139325085601178534269 - 1.97725033192085410977458im,
            3.0666201126826972102007 - 1.59002082593259971758095im,
            2.307327490410578276422 - 1.7546732543728200653674im,
            0.687200524906019065672977 - 2.040288525975844018682im,
        ]
    elseif J == 4
        bzj = T[
            0.546796859834032 + 0.037196505239277im,
            -1.046796859834027 + 2.101852568038518im,
        ]
        czj = T[
            1.23588765343592 - 1.21498213255731im,
            -0.378611612386277 - 1.35094358543273im,
        ]
    elseif J == 6
        bzj = T[
            -2.0348129516220801032 - 3.696953928431426183im,
            -0.052066769165630992082 + 0.13185419997353193789im,
            1.5868797207877110952 - 0.012584451331815698194im,
        ]
        czj = T[
            0.31420182117748071085 - 1.6013905451342740139im,
            1.7917787676379544378 - 1.3594900979984264658im,
            -0.97975990624668776443 - 1.5226788060732615273im,
        ]
    elseif J == 10
        bzj = T[
            -11.096022667802034791 + 20.894981637405839142im,
            -2.0505302680318152519 + 2.8357523253617763706im,
            -0.17457143444095772263 + 0.33944656050059462742im,
            0.0092318520465358982968 - 0.00045136341696997654943im,
            12.811892518228271867 + 1.1703829347651599709im,
        ]
        czj = T[
            -0.24645915590151640764 - 2.1213597368144015061im,
            1.2868093343179118256 - 2.0098261125516749386im,
            -1.8934221509724099323 - 1.8942779865460054073im,
            2.6471246045959517878 - 1.7296444188113259578im,
            0.74968316099982057453 - 2.0845790675218448384im,
        ]
    elseif J == 16
        bzj = T[
            -86.416592794839804566 - 147.57960545984972964im,
            -22.962540986214500398 + 46.211318219085729914im,
            -8.8757833558787660662 - 11.561957978688249474im,
            -0.025134802434111256483 + 0.19730442150379382482im,
            -0.0056462830661756538039 - 0.0027884991898011769583im,
            0.000028262945845046458372 + 0.000026335348714810255537im,
            2.3290098166119338312 - 0.57238325918028725167im,
            115.45666014287557906 - 2.8617578808752183449im,
        ]
        czj = T[
            0.1966439744113664608458045976392 - 2.5854046363167904820930238267552im,
            1.000427687089304511157923736374 - 2.5277610669350594581215470753126im,
            1.4263380087098663428834704281261 - 2.4694803409658086505344546718783im,
            2.382753075769737513956410751299 - 2.2903917960623787648467068236658im,
            2.9566517643704010426779572638885 - 2.1658992556376956216621262314559im,
            -3.6699741330155866185481740497527 - 2.008727613312046260114172119472im,
            1.8818356204685089975461092960437 - 2.3907395820644127767937911780402im,
            0.5933003629474285223202174828712 - 2.5662607006180515205167080595386im,
        ]
    elseif J == 20
        bzj = T[
            2.370563868103888e-9 - 2.539993703144051e-7im,
            9.066371501956585e-5 + 1.467741891511021e-5im,
            -0.001548195081695735 + 0.005538234797273276im,
            -0.1187449908723126 - 0.05251079094649405im,
            0.8189033525844321 - 1.174311257342413im,
            5.980769211588994 + 6.826225337943291im,
            -33.13580595567155 + 15.31444785352592im,
            -10.96002434997827 - 98.56499306153125im,
            184.0604900158517 + 44.51786749043065im,
            -147.1441297545068 + 214.6355315521182im,
        ]
        czj = T[
            -4.289538835361067 - 1.989275288759969im,
            -3.583937986426706 - 2.139183234944734im,
            -3.01512916068484 - 2.263535529839947im,
            -2.519734806587566 - 2.36995157596592im,
            -2.072161079092204 - 2.460331650306395im,
            -1.658055996353741 - 2.53524891246743im,
            -1.267869984223448 - 2.594931353792834im,
            -0.8945422232244176 - 2.639529570725944im,
            -0.5324423359117022 - 2.669173720892624im,
            -0.1767811345348292 - 2.683966971857356im,
        ]
    elseif J == 24
        bzj = T[
            -579.77656932346560644 - 844.01436313629880827im,
            -179.52530851977905732 - 86.660002027244731382im,
            -52.107235029274485215 + 453.3246806707749413im,
            -2.1607927691932962178 + 0.63681255371973499384im,
            -0.018283386874895507814 - 0.21941582055233427677im,
            -0.00006819511737162705016 + 0.00032026091897256872621im,
            -0.0000028986123310445793648 - 0.00000099510625011385493369im,
            0.0000000023382228949223867744 - 0.0000000040404517369565098657im,
            0.01221466589423530596 + 0.00097890737323377354166im,
            7.3718296773233126912 - 12.575687057120635407im,
            44.078424019374375065 - 46.322124026599601416im,
            761.62579175738689742 + 185.11797721443392707im,
        ]
        czj = T[
            0.16167711630587375808393823760988 - 2.9424665391729649010502939606152im,
            1.1509135876493567244599398043479 - 2.8745542965490153159866506667543im,
            0.81513635269214329286824152984179 - 2.9085569383176322446978082849749im,
            2.2362950589041724110736073820844 - 2.7033607074680388479084431872604im,
            2.6403561313404041541230494846625 - 2.6228400297078984516779261304916im,
            3.5620497451197056657834990483967 - 2.4245607245823420555878190731282im,
            4.116925125710675393072860873751 - 2.3036541720854573608940600179944im,
            4.8034117493360317933109830717707 - 2.1592490859689535412501218722927im,
            3.0778922349246567316482750461458 - 2.5301774598854448463007864644617im,
            -1.8572088635240765003561090479193 - 2.7720571884094886583775397071469im,
            1.4969881322466893380396663902149 - 2.8290855580900544693059801078858im,
            -0.48636891219330428093331493852099 - 2.9311741817223824196339069754696im,
        ]
    elseif J == 28
        bzj = T[
            -6.990333597785891e-11 - 3.844801630174103e-11im,
            5.958884394492422e-8 - 6.533760978309747e-8im,
            9.644940137755459e-6 + 1.156039714667892e-5im,
            -0.0007360376178448977 + 0.0004706874270260621im,
            -0.009829228795726045 - 0.02183839440379618im,
            0.3587753042250427 - 0.09013468287145511im,
            0.1264786525009892 + 3.589874046687523im,
            -23.07993841496687 - 4.901784585465681im,
            51.42260618075652 - 97.65815979660832im,
            269.3151882535272 + 270.2406101655737im,
            -897.6534736073505 + 443.2670670715868im,
            -235.2064726575382 - 2019.283674665684im,
            3145.022973432025 + 767.9240405892624im,
            -2310.795581581225 + 3359.074813701693im,
        ]
        czj = T[
            -5.274315519672472 - 2.315931599664814im,
            -4.60385910024417 - 2.455760440449168im,
            -4.060876884075117 - 2.573378587177523im,
            -3.586314479330714 - 2.677316091557037im,
            -3.1570693599154 - 2.770199141723955im,
            -2.760754763795907 - 2.853038026706498im,
            -2.389521521899566 - 2.926259908085017im,
            -2.037841188423438 - 2.990064209207999im,
            -1.701527749822992 - 3.044563094147437im,
            -1.377240603887195 - 3.089837634024895im,
            -1.062202056569621 - 3.125957495198709im,
            -0.754020037194251 - 3.152984712343023im,
            -0.4505661571573089 - 3.170971180914522im,
            -0.1498842415719675 - 3.179954519589582im,
        ]
    elseif J == 32
        bzj = T[
            1.288000175155747e-12 - 6.165437822620772e-14im,
            -5.21842626011812e-10 + 2.22484909250969e-9im,
            -5.45993065604352e-7 - 2.213024344007754e-7im,
            2.482641796977981e-5 - 4.392260026870976e-5im,
            0.00159770023288028 + 0.001218924472518216im,
            -0.03236028144991428 + 0.03075247101343072im,
            -0.3356879616865236 - 0.5202286015269713im,
            5.420360688192178 - 2.019925251380004im,
            4.278400838512319 + 38.19078666688344im,
            -186.4720211995964 - 30.1728826831305im,
            311.9071201210668 - 635.2557114736077im,
            1479.153702668081 + 1443.321511690421im,
            -4278.767516753266 + 2140.147973966637im,
            -1033.23156980727 - 8828.667007433389im,
            12987.60133150033 + 3183.861482228732im,
            -9290.023381793046 + 13489.51743209337im,
        ]
        czj = T[
            -5.711386487786359 - 2.462045142200075im,
            -5.054682350875668 - 2.598001772487061im,
            -4.522048874539196 - 2.712627489029625im,
            -4.055893312027004 - 2.81462581890055im,
            -3.633837632185388 - 2.906815186530633im,
            -3.244016833346553 - 2.990356362011712im,
            -2.879018027888166 - 3.065777008834288im,
            -2.53369893303409 - 3.133332962829416im,
            -2.204224090140716 - 3.193159826311868im,
            -1.887577844101815 - 3.245341303953im,
            -1.581292110801636 - 3.289939833174147im,
            -1.283280815643552 - 3.327009004303097im,
            -0.9917312771695507 - 3.35659712135624im,
            -0.7050277254451327 - 3.378746697704547im,
            -0.421693707628128 - 3.39349241608021im,
            -0.1403458271046983 - 3.400858886393542im,
        ]
    else
        error("Unsupported J-pole number: $J. Supported: 4, 6, 8, 10, 12, 16, 20, 24, 28, 32")
    end

    # Add conjugate pairs
    Jhalf = length(bzj)
    bzj_full = vcat(bzj, conj.(bzj))
    czj_full = vcat(czj, -conj.(czj))

    return JPoleCoefficients{T}(J, bzj_full, czj_full)
end

"""
    MatrixSolverParams{T}

Parameters for the matrix eigenvalue solver.
"""
struct MatrixSolverParams{T}
    S::Int                      # Number of species
    c2::T                       # Speed of light squared
    wcs::Vector{T}              # Cyclotron frequencies
    wps2::Vector{T}             # Plasma frequencies squared
    rhocs::Vector{T}            # Cyclotron radii
    vtzs::Vector{T}             # Parallel thermal velocities
    vtps::Vector{T}             # Perpendicular thermal velocities
    vdsz::Vector{T}             # Parallel drift velocities
    ds::Vector{T}               # Ring beam drift parameters
    As::Vector{T}               # Normalization parameters
    Nss::Vector{Int}            # Maximum harmonic numbers
    aslm::Vector{Matrix{T}}     # Hermite expansion coefficients
    msmax::Vector{Int}          # Max perpendicular Hermite index
    lsmax::Vector{Int}          # Max parallel Hermite index
    jpole::JPoleCoefficients{Complex{T}}  # J-pole coefficients
end

"""
    solve_dispersion_matrix(params, kx, kz; J=8)

Solve the kinetic dispersion relation using the matrix eigenvalue method.

Returns all eigenfrequencies ω(k) for the given wave vector (kx, kz).

This method transforms the dispersion relation into a matrix eigenvalue problem
using J-pole approximation for the plasma dispersion function, allowing
simultaneous computation of all wave modes.

# Arguments
- `params`: HHSolverParams or MatrixSolverParams with plasma parameters
- `kx`: Perpendicular wave vector component (m⁻¹)
- `kz`: Parallel wave vector component (m⁻¹)
- `J`: Number of poles for Z-function approximation (default: 8)

See also: [`fDrHH`](@ref), [`get_jpole_coefficients`](@ref)
"""
function solve_dispersion_matrix(params::HHSolverParams, kx, kz; J = 8)
    jpole = get_jpole_coefficients(J)
    mparams = MatrixSolverParams(
        params.S, params.c2, params.wcs, params.wps2, params.rhocs,
        params.vtzs, params.vtps, params.vdsz, params.ds, params.As,
        params.Nss, params.aslm, params.msmax, params.lsmax, jpole
    )
    return solve_dispersion_matrix(mparams, kx, kz)
end

function solve_dispersion_matrix(params::MatrixSolverParams, kx, kz)
    (; S, c2, wcs, wps2, rhocs, vtzs, vtps, vdsz, ds, As, Nss, aslm, msmax, lsmax, jpole) = params
    (; J, bzj, czj) = jpole

    # Handle singularities
    kx_local = kx == 0.0 ? 1.0e-30 : kx
    kz_local = kz

    # Perpendicular wavenumber parameter
    as = kx_local .* rhocs .* sqrt(2)

    # Compute matrix dimensions
    Ns = 2 .* Nss .+ 1  # Number of harmonics per species
    SNJ = sum(Ns) * J
    SNJ1 = SNJ + S
    SNJ3 = 3 * SNJ1
    NN = SNJ3 + 6

    # Adjust czj for kz sign
    czjj = kz_local < 0 ? -czj : czj

    # Precompute czjj^l for all l values
    max_lsmax = maximum(lsmax)
    czjj_l = zeros(ComplexF64, J, max_lsmax + 5)
    for l in 0:(max_lsmax + 3)
        czjj_l[:, l + 2] .= czjj .^ l
    end

    # Initialize coefficient arrays
    b11s = zeros(ComplexF64, S)
    b12s = zeros(ComplexF64, S)
    b13s = zeros(ComplexF64, S)
    b21s = zeros(ComplexF64, S)
    b22s = zeros(ComplexF64, S)
    b23s = zeros(ComplexF64, S)
    b31s = zeros(ComplexF64, S)
    b32s = zeros(ComplexF64, S)
    b33s = zeros(ComplexF64, S)

    csnj = zeros(ComplexF64, SNJ)
    b11snj = zeros(ComplexF64, SNJ)
    b12snj = zeros(ComplexF64, SNJ)
    b13snj = zeros(ComplexF64, SNJ)
    b21snj = zeros(ComplexF64, SNJ)
    b22snj = zeros(ComplexF64, SNJ)
    b23snj = zeros(ComplexF64, SNJ)
    b31snj = zeros(ComplexF64, SNJ)
    b32snj = zeros(ComplexF64, SNJ)
    b33snj = zeros(ComplexF64, SNJ)

    snj = 0
    for s in 1:S
        for n in -Nss[s]:Nss[s]
            # Compute perpendicular integrals (only once per n)
            Ans = zeros(msmax[s] + 4, 2)
            Bns = zeros(msmax[s] + 4, 2)
            Cns = zeros(msmax[s] + 4, 2)

            for m in 0:(msmax[s] + 2)
                idm = m + 2
                Ans[idm, 1] = 2.0 / As[s] * funAn(n, as[s], ds[s], m, 0)
                Ans[idm, 2] = 2.0 / As[s] * funAn(n, as[s], ds[s], m, 1)
                Bns[idm, 1] = 2.0 / As[s] * funBn(n, as[s], ds[s], m, 1)
                Bns[idm, 2] = 2.0 / As[s] * funBn(n, as[s], ds[s], m, 2)
                Cns[idm, 1] = 2.0 / As[s] * funCn(n, as[s], ds[s], m, 2)
                Cns[idm, 2] = 2.0 / As[s] * funCn(n, as[s], ds[s], m, 3)
            end

            for j in 1:J
                snj += 1

                # Pole location
                csnj[snj] = czjj[j] * kz_local * vtzs[s] + kz_local * vdsz[s] + n * wcs[s]
                cnj = csnj[snj]

                # Sum over Hermite indices
                sum11tmp1, sum11tmp2, sum11tmp3 = 0.0im, 0.0im, 0.0im
                sum12tmp1, sum12tmp2 = 0.0im, 0.0im
                sum22tmp1, sum22tmp2, sum22tmp3 = 0.0im, 0.0im, 0.0im
                sum13tmp1, sum13tmp2 = 0.0im, 0.0im
                sum23tmp1, sum23tmp2 = 0.0im, 0.0im
                sum32tmp1, sum32tmp2 = 0.0im, 0.0im
                sum33tmp1, sum33tmp2, sum33tmp3 = 0.0im, 0.0im, 0.0im

                dr = vdsz[s] / vtzs[s]

                for l in 0:lsmax[s], m in 0:msmax[s]
                    coeff = aslm[s][l + 1, m + 1]
                    idx_mp1 = m + 1 + 2  # m+1+2
                    idx_mm1 = m - 1 + 2  # m-1+2
                    idx_m = m + 2        # m+2

                    dAm = 2 * Ans[idx_mp1, 1] - m * Ans[idx_mm1, 1]
                    dBm = 2 * Bns[idx_mp1, 1] - m * Bns[idx_mm1, 1]
                    dCm = 2 * Cns[idx_mp1, 1] - m * Cns[idx_mm1, 1]

                    # czjj_l indices: l+2 for l, l+1+2 for l+1, l-1+2 for l-1, etc.
                    czl = czjj_l[j, l + 2]
                    czlp1 = czjj_l[j, l + 1 + 2]
                    czlm1 = l >= 1 ? czjj_l[j, l - 1 + 2] : 0.0im
                    czlp2 = czjj_l[j, l + 2 + 2]
                    czlp3 = czjj_l[j, l + 3 + 2]

                    dZl = 2 * czlp1 - l * czlm1

                    sum11tmp1 += coeff * czl * dAm
                    sum11tmp2 += coeff * dZl * Ans[idx_m, 2]

                    sum12tmp1 += coeff * czl * dBm
                    sum12tmp2 += coeff * dZl * Bns[idx_m, 2]

                    sum22tmp1 += coeff * czl * dCm
                    sum22tmp2 += coeff * dZl * Cns[idx_m, 2]

                    sum13tmp1 += coeff * (dr * czl + czlp1) * dAm
                    sum13tmp2 += coeff * ((2 * czlp2 - l * czl) + dr * dZl) * Ans[idx_m, 2]

                    sum23tmp1 += coeff * (dr * czl + czlp1) * dBm
                    sum23tmp2 += coeff * ((2 * czlp2 - l * czl) + dr * dZl) * Bns[idx_m, 2]

                    sum32tmp1 += coeff * (dr * czl + czlp1) * dBm
                    sum32tmp2 += coeff * ((2 * czlp2 - l * czl) + dr * dZl) * Bns[idx_m, 2]

                    sum33tmp1 += coeff * (dr^2 * czl + 2 * dr * czlp1 + czlp2) * dAm
                    sum33tmp2 += coeff * (
                        dr^2 * dZl + 2 * dr * (2 * czlp2 - l * czl) +
                            (2 * czlp3 - l * czlp1)
                    ) * Ans[idx_m, 2]

                    if j == 1
                        Il = funIn(l)
                        Ilp1 = funIn(l + 1)
                        Ilp2 = funIn(l + 2)
                        Ilm1 = l >= 1 ? funIn(l - 1) : 0.0

                        sum11tmp3 += coeff * Il * dAm
                        sum22tmp3 += coeff * Il * dCm
                        sum33tmp3 += coeff * (dr * (2 * Ilp1 - l * Ilm1) + (2 * Ilp2 - l * Il)) * Ans[idx_m, 2]
                    end
                end

                # Add background term for j=1
                if j == 1
                    nwkp = n * wcs[s] / (kx_local * vtps[s])
                    p11snj_b = nwkp^2 * sum11tmp3
                    p22snj_b = sum22tmp3
                    p33snj_b = sum33tmp3
                    b11s[s] -= wps2[s] * p11snj_b
                    b22s[s] -= wps2[s] * p22snj_b
                    b33s[s] -= wps2[s] * p33snj_b
                end

                tmp = wps2[s] * bzj[j] / cnj
                vr = vtps[s] / vtzs[s]
                nwkp = n * wcs[s] / (kx_local * vtps[s])

                p11snj = nwkp^2 * (n * wcs[s] * sum11tmp1 + kz_local * vtzs[s] * vr^2 * sum11tmp2)
                p12snj = 1im * nwkp * (n * wcs[s] * sum12tmp1 + kz_local * vtzs[s] * vr^2 * sum12tmp2)
                p21snj = -p12snj
                p22snj = n * wcs[s] * sum22tmp1 + kz_local * vtzs[s] * vr^2 * sum22tmp2
                p13snj = nwkp * ((vtzs[s] / vtps[s]) * n * wcs[s] * sum13tmp1 + kz_local * vtps[s] * sum13tmp2)
                p31snj = p13snj
                p23snj = -1im * ((vtzs[s] / vtps[s]) * n * wcs[s] * sum23tmp1 + kz_local * vtps[s] * sum23tmp2)
                p32snj = 1im * ((vtzs[s] / vtps[s]) * n * wcs[s] * sum32tmp1 + kz_local * vtps[s] * sum32tmp2)
                p33snj = (vtzs[s] / vtps[s]) * (
                    (vtzs[s] / vtps[s]) * n * wcs[s] * sum33tmp1 +
                        kz_local * vtps[s] * sum33tmp2
                )

                b11snj[snj] += tmp * p11snj
                b11s[s] -= tmp * p11snj
                b12snj[snj] += tmp * p12snj
                b12s[s] -= tmp * p12snj
                b21snj[snj] += tmp * p21snj
                b21s[s] -= tmp * p21snj
                b22snj[snj] += tmp * p22snj
                b22s[s] -= tmp * p22snj
                b13snj[snj] += tmp * p13snj
                b13s[s] -= tmp * p13snj
                b31snj[snj] += tmp * p31snj
                b31s[s] -= tmp * p31snj
                b23snj[snj] += tmp * p23snj
                b23s[s] -= tmp * p23snj
                b32snj[snj] += tmp * p32snj
                b32s[s] -= tmp * p32snj
                b33snj[snj] += tmp * p33snj
                b33s[s] -= tmp * p33snj
            end
        end
    end

    # Build sparse matrix
    I_idx = Int[]
    J_idx = Int[]
    V_val = ComplexF64[]

    function add_entry!(i, j, v)
        return if v != 0
            push!(I_idx, i)
            push!(J_idx, j)
            push!(V_val, v)
        end
    end

    for snj_idx in 1:SNJ
        jjx = snj_idx + 0 * SNJ1
        jjy = snj_idx + 1 * SNJ1
        jjz = snj_idx + 2 * SNJ1

        # v_snjx equations
        add_entry!(jjx, jjx, csnj[snj_idx])
        add_entry!(jjx, SNJ3 + 1, b11snj[snj_idx])
        add_entry!(jjx, SNJ3 + 2, b12snj[snj_idx])
        add_entry!(jjx, SNJ3 + 3, b13snj[snj_idx])

        # v_snjy equations
        add_entry!(jjy, jjy, csnj[snj_idx])
        add_entry!(jjy, SNJ3 + 1, b21snj[snj_idx])
        add_entry!(jjy, SNJ3 + 2, b22snj[snj_idx])
        add_entry!(jjy, SNJ3 + 3, b23snj[snj_idx])

        # v_snjz equations
        add_entry!(jjz, jjz, csnj[snj_idx])
        add_entry!(jjz, SNJ3 + 1, b31snj[snj_idx])
        add_entry!(jjz, SNJ3 + 2, b32snj[snj_idx])
        add_entry!(jjz, SNJ3 + 3, b33snj[snj_idx])
    end

    # E(J) coupling: J_xyz = j_xyz + sum(v_snj_xyz)
    for jj in (0 * SNJ1 + 1):(1 * SNJ1)
        add_entry!(SNJ3 + 1, jj, -1.0)
    end
    for jj in (1 * SNJ1 + 1):(2 * SNJ1)
        add_entry!(SNJ3 + 2, jj, -1.0)
    end
    for jj in (2 * SNJ1 + 1):(3 * SNJ1)
        add_entry!(SNJ3 + 3, jj, -1.0)
    end

    # Species current contributions
    for s in 1:S
        add_entry!(1 * SNJ1 - S + s, SNJ3 + 1, b11s[s])
        add_entry!(1 * SNJ1 - S + s, SNJ3 + 2, b12s[s])
        add_entry!(1 * SNJ1 - S + s, SNJ3 + 3, b13s[s])
        add_entry!(2 * SNJ1 - S + s, SNJ3 + 1, b21s[s])
        add_entry!(2 * SNJ1 - S + s, SNJ3 + 2, b22s[s])
        add_entry!(2 * SNJ1 - S + s, SNJ3 + 3, b23s[s])
        add_entry!(3 * SNJ1 - S + s, SNJ3 + 1, b31s[s])
        add_entry!(3 * SNJ1 - S + s, SNJ3 + 2, b32s[s])
        add_entry!(3 * SNJ1 - S + s, SNJ3 + 3, b33s[s])
    end

    # E(B) coupling: Maxwell's equations
    add_entry!(SNJ3 + 1, SNJ3 + 5, c2 * kz_local)
    add_entry!(SNJ3 + 2, SNJ3 + 4, -c2 * kz_local)
    add_entry!(SNJ3 + 2, SNJ3 + 6, c2 * kx_local)
    add_entry!(SNJ3 + 3, SNJ3 + 5, -c2 * kx_local)

    # B(E) coupling: Faraday's law
    add_entry!(SNJ3 + 4, SNJ3 + 2, -kz_local)
    add_entry!(SNJ3 + 5, SNJ3 + 1, kz_local)
    add_entry!(SNJ3 + 5, SNJ3 + 3, -kx_local)
    add_entry!(SNJ3 + 6, SNJ3 + 2, kx_local)

    M = sparse(I_idx, J_idx, V_val, NN, NN)

    # Solve eigenvalue problem
    eigenvalues = eigen(Matrix(M)).values

    return eigenvalues
end
