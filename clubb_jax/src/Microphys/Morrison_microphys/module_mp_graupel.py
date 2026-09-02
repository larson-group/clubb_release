"""JAX port of Morrison 2-moment microphysics special functions (module_mp_graupel.F90).

Foundational, independently-verifiable pieces (port these first, KK-playbook):
  - polysvp  : saturation vapor pressure [Pa] (Flatau 1992 Table 4, the "V1.7" coeffs)
  - gamma    : the complete gamma function Γ(x) (W. J. Cody algorithm)
  - derf1    : the error function erf(x) (Morrison's DERF1 rational approximation)

Each mirrors the Fortran exactly (Horner nesting / branch structure preserved) so the
later process rates that call them are bit-faithful. Validated vs scipy + known values
in tests/test_morrison_special.py. Morrison runs FLOAT64 in the CLUBB build (the inner WRF
graupel scheme's bare `REAL` is promoted by -fdefault-real-8; see DESIGN), so these float64
ports can be bit-faithful — reproduce the single-precision literals as their decimal values.

Mirror scope: module_mp_graupel.F90's `calc_refl10cm` + `rayleigh_soak_wetgraupel` are deliberately NOT
ported — they compute a radar-reflectivity (dBZ) DIAGNOSTIC for the WRF post-processor only, with no CLUBB
coupling (they don't feed any prognostic tendency) and no oracle path. The prognostic microphysics rates +
sedimentation + the CLUBB interface ARE all mirrored.
"""
import jax.numpy as jnp

# ── POLYSVP — saturation vapor pressure, Flatau et al. 1992 Table 4 (RHS column) ──────
# module_mp_graupel.F90:5423. Distinct from the CLUBB-core saturation.py Flatau fit
# (Morrison uses its own "V1.7" coefficients).
_POLYSVP_LIQ = (6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3,
                0.302950461e-5, 0.206739458e-7, 0.640689451e-10,
                -0.952447341e-13, -0.976195544e-15)
_POLYSVP_ICE = (6.11147274, 0.503160820, 0.188439774e-1, 0.420895665e-3,
                0.615021634e-5, 0.602588177e-7, 0.385852041e-9,
                0.146898966e-11, 0.252751365e-14)


def polysvp(T, itype):
    """Saturation vapor pressure [Pa]. T in K. itype: 0=liquid, 1=ice (static int).

    Faithful to module_mp_graupel.F90:POLYSVP — `dt = max(-80, T-273.16)`, Horner-nested
    8th-order polynomial, ×100 to Pa. At the triple point (T=273.16, dt=0) returns a0·100
    (611.239921 liquid / 611.147274 ice).
    """
    a = _POLYSVP_ICE if itype == 1 else _POLYSVP_LIQ
    dt = jnp.maximum(-80.0, T - 273.16)
    # Exact Fortran nesting: a0 + dt*(a1 + dt*(a2 + ... + dt*(a7 + a8*dt)...))
    p = a[7] + a[8] * dt
    for c in (a[6], a[5], a[4], a[3], a[2], a[1], a[0]):
        p = c + dt * p
    return p * 100.0


# ── DERF1 — error function erf(x) (Ooura table approximation) ─────────────────────────
# module_mp_graupel.F90:5728. Three argument-range branches; the coefficient block is
# selected by INT(w^2) (w<2.2) or INT(w) (2.2<=w<6.9). A/B are 5 blocks of 13.
_DERF1_A = (
    # block 0
    0.00000000005958930743, -0.00000000113739022964, 0.00000001466005199839,
    -0.00000016350354461960, 0.00000164610044809620, -0.00001492559551950604,
    0.00012055331122299265, -0.00085483269811296660, 0.00522397762482322257,
    -0.02686617064507733420, 0.11283791670954881569, -0.37612638903183748117,
    1.12837916709551257377,
    # block 1
    0.00000000002372510631, -0.00000000045493253732, 0.00000000590362766598,
    -0.00000006642090827576, 0.00000067595634268133, -0.00000621188515924000,
    0.00005103883009709690, -0.00037015410692956173, 0.00233307631218880978,
    -0.01254988477182192210, 0.05657061146827041994, -0.21379664776456006580,
    0.84270079294971486929,
    # block 2
    0.00000000000949905026, -0.00000000018310229805, 0.00000000239463074000,
    -0.00000002721444369609, 0.00000028045522331686, -0.00000261830022482897,
    0.00002195455056768781, -0.00016358986921372656, 0.00107052153564110318,
    -0.00608284718113590151, 0.02986978465246258244, -0.13055593046562267625,
    0.67493323603965504676,
    # block 3
    0.00000000000382722073, -0.00000000007421598602, 0.00000000097930574080,
    -0.00000001126008898854, 0.00000011775134830784, -0.00000111992758382650,
    0.00000962023443095201, -0.00007404402135070773, 0.00050689993654144881,
    -0.00307553051439272889, 0.01668977892553165586, -0.08548534594781312114,
    0.56909076642393639985,
    # block 4
    0.00000000000155296588, -0.00000000003032205868, 0.00000000040424830707,
    -0.00000000471135111493, 0.00000005011915876293, -0.00000048722516178974,
    0.00000430683284629395, -0.00003445026145385764, 0.00024879276133931664,
    -0.00162940941748079288, 0.00988786373932350462, -0.05962426839442303805,
    0.49766113250947636708,
)
_DERF1_B = (
    # block 0 (B 0-12)
    -0.00000000029734388465, 0.00000000269776334046, -0.00000000640788827665,
    -0.00000001667820132100, -0.00000021854388148686, 0.00000266246030457984,
    0.00001612722157047886, -0.00025616361025506629, 0.00015380842432375365,
    0.00815533022524927908, -0.01402283663896319337, -0.19746892495383021487,
    0.71511720328842845913,
    # block 1 (B 13-25)
    -0.00000000001951073787, -0.00000000032302692214, 0.00000000522461866919,
    0.00000000342940918551, -0.00000035772874310272, 0.00000019999935792654,
    0.00002687044575042908, -0.00011843240273775776, -0.00080991728956032271,
    0.00661062970502241174, 0.00909530922354827295, -0.20160072778491013140,
    0.51169696718727644908,
    # block 2 (B 26-38)
    0.00000000003147682272, -0.00000000048465972408, 0.00000000063675740242,
    0.00000003377623323271, -0.00000015451139637086, -0.00000203340624738438,
    0.00001947204525295057, 0.00002854147231653228, -0.00101565063152200272,
    0.00271187003520095655, 0.02328095035422810727, -0.16725021123116877197,
    0.32490054966649436974,
    # block 3 (B 39-51)
    0.00000000002319363370, -0.00000000006303206648, -0.00000000264888267434,
    0.00000002050708040581, 0.00000011371857327578, -0.00000211211337219663,
    0.00000368797328322935, 0.00009823686253424796, -0.00065860243990455368,
    -0.00075285814895230877, 0.02585434424202960464, -0.11637092784486193258,
    0.18267336775296612024,
    # block 4 (B 52-64)
    -0.00000000000367789363, 0.00000000020876046746, -0.00000000193319027226,
    -0.00000000435953392472, 0.00000018006992266137, -0.00000078441223763969,
    -0.00000675407647949153, 0.00008428418334440096, -0.00017604388937031815,
    -0.00239729611435071610, 0.02064129023876022970, -0.06905562880005864105,
    0.09084526782065478489,
)
_DERF1_A_BLOCKS = jnp.asarray(_DERF1_A).reshape(5, 13)
_DERF1_B_BLOCKS = jnp.asarray(_DERF1_B).reshape(5, 13)


def _horner13(coeffs, t):
    """Horner eval of 13 coeffs (coeffs[..., 13]) at t: c0·t^12 + ... + c12 (Fortran order)."""
    y = coeffs[..., 0]
    for i in range(1, 13):
        y = y * t + coeffs[..., i]
    return y


def derf1(x):
    """Error function erf(x), faithful to module_mp_graupel.F90:DERF1 (Ooura approximation).
    Array-capable; the coefficient block is gathered per element by the argument range.
    """
    x = jnp.asarray(x, dtype=jnp.float64)
    w = jnp.abs(x)

    # Branch 1: w < 2.2 — k=int(w^2) in 0..4, polynomial in t=frac(w^2), ×w
    t1 = w * w
    k1 = jnp.clip(jnp.floor(t1).astype(jnp.int32), 0, 4)
    t1f = t1 - k1
    y1 = _horner13(_DERF1_A_BLOCKS[k1], t1f) * w

    # Branch 2: 2.2 <= w < 6.9 — k=int(w) in 2..6, block k-2, polynomial in t=w-k, then 1-poly^16
    k2 = jnp.clip(jnp.floor(w).astype(jnp.int32), 2, 6)
    t2 = w - k2
    y2 = _horner13(_DERF1_B_BLOCKS[k2 - 2], t2)
    y2 = y2 * y2
    y2 = y2 * y2
    y2 = y2 * y2
    y2 = 1.0 - y2 * y2

    y = jnp.where(w < 2.2, y1, jnp.where(w < 6.9, y2, 1.0))
    return jnp.where(x < 0.0, -y, y)


# ── GAMMA — complete gamma function Γ(x) (W. J. Cody rational-minimax) ─────────────────
# module_mp_graupel.F90:5494. Branches: negative-arg reflection, small-arg 1/Y, (1,12)
# rational approx with integer argument reduction, ≥12 Stirling asymptotic.
_GAMMA_PI = 3.1415926535897932384626434
_GAMMA_SQRTPI = 0.9189385332046727417803297   # = 0.5*ln(2π), the log-gamma Stirling constant
_GAMMA_EPS = 1.19e-7        # single-precision literals (real8 in the -fdefault-real-8 build)
_GAMMA_XBIG = 35.040
_GAMMA_XMININ = 1.18e-38
_GAMMA_XINF = 3.4e38
_GAMMA_P = (-1.71618513886549492533811e0, 2.47656508055759199108314e1,
            -3.79804256470945635097577e2, 6.29331155312818442661052e2,
            8.66966202790413211295064e2, -3.14512729688483675254357e4,
            -3.61444134186911729807069e4, 6.64561438202405440627855e4)
_GAMMA_Q = (-3.08402300119738975254353e1, 3.15350626979604161529144e2,
            -1.01515636749021914166146e3, -3.10777167157231109440444e3,
            2.25381184209801510330112e4, 4.75584627752788110767815e3,
            -1.34659959864969306392456e5, -1.15132259675553483497211e5)
_GAMMA_C = (-1.910444077728e-3, 8.4171387781295e-4, -5.952379913043012e-4,
            7.93650793500350248e-4, -2.777777777777681622553e-3,
            8.333333333333333331554247e-2, 5.7083835261e-3)
_GAMMA_NMAX = 11   # INT(Y)-1 for Y<12 is at most 10


def gamma(x):
    """Complete gamma function Γ(x), faithful to module_mp_graupel.F90:GAMMA (Cody).
    Array-capable; all argument-range branches masked with jnp.where. Negative integers
    (poles) and overflow return XINF=3.4e38, as in the Fortran.
    """
    x = jnp.asarray(x, dtype=jnp.float64)

    # ── negative-argument reflection (Y = |x| for x<=0, else x) ──
    is_neg = x <= 0.0
    y = jnp.where(is_neg, -x, x)
    y1n = jnp.floor(y)              # AINT(Y) for Y>0
    resf = y - y1n                  # fractional part
    neg_nonint = is_neg & (resf != 0.0)
    neg_singular = is_neg & (resf == 0.0)   # pole at a non-positive integer
    parity = neg_nonint & (jnp.floor(y1n * 0.5) * 2.0 != y1n)   # y1n odd
    fact = jnp.where(neg_nonint, -_GAMMA_PI / jnp.sin(_GAMMA_PI * resf), 1.0)
    y = jnp.where(neg_nonint, y + 1.0, y)   # shift into the positive-arg computation

    # ── branch A: y < EPS  →  1/y ──
    res_A = jnp.where(y >= _GAMMA_XMININ, 1.0 / jnp.where(y == 0.0, 1.0, y), _GAMMA_XINF)

    # ── branch B: EPS <= y < 12  →  rational approx over (1,2) with integer reduction ──
    lt1 = y < 1.0
    n = jnp.where(lt1, 0, jnp.floor(y).astype(jnp.int32) - 1)   # 0..10
    nf = n.astype(jnp.float64)
    y_red = jnp.where(lt1, y + 1.0, y - nf)    # reduced arg in [1,2)
    z = jnp.where(lt1, y, y_red - 1.0)         # in [0,1)
    xnum = jnp.zeros_like(y)
    xden = jnp.ones_like(y)
    for i in range(8):
        xnum = (xnum + _GAMMA_P[i]) * z
        xden = xden * z + _GAMMA_Q[i]
    res_B = xnum / xden + 1.0
    # adjust for y<1: divide by the original y
    res_B = jnp.where(lt1, res_B / jnp.where(y == 0.0, 1.0, y), res_B)
    # adjust for 1<=y<12: multiply by y_red, y_red+1, ... (n factors)
    y_acc = y_red
    for i in range(_GAMMA_NMAX):
        do_mul = (~lt1) & (i < n)
        res_B = jnp.where(do_mul, res_B * y_acc, res_B)
        y_acc = jnp.where(do_mul, y_acc + 1.0, y_acc)

    # ── branch C: y >= 12  →  Stirling asymptotic (XINF beyond XBIG) ──
    ysq = y * y
    ssum = jnp.full_like(y, _GAMMA_C[6])
    for i in range(6):
        ssum = ssum / ysq + _GAMMA_C[i]
    ssum = ssum / y - y + _GAMMA_SQRTPI
    ssum = ssum + (y - 0.5) * jnp.log(jnp.where(y > 0.0, y, 1.0))
    res_C = jnp.where(y <= _GAMMA_XBIG, jnp.exp(ssum), _GAMMA_XINF)

    res = jnp.where(y < _GAMMA_EPS, res_A,
                    jnp.where(y < 12.0, res_B, res_C))

    # ── final reflection adjustments ──
    res = jnp.where(parity, -res, res)
    res = jnp.where(fact != 1.0, fact / res, res)
    res = jnp.where(neg_singular, _GAMMA_XINF, res)
    return res


# ── Warm-rain process rates (Khairoutdinov-Kogan 2000 bulk, IRAIN=0 default) ──────────
# module_mp_graupel.F90: autoconversion PRC ~2053, accretion PRA ~2174. The BULK (local)
# KK forms — distinct from the CLUBB KK scheme's upscaled PDF-integrated rates. Inputs are
# per-kg mixing ratios / number (qc=rcm, nc=Ncm, qr=rrm) + air density rho. dosb_warm_rain
# defaults .false. → IRAIN=0.
_KK_RHOW = 997.0  # density of water [kg/m^3] (module_mp_graupel.F90:495)
_KK_CONS29 = 4.0 / 3.0 * _GAMMA_PI * _KK_RHOW * (25.0e-6) ** 3  # mass of a 25um drop


def kk_warm_rain_rates(qc, nc, qr, rho, dt):
    """KK(2000) bulk warm-rain autoconversion (PRC) + accretion (PRA), with the number-rate
    companions. Returns (prc, nprc, nprc1, pra, npra). Faithful to module_mp_graupel.F90.

    Autoconversion (guard qc>=1e-6): PRC = 1350·qc^2.47·(nc/1e6·rho)^(-1.79);
      NPRC1 = PRC/CONS29 (→Nr); NPRC = min(PRC/(qc/nc), nc/dt) (→Nc).
    Accretion (guard qr>=1e-8 & qc>=1e-8): PRA = 67·(qc·qr)^1.15; NPRA = PRA/(qc/nc).
    """
    qc = jnp.asarray(qc, dtype=jnp.float64)
    nc = jnp.asarray(nc, dtype=jnp.float64)
    qr = jnp.asarray(qr, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    # safe denominators (only used where the guard holds)
    nc_safe = jnp.where(nc > 0.0, nc, 1.0)
    qc_safe = jnp.where(qc > 0.0, qc, 1.0)

    # autoconversion needs droplets: guard nc>0 (no-op where nc>0 — all real cloud points — but
    # prevents the (nc)^-1.79 → inf at degenerate qc>=1e-6 & nc=0 in-cloud-conversion edge points)
    auto_on = (qc >= 1.0e-6) & (nc > 0.0)
    prc = jnp.where(auto_on, 1350.0 * qc ** 2.47 * (nc_safe / 1.0e6 * rho) ** (-1.79), 0.0)
    nprc1 = jnp.where(auto_on, prc / _KK_CONS29, 0.0)
    nprc = jnp.where(auto_on, jnp.minimum(prc / (qc_safe / nc_safe), nc / dt), 0.0)

    accr_on = (qr >= 1.0e-8) & (qc >= 1.0e-8)
    pra = jnp.where(accr_on, 67.0 * (qc * qr) ** 1.15, 0.0)
    npra = jnp.where(accr_on, pra / (qc_safe / nc_safe), 0.0)
    return prc, nprc, nprc1, pra, npra


# ── Size-distribution slopes (gamma distribution) ────────────────────────────────────
# module_mp_graupel.F90: rain LAMR/N0RR ~1873, cloud PGAM/LAMC ~1912. The slopes feed
# nearly every process rate, fall speed, and ventilation factor. Inputs per-kg.
_M_QSMALL = 1.0e-14
_M_LAMMAXR = 1.0 / 20.0e-6
_M_LAMMINR = 1.0 / 2800.0e-6
_M_CONS26 = _GAMMA_PI / 6.0 * _KK_RHOW   # = PI/6 * RHOW


# Per-species slope constants. CONS = Γ(1+D)·(ρ·π/6) = ρ·π for D=3 (all hydrometeors here use
# D=3, spherical). Clip bounds = 1/Dmax .. 1/Dmin. (RHOW=997, RHOI=500, RHOSN=100, RHOG=400;
# DCS=125µm is the ice→snow autoconversion size → LAMMINI=1/(2·DCS+100µm)=1/350µm.)
_M_D = 3.0
_M_CONS_RAIN = _GAMMA_PI * _KK_RHOW
_M_CONS_ICE = _GAMMA_PI * 500.0
_M_CONS_SNOW = _GAMMA_PI * 100.0
_M_CONS_GRAUPEL = _GAMMA_PI * 400.0
_M_LAMMINI, _M_LAMMAXI = 1.0 / (2.0 * 125.0e-6 + 100.0e-6), 1.0 / 1.0e-6
_M_LAMMINS, _M_LAMMAXS = 1.0 / 2000.0e-6, 1.0 / 10.0e-6
_M_LAMMING, _M_LAMMAXG = 1.0 / 2000.0e-6, 1.0 / 20.0e-6


def _gamma_slope(q, n, cons, lammin, lammax, d=_M_D):
    """Generic gamma-distribution slope: LAM=(cons·n/q)^(1/d), N0=n·LAM, clipped to [lammin,lammax]
    (N0=LAM^(d+1)·q/cons when clipped). Both 0 where q<QSMALL. The shared form of the Morrison
    rain/ice/snow/graupel slope blocks (module_mp_graupel.F90:1873/1958/1985/2807)."""
    q = jnp.asarray(q, dtype=jnp.float64); n = jnp.asarray(n, dtype=jnp.float64)
    on = q >= _M_QSMALL
    q_s = jnp.where(on, q, 1.0)
    lam = (cons * n / q_s) ** (1.0 / d)
    n0 = n * lam
    lo = lam < lammin; hi = lam > lammax
    lam_c = jnp.where(lo, lammin, jnp.where(hi, lammax, lam))
    n0_c = jnp.where(lo | hi, lam_c ** (d + 1.0) * q / cons, n0)
    return jnp.where(on, lam_c, 0.0), jnp.where(on, n0_c, 0.0)


def rain_slope(qr, nr):
    """Rain slope LAMR + intercept N0RR (module_mp_graupel.F90:1873). qr=rrm, nr=Nrm."""
    return _gamma_slope(qr, nr, _M_CONS_RAIN, _M_LAMMINR, _M_LAMMAXR)


def ice_slope(qi, ni):
    """Cloud-ice slope LAMI + N0I (module_mp_graupel.F90:2807). qi=rim, ni=Nim."""
    return _gamma_slope(qi, ni, _M_CONS_ICE, _M_LAMMINI, _M_LAMMAXI)


def snow_slope(qs, ns):
    """Snow slope LAMS + N0S (module_mp_graupel.F90:1958). qs=rsm, ns=Nsm."""
    return _gamma_slope(qs, ns, _M_CONS_SNOW, _M_LAMMINS, _M_LAMMAXS)


def graupel_slope(qg, ng):
    """Graupel slope LAMG + N0G (module_mp_graupel.F90:1985). qg=rgm, ng=Ngm."""
    return _gamma_slope(qg, ng, _M_CONS_GRAUPEL, _M_LAMMING, _M_LAMMAXG)


def cloud_slope(qc, nc, rho, dofix_pgam=False, pgam_fixed=5.0):
    """Cloud gamma-distribution shape PGAM (Martin et al. 1994) + slope LAMC
    (module_mp_graupel.F90:1903-1949). PGAM=1/(0.0005714·nc_cm3+0.2714)²−1 clipped [2,10];
    LAMC=(CONS26·nc·Γ(PGAM+4)/(qc·Γ(PGAM+1)))^(1/3) clipped to [(PGAM+1)/60µm, (PGAM+1)/1µm].
    Returns (pgam, lamc); both 0 where qc<QSMALL. Inputs per-kg (qc=rcm, nc=Ncm)."""
    qc = jnp.asarray(qc, dtype=jnp.float64)
    nc = jnp.asarray(nc, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    on = qc >= _M_QSMALL
    if dofix_pgam:
        pgam = jnp.full_like(qc, pgam_fixed)
    else:
        pgam = 0.0005714 * (nc / 1.0e6 * rho) + 0.2714
        pgam = 1.0 / (pgam ** 2) - 1.0
        pgam = jnp.clip(pgam, 2.0, 10.0)
    qc_s = jnp.where(on, qc, 1.0)
    lamc = (_M_CONS26 * nc * gamma(pgam + 4.0) / (qc_s * gamma(pgam + 1.0))) ** (1.0 / 3.0)
    lammin = (pgam + 1.0) / 60.0e-6
    lammax = (pgam + 1.0) / 1.0e-6
    lamc = jnp.clip(lamc, lammin, lammax)
    return jnp.where(on, pgam, 0.0), jnp.where(on, lamc, 0.0)


# ── Rain evaporation PRE (Rutledge & Hobbs 1983) ─────────────────────────────────────
# module_mp_graupel.F90: ventilation/thermo factors ~1739-1799, EPSR+PRE ~2215-2231.
# XXLV=Lv, CPM=Cp (the moist-cp correction is commented out in the Fortran). EVS=POLYSVP.
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp as _CP, Lv as _LV, Rd as _RD, Rv as _RV, grav as _M_G,
    Ls as _M_LS, T_freeze_K as _M_TMELT,
)   # module_mp_graupel.F90:78 `use clubb_api_module` (Cp/Lv/Rd/Rv/Ls/grav/T_freeze_K via constants_clubb)
_M_EP2 = _RD / _RV
_M_AR = 841.99667          # rain fall-speed coefficient (V·D^BR), BR=0.8
_M_BR = 0.8
_M_RHOSU = 85000.0 / (_RD * _M_TMELT)   # _M_TMELT = T_freeze_K (imported); module_mp_graupel.F90:481-482
_M_F1R = 0.78
_M_F2R = 0.308


def rain_evap_rate(qr, nr, qv, T, pres, rho):
    """Rain evaporation PRE [kg/kg/s] (≤0, evaporation only where subsaturated).
    Faithful to module_mp_graupel.F90:2215-2231 (Rutledge & Hobbs ventilated diffusion).
    Inputs per-kg (qr=rrm, nr=Nrm, qv=rvm), T [K], pres [Pa], rho [kg/m^3]."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    qv = jnp.asarray(qv, dtype=jnp.float64); T = jnp.asarray(T, dtype=jnp.float64)
    pres = jnp.asarray(pres, dtype=jnp.float64); rho = jnp.asarray(rho, dtype=jnp.float64)

    cons9 = gamma(jnp.asarray(2.5 + _M_BR / 2.0))   # = Γ(2.9)
    cons34 = 2.5 + _M_BR / 2.0                            # = 2.9
    # saturation (over liquid), clipped at 0.99*pres
    evs = jnp.minimum(0.99 * pres, polysvp(T, 0))
    qvs = _M_EP2 * evs / (pres - evs)
    # ventilation / thermodynamic factors
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    arn = (_M_RHOSU / rho) ** 0.54 * _M_AR
    dv = 8.794e-5 * T ** 1.81 / pres
    sc = mu / (rho * dv)
    dqsdt = _LV * qvs / (_RV * T ** 2)
    ab = 1.0 + dqsdt * _LV / _CP
    # rain slope
    lamr, n0rr = rain_slope(qr, nr)
    on = qr >= _M_QSMALL
    lamr_s = jnp.where(lamr > 0.0, lamr, 1.0)
    epsr = jnp.where(on, 2.0 * _GAMMA_PI * n0rr * rho * dv *
                     (_M_F1R / (lamr_s * lamr_s)
                      + _M_F2R * (arn * rho / mu) ** 0.5 * sc ** (1.0 / 3.0)
                      * cons9 / (lamr_s ** cons34)), 0.0)
    pre = jnp.where(qv < qvs, jnp.minimum(epsr * (qv - qvs) / ab, 0.0), 0.0)
    return pre


# ── Ice/snow/graupel vapor deposition & sublimation (Harrington et al. 1995) ─────────
# module_mp_graupel.F90:3697-3796. Cloud ice has no ventilation correction (EPSI); snow/
# graupel are ventilated (EPSS/EPSG). The ice tail < DCS deposits as ice, the rest to snow.
# A SUM_DEP/FUDGEF limiter caps total deposition to the available (super)saturation, then
# negative rates split off as sublimation (EPRD/EPRDS/EPRDG).
# _M_LS (latent heat of sublimation) = Ls from constants_clubb (imported above; constants_clubb.F90:209)
_M_AS, _M_BS = 11.72, 0.41     # snow fall speed V=AS·D^BS
_M_AG, _M_BG = 19.3, 0.37      # graupel (default, not hail)
_M_F1S, _M_F2S = 0.86, 0.28
_M_DCS = 125.0e-6


def ice_deposition(qi, ni, qs, ns, qg, ng, qv, T, pres, rho, mnuccd, dt):
    """Ice/snow/graupel deposition+sublimation. Returns (prd, eprd, prds, eprds, prdg, eprdg).
    Faithful to module_mp_graupel.F90:3697-3796. Inputs per-kg; mnuccd = deposition-nucleation
    rate (enters the SUM_DEP supersaturation limiter)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qi, ni, qs, ns, qg, ng = a(qi), a(ni), a(qs), a(ns), a(qg), a(ng)
    qv, T, pres, rho, mnuccd = a(qv), a(T), a(pres), a(rho), a(mnuccd)
    # ice saturation (EIS capped at the liquid value EVS), thermo factor ABI
    evs = jnp.minimum(0.99 * pres, polysvp(T, 0))
    eis = jnp.minimum(jnp.minimum(0.99 * pres, polysvp(T, 1)), evs)
    qvi = _M_EP2 * eis / (pres - eis)
    abi = 1.0 + (_M_LS * qvi / (_RV * T ** 2)) * _M_LS / _CP
    # ventilation factors
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    dv = 8.794e-5 * T ** 1.81 / pres
    sc = mu / (rho * dv)
    dcorr = (_M_RHOSU / rho) ** 0.54
    asn, agn = dcorr * _M_AS, dcorr * _M_AG
    cons10 = gamma(jnp.asarray(2.5 + _M_BS / 2.0)); cons35 = 2.5 + _M_BS / 2.0
    cons11 = gamma(jnp.asarray(2.5 + _M_BG / 2.0)); cons36 = 2.5 + _M_BG / 2.0
    # slopes
    lami, n0i = ice_slope(qi, ni)
    lams, n0s = snow_slope(qs, ns)
    lamg, n0g = graupel_slope(qg, ng)
    li = jnp.where(lami > 0, lami, 1.0); ls_ = jnp.where(lams > 0, lams, 1.0)
    lg = jnp.where(lamg > 0, lamg, 1.0)
    qi_on, qs_on, qg_on = qi >= _M_QSMALL, qs >= _M_QSMALL, qg >= _M_QSMALL
    epsi = jnp.where(qi_on, 2.0 * _GAMMA_PI * n0i * rho * dv / (li * li), 0.0)
    epss = jnp.where(qs_on, 2.0 * _GAMMA_PI * n0s * rho * dv *
                     (_M_F1S / (ls_ * ls_) + _M_F2S * (asn * rho / mu) ** 0.5
                      * sc ** (1.0 / 3.0) * cons10 / (ls_ ** cons35)), 0.0)
    epsg = jnp.where(qg_on, 2.0 * _GAMMA_PI * n0g * rho * dv *
                     (_M_F1S / (lg * lg) + _M_F2S * (agn * rho / mu) ** 0.5
                      * sc ** (1.0 / 3.0) * cons11 / (lg ** cons36)), 0.0)
    depf = (qv - qvi) / abi
    dum_h = jnp.where(qi_on, 1.0 - jnp.exp(-lami * _M_DCS) * (1.0 + lami * _M_DCS), 0.0)
    prd = jnp.where(qi_on, epsi * depf * dum_h, 0.0)
    ice_tail = jnp.where(qi_on, epsi * depf * (1.0 - dum_h), 0.0)
    prds = jnp.where(qs_on, epss * depf + ice_tail, 0.0)
    prd = jnp.where(qs_on, prd, prd + ice_tail)   # if no snow, the tail goes to cloud ice
    prdg = epsg * depf
    # SUM_DEP supersaturation limiter (Reisner 2)
    dum = (qv - qvi) / dt
    fudgef = 0.9999
    sum_dep = prd + prds + mnuccd + prdg
    correct = ((dum > 0.0) & (sum_dep > dum * fudgef)) | ((dum < 0.0) & (sum_dep < dum * fudgef))
    scale = jnp.where(sum_dep != 0.0, fudgef * dum / jnp.where(sum_dep != 0.0, sum_dep, 1.0), 1.0)
    prd = jnp.where(correct, prd * scale, prd)
    prds = jnp.where(correct, prds * scale, prds)
    prdg = jnp.where(correct, prdg * scale, prdg)
    # split negatives into sublimation
    eprd = jnp.where(prd < 0.0, prd, 0.0); prd = jnp.where(prd < 0.0, 0.0, prd)
    eprds = jnp.where(prds < 0.0, prds, 0.0); prds = jnp.where(prds < 0.0, 0.0, prds)
    eprdg = jnp.where(prdg < 0.0, prdg, 0.0); prdg = jnp.where(prdg < 0.0, 0.0, prdg)
    return prd, eprd, prds, eprds, prdg, eprdg


# ── Snow collection (riming + ice accretion) + cloud-ice→snow autoconversion ─────────
# module_mp_graupel.F90: PSACWS ~3176, PRAI ~3595, PRCI ~3581. EII=0.1 (ice-ice), ECI=0.7
# (snow-cloud water). Collection rates are field-driven (no deficit); PRCI is supersat-driven.
_M_EII, _M_ECI = 0.1, 0.7
_M_ECR = 1.0               # rain collection efficiency (module_mp_graupel.F90:506)
_M_AI, _M_BI = 700.0, 1.0  # cloud-ice fall speed V=AI·D^BI (:460/:464)
_M_RHOI = 500.0


def snow_collection_rates(qc, nc, qi, ni, qs, ns, rho, dt):
    """Snow accretion of cloud water (PSACWS, riming) and cloud ice (PRAI), + number rates.
    Returns (psacws, npsacws, prai, nprai). Faithful to module_mp_graupel.F90:3176/3595 — gravitational
    collection ∝ N0S/LAMS^(BS+3). Field-driven (single timing confound)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qc, nc, qi, ni, qs, ns, rho = a(qc), a(nc), a(qi), a(ni), a(qs), a(ns), a(rho)
    lams, n0s = snow_slope(qs, ns)
    ls_ = jnp.where(lams > 0.0, lams, 1.0)
    asn = (_M_RHOSU / rho) ** 0.54 * _M_AS
    g = float(gamma(jnp.asarray(_M_BS + 3.0)))   # Γ(3.41), a constant
    cons13 = g * _GAMMA_PI / 4.0 * _M_ECI
    cons23 = _GAMMA_PI / 4.0 * _M_EII * g
    base = asn * rho * n0s / ls_ ** (_M_BS + 3.0)
    snow_on = qs >= 1.0e-8
    qc_on = snow_on & (qc >= _M_QSMALL)
    psacws = jnp.where(qc_on, cons13 * qc * base, 0.0)
    npsacws = jnp.where(qc_on, cons13 * nc * base, 0.0)
    qi_on = snow_on & (qi >= _M_QSMALL)
    prai = jnp.where(qi_on, cons23 * qi * base, 0.0)
    nprai = jnp.where(qi_on, jnp.minimum(cons23 * ni * base, ni / dt), 0.0)
    return psacws, npsacws, prai, nprai


def ice_autoconv_to_snow(qi, ni, qv, T, pres, rho, dt):
    """Cloud ice → snow autoconversion PRCI (Harrington 1995; module_mp_graupel.F90:3581).
    Occurs only where ice is supersaturated (qv≥qvi). Returns (prci, nprci). Deficit-driven."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qi, ni, qv, T, pres, rho = a(qi), a(ni), a(qv), a(T), a(pres), a(rho)
    evs = jnp.minimum(0.99 * pres, polysvp(T, 0))
    eis = jnp.minimum(jnp.minimum(0.99 * pres, polysvp(T, 1)), evs)
    qvi = _M_EP2 * eis / (pres - eis)
    abi = 1.0 + (_M_LS * qvi / (_RV * T ** 2)) * _M_LS / _CP
    dv = 8.794e-5 * T ** 1.81 / pres
    lami, n0i = ice_slope(qi, ni)
    cons21 = 4.0 / (_M_DCS * _M_RHOI)
    cons22 = _GAMMA_PI * _M_RHOI * _M_DCS ** 3 / 6.0
    on = (qi >= 1.0e-8) & (qv >= qvi)
    nprci = jnp.where(on, cons21 * (qv - qvi) * rho * n0i * jnp.exp(-lami * _M_DCS) * dv / abi, 0.0)
    prci = cons22 * nprci
    nprci = jnp.where(on, jnp.minimum(nprci, ni / dt), 0.0)
    return prci, nprci


# ── Snow self-aggregation NSAGG (Passarelli 1978) ────────────────────────────────────
# module_mp_graupel.F90:3160. Number-loss rate (snow crystals aggregate). CONS15 hard-wired ~BS=0.4.
_M_RHOSN = 100.0
_M_CONS15 = (-1108.0 * _M_EII * _GAMMA_PI ** ((1.0 - _M_BS) / 3.0)
             * _M_RHOSN ** ((-2.0 - _M_BS) / 3.0) / (4.0 * 720.0))


def snow_self_aggregation(qs, ns, rho):
    """Snow self-aggregation NSAGG (module_mp_graupel.F90:3160). Returns the (negative) snow-number
    tendency. Inputs per-kg (qs=rsm, ns=Nsm)."""
    qs = jnp.asarray(qs, dtype=jnp.float64); ns = jnp.asarray(ns, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    asn = (_M_RHOSU / rho) ** 0.54 * _M_AS
    e1 = (2.0 + _M_BS) / 3.0
    e2 = (4.0 - _M_BS) / 3.0
    nsagg = _M_CONS15 * asn * rho ** e1 * qs ** e1 * (ns * rho) ** e2 / rho
    return jnp.where(qs >= 1.0e-8, nsagg, 0.0)


# ── Deposition ice nucleation MNUCCD/NNUCCD (Cooper curve, INUC=0) ────────────────────
# module_mp_graupel.F90:3645. doarcticicenucl=.false. → INUC=0 (Rasmussen/Cooper mid-latitude).
# NNUCCD_REDUCE_COEF=1.0 default (clex9_oct14 reduces it ×100). MI0 = mass of a 10µm ice sphere.
_M_MI0 = 4.0 / 3.0 * _GAMMA_PI * _M_RHOI * (10.0e-6) ** 3


# NNUCCD reduction (module_mp_graupel.F90:90 default 1.0; set to 0.01 for clex9_oct14 ONLY,
# microphys_init_cleanup.F90:378 — reduce deposition nucleation 100× so its cloud persists).
_NNUCCD_REDUCE_COEF = 1.0


def deposition_nucleation(qv, T, ni, ns, ng, pres, rho, dt, nnuccd_reduce_coef=_NNUCCD_REDUCE_COEF):
    """Deposition ice nucleation MNUCCD (mass) + NNUCCD (number), Cooper curve.
    Faithful to module_mp_graupel.F90:3645-3682 incl. the #ifdef CLUBB `NNUCCD *= NNUCCD_REDUCE_COEF`
    (default 1.0; 0.01 for clex9_oct14). Returns (mnuccd, nnuccd). Inputs per-kg."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qv, T, ni, ns, ng, pres, rho = a(qv), a(T), a(ni), a(ns), a(ng), a(pres), a(rho)
    evs = jnp.minimum(0.99 * pres, polysvp(T, 0))
    eis = jnp.minimum(jnp.minimum(0.99 * pres, polysvp(T, 1)), evs)
    qvqvs = qv / (_M_EP2 * evs / (pres - evs))
    qvqvsi = qv / (_M_EP2 * eis / (pres - eis))
    guard = ((qvqvs >= 0.999) & (T <= 265.15)) | (qvqvsi >= 1.08)
    kc2 = 0.005 * jnp.exp(0.304 * (273.15 - T)) * 1000.0   # Cooper, L^-1 → m^-3
    kc2 = jnp.minimum(kc2, 500.0e3)
    kc2 = jnp.maximum(kc2 / rho, 0.0)                       # → kg^-1
    nnuccd = jnp.where(guard & (kc2 > ni + ns + ng), (kc2 - ni - ns - ng) / dt, 0.0)
    nnuccd = nnuccd * nnuccd_reduce_coef                    # #ifdef CLUBB reduction (F90:3682)
    return nnuccd * _M_MI0, nnuccd


# ── Immersion freezing of rain → graupel MNUCCR (Bigg 1953) ──────────────────────────
# module_mp_graupel.F90:3512. Allowed below -4 C. Uses the rain slope LAMR.
_M_AIMM, _M_BIMM = 0.66, 100.0
_M_CONS20 = 20.0 * _GAMMA_PI * _GAMMA_PI * _KK_RHOW * _M_BIMM


def rain_immersion_freezing(qr, nr, T, dt):
    """Immersion freezing of rain MNUCCR (mass) + NNUCCR (number), Bigg 1953.
    Faithful to module_mp_graupel.F90:3512-3518. Returns (mnuccr, nnuccr). Guard T<269.15 & qr>=QSMALL."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    T = jnp.asarray(T, dtype=jnp.float64)
    lamr, _ = rain_slope(qr, nr)
    lr = jnp.where(lamr > 0.0, lamr, 1.0)
    on = (T < 269.15) & (qr >= _M_QSMALL)
    ex = jnp.exp(_M_AIMM * (273.15 - T))
    mnuccr = jnp.where(on, _M_CONS20 * nr * ex / lr ** 6, 0.0)
    nnuccr = jnp.where(on, jnp.minimum(_GAMMA_PI * nr * _M_BIMM * ex / lr ** 3, nr / dt), 0.0)
    return mnuccr, nnuccr


# ── Sublimation/evaporation number tendencies NSUBI/NSUBS/NSUBR ──────────────────────
# module_mp_graupel.F90:4343-4357. The number lost ∝ the fractional mass sublimated, capped
# so it can remove at most the whole number (DUM≥−1). Driven by the negative rates EPRD/EPRDS/PRE.
def sublimation_number_rates(eprd, eprds, pre, qi, ni, qs, ns, qr, nr, dt):
    """Number tendencies from ice/snow/rain sublimation+evaporation. Returns (nsubi, nsubs, nsubr).
    NSUB = max(−1, rate·dt/q)·N/dt where rate<0 (sublimation/evap), else 0."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    def _nsub(rate, q, n):
        rate, q, n = a(rate), a(q), a(n)
        q_s = jnp.where(q != 0.0, q, 1.0)
        dum = jnp.maximum(-1.0, rate * dt / q_s)
        return jnp.where(rate < 0.0, dum * n / dt, 0.0)
    return _nsub(eprd, qi, ni), _nsub(eprds, qs, ns), _nsub(pre, qr, nr)


# ── Contact + immersion freezing of cloud water MNUCCC (Meyers/Bigg) ─────────────────
# module_mp_graupel.F90:3043-3099. Contact nuclei (NACNT, Meyers 1992) Brownian-diffuse to droplets
# (DAP); plus Bigg immersion. Uses the cloud slope (PGAM/LAMC) + CDIST1=nc/Γ(PGAM+1). The Fortran
# computes the moments in log space (exp(ln a + ln b − n·ln c)) — preserved for bit-faithfulness.
_M_RIN = 0.1e-6
_M_CONS37 = 4.0 * _GAMMA_PI * 1.38e-23 / (6.0 * _GAMMA_PI * _M_RIN)
_M_CONS38 = _GAMMA_PI ** 2 / 3.0 * _KK_RHOW
_M_CONS39 = _GAMMA_PI ** 2 / 36.0 * _KK_RHOW * _M_BIMM
_M_CONS40 = _GAMMA_PI / 6.0 * _M_BIMM


def cloud_contact_immersion_freezing(qc, nc, T, pres, rho, dt):
    """Contact + immersion freezing of cloud water MNUCCC (mass) + NNUCCC (number).
    Faithful to module_mp_graupel.F90:3043-3099. Returns (mnuccc, nnuccc). Guard qc>=QSMALL & T<269.15."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qc, nc, T, pres, rho = a(qc), a(nc), a(T), a(pres), a(rho)
    pgam, lamc = cloud_slope(qc, nc, rho)
    lc = jnp.where(lamc > 0.0, lamc, 1.0)
    cd = jnp.where(qc >= _M_QSMALL, nc / gamma(pgam + 1.0), 1.0)   # CDIST1
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    nacnt = jnp.exp(-2.80 + 0.262 * (273.15 - T)) * 1000.0
    mfp = 7.37 * T / (288.0 * 10.0 * pres) / 100.0                     # mean free path
    dap = _M_CONS37 * T * (1.0 + mfp / _M_RIN) / mu
    lncd, lnlc = jnp.log(cd), jnp.log(lc)
    # contact freezing (CONS38·DAP·NACNT · CDIST1·Γ(PGAM+5)/LAMC^4)
    mnuccc = _M_CONS38 * dap * nacnt * jnp.exp(lncd + jnp.log(gamma(pgam + 5.0)) - 4.0 * lnlc)
    nnuccc = 2.0 * _GAMMA_PI * dap * nacnt * cd * gamma(pgam + 2.0) / lc
    # immersion (Bigg)
    eb = jnp.exp(_M_AIMM * (273.15 - T))
    mnuccc = mnuccc + _M_CONS39 * jnp.exp(lncd + jnp.log(gamma(7.0 + pgam)) - 6.0 * lnlc) * eb
    nnuccc = nnuccc + _M_CONS40 * jnp.exp(lncd + jnp.log(gamma(pgam + 4.0)) - 3.0 * lnlc) * eb
    on = (qc >= _M_QSMALL) & (T < 269.15)
    return jnp.where(on, mnuccc, 0.0), jnp.where(on, jnp.minimum(nnuccc, nc / dt), 0.0)


# ── Rain self-collection NRAGG (with break-up above 300µm) ───────────────────────────
# module_mp_graupel.F90:3561. Number-loss; break-up (dum<1) when drops exceed 300µm.
def rain_self_collection(qr, nr, rho):
    """Rain self-collection/break-up NRAGG (module_mp_graupel.F90:3561). Returns the (negative)
    rain-number tendency. Guard qr>=1e-8."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lamr, _ = rain_slope(qr, nr)
    lr = jnp.where(lamr > 0.0, lamr, 1.0)
    inv = 1.0 / lr
    dum = jnp.where(inv < 300.0e-6, 1.0, 2.0 - jnp.exp(2300.0 * (inv - 300.0e-6)))
    nragg = -5.78 * dum * nr * qr * rho
    return jnp.where(qr >= 1.0e-8, nragg, 0.0)


# ── Conservation limiters (scale sink rates so they don't over-deplete each species) ──
# module_mp_graupel.F90:3850-3923 (cold branch, nov11 active terms; graupel/QMULT*/PIACR*=0).
# QC has no sources (RATIO=QC/DUM); QI/QR/QNI separate sources from sinks (RATIO=(q/dt+src)/sink).
def conserve_qc(prc, pra, mnuccc, psacws, psacwi, qc, dt):
    """Scale qc sinks (PRC/PRA/MNUCCC/PSACWS/PSACWI) if their dt-depletion exceeds qc (:3850)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    prc, pra, mnuccc, psacws, psacwi, qc = a(prc), a(pra), a(mnuccc), a(psacws), a(psacwi), a(qc)
    dum = (prc + pra + mnuccc + psacws + psacwi) * dt
    over = (dum > qc) & (qc >= _M_QSMALL)
    r = jnp.where(over, qc / jnp.where(dum > 0.0, dum, 1.0), 1.0)
    return prc * r, pra * r, mnuccc * r, psacws * r, psacwi * r


def conserve_qi(prd, eprd, mnuccc, mnuccd, psacwi, prci, prai, qi, dt):
    """Scale qi sinks (PRCI/PRAI/EPRD) if net depletion exceeds qi (:3868). Sources add ice."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    prd, eprd, mnuccc, mnuccd, psacwi, prci, prai, qi = (a(v) for v in (prd, eprd, mnuccc, mnuccd, psacwi, prci, prai, qi))
    dum = (-prd - mnuccc + prci + prai - mnuccd - eprd - psacwi) * dt
    sinks = prci + prai - eprd
    over = (dum > qi) & (qi >= _M_QSMALL)
    r = jnp.where(over, (qi / dt + prd + mnuccc + mnuccd + psacwi) / jnp.where(sinks != 0.0, sinks, 1.0), 1.0)
    return prci * r, prai * r, eprd * r


def conserve_qr(pre, prc, pra, mnuccr, qr, dt):
    """Scale qr sinks (PRE/MNUCCR) if net depletion exceeds qr (:3887). PRC/PRA are sources."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    pre, prc, pra, mnuccr, qr = a(pre), a(prc), a(pra), a(mnuccr), a(qr)
    dum = (-pre - prc + mnuccr - pra) * dt
    sinks = -pre + mnuccr
    over = (dum > qr) & (qr >= _M_QSMALL)
    r = jnp.where(over, (qr / dt + prc + pra) / jnp.where(sinks != 0.0, sinks, 1.0), 1.0)
    return pre * r, mnuccr * r


def conserve_qni(prds, psacws, prai, prci, eprds, qni, dt):
    """Scale snow sink (EPRDS) if net depletion exceeds qni (:3912, IGRAUP=0). The rest are sources."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    prds, psacws, prai, prci, eprds, qni = a(prds), a(psacws), a(prai), a(prci), a(eprds), a(qni)
    dum = (-prds - psacws - prai - prci - eprds) * dt
    over = (dum > qni) & (qni >= _M_QSMALL)
    r = jnp.where(over, (qni / dt + prds + psacws + prai + prci) / jnp.where(eprds != 0.0, -eprds, 1.0), 1.0)
    return eprds * r


# ── PCC — saturation adjustment (condense vapor above liquid saturation) ─────────────
# module_mp_graupel.F90:4018-4036. Computed AFTER the process tendencies, on the updated state
# (T+dt·T_ten, qv+dt·qv_ten). Condenses the excess vapor (≥0) or evaporates cloud (≤qc available).
def saturation_adjustment_pcc(T, qv, qc, t_ten, qv_ten, qc_ten, pres, dt):
    """Cloud condensation/evaporation rate PCC (saturation adjustment). Faithful to
    module_mp_graupel.F90:4018-4036. Returns PCC [kg/kg/s]; QV3DTEN−=PCC, QC3DTEN+=PCC, T3DTEN+=PCC·Lv/Cp."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    T, qv, qc, t_ten, qv_ten, qc_ten, pres = (a(v) for v in (T, qv, qc, t_ten, qv_ten, qc_ten, pres))
    dumt = T + dt * t_ten
    dumqv = qv + dt * qv_ten
    dum = jnp.minimum(0.99 * pres, polysvp(dumt, 0))
    dumqss = _M_EP2 * dum / (pres - dum)
    dums = dumqv - dumqss
    pcc = dums / (1.0 + _LV ** 2 * dumqss / (_CP * _RV * dumt ** 2)) / dt
    # can't evaporate more cloud than will exist
    dumqc = qc + dt * qc_ten
    pcc = jnp.where(pcc * dt + dumqc < 0.0, -dumqc / dt, pcc)
    return pcc


# ── CLUBB-Morrison subgrid (in-cloud) conversion + negative-fix ──────────────────────
# module_mp_graupel.F90:1574-1607 (in-cloud) / :1864 (neg-fix) / :4385 (×CF3D back to grid-mean).
# Within cloud (CF3D>thresh): vapor is set to liquid saturation and hydrometeors are ÷CF3D; the
# resulting in-cloud tendencies are ×CF3D back to grid-mean. Gated on CF3D > cloud_frac_thresh.
_M_CF_THRESH = 0.005


def to_in_cloud(fields, qv, qvs, cf3d):
    """Grid-mean → in-cloud: ÷CF3D for each hydrometeor in `fields` (a tuple), and set vapor=QVS,
    where CF3D>thresh (else unchanged). Returns (fields_in_cloud_tuple, qv_in_cloud)."""
    cf3d = jnp.asarray(cf3d, dtype=jnp.float64)
    on = cf3d > _M_CF_THRESH
    cf = jnp.where(on, cf3d, 1.0)
    ic = tuple(jnp.where(on, jnp.asarray(f, jnp.float64) / cf, jnp.asarray(f, jnp.float64)) for f in fields)
    qv_ic = jnp.where(on, jnp.asarray(qvs, jnp.float64), jnp.asarray(qv, jnp.float64))
    return ic, qv_ic


def tendency_to_grid_mean(tendencies, cf3d):
    """In-cloud tendencies → grid-mean: ×CF3D where CF3D>thresh (module_mp_graupel.F90:4385)."""
    cf3d = jnp.asarray(cf3d, dtype=jnp.float64)
    w = jnp.where(cf3d > _M_CF_THRESH, cf3d, 1.0)
    return tuple(jnp.asarray(t, jnp.float64) * w for t in tendencies)


def neg_fix_number(*ns):
    """Clip negative number concentrations to 0 (module_mp_graupel.F90:1864)."""
    return tuple(jnp.maximum(0.0, jnp.asarray(n, jnp.float64)) for n in ns)


# ── Rain fall speeds (mass UMR + number UNR) for sedimentation ───────────────────────
# module_mp_graupel.F90:4662-4668. ARN=(rhosu/rho)^0.54·AR; capped at 9.1·(rhosu/rho)^0.54.
def rain_fall_speed(qr, nr, rho):
    """Rain mass (UMR) and number (UNR) terminal fall speeds [m/s]. Faithful to
    module_mp_graupel.F90:4662-4668. Returns (umr, unr); 0 where qr<QSMALL."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lamr, _ = rain_slope(qr, nr)
    lr = jnp.where(lamr > 0.0, lamr, 1.0)
    dcorr = (_M_RHOSU / rho) ** 0.54
    arn = dcorr * _M_AR
    cons4 = gamma(jnp.asarray(4.0 + _M_BR)) / 6.0   # Γ(4.8)/6
    cons6 = gamma(jnp.asarray(1.0 + _M_BR))         # Γ(1.8)
    cap = 9.1 * dcorr
    on = qr >= _M_QSMALL
    umr = jnp.where(on, jnp.minimum(arn * cons4 / lr ** _M_BR, cap), 0.0)
    unr = jnp.where(on, jnp.minimum(arn * cons6 / lr ** _M_BR, cap), 0.0)
    return umr, unr


# ── Sedimentation: CFL sub-stepped upwind flux divergence ────────────────────────────
# module_mp_graupel.F90:4712-4834. The CLUBB↔M2005 interface flips the vertical index
# (microphysics.F90:1944 `m = nz-k`), so M2005's KTE↔surface maps to the JAX grid as
# index 0 = surface, nzt-1 = top. Equivalent here with "above k" = k+1, top cell (-1) has
# no inflow, rain exits the surface (index 0). Mass/centroid/non-negativity verified.
def _fall_speed_propagate(fr):
    """FR(K)=FR(K+1) where FR(K)<1e-10 (module_mp_graupel.F90:4714) — fill the fall speed from the
    rain-bearing cell above into clear cells below so sedimenting rain does not stall mid-substep.
    Operates on the VERTICAL (last) axis, so it is correct for (ngrdcol, nzt) and bare (nzt,) inputs."""
    import jax.lax as lax
    frv = jnp.moveaxis(fr, -1, 0)            # vertical → axis 0 for the scan (no-op when fr is 1-D)
    def _fill(carry, fr_k):
        new = jnp.where(fr_k < 1.0e-10, carry, fr_k)
        return new, new
    _, out = lax.scan(_fill, frv[-1], frv[::-1])   # scan top (last index) → bottom (index 0)
    return jnp.moveaxis(out[::-1], 0, -1)


def _sediment(dum0, fr, dzq, dt, nstep):
    """One species/moment through the NSTEP upwind flux-divergence loop. `dum0` = ρ·field.
    Returns ρ·field after sedimentation. faltnd[k]=(F[k+1]-F[k])/dz interior, -F[top]/dz at the top.
    Operates on the VERTICAL (last) axis — correct for (ngrdcol, nzt) and bare (nzt,) inputs (the
    older axis-0 form silently destroyed mass for a 2-D (ngrdcol, nzt) column, Iter233)."""
    import jax.lax as lax
    def _body(n, d):
        falout = fr * d
        interior = (falout[..., 1:] - falout[..., :-1]) / dzq[..., :-1]
        top = (-falout[..., -1] / dzq[..., -1])[..., None]
        tend = jnp.concatenate([interior, top], axis=-1)
        return d + tend * dt / nstep
    return lax.fori_loop(0, nstep, _body, dum0)


def rain_sedimentation(qr, nr, rho, dzq, dt):
    """Rain mass + number sedimentation tendencies (qrsten, nrsten) [per s], faithful to
    module_mp_graupel.F90:4712-4834. Inputs are the POST-process-rate fields; `dzq`=layer
    thickness [m], index 0=surface. NSTEP is the shared CFL count over both moments (RGVM)."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64); dzq = jnp.asarray(dzq, dtype=jnp.float64)
    umr, unr = rain_fall_speed(qr, nr, rho)
    fmr = _fall_speed_propagate(umr)
    fnr = _fall_speed_propagate(unr)
    rgvm = jnp.maximum(jnp.max(fmr * dt / dzq), jnp.max(fnr * dt / dzq))
    nstep = int(jnp.maximum(jnp.floor(rgvm + 1.0), 1.0))
    qr_f = _sediment(qr * rho, fmr, dzq, dt, nstep) / rho
    nr_f = _sediment(nr * rho, fnr, dzq, dt, nstep) / rho
    return (qr_f - qr) / dt, (nr_f - nr) / dt


# ── M2005 cold-branch (T<273.15) tendency assembly ───────────────────────────────────
# module_mp_graupel.F90:3801-4007 — the "CONSERVATION OF WATER" over-depletion limiters
# (QC/QI/QR/QNI/QG) followed by the mass+number tendency assembly. Composes the ported
# process rates into the d/dt(field) the driver integrates. Faithful to IGRAUP=0 (graupel
# included) and IGRAUP=1 (no graupel → rain freezing feeds snow). PCC is added separately.
def m2005_cold_tendencies(rates, qc3d, qi3d, qr3d, qni3d, qg3d, dt,
                          xxlv, xxls, xlf, cpm, pcc=0.0, igraup=0):
    """Assemble the cold-branch (T<273.15) Morrison tendencies from the process-rate dict `rates`.

    `rates` keys are the Fortran rate names (PRC, PRA, PRD, PRE, MNUCCC, PSACWS, ...); any missing
    key defaults to 0 (a case without graupel/rime-splintering simply omits those). Returns a dict of
    the 12 tendencies {qv,T,qc,qi,qr,qni,qg,nc,ni,nr,ns,ng}_ten. `xxlv/xxls/xlf` = latent heats of
    vaporization/sublimation/fusion, `cpm` = moist specific heat. Conservation contract: with pcc=0,
    Σ mass tendencies (qv+qc+qi+qr+qni+qg) = 0 exactly (each rate is a +source/−sink pair)."""
    qs = _M_QSMALL

    def g(name):
        return jnp.asarray(rates.get(name, 0.0), dtype=jnp.float64)

    # gather every rate the cold branch references
    PRC, PRA, PRD, PRE = g('PRC'), g('PRA'), g('PRD'), g('PRE')
    PRDS, PRDG, EPRD, EPRDS, EPRDG = g('PRDS'), g('PRDG'), g('EPRD'), g('EPRDS'), g('EPRDG')
    PRCI, PRAI, PRACI, PRACIS, PRACS = g('PRCI'), g('PRAI'), g('PRACI'), g('PRACIS'), g('PRACS')
    MNUCCC, MNUCCR, MNUCCD = g('MNUCCC'), g('MNUCCR'), g('MNUCCD')
    PSACWS, PSACWI, PSACWG, PGSACW = g('PSACWS'), g('PSACWI'), g('PSACWG'), g('PGSACW')
    QMULTS, QMULTG, QMULTR, QMULTRG = g('QMULTS'), g('QMULTG'), g('QMULTR'), g('QMULTRG')
    PIACR, PIACRS, PGRACS, PRACG, PSACR = g('PIACR'), g('PIACRS'), g('PGRACS'), g('PRACG'), g('PSACR')
    qc3d = jnp.asarray(qc3d, dtype=jnp.float64); qi3d = jnp.asarray(qi3d, dtype=jnp.float64)
    qr3d = jnp.asarray(qr3d, dtype=jnp.float64); qni3d = jnp.asarray(qni3d, dtype=jnp.float64)
    qg3d = jnp.asarray(qg3d, dtype=jnp.float64)

    def _safe_ratio(num, den):
        return num / jnp.where(den == 0.0, 1.0, den)

    # ── SENSITIVITY: NO GRAUPEL (:3826) — zero graupel mass rates, fold PIACR into PIACRS ──
    if igraup == 1:
        _z = jnp.zeros_like(qr3d)
        PIACRS = PIACRS + PIACR
        PIACR = _z
        PRACG = PSACR = PSACWG = PGSACW = PGRACS = PRDG = EPRDG = _z

    # ── CONSERVATION OF QC (:3850) ──
    dum = (PRC + PRA + MNUCCC + PSACWS + PSACWI + QMULTS + QMULTG + PSACWG + PGSACW) * dt
    m = (dum > qc3d) & (qc3d >= qs)
    ratio = jnp.where(m, _safe_ratio(qc3d, dum), 1.0)
    PRC, PRA, MNUCCC = PRC * ratio, PRA * ratio, MNUCCC * ratio
    PSACWS, PSACWI = PSACWS * ratio, PSACWI * ratio
    QMULTS, QMULTG, PSACWG, PGSACW = QMULTS * ratio, QMULTG * ratio, PSACWG * ratio, PGSACW * ratio

    # ── CONSERVATION OF QI (:3868) ──
    dum = (-PRD - MNUCCC + PRCI + PRAI - QMULTS - QMULTG - QMULTR - QMULTRG
           - MNUCCD + PRACI + PRACIS - EPRD - PSACWI) * dt
    m = (dum > qi3d) & (qi3d >= qs)
    den = PRCI + PRAI + PRACI + PRACIS - EPRD
    ratio = jnp.where(m, _safe_ratio(qi3d / dt + PRD + MNUCCC + QMULTS + QMULTG + QMULTR
                                     + QMULTRG + MNUCCD + PSACWI, den), 1.0)
    PRCI, PRAI, PRACI = PRCI * ratio, PRAI * ratio, PRACI * ratio
    PRACIS, EPRD = PRACIS * ratio, EPRD * ratio

    # ── CONSERVATION OF QR (:3887) ──
    dum = ((PRACS - PRE) + (QMULTR + QMULTRG - PRC) + (MNUCCR - PRA)
           + PIACR + PIACRS + PGRACS + PRACG) * dt
    m = (dum > qr3d) & (qr3d >= qs)
    den = -PRE + QMULTR + QMULTRG + PRACS + MNUCCR + PIACR + PIACRS + PGRACS + PRACG
    ratio = jnp.where(m, _safe_ratio(qr3d / dt + PRC + PRA, den), 1.0)
    PRE, PRACS, QMULTR, QMULTRG = PRE * ratio, PRACS * ratio, QMULTR * ratio, QMULTRG * ratio
    MNUCCR, PIACR, PIACRS = MNUCCR * ratio, PIACR * ratio, PIACRS * ratio
    PGRACS, PRACG = PGRACS * ratio, PRACG * ratio

    # ── CONSERVATION OF QNI (:3912, IGRAUP branch) ──
    extra = MNUCCR if igraup == 1 else 0.0
    dum = (-PRDS - PSACWS - PRAI - PRCI - PRACS - EPRDS + PSACR - PIACRS - PRACIS - extra) * dt
    m = (dum > qni3d) & (qni3d >= qs)
    den = -EPRDS + PSACR
    num = (qni3d / dt + PRDS + PSACWS + PRAI + PRCI + PRACS + PIACRS + PRACIS + extra)
    ratio = jnp.where(m, _safe_ratio(num, den), 1.0)
    EPRDS, PSACR = EPRDS * ratio, PSACR * ratio

    # ── CONSERVATION OF QG (:3941) ──
    dum = (-PSACWG - PRACG - PGSACW - PGRACS - PRDG - MNUCCR - EPRDG - PIACR - PRACI - PSACR) * dt
    m = (dum > qg3d) & (qg3d >= qs)
    num = (qg3d / dt + PSACWG + PRACG + PGSACW + PGRACS + PRDG + MNUCCR + PSACR + PIACR + PRACI)
    ratio = jnp.where(m, _safe_ratio(num, -EPRDG), 1.0)
    EPRDG = EPRDG * ratio

    # ── TENDENCIES (:3963) ──
    qv_ten = -PRE - PRD - PRDS - MNUCCD - EPRD - EPRDS - PRDG - EPRDG
    T_ten = (PRE * xxlv
             + (PRD + PRDS + MNUCCD + EPRD + EPRDS + PRDG + EPRDG) * xxls
             + (PSACWS + PSACWI + MNUCCC + MNUCCR + QMULTS + QMULTG + QMULTR + QMULTRG + PRACS
                + PSACWG + PRACG + PGSACW + PGRACS + PIACR + PIACRS) * xlf) / cpm
    qc_ten = (-PRA - PRC - MNUCCC + jnp.asarray(pcc, dtype=jnp.float64)
              - PSACWS - PSACWI - QMULTS - QMULTG - PSACWG - PGSACW)
    qi_ten = (PRD + EPRD + PSACWI + MNUCCC - PRCI - PRAI + QMULTS + QMULTG + QMULTR + QMULTRG
              + MNUCCD - PRACI - PRACIS)
    qr_ten = (PRE + PRA + PRC - PRACS - MNUCCR - QMULTR - QMULTRG - PIACR - PIACRS - PRACG - PGRACS)
    if igraup == 1:
        qni_ten = (PRAI + PSACWS + PRDS + PRACS + PRCI + EPRDS - PSACR + PIACRS + PRACIS + MNUCCR)
        ns_ten = (g('NSAGG') + g('NPRCI') - g('NSCNG') - g('NGRACS') + g('NIACRS') + g('NNUCCR'))
        qg_ten = jnp.zeros_like(qr_ten)
        ng_ten = jnp.zeros_like(qr_ten)
    else:
        qni_ten = (PRAI + PSACWS + PRDS + PRACS + PRCI + EPRDS - PSACR + PIACRS + PRACIS)
        ns_ten = (g('NSAGG') + g('NPRCI') - g('NSCNG') - g('NGRACS') + g('NIACRS'))
        qg_ten = (PRACG + PSACWG + PGSACW + PGRACS + PRDG + EPRDG + MNUCCR + PIACR + PRACI + PSACR)
        ng_ten = (g('NSCNG') + g('NGRACS') + g('NNUCCR') + g('NIACR'))
    nc_ten = -g('NNUCCC') - g('NPSACWS') - g('NPRA') - g('NPRC') - g('NPSACWI') - g('NPSACWG')
    ni_ten = (g('NNUCCC') - g('NPRCI') - g('NPRAI') + g('NMULTS') + g('NMULTG') + g('NMULTR')
              + g('NMULTRG') + g('NNUCCD') - g('NIACR') - g('NIACRS'))
    nr_ten = (g('NPRC1') - g('NPRACS') - g('NNUCCR') + g('NRAGG') - g('NIACR') - g('NIACRS')
              - g('NPRACG') - g('NGRACS'))
    return {'qv': qv_ten, 'T': T_ten, 'qc': qc_ten, 'qi': qi_ten, 'qr': qr_ten,
            'qni': qni_ten, 'qg': qg_ten, 'nc': nc_ten, 'ni': ni_ten, 'nr': nr_ten,
            'ns': ns_ten, 'ng': ng_ten}


# ── M2005 warm-branch (T>=273.15) tendency assembly ──────────────────────────────────
# module_mp_graupel.F90:2318-2440 — conservation (QC/snow/graupel/QR) + the mass/number
# tendency assembly + the number melting/sublimation sub-calcs. The warm branch has no ice
# *growth*: only melting (PSMLT/PGMLT), evaporation of melting precip (EVPMS/EVPMG), warm-rain
# (PRC/PRA/PRE). In CLUBB PRACS=PRACG=0 here (the driver asserts it). PCC added separately.
def m2005_warm_tendencies(rates, qc3d, qr3d, qni3d, qg3d, nr3d, ns3d, ng3d, dt,
                          xxlv, xxls, xlf, cpm, pcc=0.0):
    """Assemble the warm-branch (T>=273.15) Morrison tendencies from the process-rate dict `rates`
    (keys PRC/PRA/PRE/PSMLT/PGMLT/EVPMS/EVPMG/PRACS/PRACG + numbers NPRA/NPRC/NPRC1/NRAGG/NPRACG).
    Returns {qv,T,qc,qi,qr,qni,qg,nc,ni,nr,ns,ng}_ten (qi/ni ≡ 0 — no cloud ice above freezing).
    Conservation contract: with pcc=0, Σ mass tendencies (qv+qc+qi+qr+qni+qg) = 0 exactly."""
    qs = _M_QSMALL

    def g(name):
        return jnp.asarray(rates.get(name, 0.0), dtype=jnp.float64)

    PRC, PRA, PRE = g('PRC'), g('PRA'), g('PRE')
    PSMLT, PGMLT, EVPMS, EVPMG = g('PSMLT'), g('PGMLT'), g('EVPMS'), g('EVPMG')
    PRACS, PRACG = g('PRACS'), g('PRACG')
    qc3d = jnp.asarray(qc3d, dtype=jnp.float64); qr3d = jnp.asarray(qr3d, dtype=jnp.float64)
    qni3d = jnp.asarray(qni3d, dtype=jnp.float64); qg3d = jnp.asarray(qg3d, dtype=jnp.float64)
    nr3d = jnp.asarray(nr3d, dtype=jnp.float64); ns3d = jnp.asarray(ns3d, dtype=jnp.float64)
    ng3d = jnp.asarray(ng3d, dtype=jnp.float64)

    def _safe_ratio(num, den):
        return num / jnp.where(den == 0.0, 1.0, den)

    # ── CONSERVATION OF QC (:2319) ──
    dum = (PRC + PRA) * dt
    m = (dum > qc3d) & (qc3d >= qs)
    ratio = jnp.where(m, _safe_ratio(qc3d, dum), 1.0)
    PRC, PRA = PRC * ratio, PRA * ratio
    # ── CONSERVATION OF SNOW (:2331) ──
    dum = (-PSMLT - EVPMS + PRACS) * dt
    m = (dum > qni3d) & (qni3d >= qs)
    ratio = jnp.where(m, _safe_ratio(qni3d, dum), 1.0)
    PSMLT, EVPMS, PRACS = PSMLT * ratio, EVPMS * ratio, PRACS * ratio
    # ── CONSERVATION OF GRAUPEL (:2346) ──
    dum = (-PGMLT - EVPMG + PRACG) * dt
    m = (dum > qg3d) & (qg3d >= qs)
    ratio = jnp.where(m, _safe_ratio(qg3d, dum), 1.0)
    PGMLT, EVPMG, PRACG = PGMLT * ratio, EVPMG * ratio, PRACG * ratio
    # ── CONSERVATION OF QR (:2361) ── (only PRE scaled; PRE is negative)
    dum = (-PRACS - PRACG - PRE - PRA - PRC + PSMLT + PGMLT) * dt
    m = (dum > qr3d) & (qr3d >= qs)
    ratio = jnp.where(m, _safe_ratio(qr3d / dt + PRACS + PRACG + PRA + PRC - PSMLT - PGMLT, -PRE), 1.0)
    PRE = PRE * ratio

    # ── MASS TENDENCIES (:2393) ──
    qv_ten = -PRE - EVPMS - EVPMG
    T_ten = (PRE * xxlv + (EVPMS + EVPMG) * xxls
             + (PSMLT + PGMLT - PRACS - PRACG) * xlf) / cpm
    qc_ten = -PRA - PRC + jnp.asarray(pcc, dtype=jnp.float64)
    qr_ten = PRE + PRA + PRC - PSMLT - PGMLT + PRACS + PRACG
    qni_ten = PSMLT + EVPMS - PRACS
    qg_ten = PGMLT + EVPMG - PRACG

    # ── NUMBER MELTING/SUBLIMATION sub-calcs (:2408), then number tendencies ──
    def _melt_num(rate_sum, q, n):
        d = jnp.maximum(-1.0, rate_sum * dt / jnp.where(q != 0.0, q, 1.0))
        return jnp.where(rate_sum < 0.0, d * n / dt, 0.0)
    NSUBR = _melt_num(PRE, qr3d, nr3d)
    NSMLTS = _melt_num(EVPMS + PSMLT, qni3d, ns3d)
    NSMLTR = _melt_num(PSMLT, qni3d, ns3d)
    NGMLTG = _melt_num(EVPMG + PGMLT, qg3d, ng3d)
    NGMLTR = _melt_num(PGMLT, qg3d, ng3d)
    nc_ten = -g('NPRA') - g('NPRC')
    nr_ten = g('NPRC1') + g('NRAGG') - g('NPRACG') + NSUBR - NSMLTR - NGMLTR
    ns_ten = NSMLTS
    ng_ten = NGMLTG
    z = jnp.zeros_like(qr_ten)
    return {'qv': qv_ten, 'T': T_ten, 'qc': qc_ten, 'qi': z, 'qr': qr_ten,
            'qni': qni_ten, 'qg': qg_ten, 'nc': nc_ten, 'ni': z, 'nr': nr_ten,
            'ns': ns_ten, 'ng': ng_ten}


# ── M2005 single-column step: per-level warm/cold select + PCC saturation adjustment ─
# The Fortran M2005MICRO_GRAUPEL runs, per level K, the warm (T>=273.15) or cold (T<273.15)
# branch (each: rates → conservation → tendency assembly), THEN the ISATADJ=0 PCC block
# (:4015 cold / :2451 warm), which recomputes PCC from the post-assembly state and feeds it
# back into qv/T/qc. In the CLUBB build XXLV=Lv, XXLS=Ls, CPM=Cp are constants.
def m2005_step_tendencies(T3d, qv3d, qc3d, qi3d, qr3d, qni3d, qg3d,
                          nc3d, ni3d, nr3d, ns3d, ng3d, rates, pres, dt, igraup=0):
    """Assemble the full M2005 single-column tendency dict: select the warm/cold branch per level by
    T>=273.15, then apply the PCC saturation adjustment. `rates` is the process-rate dict (the rate
    functions produce 0 in the branch that doesn't apply — ice rates vanish above freezing, melting
    below). Returns {qv,T,qc,qi,qr,qni,qg,nc,ni,nr,ns,ng}_ten ready for field integration (q+=ten·dt).
    Water-conservation contract preserved: Σ mass tendencies = 0 (PCC moves only qv↔qc)."""
    xxlv, xxls, cpm = _LV, _M_LS, _CP
    xlf = _M_LS - _LV
    cold = m2005_cold_tendencies(rates, qc3d, qi3d, qr3d, qni3d, qg3d, dt,
                                 xxlv, xxls, xlf, cpm, pcc=0.0, igraup=igraup)
    warm = m2005_warm_tendencies(rates, qc3d, qr3d, qni3d, qg3d, nr3d, ns3d, ng3d, dt,
                                 xxlv, xxls, xlf, cpm, pcc=0.0)
    is_warm = jnp.asarray(T3d, dtype=jnp.float64) >= 273.15
    ten = {k: jnp.where(is_warm, warm[k], cold[k]) for k in cold}
    # ── PCC saturation adjustment (post-assembly), feeds back into qv/T/qc ──
    pcc = saturation_adjustment_pcc(T3d, qv3d, qc3d, ten['T'], ten['qv'], ten['qc'], pres, dt)
    ten['qv'] = ten['qv'] - pcc
    ten['T'] = ten['T'] + pcc * xxlv / cpm
    ten['qc'] = ten['qc'] + pcc
    return ten


# ── M2005 rate orchestration: compose the ported process rates into the rates dict ───
# Calls the ported rate functions in the Fortran dependency order (deposition nucleation →
# ice deposition → evap → warm rain → collection → numbers) and assembles the dict that
# m2005_step_tendencies / m2005_cold_tendencies consume. The rate functions each gate on
# their own QSMALL thresholds and compute their slopes internally.
#
# NOT-YET-PORTED rates (default to 0/absent — the rarer ice-ice & graupel collection that
# nov11_altocu's thin altocumulus barely exercises): PSACWI, QMULT{S,G,R,RG}, PIACR(S),
# PRACI(S), PRACS, PSACR, PGRACS, PRACG, PSACWG, PGSACW, the graupel block, and the warm
# melting rates PSMLT/PGMLT/EVPMS/EVPMG. Each added here as it is ported.
def compute_m2005_rates(qc, nc, qr, nr, qi, ni, qs, ns, qg, ng, qv, T, pres, rho, dt):
    """Assemble the M2005 process-rate dict from the ported rate functions, given the IN-CLOUD fields.
    Returns {RATE_NAME: array} for every ported rate; consumers default missing keys to 0. The values
    are the cold-branch (ice) rates — the dominant set for a cold case like nov11_altocu."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qc, nc, qr, nr = a(qc), a(nc), a(qr), a(nr)
    qi, ni, qs, ns, qg, ng = a(qi), a(ni), a(qs), a(ns), a(qg), a(ng)
    qv, T, pres, rho = a(qv), a(T), a(pres), a(rho)
    r = {}
    # deposition nucleation (Cooper) FIRST — ice_deposition needs MNUCCD
    r['MNUCCD'], r['NNUCCD'] = deposition_nucleation(qv, T, ni, ns, ng, pres, rho, dt)
    # ice/snow/graupel vapor deposition+sublimation (Harrington, with the SUM_DEP limiter)
    r['PRD'], r['EPRD'], r['PRDS'], r['EPRDS'], r['PRDG'], r['EPRDG'] = ice_deposition(
        qi, ni, qs, ns, qg, ng, qv, T, pres, rho, r['MNUCCD'], dt)
    # rain evaporation (Rutledge-Hobbs)
    r['PRE'] = rain_evap_rate(qr, nr, qv, T, pres, rho)
    # warm rain: autoconversion + accretion + number companions (bulk KK)
    r['PRC'], r['NPRC'], r['NPRC1'], r['PRA'], r['NPRA'] = kk_warm_rain_rates(qc, nc, qr, rho, dt)
    # snow rimes cloud water (PSACWS) + accretes ice (PRAI) + number companions
    r['PSACWS'], r['NPSACWS'], r['PRAI'], r['NPRAI'] = snow_collection_rates(qc, nc, qi, ni, qs, ns, rho, dt)
    # cloud ice collecting droplets (PSACWI)
    r['PSACWI'], r['NPSACWI'] = cloud_ice_collect_droplets(qc, nc, qi, ni, rho)
    # rain-ice collision → snow (PIACRS/PRACIS, the qr<0.1 g/kg branch) + number
    r['PIACRS'], r['PRACIS'], r['NIACRS'] = rain_ice_collision_snow(qr, nr, qi, ni, T, rho, dt)
    # rain collects snow (PRACS, Ikawa-Saito) — WARM-branch rate, gated to T>=TMELT (0 at cold levels)
    r['PRACS'] = jnp.where(T >= _M_TMELT, rain_accrete_snow(qr, nr, qs, ns, rho), 0.0)
    # cloud ice → snow autoconversion
    r['PRCI'], r['NPRCI'] = ice_autoconv_to_snow(qi, ni, qv, T, pres, rho, dt)
    # snow self-aggregation (number)
    r['NSAGG'] = snow_self_aggregation(qs, ns, rho)
    # rain immersion freezing (Bigg) → graupel mass + number
    r['MNUCCR'], r['NNUCCR'] = rain_immersion_freezing(qr, nr, T, dt)
    # cloud-water contact + immersion freezing → ice
    r['MNUCCC'], r['NNUCCC'] = cloud_contact_immersion_freezing(qc, nc, T, pres, rho, dt)
    # rain self-collection (number)
    r['NRAGG'] = rain_self_collection(qr, nr, rho)
    # sublimation/evaporation number losses (need EPRD/EPRDS/PRE)
    r['NSUBI'], r['NSUBS'], r['NSUBR'] = sublimation_number_rates(
        r['EPRD'], r['EPRDS'], r['PRE'], qi, ni, qs, ns, qr, nr, dt)
    return r


# ── Minor ice-ice / rain-snow collection rates (cold + warm branches) ────────────────
# The remaining nonzero collection rates for nov11_altocu (a mid-level altocumulus with cold
# cloud + warm sub-cloud): PSACWI (cold), PIACRS/PRACIS (cold), PRACS (warm).
def cloud_ice_collect_droplets(qc, nc, qi, ni, rho):
    """PSACWI / NPSACWI — cloud ice collecting droplets (Thompson 2004 size-dependent efficiency),
    module_mp_graupel.F90:3206-3217. Gate: qi>=1e-8 & qc>=QSMALL & 1/LAMI>=100µm."""
    qc = jnp.asarray(qc, dtype=jnp.float64); nc = jnp.asarray(nc, dtype=jnp.float64)
    qi = jnp.asarray(qi, dtype=jnp.float64); ni = jnp.asarray(ni, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lami, n0i = ice_slope(qi, ni)
    li = jnp.where(lami > 0.0, lami, 1.0)
    ain = (_M_RHOSU / rho) ** 0.35 * _M_AI
    cons16 = gamma(jnp.asarray(_M_BI + 3.0)) * _GAMMA_PI / 4.0 * _M_ECI  # Γ(4)·π/4·ECI
    base = cons16 * ain * rho * n0i / li ** (_M_BI + 3.0)
    on = (qi >= 1.0e-8) & (qc >= _M_QSMALL) & (1.0 / li >= 100.0e-6)
    return jnp.where(on, base * qc, 0.0), jnp.where(on, base * nc, 0.0)


def rain_ice_collision_snow(qr, nr, qi, ni, T, rho, dt):
    """PIACRS / PRACIS / NIACRS — rain-ice collision producing SNOW (the qr<0.1 g/kg branch; ≥0.1 g/kg
    routes to graupel as PIACR/PRACI, which nov11 doesn't reach), module_mp_graupel.F90:3622-3633
    (Reisner 1998). Gate: qr>=1e-8 & qi>=1e-8 & T<=TMELT & qr<0.1e-3. Returns (PIACRS, PRACIS, NIACRS)."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    qi = jnp.asarray(qi, dtype=jnp.float64); ni = jnp.asarray(ni, dtype=jnp.float64)
    T = jnp.asarray(T, dtype=jnp.float64); rho = jnp.asarray(rho, dtype=jnp.float64)
    lamr, n0rr = rain_slope(qr, nr)
    lr = jnp.where(lamr > 0.0, lamr, 1.0)
    arn = (_M_RHOSU / rho) ** 0.54 * _M_AR
    cons24 = _GAMMA_PI / 4.0 * _M_ECR * gamma(jnp.asarray(_M_BR + 3.0))          # π/4·Γ(3.8)
    cons25 = _GAMMA_PI ** 2 / 24.0 * _KK_RHOW * _M_ECR * gamma(jnp.asarray(_M_BR + 6.0))  # π²/24·ρw·Γ(6.8)
    niacrs = cons24 * ni * n0rr * arn / lr ** (_M_BR + 3.0) * rho
    piacrs = cons25 * ni * n0rr * arn / lr ** (_M_BR + 3.0) / lr ** 3 * rho
    pracis = cons24 * qi * n0rr * arn / lr ** (_M_BR + 3.0) * rho
    niacrs = jnp.minimum(jnp.minimum(niacrs, nr / dt), ni / dt)
    on = (qr >= 1.0e-8) & (qi >= 1.0e-8) & (T <= _M_TMELT) & (qr < 0.1e-3)
    return (jnp.where(on, piacrs, 0.0), jnp.where(on, pracis, 0.0), jnp.where(on, niacrs, 0.0))


def rain_accrete_snow(qr, nr, qs, ns, rho):
    """PRACS — collection of snow by rain (Ikawa-Saito 1991), module_mp_graupel.F90:2088-2108.
    Gate: qr>=1e-8 & qni>=1e-8. Uses the mass-weighted fall speeds UMS/UMR (capped). NPRACS is
    no longer subtracted (fix 053011) so only the mass rate is returned."""
    qr = jnp.asarray(qr, dtype=jnp.float64); nr = jnp.asarray(nr, dtype=jnp.float64)
    qs = jnp.asarray(qs, dtype=jnp.float64); ns = jnp.asarray(ns, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lamr, n0rr = rain_slope(qr, nr)
    lams, n0s = snow_slope(qs, ns)
    lr = jnp.where(lamr > 0.0, lamr, 1.0); ls = jnp.where(lams > 0.0, lams, 1.0)
    dcorr = (_M_RHOSU / rho) ** 0.54
    arn, asn = dcorr * _M_AR, dcorr * _M_AS
    cons3 = gamma(jnp.asarray(4.0 + _M_BS)) / 6.0   # Γ(4.41)/6
    cons4 = gamma(jnp.asarray(4.0 + _M_BR)) / 6.0   # Γ(4.8)/6
    cons31 = _GAMMA_PI ** 2 * _M_ECR * _M_RHOSN
    ums = jnp.minimum(asn * cons3 / ls ** _M_BS, 1.2 * dcorr)
    umr = jnp.minimum(arn * cons4 / lr ** _M_BR, 9.1 * dcorr)
    pracs = cons31 * (((1.2 * umr - 0.95 * ums) ** 2 + 0.08 * ums * umr) ** 0.5 * rho
                      * n0rr * n0s / ls ** 3
                      * (5.0 / (ls ** 3 * lr) + 2.0 / (ls ** 2 * lr ** 2) + 0.5 / (ls * lr ** 3)))
    on = (qr >= 1.0e-8) & (qs >= 1.0e-8)
    return jnp.where(on, pracs, 0.0)


# ── M2005 full single-column driver: grid-mean fields → grid-mean tendencies ─────────
# Composes the whole micro step: saturation → in-cloud conversion (÷CF3D) → rate
# orchestration → per-level warm/cold tendency assembly + PCC → grid-mean conversion
# (×CF3D). This is the callable the CLUBB loop wires in (mirrors the M2005MICRO_GRAUPEL
# entry/exit, module_mp_graupel.F90:1574/4385). Sedimentation is applied separately.
def m2005_driver(qc, nc, qr, nr, qi, ni, qs, ns, qg, ng, qv, T, pres, rho, cf3d, dt, igraup=0):
    """Full M2005 single-column microphysics step. Inputs are GRID-MEAN fields; returns the GRID-MEAN
    tendency dict {qv,T,qc,qi,qr,qni,qg,nc,ni,nr,ns,ng}_ten. `cf3d` = SGS cloud fraction (the M2005
    in-cloud weighting; ÷CF3D at entry, ×CF3D at exit where CF3D>cloud_frac_thresh). Water-conserving
    (Σ mass tendencies = 0 up to PCC, which moves only qv↔qc). Note `qni`≡snow."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qc, nc, qr, nr = a(qc), a(nc), a(qr), a(nr)
    qi, ni, qs, ns, qg, ng = a(qi), a(ni), a(qs), a(ns), a(qg), a(ng)
    qv, T, pres, rho = a(qv), a(T), a(pres), a(rho)
    # saturation mixing ratio w.r.t. liquid (for the in-cloud vapor set), low-pressure-capped
    evs = jnp.minimum(0.99 * pres, polysvp(T, 0))
    qvs = _M_EP2 * evs / (pres - evs)
    # grid-mean → in-cloud (÷CF3D, qv→QVS where CF3D>thresh)
    (qc_i, nc_i, qr_i, nr_i, qi_i, ni_i, qs_i, ns_i, qg_i, ng_i), qv_i = to_in_cloud(
        (qc, nc, qr, nr, qi, ni, qs, ns, qg, ng), qv, qvs, cf3d)
    # pre-rate slope clamps (F90:1881-2002): reset the in-cloud numbers so every gamma slope is in
    # bounds BEFORE the rates use them; the reset is an effective number tendency folded into ten_i.
    nc_c, nr_c, ni_c, ns_c, ng_c = _size_clamp_numbers(
        qc_i, nc_i, qr_i, nr_i, qi_i, ni_i, qs_i, ns_i, qg_i, ng_i, rho)
    sf = {'nc': (nc_c - nc_i) / dt, 'nr': (nr_c - nr_i) / dt, 'ni': (ni_c - ni_i) / dt,
          'ns': (ns_c - ns_i) / dt, 'ng': (ng_c - ng_i) / dt}
    nc_i, nr_i, ni_i, ns_i, ng_i = nc_c, nr_c, ni_c, ns_c, ng_c
    # process rates on the clamped in-cloud fields
    rates = compute_m2005_rates(qc_i, nc_i, qr_i, nr_i, qi_i, ni_i, qs_i, ns_i, qg_i, ng_i,
                                qv_i, T, pres, rho, dt)
    # per-level warm/cold tendency assembly + PCC (in-cloud)
    ten_i = m2005_step_tendencies(T, qv_i, qc_i, qi_i, qr_i, qs_i, qg_i,
                                  nc_i, ni_i, nr_i, ns_i, ng_i, rates, pres, dt, igraup=igraup)
    for k, v in sf.items():
        ten_i[k] = ten_i[k] + v
    # in-cloud tendencies → grid-mean (×CF3D where CF3D>thresh)
    w = jnp.where(cf3d > _M_CF_THRESH, cf3d, 1.0)
    return {k: ten_i[k] * w for k in ten_i}


# ── Cloud + ice + snow terminal fall speeds + multi-species sedimentation ────────────
# module_mp_graupel.F90:4648-4672 (fall speeds) + the shared-NSTEP sub-step loop (:4647-4834).
# _M_G (gravitational acceleration) = grav from constants_clubb (module_mp_graupel.F90:519 `G = grav`), imported above.


def cloud_fall_speed(qc, nc, T, rho):
    """Cloud-water mass UMC + number UNC fall speeds, Stokes/viscosity regime (module_mp_graupel.F90:
    4648, ACN=G·ρw/(18·MU), MU=1.496e-6·T^1.5/(T+120), BC=2). PGAM/DLAMC from the cloud slope; no cap.
    Cloud-water sedimentation is folded into rcm_mc (QC3DTEN += QCSTEN, :4885)."""
    qc = jnp.asarray(qc, dtype=jnp.float64); nc = jnp.asarray(nc, dtype=jnp.float64)
    T = jnp.asarray(T, dtype=jnp.float64); rho = jnp.asarray(rho, dtype=jnp.float64)
    pgam, lamc = cloud_slope(qc, nc, rho)
    lc = jnp.where(lamc > 0.0, lamc, 1.0)
    mu = 1.496e-6 * T ** 1.5 / (T + 120.0)
    acn = _M_G * _KK_RHOW / (18.0 * mu)
    bc = 2.0
    # no droplets → no defined size → no sedimentation (nc>0 guard, like PRC). The oracle has
    # rcm>0 & Ncm=0 at ~40% of cloudy points (phantom cloud water) and rcm_sd_mg_morr=0 there.
    on = (qc >= _M_QSMALL) & (nc > 0.0)
    umc = jnp.where(on, acn * gamma(4.0 + bc + pgam) / (lc ** bc * gamma(pgam + 4.0)), 0.0)
    unc = jnp.where(on, acn * gamma(1.0 + bc + pgam) / (lc ** bc * gamma(pgam + 1.0)), 0.0)
    return umc, unc


def ice_fall_speed(qi, ni, rho):
    """Ice mass UMI + number UNI terminal fall speeds (module_mp_graupel.F90:4655-4656),
    capped at 1.2·(rhosu/rho)^0.35. CONS27=Γ(1+BI)=Γ(2)=1, CONS28=Γ(4+BI)/6=Γ(5)/6=4 (BI=1)."""
    qi = jnp.asarray(qi, dtype=jnp.float64); ni = jnp.asarray(ni, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lami, _ = ice_slope(qi, ni)
    li = jnp.where(lami > 0.0, lami, 1.0)
    ain = (_M_RHOSU / rho) ** 0.35 * _M_AI
    cons27 = gamma(jnp.asarray(1.0 + _M_BI)); cons28 = gamma(jnp.asarray(4.0 + _M_BI)) / 6.0
    cap = 1.2 * (_M_RHOSU / rho) ** 0.35
    on = qi >= _M_QSMALL
    umi = jnp.where(on, jnp.minimum(ain * cons28 / li ** _M_BI, cap), 0.0)
    uni = jnp.where(on, jnp.minimum(ain * cons27 / li ** _M_BI, cap), 0.0)
    return umi, uni


def snow_fall_speed(qs, ns, rho):
    """Snow mass UMS + number UNS fall speeds (module_mp_graupel.F90:4671-4672),
    capped at 1.2·(rhosu/rho)^0.54. CONS3=Γ(4+BS)/6, CONS5=Γ(1+BS)."""
    qs = jnp.asarray(qs, dtype=jnp.float64); ns = jnp.asarray(ns, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    lams, _ = snow_slope(qs, ns)
    ls = jnp.where(lams > 0.0, lams, 1.0)
    dcorr = (_M_RHOSU / rho) ** 0.54
    asn = dcorr * _M_AS
    cons3 = gamma(jnp.asarray(4.0 + _M_BS)) / 6.0; cons5 = gamma(jnp.asarray(1.0 + _M_BS))
    cap = 1.2 * dcorr
    on = qs >= _M_QSMALL
    ums = jnp.where(on, jnp.minimum(asn * cons3 / ls ** _M_BS, cap), 0.0)
    uns = jnp.where(on, jnp.minimum(asn * cons5 / ls ** _M_BS, cap), 0.0)
    return ums, uns


def morrison_sedimentation(qr, nr, qi, ni, qs, ns, rho, dzq, dt, qc=None, nc=None, T=None):
    """Multi-species sedimentation (rain + ice + snow [+ optional cloud], mass + number) with the
    SHARED CFL sub-step count NSTEP = max(int(RGVM·dt/dz + 1), 1), RGVM = max fall speed over ALL
    species/moments (module_mp_graupel.F90:4647-4834). Each species' fall speed is propagated downward
    through clear cells. Returns {qr,nr,qi,ni,qs,ns} (+ {qc,nc} if cloud is supplied) sedimentation
    tendencies [per s]. Cloud sed is folded into rcm_mc/Ncm_mc in the interface. (Graupel = 0 for nov11.)"""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    qr, nr, qi, ni, qs, ns, rho, dzq = a(qr), a(nr), a(qi), a(ni), a(qs), a(ns), a(rho), a(dzq)
    fmr, fnr = rain_fall_speed(qr, nr, rho)
    fmi, fni = ice_fall_speed(qi, ni, rho)
    fms, fns = snow_fall_speed(qs, ns, rho)
    falls = [fmr, fnr, fmi, fni, fms, fns]
    l_cloud = qc is not None
    if l_cloud:
        qc, nc, T = a(qc), a(nc), a(T)
        fmc, fnc = cloud_fall_speed(qc, nc, T, rho)
        falls += [fmc, fnc]
    speeds = [_fall_speed_propagate(f) for f in falls]
    rgvm = jnp.max(jnp.stack([jnp.max(f * dt / dzq) for f in speeds]))
    nstep = int(jnp.maximum(jnp.floor(rgvm + 1.0), 1.0))

    def sed(q, fall):
        return (_sediment(q * rho, fall, dzq, dt, nstep) / rho - q) / dt
    out = {'qr': sed(qr, speeds[0]), 'nr': sed(nr, speeds[1]), 'qi': sed(qi, speeds[2]),
           'ni': sed(ni, speeds[3]), 'qs': sed(qs, speeds[4]), 'ns': sed(ns, speeds[5])}
    if l_cloud:
        out['qc'] = sed(qc, speeds[6]); out['nc'] = sed(nc, speeds[7])
    return out


# ── CLUBB-Morrison interface: CLUBB state → CLUBB-form tendencies (*_mc) ──────────────
# morrison_microphys_module.F90:morrison_microphys_driver. Calls m2005_driver for the process
# tendencies, integrates the hydrometeor fields (post-process field + sedimentation, clipped),
# and forms the CLUBB tendencies: hydromet_mc = (field_final − field_initial)/dt (:786, the net
# change incl. sedimentation + clipping), rcm_mc/rvm_mc are the raw process tendencies, and
# thlm_mc = (T_in_K2thlm(T_final, exner, rcm_final) − thlm)/dt with thlm = (T − Lv/Cp·rcm)/exner.
# ── Size (slope) clamps: reset the working number where the gamma-slope is out of bounds ──
# All five clamps are applied PRE-RATE in the Fortran (warm block module_mp_graupel.F90:1881-2002,
# cold 2816-2968) — the reset NR3D/NC3D/NI3D/NS3D/NG3D feeds the rate computation AND is integrated;
# there is NO post-sed clamp (line 4509 only scales the SIZEFIX stat by CF3D). For the exponential
# species (rain/ice/snow/graupel, N(D)=N0·exp(−LAM·D)) the slope is LAM=(CONS·N/Q)^(1/3) and the reset
# number is N=LAM_clamped³·Q/CONS (N0=N·LAM ⇒ N0=LAM⁴·Q/CONS, N=N0/LAM). Cloud uses the gamma-shape form.
def _sizefix_exp_number(q, n, cons, lammin, lammax):
    """Reset an exponential-species number so its slope LAM=(cons·n/q)^(1/3) stays in [lammin,lammax]
    (module_mp_graupel.F90:1881 rain / 1966 snow / 1991 graupel / 2816 ice). N=LAM_clamped³·q/cons."""
    q = jnp.asarray(q, dtype=jnp.float64); n = jnp.asarray(n, dtype=jnp.float64)
    on = q >= _M_QSMALL
    lam = jnp.where(on, (cons * n / jnp.where(on, q, 1.0)) ** (1.0 / 3.0), 0.0)
    oob = on & ((lam > lammax) | (lam < lammin))
    lam_c = jnp.clip(lam, lammin, lammax)
    return jnp.where(oob, lam_c ** 3 * q / cons, n)


def _sizefix_cloud_number(qc, nc, rho, dofix_pgam=False, pgam_fixed=5.0):
    """Reset the cloud droplet number so LAMC stays in [(pgam+1)/60µm, (pgam+1)/1µm]
    (module_mp_graupel.F90:1938-1949). NC = LAMC_clamped³·qc·Γ(pgam+1)/(Γ(pgam+4)·CONS26)."""
    qc = jnp.asarray(qc, dtype=jnp.float64); nc = jnp.asarray(nc, dtype=jnp.float64)
    rho = jnp.asarray(rho, dtype=jnp.float64)
    on = qc >= _M_QSMALL
    if dofix_pgam:
        pgam = jnp.full_like(qc, pgam_fixed)
    else:
        pgam = 0.0005714 * (nc / 1.0e6 * rho) + 0.2714
        pgam = jnp.clip(1.0 / (pgam ** 2) - 1.0, 2.0, 10.0)
    g1, g4 = gamma(pgam + 1.0), gamma(pgam + 4.0)
    qc_s = jnp.where(on, qc, 1.0)
    lamc = (_M_CONS26 * nc * g4 / (qc_s * g1)) ** (1.0 / 3.0)
    lammin, lammax = (pgam + 1.0) / 60.0e-6, (pgam + 1.0) / 1.0e-6
    oob = on & ((lamc > lammax) | (lamc < lammin))
    lamc_c = jnp.clip(lamc, lammin, lammax)
    return jnp.where(oob, lamc_c ** 3 * qc * g1 / (g4 * _M_CONS26), nc)


def _size_clamp_numbers(qc, nc, qr, nr, qi, ni, qs, ns, qg, ng, rho):
    """Apply all five Morrison pre-rate slope clamps, returning the reset (nc, nr, ni, ns, ng).
    Faithful to module_mp_graupel.F90 (warm 1881-2002 / cold 2816-2968 — identical formulas)."""
    return (_sizefix_cloud_number(qc, nc, rho),
            _sizefix_exp_number(qr, nr, _M_CONS_RAIN, _M_LAMMINR, _M_LAMMAXR),
            _sizefix_exp_number(qi, ni, _M_CONS_ICE, _M_LAMMINI, _M_LAMMAXI),
            _sizefix_exp_number(qs, ns, _M_CONS_SNOW, _M_LAMMINS, _M_LAMMAXS),
            _sizefix_exp_number(qg, ng, _M_CONS_GRAUPEL, _M_LAMMING, _M_LAMMAXG))


# morrison_microphys_driver now lives in its Fortran-home module
# Microphys/morrison_microphys_module.py (mirror-refactor iter 114).
