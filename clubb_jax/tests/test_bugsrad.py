"""Validate the JAX BUGSrad components against direct Fortran-formula replicas (bit-exact).

Covers the self-contained first pieces of the BUGSrad correlated-k radiation port:
  - planck  (bugsrad_planck.F90) : blackbody band emission, 12-band Horner polynomial
  - newexp  (newexp.F90)         : the fast-exp approximation the two-stream solvers use
  - rayle   (rayle.F)            : Rayleigh optical depth + single-scattering albedo
"""
import numpy as np
import jax

jax.config.update("jax_enable_x64", True)

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Radiation.BUGSrad.bugsrad_planck import planck, _B, _MBIR
from clubb_jax.src.Radiation.BUGSrad.newexp import newexp
from clubb_jax.src.Radiation.BUGSrad.rayle import rayle
from clubb_jax.src.Radiation.BUGSrad.gascon import gascon, parm_ckd24, _CK24, _H2OBND, _IFLB
from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import GRAVITY, R_D, R_STAR, MW_H2O, F_VIRT
from clubb_jax.src.Radiation.BUGSrad.cloudg import cloudg, _PI as _CLOUDG_PI
from clubb_jax.src.Radiation.BUGSrad.comscp1 import comscp1
from clubb_jax.src.Radiation.BUGSrad.comscp2 import comscp2
from clubb_jax.src.Radiation.BUGSrad.two_rt_lw import two_rt_lw
from clubb_jax.src.Radiation.BUGSrad.two_rt_sw import two_rt_sw
from clubb_jax.src.Radiation.BUGSrad.gases_ckd import (pscale, qk, qkio3, qoph2o, qopo3s, qophc,
                                                       qopch4, qopn2o)


def _planck_replica(tt, ts, nbir):
    b = _B[nbir]
    ncol, nlm = tt.shape
    poly = lambda T: b[0] + T * (b[1] + T * (b[2] + T * (b[3] + T * (b[4] + T * b[5]))))
    bf = np.zeros((ncol, nlm + 1), dtype=np.float64)
    bf[:, 0] = poly(tt[:, 0])
    for i in range(1, nlm):
        bf[:, i] = poly(0.5 * (tt[:, i - 1] + tt[:, i]))
    bf[:, nlm] = poly(ts)
    return bf


def test_planck_vs_fortran_replica():
    rng = np.random.default_rng(0)
    ncol, nlm = 3, 40
    tt = 180.0 + 140.0 * rng.random((ncol, nlm))
    ts = 250.0 + 60.0 * rng.random(ncol)
    worst = max(np.abs(np.asarray(planck(tt, ts, n)) - _planck_replica(tt, ts, n)).max()
                / (np.abs(_planck_replica(tt, ts, n)).max() + 1e-300) for n in range(_MBIR))
    assert worst == 0.0, f"Planck JAX vs Fortran-replica worst rel {worst:.2e} (expected 0)"
    print(f"  Planck (12 IR bands) vs Fortran-formula replica: BIT-EXACT (rel {worst:.1e})  PASS")


def test_newexp_vs_fortran_replica():
    rng = np.random.default_rng(1)
    x = -8.0 * rng.random(2000)   # optical-depth-like negative args
    def ref(x):
        y1 = 1.0 - x * (0.2507213 - x * (0.0292732 - x * 0.0038278))
        return 1.0 / (y1 * y1) ** 2
    jx = np.asarray(newexp(x)); rf = ref(x)
    rel = np.abs(jx - rf).max() / np.abs(rf).max()
    assert rel == 0.0, f"newexp JAX vs Fortran-replica rel {rel:.2e} (expected 0)"
    # sanity: it is an APPROXIMATION of exp (must NOT equal np.exp — that's the whole point of faithfulness)
    approx_err = np.abs(jx - np.exp(x)).max() / np.abs(np.exp(x)).max()
    assert approx_err > 1e-6, "newexp should differ from true exp (it is the Fortran approximation)"
    print(f"  newexp vs Fortran-formula replica: BIT-EXACT (rel {rel:.1e}); differs from true exp by {approx_err:.1e}  PASS")


def test_rayle_vs_fortran_replica():
    rng = np.random.default_rng(2)
    ncol, nlm, mbs = 4, 30, 6
    amu0 = 0.2 + 0.7 * rng.random(ncol)
    ri = rng.random(mbs) * 1e-5
    pp = np.cumsum(20.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 100.0   # increasing level pressure
    ppl = 0.5 * (pp[:, :-1] + pp[:, 1:])
    def ref(ib):
        u0 = amu0
        x = (-3.902860e-6 * u0 * u0 + 6.120070e-6 * u0 + 4.177440e-6) if ib == 0 else np.full(ncol, ri[ib])
        trp = 29.267 * ppl * np.log(pp[:, 1:] / pp[:, :-1])
        return x[:, None] * trp
    worst = 0.0
    for ib in range(mbs):
        tray, wray = rayle(ib, amu0, ri, pp, ppl)
        rf = ref(ib)
        worst = max(worst, np.abs(np.asarray(tray) - rf).max() / (np.abs(rf).max() + 1e-300))
        assert np.all(np.asarray(wray) == 1.0)
    assert worst == 0.0, f"rayle JAX vs Fortran-replica worst rel {worst:.2e} (expected 0)"
    print(f"  rayle (Rayleigh τ + ω, 6 bands incl. band-1 quadratic): BIT-EXACT (rel {worst:.1e})  PASS")


def _gascon_replica(ib, pp, ppl, dp, tt, rmix):
    """Direct transcription of gascon.F + parm_ckd24 (NumPy)."""
    cb = _IFLB[ib]
    if cb == 0:
        return np.zeros_like(ppl)
    iband0 = cb - 1
    c = _CK24[iband0]; h2b = _H2OBND[iband0]
    out = np.zeros_like(ppl)
    for i in range(ppl.shape[0]):
        for l in range(ppl.shape[1]):
            if rmix[i, l] > 0.0:
                amnt = 10.0 * dp[i, l] * rmix[i, l] / GRAVITY
                patm = ppl[i, l] / 1013.25
                tv = tt[i, l] * (1.0 + F_VIRT * rmix[i, l])
                dz = (R_D / GRAVITY) * tv * np.log(pp[i, l + 1] / pp[i, l]) * 0.001
                factor = 0.25 / dz if dz < 0.25 else (1.50 / dz if dz > 1.50 else 1.0)
                dz1 = 0.25 if dz < 0.25 else (1.50 if dz > 1.50 else dz)
                amnt1 = amnt * factor
                ireg = 1 if np.log(amnt1) > h2b else 0
                ph2o = amnt1 * (R_STAR * 1.e4 * tt[i, l]) / (dz1 * 1.0e5 * MW_H2O * 1.01325e6)
                patmx = np.log(patm)
                tau = (c[ireg, 0] + c[ireg, 1] * np.log(amnt1) + c[ireg, 2] * tt[i, l]
                       + c[ireg, 3] * patmx + c[ireg, 4] * ph2o + c[ireg, 5] * amnt1
                       + c[ireg, 6] * np.log(ph2o))
                out[i, l] = np.exp(tau) / factor
    return out


def test_gascon_vs_fortran_replica():
    rng = np.random.default_rng(3)
    ncol, nlm = 3, 30
    pp = np.cumsum(20.0 + 15.0 * rng.random((ncol, nlm + 1)), axis=1) + 100.0   # increasing level pres [hPa]
    ppl = 0.5 * (pp[:, :-1] + pp[:, 1:])
    dp = pp[:, 1:] - pp[:, :-1]
    tt = 200.0 + 100.0 * rng.random((ncol, nlm))
    rmix = np.where(rng.random((ncol, nlm)) > 0.2, 1e-4 + 1e-2 * rng.random((ncol, nlm)), 0.0)
    worst = 0.0
    fired = 0
    for ib in range(18):                 # all global bands, incl. SW (iflb=0 → zeros) and LW
        jx = np.asarray(gascon(ib, pp, ppl, dp, tt, rmix))
        ref = _gascon_replica(ib, pp, ppl, dp, tt, rmix)
        den = np.abs(ref).max()
        if den > 0:
            fired += 1
            worst = max(worst, np.abs(jx - ref).max() / den)
        else:
            assert np.all(jx == 0.0), f"band {ib} should be all-zero (no continuum)"
    # uses log/exp (transcendentals) → match to libm precision, not necessarily 0 ULP
    assert worst < 1e-12, f"gascon JAX vs Fortran-replica worst rel {worst:.2e}"
    assert fired == 12, f"expected 12 LW continuum bands to fire, got {fired}"
    print(f"  gascon (CKD2.4 H2O continuum, 12 LW bands + 6 zero SW): vs Fortran replica rel {worst:.1e}  PASS")


def _cloudg_replica(ib, pp, tt, wcont, re, pdist, cnrw, cniw, cnri, cnii, xlam, flag):
    pi = float(np.pi)   # REFACTOR A3: float64 π (was the float32 value, for bit-faithfulness)
    eps = 1.0e-5
    cnr = cnri[ib] if flag else cnrw[ib]
    cni = cnii[ib] if flag else cniw[ib]
    xl = xlam[ib]
    p0 = pdist; p1 = p0 + 1.0; p2 = p0 + 2.0; f2 = p1 * p0; f3 = p2 * f2
    ncol, nlm = wcont.shape
    tcld = np.zeros_like(wcont); wcld = np.zeros_like(wcont); asycld = np.ones_like(wcont)
    for i in range(ncol):
        for l in range(nlm):
            if wcont[i, l] > eps:
                dz = 29.286 * np.log(pp[i, l + 1] / pp[i, l]) * tt[i, l]
                rm = re[i, l] / p2
                no = wcont[i, l] / ((4.0 * pi / 3.0) * f3 * 1.0e-6 * rm ** 3)
                area = 1.0e-6 * pi * f2 * no * rm ** 2
                c0 = 2.0 * area; c1 = c0 / f2
                xm = 2.0 * pi * rm / xl
                if ib == 0:
                    um = 2.0 * xm * (cnr - 1.0) * 1j
                    ext = c0 + 2.0 * c1 * np.real(p0 / (um * (um + 1.0) ** p1)
                                                  + 1.0 / (um ** 2 * (um + 1.0) ** p0) - 1.0 / um ** 2)
                    tcld[i, l] = ext * dz; wcld[i, l] = 0.999999; asycld[i, l] = 0.85
                else:
                    cm = complex(cnr, -cni)
                    um = 2.0 * xm * (cm - 1.0) * 1j
                    ext = c0 + 2.0 * c1 * np.real(p0 / (um * (um + 1.0) ** p1)
                                                  + 1.0 / (um ** 2 * (um + 1.0) ** p0) - 1.0 / um ** 2)
                    vm = 4.0 * xm * cni
                    expr = p0 / (vm * (vm + 1.0) ** p1) + 1.0 / (vm ** 2 * (vm + 1.0) ** p0) - 1.0 / vm ** 2
                    abs_ = area + c1 * expr   # REFACTOR A3: float64 (was sngl float32 truncation)
                    tcld[i, l] = ext * dz
                    if ext < abs_:
                        ext = abs_
                    wcld[i, l] = (ext - abs_) / ext; asycld[i, l] = 0.85
    return tcld, wcld, asycld


def test_cloudg_vs_fortran_replica():
    rng = np.random.default_rng(4)
    ncol, nlm, mb = 3, 25, 18
    pp = np.cumsum(20.0 + 15.0 * rng.random((ncol, nlm + 1)), axis=1) + 100.0
    tt = 230.0 + 60.0 * rng.random((ncol, nlm))
    wcont = np.where(rng.random((ncol, nlm)) > 0.3, 1e-4 + 0.5 * rng.random((ncol, nlm)), 0.0)  # g/m^3
    re = 5.0 + 20.0 * rng.random((ncol, nlm))                                                    # µm
    pdist = 2.0
    cnrw = 1.3 + 0.1 * rng.random(mb); cniw = 0.01 + 0.2 * rng.random(mb)
    cnri = 1.3 + 0.1 * rng.random(mb); cnii = 0.01 + 0.2 * rng.random(mb)
    xlam = 0.5 + 80.0 * rng.random(mb)
    # REFACTOR A3: cloudg now uses float64 π (the float32-π bit-faithfulness contrivance was removed).
    assert _CLOUDG_PI == float(np.pi), "cloudg π must be float64 π after REFACTOR A3"
    worst = 0.0
    for flag in (False, True):
        for ib in (0, 5, 11):           # band 1 (extinction-only) + two general bands
            jt, jw, ja = (np.asarray(x) for x in cloudg(ib, pp, tt, wcont, re, pdist,
                                                        cnrw, cniw, cnri, cnii, xlam, flag))
            rt, rw, ra = _cloudg_replica(ib, pp, tt, wcont, re, pdist, cnrw, cniw, cnri, cnii, xlam, flag)
            for j, r in ((jt, rt), (jw, rw), (ja, ra)):
                den = np.abs(r).max()
                worst = max(worst, np.abs(j - r).max() / (den + 1e-300))
    assert worst < 1e-12, f"cloudg JAX vs Fortran-replica worst rel {worst:.2e}"
    print(f"  cloudg (ADT cloud optics; complex ext, float64 abs+π; water+ice): vs replica rel {worst:.1e}  PASS")


def test_comscp_vs_fortran_replica():
    rng = np.random.default_rng(5)
    sh = (4, 20)
    R = lambda: rng.random(sh)
    Z = lambda: np.where(rng.random(sh) > 0.4, R(), 0.0)   # some zeros to exercise the else-branches
    taer, tcldi, tcldw, tgm, tray = Z(), Z(), Z(), R(), R()
    waer, wcldi, wcldw, wray = R(), R(), R(), R()
    asyaer, asycldi, asycldw, tg = R(), R(), R(), R()
    # comscp1
    j1 = [np.asarray(x) for x in comscp1(taer, tcldi, tcldw, tgm, tray, waer, wcldi, wcldw, wray,
                                         asyaer, asycldi, asycldw)]
    tcclr1 = tgm + tray + taer
    tccld1 = tcclr1 + tcldi + tcldw
    wwray = wray * tray; wwaer = waer * taer; wwcldi = wcldi * tcldi; wwcldw = wcldw * tcldw
    fwclr = wwray + wwaer; fwcld = fwclr + wwcldi + wwcldw
    asyclr = np.where(fwclr > 1e-10, asyaer * wwaer / np.where(fwclr > 1e-10, fwclr, 1.0), 1.0)
    asycld = np.where(fwcld > 1e-10, (asyaer * wwaer + asycldi * wwcldi + asycldw * wwcldw)
                      / np.where(fwcld > 1e-10, fwcld, 1.0), 1.0)
    r1 = [tccld1, tcclr1, asycld, asyclr, fwcld, fwclr]
    w1 = max(np.abs(j - r).max() / (np.abs(r).max() + 1e-300) for j, r in zip(j1, r1))
    # comscp2
    j2 = [np.asarray(x) for x in comscp2(tg, fwcld, fwclr, tccld1, tcclr1)]
    tcclr = tcclr1 + tg; tccld = tccld1 + tg
    wcclr = np.minimum(0.999999, np.where(tcclr > 0., fwclr / np.where(tcclr > 0., tcclr, 1.0), 0.0))
    wccld = np.minimum(0.999999, np.where(tccld > 0., fwcld / np.where(tccld > 0., tccld, 1.0), 0.0))
    r2 = [tccld, tcclr, wccld, wcclr]
    w2 = max(np.abs(j - r).max() / (np.abs(r).max() + 1e-300) for j, r in zip(j2, r2))
    assert w1 == 0.0 and w2 == 0.0, f"comscp1/2 not bit-exact: w1={w1:.2e} w2={w2:.2e}"
    print(f"  comscp1/2 (combine cloud+aerosol+Rayleigh+gas optical props): BIT-EXACT (rel {max(w1, w2):.1e})  PASS")


def _two_rt_lw_replica(ib, wc, asym, tau, es, bf, mbs):
    diffac = 2.0
    ncol, nlm = wc.shape
    ibms = ib - mbs
    fu = np.zeros((ncol, nlm + 1)); fd = np.zeros((ncol, nlm + 1))
    for i in range(ncol):
        re = np.zeros(nlm + 1); vd = np.zeros(nlm + 1)
        rr = np.zeros(nlm); tr = np.zeros(nlm); sigu = np.zeros(nlm); sigd = np.zeros(nlm)
        vu = np.zeros(nlm); td = np.zeros(nlm)
        for l in range(nlm):
            fact = asym[i, l] ** 2
            oms = ((1 - fact) * wc[i, l]) / (1 - fact * wc[i, l])
            taus = (1 - fact * wc[i, l]) * tau[i, l]
            beta0 = (4 + asym[i, l]) / (8 * (1 + asym[i, l]))
            t = diffac * (1 - oms * (1 - beta0)); r = diffac * oms * beta0
            kappa = np.sqrt(t ** 2 - r ** 2); rinf = r / (kappa + t)
            eggtau = np.exp(-kappa * taus)
            denom = 1 - rinf ** 2 * eggtau ** 2
            tr[l] = (1 - rinf ** 2) * eggtau / denom
            rr[l] = rinf * (1 - eggtau ** 2) / denom
            if taus < 0.1481e-2:
                sigu[l] = 0.5 * diffac * (bf[i, l] + bf[i, l + 1]) * taus; sigd[l] = sigu[l]
            else:
                aa = (t + r) * (1 - rr[l]) - (1 + rr[l] - tr[l]) / taus
                bb = -(t + r) * tr[l] + (1 + rr[l] - tr[l]) / taus
                cc = diffac * (1 - oms) / kappa ** 2
                sigu[l] = cc * (aa * bf[i, l] + bb * bf[i, l + 1])
                sigd[l] = cc * (bb * bf[i, l] + aa * bf[i, l + 1])
        for l in range(nlm):
            prop = 1 / (1 - re[l] * rr[l])
            re[l + 1] = rr[l] + tr[l] ** 2 * re[l] * prop
            vd[l + 1] = sigd[l] + (tr[l] * vd[l] + tr[l] * re[l] * sigu[l]) * prop
            vu[l] = (rr[l] * vd[l] + sigu[l]) * prop; td[l] = prop
        fu[i, nlm] = es[i, ibms] * bf[i, nlm]
        for l in range(nlm, 0, -1):
            fd[i, l] = re[l] * fu[i, l] + vd[l]
            fu[i, l - 1] = tr[l - 1] * fu[i, l] * td[l - 1] + vu[l - 1]
    return fu, fd


def test_two_rt_lw_vs_fortran_replica():
    rng = np.random.default_rng(6)
    ncol, nlm, mbs, mbir = 3, 30, 6, 12
    wc = 0.3 * rng.random((ncol, nlm))                 # SSA (LW: small)
    asym = 0.85 * rng.random((ncol, nlm))
    tau = 0.01 + 2.0 * rng.random((ncol, nlm))
    es = 0.95 + 0.05 * rng.random((ncol, mbir))
    bf = 5.0 + 20.0 * rng.random((ncol, nlm + 1))
    worst = 0.0
    for ib in (6, 11, 17):                              # LW bands (ib >= mbs)
        ju, jd = (np.asarray(x) for x in two_rt_lw(ib, wc, asym, tau, es, bf, mbs))
        ru, rd = _two_rt_lw_replica(ib, wc, asym, tau, es, bf, mbs)
        for j, r in ((ju, ru), (jd, rd)):
            worst = max(worst, np.abs(j - r).max() / (np.abs(r).max() + 1e-300))
    assert worst < 1e-12, f"two_rt_lw JAX vs Fortran-replica worst rel {worst:.2e}"
    print(f"  two_rt_lw (delta-Eddington 2-stream + adding, full-calc branch): vs replica rel {worst:.1e}  PASS")


def _two_rt_sw_replica(ib, slr, amu0, wc, asym, tau, asdir, asdif, mbs):
    eps = 1.0e-2
    ncol, nlm = wc.shape
    fudif = np.zeros((ncol, nlm + 1)); fddir = np.zeros((ncol, nlm + 1)); fddif = np.zeros((ncol, nlm + 1))
    for i in range(ncol):
        direct = np.zeros(nlm + 1); direct[0] = 1.0
        re = np.zeros(nlm + 1); vd = np.zeros(nlm + 1)
        rr = np.zeros(nlm); tr = np.zeros(nlm); sigu = np.zeros(nlm); sigd = np.zeros(nlm)
        vu = np.zeros(nlm); td = np.zeros(nlm); exptau = np.zeros(nlm)
        for l in range(nlm):
            fact = asym[i, l] ** 2
            oms = ((1 - fact) * wc[i, l]) / (1 - fact * wc[i, l])
            taus = (1 - fact * wc[i, l]) * tau[i, l]
            asy = asym[i, l] / (1 + asym[i, l])
            exptau[l] = np.exp(-taus / amu0[i]); direct[l + 1] = exptau[l] * direct[l]
            t = 0.25 * (7 - oms * (4 + 3 * asy)); r = -0.25 * (1 - oms * (4 - 3 * asy))
            kappa = np.sqrt(t ** 2 - r ** 2); rinf = r / (kappa + t)
            eggtau = np.exp(-kappa * taus); denom = 1 - rinf ** 2 * eggtau ** 2
            tr[l] = (1 - rinf ** 2) * eggtau / denom; rr[l] = rinf * (1 - eggtau ** 2) / denom
            d = kappa ** 2 - 1 / amu0[i] ** 2
            factb = 1 / eps if abs(d) < eps else 1 / d
            cc = oms * slr[i] * factb
            g3 = 0.5 - 0.75 * asy * amu0[i]; g4 = 1 - g3
            aa = g3 * (t - 1 / amu0[i]) + g4 * r; bb = g4 * (t + 1 / amu0[i]) + g3 * r
            sigu[l] = cc * ((aa - rr[l] * bb) - aa * tr[l] * exptau[l]) * direct[l]
            sigd[l] = cc * (-bb * tr[l] + (bb - rr[l] * aa) * exptau[l]) * direct[l]
        for l in range(nlm):
            prop = 1 / (1 - re[l] * rr[l])
            re[l + 1] = rr[l] + tr[l] ** 2 * re[l] * prop
            vd[l + 1] = sigd[l] + (tr[l] * vd[l] + tr[l] * re[l] * sigu[l]) * prop
            vu[l] = (rr[l] * vd[l] + sigu[l]) * prop; td[l] = prop
        fudif[i, nlm] = (asdif[i, ib] * vd[nlm] + asdir[i, ib] * slr[i] * amu0[i] * direct[nlm]) \
            / (1 - asdif[i, ib] * re[nlm])
        for l in range(nlm, 0, -1):
            fddif[i, l] = re[l] * fudif[i, l] + vd[l]
            fudif[i, l - 1] = tr[l - 1] * fudif[i, l] * td[l - 1] + vu[l - 1]
        for l in range(nlm + 1):
            fddir[i, l] = amu0[i] * slr[i] * direct[l]
    return fudif, fddir, fddif


def test_two_rt_sw_vs_fortran_replica():
    rng = np.random.default_rng(7)
    ncol, nlm, mbs = 3, 30, 6
    slr = 0.5 + 0.5 * rng.random(ncol)
    amu0 = 0.2 + 0.7 * rng.random(ncol)
    wc = 0.5 + 0.49 * rng.random((ncol, nlm))          # SW: high SSA
    asym = 0.85 * rng.random((ncol, nlm))
    tau = 0.01 + 1.5 * rng.random((ncol, nlm))
    asdir = 0.1 + 0.5 * rng.random((ncol, mbs)); asdif = 0.1 + 0.5 * rng.random((ncol, mbs))
    worst = 0.0
    for ib in (0, 3, 5):
        jf = [np.asarray(x) for x in two_rt_sw(ib, slr, amu0, wc, asym, tau, asdir, asdif, mbs)]
        rf = _two_rt_sw_replica(ib, slr, amu0, wc, asym, tau, asdir, asdif, mbs)
        for j, r in zip(jf, rf):
            worst = max(worst, np.abs(j - r).max() / (np.abs(r).max() + 1e-300))
    assert worst < 1e-12, f"two_rt_sw JAX vs Fortran-replica worst rel {worst:.2e}"
    print(f"  two_rt_sw (delta-Eddington 2-stream + direct beam + adding): vs replica rel {worst:.1e}  PASS")


def test_gases_ckd_data_parser():
    """The gases_ckd_data.h parser builds the correlated-k tables with correct shapes, fully filled
    (no missing T-coefficient slices), and the right spot-checked values."""
    from clubb_jax.src.Radiation.BUGSrad.gases_ckd_data import tables, KG, NUMPS, NUMPIR, NUMTS, NUMTIR
    t = tables()
    # 1D weights hk1..hk18 (length KG[band]) + fk1o3 (KG[0])
    for b in range(18):
        assert t[f'hk{b + 1}'].shape == (KG[b],), f'hk{b+1} shape'
    assert t['fk1o3'].shape == (KG[0],)
    # spot-check known values from the file
    assert t['hk1'][0] == 0.24 and t['hk1'][-1] == 0.008
    assert t['hk7'][0] == 0.7 and t['hk7'][1] == 0.3
    # SW H2O 3D arrays (NUMTS, NUMPS, KG_band)
    for b, name in [(2, 'c2h2o'), (3, 'c3h2o'), (4, 'c4h2o'), (5, 'c5h2o'), (6, 'c6h2o')]:
        assert t[name].shape == (NUMTS, NUMPS, KG[b - 1]), f'{name} shape {t[name].shape}'
        assert not np.any([np.all(t[name][k] == 0) for k in range(NUMTS)]), f'{name} has empty T-slice'
    # LW 3D arrays (NUMTIR, NUMPIR, KG_band) + 2D (NUMTIR, NUMPIR)
    for name in ('c7h2o', 'c8h2o', 'c9h2o', 'c10h2o', 'c11h2o', 'c12o3'):
        assert t[name].shape[:2] == (NUMTIR, NUMPIR), f'{name} shape {t[name].shape}'
    for name in ('c10ch4', 'c10n2o', 'c11ch4', 'c11n2o'):
        assert t[name].shape == (NUMTIR, NUMPIR), f'{name} shape {t[name].shape}'
    assert t['c2h2o'][0, 0, 0] == -17.35 and t['c2h2o'][0, 0, 1] == -14.07   # first T_1 values
    assert t['c10ch4'][0, 0] == -8.909
    print(f"  gases_ckd_data.h parser ({len(t)} arrays: hk*/fk1o3/c*h2o/c*o3/c*ch4/c*n2o, correct shapes "
          f"+ values)  PASS")


def test_gases_dispatch_vs_fortran_replica():
    """The `gases` 18-band dispatch (table indexing + overlap formulas + hk weights) vs a NumPy replica."""
    from clubb_jax.src.Radiation.BUGSrad.gases_ckd import gases, pscale
    from clubb_jax.src.Radiation.BUGSrad.gases_ckd_data import tables, KG, STANPS, STANPIR
    from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import MOLAR_VOLUME as MV, GRAVITY as G, \
        MW_H2O as MWh, MW_O3 as MWo, MW_DRY_AIR as MWd
    T = tables()
    rng = np.random.default_rng(9)
    ncol, nlm = 2, 20
    pp = np.cumsum(15.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 1.0
    dp = pp[:, 1:] - pp[:, :-1]; ppl = 0.5 * (pp[:, :-1] + pp[:, 1:])
    tt = 200.0 + 80.0 * rng.random((ncol, nlm))
    rmix = 1e-3 * rng.random((ncol, nlm)); o3mix = 1e-6 * rng.random((ncol, nlm))
    umco2, umch4, umn2o = 330.0, 1.6, 0.28
    sol = {1: 619.618, 2: 484.295, 3: 149.845, 4: 48.7302, 5: 31.6576, 6: 5.79927}

    def qk_np(coefk, t0=245.0, ylo=-65.0, yhi=75.0, yset=75.0, sw=True):
        pkd, ip1, ip2 = (np.asarray(x) for x in pscale(ppl, STANPS if sw else STANPIR))
        y = tt - t0
        y = np.where(y < ylo, ylo, y)
        y = np.where(y > yhi, yset, y)
        kf = lambda ip: np.exp(coefk[0][ip] + coefk[1][ip] * y + coefk[2][ip] * y ** 2)
        x1, x2 = kf(ip1), kf(ip2)
        return np.where(ip1 == ip2, x1 * pkd, x1 + (x2 - x1) * pkd)

    worst = 0.0
    for ib in range(1, 19):
        sw = ib <= 6
        pkd, ip1, ip2 = (np.asarray(x) for x in pscale(ppl, STANPS if sw else STANPIR))
        for ig in range(KG[ib - 1]):
            tg_j, hk_j = gases(ib, ig, dp, tt, rmix, o3mix, umco2, umch4, umn2o, pkd, ip1, ip2, pp)
            tg_j = np.asarray(tg_j)
            # replica
            if ib == 1:
                tg_r = o3mix * dp * float(T['fk1o3'][ig]) * MV / MWo * 10.0 / G; hk_r = sol[1] * T['hk1'][ig]
            elif ib in (2, 3, 4, 5, 6):
                tg_r = qk_np(T[f'c{ib}h2o'][:, :, ig], sw=True) * rmix * dp * MV / MWh * 10.0 / G
                hk_r = sol[ib] * T[f'hk{ib}'][ig]
            elif ib in (7, 8, 9, 13, 16, 17, 18):
                tg_r = qk_np(T[f'c{ib}h2o'][:, :, ig], sw=False) * rmix * dp * MV / MWh * 10.0 / G
                hk_r = T[f'hk{ib}'][ig]
            elif ib in (10, 11):
                tg1 = qk_np(T[f'c{ib}h2o'][:, :, ig], sw=False) * rmix * dp * MV / MWh * 10.0 / G
                tg2 = qk_np(T[f'c{ib}ch4'], sw=False) * dp * MV * 10.0 / G * 1.6e-6 / MWd
                tg3 = qk_np(T[f'c{ib}n2o'], sw=False) * dp * MV * 10.0 / G * 0.28e-6 / MWd
                tg_r = tg1 + tg2 / 1.6 * umch4 + tg3 / 0.28 * umn2o; hk_r = T[f'hk{ib}'][ig]
            elif ib == 12:
                tg1 = qk_np(T['c12o3'][:, :, ig], 250.0, -70.0, 75.0, 70.0, sw=False) * o3mix * dp * MV / MWo * 10.0 / G
                tg2 = qk_np(T['c12h2o'], sw=False) * rmix * dp * MV / MWh * 10.0 / G
                tg_r = tg1 + tg2; hk_r = T['hk12'][ig]
            else:  # 14, 15
                pq = np.where(pp[:, :-1] >= 63.1, rmix, 0.0)
                fkg = qk_np(T[f'c{ib}hca'][:, :, ig], sw=False) / 330.0 * umco2 + pq * qk_np(T[f'c{ib}hcb'][:, :, ig], sw=False)
                tg_r = fkg * dp; hk_r = T[f'hk{ib}'][ig]
            assert abs(hk_j - hk_r) < 1e-13 * (abs(hk_r) + 1e-30), f'band {ib} hk'
            worst = max(worst, np.abs(tg_j - tg_r).max() / (np.abs(tg_r).max() + 1e-300))
    assert worst < 1e-12, f"gases dispatch worst rel {worst:.2e}"
    print(f"  gases 18-band dispatch (H2O/O3/CH4/N2O/CO2 overlaps + hk weights): vs replica rel {worst:.1e}  PASS")


def test_gases_ckd_helpers_vs_fortran_replica():
    rng = np.random.default_rng(8)
    ncol, nlm, n = 3, 25, 19
    stanp = np.sort(1.0 + rng.random(n)) * 60.0 + 1.0          # ascending standard pressures
    # ppl spanning below stanp[0], interior, and above stanp[-1]
    ppl = np.concatenate([np.full((ncol, 2), stanp[0] * 0.5),
                          stanp[0] + (stanp[-1] - stanp[0]) * rng.random((ncol, nlm - 4)),
                          np.full((ncol, 2), stanp[-1] * 1.2)], axis=1)
    # realistic correlated-k coef magnitudes: a~O(1), b~1e-2, c~1e-4 (so exp(a+b·y+c·y²) doesn't overflow)
    coefk = np.stack([rng.standard_normal(n) - 8.0, 1e-2 * rng.standard_normal(n), 1e-4 * rng.standard_normal(n)])
    tt = 180.0 + 150.0 * rng.random((ncol, nlm))
    dp = 5.0 + 10.0 * rng.random((ncol, nlm))
    rmix = 1e-3 * rng.random((ncol, nlm)); o3mix = 1e-6 * rng.random((ncol, nlm))

    # pscale replica
    pkd_j, ip1_j, ip2_j = (np.asarray(x) for x in pscale(ppl, stanp))
    pkd_r = np.zeros((ncol, nlm)); ip1_r = np.zeros((ncol, nlm), int); ip2_r = np.zeros((ncol, nlm), int)
    for i in range(ncol):
        for l in range(nlm):
            p = ppl[i, l]
            if p < stanp[0]:
                pkd_r[i, l] = p / stanp[0]; ip1_r[i, l] = 0; ip2_r[i, l] = 0
            elif p >= stanp[-1]:
                pkd_r[i, l] = (p - stanp[n - 2]) / (stanp[n - 1] - stanp[n - 2]); ip1_r[i, l] = n - 2; ip2_r[i, l] = n - 1
            else:
                j = 0
                while p >= stanp[j + 1]:
                    j += 1
                pkd_r[i, l] = (p - stanp[j]) / (stanp[j + 1] - stanp[j]); ip1_r[i, l] = j; ip2_r[i, l] = j + 1
    assert np.array_equal(ip1_j, ip1_r) and np.array_equal(ip2_j, ip2_r), "pscale indices mismatch"
    assert np.abs(pkd_j - pkd_r).max() == 0.0, "pscale pkd mismatch"

    # qk replica
    def qk_ref(t0, ylo, yhi, yhi_set):
        out = np.zeros((ncol, nlm))
        for i in range(ncol):
            for l in range(nlm):
                y1 = tt[i, l] - t0
                if y1 < ylo:
                    y1 = ylo
                if y1 > yhi:
                    y1 = yhi_set
                i1, i2 = ip1_r[i, l], ip2_r[i, l]
                x1 = np.exp(coefk[0, i1] + coefk[1, i1] * y1 + coefk[2, i1] * y1 ** 2)
                if i1 == i2:
                    out[i, l] = x1 * pkd_r[i, l]
                else:
                    x2 = np.exp(coefk[0, i2] + coefk[1, i2] * y1 + coefk[2, i2] * y1 ** 2)
                    out[i, l] = x1 + (x2 - x1) * pkd_r[i, l]
        return out
    fkg = np.asarray(qk(coefk, tt, pkd_j, ip1_j, ip2_j))
    rk = qk_ref(245.0, -65.0, 75.0, 75.0)
    rel_qk = np.abs(fkg - rk).max() / np.abs(rk).max()
    fkg3 = np.asarray(qkio3(coefk, tt, pkd_j, ip1_j, ip2_j))
    rel_qk3 = np.abs(fkg3 - qk_ref(250.0, -70.0, 75.0, 70.0)).max() / np.abs(qk_ref(250.0, -70.0, 75.0, 70.0)).max()

    # qop transmission functions (exact formulas)
    from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import MOLAR_VOLUME, GRAVITY, MW_H2O, MW_O3, MW_DRY_AIR
    h2o_r = fkg * rmix * dp * MOLAR_VOLUME / MW_H2O * 10.0 / GRAVITY
    rel_h2o = np.abs(np.asarray(qoph2o(fkg, dp, rmix)) - h2o_r).max() / np.abs(h2o_r).max()
    o3s_r = o3mix * dp * 0.7 * MOLAR_VOLUME / MW_O3 * 10.0 / GRAVITY
    rel_o3 = np.abs(np.asarray(qopo3s(0.7, dp, o3mix)) - o3s_r).max() / np.abs(o3s_r).max()
    ch4_r = fkg * dp * MOLAR_VOLUME * 10.0 / GRAVITY * 1.6e-6 / MW_DRY_AIR
    rel_ch4 = np.abs(np.asarray(qopch4(fkg, dp)) - ch4_r).max() / np.abs(ch4_r).max()
    rel_hc = np.abs(np.asarray(qophc(fkg, dp)) - fkg * dp).max() / np.abs(fkg * dp).max()

    worst = max(rel_qk, rel_qk3, rel_h2o, rel_o3, rel_ch4, rel_hc)
    assert worst < 1e-13, f"gases helpers worst rel {worst:.2e}"
    print(f"  gases_ckd helpers (pscale exact + qk/qkio3 interp + qop* transmission): vs replica rel {worst:.1e}  PASS")


def test_bugs_lwr():
    """LW band driver: the no-cloud invariant (all-sky == clear-sky, bit-exact), finiteness, physical
    flux ranges, and clouds adding downward LW."""
    from clubb_jax.src.Radiation.BUGSrad.bugs_lwr import bugs_lwr
    rng = np.random.default_rng(11)
    ncol, nlm = 2, 25
    pp = np.cumsum(15.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 1.0
    dp = pp[:, 1:] - pp[:, :-1]; ppl = 0.5 * (pp[:, :-1] + pp[:, 1:])
    tt = 220.0 + 70.0 * rng.random((ncol, nlm)); ts = 290.0 + 10.0 * rng.random(ncol)
    rmix = 1e-3 * rng.random((ncol, nlm)); o3mix = 1e-6 * rng.random((ncol, nlm))
    cw0 = np.zeros((ncol, nlm)); ci0 = np.zeros((ncol, nlm))
    fd, fu, fdcl, fucl = (np.asarray(x) for x in
                          bugs_lwr(ts, ppl, dp, tt, rmix, cw0, ci0, o3mix, pp, 330.0, 1.6, 0.28))
    assert np.array_equal(fd, fdcl) and np.array_equal(fu, fucl), "no-cloud all-sky != clear-sky"
    assert np.all(np.isfinite(fd)) and np.all(np.isfinite(fu))
    assert np.all((fu[:, 0] > 100) & (fu[:, 0] < 600)), f"OLR unphysical {fu[:, 0]}"   # outgoing LW
    cw = np.where(rng.random((ncol, nlm)) > 0.5, 0.1, 0.0)
    ci = np.where(rng.random((ncol, nlm)) > 0.7, 0.05, 0.0)
    fd2, fu2, fdc2, fuc2 = (np.asarray(x) for x in
                            bugs_lwr(ts, ppl, dp, tt, rmix, cw, ci, o3mix, pp, 330.0, 1.6, 0.28))
    assert np.all(np.isfinite(fd2)) and np.all(fd2 >= fdc2 - 1e-9), "clouds should add downward LW"
    print(f"  bugs_lwr (LW band driver): no-cloud all-sky==clear-sky (bit-exact), OLR≈{fu[:, 0].mean():.0f} W/m²  PASS")


def test_bugs_swr():
    """SW band driver: no-cloud invariant (all-sky==clear-sky bit-exact), the TOA incident solar
    (fdsw[TOA] = amu0·slr·Σband_solar = amu0·slr·1339.946) validating the accumulation, and finiteness."""
    from clubb_jax.src.Radiation.BUGSrad.bugs_swr import bugs_swr
    rng = np.random.default_rng(12)
    ncol, nlm = 2, 25
    pp = np.cumsum(15.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 1.0
    dp = pp[:, 1:] - pp[:, :-1]; ppl = 0.5 * (pp[:, :-1] + pp[:, 1:])
    tt = 230.0 + 60.0 * rng.random((ncol, nlm)); ts = 290.0 * np.ones(ncol)
    rmix = 1e-3 * rng.random((ncol, nlm)); o3mix = 1e-6 * rng.random((ncol, nlm))
    amu0 = np.array([0.5, 0.8]); slr = np.ones(ncol); alb = 0.2 * np.ones(ncol)
    cw0 = np.zeros((ncol, nlm)); ci0 = np.zeros((ncol, nlm))
    out = [np.asarray(x) for x in
           bugs_swr(ts, amu0, slr, alb, alb, alb, alb, ppl, dp, tt, rmix, cw0, ci0, o3mix, pp, 330.0, 1.6, 0.28)]
    fdsw, fusw, fdswcl, fuswcl = out[:4]
    assert np.array_equal(fdsw, fdswcl) and np.array_equal(fusw, fuswcl), "no-cloud all-sky != clear-sky"
    assert all(np.all(np.isfinite(x)) for x in out)
    sigsol = 1339.946   # Σ of the 6 SW band solar energies (619.618+484.295+149.845+48.7302+31.6576+5.79927)
    assert np.allclose(fdsw[:, 0], amu0 * slr * sigsol, rtol=1e-4), f"TOA solar wrong {fdsw[:, 0]}"
    assert np.all(fdsw[:, 0] >= fdsw[:, -1]), "SW flux should decrease downward"
    print(f"  bugs_swr (SW band driver): no-cloud all-sky==clear-sky, TOA solar = amu0·slr·{sigsol:.0f}  PASS")


def test_bugs_rad():
    """Main orchestration: LW heating conserves (telescoping flux divergence), night columns get zero
    SW heating, the no-cloud invariant (all-sky heating == clear-sky), and physical heating rates."""
    from clubb_jax.src.Radiation.BUGSrad.bugs_rad import bugs_rad
    from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import GRAVITY, CP_DRY_AIR
    rng = np.random.default_rng(13)
    ncol, nlm = 2, 30
    pl2 = np.cumsum(15.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 1.0
    dpl = pl2[:, 1:] - pl2[:, :-1]; ppl = 0.5 * (pl2[:, :-1] + pl2[:, 1:])
    tl = 220.0 + 70.0 * rng.random((ncol, nlm)); ts = 290.0 * np.ones(ncol)
    ql = 1e-3 * rng.random((ncol, nlm)); o3 = 1e-6 * rng.random((ncol, nlm))      # specific humidity, ozone
    amu0 = np.array([0.5, 0.005]); slr = np.ones(ncol); alb = 0.2 * np.ones(ncol)   # col 0 day, col 1 night
    qcwl = np.where(rng.random((ncol, nlm)) > 0.6, 1e-4, 0.0)                      # cloud-water mixing ratio [kg/kg]
    qcil = np.zeros((ncol, nlm)); acld = np.where(qcwl > 0, 1.0, 0.0)             # cloud ice, cloud fraction
    r = bugs_rad(pl2, ppl, dpl, tl, ql, qcwl, qcil, o3, ts, amu0, slr, alb, alb, alb, alb, acld, s0=1367.0)
    atl, asl = np.asarray(r['atl']), np.asarray(r['asl'])
    hf = GRAVITY * 0.01 / CP_DRY_AIR
    fnet = np.asarray(r['fulw'] - r['fdlw'])
    assert np.allclose((atl * dpl).sum(1), -hf * (fnet[:, 0] - fnet[:, nlm])), "LW heating not conservative"
    assert np.all(asl[1] == 0.0), "night column should have zero SW heating"
    assert np.any(asl[0] != 0.0), "day column should have nonzero SW heating"
    z = np.zeros((ncol, nlm))
    r2 = bugs_rad(pl2, ppl, dpl, tl, ql, z, z, o3, ts, amu0, slr, alb, alb, alb, alb, z, s0=1367.0)
    assert np.array_equal(np.asarray(r2['atl']), np.asarray(r2['atlcl'])), "no-cloud LW all-sky != clear-sky"
    assert np.all(np.isfinite(atl)) and np.all(np.isfinite(asl))
    print(f"  bugs_rad (preprocessing + orchestration): LW conserves, night SW=0, no-cloud invariant  PASS")


def test_bugsrad_driver():
    """CLUBB↔BUGSrad interface: std-atmosphere loads (50 levels), determine_extended_atmos_bounds gives
    sane bounds for a gabls3-like grid (model top 5 km, radiation_top 30 km), and compute_bugsrad_radiation
    produces a strictly-increasing radiation-grid pressure, finite CLUBB-grid heating, and zero SW at night."""
    import numpy as np
    from clubb_jax.src.Input_fields.sounding import load_extended_std_atm
    from clubb_jax.src.Radiation.bugsrad_driver import (
        determine_extended_atmos_bounds, build_rad_grid_setup,
        compute_bugsrad_radiation)
    ext = load_extended_std_atm()
    assert ext['alt'].shape == (50,) and ext['alt'][0] == 1000.0 and ext['alt'][-1] == 50000.0
    nzm, dz, ncol = 51, 100.0, 2
    zm = np.arange(nzm) * dz                                  # 0..5000 m
    zm_spacing = np.full(nzm - 1, dz)
    p_zm = 1013e2 * np.exp(-zm / 8000.0)                      # Pa, hydrostatic (top ≈ 542 mb)
    rad_top = 30000.0
    bottom, top, rsize, linbuf = determine_extended_atmos_bounds(zm, zm_spacing, p_zm, rad_top, ext)
    assert rsize >= 1 and linbuf >= 0 and top >= bottom, "degenerate extension bounds"
    assert ext['alt'][bottom - 1] >= 5000.0 and ext['alt'][top - 1] <= 30000.0, "bounds outside [5,30] km"
    setup = build_rad_grid_setup(zm, zm_spacing, p_zm, rad_top, ext)
    assert setup['buffer'] == linbuf + rsize
    # synthetic state on the nzt grid
    nzt = nzm - 1; rng = np.random.default_rng(7)
    zt = 0.5 * (zm[1:] + zm[:-1])
    T = np.tile(288.0 - 6.5e-3 * zt, (ncol, 1))
    p_zt = np.tile(1013e2 * np.exp(-zt / 8000.0), (ncol, 1))
    exner = (p_zt / 1.0e5) ** (287.0 / 1004.0)
    rho = np.tile(p_zm / (287.0 * (288.0 - 6.5e-3 * zm)), (ncol, 1))
    rtm = 5e-3 * np.ones((ncol, nzt)); rcm = np.zeros((ncol, nzt)); rcm[:, 5:8] = 1e-4
    cf = np.where(rcm > 0, 1.0, 0.0); rim = np.zeros_like(rcm); rsm = np.zeros_like(rcm); isf = np.zeros_like(rcm)
    amu0 = np.array([0.5, 0.0]); slr = np.ones(ncol); alb = 0.23 * np.ones(ncol)
    p_zm2 = np.tile(p_zm, (ncol, 1))
    res = compute_bugsrad_radiation(setup, T, rcm, rtm, rsm, rim, cf, isf, p_zt, p_zm2, rho, exner,
                                    amu0, slr, alb, alb, alb, alb, sol_const=1321.0)
    assert res['radht'].shape == (ncol, nzt), "radht wrong shape"
    assert np.all(np.isfinite(np.asarray(res['radht']))), "non-finite heating rate"
    assert np.all(np.asarray(res['radht_SW'])[1] == 0.0), "night column SW heating should be zero"
    assert np.any(np.asarray(res['radht_SW'])[0] != 0.0), "day column should have SW heating"
    print(f"  bugsrad_driver (std-atm + bounds + grid map): bounds=({bottom},{top},rsize={rsize},"
          f"buf={linbuf}), heating finite, night SW=0  PASS")


def test_bugsrad_radiation_dispatch():
    """End-to-end wiring: advance_clubb_radiation(rad_scheme='bugsrad') builds the grid setup, maps a gabls3-like
    CLUBB state through bugs_rad, and writes a finite K/s heating rate into state['radht'] (shape (1,nzt))."""
    import numpy as np
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.Radiation.radiation_module import advance_clubb_radiation
    gr = setup_grid(ngrdcol=1, deltaz=100.0, zm_init=0.0, zm_top=5000.0)
    nzt, nzm = gr.nzt, gr.nzm
    zt = np.asarray(gr.zt)[0]
    p_zt = (1013e2 * np.exp(-zt / 8000.0))[None, :]
    exner = (p_zt / 1.0e5) ** (287.0 / 1004.67)
    thlm = np.full((1, nzt), 295.0); rtm = np.full((1, nzt), 5e-3)
    rcm = np.zeros((1, nzt)); rcm[:, 5:9] = 2e-4
    cloud_frac = np.where(rcm > 0, 1.0, 0.0)
    rho_zm = (np.asarray(p_zt) / (287.0 * 290.0))
    rho_zm = np.concatenate([rho_zm, rho_zm[:, -1:]], axis=1)        # pad to nzm
    cfg = dict(rad_scheme="bugsrad", radiation_top=30000.0, sol_const=1321.0, slr=1.0,
               alvdr=0.23, alvdf=0.23, alndr=0.23, alndf=0.23,
               l_fix_cos_solar_zen=True, cos_solar_zen_values=[0.5])
    state = dict(cfg=cfg, gr=gr, ngrdcol=1, nzt=nzt, nzm=nzm, stats_writer=None,
                 rad_scheme="bugsrad", p_in_Pa=p_zt, exner=exner, thlm=thlm, rtm=rtm, rcm=rcm,
                 cloud_frac=cloud_frac, ice_supersat_frac=np.zeros((1, nzt)), rho_zm=rho_zm,
                 radht=np.zeros((1, nzt)))
    advance_clubb_radiation(state, time_current=43200.0)                  # local noon-ish
    assert state['radht'].shape == (1, nzt), "radht wrong shape"
    assert np.all(np.isfinite(state['radht'])), "non-finite radht from bugsrad dispatch"
    assert '_bugsrad_setup' in state, "rad-grid setup not cached"
    assert np.any(state['radht_SW'] != 0.0), "daytime SW heating should be nonzero"
    advance_clubb_radiation(state, time_current=43500.0)                  # second call reuses cached setup
    assert np.all(np.isfinite(state['radht']))
    print(f"  bugsrad dispatch (advance_clubb_radiation end-to-end): radht finite {state['radht'].shape}, "
          f"setup cached, day SW≠0  PASS")


def test_bugs_rad_differentiable():
    """The BUGSrad radiative transfer is DIFFERENTIABLE (a stated project goal): jax.grad of a radiative
    loss (TOA outgoing LW + total SW heating) w.r.t. temperature and cloud water is finite AND nonzero
    (the radiation responds to its physical inputs). cloudg is now fully float64 (REFACTOR A3)."""
    import jax, jax.numpy as jnp
    from clubb_jax.src.Radiation.BUGSrad.bugs_rad import bugs_rad
    rng = np.random.default_rng(1); ncol, nlm = 1, 40
    pl2 = jnp.asarray(np.cumsum(15.0 + 10.0 * rng.random((ncol, nlm + 1)), axis=1) + 1.0)
    dpl = pl2[:, 1:] - pl2[:, :-1]; ppl = 0.5 * (pl2[:, :-1] + pl2[:, 1:])
    tl0 = jnp.asarray(220.0 + 70.0 * rng.random((ncol, nlm))); ts = jnp.asarray(290.0 * np.ones(ncol))
    ql = jnp.asarray(1e-3 * rng.random((ncol, nlm))); o3 = jnp.asarray(1e-6 * rng.random((ncol, nlm)))
    qcwl0 = jnp.asarray(np.where(rng.random((ncol, nlm)) > 0.5, 1e-4, 0.0)); qcil = jnp.zeros((ncol, nlm))
    acld = jnp.asarray(np.where(np.asarray(qcwl0) > 0, 1.0, 0.0))
    amu0 = jnp.asarray([0.5]); slr = jnp.ones(ncol); alb = jnp.asarray(0.2 * np.ones(ncol))

    def loss(tl, qcwl):
        r = bugs_rad(pl2, ppl, dpl, tl, ql, qcwl, qcil, o3, ts, amu0, slr, alb, alb, alb, alb, acld, s0=1367.0)
        return jnp.sum(r['fulw'][:, 0]) + jnp.sum(r['asl'])

    g_t, g_q = jax.grad(loss, argnums=(0, 1))(tl0, qcwl0)
    assert jnp.all(jnp.isfinite(g_t)) and jnp.all(jnp.isfinite(g_q)), "non-finite BUGSrad gradient"
    assert float(jnp.max(jnp.abs(g_t))) > 0 and float(jnp.max(jnp.abs(g_q))) > 0, "BUGSrad gradient is zero"
    print(f"  bugs_rad differentiable: grad finite+nonzero (|dL/dT|max={float(jnp.max(jnp.abs(g_t))):.2e}, "
          f"|dL/dqc|max={float(jnp.max(jnp.abs(g_q))):.2e})  PASS")


if __name__ == "__main__":
    import gc
    print("BUGSrad component validation:")
    # The band-driver tests (bugs_lwr/swr/rad/driver) call BUGSrad EAGERLY, which accumulates JAX device
    # buffers (~hundreds of MB/call, the same leak the product path fixes by jitting); clear caches + gc
    # between tests so the whole suite stays within memory (else the cumulative eager footprint OOMs).
    tests = [
        test_planck_vs_fortran_replica, test_newexp_vs_fortran_replica, test_rayle_vs_fortran_replica,
        test_gascon_vs_fortran_replica, test_cloudg_vs_fortran_replica, test_comscp_vs_fortran_replica,
        test_two_rt_lw_vs_fortran_replica, test_two_rt_sw_vs_fortran_replica, test_gases_ckd_data_parser,
        test_gases_ckd_helpers_vs_fortran_replica, test_gases_dispatch_vs_fortran_replica,
        test_bugs_lwr, test_bugs_swr, test_bugs_rad, test_bugsrad_driver, test_bugsrad_radiation_dispatch,
        test_bugs_rad_differentiable,
    ]
    for t in tests:
        t()
        jax.clear_caches(); gc.collect()
    print("All BUGSrad component tests PASSED.")
