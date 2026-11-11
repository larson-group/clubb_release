// Native Morrison numerical engine in call-local C++ host storage.
// Default-real arithmetic/coefficients remain float; core constants retain their
// imported kind. Array bounds and lookup slots are zero-based C++ coordinates.
// All saved module fields live in the appended case record.
// GRAUPEL_INIT initializes that record before M2005MICRO_GRAUPEL advances a column.
// References: Morrison et al.(2005),JAS62,1665–1677; source V3.3 change history.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/Morrison_microphys/engine_contracts.hpp"
#include <complex>
#include <iostream>
namespace clubb_kokkos::morrison_engine {
inline void M2005MICRO_GRAUPEL(
    ScratchWorkspace& workspace, MorrisonColumn QC3DTEN, MorrisonColumn QI3DTEN,
    MorrisonColumn QNI3DTEN, MorrisonColumn QR3DTEN, MorrisonColumn NC3DTEN, MorrisonColumn NI3DTEN,
    MorrisonColumn NS3DTEN, MorrisonColumn NR3DTEN, MorrisonColumn QC3D, MorrisonColumn QI3D,
    MorrisonColumn QNI3D, MorrisonColumn QR3D, MorrisonColumn NC3D, MorrisonColumn NI3D,
    MorrisonColumn NS3D, MorrisonColumn NR3D, MorrisonColumn T3DTEN, MorrisonColumn QV3DTEN,
    MorrisonColumn T3D, MorrisonColumn QV3D, MorrisonColumn PRES, MorrisonColumn RHO,
    MorrisonColumn DZQ, MorrisonColumn W3D, MorrisonColumn WVAR, MorrisonColumn FR, float& PRECRT,
    float& SNOWRT, MorrisonColumn EFFC, MorrisonColumn EFFI, MorrisonColumn EFFS,
    MorrisonColumn EFFR, float DT, int IMS, int IME, int JMS, int JME, int KMS, int KME, int ITS,
    int ITE, int JTS, int JTE, int KTS, int KTE, MorrisonColumn QG3DTEN, MorrisonColumn NG3DTEN,
    MorrisonColumn QG3D, MorrisonColumn NG3D, MorrisonColumn EFFG, MorrisonColumn QGSTEN,
    MorrisonColumn QRSTEN, MorrisonColumn QISTEN, MorrisonColumn QNISTEN, MorrisonColumn QCSTEN,
    MorrisonColumn NGSTEN, MorrisonColumn NRSTEN, MorrisonColumn NISTEN, MorrisonColumn NSSTEN,
    MorrisonColumn NCSTEN, MorrisonColumn CF3D, MorrisonColumn PRC, MorrisonColumn PRA,
    MorrisonColumn PSMLT, MorrisonColumn EVPMS, MorrisonColumn PRACS, MorrisonColumn EVPMG,
    MorrisonColumn PRACG, MorrisonColumn PRE, MorrisonColumn PGMLT, MorrisonColumn MNUCCC,
    MorrisonColumn PSACWS, MorrisonColumn PSACWI, MorrisonColumn QMULTS, MorrisonColumn QMULTG,
    MorrisonColumn PSACWG, MorrisonColumn PGSACW, MorrisonColumn PRD, MorrisonColumn PRCI,
    MorrisonColumn PRAI, MorrisonColumn QMULTR, MorrisonColumn QMULTRG, MorrisonColumn MNUCCD,
    MorrisonColumn PRACI, MorrisonColumn PRACIS, MorrisonColumn EPRD, MorrisonColumn MNUCCR,
    MorrisonColumn PIACR, MorrisonColumn PIACRS, MorrisonColumn PGRACS, MorrisonColumn PRDS,
    MorrisonColumn EPRDS, MorrisonColumn PSACR, MorrisonColumn PRDG, MorrisonColumn EPRDG,
    MorrisonColumn NPRC1, MorrisonColumn NRAGG, MorrisonColumn NPRACG, MorrisonColumn NSUBR,
    MorrisonColumn NSMLTR, MorrisonColumn NGMLTR, MorrisonColumn NPRACS, MorrisonColumn NNUCCR,
    MorrisonColumn NIACR, MorrisonColumn NIACRS, MorrisonColumn NGRACS, MorrisonColumn NSMLTS,
    MorrisonColumn NSAGG, MorrisonColumn NPRCI, MorrisonColumn NSCNG, MorrisonColumn NSUBS,
    MorrisonColumn PCC, MorrisonColumn NNUCCC, MorrisonColumn NPSACWS, MorrisonColumn NPRA,
    MorrisonColumn NPRC, MorrisonColumn NPSACWI, MorrisonColumn NPSACWG, MorrisonColumn NPRAI,
    MorrisonColumn NMULTS, MorrisonColumn NMULTG, MorrisonColumn NMULTR, MorrisonColumn NMULTRG,
    MorrisonColumn NNUCCD, MorrisonColumn NSUBI, MorrisonColumn NGMLTG, MorrisonColumn NSUBG,
    MorrisonColumn NACT, MorrisonColumn SIZEFIX_NR, MorrisonColumn SIZEFIX_NC,
    MorrisonColumn SIZEFIX_NI, MorrisonColumn SIZEFIX_NS, MorrisonColumn SIZEFIX_NG,
    MorrisonColumn NEGFIX_NR, MorrisonColumn NEGFIX_NC, MorrisonColumn NEGFIX_NI,
    MorrisonColumn NEGFIX_NS, MorrisonColumn NEGFIX_NG, MorrisonColumn NIM_MORR_CL,
    MorrisonColumn QC_INST, MorrisonColumn QR_INST, MorrisonColumn QI_INST, MorrisonColumn QS_INST,
    MorrisonColumn QG_INST, MorrisonColumn NC_INST, MorrisonColumn NR_INST, MorrisonColumn NI_INST,
    MorrisonColumn NS_INST, MorrisonColumn NG_INST, morrison_microphysics_state& morrison);
inline float POLYSVP(float T, int TYPE, morrison_microphysics_state& morrison);
inline float GAMMA(ScratchWorkspace& workspace, float X, morrison_microphysics_state& morrison);
inline float DERF1(ScratchWorkspace& workspace, float X, morrison_microphysics_state& morrison);
inline void radar_init(morrison_microphysics_state& morrison);
inline std::complex<double> m_complex_water_ray(double lambda, double T,
                                                morrison_microphysics_state& morrison);
inline std::complex<double> m_complex_ice_maetzler(double lambda, double T,
                                                   morrison_microphysics_state& morrison);
inline void POSITIVE_QV_ADJ(float& QV, float& QC, float& QR, float& QI, float& QS, float& QG,
                            float& T_IN_K, morrison_microphysics_state& morrison);
} // namespace clubb_kokkos::morrison_engine
namespace clubb_kokkos {
inline void GRAUPEL_INIT(ScratchWorkspace& workspace, morrison_microphysics_state& morrison) {
  // Initialize source configuration, coefficients and radar tables.
  auto workspace_scope = workspace.scope();
  int n;
  int i;

  morrison.INUM = 1;
  if (morrison.dopredictNc) {
    morrison.INUM = 0;
  }
  morrison.NDCNST = morrison.Nc0;
  if ((morrison.aerosol_mode == 2)) {
    morrison.IACT = 2;
  } else {
    if ((morrison.aerosol_mode == 1)) {
      morrison.IACT = 1;
    } else {
      morrison.IACT = 0;
    }
  }
  if (morrison.docloudedgeactivation) {
    morrison.IBASE = 2;
  } else {
    morrison.IBASE = 1;
  }
  if (morrison.dosubgridw) {
    morrison.ISUB = 0;
  } else {
    morrison.ISUB = 1;
  }
  if (morrison.doicemicro) {
    morrison.ILIQ = 0;
  } else {
    morrison.ILIQ = 1;
  }
  if (morrison.doarcticicenucl) {
    morrison.INUC = 1;
  } else {
    morrison.INUC = 0;
  }
  if (morrison.dograupel) {
    morrison.IGRAUP = 0;
  } else {
    morrison.IGRAUP = 1;
  }
  if (morrison.dohail) {
    morrison.IHAIL = 1;
  } else {
    morrison.IHAIL = 0;
  }
  if (morrison.dosb_warm_rain) {
    morrison.IRAIN = 1;
  } else {
    morrison.IRAIN = 0;
  }
  morrison.ISATADJ = 0;
  morrison.AI      = 700.f;
  morrison.AC      = 3.E7f;
  morrison.AS      = 11.72f;
  morrison.AR      = 841.99667f;
  morrison.BI      = 1.f;
  morrison.BC      = 2.f;
  morrison.BS      = 0.41f;
  morrison.BR      = 0.8f;
  if ((morrison.IHAIL == 0)) {
    morrison.AG = 19.3f;
    morrison.BG = 0.37f;
  } else {
    morrison.AG = 114.5f;
    morrison.BG = 0.5f;
  }
  morrison.R     = Rd;
  morrison.RHOW  = rho_lw;
  morrison.TMELT = T_freeze_K;
  morrison.RHOSU = (85000.f / (morrison.R * morrison.TMELT));
  morrison.RHOSU = (85000.f / (morrison.R * morrison.TMELT));
  morrison.RHOW  = 997.f;
  morrison.RHOI  = 500.f;
  morrison.RHOSN = 100.f;
  if ((morrison.IHAIL == 0)) {
    morrison.RHOG = 400.f;
  } else {
    morrison.RHOG = 900.f;
  }
  morrison.AIMM = 0.66f;
  morrison.BIMM = 100.f;
  morrison.ECR  = 1.f;
  morrison.DCS  = 125.E-6f;
  morrison.MI0  = ((((4.f / 3.f) * morrison.PI) * morrison.RHOI) * ((10.E-6f * 10.E-6f) * 10.E-6f));
  morrison.MG0  = 1.6E-10f;
  morrison.F1S  = 0.86f;
  morrison.F2S  = 0.28f;
  morrison.F1R  = 0.78f;
  morrison.F2R  = 0.308f;
  morrison.G    = grav;
  morrison.QSMALL  = 1.E-14f;
  morrison.EII     = 0.1f;
  morrison.ECI     = 0.7f;
  morrison.CPW     = 4218.f;
  morrison.CI      = ((morrison.RHOI * morrison.PI) / 6.f);
  morrison.DI      = 3.f;
  morrison.CS      = ((morrison.RHOSN * morrison.PI) / 6.f);
  morrison.DS      = 3.f;
  morrison.CG      = ((morrison.RHOG * morrison.PI) / 6.f);
  morrison.DG      = 3.f;
  morrison.RIN     = 0.1E-6f;
  morrison.MMULT   = ((((4.f / 3.f) * morrison.PI) * morrison.RHOI) * ((5.E-6f * 5.E-6f) * 5.E-6f));
  morrison.LAMMAXI = (1.f / 1.E-6f);
  morrison.LAMMINI = (1.f / ((2.f * morrison.DCS) + 100.E-6f));
  morrison.LAMMAXR = (1.f / 20.E-6f);
  morrison.LAMMINR = (1.f / 2800.E-6f);
  morrison.LAMMAXS = (1.f / 10.E-6f);
  morrison.LAMMINS = (1.f / 2000.E-6f);
  morrison.LAMMAXG = (1.f / 20.E-6f);
  morrison.LAMMING = (1.f / 2000.E-6f);
  morrison.K1      = morrison.ccnexpnt;
  morrison.C1      = morrison.ccnconst;
  morrison.MW      = 0.018f;
  morrison.OSM     = 1.f;
  morrison.VI      = 3.f;
  morrison.EPSM    = 0.7f;
  morrison.RHOA    = 1777.f;
  morrison.MAP     = 0.132f;
  morrison.MA      = 0.0284f;
  morrison.RR      = 8.3187f;
  morrison.BACT =
      (((((morrison.VI * morrison.OSM) * morrison.EPSM) * morrison.MW) * morrison.RHOA) /
       (morrison.MAP * morrison.RHOW));
  morrison.RM1    = morrison.aer_rm1;
  morrison.SIG1   = morrison.aer_sig1;
  morrison.NANEW1 = morrison.aer_n1;
  morrison.F11 =
      (0.5f * Kokkos::exp((2.5f * (Kokkos::log(morrison.SIG1) * Kokkos::log(morrison.SIG1)))));
  morrison.F21    = (1.f + (0.25f * Kokkos::log(morrison.SIG1)));
  morrison.RM2    = morrison.aer_rm2;
  morrison.SIG2   = morrison.aer_sig2;
  morrison.NANEW2 = morrison.aer_n2;
  morrison.F12 =
      (0.5f * Kokkos::exp((2.5f * (Kokkos::log(morrison.SIG2) * Kokkos::log(morrison.SIG2)))));
  morrison.F22    = (1.f + (0.25f * Kokkos::log(morrison.SIG2)));
  morrison.CONS1 = (morrison_engine::GAMMA(workspace, (1.f + morrison.DS), morrison) * morrison.CS);
  morrison.CONS2 = (morrison_engine::GAMMA(workspace, (1.f + morrison.DG), morrison) * morrison.CG);
  morrison.CONS3 = (morrison_engine::GAMMA(workspace, (4.f + morrison.BS), morrison) / 6.f);
  morrison.CONS4 = (morrison_engine::GAMMA(workspace, (4.f + morrison.BR), morrison) / 6.f);
  morrison.CONS5 = morrison_engine::GAMMA(workspace, (1.f + morrison.BS), morrison);
  morrison.CONS6 = morrison_engine::GAMMA(workspace, (1.f + morrison.BR), morrison);
  morrison.CONS7 = (morrison_engine::GAMMA(workspace, (4.f + morrison.BG), morrison) / 6.f);
  morrison.CONS8 = morrison_engine::GAMMA(workspace, (1.f + morrison.BG), morrison);
  morrison.CONS9 = morrison_engine::GAMMA(workspace, ((5.f / 2.f) + (morrison.BR / 2.f)), morrison);
  morrison.CONS10 =
      morrison_engine::GAMMA(workspace, ((5.f / 2.f) + (morrison.BS / 2.f)), morrison);
  morrison.CONS11 =
      morrison_engine::GAMMA(workspace, ((5.f / 2.f) + (morrison.BG / 2.f)), morrison);
  morrison.CONS12 =
      (morrison_engine::GAMMA(workspace, (1.f + morrison.DI), morrison) * morrison.CI);
  morrison.CONS13 =
      (((morrison_engine::GAMMA(workspace, (morrison.BS + 3.f), morrison) * morrison.PI) / 4.f) *
       morrison.ECI);
  morrison.CONS14 =
      (((morrison_engine::GAMMA(workspace, (morrison.BG + 3.f), morrison) * morrison.PI) / 4.f) *
       morrison.ECI);
  morrison.CONS15 =
      (((((-1108.f) * morrison.EII) * Kokkos::pow(morrison.PI, ((1.f - morrison.BS) / 3.f))) *
        Kokkos::pow(morrison.RHOSN, (((-2.f) - morrison.BS) / 3.f))) /
       (4.f * 720.f));
  morrison.CONS16 =
      (((morrison_engine::GAMMA(workspace, (morrison.BI + 3.f), morrison) * morrison.PI) / 4.f) *
       morrison.ECI);
  morrison.CONS17 =
      ((((((((4.f * 2.f) * 3.f) * morrison.RHOSU) * morrison.PI) * morrison.ECI) * morrison.ECI) *
        morrison_engine::GAMMA(workspace, ((2.f * morrison.BS) + 2.f), morrison)) /
       (8.f * (morrison.RHOG - morrison.RHOSN)));
  morrison.CONS18 = (morrison.RHOSN * morrison.RHOSN);
  morrison.CONS19 = (morrison.RHOW * morrison.RHOW);
  morrison.CONS20 = ((((20.f * morrison.PI) * morrison.PI) * morrison.RHOW) * morrison.BIMM);
  morrison.CONS21 = (4.f / (morrison.DCS * morrison.RHOI));
  morrison.CONS22 =
      (((morrison.PI * morrison.RHOI) * ((morrison.DCS * morrison.DCS) * morrison.DCS)) / 6.f);
  morrison.CONS23 = (((morrison.PI / 4.f) * morrison.EII) *
                     morrison_engine::GAMMA(workspace, (morrison.BS + 3.f), morrison));
  morrison.CONS24 = (((morrison.PI / 4.f) * morrison.ECR) *
                     morrison_engine::GAMMA(workspace, (morrison.BR + 3.f), morrison));
  morrison.CONS25 = (((((morrison.PI * morrison.PI) / 24.f) * morrison.RHOW) * morrison.ECR) *
                     morrison_engine::GAMMA(workspace, (morrison.BR + 6.f), morrison));
  morrison.CONS26 = ((morrison.PI / 6.f) * morrison.RHOW);
  morrison.CONS27 = morrison_engine::GAMMA(workspace, (1.f + morrison.BI), morrison);
  morrison.CONS28 = (morrison_engine::GAMMA(workspace, (4.f + morrison.BI), morrison) / 6.f);
  morrison.CONS29 =
      ((((4.f / 3.f) * morrison.PI) * morrison.RHOW) * ((25.E-6f * 25.E-6f) * 25.E-6f));
  morrison.CONS30 = (((4.f / 3.f) * morrison.PI) * morrison.RHOW);
  morrison.CONS31 = (((morrison.PI * morrison.PI) * morrison.ECR) * morrison.RHOSN);
  morrison.CONS32 = ((morrison.PI / 2.f) * morrison.ECR);
  morrison.CONS33 = (((morrison.PI * morrison.PI) * morrison.ECR) * morrison.RHOG);
  morrison.CONS34 = ((5.f / 2.f) + (morrison.BR / 2.f));
  morrison.CONS35 = ((5.f / 2.f) + (morrison.BS / 2.f));
  morrison.CONS36 = ((5.f / 2.f) + (morrison.BG / 2.f));
  morrison.CONS37 = (((4.f * morrison.PI) * 1.38E-23f) / ((6.f * morrison.PI) * morrison.RIN));
  morrison.CONS38 = (((morrison.PI * morrison.PI) / 3.f) * morrison.RHOW);
  morrison.CONS39 = ((((morrison.PI * morrison.PI) / 36.f) * morrison.RHOW) * morrison.BIMM);
  morrison.CONS40 = ((morrison.PI / 6.f) * morrison.BIMM);
  morrison.CONS41 = (((morrison.PI * morrison.PI) * morrison.ECR) * morrison.RHOW);
  morrison.dnu[0]            = (-0.557f);
  morrison.dnu[1]            = (-0.557f);
  morrison.dnu[2]            = (-0.430f);
  morrison.dnu[3]            = (-0.307f);
  morrison.dnu[4]            = (-0.186f);
  morrison.dnu[5]            = (-0.067f);
  morrison.dnu[6]            = 0.050f;
  morrison.dnu[7]            = 0.167f;
  morrison.dnu[8]            = 0.282f;
  morrison.dnu[9]            = 0.397f;
  morrison.dnu[10]           = 0.512f;
  morrison.dnu[11]           = 0.626f;
  morrison.dnu[12]           = 0.739f;
  morrison.dnu[13]           = 0.853f;
  morrison.dnu[14]           = 0.966f;
  morrison.dnu[15]           = 0.966f;
  morrison.ddx[0]            = (morrison.D0r * 1.0e0);
  morrison.ddx[morrison.nbr] = 0.005e0;
  for (n = 1; n < morrison.nbr; n += 1) {
    morrison.ddx[n] = Kokkos::exp((((static_cast<double>(n) / static_cast<double>(morrison.nbr)) *
                                    Kokkos::log((morrison.ddx[morrison.nbr] / morrison.ddx[0]))) +
                                   Kokkos::log(morrison.ddx[0])));
  }
  for (n = 0; n < morrison.nbr; n += 1) {
    morrison.Dr[n]  = Kokkos::sqrt((morrison.ddx[n] * morrison.ddx[((n + 1))]));
    morrison.dtr[n] = (morrison.ddx[((n + 1))] - morrison.ddx[n]);
  }
  morrison.ddx[0]            = (morrison.D0s * 1.0e0);
  morrison.ddx[morrison.nbs] = 0.02e0;
  for (n = 1; n < morrison.nbs; n += 1) {
    morrison.ddx[n] = Kokkos::exp((((static_cast<double>(n) / static_cast<double>(morrison.nbs)) *
                                    Kokkos::log((morrison.ddx[morrison.nbs] / morrison.ddx[0]))) +
                                   Kokkos::log(morrison.ddx[0])));
  }
  for (n = 0; n < morrison.nbs; n += 1) {
    morrison.Dds[n] = Kokkos::sqrt((morrison.ddx[n] * morrison.ddx[((n + 1))]));
    morrison.dts[n] = (morrison.ddx[((n + 1))] - morrison.ddx[n]);
  }
  morrison.ddx[0]            = (morrison.D0g * 1.0e0);
  morrison.ddx[morrison.nbg] = 0.05e0;
  for (n = 1; n < morrison.nbg; n += 1) {
    morrison.ddx[n] = Kokkos::exp((((static_cast<double>(n) / static_cast<double>(morrison.nbg)) *
                                    Kokkos::log((morrison.ddx[morrison.nbg] / morrison.ddx[0]))) +
                                   Kokkos::log(morrison.ddx[0])));
  }
  for (n = 0; n < morrison.nbg; n += 1) {
    morrison.Ddg[n] = Kokkos::sqrt((morrison.ddx[n] * morrison.ddx[((n + 1))]));
    morrison.dtg[n] = (morrison.ddx[((n + 1))] - morrison.ddx[n]);
  }
  for (i = 0; i < 256; i += 1) {
    morrison.mp_debug.replace(i, 1, std::string(1, static_cast<char>(0)));
  }
  morrison_engine::radar_init(morrison);
}
} // namespace clubb_kokkos
namespace clubb_kokkos::morrison_engine {
inline void M2005MICRO_GRAUPEL(
    ScratchWorkspace& workspace, MorrisonColumn QC3DTEN, MorrisonColumn QI3DTEN,
    MorrisonColumn QNI3DTEN, MorrisonColumn QR3DTEN, MorrisonColumn NC3DTEN, MorrisonColumn NI3DTEN,
    MorrisonColumn NS3DTEN, MorrisonColumn NR3DTEN, MorrisonColumn QC3D, MorrisonColumn QI3D,
    MorrisonColumn QNI3D, MorrisonColumn QR3D, MorrisonColumn NC3D, MorrisonColumn NI3D,
    MorrisonColumn NS3D, MorrisonColumn NR3D, MorrisonColumn T3DTEN, MorrisonColumn QV3DTEN,
    MorrisonColumn T3D, MorrisonColumn QV3D, MorrisonColumn PRES, MorrisonColumn RHO,
    MorrisonColumn DZQ, MorrisonColumn W3D, MorrisonColumn WVAR, MorrisonColumn FR, float& PRECRT,
    float& SNOWRT, MorrisonColumn EFFC, MorrisonColumn EFFI, MorrisonColumn EFFS,
    MorrisonColumn EFFR, float DT, int IMS, int IME, int JMS, int JME, int KMS, int KME, int ITS,
    int ITE, int JTS, int JTE, int KTS, int KTE, MorrisonColumn QG3DTEN, MorrisonColumn NG3DTEN,
    MorrisonColumn QG3D, MorrisonColumn NG3D, MorrisonColumn EFFG, MorrisonColumn QGSTEN,
    MorrisonColumn QRSTEN, MorrisonColumn QISTEN, MorrisonColumn QNISTEN, MorrisonColumn QCSTEN,
    MorrisonColumn NGSTEN, MorrisonColumn NRSTEN, MorrisonColumn NISTEN, MorrisonColumn NSSTEN,
    MorrisonColumn NCSTEN, MorrisonColumn CF3D, MorrisonColumn PRC, MorrisonColumn PRA,
    MorrisonColumn PSMLT, MorrisonColumn EVPMS, MorrisonColumn PRACS, MorrisonColumn EVPMG,
    MorrisonColumn PRACG, MorrisonColumn PRE, MorrisonColumn PGMLT, MorrisonColumn MNUCCC,
    MorrisonColumn PSACWS, MorrisonColumn PSACWI, MorrisonColumn QMULTS, MorrisonColumn QMULTG,
    MorrisonColumn PSACWG, MorrisonColumn PGSACW, MorrisonColumn PRD, MorrisonColumn PRCI,
    MorrisonColumn PRAI, MorrisonColumn QMULTR, MorrisonColumn QMULTRG, MorrisonColumn MNUCCD,
    MorrisonColumn PRACI, MorrisonColumn PRACIS, MorrisonColumn EPRD, MorrisonColumn MNUCCR,
    MorrisonColumn PIACR, MorrisonColumn PIACRS, MorrisonColumn PGRACS, MorrisonColumn PRDS,
    MorrisonColumn EPRDS, MorrisonColumn PSACR, MorrisonColumn PRDG, MorrisonColumn EPRDG,
    MorrisonColumn NPRC1, MorrisonColumn NRAGG, MorrisonColumn NPRACG, MorrisonColumn NSUBR,
    MorrisonColumn NSMLTR, MorrisonColumn NGMLTR, MorrisonColumn NPRACS, MorrisonColumn NNUCCR,
    MorrisonColumn NIACR, MorrisonColumn NIACRS, MorrisonColumn NGRACS, MorrisonColumn NSMLTS,
    MorrisonColumn NSAGG, MorrisonColumn NPRCI, MorrisonColumn NSCNG, MorrisonColumn NSUBS,
    MorrisonColumn PCC, MorrisonColumn NNUCCC, MorrisonColumn NPSACWS, MorrisonColumn NPRA,
    MorrisonColumn NPRC, MorrisonColumn NPSACWI, MorrisonColumn NPSACWG, MorrisonColumn NPRAI,
    MorrisonColumn NMULTS, MorrisonColumn NMULTG, MorrisonColumn NMULTR, MorrisonColumn NMULTRG,
    MorrisonColumn NNUCCD, MorrisonColumn NSUBI, MorrisonColumn NGMLTG, MorrisonColumn NSUBG,
    MorrisonColumn NACT, MorrisonColumn SIZEFIX_NR, MorrisonColumn SIZEFIX_NC,
    MorrisonColumn SIZEFIX_NI, MorrisonColumn SIZEFIX_NS, MorrisonColumn SIZEFIX_NG,
    MorrisonColumn NEGFIX_NR, MorrisonColumn NEGFIX_NC, MorrisonColumn NEGFIX_NI,
    MorrisonColumn NEGFIX_NS, MorrisonColumn NEGFIX_NG, MorrisonColumn NIM_MORR_CL,
    MorrisonColumn QC_INST, MorrisonColumn QR_INST, MorrisonColumn QI_INST, MorrisonColumn QS_INST,
    MorrisonColumn QG_INST, MorrisonColumn NC_INST, MorrisonColumn NR_INST, MorrisonColumn NI_INST,
    MorrisonColumn NS_INST, MorrisonColumn NG_INST, morrison_microphysics_state& morrison) {
  // Execute the Native M2005MICRO_GRAUPEL body and update its declared outputs.
  auto workspace_scope = workspace.scope();
  HostFloat1D LAMC = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::LAMC", (KME - KMS + 1));
  HostFloat1D LAMI = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::LAMI", (KME - KMS + 1));
  HostFloat1D LAMS = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::LAMS", (KME - KMS + 1));
  HostFloat1D LAMR = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::LAMR", (KME - KMS + 1));
  HostFloat1D LAMG = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::LAMG", (KME - KMS + 1));
  HostFloat1D CDIST1 =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::CDIST1", (KME - KMS + 1));
  HostFloat1D N0I  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::N0I", (KME - KMS + 1));
  HostFloat1D N0S  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::N0S", (KME - KMS + 1));
  HostFloat1D N0RR = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::N0RR", (KME - KMS + 1));
  HostFloat1D N0G  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::N0G", (KME - KMS + 1));
  HostFloat1D PGAM = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::PGAM", (KME - KMS + 1));
  HostFloat1D NSUBC =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NSUBC", (KME - KMS + 1));
  HostFloat1D PCCN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::PCCN", (KME - KMS + 1));
  HostFloat1D KAP  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::KAP", (KME - KMS + 1));
  HostFloat1D EVS  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::EVS", (KME - KMS + 1));
  HostFloat1D EIS  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::EIS", (KME - KMS + 1));
  HostFloat1D QVS  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QVS", (KME - KMS + 1));
  HostFloat1D QVI  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QVI", (KME - KMS + 1));
  HostFloat1D QVQVS =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QVQVS", (KME - KMS + 1));
  HostFloat1D QVQVSI =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QVQVSI", (KME - KMS + 1));
  HostFloat1D DV   = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DV", (KME - KMS + 1));
  HostFloat1D XXLS = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::XXLS", (KME - KMS + 1));
  HostFloat1D XXLV = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::XXLV", (KME - KMS + 1));
  HostFloat1D CPM  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::CPM", (KME - KMS + 1));
  HostFloat1D MU   = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::MU", (KME - KMS + 1));
  HostFloat1D SC   = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::SC", (KME - KMS + 1));
  HostFloat1D XLF  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::XLF", (KME - KMS + 1));
  HostFloat1D AB   = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::AB", (KME - KMS + 1));
  HostFloat1D ABI  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::ABI", (KME - KMS + 1));
  HostFloat1D DAP  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DAP", (KME - KMS + 1));
  float NACNT;
  float FMULT;
  float COFFI;
  HostFloat1D DUMI = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMI", (KME - KMS + 1));
  HostFloat1D DUMR = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMR", (KME - KMS + 1));
  HostFloat1D DUMFNI =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMFNI", (KME - KMS + 1));
  HostFloat1D DUMG = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMG", (KME - KMS + 1));
  HostFloat1D DUMFNG =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMFNG", (KME - KMS + 1));
  float UNI;
  float UMI;
  float UMR;
  HostFloat1D FI  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FI", (KME - KMS + 1));
  HostFloat1D FNI = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FNI", (KME - KMS + 1));
  HostFloat1D FG  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FG", (KME - KMS + 1));
  HostFloat1D FNG = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FNG", (KME - KMS + 1));
  float RGVM;
  HostFloat1D FALOUTR =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTR", (KME - KMS + 1));
  HostFloat1D FALOUTI =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTI", (KME - KMS + 1));
  HostFloat1D FALOUTNI =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTNI", (KME - KMS + 1));
  float FALTNDR;
  float FALTNDI;
  float FALTNDNI;
  float RHO2;
  HostFloat1D DUMQS =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMQS", (KME - KMS + 1));
  HostFloat1D DUMFNS =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMFNS", (KME - KMS + 1));
  float UMS;
  float UNS;
  HostFloat1D FS  = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FS", (KME - KMS + 1));
  HostFloat1D FNS = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FNS", (KME - KMS + 1));
  HostFloat1D FALOUTS =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTS", (KME - KMS + 1));
  HostFloat1D FALOUTNS =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTNS", (KME - KMS + 1));
  HostFloat1D FALOUTG =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTG", (KME - KMS + 1));
  HostFloat1D FALOUTNG =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTNG", (KME - KMS + 1));
  float FALTNDS;
  float FALTNDNS;
  float UNR;
  float FALTNDG;
  float FALTNDNG;
  HostFloat1D DUMC = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMC", (KME - KMS + 1));
  HostFloat1D DUMFNC =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMFNC", (KME - KMS + 1));
  float UNC;
  float UMC;
  float UNG;
  float UMG;
  HostFloat1D FC = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FC", (KME - KMS + 1));
  HostFloat1D FALOUTC =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTC", (KME - KMS + 1));
  HostFloat1D FALOUTNC =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTNC", (KME - KMS + 1));
  float FALTNDC;
  float FALTNDNC;
  HostFloat1D FNC = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FNC", (KME - KMS + 1));
  HostFloat1D DUMFNR =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::DUMFNR", (KME - KMS + 1));
  HostFloat1D FALOUTNR =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FALOUTNR", (KME - KMS + 1));
  float FALTNDNR;
  HostFloat1D FNR = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::FNR", (KME - KMS + 1));
  HostFloat1D AIN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::AIN", (KME - KMS + 1));
  HostFloat1D ARN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::ARN", (KME - KMS + 1));
  HostFloat1D ASN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::ASN", (KME - KMS + 1));
  HostFloat1D ACN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::ACN", (KME - KMS + 1));
  HostFloat1D AGN = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::AGN", (KME - KMS + 1));
  float DUM;
  float DUM1;
  float DUM2;
  float DUMT;
  float DUMQV;
  float DUMQSS;
  float DUMQSI;
  float DUMS;
  float TMPNUM;
  float DQSDT;
  float DQSIDT;
  float EPSI;
  float EPSS;
  float EPSR;
  float EPSG;
  float TAUC;
  float TAUR;
  float TAUI;
  float TAUS;
  float TAUG;
  float DUMACT;
  float DUM3;
  int K;
  int NSTEP;
  int N;
  int LTRUE;
  float CT;
  float TEMP1;
  float SAT1;
  float SIGVL;
  float KEL;
  float KC2;
  float CRY;
  float KRY;
  float DUMQI;
  float DUMNI;
  float DC0;
  float DS0;
  float DG0;
  float DUMQC;
  float DUMQR;
  float RATIO;
  float SUM_DEP;
  float FUDGEF;
  float WEF;
  float ANUC;
  float BNUC;
  float AACT;
  float GAMM;
  float GG;
  float PSI;
  float ETA1;
  float ETA2;
  float SM1;
  float SM2;
  float SMAX;
  float UU1;
  float UU2;
  float ALPHA;
  float DLAMS;
  float DLAMR;
  float DLAMI;
  float DLAMC;
  float DLAMG;
  float LAMMAX;
  float LAMMIN;
  int IDROP;
  HostFloat1D nu = workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::nu", (KME - KMS + 1));
  int dumii;
  HostFloat1D QC3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QC3D_INIT", (KME - KMS + 1));
  HostFloat1D QI3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QI3D_INIT", (KME - KMS + 1));
  HostFloat1D QNI3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QNI3D_INIT", (KME - KMS + 1));
  HostFloat1D QR3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QR3D_INIT", (KME - KMS + 1));
  HostFloat1D NC3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NC3D_INIT", (KME - KMS + 1));
  HostFloat1D NI3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NI3D_INIT", (KME - KMS + 1));
  HostFloat1D NS3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NS3D_INIT", (KME - KMS + 1));
  HostFloat1D NR3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NR3D_INIT", (KME - KMS + 1));
  HostFloat1D T3D_TEMP =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::T3D_TEMP", (KME - KMS + 1));
  HostFloat1D QV3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QV3D_INIT", (KME - KMS + 1));
  HostFloat1D QG3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::QG3D_INIT", (KME - KMS + 1));
  HostFloat1D NG3D_INIT =
      workspace_scope.view<HostFloat1D>("M2005MICRO_GRAUPEL::NG3D_INIT", (KME - KMS + 1));
  float QV_INIT;
  float QSAT_INIT;
  float TMPQSMALL;
  float T3D_INIT;

  LTRUE = 0;
  for (int index_0 = 0; index_0 < KTE - KTS + 1; ++index_0) {
    EFFC(((KTS + index_0))) = 25.f;
  }
  for (int index_0 = 0; index_0 < KTE - KTS + 1; ++index_0) {
    EFFI(((KTS + index_0))) = 25.f;
  }
  for (int index_0 = 0; index_0 < KTE - KTS + 1; ++index_0) {
    EFFS(((KTS + index_0))) = 25.f;
  }
  for (int index_0 = 0; index_0 < KTE - KTS + 1; ++index_0) {
    EFFR(((KTS + index_0))) = 25.f;
  }
  for (int index_0 = 0; index_0 < KTE - KTS + 1; ++index_0) {
    EFFG(((KTS + index_0))) = 25.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PRC(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NPRC(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NPRC1(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PRA(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NPRA(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NRAGG(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PSMLT(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NSMLTS(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NSMLTR(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PCC(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PRE(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NSUBC(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NSUBR(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PRACG(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NPRACG(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PSMLT(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    EVPMS(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PGMLT(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    EVPMG(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    PRACS(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NPRACS(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NGMLTG(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    NGMLTR(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < KME - KMS + 1; ++index_0) {
    FR(((KMS + index_0))) = 0.f;
  }
  for (int index_0 = 0; index_0 < NPRC1.extent_int(0); ++index_0) {
    NPRC1(index_0) = 0.0f;
  }
  for (int index_0 = 0; index_0 < NRAGG.extent_int(0); ++index_0) {
    NRAGG(index_0) = 0.0f;
  }
  for (int index_0 = 0; index_0 < NPRACG.extent_int(0); ++index_0) {
    NPRACG(index_0) = 0.0f;
  }
  for (int index_0 = 0; index_0 < NSUBR.extent_int(0); ++index_0) {
    NSUBR(index_0) = 0.0f;
  }
  for (int index_0 = 0; index_0 < NSMLTR.extent_int(0); ++index_0) {
    NSMLTR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NGMLTR.extent_int(0); ++index_0) {
    NGMLTR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NPRACS.extent_int(0); ++index_0) {
    NPRACS(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NNUCCR.extent_int(0); ++index_0) {
    NNUCCR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NIACR.extent_int(0); ++index_0) {
    NIACR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NIACRS.extent_int(0); ++index_0) {
    NIACRS(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NGRACS.extent_int(0); ++index_0) {
    NGRACS(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NACT.extent_int(0); ++index_0) {
    NACT(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < SIZEFIX_NR.extent_int(0); ++index_0) {
    SIZEFIX_NR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < SIZEFIX_NC.extent_int(0); ++index_0) {
    SIZEFIX_NC(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < SIZEFIX_NI.extent_int(0); ++index_0) {
    SIZEFIX_NI(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < SIZEFIX_NS.extent_int(0); ++index_0) {
    SIZEFIX_NS(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < SIZEFIX_NG.extent_int(0); ++index_0) {
    SIZEFIX_NG(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NEGFIX_NI.extent_int(0); ++index_0) {
    NEGFIX_NI(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NEGFIX_NS.extent_int(0); ++index_0) {
    NEGFIX_NS(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NEGFIX_NC.extent_int(0); ++index_0) {
    NEGFIX_NC(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NEGFIX_NR.extent_int(0); ++index_0) {
    NEGFIX_NR(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NEGFIX_NG.extent_int(0); ++index_0) {
    NEGFIX_NG(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NIM_MORR_CL.extent_int(0); ++index_0) {
    NIM_MORR_CL(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QC_INST.extent_int(0); ++index_0) {
    QC_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QR_INST.extent_int(0); ++index_0) {
    QR_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QI_INST.extent_int(0); ++index_0) {
    QI_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QS_INST.extent_int(0); ++index_0) {
    QS_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QG_INST.extent_int(0); ++index_0) {
    QG_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NC_INST.extent_int(0); ++index_0) {
    NC_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NR_INST.extent_int(0); ++index_0) {
    NR_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NI_INST.extent_int(0); ++index_0) {
    NI_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NS_INST.extent_int(0); ++index_0) {
    NS_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < NG_INST.extent_int(0); ++index_0) {
    NG_INST(index_0) = 0.f;
  }
  for (int index_0 = 0; index_0 < QC3D_INIT.extent_int(0); ++index_0) {
    QC3D_INIT(index_0) = QC3D(index_0);
  }
  for (int index_0 = 0; index_0 < QI3D_INIT.extent_int(0); ++index_0) {
    QI3D_INIT(index_0) = QI3D(index_0);
  }
  for (int index_0 = 0; index_0 < QNI3D_INIT.extent_int(0); ++index_0) {
    QNI3D_INIT(index_0) = QNI3D(index_0);
  }
  for (int index_0 = 0; index_0 < QR3D_INIT.extent_int(0); ++index_0) {
    QR3D_INIT(index_0) = QR3D(index_0);
  }
  for (int index_0 = 0; index_0 < NC3D_INIT.extent_int(0); ++index_0) {
    NC3D_INIT(index_0) = NC3D(index_0);
  }
  for (int index_0 = 0; index_0 < NI3D_INIT.extent_int(0); ++index_0) {
    NI3D_INIT(index_0) = NI3D(index_0);
  }
  for (int index_0 = 0; index_0 < NS3D_INIT.extent_int(0); ++index_0) {
    NS3D_INIT(index_0) = NS3D(index_0);
  }
  for (int index_0 = 0; index_0 < NR3D_INIT.extent_int(0); ++index_0) {
    NR3D_INIT(index_0) = NR3D(index_0);
  }
  for (int index_0 = 0; index_0 < QG3D_INIT.extent_int(0); ++index_0) {
    QG3D_INIT(index_0) = QG3D(index_0);
  }
  for (int index_0 = 0; index_0 < NG3D_INIT.extent_int(0); ++index_0) {
    NG3D_INIT(index_0) = NG3D(index_0);
  }
  for (int index_0 = 0; index_0 < T3D_TEMP.extent_int(0); ++index_0) {
    T3D_TEMP(index_0) = T3D(index_0);
  }
  for (int index_0 = 0; index_0 < QV3D_INIT.extent_int(0); ++index_0) {
    QV3D_INIT(index_0) = QV3D(index_0);
  }
  for (K = KTS; K < KTE + 1; K += 1) {
    for (int index_0 = 0; index_0 < XXLV.extent_int(0); ++index_0) {
      XXLV(index_0) = Lv;
    }
    XXLS(((K - KMS))) = Ls;
    CPM(((K - KMS)))  = Cp;
    EVS(((K - KMS))) =
        (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(T3D(((K - KMS))), 0, morrison))
             ? ((0.99f * PRES(((K - KMS)))))
             : (morrison_engine::POLYSVP(T3D(((K - KMS))), 0, morrison)));
    EIS(((K - KMS))) =
        (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(T3D(((K - KMS))), 1, morrison))
             ? ((0.99f * PRES(((K - KMS)))))
             : (morrison_engine::POLYSVP(T3D(((K - KMS))), 1, morrison)));
    if ((EIS(((K - KMS))) > EVS(((K - KMS))))) {
      EIS(((K - KMS))) = EVS(((K - KMS)));
    }
    QVS(((K - KMS))) = ((ep * EVS(((K - KMS)))) / (PRES(((K - KMS))) - EVS(((K - KMS)))));
    QVI(((K - KMS))) = ((ep * EIS(((K - KMS)))) / (PRES(((K - KMS))) - EIS(((K - KMS)))));
    if ((CF3D(((K - KMS))) > morrison.cloud_frac_thresh)) {
      T3D_INIT          = T3D(((K - KMS)));
      QV_INIT           = QV3D(((K - KMS)));
      QV3D(((K - KMS))) = QVS(((K - KMS)));
      QSAT_INIT         = QVS(((K - KMS)));
      QC3D(((K - KMS))) = (QC3D(((K - KMS))) / CF3D(((K - KMS))));
      if ((morrison.INUM == 0)) {
        NC3D(((K - KMS))) = (NC3D(((K - KMS))) / CF3D(((K - KMS))));
      }
      QR3D(((K - KMS))) = (QR3D(((K - KMS))) / CF3D(((K - KMS))));
      NR3D(((K - KMS))) = (NR3D(((K - KMS))) / CF3D(((K - KMS))));
      if ((morrison.ILIQ == 0)) {
        QI3D(((K - KMS)))  = (QI3D(((K - KMS))) / CF3D(((K - KMS))));
        NI3D(((K - KMS)))  = (NI3D(((K - KMS))) / CF3D(((K - KMS))));
        QNI3D(((K - KMS))) = (QNI3D(((K - KMS))) / CF3D(((K - KMS))));
        NS3D(((K - KMS)))  = (NS3D(((K - KMS))) / CF3D(((K - KMS))));
      }
      if ((morrison.IGRAUP == 0)) {
        QG3D(((K - KMS))) = (QG3D(((K - KMS))) / CF3D(((K - KMS))));
        NG3D(((K - KMS))) = (NG3D(((K - KMS))) / CF3D(((K - KMS))));
      }
    }
    QVQVS(((K - KMS)))  = (QV3D(((K - KMS))) / QVS(((K - KMS))));
    QVQVSI(((K - KMS))) = (QV3D(((K - KMS))) / QVI(((K - KMS))));
    if ((QVQVS(((K - KMS))) < 0.9f)) {
      if ((QR3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QR3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QR3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
        QR3D(((K - KMS)))    = 0.f;
      }
      if ((QC3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QC3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QC3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QC_INST(((K - KMS))) = (QC_INST(((K - KMS))) - (QC3D(((K - KMS))) / DT));
        QC3D(((K - KMS)))    = 0.f;
      }
    }
    if ((QVQVSI(((K - KMS))) < 0.9f)) {
      if ((QI3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QI3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) - (QI3D(((K - KMS))) / DT));
        QI3D(((K - KMS)))    = 0.f;
      }
      if ((QNI3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QNI3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QNI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) - (QNI3D(((K - KMS))) / DT));
        QNI3D(((K - KMS)))   = 0.f;
      }
      if ((QG3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QG3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QG3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) - (QG3D(((K - KMS))) / DT));
        QG3D(((K - KMS)))    = 0.f;
      }
    }
    XLF(((K - KMS))) = (XXLS(((K - KMS))) - XXLV(((K - KMS))));
    if ((QC3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QC3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QC3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
      NC_INST(((K - KMS))) = (NC_INST(((K - KMS))) - (NC3D(((K - KMS))) / DT));
      QC_INST(((K - KMS))) = (QC_INST(((K - KMS))) - (QC3D(((K - KMS))) / DT));
      QC3D(((K - KMS)))    = 0.f;
      NC3D(((K - KMS)))    = 0.f;
      EFFC(((K - KMS)))    = 0.f;
    }
    if ((QR3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QR3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QR3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
      NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) - (NR3D(((K - KMS))) / DT));
      QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
      QR3D(((K - KMS)))    = 0.f;
      NR3D(((K - KMS)))    = 0.f;
      EFFR(((K - KMS)))    = 0.f;
    }
    if ((QI3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QI3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NI_INST(((K - KMS))) = (NI_INST(((K - KMS))) - (NI3D(((K - KMS))) / DT));
      QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) - (QI3D(((K - KMS))) / DT));
      QI3D(((K - KMS)))    = 0.f;
      NI3D(((K - KMS)))    = 0.f;
      EFFI(((K - KMS)))    = 0.f;
    }
    if ((QNI3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QNI3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QNI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NS_INST(((K - KMS))) = (NS_INST(((K - KMS))) - (NS3D(((K - KMS))) / DT));
      QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) - (QNI3D(((K - KMS))) / DT));
      QNI3D(((K - KMS)))   = 0.f;
      NS3D(((K - KMS)))    = 0.f;
      EFFS(((K - KMS)))    = 0.f;
    }
    if ((QG3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QG3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QG3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NG_INST(((K - KMS))) = (NG_INST(((K - KMS))) - (NG3D(((K - KMS))) / DT));
      QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) - (QG3D(((K - KMS))) / DT));
      QG3D(((K - KMS)))    = 0.f;
      NG3D(((K - KMS)))    = 0.f;
      EFFG(((K - KMS)))    = 0.f;
    }
    QRSTEN(((K - KMS)))  = 0.f;
    QISTEN(((K - KMS)))  = 0.f;
    QNISTEN(((K - KMS))) = 0.f;
    QCSTEN(((K - KMS)))  = 0.f;
    QGSTEN(((K - KMS)))  = 0.f;
    NRSTEN(((K - KMS)))  = 0.f;
    NISTEN(((K - KMS)))  = 0.f;
    NSSTEN(((K - KMS)))  = 0.f;
    NCSTEN(((K - KMS)))  = 0.f;
    NGSTEN(((K - KMS)))  = 0.f;
    MU(((K - KMS))) =
        ((1.496E-6f * Kokkos::pow(T3D(((K - KMS))), 1.5f)) / (T3D(((K - KMS))) + 120.f));
    DUM               = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
    AIN(((K - KMS)))  = (Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.35f) * morrison.AI);
    ARN(((K - KMS)))  = (DUM * morrison.AR);
    ASN(((K - KMS)))  = (DUM * morrison.AS);
    ACN(((K - KMS)))  = ((morrison.G * morrison.RHOW) / (18.f * MU(((K - KMS)))));
    AGN(((K - KMS)))  = (DUM * morrison.AG);
    LAMI(((K - KMS))) = 0.f;
    if ((((((QC3D(((K - KMS))) < morrison.QSMALL) && (QI3D(((K - KMS))) < morrison.QSMALL)) &&
           (QNI3D(((K - KMS))) < morrison.QSMALL)) &&
          (QR3D(((K - KMS))) < morrison.QSMALL)) &&
         (QG3D(((K - KMS))) < morrison.QSMALL))) {
      if (((T3D(((K - KMS))) < morrison.TMELT) && (QVQVSI(((K - KMS))) < 0.999f))) {
        goto label_200;
      }
      if (((T3D(((K - KMS))) >= morrison.TMELT) && (QVQVS(((K - KMS))) < 0.999f))) {
        goto label_200;
      }
    }
    KAP(((K - KMS))) = (1.414E3f * MU(((K - KMS))));
    DV(((K - KMS)))  = ((8.794E-5f * Kokkos::pow(T3D(((K - KMS))), 1.81f)) / PRES(((K - KMS))));
    SC(((K - KMS)))  = (MU(((K - KMS))) / (RHO(((K - KMS))) * DV(((K - KMS)))));
    DUM              = (Rv * (T3D(((K - KMS))) * T3D(((K - KMS)))));
    DQSDT            = ((XXLV(((K - KMS))) * QVS(((K - KMS)))) / DUM);
    DQSIDT           = ((XXLS(((K - KMS))) * QVI(((K - KMS)))) / DUM);
    ABI(((K - KMS))) = (1.f + ((DQSIDT * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
    AB(((K - KMS)))  = (1.f + ((DQSDT * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
    if ((T3D(((K - KMS))) >= morrison.TMELT)) {
      if ((morrison.INUM == 1)) {
        NC3D(((K - KMS))) = ((morrison.NDCNST * 1.E6f) / RHO(((K - KMS))));
      }
      if ((QNI3D(((K - KMS))) < 1.E-6f)) {
        QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) + (QNI3D(((K - KMS))) / DT));
        QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) - (QNI3D(((K - KMS))) / DT));
        NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) + (NS3D(((K - KMS))) / DT));
        NS_INST(((K - KMS))) = (NS_INST(((K - KMS))) - (NS3D(((K - KMS))) / DT));
        QR3D(((K - KMS)))    = (QR3D(((K - KMS))) + QNI3D(((K - KMS))));
        NR3D(((K - KMS)))    = (NR3D(((K - KMS))) + NS3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QNI3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
        QNI3D(((K - KMS))) = 0.f;
        NS3D(((K - KMS)))  = 0.f;
      }
      if ((QG3D(((K - KMS))) < 1.E-6f)) {
        QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) + (QG3D(((K - KMS))) / DT));
        QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) - (QG3D(((K - KMS))) / DT));
        NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) + (NG3D(((K - KMS))) / DT));
        NG_INST(((K - KMS))) = (NG_INST(((K - KMS))) - (NG3D(((K - KMS))) / DT));
        QR3D(((K - KMS)))    = (QR3D(((K - KMS))) + QG3D(((K - KMS))));
        NR3D(((K - KMS)))    = (NR3D(((K - KMS))) + NG3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QG3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
        QG3D(((K - KMS))) = 0.f;
        NG3D(((K - KMS))) = 0.f;
      }
      if (((((QC3D(((K - KMS))) < morrison.QSMALL) && (QNI3D(((K - KMS))) < 1.E-8f)) &&
            (QR3D(((K - KMS))) < morrison.QSMALL)) &&
           (QG3D(((K - KMS))) < 1.E-8f))) {
        goto label_300;
      }
      if ((NI3D(((K - KMS))) < 0.f)) {
        NEGFIX_NI(((K - KMS))) = (NEGFIX_NI(((K - KMS))) + (NI3D(((K - KMS))) / DT));
      }
      if ((NS3D(((K - KMS))) < 0.f)) {
        NEGFIX_NS(((K - KMS))) = (NEGFIX_NS(((K - KMS))) + (NS3D(((K - KMS))) / DT));
      }
      if ((NC3D(((K - KMS))) < 0.f)) {
        NEGFIX_NC(((K - KMS))) = (NEGFIX_NC(((K - KMS))) + (NC3D(((K - KMS))) / DT));
      }
      if ((NR3D(((K - KMS))) < 0.f)) {
        NEGFIX_NR(((K - KMS))) = (NEGFIX_NR(((K - KMS))) + (NR3D(((K - KMS))) / DT));
      }
      if ((NG3D(((K - KMS))) < 0.f)) {
        NEGFIX_NG(((K - KMS))) = (NEGFIX_NG(((K - KMS))) + (NG3D(((K - KMS))) / DT));
      }
      NS3D(((K - KMS))) = ((0.f) > (NS3D(((K - KMS)))) ? (0.f) : (NS3D(((K - KMS)))));
      NC3D(((K - KMS))) = ((0.f) > (NC3D(((K - KMS)))) ? (0.f) : (NC3D(((K - KMS)))));
      NR3D(((K - KMS))) = ((0.f) > (NR3D(((K - KMS)))) ? (0.f) : (NR3D(((K - KMS)))));
      NG3D(((K - KMS))) = ((0.f) > (NG3D(((K - KMS)))) ? (0.f) : (NG3D(((K - KMS)))));
      if ((QR3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMR(((K - KMS))) = Kokkos::pow(
            (((morrison.PI * morrison.RHOW) * NR3D(((K - KMS)))) / QR3D(((K - KMS)))), (1.f / 3.f));
        N0RR(((K - KMS))) = (NR3D(((K - KMS))) * LAMR(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMR(((K - KMS))) < morrison.LAMMINR)) {
          TMPNUM            = NR3D(((K - KMS)));
          LAMR(((K - KMS))) = morrison.LAMMINR;
          // BFB-030: Native integer fourth power keeps both squares in default REAL.
          const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
          N0RR(((K - KMS))) =
              (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
          NR3D(((K - KMS)))       = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
          SIZEFIX_NR(((K - KMS))) = (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMR(((K - KMS))) > morrison.LAMMAXR)) {
            TMPNUM            = NR3D(((K - KMS)));
            LAMR(((K - KMS))) = morrison.LAMMAXR;
            // BFB-030: Native integer fourth power keeps both squares in default REAL.
            const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
            N0RR(((K - KMS))) =
                (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
            NR3D(((K - KMS))) = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
            SIZEFIX_NR(((K - KMS))) =
                (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QC3D(((K - KMS))) >= morrison.QSMALL)) {
        if (morrison.dofix_pgam) {
          PGAM(((K - KMS))) = morrison.pgam_fixed;
        } else {
          PGAM(((K - KMS))) =
              ((0.0005714f * ((NC3D(((K - KMS))) / 1.E6f) * RHO(((K - KMS))))) + 0.2714f);
          PGAM(((K - KMS))) = ((1.f / (PGAM(((K - KMS))) * PGAM(((K - KMS))))) - 1.f);
          PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) > (2.f) ? (PGAM(((K - KMS)))) : (2.f));
          PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) < (10.f) ? (PGAM(((K - KMS)))) : (10.f));
        }
        dumii = (static_cast<int>(PGAM(((K - KMS))))) - 1;
        nu(((K - KMS))) =
            (morrison.dnu[dumii] + ((morrison.dnu[((dumii + 1))] - morrison.dnu[dumii]) *
                                    (PGAM(((K - KMS))) - static_cast<float>(((dumii + 1))))));
        LAMC(((K - KMS))) =
            Kokkos::pow((((morrison.CONS26 * NC3D(((K - KMS)))) *
                          morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)) /
                         (QC3D(((K - KMS))) *
                          morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))),
                        (1.f / 3.f));
        LAMMIN = ((PGAM(((K - KMS))) + 1.f) / 60.E-6f);
        LAMMAX = ((PGAM(((K - KMS))) + 1.f) / 1.E-6f);
        TMPNUM = 0.f;
        if ((LAMC(((K - KMS))) < LAMMIN)) {
          LAMC(((K - KMS))) = LAMMIN;
          TMPNUM            = NC3D(((K - KMS)));
          NC3D(((K - KMS))) =
              (Kokkos::exp(
                   ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                     Kokkos::log(
                         morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))) -
                    Kokkos::log(
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)))) /
               morrison.CONS26);
          SIZEFIX_NC(((K - KMS))) = (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMC(((K - KMS))) > LAMMAX)) {
            LAMC(((K - KMS))) = LAMMAX;
            TMPNUM            = NC3D(((K - KMS)));
            NC3D(((K - KMS))) =
                (Kokkos::exp(
                     ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                       Kokkos::log(morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f),
                                                          morrison))) -
                      Kokkos::log(morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f),
                                                         morrison)))) /
                 morrison.CONS26);
            SIZEFIX_NC(((K - KMS))) =
                (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QNI3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMS(((K - KMS))) = Kokkos::pow(((morrison.CONS1 * NS3D(((K - KMS)))) / QNI3D(((K - KMS)))),
                                        (1.f / morrison.DS));
        N0S(((K - KMS)))  = (NS3D(((K - KMS))) * LAMS(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMS(((K - KMS))) < morrison.LAMMINS)) {
          LAMS(((K - KMS))) = morrison.LAMMINS;
          N0S(((K - KMS))) =
              ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
               morrison.CONS1);
          TMPNUM                  = NS3D(((K - KMS)));
          NS3D(((K - KMS)))       = (N0S(((K - KMS))) / LAMS(((K - KMS))));
          SIZEFIX_NS(((K - KMS))) = (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMS(((K - KMS))) > morrison.LAMMAXS)) {
            LAMS(((K - KMS))) = morrison.LAMMAXS;
            N0S(((K - KMS))) =
                ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
                 morrison.CONS1);
            TMPNUM            = NS3D(((K - KMS)));
            NS3D(((K - KMS))) = (N0S(((K - KMS))) / LAMS(((K - KMS))));
            SIZEFIX_NS(((K - KMS))) =
                (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QG3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMG(((K - KMS))) = Kokkos::pow(((morrison.CONS2 * NG3D(((K - KMS)))) / QG3D(((K - KMS)))),
                                        (1.f / morrison.DG));
        N0G(((K - KMS)))  = (NG3D(((K - KMS))) * LAMG(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMG(((K - KMS))) < morrison.LAMMING)) {
          LAMG(((K - KMS))) = morrison.LAMMING;
          N0G(((K - KMS))) =
              ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
               morrison.CONS2);
          TMPNUM                  = NG3D(((K - KMS)));
          NG3D(((K - KMS)))       = (N0G(((K - KMS))) / LAMG(((K - KMS))));
          SIZEFIX_NG(((K - KMS))) = (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMG(((K - KMS))) > morrison.LAMMAXG)) {
            LAMG(((K - KMS))) = morrison.LAMMAXG;
            N0G(((K - KMS))) =
                ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
                 morrison.CONS2);
            TMPNUM            = NG3D(((K - KMS)));
            NG3D(((K - KMS))) = (N0G(((K - KMS))) / LAMG(((K - KMS))));
            SIZEFIX_NG(((K - KMS))) =
                (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      PRC(((K - KMS)))    = 0.f;
      NPRC(((K - KMS)))   = 0.f;
      NPRC1(((K - KMS)))  = 0.f;
      PRA(((K - KMS)))    = 0.f;
      NPRA(((K - KMS)))   = 0.f;
      NRAGG(((K - KMS)))  = 0.f;
      PSMLT(((K - KMS)))  = 0.f;
      NSMLTS(((K - KMS))) = 0.f;
      NSMLTR(((K - KMS))) = 0.f;
      PCC(((K - KMS)))    = 0.f;
      PRE(((K - KMS)))    = 0.f;
      NSUBC(((K - KMS)))  = 0.f;
      NSUBR(((K - KMS)))  = 0.f;
      PRACG(((K - KMS)))  = 0.f;
      NPRACG(((K - KMS))) = 0.f;
      PSMLT(((K - KMS)))  = 0.f;
      EVPMS(((K - KMS)))  = 0.f;
      PGMLT(((K - KMS)))  = 0.f;
      EVPMG(((K - KMS)))  = 0.f;
      PRACS(((K - KMS)))  = 0.f;
      NPRACS(((K - KMS))) = 0.f;
      NGMLTG(((K - KMS))) = 0.f;
      NGMLTR(((K - KMS))) = 0.f;
      if ((QC3D(((K - KMS))) >= 1.E-6f)) {
        if ((morrison.IRAIN == 0)) {
          PRC(((K - KMS))) =
              ((1350.f * Kokkos::pow(QC3D(((K - KMS))), 2.47f)) *
               Kokkos::pow(((NC3D(((K - KMS))) / 1.e6f) * RHO(((K - KMS)))), (-1.79f)));
          NPRC1(((K - KMS))) = (PRC(((K - KMS))) / morrison.CONS29);
          NPRC(((K - KMS)))  = (PRC(((K - KMS))) / (QC3D(((K - KMS))) / NC3D(((K - KMS)))));
          NPRC(((K - KMS))) =
              ((NPRC(((K - KMS)))) < ((NC3D(((K - KMS))) / DT)) ? (NPRC(((K - KMS))))
                                                                : ((NC3D(((K - KMS))) / DT)));
        } else {
          if ((morrison.IRAIN == 1)) {
            DUM  = (1.f - (QC3D(((K - KMS))) / (QC3D(((K - KMS))) + QR3D(((K - KMS))))));
            DUM1 = ((600.f * Kokkos::pow(DUM, 0.68f)) *
                    (((1.f - Kokkos::pow(DUM, 0.68f)) * (1.f - Kokkos::pow(DUM, 0.68f))) *
                     (1.f - Kokkos::pow(DUM, 0.68f))));
            PRC(((K - KMS))) =
                (((((((((9.44e9f / (20.f * 2.6e-7f)) * (nu(((K - KMS))) + 2.f)) *
                       (nu(((K - KMS))) + 4.f)) /
                      ((nu(((K - KMS))) + 1.f) * (nu(((K - KMS))) + 1.f))) *
                     Kokkos::pow(((RHO(((K - KMS))) * QC3D(((K - KMS)))) / 1000.f), 4)) /
                    (((RHO(((K - KMS))) * NC3D(((K - KMS)))) / 1.e6f) *
                     ((RHO(((K - KMS))) * NC3D(((K - KMS)))) / 1.e6f))) *
                   (1.f + (DUM1 / ((1.f - DUM) * (1.f - DUM))))) *
                  1000.f) /
                 RHO(((K - KMS))));
            NPRC(((K - KMS)))  = (((PRC(((K - KMS))) * 2.f) / 2.6e-7f) * 1000.f);
            NPRC1(((K - KMS))) = (0.5f * NPRC(((K - KMS))));
          }
        }
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QNI3D(((K - KMS))) >= 1.E-8f))) {
        UMS = ((ASN(((K - KMS))) * morrison.CONS3) / Kokkos::pow(LAMS(((K - KMS))), morrison.BS));
        UMR = ((ARN(((K - KMS))) * morrison.CONS4) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        UNS = ((ASN(((K - KMS))) * morrison.CONS5) / Kokkos::pow(LAMS(((K - KMS))), morrison.BS));
        UNR = ((ARN(((K - KMS))) * morrison.CONS6) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        DUM = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
        UMS = ((UMS) < ((1.2f * DUM)) ? (UMS) : ((1.2f * DUM)));
        UNS = ((UNS) < ((1.2f * DUM)) ? (UNS) : ((1.2f * DUM)));
        UMR = ((UMR) < ((9.1f * DUM)) ? (UMR) : ((9.1f * DUM)));
        UNR = ((UNR) < ((9.1f * DUM)) ? (UNR) : ((9.1f * DUM)));
        PRACS(((K - KMS))) =
            (morrison.CONS31 *
             (((((Kokkos::pow(((((1.2f * UMR) - (0.95f * UMS)) * ((1.2f * UMR) - (0.95f * UMS))) +
                               ((0.08f * UMS) * UMR)),
                              0.5f) *
                  RHO(((K - KMS)))) *
                 N0RR(((K - KMS)))) *
                N0S(((K - KMS)))) /
               ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS))))) *
              (((5.f / (((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS)))) *
                        LAMR(((K - KMS))))) +
                (2.f / ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) *
                        (LAMR(((K - KMS))) * LAMR(((K - KMS))))))) +
               (0.5f / (LAMS(((K - KMS))) *
                        ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))))))));
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QG3D(((K - KMS))) >= 1.E-8f))) {
        UMG = ((AGN(((K - KMS))) * morrison.CONS7) / Kokkos::pow(LAMG(((K - KMS))), morrison.BG));
        UMR = ((ARN(((K - KMS))) * morrison.CONS4) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        UNG = ((AGN(((K - KMS))) * morrison.CONS8) / Kokkos::pow(LAMG(((K - KMS))), morrison.BG));
        UNR = ((ARN(((K - KMS))) * morrison.CONS6) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        DUM = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
        UMG = ((UMG) < ((20.f * DUM)) ? (UMG) : ((20.f * DUM)));
        UNG = ((UNG) < ((20.f * DUM)) ? (UNG) : ((20.f * DUM)));
        UMR = ((UMR) < ((9.1f * DUM)) ? (UMR) : ((9.1f * DUM)));
        UNR = ((UNR) < ((9.1f * DUM)) ? (UNR) : ((9.1f * DUM)));
        PRACG(((K - KMS))) =
            (morrison.CONS41 *
             (((((Kokkos::pow(((((1.2f * UMR) - (0.95f * UMG)) * ((1.2f * UMR) - (0.95f * UMG))) +
                               ((0.08f * UMG) * UMR)),
                              0.5f) *
                  RHO(((K - KMS)))) *
                 N0RR(((K - KMS)))) *
                N0G(((K - KMS)))) /
               ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) *
              (((5.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                        LAMG(((K - KMS))))) +
                (2.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                        (LAMG(((K - KMS))) * LAMG(((K - KMS))))))) +
               (0.5f / (LAMR(((K - KMS))) *
                        ((LAMG(((K - KMS))) * LAMG(((K - KMS)))) * LAMG(((K - KMS)))))))));
        DUM = (PRACG(((K - KMS))) / 5.2E-7f);
        NPRACG(((K - KMS))) =
            (((((morrison.CONS32 * RHO(((K - KMS)))) *
                Kokkos::pow(((1.7f * ((UNR - UNG) * (UNR - UNG))) + ((0.3f * UNR) * UNG)), 0.5f)) *
               N0RR(((K - KMS)))) *
              N0G(((K - KMS)))) *
             (((1.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                       LAMG(((K - KMS))))) +
               (1.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                       (LAMG(((K - KMS))) * LAMG(((K - KMS))))))) +
              (1.f / (LAMR(((K - KMS))) *
                      ((LAMG(((K - KMS))) * LAMG(((K - KMS)))) * LAMG(((K - KMS))))))));
        NPRACG(((K - KMS))) =
            (((NPRACG(((K - KMS))) - DUM)) > (0.f) ? ((NPRACG(((K - KMS))) - DUM)) : (0.f));
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QC3D(((K - KMS))) >= 1.E-8f))) {
        if ((morrison.IRAIN == 0)) {
          DUM               = (QC3D(((K - KMS))) * QR3D(((K - KMS))));
          PRA(((K - KMS)))  = (67.f * Kokkos::pow(DUM, 1.15f));
          NPRA(((K - KMS))) = (PRA(((K - KMS))) / (QC3D(((K - KMS))) / NC3D(((K - KMS)))));
        } else {
          if ((morrison.IRAIN == 1)) {
            DUM  = (1.f - (QC3D(((K - KMS))) / (QC3D(((K - KMS))) + QR3D(((K - KMS))))));
            DUM1 = Kokkos::pow((DUM / (DUM + 5.e-4f)), 4);
            PRA(((K - KMS)))  = (((((5.78e3f * RHO(((K - KMS)))) / 1000.f) * QC3D(((K - KMS)))) *
                                  QR3D(((K - KMS)))) *
                                 DUM1);
            NPRA(((K - KMS))) = ((((((PRA(((K - KMS))) * RHO(((K - KMS)))) / 1000.f) *
                                    ((NC3D(((K - KMS))) * RHO(((K - KMS)))) / 1.e6f)) /
                                   ((QC3D(((K - KMS))) * RHO(((K - KMS)))) / 1000.f)) *
                                  1.e6f) /
                                 RHO(((K - KMS))));
          }
        }
      }
      if ((QR3D(((K - KMS))) >= 1.E-8f)) {
        DUM1 = 300.e-6f;
        if (((1.f / LAMR(((K - KMS)))) < DUM1)) {
          DUM = 1.f;
        } else {
          if (((1.f / LAMR(((K - KMS)))) >= DUM1)) {
            DUM = (2.f - Kokkos::exp((2300.f * ((1.f / LAMR(((K - KMS)))) - DUM1))));
          }
        }
        NRAGG(((K - KMS))) =
            (((((-5.78f) * DUM) * NR3D(((K - KMS)))) * QR3D(((K - KMS)))) * RHO(((K - KMS))));
      }
      if ((QR3D(((K - KMS))) >= morrison.QSMALL)) {
        EPSR = (((((2.f * morrison.PI) * N0RR(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                ((morrison.F1R / (LAMR(((K - KMS))) * LAMR(((K - KMS))))) +
                 ((((morrison.F2R *
                     Kokkos::pow(((ARN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                    Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                   morrison.CONS9) /
                  Kokkos::pow(LAMR(((K - KMS))), morrison.CONS34))));
      } else {
        EPSR = 0.f;
      }
      if ((QV3D(((K - KMS))) < QVS(((K - KMS))))) {
        PRE(((K - KMS))) = ((EPSR * (QV3D(((K - KMS))) - QVS(((K - KMS))))) / AB(((K - KMS))));
        PRE(((K - KMS))) = ((PRE(((K - KMS)))) < (0.f) ? (PRE(((K - KMS)))) : (0.f));
      } else {
        PRE(((K - KMS))) = 0.f;
      }
      if ((QNI3D(((K - KMS))) >= 1.E-8f)) {
        DUM =
            ((((-morrison.CPW) / XLF(((K - KMS)))) * (((T3D(((K - KMS))) - morrison.TMELT)) > (0.0f)
                                                          ? ((T3D(((K - KMS))) - morrison.TMELT))
                                                          : (0.0f))) *
             PRACS(((K - KMS))));
        PSMLT(((K - KMS))) =
            ((((((((2.f * morrison.PI) * N0S(((K - KMS)))) * KAP(((K - KMS)))) *
                 (((morrison.TMELT - T3D(((K - KMS))))) < (0.0f)
                      ? ((morrison.TMELT - T3D(((K - KMS)))))
                      : (0.0f))) /
                XLF(((K - KMS)))) *
               RHO(((K - KMS)))) *
              ((morrison.F1S / (LAMS(((K - KMS))) * LAMS(((K - KMS))))) +
               ((((morrison.F2S *
                   Kokkos::pow(((ASN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                  Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                 morrison.CONS10) /
                Kokkos::pow(LAMS(((K - KMS))), morrison.CONS35)))) +
             DUM);
        if ((QVQVS(((K - KMS))) < 1.f)) {
          EPSS =
              (((((2.f * morrison.PI) * N0S(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
               ((morrison.F1S / (LAMS(((K - KMS))) * LAMS(((K - KMS))))) +
                ((((morrison.F2S *
                    Kokkos::pow(((ASN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                   Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                  morrison.CONS10) /
                 Kokkos::pow(LAMS(((K - KMS))), morrison.CONS35))));
          EVPMS(((K - KMS))) = (((QV3D(((K - KMS))) - QVS(((K - KMS)))) * EPSS) / AB(((K - KMS))));
          EVPMS(((K - KMS))) = ((EVPMS(((K - KMS)))) > (PSMLT(((K - KMS)))) ? (EVPMS(((K - KMS))))
                                                                            : (PSMLT(((K - KMS)))));
          PSMLT(((K - KMS))) = (PSMLT(((K - KMS))) - EVPMS(((K - KMS))));
        }
      }
      if ((QG3D(((K - KMS))) >= 1.E-8f)) {
        DUM =
            ((((-morrison.CPW) / XLF(((K - KMS)))) * (((T3D(((K - KMS))) - morrison.TMELT)) > (0.0f)
                                                          ? ((T3D(((K - KMS))) - morrison.TMELT))
                                                          : (0.0f))) *
             PRACG(((K - KMS))));
        PGMLT(((K - KMS))) =
            ((((((((2.f * morrison.PI) * N0G(((K - KMS)))) * KAP(((K - KMS)))) *
                 (morrison.TMELT - T3D(((K - KMS))))) /
                XLF(((K - KMS)))) *
               RHO(((K - KMS)))) *
              ((morrison.F1S / (LAMG(((K - KMS))) * LAMG(((K - KMS))))) +
               ((((morrison.F2S *
                   Kokkos::pow(((AGN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                  Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                 morrison.CONS11) /
                Kokkos::pow(LAMG(((K - KMS))), morrison.CONS36)))) +
             DUM);
        if ((QVQVS(((K - KMS))) < 1.f)) {
          EPSG =
              (((((2.f * morrison.PI) * N0G(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
               ((morrison.F1S / (LAMG(((K - KMS))) * LAMG(((K - KMS))))) +
                ((((morrison.F2S *
                    Kokkos::pow(((AGN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                   Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                  morrison.CONS11) /
                 Kokkos::pow(LAMG(((K - KMS))), morrison.CONS36))));
          EVPMG(((K - KMS))) = (((QV3D(((K - KMS))) - QVS(((K - KMS)))) * EPSG) / AB(((K - KMS))));
          EVPMG(((K - KMS))) = ((EVPMG(((K - KMS)))) > (PGMLT(((K - KMS)))) ? (EVPMG(((K - KMS))))
                                                                            : (PGMLT(((K - KMS)))));
          PGMLT(((K - KMS))) = (PGMLT(((K - KMS))) - EVPMG(((K - KMS))));
        }
      }
      PRACG(((K - KMS))) = 0.f;
      PRACS(((K - KMS))) = 0.f;
      DUM                = ((PRC(((K - KMS))) + PRA(((K - KMS)))) * DT);
      if (((DUM > QC3D(((K - KMS)))) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO            = (QC3D(((K - KMS))) / DUM);
        PRC(((K - KMS))) = (PRC(((K - KMS))) * RATIO);
        PRA(((K - KMS))) = (PRA(((K - KMS))) * RATIO);
      }
      DUM = ((((-PSMLT(((K - KMS)))) - EVPMS(((K - KMS)))) + PRACS(((K - KMS)))) * DT);
      if (((DUM > QNI3D(((K - KMS)))) && (QNI3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO              = (QNI3D(((K - KMS))) / DUM);
        PSMLT(((K - KMS))) = (PSMLT(((K - KMS))) * RATIO);
        EVPMS(((K - KMS))) = (EVPMS(((K - KMS))) * RATIO);
        PRACS(((K - KMS))) = (PRACS(((K - KMS))) * RATIO);
      }
      DUM = ((((-PGMLT(((K - KMS)))) - EVPMG(((K - KMS)))) + PRACG(((K - KMS)))) * DT);
      if (((DUM > QG3D(((K - KMS)))) && (QG3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO              = (QG3D(((K - KMS))) / DUM);
        PGMLT(((K - KMS))) = (PGMLT(((K - KMS))) * RATIO);
        EVPMG(((K - KMS))) = (EVPMG(((K - KMS))) * RATIO);
        PRACG(((K - KMS))) = (PRACG(((K - KMS))) * RATIO);
      }
      DUM = ((((((((-PRACS(((K - KMS)))) - PRACG(((K - KMS)))) - PRE(((K - KMS)))) -
                 PRA(((K - KMS)))) -
                PRC(((K - KMS)))) +
               PSMLT(((K - KMS)))) +
              PGMLT(((K - KMS)))) *
             DT);
      if (((DUM > QR3D(((K - KMS)))) && (QR3D(((K - KMS))) >= morrison.QSMALL))) {
        // BFB-033: Match Native's balanced source sum before rain-evaporation limiting.
        RATIO            = ((((PRC(((K - KMS))) + (PRACG(((K - KMS))) + PRA(((K - KMS))))) +
                              ((QR3D(((K - KMS))) / DT) + PRACS(((K - KMS))))) -
                             PSMLT(((K - KMS)))) -
                            PGMLT(((K - KMS)))) /
                           (-PRE(((K - KMS))));
        PRE(((K - KMS))) = (PRE(((K - KMS))) * RATIO);
      }
      if (((PRE(((K - KMS))) < 0.0f) && (CF3D(((K - KMS))) > morrison.cloud_frac_thresh))) {
        std::cerr << std::string("Rain evaporation occurring within cloud") << "\n";
      }
      if (((PRACS(((K - KMS))) != 0.f) || (PRACG(((K - KMS))) != 0.f))) {
        std::cerr << std::string("PRACS or PRACG is violating expectation of zero value") << "\n";
        throw std::runtime_error(std::string("Error in Morrison microphysics."));
      }
      QV3DTEN(((K - KMS))) = (QV3DTEN(((K - KMS))) +
                              (((-PRE(((K - KMS)))) - EVPMS(((K - KMS)))) - EVPMG(((K - KMS)))));
      T3DTEN(((K - KMS)))  = (T3DTEN(((K - KMS))) +
                              ((((PRE(((K - KMS))) * XXLV(((K - KMS)))) +
                                 ((EVPMS(((K - KMS))) + EVPMG(((K - KMS)))) * XXLS(((K - KMS))))) +
                                ((((PSMLT(((K - KMS))) + PGMLT(((K - KMS)))) - PRACS(((K - KMS)))) -
                                  PRACG(((K - KMS)))) *
                                 XLF(((K - KMS))))) /
                               CPM(((K - KMS)))));
      QC3DTEN(((K - KMS))) = (QC3DTEN(((K - KMS))) + ((-PRA(((K - KMS)))) - PRC(((K - KMS)))));
      QR3DTEN(((K - KMS))) =
          (QR3DTEN(((K - KMS))) +
           ((((((PRE(((K - KMS))) + PRA(((K - KMS)))) + PRC(((K - KMS)))) - PSMLT(((K - KMS)))) -
              PGMLT(((K - KMS)))) +
             PRACS(((K - KMS)))) +
            PRACG(((K - KMS)))));
      QNI3DTEN(((K - KMS))) = (QNI3DTEN(((K - KMS))) +
                               ((PSMLT(((K - KMS))) + EVPMS(((K - KMS)))) - PRACS(((K - KMS)))));
      QG3DTEN(((K - KMS))) =
          (QG3DTEN(((K - KMS))) + ((PGMLT(((K - KMS))) + EVPMG(((K - KMS)))) - PRACG(((K - KMS)))));
      NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + ((-NPRA(((K - KMS)))) - NPRC(((K - KMS)))));
      NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) +
                              ((NPRC1(((K - KMS))) + NRAGG(((K - KMS)))) - NPRACG(((K - KMS)))));
      if ((PRE(((K - KMS))) < 0.f)) {
        DUM                = ((PRE(((K - KMS))) * DT) / QR3D(((K - KMS))));
        DUM                = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSUBR(((K - KMS))) = ((DUM * NR3D(((K - KMS)))) / DT);
      }
      if (((EVPMS(((K - KMS))) + PSMLT(((K - KMS)))) < 0.f)) {
        DUM = (((EVPMS(((K - KMS))) + PSMLT(((K - KMS)))) * DT) / QNI3D(((K - KMS))));
        DUM = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSMLTS(((K - KMS))) = ((DUM * NS3D(((K - KMS)))) / DT);
      }
      if ((PSMLT(((K - KMS))) < 0.f)) {
        DUM                 = ((PSMLT(((K - KMS))) * DT) / QNI3D(((K - KMS))));
        DUM                 = (((-1.0f)) > (DUM) ? ((-1.0f)) : (DUM));
        NSMLTR(((K - KMS))) = ((DUM * NS3D(((K - KMS)))) / DT);
      }
      if (((EVPMG(((K - KMS))) + PGMLT(((K - KMS)))) < 0.f)) {
        DUM = (((EVPMG(((K - KMS))) + PGMLT(((K - KMS)))) * DT) / QG3D(((K - KMS))));
        DUM = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NGMLTG(((K - KMS))) = ((DUM * NG3D(((K - KMS)))) / DT);
      }
      if ((PGMLT(((K - KMS))) < 0.f)) {
        DUM                 = ((PGMLT(((K - KMS))) * DT) / QG3D(((K - KMS))));
        DUM                 = (((-1.0f)) > (DUM) ? ((-1.0f)) : (DUM));
        NGMLTR(((K - KMS))) = ((DUM * NG3D(((K - KMS)))) / DT);
      }
      NS3DTEN(((K - KMS))) = (NS3DTEN(((K - KMS))) + NSMLTS(((K - KMS))));
      NG3DTEN(((K - KMS))) = (NG3DTEN(((K - KMS))) + NGMLTG(((K - KMS))));
      NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) +
                              ((NSUBR(((K - KMS))) - NSMLTR(((K - KMS)))) - NGMLTR(((K - KMS)))));
    label_300:;
      if ((morrison.ISATADJ == 0)) {
        DUMT   = (T3D(((K - KMS))) + (DT * T3DTEN(((K - KMS)))));
        DUMQV  = (QV3D(((K - KMS))) + (DT * QV3DTEN(((K - KMS)))));
        DUM    = (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(DUMT, 0, morrison))
                      ? ((0.99f * PRES(((K - KMS)))))
                      : (morrison_engine::POLYSVP(DUMT, 0, morrison)));
        DUMQSS = ((ep * DUM) / (PRES(((K - KMS))) - DUM));
        DUMQC  = (QC3D(((K - KMS))) + (DT * QC3DTEN(((K - KMS)))));
        DUMQC  = ((DUMQC) > (0.f) ? (DUMQC) : (0.f));
        DUMS   = (DUMQV - DUMQSS);
        PCC(((K - KMS))) = ((DUMS / (1.f + (((XXLV(((K - KMS))) * XXLV(((K - KMS)))) * DUMQSS) /
                                            ((CPM(((K - KMS))) * Rv) * (DUMT * DUMT))))) /
                            DT);
        if (((((PCC(((K - KMS))) * DT) + QC3D(((K - KMS)))) + (DT * QC3DTEN(((K - KMS))))) < 0.f)) {
          PCC(((K - KMS))) = ((-(QC3D(((K - KMS))) + (DT * QC3DTEN(((K - KMS)))))) / DT);
        }
        QV3DTEN(((K - KMS))) = (QV3DTEN(((K - KMS))) - PCC(((K - KMS))));
        T3DTEN(((K - KMS))) =
            (T3DTEN(((K - KMS))) + ((PCC(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QC3DTEN(((K - KMS))) = (QC3DTEN(((K - KMS))) + PCC(((K - KMS))));
      }
      if ((((QC3D(((K - KMS))) + (QC3DTEN(((K - KMS))) * DT)) >= morrison.QSMALL) &&
           (morrison.INUM == 0))) {
        if ((morrison.ISUB == 0)) {
          DUM = (W3D(((K - KMS))) + WVAR(((K - KMS))));
          DUM = ((DUM) > (0.10f) ? (DUM) : (0.10f));
        } else {
          if ((morrison.ISUB == 1)) {
            DUM = W3D(((K - KMS)));
          }
        }
        if ((DUM >= 0.001f)) {
          if ((morrison.IBASE == 1)) {
            IDROP = 0;
            if ((QC3D(((K - KMS))) <= (0.05E-3f / RHO(((K - KMS)))))) {
              IDROP = 1;
            }
            if ((K == 0)) {
              IDROP = 1;
            } else {
              if ((K >= 1)) {
                if (((QC3D(((K - KMS))) > (0.05E-3f / RHO(((K - KMS))))) &&
                     (QC3D(((K - KMS - 1))) <= (0.05E-3f / RHO(((K - KMS - 1))))))) {
                  IDROP = 1;
                }
              }
            }
            if ((IDROP == 1)) {
              if ((morrison.IACT == 1)) {
                DUM  = (DUM * 100.f);
                DUM2 = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                        Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                    (morrison.K1 / (morrison.K1 + 2.f))));
                DUM2 = (DUM2 * 1.E6f);
                DUM2                 = (DUM2 / RHO(((K - KMS))));
                DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
              } else {
                if ((morrison.IACT == 2)) {
                  SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                  AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                           T3D(((K - KMS))));
                  ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                            ((CPM(((K - KMS))) * morrison.RR) *
                             (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                           ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                  GAMM  = (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                  GG  = (1.f /
                         ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                           ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                          (((XXLV(((K - KMS))) * morrison.RHOW) /
                            (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                           (((XXLV(((K - KMS))) * morrison.MW) / (T3D(((K - KMS))) * morrison.RR)) -
                            1.f))));
                  PSI = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                  ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                  ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                  SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                  SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                  DUM1 =
                      ((1.f / (SM1 * SM1)) *
                       ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                        (morrison.F21 * Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                  DUM2 =
                      ((1.f / (SM2 * SM2)) *
                       ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                        (morrison.F22 * Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                  SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                  UU1 = ((2.f * Kokkos::log((SM1 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG1)));
                  UU2 = ((2.f * Kokkos::log((SM2 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG2)));
                  DUM1 = ((morrison.NANEW1 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                  DUM2 = ((morrison.NANEW2 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                  DUM2 = ((DUM1 + DUM2) / RHO(((K - KMS))));
                  DUM2 = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                              ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                              : (DUM2));
                  DUM2 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                  DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                  NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                  NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                }
              }
            } else {
              if ((IDROP == 0)) {
                TAUC = (1.f / ((((((2.f * morrison.PI) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                                 NC3D(((K - KMS)))) *
                                (PGAM(((K - KMS))) + 1.f)) /
                               LAMC(((K - KMS)))));
                if ((EPSR > 1.E-8f)) {
                  TAUR = (1.f / EPSR);
                } else {
                  TAUR = 1.E8f;
                }
                DUM3 = (((((QVS(((K - KMS))) * RHO(((K - KMS)))) /
                           (PRES(((K - KMS))) - EVS(((K - KMS))))) +
                          (DQSDT / Cp)) *
                         morrison.G) *
                        DUM);
                DUM3 = (((DUM3 * TAUC) * TAUR) / (TAUC + TAUR));
                if (((DUM3 / QVS(((K - KMS)))) >= 1.E-6f)) {
                  if ((morrison.IACT == 1)) {
                    DUM    = (DUM * 100.f);
                    DUMACT = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                              Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                          (morrison.K1 / (morrison.K1 + 2.f))));
                    DUM3                 = ((DUM3 / QVS(((K - KMS)))) * 100.f);
                    DUM2   = (morrison.C1 * Kokkos::pow(DUM3, morrison.K1));
                    DUM2   = ((DUM2) < (DUMACT) ? (DUM2) : (DUMACT));
                    DUM2   = (DUM2 * 1.E6f);
                    DUM2                 = (DUM2 / RHO(((K - KMS))));
                    DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                    DUM2   = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                    NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                    NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                  } else {
                    if ((morrison.IACT == 2)) {
                      SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                      AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                               T3D(((K - KMS))));
                      ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                                ((CPM(((K - KMS))) * morrison.RR) *
                                 (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                               ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                      GAMM =
                          (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                      GG     = (1.f / ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                                        ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                                       (((XXLV(((K - KMS))) * morrison.RHOW) /
                                         (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                                        (((XXLV(((K - KMS))) * morrison.MW) /
                                          (T3D(((K - KMS))) * morrison.RR)) -
                                         1.f))));
                      PSI  = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                      ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                              ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                      ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                              ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                      SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                              Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                      SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                              Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                      DUM1 = ((1.f / (SM1 * SM1)) *
                              ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                               (morrison.F21 *
                                Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                      DUM2 = ((1.f / (SM2 * SM2)) *
                              ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                               (morrison.F22 *
                                Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                      SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                      UU1  = ((2.f * Kokkos::log((SM1 / SMAX))) /
                              (4.242f * Kokkos::log(morrison.SIG1)));
                      UU2  = ((2.f * Kokkos::log((SM2 / SMAX))) /
                              (4.242f * Kokkos::log(morrison.SIG2)));
                      DUM1   = ((morrison.NANEW1 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                      DUM2   = ((morrison.NANEW2 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                      DUM2   = ((DUM1 + DUM2) / RHO(((K - KMS))));
                      DUMACT = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                                    ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                                    : (DUM2));
                      SIGVL  = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                      AACT   = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                                T3D(((K - KMS))));
                      SM1    = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                                Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                      SM2    = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                                Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                      SMAX   = (DUM3 / QVS(((K - KMS))));
                      UU1    = ((2.f * Kokkos::log((SM1 / SMAX))) /
                                (4.242f * Kokkos::log(morrison.SIG1)));
                      UU2    = ((2.f * Kokkos::log((SM2 / SMAX))) /
                                (4.242f * Kokkos::log(morrison.SIG2)));
                      DUM1   = ((morrison.NANEW1 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                      DUM2   = ((morrison.NANEW2 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                      DUM2   = ((DUM1 + DUM2) / RHO(((K - KMS))));
                      DUM2   = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                                    ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                                    : (DUM2));
                      DUM2 = ((DUM2) < (DUMACT) ? (DUM2) : (DUMACT));
                      DUM2   = ((DUM2 - NC3D(((K - KMS)))) / DT);
                      DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                      NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                      NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                    }
                  }
                }
              }
            }
          } else {
            if ((morrison.IBASE == 2)) {
              if ((morrison.IACT == 1)) {
                DUM  = (DUM * 100.f);
                DUM2 = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                        Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                    (morrison.K1 / (morrison.K1 + 2.f))));
                DUM2 = (DUM2 * 1.E6f);
                DUM2                 = (DUM2 / RHO(((K - KMS))));
                DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
              } else {
                if ((morrison.IACT == 2)) {
                  SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                  AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                           T3D(((K - KMS))));
                  ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                            ((CPM(((K - KMS))) * morrison.RR) *
                             (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                           ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                  GAMM  = (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                  GG  = (1.f /
                         ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                           ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                          (((XXLV(((K - KMS))) * morrison.RHOW) /
                            (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                           (((XXLV(((K - KMS))) * morrison.MW) / (T3D(((K - KMS))) * morrison.RR)) -
                            1.f))));
                  PSI = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                  ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                  ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                  SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                  SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                  DUM1 =
                      ((1.f / (SM1 * SM1)) *
                       ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                        (morrison.F21 * Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                  DUM2 =
                      ((1.f / (SM2 * SM2)) *
                       ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                        (morrison.F22 * Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                  SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                  UU1 = ((2.f * Kokkos::log((SM1 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG1)));
                  UU2 = ((2.f * Kokkos::log((SM2 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG2)));
                  DUM1 = ((morrison.NANEW1 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                  DUM2 = ((morrison.NANEW2 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                  DUM2 = ((DUM1 + DUM2) / RHO(((K - KMS))));
                  DUM2 = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                              ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                              : (DUM2));
                  DUM2 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                  DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                  NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                  NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                }
              }
            }
          }
        }
      }
    } else {
      if ((morrison.INUM == 1)) {
        NC3D(((K - KMS))) = ((morrison.NDCNST * 1.E6f) / RHO(((K - KMS))));
      }
      if ((NI3D(((K - KMS))) < 0.f)) {
        NEGFIX_NI(((K - KMS))) = (NEGFIX_NI(((K - KMS))) + (NI3D(((K - KMS))) / DT));
      }
      if ((NS3D(((K - KMS))) < 0.f)) {
        NEGFIX_NS(((K - KMS))) = (NEGFIX_NS(((K - KMS))) + (NS3D(((K - KMS))) / DT));
      }
      if ((NC3D(((K - KMS))) < 0.f)) {
        NEGFIX_NC(((K - KMS))) = (NEGFIX_NC(((K - KMS))) + (NC3D(((K - KMS))) / DT));
      }
      if ((NR3D(((K - KMS))) < 0.f)) {
        NEGFIX_NR(((K - KMS))) = (NEGFIX_NR(((K - KMS))) + (NR3D(((K - KMS))) / DT));
      }
      if ((NG3D(((K - KMS))) < 0.f)) {
        NEGFIX_NG(((K - KMS))) = (NEGFIX_NG(((K - KMS))) + (NG3D(((K - KMS))) / DT));
      }
      NI3D(((K - KMS))) = ((0.f) > (NI3D(((K - KMS)))) ? (0.f) : (NI3D(((K - KMS)))));
      NS3D(((K - KMS))) = ((0.f) > (NS3D(((K - KMS)))) ? (0.f) : (NS3D(((K - KMS)))));
      NC3D(((K - KMS))) = ((0.f) > (NC3D(((K - KMS)))) ? (0.f) : (NC3D(((K - KMS)))));
      NR3D(((K - KMS))) = ((0.f) > (NR3D(((K - KMS)))) ? (0.f) : (NR3D(((K - KMS)))));
      NG3D(((K - KMS))) = ((0.f) > (NG3D(((K - KMS)))) ? (0.f) : (NG3D(((K - KMS)))));
      if ((QI3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMI(((K - KMS))) = Kokkos::pow(((morrison.CONS12 * NI3D(((K - KMS)))) / QI3D(((K - KMS)))),
                                        (1.f / morrison.DI));
        N0I(((K - KMS)))  = (NI3D(((K - KMS))) * LAMI(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMI(((K - KMS))) < morrison.LAMMINI)) {
          LAMI(((K - KMS))) = morrison.LAMMINI;
          TMPNUM            = NI3D(((K - KMS)));
          N0I(((K - KMS))) =
              ((Kokkos::pow(LAMI(((K - KMS))), (morrison.DI + 1.f)) * QI3D(((K - KMS)))) /
               morrison.CONS12);
          NI3D(((K - KMS)))       = (N0I(((K - KMS))) / LAMI(((K - KMS))));
          SIZEFIX_NI(((K - KMS))) = (SIZEFIX_NI(((K - KMS))) + ((NI3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMI(((K - KMS))) > morrison.LAMMAXI)) {
            LAMI(((K - KMS))) = morrison.LAMMAXI;
            N0I(((K - KMS))) =
                ((Kokkos::pow(LAMI(((K - KMS))), (morrison.DI + 1.f)) * QI3D(((K - KMS)))) /
                 morrison.CONS12);
            TMPNUM            = NI3D(((K - KMS)));
            NI3D(((K - KMS))) = (N0I(((K - KMS))) / LAMI(((K - KMS))));
            SIZEFIX_NI(((K - KMS))) =
                (SIZEFIX_NI(((K - KMS))) + ((NI3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QR3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMR(((K - KMS))) = Kokkos::pow(
            (((morrison.PI * morrison.RHOW) * NR3D(((K - KMS)))) / QR3D(((K - KMS)))), (1.f / 3.f));
        N0RR(((K - KMS))) = (NR3D(((K - KMS))) * LAMR(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMR(((K - KMS))) < morrison.LAMMINR)) {
          TMPNUM            = NR3D(((K - KMS)));
          LAMR(((K - KMS))) = morrison.LAMMINR;
          // BFB-030: Native integer fourth power keeps both squares in default REAL.
          const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
          N0RR(((K - KMS))) =
              (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
          NR3D(((K - KMS)))       = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
          SIZEFIX_NR(((K - KMS))) = (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMR(((K - KMS))) > morrison.LAMMAXR)) {
            TMPNUM            = NR3D(((K - KMS)));
            LAMR(((K - KMS))) = morrison.LAMMAXR;
            // BFB-030: Native integer fourth power keeps both squares in default REAL.
            const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
            N0RR(((K - KMS))) =
                (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
            NR3D(((K - KMS))) = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
            SIZEFIX_NR(((K - KMS))) =
                (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QC3D(((K - KMS))) >= morrison.QSMALL)) {
        if (morrison.dofix_pgam) {
          PGAM(((K - KMS))) = morrison.pgam_fixed;
        } else {
          PGAM(((K - KMS))) =
              ((0.0005714f * ((NC3D(((K - KMS))) / 1.E6f) * RHO(((K - KMS))))) + 0.2714f);
          PGAM(((K - KMS))) = ((1.f / (PGAM(((K - KMS))) * PGAM(((K - KMS))))) - 1.f);
          PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) > (2.f) ? (PGAM(((K - KMS)))) : (2.f));
          PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) < (10.f) ? (PGAM(((K - KMS)))) : (10.f));
        }
        dumii = (static_cast<int>(PGAM(((K - KMS))))) - 1;
        nu(((K - KMS))) =
            (morrison.dnu[dumii] + ((morrison.dnu[((dumii + 1))] - morrison.dnu[dumii]) *
                                    (PGAM(((K - KMS))) - static_cast<float>(((dumii + 1))))));
        LAMC(((K - KMS))) =
            Kokkos::pow((((morrison.CONS26 * NC3D(((K - KMS)))) *
                          morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)) /
                         (QC3D(((K - KMS))) *
                          morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))),
                        (1.f / 3.f));
        LAMMIN = ((PGAM(((K - KMS))) + 1.f) / 60.E-6f);
        LAMMAX = ((PGAM(((K - KMS))) + 1.f) / 1.E-6f);
        TMPNUM = 0.f;
        if ((LAMC(((K - KMS))) < LAMMIN)) {
          LAMC(((K - KMS))) = LAMMIN;
          TMPNUM            = NC3D(((K - KMS)));
          NC3D(((K - KMS))) =
              (Kokkos::exp(
                   ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                     Kokkos::log(
                         morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))) -
                    Kokkos::log(
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)))) /
               morrison.CONS26);
          SIZEFIX_NC(((K - KMS))) = (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMC(((K - KMS))) > LAMMAX)) {
            LAMC(((K - KMS))) = LAMMAX;
            TMPNUM            = NC3D(((K - KMS)));
            NC3D(((K - KMS))) =
                (Kokkos::exp(
                     ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                       Kokkos::log(morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f),
                                                          morrison))) -
                      Kokkos::log(morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f),
                                                         morrison)))) /
                 morrison.CONS26);
            SIZEFIX_NC(((K - KMS))) =
                (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
        CDIST1(((K - KMS))) =
            (NC3D(((K - KMS))) /
             morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison));
      }
      if ((QNI3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMS(((K - KMS))) = Kokkos::pow(((morrison.CONS1 * NS3D(((K - KMS)))) / QNI3D(((K - KMS)))),
                                        (1.f / morrison.DS));
        N0S(((K - KMS)))  = (NS3D(((K - KMS))) * LAMS(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMS(((K - KMS))) < morrison.LAMMINS)) {
          LAMS(((K - KMS))) = morrison.LAMMINS;
          N0S(((K - KMS))) =
              ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
               morrison.CONS1);
          TMPNUM                  = NS3D(((K - KMS)));
          NS3D(((K - KMS)))       = (N0S(((K - KMS))) / LAMS(((K - KMS))));
          SIZEFIX_NS(((K - KMS))) = (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMS(((K - KMS))) > morrison.LAMMAXS)) {
            LAMS(((K - KMS))) = morrison.LAMMAXS;
            N0S(((K - KMS))) =
                ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
                 morrison.CONS1);
            TMPNUM            = NS3D(((K - KMS)));
            NS3D(((K - KMS))) = (N0S(((K - KMS))) / LAMS(((K - KMS))));
            SIZEFIX_NS(((K - KMS))) =
                (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      if ((QG3D(((K - KMS))) >= morrison.QSMALL)) {
        LAMG(((K - KMS))) = Kokkos::pow(((morrison.CONS2 * NG3D(((K - KMS)))) / QG3D(((K - KMS)))),
                                        (1.f / morrison.DG));
        N0G(((K - KMS)))  = (NG3D(((K - KMS))) * LAMG(((K - KMS))));
        TMPNUM            = 0.f;
        if ((LAMG(((K - KMS))) < morrison.LAMMING)) {
          LAMG(((K - KMS))) = morrison.LAMMING;
          N0G(((K - KMS))) =
              ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
               morrison.CONS2);
          TMPNUM                  = NG3D(((K - KMS)));
          NG3D(((K - KMS)))       = (N0G(((K - KMS))) / LAMG(((K - KMS))));
          SIZEFIX_NG(((K - KMS))) = (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
        } else {
          if ((LAMG(((K - KMS))) > morrison.LAMMAXG)) {
            LAMG(((K - KMS))) = morrison.LAMMAXG;
            N0G(((K - KMS))) =
                ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
                 morrison.CONS2);
            TMPNUM            = NG3D(((K - KMS)));
            NG3D(((K - KMS))) = (N0G(((K - KMS))) / LAMG(((K - KMS))));
            SIZEFIX_NG(((K - KMS))) =
                (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
          }
        }
      }
      MNUCCC(((K - KMS)))  = 0.f;
      NNUCCC(((K - KMS)))  = 0.f;
      PRC(((K - KMS)))     = 0.f;
      NPRC(((K - KMS)))    = 0.f;
      NPRC1(((K - KMS)))   = 0.f;
      NSAGG(((K - KMS)))   = 0.f;
      PSACWS(((K - KMS)))  = 0.f;
      NPSACWS(((K - KMS))) = 0.f;
      PSACWI(((K - KMS)))  = 0.f;
      NPSACWI(((K - KMS))) = 0.f;
      PRACS(((K - KMS)))   = 0.f;
      NPRACS(((K - KMS)))  = 0.f;
      NMULTS(((K - KMS)))  = 0.f;
      QMULTS(((K - KMS)))  = 0.f;
      NMULTR(((K - KMS)))  = 0.f;
      QMULTR(((K - KMS)))  = 0.f;
      NMULTG(((K - KMS)))  = 0.f;
      QMULTG(((K - KMS)))  = 0.f;
      NMULTRG(((K - KMS))) = 0.f;
      QMULTRG(((K - KMS))) = 0.f;
      MNUCCR(((K - KMS)))  = 0.f;
      NNUCCR(((K - KMS)))  = 0.f;
      PRA(((K - KMS)))     = 0.f;
      NPRA(((K - KMS)))    = 0.f;
      NRAGG(((K - KMS)))   = 0.f;
      PRCI(((K - KMS)))    = 0.f;
      NPRCI(((K - KMS)))   = 0.f;
      PRAI(((K - KMS)))    = 0.f;
      NPRAI(((K - KMS)))   = 0.f;
      NNUCCD(((K - KMS)))  = 0.f;
      MNUCCD(((K - KMS)))  = 0.f;
      PCC(((K - KMS)))     = 0.f;
      PRE(((K - KMS)))     = 0.f;
      PRD(((K - KMS)))     = 0.f;
      PRDS(((K - KMS)))    = 0.f;
      EPRD(((K - KMS)))    = 0.f;
      EPRDS(((K - KMS)))   = 0.f;
      NSUBC(((K - KMS)))   = 0.f;
      NSUBI(((K - KMS)))   = 0.f;
      NSUBS(((K - KMS)))   = 0.f;
      NSUBR(((K - KMS)))   = 0.f;
      PIACR(((K - KMS)))   = 0.f;
      NIACR(((K - KMS)))   = 0.f;
      PRACI(((K - KMS)))   = 0.f;
      PIACRS(((K - KMS)))  = 0.f;
      NIACRS(((K - KMS)))  = 0.f;
      PRACIS(((K - KMS)))  = 0.f;
      PRACG(((K - KMS)))   = 0.f;
      PSACR(((K - KMS)))   = 0.f;
      PSACWG(((K - KMS)))  = 0.f;
      PGSACW(((K - KMS)))  = 0.f;
      PGRACS(((K - KMS)))  = 0.f;
      PRDG(((K - KMS)))    = 0.f;
      EPRDG(((K - KMS)))   = 0.f;
      NPRACG(((K - KMS)))  = 0.f;
      NPSACWG(((K - KMS))) = 0.f;
      NSCNG(((K - KMS)))   = 0.f;
      NGRACS(((K - KMS)))  = 0.f;
      NSUBG(((K - KMS)))   = 0.f;
      if (((QC3D(((K - KMS))) >= morrison.QSMALL) && (T3D(((K - KMS))) < 269.15f))) {
        NACNT = (Kokkos::exp(((-2.80f) + (0.262f * (morrison.TMELT - T3D(((K - KMS))))))) * 1000.f);
        DUM   = (((7.37f * T3D(((K - KMS)))) / ((288.f * 10.f) * PRES(((K - KMS))))) / 100.f);
        DAP(((K - KMS))) = (((morrison.CONS37 * T3D(((K - KMS)))) * (1.f + (DUM / morrison.RIN))) /
                            MU(((K - KMS))));
        MNUCCC(((K - KMS))) = (((morrison.CONS38 * DAP(((K - KMS)))) * NACNT) *
                               Kokkos::exp(((Kokkos::log(CDIST1(((K - KMS)))) +
                                             Kokkos::log(morrison_engine::GAMMA(
                                                 workspace, (PGAM(((K - KMS))) + 5.f), morrison))) -
                                            (4.f * Kokkos::log(LAMC(((K - KMS))))))));
        NNUCCC(((K - KMS))) =
            ((((((2.f * morrison.PI) * DAP(((K - KMS)))) * NACNT) * CDIST1(((K - KMS)))) *
              morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 2.f), morrison)) /
             LAMC(((K - KMS))));
        MNUCCC(((K - KMS))) =
            (MNUCCC(((K - KMS))) +
             ((morrison.CONS39 *
               Kokkos::exp(((Kokkos::log(CDIST1(((K - KMS)))) +
                             Kokkos::log(morrison_engine::GAMMA(
                                 workspace, (7.f + PGAM(((K - KMS)))), morrison))) -
                            (6.f * Kokkos::log(LAMC(((K - KMS)))))))) *
              Kokkos::exp((morrison.AIMM * (morrison.TMELT - T3D(((K - KMS))))))));
        NNUCCC(((K - KMS))) =
            (NNUCCC(((K - KMS))) +
             ((morrison.CONS40 *
               Kokkos::exp(((Kokkos::log(CDIST1(((K - KMS)))) +
                             Kokkos::log(morrison_engine::GAMMA(
                                 workspace, (PGAM(((K - KMS))) + 4.f), morrison))) -
                            (3.f * Kokkos::log(LAMC(((K - KMS)))))))) *
              Kokkos::exp((morrison.AIMM * (morrison.TMELT - T3D(((K - KMS))))))));
        NNUCCC(((K - KMS))) =
            ((NNUCCC(((K - KMS)))) < ((NC3D(((K - KMS))) / DT)) ? (NNUCCC(((K - KMS))))
                                                                : ((NC3D(((K - KMS))) / DT)));
      }
      NNUCCC(((K - KMS))) = (NNUCCC(((K - KMS))) * morrison.NNUCCC_REDUCE_COEF);
      if ((QC3D(((K - KMS))) >= 1.E-6f)) {
        if ((morrison.IRAIN == 0)) {
          PRC(((K - KMS))) =
              ((1350.f * Kokkos::pow(QC3D(((K - KMS))), 2.47f)) *
               Kokkos::pow(((NC3D(((K - KMS))) / 1.e6f) * RHO(((K - KMS)))), (-1.79f)));
          NPRC1(((K - KMS))) = (PRC(((K - KMS))) / morrison.CONS29);
          NPRC(((K - KMS)))  = (PRC(((K - KMS))) / (QC3D(((K - KMS))) / NC3D(((K - KMS)))));
          NPRC(((K - KMS))) =
              ((NPRC(((K - KMS)))) < ((NC3D(((K - KMS))) / DT)) ? (NPRC(((K - KMS))))
                                                                : ((NC3D(((K - KMS))) / DT)));
        } else {
          if ((morrison.IRAIN == 1)) {
            DUM  = (1.f - (QC3D(((K - KMS))) / (QC3D(((K - KMS))) + QR3D(((K - KMS))))));
            DUM1 = ((600.f * Kokkos::pow(DUM, 0.68f)) *
                    (((1.f - Kokkos::pow(DUM, 0.68f)) * (1.f - Kokkos::pow(DUM, 0.68f))) *
                     (1.f - Kokkos::pow(DUM, 0.68f))));
            PRC(((K - KMS))) =
                (((((((((9.44e9f / (20.f * 2.6e-7f)) * (nu(((K - KMS))) + 2.f)) *
                       (nu(((K - KMS))) + 4.f)) /
                      ((nu(((K - KMS))) + 1.f) * (nu(((K - KMS))) + 1.f))) *
                     Kokkos::pow(((RHO(((K - KMS))) * QC3D(((K - KMS)))) / 1000.f), 4)) /
                    (((RHO(((K - KMS))) * NC3D(((K - KMS)))) / 1.e6f) *
                     ((RHO(((K - KMS))) * NC3D(((K - KMS)))) / 1.e6f))) *
                   (1.f + (DUM1 / ((1.f - DUM) * (1.f - DUM))))) *
                  1000.f) /
                 RHO(((K - KMS))));
            NPRC(((K - KMS)))  = (((PRC(((K - KMS))) * 2.f) / 2.6e-7f) * 1000.f);
            NPRC1(((K - KMS))) = (0.5f * NPRC(((K - KMS))));
          }
        }
      }
      if ((QNI3D(((K - KMS))) >= 1.E-8f)) {
        NSAGG(((K - KMS))) =
            (((((morrison.CONS15 * ASN(((K - KMS)))) *
                Kokkos::pow(RHO(((K - KMS))), ((2.f + morrison.BS) / 3.f))) *
               Kokkos::pow(QNI3D(((K - KMS))), ((2.f + morrison.BS) / 3.f))) *
              Kokkos::pow((NS3D(((K - KMS))) * RHO(((K - KMS)))), ((4.f - morrison.BS) / 3.f))) /
             RHO(((K - KMS))));
      }
      if (((QNI3D(((K - KMS))) >= 1.E-8f) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
        PSACWS(((K - KMS))) =
            (((((morrison.CONS13 * ASN(((K - KMS)))) * QC3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0S(((K - KMS)))) /
             Kokkos::pow(LAMS(((K - KMS))), (morrison.BS + 3.f)));
        NPSACWS(((K - KMS))) =
            (((((morrison.CONS13 * ASN(((K - KMS)))) * NC3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0S(((K - KMS)))) /
             Kokkos::pow(LAMS(((K - KMS))), (morrison.BS + 3.f)));
      }
      if (((QG3D(((K - KMS))) >= 1.E-8f) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
        PSACWG(((K - KMS))) =
            (((((morrison.CONS14 * AGN(((K - KMS)))) * QC3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0G(((K - KMS)))) /
             Kokkos::pow(LAMG(((K - KMS))), (morrison.BG + 3.f)));
        NPSACWG(((K - KMS))) =
            (((((morrison.CONS14 * AGN(((K - KMS)))) * NC3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0G(((K - KMS)))) /
             Kokkos::pow(LAMG(((K - KMS))), (morrison.BG + 3.f)));
      }
      if (((QI3D(((K - KMS))) >= 1.E-8f) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
        if (((1.f / LAMI(((K - KMS)))) >= 100.E-6f)) {
          PSACWI(((K - KMS))) =
              (((((morrison.CONS16 * AIN(((K - KMS)))) * QC3D(((K - KMS)))) * RHO(((K - KMS)))) *
                N0I(((K - KMS)))) /
               Kokkos::pow(LAMI(((K - KMS))), (morrison.BI + 3.f)));
          NPSACWI(((K - KMS))) =
              (((((morrison.CONS16 * AIN(((K - KMS)))) * NC3D(((K - KMS)))) * RHO(((K - KMS)))) *
                N0I(((K - KMS)))) /
               Kokkos::pow(LAMI(((K - KMS))), (morrison.BI + 3.f)));
        }
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QNI3D(((K - KMS))) >= 1.E-8f))) {
        UMS = ((ASN(((K - KMS))) * morrison.CONS3) / Kokkos::pow(LAMS(((K - KMS))), morrison.BS));
        UMR = ((ARN(((K - KMS))) * morrison.CONS4) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        UNS = ((ASN(((K - KMS))) * morrison.CONS5) / Kokkos::pow(LAMS(((K - KMS))), morrison.BS));
        UNR = ((ARN(((K - KMS))) * morrison.CONS6) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        DUM = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
        UMS = ((UMS) < ((1.2f * DUM)) ? (UMS) : ((1.2f * DUM)));
        UNS = ((UNS) < ((1.2f * DUM)) ? (UNS) : ((1.2f * DUM)));
        UMR = ((UMR) < ((9.1f * DUM)) ? (UMR) : ((9.1f * DUM)));
        UNR = ((UNR) < ((9.1f * DUM)) ? (UNR) : ((9.1f * DUM)));
        PRACS(((K - KMS))) =
            (morrison.CONS41 *
             (((((Kokkos::pow(((((1.2f * UMR) - (0.95f * UMS)) * ((1.2f * UMR) - (0.95f * UMS))) +
                               ((0.08f * UMS) * UMR)),
                              0.5f) *
                  RHO(((K - KMS)))) *
                 N0RR(((K - KMS)))) *
                N0S(((K - KMS)))) /
               ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) *
              (((5.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                        LAMS(((K - KMS))))) +
                (2.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                        (LAMS(((K - KMS))) * LAMS(((K - KMS))))))) +
               (0.5f / (LAMR(((K - KMS))) *
                        ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS)))))))));
        NPRACS(((K - KMS))) =
            (((((morrison.CONS32 * RHO(((K - KMS)))) *
                Kokkos::pow(((1.7f * ((UNR - UNS) * (UNR - UNS))) + ((0.3f * UNR) * UNS)), 0.5f)) *
               N0RR(((K - KMS)))) *
              N0S(((K - KMS)))) *
             (((1.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                       LAMS(((K - KMS))))) +
               (1.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                       (LAMS(((K - KMS))) * LAMS(((K - KMS))))))) +
              (1.f / (LAMR(((K - KMS))) *
                      ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS))))))));
        PRACS(((K - KMS))) =
            ((PRACS(((K - KMS)))) < ((QR3D(((K - KMS))) / DT)) ? (PRACS(((K - KMS))))
                                                               : ((QR3D(((K - KMS))) / DT)));
        if (((QNI3D(((K - KMS))) >= 0.1E-3f) && (QR3D(((K - KMS))) >= 0.1E-3f))) {
          PSACR(((K - KMS))) =
              (morrison.CONS31 *
               (((((Kokkos::pow(((((1.2f * UMR) - (0.95f * UMS)) * ((1.2f * UMR) - (0.95f * UMS))) +
                                 ((0.08f * UMS) * UMR)),
                                0.5f) *
                    RHO(((K - KMS)))) *
                   N0RR(((K - KMS)))) *
                  N0S(((K - KMS)))) /
                 ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS))))) *
                (((5.f / (((LAMS(((K - KMS))) * LAMS(((K - KMS)))) * LAMS(((K - KMS)))) *
                          LAMR(((K - KMS))))) +
                  (2.f / ((LAMS(((K - KMS))) * LAMS(((K - KMS)))) *
                          (LAMR(((K - KMS))) * LAMR(((K - KMS))))))) +
                 (0.5f / (LAMS(((K - KMS))) *
                          ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))))))));
        }
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QG3D(((K - KMS))) >= 1.E-8f))) {
        UMG = ((AGN(((K - KMS))) * morrison.CONS7) / Kokkos::pow(LAMG(((K - KMS))), morrison.BG));
        UMR = ((ARN(((K - KMS))) * morrison.CONS4) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        UNG = ((AGN(((K - KMS))) * morrison.CONS8) / Kokkos::pow(LAMG(((K - KMS))), morrison.BG));
        UNR = ((ARN(((K - KMS))) * morrison.CONS6) / Kokkos::pow(LAMR(((K - KMS))), morrison.BR));
        DUM = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
        UMG = ((UMG) < ((20.f * DUM)) ? (UMG) : ((20.f * DUM)));
        UNG = ((UNG) < ((20.f * DUM)) ? (UNG) : ((20.f * DUM)));
        UMR = ((UMR) < ((9.1f * DUM)) ? (UMR) : ((9.1f * DUM)));
        UNR = ((UNR) < ((9.1f * DUM)) ? (UNR) : ((9.1f * DUM)));
        PRACG(((K - KMS))) =
            (morrison.CONS41 *
             (((((Kokkos::pow(((((1.2f * UMR) - (0.95f * UMG)) * ((1.2f * UMR) - (0.95f * UMG))) +
                               ((0.08f * UMG) * UMR)),
                              0.5f) *
                  RHO(((K - KMS)))) *
                 N0RR(((K - KMS)))) *
                N0G(((K - KMS)))) /
               ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) *
              (((5.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                        LAMG(((K - KMS))))) +
                (2.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                        (LAMG(((K - KMS))) * LAMG(((K - KMS))))))) +
               (0.5f / (LAMR(((K - KMS))) *
                        ((LAMG(((K - KMS))) * LAMG(((K - KMS)))) * LAMG(((K - KMS)))))))));
        NPRACG(((K - KMS))) =
            (((((morrison.CONS32 * RHO(((K - KMS)))) *
                Kokkos::pow(((1.7f * ((UNR - UNG) * (UNR - UNG))) + ((0.3f * UNR) * UNG)), 0.5f)) *
               N0RR(((K - KMS)))) *
              N0G(((K - KMS)))) *
             (((1.f / (((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))) *
                       LAMG(((K - KMS))))) +
               (1.f / ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) *
                       (LAMG(((K - KMS))) * LAMG(((K - KMS))))))) +
              (1.f / (LAMR(((K - KMS))) *
                      ((LAMG(((K - KMS))) * LAMG(((K - KMS)))) * LAMG(((K - KMS))))))));
        PRACG(((K - KMS))) =
            ((PRACG(((K - KMS)))) < ((QR3D(((K - KMS))) / DT)) ? (PRACG(((K - KMS))))
                                                               : ((QR3D(((K - KMS))) / DT)));
      }
      if ((QNI3D(((K - KMS))) >= 0.1E-3f)) {
        if (((QC3D(((K - KMS))) >= 0.5E-3f) || (QR3D(((K - KMS))) >= 0.1E-3f))) {
          if (((PSACWS(((K - KMS))) > 0.f) || (PRACS(((K - KMS))) > 0.f))) {
            if (((T3D(((K - KMS))) < 270.16f) && (T3D(((K - KMS))) > 265.16f))) {
              if ((T3D(((K - KMS))) > 270.16f)) {
                FMULT = 0.f;
              } else {
                if (((T3D(((K - KMS))) <= 270.16f) && (T3D(((K - KMS))) > 268.16f))) {
                  FMULT = ((270.16f - T3D(((K - KMS)))) / 2.f);
                } else {
                  if (((T3D(((K - KMS))) >= 265.16f) && (T3D(((K - KMS))) <= 268.16f))) {
                    FMULT = ((T3D(((K - KMS))) - 265.16f) / 3.f);
                  } else {
                    if ((T3D(((K - KMS))) < 265.16f)) {
                      FMULT = 0.f;
                    }
                  }
                }
              }
              if ((PSACWS(((K - KMS))) > 0.f)) {
                NMULTS(((K - KMS))) = (((35.E4f * PSACWS(((K - KMS)))) * FMULT) * 1000.f);
                QMULTS(((K - KMS))) = (NMULTS(((K - KMS))) * morrison.MMULT);
                QMULTS(((K - KMS))) =
                    ((QMULTS(((K - KMS)))) < (PSACWS(((K - KMS)))) ? (QMULTS(((K - KMS))))
                                                                   : (PSACWS(((K - KMS)))));
                PSACWS(((K - KMS))) = (PSACWS(((K - KMS))) - QMULTS(((K - KMS))));
              }
              if ((PRACS(((K - KMS))) > 0.f)) {
                NMULTR(((K - KMS))) = (((35.E4f * PRACS(((K - KMS)))) * FMULT) * 1000.f);
                QMULTR(((K - KMS))) = (NMULTR(((K - KMS))) * morrison.MMULT);
                QMULTR(((K - KMS))) =
                    ((QMULTR(((K - KMS)))) < (PRACS(((K - KMS)))) ? (QMULTR(((K - KMS))))
                                                                  : (PRACS(((K - KMS)))));
                PRACS(((K - KMS))) = (PRACS(((K - KMS))) - QMULTR(((K - KMS))));
              }
            }
          }
        }
      }
      if ((QG3D(((K - KMS))) >= 0.1E-3f)) {
        if (((QC3D(((K - KMS))) >= 0.5E-3f) || (QR3D(((K - KMS))) >= 0.1E-3f))) {
          if (((PSACWG(((K - KMS))) > 0.f) || (PRACG(((K - KMS))) > 0.f))) {
            if (((T3D(((K - KMS))) < 270.16f) && (T3D(((K - KMS))) > 265.16f))) {
              if ((T3D(((K - KMS))) > 270.16f)) {
                FMULT = 0.f;
              } else {
                if (((T3D(((K - KMS))) <= 270.16f) && (T3D(((K - KMS))) > 268.16f))) {
                  FMULT = ((270.16f - T3D(((K - KMS)))) / 2.f);
                } else {
                  if (((T3D(((K - KMS))) >= 265.16f) && (T3D(((K - KMS))) <= 268.16f))) {
                    FMULT = ((T3D(((K - KMS))) - 265.16f) / 3.f);
                  } else {
                    if ((T3D(((K - KMS))) < 265.16f)) {
                      FMULT = 0.f;
                    }
                  }
                }
              }
              if ((PSACWG(((K - KMS))) > 0.f)) {
                NMULTG(((K - KMS))) = (((35.E4f * PSACWG(((K - KMS)))) * FMULT) * 1000.f);
                QMULTG(((K - KMS))) = (NMULTG(((K - KMS))) * morrison.MMULT);
                QMULTG(((K - KMS))) =
                    ((QMULTG(((K - KMS)))) < (PSACWG(((K - KMS)))) ? (QMULTG(((K - KMS))))
                                                                   : (PSACWG(((K - KMS)))));
                PSACWG(((K - KMS))) = (PSACWG(((K - KMS))) - QMULTG(((K - KMS))));
              }
              if ((PRACG(((K - KMS))) > 0.f)) {
                NMULTRG(((K - KMS))) = (((35.E4f * PRACG(((K - KMS)))) * FMULT) * 1000.f);
                QMULTRG(((K - KMS))) = (NMULTRG(((K - KMS))) * morrison.MMULT);
                QMULTRG(((K - KMS))) =
                    ((QMULTRG(((K - KMS)))) < (PRACG(((K - KMS)))) ? (QMULTRG(((K - KMS))))
                                                                   : (PRACG(((K - KMS)))));
                PRACG(((K - KMS))) = (PRACG(((K - KMS))) - QMULTRG(((K - KMS))));
              }
            }
          }
        }
      }
      if ((PSACWS(((K - KMS))) > 0.f)) {
        if (((QNI3D(((K - KMS))) >= 0.1E-3f) && (QC3D(((K - KMS))) >= 0.5E-3f))) {
          PGSACW(((K - KMS))) =
              ((PSACWS(((K - KMS)))) <
                       ((((((((morrison.CONS17 * DT) * N0S(((K - KMS)))) * QC3D(((K - KMS)))) *
                            QC3D(((K - KMS)))) *
                           ASN(((K - KMS)))) *
                          ASN(((K - KMS)))) /
                         (RHO(((K - KMS))) *
                          Kokkos::pow(LAMS(((K - KMS))), ((2.f * morrison.BS) + 2.f)))))
                   ? (PSACWS(((K - KMS))))
                   : ((((((((morrison.CONS17 * DT) * N0S(((K - KMS)))) * QC3D(((K - KMS)))) *
                          QC3D(((K - KMS)))) *
                         ASN(((K - KMS)))) *
                        ASN(((K - KMS)))) /
                       (RHO(((K - KMS))) *
                        Kokkos::pow(LAMS(((K - KMS))), ((2.f * morrison.BS) + 2.f))))));
          DUM =
              ((((morrison.RHOSN / (morrison.RHOG - morrison.RHOSN)) * PGSACW(((K - KMS))))) > (0.f)
                   ? (((morrison.RHOSN / (morrison.RHOG - morrison.RHOSN)) * PGSACW(((K - KMS)))))
                   : (0.f));
          NSCNG(((K - KMS))) = ((DUM / morrison.MG0) * RHO(((K - KMS))));
          NSCNG(((K - KMS))) =
              ((NSCNG(((K - KMS)))) < ((NS3D(((K - KMS))) / DT)) ? (NSCNG(((K - KMS))))
                                                                 : ((NS3D(((K - KMS))) / DT)));
          PSACWS(((K - KMS))) = (PSACWS(((K - KMS))) - PGSACW(((K - KMS))));
        }
      }
      if ((PRACS(((K - KMS))) > 0.f)) {
        if (((QNI3D(((K - KMS))) >= 0.1E-3f) && (QR3D(((K - KMS))) >= 0.1E-3f))) {
          DUM = (((morrison.CONS18 * (((4.f / LAMS(((K - KMS)))) * (4.f / LAMS(((K - KMS))))) *
                                      (4.f / LAMS(((K - KMS)))))) *
                  (((4.f / LAMS(((K - KMS)))) * (4.f / LAMS(((K - KMS))))) *
                   (4.f / LAMS(((K - KMS)))))) /
                 (((morrison.CONS18 * (((4.f / LAMS(((K - KMS)))) * (4.f / LAMS(((K - KMS))))) *
                                       (4.f / LAMS(((K - KMS)))))) *
                   (((4.f / LAMS(((K - KMS)))) * (4.f / LAMS(((K - KMS))))) *
                    (4.f / LAMS(((K - KMS)))))) +
                  ((morrison.CONS19 * (((4.f / LAMR(((K - KMS)))) * (4.f / LAMR(((K - KMS))))) *
                                       (4.f / LAMR(((K - KMS)))))) *
                   (((4.f / LAMR(((K - KMS)))) * (4.f / LAMR(((K - KMS))))) *
                    (4.f / LAMR(((K - KMS))))))));
          DUM = ((DUM) < (1.f) ? (DUM) : (1.f));
          DUM = ((DUM) > (0.f) ? (DUM) : (0.f));
          PGRACS(((K - KMS))) = ((1.f - DUM) * PRACS(((K - KMS))));
          NGRACS(((K - KMS))) = ((1.f - DUM) * NPRACS(((K - KMS))));
          NGRACS(((K - KMS))) =
              ((NGRACS(((K - KMS)))) < ((NR3D(((K - KMS))) / DT)) ? (NGRACS(((K - KMS))))
                                                                  : ((NR3D(((K - KMS))) / DT)));
          NGRACS(((K - KMS))) =
              ((NGRACS(((K - KMS)))) < ((NS3D(((K - KMS))) / DT)) ? (NGRACS(((K - KMS))))
                                                                  : ((NS3D(((K - KMS))) / DT)));
          PRACS(((K - KMS)))  = (PRACS(((K - KMS))) - PGRACS(((K - KMS))));
          NPRACS(((K - KMS))) = (NPRACS(((K - KMS))) - NGRACS(((K - KMS))));
          PSACR(((K - KMS)))  = (PSACR(((K - KMS))) * (1.f - DUM));
        }
      }
      if (((T3D(((K - KMS))) < 269.15f) && (QR3D(((K - KMS))) >= morrison.QSMALL))) {
        MNUCCR(((K - KMS))) =
            ((((morrison.CONS20 * NR3D(((K - KMS)))) *
               Kokkos::exp((morrison.AIMM * (morrison.TMELT - T3D(((K - KMS))))))) /
              ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) /
             ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))));
        NNUCCR(((K - KMS))) =
            ((((morrison.PI * NR3D(((K - KMS)))) * morrison.BIMM) *
              Kokkos::exp((morrison.AIMM * (morrison.TMELT - T3D(((K - KMS))))))) /
             ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS)))));
        NNUCCR(((K - KMS))) =
            ((NNUCCR(((K - KMS)))) < ((NR3D(((K - KMS))) / DT)) ? (NNUCCR(((K - KMS))))
                                                                : ((NR3D(((K - KMS))) / DT)));
      }
      if (((QR3D(((K - KMS))) >= 1.E-8f) && (QC3D(((K - KMS))) >= 1.E-8f))) {
        if ((morrison.IRAIN == 0)) {
          DUM               = (QC3D(((K - KMS))) * QR3D(((K - KMS))));
          PRA(((K - KMS)))  = (67.f * Kokkos::pow(DUM, 1.15f));
          NPRA(((K - KMS))) = (PRA(((K - KMS))) / (QC3D(((K - KMS))) / NC3D(((K - KMS)))));
        } else {
          if ((morrison.IRAIN == 1)) {
            DUM  = (1.f - (QC3D(((K - KMS))) / (QC3D(((K - KMS))) + QR3D(((K - KMS))))));
            DUM1 = Kokkos::pow((DUM / (DUM + 5.e-4f)), 4);
            PRA(((K - KMS)))  = (((((5.78e3f * RHO(((K - KMS)))) / 1000.f) * QC3D(((K - KMS)))) *
                                  QR3D(((K - KMS)))) *
                                 DUM1);
            NPRA(((K - KMS))) = ((((((PRA(((K - KMS))) * RHO(((K - KMS)))) / 1000.f) *
                                    ((NC3D(((K - KMS))) * RHO(((K - KMS)))) / 1.e6f)) /
                                   ((QC3D(((K - KMS))) * RHO(((K - KMS)))) / 1000.f)) *
                                  1.e6f) /
                                 RHO(((K - KMS))));
          }
        }
      }
      if ((QR3D(((K - KMS))) >= 1.E-8f)) {
        DUM1 = 300.e-6f;
        if (((1.f / LAMR(((K - KMS)))) < DUM1)) {
          DUM = 1.f;
        } else {
          if (((1.f / LAMR(((K - KMS)))) >= DUM1)) {
            DUM = (2.f - Kokkos::exp((2300.f * ((1.f / LAMR(((K - KMS)))) - DUM1))));
          }
        }
        NRAGG(((K - KMS))) =
            (((((-5.78f) * DUM) * NR3D(((K - KMS)))) * QR3D(((K - KMS)))) * RHO(((K - KMS))));
      }
      if (((QI3D(((K - KMS))) >= 1.E-8f) && (QVQVSI(((K - KMS))) >= 1.f))) {
        NPRCI(((K - KMS))) =
            ((((((morrison.CONS21 * (QV3D(((K - KMS))) - QVI(((K - KMS))))) * RHO(((K - KMS)))) *
                N0I(((K - KMS)))) *
               Kokkos::exp(((-LAMI(((K - KMS)))) * morrison.DCS))) *
              DV(((K - KMS)))) /
             ABI(((K - KMS))));
        PRCI(((K - KMS))) = (morrison.CONS22 * NPRCI(((K - KMS))));
        NPRCI(((K - KMS))) =
            ((NPRCI(((K - KMS)))) < ((NI3D(((K - KMS))) / DT)) ? (NPRCI(((K - KMS))))
                                                               : ((NI3D(((K - KMS))) / DT)));
      }
      if (((QNI3D(((K - KMS))) >= 1.E-8f) && (QI3D(((K - KMS))) >= morrison.QSMALL))) {
        PRAI(((K - KMS))) =
            (((((morrison.CONS23 * ASN(((K - KMS)))) * QI3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0S(((K - KMS)))) /
             Kokkos::pow(LAMS(((K - KMS))), (morrison.BS + 3.f)));
        NPRAI(((K - KMS))) =
            (((((morrison.CONS23 * ASN(((K - KMS)))) * NI3D(((K - KMS)))) * RHO(((K - KMS)))) *
              N0S(((K - KMS)))) /
             Kokkos::pow(LAMS(((K - KMS))), (morrison.BS + 3.f)));
        NPRAI(((K - KMS))) =
            ((NPRAI(((K - KMS)))) < ((NI3D(((K - KMS))) / DT)) ? (NPRAI(((K - KMS))))
                                                               : ((NI3D(((K - KMS))) / DT)));
      }
      if ((((QR3D(((K - KMS))) >= 1.E-8f) && (QI3D(((K - KMS))) >= 1.E-8f)) &&
           (T3D(((K - KMS))) <= morrison.TMELT))) {
        if ((QR3D(((K - KMS))) >= 0.1E-3f)) {
          NIACR(((K - KMS))) =
              (((((morrison.CONS24 * NI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) *
               RHO(((K - KMS))));
          PIACR(((K - KMS))) =
              ((((((morrison.CONS25 * NI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                 Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) /
                ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) *
               RHO(((K - KMS))));
          PRACI(((K - KMS))) =
              (((((morrison.CONS24 * QI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) *
               RHO(((K - KMS))));
          NIACR(((K - KMS))) =
              ((NIACR(((K - KMS)))) < ((NR3D(((K - KMS))) / DT)) ? (NIACR(((K - KMS))))
                                                                 : ((NR3D(((K - KMS))) / DT)));
          NIACR(((K - KMS))) =
              ((NIACR(((K - KMS)))) < ((NI3D(((K - KMS))) / DT)) ? (NIACR(((K - KMS))))
                                                                 : ((NI3D(((K - KMS))) / DT)));
        } else {
          NIACRS(((K - KMS))) =
              (((((morrison.CONS24 * NI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) *
               RHO(((K - KMS))));
          PIACRS(((K - KMS))) =
              ((((((morrison.CONS25 * NI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                 Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) /
                ((LAMR(((K - KMS))) * LAMR(((K - KMS)))) * LAMR(((K - KMS))))) *
               RHO(((K - KMS))));
          PRACIS(((K - KMS))) =
              (((((morrison.CONS24 * QI3D(((K - KMS)))) * N0RR(((K - KMS)))) * ARN(((K - KMS)))) /
                Kokkos::pow(LAMR(((K - KMS))), (morrison.BR + 3.f))) *
               RHO(((K - KMS))));
          NIACRS(((K - KMS))) =
              ((NIACRS(((K - KMS)))) < ((NR3D(((K - KMS))) / DT)) ? (NIACRS(((K - KMS))))
                                                                  : ((NR3D(((K - KMS))) / DT)));
          NIACRS(((K - KMS))) =
              ((NIACRS(((K - KMS)))) < ((NI3D(((K - KMS))) / DT)) ? (NIACRS(((K - KMS))))
                                                                  : ((NI3D(((K - KMS))) / DT)));
        }
      }
      if ((morrison.INUC == 0)) {
        if ((((QVQVS(((K - KMS))) >= 0.999f) && (T3D(((K - KMS))) <= 265.15f)) ||
             (QVQVSI(((K - KMS))) >= 1.08f))) {
          KC2 = ((0.005f * Kokkos::exp((0.304f * (morrison.TMELT - T3D(((K - KMS))))))) * 1000.f);
          KC2 = ((KC2) < (500.e3f) ? (KC2) : (500.e3f));
          KC2 = (((KC2 / RHO(((K - KMS))))) > (0.f) ? ((KC2 / RHO(((K - KMS))))) : (0.f));
          if ((KC2 > ((NI3D(((K - KMS))) + NS3D(((K - KMS)))) + NG3D(((K - KMS)))))) {
            NNUCCD(((K - KMS))) =
                ((((KC2 - NI3D(((K - KMS)))) - NS3D(((K - KMS)))) - NG3D(((K - KMS)))) / DT);
            MNUCCD(((K - KMS))) = (NNUCCD(((K - KMS))) * morrison.MI0);
          }
        }
      } else {
        if ((morrison.INUC == 1)) {
          if (((T3D(((K - KMS))) < morrison.TMELT) && (QVQVSI(((K - KMS))) > 1.f))) {
            KC2 = ((0.16f * 1000.f) / RHO(((K - KMS))));
            if ((KC2 > ((NI3D(((K - KMS))) + NS3D(((K - KMS)))) + NG3D(((K - KMS)))))) {
              NNUCCD(((K - KMS))) =
                  ((((KC2 - NI3D(((K - KMS)))) - NS3D(((K - KMS)))) - NG3D(((K - KMS)))) / DT);
              MNUCCD(((K - KMS))) = (NNUCCD(((K - KMS))) * morrison.MI0);
            }
          }
        }
      }
      NNUCCD(((K - KMS))) = (NNUCCD(((K - KMS))) * morrison.NNUCCD_REDUCE_COEF);
    label_101:;
      if ((QI3D(((K - KMS))) >= morrison.QSMALL)) {
        EPSI = (((((2.f * morrison.PI) * N0I(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) /
                (LAMI(((K - KMS))) * LAMI(((K - KMS)))));
      } else {
        EPSI = 0.f;
      }
      if ((QNI3D(((K - KMS))) >= morrison.QSMALL)) {
        EPSS = (((((2.f * morrison.PI) * N0S(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                ((morrison.F1S / (LAMS(((K - KMS))) * LAMS(((K - KMS))))) +
                 ((((morrison.F2S *
                     Kokkos::pow(((ASN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                    Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                   morrison.CONS10) /
                  Kokkos::pow(LAMS(((K - KMS))), morrison.CONS35))));
      } else {
        EPSS = 0.f;
      }
      if ((QG3D(((K - KMS))) >= morrison.QSMALL)) {
        EPSG = (((((2.f * morrison.PI) * N0G(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                ((morrison.F1S / (LAMG(((K - KMS))) * LAMG(((K - KMS))))) +
                 ((((morrison.F2S *
                     Kokkos::pow(((AGN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                    Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                   morrison.CONS11) /
                  Kokkos::pow(LAMG(((K - KMS))), morrison.CONS36))));
      } else {
        EPSG = 0.f;
      }
      if ((QR3D(((K - KMS))) >= morrison.QSMALL)) {
        EPSR = (((((2.f * morrison.PI) * N0RR(((K - KMS)))) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                ((morrison.F1R / (LAMR(((K - KMS))) * LAMR(((K - KMS))))) +
                 ((((morrison.F2R *
                     Kokkos::pow(((ARN(((K - KMS))) * RHO(((K - KMS)))) / MU(((K - KMS)))), 0.5f)) *
                    Kokkos::pow(SC(((K - KMS))), (1.f / 3.f))) *
                   morrison.CONS9) /
                  Kokkos::pow(LAMR(((K - KMS))), morrison.CONS34))));
      } else {
        EPSR = 0.f;
      }
      if ((QI3D(((K - KMS))) >= morrison.QSMALL)) {
        DUM = (1.f - (Kokkos::exp(((-LAMI(((K - KMS)))) * morrison.DCS)) *
                      (1.f + (LAMI(((K - KMS))) * morrison.DCS))));
        PRD(((K - KMS))) =
            (((EPSI * (QV3D(((K - KMS))) - QVI(((K - KMS))))) / ABI(((K - KMS)))) * DUM);
      } else {
        DUM = 0.f;
      }
      if ((QNI3D(((K - KMS))) >= morrison.QSMALL)) {
        PRDS(((K - KMS))) =
            (((EPSS * (QV3D(((K - KMS))) - QVI(((K - KMS))))) / ABI(((K - KMS)))) +
             (((EPSI * (QV3D(((K - KMS))) - QVI(((K - KMS))))) / ABI(((K - KMS)))) * (1.f - DUM)));
      } else {
        PRD(((K - KMS))) =
            (PRD(((K - KMS))) +
             (((EPSI * (QV3D(((K - KMS))) - QVI(((K - KMS))))) / ABI(((K - KMS)))) * (1.f - DUM)));
      }
      PRDG(((K - KMS))) = ((EPSG * (QV3D(((K - KMS))) - QVI(((K - KMS))))) / ABI(((K - KMS))));
      if ((QV3D(((K - KMS))) < QVS(((K - KMS))))) {
        PRE(((K - KMS))) = ((EPSR * (QV3D(((K - KMS))) - QVS(((K - KMS))))) / AB(((K - KMS))));
        PRE(((K - KMS))) = ((PRE(((K - KMS)))) < (0.f) ? (PRE(((K - KMS)))) : (0.f));
      } else {
        PRE(((K - KMS))) = 0.f;
      }
      DUM    = ((QV3D(((K - KMS))) - QVI(((K - KMS)))) / DT);
      FUDGEF = 0.9999f;
      SUM_DEP =
          (((PRD(((K - KMS))) + PRDS(((K - KMS)))) + MNUCCD(((K - KMS)))) + PRDG(((K - KMS))));
      if ((((DUM > 0.f) && (SUM_DEP > (DUM * FUDGEF))) ||
           ((DUM < 0.f) && (SUM_DEP < (DUM * FUDGEF))))) {
        MNUCCD(((K - KMS))) = (((FUDGEF * MNUCCD(((K - KMS)))) * DUM) / SUM_DEP);
        PRD(((K - KMS)))    = (((FUDGEF * PRD(((K - KMS)))) * DUM) / SUM_DEP);
        PRDS(((K - KMS)))   = (((FUDGEF * PRDS(((K - KMS)))) * DUM) / SUM_DEP);
        PRDG(((K - KMS)))   = (((FUDGEF * PRDG(((K - KMS)))) * DUM) / SUM_DEP);
      }
      if ((PRD(((K - KMS))) < 0.f)) {
        EPRD(((K - KMS))) = PRD(((K - KMS)));
        PRD(((K - KMS)))  = 0.f;
      }
      if ((PRDS(((K - KMS))) < 0.f)) {
        EPRDS(((K - KMS))) = PRDS(((K - KMS)));
        PRDS(((K - KMS)))  = 0.f;
      }
      if ((PRDG(((K - KMS))) < 0.f)) {
        EPRDG(((K - KMS))) = PRDG(((K - KMS)));
        PRDG(((K - KMS)))  = 0.f;
      }
      if ((morrison.ILIQ == 1)) {
        MNUCCC(((K - KMS))) = 0.f;
        NNUCCC(((K - KMS))) = 0.f;
        MNUCCR(((K - KMS))) = 0.f;
        NNUCCR(((K - KMS))) = 0.f;
        MNUCCD(((K - KMS))) = 0.f;
        NNUCCD(((K - KMS))) = 0.f;
      }
      if ((morrison.IGRAUP == 1)) {
        PRACG(((K - KMS)))   = 0.f;
        PSACR(((K - KMS)))   = 0.f;
        PSACWG(((K - KMS)))  = 0.f;
        PGSACW(((K - KMS)))  = 0.f;
        PGRACS(((K - KMS)))  = 0.f;
        PRDG(((K - KMS)))    = 0.f;
        EPRDG(((K - KMS)))   = 0.f;
        EVPMG(((K - KMS)))   = 0.f;
        PGMLT(((K - KMS)))   = 0.f;
        NPRACG(((K - KMS)))  = 0.f;
        NPSACWG(((K - KMS))) = 0.f;
        NSCNG(((K - KMS)))   = 0.f;
        NGRACS(((K - KMS)))  = 0.f;
        NSUBG(((K - KMS)))   = 0.f;
        NGMLTG(((K - KMS)))  = 0.f;
        NGMLTR(((K - KMS)))  = 0.f;
        PIACRS(((K - KMS)))  = (PIACRS(((K - KMS))) + PIACR(((K - KMS))));
        PIACR(((K - KMS)))   = 0.f;
      }
      DUM = (((((((((PRC(((K - KMS))) + PRA(((K - KMS)))) + MNUCCC(((K - KMS)))) +
                   PSACWS(((K - KMS)))) +
                  PSACWI(((K - KMS)))) +
                 QMULTS(((K - KMS)))) +
                PSACWG(((K - KMS)))) +
               PGSACW(((K - KMS)))) +
              QMULTG(((K - KMS)))) *
             DT);
      if (((DUM > QC3D(((K - KMS)))) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO               = (QC3D(((K - KMS))) / DUM);
        PRC(((K - KMS)))    = (PRC(((K - KMS))) * RATIO);
        PRA(((K - KMS)))    = (PRA(((K - KMS))) * RATIO);
        MNUCCC(((K - KMS))) = (MNUCCC(((K - KMS))) * RATIO);
        PSACWS(((K - KMS))) = (PSACWS(((K - KMS))) * RATIO);
        PSACWI(((K - KMS))) = (PSACWI(((K - KMS))) * RATIO);
        QMULTS(((K - KMS))) = (QMULTS(((K - KMS))) * RATIO);
        QMULTG(((K - KMS))) = (QMULTG(((K - KMS))) * RATIO);
        PSACWG(((K - KMS))) = (PSACWG(((K - KMS))) * RATIO);
        PGSACW(((K - KMS))) = (PGSACW(((K - KMS))) * RATIO);
      }
      DUM = ((((((((((((((-PRD(((K - KMS)))) - MNUCCC(((K - KMS)))) + PRCI(((K - KMS)))) +
                       PRAI(((K - KMS)))) -
                      QMULTS(((K - KMS)))) -
                     QMULTG(((K - KMS)))) -
                    QMULTR(((K - KMS)))) -
                   QMULTRG(((K - KMS)))) -
                  MNUCCD(((K - KMS)))) +
                 PRACI(((K - KMS)))) +
                PRACIS(((K - KMS)))) -
               EPRD(((K - KMS)))) -
              PSACWI(((K - KMS)))) *
             DT);
      if (((DUM > QI3D(((K - KMS)))) && (QI3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO = ((((((((((QI3D(((K - KMS))) / DT) + PRD(((K - KMS)))) + MNUCCC(((K - KMS)))) +
                       QMULTS(((K - KMS)))) +
                      QMULTG(((K - KMS)))) +
                     QMULTR(((K - KMS)))) +
                    QMULTRG(((K - KMS)))) +
                   MNUCCD(((K - KMS)))) +
                  PSACWI(((K - KMS)))) /
                 ((((PRCI(((K - KMS))) + PRAI(((K - KMS)))) + PRACI(((K - KMS)))) +
                   PRACIS(((K - KMS)))) -
                  EPRD(((K - KMS)))));
        PRCI(((K - KMS)))   = (PRCI(((K - KMS))) * RATIO);
        PRAI(((K - KMS)))   = (PRAI(((K - KMS))) * RATIO);
        PRACI(((K - KMS)))  = (PRACI(((K - KMS))) * RATIO);
        PRACIS(((K - KMS))) = (PRACIS(((K - KMS))) * RATIO);
        EPRD(((K - KMS)))   = (EPRD(((K - KMS))) * RATIO);
      }
      DUM = ((((((((PRACS(((K - KMS))) - PRE(((K - KMS)))) +
                   ((QMULTR(((K - KMS))) + QMULTRG(((K - KMS)))) - PRC(((K - KMS))))) +
                  (MNUCCR(((K - KMS))) - PRA(((K - KMS))))) +
                 PIACR(((K - KMS)))) +
                PIACRS(((K - KMS)))) +
               PGRACS(((K - KMS)))) +
              PRACG(((K - KMS)))) *
             DT);
      if (((DUM > QR3D(((K - KMS)))) && (QR3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO = ((((QR3D(((K - KMS))) / DT) + PRC(((K - KMS)))) + PRA(((K - KMS)))) /
                 (((((((((-PRE(((K - KMS)))) + QMULTR(((K - KMS)))) + QMULTRG(((K - KMS)))) +
                       PRACS(((K - KMS)))) +
                      MNUCCR(((K - KMS)))) +
                     PIACR(((K - KMS)))) +
                    PIACRS(((K - KMS)))) +
                   PGRACS(((K - KMS)))) +
                  PRACG(((K - KMS)))));
        PRE(((K - KMS)))     = (PRE(((K - KMS))) * RATIO);
        PRACS(((K - KMS)))   = (PRACS(((K - KMS))) * RATIO);
        QMULTR(((K - KMS)))  = (QMULTR(((K - KMS))) * RATIO);
        QMULTRG(((K - KMS))) = (QMULTRG(((K - KMS))) * RATIO);
        MNUCCR(((K - KMS)))  = (MNUCCR(((K - KMS))) * RATIO);
        PIACR(((K - KMS)))   = (PIACR(((K - KMS))) * RATIO);
        PIACRS(((K - KMS)))  = (PIACRS(((K - KMS))) * RATIO);
        PGRACS(((K - KMS)))  = (PGRACS(((K - KMS))) * RATIO);
        PRACG(((K - KMS)))   = (PRACG(((K - KMS))) * RATIO);
      }
      if ((morrison.IGRAUP == 0)) {
        DUM = ((((((((((-PRDS(((K - KMS)))) - PSACWS(((K - KMS)))) - PRAI(((K - KMS)))) -
                     PRCI(((K - KMS)))) -
                    PRACS(((K - KMS)))) -
                   EPRDS(((K - KMS)))) +
                  PSACR(((K - KMS)))) -
                 PIACRS(((K - KMS)))) -
                PRACIS(((K - KMS)))) *
               DT);
        if (((DUM > QNI3D(((K - KMS)))) && (QNI3D(((K - KMS))) >= morrison.QSMALL))) {
          RATIO = (((((((((QNI3D(((K - KMS))) / DT) + PRDS(((K - KMS)))) + PSACWS(((K - KMS)))) +
                        PRAI(((K - KMS)))) +
                       PRCI(((K - KMS)))) +
                      PRACS(((K - KMS)))) +
                     PIACRS(((K - KMS)))) +
                    PRACIS(((K - KMS)))) /
                   ((-EPRDS(((K - KMS)))) + PSACR(((K - KMS)))));
          EPRDS(((K - KMS))) = (EPRDS(((K - KMS))) * RATIO);
          PSACR(((K - KMS))) = (PSACR(((K - KMS))) * RATIO);
        }
      } else {
        if ((morrison.IGRAUP == 1)) {
          DUM = (((((((((((-PRDS(((K - KMS)))) - PSACWS(((K - KMS)))) - PRAI(((K - KMS)))) -
                        PRCI(((K - KMS)))) -
                       PRACS(((K - KMS)))) -
                      EPRDS(((K - KMS)))) +
                     PSACR(((K - KMS)))) -
                    PIACRS(((K - KMS)))) -
                   PRACIS(((K - KMS)))) -
                  MNUCCR(((K - KMS)))) *
                 DT);
          if (((DUM > QNI3D(((K - KMS)))) && (QNI3D(((K - KMS))) >= morrison.QSMALL))) {
            RATIO = ((((((((((QNI3D(((K - KMS))) / DT) + PRDS(((K - KMS)))) + PSACWS(((K - KMS)))) +
                           PRAI(((K - KMS)))) +
                          PRCI(((K - KMS)))) +
                         PRACS(((K - KMS)))) +
                        PIACRS(((K - KMS)))) +
                       PRACIS(((K - KMS)))) +
                      MNUCCR(((K - KMS)))) /
                     ((-EPRDS(((K - KMS)))) + PSACR(((K - KMS)))));
            EPRDS(((K - KMS))) = (EPRDS(((K - KMS))) * RATIO);
            PSACR(((K - KMS))) = (PSACR(((K - KMS))) * RATIO);
          }
        }
      }
      DUM = (((((((((((-PSACWG(((K - KMS)))) - PRACG(((K - KMS)))) - PGSACW(((K - KMS)))) -
                    PGRACS(((K - KMS)))) -
                   PRDG(((K - KMS)))) -
                  MNUCCR(((K - KMS)))) -
                 EPRDG(((K - KMS)))) -
                PIACR(((K - KMS)))) -
               PRACI(((K - KMS)))) -
              PSACR(((K - KMS)))) *
             DT);
      if (((DUM > QG3D(((K - KMS)))) && (QG3D(((K - KMS))) >= morrison.QSMALL))) {
        RATIO = (((((((((((QG3D(((K - KMS))) / DT) + PSACWG(((K - KMS)))) + PRACG(((K - KMS)))) +
                        PGSACW(((K - KMS)))) +
                       PGRACS(((K - KMS)))) +
                      PRDG(((K - KMS)))) +
                     MNUCCR(((K - KMS)))) +
                    PSACR(((K - KMS)))) +
                   PIACR(((K - KMS)))) +
                  PRACI(((K - KMS)))) /
                 (-EPRDG(((K - KMS)))));
        EPRDG(((K - KMS))) = (EPRDG(((K - KMS))) * RATIO);
      }
      QV3DTEN(((K - KMS))) = (QV3DTEN(((K - KMS))) +
                              ((((((((-PRE(((K - KMS)))) - PRD(((K - KMS)))) - PRDS(((K - KMS)))) -
                                   MNUCCD(((K - KMS)))) -
                                  EPRD(((K - KMS)))) -
                                 EPRDS(((K - KMS)))) -
                                PRDG(((K - KMS)))) -
                               EPRDG(((K - KMS)))));
      T3DTEN(((K - KMS))) =
          (T3DTEN(((K - KMS))) +
           ((((PRE(((K - KMS))) * XXLV(((K - KMS)))) +
              (((((((PRD(((K - KMS))) + PRDS(((K - KMS)))) + MNUCCD(((K - KMS)))) +
                   EPRD(((K - KMS)))) +
                  EPRDS(((K - KMS)))) +
                 PRDG(((K - KMS)))) +
                EPRDG(((K - KMS)))) *
               XXLS(((K - KMS))))) +
             (((((((((((((((PSACWS(((K - KMS))) + PSACWI(((K - KMS)))) + MNUCCC(((K - KMS)))) +
                          MNUCCR(((K - KMS)))) +
                         QMULTS(((K - KMS)))) +
                        QMULTG(((K - KMS)))) +
                       QMULTR(((K - KMS)))) +
                      QMULTRG(((K - KMS)))) +
                     PRACS(((K - KMS)))) +
                    PSACWG(((K - KMS)))) +
                   PRACG(((K - KMS)))) +
                  PGSACW(((K - KMS)))) +
                 PGRACS(((K - KMS)))) +
                PIACR(((K - KMS)))) +
               PIACRS(((K - KMS)))) *
              XLF(((K - KMS))))) /
            CPM(((K - KMS)))));
      if ((PCC(((K - KMS))) != 0.f)) {
        std::cerr << std::string("PCC has a nonzero value and will be added to QC3DTEN twice.")
                  << "\n";
        throw std::runtime_error(std::string("Error in Morrison microphysics."));
      }
      QC3DTEN(((K - KMS))) =
          (QC3DTEN(((K - KMS))) +
           ((((((((((-PRA(((K - KMS)))) - PRC(((K - KMS)))) - MNUCCC(((K - KMS)))) +
                  PCC(((K - KMS)))) -
                 PSACWS(((K - KMS)))) -
                PSACWI(((K - KMS)))) -
               QMULTS(((K - KMS)))) -
              QMULTG(((K - KMS)))) -
             PSACWG(((K - KMS)))) -
            PGSACW(((K - KMS)))));
      QI3DTEN(((K - KMS))) =
          (QI3DTEN(((K - KMS))) +
           ((((((((((((PRD(((K - KMS))) + EPRD(((K - KMS)))) + PSACWI(((K - KMS)))) +
                     MNUCCC(((K - KMS)))) -
                    PRCI(((K - KMS)))) -
                   PRAI(((K - KMS)))) +
                  QMULTS(((K - KMS)))) +
                 QMULTG(((K - KMS)))) +
                QMULTR(((K - KMS)))) +
               QMULTRG(((K - KMS)))) +
              MNUCCD(((K - KMS)))) -
             PRACI(((K - KMS)))) -
            PRACIS(((K - KMS)))));
      QR3DTEN(((K - KMS))) = (QR3DTEN(((K - KMS))) +
                              ((((((((((PRE(((K - KMS))) + PRA(((K - KMS)))) + PRC(((K - KMS)))) -
                                      PRACS(((K - KMS)))) -
                                     MNUCCR(((K - KMS)))) -
                                    QMULTR(((K - KMS)))) -
                                   QMULTRG(((K - KMS)))) -
                                  PIACR(((K - KMS)))) -
                                 PIACRS(((K - KMS)))) -
                                PRACG(((K - KMS)))) -
                               PGRACS(((K - KMS)))));
      if ((morrison.IGRAUP == 0)) {
        QNI3DTEN(((K - KMS))) =
            (QNI3DTEN(((K - KMS))) +
             ((((((((PRAI(((K - KMS))) + PSACWS(((K - KMS)))) + PRDS(((K - KMS)))) +
                   PRACS(((K - KMS)))) +
                  PRCI(((K - KMS)))) +
                 EPRDS(((K - KMS)))) -
                PSACR(((K - KMS)))) +
               PIACRS(((K - KMS)))) +
              PRACIS(((K - KMS)))));
        NS3DTEN(((K - KMS))) = (NS3DTEN(((K - KMS))) +
                                ((((NSAGG(((K - KMS))) + NPRCI(((K - KMS)))) - NSCNG(((K - KMS)))) -
                                  NGRACS(((K - KMS)))) +
                                 NIACRS(((K - KMS)))));
        QG3DTEN(((K - KMS))) =
            (QG3DTEN(((K - KMS))) +
             (((((((((PRACG(((K - KMS))) + PSACWG(((K - KMS)))) + PGSACW(((K - KMS)))) +
                    PGRACS(((K - KMS)))) +
                   PRDG(((K - KMS)))) +
                  EPRDG(((K - KMS)))) +
                 MNUCCR(((K - KMS)))) +
                PIACR(((K - KMS)))) +
               PRACI(((K - KMS)))) +
              PSACR(((K - KMS)))));
        NG3DTEN(((K - KMS))) =
            (NG3DTEN(((K - KMS))) +
             (((NSCNG(((K - KMS))) + NGRACS(((K - KMS)))) + NNUCCR(((K - KMS)))) +
              NIACR(((K - KMS)))));
      } else {
        if ((morrison.IGRAUP == 1)) {
          QNI3DTEN(((K - KMS))) =
              (QNI3DTEN(((K - KMS))) +
               (((((((((PRAI(((K - KMS))) + PSACWS(((K - KMS)))) + PRDS(((K - KMS)))) +
                      PRACS(((K - KMS)))) +
                     PRCI(((K - KMS)))) +
                    EPRDS(((K - KMS)))) -
                   PSACR(((K - KMS)))) +
                  PIACRS(((K - KMS)))) +
                 PRACIS(((K - KMS)))) +
                MNUCCR(((K - KMS)))));
          NS3DTEN(((K - KMS))) =
              (NS3DTEN(((K - KMS))) +
               (((((NSAGG(((K - KMS))) + NPRCI(((K - KMS)))) - NSCNG(((K - KMS)))) -
                  NGRACS(((K - KMS)))) +
                 NIACRS(((K - KMS)))) +
                NNUCCR(((K - KMS)))));
        }
      }
      NC3DTEN(((K - KMS))) =
          (NC3DTEN(((K - KMS))) +
           ((((((-NNUCCC(((K - KMS)))) - NPSACWS(((K - KMS)))) - NPRA(((K - KMS)))) -
              NPRC(((K - KMS)))) -
             NPSACWI(((K - KMS)))) -
            NPSACWG(((K - KMS)))));
      NI3DTEN(((K - KMS))) =
          (NI3DTEN(((K - KMS))) +
           (((((((((NNUCCC(((K - KMS))) - NPRCI(((K - KMS)))) - NPRAI(((K - KMS)))) +
                  NMULTS(((K - KMS)))) +
                 NMULTG(((K - KMS)))) +
                NMULTR(((K - KMS)))) +
               NMULTRG(((K - KMS)))) +
              NNUCCD(((K - KMS)))) -
             NIACR(((K - KMS)))) -
            NIACRS(((K - KMS)))));
      NR3DTEN(((K - KMS))) =
          (NR3DTEN(((K - KMS))) +
           (((((((NPRC1(((K - KMS))) - NPRACS(((K - KMS)))) - NNUCCR(((K - KMS)))) +
                NRAGG(((K - KMS)))) -
               NIACR(((K - KMS)))) -
              NIACRS(((K - KMS)))) -
             NPRACG(((K - KMS)))) -
            NGRACS(((K - KMS)))));
      if ((morrison.ISATADJ == 0)) {
        DUMT   = (T3D(((K - KMS))) + (DT * T3DTEN(((K - KMS)))));
        DUMQV  = (QV3D(((K - KMS))) + (DT * QV3DTEN(((K - KMS)))));
        DUM    = (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(DUMT, 0, morrison))
                      ? ((0.99f * PRES(((K - KMS)))))
                      : (morrison_engine::POLYSVP(DUMT, 0, morrison)));
        DUMQSS = ((ep * DUM) / (PRES(((K - KMS))) - DUM));
        DUMQC  = (QC3D(((K - KMS))) + (DT * QC3DTEN(((K - KMS)))));
        DUMQC  = ((DUMQC) > (0.f) ? (DUMQC) : (0.f));
        DUMS   = (DUMQV - DUMQSS);
        PCC(((K - KMS))) = ((DUMS / (1.f + (((XXLV(((K - KMS))) * XXLV(((K - KMS)))) * DUMQSS) /
                                            ((CPM(((K - KMS))) * Rv) * (DUMT * DUMT))))) /
                            DT);
        if (((((PCC(((K - KMS))) * DT) + QC3D(((K - KMS)))) + (DT * QC3DTEN(((K - KMS))))) < 0.f)) {
          PCC(((K - KMS))) = ((-(QC3D(((K - KMS))) + (DT * QC3DTEN(((K - KMS)))))) / DT);
        }
        QV3DTEN(((K - KMS))) = (QV3DTEN(((K - KMS))) - PCC(((K - KMS))));
        T3DTEN(((K - KMS))) =
            (T3DTEN(((K - KMS))) + ((PCC(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QC3DTEN(((K - KMS))) = (QC3DTEN(((K - KMS))) + PCC(((K - KMS))));
      }
      if ((((QC3D(((K - KMS))) + (QC3DTEN(((K - KMS))) * DT)) >= morrison.QSMALL) &&
           (morrison.INUM == 0))) {
        if ((morrison.ISUB == 0)) {
          DUM = (W3D(((K - KMS))) + WVAR(((K - KMS))));
          DUM = ((DUM) > (0.10f) ? (DUM) : (0.10f));
        } else {
          if ((morrison.ISUB == 1)) {
            DUM = W3D(((K - KMS)));
          }
        }
        if ((DUM >= 0.001f)) {
          if ((morrison.IBASE == 1)) {
            IDROP = 0;
            if ((QC3D(((K - KMS))) <= (0.05E-3f / RHO(((K - KMS)))))) {
              IDROP = 1;
            }
            if ((K == 0)) {
              IDROP = 1;
            } else {
              if ((K >= 1)) {
                if (((QC3D(((K - KMS))) > (0.05E-3f / RHO(((K - KMS))))) &&
                     (QC3D(((K - KMS - 1))) <= (0.05E-3f / RHO(((K - KMS - 1))))))) {
                  IDROP = 1;
                }
              }
            }
            if ((IDROP == 1)) {
              if ((morrison.IACT == 1)) {
                DUM  = (DUM * 100.f);
                DUM2 = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                        Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                    (morrison.K1 / (morrison.K1 + 2.f))));
                DUM2 = (DUM2 * 1.E6f);
                DUM2                 = (DUM2 / RHO(((K - KMS))));
                DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
              } else {
                if ((morrison.IACT == 2)) {
                  SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                  AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                           T3D(((K - KMS))));
                  ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                            ((CPM(((K - KMS))) * morrison.RR) *
                             (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                           ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                  GAMM  = (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                  GG  = (1.f /
                         ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                           ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                          (((XXLV(((K - KMS))) * morrison.RHOW) /
                            (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                           (((XXLV(((K - KMS))) * morrison.MW) / (T3D(((K - KMS))) * morrison.RR)) -
                            1.f))));
                  PSI = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                  ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                  ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                  SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                  SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                  DUM1 =
                      ((1.f / (SM1 * SM1)) *
                       ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                        (morrison.F21 * Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                  DUM2 =
                      ((1.f / (SM2 * SM2)) *
                       ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                        (morrison.F22 * Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                  SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                  UU1 = ((2.f * Kokkos::log((SM1 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG1)));
                  UU2 = ((2.f * Kokkos::log((SM2 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG2)));
                  DUM1 = ((morrison.NANEW1 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                  DUM2 = ((morrison.NANEW2 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                  DUM2 = ((DUM1 + DUM2) / RHO(((K - KMS))));
                  DUM2 = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                              ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                              : (DUM2));
                  DUM2 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                  DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                  NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                  NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                }
              }
            } else {
              if ((IDROP == 0)) {
                TAUC = (1.f / ((((((2.f * morrison.PI) * RHO(((K - KMS)))) * DV(((K - KMS)))) *
                                 NC3D(((K - KMS)))) *
                                (PGAM(((K - KMS))) + 1.f)) /
                               LAMC(((K - KMS)))));
                if ((EPSR > 1.E-8f)) {
                  TAUR = (1.f / EPSR);
                } else {
                  TAUR = 1.E8f;
                }
                if ((EPSI > 1.E-8f)) {
                  TAUI = (1.f / EPSI);
                } else {
                  TAUI = 1.E8f;
                }
                if ((EPSS > 1.E-8f)) {
                  TAUS = (1.f / EPSS);
                } else {
                  TAUS = 1.E8f;
                }
                if ((EPSG > 1.E-8f)) {
                  TAUG = (1.f / EPSG);
                } else {
                  TAUG = 1.E8f;
                }
                DUM3 = (((((QVS(((K - KMS))) * RHO(((K - KMS)))) /
                           (PRES(((K - KMS))) - EVS(((K - KMS))))) +
                          (DQSDT / Cp)) *
                         morrison.G) *
                        DUM);
                DUM3 = (((((((DUM3 * TAUC) * TAUR) * TAUI) * TAUS) * TAUG) -
                         ((QVS(((K - KMS))) - QVI(((K - KMS)))) *
                          (((((TAUC * TAUR) * TAUI) * TAUG) + (((TAUC * TAUR) * TAUS) * TAUG)) +
                           (((TAUC * TAUR) * TAUI) * TAUS)))) /
                        (((((((TAUC * TAUR) * TAUI) * TAUG) + (((TAUC * TAUR) * TAUS) * TAUG)) +
                           (((TAUC * TAUR) * TAUI) * TAUS)) +
                          (((TAUR * TAUI) * TAUS) * TAUG)) +
                         (((TAUC * TAUI) * TAUS) * TAUG)));
                if (((DUM3 / QVS(((K - KMS)))) >= 1.E-6f)) {
                  if ((morrison.IACT == 1)) {
                    DUM    = (DUM * 100.f);
                    DUMACT = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                              Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                          (morrison.K1 / (morrison.K1 + 2.f))));
                    DUM3                 = ((DUM3 / QVS(((K - KMS)))) * 100.f);
                    DUM2   = (morrison.C1 * Kokkos::pow(DUM3, morrison.K1));
                    DUM2   = ((DUM2) < (DUMACT) ? (DUM2) : (DUMACT));
                    DUM2   = (DUM2 * 1.E6f);
                    DUM2                 = (DUM2 / RHO(((K - KMS))));
                    DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                    DUM2   = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                    NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                    NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                  } else {
                    if ((morrison.IACT == 2)) {
                      SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                      AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                               T3D(((K - KMS))));
                      ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                                ((CPM(((K - KMS))) * morrison.RR) *
                                 (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                               ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                      GAMM =
                          (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                      GG     = (1.f / ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                                        ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                                       (((XXLV(((K - KMS))) * morrison.RHOW) /
                                         (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                                        (((XXLV(((K - KMS))) * morrison.MW) /
                                          (T3D(((K - KMS))) * morrison.RR)) -
                                         1.f))));
                      PSI  = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                      ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                              ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                      ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                              ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                      SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                              Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                      SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                              Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                      DUM1 = ((1.f / (SM1 * SM1)) *
                              ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                               (morrison.F21 *
                                Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                      DUM2 = ((1.f / (SM2 * SM2)) *
                              ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                               (morrison.F22 *
                                Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                      SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                      UU1  = ((2.f * Kokkos::log((SM1 / SMAX))) /
                              (4.242f * Kokkos::log(morrison.SIG1)));
                      UU2  = ((2.f * Kokkos::log((SM2 / SMAX))) /
                              (4.242f * Kokkos::log(morrison.SIG2)));
                      DUM1   = ((morrison.NANEW1 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                      DUM2   = ((morrison.NANEW2 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                      DUM2   = ((DUM1 + DUM2) / RHO(((K - KMS))));
                      DUMACT = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                                    ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                                    : (DUM2));
                      SIGVL  = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                      AACT   = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                                T3D(((K - KMS))));
                      SM1    = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                                Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                      SM2    = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                                Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                      SMAX   = (DUM3 / QVS(((K - KMS))));
                      UU1    = ((2.f * Kokkos::log((SM1 / SMAX))) /
                                (4.242f * Kokkos::log(morrison.SIG1)));
                      UU2    = ((2.f * Kokkos::log((SM2 / SMAX))) /
                                (4.242f * Kokkos::log(morrison.SIG2)));
                      DUM1   = ((morrison.NANEW1 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                      DUM2   = ((morrison.NANEW2 / 2.f) *
                                (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                      DUM2   = ((DUM1 + DUM2) / RHO(((K - KMS))));
                      DUM2   = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                                    ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                                    : (DUM2));
                      DUM2 = ((DUM2) < (DUMACT) ? (DUM2) : (DUMACT));
                      DUM2   = ((DUM2 - NC3D(((K - KMS)))) / DT);
                      DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                      NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                      NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                    }
                  }
                }
              }
            }
          } else {
            if ((morrison.IBASE == 2)) {
              if ((morrison.IACT == 1)) {
                DUM  = (DUM * 100.f);
                DUM2 = ((0.88f * Kokkos::pow(morrison.C1, (2.f / (morrison.K1 + 2.f)))) *
                        Kokkos::pow((7.E-2f * Kokkos::pow(DUM, 1.5f)),
                                    (morrison.K1 / (morrison.K1 + 2.f))));
                DUM2 = (DUM2 * 1.E6f);
                DUM2                 = (DUM2 / RHO(((K - KMS))));
                DUM2                 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
              } else {
                if ((morrison.IACT == 2)) {
                  SIGVL = (0.0761f - (1.55E-4f * (T3D(((K - KMS))) - morrison.TMELT)));
                  AACT  = ((((2.f * morrison.MW) / (morrison.RHOW * morrison.RR)) * SIGVL) /
                           T3D(((K - KMS))));
                  ALPHA = ((((morrison.G * morrison.MW) * XXLV(((K - KMS)))) /
                            ((CPM(((K - KMS))) * morrison.RR) *
                             (T3D(((K - KMS))) * T3D(((K - KMS)))))) -
                           ((morrison.G * morrison.MA) / (morrison.RR * T3D(((K - KMS))))));
                  GAMM  = (((morrison.RR * T3D(((K - KMS)))) / (EVS(((K - KMS))) * morrison.MW)) +
                           ((morrison.MW * (XXLV(((K - KMS))) * XXLV(((K - KMS))))) /
                            (((CPM(((K - KMS))) * PRES(((K - KMS)))) * morrison.MA) *
                             T3D(((K - KMS))))));
                  GG  = (1.f /
                         ((((morrison.RHOW * morrison.RR) * T3D(((K - KMS)))) /
                           ((EVS(((K - KMS))) * DV(((K - KMS)))) * morrison.MW)) +
                          (((XXLV(((K - KMS))) * morrison.RHOW) /
                            (KAP(((K - KMS))) * T3D(((K - KMS))))) *
                           (((XXLV(((K - KMS))) * morrison.MW) / (T3D(((K - KMS))) * morrison.RR)) -
                            1.f))));
                  PSI = (((2.f / 3.f) * Kokkos::pow(((ALPHA * DUM) / GG), 0.5f)) * AACT);
                  ETA1 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW1));
                  ETA2 = (Kokkos::pow(((ALPHA * DUM) / GG), 1.5f) /
                          ((((2.f * morrison.PI) * morrison.RHOW) * GAMM) * morrison.NANEW2));
                  SM1  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM1)), 1.5f));
                  SM2  = ((2.f / Kokkos::pow(morrison.BACT, 0.5f)) *
                          Kokkos::pow((AACT / (3.f * morrison.RM2)), 1.5f));
                  DUM1 =
                      ((1.f / (SM1 * SM1)) *
                       ((morrison.F11 * Kokkos::pow((PSI / ETA1), 1.5f)) +
                        (morrison.F21 * Kokkos::pow(((SM1 * SM1) / (ETA1 + (3.f * PSI))), 0.75f))));
                  DUM2 =
                      ((1.f / (SM2 * SM2)) *
                       ((morrison.F12 * Kokkos::pow((PSI / ETA2), 1.5f)) +
                        (morrison.F22 * Kokkos::pow(((SM2 * SM2) / (ETA2 + (3.f * PSI))), 0.75f))));
                  SMAX = (1.f / Kokkos::pow((DUM1 + DUM2), 0.5f));
                  UU1 = ((2.f * Kokkos::log((SM1 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG1)));
                  UU2 = ((2.f * Kokkos::log((SM2 / SMAX))) / (4.242f * Kokkos::log(morrison.SIG2)));
                  DUM1 = ((morrison.NANEW1 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU1, morrison)));
                  DUM2 = ((morrison.NANEW2 / 2.f) *
                          (1.f - morrison_engine::DERF1(workspace, UU2, morrison)));
                  DUM2 = ((DUM1 + DUM2) / RHO(((K - KMS))));
                  DUM2 = ((((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))) < (DUM2)
                              ? (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
                              : (DUM2));
                  DUM2 = ((DUM2 - NC3D(((K - KMS)))) / DT);
                  DUM2 = ((0.f) > (DUM2) ? (0.f) : (DUM2));
                  NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + DUM2);
                  NACT(((K - KMS)))    = (NACT(((K - KMS))) + DUM2);
                }
              }
            }
          }
        }
      }
      if ((EPRD(((K - KMS))) < 0.f)) {
        DUM                = ((EPRD(((K - KMS))) * DT) / QI3D(((K - KMS))));
        DUM                = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSUBI(((K - KMS))) = ((DUM * NI3D(((K - KMS)))) / DT);
      }
      if ((EPRDS(((K - KMS))) < 0.f)) {
        DUM                = ((EPRDS(((K - KMS))) * DT) / QNI3D(((K - KMS))));
        DUM                = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSUBS(((K - KMS))) = ((DUM * NS3D(((K - KMS)))) / DT);
      }
      if ((PRE(((K - KMS))) < 0.f)) {
        DUM                = ((PRE(((K - KMS))) * DT) / QR3D(((K - KMS))));
        DUM                = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSUBR(((K - KMS))) = ((DUM * NR3D(((K - KMS)))) / DT);
      }
      if ((EPRDG(((K - KMS))) < 0.f)) {
        DUM                = ((EPRDG(((K - KMS))) * DT) / QG3D(((K - KMS))));
        DUM                = (((-1.f)) > (DUM) ? ((-1.f)) : (DUM));
        NSUBG(((K - KMS))) = ((DUM * NG3D(((K - KMS)))) / DT);
      }
      NI3DTEN(((K - KMS))) = (NI3DTEN(((K - KMS))) + NSUBI(((K - KMS))));
      NS3DTEN(((K - KMS))) = (NS3DTEN(((K - KMS))) + NSUBS(((K - KMS))));
      NG3DTEN(((K - KMS))) = (NG3DTEN(((K - KMS))) + NSUBG(((K - KMS))));
      NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) + NSUBR(((K - KMS))));
    }
    LTRUE = 1;
  label_200:;
    if ((CF3D(((K - KMS))) > morrison.cloud_frac_thresh)) {
      T3D(((K - KMS)))     = (T3D_INIT + ((T3D(((K - KMS))) - T3D_INIT) * CF3D(((K - KMS)))));
      T3DTEN(((K - KMS)))  = (T3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      QV3D(((K - KMS)))    = (QV_INIT + ((QV3D(((K - KMS))) - QSAT_INIT) * CF3D(((K - KMS)))));
      QV3DTEN(((K - KMS))) = (QV3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      QC3D(((K - KMS)))    = (QC3D(((K - KMS))) * CF3D(((K - KMS))));
      QC3DTEN(((K - KMS))) = (QC3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      if ((morrison.INUM == 0)) {
        NC3D(((K - KMS)))    = (NC3D(((K - KMS))) * CF3D(((K - KMS))));
        NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      }
      QR3D(((K - KMS)))    = (QR3D(((K - KMS))) * CF3D(((K - KMS))));
      QR3DTEN(((K - KMS))) = (QR3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      NR3D(((K - KMS)))    = (NR3D(((K - KMS))) * CF3D(((K - KMS))));
      NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      if ((morrison.ILIQ == 0)) {
        QI3D(((K - KMS)))     = (QI3D(((K - KMS))) * CF3D(((K - KMS))));
        QI3DTEN(((K - KMS)))  = (QI3DTEN(((K - KMS))) * CF3D(((K - KMS))));
        NI3D(((K - KMS)))     = (NI3D(((K - KMS))) * CF3D(((K - KMS))));
        NI3DTEN(((K - KMS)))  = (NI3DTEN(((K - KMS))) * CF3D(((K - KMS))));
        QNI3D(((K - KMS)))    = (QNI3D(((K - KMS))) * CF3D(((K - KMS))));
        QNI3DTEN(((K - KMS))) = (QNI3DTEN(((K - KMS))) * CF3D(((K - KMS))));
        NS3D(((K - KMS)))     = (NS3D(((K - KMS))) * CF3D(((K - KMS))));
        NS3DTEN(((K - KMS)))  = (NS3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      }
      if ((morrison.IGRAUP == 0)) {
        QG3D(((K - KMS)))    = (QG3D(((K - KMS))) * CF3D(((K - KMS))));
        QG3DTEN(((K - KMS))) = (QG3DTEN(((K - KMS))) * CF3D(((K - KMS))));
        NG3D(((K - KMS)))    = (NG3D(((K - KMS))) * CF3D(((K - KMS))));
        NG3DTEN(((K - KMS))) = (NG3DTEN(((K - KMS))) * CF3D(((K - KMS))));
      }
      PRC(((K - KMS)))        = (PRC(((K - KMS))) * CF3D(((K - KMS))));
      PRA(((K - KMS)))        = (PRA(((K - KMS))) * CF3D(((K - KMS))));
      PSMLT(((K - KMS)))      = (PSMLT(((K - KMS))) * CF3D(((K - KMS))));
      EVPMS(((K - KMS)))      = (EVPMS(((K - KMS))) * CF3D(((K - KMS))));
      PRACS(((K - KMS)))      = (PRACS(((K - KMS))) * CF3D(((K - KMS))));
      EVPMG(((K - KMS)))      = (EVPMG(((K - KMS))) * CF3D(((K - KMS))));
      PRACG(((K - KMS)))      = (PRACG(((K - KMS))) * CF3D(((K - KMS))));
      PRE(((K - KMS)))        = (PRE(((K - KMS))) * CF3D(((K - KMS))));
      PGMLT(((K - KMS)))      = (PGMLT(((K - KMS))) * CF3D(((K - KMS))));
      MNUCCC(((K - KMS)))     = (MNUCCC(((K - KMS))) * CF3D(((K - KMS))));
      PSACWS(((K - KMS)))     = (PSACWS(((K - KMS))) * CF3D(((K - KMS))));
      PSACWI(((K - KMS)))     = (PSACWI(((K - KMS))) * CF3D(((K - KMS))));
      QMULTS(((K - KMS)))     = (QMULTS(((K - KMS))) * CF3D(((K - KMS))));
      QMULTG(((K - KMS)))     = (QMULTG(((K - KMS))) * CF3D(((K - KMS))));
      PSACWG(((K - KMS)))     = (PSACWG(((K - KMS))) * CF3D(((K - KMS))));
      PGSACW(((K - KMS)))     = (PGSACW(((K - KMS))) * CF3D(((K - KMS))));
      PRD(((K - KMS)))        = (PRD(((K - KMS))) * CF3D(((K - KMS))));
      PRCI(((K - KMS)))       = (PRCI(((K - KMS))) * CF3D(((K - KMS))));
      PRAI(((K - KMS)))       = (PRAI(((K - KMS))) * CF3D(((K - KMS))));
      QMULTR(((K - KMS)))     = (QMULTR(((K - KMS))) * CF3D(((K - KMS))));
      QMULTRG(((K - KMS)))    = (QMULTRG(((K - KMS))) * CF3D(((K - KMS))));
      MNUCCD(((K - KMS)))     = (MNUCCD(((K - KMS))) * CF3D(((K - KMS))));
      PRACI(((K - KMS)))      = (PRACI(((K - KMS))) * CF3D(((K - KMS))));
      PRACIS(((K - KMS)))     = (PRACIS(((K - KMS))) * CF3D(((K - KMS))));
      EPRD(((K - KMS)))       = (EPRD(((K - KMS))) * CF3D(((K - KMS))));
      MNUCCR(((K - KMS)))     = (MNUCCR(((K - KMS))) * CF3D(((K - KMS))));
      PIACR(((K - KMS)))      = (PIACR(((K - KMS))) * CF3D(((K - KMS))));
      PIACRS(((K - KMS)))     = (PIACRS(((K - KMS))) * CF3D(((K - KMS))));
      PGRACS(((K - KMS)))     = (PGRACS(((K - KMS))) * CF3D(((K - KMS))));
      PRDS(((K - KMS)))       = (PRDS(((K - KMS))) * CF3D(((K - KMS))));
      EPRDS(((K - KMS)))      = (EPRDS(((K - KMS))) * CF3D(((K - KMS))));
      PSACR(((K - KMS)))      = (PSACR(((K - KMS))) * CF3D(((K - KMS))));
      PRDG(((K - KMS)))       = (PRDG(((K - KMS))) * CF3D(((K - KMS))));
      EPRDG(((K - KMS)))      = (EPRDG(((K - KMS))) * CF3D(((K - KMS))));
      NPRC1(((K - KMS)))      = (NPRC1(((K - KMS))) * CF3D(((K - KMS))));
      NRAGG(((K - KMS)))      = (NRAGG(((K - KMS))) * CF3D(((K - KMS))));
      NPRACG(((K - KMS)))     = (NPRACG(((K - KMS))) * CF3D(((K - KMS))));
      NSUBR(((K - KMS)))      = (NSUBR(((K - KMS))) * CF3D(((K - KMS))));
      NSMLTR(((K - KMS)))     = (NSMLTR(((K - KMS))) * CF3D(((K - KMS))));
      NGMLTR(((K - KMS)))     = (NGMLTR(((K - KMS))) * CF3D(((K - KMS))));
      NPRACS(((K - KMS)))     = (NPRACS(((K - KMS))) * CF3D(((K - KMS))));
      NNUCCR(((K - KMS)))     = (NNUCCR(((K - KMS))) * CF3D(((K - KMS))));
      NIACR(((K - KMS)))      = (NIACR(((K - KMS))) * CF3D(((K - KMS))));
      NIACRS(((K - KMS)))     = (NIACRS(((K - KMS))) * CF3D(((K - KMS))));
      NGRACS(((K - KMS)))     = (NGRACS(((K - KMS))) * CF3D(((K - KMS))));
      NSMLTS(((K - KMS)))     = (NSMLTS(((K - KMS))) * CF3D(((K - KMS))));
      NSAGG(((K - KMS)))      = (NSAGG(((K - KMS))) * CF3D(((K - KMS))));
      NPRCI(((K - KMS)))      = (NPRCI(((K - KMS))) * CF3D(((K - KMS))));
      NSCNG(((K - KMS)))      = (NSCNG(((K - KMS))) * CF3D(((K - KMS))));
      NSUBS(((K - KMS)))      = (NSUBS(((K - KMS))) * CF3D(((K - KMS))));
      PCC(((K - KMS)))        = (PCC(((K - KMS))) * CF3D(((K - KMS))));
      NNUCCC(((K - KMS)))     = (NNUCCC(((K - KMS))) * CF3D(((K - KMS))));
      NPSACWS(((K - KMS)))    = (NPSACWS(((K - KMS))) * CF3D(((K - KMS))));
      NPRA(((K - KMS)))       = (NPRA(((K - KMS))) * CF3D(((K - KMS))));
      NPRC(((K - KMS)))       = (NPRC(((K - KMS))) * CF3D(((K - KMS))));
      NPSACWI(((K - KMS)))    = (NPSACWI(((K - KMS))) * CF3D(((K - KMS))));
      NPSACWG(((K - KMS)))    = (NPSACWG(((K - KMS))) * CF3D(((K - KMS))));
      NPRAI(((K - KMS)))      = (NPRAI(((K - KMS))) * CF3D(((K - KMS))));
      NMULTS(((K - KMS)))     = (NMULTS(((K - KMS))) * CF3D(((K - KMS))));
      NMULTG(((K - KMS)))     = (NMULTG(((K - KMS))) * CF3D(((K - KMS))));
      NMULTR(((K - KMS)))     = (NMULTR(((K - KMS))) * CF3D(((K - KMS))));
      NMULTRG(((K - KMS)))    = (NMULTRG(((K - KMS))) * CF3D(((K - KMS))));
      NNUCCD(((K - KMS)))     = (NNUCCD(((K - KMS))) * CF3D(((K - KMS))));
      NSUBI(((K - KMS)))      = (NSUBI(((K - KMS))) * CF3D(((K - KMS))));
      NGMLTG(((K - KMS)))     = (NGMLTG(((K - KMS))) * CF3D(((K - KMS))));
      NSUBG(((K - KMS)))      = (NSUBG(((K - KMS))) * CF3D(((K - KMS))));
      NACT(((K - KMS)))       = (NACT(((K - KMS))) * CF3D(((K - KMS))));
      NEGFIX_NI(((K - KMS)))  = (NEGFIX_NI(((K - KMS))) * CF3D(((K - KMS))));
      NEGFIX_NS(((K - KMS)))  = (NEGFIX_NS(((K - KMS))) * CF3D(((K - KMS))));
      NEGFIX_NC(((K - KMS)))  = (NEGFIX_NC(((K - KMS))) * CF3D(((K - KMS))));
      NEGFIX_NR(((K - KMS)))  = (NEGFIX_NR(((K - KMS))) * CF3D(((K - KMS))));
      NEGFIX_NG(((K - KMS)))  = (NEGFIX_NG(((K - KMS))) * CF3D(((K - KMS))));
      SIZEFIX_NR(((K - KMS))) = (SIZEFIX_NR(((K - KMS))) * CF3D(((K - KMS))));
      SIZEFIX_NC(((K - KMS))) = (SIZEFIX_NC(((K - KMS))) * CF3D(((K - KMS))));
      SIZEFIX_NI(((K - KMS))) = (SIZEFIX_NI(((K - KMS))) * CF3D(((K - KMS))));
      SIZEFIX_NS(((K - KMS))) = (SIZEFIX_NS(((K - KMS))) * CF3D(((K - KMS))));
      SIZEFIX_NG(((K - KMS))) = (SIZEFIX_NG(((K - KMS))) * CF3D(((K - KMS))));
      QC_INST(((K - KMS)))    = (QC_INST(((K - KMS))) * CF3D(((K - KMS))));
      QR_INST(((K - KMS)))    = (QR_INST(((K - KMS))) * CF3D(((K - KMS))));
      QI_INST(((K - KMS)))    = (QI_INST(((K - KMS))) * CF3D(((K - KMS))));
      QS_INST(((K - KMS)))    = (QS_INST(((K - KMS))) * CF3D(((K - KMS))));
      QG_INST(((K - KMS)))    = (QG_INST(((K - KMS))) * CF3D(((K - KMS))));
      NC_INST(((K - KMS)))    = (NC_INST(((K - KMS))) * CF3D(((K - KMS))));
      NR_INST(((K - KMS)))    = (NR_INST(((K - KMS))) * CF3D(((K - KMS))));
      NI_INST(((K - KMS)))    = (NI_INST(((K - KMS))) * CF3D(((K - KMS))));
      NS_INST(((K - KMS)))    = (NS_INST(((K - KMS))) * CF3D(((K - KMS))));
      NG_INST(((K - KMS)))    = (NG_INST(((K - KMS))) * CF3D(((K - KMS))));
    }
    NACT(((K - KMS))) = (NACT(((K - KMS))) / DT);
  }
  PRECRT = 0.f;
  SNOWRT = 0.f;
  if ((LTRUE == 0)) {
    goto label_400;
  }
  NSTEP = 1;
  for (K = KTE; K >= KTS; K += (-1)) {
    DUMI(((K - KMS)))   = (QI3D(((K - KMS))) + (QI3DTEN(((K - KMS))) * DT));
    DUMQS(((K - KMS)))  = (QNI3D(((K - KMS))) + (QNI3DTEN(((K - KMS))) * DT));
    DUMR(((K - KMS)))   = (QR3D(((K - KMS))) + (QR3DTEN(((K - KMS))) * DT));
    DUMFNI(((K - KMS))) = (NI3D(((K - KMS))) + (NI3DTEN(((K - KMS))) * DT));
    DUMFNS(((K - KMS))) = (NS3D(((K - KMS))) + (NS3DTEN(((K - KMS))) * DT));
    DUMFNR(((K - KMS))) = (NR3D(((K - KMS))) + (NR3DTEN(((K - KMS))) * DT));
    DUMC(((K - KMS)))   = (QC3D(((K - KMS))) + (QC3DTEN(((K - KMS))) * DT));
    DUMFNC(((K - KMS))) = (NC3D(((K - KMS))) + (NC3DTEN(((K - KMS))) * DT));
    DUMG(((K - KMS)))   = (QG3D(((K - KMS))) + (QG3DTEN(((K - KMS))) * DT));
    DUMFNG(((K - KMS))) = (NG3D(((K - KMS))) + (NG3DTEN(((K - KMS))) * DT));
    if ((morrison.INUM == 1)) {
      DUMFNC(((K - KMS))) = NC3D(((K - KMS)));
    }
    DUMFNI(((K - KMS))) = ((0.f) > (DUMFNI(((K - KMS)))) ? (0.f) : (DUMFNI(((K - KMS)))));
    DUMFNS(((K - KMS))) = ((0.f) > (DUMFNS(((K - KMS)))) ? (0.f) : (DUMFNS(((K - KMS)))));
    DUMFNC(((K - KMS))) = ((0.f) > (DUMFNC(((K - KMS)))) ? (0.f) : (DUMFNC(((K - KMS)))));
    DUMFNR(((K - KMS))) = ((0.f) > (DUMFNR(((K - KMS)))) ? (0.f) : (DUMFNR(((K - KMS)))));
    DUMFNG(((K - KMS))) = ((0.f) > (DUMFNG(((K - KMS)))) ? (0.f) : (DUMFNG(((K - KMS)))));
    if ((DUMI(((K - KMS))) >= morrison.QSMALL)) {
      DLAMI = Kokkos::pow(((morrison.CONS12 * DUMFNI(((K - KMS)))) / DUMI(((K - KMS)))),
                          (1.f / morrison.DI));
      DLAMI = ((DLAMI) > (morrison.LAMMINI) ? (DLAMI) : (morrison.LAMMINI));
      DLAMI = ((DLAMI) < (morrison.LAMMAXI) ? (DLAMI) : (morrison.LAMMAXI));
    }
    if ((DUMR(((K - KMS))) >= morrison.QSMALL)) {
      DLAMR = Kokkos::pow(
          (((morrison.PI * morrison.RHOW) * DUMFNR(((K - KMS)))) / DUMR(((K - KMS)))), (1.f / 3.f));
      DLAMR = ((DLAMR) > (morrison.LAMMINR) ? (DLAMR) : (morrison.LAMMINR));
      DLAMR = ((DLAMR) < (morrison.LAMMAXR) ? (DLAMR) : (morrison.LAMMAXR));
    }
    if ((DUMC(((K - KMS))) >= morrison.QSMALL)) {
      if (morrison.dofix_pgam) {
        PGAM(((K - KMS))) = morrison.pgam_fixed;
      } else {
        PGAM(((K - KMS))) =
            ((0.0005714f * ((NC3D(((K - KMS))) / 1.E6f) * RHO(((K - KMS))))) + 0.2714f);
        PGAM(((K - KMS))) = ((1.f / (PGAM(((K - KMS))) * PGAM(((K - KMS))))) - 1.f);
        PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) > (2.f) ? (PGAM(((K - KMS)))) : (2.f));
        PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) < (10.f) ? (PGAM(((K - KMS)))) : (10.f));
      }
      DLAMC =
          Kokkos::pow((((morrison.CONS26 * DUMFNC(((K - KMS)))) *
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)) /
                       (DUMC(((K - KMS))) *
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))),
                      (1.f / 3.f));
      LAMMIN = ((PGAM(((K - KMS))) + 1.f) / 60.E-6f);
      LAMMAX = ((PGAM(((K - KMS))) + 1.f) / 1.E-6f);
      DLAMC  = ((DLAMC) > (LAMMIN) ? (DLAMC) : (LAMMIN));
      DLAMC  = ((DLAMC) < (LAMMAX) ? (DLAMC) : (LAMMAX));
    }
    if ((DUMQS(((K - KMS))) >= morrison.QSMALL)) {
      DLAMS = Kokkos::pow(((morrison.CONS1 * DUMFNS(((K - KMS)))) / DUMQS(((K - KMS)))),
                          (1.f / morrison.DS));
      DLAMS = ((DLAMS) > (morrison.LAMMINS) ? (DLAMS) : (morrison.LAMMINS));
      DLAMS = ((DLAMS) < (morrison.LAMMAXS) ? (DLAMS) : (morrison.LAMMAXS));
    }
    if ((DUMG(((K - KMS))) >= morrison.QSMALL)) {
      DLAMG = Kokkos::pow(((morrison.CONS2 * DUMFNG(((K - KMS)))) / DUMG(((K - KMS)))),
                          (1.f / morrison.DG));
      DLAMG = ((DLAMG) > (morrison.LAMMING) ? (DLAMG) : (morrison.LAMMING));
      DLAMG = ((DLAMG) < (morrison.LAMMAXG) ? (DLAMG) : (morrison.LAMMAXG));
    }
    if ((DUMC(((K - KMS))) >= morrison.QSMALL)) {
      UNC = ((ACN(((K - KMS))) * morrison_engine::GAMMA(workspace,
                                                        ((1.f + morrison.BC) + PGAM(((K - KMS)))),
                                                        morrison)) /
             (Kokkos::pow(DLAMC, morrison.BC) *
              morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison)));
      UMC = ((ACN(((K - KMS))) * morrison_engine::GAMMA(workspace,
                                                        ((4.f + morrison.BC) + PGAM(((K - KMS)))),
                                                        morrison)) /
             (Kokkos::pow(DLAMC, morrison.BC) *
              morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)));
    } else {
      UMC = 0.f;
      UNC = 0.f;
    }
    if ((DUMI(((K - KMS))) >= morrison.QSMALL)) {
      UNI = ((AIN(((K - KMS))) * morrison.CONS27) / Kokkos::pow(DLAMI, morrison.BI));
      UMI = ((AIN(((K - KMS))) * morrison.CONS28) / Kokkos::pow(DLAMI, morrison.BI));
    } else {
      UMI = 0.f;
      UNI = 0.f;
    }
    if ((DUMR(((K - KMS))) >= morrison.QSMALL)) {
      UNR = ((ARN(((K - KMS))) * morrison.CONS6) / Kokkos::pow(DLAMR, morrison.BR));
      UMR = ((ARN(((K - KMS))) * morrison.CONS4) / Kokkos::pow(DLAMR, morrison.BR));
    } else {
      UMR = 0.f;
      UNR = 0.f;
    }
    if ((DUMQS(((K - KMS))) >= morrison.QSMALL)) {
      UMS = ((ASN(((K - KMS))) * morrison.CONS3) / Kokkos::pow(DLAMS, morrison.BS));
      UNS = ((ASN(((K - KMS))) * morrison.CONS5) / Kokkos::pow(DLAMS, morrison.BS));
    } else {
      UMS = 0.f;
      UNS = 0.f;
    }
    if ((DUMG(((K - KMS))) >= morrison.QSMALL)) {
      UMG = ((AGN(((K - KMS))) * morrison.CONS7) / Kokkos::pow(DLAMG, morrison.BG));
      UNG = ((AGN(((K - KMS))) * morrison.CONS8) / Kokkos::pow(DLAMG, morrison.BG));
    } else {
      UMG = 0.f;
      UNG = 0.f;
    }
    DUM              = Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.54f);
    UMS              = ((UMS) < ((1.2f * DUM)) ? (UMS) : ((1.2f * DUM)));
    UNS              = ((UNS) < ((1.2f * DUM)) ? (UNS) : ((1.2f * DUM)));
    UMI              = ((UMI) < ((1.2f * Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.35f)))
                            ? (UMI)
                            : ((1.2f * Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.35f))));
    UNI              = ((UNI) < ((1.2f * Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.35f)))
                            ? (UNI)
                            : ((1.2f * Kokkos::pow((morrison.RHOSU / RHO(((K - KMS)))), 0.35f))));
    UMR              = ((UMR) < ((9.1f * DUM)) ? (UMR) : ((9.1f * DUM)));
    UNR              = ((UNR) < ((9.1f * DUM)) ? (UNR) : ((9.1f * DUM)));
    UMG              = ((UMG) < ((20.f * DUM)) ? (UMG) : ((20.f * DUM)));
    UNG              = ((UNG) < ((20.f * DUM)) ? (UNG) : ((20.f * DUM)));
    FR(((K - KMS)))  = UMR;
    FI(((K - KMS)))  = UMI;
    FNI(((K - KMS))) = UNI;
    FS(((K - KMS)))  = UMS;
    FNS(((K - KMS))) = UNS;
    FNR(((K - KMS))) = UNR;
    FC(((K - KMS)))  = UMC;
    FNC(((K - KMS))) = UNC;
    FG(((K - KMS)))  = UMG;
    FNG(((K - KMS))) = UNG;
    if ((K < KTE)) {
      if ((FR(((K - KMS))) < 1.E-10f)) {
        FR(((K - KMS))) = FR(((K - KMS + 1)));
      }
      if ((FI(((K - KMS))) < 1.E-10f)) {
        FI(((K - KMS))) = FI(((K - KMS + 1)));
      }
      if ((FNI(((K - KMS))) < 1.E-10f)) {
        FNI(((K - KMS))) = FNI(((K - KMS + 1)));
      }
      if ((FS(((K - KMS))) < 1.E-10f)) {
        FS(((K - KMS))) = FS(((K - KMS + 1)));
      }
      if ((FNS(((K - KMS))) < 1.E-10f)) {
        FNS(((K - KMS))) = FNS(((K - KMS + 1)));
      }
      if ((FNR(((K - KMS))) < 1.E-10f)) {
        FNR(((K - KMS))) = FNR(((K - KMS + 1)));
      }
      if ((FC(((K - KMS))) < 1.E-10f)) {
        FC(((K - KMS))) = FC(((K - KMS + 1)));
      }
      if ((FNC(((K - KMS))) < 1.E-10f)) {
        FNC(((K - KMS))) = FNC(((K - KMS + 1)));
      }
      if ((FG(((K - KMS))) < 1.E-10f)) {
        FG(((K - KMS))) = FG(((K - KMS + 1)));
      }
      if ((FNG(((K - KMS))) < 1.E-10f)) {
        FNG(((K - KMS))) = FNG(((K - KMS + 1)));
      }
    }
    RGVM =
        ((((((((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                               : (FI(((K - KMS)))))) >
                               (FS(((K - KMS))))
                           ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS))))))
                           : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                         ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS))))))
                         : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                       ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS))))))
                       : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                     ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                             ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS))))))
                             : (FNI(((K - KMS))))))
                     : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                   ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                             ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS))))))
                             : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                           ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS))))))
                           : (FNR(((K - KMS))))))
                   : (FNS(((K - KMS)))))) > (FNC(((K - KMS))))
                 ? (((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                             ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS))))))
                             : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                           ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS))))))
                           : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                         ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS))))))
                         : (FNS(((K - KMS))))))
                 : (FNC(((K - KMS)))))) > (FG(((K - KMS))))
               ? (((((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                             ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS))))))
                             : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                           ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS))))))
                           : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                         ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS))))))
                         : (FNS(((K - KMS)))))) > (FNC(((K - KMS))))
                       ? (((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                               ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                       ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                               ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS))))))
                                               : (FNI(((K - KMS))))))
                                       : (FNR(((K - KMS))))))
                               : (FNS(((K - KMS))))))
                       : (FNC(((K - KMS))))))
               : (FG(((K - KMS)))))) > (FNG(((K - KMS))))
             ? (((((((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                     : (FI(((K - KMS)))))) >
                                     (FS(((K - KMS))))
                                 ? (((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS))))))
                                 : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                               ? (((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS))))))
                               : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                             ? (((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS))))))
                             : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                           ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS))))))
                           : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                         ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS))))))
                         : (FNS(((K - KMS)))))) > (FNC(((K - KMS))))
                       ? (((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                               ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                       ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                               ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS))))))
                                               : (FNI(((K - KMS))))))
                                       : (FNR(((K - KMS))))))
                               : (FNS(((K - KMS))))))
                       : (FNC(((K - KMS)))))) > (FG(((K - KMS))))
                     ? (((((((((((((((FR(((K - KMS)))) > (FI(((K - KMS)))) ? (FR(((K - KMS))))
                                                                           : (FI(((K - KMS)))))) >
                                           (FS(((K - KMS))))
                                       ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS))))))
                                       : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                     ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS))))))
                                     : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                   ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS))))))
                                   : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                 ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS))))))
                                 : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                               ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                       ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                               ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS))))))
                                               : (FNI(((K - KMS))))))
                                       : (FNR(((K - KMS))))))
                               : (FNS(((K - KMS)))))) > (FNC(((K - KMS))))
                             ? (((((((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                               ? (FR(((K - KMS))))
                                               : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                             ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS))))))
                                             : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                           ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS))))))
                                           : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                         ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS))))))
                                         : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                       ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                               ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS))))))
                                               : (FNI(((K - KMS))))))
                                       : (FNR(((K - KMS)))))) > (FNS(((K - KMS))))
                                     ? (((((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                     ? (FR(((K - KMS))))
                                                     : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                   ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS))))))
                                                   : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                 ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS))))))
                                                 : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                               ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS))))))
                                               : (FNI(((K - KMS)))))) > (FNR(((K - KMS))))
                                             ? (((((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                           ? (FR(((K - KMS))))
                                                           : (FI(((K - KMS)))))) > (FS(((K - KMS))))
                                                         ? (((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS))))))
                                                         : (FS(((K - KMS)))))) > (FC(((K - KMS))))
                                                       ? (((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS))))))
                                                       : (FC(((K - KMS)))))) > (FNI(((K - KMS))))
                                                     ? (((((((FR(((K - KMS)))) > (FI(((K - KMS))))
                                                                 ? (FR(((K - KMS))))
                                                                 : (FI(((K - KMS)))))) >
                                                                   (FS(((K - KMS))))
                                                               ? (((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS))))))
                                                               : (FS(((K - KMS)))))) >
                                                                 (FC(((K - KMS))))
                                                             ? (((((FR(((K - KMS)))) >
                                                                           (FI(((K - KMS))))
                                                                       ? (FR(((K - KMS))))
                                                                       : (FI(((K - KMS)))))) >
                                                                         (FS(((K - KMS))))
                                                                     ? (((FR(((K - KMS)))) >
                                                                                 (FI(((K - KMS))))
                                                                             ? (FR(((K - KMS))))
                                                                             : (FI(((K - KMS))))))
                                                                     : (FS(((K - KMS))))))
                                                             : (FC(((K - KMS))))))
                                                     : (FNI(((K - KMS))))))
                                             : (FNR(((K - KMS))))))
                                     : (FNS(((K - KMS))))))
                             : (FNC(((K - KMS))))))
                     : (FG(((K - KMS))))))
             : (FNG(((K - KMS)))));
    NSTEP               = ((static_cast<int>((((RGVM * DT) / DZQ(((K - KMS)))) + 1.f))) > NSTEP
                               ? (static_cast<int>((((RGVM * DT) / DZQ(((K - KMS)))) + 1.f)))
                               : NSTEP);
    DUMR(((K - KMS)))   = (DUMR(((K - KMS))) * RHO(((K - KMS))));
    DUMI(((K - KMS)))   = (DUMI(((K - KMS))) * RHO(((K - KMS))));
    DUMFNI(((K - KMS))) = (DUMFNI(((K - KMS))) * RHO(((K - KMS))));
    DUMQS(((K - KMS)))  = (DUMQS(((K - KMS))) * RHO(((K - KMS))));
    DUMFNS(((K - KMS))) = (DUMFNS(((K - KMS))) * RHO(((K - KMS))));
    DUMFNR(((K - KMS))) = (DUMFNR(((K - KMS))) * RHO(((K - KMS))));
    DUMC(((K - KMS)))   = (DUMC(((K - KMS))) * RHO(((K - KMS))));
    DUMFNC(((K - KMS))) = (DUMFNC(((K - KMS))) * RHO(((K - KMS))));
    DUMG(((K - KMS)))   = (DUMG(((K - KMS))) * RHO(((K - KMS))));
    DUMFNG(((K - KMS))) = (DUMFNG(((K - KMS))) * RHO(((K - KMS))));
  }
  for (N = 1; N < NSTEP + 1; N += 1) {
    for (K = KTS; K < KTE + 1; K += 1) {
      FALOUTR(((K - KMS)))  = (FR(((K - KMS))) * DUMR(((K - KMS))));
      FALOUTI(((K - KMS)))  = (FI(((K - KMS))) * DUMI(((K - KMS))));
      FALOUTNI(((K - KMS))) = (FNI(((K - KMS))) * DUMFNI(((K - KMS))));
      FALOUTS(((K - KMS)))  = (FS(((K - KMS))) * DUMQS(((K - KMS))));
      FALOUTNS(((K - KMS))) = (FNS(((K - KMS))) * DUMFNS(((K - KMS))));
      FALOUTNR(((K - KMS))) = (FNR(((K - KMS))) * DUMFNR(((K - KMS))));
      FALOUTC(((K - KMS)))  = (FC(((K - KMS))) * DUMC(((K - KMS))));
      FALOUTNC(((K - KMS))) = (FNC(((K - KMS))) * DUMFNC(((K - KMS))));
      FALOUTG(((K - KMS)))  = (FG(((K - KMS))) * DUMG(((K - KMS))));
      FALOUTNG(((K - KMS))) = (FNG(((K - KMS))) * DUMFNG(((K - KMS))));
    }
    K                    = KTE;
    FALTNDR              = (FALOUTR(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDI              = (FALOUTI(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDNI             = (FALOUTNI(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDS              = (FALOUTS(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDNS             = (FALOUTNS(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDNR             = (FALOUTNR(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDC              = (FALOUTC(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDNC             = (FALOUTNC(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDG              = (FALOUTG(((K - KMS))) / DZQ(((K - KMS))));
    FALTNDNG             = (FALOUTNG(((K - KMS))) / DZQ(((K - KMS))));
    QRSTEN(((K - KMS)))  = (QRSTEN(((K - KMS))) - ((FALTNDR / NSTEP) / RHO(((K - KMS)))));
    QISTEN(((K - KMS)))  = (QISTEN(((K - KMS))) - ((FALTNDI / NSTEP) / RHO(((K - KMS)))));
    NI3DTEN(((K - KMS))) = (NI3DTEN(((K - KMS))) - ((FALTNDNI / NSTEP) / RHO(((K - KMS)))));
    QNISTEN(((K - KMS))) = (QNISTEN(((K - KMS))) - ((FALTNDS / NSTEP) / RHO(((K - KMS)))));
    NS3DTEN(((K - KMS))) = (NS3DTEN(((K - KMS))) - ((FALTNDNS / NSTEP) / RHO(((K - KMS)))));
    NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) - ((FALTNDNR / NSTEP) / RHO(((K - KMS)))));
    QCSTEN(((K - KMS)))  = (QCSTEN(((K - KMS))) - ((FALTNDC / NSTEP) / RHO(((K - KMS)))));
    NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) - ((FALTNDNC / NSTEP) / RHO(((K - KMS)))));
    QGSTEN(((K - KMS)))  = (QGSTEN(((K - KMS))) - ((FALTNDG / NSTEP) / RHO(((K - KMS)))));
    NG3DTEN(((K - KMS))) = (NG3DTEN(((K - KMS))) - ((FALTNDNG / NSTEP) / RHO(((K - KMS)))));
    NISTEN(((K - KMS)))  = (NISTEN(((K - KMS))) - ((FALTNDNI / NSTEP) / RHO(((K - KMS)))));
    NSSTEN(((K - KMS)))  = (NSSTEN(((K - KMS))) - ((FALTNDNS / NSTEP) / RHO(((K - KMS)))));
    NRSTEN(((K - KMS)))  = (NRSTEN(((K - KMS))) - ((FALTNDNR / NSTEP) / RHO(((K - KMS)))));
    NCSTEN(((K - KMS)))  = (NCSTEN(((K - KMS))) - ((FALTNDNC / NSTEP) / RHO(((K - KMS)))));
    NGSTEN(((K - KMS)))  = (NGSTEN(((K - KMS))) - ((FALTNDNG / NSTEP) / RHO(((K - KMS)))));
    DUMR(((K - KMS)))    = (DUMR(((K - KMS))) - ((FALTNDR * DT) / NSTEP));
    DUMI(((K - KMS)))    = (DUMI(((K - KMS))) - ((FALTNDI * DT) / NSTEP));
    DUMFNI(((K - KMS)))  = (DUMFNI(((K - KMS))) - ((FALTNDNI * DT) / NSTEP));
    DUMQS(((K - KMS)))   = (DUMQS(((K - KMS))) - ((FALTNDS * DT) / NSTEP));
    DUMFNS(((K - KMS)))  = (DUMFNS(((K - KMS))) - ((FALTNDNS * DT) / NSTEP));
    DUMFNR(((K - KMS)))  = (DUMFNR(((K - KMS))) - ((FALTNDNR * DT) / NSTEP));
    DUMC(((K - KMS)))    = (DUMC(((K - KMS))) - ((FALTNDC * DT) / NSTEP));
    DUMFNC(((K - KMS)))  = (DUMFNC(((K - KMS))) - ((FALTNDNC * DT) / NSTEP));
    DUMG(((K - KMS)))    = (DUMG(((K - KMS))) - ((FALTNDG * DT) / NSTEP));
    DUMFNG(((K - KMS)))  = (DUMFNG(((K - KMS))) - ((FALTNDNG * DT) / NSTEP));
    for (K = (KTE - 1); K >= KTS; K += (-1)) {
      FALTNDR  = ((FALOUTR(((K - KMS + 1))) - FALOUTR(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDI  = ((FALOUTI(((K - KMS + 1))) - FALOUTI(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDNI = ((FALOUTNI(((K - KMS + 1))) - FALOUTNI(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDS  = ((FALOUTS(((K - KMS + 1))) - FALOUTS(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDNS = ((FALOUTNS(((K - KMS + 1))) - FALOUTNS(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDNR = ((FALOUTNR(((K - KMS + 1))) - FALOUTNR(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDC  = ((FALOUTC(((K - KMS + 1))) - FALOUTC(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDNC = ((FALOUTNC(((K - KMS + 1))) - FALOUTNC(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDG  = ((FALOUTG(((K - KMS + 1))) - FALOUTG(((K - KMS)))) / DZQ(((K - KMS))));
      FALTNDNG = ((FALOUTNG(((K - KMS + 1))) - FALOUTNG(((K - KMS)))) / DZQ(((K - KMS))));
      QRSTEN(((K - KMS)))  = (QRSTEN(((K - KMS))) + ((FALTNDR / NSTEP) / RHO(((K - KMS)))));
      QISTEN(((K - KMS)))  = (QISTEN(((K - KMS))) + ((FALTNDI / NSTEP) / RHO(((K - KMS)))));
      NI3DTEN(((K - KMS))) = (NI3DTEN(((K - KMS))) + ((FALTNDNI / NSTEP) / RHO(((K - KMS)))));
      QNISTEN(((K - KMS))) = (QNISTEN(((K - KMS))) + ((FALTNDS / NSTEP) / RHO(((K - KMS)))));
      NS3DTEN(((K - KMS))) = (NS3DTEN(((K - KMS))) + ((FALTNDNS / NSTEP) / RHO(((K - KMS)))));
      NR3DTEN(((K - KMS))) = (NR3DTEN(((K - KMS))) + ((FALTNDNR / NSTEP) / RHO(((K - KMS)))));
      QCSTEN(((K - KMS)))  = (QCSTEN(((K - KMS))) + ((FALTNDC / NSTEP) / RHO(((K - KMS)))));
      NC3DTEN(((K - KMS))) = (NC3DTEN(((K - KMS))) + ((FALTNDNC / NSTEP) / RHO(((K - KMS)))));
      QGSTEN(((K - KMS)))  = (QGSTEN(((K - KMS))) + ((FALTNDG / NSTEP) / RHO(((K - KMS)))));
      NG3DTEN(((K - KMS))) = (NG3DTEN(((K - KMS))) + ((FALTNDNG / NSTEP) / RHO(((K - KMS)))));
      NISTEN(((K - KMS)))  = (NISTEN(((K - KMS))) + ((FALTNDNI / NSTEP) / RHO(((K - KMS)))));
      NSSTEN(((K - KMS)))  = (NSSTEN(((K - KMS))) + ((FALTNDNS / NSTEP) / RHO(((K - KMS)))));
      NRSTEN(((K - KMS)))  = (NRSTEN(((K - KMS))) + ((FALTNDNR / NSTEP) / RHO(((K - KMS)))));
      NCSTEN(((K - KMS)))  = (NCSTEN(((K - KMS))) + ((FALTNDNC / NSTEP) / RHO(((K - KMS)))));
      NGSTEN(((K - KMS)))  = (NGSTEN(((K - KMS))) + ((FALTNDNG / NSTEP) / RHO(((K - KMS)))));
      DUMR(((K - KMS)))    = (DUMR(((K - KMS))) + ((FALTNDR * DT) / NSTEP));
      DUMI(((K - KMS)))    = (DUMI(((K - KMS))) + ((FALTNDI * DT) / NSTEP));
      DUMFNI(((K - KMS)))  = (DUMFNI(((K - KMS))) + ((FALTNDNI * DT) / NSTEP));
      DUMQS(((K - KMS)))   = (DUMQS(((K - KMS))) + ((FALTNDS * DT) / NSTEP));
      DUMFNS(((K - KMS)))  = (DUMFNS(((K - KMS))) + ((FALTNDNS * DT) / NSTEP));
      DUMFNR(((K - KMS)))  = (DUMFNR(((K - KMS))) + ((FALTNDNR * DT) / NSTEP));
      DUMC(((K - KMS)))    = (DUMC(((K - KMS))) + ((FALTNDC * DT) / NSTEP));
      DUMFNC(((K - KMS)))  = (DUMFNC(((K - KMS))) + ((FALTNDNC * DT) / NSTEP));
      DUMG(((K - KMS)))    = (DUMG(((K - KMS))) + ((FALTNDG * DT) / NSTEP));
      DUMFNG(((K - KMS)))  = (DUMFNG(((K - KMS))) + ((FALTNDNG * DT) / NSTEP));
    }
    PRECRT = (PRECRT +
              ((((((FALOUTR(((-KMS + KTS))) + FALOUTC(((-KMS + KTS)))) + FALOUTS(((-KMS + KTS)))) +
                  FALOUTI(((-KMS + KTS)))) +
                 FALOUTG(((-KMS + KTS)))) *
                DT) /
               NSTEP));
    SNOWRT =
        (SNOWRT +
         ((((FALOUTS(((-KMS + KTS))) + FALOUTI(((-KMS + KTS)))) + FALOUTG(((-KMS + KTS)))) * DT) /
          NSTEP));
  }
  for (K = KTS; K < KTE + 1; K += 1) {
    QR3DTEN(((K - KMS)))  = (QR3DTEN(((K - KMS))) + QRSTEN(((K - KMS))));
    QI3DTEN(((K - KMS)))  = (QI3DTEN(((K - KMS))) + QISTEN(((K - KMS))));
    QC3DTEN(((K - KMS)))  = (QC3DTEN(((K - KMS))) + QCSTEN(((K - KMS))));
    QG3DTEN(((K - KMS)))  = (QG3DTEN(((K - KMS))) + QGSTEN(((K - KMS))));
    QNI3DTEN(((K - KMS))) = (QNI3DTEN(((K - KMS))) + QNISTEN(((K - KMS))));
    QC3D(((K - KMS)))     = (QC3D(((K - KMS))) + (QC3DTEN(((K - KMS))) * DT));
    QI3D(((K - KMS)))     = (QI3D(((K - KMS))) + (QI3DTEN(((K - KMS))) * DT));
    QNI3D(((K - KMS)))    = (QNI3D(((K - KMS))) + (QNI3DTEN(((K - KMS))) * DT));
    QR3D(((K - KMS)))     = (QR3D(((K - KMS))) + (QR3DTEN(((K - KMS))) * DT));
    NC3D(((K - KMS)))     = (NC3D(((K - KMS))) + (NC3DTEN(((K - KMS))) * DT));
    NI3D(((K - KMS)))     = (NI3D(((K - KMS))) + (NI3DTEN(((K - KMS))) * DT));
    NS3D(((K - KMS)))     = (NS3D(((K - KMS))) + (NS3DTEN(((K - KMS))) * DT));
    NR3D(((K - KMS)))     = (NR3D(((K - KMS))) + (NR3DTEN(((K - KMS))) * DT));
    if ((morrison.IGRAUP == 0)) {
      QG3D(((K - KMS))) = (QG3D(((K - KMS))) + (QG3DTEN(((K - KMS))) * DT));
      NG3D(((K - KMS))) = (NG3D(((K - KMS))) + (NG3DTEN(((K - KMS))) * DT));
    }
    T3D(((K - KMS)))  = (T3D(((K - KMS))) + (T3DTEN(((K - KMS))) * DT));
    QV3D(((K - KMS))) = (QV3D(((K - KMS))) + (QV3DTEN(((K - KMS))) * DT));
    EVS(((K - KMS))) =
        (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(T3D(((K - KMS))), 0, morrison))
             ? ((0.99f * PRES(((K - KMS)))))
             : (morrison_engine::POLYSVP(T3D(((K - KMS))), 0, morrison)));
    EIS(((K - KMS))) =
        (((0.99f * PRES(((K - KMS))))) < (morrison_engine::POLYSVP(T3D(((K - KMS))), 1, morrison))
             ? ((0.99f * PRES(((K - KMS)))))
             : (morrison_engine::POLYSVP(T3D(((K - KMS))), 1, morrison)));
    if ((EIS(((K - KMS))) > EVS(((K - KMS))))) {
      EIS(((K - KMS))) = EVS(((K - KMS)));
    }
    QVS(((K - KMS)))    = ((ep * EVS(((K - KMS)))) / (PRES(((K - KMS))) - EVS(((K - KMS)))));
    QVI(((K - KMS)))    = ((ep * EIS(((K - KMS)))) / (PRES(((K - KMS))) - EIS(((K - KMS)))));
    QVQVS(((K - KMS)))  = (QV3D(((K - KMS))) / QVS(((K - KMS))));
    QVQVSI(((K - KMS))) = (QV3D(((K - KMS))) / QVI(((K - KMS))));
    if ((QVQVS(((K - KMS))) < 0.9f)) {
      if ((QR3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QR3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QR3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
        QR3D(((K - KMS)))    = 0.f;
      }
      if ((QC3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QC3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QC3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
        QC_INST(((K - KMS))) = (QC_INST(((K - KMS))) - (QC3D(((K - KMS))) / DT));
        QC3D(((K - KMS)))    = 0.f;
      }
    }
    if ((QVQVSI(((K - KMS))) < 0.9f)) {
      if ((QI3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QI3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) - (QI3D(((K - KMS))) / DT));
        QI3D(((K - KMS)))    = 0.f;
      }
      if ((QNI3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QNI3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QNI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) - (QNI3D(((K - KMS))) / DT));
        QNI3D(((K - KMS)))   = 0.f;
      }
      if ((QG3D(((K - KMS))) < 1.E-8f)) {
        QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QG3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) - ((QG3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
        QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) - (QG3D(((K - KMS))) / DT));
        QG3D(((K - KMS)))    = 0.f;
      }
    }
    if ((QC3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QC3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QC3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
      NC_INST(((K - KMS))) = (NC_INST(((K - KMS))) - (NC3D(((K - KMS))) / DT));
      QC_INST(((K - KMS))) = (QC_INST(((K - KMS))) - (QC3D(((K - KMS))) / DT));
      QC3D(((K - KMS)))    = 0.f;
      NC3D(((K - KMS)))    = 0.f;
      EFFC(((K - KMS)))    = 0.f;
    }
    if ((QR3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QR3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QR3D(((K - KMS))) * XXLV(((K - KMS)))) / CPM(((K - KMS)))));
      NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) - (NR3D(((K - KMS))) / DT));
      QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
      QR3D(((K - KMS)))    = 0.f;
      NR3D(((K - KMS)))    = 0.f;
      EFFR(((K - KMS)))    = 0.f;
    }
    if ((QI3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QI3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NI_INST(((K - KMS))) = (NI_INST(((K - KMS))) - (NI3D(((K - KMS))) / DT));
      QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) - (QI3D(((K - KMS))) / DT));
      QI3D(((K - KMS)))    = 0.f;
      NI3D(((K - KMS)))    = 0.f;
      EFFI(((K - KMS)))    = 0.f;
    }
    if ((QNI3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QNI3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QNI3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NS_INST(((K - KMS))) = (NS_INST(((K - KMS))) - (NS3D(((K - KMS))) / DT));
      QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) - (QNI3D(((K - KMS))) / DT));
      QNI3D(((K - KMS)))   = 0.f;
      NS3D(((K - KMS)))    = 0.f;
      EFFS(((K - KMS)))    = 0.f;
    }
    if ((QG3D(((K - KMS))) < morrison.QSMALL)) {
      QV3D(((K - KMS))) = (QV3D(((K - KMS))) + QG3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QG3D(((K - KMS))) * XXLS(((K - KMS)))) / CPM(((K - KMS)))));
      NG_INST(((K - KMS))) = (NG_INST(((K - KMS))) - (NG3D(((K - KMS))) / DT));
      QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) - (QG3D(((K - KMS))) / DT));
      QG3D(((K - KMS)))    = 0.f;
      NG3D(((K - KMS)))    = 0.f;
      EFFG(((K - KMS)))    = 0.f;
    }
    if ((((((QC3D(((K - KMS))) < morrison.QSMALL) && (QI3D(((K - KMS))) < morrison.QSMALL)) &&
           (QNI3D(((K - KMS))) < morrison.QSMALL)) &&
          (QR3D(((K - KMS))) < morrison.QSMALL)) &&
         (QG3D(((K - KMS))) < morrison.QSMALL))) {
      goto label_500;
    }
    if (((QI3D(((K - KMS))) >= morrison.QSMALL) && (T3D(((K - KMS))) >= morrison.TMELT))) {
      QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) + (QI3D(((K - KMS))) / DT));
      QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) - (QI3D(((K - KMS))) / DT));
      NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) + (NI3D(((K - KMS))) / DT));
      NI_INST(((K - KMS))) = (NI_INST(((K - KMS))) - (NI3D(((K - KMS))) / DT));
      QR3D(((K - KMS)))    = (QR3D(((K - KMS))) + QI3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) - ((QI3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
      QI3D(((K - KMS))) = 0.f;
      NR3D(((K - KMS))) = (NR3D(((K - KMS))) + NI3D(((K - KMS))));
      NI3D(((K - KMS))) = 0.f;
    }
    if ((morrison.ILIQ == 1)) {
      goto label_778;
    }
    if (((T3D(((K - KMS))) <= 233.15f) && (QC3D(((K - KMS))) >= morrison.QSMALL))) {
      QI_INST(((K - KMS))) = (QI_INST(((K - KMS))) + (QC3D(((K - KMS))) / DT));
      QC_INST(((K - KMS))) = (QC_INST(((K - KMS))) - (QC3D(((K - KMS))) / DT));
      NI_INST(((K - KMS))) = (NI_INST(((K - KMS))) + (NC3D(((K - KMS))) / DT));
      NC_INST(((K - KMS))) = (NC_INST(((K - KMS))) - (NC3D(((K - KMS))) / DT));
      QI3D(((K - KMS)))    = (QI3D(((K - KMS))) + QC3D(((K - KMS))));
      T3D(((K - KMS))) =
          (T3D(((K - KMS))) + ((QC3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
      QC3D(((K - KMS))) = 0.f;
      NI3D(((K - KMS))) = (NI3D(((K - KMS))) + NC3D(((K - KMS))));
      NC3D(((K - KMS))) = 0.f;
    }
    if ((morrison.IGRAUP == 0)) {
      if (((T3D(((K - KMS))) <= 233.15f) && (QR3D(((K - KMS))) >= morrison.QSMALL))) {
        QG_INST(((K - KMS))) = (QG_INST(((K - KMS))) + (QR3D(((K - KMS))) / DT));
        QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
        NG_INST(((K - KMS))) = (NG_INST(((K - KMS))) + (NR3D(((K - KMS))) / DT));
        NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) - (NR3D(((K - KMS))) / DT));
        QG3D(((K - KMS)))    = (QG3D(((K - KMS))) + QR3D(((K - KMS))));
        T3D(((K - KMS))) =
            (T3D(((K - KMS))) + ((QR3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
        QR3D(((K - KMS))) = 0.f;
        NG3D(((K - KMS))) = (NG3D(((K - KMS))) + NR3D(((K - KMS))));
        NR3D(((K - KMS))) = 0.f;
      }
    } else {
      if ((morrison.IGRAUP == 1)) {
        if (((T3D(((K - KMS))) <= 233.15f) && (QR3D(((K - KMS))) >= morrison.QSMALL))) {
          QS_INST(((K - KMS))) = (QS_INST(((K - KMS))) + (QR3D(((K - KMS))) / DT));
          QR_INST(((K - KMS))) = (QR_INST(((K - KMS))) - (QR3D(((K - KMS))) / DT));
          NS_INST(((K - KMS))) = (NS_INST(((K - KMS))) + (NR3D(((K - KMS))) / DT));
          NR_INST(((K - KMS))) = (NR_INST(((K - KMS))) - (NR3D(((K - KMS))) / DT));
          QNI3D(((K - KMS)))   = (QNI3D(((K - KMS))) + QR3D(((K - KMS))));
          T3D(((K - KMS))) =
              (T3D(((K - KMS))) + ((QR3D(((K - KMS))) * XLF(((K - KMS)))) / CPM(((K - KMS)))));
          QR3D(((K - KMS))) = 0.f;
          NS3D(((K - KMS))) = (NS3D(((K - KMS))) + NR3D(((K - KMS))));
          NR3D(((K - KMS))) = 0.f;
        }
      }
    }
  label_778:;
    if ((NI3D(((K - KMS))) < 0.f)) {
      NEGFIX_NI(((K - KMS))) = (NEGFIX_NI(((K - KMS))) + (NI3D(((K - KMS))) / DT));
    }
    if ((NS3D(((K - KMS))) < 0.f)) {
      NEGFIX_NS(((K - KMS))) = (NEGFIX_NS(((K - KMS))) + (NS3D(((K - KMS))) / DT));
    }
    if ((NC3D(((K - KMS))) < 0.f)) {
      NEGFIX_NC(((K - KMS))) = (NEGFIX_NC(((K - KMS))) + (NC3D(((K - KMS))) / DT));
    }
    if ((NR3D(((K - KMS))) < 0.f)) {
      NEGFIX_NR(((K - KMS))) = (NEGFIX_NR(((K - KMS))) + (NR3D(((K - KMS))) / DT));
    }
    if ((NG3D(((K - KMS))) < 0.f)) {
      NEGFIX_NG(((K - KMS))) = (NEGFIX_NG(((K - KMS))) + (NG3D(((K - KMS))) / DT));
    }
    NI3D(((K - KMS))) = ((0.f) > (NI3D(((K - KMS)))) ? (0.f) : (NI3D(((K - KMS)))));
    NS3D(((K - KMS))) = ((0.f) > (NS3D(((K - KMS)))) ? (0.f) : (NS3D(((K - KMS)))));
    NC3D(((K - KMS))) = ((0.f) > (NC3D(((K - KMS)))) ? (0.f) : (NC3D(((K - KMS)))));
    NR3D(((K - KMS))) = ((0.f) > (NR3D(((K - KMS)))) ? (0.f) : (NR3D(((K - KMS)))));
    NG3D(((K - KMS))) = ((0.f) > (NG3D(((K - KMS)))) ? (0.f) : (NG3D(((K - KMS)))));
    if ((QI3D(((K - KMS))) >= morrison.QSMALL)) {
      LAMI(((K - KMS))) = Kokkos::pow(((morrison.CONS12 * NI3D(((K - KMS)))) / QI3D(((K - KMS)))),
                                      (1.f / morrison.DI));
      TMPNUM            = 0.f;
      if ((LAMI(((K - KMS))) < morrison.LAMMINI)) {
        LAMI(((K - KMS))) = morrison.LAMMINI;
        N0I(((K - KMS))) =
            ((Kokkos::pow(LAMI(((K - KMS))), (morrison.DI + 1.f)) * QI3D(((K - KMS)))) /
             morrison.CONS12);
        TMPNUM                  = NI3D(((K - KMS)));
        NI3D(((K - KMS)))       = (N0I(((K - KMS))) / LAMI(((K - KMS))));
        SIZEFIX_NI(((K - KMS))) = (SIZEFIX_NI(((K - KMS))) + ((NI3D(((K - KMS))) - TMPNUM) / DT));
      } else {
        if ((LAMI(((K - KMS))) > morrison.LAMMAXI)) {
          LAMI(((K - KMS))) = morrison.LAMMAXI;
          N0I(((K - KMS))) =
              ((Kokkos::pow(LAMI(((K - KMS))), (morrison.DI + 1.f)) * QI3D(((K - KMS)))) /
               morrison.CONS12);
          TMPNUM                  = NI3D(((K - KMS)));
          NI3D(((K - KMS)))       = (N0I(((K - KMS))) / LAMI(((K - KMS))));
          SIZEFIX_NI(((K - KMS))) = (SIZEFIX_NI(((K - KMS))) + ((NI3D(((K - KMS))) - TMPNUM) / DT));
        }
      }
    }
    if ((QR3D(((K - KMS))) >= morrison.QSMALL)) {
      LAMR(((K - KMS))) = Kokkos::pow(
          (((morrison.PI * morrison.RHOW) * NR3D(((K - KMS)))) / QR3D(((K - KMS)))), (1.f / 3.f));
      TMPNUM = 0.f;
      if ((LAMR(((K - KMS))) < morrison.LAMMINR)) {
        TMPNUM            = NR3D(((K - KMS)));
        LAMR(((K - KMS))) = morrison.LAMMINR;
        // BFB-030: Native integer fourth power keeps both squares in default REAL.
        const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
        N0RR(((K - KMS))) =
            (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
        NR3D(((K - KMS)))       = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
        SIZEFIX_NR(((K - KMS))) = (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
      } else {
        if ((LAMR(((K - KMS))) > morrison.LAMMAXR)) {
          TMPNUM            = NR3D(((K - KMS)));
          LAMR(((K - KMS))) = morrison.LAMMAXR;
          // BFB-030: Native integer fourth power keeps both squares in default REAL.
          const float LAMR_sqd = LAMR(((K - KMS))) * LAMR(((K - KMS)));
          N0RR(((K - KMS))) =
              (((LAMR_sqd * LAMR_sqd) * QR3D(((K - KMS)))) / (morrison.PI * morrison.RHOW));
          NR3D(((K - KMS)))       = (N0RR(((K - KMS))) / LAMR(((K - KMS))));
          SIZEFIX_NR(((K - KMS))) = (SIZEFIX_NR(((K - KMS))) + ((NR3D(((K - KMS))) - TMPNUM) / DT));
        }
      }
    }
    if ((QC3D(((K - KMS))) >= morrison.QSMALL)) {
      if (morrison.dofix_pgam) {
        PGAM(((K - KMS))) = morrison.pgam_fixed;
      } else {
        PGAM(((K - KMS))) =
            ((0.0005714f * ((NC3D(((K - KMS))) / 1.E6f) * RHO(((K - KMS))))) + 0.2714f);
        PGAM(((K - KMS))) = ((1.f / (PGAM(((K - KMS))) * PGAM(((K - KMS))))) - 1.f);
        PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) > (2.f) ? (PGAM(((K - KMS)))) : (2.f));
        PGAM(((K - KMS))) = ((PGAM(((K - KMS)))) < (10.f) ? (PGAM(((K - KMS)))) : (10.f));
      }
      LAMC(((K - KMS))) =
          Kokkos::pow((((morrison.CONS26 * NC3D(((K - KMS)))) *
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)) /
                       (QC3D(((K - KMS))) *
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))),
                      (1.f / 3.f));
      LAMMIN = ((PGAM(((K - KMS))) + 1.f) / 60.E-6f);
      LAMMAX = ((PGAM(((K - KMS))) + 1.f) / 1.E-6f);
      TMPNUM = 0.f;
      if ((LAMC(((K - KMS))) < LAMMIN)) {
        LAMC(((K - KMS))) = LAMMIN;
        TMPNUM            = NC3D(((K - KMS)));
        NC3D(((K - KMS))) =
            (Kokkos::exp(
                 ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                   Kokkos::log(
                       morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))) -
                  Kokkos::log(
                      morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)))) /
             morrison.CONS26);
        SIZEFIX_NC(((K - KMS))) = (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
      } else {
        if ((LAMC(((K - KMS))) > LAMMAX)) {
          LAMC(((K - KMS))) = LAMMAX;
          TMPNUM            = NC3D(((K - KMS)));
          NC3D(((K - KMS))) =
              (Kokkos::exp(
                   ((((3.f * Kokkos::log(LAMC(((K - KMS))))) + Kokkos::log(QC3D(((K - KMS))))) +
                     Kokkos::log(
                         morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 1.f), morrison))) -
                    Kokkos::log(
                        morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison)))) /
               morrison.CONS26);
          SIZEFIX_NC(((K - KMS))) = (SIZEFIX_NC(((K - KMS))) + ((NC3D(((K - KMS))) - TMPNUM) / DT));
        }
      }
    }
    if ((QNI3D(((K - KMS))) >= morrison.QSMALL)) {
      LAMS(((K - KMS))) = Kokkos::pow(((morrison.CONS1 * NS3D(((K - KMS)))) / QNI3D(((K - KMS)))),
                                      (1.f / morrison.DS));
      TMPNUM            = 0.f;
      if ((LAMS(((K - KMS))) < morrison.LAMMINS)) {
        LAMS(((K - KMS))) = morrison.LAMMINS;
        N0S(((K - KMS))) =
            ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
             morrison.CONS1);
        TMPNUM                  = NS3D(((K - KMS)));
        NS3D(((K - KMS)))       = (N0S(((K - KMS))) / LAMS(((K - KMS))));
        SIZEFIX_NS(((K - KMS))) = (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
      } else {
        if ((LAMS(((K - KMS))) > morrison.LAMMAXS)) {
          LAMS(((K - KMS))) = morrison.LAMMAXS;
          N0S(((K - KMS))) =
              ((Kokkos::pow(LAMS(((K - KMS))), (morrison.DS + 1.f)) * QNI3D(((K - KMS)))) /
               morrison.CONS1);
          TMPNUM                  = NS3D(((K - KMS)));
          NS3D(((K - KMS)))       = (N0S(((K - KMS))) / LAMS(((K - KMS))));
          SIZEFIX_NS(((K - KMS))) = (SIZEFIX_NS(((K - KMS))) + ((NS3D(((K - KMS))) - TMPNUM) / DT));
        }
      }
    }
    if ((QG3D(((K - KMS))) >= morrison.QSMALL)) {
      LAMG(((K - KMS))) = Kokkos::pow(((morrison.CONS2 * NG3D(((K - KMS)))) / QG3D(((K - KMS)))),
                                      (1.f / morrison.DG));
      TMPNUM            = 0.f;
      if ((LAMG(((K - KMS))) < morrison.LAMMING)) {
        LAMG(((K - KMS))) = morrison.LAMMING;
        N0G(((K - KMS))) =
            ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
             morrison.CONS2);
        TMPNUM                  = NG3D(((K - KMS)));
        NG3D(((K - KMS)))       = (N0G(((K - KMS))) / LAMG(((K - KMS))));
        SIZEFIX_NG(((K - KMS))) = (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
      } else {
        if ((LAMG(((K - KMS))) > morrison.LAMMAXG)) {
          LAMG(((K - KMS))) = morrison.LAMMAXG;
          N0G(((K - KMS))) =
              ((Kokkos::pow(LAMG(((K - KMS))), (morrison.DG + 1.f)) * QG3D(((K - KMS)))) /
               morrison.CONS2);
          TMPNUM                  = NG3D(((K - KMS)));
          NG3D(((K - KMS)))       = (N0G(((K - KMS))) / LAMG(((K - KMS))));
          SIZEFIX_NG(((K - KMS))) = (SIZEFIX_NG(((K - KMS))) + ((NG3D(((K - KMS))) - TMPNUM) / DT));
        }
      }
    }
  label_500:;
    if ((CF3D(((K - KMS))) > morrison.cloud_frac_thresh)) {
      TMPQSMALL = (morrison.QSMALL / CF3D(((K - KMS))));
    } else {
      TMPQSMALL = morrison.QSMALL;
    }
    if ((QI3D(((K - KMS))) >= TMPQSMALL)) {
      EFFI(((K - KMS))) = (((3.f / LAMI(((K - KMS)))) / 2.f) * 1.E6f);
    } else {
      EFFI(((K - KMS))) = 25.f;
    }
    if ((QNI3D(((K - KMS))) >= TMPQSMALL)) {
      EFFS(((K - KMS))) = (((3.f / LAMS(((K - KMS)))) / 2.f) * 1.E6f);
    } else {
      EFFS(((K - KMS))) = 25.f;
    }
    if ((QR3D(((K - KMS))) >= TMPQSMALL)) {
      EFFR(((K - KMS))) = (((3.f / LAMR(((K - KMS)))) / 2.f) * 1.E6f);
    } else {
      EFFR(((K - KMS))) = 25.f;
    }
    if ((QC3D(((K - KMS))) >= TMPQSMALL)) {
      EFFC(((K - KMS))) =
          ((((morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 4.f), morrison) /
              morrison_engine::GAMMA(workspace, (PGAM(((K - KMS))) + 3.f), morrison)) /
             LAMC(((K - KMS)))) /
            2.f) *
           1.E6f);
    } else {
      EFFC(((K - KMS))) = 25.f;
    }
    if ((QG3D(((K - KMS))) >= TMPQSMALL)) {
      EFFG(((K - KMS))) = (((3.f / LAMG(((K - KMS)))) / 2.f) * 1.E6f);
    } else {
      EFFG(((K - KMS))) = 25.f;
    }
    if ((NI3D(((K - KMS))) > (10.E6f / RHO(((K - KMS)))))) {
      NIM_MORR_CL(((K - KMS))) = ((NI3D(((K - KMS))) - (10.E6f / RHO(((K - KMS))))) / DT);
    }
    NI3D(((K - KMS))) =
        ((NI3D(((K - KMS)))) < ((10.E6f / RHO(((K - KMS))))) ? (NI3D(((K - KMS))))
                                                             : ((10.E6f / RHO(((K - KMS))))));
    if (((morrison.INUM == 0) && (morrison.IACT == 2))) {
      NC3D(((K - KMS))) =
          ((NC3D(((K - KMS)))) < (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS)))))
               ? (NC3D(((K - KMS))))
               : (((morrison.NANEW1 + morrison.NANEW2) / RHO(((K - KMS))))));
    }
    if ((morrison.INUM == 1)) {
      NC3D(((K - KMS))) = ((morrison.NDCNST * 1.E6f) / RHO(((K - KMS))));
    }
    morrison_engine::POSITIVE_QV_ADJ(QV3D(((K - KMS))), QC3D(((K - KMS))), QR3D(((K - KMS))),
                                     QI3D(((K - KMS))), QNI3D(((K - KMS))), QG3D(((K - KMS))),
                                     T3D(((K - KMS))), morrison);
  }
label_400:;
  return;
}
inline float POLYSVP(float T, int TYPE, morrison_microphysics_state& morrison) {
  // Execute the Native POLYSVP body and update its declared outputs.
  float DUM;
  float a0i;
  float a1i;
  float a2i;
  float a3i;
  float a4i;
  float a5i;
  float a6i;
  float a7i;
  float a8i;
  float a0;
  float a1;
  float a2;
  float a3;
  float a4;
  float a5;
  float a6;
  float a7;
  float a8;
  float dt;
  float POLYSVP;
  a0i = 6.11147274f;
  a1i = 0.503160820f;
  a2i = 0.188439774e-1f;
  a3i = 0.420895665e-3f;
  a4i = 0.615021634e-5f;
  a5i = 0.602588177e-7f;
  a6i = 0.385852041e-9f;
  a7i = 0.146898966e-11f;
  a8i = 0.252751365e-14f;
  a0  = 6.11239921f;
  a1  = 0.443987641f;
  a2  = 0.142986287e-1f;
  a3  = 0.264847430e-3f;
  a4  = 0.302950461e-5f;
  a5  = 0.206739458e-7f;
  a6  = 0.640689451e-10f;
  a7  = (-0.952447341e-13f);
  a8  = (-0.976195544e-15f);
  if ((TYPE == 1)) {
    dt = (((-80.f)) > ((T - 273.16f)) ? ((-80.f)) : ((T - 273.16f)));
    POLYSVP =
        (a0i +
         (dt *
          (a1i +
           (dt *
            (a2i +
             (dt *
              (a3i + (dt * (a4i + (dt * (a5i + (dt * (a6i + (dt * (a7i + (a8i * dt))))))))))))))));
    POLYSVP = (POLYSVP * 100.f);
  }
  if ((TYPE == 0)) {
    dt = (((-80.f)) > ((T - 273.16f)) ? ((-80.f)) : ((T - 273.16f)));
    POLYSVP =
        (a0 +
         (dt *
          (a1 +
           (dt *
            (a2 +
             (dt * (a3 + (dt * (a4 + (dt * (a5 + (dt * (a6 + (dt * (a7 + (a8 * dt))))))))))))))));
    POLYSVP = (POLYSVP * 100.f);
  }
  return POLYSVP;
}
inline float GAMMA(ScratchWorkspace& workspace, float X, morrison_microphysics_state& morrison) {
  // Execute the Native GAMMA body and update its declared outputs.
  auto workspace_scope = workspace.scope();
  int I;
  int N;
  bool PARITY;
  float CONV;
  float EPS;
  float FACT;
  float HALF;
  float ONE;
  float RES;
  float SUM;
  float TWELVE;
  float TWO;
  float XBIG;
  float XDEN;
  float XINF;
  float XMININ;
  float XNUM;
  float Y;
  float Y1;
  float YSQ;
  float Z;
  float ZERO;
  HostFloat1D C = workspace_scope.view<HostFloat1D>("GAMMA::C", 7);
  HostFloat1D P = workspace_scope.view<HostFloat1D>("GAMMA::P", 8);
  HostFloat1D Q = workspace_scope.view<HostFloat1D>("GAMMA::Q", 8);
  float GAMMA;
  ONE          = 1.0E0f;
  HALF         = 0.5E0f;
  TWELVE       = 12.0E0f;
  TWO          = 2.0E0f;
  ZERO         = 0.0E0f;
  XBIG         = 35.040E0f;
  XMININ       = 1.18E-38f;
  EPS          = 1.19E-7f;
  XINF         = 3.4E38f;
  P(0)         = (-1.71618513886549492533811E+0f);
  P(1)         = 2.47656508055759199108314E+1f;
  P(2)         = (-3.79804256470945635097577E+2f);
  P(3)         = 6.29331155312818442661052E+2f;
  P(4)         = 8.66966202790413211295064E+2f;
  P(5)         = (-3.14512729688483675254357E+4f);
  P(6)         = (-3.61444134186911729807069E+4f);
  P(7)         = 6.64561438202405440627855E+4f;
  Q(0)         = (-3.08402300119738975254353E+1f);
  Q(1)         = 3.15350626979604161529144E+2f;
  Q(2)         = (-1.01515636749021914166146E+3f);
  Q(3)         = (-3.10777167157231109440444E+3f);
  Q(4)         = 2.25381184209801510330112E+4f;
  Q(5)         = 4.75584627752788110767815E+3f;
  Q(6)         = (-1.34659959864969306392456E+5f);
  Q(7)         = (-1.15132259675553483497211E+5f);
  C(0)         = (-1.910444077728E-03f);
  C(1)         = 8.4171387781295E-04f;
  C(2)         = (-5.952379913043012E-04f);
  C(3)         = 7.93650793500350248E-04f;
  C(4)         = (-2.777777777777681622553E-03f);
  C(5)         = 8.333333333333333331554247E-02f;
  C(6)         = 5.7083835261E-03f;
  // Native statement function CONV(I)=REAL(I) is expanded at its call sites.
  PARITY = false;
  FACT   = ONE;
  N      = 0;
  Y      = X;
  if ((Y <= ZERO)) {
    Y   = (-X);
    Y1  = Kokkos::trunc(Y);
    RES = (Y - Y1);
    if ((RES != ZERO)) {
      if ((Y1 != (Kokkos::trunc((Y1 * HALF)) * TWO))) {
        PARITY = true;
      }
      FACT = ((-morrison.PI) / Kokkos::sin((morrison.PI * RES)));
      Y    = (Y + ONE);
    } else {
      RES = XINF;
      goto label_900;
    }
  }
  if ((Y < EPS)) {
    if ((Y >= XMININ)) {
      RES = (ONE / Y);
    } else {
      RES = XINF;
      goto label_900;
    }
  } else {
    if ((Y < TWELVE)) {
      Y1 = Y;
      if ((Y < ONE)) {
        Z = Y;
        Y = (Y + ONE);
      } else {
        N = (static_cast<int>(Y) - 1);
        Y = (Y - static_cast<float>(N));
        Z = (Y - ONE);
      }
      XNUM = ZERO;
      XDEN = ONE;
      for (I = 0; I < 8; I += 1) {
        XNUM = ((XNUM + P(I)) * Z);
        XDEN = ((XDEN * Z) + Q(I));
      }
      RES = ((XNUM / XDEN) + ONE);
      if ((Y1 < Y)) {
        RES = (RES / Y1);
      } else {
        if ((Y1 > Y)) {
          for (I = 0; I < N; I += 1) {
            RES = (RES * Y);
            Y   = (Y + ONE);
          }
        }
      }
    } else {
      if ((Y <= XBIG)) {
        YSQ = (Y * Y);
        SUM = C(6);
        for (I = 0; I < 6; I += 1) {
          SUM = ((SUM / YSQ) + C(I));
        }
        SUM = (((SUM / Y) - Y) + morrison.SQRTPI);
        SUM = (SUM + ((Y - HALF) * Kokkos::log(Y)));
        RES = Kokkos::exp(SUM);
      } else {
        RES = XINF;
        goto label_900;
      }
    }
  }
  if (PARITY) {
    RES = (-RES);
  }
  if ((FACT != ONE)) {
    RES = (FACT / RES);
  }
label_900:
  GAMMA = RES;
  return GAMMA;
}
inline float DERF1(ScratchWorkspace& workspace, float X, morrison_microphysics_state& morrison) {
  // Execute the Native DERF1 body and update its declared outputs.
  auto workspace_scope = workspace.scope();
  HostFloat1D A        = workspace_scope.view<HostFloat1D>("DERF1::A", 65);
  HostFloat1D B        = workspace_scope.view<HostFloat1D>("DERF1::B", 65);
  float W;
  float T;
  float Y;
  int K;
  int I;
  float DERF1;
  A(0)          = 0.00000000005958930743E0f;
  A(1)          = (-0.00000000113739022964E0f);
  A(2)          = 0.00000001466005199839E0f;
  A(3)          = (-0.00000016350354461960E0f);
  A(4)          = 0.00000164610044809620E0f;
  A(5)          = (-0.00001492559551950604E0f);
  A(6)          = 0.00012055331122299265E0f;
  A(7)          = (-0.00085483269811296660E0f);
  A(8)          = 0.00522397762482322257E0f;
  A(9)          = (-0.02686617064507733420E0f);
  A(10)         = 0.11283791670954881569E0f;
  A(11)         = (-0.37612638903183748117E0f);
  A(12)         = 1.12837916709551257377E0f;
  A(13)         = 0.00000000002372510631E0f;
  A(14)         = (-0.00000000045493253732E0f);
  A(15)         = 0.00000000590362766598E0f;
  A(16)         = (-0.00000006642090827576E0f);
  A(17)         = 0.00000067595634268133E0f;
  A(18)         = (-0.00000621188515924000E0f);
  A(19)         = 0.00005103883009709690E0f;
  A(20)         = (-0.00037015410692956173E0f);
  A(21)         = 0.00233307631218880978E0f;
  A(22)         = (-0.01254988477182192210E0f);
  A(23)         = 0.05657061146827041994E0f;
  A(24)         = (-0.21379664776456006580E0f);
  A(25)         = 0.84270079294971486929E0f;
  A(26)         = 0.00000000000949905026E0f;
  A(27)         = (-0.00000000018310229805E0f);
  A(28)         = 0.00000000239463074000E0f;
  A(29)         = (-0.00000002721444369609E0f);
  A(30)         = 0.00000028045522331686E0f;
  A(31)         = (-0.00000261830022482897E0f);
  A(32)         = 0.00002195455056768781E0f;
  A(33)         = (-0.00016358986921372656E0f);
  A(34)         = 0.00107052153564110318E0f;
  A(35)         = (-0.00608284718113590151E0f);
  A(36)         = 0.02986978465246258244E0f;
  A(37)         = (-0.13055593046562267625E0f);
  A(38)         = 0.67493323603965504676E0f;
  A(39)         = 0.00000000000382722073E0f;
  A(40)         = (-0.00000000007421598602E0f);
  A(41)         = 0.00000000097930574080E0f;
  A(42)         = (-0.00000001126008898854E0f);
  A(43)         = 0.00000011775134830784E0f;
  A(44)         = (-0.00000111992758382650E0f);
  A(45)         = 0.00000962023443095201E0f;
  A(46)         = (-0.00007404402135070773E0f);
  A(47)         = 0.00050689993654144881E0f;
  A(48)         = (-0.00307553051439272889E0f);
  A(49)         = 0.01668977892553165586E0f;
  A(50)         = (-0.08548534594781312114E0f);
  A(51)         = 0.56909076642393639985E0f;
  A(52)         = 0.00000000000155296588E0f;
  A(53)         = (-0.00000000003032205868E0f);
  A(54)         = 0.00000000040424830707E0f;
  A(55)         = (-0.00000000471135111493E0f);
  A(56)         = 0.00000005011915876293E0f;
  A(57)         = (-0.00000048722516178974E0f);
  A(58)         = 0.00000430683284629395E0f;
  A(59)         = (-0.00003445026145385764E0f);
  A(60)         = 0.00024879276133931664E0f;
  A(61)         = (-0.00162940941748079288E0f);
  A(62)         = 0.00988786373932350462E0f;
  A(63)         = (-0.05962426839442303805E0f);
  A(64)         = 0.49766113250947636708E0f;
  B(0)          = (-0.00000000029734388465E0f);
  B(1)          = 0.00000000269776334046E0f;
  B(2)          = (-0.00000000640788827665E0f);
  B(3)          = (-0.00000001667820132100E0f);
  B(4)          = (-0.00000021854388148686E0f);
  B(5)          = 0.00000266246030457984E0f;
  B(6)          = 0.00001612722157047886E0f;
  B(7)          = (-0.00025616361025506629E0f);
  B(8)          = 0.00015380842432375365E0f;
  B(9)          = 0.00815533022524927908E0f;
  B(10)         = (-0.01402283663896319337E0f);
  B(11)         = (-0.19746892495383021487E0f);
  B(12)         = 0.71511720328842845913E0f;
  B(13)         = (-0.00000000001951073787E0f);
  B(14)         = (-0.00000000032302692214E0f);
  B(15)         = 0.00000000522461866919E0f;
  B(16)         = 0.00000000342940918551E0f;
  B(17)         = (-0.00000035772874310272E0f);
  B(18)         = 0.00000019999935792654E0f;
  B(19)         = 0.00002687044575042908E0f;
  B(20)         = (-0.00011843240273775776E0f);
  B(21)         = (-0.00080991728956032271E0f);
  B(22)         = 0.00661062970502241174E0f;
  B(23)         = 0.00909530922354827295E0f;
  B(24)         = (-0.20160072778491013140E0f);
  B(25)         = 0.51169696718727644908E0f;
  B(26)         = 0.00000000003147682272E0f;
  B(27)         = (-0.00000000048465972408E0f);
  B(28)         = 0.00000000063675740242E0f;
  B(29)         = 0.00000003377623323271E0f;
  B(30)         = (-0.00000015451139637086E0f);
  B(31)         = (-0.00000203340624738438E0f);
  B(32)         = 0.00001947204525295057E0f;
  B(33)         = 0.00002854147231653228E0f;
  B(34)         = (-0.00101565063152200272E0f);
  B(35)         = 0.00271187003520095655E0f;
  B(36)         = 0.02328095035422810727E0f;
  B(37)         = (-0.16725021123116877197E0f);
  B(38)         = 0.32490054966649436974E0f;
  B(39)         = 0.00000000002319363370E0f;
  B(40)         = (-0.00000000006303206648E0f);
  B(41)         = (-0.00000000264888267434E0f);
  B(42)         = 0.00000002050708040581E0f;
  B(43)         = 0.00000011371857327578E0f;
  B(44)         = (-0.00000211211337219663E0f);
  B(45)         = 0.00000368797328322935E0f;
  B(46)         = 0.00009823686253424796E0f;
  B(47)         = (-0.00065860243990455368E0f);
  B(48)         = (-0.00075285814895230877E0f);
  B(49)         = 0.02585434424202960464E0f;
  B(50)         = (-0.11637092784486193258E0f);
  B(51)         = 0.18267336775296612024E0f;
  B(52)         = (-0.00000000000367789363E0f);
  B(53)         = 0.00000000020876046746E0f;
  B(54)         = (-0.00000000193319027226E0f);
  B(55)         = (-0.00000000435953392472E0f);
  B(56)         = 0.00000018006992266137E0f;
  B(57)         = (-0.00000078441223763969E0f);
  B(58)         = (-0.00000675407647949153E0f);
  B(59)         = 0.00008428418334440096E0f;
  B(60)         = (-0.00017604388937031815E0f);
  B(61)         = (-0.00239729611435071610E0f);
  B(62)         = 0.02064129023876022970E0f;
  B(63)         = (-0.06905562880005864105E0f);
  B(64)         = 0.09084526782065478489E0f;
  W             = Kokkos::abs(X);
  if ((W < 2.2e0)) {
    T = (W * W);
    K = static_cast<int>(T);
    T = (T - K);
    K = ((13 * K));
    Y = (((((((((((((((((((((((((A(K) * T) + A(((K + 1)))) * T) + A(((K + 2)))) * T) +
                            A(((K + 3)))) *
                           T) +
                          A(((K + 4)))) *
                         T) +
                        A(((K + 5)))) *
                       T) +
                      A(((K + 6)))) *
                     T) +
                    A(((K + 7)))) *
                   T) +
                  A(((K + 8)))) *
                 T) +
                A(((K + 9)))) *
               T) +
              A(((K + 10)))) *
             T) +
            A(((K + 11)))) *
           T) +
          A(((K + 12)))) *
         W);
  } else {
    if ((W < 6.9e0)) {
      K = static_cast<int>(W);
      T = (W - K);
      K = ((13 * K - 26));
      Y = ((((((((((((((((((((((((B(K) * T) + B(((K + 1)))) * T) + B(((K + 2)))) * T) +
                             B(((K + 3)))) *
                            T) +
                           B(((K + 4)))) *
                          T) +
                         B(((K + 5)))) *
                        T) +
                       B(((K + 6)))) *
                      T) +
                     B(((K + 7)))) *
                    T) +
                   B(((K + 8)))) *
                  T) +
                 B(((K + 9)))) *
                T) +
               B(((K + 10)))) *
              T) +
             B(((K + 11)))) *
            T) +
           B(((K + 12))));
      Y = (Y * Y);
      Y = (Y * Y);
      Y = (Y * Y);
      Y = (1 - (Y * Y));
    } else {
      Y = 1;
    }
  }
  if ((X < 0)) {
    Y = (-Y);
  }
  DERF1 = Y;
  return DERF1;
}
inline void radar_init(morrison_microphysics_state& morrison) {
  // Execute the Native radar_init body and update its declared outputs.
  int n;

  morrison.PI5    = ((((morrison.PI * morrison.PI) * morrison.PI) * morrison.PI) * morrison.PI);
  morrison.lamda4 = (((morrison.lamda_radar * morrison.lamda_radar) * morrison.lamda_radar) *
                     morrison.lamda_radar);
  morrison.m_w_0  = morrison_engine::m_complex_water_ray(morrison.lamda_radar, 0.0e0, morrison);
  morrison.m_i_0  = morrison_engine::m_complex_ice_maetzler(morrison.lamda_radar, 0.0e0, morrison);
  morrison.K_w =
      (std::abs(
           (std::complex<double>((std::complex<double>((std::complex<double>(morrison.m_w_0) *
                                                        std::complex<double>(morrison.m_w_0))) -
                                  std::complex<double>(1.0f))) /
            std::complex<double>((std::complex<double>((std::complex<double>(morrison.m_w_0) *
                                                        std::complex<double>(morrison.m_w_0))) +
                                  std::complex<double>(2.0f))))) *
       std::abs(
           (std::complex<double>((std::complex<double>((std::complex<double>(morrison.m_w_0) *
                                                        std::complex<double>(morrison.m_w_0))) -
                                  std::complex<double>(1.0f))) /
            std::complex<double>((std::complex<double>((std::complex<double>(morrison.m_w_0) *
                                                        std::complex<double>(morrison.m_w_0))) +
                                  std::complex<double>(2.0f))))));
  for (n = 0; n < morrison.nbins + 1; n += 1) {
    morrison.simpson[n] = 0.0e0;
  }
  for (n = 0; n < morrison.nbins - 1; n += 2) {
    morrison.simpson[n]         = (morrison.simpson[n] + morrison.basis[0]);
    morrison.simpson[((n + 1))] = (morrison.simpson[((n + 1))] + morrison.basis[1]);
    morrison.simpson[((n + 2))] = (morrison.simpson[((n + 2))] + morrison.basis[2]);
  }
  for (n = 0; n < morrison.slen; n += 1) {
    morrison.mixingrulestring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.matrixstring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.inclusionstring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hoststring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hostmatrixstring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hostinclusionstring_s.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.mixingrulestring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.matrixstring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.inclusionstring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hoststring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hostmatrixstring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
    morrison.hostinclusionstring_g.replace(n, 1, std::string(1, static_cast<char>(0)));
  }
  morrison.mixingrulestring_s = std::string("maxwellgarnett");
  morrison.mixingrulestring_s.resize(morrison.slen, ' ');
  morrison.hoststring_s = std::string("air");
  morrison.hoststring_s.resize(morrison.slen, ' ');
  morrison.matrixstring_s = std::string("water");
  morrison.matrixstring_s.resize(morrison.slen, ' ');
  morrison.inclusionstring_s = std::string("spheroidal");
  morrison.inclusionstring_s.resize(morrison.slen, ' ');
  morrison.hostmatrixstring_s = std::string("icewater");
  morrison.hostmatrixstring_s.resize(morrison.slen, ' ');
  morrison.hostinclusionstring_s = std::string("spheroidal");
  morrison.hostinclusionstring_s.resize(morrison.slen, ' ');
  morrison.mixingrulestring_g = std::string("maxwellgarnett");
  morrison.mixingrulestring_g.resize(morrison.slen, ' ');
  morrison.hoststring_g = std::string("air");
  morrison.hoststring_g.resize(morrison.slen, ' ');
  morrison.matrixstring_g = std::string("water");
  morrison.matrixstring_g.resize(morrison.slen, ' ');
  morrison.inclusionstring_g = std::string("spheroidal");
  morrison.inclusionstring_g.resize(morrison.slen, ' ');
  morrison.hostmatrixstring_g = std::string("icewater");
  morrison.hostmatrixstring_g.resize(morrison.slen, ' ');
  morrison.hostinclusionstring_g = std::string("spheroidal");
  morrison.hostinclusionstring_g.resize(morrison.slen, ' ');
}
inline std::complex<double> m_complex_water_ray(double lambda, double T,
                                                morrison_microphysics_state& morrison) {
  // Execute the Native m_complex_water_ray body and update its declared outputs.
  double epsinf;
  double epss;
  double epsr;
  double epsi;
  double alpha;
  double lambdas;
  double sigma;
  double nenner;
  const std::complex<double> i = std::complex<double>(0., 1.);
  std::complex<double> m_complex_water_ray;

  epsinf = ((5.27137e0 + (0.02164740e0 * T)) - ((0.00131198e0 * T) * T));
  epss =
      (78.54e+0 * (((1.0f - (4.579e-3 * (T - 25.0f))) + ((1.190e-5 * (T - 25.0f)) * (T - 25.0f))) -
                   (((2.800e-8 * (T - 25.0f)) * (T - 25.0f)) * (T - 25.0f))));
  alpha   = (((-16.8129e0) / (T + 273.16f)) + 0.0609265e0);
  lambdas = ((0.00033836e0 * Kokkos::exp((2513.98e0 / (T + 273.16f)))) * 1e-2f);
  nenner  = ((1.e0 + ((2.e0 * Kokkos::pow((lambdas / lambda), (1e0 - alpha))) *
                      Kokkos::sin(((alpha * morrison.PI) * 0.5f)))) +
             Kokkos::pow((lambdas / lambda), (2e0 - (2e0 * alpha))));
  epsr    = (epsinf + (((epss - epsinf) * ((Kokkos::pow((lambdas / lambda), (1e0 - alpha)) *
                                            Kokkos::sin(((alpha * morrison.PI) * 0.5f))) +
                                           1e0)) /
                       nenner));
  epsi    = ((((epss - epsinf) * ((Kokkos::pow((lambdas / lambda), (1e0 - alpha)) *
                                   Kokkos::cos(((alpha * morrison.PI) * 0.5f))) +
                                  0e0)) /
              nenner) +
             ((lambda * 1.25664f) / 1.88496f));
  m_complex_water_ray =
      std::sqrt(std::complex<float>(static_cast<float>(epsr), static_cast<float>((-epsi))));
  return m_complex_water_ray;
}
inline std::complex<double> m_complex_ice_maetzler(double lambda, double T,
                                                   morrison_microphysics_state& morrison) {
  // Execute the Native m_complex_ice_maetzler body and update its declared outputs.
  double f;
  double c;
  double TK;
  double B1;
  double B2;
  double b;
  double deltabeta;
  double betam;
  double beta;
  double theta;
  double alfa;
  std::complex<double> m_complex_ice_maetzler;

  c         = 2.99e8;
  TK        = (T + 273.16f);
  f         = ((c / lambda) * 1e-9);
  B1        = 0.0207f;
  B2        = 1.16e-11;
  b         = 335.0e0;
  deltabeta = Kokkos::exp(((-10.02f) + (0.0364f * (TK - 273.16f))));
  betam = (((B1 / TK) *
            (Kokkos::exp((b / TK)) / ((Kokkos::exp((b / TK)) - 1) * (Kokkos::exp((b / TK)) - 1)))) +
           ((B2 * f) * f));
  beta  = (betam + deltabeta);
  theta = ((300.f / TK) - 1.f);
  alfa  = ((0.00504e0 + (0.0062e0 * theta)) * Kokkos::exp(((-22.1e0) * theta)));
  m_complex_ice_maetzler = (3.1884f + (9.1e-4f * (TK - 273.16f)));
  m_complex_ice_maetzler =
      (std::complex<double>(m_complex_ice_maetzler) +
       std::complex<double>(std::complex<float>(static_cast<float>(0.0e0),
                                                static_cast<float>(((alfa / f) + (beta * f))))));
  m_complex_ice_maetzler = std::sqrt(std::conj(m_complex_ice_maetzler));
  return m_complex_ice_maetzler;
}
inline void POSITIVE_QV_ADJ(float& QV, float& QC, float& QR, float& QI, float& QS, float& QG,
                            float& T_IN_K, morrison_microphysics_state& morrison) {
  // Execute the Native POSITIVE_QV_ADJ body and update its declared outputs.
  const float EPS = 1.E-12f;
  float QT_COND_LIQ;
  float QT_COND_ICE;
  float QT_TOTAL;
  float DELTA_QV;
  float DELTA_QT_COND_LIQ;
  float DELTA_QT_COND_ICE;
  float REDUCE_COEF;

  if ((QV >= EPS)) {
    return;
  }
  QT_COND_LIQ = (QC + QR);
  QT_COND_ICE = 0.0f;
  if ((morrison.ILIQ == 0)) {
    QT_COND_ICE = ((QT_COND_ICE + QS) + QI);
  }
  if ((morrison.IGRAUP == 0)) {
    QT_COND_ICE = (QT_COND_ICE + QG);
  }
  QT_TOTAL = ((QV + QT_COND_LIQ) + QT_COND_ICE);
  if ((QT_TOTAL < (2 * EPS))) {
    return;
  }
  DELTA_QV          = (EPS - QV);
  QV                = EPS;
  REDUCE_COEF       = (1.f - (DELTA_QV / (QT_COND_LIQ + QT_COND_ICE)));
  QC                = (QC * REDUCE_COEF);
  QR                = (QR * REDUCE_COEF);
  DELTA_QT_COND_LIQ = (QT_COND_LIQ - (QC + QR));
  DELTA_QT_COND_ICE = 0.0f;
  if ((morrison.ILIQ == 0)) {
    QI = (QI * REDUCE_COEF);
    QS = (QS * REDUCE_COEF);
    if ((morrison.IGRAUP != 0)) {
      DELTA_QT_COND_ICE = (QT_COND_ICE - (QI + QS));
    }
  }
  if ((morrison.IGRAUP == 0)) {
    QG                = (QG * REDUCE_COEF);
    DELTA_QT_COND_ICE = (QT_COND_ICE - ((QI + QS) + QG));
  }
  T_IN_K = ((T_IN_K - ((Lv / Cp) * DELTA_QT_COND_LIQ)) - ((Ls / Cp) * DELTA_QT_COND_ICE));
  return;
}
} // namespace clubb_kokkos::morrison_engine
