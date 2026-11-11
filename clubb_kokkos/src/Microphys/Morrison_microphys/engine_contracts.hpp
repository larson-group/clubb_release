// Morrison column bounds are inclusive zero-based C++ storage coordinates.
// Counts and physical recurrence/iteration parameters retain their Native values.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/microphys_host_types.hpp"
namespace clubb_kokkos::morrison_engine {
using MorrisonColumn = Kokkos::View<float*, Kokkos::LayoutStride, Kokkos::HostSpace>;
void M2005MICRO_GRAUPEL(
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
} // namespace clubb_kokkos::morrison_engine
