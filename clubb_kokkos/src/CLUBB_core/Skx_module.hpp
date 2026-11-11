// Skewness and third-moment closures from src/CLUBB_core/Skx_module.F90.
//
// Skx_func computes skewness with the native small-variance regularization.
// compute_gamma_Skw gives the PDF width coefficient as a function of skewness.
// LG_2005_ansatz diagnoses a responder variable's skewness, and
// xp3_LG_2005_ansatz converts that diagnosis to a third moment.
// Inputs and outputs are execution-space Views, indexed (column, level).
// Reference: Larson and Golaz (2005), Mon. Wea. Rev., 133, 1023–1042.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"

namespace clubb_kokkos {

void Skx_func(int nz, int ngrdcol, Real2D xp2, Real2D xp3, // In
              real_t x_tol, Real2D clubb_params,           // In
              Real2D Skx);                                 // Out

void compute_gamma_Skw(int nz, int ngrdcol, Real2D Skw, Real2D clubb_params, // In
                       Real2D gamma_Skw_fnc);                                // Out

void LG_2005_ansatz(int nz, int ngrdcol, Real2D Skw, Real2D wpxp, Real2D wp2,  // In
                    Real2D xp2, Real1D beta, Real2D sigma_sqd_w, real_t x_tol, // In
                    Real2D Skx);                                               // Out

void xp3_LG_2005_ansatz(ScratchWorkspace& workspace,                                        // InOut
                        int nzt, int ngrdcol, Real2D Skw_zt, Real2D wpxp_zt, Real2D wp2_zt, // In
                        Real2D xp2_zt, Real2D sigma_sqd_w_zt,                               // In
                        Real2D clubb_params, real_t x_tol,                                  // In
                        Real2D xp3);                                                        // Out

} // namespace clubb_kokkos
