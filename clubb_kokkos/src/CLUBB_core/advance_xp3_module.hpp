// Predict and diagnose third-order scalar and horizontal-velocity moments.
//
// advance_xp3 advances rt/thl/passive-scalar third moments with the Native
// simplified equation, then diagnoses u/v third moments. diagnose_xp3 applies
// the Larson–Golaz ansatz to all those fields with the selected PDF width.
// Inputs/outputs are execution-space LayoutLeft Views indexed (column,level)
// or (column,level,scalar). Native units are listed in the interfaces.
// Private production/accumulation expressions execute within independent cell
// kernels. Stats retain the existing Fortran backend; no state is kept here.
// Native predictive/LMM and steady-state branches are retained. The last
// storage level is the Native upper boundary in either grid direction.
// Reference: Larson and Golaz (2005) for the diagnostic ansatz.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
#include "stats_netcdf.hpp"
namespace clubb_kokkos {
void advance_xp3(
    ScratchWorkspace& workspace, int nzm, // In: Number of momentum vertical levels
    int nzt,                              // In: Number of thermodynamic vertical levels
    int ngrdcol,                          // In: Number of grid columns
    int sclr_dim,                         // In: Number of passive scalars
    Real1D sclr_tol,                      // In: Threshold(s) on the passive scalars  [units vary]
    const grid& gr,                       // In
    real_t dt,                            // In: Model timestep                            [s]
    Real2D rtm,    // In: Mean (overall) of rt (thermo. levels)                [kg/kg]
    Real2D thlm,   // In: Mean (overall) of thl (thermo. levels)               [K]
    Real2D rtp2,   // In: Variance (overall) of rt (m-levs.)                   [kg^2/kg^2]
    Real2D thlp2,  // In: Variance (overall) of thl (m-levs.)                  [K^2]
    Real2D wprtp,  // In: Turbulent flux of rt (momentum levs.)                [m/s kg/kg]
    Real2D wpthlp, // In: Turbulent flux of thl (momentum levs.)               [m/s K]
    // In: <w'rt'^2> (thermodynamic levels)                     [m/s(kg/kg)^2]
    Real2D wprtp2,
    Real2D wpthlp2,   // In: <w'thl'^2> (thermodynamic levels)                    [m/s K^2]
    Real2D rho_ds_zm, // In: Dry, static density on momentum levels               [kg/m^3]
    // In: Inv. dry, static density at thermo. levels           [m^3/kg]
    Real2D invrs_rho_ds_zt,
    Real2D invrs_tau_zt, // In: Inverse time-scale tau on thermodynamic levels       [1/s]
    Real2D tau_max_zt,   // In: Max. allowable eddy dissipation time scale on t-levs [s]
    Real3D sclrm,        // In: Mean (overall) of sclr (thermo. levels) [sclr units]
    Real3D sclrp2,       // In: Variance (overall) of sclr (m-levs.)    [(sclr units)^2]
    Real3D wpsclrp,      // In: Turbulent flux of sclr (momentum levs.) [m/s(sclr units)]
    Real3D wpsclrp2,     // In: <w'sclr'^2> (thermodynamic levels)      [m/s(sclr units)^2]
    Real2D wp2,          // In: w'^2 (momentum levels)                              [m^2/s^2]
    Real2D wp3,          // In: w'^3 on thermo. grid     [m^3/s^3]
    Real2D upwp,         // In: u'w' (momentum levels)                              [m^2/s^2]
    Real2D vpwp,         // In: v'w' (momentum levels)                              [m^2/s^2]
    Real2D up2,          // In: u'^2 (momentum levels)                              [m^2/s^2]
    Real2D vp2,          // In: v'^2 (momentum levels)                              [m^2/s^2]
    Real2D thvm,         // In: Virtual potential temperature on thermodynamic levels [K]
    Real2D clubb_params, // In: Array of CLUBB's tunable parameters    [units vary]
    bool l_lmm_stepping, // In: Apply Linear Multistep Method (LMM) Stepping
    stats_type& stats,   // Inout
    Real2D rtp3,         // Inout: <rt'^3> (thermodynamic levels)     [kg^3/kg^3]
    Real2D thlp3,        // Inout: <thl'^3> (thermodynamic levels)    [K^3]
    Real3D sclrp3,       // Inout: <sclr'^3> (thermodynamic levels)    [(sclr units)^3]
    Real2D up3,          // Inout: u'^3 (thermodynamic levels)        [m^3/s^3]
    // Inout: v'^3 (thermodynamic levels)        [m^3/s^3]
    Real2D vp3);

void diagnose_xp3(ScratchWorkspace& workspace, int nzm, // In: Number of momentum levels
                  int nzt,                              // In: Number of thermodynamic levels
                  int ngrdcol,                          // In: Number of grid columns
                  int sclr_dim,                         // In: Number of passive scalars
                  Real1D sclr_tol,     // In: Thresholds for passive scalars [units vary]
                  const grid& gr,      // In: Grid structure
                  int iiPDF_type,      // In: PDF type selector
                  Real2D clubb_params, // In: Array of CLUBB tunable parameters [units vary]
                  Real2D wp2,          // In: Variance of vertical velocity [m^2/s^2]
                  Real2D wp3,          // In: Third moment of vertical velocity [m^3/s^3]
                  Real2D thvm,         // In: Virtual potential temperature [K]
                  Real2D wprtp,        // In: Turbulent flux of rt [m kg/(kg s)]
                  Real2D wpthlp,       // In: Turbulent flux of thl [m K/s]
                  Real2D rtp2,         // In: Variance of rt [(kg/kg)^2]
                  Real2D thlp2,        // In: Variance of thl [K^2]
                  Real2D upwp,         // In: Turbulent flux of u [m^2/s^2]
                  Real2D vpwp,         // In: Turbulent flux of v [m^2/s^2]
                  Real2D up2,          // In: Variance of u wind [m^2/s^2]
                  Real2D vp2,          // In: Variance of v wind [m^2/s^2]
                  Real2D sigma_sqd_w,  // In: PDF width parameter [-]
                  Real3D wpsclrp,      // In: Turbulent flux of passive scalars [m sclr/s]
                  Real3D sclrp2,       // In: Variance of passive scalars [units vary]
                  stats_type& stats,   // Inout: Statistics accumulator
                  Real2D rtp3,         // Inout: Third moment of rt [(kg/kg)^3]
                  Real2D thlp3,        // Inout: Third moment of thl [K^3]
                  Real2D up3,          // Inout: Third moment of u wind [m^3/s^3]
                  Real2D vp3,          // Inout: Third moment of v wind [m^3/s^3]
                  // Inout: Third moments of passive scalars [units vary]
                  Real3D sclrp3);

} // namespace clubb_kokkos
