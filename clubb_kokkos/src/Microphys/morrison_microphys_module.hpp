// Coordinate the Morrison two-moment scheme and its CLUBB diagnostics.
//
// morrison_microphys_driver stages each column in the engine's Native default
// real kind, calls M2005MICRO_GRAUPEL, restores core-real tendencies/velocities
// and updates process budgets and surface precipitation statistics. Internal
// helpers keep SILHS weighting and detailed nonfinite diagnostics in source order.
// Arrays are column-first LayoutLeft (column,level[,species]); the appended
// case records own configuration and initialized Morrison coefficients/tables.
// Each host engine call uses local float arrays, with no shared model storage.
// Reference: Morrison et al.(2005), JAS62,1665–1677.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/error_code.hpp"
#include "Microphys/microphys_host_types.hpp"
#include "Microphys/KK_microphys/KK_microphys.hpp"
#include "Microphys/advance_microphys_module.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel.hpp"
namespace clubb_kokkos::microphys_host {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void morrison_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel_zt,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, morrison_microphysics_state& morrison);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void morrison_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel_zt,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, morrison_microphysics_state& morrison) {
  // Description:
  // Wrapper for the Morrison microphysics, applied to all columns.
  //
  // References:
  // None
  //-----------------------------------------------------------------------
  // The version of the Morrison 2005 microphysics that is in SAM.
  // External
  // Input Variables
  // Output Variables
  // Local Variables
  // In the comments below, by "adds to" we mean that if the quantity is
  // positive, it adds positively to the prognostic variable, but if the
  // quantity is negative, it subtracts from the prognostic variable.
  // Local Variables
  // Just to avoid typing hm_metadata%iixx everywhere
  // In: gr
  // In: ngrdcol — Number of model columns       [-]
  // In: dt — Model timestep        [s]
  // In: nzt — Points in the Vertical      [-]
  // In: hydromet_dim — Number of hydrometeor species [-]
  // In: hm_metadata
  // In: l_latin_hypercube — Whether we're using latin hypercube sampling
  // In: thlm (ngrdcol,nzt) — Liquid potential temperature       [K]
  // In: wm_zt (ngrdcol,nzt) — Mean vertical velocity on the thermo grid     [m/s]
  // In: p_in_Pa (ngrdcol,nzt) — Pressure                           [Pa]
  // In: exner (ngrdcol,nzt) — Exner function                     [-]
  // In: rho (ngrdcol,nzt) — Density on thermo. grid            [kg/m^3]
  // In: cloud_frac (ngrdcol,nzt) — Cloud fraction                     [-]
  // In: w_std_dev (ngrdcol,nzt) — Standard deviation of vertical vel. [m/s]
  // In: dzq (ngrdcol,nzt) — Change in altitude                  [m]
  // In: rcm (ngrdcol,nzt) — Cloud water mixing ratio                  [kg/kg]
  // In: Ncm (ngrdcol,nzt) — Grid mean value for cloud droplet conc.    [#/kg]
  // In: chi (ngrdcol,nzt) — The variable 's' from Mellor              [kg/kg]
  // In: rvm (ngrdcol,nzt) — Vapor water mixing ratio                  [kg/kg]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor species    [units vary]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // In: sample_weight (ngrdcol,nzt) — SILHS sample weight; unity for an ordinary call
  // Inout: stats
  // Out: hydromet_mc (ngrdcol,nzt,hydromet_dim) — Hydrometeor time tendency          [(units vary)/s]
  // Out: hydromet_vel_zt (ngrdcol,nzt,hydromet_dim) — Hydrometeor sedimentation velocity [m/s]
  // Out: Ncm_mc (ngrdcol,nzt) — Cloud droplet concentration time tendency    [num/kg/s]
  // Out: rcm_mc (ngrdcol,nzt) — Time tendency of liquid water mixing ratio    [kg/kg/s]
  // Out: rvm_mc (ngrdcol,nzt) — Time tendency of vapor water mixing ratio     [kg/kg/s]
  // Out: thlm_mc (ngrdcol,nzt) — Time tendency of liquid potential temperature [K/s]
  // Out: rrm_auto_diag (ngrdcol,nzt)
  // Out: rrm_accr_diag (ngrdcol,nzt)
  // Out: rrm_evap_diag (ngrdcol,nzt)
  // Out: Nrm_auto_diag (ngrdcol,nzt)
  // Out: Nrm_evap_diag (ngrdcol,nzt)
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: morrison replaces the Native morrison_microphysics_state module.

  auto workspace_scope        = workspace.scope();
  const auto gr_zt            = to_host_space(workspace, gr.zt);
  const auto gr_dzt           = to_host_space(workspace, gr.dzt);
  const auto hm_metadata_iirr = hm_metadata.iirr;
  const auto hm_metadata_iirs = hm_metadata.iirs;
  const auto hm_metadata_iiri = hm_metadata.iiri;
  const auto hm_metadata_iirg = hm_metadata.iirg;
  const auto hm_metadata_iiNr = hm_metadata.iiNr;
  const auto hm_metadata_iiNs = hm_metadata.iiNs;
  const auto hm_metadata_iiNi = hm_metadata.iiNi;
  const auto hm_metadata_iiNg = hm_metadata.iiNg;
  constexpr real_t w_thresh   = real_t(0.1);
  HostReal2D rrm_auto =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rrm_auto", ngrdcol,
                                       nzt); // Autoconversion rate                 [kg/kg/s]
  HostReal2D rrm_accr =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rrm_accr", ngrdcol,
                                       nzt); // Accretion rate                      [kg/kg/s]
  HostReal2D rrm_evap =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rrm_evap", ngrdcol,
                                       nzt); // Rain evaporation rate               [kg/kg/s]
  HostReal2D T_in_K_core =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::T_in_K_core", ngrdcol,
                                       nzt); // Temperature at CLUBB precision      [K]
  HostFloat2D effc   = workspace_scope.view<HostFloat2D>("morrison_microphys_driver::effc", ngrdcol,
                                                         nzt); // Effective droplet radii [μ]
  HostFloat2D effi   = workspace_scope.view<HostFloat2D>("morrison_microphys_driver::effi", ngrdcol,
                                                         nzt); // Effective droplet radii [μ]
  HostFloat2D effg   = workspace_scope.view<HostFloat2D>("morrison_microphys_driver::effg", ngrdcol,
                                                         nzt); // Effective droplet radii [μ]
  HostFloat2D effs   = workspace_scope.view<HostFloat2D>("morrison_microphys_driver::effs", ngrdcol,
                                                         nzt); // Effective droplet radii [μ]
  HostFloat2D effr   = workspace_scope.view<HostFloat2D>("morrison_microphys_driver::effr", ngrdcol,
                                                         nzt); // Effective droplet radii [μ]
  HostFloat2D T_in_K = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::T_in_K", ngrdcol,
      nzt); // Temperature                                      [K]
  HostFloat2D T_in_K_mc = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::T_in_K_mc", ngrdcol,
      nzt); // Temperature tendency                           [K/s]
  HostFloat2D rcm_r4 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rcm_r4", ngrdcol,
      nzt); // Temporary array for cloud water mixing ratio [kg/kg]
  HostFloat2D Ncm_r4 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::Ncm_r4", ngrdcol,
      nzt); // Temporary array for cloud number conc.        [#/kg]
  HostFloat2D Ncm_mc_r4 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::Ncm_mc_r4", ngrdcol,
      nzt); // Temporary array for cloud number conc.      [#/kg/s]
  HostFloat2D rvm_r4 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rvm_r4", ngrdcol,
      nzt); // Temporary array for vapor water mixing ratio [kg/kg]
  HostFloat2D rcm_sten = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rcm_sten", ngrdcol,
      nzt); // Mean rc sedimentation tendency             [kg/kg/s]
  HostFloat2D rrm_sten = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rrm_sten", ngrdcol,
      nzt); // Mean rr sedimentation tendency             [kg/kg/s]
  HostFloat2D rim_sten = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rim_sten", ngrdcol,
      nzt); // Mean ri sedimentation tendency             [kg/kg/s]
  HostFloat2D rsm_sten = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rsm_sten", ngrdcol,
      nzt); // Mean rs sedimentation tendency             [kg/kg/s]
  HostFloat2D rgm_sten = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::rgm_sten", ngrdcol,
      nzt); // Mean rg sedimentation tendency             [kg/kg/s]
  HostFloat2D morr_rain_vel_r4 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::morr_rain_vel_r4", ngrdcol,
      nzt); // Rain fall velocity from Morrison microphysics  [m/s]
  HostFloat2D cloud_frac_in = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::cloud_frac_in", ngrdcol,
      nzt); // Cloud frac used as input for the Morrison scheme [-]
  HostFloat2D PSMLT = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PSMLT", ngrdcol,
      nzt); // Freezing of rain to form snow. Adds to rsm, subtracts from rrm [kg/kg/s]
  HostFloat2D EVPMS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::EVPMS", ngrdcol,
      nzt); // Evaporation of melted snow. Adds to rsm, subtracts from rvm [kg/kg/s]
  HostFloat2D PRACS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRACS", ngrdcol,
      nzt); // Collection of rain by snow. Adds to rsm, subtracts from rrm [kg/kg/s]
  HostFloat2D EVPMG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::EVPMG", ngrdcol,
      nzt); // Evaporation of melted graupel. Adds to rgm, subtracts from rvm [kg/kg/s]
  HostFloat2D PRACG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRACG", ngrdcol,
      nzt); // Negative of collection of rain by graupel. Adds to rrm, subtracts from rgm [kg/kg/s]
  HostFloat2D PGMLT = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PGMLT", ngrdcol,
      nzt); // Negative of melting of graupel. Adds to rgm, subtracts from rrm [kg/kg/s]
  HostFloat2D MNUCCC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::MNUCCC", ngrdcol,
      nzt); // Contact freezing of cloud droplets. Adds to rim, subtracts from rcm [kg/kg/s]
  HostFloat2D PSACWS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PSACWS", ngrdcol,
      nzt); // Collection of cloud water by snow. Adds to rsm, subtracts from rcm [kg/kg/s]
  HostFloat2D PSACWI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PSACWI", ngrdcol,
      nzt); // Collection of cloud water by cloud ice. Adds to rim, subtracts from rcm [kg/kg/s]
  HostFloat2D QMULTS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QMULTS", ngrdcol,
      nzt); // Splintering from cloud droplets accreted onto snow. Adds to rim, subtracts from rcm [kg/kg/s]
  HostFloat2D QMULTG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QMULTG", ngrdcol,
      nzt); // Splintering from droplets accreted onto graupel. Adds to rim, subtracts from rcm [kg/kg/s]
  HostFloat2D PSACWG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PSACWG", ngrdcol,
      nzt); // Collection of cloud water by graupel. Adds to rgm, subtracts from rcm [kg/kg/s]
  HostFloat2D PGSACW = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PGSACW", ngrdcol,
      nzt); // Reclassification of rimed snow as graupel. Adds to rgm, subtracts from rcm [kg/kg/s]
  HostFloat2D PRD = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRD", ngrdcol,
      nzt); // Depositional growth of cloud ice. Adds to rim, subtracts from rcm [kg/kg/s]
  HostFloat2D PRCI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRCI", ngrdcol,
      nzt); // Autoconversion of cloud ice to snow. Adds to rsm, subtracts from rim [kg/kg/s]
  HostFloat2D PRAI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRAI", ngrdcol,
      nzt); // Collection of cloud ice by snow. Adds to rsm, subtracts from rim [kg/kg/s]
  HostFloat2D QMULTR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QMULTR", ngrdcol,
      nzt); // Splintering from rain droplets accreted onto snow. Adds to rim, subtracts from rrm [kg/kg/s]
  HostFloat2D QMULTRG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QMULTRG", ngrdcol,
      nzt); // Splintering from rain droplets accreted onto graupel. Adds to rim, subtracts from rrm [kg/kg/s]
  HostFloat2D MNUCCD = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::MNUCCD", ngrdcol,
      nzt); // Freezing of aerosol. Adds to rim, subtracts from rvm [kg/kg/s]
  HostFloat2D PRACI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRACI", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to rgm, subtracts from rim [kg/kg/s]
  HostFloat2D PRACIS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRACIS", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to rsm, subtracts from rim [kg/kg/s]
  HostFloat2D EPRD = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::EPRD", ngrdcol,
      nzt); // Negative of sublimation of cloud ice. Adds to rim, subtracts from rvm [kg/kg/s]
  HostFloat2D MNUCCR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::MNUCCR", ngrdcol,
      nzt); // Contact freezing of rain droplets. Adds to rgm, subtracts from rrm [kg/kg/s]
  HostFloat2D PIACR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PIACR", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to rgm, subtracts from rrm [kg/kg/s]
  HostFloat2D PIACRS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PIACRS", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to rsm, subtracts from rrm [kg/kg/s]
  HostFloat2D PGRACS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PGRACS", ngrdcol,
      nzt); // Collection of rain by snow. Adds to rgm, subtracts from rrm [kg/kg/s]
  HostFloat2D PRDS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRDS", ngrdcol,
      nzt); // Depositional growth of snow. Adds to rsm, subtracts from rvm [kg/kg/s]
  HostFloat2D EPRDS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::EPRDS", ngrdcol,
      nzt); // Negative of sublimation of snow. Adds to rsm, subtracts from rvm [kg/kg/s]
  HostFloat2D PSACR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PSACR", ngrdcol,
      nzt); // Collection of snow by rain. Adds to rgm, subtracts from rsm [kg/kg/s]
  HostFloat2D PRDG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRDG", ngrdcol,
      nzt); // Depositional growth of graupel. Adds to rgm, subtracts from rvm [kg/kg/s]
  HostFloat2D EPRDG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::EPRDG", ngrdcol,
      nzt); // Negative of sublimation of graupel. Adds to rgm, subtracts from rvm [kg/kg/s]
  HostFloat2D NGSTEN = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NGSTEN", ngrdcol, nzt); // Graupel sedimentation tendency [#/kg/s]
  HostFloat2D NRSTEN = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NRSTEN", ngrdcol, nzt); // Rain sedimentation tendency [#/kg/s]
  HostFloat2D NISTEN =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NISTEN", ngrdcol,
                                        nzt); // Cloud ice sedimentation tendency [#/kg/s]
  HostFloat2D NSSTEN = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSSTEN", ngrdcol, nzt); // Snow sedimentation tendency [#/kg/s]
  HostFloat2D NCSTEN =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NCSTEN", ngrdcol,
                                        nzt); // Cloud water sedimentation tendency [#/kg/s]
  HostFloat2D NPRC1 = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRC1", ngrdcol,
      nzt); // Change in Nrm due to autoconversion of droplets. Adds to Nrm [#/kg/s]
  HostFloat2D NRAGG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NRAGG", ngrdcol,
      nzt); // Change in Nrm due to self-collection of raindrops. Adds to Nrm [#/kg/s]
  HostFloat2D NPRACG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRACG", ngrdcol,
      nzt); // Collection of rainwater by graupel. Subtracts from Nrm [#/kg/s]
  HostFloat2D NSUBR =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NSUBR", ngrdcol,
                                        nzt); // Loss of Nrm by evaporation. Adds to Nrm [#/kg/s]
  HostFloat2D NSMLTR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSMLTR", ngrdcol,
      nzt); // Melting of snow to form rain. Subtracts from Nrm [#/kg/s]
  HostFloat2D NGMLTR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NGMLTR", ngrdcol,
      nzt); // Melting of graupel to form rain. Subtracts from Nrm [#/kg/s]
  HostFloat2D NPRACS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRACS", ngrdcol,
      nzt); // Collection of rainwater by snow. Subtracts from Nrm [#/kg/s]
  HostFloat2D NNUCCR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NNUCCR", ngrdcol,
      nzt); // Contact freezing of rain. Adds to Ngm, subtracts from Nrm [#/kg/s]
  HostFloat2D NIACR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NIACR", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to Ngm, subtracts from Nrm and Nim [#/kg/s]
  HostFloat2D NIACRS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NIACRS", ngrdcol,
      nzt); // Collection of cloud ice by rain. Adds to Nsm, subtracts from Nrm and Nim [#/kg/s]
  HostFloat2D NGRACS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NGRACS", ngrdcol,
      nzt); // Collection of rain by snow. Adds to Ngm, subtracts from Nrm and Nsm [#/kg/s]
  HostFloat2D NSMLTS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSMLTS", ngrdcol, nzt); // Melting of snow Adds to Nsm [#/kg/s]
  HostFloat2D NSAGG =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NSAGG", ngrdcol,
                                        nzt); // Self collection of snow Adds to Nsm [#/kg/s]
  HostFloat2D NPRCI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRCI", ngrdcol,
      nzt); // Autoconversion of cloud ice to snow Adds to Nsm, subtracts from Nim [#/kg/s]
  HostFloat2D NSCNG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSCNG", ngrdcol,
      nzt); // Conversion of snow to graupel Adds to Ngm, subtracts from Nsm [#/kg/s]
  HostFloat2D NSUBS =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NSUBS", ngrdcol,
                                        nzt); // Loss of Nsm due to sublimation Adds to Nsm [#/kg/s]
  HostFloat2D PRA = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRA", ngrdcol,
      nzt); // Accretion. Adds to rrm, subtracts from rcm [kg/kg/s]
  HostFloat2D PRC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PRC", ngrdcol,
      nzt); // Autoconversion. Adds to rrm, subtracts from rcm [kg/kg/s]
  HostFloat2D PRE =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::PRE", ngrdcol,
                                        nzt); // Rain evaporation. Subtracts from rrm [kg/kg/s]
  HostFloat2D PCC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::PCC", ngrdcol,
      nzt); // Saturation adjustment Adds to rcm, substracts from rvm [kg/kg/s]
  HostFloat2D NNUCCC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NNUCCC", ngrdcol,
      nzt); // Contact freezing of cloud drops Adds to Nim, subtracts from Ncm [#/kg/s]
  HostFloat2D NPSACWS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPSACWS", ngrdcol,
      nzt); // Droplet accretion by snow. Subtracts from Ncm [#/kg/s]
  HostFloat2D NPRA = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRA", ngrdcol,
      nzt); // Droplet accretion by rain. Subtracts from Ncm [#/kg/s]
  HostFloat2D NPRC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRC", ngrdcol,
      nzt); // Autoconversion of cloud drops. Subtracts from Ncm [#/kg/s]
  HostFloat2D NPSACWI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPSACWI", ngrdcol,
      nzt); // Droplet accretion by cloud ice. Subtracts from Ncm [#/kg/s]
  HostFloat2D NPSACWG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPSACWG", ngrdcol,
      nzt); // Collection of cloud drops by graupel. Subtracts from Ncm [#/kg/s]
  HostFloat2D NPRAI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NPRAI", ngrdcol,
      nzt); // Accretion of cloud ice by snow. Subtracts from Nim [#/kg/s]
  HostFloat2D NMULTS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NMULTS", ngrdcol,
      nzt); // Ice mult. due to riming of cloud drops by snow. Adds to Nim [#/kg/s]
  HostFloat2D NMULTG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NMULTG", ngrdcol,
      nzt); // Ice mult. due to accretion of cloud drops by graupel. Adds to Nim [#/kg/s]
  HostFloat2D NMULTR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NMULTR", ngrdcol,
      nzt); // Ice mult. due to riming of rain by snow. Adds to Nim [#/kg/s]
  HostFloat2D NMULTRG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NMULTRG", ngrdcol,
      nzt); // Ice mult. due to accretion of rain by graupel. Adds to Nim [#/kg/s]
  HostFloat2D NNUCCD = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NNUCCD", ngrdcol,
      nzt); // Primary ice nucleation, freezing of aerosol. Adds to Nim [#/kg/s]
  HostFloat2D NSUBI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSUBI", ngrdcol,
      nzt); // Loss of ice due to sublimation. Subtracts from Nim [#/kg/s]
  HostFloat2D NGMLTG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NGMLTG", ngrdcol,
      nzt); // Loss of graupel due to melting. Subtracts from Ngm [#/kg/s]
  HostFloat2D NSUBG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NSUBG", ngrdcol,
      nzt); // Loss of graupel due to sublimation. Subtracts from Ngm [#/kg/s]
  HostFloat2D NACT = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NACT", ngrdcol,
      nzt); // Cloud droplet formation by aerosol activation. Adds to Ncm [#/kg/s]
  HostFloat2D SIZEFIX_NR = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::SIZEFIX_NR", ngrdcol,
      nzt); // Adjustment to rain drop number concentration for large/small drops
  HostFloat2D SIZEFIX_NC = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::SIZEFIX_NC", ngrdcol,
      nzt); // Adjustment to cloud drop number concentration for large/small drops
  HostFloat2D SIZEFIX_NI = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::SIZEFIX_NI", ngrdcol,
      nzt); // Adjustment to ice number concentration for large/small drops
  HostFloat2D SIZEFIX_NS = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::SIZEFIX_NS", ngrdcol,
      nzt); // Adjustment to snow number concentration for large/small drops
  HostFloat2D SIZEFIX_NG = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::SIZEFIX_NG", ngrdcol,
      nzt); // Adjustment to graupel number concentration for large/small drops
  HostFloat2D NEGFIX_NI =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NEGFIX_NI", ngrdcol,
                                        nzt); // Removal of negative ice number concentration
  HostFloat2D NEGFIX_NS =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NEGFIX_NS", ngrdcol,
                                        nzt); // Removal of negative snow number concentration
  HostFloat2D NEGFIX_NC =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NEGFIX_NC", ngrdcol,
                                        nzt); // Removal of negative cloud drop number concentration
  HostFloat2D NEGFIX_NR =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NEGFIX_NR", ngrdcol,
                                        nzt); // Removal of negative rain drop number concentration
  HostFloat2D NEGFIX_NG =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NEGFIX_NG", ngrdcol,
                                        nzt); // Removal of negative graupel number concentration
  HostFloat2D NIM_MORR_CL =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::NIM_MORR_CL", ngrdcol,
                                        nzt); // Clipping of large ice number concentrations
  HostFloat2D QC_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QC_INST", ngrdcol,
      nzt); // Change in cloud mixing ratio due to instantaneous processes
  HostFloat2D QR_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QR_INST", ngrdcol,
      nzt); // Change in rain mixing ratio due to instantaneous processes
  HostFloat2D QI_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QI_INST", ngrdcol,
      nzt); // Change in ice mixing ratio due to instantaneous processes
  HostFloat2D QS_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QS_INST", ngrdcol,
      nzt); // Change in snow mixing ratio due to instantaneous processes
  HostFloat2D QG_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::QG_INST", ngrdcol,
      nzt); // Change in graupel mixing ratio due to instantaneous processes
  HostFloat2D NC_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NC_INST", ngrdcol,
      nzt); // Change in cloud number concentration due to instantaneous processes
  HostFloat2D NR_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NR_INST", ngrdcol,
      nzt); // Change in rain number concentration due to instantaneous processes
  HostFloat2D NI_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NI_INST", ngrdcol,
      nzt); // Change in ice number concentration due to instantaneous processes
  HostFloat2D NS_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NS_INST", ngrdcol,
      nzt); // Change in snow number concentration due to instantaneous processes
  HostFloat2D NG_INST = workspace_scope.view<HostFloat2D>(
      "morrison_microphys_driver::NG_INST", ngrdcol,
      nzt); // Change in graupel number concentration due to instantaneous processes
  HostReal2D rrm =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rrm", ngrdcol,
                                       nzt); // Mean rain water mixing ratio            [kg/kg]
  HostReal2D rim =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rim", ngrdcol,
                                       nzt); // Mean ice mixing ratio                   [kg/kg]
  HostReal2D rsm =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rsm", ngrdcol,
                                       nzt); // Mean snow mixing ratio                  [kg/kg]
  HostReal2D rgm =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::rgm", ngrdcol,
                                       nzt); // Mean graupel mixing ratio               [kg/kg]
  HostFloat1D Morr_snow_rate =
      workspace_scope.view<HostFloat1D>("morrison_microphys_driver::Morr_snow_rate", ngrdcol);
  HostFloat1D Morr_precip_rate =
      workspace_scope.view<HostFloat1D>("morrison_microphys_driver::Morr_precip_rate", ngrdcol);
  HostFloat3D hydromet_r4 = workspace_scope.view<HostFloat3D>(
      "morrison_microphys_driver::hydromet_r4", ngrdcol, nzt, hydromet_dim); // Temporary variable
  HostFloat3D hydromet_mc_r4 =
      workspace_scope.view<HostFloat3D>("morrison_microphys_driver::hydromet_mc_r4", ngrdcol, nzt,
                                        hydromet_dim); // Temporary variables
  HostFloat2D rrm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rrm_r4", ngrdcol,
                                        nzt); // Mean rain water mixing ratio            [kg/kg]
  HostFloat2D Nrm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nrm_r4", ngrdcol,
                                        nzt); // Mean rain drop concentration            [num/kg]
  HostFloat2D rim_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rim_r4", ngrdcol,
                                        nzt); // Mean ice mixing ratio                   [kg/kg]
  HostFloat2D Nim_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nim_r4", ngrdcol,
                                        nzt); // Mean ice crystal concentration          [num/kg]
  HostFloat2D rsm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rsm_r4", ngrdcol,
                                        nzt); // Mean snow mixing ratio                  [kg/kg]
  HostFloat2D Nsm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nsm_r4", ngrdcol,
                                        nzt); // Mean snow flake concentration           [num/kg]
  HostFloat2D rgm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rgm_r4", ngrdcol,
                                        nzt); // Mean graupel mixing ratio               [kg/kg]
  HostFloat2D Ngm_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Ngm_r4", ngrdcol,
                                        nzt); // Mean graupel concentration              [num/kg]
  HostFloat2D rrm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rrm_mc_r4", ngrdcol,
                                        nzt); // Mean rain water mixing ratio tendency   [kg/kg/s]
  HostFloat2D Nrm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nrm_mc_r4", ngrdcol,
                                        nzt); // Mean rain drop concentration tendency   [num/kg/s]
  HostFloat2D rim_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rim_mc_r4", ngrdcol,
                                        nzt); // Mean ice mixing ratio tendency          [kg/kg/s]
  HostFloat2D Nim_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nim_mc_r4", ngrdcol,
                                        nzt); // Mean ice crystal concentration tendency [num/kg/s]
  HostFloat2D rsm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rsm_mc_r4", ngrdcol,
                                        nzt); // Mean snow mixing ratio tendency         [kg/kg/s]
  HostFloat2D Nsm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Nsm_mc_r4", ngrdcol,
                                        nzt); // Mean snow flake concentration tendency  [num/kg/s]
  HostFloat2D rgm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rgm_mc_r4", ngrdcol,
                                        nzt); // Mean graupel mixing ratio tendency      [kg/kg/s]
  HostFloat2D Ngm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::Ngm_mc_r4", ngrdcol,
                                        nzt); // Mean graupel concentration tendency     [num/kg/s]
  HostFloat2D rcm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rcm_mc_r4", ngrdcol, nzt);
  HostFloat2D rvm_mc_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rvm_mc_r4", ngrdcol, nzt);
  HostFloat2D P_in_pa_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::P_in_pa_r4", ngrdcol, nzt);
  HostFloat2D rho_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::rho_r4", ngrdcol, nzt);
  HostFloat2D dzq_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::dzq_r4", ngrdcol, nzt);
  HostFloat2D wm_zt_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::wm_zt_r4", ngrdcol,
                                        nzt); // Mean vertical velocity on the thermo grid  [m/s]
  HostFloat2D w_std_dev_r4 =
      workspace_scope.view<HostFloat2D>("morrison_microphys_driver::w_std_dev_r4", ngrdcol,
                                        nzt); // Standard deviation of w                    [m/s]
  int i;                                      // Column and vertical level indices [-]
  int k;                                      // Column and vertical level indices [-]
  HostReal2D Nrm_auto = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::Nrm_auto", ngrdcol,
      nzt); // Change in Nrm due to autoconversion               [num/kg/s]
  HostReal2D Nrm_evap = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::Nrm_evap", ngrdcol,
      nzt); // Change in Nrm due to evaporation                  [num/kg/s]
  HostReal1D rsm_sd_morr_int =
      workspace_scope.view<HostReal1D>("morrison_microphys_driver::rsm_sd_morr_int", ngrdcol);
  HostReal2D hl_before =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::hl_before", ngrdcol, nzt);
  HostReal2D qto_before =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::qto_before", ngrdcol, nzt);
  HostReal2D hl_after =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::hl_after", ngrdcol, nzt);
  HostReal2D hl_on_Cp_residual = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::hl_on_Cp_residual", ngrdcol, nzt);
  HostReal2D qto_after =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::qto_after", ngrdcol, nzt);
  HostReal2D qto_residual =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::qto_residual", ngrdcol, nzt);
  int iirr;
  int iirs;
  int iiri;
  int iirg;
  int iiNr;
  int iiNs;
  int iiNi;
  int iiNg;
  // Update weighted SILHS diagnostics under their source names.
  auto update_microphys_stat = [&](std::string name, auto value) {
    auto stats_scope = workspace.scope();
    std::string output_name;

    if ((!stats.l_sample)) {
      return;
    }

    output_name = name;
    if ((l_latin_hypercube && (microphys.lh_microphys_type == lh_microphys_non_interactive))) {
      if (microphys_detail::trim(name) == std::string("rrm_auto") ||
          microphys_detail::trim(name) == std::string("rrm_accr") ||
          microphys_detail::trim(name) == std::string("rrm_evap") ||
          microphys_detail::trim(name) == std::string("Nrm_auto") ||
          microphys_detail::trim(name) == std::string("Nrm_evap")) {
        output_name = (std::string("lh_") + microphys_detail::trim(name));
      } else {
        return;
      }
    }

    if (l_latin_hypercube) {
      HostReal2D array_actual_1 = stats_scope.view<HostReal2D>(
          Kokkos::WithoutInitializing,
          "morrison_microphys_driver/update_microphys_stat::array_actual_1", value.extent_int(0),
          value.extent_int(1));
      for (int index_1 = 0; index_1 < value.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < value.extent_int(0); ++index_0) {
          array_actual_1(index_0, index_1) =
              (value(index_0, index_1) * sample_weight(index_0, index_1));
        }
      }
      stats_update(microphys_detail::trim(output_name), array_actual_1, stats, true);
    } else {
      stats_update(microphys_detail::trim(output_name), value, stats);
    }
  };
  // Update the source surface precipitation diagnostic for this sample.
  auto update_microphys_stat_sfc = [&](std::string name, auto value) {
    auto stats_scope = workspace.scope();
    if (((!stats.l_sample) ||
         (l_latin_hypercube && (microphys.lh_microphys_type == lh_microphys_non_interactive)))) {
      return;
    }

    if (l_latin_hypercube) {
      HostReal1D array_actual_1 = stats_scope.view<HostReal1D>(
          Kokkos::WithoutInitializing,
          "morrison_microphys_driver/update_microphys_stat_sfc::array_actual_1",
          value.extent_int(0));
      for (int index_0 = 0; index_0 < value.extent_int(0); ++index_0) {
        array_actual_1(index_0) = (value(index_0) * sample_weight(index_0, 0));
      }
      stats_update(name, array_actual_1, stats, std::nullopt, true);
    } else {
      stats_update(name, value, stats);
    }
  };
  // Print the complete source input/process fields for a nonfinite result.
  auto print_morr_error_output = [&]() {
    bool nonfinite_at_lev;
    int i;
    int k;

    //--------------------------- Begin Code ---------------------------

    // Check for a non-finite value in any of the major microphysics tendency variables
    // that are output from Morrison microphysics for use in the code.
    auto diagnostic_scope    = workspace.scope();
    const auto any_nonfinite = [](auto field) {
      for (int k = 0; k < field.extent_int(1); ++k)
        for (int i = 0; i < field.extent_int(0); ++i)
          if (!std::isfinite(static_cast<real_t>(field(i, k))))
            return true;
      return false;
    };
    if (any_nonfinite(rcm_mc_r4) || any_nonfinite(rim_mc_r4) || any_nonfinite(rsm_mc_r4) ||
        any_nonfinite(rrm_mc_r4) || any_nonfinite(rgm_mc_r4) || any_nonfinite(Ncm_mc_r4) ||
        any_nonfinite(Nim_mc_r4) || any_nonfinite(Nsm_mc_r4) || any_nonfinite(Nrm_mc_r4) ||
        any_nonfinite(Ngm_mc_r4) || any_nonfinite(rvm_mc_r4) || any_nonfinite(T_in_K_mc)) {
      const auto diagnostic_gr_zt         = to_host_space(workspace, gr_zt);
      const auto diagnostic_thlm          = to_host_space(workspace, thlm);
      const auto diagnostic_rvm           = to_host_space(workspace, rvm);
      const auto diagnostic_rcm           = to_host_space(workspace, rcm);
      const auto diagnostic_hydromet      = to_host_space(workspace, hydromet);
      const auto diagnostic_Ncm           = to_host_space(workspace, Ncm);
      const auto diagnostic_rgm_sten      = to_host_space(workspace, rgm_sten);
      const auto diagnostic_rrm_sten      = to_host_space(workspace, rrm_sten);
      const auto diagnostic_rim_sten      = to_host_space(workspace, rim_sten);
      const auto diagnostic_rsm_sten      = to_host_space(workspace, rsm_sten);
      const auto diagnostic_rcm_sten      = to_host_space(workspace, rcm_sten);
      const auto diagnostic_NGSTEN        = to_host_space(workspace, NGSTEN);
      const auto diagnostic_NRSTEN        = to_host_space(workspace, NRSTEN);
      const auto diagnostic_NISTEN        = to_host_space(workspace, NISTEN);
      const auto diagnostic_NSSTEN        = to_host_space(workspace, NSSTEN);
      const auto diagnostic_NCSTEN        = to_host_space(workspace, NCSTEN);
      const auto diagnostic_cloud_frac_in = to_host_space(workspace, cloud_frac_in);
      const auto diagnostic_PRC           = to_host_space(workspace, PRC);
      const auto diagnostic_PRA           = to_host_space(workspace, PRA);
      const auto diagnostic_PSMLT         = to_host_space(workspace, PSMLT);
      const auto diagnostic_EVPMS         = to_host_space(workspace, EVPMS);
      const auto diagnostic_PRACS         = to_host_space(workspace, PRACS);
      const auto diagnostic_EVPMG         = to_host_space(workspace, EVPMG);
      const auto diagnostic_PRACG         = to_host_space(workspace, PRACG);
      const auto diagnostic_PRE           = to_host_space(workspace, PRE);
      const auto diagnostic_PGMLT         = to_host_space(workspace, PGMLT);
      const auto diagnostic_MNUCCC        = to_host_space(workspace, MNUCCC);
      const auto diagnostic_PSACWS        = to_host_space(workspace, PSACWS);
      const auto diagnostic_PSACWI        = to_host_space(workspace, PSACWI);
      const auto diagnostic_QMULTS        = to_host_space(workspace, QMULTS);
      const auto diagnostic_QMULTG        = to_host_space(workspace, QMULTG);
      const auto diagnostic_PSACWG        = to_host_space(workspace, PSACWG);
      const auto diagnostic_PGSACW        = to_host_space(workspace, PGSACW);
      const auto diagnostic_PRD           = to_host_space(workspace, PRD);
      const auto diagnostic_PRCI          = to_host_space(workspace, PRCI);
      const auto diagnostic_PRAI          = to_host_space(workspace, PRAI);
      const auto diagnostic_QMULTR        = to_host_space(workspace, QMULTR);
      const auto diagnostic_QMULTRG       = to_host_space(workspace, QMULTRG);
      const auto diagnostic_MNUCCD        = to_host_space(workspace, MNUCCD);
      const auto diagnostic_PRACI         = to_host_space(workspace, PRACI);
      const auto diagnostic_PRACIS        = to_host_space(workspace, PRACIS);
      const auto diagnostic_EPRD          = to_host_space(workspace, EPRD);
      const auto diagnostic_MNUCCR        = to_host_space(workspace, MNUCCR);
      const auto diagnostic_PIACR         = to_host_space(workspace, PIACR);
      const auto diagnostic_PIACRS        = to_host_space(workspace, PIACRS);
      const auto diagnostic_PGRACS        = to_host_space(workspace, PGRACS);
      const auto diagnostic_PRDS          = to_host_space(workspace, PRDS);
      const auto diagnostic_EPRDS         = to_host_space(workspace, EPRDS);
      const auto diagnostic_PSACR         = to_host_space(workspace, PSACR);
      const auto diagnostic_PRDG          = to_host_space(workspace, PRDG);
      const auto diagnostic_EPRDG         = to_host_space(workspace, EPRDG);
      const auto diagnostic_NPRC1         = to_host_space(workspace, NPRC1);
      const auto diagnostic_NRAGG         = to_host_space(workspace, NRAGG);
      const auto diagnostic_NPRACG        = to_host_space(workspace, NPRACG);
      const auto diagnostic_NSUBR         = to_host_space(workspace, NSUBR);
      const auto diagnostic_NSMLTR        = to_host_space(workspace, NSMLTR);
      const auto diagnostic_NGMLTR        = to_host_space(workspace, NGMLTR);
      const auto diagnostic_NPRACS        = to_host_space(workspace, NPRACS);
      const auto diagnostic_NNUCCR        = to_host_space(workspace, NNUCCR);
      const auto diagnostic_NIACR         = to_host_space(workspace, NIACR);
      const auto diagnostic_NIACRS        = to_host_space(workspace, NIACRS);
      const auto diagnostic_NGRACS        = to_host_space(workspace, NGRACS);
      const auto diagnostic_NSMLTS        = to_host_space(workspace, NSMLTS);
      const auto diagnostic_NSAGG         = to_host_space(workspace, NSAGG);
      const auto diagnostic_NPRCI         = to_host_space(workspace, NPRCI);
      const auto diagnostic_NSCNG         = to_host_space(workspace, NSCNG);
      const auto diagnostic_NSUBS         = to_host_space(workspace, NSUBS);
      const auto diagnostic_PCC           = to_host_space(workspace, PCC);
      const auto diagnostic_NNUCCC        = to_host_space(workspace, NNUCCC);
      const auto diagnostic_NPSACWS       = to_host_space(workspace, NPSACWS);
      const auto diagnostic_NPRA          = to_host_space(workspace, NPRA);
      const auto diagnostic_NPRC          = to_host_space(workspace, NPRC);
      const auto diagnostic_NPSACWI       = to_host_space(workspace, NPSACWI);
      const auto diagnostic_NPSACWG       = to_host_space(workspace, NPSACWG);
      const auto diagnostic_NPRAI         = to_host_space(workspace, NPRAI);
      const auto diagnostic_NMULTS        = to_host_space(workspace, NMULTS);
      const auto diagnostic_NMULTG        = to_host_space(workspace, NMULTG);
      const auto diagnostic_NMULTR        = to_host_space(workspace, NMULTR);
      const auto diagnostic_NMULTRG       = to_host_space(workspace, NMULTRG);
      const auto diagnostic_NNUCCD        = to_host_space(workspace, NNUCCD);
      const auto diagnostic_NSUBI         = to_host_space(workspace, NSUBI);
      const auto diagnostic_NGMLTG        = to_host_space(workspace, NGMLTG);
      const auto diagnostic_NSUBG         = to_host_space(workspace, NSUBG);
      const auto diagnostic_NACT          = to_host_space(workspace, NACT);
      const auto diagnostic_SIZEFIX_NR    = to_host_space(workspace, SIZEFIX_NR);
      const auto diagnostic_SIZEFIX_NC    = to_host_space(workspace, SIZEFIX_NC);
      const auto diagnostic_SIZEFIX_NI    = to_host_space(workspace, SIZEFIX_NI);
      const auto diagnostic_SIZEFIX_NS    = to_host_space(workspace, SIZEFIX_NS);
      const auto diagnostic_SIZEFIX_NG    = to_host_space(workspace, SIZEFIX_NG);
      const auto diagnostic_NEGFIX_NR     = to_host_space(workspace, NEGFIX_NR);
      const auto diagnostic_NEGFIX_NC     = to_host_space(workspace, NEGFIX_NC);
      const auto diagnostic_NEGFIX_NI     = to_host_space(workspace, NEGFIX_NI);
      const auto diagnostic_NEGFIX_NS     = to_host_space(workspace, NEGFIX_NS);
      const auto diagnostic_NEGFIX_NG     = to_host_space(workspace, NEGFIX_NG);
      const auto diagnostic_NIM_MORR_CL   = to_host_space(workspace, NIM_MORR_CL);
      const auto diagnostic_QC_INST       = to_host_space(workspace, QC_INST);
      const auto diagnostic_QR_INST       = to_host_space(workspace, QR_INST);
      const auto diagnostic_QI_INST       = to_host_space(workspace, QI_INST);
      const auto diagnostic_QS_INST       = to_host_space(workspace, QS_INST);
      const auto diagnostic_QG_INST       = to_host_space(workspace, QG_INST);
      const auto diagnostic_NC_INST       = to_host_space(workspace, NC_INST);
      const auto diagnostic_NR_INST       = to_host_space(workspace, NR_INST);
      const auto diagnostic_NI_INST       = to_host_space(workspace, NI_INST);
      const auto diagnostic_NS_INST       = to_host_space(workspace, NS_INST);
      const auto diagnostic_NG_INST       = to_host_space(workspace, NG_INST);
      std::cerr << std::string("non-finite detected in a Morrison microphysics tendency") << "\n";
      for (k = 0; k < nzt; k += 1) {
        for (i = 0; i < ngrdcol; i += 1) {
          nonfinite_at_lev = false;
          if ((!std::isfinite(static_cast<real_t>(rcm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rcm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(rim_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rim_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(rsm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rsm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(rrm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rrm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(rgm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rgm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(Ncm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in Ncm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(Nim_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in Nim_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(Nsm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in Nsm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(Nrm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in Nrm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(Ngm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in Ngm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(rvm_mc_r4(i, k))))) {
            std::cerr << std::string("non-finite detected in rvm_mc_r4 at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if ((!std::isfinite(static_cast<real_t>(T_in_K_mc(i, k))))) {
            std::cerr << std::string("non-finite detected in T_in_K_mc at column = ") << " "
                      << (((i + 1))) << " " << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            nonfinite_at_lev = true;
          }
          if (nonfinite_at_lev) {
            std::cerr << std::string("At column = ") << " " << (((i + 1))) << " "
                      << std::string(", k = ") << " " << (((k + 1))) << " "
                      << std::string("Altitude (m) = ") << " " << diagnostic_gr_zt(i, k) << "\n";
            std::cerr << std::string("thlm (in) = ") << " " << diagnostic_thlm(i, k) << "\n";
            std::cerr << std::string("rvm (in) = ") << " " << diagnostic_rvm(i, k) << "\n";
            std::cerr << std::string("rcm (in) = ") << " " << diagnostic_rcm(i, k) << "\n";
            std::cerr << std::string("rrm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iirr) << "\n";
            std::cerr << std::string("rim (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iiri) << "\n";
            std::cerr << std::string("rsm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iirs) << "\n";
            std::cerr << std::string("rgm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iirg) << "\n";
            std::cerr << std::string("Ncm (in) = ") << " " << diagnostic_Ncm(i, k) << "\n";
            std::cerr << std::string("Nrm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iiNr) << "\n";
            std::cerr << std::string("Nim (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iiNi) << "\n";
            std::cerr << std::string("Nsm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iiNs) << "\n";
            std::cerr << std::string("Ngm (in) = ") << " "
                      << diagnostic_hydromet(i, k, hm_metadata_iiNg) << "\n";
            std::cerr << std::string("rgm_sten = ") << " " << diagnostic_rgm_sten(i, k) << "\n";
            std::cerr << std::string("rrm_sten = ") << " " << diagnostic_rrm_sten(i, k) << "\n";
            std::cerr << std::string("rim_sten = ") << " " << diagnostic_rim_sten(i, k) << "\n";
            std::cerr << std::string("rsm_sten = ") << " " << diagnostic_rsm_sten(i, k) << "\n";
            std::cerr << std::string("rcm_sten = ") << " " << diagnostic_rcm_sten(i, k) << "\n";
            std::cerr << std::string("NGSTEN = ") << " " << diagnostic_NGSTEN(i, k) << "\n";
            std::cerr << std::string("NRSTEN = ") << " " << diagnostic_NRSTEN(i, k) << "\n";
            std::cerr << std::string("NISTEN = ") << " " << diagnostic_NISTEN(i, k) << "\n";
            std::cerr << std::string("NSSTEN = ") << " " << diagnostic_NSSTEN(i, k) << "\n";
            std::cerr << std::string("NCSTEN = ") << " " << diagnostic_NCSTEN(i, k) << "\n";
            std::cerr << std::string("cloud_frac_in = ") << " " << diagnostic_cloud_frac_in(i, k)
                      << "\n";
            std::cerr << std::string("PRC = ") << " " << diagnostic_PRC(i, k) << "\n";
            std::cerr << std::string("PRA = ") << " " << diagnostic_PRA(i, k) << "\n";
            std::cerr << std::string("PSMLT = ") << " " << diagnostic_PSMLT(i, k) << "\n";
            std::cerr << std::string("EVPMS = ") << " " << diagnostic_EVPMS(i, k) << "\n";
            std::cerr << std::string("PRACS = ") << " " << diagnostic_PRACS(i, k) << "\n";
            std::cerr << std::string("EVPMG = ") << " " << diagnostic_EVPMG(i, k) << "\n";
            std::cerr << std::string("PRACG = ") << " " << diagnostic_PRACG(i, k) << "\n";
            std::cerr << std::string("PRE = ") << " " << diagnostic_PRE(i, k) << "\n";
            std::cerr << std::string("PGMLT = ") << " " << diagnostic_PGMLT(i, k) << "\n";
            std::cerr << std::string("MNUCCC = ") << " " << diagnostic_MNUCCC(i, k) << "\n";
            std::cerr << std::string("PSACWS = ") << " " << diagnostic_PSACWS(i, k) << "\n";
            std::cerr << std::string("PSACWI = ") << " " << diagnostic_PSACWI(i, k) << "\n";
            std::cerr << std::string("QMULTS = ") << " " << diagnostic_QMULTS(i, k) << "\n";
            std::cerr << std::string("QMULTG = ") << " " << diagnostic_QMULTG(i, k) << "\n";
            std::cerr << std::string("PSACWG = ") << " " << diagnostic_PSACWG(i, k) << "\n";
            std::cerr << std::string("PGSACW = ") << " " << diagnostic_PGSACW(i, k) << "\n";
            std::cerr << std::string("PRD = ") << " " << diagnostic_PRD(i, k) << "\n";
            std::cerr << std::string("PRCI = ") << " " << diagnostic_PRCI(i, k) << "\n";
            std::cerr << std::string("PRAI = ") << " " << diagnostic_PRAI(i, k) << "\n";
            std::cerr << std::string("QMULTR = ") << " " << diagnostic_QMULTR(i, k) << "\n";
            std::cerr << std::string("QMULTRG = ") << " " << diagnostic_QMULTRG(i, k) << "\n";
            std::cerr << std::string("MNUCCD = ") << " " << diagnostic_MNUCCD(i, k) << "\n";
            std::cerr << std::string("PRACI = ") << " " << diagnostic_PRACI(i, k) << "\n";
            std::cerr << std::string("PRACIS = ") << " " << diagnostic_PRACIS(i, k) << "\n";
            std::cerr << std::string("EPRD = ") << " " << diagnostic_EPRD(i, k) << "\n";
            std::cerr << std::string("MNUCCR = ") << " " << diagnostic_MNUCCR(i, k) << "\n";
            std::cerr << std::string("PIACR = ") << " " << diagnostic_PIACR(i, k) << "\n";
            std::cerr << std::string("PIACRS = ") << " " << diagnostic_PIACRS(i, k) << "\n";
            std::cerr << std::string("PGRACS = ") << " " << diagnostic_PGRACS(i, k) << "\n";
            std::cerr << std::string("PRDS = ") << " " << diagnostic_PRDS(i, k) << "\n";
            std::cerr << std::string("EPRDS = ") << " " << diagnostic_EPRDS(i, k) << "\n";
            std::cerr << std::string("PSACR = ") << " " << diagnostic_PSACR(i, k) << "\n";
            std::cerr << std::string("PRDG = ") << " " << diagnostic_PRDG(i, k) << "\n";
            std::cerr << std::string("EPRDG = ") << " " << diagnostic_EPRDG(i, k) << "\n";
            std::cerr << std::string("NPRC1 = ") << " " << diagnostic_NPRC1(i, k) << "\n";
            std::cerr << std::string("NRAGG = ") << " " << diagnostic_NRAGG(i, k) << "\n";
            std::cerr << std::string("NPRACG = ") << " " << diagnostic_NPRACG(i, k) << "\n";
            std::cerr << std::string("NSUBR = ") << " " << diagnostic_NSUBR(i, k) << "\n";
            std::cerr << std::string("NSMLTR = ") << " " << diagnostic_NSMLTR(i, k) << "\n";
            std::cerr << std::string("NGMLTR = ") << " " << diagnostic_NGMLTR(i, k) << "\n";
            std::cerr << std::string("NPRACS = ") << " " << diagnostic_NPRACS(i, k) << "\n";
            std::cerr << std::string("NNUCCR = ") << " " << diagnostic_NNUCCR(i, k) << "\n";
            std::cerr << std::string("NIACR = ") << " " << diagnostic_NIACR(i, k) << "\n";
            std::cerr << std::string("NIACRS = ") << " " << diagnostic_NIACRS(i, k) << "\n";
            std::cerr << std::string("NGRACS = ") << " " << diagnostic_NGRACS(i, k) << "\n";
            std::cerr << std::string("NSMLTS = ") << " " << diagnostic_NSMLTS(i, k) << "\n";
            std::cerr << std::string("NSAGG = ") << " " << diagnostic_NSAGG(i, k) << "\n";
            std::cerr << std::string("NPRCI = ") << " " << diagnostic_NPRCI(i, k) << "\n";
            std::cerr << std::string("NSCNG = ") << " " << diagnostic_NSCNG(i, k) << "\n";
            std::cerr << std::string("NSUBS = ") << " " << diagnostic_NSUBS(i, k) << "\n";
            std::cerr << std::string("PCC = ") << " " << diagnostic_PCC(i, k) << "\n";
            std::cerr << std::string("NNUCCC = ") << " " << diagnostic_NNUCCC(i, k) << "\n";
            std::cerr << std::string("NPSACWS = ") << " " << diagnostic_NPSACWS(i, k) << "\n";
            std::cerr << std::string("NPRA = ") << " " << diagnostic_NPRA(i, k) << "\n";
            std::cerr << std::string("NPRC = ") << " " << diagnostic_NPRC(i, k) << "\n";
            std::cerr << std::string("NPSACWI = ") << " " << diagnostic_NPSACWI(i, k) << "\n";
            std::cerr << std::string("NPSACWG = ") << " " << diagnostic_NPSACWG(i, k) << "\n";
            std::cerr << std::string("NPRAI = ") << " " << diagnostic_NPRAI(i, k) << "\n";
            std::cerr << std::string("NMULTS = ") << " " << diagnostic_NMULTS(i, k) << "\n";
            std::cerr << std::string("NMULTG = ") << " " << diagnostic_NMULTG(i, k) << "\n";
            std::cerr << std::string("NMULTR = ") << " " << diagnostic_NMULTR(i, k) << "\n";
            std::cerr << std::string("NMULTRG = ") << " " << diagnostic_NMULTRG(i, k) << "\n";
            std::cerr << std::string("NNUCCD = ") << " " << diagnostic_NNUCCD(i, k) << "\n";
            std::cerr << std::string("NSUBI = ") << " " << diagnostic_NSUBI(i, k) << "\n";
            std::cerr << std::string("NGMLTG = ") << " " << diagnostic_NGMLTG(i, k) << "\n";
            std::cerr << std::string("NSUBG = ") << " " << diagnostic_NSUBG(i, k) << "\n";
            std::cerr << std::string("NACT = ") << " " << diagnostic_NACT(i, k) << "\n";
            std::cerr << std::string("SIZEFIX_NR = ") << " " << diagnostic_SIZEFIX_NR(i, k) << "\n";
            std::cerr << std::string("SIZEFIX_NC = ") << " " << diagnostic_SIZEFIX_NC(i, k) << "\n";
            std::cerr << std::string("SIZEFIX_NI = ") << " " << diagnostic_SIZEFIX_NI(i, k) << "\n";
            std::cerr << std::string("SIZEFIX_NS = ") << " " << diagnostic_SIZEFIX_NS(i, k) << "\n";
            std::cerr << std::string("SIZEFIX_NG = ") << " " << diagnostic_SIZEFIX_NG(i, k) << "\n";
            std::cerr << std::string("NEGFIX_NR = ") << " " << diagnostic_NEGFIX_NR(i, k) << "\n";
            std::cerr << std::string("NEGFIX_NC = ") << " " << diagnostic_NEGFIX_NC(i, k) << "\n";
            std::cerr << std::string("NEGFIX_NI = ") << " " << diagnostic_NEGFIX_NI(i, k) << "\n";
            std::cerr << std::string("NEGFIX_NS = ") << " " << diagnostic_NEGFIX_NS(i, k) << "\n";
            std::cerr << std::string("NEGFIX_NG = ") << " " << diagnostic_NEGFIX_NG(i, k) << "\n";
            std::cerr << std::string("NIM_MORR_CL = ") << " " << diagnostic_NIM_MORR_CL(i, k)
                      << "\n";
            std::cerr << std::string("QC_INST = ") << " " << diagnostic_QC_INST(i, k) << "\n";
            std::cerr << std::string("QR_INST = ") << " " << diagnostic_QR_INST(i, k) << "\n";
            std::cerr << std::string("QI_INST = ") << " " << diagnostic_QI_INST(i, k) << "\n";
            std::cerr << std::string("QS_INST = ") << " " << diagnostic_QS_INST(i, k) << "\n";
            std::cerr << std::string("QG_INST = ") << " " << diagnostic_QG_INST(i, k) << "\n";
            std::cerr << std::string("NC_INST = ") << " " << diagnostic_NC_INST(i, k) << "\n";
            std::cerr << std::string("NR_INST = ") << " " << diagnostic_NR_INST(i, k) << "\n";
            std::cerr << std::string("NI_INST = ") << " " << diagnostic_NI_INST(i, k) << "\n";
            std::cerr << std::string("NS_INST = ") << " " << diagnostic_NS_INST(i, k) << "\n";
            std::cerr << std::string("NG_INST = ") << " " << diagnostic_NG_INST(i, k) << "\n";
            std::cerr << std::string("---------------------------------------------") << "\n";
            // nonfinite_at_lev
          }
          // i = 1, ngrdcol
        }
        // k = 1, nzt, 1
      }
    }

    return;
  };

  // ---- Begin Code ----

  iirr = hm_metadata_iirr;
  iirs = hm_metadata_iirs;
  iiri = hm_metadata_iiri;
  iirg = hm_metadata_iirg;
  iiNr = hm_metadata_iiNr;
  iiNs = hm_metadata_iiNs;
  iiNi = hm_metadata_iiNi;
  iiNg = hm_metadata_iiNg;

  // Determine temperature
  for (int index_1 = 0; index_1 < T_in_K_core.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < T_in_K_core.extent_int(0); ++index_0) {
      T_in_K_core(index_0, index_1) =
          thlm2T_in_K_api(thlm(index_0, index_1), exner(index_0, index_1), rcm(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < T_in_K.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < T_in_K.extent_int(0); ++index_0) {
      T_in_K(index_0, index_1) = static_cast<float>(T_in_K_core(index_0, index_1));
    }
  }

  if (l_latin_hypercube) {
    // Don't use sgs cloud fraction to weight the tendencies
    for (int index_1 = 0; index_1 < cloud_frac_in.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < cloud_frac_in.extent_int(0); ++index_0) {
        cloud_frac_in(index_0, index_1) = 0.0f;
      }
    }

    // Impose a minimum value on w
    for (int index_1 = 0; index_1 < wm_zt_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < wm_zt_r4.extent_int(0); ++index_0) {
        wm_zt_r4(index_0, index_1) = static_cast<float>(
            ((wm_zt(index_0, index_1)) > (w_thresh) ? (wm_zt(index_0, index_1)) : (w_thresh)));
      }
    }
    // Don't add in a standard deviation for aerosol activation
    for (int index_1 = 0; index_1 < w_std_dev_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < w_std_dev_r4.extent_int(0); ++index_0) {
        w_std_dev_r4(index_0, index_1) = 0.f;
      }
    }

  } else {
    // Use sgs cloud fraction to weight tendencies
    for (int index_1 = 0; index_1 < cloud_frac_in.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < cloud_frac_in.extent_int(0); ++index_0) {
        cloud_frac_in(index_0, index_1) = static_cast<float>(cloud_frac(index_0, index_1));
      }
    }

    // Use the mean value without a threshold
    for (int index_1 = 0; index_1 < wm_zt_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < wm_zt_r4.extent_int(0); ++index_0) {
        wm_zt_r4(index_0, index_1) = static_cast<float>(wm_zt(index_0, index_1));
      }
    }
    // Add in a standard deviation
    for (int index_1 = 0; index_1 < w_std_dev_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < w_std_dev_r4.extent_int(0); ++index_0) {
        w_std_dev_r4(index_0, index_1) = static_cast<float>(w_std_dev(index_0, index_1));
      }
    }
  }

  for (int index_1 = 0; index_1 < rcm_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_r4.extent_int(0); ++index_0) {
      rcm_r4(index_0, index_1) = static_cast<float>(rcm(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rvm_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rvm_r4.extent_int(0); ++index_0) {
      rvm_r4(index_0, index_1) = static_cast<float>(rvm(index_0, index_1));
    }
  }

  // Note: The Ncm_r4 variable is only used if INUM = 0 in the Morrison code;
  // otherwise the NDCNST variable is used as the fixed value.
  for (int index_1 = 0; index_1 < Ncm_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_r4.extent_int(0); ++index_0) {
      Ncm_r4(index_0, index_1) = static_cast<float>(Ncm(index_0, index_1));
    }
  }

  for (int index_2 = 0; index_2 < hydromet_r4.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
        hydromet_r4(index_0, index_1, index_2) =
            static_cast<float>(hydromet(index_0, index_1, index_2));
      }
    }
  }

  if (l_evaporate_cold_rcm) {
    // Convert liquid to vapor at temperatures colder than -37C
    for (k = 0; k < nzt; k += 1) {
      for (i = 0; i < ngrdcol; i += 1) {
        if ((T_in_K(i, k) < 236.15f)) {
          rcm_r4(i, k)        = 0.0f;
          cloud_frac_in(i, k) = 0.0f;
          Ncm_r4(i, k)        = 0.0f;
        }
      }
    }
  }

  // Initialize tendencies to zero
  for (int index_1 = 0; index_1 < T_in_K_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < T_in_K_mc.extent_int(0); ++index_0) {
      T_in_K_mc(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < rcm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_mc.extent_int(0); ++index_0) {
      rcm_mc(index_0, index_1) = real_t(0.0);
    }
  }
  for (int index_1 = 0; index_1 < rvm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rvm_mc.extent_int(0); ++index_0) {
      rvm_mc(index_0, index_1) = real_t(0.0);
    }
  }
  for (int index_2 = 0; index_2 < hydromet_mc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
        hydromet_mc(index_0, index_1, index_2) = real_t(0.0);
      }
    }
  }
  for (int index_1 = 0; index_1 < Ncm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_mc.extent_int(0); ++index_0) {
      Ncm_mc(index_0, index_1) = real_t(0.0);
    }
  }

  // Initialize process rates not guaranteed to be set on Morrison's early-exit path.
  for (int index_1 = 0; index_1 < MNUCCC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCC.extent_int(0); ++index_0) {
      MNUCCC(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PSACWS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWS.extent_int(0); ++index_0) {
      PSACWS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PSACWI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWI.extent_int(0); ++index_0) {
      PSACWI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < QMULTS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTS.extent_int(0); ++index_0) {
      QMULTS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < QMULTG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTG.extent_int(0); ++index_0) {
      QMULTG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PSACWG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWG.extent_int(0); ++index_0) {
      PSACWG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PGSACW.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PGSACW.extent_int(0); ++index_0) {
      PGSACW(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRD.extent_int(0); ++index_0) {
      PRD(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRCI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRCI.extent_int(0); ++index_0) {
      PRCI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRAI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRAI.extent_int(0); ++index_0) {
      PRAI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < QMULTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTR.extent_int(0); ++index_0) {
      QMULTR(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < QMULTRG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTRG.extent_int(0); ++index_0) {
      QMULTRG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < MNUCCD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCD.extent_int(0); ++index_0) {
      MNUCCD(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRACI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACI.extent_int(0); ++index_0) {
      PRACI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRACIS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACIS.extent_int(0); ++index_0) {
      PRACIS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < EPRD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRD.extent_int(0); ++index_0) {
      EPRD(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < MNUCCR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCR.extent_int(0); ++index_0) {
      MNUCCR(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PIACR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PIACR.extent_int(0); ++index_0) {
      PIACR(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PIACRS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PIACRS.extent_int(0); ++index_0) {
      PIACRS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PGRACS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PGRACS.extent_int(0); ++index_0) {
      PGRACS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRDS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRDS.extent_int(0); ++index_0) {
      PRDS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < EPRDS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRDS.extent_int(0); ++index_0) {
      EPRDS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PSACR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACR.extent_int(0); ++index_0) {
      PSACR(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < PRDG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRDG.extent_int(0); ++index_0) {
      PRDG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < EPRDG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRDG.extent_int(0); ++index_0) {
      EPRDG(index_0, index_1) = 0.0f;
    }
  }

  // Initialize these Morrison budgets to zero too
  for (int index_1 = 0; index_1 < NSAGG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSAGG.extent_int(0); ++index_0) {
      NSAGG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NPRCI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRCI.extent_int(0); ++index_0) {
      NPRCI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NSCNG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSCNG.extent_int(0); ++index_0) {
      NSCNG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NSUBS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBS.extent_int(0); ++index_0) {
      NSUBS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NNUCCC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NNUCCC.extent_int(0); ++index_0) {
      NNUCCC(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NPSACWS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWS.extent_int(0); ++index_0) {
      NPSACWS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NPSACWI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWI.extent_int(0); ++index_0) {
      NPSACWI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NPSACWG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWG.extent_int(0); ++index_0) {
      NPSACWG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NPRAI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRAI.extent_int(0); ++index_0) {
      NPRAI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NMULTS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTS.extent_int(0); ++index_0) {
      NMULTS(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NMULTG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTG.extent_int(0); ++index_0) {
      NMULTG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NMULTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTR.extent_int(0); ++index_0) {
      NMULTR(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NMULTRG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTRG.extent_int(0); ++index_0) {
      NMULTRG(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NNUCCD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NNUCCD.extent_int(0); ++index_0) {
      NNUCCD(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NSUBI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBI.extent_int(0); ++index_0) {
      NSUBI(index_0, index_1) = 0.0f;
    }
  }
  for (int index_1 = 0; index_1 < NSUBG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBG.extent_int(0); ++index_0) {
      NSUBG(index_0, index_1) = 0.0f;
    }
  }

  for (int index_2 = 0; index_2 < hydromet_mc_r4.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
        hydromet_mc_r4(index_0, index_1, index_2) =
            static_cast<float>(hydromet_mc(index_0, index_1, index_2));
      }
    }
  }
  for (int index_1 = 0; index_1 < Ncm_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_mc_r4.extent_int(0); ++index_0) {
      Ncm_mc_r4(index_0, index_1) = static_cast<float>(Ncm_mc(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rcm_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_mc_r4.extent_int(0); ++index_0) {
      rcm_mc_r4(index_0, index_1) = static_cast<float>(rcm_mc(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rvm_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rvm_mc_r4.extent_int(0); ++index_0) {
      rvm_mc_r4(index_0, index_1) = static_cast<float>(rvm_mc(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < P_in_pa_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < P_in_pa_r4.extent_int(0); ++index_0) {
      P_in_pa_r4(index_0, index_1) = static_cast<float>(p_in_Pa(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rho_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rho_r4.extent_int(0); ++index_0) {
      rho_r4(index_0, index_1) = static_cast<float>(rho(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < dzq_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < dzq_r4.extent_int(0); ++index_0) {
      dzq_r4(index_0, index_1) = static_cast<float>(dzq(index_0, index_1));
    }
  }

  // Unpack hydrometeor arrays.
  for (int index_1 = 0; index_1 < rrm.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm.extent_int(0); ++index_0) {
      rrm(index_0, index_1) = hydromet(index_0, index_1, iirr);
    }
  }

  for (int index_1 = 0; index_1 < rrm_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_r4.extent_int(0); ++index_0) {
      rrm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iirr);
    }
  }
  for (int index_1 = 0; index_1 < Nrm_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_r4.extent_int(0); ++index_0) {
      Nrm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iiNr);
    }
  }

  for (int index_1 = 0; index_1 < rrm_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_mc_r4.extent_int(0); ++index_0) {
      rrm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iirr);
    }
  }
  for (int index_1 = 0; index_1 < Nrm_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_mc_r4.extent_int(0); ++index_0) {
      Nrm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iiNr);
    }
  }

  if (microphys.l_ice_microphys) {

    for (int index_1 = 0; index_1 < rim.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim.extent_int(0); ++index_0) {
        rim(index_0, index_1) = hydromet(index_0, index_1, iiri);
      }
    }
    for (int index_1 = 0; index_1 < rsm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm.extent_int(0); ++index_0) {
        rsm(index_0, index_1) = hydromet(index_0, index_1, iirs);
      }
    }

    for (int index_1 = 0; index_1 < rim_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim_r4.extent_int(0); ++index_0) {
        rim_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iiri);
      }
    }
    for (int index_1 = 0; index_1 < Nim_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nim_r4.extent_int(0); ++index_0) {
        Nim_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iiNi);
      }
    }
    for (int index_1 = 0; index_1 < rsm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm_r4.extent_int(0); ++index_0) {
        rsm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iirs);
      }
    }
    for (int index_1 = 0; index_1 < Nsm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nsm_r4.extent_int(0); ++index_0) {
        Nsm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iiNs);
      }
    }

    for (int index_1 = 0; index_1 < rim_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim_mc_r4.extent_int(0); ++index_0) {
        rim_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iiri);
      }
    }
    for (int index_1 = 0; index_1 < Nim_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nim_mc_r4.extent_int(0); ++index_0) {
        Nim_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iiNi);
      }
    }
    for (int index_1 = 0; index_1 < rsm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm_mc_r4.extent_int(0); ++index_0) {
        rsm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iirs);
      }
    }
    for (int index_1 = 0; index_1 < Nsm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nsm_mc_r4.extent_int(0); ++index_0) {
        Nsm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iiNs);
      }
    }

    if (microphys.l_graupel) {

      for (int index_1 = 0; index_1 < rgm.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm.extent_int(0); ++index_0) {
          rgm(index_0, index_1) = hydromet(index_0, index_1, iirg);
        }
      }

      for (int index_1 = 0; index_1 < rgm_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm_r4.extent_int(0); ++index_0) {
          rgm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iirg);
        }
      }
      for (int index_1 = 0; index_1 < Ngm_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < Ngm_r4.extent_int(0); ++index_0) {
          Ngm_r4(index_0, index_1) = hydromet_r4(index_0, index_1, iiNg);
        }
      }

      for (int index_1 = 0; index_1 < rgm_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm_mc_r4.extent_int(0); ++index_0) {
          rgm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iirg);
        }
      }
      for (int index_1 = 0; index_1 < Ngm_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < Ngm_mc_r4.extent_int(0); ++index_0) {
          Ngm_mc_r4(index_0, index_1) = hydromet_mc_r4(index_0, index_1, iiNg);
        }
      }

    } else {
      // l_graupel disabled

      for (int index_1 = 0; index_1 < rgm.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm.extent_int(0); ++index_0) {
          rgm(index_0, index_1) = zero;
        }
      }

      for (int index_1 = 0; index_1 < rgm_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm_r4.extent_int(0); ++index_0) {
          rgm_r4(index_0, index_1) = 0.0f;
        }
      }
      for (int index_1 = 0; index_1 < Ngm_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < Ngm_r4.extent_int(0); ++index_0) {
          Ngm_r4(index_0, index_1) = 0.0f;
        }
      }

      for (int index_1 = 0; index_1 < rgm_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < rgm_mc_r4.extent_int(0); ++index_0) {
          rgm_mc_r4(index_0, index_1) = 0.0f;
        }
      }
      for (int index_1 = 0; index_1 < Ngm_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < Ngm_mc_r4.extent_int(0); ++index_0) {
          Ngm_mc_r4(index_0, index_1) = 0.0f;
        }
      }

      // l_graupel
    }

  } else {
    // l_ice_microphys disabled

    for (int index_1 = 0; index_1 < rim.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim.extent_int(0); ++index_0) {
        rim(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < rsm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm.extent_int(0); ++index_0) {
        rsm(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < rgm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rgm.extent_int(0); ++index_0) {
        rgm(index_0, index_1) = zero;
      }
    }

    for (int index_1 = 0; index_1 < rim_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim_r4.extent_int(0); ++index_0) {
        rim_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Nim_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nim_r4.extent_int(0); ++index_0) {
        Nim_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < rsm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm_r4.extent_int(0); ++index_0) {
        rsm_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Nsm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nsm_r4.extent_int(0); ++index_0) {
        Nsm_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < rgm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rgm_r4.extent_int(0); ++index_0) {
        rgm_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Ngm_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Ngm_r4.extent_int(0); ++index_0) {
        Ngm_r4(index_0, index_1) = 0.0f;
      }
    }

    for (int index_1 = 0; index_1 < rim_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rim_mc_r4.extent_int(0); ++index_0) {
        rim_mc_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Nim_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nim_mc_r4.extent_int(0); ++index_0) {
        Nim_mc_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < rsm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rsm_mc_r4.extent_int(0); ++index_0) {
        rsm_mc_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Nsm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Nsm_mc_r4.extent_int(0); ++index_0) {
        Nsm_mc_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < rgm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < rgm_mc_r4.extent_int(0); ++index_0) {
        rgm_mc_r4(index_0, index_1) = 0.0f;
      }
    }
    for (int index_1 = 0; index_1 < Ngm_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < Ngm_mc_r4.extent_int(0); ++index_0) {
        Ngm_mc_r4(index_0, index_1) = 0.0f;
      }
    }

    // l_ice_microphys
  }

  for (k = 0; k < nzt; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      // BFB-032: Native rounds grav*z before fusing Cp*T into the energy sum.
      hl_before(i, k) = ((Kokkos::fma(Cp, static_cast<real_t>(T_in_K(i, k)), grav * gr_zt(i, k)) -
                          (Lv * (rcm(i, k) + rrm(i, k)))) -
                         (Ls * ((rim(i, k) + rsm(i, k)) + rgm(i, k))));
    }
  }

  for (int index_1 = 0; index_1 < qto_before.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < qto_before.extent_int(0); ++index_0) {
      qto_before(index_0, index_1) =
          (((((rvm(index_0, index_1) + rcm(index_0, index_1)) + rrm(index_0, index_1)) +
             rim(index_0, index_1)) +
            rsm(index_0, index_1)) +
           rgm(index_0, index_1));
    }
  }

  // Call the one-column Morrison microphysics core.
  for (i = 0; i < ngrdcol; i += 1) {
    // Call the Morrison microphysics
    morrison_engine::M2005MICRO_GRAUPEL(
        workspace, Kokkos::subview(rcm_mc_r4, i, Kokkos::ALL),
        Kokkos::subview(rim_mc_r4, i, Kokkos::ALL), Kokkos::subview(rsm_mc_r4, i, Kokkos::ALL),
        Kokkos::subview(rrm_mc_r4, i, Kokkos::ALL), Kokkos::subview(Ncm_mc_r4, i, Kokkos::ALL),
        Kokkos::subview(Nim_mc_r4, i, Kokkos::ALL), Kokkos::subview(Nsm_mc_r4, i, Kokkos::ALL),
        Kokkos::subview(Nrm_mc_r4, i, Kokkos::ALL), Kokkos::subview(rcm_r4, i, Kokkos::ALL),
        Kokkos::subview(rim_r4, i, Kokkos::ALL), Kokkos::subview(rsm_r4, i, Kokkos::ALL),
        Kokkos::subview(rrm_r4, i, Kokkos::ALL), Kokkos::subview(Ncm_r4, i, Kokkos::ALL),
        Kokkos::subview(Nim_r4, i, Kokkos::ALL), Kokkos::subview(Nsm_r4, i, Kokkos::ALL),
        Kokkos::subview(Nrm_r4, i, Kokkos::ALL), Kokkos::subview(T_in_K_mc, i, Kokkos::ALL),
        Kokkos::subview(rvm_mc_r4, i, Kokkos::ALL), Kokkos::subview(T_in_K, i, Kokkos::ALL),
        Kokkos::subview(rvm_r4, i, Kokkos::ALL), Kokkos::subview(P_in_pa_r4, i, Kokkos::ALL),
        Kokkos::subview(rho_r4, i, Kokkos::ALL), Kokkos::subview(dzq_r4, i, Kokkos::ALL),
        Kokkos::subview(wm_zt_r4, i, Kokkos::ALL), Kokkos::subview(w_std_dev_r4, i, Kokkos::ALL),
        Kokkos::subview(morr_rain_vel_r4, i, Kokkos::ALL), Morr_precip_rate(i), Morr_snow_rate(i),
        Kokkos::subview(effc, i, Kokkos::ALL), Kokkos::subview(effi, i, Kokkos::ALL),
        Kokkos::subview(effs, i, Kokkos::ALL), Kokkos::subview(effr, i, Kokkos::ALL),
        static_cast<float>(dt), 0, 0, 0, 0, 0, (nzt - 1), 0, 0, 0, 0, 0, (nzt - 1),
        Kokkos::subview(rgm_mc_r4, i, Kokkos::ALL), Kokkos::subview(Ngm_mc_r4, i, Kokkos::ALL),
        Kokkos::subview(rgm_r4, i, Kokkos::ALL), Kokkos::subview(Ngm_r4, i, Kokkos::ALL),
        Kokkos::subview(effg, i, Kokkos::ALL), Kokkos::subview(rgm_sten, i, Kokkos::ALL),
        Kokkos::subview(rrm_sten, i, Kokkos::ALL), Kokkos::subview(rim_sten, i, Kokkos::ALL),
        Kokkos::subview(rsm_sten, i, Kokkos::ALL), Kokkos::subview(rcm_sten, i, Kokkos::ALL),
        Kokkos::subview(NGSTEN, i, Kokkos::ALL), Kokkos::subview(NRSTEN, i, Kokkos::ALL),
        Kokkos::subview(NISTEN, i, Kokkos::ALL), Kokkos::subview(NSSTEN, i, Kokkos::ALL),
        Kokkos::subview(NCSTEN, i, Kokkos::ALL), Kokkos::subview(cloud_frac_in, i, Kokkos::ALL),
        Kokkos::subview(PRC, i, Kokkos::ALL), Kokkos::subview(PRA, i, Kokkos::ALL),
        Kokkos::subview(PSMLT, i, Kokkos::ALL), Kokkos::subview(EVPMS, i, Kokkos::ALL),
        Kokkos::subview(PRACS, i, Kokkos::ALL), Kokkos::subview(EVPMG, i, Kokkos::ALL),
        Kokkos::subview(PRACG, i, Kokkos::ALL), Kokkos::subview(PRE, i, Kokkos::ALL),
        Kokkos::subview(PGMLT, i, Kokkos::ALL), Kokkos::subview(MNUCCC, i, Kokkos::ALL),
        Kokkos::subview(PSACWS, i, Kokkos::ALL), Kokkos::subview(PSACWI, i, Kokkos::ALL),
        Kokkos::subview(QMULTS, i, Kokkos::ALL), Kokkos::subview(QMULTG, i, Kokkos::ALL),
        Kokkos::subview(PSACWG, i, Kokkos::ALL), Kokkos::subview(PGSACW, i, Kokkos::ALL),
        Kokkos::subview(PRD, i, Kokkos::ALL), Kokkos::subview(PRCI, i, Kokkos::ALL),
        Kokkos::subview(PRAI, i, Kokkos::ALL), Kokkos::subview(QMULTR, i, Kokkos::ALL),
        Kokkos::subview(QMULTRG, i, Kokkos::ALL), Kokkos::subview(MNUCCD, i, Kokkos::ALL),
        Kokkos::subview(PRACI, i, Kokkos::ALL), Kokkos::subview(PRACIS, i, Kokkos::ALL),
        Kokkos::subview(EPRD, i, Kokkos::ALL), Kokkos::subview(MNUCCR, i, Kokkos::ALL),
        Kokkos::subview(PIACR, i, Kokkos::ALL), Kokkos::subview(PIACRS, i, Kokkos::ALL),
        Kokkos::subview(PGRACS, i, Kokkos::ALL), Kokkos::subview(PRDS, i, Kokkos::ALL),
        Kokkos::subview(EPRDS, i, Kokkos::ALL), Kokkos::subview(PSACR, i, Kokkos::ALL),
        Kokkos::subview(PRDG, i, Kokkos::ALL), Kokkos::subview(EPRDG, i, Kokkos::ALL),
        Kokkos::subview(NPRC1, i, Kokkos::ALL), Kokkos::subview(NRAGG, i, Kokkos::ALL),
        Kokkos::subview(NPRACG, i, Kokkos::ALL), Kokkos::subview(NSUBR, i, Kokkos::ALL),
        Kokkos::subview(NSMLTR, i, Kokkos::ALL), Kokkos::subview(NGMLTR, i, Kokkos::ALL),
        Kokkos::subview(NPRACS, i, Kokkos::ALL), Kokkos::subview(NNUCCR, i, Kokkos::ALL),
        Kokkos::subview(NIACR, i, Kokkos::ALL), Kokkos::subview(NIACRS, i, Kokkos::ALL),
        Kokkos::subview(NGRACS, i, Kokkos::ALL), Kokkos::subview(NSMLTS, i, Kokkos::ALL),
        Kokkos::subview(NSAGG, i, Kokkos::ALL), Kokkos::subview(NPRCI, i, Kokkos::ALL),
        Kokkos::subview(NSCNG, i, Kokkos::ALL), Kokkos::subview(NSUBS, i, Kokkos::ALL),
        Kokkos::subview(PCC, i, Kokkos::ALL), Kokkos::subview(NNUCCC, i, Kokkos::ALL),
        Kokkos::subview(NPSACWS, i, Kokkos::ALL), Kokkos::subview(NPRA, i, Kokkos::ALL),
        Kokkos::subview(NPRC, i, Kokkos::ALL), Kokkos::subview(NPSACWI, i, Kokkos::ALL),
        Kokkos::subview(NPSACWG, i, Kokkos::ALL), Kokkos::subview(NPRAI, i, Kokkos::ALL),
        Kokkos::subview(NMULTS, i, Kokkos::ALL), Kokkos::subview(NMULTG, i, Kokkos::ALL),
        Kokkos::subview(NMULTR, i, Kokkos::ALL), Kokkos::subview(NMULTRG, i, Kokkos::ALL),
        Kokkos::subview(NNUCCD, i, Kokkos::ALL), Kokkos::subview(NSUBI, i, Kokkos::ALL),
        Kokkos::subview(NGMLTG, i, Kokkos::ALL), Kokkos::subview(NSUBG, i, Kokkos::ALL),
        Kokkos::subview(NACT, i, Kokkos::ALL), Kokkos::subview(SIZEFIX_NR, i, Kokkos::ALL),
        Kokkos::subview(SIZEFIX_NC, i, Kokkos::ALL), Kokkos::subview(SIZEFIX_NI, i, Kokkos::ALL),
        Kokkos::subview(SIZEFIX_NS, i, Kokkos::ALL), Kokkos::subview(SIZEFIX_NG, i, Kokkos::ALL),
        Kokkos::subview(NEGFIX_NR, i, Kokkos::ALL), Kokkos::subview(NEGFIX_NC, i, Kokkos::ALL),
        Kokkos::subview(NEGFIX_NI, i, Kokkos::ALL), Kokkos::subview(NEGFIX_NS, i, Kokkos::ALL),
        Kokkos::subview(NEGFIX_NG, i, Kokkos::ALL), Kokkos::subview(NIM_MORR_CL, i, Kokkos::ALL),
        Kokkos::subview(QC_INST, i, Kokkos::ALL), Kokkos::subview(QR_INST, i, Kokkos::ALL),
        Kokkos::subview(QI_INST, i, Kokkos::ALL), Kokkos::subview(QS_INST, i, Kokkos::ALL),
        Kokkos::subview(QG_INST, i, Kokkos::ALL), Kokkos::subview(NC_INST, i, Kokkos::ALL),
        Kokkos::subview(NR_INST, i, Kokkos::ALL), Kokkos::subview(NI_INST, i, Kokkos::ALL),
        Kokkos::subview(NS_INST, i, Kokkos::ALL), Kokkos::subview(NG_INST, i, Kokkos::ALL),
        morrison);
  }

  if (clubb_at_least_debug_level_api(2)) {
    print_morr_error_output();
    // clubb_at_least_debug_level_api( 2 )
  }

  for (k = 0; k < nzt; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      // BFB-032: Native rounds grav*z before fusing Cp*T into the energy sum.
      hl_after(i, k) =
          ((Kokkos::fma(Cp, static_cast<real_t>(T_in_K(i, k)), grav * gr_zt(i, k)) -
            (Lv * (static_cast<real_t>(rcm_r4(i, k)) + static_cast<real_t>(rrm_r4(i, k))))) -
           (Ls * ((static_cast<real_t>(rim_r4(i, k)) + static_cast<real_t>(rsm_r4(i, k))) +
                  static_cast<real_t>(rgm_r4(i, k)))));
    }
  }

  for (int index_1 = 0; index_1 < hl_on_Cp_residual.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hl_on_Cp_residual.extent_int(0); ++index_0) {
      hl_on_Cp_residual(index_0, index_1) =
          ((((hl_after(index_0, index_1) - hl_before(index_0, index_1)) -
             ((dt * Lv) * (static_cast<real_t>(rcm_sten(index_0, index_1)) +
                           static_cast<real_t>(rrm_sten(index_0, index_1))))) -
            ((dt * Ls) * ((static_cast<real_t>(rim_sten(index_0, index_1)) +
                           static_cast<real_t>(rsm_sten(index_0, index_1))) +
                          static_cast<real_t>(rgm_sten(index_0, index_1))))) /
           Cp);
    }
  }

  for (int index_1 = 0; index_1 < qto_after.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < qto_after.extent_int(0); ++index_0) {
      qto_after(index_0, index_1) = (((((static_cast<real_t>(rvm_r4(index_0, index_1)) +
                                         static_cast<real_t>(rcm_r4(index_0, index_1))) +
                                        static_cast<real_t>(rrm_r4(index_0, index_1))) +
                                       static_cast<real_t>(rim_r4(index_0, index_1))) +
                                      static_cast<real_t>(rsm_r4(index_0, index_1))) +
                                     static_cast<real_t>(rgm_r4(index_0, index_1)));
    }
  }

  for (int index_1 = 0; index_1 < qto_residual.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < qto_residual.extent_int(0); ++index_0) {
      qto_residual(index_0, index_1) =
          ((qto_after(index_0, index_1) - qto_before(index_0, index_1)) -
           (dt * ((((static_cast<real_t>(rcm_sten(index_0, index_1)) +
                     static_cast<real_t>(rrm_sten(index_0, index_1))) +
                    static_cast<real_t>(rim_sten(index_0, index_1))) +
                   static_cast<real_t>(rsm_sten(index_0, index_1))) +
                  static_cast<real_t>(rgm_sten(index_0, index_1)))));
    }
  }

  // Pack hydrometeor arrays.
  for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
      hydromet_r4(index_0, index_1, iirr) = rrm_r4(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
      hydromet_r4(index_0, index_1, iiNr) = Nrm_r4(index_0, index_1);
    }
  }

  for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
      hydromet_mc_r4(index_0, index_1, iirr) = rrm_mc_r4(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
      hydromet_mc_r4(index_0, index_1, iiNr) = Nrm_mc_r4(index_0, index_1);
    }
  }

  if (microphys.l_ice_microphys) {

    for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
        hydromet_r4(index_0, index_1, iiri) = rim_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
        hydromet_r4(index_0, index_1, iiNi) = Nim_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
        hydromet_r4(index_0, index_1, iirs) = rsm_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
        hydromet_r4(index_0, index_1, iiNs) = Nsm_r4(index_0, index_1);
      }
    }

    for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
        hydromet_mc_r4(index_0, index_1, iiri) = rim_mc_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
        hydromet_mc_r4(index_0, index_1, iiNi) = Nim_mc_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
        hydromet_mc_r4(index_0, index_1, iirs) = rsm_mc_r4(index_0, index_1);
      }
    }
    for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
        hydromet_mc_r4(index_0, index_1, iiNs) = Nsm_mc_r4(index_0, index_1);
      }
    }

    if (microphys.l_graupel) {

      for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
          hydromet_r4(index_0, index_1, iirg) = rgm_r4(index_0, index_1);
        }
      }
      for (int index_1 = 0; index_1 < hydromet_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_r4.extent_int(0); ++index_0) {
          hydromet_r4(index_0, index_1, iiNg) = Ngm_r4(index_0, index_1);
        }
      }

      for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
          hydromet_mc_r4(index_0, index_1, iirg) = rgm_mc_r4(index_0, index_1);
        }
      }
      for (int index_1 = 0; index_1 < hydromet_mc_r4.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_mc_r4.extent_int(0); ++index_0) {
          hydromet_mc_r4(index_0, index_1, iiNg) = Ngm_mc_r4(index_0, index_1);
        }
      }
    }
  }

  //hydromet_mc = real( hydromet_mc_r4, kind = core_rknd )
  for (int index_1 = 0; index_1 < rcm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_mc.extent_int(0); ++index_0) {
      rcm_mc(index_0, index_1) = static_cast<real_t>(rcm_mc_r4(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rvm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rvm_mc.extent_int(0); ++index_0) {
      rvm_mc(index_0, index_1) = static_cast<real_t>(rvm_mc_r4(index_0, index_1));
    }
  }

  for (int index_1 = 0; index_1 < rrm_auto.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_auto.extent_int(0); ++index_0) {
      rrm_auto(index_0, index_1) = static_cast<real_t>(PRC(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rrm_accr.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_accr.extent_int(0); ++index_0) {
      rrm_accr(index_0, index_1) = static_cast<real_t>(PRA(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < rrm_evap.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_evap.extent_int(0); ++index_0) {
      rrm_evap(index_0, index_1) = static_cast<real_t>(PRE(index_0, index_1));
    }
  }

  for (int index_1 = 0; index_1 < Nrm_auto.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_auto.extent_int(0); ++index_0) {
      Nrm_auto(index_0, index_1) = static_cast<real_t>(NPRC1(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < Nrm_evap.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_evap.extent_int(0); ++index_0) {
      Nrm_evap(index_0, index_1) = static_cast<real_t>(NSUBR(index_0, index_1));
    }
  }
  // Update hydrometeor tendencies
  // This done because the hydromet_mc arrays that are produced by
  // M2005MICRO_GRAUPEL don't include the clipping term.
  for (int index_2 = 0; index_2 < hydromet_mc.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
        hydromet_mc(index_0, index_1, index_2) =
            ((static_cast<real_t>(hydromet_r4(index_0, index_1, index_2)) -
              hydromet(index_0, index_1, index_2)) /
             dt);
      }
    }
  }

  for (int index_1 = 0; index_1 < Ncm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_mc.extent_int(0); ++index_0) {
      Ncm_mc(index_0, index_1) =
          ((static_cast<real_t>(Ncm_r4(index_0, index_1)) - Ncm(index_0, index_1)) / dt);
    }
  }

  // Update thetal based on absolute temperature
  for (int index_1 = 0; index_1 < thlm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < thlm_mc.extent_int(0); ++index_0) {
      thlm_mc(index_0, index_1) =
          ((T_in_K2thlm_api(static_cast<real_t>(T_in_K(index_0, index_1)), exner(index_0, index_1),
                            static_cast<real_t>(rcm_r4(index_0, index_1))) -
            thlm(index_0, index_1)) /
           dt);
    }
  }

  // Sedimentation is handled within the Morrison microphysics
  for (int index_2 = 0; index_2 < hydromet_vel_zt.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < hydromet_vel_zt.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < hydromet_vel_zt.extent_int(0); ++index_0) {
        hydromet_vel_zt(index_0, index_1, index_2) = real_t(0.0);
      }
    }
  }

  // Output rain sedimentation velocity
  // Multiply by -1 so that negative is associated with falling precip
  for (int index_1 = 0; index_1 < morr_rain_vel_r4.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < morr_rain_vel_r4.extent_int(0); ++index_0) {
      morr_rain_vel_r4(index_0, index_1) = (morr_rain_vel_r4(index_0, index_1) * (-1.0f));
    }
  }
  for (int index_1 = 0; index_1 < hydromet_vel_zt.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < hydromet_vel_zt.extent_int(0); ++index_0) {
      hydromet_vel_zt(index_0, index_1, iirr) =
          static_cast<real_t>(morr_rain_vel_r4(index_0, index_1));
    }
  }

  for (int index_0 = 0; index_0 < rsm_sd_morr_int.extent_int(0); ++index_0) {
    rsm_sd_morr_int(index_0) = zero;
  }
  for (k = 0; k < nzt; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      rsm_sd_morr_int(i) =
          (rsm_sd_morr_int(i) + ((rho(i, k) * static_cast<real_t>(rsm_sten(i, k))) * gr_dzt(i, k)));
    }
  }

  update_microphys_stat_sfc(std::string("rs_sd_morr_int"), rsm_sd_morr_int);

  HostReal1D array_actual_1 = workspace_scope.view<HostReal1D>(
      "morrison_microphys_driver::array_actual_1", rsm_sd_morr_int.extent_int(0));
  HostReal2D array_actual_2 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_2", rsm_sten.extent_int(0), rsm_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rsm_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rsm_sten.extent_int(0); ++index_0) {
      array_actual_2(index_0, index_1) = static_cast<real_t>(rsm_sten(index_0, index_1));
    }
  }
  for (int index_0 = 0; index_0 < rsm_sd_morr_int.extent_int(0); ++index_0) {
    array_actual_1(index_0) =
        (rsm_sd_morr_int(index_0) > microphys_host_detail::maxval_levels(array_actual_2)(index_0));
  }
  if ((clubb_at_least_debug_level_api(1) && microphys_host_detail::any_true(array_actual_1))) {
    std::cerr << (std::string("Warning: rsm_sd_morr was not conservative!") +
                  std::string(" rsm_sd_morr_verical_integr = "))
              << " " << microphys_detail::dump_values(rsm_sd_morr_int) << "\n";
  }

  update_microphys_stat(std::string("rrm_auto"), rrm_auto);
  update_microphys_stat(std::string("rrm_accr"), rrm_accr);
  update_microphys_stat(std::string("rrm_evap"), rrm_evap);

  update_microphys_stat(std::string("Nrm_auto"), Nrm_auto);
  update_microphys_stat(std::string("Nrm_evap"), Nrm_evap);

  // Update Morrison budgets
  update_microphys_stat(std::string("hl_on_Cp_residual"), hl_on_Cp_residual);
  update_microphys_stat(std::string("qto_residual"), qto_residual);
  HostReal2D array_actual_3 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_3", rgm_sten.extent_int(0), rgm_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rgm_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rgm_sten.extent_int(0); ++index_0) {
      array_actual_3(index_0, index_1) = static_cast<real_t>(rgm_sten(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("rgm_sd_morr"), array_actual_3);
  HostReal2D array_actual_4 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_4", rrm_sten.extent_int(0), rrm_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rrm_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_sten.extent_int(0); ++index_0) {
      array_actual_4(index_0, index_1) = static_cast<real_t>(rrm_sten(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("rrm_sd_morr"), array_actual_4);
  HostReal2D array_actual_5 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_5", rsm_sten.extent_int(0), rsm_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rsm_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rsm_sten.extent_int(0); ++index_0) {
      array_actual_5(index_0, index_1) = static_cast<real_t>(rsm_sten(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("rsm_sd_morr"), array_actual_5);
  HostReal2D array_actual_6 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_6", rim_sten.extent_int(0), rim_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rim_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rim_sten.extent_int(0); ++index_0) {
      array_actual_6(index_0, index_1) = static_cast<real_t>(rim_sten(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("rim_sd_mg_morr"), array_actual_6);
  HostReal2D array_actual_7 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_7", rcm_sten.extent_int(0), rcm_sten.extent_int(1));
  for (int index_1 = 0; index_1 < rcm_sten.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_sten.extent_int(0); ++index_0) {
      array_actual_7(index_0, index_1) = static_cast<real_t>(rcm_sten(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("rcm_sd_mg_morr"), array_actual_7);
  HostReal2D array_actual_8 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_8", PRC.extent_int(0), PRC.extent_int(1));
  for (int index_1 = 0; index_1 < PRC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRC.extent_int(0); ++index_0) {
      array_actual_8(index_0, index_1) = static_cast<real_t>(PRC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRC"), array_actual_8);
  HostReal2D array_actual_9 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_9", PRA.extent_int(0), PRA.extent_int(1));
  for (int index_1 = 0; index_1 < PRA.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRA.extent_int(0); ++index_0) {
      array_actual_9(index_0, index_1) = static_cast<real_t>(PRA(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRA"), array_actual_9);
  HostReal2D array_actual_10 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_10", PRE.extent_int(0), PRE.extent_int(1));
  for (int index_1 = 0; index_1 < PRE.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRE.extent_int(0); ++index_0) {
      array_actual_10(index_0, index_1) = static_cast<real_t>(PRE(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRE"), array_actual_10);
  HostReal2D array_actual_11 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_11", PSMLT.extent_int(0), PSMLT.extent_int(1));
  for (int index_1 = 0; index_1 < PSMLT.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSMLT.extent_int(0); ++index_0) {
      array_actual_11(index_0, index_1) = static_cast<real_t>(PSMLT(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PSMLT"), array_actual_11);
  HostReal2D array_actual_12 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_12", EVPMS.extent_int(0), EVPMS.extent_int(1));
  for (int index_1 = 0; index_1 < EVPMS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EVPMS.extent_int(0); ++index_0) {
      array_actual_12(index_0, index_1) = static_cast<real_t>(EVPMS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("EVPMS"), array_actual_12);
  HostReal2D array_actual_13 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_13", PRACS.extent_int(0), PRACS.extent_int(1));
  for (int index_1 = 0; index_1 < PRACS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACS.extent_int(0); ++index_0) {
      array_actual_13(index_0, index_1) = static_cast<real_t>(PRACS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRACS"), array_actual_13);
  HostReal2D array_actual_14 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_14", EVPMG.extent_int(0), EVPMG.extent_int(1));
  for (int index_1 = 0; index_1 < EVPMG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EVPMG.extent_int(0); ++index_0) {
      array_actual_14(index_0, index_1) = static_cast<real_t>(EVPMG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("EVPMG"), array_actual_14);
  HostReal2D array_actual_15 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_15", PRACG.extent_int(0), PRACG.extent_int(1));
  for (int index_1 = 0; index_1 < PRACG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACG.extent_int(0); ++index_0) {
      array_actual_15(index_0, index_1) = static_cast<real_t>(PRACG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRACG"), array_actual_15);
  HostReal2D array_actual_16 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_16", PGMLT.extent_int(0), PGMLT.extent_int(1));
  for (int index_1 = 0; index_1 < PGMLT.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PGMLT.extent_int(0); ++index_0) {
      array_actual_16(index_0, index_1) = static_cast<real_t>(PGMLT(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PGMLT"), array_actual_16);
  HostReal2D array_actual_17 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_17", MNUCCC.extent_int(0), MNUCCC.extent_int(1));
  for (int index_1 = 0; index_1 < MNUCCC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCC.extent_int(0); ++index_0) {
      array_actual_17(index_0, index_1) = static_cast<real_t>(MNUCCC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("MNUCCC"), array_actual_17);
  HostReal2D array_actual_18 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_18", PSACWS.extent_int(0), PSACWS.extent_int(1));
  for (int index_1 = 0; index_1 < PSACWS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWS.extent_int(0); ++index_0) {
      array_actual_18(index_0, index_1) = static_cast<real_t>(PSACWS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PSACWS"), array_actual_18);
  HostReal2D array_actual_19 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_19", PSACWI.extent_int(0), PSACWI.extent_int(1));
  for (int index_1 = 0; index_1 < PSACWI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWI.extent_int(0); ++index_0) {
      array_actual_19(index_0, index_1) = static_cast<real_t>(PSACWI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PSACWI"), array_actual_19);
  HostReal2D array_actual_20 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_20", QMULTS.extent_int(0), QMULTS.extent_int(1));
  for (int index_1 = 0; index_1 < QMULTS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTS.extent_int(0); ++index_0) {
      array_actual_20(index_0, index_1) = static_cast<real_t>(QMULTS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QMULTS"), array_actual_20);
  HostReal2D array_actual_21 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_21", QMULTG.extent_int(0), QMULTG.extent_int(1));
  for (int index_1 = 0; index_1 < QMULTG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTG.extent_int(0); ++index_0) {
      array_actual_21(index_0, index_1) = static_cast<real_t>(QMULTG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QMULTG"), array_actual_21);
  HostReal2D array_actual_22 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_22", PSACWG.extent_int(0), PSACWG.extent_int(1));
  for (int index_1 = 0; index_1 < PSACWG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACWG.extent_int(0); ++index_0) {
      array_actual_22(index_0, index_1) = static_cast<real_t>(PSACWG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PSACWG"), array_actual_22);
  HostReal2D array_actual_23 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_23", PGSACW.extent_int(0), PGSACW.extent_int(1));
  for (int index_1 = 0; index_1 < PGSACW.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PGSACW.extent_int(0); ++index_0) {
      array_actual_23(index_0, index_1) = static_cast<real_t>(PGSACW(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PGSACW"), array_actual_23);
  HostReal2D array_actual_24 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_24", PRD.extent_int(0), PRD.extent_int(1));
  for (int index_1 = 0; index_1 < PRD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRD.extent_int(0); ++index_0) {
      array_actual_24(index_0, index_1) = static_cast<real_t>(PRD(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRD"), array_actual_24);
  HostReal2D array_actual_25 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_25", PRCI.extent_int(0), PRCI.extent_int(1));
  for (int index_1 = 0; index_1 < PRCI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRCI.extent_int(0); ++index_0) {
      array_actual_25(index_0, index_1) = static_cast<real_t>(PRCI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRCI"), array_actual_25);
  HostReal2D array_actual_26 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_26", PRAI.extent_int(0), PRAI.extent_int(1));
  for (int index_1 = 0; index_1 < PRAI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRAI.extent_int(0); ++index_0) {
      array_actual_26(index_0, index_1) = static_cast<real_t>(PRAI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRAI"), array_actual_26);
  HostReal2D array_actual_27 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_27", QMULTR.extent_int(0), QMULTR.extent_int(1));
  for (int index_1 = 0; index_1 < QMULTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTR.extent_int(0); ++index_0) {
      array_actual_27(index_0, index_1) = static_cast<real_t>(QMULTR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QMULTR"), array_actual_27);
  HostReal2D array_actual_28 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_28", QMULTRG.extent_int(0), QMULTRG.extent_int(1));
  for (int index_1 = 0; index_1 < QMULTRG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QMULTRG.extent_int(0); ++index_0) {
      array_actual_28(index_0, index_1) = static_cast<real_t>(QMULTRG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QMULTRG"), array_actual_28);
  HostReal2D array_actual_29 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_29", MNUCCD.extent_int(0), MNUCCD.extent_int(1));
  for (int index_1 = 0; index_1 < MNUCCD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCD.extent_int(0); ++index_0) {
      array_actual_29(index_0, index_1) = static_cast<real_t>(MNUCCD(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("MNUCCD"), array_actual_29);
  HostReal2D array_actual_30 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_30", PRACI.extent_int(0), PRACI.extent_int(1));
  for (int index_1 = 0; index_1 < PRACI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACI.extent_int(0); ++index_0) {
      array_actual_30(index_0, index_1) = static_cast<real_t>(PRACI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRACI"), array_actual_30);
  HostReal2D array_actual_31 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_31", PRACIS.extent_int(0), PRACIS.extent_int(1));
  for (int index_1 = 0; index_1 < PRACIS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRACIS.extent_int(0); ++index_0) {
      array_actual_31(index_0, index_1) = static_cast<real_t>(PRACIS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRACIS"), array_actual_31);
  HostReal2D array_actual_32 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_32", EPRD.extent_int(0), EPRD.extent_int(1));
  for (int index_1 = 0; index_1 < EPRD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRD.extent_int(0); ++index_0) {
      array_actual_32(index_0, index_1) = static_cast<real_t>(EPRD(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("EPRD"), array_actual_32);
  HostReal2D array_actual_33 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_33", MNUCCR.extent_int(0), MNUCCR.extent_int(1));
  for (int index_1 = 0; index_1 < MNUCCR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < MNUCCR.extent_int(0); ++index_0) {
      array_actual_33(index_0, index_1) = static_cast<real_t>(MNUCCR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("MNUCCR"), array_actual_33);
  HostReal2D array_actual_34 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_34", PIACR.extent_int(0), PIACR.extent_int(1));
  for (int index_1 = 0; index_1 < PIACR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PIACR.extent_int(0); ++index_0) {
      array_actual_34(index_0, index_1) = static_cast<real_t>(PIACR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PIACR"), array_actual_34);
  HostReal2D array_actual_35 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_35", PIACRS.extent_int(0), PIACRS.extent_int(1));
  for (int index_1 = 0; index_1 < PIACRS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PIACRS.extent_int(0); ++index_0) {
      array_actual_35(index_0, index_1) = static_cast<real_t>(PIACRS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PIACRS"), array_actual_35);
  HostReal2D array_actual_36 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_36", PGRACS.extent_int(0), PGRACS.extent_int(1));
  for (int index_1 = 0; index_1 < PGRACS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PGRACS.extent_int(0); ++index_0) {
      array_actual_36(index_0, index_1) = static_cast<real_t>(PGRACS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PGRACS"), array_actual_36);
  HostReal2D array_actual_37 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_37", PRDS.extent_int(0), PRDS.extent_int(1));
  for (int index_1 = 0; index_1 < PRDS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRDS.extent_int(0); ++index_0) {
      array_actual_37(index_0, index_1) = static_cast<real_t>(PRDS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRDS"), array_actual_37);
  HostReal2D array_actual_38 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_38", EPRDS.extent_int(0), EPRDS.extent_int(1));
  for (int index_1 = 0; index_1 < EPRDS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRDS.extent_int(0); ++index_0) {
      array_actual_38(index_0, index_1) = static_cast<real_t>(EPRDS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("EPRDS"), array_actual_38);
  HostReal2D array_actual_39 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_39", PSACR.extent_int(0), PSACR.extent_int(1));
  for (int index_1 = 0; index_1 < PSACR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PSACR.extent_int(0); ++index_0) {
      array_actual_39(index_0, index_1) = static_cast<real_t>(PSACR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PSACR"), array_actual_39);
  HostReal2D array_actual_40 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_40", PRDG.extent_int(0), PRDG.extent_int(1));
  for (int index_1 = 0; index_1 < PRDG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PRDG.extent_int(0); ++index_0) {
      array_actual_40(index_0, index_1) = static_cast<real_t>(PRDG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PRDG"), array_actual_40);
  HostReal2D array_actual_41 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_41", EPRDG.extent_int(0), EPRDG.extent_int(1));
  for (int index_1 = 0; index_1 < EPRDG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < EPRDG.extent_int(0); ++index_0) {
      array_actual_41(index_0, index_1) = static_cast<real_t>(EPRDG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("EPRDG"), array_actual_41);

  // Update more Morrison budgets
  HostReal2D array_actual_42 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_42", NGSTEN.extent_int(0), NGSTEN.extent_int(1));
  for (int index_1 = 0; index_1 < NGSTEN.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NGSTEN.extent_int(0); ++index_0) {
      array_actual_42(index_0, index_1) = static_cast<real_t>(NGSTEN(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NGSTEN"), array_actual_42);
  HostReal2D array_actual_43 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_43", NRSTEN.extent_int(0), NRSTEN.extent_int(1));
  for (int index_1 = 0; index_1 < NRSTEN.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NRSTEN.extent_int(0); ++index_0) {
      array_actual_43(index_0, index_1) = static_cast<real_t>(NRSTEN(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NRSTEN"), array_actual_43);
  HostReal2D array_actual_44 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_44", NISTEN.extent_int(0), NISTEN.extent_int(1));
  for (int index_1 = 0; index_1 < NISTEN.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NISTEN.extent_int(0); ++index_0) {
      array_actual_44(index_0, index_1) = static_cast<real_t>(NISTEN(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NISTEN"), array_actual_44);
  HostReal2D array_actual_45 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_45", NSSTEN.extent_int(0), NSSTEN.extent_int(1));
  for (int index_1 = 0; index_1 < NSSTEN.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSSTEN.extent_int(0); ++index_0) {
      array_actual_45(index_0, index_1) = static_cast<real_t>(NSSTEN(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSSTEN"), array_actual_45);
  HostReal2D array_actual_46 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_46", NCSTEN.extent_int(0), NCSTEN.extent_int(1));
  for (int index_1 = 0; index_1 < NCSTEN.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NCSTEN.extent_int(0); ++index_0) {
      array_actual_46(index_0, index_1) = static_cast<real_t>(NCSTEN(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NCSTEN"), array_actual_46);
  HostReal2D array_actual_47 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_47", NPRC1.extent_int(0), NPRC1.extent_int(1));
  for (int index_1 = 0; index_1 < NPRC1.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRC1.extent_int(0); ++index_0) {
      array_actual_47(index_0, index_1) = static_cast<real_t>(NPRC1(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRC1"), array_actual_47);
  HostReal2D array_actual_48 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_48", NRAGG.extent_int(0), NRAGG.extent_int(1));
  for (int index_1 = 0; index_1 < NRAGG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NRAGG.extent_int(0); ++index_0) {
      array_actual_48(index_0, index_1) = static_cast<real_t>(NRAGG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NRAGG"), array_actual_48);
  HostReal2D array_actual_49 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_49", NPRACG.extent_int(0), NPRACG.extent_int(1));
  for (int index_1 = 0; index_1 < NPRACG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRACG.extent_int(0); ++index_0) {
      array_actual_49(index_0, index_1) = static_cast<real_t>(NPRACG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRACG"), array_actual_49);
  HostReal2D array_actual_50 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_50", NSUBR.extent_int(0), NSUBR.extent_int(1));
  for (int index_1 = 0; index_1 < NSUBR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBR.extent_int(0); ++index_0) {
      array_actual_50(index_0, index_1) = static_cast<real_t>(NSUBR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSUBR"), array_actual_50);
  HostReal2D array_actual_51 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_51", NSMLTR.extent_int(0), NSMLTR.extent_int(1));
  for (int index_1 = 0; index_1 < NSMLTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSMLTR.extent_int(0); ++index_0) {
      array_actual_51(index_0, index_1) = static_cast<real_t>(NSMLTR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSMLTR"), array_actual_51);
  HostReal2D array_actual_52 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_52", NGMLTR.extent_int(0), NGMLTR.extent_int(1));
  for (int index_1 = 0; index_1 < NGMLTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NGMLTR.extent_int(0); ++index_0) {
      array_actual_52(index_0, index_1) = static_cast<real_t>(NGMLTR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NGMLTR"), array_actual_52);
  HostReal2D array_actual_53 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_53", NPRACS.extent_int(0), NPRACS.extent_int(1));
  for (int index_1 = 0; index_1 < NPRACS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRACS.extent_int(0); ++index_0) {
      array_actual_53(index_0, index_1) = static_cast<real_t>(NPRACS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRACS"), array_actual_53);
  HostReal2D array_actual_54 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_54", NNUCCR.extent_int(0), NNUCCR.extent_int(1));
  for (int index_1 = 0; index_1 < NNUCCR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NNUCCR.extent_int(0); ++index_0) {
      array_actual_54(index_0, index_1) = static_cast<real_t>(NNUCCR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NNUCCR"), array_actual_54);
  HostReal2D array_actual_55 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_55", NIACR.extent_int(0), NIACR.extent_int(1));
  for (int index_1 = 0; index_1 < NIACR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NIACR.extent_int(0); ++index_0) {
      array_actual_55(index_0, index_1) = static_cast<real_t>(NIACR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NIACR"), array_actual_55);
  HostReal2D array_actual_56 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_56", NIACRS.extent_int(0), NIACRS.extent_int(1));
  for (int index_1 = 0; index_1 < NIACRS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NIACRS.extent_int(0); ++index_0) {
      array_actual_56(index_0, index_1) = static_cast<real_t>(NIACRS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NIACRS"), array_actual_56);
  HostReal2D array_actual_57 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_57", NGRACS.extent_int(0), NGRACS.extent_int(1));
  for (int index_1 = 0; index_1 < NGRACS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NGRACS.extent_int(0); ++index_0) {
      array_actual_57(index_0, index_1) = static_cast<real_t>(NGRACS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NGRACS"), array_actual_57);
  HostReal2D array_actual_58 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_58", NSMLTS.extent_int(0), NSMLTS.extent_int(1));
  for (int index_1 = 0; index_1 < NSMLTS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSMLTS.extent_int(0); ++index_0) {
      array_actual_58(index_0, index_1) = static_cast<real_t>(NSMLTS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSMLTS"), array_actual_58);
  HostReal2D array_actual_59 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_59", NSAGG.extent_int(0), NSAGG.extent_int(1));
  for (int index_1 = 0; index_1 < NSAGG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSAGG.extent_int(0); ++index_0) {
      array_actual_59(index_0, index_1) = static_cast<real_t>(NSAGG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSAGG"), array_actual_59);
  HostReal2D array_actual_60 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_60", NPRCI.extent_int(0), NPRCI.extent_int(1));
  for (int index_1 = 0; index_1 < NPRCI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRCI.extent_int(0); ++index_0) {
      array_actual_60(index_0, index_1) = static_cast<real_t>(NPRCI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRCI"), array_actual_60);
  HostReal2D array_actual_61 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_61", NSCNG.extent_int(0), NSCNG.extent_int(1));
  for (int index_1 = 0; index_1 < NSCNG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSCNG.extent_int(0); ++index_0) {
      array_actual_61(index_0, index_1) = static_cast<real_t>(NSCNG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSCNG"), array_actual_61);
  HostReal2D array_actual_62 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_62", NSUBS.extent_int(0), NSUBS.extent_int(1));
  for (int index_1 = 0; index_1 < NSUBS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBS.extent_int(0); ++index_0) {
      array_actual_62(index_0, index_1) = static_cast<real_t>(NSUBS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSUBS"), array_actual_62);

  HostReal2D array_actual_63 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_63", PCC.extent_int(0), PCC.extent_int(1));
  for (int index_1 = 0; index_1 < PCC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < PCC.extent_int(0); ++index_0) {
      array_actual_63(index_0, index_1) = static_cast<real_t>(PCC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("PCC"), array_actual_63);
  HostReal2D array_actual_64 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_64", NNUCCC.extent_int(0), NNUCCC.extent_int(1));
  for (int index_1 = 0; index_1 < NNUCCC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NNUCCC.extent_int(0); ++index_0) {
      array_actual_64(index_0, index_1) = static_cast<real_t>(NNUCCC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NNUCCC"), array_actual_64);
  HostReal2D array_actual_65 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_65", NPSACWS.extent_int(0), NPSACWS.extent_int(1));
  for (int index_1 = 0; index_1 < NPSACWS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWS.extent_int(0); ++index_0) {
      array_actual_65(index_0, index_1) = static_cast<real_t>(NPSACWS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPSACWS"), array_actual_65);
  HostReal2D array_actual_66 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_66", NPRA.extent_int(0), NPRA.extent_int(1));
  for (int index_1 = 0; index_1 < NPRA.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRA.extent_int(0); ++index_0) {
      array_actual_66(index_0, index_1) = static_cast<real_t>(NPRA(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRA"), array_actual_66);
  HostReal2D array_actual_67 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_67", NPRC.extent_int(0), NPRC.extent_int(1));
  for (int index_1 = 0; index_1 < NPRC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRC.extent_int(0); ++index_0) {
      array_actual_67(index_0, index_1) = static_cast<real_t>(NPRC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRC"), array_actual_67);
  HostReal2D array_actual_68 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_68", NPSACWI.extent_int(0), NPSACWI.extent_int(1));
  for (int index_1 = 0; index_1 < NPSACWI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWI.extent_int(0); ++index_0) {
      array_actual_68(index_0, index_1) = static_cast<real_t>(NPSACWI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPSACWI"), array_actual_68);
  HostReal2D array_actual_69 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_69", NPSACWG.extent_int(0), NPSACWG.extent_int(1));
  for (int index_1 = 0; index_1 < NPSACWG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPSACWG.extent_int(0); ++index_0) {
      array_actual_69(index_0, index_1) = static_cast<real_t>(NPSACWG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPSACWG"), array_actual_69);
  HostReal2D array_actual_70 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_70", NPRAI.extent_int(0), NPRAI.extent_int(1));
  for (int index_1 = 0; index_1 < NPRAI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NPRAI.extent_int(0); ++index_0) {
      array_actual_70(index_0, index_1) = static_cast<real_t>(NPRAI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NPRAI"), array_actual_70);
  HostReal2D array_actual_71 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_71", NMULTS.extent_int(0), NMULTS.extent_int(1));
  for (int index_1 = 0; index_1 < NMULTS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTS.extent_int(0); ++index_0) {
      array_actual_71(index_0, index_1) = static_cast<real_t>(NMULTS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NMULTS"), array_actual_71);
  HostReal2D array_actual_72 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_72", NMULTG.extent_int(0), NMULTG.extent_int(1));
  for (int index_1 = 0; index_1 < NMULTG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTG.extent_int(0); ++index_0) {
      array_actual_72(index_0, index_1) = static_cast<real_t>(NMULTG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NMULTG"), array_actual_72);
  HostReal2D array_actual_73 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_73", NMULTR.extent_int(0), NMULTR.extent_int(1));
  for (int index_1 = 0; index_1 < NMULTR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTR.extent_int(0); ++index_0) {
      array_actual_73(index_0, index_1) = static_cast<real_t>(NMULTR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NMULTR"), array_actual_73);
  HostReal2D array_actual_74 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_74", NMULTRG.extent_int(0), NMULTRG.extent_int(1));
  for (int index_1 = 0; index_1 < NMULTRG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NMULTRG.extent_int(0); ++index_0) {
      array_actual_74(index_0, index_1) = static_cast<real_t>(NMULTRG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NMULTRG"), array_actual_74);
  HostReal2D array_actual_75 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_75", NNUCCD.extent_int(0), NNUCCD.extent_int(1));
  for (int index_1 = 0; index_1 < NNUCCD.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NNUCCD.extent_int(0); ++index_0) {
      array_actual_75(index_0, index_1) = static_cast<real_t>(NNUCCD(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NNUCCD"), array_actual_75);
  HostReal2D array_actual_76 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_76", NSUBI.extent_int(0), NSUBI.extent_int(1));
  for (int index_1 = 0; index_1 < NSUBI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBI.extent_int(0); ++index_0) {
      array_actual_76(index_0, index_1) = static_cast<real_t>(NSUBI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSUBI"), array_actual_76);
  HostReal2D array_actual_77 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_77", NGMLTG.extent_int(0), NGMLTG.extent_int(1));
  for (int index_1 = 0; index_1 < NGMLTG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NGMLTG.extent_int(0); ++index_0) {
      array_actual_77(index_0, index_1) = static_cast<real_t>(NGMLTG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NGMLTG"), array_actual_77);
  HostReal2D array_actual_78 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_78", NSUBG.extent_int(0), NSUBG.extent_int(1));
  for (int index_1 = 0; index_1 < NSUBG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NSUBG.extent_int(0); ++index_0) {
      array_actual_78(index_0, index_1) = static_cast<real_t>(NSUBG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NSUBG"), array_actual_78);
  HostReal2D array_actual_79 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_79", NACT.extent_int(0), NACT.extent_int(1));
  for (int index_1 = 0; index_1 < NACT.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NACT.extent_int(0); ++index_0) {
      array_actual_79(index_0, index_1) = static_cast<real_t>(NACT(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NACT"), array_actual_79);
  HostReal2D array_actual_80 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_80",
                                       SIZEFIX_NR.extent_int(0), SIZEFIX_NR.extent_int(1));
  for (int index_1 = 0; index_1 < SIZEFIX_NR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < SIZEFIX_NR.extent_int(0); ++index_0) {
      array_actual_80(index_0, index_1) = static_cast<real_t>(SIZEFIX_NR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("SIZEFIX_NR"), array_actual_80);
  HostReal2D array_actual_81 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_81",
                                       SIZEFIX_NC.extent_int(0), SIZEFIX_NC.extent_int(1));
  for (int index_1 = 0; index_1 < SIZEFIX_NC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < SIZEFIX_NC.extent_int(0); ++index_0) {
      array_actual_81(index_0, index_1) = static_cast<real_t>(SIZEFIX_NC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("SIZEFIX_NC"), array_actual_81);
  HostReal2D array_actual_82 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_82",
                                       SIZEFIX_NI.extent_int(0), SIZEFIX_NI.extent_int(1));
  for (int index_1 = 0; index_1 < SIZEFIX_NI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < SIZEFIX_NI.extent_int(0); ++index_0) {
      array_actual_82(index_0, index_1) = static_cast<real_t>(SIZEFIX_NI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("SIZEFIX_NI"), array_actual_82);
  HostReal2D array_actual_83 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_83",
                                       SIZEFIX_NS.extent_int(0), SIZEFIX_NS.extent_int(1));
  for (int index_1 = 0; index_1 < SIZEFIX_NS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < SIZEFIX_NS.extent_int(0); ++index_0) {
      array_actual_83(index_0, index_1) = static_cast<real_t>(SIZEFIX_NS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("SIZEFIX_NS"), array_actual_83);
  HostReal2D array_actual_84 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_84",
                                       SIZEFIX_NG.extent_int(0), SIZEFIX_NG.extent_int(1));
  for (int index_1 = 0; index_1 < SIZEFIX_NG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < SIZEFIX_NG.extent_int(0); ++index_0) {
      array_actual_84(index_0, index_1) = static_cast<real_t>(SIZEFIX_NG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("SIZEFIX_NG"), array_actual_84);
  HostReal2D array_actual_85 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_85",
                                       NEGFIX_NR.extent_int(0), NEGFIX_NR.extent_int(1));
  for (int index_1 = 0; index_1 < NEGFIX_NR.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NEGFIX_NR.extent_int(0); ++index_0) {
      array_actual_85(index_0, index_1) = static_cast<real_t>(NEGFIX_NR(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NEGFIX_NR"), array_actual_85);
  HostReal2D array_actual_86 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_86",
                                       NEGFIX_NC.extent_int(0), NEGFIX_NC.extent_int(1));
  for (int index_1 = 0; index_1 < NEGFIX_NC.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NEGFIX_NC.extent_int(0); ++index_0) {
      array_actual_86(index_0, index_1) = static_cast<real_t>(NEGFIX_NC(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NEGFIX_NC"), array_actual_86);
  HostReal2D array_actual_87 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_87",
                                       NEGFIX_NI.extent_int(0), NEGFIX_NI.extent_int(1));
  for (int index_1 = 0; index_1 < NEGFIX_NI.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NEGFIX_NI.extent_int(0); ++index_0) {
      array_actual_87(index_0, index_1) = static_cast<real_t>(NEGFIX_NI(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NEGFIX_NI"), array_actual_87);
  HostReal2D array_actual_88 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_88",
                                       NEGFIX_NS.extent_int(0), NEGFIX_NS.extent_int(1));
  for (int index_1 = 0; index_1 < NEGFIX_NS.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NEGFIX_NS.extent_int(0); ++index_0) {
      array_actual_88(index_0, index_1) = static_cast<real_t>(NEGFIX_NS(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NEGFIX_NS"), array_actual_88);
  HostReal2D array_actual_89 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_89",
                                       NEGFIX_NG.extent_int(0), NEGFIX_NG.extent_int(1));
  for (int index_1 = 0; index_1 < NEGFIX_NG.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NEGFIX_NG.extent_int(0); ++index_0) {
      array_actual_89(index_0, index_1) = static_cast<real_t>(NEGFIX_NG(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NEGFIX_NG"), array_actual_89);
  HostReal2D array_actual_90 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_90",
                                       NIM_MORR_CL.extent_int(0), NIM_MORR_CL.extent_int(1));
  for (int index_1 = 0; index_1 < NIM_MORR_CL.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NIM_MORR_CL.extent_int(0); ++index_0) {
      array_actual_90(index_0, index_1) = static_cast<real_t>(NIM_MORR_CL(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NIM_MORR_CL"), array_actual_90);
  HostReal2D array_actual_91 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_91", QC_INST.extent_int(0), QC_INST.extent_int(1));
  for (int index_1 = 0; index_1 < QC_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QC_INST.extent_int(0); ++index_0) {
      array_actual_91(index_0, index_1) = static_cast<real_t>(QC_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QC_INST"), array_actual_91);
  HostReal2D array_actual_92 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_92", QR_INST.extent_int(0), QR_INST.extent_int(1));
  for (int index_1 = 0; index_1 < QR_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QR_INST.extent_int(0); ++index_0) {
      array_actual_92(index_0, index_1) = static_cast<real_t>(QR_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QR_INST"), array_actual_92);
  HostReal2D array_actual_93 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_93", QI_INST.extent_int(0), QI_INST.extent_int(1));
  for (int index_1 = 0; index_1 < QI_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QI_INST.extent_int(0); ++index_0) {
      array_actual_93(index_0, index_1) = static_cast<real_t>(QI_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QI_INST"), array_actual_93);
  HostReal2D array_actual_94 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_94", QS_INST.extent_int(0), QS_INST.extent_int(1));
  for (int index_1 = 0; index_1 < QS_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QS_INST.extent_int(0); ++index_0) {
      array_actual_94(index_0, index_1) = static_cast<real_t>(QS_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QS_INST"), array_actual_94);
  HostReal2D array_actual_95 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_95", QG_INST.extent_int(0), QG_INST.extent_int(1));
  for (int index_1 = 0; index_1 < QG_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < QG_INST.extent_int(0); ++index_0) {
      array_actual_95(index_0, index_1) = static_cast<real_t>(QG_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("QG_INST"), array_actual_95);
  HostReal2D array_actual_96 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_96", NC_INST.extent_int(0), NC_INST.extent_int(1));
  for (int index_1 = 0; index_1 < NC_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NC_INST.extent_int(0); ++index_0) {
      array_actual_96(index_0, index_1) = static_cast<real_t>(NC_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NC_INST"), array_actual_96);
  HostReal2D array_actual_97 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_97", NR_INST.extent_int(0), NR_INST.extent_int(1));
  for (int index_1 = 0; index_1 < NR_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NR_INST.extent_int(0); ++index_0) {
      array_actual_97(index_0, index_1) = static_cast<real_t>(NR_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NR_INST"), array_actual_97);
  HostReal2D array_actual_98 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_98", NI_INST.extent_int(0), NI_INST.extent_int(1));
  for (int index_1 = 0; index_1 < NI_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NI_INST.extent_int(0); ++index_0) {
      array_actual_98(index_0, index_1) = static_cast<real_t>(NI_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NI_INST"), array_actual_98);
  HostReal2D array_actual_99 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_99", NS_INST.extent_int(0), NS_INST.extent_int(1));
  for (int index_1 = 0; index_1 < NS_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NS_INST.extent_int(0); ++index_0) {
      array_actual_99(index_0, index_1) = static_cast<real_t>(NS_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NS_INST"), array_actual_99);
  HostReal2D array_actual_100 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_100", NG_INST.extent_int(0), NG_INST.extent_int(1));
  for (int index_1 = 0; index_1 < NG_INST.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < NG_INST.extent_int(0); ++index_0) {
      array_actual_100(index_0, index_1) = static_cast<real_t>(NG_INST(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("NG_INST"), array_actual_100);

  HostReal2D array_actual_101 =
      workspace_scope.view<HostReal2D>("morrison_microphys_driver::array_actual_101",
                                       T_in_K_mc.extent_int(0), T_in_K_mc.extent_int(1));
  for (int index_1 = 0; index_1 < T_in_K_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < T_in_K_mc.extent_int(0); ++index_0) {
      array_actual_101(index_0, index_1) = static_cast<real_t>(T_in_K_mc(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("T_in_K_mc"), array_actual_101);

  // --- Number concentrations ---
  // No budgets for sedimentation are output
  // Effective radii of hydrometeor species
  HostReal2D array_actual_102 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_102", effc.extent_int(0), effc.extent_int(1));
  for (int index_1 = 0; index_1 < effc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < effc.extent_int(0); ++index_0) {
      array_actual_102(index_0, index_1) = static_cast<real_t>(effc(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("eff_rad_cloud"), array_actual_102);
  HostReal2D array_actual_103 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_103", effi.extent_int(0), effi.extent_int(1));
  for (int index_1 = 0; index_1 < effi.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < effi.extent_int(0); ++index_0) {
      array_actual_103(index_0, index_1) = static_cast<real_t>(effi(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("eff_rad_ice"), array_actual_103);
  HostReal2D array_actual_104 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_104", effs.extent_int(0), effs.extent_int(1));
  for (int index_1 = 0; index_1 < effs.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < effs.extent_int(0); ++index_0) {
      array_actual_104(index_0, index_1) = static_cast<real_t>(effs(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("eff_rad_snow"), array_actual_104);
  HostReal2D array_actual_105 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_105", effr.extent_int(0), effr.extent_int(1));
  for (int index_1 = 0; index_1 < effr.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < effr.extent_int(0); ++index_0) {
      array_actual_105(index_0, index_1) = static_cast<real_t>(effr(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("eff_rad_rain"), array_actual_105);
  HostReal2D array_actual_106 = workspace_scope.view<HostReal2D>(
      "morrison_microphys_driver::array_actual_106", effg.extent_int(0), effg.extent_int(1));
  for (int index_1 = 0; index_1 < effg.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < effg.extent_int(0); ++index_0) {
      array_actual_106(index_0, index_1) = static_cast<real_t>(effg(index_0, index_1));
    }
  }
  update_microphys_stat(std::string("eff_rad_graupel"), array_actual_106);

  // Snow and Rain rates at the bottom of the domain, in mm/day
  HostReal1D array_actual_107 = workspace_scope.view<HostReal1D>(
      "morrison_microphys_driver::array_actual_107", Morr_precip_rate.extent_int(0));
  for (int index_0 = 0; index_0 < Morr_precip_rate.extent_int(0); ++index_0) {
    array_actual_107(index_0) =
        ((static_cast<real_t>(Morr_precip_rate(index_0)) * sec_per_day) / dt);
  }
  update_microphys_stat_sfc(std::string("precip_rate_sfc"), array_actual_107);

  HostReal1D array_actual_108 = workspace_scope.view<HostReal1D>(
      "morrison_microphys_driver::array_actual_108", Morr_snow_rate.extent_int(0));
  for (int index_0 = 0; index_0 < Morr_snow_rate.extent_int(0); ++index_0) {
    array_actual_108(index_0) = ((static_cast<real_t>(Morr_snow_rate(index_0)) * sec_per_day) / dt);
  }
  update_microphys_stat_sfc(std::string("morr_snow_rate"), array_actual_108);

  for (int index_1 = 0; index_1 < rrm_auto_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_auto_diag.extent_int(0); ++index_0) {
      rrm_auto_diag(index_0, index_1) = rrm_auto(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < rrm_accr_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_accr_diag.extent_int(0); ++index_0) {
      rrm_accr_diag(index_0, index_1) = rrm_accr(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < rrm_evap_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_evap_diag.extent_int(0); ++index_0) {
      rrm_evap_diag(index_0, index_1) = rrm_evap(index_0, index_1);
    }
  }
  for (int index_1 = 0; index_1 < Nrm_auto_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_auto_diag.extent_int(0); ++index_0) {
      Nrm_auto_diag(index_0, index_1) = static_cast<real_t>(Nrm_auto(index_0, index_1));
    }
  }
  for (int index_1 = 0; index_1 < Nrm_evap_diag.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_evap_diag.extent_int(0); ++index_0) {
      Nrm_evap_diag(index_0, index_1) = static_cast<real_t>(Nrm_evap(index_0, index_1));
    }
  }

  return;
}
} // namespace clubb_kokkos::microphys_host

namespace clubb_kokkos {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24>
inline void morrison_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, bool l_latin_hypercube, Array0 thlm, Array1 wm_zt,
    Array2 p_in_Pa, Array3 exner, Array4 rho, Array5 cloud_frac, Array6 w_std_dev, Array7 dzq,
    Array8 rcm, Array9 Ncm, Array10 chi, Array11 rvm, Array12 hydromet, int saturation_formula,
    Array13 sample_weight, stats_type& stats, Array14 hydromet_mc, Array15 hydromet_vel_zt,
    Array16 Ncm_mc, Array17 rcm_mc, Array18 rvm_mc, Array19 thlm_mc, Array20 rrm_auto_diag,
    Array21 rrm_accr_diag, Array22 rrm_evap_diag, Array23 Nrm_auto_diag, Array24 Nrm_evap_diag,
    const parameters_microphys& microphys, morrison_microphysics_state& morrison) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  auto workspace_scope      = workspace.scope();
  auto thlm_host            = to_host_space(workspace, thlm);
  auto wm_zt_host           = to_host_space(workspace, wm_zt);
  auto p_in_Pa_host         = to_host_space(workspace, p_in_Pa);
  auto exner_host           = to_host_space(workspace, exner);
  auto rho_host             = to_host_space(workspace, rho);
  auto cloud_frac_host      = to_host_space(workspace, cloud_frac);
  auto w_std_dev_host       = to_host_space(workspace, w_std_dev);
  auto dzq_host             = to_host_space(workspace, dzq);
  auto rcm_host             = to_host_space(workspace, rcm);
  auto Ncm_host             = to_host_space(workspace, Ncm);
  auto chi_host             = to_host_space(workspace, chi);
  auto rvm_host             = to_host_space(workspace, rvm);
  auto hydromet_host        = to_host_space(workspace, hydromet);
  auto sample_weight_host   = to_host_space(workspace, sample_weight);
  auto hydromet_mc_host     = to_host_space(workspace, hydromet_mc);
  auto hydromet_vel_zt_host = to_host_space(workspace, hydromet_vel_zt);
  auto Ncm_mc_host          = to_host_space(workspace, Ncm_mc);
  auto rcm_mc_host          = to_host_space(workspace, rcm_mc);
  auto rvm_mc_host          = to_host_space(workspace, rvm_mc);
  auto thlm_mc_host         = to_host_space(workspace, thlm_mc);
  auto rrm_auto_diag_host   = to_host_space(workspace, rrm_auto_diag);
  auto rrm_accr_diag_host   = to_host_space(workspace, rrm_accr_diag);
  auto rrm_evap_diag_host   = to_host_space(workspace, rrm_evap_diag);
  auto Nrm_auto_diag_host   = to_host_space(workspace, Nrm_auto_diag);
  auto Nrm_evap_diag_host   = to_host_space(workspace, Nrm_evap_diag);
  microphys_host::morrison_microphys_driver(
      workspace, gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube, thlm_host,
      wm_zt_host, p_in_Pa_host, exner_host, rho_host, cloud_frac_host, w_std_dev_host, dzq_host,
      rcm_host, Ncm_host, chi_host, rvm_host, hydromet_host, saturation_formula, sample_weight_host,
      stats, hydromet_mc_host, hydromet_vel_zt_host, Ncm_mc_host, rcm_mc_host, rvm_mc_host,
      thlm_mc_host, rrm_auto_diag_host, rrm_accr_diag_host, rrm_evap_diag_host, Nrm_auto_diag_host,
      Nrm_evap_diag_host, microphys, morrison);
  from_execution_space(workspace, hydromet_mc, hydromet_mc_host);
  from_execution_space(workspace, hydromet_vel_zt, hydromet_vel_zt_host);
  from_execution_space(workspace, Ncm_mc, Ncm_mc_host);
  from_execution_space(workspace, rcm_mc, rcm_mc_host);
  from_execution_space(workspace, rvm_mc, rvm_mc_host);
  from_execution_space(workspace, thlm_mc, thlm_mc_host);
  from_execution_space(workspace, rrm_auto_diag, rrm_auto_diag_host);
  from_execution_space(workspace, rrm_accr_diag, rrm_accr_diag_host);
  from_execution_space(workspace, rrm_evap_diag, rrm_evap_diag_host);
  from_execution_space(workspace, Nrm_auto_diag, Nrm_auto_diag_host);
  from_execution_space(workspace, Nrm_evap_diag, Nrm_evap_diag_host);
}
} // namespace clubb_kokkos
