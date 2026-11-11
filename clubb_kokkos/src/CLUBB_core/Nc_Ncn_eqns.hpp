// Convert cloud droplet and simplified cloud nuclei concentration.
// Four public elemental formulas relate in-cloud/overall Nc and overall Ncn;
// two private PDF integrals retain the Native ERFC branches and thresholds.
// Scalar functions are device callable. Homogeneous rank1/rank2 View overloads
// preserve caller memory space and Native argument lists, without bridge calls.
#pragma once
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <type_traits>
namespace clubb_kokkos {
// Description:
// Equations are provided to perform calculations back-and-forth between Nc and
// Ncn, where Nc is cloud droplet concentration and Ncn is simplified cloud
// nuclei concentration.  The equation that relates the two is:
//
// Nc = Ncn * H(chi);
//
// where chi is extended liquid water mixing ratio, which is equal to cloud
// water mixing ratio, rc, when both are positive.  However, chi is negative in
// subsaturated air.
//
// Equation are provided relating mean cloud droplet concentration (overall),
// Ncm, and/or mean cloud droplet concentration (in-cloud), Nc_in_cloud, to
// mean simplified cloud nuclei concentration, Ncnm.
// Notes:
//
// Meaning of Nc flag combinations:
//
// l_const_Nc_in_cloud:
// When this flag is enabled, cloud droplet concentration (in-cloud) is
// constant (spatially) at a grid level (it is constant over the subgrid
// domain, but could vary over time depending on the value of l_predict_Nc).
// The value of in-cloud Nc does not vary at a grid level.  This also means
// that Ncn is constant across the entire grid level.  When this flag is turned
// off, both in-cloud Nc and Ncn vary at a grid level.
//
// l_predict_Nc:
// When this flag is enabled, Nc_in_cloud (or alternatively Ncm) is predicted.
// It is advanced every time step by a predictive equation, and can change
// at every time step at a grid level.  When this flag is turned off,
// Nc_in_cloud does not change at a grid level over the course of a model run.
//
// 1) l_predict_Nc turned on and l_const_Nc_in_cloud turned on:
//    The value of Nc_in_cloud (mean in-cloud Nc) is predicted and can change
//    at every timestep at a grid level.  However, the value of in-cloud Nc is
//    constant (spatially) at a grid level (no subgrid variability).
//
// 2) l_predict_Nc turned on and l_const_Nc_in_cloud turned off:
//    The value of Nc_in_cloud (mean in-cloud Nc) is predicted and can change
//    at every timestep at a grid level.  The value of in-cloud Nc also varies
//    (spatially) at a grid level (subgrid variability around mean
//    in-cloud Nc).
//
// 3) l_predict_Nc turned off and l_const_Nc_in_cloud turned on:
//    The value of Nc_in_cloud (mean in-cloud Nc) is constant over time at a
//    grid level.  It retains its initial value.  Additionally, the value of
//    in-cloud Nc is constant (spatially) at a grid level (no subgrid
//    variability).  This configuration is used most often in idealized cases.
//
// 4) l_predict_Nc turned off and l_const_Nc_in_cloud turned off:
//    The value of Nc_in_cloud (mean in-cloud Nc) is constant over time at a
//    grid level.  It retains its initial value.  However, the value of
//    in-cloud Nc varies (spatially) at a grid level (subgrid variability
//    around mean in-cloud Nc).
//
//
//
// Nc_in_cloud/Nc - Ncn flow chart of CLUBB code:
//
// (Please update when warranted).
//
//
// Ncm/Nc-in-cloud                                         Ncnm/Ncn PDF params.
// --->
// |  |                    Start of CLUBB main time step loop
// |  |
// |  |                    advance_clubb_core
// |  |
// |  |
// |  |\
// |  | \
// |  |  (intent in)-------setup_pdf_parameters-------->calc. Ncnm (local)
// |  |                                                       |
// |  |                                                      \ /
// |  |                                             mu_Ncn_i, sigma_Ncn_i,
// |  |                                                  corr_xNcn_i
// |  |                                                       |
// |  |                                                      \ /
// |  |                                               PDF param. arrays:
// |  |                                             mu_x_i_n, sigma_x_i_n,
// |  |                                                 corr_array_i_n
// |  |                                                  (intent out)
// |  |                                                       |
// |  |                                                       |
// |  |                                                       |
// |  |                                                       |
// |  |                                                       |
// |  |                                                       |
// |  |--(intent in)---calc_microphys_scheme_tendcies----(intent in)
// |  |                            |
// |  |                            |
// |  |                call a microphysics scheme
// |  |                            |
// |  |  Local micro. scheme-----------Latin Hypercube-----------Upscaled KK
// |  |         |                            |                        |
// |  |  Ncm/Nc-in-cloud:          Populate sample points      Use PDF params.
// |  |  used to find micro.       using PDF params (Ncn).         of Ncn
// |  |    tendencies.             At every sample point:      (mu_Ncn_i, etc.)
// |  |         |                    Nc = Ncn * H(chi).         to find micro.
// |  |         |                  Use sample-point Nc to         tendencies.
// |  |         |                  find micro. tendencies             |
// |  |         |                when calling micro. scheme.          |
// |  |         |                            |                        |
// |  |    hydromet_mc/-----------------hydromet_mc/-------------hydromet_mc
// |  |    Ncm_mc (intent out)     |    Ncm_mc (intent out)      (intent out)
// |  |                            |
// |  |                            |
// |  |                            |
// |  |                            |
// |  |                            |
// |  |                            |
// |  |                            |
// |  |                        (intent in)
// |  |                            |
// |  |--(intent inout)----advance_microphys
// |  |
// |  |
// |  | advance microphysics variables (hydromet, Nc_in_cloud/Ncm) one timestep
// |  |
// |  |                   l_predict_Nc = true:
// |  |            Nc_in_cloud/Ncm necessary for starting
// |  |            value of Nc_in_cloud/Ncm when advancing
// |  |            one timestep using predictive equation.
// |  |
// |  |
// |  |                    End of CLUBB main time step loop
// <---
// References:
//-------------------------------------------------------------------------
// mu_chi_1: Mean of chi (old s) (1st PDF component)      [kg/kg]
// mu_chi_2: Mean of chi (old s) (2nd PDF component)      [kg/kg]
// mu_Ncn_1: Mean of Ncn (1st PDF component)             [num/kg]
// mu_Ncn_2: Mean of Ncn (2nd PDF component)             [num/kg]
// sigma_chi_1: Standard deviation of chi (1st PDF comp.)    [kg/kg]
// sigma_chi_2: Standard deviation of chi (2nd PDF comp.)    [kg/kg]
// sigma_Ncn_1: Standard deviation of Ncn (1st PDF comp.)   [num/kg]
// sigma_Ncn_2: Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
// sigma_Ncn_1_n: Standard deviation of ln Ncn (1st PDF component) [-]
// sigma_Ncn_2_n: Standard deviation of ln Ncn (2nd PDF component) [-]
// corr_chi_Ncn_1_n: Correlation of chi and ln Ncn (1st PDF comp.)    [-]
// corr_chi_Ncn_2_n: Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
// mixt_frac: Mixture fraction                                 [-]
// cloud_frac_1: Cloud fraction (1st PDF component)               [-]
// cloud_frac_2: Cloud fraction (2nd PDF component)               [-]
KOKKOS_INLINE_FUNCTION real_t Ncnm_to_Nc_in_cloud(real_t mu_chi_1, real_t mu_chi_2,
                                                  real_t mu_Ncn_1, // In
                                                  real_t mu_Ncn_2, real_t sigma_chi_1,
                                                  real_t sigma_chi_2,                         // In
                                                  real_t sigma_Ncn_1, real_t sigma_Ncn_2,     // In
                                                  real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n, // In
                                                  real_t corr_chi_Ncn_1_n,
                                                  real_t corr_chi_Ncn_2_n,               // In
                                                  real_t mixt_frac, real_t cloud_frac_1, // In
                                                  real_t cloud_frac_2                    // In
);
// mu_chi_1: Mean of chi (old s) (1st PDF component)           [kg/kg]
// mu_chi_2: Mean of chi (old s) (2nd PDF component)           [kg/kg]
// sigma_chi_1: Standard deviation of chi (1st PDF component)     [kg/kg]
// sigma_chi_2: Standard deviation of chi (2nd PDF component)     [kg/kg]
// mixt_frac: Mixture fraction                                      [-]
// Nc_in_cloud: Mean cloud droplet conc. (in-cloud)     [num/kg]
// cloud_frac_1: Cloud fraction (1st PDF component)           [-]
// cloud_frac_2: Cloud fraction (2nd PDF component)           [-]
// const_Ncnp2_on_Ncnm2: Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
// const_corr_chi_Ncn: Prescribed correlation of chi and Ncn        [-]
KOKKOS_INLINE_FUNCTION real_t Nc_in_cloud_to_Ncnm(real_t mu_chi_1, real_t mu_chi_2,
                                                  real_t sigma_chi_1, // In
                                                  real_t sigma_chi_2, real_t mixt_frac,
                                                  real_t Nc_in_cloud,                       // In
                                                  real_t cloud_frac_1, real_t cloud_frac_2, // In
                                                  real_t const_Ncnp2_on_Ncnm2,              // In
                                                  real_t const_corr_chi_Ncn                 // In
);
// mu_chi_1: Mean of chi (old s) (1st PDF component)      [kg/kg]
// mu_chi_2: Mean of chi (old s) (2nd PDF component)      [kg/kg]
// mu_Ncn_1: Mean of Ncn (1st PDF component)             [num/kg]
// mu_Ncn_2: Mean of Ncn (2nd PDF component)             [num/kg]
// sigma_chi_1: Standard deviation of chi (1st PDF comp.)    [kg/kg]
// sigma_chi_2: Standard deviation of chi (2nd PDF comp.)    [kg/kg]
// sigma_Ncn_1: Standard deviation of Ncn (1st PDF comp.)   [num/kg]
// sigma_Ncn_2: Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
// sigma_Ncn_1_n: Standard deviation of ln Ncn (1st PDF component) [-]
// sigma_Ncn_2_n: Standard deviation of ln Ncn (2nd PDF component) [-]
// corr_chi_Ncn_1_n: Correlation of chi and ln Ncn (1st PDF comp.)    [-]
// corr_chi_Ncn_2_n: Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
// mixt_frac: Mixture fraction                                 [-]
KOKKOS_INLINE_FUNCTION real_t Ncnm_to_Ncm(real_t mu_chi_1, real_t mu_chi_2, real_t mu_Ncn_1,
                                          real_t mu_Ncn_2, // In
                                          real_t sigma_chi_1, real_t sigma_chi_2,
                                          real_t sigma_Ncn_1, // In
                                          real_t sigma_Ncn_2, real_t sigma_Ncn_1_n,
                                          real_t sigma_Ncn_2_n,                             // In
                                          real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, // In
                                          real_t mixt_frac                                  // In
);
// mu_chi_1: Mean of chi (old s) (1st PDF component)           [kg/kg]
// mu_chi_2: Mean of chi (old s) (2nd PDF component)           [kg/kg]
// sigma_chi_1: Standard deviation of chi (1st PDF component)     [kg/kg]
// sigma_chi_2: Standard deviation of chi (2nd PDF component)     [kg/kg]
// mixt_frac: Mixture fraction                                      [-]
// Ncm: Mean cloud droplet conc. (overall)      [num/kg]
// const_Ncnp2_on_Ncnm2: Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
// const_corr_chi_Ncn: Prescribed correlation of chi and Ncn        [-]
// Ncnm_val_denom_0: Ncnm value -- denominator in eqn. is 0  [num/kg]
KOKKOS_INLINE_FUNCTION real_t Ncm_to_Ncnm(real_t mu_chi_1, real_t mu_chi_2,
                                          real_t sigma_chi_1,                               // In
                                          real_t sigma_chi_2, real_t mixt_frac, real_t Ncm, // In
                                          real_t const_Ncnp2_on_Ncnm2,
                                          real_t const_corr_chi_Ncn, // In
                                          real_t Ncnm_val_denom_0    // In
);
namespace Nc_Ncn_eqns_detail {
KOKKOS_INLINE_FUNCTION real_t bivar_NL_chi_Ncn_mean(real_t mu_chi_i, real_t mu_Ncn_i,
                                                    real_t sigma_chi_i,                       // In
                                                    real_t sigma_Ncn_i, real_t sigma_Ncn_i_n, // In
                                                    real_t corr_chi_Ncn_i_n                   // In
);
KOKKOS_INLINE_FUNCTION real_t bivar_Ncnm_eqn_comp(real_t mu_chi_i, real_t sigma_chi_i, // In
                                                  real_t const_Ncnp2_on_Ncnm2,         // In
                                                  real_t const_corr_chi_Ncn            // In
);
} // namespace Nc_Ncn_eqns_detail
// =============================================================================
KOKKOS_INLINE_FUNCTION real_t Ncnm_to_Nc_in_cloud(real_t mu_chi_1, real_t mu_chi_2,
                                                  real_t mu_Ncn_1, // In
                                                  real_t mu_Ncn_2, real_t sigma_chi_1,
                                                  real_t sigma_chi_2,                         // In
                                                  real_t sigma_Ncn_1, real_t sigma_Ncn_2,     // In
                                                  real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n, // In
                                                  real_t corr_chi_Ncn_1_n,
                                                  real_t corr_chi_Ncn_2_n,               // In
                                                  real_t mixt_frac, real_t cloud_frac_1, // In
                                                  real_t cloud_frac_2                    // In
) {
  // Description:
  // The in-cloud mean of cloud droplet concentration is calculated from the
  // PDF parameters involving simplified cloud nuclei concentration, Ncn, and
  // cloud fraction.  At any point, cloud droplet concentration, Nc, is given
  // by:
  //
  // Nc = Ncn * H(chi);
  //
  // where extended liquid water mixing ratio, chi, is equal to cloud water
  // ratio, rc, when positive.  When the atmosphere is saturated at this point,
  // cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
  // and Nc = 0.
  //
  // The overall mean of cloud droplet concentration, <Nc>, is calculated from
  // the PDF parameters involving Ncn.  The in-cloud mean of cloud droplet
  // concentration is calculated from <Nc> and cloud fraction.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_1: Mean of chi (old s) (1st PDF component)      [kg/kg]
  // mu_chi_2: Mean of chi (old s) (2nd PDF component)      [kg/kg]
  // mu_Ncn_1: Mean of Ncn (1st PDF component)             [num/kg]
  // mu_Ncn_2: Mean of Ncn (2nd PDF component)             [num/kg]
  // sigma_chi_1: Standard deviation of chi (1st PDF comp.)    [kg/kg]
  // sigma_chi_2: Standard deviation of chi (2nd PDF comp.)    [kg/kg]
  // sigma_Ncn_1: Standard deviation of Ncn (1st PDF comp.)   [num/kg]
  // sigma_Ncn_2: Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
  // sigma_Ncn_1_n: Standard deviation of ln Ncn (1st PDF component) [-]
  // sigma_Ncn_2_n: Standard deviation of ln Ncn (2nd PDF component) [-]
  // corr_chi_Ncn_1_n: Correlation of chi and ln Ncn (1st PDF comp.)    [-]
  // corr_chi_Ncn_2_n: Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
  // mixt_frac: Mixture fraction                                 [-]
  // cloud_frac_1: Cloud fraction (1st PDF component)               [-]
  // cloud_frac_2: Cloud fraction (2nd PDF component)               [-]

  real_t Nc_in_cloud{};
  real_t Ncm{};
  real_t cloud_frac{};

  // Calculate overall cloud fraction as calculated by the PDF.
  // The variable cloud_frac is not used here because it is altered by factors
  // such as the trapezoidal rule calculation.
  // Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
  // as long neither of these variables are altered by any factor.  They can
  // only be calculated from PDF.
  cloud_frac = ((mixt_frac * cloud_frac_1) + ((one - mixt_frac) * cloud_frac_2));
  if ((cloud_frac > cloud_frac_min)) {
    // There is cloud found at this grid level.  Calculate Nc_in_cloud.
    Ncm = Ncnm_to_Ncm(mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2,       // In
                      sigma_chi_1, sigma_chi_2, sigma_Ncn_1,        // In
                      sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n,    // In
                      corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, mixt_frac // In
    );
    Nc_in_cloud = (Ncm / cloud_frac);
  } else {
    // cloud_frac <= cloud_frac_min
    // This level is entirely clear.  Set Nc_in_cloud to <Ncn>.
    // Since <Ncn> = mu_Ncn_1 = mu_Ncn_2, use mu_Ncn_1 here.
    Nc_in_cloud = mu_Ncn_1;
  }
  return Nc_in_cloud;
}
// =============================================================================
KOKKOS_INLINE_FUNCTION real_t Nc_in_cloud_to_Ncnm(real_t mu_chi_1, real_t mu_chi_2,
                                                  real_t sigma_chi_1, // In
                                                  real_t sigma_chi_2, real_t mixt_frac,
                                                  real_t Nc_in_cloud,                       // In
                                                  real_t cloud_frac_1, real_t cloud_frac_2, // In
                                                  real_t const_Ncnp2_on_Ncnm2,              // In
                                                  real_t const_corr_chi_Ncn                 // In
) {
  // Description:
  // The overall mean of simplified cloud nuclei concentration, <Ncn>, is
  // calculated from the in-cloud mean of cloud droplet concentration, <Nc>,
  // cloud fraction, and some of the PDF parameters.
  //
  // At any point, cloud droplet concentration, Nc, is given by:
  //
  // Nc = Ncn * H(chi);
  //
  // where extended liquid water mixing ratio, chi, is equal to cloud water
  // ratio, rc, when positive.  When the atmosphere is saturated at this point,
  // cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
  // and Nc = 0.
  //
  // The overall mean of cloud droplet concentration, <Nc>, is calculated from
  // Nc_in_cloud and cloud fraction.  The value of <Ncn> is calculated from
  // <Nc> and PDF parameters.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_1: Mean of chi (old s) (1st PDF component)           [kg/kg]
  // mu_chi_2: Mean of chi (old s) (2nd PDF component)           [kg/kg]
  // sigma_chi_1: Standard deviation of chi (1st PDF component)     [kg/kg]
  // sigma_chi_2: Standard deviation of chi (2nd PDF component)     [kg/kg]
  // mixt_frac: Mixture fraction                                      [-]
  // Nc_in_cloud: Mean cloud droplet conc. (in-cloud)     [num/kg]
  // cloud_frac_1: Cloud fraction (1st PDF component)           [-]
  // cloud_frac_2: Cloud fraction (2nd PDF component)           [-]
  // const_Ncnp2_on_Ncnm2: Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
  // const_corr_chi_Ncn: Prescribed correlation of chi and Ncn        [-]

  real_t Ncnm{};
  real_t Ncm{};
  real_t cloud_frac{};

  // Calculate overall cloud fraction as calculated by the PDF.
  // The variable cloud_frac is not used here because it is altered by factors
  // such as the trapezoidal rule calculation.
  // Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
  // as long neither of these variables are altered by any factor.  They can
  // only be calculated from the PDF.
  cloud_frac = ((mixt_frac * cloud_frac_1) + ((one - mixt_frac) * cloud_frac_2));
  if (((cloud_frac > cloud_frac_min) &&
       (Kokkos::abs((const_corr_chi_Ncn * const_Ncnp2_on_Ncnm2)) > eps))) {
    // There is cloud found at this grid level.  Additionally, Ncn varies.
    // Calculate Nc_in_cloud.
    Ncm = (Nc_in_cloud * cloud_frac);
    Ncnm = Ncm_to_Ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, // In
                       mixt_frac, Ncm, const_Ncnp2_on_Ncnm2,         // In
                       const_corr_chi_Ncn, Nc_in_cloud               // In
    );
  } else {
    // cloud_frac <= cloud_frac_min .or. const_Ncnp2_on_Ncnm2 = 0
    // When Ncn is constant a a grid level, it is equal to Nc_in_cloud.
    // Additionally, when a level is entirely clear, <Ncn>, which is based on
    // Nc_in_cloud, here, must be set to something.  Set <Ncn> to Nc_in_cloud.
    Ncnm = Nc_in_cloud;
  }
  return Ncnm;
}
// =============================================================================
KOKKOS_INLINE_FUNCTION real_t Ncnm_to_Ncm(real_t mu_chi_1, real_t mu_chi_2, real_t mu_Ncn_1,
                                          real_t mu_Ncn_2, // In
                                          real_t sigma_chi_1, real_t sigma_chi_2,
                                          real_t sigma_Ncn_1, // In
                                          real_t sigma_Ncn_2, real_t sigma_Ncn_1_n,
                                          real_t sigma_Ncn_2_n,                             // In
                                          real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, // In
                                          real_t mixt_frac                                  // In
) {
  // Description:
  // The overall mean of cloud droplet concentration, <Nc>, is calculated from
  // the PDF parameters involving the simplified cloud nuclei concentration,
  // Ncn.  At any point, cloud droplet concentration, Nc, is given by:
  //
  // Nc = Ncn * H(chi);
  //
  // where extended liquid water mixing ratio, chi, is equal to cloud water
  // ratio, rc, when positive.  When the atmosphere is saturated at this point,
  // cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
  // and Nc = 0.
  //
  // The overall mean of cloud droplet concentration, <Nc>, is found by
  // integrating over the PDF of chi and Ncn, such that:
  //
  // <Nc> = INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P(chi,Ncn) dNcn dchi;
  //
  // which can also be written as:
  //
  // <Nc> = SUM(i=1,n) mixt_frac_i
  //        * INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
  //
  // where n is the number of multivariate joint PDF components, mixt_frac_i is
  // the weight of the ith PDF component, and P_i is the functional form of the
  // multivariate joint PDF in the ith PDF component.
  //
  // This equation is rewritten as:
  //
  // <Nc> = SUM(i=1,n) mixt_frac_i
  //        * INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi.
  //
  // When both chi and Ncn vary in the ith PDF component, the integral is
  // evaluated and the result is:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
  //
  // which can be reduced to:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * mu_Ncn_i
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
  //
  // When chi is constant, but Ncn varies, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  //
  // When chi varies, but Ncn is constant, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
  //
  // When both chi and Ncn are constant in the ith PDF component, the integral
  // is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_1: Mean of chi (old s) (1st PDF component)      [kg/kg]
  // mu_chi_2: Mean of chi (old s) (2nd PDF component)      [kg/kg]
  // mu_Ncn_1: Mean of Ncn (1st PDF component)             [num/kg]
  // mu_Ncn_2: Mean of Ncn (2nd PDF component)             [num/kg]
  // sigma_chi_1: Standard deviation of chi (1st PDF comp.)    [kg/kg]
  // sigma_chi_2: Standard deviation of chi (2nd PDF comp.)    [kg/kg]
  // sigma_Ncn_1: Standard deviation of Ncn (1st PDF comp.)   [num/kg]
  // sigma_Ncn_2: Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
  // sigma_Ncn_1_n: Standard deviation of ln Ncn (1st PDF component) [-]
  // sigma_Ncn_2_n: Standard deviation of ln Ncn (2nd PDF component) [-]
  // corr_chi_Ncn_1_n: Correlation of chi and ln Ncn (1st PDF comp.)    [-]
  // corr_chi_Ncn_2_n: Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
  // mixt_frac: Mixture fraction                                 [-]

  real_t Ncm{};

  // Calculate mean cloud droplet concentration (overall), <Nc>.
  Ncm = ((mixt_frac * Nc_Ncn_eqns_detail::bivar_NL_chi_Ncn_mean(
                          mu_chi_1, mu_Ncn_1, sigma_chi_1,             // In
                          sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n // In
                          )) +
         ((one - mixt_frac) * Nc_Ncn_eqns_detail::bivar_NL_chi_Ncn_mean(
                                  mu_chi_2, mu_Ncn_2, sigma_chi_2,             // In
                                  sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n // In
                                  )));
  return Ncm;
}
// =============================================================================
KOKKOS_INLINE_FUNCTION real_t Ncm_to_Ncnm(real_t mu_chi_1, real_t mu_chi_2,
                                          real_t sigma_chi_1,                               // In
                                          real_t sigma_chi_2, real_t mixt_frac, real_t Ncm, // In
                                          real_t const_Ncnp2_on_Ncnm2,
                                          real_t const_corr_chi_Ncn, // In
                                          real_t Ncnm_val_denom_0    // In
) {
  // Description:
  // The overall mean of simplified cloud nuclei concentration, <Ncn>, is
  // calculated from the overall mean of cloud droplet concentration, <Nc>, and
  // some of the PDF parameters.
  //
  // At any point, cloud droplet concentration, Nc, is given by:
  //
  // Nc = Ncn * H(chi);
  //
  // where extended liquid water mixing ratio, chi, is equal to cloud water
  // ratio, rc, when positive.  When the atmosphere is saturated at this point,
  // cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
  // and Nc = 0.
  //
  // The overall mean of cloud droplet concentration, <Nc>, is found by
  // integrating over the PDF of chi and Ncn, such that:
  //
  // <Nc> = INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P(chi,Ncn) dNcn dchi;
  //
  // which can also be written as:
  //
  // <Nc> = SUM(i=1,n) mixt_frac_i
  //        * INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
  //
  // where n is the number of multivariate joint PDF components, mixt_frac_i is
  // the weight of the ith PDF component, and P_i is the functional form of the
  // multivariate joint PDF in the ith PDF component.
  //
  // This equation is rewritten as:
  //
  // <Nc> = SUM(i=1,n) mixt_frac_i
  //        * INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi.
  //
  // When both chi and Ncn vary in the ith PDF component, the integral is
  // evaluated and the result is:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
  //
  // which can be reduced to:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * mu_Ncn_i
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
  //
  // When chi is constant, but Ncn varies, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  //
  // When chi varies, but Ncn is constant, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
  //
  // When both chi and Ncn are constant in the ith PDF component, the integral
  // is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  //
  //
  // Solving for <Ncn>
  // =================
  //
  // The individual marginal for simplified cloud nuclei concentration, Ncn, is
  // a single lognormal distribution over the entire horizontal domain.  In
  // order to accomplish this in a two-component PDF structure, the PDF
  // parameters involving Ncn are set equal between the two components.  This
  // results in:
  //
  // mu_Ncn_1 = mu_Ncn_2 = mu_Ncn_i = <Ncn>;
  // mu_Ncn_1_n = mu_Ncn_2_n = mu_Ncn_i_n;
  // sigma_Ncn_1 = sigma_Ncn_2 = sigma_Ncn_i = sqrt( <Ncn'^2> );
  // sigma_Ncn_1_n = sigma_Ncn_2_n = sigma_Ncn_i_n;
  // rho_chi_Ncn_1 = rho_chi_Ncn_2 = rho_chi_Ncn_i = rho_chi_Ncn; and
  // rho_chi_Ncn_1_n = rho_chi_Ncn_2_n = rho_chi_Ncn_i_n.
  //
  // Additionally, the equation for sigma_Ncn_i_n is:
  //
  // sigma_Ncn_i_n = sqrt( ln( 1 + ( sigma_Ncn_i^2 / mu_Ncn_i^2 ) ) );
  //
  // and the equation for rho_chi_Ncn_i_n is:
  //
  // rho_chi_Ncn_i_n
  // = rho_chi_Ncn_i * sqrt( exp{ sigma_Ncn_i_n^2 } - 1 ) / sigma_Ncn_i_n.
  //
  // The product of rho_chi_Ncn_i_n and sigma_Ncn_i_n is:
  //
  // rho_chi_Ncn_i_n * sigma_Ncn_i_n
  // = rho_chi_Ncn_i * sqrt( exp{ sigma_Ncn_i_n^2 } - 1 ).
  //
  // After substituting for sigma_Ncn_i_n^2, the equation for the product of
  // rho_chi_Ncn_i_n and sigma_Ncn_i_n is:
  //
  // rho_chi_Ncn_i_n * sigma_Ncn_i_n
  // = rho_chi_Ncn_i * sqrt( sigma_Ncn_i^2 / mu_Ncn_i^2 );
  //
  // which can be rewritten as:
  //
  // rho_chi_Ncn_i_n * sigma_Ncn_i_n
  // = rho_chi_Ncn * sqrt( <Ncn'^2> / <Ncn>^2 ).
  //
  // Substituting all of this into the equation for <Nc>, the equation for <Nc>
  // becomes:
  //
  // <Nc> = <Ncn>
  //        * SUM(i=1,n) mixt_frac_i
  //          ---
  //          | (1/2) * erfc( - ( 1 / sqrt(2) )
  //          |                 * ( ( mu_chi_i / sigma_chi_i )
  //          |                     + rho_chi_Ncn * sqrt(<Ncn'^2>/<Ncn>^2) ) );
  //          | where sigma_chi_i > 0 and <Ncn'^2> > 0;
  //          |
  //        * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
  //          | where sigma_chi_i > 0 and <Ncn'^2> = 0;
  //          |
  //          | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
  //          |
  //          | 0; where sigma_chi_i = 0 and mu_chi_i <= 0.
  //          ---
  //
  // In order to isolate <Ncn>, the value of <Ncn'^2>/<Ncn>^2 is set to a
  // constant value, const_Ncn.  The value of this constant does not depend on
  // <Ncn>.  Likewise, the value of rho_chi_Ncn does not depend on <Ncn>.
  // Solving for <Ncn>, the equation becomes:
  //
  // <Ncn>
  // = <Nc> / ( SUM(i=1,n) mixt_frac_i
  //              ---
  //              | (1/2) * erfc( - ( 1 / sqrt(2) )
  //              |                 * ( ( mu_chi_i / sigma_chi_i )
  //              |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
  //              | where sigma_chi_i > 0 and const_Ncn > 0;
  //              |
  //            * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
  //              | where sigma_chi_i > 0 and const_Ncn = 0;
  //              |
  //              | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
  //              |
  //              | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
  //              ---
  //
  // When the denominator term is 0, there is only clear air.  Both the
  // numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
  // to an appropriate value.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_1: Mean of chi (old s) (1st PDF component)           [kg/kg]
  // mu_chi_2: Mean of chi (old s) (2nd PDF component)           [kg/kg]
  // sigma_chi_1: Standard deviation of chi (1st PDF component)     [kg/kg]
  // sigma_chi_2: Standard deviation of chi (2nd PDF component)     [kg/kg]
  // mixt_frac: Mixture fraction                                      [-]
  // Ncm: Mean cloud droplet conc. (overall)      [num/kg]
  // const_Ncnp2_on_Ncnm2: Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
  // const_corr_chi_Ncn: Prescribed correlation of chi and Ncn        [-]
  // Ncnm_val_denom_0: Ncnm value -- denominator in eqn. is 0  [num/kg]

  real_t Ncnm{};
  real_t denominator_term{};

  denominator_term = ((mixt_frac * Nc_Ncn_eqns_detail::bivar_Ncnm_eqn_comp(
                                       mu_chi_1, sigma_chi_1,                   // In
                                       const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn // In
                                       )) +
                      ((one - mixt_frac) * Nc_Ncn_eqns_detail::bivar_Ncnm_eqn_comp(
                                               mu_chi_2, sigma_chi_2,                   // In
                                               const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn // In
                                               )));
  if ((denominator_term > zero)) {
    Ncnm = (Ncm / denominator_term);
  } else {
    // denominator_term = 0
    // When the denominator is 0, it is usually because there is only clear
    // air.  In that scenario, Ncm should also be 0.  Set Ncnm to a value that
    // is usual or typical
    Ncnm = Ncnm_val_denom_0;
  }
  // denominator_term > 0
  return Ncnm;
}
// =============================================================================
namespace Nc_Ncn_eqns_detail {
KOKKOS_INLINE_FUNCTION real_t bivar_NL_chi_Ncn_mean(real_t mu_chi_i, real_t mu_Ncn_i,
                                                    real_t sigma_chi_i,                       // In
                                                    real_t sigma_Ncn_i, real_t sigma_Ncn_i_n, // In
                                                    real_t corr_chi_Ncn_i_n                   // In
) {
  // Description:
  // The double integral over Ncn * H(chi) multiplied by the
  // bivariate normal-lognormal joint PDF of chi and Ncn is evaluated.  The
  // integral is given by:
  //
  // INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
  //
  // which reduces to:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi;
  //
  // where the individual marginal distribution of chi is normal in the ith PDF
  // component and the individual marginal distribution of Ncn is lognormal in
  // the ith PDF component.
  //
  // When both chi and Ncn vary in the ith PDF component, the integral is
  // evaluated and the result is:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
  //
  // which can be reduced to:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = (1/2) * mu_Ncn_i
  //   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
  //                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
  //
  // When chi is constant, but Ncn varies, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  //
  // When chi varies, but Ncn is constant, in the ith PDF component, the
  // integral is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
  // = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
  //
  // When both chi and Ncn are constant in the ith PDF component, the integral
  // is evaluated and results in:
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
  //
  // when mu_chi_i > 0; and
  //
  // INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
  //
  // when mu_chi_i <= 0.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_i: Mean of chi (old s) (ith PDF component)      [kg/kg]
  // mu_Ncn_i: Mean of Ncn (ith PDF component)             [num/kg]
  // sigma_chi_i: Standard deviation of chi (ith PDF comp.)    [kg/kg]
  // sigma_Ncn_i: Standard deviation of Ncn (ith PDF comp.)   [num/kg]
  // sigma_Ncn_i_n: Standard deviation of ln Ncn (ith PDF component) [-]
  // corr_chi_Ncn_i_n: Correlation of chi and ln Ncn (ith PDF comp.)    [-]

  real_t bivar_NL_chi_Ncn_mean{};

  if (((sigma_chi_i <= chi_tol) && (sigma_Ncn_i <= Ncn_tol))) {
    // The ith PDF component variances of both chi and Ncn are 0.
    if ((mu_chi_i > zero)) {
      bivar_NL_chi_Ncn_mean = mu_Ncn_i;
    } else {
      // mu_chi_i <= 0
      bivar_NL_chi_Ncn_mean = zero;
    }
  } else if ((sigma_chi_i <= chi_tol)) {
    // The ith PDF component variance of chi is 0.
    if ((mu_chi_i > zero)) {
      bivar_NL_chi_Ncn_mean = mu_Ncn_i;
    } else {
      // mu_chi_i <= 0
      bivar_NL_chi_Ncn_mean = zero;
    }
  } else if ((sigma_Ncn_i <= Ncn_tol)) {
    // The ith PDF component variance of Ncn is 0.
    bivar_NL_chi_Ncn_mean =
        ((mu_Ncn_i * one_half) * Kokkos::erfc(-((mu_chi_i / (sqrt_2 * sigma_chi_i)))));
  } else {
    // Both chi and Ncn vary in the ith PDF component.
    bivar_NL_chi_Ncn_mean =
        ((one_half * mu_Ncn_i) *
         Kokkos::erfc((-((one / sqrt_2)) *
                       ((mu_chi_i / sigma_chi_i) + (corr_chi_Ncn_i_n * sigma_Ncn_i_n)))));
  }
  return bivar_NL_chi_Ncn_mean;
}
} // namespace Nc_Ncn_eqns_detail
// =============================================================================
namespace Nc_Ncn_eqns_detail {
KOKKOS_INLINE_FUNCTION real_t bivar_Ncnm_eqn_comp(real_t mu_chi_i, real_t sigma_chi_i, // In
                                                  real_t const_Ncnp2_on_Ncnm2,         // In
                                                  real_t const_corr_chi_Ncn            // In
) {
  // Description:
  // When <Ncn> is found based on the value of <Nc>, the following equation is
  // used:
  //
  // <Ncn>
  // = <Nc> / ( SUM(i=1,n) mixt_frac_i
  //              ---
  //              | (1/2) * erfc( - ( 1 / sqrt(2) )
  //              |                 * ( ( mu_chi_i / sigma_chi_i )
  //              |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
  //              | where sigma_chi_i > 0 and const_Ncn > 0;
  //              |
  //            * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
  //              | where sigma_chi_i > 0 and const_Ncn = 0;
  //              |
  //              | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
  //              |
  //              | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
  //              ---
  //
  // In the above equation, const_Ncn = <Ncn'^2> / <Ncn>^2.  It is a constant,
  // prescribed parameter.  Likewise, rho_chi_Ncn is a parameter that is not
  // based on the value of <Ncn>.
  //
  // When the denominator term is 0, there is only clear air.  Both the
  // numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
  // to an appropriate value.
  //
  // The contribution of the ith PDF component to the denominator term in the
  // equation is calculated here.
  // References:
  //-----------------------------------------------------------------------
  // mu_chi_i: Mean of chi (old s) (ith PDF component)           [kg/kg]
  // sigma_chi_i: Standard deviation of chi (ith PDF component)     [kg/kg]
  // const_Ncnp2_on_Ncnm2: Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
  // const_corr_chi_Ncn: Prescribed correlation of chi and Ncn        [-]

  real_t bivar_Ncnm_eqn_comp{};

  if ((sigma_chi_i <= chi_tol)) {
    // The ith PDF component variances of chi is 0.  The value of the ith PDF
    // component variance of Ncn does not matter in this scenario.
    if ((mu_chi_i > zero)) {
      bivar_Ncnm_eqn_comp = one;
    } else {
      // mu_chi_i <= 0
      bivar_Ncnm_eqn_comp = zero;
    }
  } else {
    // Both chi and Ncn vary in the ith PDF component.
    bivar_Ncnm_eqn_comp =
        (one_half * Kokkos::erfc((-((one / sqrt_2)) *
                                  ((mu_chi_i / sigma_chi_i) +
                                   (const_corr_chi_Ncn * Kokkos::sqrt(const_Ncnp2_on_Ncnm2))))));
  }
  return bivar_Ncnm_eqn_comp;
}
} // namespace Nc_Ncn_eqns_detail
// =============================================================================
// Native elemental dispatch for the actual rank1/profile and rank2/column-field
// boundaries. All inputs are conformable Views; scalar-only calls use the
// Native scalar overload above. Owned results follow the first input's memory
// space; direct outputs use matching scalar kind/rank and accessible memory.
// CUDA extended lambdas require every enclosing template parameter to be named.
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Output,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
void Ncnm_to_Nc_in_cloud(Field0 mu_chi_1, Field1 mu_chi_2, Field2 mu_Ncn_1,       // In
                         Field3 mu_Ncn_2, Field4 sigma_chi_1, Field5 sigma_chi_2, // In
                         Field6 sigma_Ncn_1, Field7 sigma_Ncn_2,                  // In
                         Field8 sigma_Ncn_1_n, Field9 sigma_Ncn_2_n,              // In
                         Field10 corr_chi_Ncn_1_n, Field11 corr_chi_Ncn_2_n,      // In
                         Field12 mixt_frac, Field13 cloud_frac_1,                 // In
                         Field14 cloud_frac_2,                                    // In
                         Output result                                            // Out
) {
  // Fill a conformable profile or column/level result in the first input's
  // execution space. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == Field0::rank);
  static_assert(
      std::is_same_v<typename Output::non_const_value_type, typename Field0::non_const_value_type>);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 1 || Field0::rank == 2);
  static_assert(Field1::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field1::memory_space>::accessible);
  static_assert(Field2::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field2::memory_space>::accessible);
  static_assert(Field3::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field3::memory_space>::accessible);
  static_assert(Field4::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field4::memory_space>::accessible);
  static_assert(Field5::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field5::memory_space>::accessible);
  static_assert(Field6::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field6::memory_space>::accessible);
  static_assert(Field7::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field7::memory_space>::accessible);
  static_assert(Field8::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field8::memory_space>::accessible);
  static_assert(Field9::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field9::memory_space>::accessible);
  static_assert(Field10::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field10::memory_space>::accessible);
  static_assert(Field11::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field11::memory_space>::accessible);
  static_assert(Field12::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field12::memory_space>::accessible);
  static_assert(Field13::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field13::memory_space>::accessible);
  static_assert(Field14::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field14::memory_space>::accessible);
  using Policy = Kokkos::RangePolicy<typename Field0::execution_space>;
  Kokkos::parallel_for(
      "Ncnm_to_Nc_in_cloud", Policy(0, mu_chi_1.size()), KOKKOS_LAMBDA(const int cell) {
        // CUDA requires each View capture before the rank-specific constexpr branch.
        (void)result;
        (void)mu_chi_1;
        (void)mu_chi_2;
        (void)mu_Ncn_1;
        (void)mu_Ncn_2;
        (void)sigma_chi_1;
        (void)sigma_chi_2;
        (void)sigma_Ncn_1;
        (void)sigma_Ncn_2;
        (void)sigma_Ncn_1_n;
        (void)sigma_Ncn_2_n;
        (void)corr_chi_Ncn_1_n;
        (void)corr_chi_Ncn_2_n;
        (void)mixt_frac;
        (void)cloud_frac_1;
        (void)cloud_frac_2;

        if constexpr (Field0::rank == 1)
          result(cell) = Ncnm_to_Nc_in_cloud(
              mu_chi_1(cell), mu_chi_2(cell), mu_Ncn_1(cell), mu_Ncn_2(cell), sigma_chi_1(cell),
              sigma_chi_2(cell), sigma_Ncn_1(cell), sigma_Ncn_2(cell), sigma_Ncn_1_n(cell),
              sigma_Ncn_2_n(cell), corr_chi_Ncn_1_n(cell), corr_chi_Ncn_2_n(cell), mixt_frac(cell),
              cloud_frac_1(cell), cloud_frac_2(cell));
        else {
          const int column      = cell % mu_chi_1.extent_int(0);
          const int level       = cell / mu_chi_1.extent_int(0);
          result(column, level) = Ncnm_to_Nc_in_cloud(
              mu_chi_1(column, level), mu_chi_2(column, level), mu_Ncn_1(column, level),
              mu_Ncn_2(column, level), sigma_chi_1(column, level), sigma_chi_2(column, level),
              sigma_Ncn_1(column, level), sigma_Ncn_2(column, level), sigma_Ncn_1_n(column, level),
              sigma_Ncn_2_n(column, level), corr_chi_Ncn_1_n(column, level),
              corr_chi_Ncn_2_n(column, level), mixt_frac(column, level),
              cloud_frac_1(column, level), cloud_frac_2(column, level));
        }
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
             typename Field0::memory_space>
Ncnm_to_Nc_in_cloud(Field0 mu_chi_1, Field1 mu_chi_2, Field2 mu_Ncn_1,       // In
                    Field3 mu_Ncn_2, Field4 sigma_chi_1, Field5 sigma_chi_2, // In
                    Field6 sigma_Ncn_1, Field7 sigma_Ncn_2,                  // In
                    Field8 sigma_Ncn_1_n, Field9 sigma_Ncn_2_n,              // In
                    Field10 corr_chi_Ncn_1_n, Field11 corr_chi_Ncn_2_n,      // In
                    Field12 mixt_frac, Field13 cloud_frac_1,                 // In
                    Field14 cloud_frac_2                                     // In
) {
  // Return independent owned storage through the same direct-output kernel.
  using Result = Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
                              typename Field0::memory_space>;
  Result result;
  if constexpr (Field0::rank == 1)
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncnm_to_Nc_in_cloud::result"),
                    mu_chi_1.extent(0));
  else
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncnm_to_Nc_in_cloud::result"),
                    mu_chi_1.extent(0), mu_chi_1.extent(1));
  Ncnm_to_Nc_in_cloud(mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1,
                      sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n,
                      mixt_frac, cloud_frac_1, cloud_frac_2, result);
  return result;
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Output,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
void Nc_in_cloud_to_Ncnm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 sigma_chi_1,     // In
                         Field3 sigma_chi_2, Field4 mixt_frac, Field5 Nc_in_cloud, // In
                         Field6 cloud_frac_1, Field7 cloud_frac_2,                 // In
                         Field8 const_Ncnp2_on_Ncnm2,                              // In
                         Field9 const_corr_chi_Ncn,                                // In
                         Output result                                             // Out
) {
  // Fill a conformable profile or column/level result in the first input's
  // execution space. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == Field0::rank);
  static_assert(
      std::is_same_v<typename Output::non_const_value_type, typename Field0::non_const_value_type>);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 1 || Field0::rank == 2);
  static_assert(Field1::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field1::memory_space>::accessible);
  static_assert(Field2::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field2::memory_space>::accessible);
  static_assert(Field3::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field3::memory_space>::accessible);
  static_assert(Field4::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field4::memory_space>::accessible);
  static_assert(Field5::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field5::memory_space>::accessible);
  static_assert(Field6::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field6::memory_space>::accessible);
  static_assert(Field7::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field7::memory_space>::accessible);
  static_assert(Field8::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field8::memory_space>::accessible);
  static_assert(Field9::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field9::memory_space>::accessible);
  using Policy = Kokkos::RangePolicy<typename Field0::execution_space>;
  Kokkos::parallel_for(
      "Nc_in_cloud_to_Ncnm", Policy(0, mu_chi_1.size()), KOKKOS_LAMBDA(const int cell) {
        // CUDA requires each View capture before the rank-specific constexpr branch.
        (void)result;
        (void)mu_chi_1;
        (void)mu_chi_2;
        (void)sigma_chi_1;
        (void)sigma_chi_2;
        (void)mixt_frac;
        (void)Nc_in_cloud;
        (void)cloud_frac_1;
        (void)cloud_frac_2;
        (void)const_Ncnp2_on_Ncnm2;
        (void)const_corr_chi_Ncn;

        if constexpr (Field0::rank == 1)
          result(cell) = Nc_in_cloud_to_Ncnm(mu_chi_1(cell), mu_chi_2(cell), sigma_chi_1(cell),
                                             sigma_chi_2(cell), mixt_frac(cell), Nc_in_cloud(cell),
                                             cloud_frac_1(cell), cloud_frac_2(cell),
                                             const_Ncnp2_on_Ncnm2(cell), const_corr_chi_Ncn(cell));
        else {
          const int column      = cell % mu_chi_1.extent_int(0);
          const int level       = cell / mu_chi_1.extent_int(0);
          result(column, level) = Nc_in_cloud_to_Ncnm(
              mu_chi_1(column, level), mu_chi_2(column, level), sigma_chi_1(column, level),
              sigma_chi_2(column, level), mixt_frac(column, level), Nc_in_cloud(column, level),
              cloud_frac_1(column, level), cloud_frac_2(column, level),
              const_Ncnp2_on_Ncnm2(column, level), const_corr_chi_Ncn(column, level));
        }
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
             typename Field0::memory_space>
Nc_in_cloud_to_Ncnm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 sigma_chi_1,     // In
                    Field3 sigma_chi_2, Field4 mixt_frac, Field5 Nc_in_cloud, // In
                    Field6 cloud_frac_1, Field7 cloud_frac_2,                 // In
                    Field8 const_Ncnp2_on_Ncnm2,                              // In
                    Field9 const_corr_chi_Ncn                                 // In
) {
  // Return independent owned storage through the same direct-output kernel.
  using Result = Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
                              typename Field0::memory_space>;
  Result result;
  if constexpr (Field0::rank == 1)
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Nc_in_cloud_to_Ncnm::result"),
                    mu_chi_1.extent(0));
  else
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Nc_in_cloud_to_Ncnm::result"),
                    mu_chi_1.extent(0), mu_chi_1.extent(1));
  Nc_in_cloud_to_Ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac, Nc_in_cloud,
                      cloud_frac_1, cloud_frac_2, const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn, result);
  return result;
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Output,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
void Ncnm_to_Ncm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 mu_Ncn_1, Field3 mu_Ncn_2, // In
                 Field4 sigma_chi_1, Field5 sigma_chi_2, Field6 sigma_Ncn_1,         // In
                 Field7 sigma_Ncn_2, Field8 sigma_Ncn_1_n, Field9 sigma_Ncn_2_n,     // In
                 Field10 corr_chi_Ncn_1_n, Field11 corr_chi_Ncn_2_n,                 // In
                 Field12 mixt_frac,                                                  // In
                 Output result                                                       // Out
) {
  // Fill a conformable profile or column/level result in the first input's
  // execution space. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == Field0::rank);
  static_assert(
      std::is_same_v<typename Output::non_const_value_type, typename Field0::non_const_value_type>);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 1 || Field0::rank == 2);
  static_assert(Field1::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field1::memory_space>::accessible);
  static_assert(Field2::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field2::memory_space>::accessible);
  static_assert(Field3::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field3::memory_space>::accessible);
  static_assert(Field4::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field4::memory_space>::accessible);
  static_assert(Field5::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field5::memory_space>::accessible);
  static_assert(Field6::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field6::memory_space>::accessible);
  static_assert(Field7::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field7::memory_space>::accessible);
  static_assert(Field8::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field8::memory_space>::accessible);
  static_assert(Field9::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field9::memory_space>::accessible);
  static_assert(Field10::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field10::memory_space>::accessible);
  static_assert(Field11::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field11::memory_space>::accessible);
  static_assert(Field12::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field12::memory_space>::accessible);
  using Policy = Kokkos::RangePolicy<typename Field0::execution_space>;
  Kokkos::parallel_for(
      "Ncnm_to_Ncm", Policy(0, mu_chi_1.size()), KOKKOS_LAMBDA(const int cell) {
        // CUDA requires each View capture before the rank-specific constexpr branch.
        (void)result;
        (void)mu_chi_1;
        (void)mu_chi_2;
        (void)mu_Ncn_1;
        (void)mu_Ncn_2;
        (void)sigma_chi_1;
        (void)sigma_chi_2;
        (void)sigma_Ncn_1;
        (void)sigma_Ncn_2;
        (void)sigma_Ncn_1_n;
        (void)sigma_Ncn_2_n;
        (void)corr_chi_Ncn_1_n;
        (void)corr_chi_Ncn_2_n;
        (void)mixt_frac;

        if constexpr (Field0::rank == 1)
          result(cell) = Ncnm_to_Ncm(
              mu_chi_1(cell), mu_chi_2(cell), mu_Ncn_1(cell), mu_Ncn_2(cell), sigma_chi_1(cell),
              sigma_chi_2(cell), sigma_Ncn_1(cell), sigma_Ncn_2(cell), sigma_Ncn_1_n(cell),
              sigma_Ncn_2_n(cell), corr_chi_Ncn_1_n(cell), corr_chi_Ncn_2_n(cell), mixt_frac(cell));
        else {
          const int column      = cell % mu_chi_1.extent_int(0);
          const int level       = cell / mu_chi_1.extent_int(0);
          result(column, level) = Ncnm_to_Ncm(
              mu_chi_1(column, level), mu_chi_2(column, level), mu_Ncn_1(column, level),
              mu_Ncn_2(column, level), sigma_chi_1(column, level), sigma_chi_2(column, level),
              sigma_Ncn_1(column, level), sigma_Ncn_2(column, level), sigma_Ncn_1_n(column, level),
              sigma_Ncn_2_n(column, level), corr_chi_Ncn_1_n(column, level),
              corr_chi_Ncn_2_n(column, level), mixt_frac(column, level));
        }
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
             typename Field0::memory_space>
Ncnm_to_Ncm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 mu_Ncn_1, Field3 mu_Ncn_2, // In
            Field4 sigma_chi_1, Field5 sigma_chi_2, Field6 sigma_Ncn_1,         // In
            Field7 sigma_Ncn_2, Field8 sigma_Ncn_1_n, Field9 sigma_Ncn_2_n,     // In
            Field10 corr_chi_Ncn_1_n, Field11 corr_chi_Ncn_2_n,                 // In
            Field12 mixt_frac                                                   // In
) {
  // Return independent owned storage through the same direct-output kernel.
  using Result = Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
                              typename Field0::memory_space>;
  Result result;
  if constexpr (Field0::rank == 1)
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncnm_to_Ncm::result"),
                    mu_chi_1.extent(0));
  else
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncnm_to_Ncm::result"),
                    mu_chi_1.extent(0), mu_chi_1.extent(1));
  Ncnm_to_Ncm(mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1,
              sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n,
              mixt_frac, result);
  return result;
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Output,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
void Ncm_to_Ncnm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 sigma_chi_1,   // In
                 Field3 sigma_chi_2, Field4 mixt_frac, Field5 Ncm,       // In
                 Field6 const_Ncnp2_on_Ncnm2, Field7 const_corr_chi_Ncn, // In
                 Field8 Ncnm_val_denom_0,                                // In
                 Output result                                           // Out
) {
  // Fill a conformable profile or column/level result in the first input's
  // execution space. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == Field0::rank);
  static_assert(
      std::is_same_v<typename Output::non_const_value_type, typename Field0::non_const_value_type>);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 1 || Field0::rank == 2);
  static_assert(Field1::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field1::memory_space>::accessible);
  static_assert(Field2::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field2::memory_space>::accessible);
  static_assert(Field3::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field3::memory_space>::accessible);
  static_assert(Field4::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field4::memory_space>::accessible);
  static_assert(Field5::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field5::memory_space>::accessible);
  static_assert(Field6::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field6::memory_space>::accessible);
  static_assert(Field7::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field7::memory_space>::accessible);
  static_assert(Field8::rank == Field0::rank);
  static_assert(Kokkos::SpaceAccessibility<typename Field0::execution_space,
                                           typename Field8::memory_space>::accessible);
  using Policy = Kokkos::RangePolicy<typename Field0::execution_space>;
  Kokkos::parallel_for(
      "Ncm_to_Ncnm", Policy(0, mu_chi_1.size()), KOKKOS_LAMBDA(const int cell) {
        // CUDA requires each View capture before the rank-specific constexpr branch.
        (void)result;
        (void)mu_chi_1;
        (void)mu_chi_2;
        (void)sigma_chi_1;
        (void)sigma_chi_2;
        (void)mixt_frac;
        (void)Ncm;
        (void)const_Ncnp2_on_Ncnm2;
        (void)const_corr_chi_Ncn;
        (void)Ncnm_val_denom_0;

        if constexpr (Field0::rank == 1)
          result(cell) =
              Ncm_to_Ncnm(mu_chi_1(cell), mu_chi_2(cell), sigma_chi_1(cell), sigma_chi_2(cell),
                          mixt_frac(cell), Ncm(cell), const_Ncnp2_on_Ncnm2(cell),
                          const_corr_chi_Ncn(cell), Ncnm_val_denom_0(cell));
        else {
          const int column      = cell % mu_chi_1.extent_int(0);
          const int level       = cell / mu_chi_1.extent_int(0);
          result(column, level) = Ncm_to_Ncnm(
              mu_chi_1(column, level), mu_chi_2(column, level), sigma_chi_1(column, level),
              sigma_chi_2(column, level), mixt_frac(column, level), Ncm(column, level),
              const_Ncnp2_on_Ncnm2(column, level), const_corr_chi_Ncn(column, level),
              Ncnm_val_denom_0(column, level));
        }
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8,
          std::enable_if_t<Kokkos::is_view_v<Field0>, int> view_input = 0>
Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
             typename Field0::memory_space>
Ncm_to_Ncnm(Field0 mu_chi_1, Field1 mu_chi_2, Field2 sigma_chi_1,   // In
            Field3 sigma_chi_2, Field4 mixt_frac, Field5 Ncm,       // In
            Field6 const_Ncnp2_on_Ncnm2, Field7 const_corr_chi_Ncn, // In
            Field8 Ncnm_val_denom_0                                 // In
) {
  // Return independent owned storage through the same direct-output kernel.
  using Result = Kokkos::View<typename Field0::non_const_data_type, Kokkos::LayoutLeft,
                              typename Field0::memory_space>;
  Result result;
  if constexpr (Field0::rank == 1)
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncm_to_Ncnm::result"),
                    mu_chi_1.extent(0));
  else
    result = Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Ncm_to_Ncnm::result"),
                    mu_chi_1.extent(0), mu_chi_1.extent(1));
  Ncm_to_Ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac, Ncm, const_Ncnp2_on_Ncnm2,
              const_corr_chi_Ncn, Ncnm_val_denom_0, result);
  return result;
}
} // namespace clubb_kokkos
