// Frequently used model constants from src/CLUBB_core/constants_clubb.F90.
// Standalone values preserve native precision and expression order. CLUBB_CAM
// obtains its host constants separately and is outside this standalone port.
#pragma once

#include "clubb_real.hpp"
#include <algorithm>
#include <limits>

namespace clubb_kokkos {

//-----------------------------------------------------------------------------
// Numerical/Arbitrary Constants
//-----------------------------------------------------------------------------

// Fortran file unit I/O constants
inline constexpr int fstderr = 0;
inline constexpr int fstdin  = 5;
inline constexpr int fstdout = 6;

// Maximum variable name length in CLUBB GrADS or netCDF output
inline constexpr int var_length = 30;

// Number of neighboring points to draw from in the hole filling algorithm
inline constexpr int num_hf_draw_points = 2;

// The parameter parab_cyl_max_input is the largest magnitude that the input to
// the parabolic cylinder function is allowed to have.  When the value of the
// input to the parabolic cylinder function is too large in magnitude
// (depending on the order of the parabolic cylinder function), overflow
// occurs, and the output of the parabolic cylinder function is +/-Inf.  The
// parameter parab_cyl_max_input places a limit on the absolute value of the
// input to the parabolic cylinder function.  When the value of the potential
// input exceeds this parameter (usually due to a very large ratio of ith PDF
// component mean of x to ith PDF component standard deviation of x), the
// variable x is considered to be constant and a different version of the
// equation called.
//
// The largest allowable magnitude of the input to the parabolic cylinder
// function (before overflow occurs) is dependent on the order of parabolic
// cylinder function.  However, after a lot of testing, it was determined that
// an absolute value of 49 works well for an order of 12 or less.
inline constexpr real_t parab_cyl_max_input =
    real_t(49.0); // Largest allowable input to parab. cyl. fnct.

// "Over-implicit" weighted time step.
//
// The weight of the implicit portion of a term is controlled by the factor
// gamma_over_implicit_ts (abbreviated "gamma" in the expression below).  A
// factor is added to the right-hand side of the equation in order to balance a
// weight that is not equal to 1, such that:
//
//      -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
//
// where X is the variable that is being solved for in a predictive equation
// (such as w'^3, w'th_l', r_t'^2, etc), y(t) is the linearized portion of the
// term that gets treated implicitly, and RHS is the portion of the term that
// is always treated explicitly.  A weight of greater than 1 can be applied to
// make the term more numerically stable.
//
//    gamma_over_implicit_ts          Effect on term
//
//            0.0               Term becomes completely explicit
//
//            1.0               Standard implicit portion of the term;
//                              as it was without the weighting factor.
//
//            1.5               Strongly weighted implicit portion of the term;
//                              increased numerical stability.
//
//            2.0               More strongly weighted implicit portion of the
//                              term; increased numerical stability.
//
// Note:  The "over-implicit" weighted time step is only applied to terms that
//        tend to significantly decrease the amount of numerical stability for
//        variable X.
//        The "over-implicit" weighted time step is applied to the turbulent
//        advection term for the following variables:
//           w'^3 (also applied to the turbulent production term), found in
//           module advance_wp2_wp3_module;
//           w'r_t', w'th_l', and w'sclr', found in
//           module advance_xm_wpxp_module; and
//           r_t'^2, th_l'^2, r_t'th_l', u'^2, v'^2, sclr'^2, sclr'r_t',
//           and sclr'th_l', found in module advance_xp2_xpyp_module.
inline constexpr real_t gamma_over_implicit_ts = real_t(1.50);

//-----------------------------------------------------------------------------
// Mathematical Constants
//-----------------------------------------------------------------------------
inline constexpr double pi_dp = 3.14159265358979323846;

inline constexpr real_t pi       = real_t(3.141592654); // The ratio of radii to their circumference
inline constexpr real_t invrs_pi = real_t(0.31830988618); // 1 / pi
inline constexpr real_t radians_per_deg = pi / real_t(180.);

inline constexpr double radians_per_deg_dp = pi_dp / 180.;

inline constexpr double sqrt_2pi_dp = 2.5066282746310005024; // sqrt(2*pi)
inline constexpr double sqrt_2_dp   = 1.4142135623730950488; // sqrt(2)

inline constexpr real_t sqrt_2pi = real_t(2.5066282746310005024); // sqrt(2*pi)
inline constexpr real_t sqrt_2   = real_t(1.4142135623730950488); // sqrt(2)

inline constexpr double two_dp        = 2.0;  // 2
inline constexpr double one_dp        = 1.0;  // 1
inline constexpr double one_half_dp   = 0.5;  // 1/2
inline constexpr double one_fourth_dp = 0.25; // 1/4
inline constexpr double zero_dp       = 0.0;  // 0

inline constexpr real_t one_hundred     = real_t(100.0);             // 100
inline constexpr real_t fifty           = real_t(50.0);              // 50
inline constexpr real_t thirty_six      = real_t(36.0);              // 36
inline constexpr real_t twenty          = real_t(20.0);              // 20
inline constexpr real_t eighteen        = real_t(18.0);              // 18
inline constexpr real_t twelve          = real_t(12.0);              // 12
inline constexpr real_t ten             = real_t(10.0);              // 10
inline constexpr real_t nine            = real_t(9.0);               // 9
inline constexpr real_t eight           = real_t(8.0);               // 8
inline constexpr real_t seven           = real_t(7.0);               // 7
inline constexpr real_t six             = real_t(6.0);               // 6
inline constexpr real_t five            = real_t(5.0);               // 5
inline constexpr real_t four            = real_t(4.0);               // 4
inline constexpr real_t three           = real_t(3.0);               // 3
inline constexpr real_t two             = real_t(2.0);               // 2
inline constexpr real_t three_halves    = real_t(3.0) / real_t(2.0); // 3/2
inline constexpr real_t four_thirds     = real_t(4.0) / real_t(3.0); // 4/3
inline constexpr real_t one             = real_t(1.0);               // 1
inline constexpr real_t three_fourths   = real_t(0.75);              // 3/4
inline constexpr real_t two_thirds      = real_t(2.0) / real_t(3.0); // 2/3
inline constexpr real_t one_half        = real_t(0.5);               // 1/2
inline constexpr real_t one_third       = real_t(1.0) / real_t(3.0); // 1/3
inline constexpr real_t one_fourth      = real_t(0.25);              // 1/4
inline constexpr real_t three_sixteenth = real_t(0.1875);            // 3/16
inline constexpr real_t one_eighth      = real_t(0.125);             // 1/8
inline constexpr real_t zero            = real_t(0.0);               // 0

//-----------------------------------------------------------------------------
// Physical constants
//-----------------------------------------------------------------------------

inline constexpr real_t Cp = real_t(1004.67); // Dry air specific heat at constant p [J/kg/K]
inline constexpr real_t Lv = real_t(2.5e6);   // Latent heat of vaporization         [J/kg]
inline constexpr real_t Ls = real_t(2.834e6); // Latent heat of sublimation          [J/kg]
inline constexpr real_t Lf = real_t(3.33e5);  // Latent heat of fusion               [J/kg]
inline constexpr real_t Rd = real_t(287.04);  // Dry air gas constant                [J/kg/K]
inline constexpr real_t Rv = real_t(461.5);   // Water vapor gas constant            [J/kg/K]

inline constexpr real_t stefan_boltzmann =
    real_t(5.6704e-8); // Stefan-Boltzmann constant [W/(m^2 K^4)]

inline constexpr real_t T_freeze_K = real_t(273.15); // Freezing point of water [K]

// Useful combinations of Rd and Rv
inline constexpr real_t ep  = Rd / Rv;                 // ep  = 0.622_core_rknd  [-]
inline constexpr real_t ep1 = (real_t(1.0) - ep) / ep; // ep1 = 0.61_core_rknd   [-]
inline constexpr real_t ep2 = real_t(1.0) / ep;        // ep2 = 1.61_core_rknd   [-]

inline constexpr real_t kappa = Rd / Cp; // kappa        [-]

// Changed g to grav to make it easier to find in the code 5/25/05
// real, parameter, public :: grav  = 9.80665_core_rknd ! Gravitational acceleration [m/s^2]
inline constexpr real_t omega_planet = real_t(7.292e-5); // Planetary rotation rate [s^-1]
inline constexpr real_t grav         = real_t(9.81);     // Gravitational acceleration     [m/s^2]
inline constexpr real_t p0           = real_t(1.0e5);    // Reference pressure             [Pa]

// Von Karman's constant
// Constant of the logarithmic wind profile in the surface layer
inline constexpr real_t vonk   = real_t(0.4);    // Accepted value is 0.40 (+/-) 0.01      [-]
inline constexpr real_t rho_lw = real_t(1000.0); // Density of liquid water                [kg/m^3]

inline constexpr real_t rho_ice = real_t(917.0); // Density of ice      [kg/m^3]

// Tolerances below which we consider moments to be zero
inline constexpr real_t w_tol   = real_t(2.e-2); // [m/s]
inline constexpr real_t thl_tol = real_t(1.e-2); // [K]
inline constexpr real_t rt_tol =
    std::max(real_t(1.e-8), std::numeric_limits<real_t>::epsilon()); // [kg/kg]
inline constexpr real_t chi_tol =
    std::max(real_t(1.e-8), std::numeric_limits<real_t>::epsilon()); // [kg/kg]
inline constexpr real_t eta_tol = chi_tol;                           // [kg/kg]

// Tolerances for use by the monatonic flux limiter.
// rt_tol_mfl is larger than rt_tol. rt_tol is extremely small
// (1e-8) to prevent spurious cloud formation aloft in LBA.
// rt_tol_mfl is larger (1e-4) to prevent the mfl from
// depositing moisture at the top of the domain.
inline constexpr real_t thl_tol_mfl = real_t(0.2);   // [K]
inline constexpr real_t rt_tol_mfl  = real_t(1.e-4); // [kg/kg]

// The tolerance for w'^2 is the square of the tolerance for w.
inline constexpr real_t w_tol_sqd = (w_tol * w_tol); // [m^2/s^2]

// Set tolerances for Khairoutdinov and Kogan rain microphysics to insure
// against numerical errors.  The tolerance values for Nc, rr, and Nr insure
// against underflow errors in computing the PDF for l_kk_rain.  Basically,
// they insure that those values squared won't be less then 10^-38, which is
// the lowest number that can be numerically represented.  However, the
// tolerance value for rc doubles as the lowest mixing ratio there can be to
// still officially have a cloud at that level.  This is figured to be about
// 1.0_core_rknd x 10^-7 kg/kg.  Brian; February 10, 2007.
inline constexpr real_t rc_tol  = real_t(1.0E-6); // Tolerance value for r_c  [kg/kg]
inline constexpr real_t Nc_tol  = real_t(1.0E+2); // Tolerance value for N_c  [#/kg]
inline constexpr real_t Ncn_tol = real_t(1.0E+2); // Tolerance value for N_cn [#/kg]

inline constexpr real_t mvr_cloud_max = real_t(1.6E-5); // Max. avg. mean vol. rad. cloud    [m]

inline constexpr real_t Nc_in_cloud_min = real_t(2.0e+4);

// Precipitating hydrometeor tolerances for mixing ratios.
inline constexpr real_t rr_tol = real_t(1.0E-10); // Tolerance value for r_r [kg/kg]
inline constexpr real_t ri_tol = real_t(1.0E-10); // Tolerance value for r_i [kg/kg]
inline constexpr real_t rs_tol = real_t(1.0E-10); // Tolerance value for r_s [kg/kg]
inline constexpr real_t rg_tol = real_t(1.0E-10); // Tolerance value for r_g [kg/kg]

// Maximum allowable values for the average mean volume radius of the various
// hydrometeor species.
inline constexpr real_t mvr_rain_max    = real_t(5.0E-3); // Max. avg. mean vol. rad. rain    [m]
inline constexpr real_t mvr_ice_max     = real_t(1.3E-4); // Max. avg. mean vol. rad. ice     [m]
inline constexpr real_t mvr_snow_max    = real_t(1.0E-2); // Max. avg. mean vol. rad. snow    [m]
inline constexpr real_t mvr_graupel_max = real_t(2.0E-2); // Max. avg. mean vol. rad. graupel [m]

// Precipitating hydrometeor tolerances for concentrations.
// Tolerance value for N_r [#/kg]
inline constexpr real_t Nr_tol =
    (one / (four_thirds * pi * rho_lw * (mvr_rain_max * mvr_rain_max * mvr_rain_max))) * rr_tol;

// Tolerance value for N_i [#/kg]
inline constexpr real_t Ni_tol =
    (one / (four_thirds * pi * rho_ice * (mvr_ice_max * mvr_ice_max * mvr_ice_max))) * ri_tol;

// Tolerance value for N_s [#/kg]
inline constexpr real_t Ns_tol =
    (one / (four_thirds * pi * rho_ice * (mvr_snow_max * mvr_snow_max * mvr_snow_max))) * rs_tol;

// Tolerance value for N_s [#/kg]
inline constexpr real_t Ng_tol =
    (one / (four_thirds * pi * rho_ice * (mvr_graupel_max * mvr_graupel_max * mvr_graupel_max))) *
    rg_tol;

// Minimum value for em (turbulence kinetic energy)
// With anisotropic TKE, em = (1/2) * ( up2 + vp2 + wp2 ).  Since
// up2, vp2, and wp2 all have the same minimum threshold value of
// w_tol_sqd, em cannot be less than (3/2) * w_tol_sqd.
inline constexpr real_t em_min = real_t(1.5) * w_tol_sqd; // [m^2/s^2]

inline constexpr real_t eps =
    std::max(real_t(1.0e-10),
             std::numeric_limits<real_t>::epsilon()); // Small value to prevent a divide by zero
// max() and epsilon() used here to raise
// threshold when using single precision

inline constexpr real_t max_num_stdevs =
    real_t(5.0); // Range of standard deviations for statistical significance

inline constexpr real_t zero_threshold =
    real_t(0.0); // Defining a threshold on a physical quantity to be 0.

// The maximum absolute value (or magnitude) that a correlation is allowed to
// have.  Statistically, a correlation is not allowed to be less than -1 or
// greater than 1, so the maximum magnitude would be 1.
inline constexpr real_t max_mag_correlation      = real_t(0.99); // Most correlations use this
inline constexpr real_t max_mag_correlation_flux = real_t(0.99); // Special for wprtp and wpthlp

// "base" smoothing magnitude before scaling for the respective data structure.
// See https://github.com/larson-group/clubb/issues/965#issuecomment-1119816722
// for a plot on how output behaves with varying min_max_smth_mag
inline constexpr real_t min_max_smth_mag = real_t(1.0e-9);

inline constexpr real_t cloud_frac_min = real_t(0.005); // Threshold for cloud fractions

inline constexpr real_t wp2_max = real_t(1000.); // Max value for variance clipping

//-----------------------------------------------------------------------------
// Useful conversion factors.
//-----------------------------------------------------------------------------
inline constexpr real_t sec_per_day = real_t(86400.0); // Seconds in a day.
inline constexpr real_t sec_per_hr  = real_t(3600.0);  // Seconds in an hour.
inline constexpr real_t sec_per_min = real_t(60.0);    // Seconds in a minute.
inline constexpr real_t min_per_hr  = real_t(60.0);    // Minutes in an hour.

inline constexpr real_t g_per_kg = real_t(1000.0); // Grams in a kilogram.

inline constexpr real_t pascal_per_mb = real_t(100.0); // Pascals per Millibar

inline constexpr real_t cm3_per_m3   = real_t(1.e6);  // Cubic centimeters per cubic meter
inline constexpr real_t micron_per_m = real_t(1.e6);  // Micrometers per meter
inline constexpr real_t cm_per_m     = real_t(100.);  // Centimeters per meter
inline constexpr real_t mm_per_m     = real_t(1000.); // Millimeters per meter

//-----------------------------------------------------------------------------
// Unused variable
//-----------------------------------------------------------------------------
inline constexpr real_t unused_var = -real_t(999.); // The standard value for unused variables

} // namespace clubb_kokkos
