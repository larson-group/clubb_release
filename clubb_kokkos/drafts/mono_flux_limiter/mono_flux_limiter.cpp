// Bound turbulent fluxes and compensate the affected mean fields.
// Native derivations/routine order are retained; see mono_flux_limiter.hpp.
#include "mono_flux_limiter.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp" // l_force_descending_solves must come from its existing owner.
#include "mean_adv.hpp"    // Reuse stable owner; promotion pending.
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
namespace clubb_kokkos {
namespace {
// These selector values must match the Native advance_xm_wpxp_module owner.
constexpr int mono_flux_thlm = 1, mono_flux_rtm = 2, mono_flux_um = 4, mono_flux_vm = 5;
constexpr int ndiags3 = 3;
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// dt: Model timestep length                      [s]
// wm_zt: w wind component on thermodynamic levels   [m/s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// grid_dir: Grid direction multiplier (+1 = ascending; -1 = descending)
// lhs: Left hand side of tridiagonal matrix
void mfl_xm_lhs(int nzm, int nzt, int ngrdcol, real_t dt, Real3D weights_zt2zm,         // In
                Real2D invrs_dzt, Real2D invrs_dzm,                                     // In
                Real2D wm_zt, bool l_implemented, bool l_upwind_xm_ma, real_t grid_dir, // In
                Real3D lhs)                                                             // Out
    ;
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// dt: Model timestep length                    [s]
// xm_old: xm; timestep (t) (thermodynamic levels)  [units vary]
// wpxp: w'x'; timestep (t+1); limited (m-levs.)  [units vary]
// xm_forcing: xm forcings (thermodynamic levels)       [units vary]
// invrs_dzt: The inverse spacing between momentum grid levels; centered over thermodynamic grid
// levels.
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// rhs: Right hand side of tridiagonal matrix equation
void mfl_xm_rhs(int nzm, int nzt, int ngrdcol, real_t dt, Real2D xm_old, Real2D wpxp, // In
                Real2D xm_forcing,                                                    // In
                Real2D invrs_dzt, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,           // In
                Real2D rhs)                                                           // Out
    ;
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// solve_type: Variables being solved for.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs: Left hand side of tridiagonal matrix
// rhs: Right hand side of tridiagonal matrix equation
// err_info: err_info struct containing err_code and err_header
// xm: Value of variable being solved for at timestep (t+1)   [units vary]
void mfl_xm_solve(int nzt, int ngrdcol, const grid &gr, int solve_type, // In
                  int tridiag_solve_method,                             // In
                  bool l_implemented,                                   // In
                  Real3D lhs, Real2D rhs, err_info_type &err_info,      // InOut
                  Real2D xm)                                            // Out
    ;
// nzm: Number of momentum levels [-]
// ngrdcol: Number of columns [-]
// w_1_zm: Mean w (1st PDF component)                   [m/s]
// w_2_zm: Mean w (2nd PDF component)                   [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)            [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)            [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// w_ref: Reference velocity, w|_ref (normally = 0)   [m/s]
// w_min: Minimum velocity to affect adjacent level    [m/s]
// stats: Native statistics context, mutated by sampled budget calls [-]
// mean_w_down: Overall mean w (<= w|_ref)                  [m/s]
// mean_w_up: Overall mean w (>= w|_ref)                  [m/s]
void mean_vert_vel_up_down(int nzm, int ngrdcol, // In
                           Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm,
                           Real2D varnce_w_2_zm,                            // In
                           Real2D mixt_frac_zm, real_t w_ref, Real2D w_min, // In
                           f2k::stats_type &stats,                          // InOut
                           Real2D mean_w_down, Real2D mean_w_up)            // Out
    ;
// nzm: Number of momentum levels [-]
// ngrdcol: Number of columns [-]
// w_i_zm: Mean of w                                   [m/s]
// varnce_w_i_zm: Variance of w                               [m^2/s^2]
// w_ref: Reference velocity, w|_ref (normally = 0)   [m/s]
// w_min: Min. velocity required to affect adj. level [m/s]
// mean_w_down_i: Mean w (<= w|_ref) from normal distribution [m/s]
// mean_w_up_i: Mean w (>= w|_ref) from normal distribution [m/s]
void calc_mean_w_up_down_component(int nzm, int ngrdcol,                     // In
                                   Real2D w_i_zm, Real2D varnce_w_i_zm,      // In
                                   real_t w_ref, Real2D w_min,               // In
                                   Real2D mean_w_down_i, Real2D mean_w_up_i) // Out
    ;
} // namespace

//=============================================================================
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// solve_type: Variables being solved for.
// dt: Model timestep length                           [s]
// xm_old: xm at previous time step (thermo. levs.) [units vary]
// xp2: x'^2 (momentum levels)                   [units vary]
// wm_zt: w wind component on thermodynamic levels [m/s]
// xm_forcing: xm forcings (thermodynamic levels)       [units vary]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// xp2_threshold: Lower limit of x'^2                      [units vary]
// xm_tol: Lower limit of maxdev                    [units vary]
// l_implemented: True if CLUBB is being implemented and run in a host model
// low_lev_effect: Index of lowest level that has an effect (for lev. k)
// high_lev_effect: Index of highest level that has an effect (for lev. k)
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// stats: Native statistics context, mutated by sampled budget calls [-]
// xm: xm at current time step (thermodynamic levels)  [units vary]
// wpxp: w'x' (momentum levels)                          [units vary]
// err_info: err_info struct containing err_code and err_header
void monotonic_turbulent_flux_limit(int nzm, int nzt, int ngrdcol, const grid &gr, int solve_type,
                                    real_t dt, // In
                                    Real2D xm_old, Real2D xp2, Real2D wm_zt,
                                    Real2D xm_forcing,                                        // In
                                    Real2D rho_ds_zm, Real2D rho_ds_zt,                       // In
                                    Real2D invrs_rho_ds_zm, Real2D invrs_rho_ds_zt,           // In
                                    real_t xp2_threshold, real_t xm_tol, bool l_implemented,  // In
                                    Kokkos::View<int **, Kokkos::LayoutLeft> low_lev_effect,  // In
                                    Kokkos::View<int **, Kokkos::LayoutLeft> high_lev_effect, // In
                                    int tridiag_solve_method,                                 // In
                                    bool l_upwind_xm_ma,                                      // In
                                    bool l_mono_flux_lim_spikefix,                            // In
                                    f2k::stats_type &stats,                          // InOut
                                    Real2D xm, Real2D wpxp, err_info_type &err_info) // InOut
{
  // Description:
  // Limits the value of w'x' and corrects the value of xm when the xm turbulent
  // advection term is not monotonic.  A monotonic turbulent advection scheme
  // will not create new extrema for variable x, based only on turbulent
  // advection (not considering mean advection and xm forcings).
  //
  // Montonic turbulent advection
  // ----------------------------
  //
  // A monotonic turbulent advection scheme does not allow new extrema for
  // variable x to be created (by means of turbulent advection).  In a
  // monotonic turbulent advection scheme, when only the effects of turbulent
  // advection are considered (neglecting forcings and mean advection), the
  // value of variable x at a given point should not increase above the
  // greatest value of variable x at nearby points, nor decrease below the
  // smallest value of variable x at nearby points.  Nearby points are points
  // that are close enough to the given point so that the value of variable x
  // at the given point is effected by the values of variable x at the nearby
  // points by means of transfer by turbulent winds during a time step.  Again,
  // a monotonic scheme insures that advection only transfers around values of
  // variable x and does not create new extrema for variable x.  A monotonic
  // turbulent advection scheme is useful because the turbulent advection term
  // (w'x') may go numerically unstable, resulting in large instabilities in
  // the mean field (xm).  A monotonic turbulent advection scheme will limit
  // the change in xm, and also in w'x'.
  //
  // The following example illustrates the concept of monotonic turbulent
  // advection.  Three successive vertical grid levels are shown (k-1, k, and
  // k+1).  Three point values of theta-l are listed at every vertical grid
  // level.  All three vertical levels have a mean theta-l (thlm) of 288.0 K.
  // A circulation is occuring (in the direction of the arrows) in the vertical
  // (w wind component) and in the horizontal (u and/or v wind components),
  // such that the mean value of vertical velocity (wmm) is 0, but there is a
  // turbulent component such that w'^2 > 0.
  //
  // level = k+1 || --- 287.0 K --- 288.0 K --- 289.0 K --- || thlm = 288.0
  //             ||      / \--------------------->|         ||
  //             ||       |                       |         || wmm = 0; wp2 > 0
  //             ||       |<---------------------\ /        ||
  // level = k   || --- 288.0 K --- 288.0 K --- 288.0 K --- || thlm = 288.0
  //             ||       |<---------------------/ \        ||
  //             ||       |                       |         || wmm = 0; wp2 > 0
  //             ||      \ /--------------------->|         ||
  // level = k-1 || --- 287.5 K --- 288.0 K --- 288.5 K --- || thlm = 288.0
  //
  // Neglecting any contributions from thlm forcings (effects of radiation,
  // microphysics, large-scale horizontal advection, etc.), the values of
  // theta-l as shown will be altered by only turbulent advection.  As a side
  // note, the contribution of mean advection will be 0 since wmm = 0.  The
  // diagram shows that the value of theta-l at the point on the right at level
  // k will increase.  However, the values of theta-l at the other two points
  // at level k will remain the same.  Thus, the value of thlm at level k will
  // become greater than 288.0 K.  In the same manner, the values of thlm at
  // the other two vertical levels (k-1 and k+1) will become smaller than
  // 288.0 K.  However, the monotonic turbulent advection scheme insures that
  // any theta-l point value cannot become smaller than the smallest theta-l
  // point value (287.0 K) or larger than the largest theta-l point value
  // (289.0 K).  Since all theta-l point values must fall between 287.0 K and
  // 289.0 K, the level averages of theta-l (thlm) must fall between 287.0 K
  // and 289.0 K.  Thus, any values of the turbulent flux, w'th_l', that would
  // cause thlm to rise above 289.0 K or fall below 287.0 K, not considering
  // the effect of other terms on thlm (such as forcings), are faulty and need
  // to be limited appropriately.  The values of thlm also need to be corrected
  // appropriately.
  //
  // Formula for the limitation of w'x' and xm
  // -----------------------------------------
  //
  // The equation for change in the mean field, xm, over time is:
  //
  // d(xm)/dt = -w*d(xm)/dz - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing;
  //
  // where w*d(xm)/dz is the mean advection component,
  // (1/rho_ds) * d( rho_ds * w'x' )/dz is the turbulent advection component,
  // and xm_forcing is the xm forcing component.  The d(xm)/dt time tendency
  // component is discretized as:
  //
  // xm(k,<t+1>)/dt = xm(k,<t>)/dt - w*d(xm)/dz
  //                  - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing.
  //
  // The value of xm after it has been advanced to timestep (t+1) must be in an
  // appropriate range based on the values of xm at timestep (t), the amount of
  // xm forcings applied over the ensuing time step, and the amount of mean
  // advection applied over the ensuing time step.  This is exactly the same
  // thing as saying that the value of xm(k,<t+1>), with the contribution of
  // turbulent advection included, must fall into a certain range based on the
  // value of xm(k,<t+1>) without the contribution of the turbulent advection
  // component over the last time step.  The following inequality is used to
  // limit the value of xm(k,<t+1>):
  //
  // MIN{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
  //         - x_max_dev_low(k-1,<t>),
  //      xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //         - x_max_dev_low(k,<t>),
  //      xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
  //         - x_max_dev_low(k+1,<t>) }
  // <= xm(k,<t+1>) <=
  // MAX{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
  //         + x_max_dev_high(k-1,<t>),
  //      xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //         + x_max_dev_high(k,<t>),
  //      xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
  //         + x_max_dev_high(k+1,<t>) };
  //
  // where x_max_dev_low is the absolute value of the deviation from the mean
  // of the smallest point value of variable x at the given vertical level and
  // timestep; and where x_max_dev_high is the deviation from the mean of the
  // largest point value of variable x at the given vertical level and
  // timestep.  For example, at vertical level (k+1) and timestep (t):
  //
  // x_max_dev_low(k+1,<t>)  = | MIN( x(k+1,<t>) ) - xm(k+1,<t>) |;
  // x_max_dev_high(k+1,<t>) = MAX( x(k+1,<t>) ) - xm(k+1,<t>).
  //
  // The inequality shown above only takes into account values from the central
  // level, one-level-below the central level, and one-level-above the central
  // level.  This is the minimal amount of vertical levels that can have their
  // values taken into consideration.  Any vertical level that can have it's
  // properties advect to the given level during the course of a single time
  // step can be taken into consideration.  However, only three levels will be
  // considered in this example for the sake of simplicity.
  //
  // The inequality will be written in more simple terms:
  //
  // xm_lower_lim_allowable(k) <= xm(k,<t+1>) <= xm_upper_lim_allowable(k).
  //
  // The inequality can now be related to the turbulent flux, w'x'(k,<t+1>),
  // through a substitution that is made for xm(k,<t+1>), such that:
  //
  // xm(k,<t+1>) = xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //               - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k).
  //
  // The inequality becomes:
  //
  // xm_lower_lim_allowable(k)
  // <=
  //    xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //    - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k)
  // <=
  // xm_upper_lim_allowable(k).
  //
  // The inequality is rearranged, and the turbulent advection term,
  // d(w'x')/dz, is discretized:
  //
  // xm_lower_lim_allowable(k)
  // - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ]
  // <=
  //    - dt * (1/rho_ds_zt(k))
  //           * invrs_dzt(k)
  //             * [   rho_ds_zm(k+1) * w'x'(k+1,<t+1>)
  //                 - rho_ds_zm(k) * w'x'(k,<t+1>) ]
  // <=
  // xm_upper_lim_allowable(k)
  // - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ];
  //
  // where invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
  //
  // Multiplying the inequality by -rho_ds_zt(k)/(dt*invrs_dzt(k)):
  //
  // rho_ds_zt(k)/(dt*invrs_dzt(k))
  // * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //     - xm_lower_lim_allowable(k) ]
  // >=
  //    rho_ds_zm(k+1) * w'x'(k+1,<t+1>) - rho_ds_zm(k) * w'x'(k,<t+1>)
  // >=
  // rho_ds_zt(k)/(dt*invrs_dzt(k))
  // * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //     - xm_upper_lim_allowable(k) ].
  //
  // Note:  The inequality symbols have been flipped due to multiplication
  //        involving a (-) sign.
  //
  // Adding rho_ds_zm(k) * w'x'(k,<t+1>) to the inequality:
  //
  // rho_ds_zt(k)/(dt*invrs_dzt(k))
  // * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //     - xm_lower_lim_allowable(k) ]
  // + rho_ds_zm(k) * w'x'(k,<t+1>)
  // >= rho_ds_zm(k+1) * w'x'(k+1,<t+1>) >=
  // rho_ds_zt(k)/(dt*invrs_dzt(k))
  // * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //     - xm_upper_lim_allowable(k) ]
  // + rho_ds_zm(k) * w'x'(k,<t+1>).
  //
  // The inequality is then rearranged to be based around w'x'(k+1,<t+1>):
  //
  // (1/rho_ds_zm(k+1))
  // * [ rho_ds_zt(k)/(dt*invrs_dzt(k))
  //     * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //         - xm_lower_lim_allowable(k) }
  //     + rho_ds_zm(k) * w'x'(k,<t+1>) ]
  // >=   w'x'(k+1,<t+1>)   >=
  // (1/rho_ds_zm(k+1))
  // * [ rho_ds_zt(k)/(dt*invrs_dzt(k))
  //     * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
  //         - xm_upper_lim_allowable(k) }
  //     + rho_ds_zm(k) * w'x'(k,<t+1>) ].
  //
  // The values of w'x' are found on the momentum levels, while the values of
  // xm are found on the thermodynamic levels.  Additionally, the values of
  // rho_ds_zm are found on the momentum levels, and the values of rho_ds_zt
  // are found on the thermodynamic levels.  The inequality is applied to
  // w'x'(k+1,<t+1>) from vertical levels (k+1=) 2 through the second-highest
  // level, (k+1=) gr%nzm-1.  The value of w'x' at level 1 is a set surface
  // (or lowest level) flux.  The value of w'x' at the highest level is also a
  // set value, and therefore is not altered.
  //
  // Approximating maximum and minimum values of x at any given vertical level
  // -------------------------------------------------------------------------
  //
  // The CLUBB code provides means, variances, and covariances for certain
  // variables at all vertical levels.  However, there is no way to find the
  // maximum or minimum point value of any variable on any vertical level.
  // Without that information, x_max_dev_low and x_max_dev_high can't be found,
  // and the inequality above is useless.  However, there is a way to
  // approximate the maximum and minimum point values at any given vertical
  // level.  The maximum and minimum point values can be approximated through
  // the use of the variance, x'^2.
  //
  // Just as the mean value of x, which is xm, and the turbulent flux of x,
  // which is w'x', are known, so is the variance of x, which is x'^2.  The
  // standard deviation of x is the square root of the variance of x.  The
  // distribution of x along the horizontal plane (at vertical level k) is
  // approximated to be the sum of two normal (or Gaussian) distributions.
  // Most of the values in a normal distribution are found within 2 standard
  // deviations from the mean.  Thus, the maximum point value of x along the
  // horizontal plance at any vertical level can be approximated as:
  // xm + 2*sqrt(x'^2).  Likewise, the minimum value of x along the horizontal
  // plane at any vertical level can be approximated as:  xm - 2*sqrt(x'^2).
  //
  // The values of x'^2 are found on the momentum levels.  The values of xm
  // are found on the thermodynamic levels.  Thus, the values of x'^2 are
  // interpolated to the thermodynamic levels in order to find the maximum
  // and minimum point values of variable x.
  //
  // The one downfall of this method is that instabilities can arise in the
  // model where unphysically large values of x'^2 are produced.  Thus, this
  // allows for an unphysically large deviation of xm from its values at the
  // previous time step due to turbulent advection.  Thus, for purposes of
  // determining the maximum and minimum point values of x, a upper limit
  // is placed on x'^2, in order to limit the standard deviation of x.  This
  // limit is only applied in this subroutine, and is not applied to x'^2
  // elsewhere in the model code.
  // References:
  //-----------------------------------------------------------------------
  // Constant Parameters
  // Flag for using a semi-implicit, tridiagonal method to solve for xm(t+1)
  // when xm(t+1) needs to be changed.
  //----------------------- Input Variables -----------------------
  // eliminates spurious drying tendencies at model top
  //----------------------- Input/Output Variables -----------------------
  //------------------------- Inout Variables -------------------------
  //----------------------- Local Variables -----------------------
  //    integer, parameter :: &
  //      num_levs = 10  ! Number of levels above and below level k to look for
  //                     ! maxima and minima of variable x.
  constexpr bool l_mfl_xm_imp_adj = true;
  Real2D xp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "monotonic_turbulent_flux_limit::xp2_zt"),
      ngrdcol, nzt); // x'^2 interpolated to thermodynamic levels  [units vary]
  Real2D xm_enter_mfl(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                         "monotonic_turbulent_flux_limit::xm_enter_mfl"),
                      ngrdcol, nzt); // xm as it enters the MFL                    [units vary]
  Real2D xm_without_ta(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                          "monotonic_turbulent_flux_limit::xm_without_ta"),
                       ngrdcol, nzt); // Value of xm without turb. adv. contrib.    [units vary]
  Real2D wpxp_net_adjust(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                            "monotonic_turbulent_flux_limit::wpxp_net_adjust"),
                         ngrdcol, nzm); // Net amount of adjustment needed on w'x'    [units vary]
  Real2D min_x_allowable_lev(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::min_x_allowable_lev"),
      ngrdcol, nzt); // Smallest usuable value of x at lev k [units vary]
  Real2D max_x_allowable_lev(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::max_x_allowable_lev"),
      ngrdcol, nzt); // Largest usuable value of x at lev k  [units vary]
  Real2D min_x_allowable(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                            "monotonic_turbulent_flux_limit::min_x_allowable"),
                         ngrdcol, nzt); // Smallest usuable x within k +/- num_levs [units vary]
  Real2D max_x_allowable(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                            "monotonic_turbulent_flux_limit::max_x_allowable"),
                         ngrdcol, nzt); // Largest usuable x within k +/- num_levs  [units vary]
  Real2D wpxp_mfl_max(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                         "monotonic_turbulent_flux_limit::wpxp_mfl_max"),
                      ngrdcol, nzm); // Upper limit on w'x'(k)                [units vary]
  Real2D wpxp_mfl_min(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                         "monotonic_turbulent_flux_limit::wpxp_mfl_min"),
                      ngrdcol, nzm); // Lower limit on w'x'(k)                [units vary]
  real_t max_xp2;                    // Maximum allowable x'^2                        [units vary]
  real_t max_dev;                    // Determines approximate upper/lower limit of x [units vary]
  real_t m_adv_term;                 // Contribution of mean advection to d(xm)/dt    [units vary]
  real_t xm_density_weighted;        // Density weighted xm at domain top             [units vary]
  real_t xm_adj_coef;                // Coeffecient to eliminate spikes at domain top [units vary]
  real_t xm_vert_integral;           // Vertical integral of xm                       [units_vary]
  real_t dxm_dt_mfl_adjust;          // Rate of change of adjustment to xm            [units vary]
  real_t dz;                         // zm grid spacing at top of domain              [m]
  Real3D lhs_mfl_xm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "monotonic_turbulent_flux_limit::lhs_mfl_xm"),
      ndiags3, ngrdcol, nzt); // Left hand side of tridiagonal matrix
  Real2D rhs_mfl_xm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "monotonic_turbulent_flux_limit::rhs_mfl_xm"),
      ngrdcol, nzt); // Right hand side of tridiagonal matrix equation
  int k;
  int i;
  int j;        // Array indices
  int low_lev;  // Lowest level (from level k) to look for x minima and maxima
  int high_lev; // Highest level (from level k) to look for x minima and maxima
  Bool1D l_xm_adjustment_needed(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::l_xm_adjustment_needed"),
      ngrdcol); // Indicates if we need an adjustment for a column
  bool l_any_adjustment_needed;
  Real2D xm_mfl(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "monotonic_turbulent_flux_limit::xm_mfl"),
      ngrdcol, nzt);
  real_t max_tmp;
  real_t min_tmp;
  real_t invrs_dt;
  Real2D wpxp_mfl_max_term_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::wpxp_mfl_max_term_zt"),
      ngrdcol, nzt);
  Real2D wpxp_mfl_min_term_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::wpxp_mfl_min_term_zt"),
      ngrdcol, nzt);
  Real2D wpxp_thresh_term_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "monotonic_turbulent_flux_limit::wpxp_thresh_term_zt"),
      ngrdcol, nzt);
  Real2D wpxp_thresh_term(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                             "monotonic_turbulent_flux_limit::wpxp_thresh_term"),
                          ngrdcol, nzm);
  int k_zt;
  std::string name_wpxp_mfl;
  std::string name_xm_mfl;
  std::string name_xm_enter_mfl;
  std::string name_xm_old;
  std::string name_wpxp_enter_mfl;
  std::string name_xm_without_ta;
  std::string name_xm_mfl_min;
  std::string name_xm_mfl_max;
  std::string name_wpxp_mfl_min;
  std::string name_wpxp_mfl_max;
  std::string name_xm_exit_mfl;
  std::string name_wpxp_exit_mfl;
  int adjustment_count = 0;
  const auto gr_dzt = gr.dzt;
  const auto gr_grid_dir = gr.grid_dir;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_invrs_dzt = gr.invrs_dzt;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_weights_zt2zm = gr.weights_zt2zm;
  const auto gr_zm = gr.zm;
  //---------------------------- Begin Code ----------------------------
  // These lines synchonize the GPU argument values to the CPU argument values, allowing us
  // to test the divergence of this procedure only, rather than the divergence of the whole code
  // MONOFLUX TEST !$acc update device( wpxp, xm, xm_old, wm_zt, &
  // MONOFLUX TEST !$acc  xm_forcing, rho_ds_zt, invrs_rho_ds_zt, &
  // MONOFLUX TEST !$acc  xp2, rho_ds_zm, invrs_rho_ds_zm, &
  // MONOFLUX TEST !$acc  low_lev_effect, high_lev_effect )
  name_wpxp_mfl = "";
  name_xm_mfl = "";
  name_xm_enter_mfl = "";
  name_xm_old = "";
  name_wpxp_enter_mfl = "";
  name_xm_without_ta = "";
  name_xm_mfl_min = "";
  name_xm_mfl_max = "";
  name_wpxp_mfl_min = "";
  name_wpxp_mfl_max = "";
  name_xm_exit_mfl = "";
  name_wpxp_exit_mfl = "";

  switch (solve_type) {
  case mono_flux_rtm: {
    max_xp2 = real_t(5.0e-6);
    name_wpxp_mfl = "wprtp_mfl";
    name_xm_mfl = "rtm_mfl";
    name_xm_enter_mfl = "rtm_enter_mfl";
    name_xm_old = "rtm_old";
    name_wpxp_enter_mfl = "wprtp_enter_mfl";
    name_xm_without_ta = "rtm_without_ta";
    name_xm_mfl_min = "rtm_mfl_min";
    name_xm_mfl_max = "rtm_mfl_max";
    name_wpxp_mfl_min = "wprtp_mfl_min";
    name_wpxp_mfl_max = "wprtp_mfl_max";
    name_xm_exit_mfl = "rtm_exit_mfl";
    name_wpxp_exit_mfl = "wprtp_exit_mfl";
    break;
  }
  case mono_flux_thlm: {
    max_xp2 = real_t(5.0);
    name_wpxp_mfl = "wpthlp_mfl";
    name_xm_mfl = "thlm_mfl";
    name_xm_enter_mfl = "thlm_enter_mfl";
    name_xm_old = "thlm_old";
    name_wpxp_enter_mfl = "wpthlp_enter_mfl";
    name_xm_without_ta = "thlm_without_ta";
    name_xm_mfl_min = "thlm_mfl_min";
    name_xm_mfl_max = "thlm_mfl_max";
    name_wpxp_mfl_min = "wpthlp_mfl_min";
    name_wpxp_mfl_max = "wpthlp_mfl_max";
    name_xm_exit_mfl = "thlm_exit_mfl";
    name_wpxp_exit_mfl = "wpthlp_exit_mfl";
    break;
  }
  case mono_flux_um: {
    max_xp2 = real_t(10.0);
    name_wpxp_mfl = "upwp_mfl";
    name_xm_mfl = "um_mfl";
    break;
  }
  case mono_flux_vm: {
    max_xp2 = real_t(10.0);
    name_wpxp_mfl = "vpwp_mfl";
    name_xm_mfl = "vm_mfl";
    break;
  }
  default: {
    max_xp2 = real_t(5.0);
    break;
  }
  }

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_1(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                               "monotonic_turbulent_flux_limit::stats_values_1"),
                            wpxp.extent(0), wpxp.extent(1));
      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::stats_values", wpxp.size(), KOKKOS_LAMBDA(int cell) {
            const int i = cell % wpxp.extent(0), k = cell / wpxp.extent(0);
            stats_values_1(i, k) = wpxp(i, k) / dt;
          });
      f2k::stats_begin_budget(name_wpxp_mfl, stats_values_1, stats);
    }
    {
      Real2D stats_values_2(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                               "monotonic_turbulent_flux_limit::stats_values_2"),
                            xm.extent(0), xm.extent(1));
      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::stats_values", xm.size(), KOKKOS_LAMBDA(int cell) {
            const int i = cell % xm.extent(0), k = cell / xm.extent(0);
            stats_values_2(i, k) = xm(i, k) / dt;
          });
      f2k::stats_begin_budget(name_xm_mfl, stats_values_2, stats);
    }
    f2k::stats_update(name_xm_enter_mfl, xm, stats);
    f2k::stats_update(name_xm_old, xm_old, stats);
    f2k::stats_update(name_wpxp_enter_mfl, wpxp, stats);
  }
  invrs_dt = one / dt;
  const int loop_3_extent_0 = std::max(0, (nzm - (1)) / (1) + 1);
  const int loop_3_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "monotonic_turbulent_flux_limit::loop3", loop_3_extent_0 * loop_3_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_3_extent_1)) * (1);
        remaining /= (loop_3_extent_1);
        const int k = 1 + (remaining % (loop_3_extent_0)) * (1);
        remaining /= (loop_3_extent_0);

        // Initialize arrays
        wpxp_net_adjust((i)-1, (k)-1) = real_t(0.0);
        wpxp_mfl_min((i)-1, (k)-1) = real_t(0.0);
        wpxp_mfl_max((i)-1, (k)-1) = real_t(0.0);
      });
  // Interpolate x'^2 to thermodynamic levels.
  Kokkos::deep_copy(xp2_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, xp2));
  // Place an upper limit on xp2_zt.
  // For purposes of this subroutine, an upper limit has been placed on the
  // variance, x'^2.  This does not effect the value of x'^2 anywhere else in
  // the model code.  The upper limit is a reasonable upper limit.  This is
  // done to prevent unphysically large standard deviations caused by numerical
  // instabilities in the x'^2 profile.
  const int loop_4_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_4_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "monotonic_turbulent_flux_limit::loop4", loop_4_extent_0 * loop_4_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_4_extent_1)) * (1);
        remaining /= (loop_4_extent_1);
        const int k = 1 + (remaining % (loop_4_extent_0)) * (1);
        remaining /= (loop_4_extent_0);

        xp2_zt((i)-1, (k)-1) =
            ((((xp2_zt((i)-1, (k)-1)) > (xp2_threshold) ? (xp2_zt((i)-1, (k)-1))
                                                        : (xp2_threshold))) < (max_xp2)
                 ? (((xp2_zt((i)-1, (k)-1)) > (xp2_threshold) ? (xp2_zt((i)-1, (k)-1))
                                                              : (xp2_threshold)))
                 : (max_xp2));
        // Store the value of xm as it enters the mfl
        xm_enter_mfl((i)-1, (k)-1) = xm((i)-1, (k)-1);
      });
  // Find the maximum and minimum usuable values of variable x at each
  // vertical level.  This computation needs to be performed for all
  // vertical levels above the ground (or model surface).
  const int loop_5_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_5_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "monotonic_turbulent_flux_limit::loop5", loop_5_extent_0 * loop_5_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_5_extent_1)) * (1);
        remaining /= (loop_5_extent_1);
        const int k = 1 + (remaining % (loop_5_extent_0)) * (1);
        remaining /= (loop_5_extent_0);
        real_t max_dev;
        real_t m_adv_term;
        // Most values are found within +/- 2 standard deviations from the mean.
        // Use +/- 2 standard deviations from the mean as the maximum/minimum
        // values.
        // max_dev = 2.0_core_rknd*stnd_dev_x
        // Set a minimum on max_dev
        max_dev = ((real_t(2.0) * Kokkos::sqrt(xp2_zt((i)-1, (k)-1))) > (xm_tol)
                       ? (real_t(2.0) * Kokkos::sqrt(xp2_zt((i)-1, (k)-1)))
                       : (xm_tol));
        // Shut off to avoid using new, possibly corrupt mean advection term
        m_adv_term = real_t(0.0);
        // Find the value of xm without the contribution from the turbulent
        // advection term.
        // Note:  the contribution of xm_forcing at level gr%nzt should be 0.
        xm_without_ta((i)-1, (k)-1) =
            xm_old((i)-1, (k)-1) + dt * xm_forcing((i)-1, (k)-1) + dt * m_adv_term;
        // Find the minimum usuable value of variable x at each vertical level.

        if (solve_type != mono_flux_um && solve_type != mono_flux_vm) {
          // Since variable x must be one of theta_l, r_t, or a scalar, all of
          // which are positive definite quantities, the value must be >= 0.
          min_x_allowable_lev((i)-1, (k)-1) =
              ((xm_without_ta((i)-1, (k)-1) - max_dev) > (zero_threshold)
                   ? (xm_without_ta((i)-1, (k)-1) - max_dev)
                   : (zero_threshold));
        } else {
          // Variable x must be one of u or v.
          min_x_allowable_lev((i)-1, (k)-1) = xm_without_ta((i)-1, (k)-1) - max_dev;
        }
        // Find the maximum usuable value of variable x at each vertical level.
        max_x_allowable_lev((i)-1, (k)-1) = xm_without_ta((i)-1, (k)-1) + max_dev;
      });
  // Find the maximum and minimum usuable values of x that can effect the value
  // of x at level k.  Then, find the upper and lower limits of w'x'.  Reset
  // the value of w'x' if it is outside of those limits, and store the amount
  // of adjustment that was needed to w'x'.
  // The values of w'x' at level 1 and at level gr%nz are set values and
  // are not altered.
  // Find the smallest and largest value of all relevant levels for variable x.
  const int loop_6_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_6_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "monotonic_turbulent_flux_limit::loop6", loop_6_extent_0 * loop_6_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_6_extent_1)) * (1);
        remaining /= (loop_6_extent_1);
        const int k = 1 + (remaining % (loop_6_extent_0)) * (1);
        remaining /= (loop_6_extent_0);
        int j;
        int low_lev;
        int high_lev;
        real_t max_tmp;
        real_t min_tmp;

        if (gr_grid_dir_indx > 0) {
          // Ascending grid
          low_lev = ((low_lev_effect((i)-1, (k)-1)) > (gr_k_lb_zt) ? (low_lev_effect((i)-1, (k)-1))
                                                                   : (gr_k_lb_zt));
          high_lev =
              ((high_lev_effect((i)-1, (k)-1)) < (gr_k_ub_zt) ? (high_lev_effect((i)-1, (k)-1))
                                                              : (gr_k_ub_zt));
        } else {
          // Descending grid
          low_lev = ((low_lev_effect((i)-1, (k)-1)) < (gr_k_lb_zt) ? (low_lev_effect((i)-1, (k)-1))
                                                                   : (gr_k_lb_zt));
          high_lev =
              ((high_lev_effect((i)-1, (k)-1)) > (gr_k_ub_zt) ? (high_lev_effect((i)-1, (k)-1))
                                                              : (gr_k_ub_zt));
        }
        min_tmp = min_x_allowable_lev((i)-1, (low_lev)-1);
        max_tmp = max_x_allowable_lev((i)-1, (low_lev)-1);

        for (j = low_lev + gr_grid_dir_indx; (gr_grid_dir_indx) > 0 ? j <= high_lev : j >= high_lev;
             j += gr_grid_dir_indx) {
          min_tmp = ((min_tmp) < (min_x_allowable_lev((i)-1, (j)-1))
                         ? (min_tmp)
                         : (min_x_allowable_lev((i)-1, (j)-1)));
          max_tmp = ((max_tmp) > (max_x_allowable_lev((i)-1, (j)-1))
                         ? (max_tmp)
                         : (max_x_allowable_lev((i)-1, (j)-1)));
        }
        min_x_allowable((i)-1, (k)-1) = min_tmp;
        max_x_allowable((i)-1, (k)-1) = max_tmp;
      });
  const int loop_7_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_7_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "monotonic_turbulent_flux_limit::loop7", loop_7_extent_0 * loop_7_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_7_extent_1)) * (1);
        remaining /= (loop_7_extent_1);
        const int k = 1 + (remaining % (loop_7_extent_0)) * (1);
        remaining /= (loop_7_extent_0);

        wpxp_thresh_term_zt((i)-1, (k)-1) =
            invrs_dt * gr_grid_dir * gr_dzt((i)-1, (k)-1) *
            (xm_without_ta((i)-1, (k)-1) - min_x_allowable((i)-1, (k)-1));
        wpxp_mfl_max_term_zt((i)-1, (k)-1) =
            rho_ds_zt((i)-1, (k)-1) * wpxp_thresh_term_zt((i)-1, (k)-1);
        wpxp_mfl_min_term_zt((i)-1, (k)-1) =
            rho_ds_zt((i)-1, (k)-1) * invrs_dt * gr_grid_dir * gr_dzt((i)-1, (k)-1) *
            (xm_without_ta((i)-1, (k)-1) - max_x_allowable((i)-1, (k)-1));
      });
  // Interpolate wpxp_thresh_term_zt to momentum levels
  Kokkos::deep_copy(wpxp_thresh_term, zt2zm_api(nzm, nzt, ngrdcol, gr, wpxp_thresh_term_zt));
  l_any_adjustment_needed = false;
  // NOTE: This loop first checks if wpxp is outside of the range determined by wpxp_mfl_min
  //       and wpxp_mfl_max, and if even a single value is outside that range, we set
  //       l_any_adjustment_needed = .true. then perform the calculation later. This is done
  //       because the calculation contains a vertical dependency, so it cannot
  //       be parallelized, which is slow on GPUs.
  //       This initial check CAN be done in parallel though, hence the "collapse(2), so
  //       we run this quickly and only perform the slow version if we have to.
  const int loop_8_extent_0 = std::max(0, (ngrdcol - (1)) / (1) + 1);
  const int loop_8_extent_1 = std::max(
      0,
      (gr_k_ub_zm - gr_grid_dir_indx - (gr_k_lb_zm + gr_grid_dir_indx)) / (gr_grid_dir_indx) + 1);

  Kokkos::parallel_reduce(
      "monotonic_turbulent_flux_limit::loop8", loop_8_extent_0 * loop_8_extent_1,
      KOKKOS_LAMBDA(const int cell, int &adjustment_hits) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_8_extent_0)) * (1);
        remaining /= (loop_8_extent_0);
        const int k =
            gr_k_lb_zm + gr_grid_dir_indx + (remaining % (loop_8_extent_1)) * (gr_grid_dir_indx);
        remaining /= (loop_8_extent_1);
        int k_zt;
        // Calculate index needed for zt level variables based on grid direction.

        if (gr_grid_dir_indx > 0) {
          // Ascending grid
          k_zt = k - 1;
        } else {
          // Descending grid
          k_zt = k;
        }
        // Find the upper limit for w'x' for a monotonic turbulent flux.
        // The following "if" statement ensures there are no "spikes" at the top of the column,
        // which can cause unphysical rtm and thlm tendencies over the height of the column.
        // The fix essentially turns off the monotonic flux limiter for these special cases,
        // but tests show that it still performs well otherwise and runs stably.

        if (l_mono_flux_lim_spikefix && solve_type == mono_flux_rtm &&
            Kokkos::abs(wpxp((i)-1, (k - gr_grid_dir_indx) - 1)) >
                wpxp_thresh_term((i)-1, (k - gr_grid_dir_indx) - 1) &&
            wpxp((i)-1, (k - gr_grid_dir_indx) - 1) < real_t(0.0)) {
          wpxp_mfl_max((i)-1, (k)-1) = zero;
        } else {
          wpxp_mfl_max((i)-1, (k)-1) =
              invrs_rho_ds_zm((i)-1, (k)-1) * (wpxp_mfl_max_term_zt((i)-1, (k_zt)-1) +
                                               rho_ds_zm((i)-1, (k - gr_grid_dir_indx) - 1) *
                                                   wpxp((i)-1, (k - gr_grid_dir_indx) - 1));
        }

        if (wpxp((i)-1, (k)-1) > wpxp_mfl_max((i)-1, (k)-1)) {
          // wpxp is above wpxp_mfl_max, we will need to perform the slow vertically
          // dependent calculation to ensure correct values
          adjustment_hits = 1;
        } else {
          // Find the lower limit for w'x' for a monotonic turbulent flux.
          wpxp_mfl_min((i)-1, (k)-1) =
              invrs_rho_ds_zm((i)-1, (k)-1) * (wpxp_mfl_min_term_zt((i)-1, (k_zt)-1) +
                                               rho_ds_zm((i)-1, (k - gr_grid_dir_indx) - 1) *
                                                   wpxp((i)-1, (k - gr_grid_dir_indx) - 1));

          if (wpxp((i)-1, (k)-1) < wpxp_mfl_min((i)-1, (k)-1)) {
            // wpxp is below wpxp_mfl_min, we will need to perform the slow vertically
            // dependent calculation to ensure correct values
            adjustment_hits = 1;
          }
        }
      },
      adjustment_count);
  l_any_adjustment_needed = adjustment_count != 0;
  // If l_any_adjustment_needed = .true. then wpxp must be modified somewhere, and requires
  // a slow vertically dependent loop to compute properly

  if (l_any_adjustment_needed) {
    // Compute wpxp with a vertically dependent calculation
    const int loop_9_extent_0 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "monotonic_turbulent_flux_limit::loop9", loop_9_extent_0, KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_9_extent_0)) * (1);
          remaining /= (loop_9_extent_0);
          int k;
          int k_zt;
          l_xm_adjustment_needed((i)-1) = false;

          for (k = gr_k_lb_zm + gr_grid_dir_indx;
               (gr_grid_dir_indx) > 0 ? k <= gr_k_ub_zm - gr_grid_dir_indx
                                      : k >= gr_k_ub_zm - gr_grid_dir_indx;
               k += gr_grid_dir_indx) {
            // Calculate index needed for zt level variables based on grid direction.

            if (gr_grid_dir_indx > 0) {
              // Ascending grid
              k_zt = k - 1;
            } else {
              // Descending grid
              k_zt = k;
            }
            // Find the upper limit for w'x' for a monotonic turbulent flux.
            // The following "if" statement ensures there are no "spikes" at the top of the column,
            // which can cause unphysical rtm and thlm tendencies over the height of the column.
            // The fix essentially turns off the monotonic flux limiter for these special cases,
            // but tests show that it still performs well otherwise and runs stably.

            if (l_mono_flux_lim_spikefix && solve_type == mono_flux_rtm &&
                Kokkos::abs(wpxp((i)-1, (k - gr_grid_dir_indx) - 1)) >
                    wpxp_thresh_term((i)-1, (k - gr_grid_dir_indx) - 1) &&
                wpxp((i)-1, (k - gr_grid_dir_indx) - 1) < real_t(0.0)) {
              wpxp_mfl_max((i)-1, (k)-1) = zero;
            } else {
              wpxp_mfl_max((i)-1, (k)-1) =
                  invrs_rho_ds_zm((i)-1, (k)-1) * (wpxp_mfl_max_term_zt((i)-1, (k_zt)-1) +
                                                   rho_ds_zm((i)-1, (k - gr_grid_dir_indx) - 1) *
                                                       wpxp((i)-1, (k - gr_grid_dir_indx) - 1));
            }

            if (wpxp((i)-1, (k)-1) > wpxp_mfl_max((i)-1, (k)-1)) {
              // Determine the net amount of adjustment needed for w'x'.
              wpxp_net_adjust((i)-1, (k)-1) = wpxp_mfl_max((i)-1, (k)-1) - wpxp((i)-1, (k)-1);
              // Reset the value of w'x' to the upper limit allowed by the
              // monotonic flux limiter.
              wpxp((i)-1, (k)-1) = wpxp_mfl_max((i)-1, (k)-1);
              // This value is unneeded for calculation, set to unused_var for stats call
              wpxp_mfl_min((i)-1, (k)-1) = unused_var;
              // If the amount to adjust wpxp by is greater than eps, we will adjust the xm
              // values in the ith column as well

              if (Kokkos::abs(wpxp_net_adjust((i)-1, (k)-1)) > eps) {
                l_xm_adjustment_needed((i)-1) = true;
              }
            } else {
              // Find the lower limit for w'x' for a monotonic turbulent flux.
              wpxp_mfl_min((i)-1, (k)-1) =
                  invrs_rho_ds_zm((i)-1, (k)-1) * (wpxp_mfl_min_term_zt((i)-1, (k_zt)-1) +
                                                   rho_ds_zm((i)-1, (k - gr_grid_dir_indx) - 1) *
                                                       wpxp((i)-1, (k - gr_grid_dir_indx) - 1));

              if (wpxp((i)-1, (k)-1) < wpxp_mfl_min((i)-1, (k)-1)) {
                // Determine the net amount of adjustment needed for w'x'.
                wpxp_net_adjust((i)-1, (k)-1) = wpxp_mfl_min((i)-1, (k)-1) - wpxp((i)-1, (k)-1);
                // Reset the value of w'x' to the lower limit allowed by the
                // monotonic flux limiter.
                wpxp((i)-1, (k)-1) = wpxp_mfl_min((i)-1, (k)-1);
                // If the amount to adjust wpxp by is greater than eps, we will adjust the xm
                // values in the ith column as well

                if (Kokkos::abs(wpxp_net_adjust((i)-1, (k)-1)) > eps) {
                  l_xm_adjustment_needed((i)-1) = true;
                }
              }
            }
          }
        });
    // DO NOT REMOVE COMMENTS - they may be uncommented for testing purposes
    // MONOFLUX TEST !$acc compare( wpxp )
    // MONOFLUX TEST if( any(wpxp_net_adjust /= 0.0) ) write(fstderr,*) "MONOFLUX: wpxp adjusted"
    // Reset the value of xm to compensate for the change to w'x'.

    if (l_mfl_xm_imp_adj) {
      // A tridiagonal matrix is used to semi-implicitly re-solve for the
      // values of xm at timestep index (t+1).
      // Set up the left-hand side of the tridiagonal matrix equation.
      mfl_xm_lhs(nzm, nzt, ngrdcol, dt, gr_weights_zt2zm,    // In
                 gr_invrs_dzt, gr_invrs_dzm, wm_zt,          // In
                 l_implemented, l_upwind_xm_ma, gr_grid_dir, // In
                 lhs_mfl_xm);                                // Out
      // Set up the right-hand side of tridiagonal matrix equation.
      mfl_xm_rhs(nzm, nzt, ngrdcol, dt, xm_old, wpxp, xm_forcing, // In
                 gr_invrs_dzt, rho_ds_zm, invrs_rho_ds_zt,        // In
                 rhs_mfl_xm);                                     // Out
      // Solve the tridiagonal matrix equation.
      mfl_xm_solve(nzt, ngrdcol, gr, solve_type, tridiag_solve_method, // In
                   l_implemented,                                      // In
                   lhs_mfl_xm, rhs_mfl_xm, err_info,                   // InOut
                   xm_mfl);                                            // Out
      // If an adjustment is for a column
      const int loop_10_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
      const int loop_10_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::loop10", loop_10_extent_0 * loop_10_extent_1,
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i = 1 + (remaining % (loop_10_extent_1)) * (1);
            remaining /= (loop_10_extent_1);
            const int k = 1 + (remaining % (loop_10_extent_0)) * (1);
            remaining /= (loop_10_extent_0);

            if (l_xm_adjustment_needed((i)-1)) {
              xm((i)-1, (k)-1) = xm_mfl((i)-1, (k)-1);
            }
          });
      // Check for errors

      if (f2k::clubb_at_least_debug_level_api(0)) {

        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << "\n";
          std::cerr << "Error in monotonic_turbulent_flux_limit after calling mfl_xm_solve" << "\n";
          return;
        }
      }
    } else {
      // An explicit adjustment is made to the values of xm at timestep
      // index (t+1), which is based upon the array of the amounts of w'x'
      // adjustments.
      const int loop_11_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
      const int loop_11_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::loop11", loop_11_extent_0 * loop_11_extent_1,
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i = 1 + (remaining % (loop_11_extent_1)) * (1);
            remaining /= (loop_11_extent_1);
            const int k = 1 + (remaining % (loop_11_extent_0)) * (1);
            remaining /= (loop_11_extent_0);
            real_t dxm_dt_mfl_adjust;

            if (l_xm_adjustment_needed((i)-1)) {
              // The rate of change of the adjustment to xm due to the monotonic
              // flux limiter.
              dxm_dt_mfl_adjust =
                  -invrs_rho_ds_zt((i)-1, (k)-1) * gr_invrs_dzt((i)-1, (k)-1) *
                  (rho_ds_zm((i)-1, (k + 1) - 1) * wpxp_net_adjust((i)-1, (k + 1) - 1) -
                   rho_ds_zm((i)-1, (k)-1) * wpxp_net_adjust((i)-1, (k)-1));
              // The net change to xm due to the monotonic flux limiter is the
              // rate of change multiplied by the time step length.  Add the
              // product to xm to find the new xm resulting from the monotonic
              // flux limiter.
              xm((i)-1, (k)-1) = xm((i)-1, (k)-1) + dxm_dt_mfl_adjust * dt;
            }
          });
    }
    //Ensure there are no spikes at the top of the domain

    // Ensure there are no spikes at the top of the domain. Each integral and
    // subsequent whole-column scaling retains Native sequential level order.
    Kokkos::View<int *, Kokkos::LayoutLeft> spike_warning(
        Kokkos::view_alloc(Kokkos::WithoutInitializing,
                           "monotonic_turbulent_flux_limit::spike_warning"),
        ngrdcol);
    Kokkos::parallel_for(
        "monotonic_turbulent_flux_limit::spike", ngrdcol, KOKKOS_LAMBDA(int column) {
          const int i = column + 1;
          spike_warning(column) = 0;
          if (Kokkos::abs(xm(i - 1, gr_k_ub_zt - 1) - xm_enter_mfl(i - 1, gr_k_ub_zt - 1)) >
              real_t(10.) * xm_tol) {
            const real_t dz =
                gr_zm(i - 1, gr_k_ub_zm - 1) - gr_zm(i - 1, gr_k_ub_zm - gr_grid_dir_indx - 1);
            const real_t xm_density_weighted =
                rho_ds_zt(i - 1, gr_k_ub_zt - 1) *
                (xm(i - 1, gr_k_ub_zt - 1) - xm_enter_mfl(i - 1, gr_k_ub_zt - 1)) * dz;
            real_t xm_vert_integral = real_t(0.0);
            for (int k = gr_k_lb_zt; gr_grid_dir_indx * (k - (gr_k_ub_zt - gr_grid_dir_indx)) <= 0;
                 k += gr_grid_dir_indx)
              xm_vert_integral = xm_vert_integral + rho_ds_zt(i - 1, k - 1) * xm(i - 1, k - 1) *
                                                        gr_grid_dir * gr_dzt(i - 1, k - 1);
            // Check the integral before dividing, without a new numerical denominator guard.
            if (Kokkos::abs(xm_vert_integral) < eps) {
              spike_warning(column) = 1;
              xm(i - 1, gr_k_ub_zt - 1) = xm_enter_mfl(i - 1, gr_k_ub_zt - 1);
            } else {
              real_t xm_adj_coef = xm_density_weighted / xm_vert_integral;
              // The source bound is -0.99, rather than the theoretical -1.
              if (xm_adj_coef < -real_t(0.99)) {
                spike_warning(column) = 2;
                xm_adj_coef = -real_t(0.99);
              }
              for (int k = 1; k <= nzt; ++k)
                xm(i - 1, k - 1) = xm(i - 1, k - 1) * (real_t(1.) + xm_adj_coef);
              xm(i - 1, gr_k_ub_zt - 1) = xm_enter_mfl(i - 1, gr_k_ub_zt - 1);
            }
          }
        });
    // Native CPU warnings are preserved at a host boundary; device kernels never print.
    auto spike_warning_host = to_host_space(spike_warning);
    for (int i = 1; i <= ngrdcol; ++i) {
      if (spike_warning_host(i - 1) == 1)
        std::cerr << "Vertical integral of xm is zero;"
                     "mfl will remove spike at top of domain,"
                     "but it will not conserve xm.\n";
      if (spike_warning_host(i - 1) == 2)
        std::cerr << "xm_adj_coef in mfl less than -0.99, mx_adj_coef set to -0.99\n";
    }

    // DO NOT REMOVE COMMENTS - they may be uncommented for testing purposes
    // MONOFLUX TEST !$acc compare( xm )
    // MONOFLUX TEST write(fstderr,*) "MONOFLUX: xm adjusted"
  }

  if (stats.flags().l_sample) {
    // Boundary conditions
    const int loop_12_extent_0 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "monotonic_turbulent_flux_limit::loop12", loop_12_extent_0, KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_12_extent_0)) * (1);
          remaining /= (loop_12_extent_0);

          wpxp_mfl_min((i)-1, (gr_k_lb_zm)-1) = real_t(0.);
          wpxp_mfl_max((i)-1, (gr_k_lb_zm)-1) = real_t(0.);
          wpxp_mfl_min((i)-1, (gr_k_ub_zm)-1) = real_t(0.);
          wpxp_mfl_max((i)-1, (gr_k_ub_zm)-1) = real_t(0.);
        });
    f2k::stats_update(name_xm_without_ta, xm_without_ta, stats);
    f2k::stats_update(name_xm_mfl_min, min_x_allowable, stats);
    f2k::stats_update(name_xm_mfl_max, max_x_allowable, stats);
    f2k::stats_update(name_wpxp_mfl_min, wpxp_mfl_min, stats);
    f2k::stats_update(name_wpxp_mfl_max, wpxp_mfl_max, stats);
    {
      Real2D stats_values_13(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                "monotonic_turbulent_flux_limit::stats_values_13"),
                             wpxp.extent(0), wpxp.extent(1));
      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::stats_values", wpxp.size(), KOKKOS_LAMBDA(int cell) {
            const int i = cell % wpxp.extent(0), k = cell / wpxp.extent(0);
            stats_values_13(i, k) = wpxp(i, k) / dt;
          });
      f2k::stats_finalize_budget(name_wpxp_mfl, stats_values_13, stats);
    }
    {
      Real2D stats_values_14(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                "monotonic_turbulent_flux_limit::stats_values_14"),
                             xm.extent(0), xm.extent(1));
      Kokkos::parallel_for(
          "monotonic_turbulent_flux_limit::stats_values", xm.size(), KOKKOS_LAMBDA(int cell) {
            const int i = cell % xm.extent(0), k = cell / xm.extent(0);
            stats_values_14(i, k) = xm(i, k) / dt;
          });
      f2k::stats_finalize_budget(name_xm_mfl, stats_values_14, stats);
    }
    f2k::stats_update(name_xm_exit_mfl, xm, stats);
    f2k::stats_update(name_wpxp_exit_mfl, wpxp, stats);
  }
  return;
}

namespace {

//=============================================================================
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// dt: Model timestep length                      [s]
// wm_zt: w wind component on thermodynamic levels   [m/s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// grid_dir: Grid direction multiplier (+1 = ascending; -1 = descending)
// lhs: Left hand side of tridiagonal matrix
void mfl_xm_lhs(int nzm, int nzt, int ngrdcol, real_t dt, Real3D weights_zt2zm,         // In
                Real2D invrs_dzt, Real2D invrs_dzm,                                     // In
                Real2D wm_zt, bool l_implemented, bool l_upwind_xm_ma, real_t grid_dir, // In
                Real3D lhs)                                                             // Out
{
  // Description:
  // This subroutine is part of the process of re-solving for xm at timestep
  // index (t+1).  This is done because the original solving process produced
  // values outside of what is deemed acceptable by the monotonic flux limiter.
  // Unlike the original formulation for advancing xm one timestep, which
  // combines w'x' and xm in a band-diagonal solver, this formulation uses a
  // tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
  // is known.
  //
  // Subroutine mfl_xm_lhs sets up the left-hand side of the matrix equation.
  // Constant parameters
  //---------------------------- Input Variables ----------------------------
  //---------------------------- Output Variables ----------------------------
  //---------------------------- Local Variables ----------------------------
  constexpr int t_above = 1;
  constexpr int t_below = 2;
  constexpr int k_tdiag = 2;
  int i;
  int k;
  int b; // Array index
  //---------------------------- Begin Code ----------------------------
  // Setup LHS of the tridiagonal system
  // LHS xm mean advection (ma) term.

  if (!l_implemented) {
    term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, weights_zt2zm, // In
                   invrs_dzt, invrs_dzm,                    // In
                   l_upwind_xm_ma, grid_dir,                // In
                   lhs);                                    // Out
  } else {
    const int loop_1_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
    const int loop_1_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);
    const int loop_1_extent_2 = std::max(0, (ndiags3 - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "mfl_xm_lhs::loop1", loop_1_extent_0 * loop_1_extent_1 * loop_1_extent_2,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_1_extent_1)) * (1);
          remaining /= (loop_1_extent_1);
          const int k = 1 + (remaining % (loop_1_extent_0)) * (1);
          remaining /= (loop_1_extent_0);
          const int b = 1 + (remaining % (loop_1_extent_2)) * (1);
          remaining /= (loop_1_extent_2);

          lhs((b)-1, (i)-1, (k)-1) = real_t(0.0);
        });
  }
  const int loop_2_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_2_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "mfl_xm_lhs::loop2", loop_2_extent_0 * loop_2_extent_1, KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_2_extent_1)) * (1);
        remaining /= (loop_2_extent_1);
        const int k = 1 + (remaining % (loop_2_extent_0)) * (1);
        remaining /= (loop_2_extent_0);

        // LHS xm time tendency.
        lhs((k_tdiag)-1, (i)-1, (k)-1) = lhs((k_tdiag)-1, (i)-1, (k)-1) + real_t(1.0) / dt;
      });
  return;
}

//=============================================================================
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// dt: Model timestep length                    [s]
// xm_old: xm; timestep (t) (thermodynamic levels)  [units vary]
// wpxp: w'x'; timestep (t+1); limited (m-levs.)  [units vary]
// xm_forcing: xm forcings (thermodynamic levels)       [units vary]
// invrs_dzt: The inverse spacing between momentum grid levels; centered over thermodynamic grid
// levels.
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// rhs: Right hand side of tridiagonal matrix equation
void mfl_xm_rhs(int nzm, int nzt, int ngrdcol, real_t dt, Real2D xm_old, Real2D wpxp, // In
                Real2D xm_forcing,                                                    // In
                Real2D invrs_dzt, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,           // In
                Real2D rhs)                                                           // Out
{
  // Description:
  // This subroutine is part of the process of re-solving for xm at timestep
  // index (t+1).  This is done because the original solving process produced
  // values outside of what is deemed acceptable by the monotonic flux limiter.
  // Unlike the original formulation for advancing xm one timestep, which
  // combines w'x' and xm in a band-diagonal solver, this formulation uses a
  // tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
  // is known.
  //
  // Subroutine mfl_xm_rhs sets up the right-hand side of the matrix equation.
  //---------------------------- Input Variables ----------------------------
  //---------------------------- Output Variable ----------------------------
  //---------------------------- Local Variables ----------------------------
  int i;
  int k; // Array indices
  //---------------------------- Begin Code ----------------------------
  const int loop_1_extent_0 = std::max(0, (nzt - (1)) / (1) + 1);
  const int loop_1_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "mfl_xm_rhs::loop1", loop_1_extent_0 * loop_1_extent_1, KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_1_extent_1)) * (1);
        remaining /= (loop_1_extent_1);
        const int k = 1 + (remaining % (loop_1_extent_0)) * (1);
        remaining /= (loop_1_extent_0);

        // RHS xm time tendency.
        rhs((i)-1, (k)-1) = xm_old((i)-1, (k)-1) / dt;
        // RHS xm turbulent advection (ta) term.
        // Note:  Normally, the turbulent advection (ta) term is treated
        //        implicitly when advancing xm one timestep, as both xm and w'x'
        //        are advanced together from timestep index (t) to timestep
        //        index (t+1).  However, in this case, both xm and w'x' have
        //        already been advanced one timestep.  However, w'x'(t+1) has been
        //        limited after the fact, and therefore it's values at timestep
        //        index (t+1) are known.  Thus, in re-solving for xm(t+1), the
        //        derivative of w'x'(t+1) can be placed on the right-hand side of
        //        the d(xm)/dt equation.
        rhs((i)-1, (k)-1) =
            rhs((i)-1, (k)-1) - invrs_rho_ds_zt((i)-1, (k)-1) * invrs_dzt((i)-1, (k)-1) *
                                    (rho_ds_zm((i)-1, (k + 1) - 1) * wpxp((i)-1, (k + 1) - 1) -
                                     rho_ds_zm((i)-1, (k)-1) * wpxp((i)-1, (k)-1));
        // RHS xm forcings.
        // Note: xm forcings include the effects of microphysics,
        //       cloud water sedimentation, radiation, and any
        //       imposed forcings on xm.
        rhs((i)-1, (k)-1) = rhs((i)-1, (k)-1) + xm_forcing((i)-1, (k)-1);
      });
  return;
}

//=============================================================================
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// solve_type: Variables being solved for.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs: Left hand side of tridiagonal matrix
// rhs: Right hand side of tridiagonal matrix equation
// err_info: err_info struct containing err_code and err_header
// xm: Value of variable being solved for at timestep (t+1)   [units vary]
void mfl_xm_solve(int nzt, int ngrdcol, const grid &gr, int solve_type, // In
                  int tridiag_solve_method,                             // In
                  bool l_implemented,                                   // In
                  Real3D lhs, Real2D rhs, err_info_type &err_info,      // InOut
                  Real2D xm)                                            // Out
{
  // Description:
  // This subroutine is part of the process of re-solving for xm at timestep
  // index (t+1).  This is done because the original solving process produced
  // values outside of what is deemed acceptable by the monotonic flux limiter.
  // Unlike the original formulation for advancing xm one timestep, which
  // combines w'x' and xm in a band-diagonal solver, this formulation uses a
  // tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
  // is known.
  //
  // Subroutine mfl_xm_solve solves the tridiagonal matrix equation for xm at
  // timestep index (t+1).
  //---------------------------- Input Variables ----------------------------
  //---------------------------- InOut Variables ----------------------------
  //---------------------------- Output Variables ----------------------------
  //---------------------------- Local Variable ----------------------------
  std::string solve_type_str; // solve_type as a string for debug output purposes
  const auto gr_grid_dir_indx = gr.grid_dir_indx;

  switch (solve_type) {
  case mono_flux_rtm:
    solve_type_str = "rtm";
    break;
  case mono_flux_thlm:
    solve_type_str = "thlm";
    break;
  default:
    solve_type_str = "scalars";
    break;
  }
  // Matrix solves differ by grid direction. Native forces both vertical rows
  // and bands to descending order when requested, before solving.
  const bool reverse = l_force_descending_solves && gr.grid_dir_indx > 0;
  Real3D solver_lhs(Kokkos::view_alloc(Kokkos::WithoutInitializing, "mfl_xm_solve::solver_lhs"),
                    ngrdcol, ndiags3, nzt);
  Real2D solver_rhs(Kokkos::view_alloc(Kokkos::WithoutInitializing, "mfl_xm_solve::solver_rhs"),
                    ngrdcol, nzt);
  Real2D solver_xm(Kokkos::view_alloc(Kokkos::WithoutInitializing, "mfl_xm_solve::solver_xm"),
                   ngrdcol, nzt);
  // This is the actual matrix API boundary: Native bands-first -> column-first.
  Kokkos::parallel_for(
      "mfl_xm_solve::pack", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        const int row = reverse ? nzt - 1 - k : k;
        solver_rhs(i, k) = rhs(i, row);
        for (int b = 0; b < ndiags3; ++b)
          solver_lhs(i, b, k) = lhs(reverse ? ndiags3 - 1 - b : b, i, row);
      });
  tridiag_solve(solve_type_str, tridiag_solve_method, ngrdcol, nzt, l_implemented, // In
                solver_lhs, solver_rhs, err_info,                                  // InOut
                solver_xm);                                                        // Out
  // Restore Native InOut arrays in their post-flip/post-solver order. They do
  // not flip back. Restore the still-unflipped solution before the fatal return.
  Kokkos::parallel_for(
      "mfl_xm_solve::restore", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        rhs(i, k) = solver_rhs(i, k);
        xm(i, k) = solver_xm(i, k);
        for (int b = 0; b < ndiags3; ++b)
          lhs(b, i, k) = solver_lhs(i, b, k);
      });
  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n"
                << "calling tridiag_solve in mfl_xm_solve\n";
      return;
    }
  }
  // Only a successful Native return reaches the final solution direction reversal.
  if (reverse) {
    Kokkos::parallel_for(
        "mfl_xm_solve::flip_solution", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          xm(i, k) = solver_xm(i, nzt - 1 - k);
        });
  }
  return;
}
} // namespace

//=============================================================================
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// dt: Model timestep length                       [s]
// w_1_zm: Mean w (1st PDF component)                   [m/s]
// w_2_zm: Mean w (2nd PDF component)                   [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)            [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)            [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// stats: Native statistics context, mutated by sampled budget calls [-]
// low_lev_effect: Index of lowest level that has an effect (for lev. k)
// high_lev_effect: Index of highest level that has an effect (for lev. k)
void calc_turb_adv_range(int nzm, int nzt, int ngrdcol, const grid &gr, real_t dt, // In
                         Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm,
                         Real2D varnce_w_2_zm,                                     // In
                         Real2D mixt_frac_zm,                                      // In
                         f2k::stats_type &stats,                                   // InOut
                         Kokkos::View<int **, Kokkos::LayoutLeft> low_lev_effect,  // Out
                         Kokkos::View<int **, Kokkos::LayoutLeft> high_lev_effect) // Out
{
  // Description:
  // Calculates the lowermost and uppermost thermodynamic grid levels that can
  // effect the base (or central) thermodynamic level through the effects of
  // turbulent advection over the course of one time step.  This is used as
  // part of the monotonic turbulent advection scheme.
  //
  // One method is to use the vertical velocity at each level to determine the
  // amount of time that it takes to travel across that particular grid level.
  // The method is to keep on advancing one grid level until either (a) the
  // total sum of time taken reaches or exceeds the model time step length,
  // (b) the top or bottom of the model is reached, or (c) a level is reached
  // where the vertical velocity component (with turbulence included) is
  // oriented completely opposite of the direction of travel towards the base
  // (or central) thermodynamic level.  An example of situation (c) would be,
  // while starting from a higher altitude and searching downward for all
  // upward vertical velocity components, encountering a strong downdraft
  // where the vertical velocity at every single point is oriented downward.
  // Such a situation would occur when the mean vertical velocity (wm_zm)
  // exceeds any turbulent component (w') that would be oriented upwards.
  //
  // Another method is to simply set the thickness (in meters) of the layer
  // that turbulent advection is allowed to act over, for purposes of the
  // monotonic turbulent advection scheme.  The lowermost and uppermost
  // grid level that can effect the base (or central) thermodynamic level
  // is computed based on the thickness and altitude of each level.
  //
  // References:
  //-----------------------------------------------------------------------
  // Constant parameters
  //------------------------- Input Variables -------------------------
  //------------------------- Output Variables -------------------------
  //------------------------- Local Variables -------------------------
  constexpr bool l_constant_thickness = false;
  constexpr real_t const_thick = real_t(150.0);
  Real2D vert_vel_up_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_turb_adv_range::vert_vel_up_zm"),
      ngrdcol, nzm); // Average upwards vertical velocity component   [m/s]
  Real2D vert_vel_down_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_turb_adv_range::vert_vel_down_zm"),
      ngrdcol, nzm); // Average downwards vertical velocity component [m/s]
  Real2D w_min(Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_turb_adv_range::w_min"),
               ngrdcol, nzm); // Minimum velocity to affect adjacent levels    [m/s]
  real_t dt_one_grid_lev;     // Amount of time to travel one grid box           [s]
  real_t dt_all_grid_levs;    // Running count of amount of time taken to travel [s]
  real_t invrs_dt;            // Inverse of timestep, used to reduce divides     [1/s]
  int k;
  int i;
  int j;
  int j_adj;
  const auto gr_dzm = gr.dzm;
  const auto gr_grid_dir = gr.grid_dir;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_zt = gr.zt;
  //------------------------- Begin Code -------------------------

  if (l_constant_thickness) {
    // The value of w'x' may only be altered between (momentum)
    // levels 3 and gr%nzm-2.
    const int loop_1_extent_0 =
        std::max(0, (gr_k_ub_zt - 2 * gr_grid_dir_indx - (gr_k_lb_zt + gr_grid_dir_indx)) /
                            (gr_grid_dir_indx) +
                        1);
    const int loop_1_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "calc_turb_adv_range::loop1", loop_1_extent_0 * loop_1_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_1_extent_1)) * (1);
          remaining /= (loop_1_extent_1);
          const int k =
              gr_k_lb_zt + gr_grid_dir_indx + (remaining % (loop_1_extent_0)) * (gr_grid_dir_indx);
          remaining /= (loop_1_extent_0);
          int j;
          // Compute the number of levels that effect the central thermodynamic
          // level through upwards motion (traveling from lower thermodynamic
          // levels to reach the central thermodynamic level).
          // Start with the index of the thermodynamic level immediately below
          // the central thermodynamic level.
          j = k - gr_grid_dir_indx;

          while (true) {

            if (gr_zt((i)-1, (k)-1) - gr_zt((i)-1, (j)-1) >= const_thick) {
              // Stop, the current grid level is the lowest level that can
              // be considered.
              low_lev_effect((i)-1, (k)-1) = j;
              break;
            } else {
              // The model lower boundary has been reached.

              if (j == gr_k_lb_zt) {
                // The current thermodynamic level (lower boundary) is the
                // lowest level that can be considered.
                low_lev_effect((i)-1, (k)-1) = j;
                break;
              } else {
                // Increment to the next vertical level down.
                j = j - gr_grid_dir_indx;
              }
            }
          }
        });
    // Compute the number of levels that effect the central thermodynamic
    // level through downwards motion (traveling from higher thermodynamic
    // levels to reach the central thermodynamic level).
    const int loop_2_extent_0 =
        std::max(0, (gr_k_ub_zt - 2 * gr_grid_dir_indx - (gr_k_lb_zt + gr_grid_dir_indx)) /
                            (gr_grid_dir_indx) +
                        1);
    const int loop_2_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "calc_turb_adv_range::loop2", loop_2_extent_0 * loop_2_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_2_extent_1)) * (1);
          remaining /= (loop_2_extent_1);
          const int k =
              gr_k_lb_zt + gr_grid_dir_indx + (remaining % (loop_2_extent_0)) * (gr_grid_dir_indx);
          remaining /= (loop_2_extent_0);
          int j;
          // Start with the index of the thermodynamic level immediately above
          // the central thermodynamic level.
          j = k + gr_grid_dir_indx;

          while (true) {

            if (gr_zt((i)-1, (j)-1) - gr_zt((i)-1, (k)-1) >= const_thick) {
              // Stop, the current grid level is the highest level that can
              // be considered.
              high_lev_effect((i)-1, (k)-1) = j;
              break;
            } else {
              // The highest level that can be considered is thermodynamic
              // level gr%nzt.

              if (j == gr_k_ub_zt) {
                // The current level (upper boundary) is the highest level
                // that can be considered.
                high_lev_effect((i)-1, (k)-1) = j;
                break;
              } else {
                // Increment to the next vertical level up.
                j = j + gr_grid_dir_indx;
              }
            }
          }
        });
  } else {
    invrs_dt = real_t(1.0) / dt;
    const int loop_3_extent_0 = std::max(0, (nzm - (1)) / (1) + 1);
    const int loop_3_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "calc_turb_adv_range::loop3", loop_3_extent_0 * loop_3_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_3_extent_1)) * (1);
          remaining /= (loop_3_extent_1);
          const int k = 1 + (remaining % (loop_3_extent_0)) * (1);
          remaining /= (loop_3_extent_0);

          w_min((i)-1, (k)-1) = gr_grid_dir * gr_dzm((i)-1, (k)-1) * invrs_dt;
        });
    // Find the average upwards vertical velocity and the average downwards
    // vertical velocity.
    // Note:  A level that has all vertical wind moving downwards will have a
    //        vert_vel_up value that is 0, and vice versa.
    mean_vert_vel_up_down(nzm, ngrdcol,                                 // In
                          w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, // In
                          mixt_frac_zm, real_t(0.0), w_min,             // In
                          stats,                                        // InOut
                          vert_vel_down_zm, vert_vel_up_zm);            // Out
    const int loop_4_extent_0 =
        std::max(0, (gr_k_ub_zt - 2 * gr_grid_dir_indx - (gr_k_lb_zt + gr_grid_dir_indx)) /
                            (gr_grid_dir_indx) +
                        1);
    const int loop_4_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "calc_turb_adv_range::loop4", loop_4_extent_0 * loop_4_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_4_extent_1)) * (1);
          remaining /= (loop_4_extent_1);
          const int k =
              gr_k_lb_zt + gr_grid_dir_indx + (remaining % (loop_4_extent_0)) * (gr_grid_dir_indx);
          remaining /= (loop_4_extent_0);
          real_t dt_one_grid_lev;
          real_t dt_all_grid_levs;
          int j;
          int j_adj;
          // Compute the number of levels that effect the central thermodynamic
          // level through upwards motion (traveling from lower thermodynamic
          // levels to reach the central thermodynamic level).
          // Initialize the overall delta t counter to 0.
          dt_all_grid_levs = real_t(0.0);
          // Start with the index of the thermodynamic level immediately below
          // the central thermodynamic level.

          for (j = k - gr_grid_dir_indx;
               (-gr_grid_dir_indx) > 0 ? j <= gr_k_lb_zt : j >= gr_k_lb_zt;
               j += -gr_grid_dir_indx) {
            low_lev_effect((i)-1, (k)-1) = j;
            // Accounting for the offset between thermodynamic and momentum
            // level indices depending on an ascending or descending grid
            // direction. This is for use with momentum-level variables in
            // the following calculation.

            if (gr_grid_dir_indx > 0) {
              // Ascending grid
              j_adj = j + 1;
            } else {
              // Descending grid
              j_adj = j;
            }
            // Continue if there is some component of upwards vertical velocity.

            if (vert_vel_up_zm((i)-1, (j_adj)-1) > real_t(0.0)) {
              // Compute the amount of time it takes to travel one grid level
              // upwards:  delta_t = delta_z / vert_vel_up_zm.
              dt_one_grid_lev =
                  gr_grid_dir * gr_dzm((i)-1, (j_adj)-1) / vert_vel_up_zm((i)-1, (j_adj)-1);
              // Total time elapsed for crossing all grid levels that have been
              // passed, thus far.
              dt_all_grid_levs = dt_all_grid_levs + dt_one_grid_lev;
              // Stop if has taken more than one model time step (overall) to
              // travel the entire extent of the current vertical grid level.

              if (dt_all_grid_levs >= dt) {
                // The current level is the lowest level that can be
                // considered.
                break;
              }
              // Stop if there isn't a component of upwards vertical velocity.
            } else {
              // The current level cannot be considered.  The lowest level that
              // can be considered is one-level-above the current level.
              low_lev_effect((i)-1, (k)-1) = j + gr_grid_dir_indx;
              break;
            }
          }
        });
    // Compute the number of levels that effect the central thermodynamic
    // level through downwards motion (traveling from higher thermodynamic
    // levels to reach the central thermodynamic level).
    const int loop_5_extent_0 =
        std::max(0, (gr_k_ub_zt - 2 * gr_grid_dir_indx - (gr_k_lb_zt + gr_grid_dir_indx)) /
                            (gr_grid_dir_indx) +
                        1);
    const int loop_5_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

    Kokkos::parallel_for(
        "calc_turb_adv_range::loop5", loop_5_extent_0 * loop_5_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i = 1 + (remaining % (loop_5_extent_1)) * (1);
          remaining /= (loop_5_extent_1);
          const int k =
              gr_k_lb_zt + gr_grid_dir_indx + (remaining % (loop_5_extent_0)) * (gr_grid_dir_indx);
          remaining /= (loop_5_extent_0);
          real_t dt_one_grid_lev;
          real_t dt_all_grid_levs;
          int j;
          int j_adj;
          // Initialize the overall delta t counter to 0.
          dt_all_grid_levs = real_t(0.0);
          // Start with the index of the thermodynamic level immediately above
          // the central thermodynamic level.

          for (j = k + gr_grid_dir_indx; (gr_grid_dir_indx) > 0 ? j <= gr_k_ub_zt : j >= gr_k_ub_zt;
               j += gr_grid_dir_indx) {
            high_lev_effect((i)-1, (k)-1) = j;
            // Accounting for the offset between thermodynamic and momentum
            // level indices depending on an ascending or descending grid
            // direction. This is for use with momentum-level variables in
            // the following calculation.

            if (gr_grid_dir_indx > 0) {
              // Ascending grid
              j_adj = j;
            } else {
              // Descending grid
              j_adj = j + 1;
            }
            // Continue if there is some component of downwards vertical velocity.

            if (vert_vel_down_zm((i)-1, (j_adj)-1) < real_t(0.0)) {
              // Compute the amount of time it takes to travel one grid level
              // downwards:  delta_t = - delta_z / vert_vel_down_zm.
              // Note:  There is a (-) sign in front of delta_z because the
              //        distance traveled is downwards.  Since vert_vel_down
              //        has a negative value, dt_one_grid_lev will be a
              //        positive value.
              dt_one_grid_lev =
                  -gr_grid_dir * gr_dzm((i)-1, (j_adj)-1) / vert_vel_down_zm((i)-1, (j_adj)-1);
              // Total time elapsed for crossing all grid levels that have been
              // passed, thus far.
              dt_all_grid_levs = dt_all_grid_levs + dt_one_grid_lev;
              // Stop if has taken more than one model time step (overall) to
              // travel the entire extent of the current vertical grid level.

              if (dt_all_grid_levs >= dt) {
                // The current level is the highest level that can be
                // considered.
                break;
              }
              // Stop if there isn't a component of downwards vertical velocity.
            } else {
              // The current level cannot be considered.  The highest level
              // that can be considered is one-level-below the current level.
              high_lev_effect((i)-1, (k)-1) = j - gr_grid_dir_indx;
              break;
            }
          }
        });
  }
  // Information for levels gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, and
  // gr%k_ub_zt is not needed. However, set the values at these levels for
  // purposes of not having odd values in the arrays.
  const int loop_6_extent_0 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "calc_turb_adv_range::loop6", loop_6_extent_0, KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_6_extent_0)) * (1);
        remaining /= (loop_6_extent_0);

        low_lev_effect((i)-1, (gr_k_lb_zt)-1) = gr_k_lb_zt;
        high_lev_effect((i)-1, (gr_k_lb_zt)-1) = gr_k_lb_zt;
        low_lev_effect((i)-1, (gr_k_ub_zt - gr_grid_dir_indx) - 1) = gr_k_ub_zt - gr_grid_dir_indx;
        high_lev_effect((i)-1, (gr_k_ub_zt - gr_grid_dir_indx) - 1) = gr_k_ub_zt;
        low_lev_effect((i)-1, (gr_k_ub_zt)-1) = gr_k_ub_zt;
        high_lev_effect((i)-1, (gr_k_ub_zt)-1) = gr_k_ub_zt;
      });
  return;
}

namespace {

//=============================================================================
// nzm: Number of momentum levels [-]
// ngrdcol: Number of columns [-]
// w_1_zm: Mean w (1st PDF component)                   [m/s]
// w_2_zm: Mean w (2nd PDF component)                   [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)            [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)            [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// w_ref: Reference velocity, w|_ref (normally = 0)   [m/s]
// w_min: Minimum velocity to affect adjacent level    [m/s]
// stats: Native statistics context, mutated by sampled budget calls [-]
// mean_w_down: Overall mean w (<= w|_ref)                  [m/s]
// mean_w_up: Overall mean w (>= w|_ref)                  [m/s]
void mean_vert_vel_up_down(int nzm, int ngrdcol, // In
                           Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm,
                           Real2D varnce_w_2_zm,                            // In
                           Real2D mixt_frac_zm, real_t w_ref, Real2D w_min, // In
                           f2k::stats_type &stats,                          // InOut
                           Real2D mean_w_down, Real2D mean_w_up)            // Out
{
  // Description
  // The values of vertical velocity, along a horizontal plane at any given
  // vertical level, are not allowed by CLUBB to be uniform.  In other words,
  // there must be some variance in vertical velocity.  This subroutine
  // calculates the mean of all values of vertical velocity, at any given
  // vertical level, that are greater than a certain reference velocity.  This
  // subroutine also calculates the mean of all values of vertical velocity, at
  // any given vertical level, that are less than a certain reference velocity.
  // The reference velocity is usually 0 m/s, in which case this subroutine
  // calculates the average positive (upward) velocity and the average negative
  // (downward) velocity.  However, the reference velocity may be other values,
  // such as wm_zm, which is the overall mean vertical velocity.  If the
  // reference velocity is wm_zm, this subroutine calculates the average of all
  // values of w that are on the positive ("upward") side of the mean and the
  // average of all values of w that are on the negative ("downward") side of
  // the mean.  These mean positive and negative vertical velocities are useful
  // in determining how long, on average, it takes a parcel of air, being
  // driven by subgrid updrafts or downdrafts, to traverse the length of the
  // vertical grid level.
  //
  // Method
  // ------
  //
  // The CLUBB model uses a joint PDF of vertical velocity, liquid water
  // potential temperature, and total water mixing ratio to determine subgrid
  // variability.
  //
  // The values of vertical velocity, w, along an undefined horizontal plane
  // at any vertical level, are considered to approximately follow a
  // distribution that is a mixture of two normal (or Gaussian) distributions.
  // The values of w that are a part of the 1st normal distribution are
  // referred to as w_1, and the values of w that are part of the 2nd normal
  // distribution are referred to as w_2.  Note that these distributions
  // overlap, and there are many values of w that are found in both w_1 and w_2.
  //
  // The probability density function (PDF) for w, P(w), is:
  //
  // P(w) = mixt_frac*P(w_1) + (1-mixt_frac)*P(w_2);
  //
  // where "mixt_frac" is the weight of the 1st normal distribution, and P(w_1) and
  // P(w_2) are the equations for the 1st and 2nd normal distributions,
  // respectively:
  //
  // P(w_1) = 1 / ( sigma_w_1 * sqrt(2*PI) )
  //         * EXP[ -(w_1-mu_w_1)^2 / (2*sigma_w_1^2) ]; and
  //
  // P(w_2) = 1 / ( sigma_w_2 * sqrt(2*PI) )
  //         * EXP[ -(w_2-mu_w_2)^2 / (2*sigma_w_2^2) ].
  //
  // The mean of the 1st normal distribution is mu_w_1, and the standard
  // deviation of the 1st normal distribution is sigma_w_1.  The mean of the
  // 2nd normal distribution is mu_w_2, and the standard deviation of the 2nd
  // normal distribution is sigma_w_2.
  //
  // The average value of w, distributed according to the probability
  // distribution, between limits alpha and beta, is:
  //
  // <w|_(alpha:beta)> = INT(alpha:beta) w P(w) dw.
  //
  // The average value of w over a certain domain is used to determine the
  // average positive and negative (as compared to the reference velocity)
  // values of w at any vertical level.
  //
  // Average Negative Vertical Velocity
  // ----------------------------------
  //
  // The average of all values of w in the distribution that are below the
  // reference velocity, w|_ref, is the mean value of w over the domain
  // -inf <= w <= w|_ref, such that:
  //
  // <w|_(-inf:w|_ref)> = INT(-inf:w|_ref) w P(w) dw.
  //                    = mixt_frac * INT(-inf:w|_ref) w_1 P(w_1) dw_1
  //                      + (1-mixt_frac) * INT(-inf:w|_ref) w_2 P(w_2) dw_2.
  //
  // For each normal distribution in the mixture of normal distribution, i
  // (where "i" can be 1 or 2):
  //
  // INT(-inf:w|_ref) wi P(wi) dwi =
  //   - ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
  //   + mu_wi * (1/2)*[ 1 + erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];
  //
  // where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
  // the standard deviations of w for the ith normal distribution, and erf( )
  // is the error function.
  //
  // The mean of all values of w <= w|_ref is:
  //
  // <w|_(-inf:w|_ref)> =
  // mixt_frac * { - ( sigma_w_1 / sqrt(2*PI) )
  //                 * EXP[ -(w|_ref-mu_w_1)^2 / (2*sigma_w_1^2) ]
  //               + mu_w_1 * (1/2)
  //                 *[1 + erf( (w|_ref-mu_w_1) / (sqrt(2)*sigma_w_1) )] }
  // + (1-mixt_frac) * { - ( sigma_w_2 / sqrt(2*PI) )
  //                       * EXP[ -(w|_ref-mu_w_2)^2 / (2*sigma_w_2^2) ]
  //                     + mu_w_2 * (1/2)
  //                       *[1 + erf( (w|_ref-mu_w_2) / (sqrt(2)*sigma_w_2) )] }.
  //
  // Average Positive Vertical Velocity
  // ----------------------------------
  //
  // The average of all values of w in the distribution that are above the
  // reference velocity, w|_ref, is the mean value of w over the domain
  // w|_ref <= w <= inf, such that:
  //
  // <w|_(w|_ref:inf)> = INT(w|_ref:inf) w P(w) dw.
  //                   = mixt_frac * INT(w|_ref:inf) w_1 P(w_1) dw_1
  //                     + (1-mixt_frac) * INT(w|_ref:inf) w_2 P(w_2) dw_2.
  //
  // For each normal distribution in the mixture of normal distribution, i
  // (where "i" can be 1 or 2):
  //
  // INT(w|_ref:inf) wi P(wi) dwi =
  //     ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
  //   + mu_wi * (1/2)*[ 1 - erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];
  //
  // where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
  // the standard deviations of w for the ith normal distribution, and erf( )
  // is the error function.
  //
  // The mean of all values of w >= w|_ref is:
  //
  // <w|_(w|_ref:inf)> =
  // mixt_frac * {   ( sigma_w_1 / sqrt(2*PI) )
  //                * EXP[ -(w|_ref-mu_w_1)^2 / (2*sigma_w_1^2) ]
  //               + mu_w_1 * (1/2)
  //                 *[1 - erf( (w|_ref-mu_w_1) / (sqrt(2)*sigma_w_1) )] }
  // + (1-mixt_frac) * {   ( sigma_w_2 / sqrt(2*PI) )
  //                      * EXP[ -(w|_ref-mu_w_2)^2 / (2*sigma_w_2^2) ]
  //                     + mu_w_2 * (1/2)
  //                       *[1 - erf( (w|_ref-mu_w_2) / (sqrt(2)*sigma_w_2) )] }.
  //
  // Special Limitations:
  // --------------------
  //
  // A normal distribution has a domain from -inf to inf.  However, the mixture
  // of normal distributions is an approximation of the distribution of values
  // of w along a horizontal plane at any given vertical level.  Vertical
  // velocity, w, has absolute minimum and maximum values (that cannot be
  // predicted by the PDF).  The absolute maximum and minimum for each normal
  // distribution is most likely found within 2 or 3 standard deviations of the
  // mean for the relevant normal distribution.  In other words, for each
  // normal distribution in the mixture of normal distributions, all the values
  // of w are found within 2 or 3 standard deviations on both sides of the
  // mean.  Therefore, if one (or both) of the normal distributions has a mean
  // that is more than 3 standard deviations away from the reference velocity,
  // then that entire w distribution is found on ONE side of the reference
  // velocity.
  //
  // Therefore:
  //
  // a) where mu_wi + 3*sigma_wi <= w|_ref:
  //
  //       The entire ith normal distribution of w is on the negative side of
  //       w|_ref; and
  //
  //       INT(-inf:w|_ref) wi P(wi) dwi = mu_wi; and
  //       INT(inf:w|_ref) wi P(wi) dwi = 0.
  //
  // b) where mu_wi - 3*sigma_wi >= w|_ref:
  //
  //       The entire ith normal distribution of w is on the positive side of
  //       w|_ref; and
  //
  //       INT(-inf:w|_ref) wi P(wi) dwi = 0; and
  //       INT(inf:w|_ref) wi P(wi) dwi = mu_wi.
  //
  // Notes: A value of 3 standard deviations above and below the mean of the
  //        ith normal distribution was chosen for the approximate maximum and
  //        minimum values of the ith normal distribution because 99.7% of
  //        values in a normal distribution are found within 3 standard
  //        deviations from the mean (compared to 95.4% for 2 standard
  //        deviations).  The value of 3 standard deviations provides for a
  //        reasonable estimate of the absolute maximum and minimum of w, while
  //        covering a great majority of the normal distribution.
  //
  //        In addition to approximating the up and down components of w
  //        by checking if the pdfs are greater than 3 standard deviations
  //        from the mean, there is now a case to approximate when w is
  //        too small in general. The input array, w_min, contains the
  //        minimum values of vertical velocity that would be required
  //        at a given grid level for that grid box to be able to affect
  //        the adjacent levels. If the magnitude of w at a given level
  //        is less than 3 standard deviations below w_min for that level,
  //        then there is no significant portion of the air from that grid
  //        box that is capable of interacting with the next level, and
  //        the upward and downward components for that pdf are set to 0.
  //
  // References:
  //-----------------------------------------------------------------------
  //------------------------- Input Variables -------------------------
  //------------------------- Output Variables -------------------------
  //------------------------- Local Variables -------------------------
  Real2D mean_w_down_1st(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mean_vert_vel_up_down::mean_w_down_1st"),
      ngrdcol, nzm); // Mean w (<= w|_ref) from 1st normal distribution [m/s]
  Real2D mean_w_down_2nd(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mean_vert_vel_up_down::mean_w_down_2nd"),
      ngrdcol, nzm); // Mean w (<= w|_ref) from 2nd normal distribution [m/s]
  Real2D mean_w_up_1st(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mean_vert_vel_up_down::mean_w_up_1st"),
      ngrdcol, nzm); // Mean w (>= w|_ref) from 1st normal distribution [m/s]
  Real2D mean_w_up_2nd(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mean_vert_vel_up_down::mean_w_up_2nd"),
      ngrdcol, nzm); // Mean w (>= w|_ref) from 2nd normal distribution [m/s]
  int i;
  int k;
  //------------------------- Begin Code -------------------------
  calc_mean_w_up_down_component(nzm, ngrdcol,                    // In
                                w_1_zm, varnce_w_1_zm,           // In
                                w_ref, w_min,                    // In
                                mean_w_down_1st, mean_w_up_1st); // Out
  calc_mean_w_up_down_component(nzm, ngrdcol,                    // In
                                w_2_zm, varnce_w_2_zm,           // In
                                w_ref, w_min,                    // In
                                mean_w_down_2nd, mean_w_up_2nd); // Out
  // Overall mean of downwards w.
  const int loop_1_extent_0 = std::max(0, (nzm - (1)) / (1) + 1);
  const int loop_1_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "mean_vert_vel_up_down::loop1", loop_1_extent_0 * loop_1_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_1_extent_1)) * (1);
        remaining /= (loop_1_extent_1);
        const int k = 1 + (remaining % (loop_1_extent_0)) * (1);
        remaining /= (loop_1_extent_0);

        mean_w_down((i)-1, (k)-1) =
            mixt_frac_zm((i)-1, (k)-1) * mean_w_down_1st((i)-1, (k)-1) +
            (real_t(1.0) - mixt_frac_zm((i)-1, (k)-1)) * mean_w_down_2nd((i)-1, (k)-1);
      });
  // Overall mean of upwards w.
  const int loop_2_extent_0 = std::max(0, (nzm - (1)) / (1) + 1);
  const int loop_2_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "mean_vert_vel_up_down::loop2", loop_2_extent_0 * loop_2_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_2_extent_1)) * (1);
        remaining /= (loop_2_extent_1);
        const int k = 1 + (remaining % (loop_2_extent_0)) * (1);
        remaining /= (loop_2_extent_0);

        mean_w_up((i)-1, (k)-1) =
            mixt_frac_zm((i)-1, (k)-1) * mean_w_up_1st((i)-1, (k)-1) +
            (real_t(1.0) - mixt_frac_zm((i)-1, (k)-1)) * mean_w_up_2nd((i)-1, (k)-1);
      });

  if (stats.flags().l_sample) {
    f2k::stats_update("mean_w_up", mean_w_up, stats);
    f2k::stats_update("mean_w_down", mean_w_down, stats);
  }
  return;
}

//=============================================================================
// nzm: Number of momentum levels [-]
// ngrdcol: Number of columns [-]
// w_i_zm: Mean of w                                   [m/s]
// varnce_w_i_zm: Variance of w                               [m^2/s^2]
// w_ref: Reference velocity, w|_ref (normally = 0)   [m/s]
// w_min: Min. velocity required to affect adj. level [m/s]
// mean_w_down_i: Mean w (<= w|_ref) from normal distribution [m/s]
// mean_w_up_i: Mean w (>= w|_ref) from normal distribution [m/s]
void calc_mean_w_up_down_component(int nzm, int ngrdcol,                     // In
                                   Real2D w_i_zm, Real2D varnce_w_i_zm,      // In
                                   real_t w_ref, Real2D w_min,               // In
                                   Real2D mean_w_down_i, Real2D mean_w_up_i) // Out
{
  // Description: This procedure is used to split the PDF component of
  //              vertical velocity into upward and downward components.
  //
  //              The method used is described in the description of
  //              mean_vert_vel_up_down, which calls this function.
  //
  // Notes: The calculation has been updated to optionally use intel's
  //        mkl_vml functions to allow vectorized calculations. Not all
  //        grid levels require expensive calculations though, so the
  //        strategy is as follows
  //           1. Keep track of which levels do need the calculation
  //           2. Store those the relavent values from those levels in
  //              a contigous array
  //           3. Perform vectorized calculation on contiguous arrays
  //              using mkl_vml functions
  //           4. Unpack results from contiguous array into output array
  //        Enabling this faster version requires compilation with MKL, by
  //        using -DMKL as a compiler flag
  //
  //-----------------------------------------------------------------------
  //------------------------- Input Variables -------------------------
  //------------------------- Output Variables -------------------------
  //------------------------- Local Variables -------------------------
  real_t erf_cache;      // erf/cdfnorm values
  real_t exp_cache;      // exp() values
  real_t sigma_w_i_zm;   // Standard deviation of w [m/s]
  real_t invrs_sqrt_2pi; // The inverse of sqrt(2*pi), calculated to save divide operations
  int i;
  int k; // Vertical loop index
  //------------------------- Begin Code -------------------------
  invrs_sqrt_2pi = one / sqrt_2pi;
  // Loop over momentum levels from 2 to nz-1.  Levels 1 and nz
  // are not needed.
  const int loop_1_extent_0 = std::max(0, (nzm - 1 - (2)) / (1) + 1);
  const int loop_1_extent_1 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "calc_mean_w_up_down_component::loop1", loop_1_extent_0 * loop_1_extent_1,
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_1_extent_1)) * (1);
        remaining /= (loop_1_extent_1);
        const int k = 2 + (remaining % (loop_1_extent_0)) * (1);
        remaining /= (loop_1_extent_0);
        real_t erf_cache;
        real_t exp_cache;
        real_t sigma_w_i_zm;
        // Standard deviation of w for the normal distribution.
        sigma_w_i_zm = Kokkos::sqrt(varnce_w_i_zm((i)-1, (k)-1));

        if (Kokkos::abs(w_i_zm((i)-1, (k)-1)) + real_t(3.0) * sigma_w_i_zm <= w_min((i)-1, (k)-1)) {
          // The entire normal is too weak to affect adjacent grid levels
          // w is considered to be 0 in both up and down directions
          mean_w_down_i((i)-1, (k)-1) = real_t(0.0);
          mean_w_up_i((i)-1, (k)-1) = real_t(0.0);
        } else {

          if (w_i_zm((i)-1, (k)-1) + real_t(3.) * sigma_w_i_zm <= w_ref) {
            // The entire normal is on the negative side of w|_ref.
            mean_w_down_i((i)-1, (k)-1) = w_i_zm((i)-1, (k)-1);
            mean_w_up_i((i)-1, (k)-1) = real_t(0.0);
          } else {

            if (w_i_zm((i)-1, (k)-1) - real_t(3.) * sigma_w_i_zm >= w_ref) {
              // The entire normal is on the positive side of w|_ref.
              mean_w_down_i((i)-1, (k)-1) = real_t(0.0);
              mean_w_up_i((i)-1, (k)-1) = w_i_zm((i)-1, (k)-1);
            } else {
              // The normal has significant values on both sides of w_ref.
              exp_cache =
                  Kokkos::exp(-((w_ref - w_i_zm((i)-1, (k)-1)) * (w_ref - w_i_zm((i)-1, (k)-1))) /
                              (real_t(2.0) * (sigma_w_i_zm * sigma_w_i_zm)));
              erf_cache = Kokkos::erf((w_ref - w_i_zm((i)-1, (k)-1)) / (sqrt_2 * sigma_w_i_zm));
              mean_w_down_i((i)-1, (k)-1) =
                  -sigma_w_i_zm * invrs_sqrt_2pi * exp_cache +
                  w_i_zm((i)-1, (k)-1) * real_t(0.5) * (real_t(1.0) + erf_cache);
              mean_w_up_i((i)-1, (k)-1) =
                  +sigma_w_i_zm * invrs_sqrt_2pi * exp_cache +
                  w_i_zm((i)-1, (k)-1) * real_t(0.5) * (real_t(1.0) - erf_cache);
            }
          }
        }
      });
  // Upper and lower levels are not used, set to 0 to besafe and avoid NaN problems
  const int loop_2_extent_0 = std::max(0, (ngrdcol - (1)) / (1) + 1);

  Kokkos::parallel_for(
      "calc_mean_w_up_down_component::loop2", loop_2_extent_0, KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i = 1 + (remaining % (loop_2_extent_0)) * (1);
        remaining /= (loop_2_extent_0);

        mean_w_down_i((i)-1, (1) - 1) = real_t(0.0);
        mean_w_up_i((i)-1, (1) - 1) = real_t(0.0);
        mean_w_down_i((i)-1, (nzm)-1) = real_t(0.0);
        mean_w_up_i((i)-1, (nzm)-1) = real_t(0.0);
      });
  return;
}
} // namespace
} // namespace clubb_kokkos
