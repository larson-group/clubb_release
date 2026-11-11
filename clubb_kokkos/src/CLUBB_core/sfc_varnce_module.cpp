// Surface velocity/scalar second moments and their Native budget updates.
// The public contract is in sfc_varnce_module.hpp; formulas and branches mirror
// src/CLUBB_core/sfc_varnce_module.F90, including its dormant Andre1978 scheme.
#include "scratch_workspace.hpp"
#include "numerical_check.hpp"
#include "sfc_varnce_module.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "stats_netcdf.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace clubb_kokkos {

void calc_sfc_varnce(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, int sclr_dim,
                     const sclr_idx_type& sclr_idx,                   // In
                     const grid& gr, real_t dt, Real1D sfc_elevation, // In
                     Real1D upwp_sfc, Real1D vpwp_sfc, Real2D wpthlp,
                     Real1D wprtp_sfc, // In
                     Real2D um, Real2D vm, Real2D Lscale_up,
                     Real2D wpsclrp_sfc,                               // In
                     Real2D lhs_splat_wp2, Real2D tau_zm,              // In
                     bool l_vary_convect_depth,                        // In
                     real_t T0,                                        // In
                     Real1D up2_sfc_coef,                              // In
                     Real1D a_const,                                   // In
                     stats_type& stats,                                // InOut
                     Real2D wp2, Real2D up2, Real2D vp2,               // InOut
                     Real2D thlp2, Real2D rtp2, Real2D rtpthlp,        // InOut
                     Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp, // InOut
                     err_info_type& err_info)                          // InOut
{
  // Description:
  // Estimate lower-boundary velocity/scalar variances and covariances from
  // surface fluxes and per-column parameters. Update only the lower momentum
  // level of the caller-owned fields, recording sampled budgets and errors.
  // References:
  // Andre, J. C., G. De Moor, P. Lacarrere, G. Therry, and R. Du Vachat, 1978.
  // Modeling the 24-Hour Evolution of the Mean and Turbulent Structures of
  // the Planetary Boundary Layer.  J. Atmos. Sci., 35, 1861 -- 1883.
  // -----------------------------------------------------------------------

  // Constant Parameters
  // Logical for Andre et al., 1978 parameterization.
  // Vince Larson increased ufmin to stabilize arm_97.  24 Jul 2007
  // ufmin = 0.0001_core_rknd, &
  // End Vince Larson's change.
  // Vince Larson changed in order to make correlations between [-1,1].  31 Jan 2008.
  // sclr_var_coef = 0.25_core_rknd, & ! This value is made up! - Vince Larson 12 Jul 2005
  // End Vince Larson's change
  // Vince Larson reduced surface spike in scalar variances associated
  // w/ Andre et al. 1978 scheme
  // -------------------------- Input Variables --------------------------
  // nzm/nzt: Momentum/thermodynamic level counts; ngrdcol: columns.
  // sclr_dim/sclr_idx: Passive scalar count and Native one-based emulation indices.
  // gr/sfc_elevation: Grid metadata and surface elevation [m].
  // dt: Main timestep [s]; l_vary_convect_depth selects diagnosed convective depth.
  // stats: Caller-owned Native sampled budget accumulator.
  // upwp_sfc: Surface u momentum flux, <u'w'>|_sfc   [m^2/s^2]
  // vpwp_sfc: Surface v momentum flux, <v'w'>|_sfc   [m^2/s^2]
  // wprtp_sfc: Surface moisture flux, <w'rt'>|_sfc    [kg/kg m/s]
  // wpthlp: Surface thetal flux, <w'thl'>|_sfc     [K m/s]
  // lhs_splat_wp2: LHS coefficient for vertical-velocity splatting [1/s].
  // tau_zm: Turbulent dissipation time at level zm(1)  [s]
  // um: Surface u wind component, <u>          [m/s]
  // vm: Surface v wind component, <v>          [m/s]
  // Lscale_up: Upward component of Lscale at surface  [m]
  // wpsclrp_sfc: Passive scalar flux, <w'sclr'>|_sfc   [units m/s]
  // T0: Reference temperature (usually 300)  [K]
  // up2_sfc_coef: CLUBB tunable parameter up2_sfc_coef   [-]
  // a_const: Coefficient in front of wp2, up2, and vp2
  // -------------------------- InOut Variables --------------------------
  // wp2: Surface variance of w, <w'^2>|_sfc            [m^2/s^2]
  // up2: Surface variance of u, <u'^2>|_sfc            [m^2/s^2]
  // vp2: Surface variance of v, <v'^2>|_sfc            [m^2/s^2]
  // thlp2: Surface variance of theta-l, <thl'^2>|_sfc    [K^2]
  // rtp2: Surface variance of rt, <rt'^2>|_sfc          [(kg/kg)^2]
  // rtpthlp: Surface covariance of rt and theta-l          [kg K/kg]
  // sclrp2: Surface variance of passive scalar            [units^2]
  // sclrprtp: Surface covariance of pssv scalar and rt  [units kg/kg]
  // sclrpthlp: Surface covariance of pssv scalar and theta-l [units K]
  // err_info: err_info struct containing err_code and err_header

  // Local storage coordinates are zero-based. Shared scalar-emulation IDs
  // normalize once at the metadata capture below.
  auto workspace_scope           = workspace.scope();
  constexpr bool l_andre_1978    = false;
  constexpr real_t z_const       = one;          // Defined height of 1 meter                [m]
  constexpr real_t ufmin         = real_t(0.01); // Minimum allowable value of u*   [m/s]
  constexpr real_t sclr_var_coef = real_t(0.4); // This value is made up! - Vince Larson 12 Jul 2005
  constexpr real_t reduce_coef   = real_t(0.2);
  Real1D uf                      = workspace_scope.view<Real1D>(
      "calc_sfc_varnce::uf", ngrdcol); // Surface friction vel., u*, in older version   [m/s]
  Real1D depth_pos_wpthlp = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "calc_sfc_varnce::depth_pos_wpthlp",
      ngrdcol); // Thickness of the layer near the surface with wpthlp > 0 [m]
  Real1D min_wp2_sfc_val =
      workspace_scope.view<Real1D>("calc_sfc_varnce::min_wp2_sfc_val",
                                   ngrdcol); // Minimum value of wp2_sfc that guarantees  [m^2/s^2]
  Real1D um_sfc_sqd = workspace_scope.view<Real1D>("calc_sfc_varnce::um_sfc_sqd",
                                                   ngrdcol); // Surface value of <u>^2 [m^2/s^2]
  Real1D vm_sfc_sqd = workspace_scope.view<Real1D>("calc_sfc_varnce::vm_sfc_sqd",
                                                   ngrdcol); // Surface value of <v>^2 [m^2/s^2]
  Real1D usp2_sfc   = workspace_scope.view<Real1D>(
      "calc_sfc_varnce::usp2_sfc",
      ngrdcol); // u_s (vector oriented w/ mean sfc. wind) variance [m^2/s^2]
  Real1D vsp2_sfc = workspace_scope.view<Real1D>(
      "calc_sfc_varnce::vsp2_sfc",
      ngrdcol); // v_s (vector perpen. to mean sfc. wind) variance  [m^2/s^2]
  Real1D ustar = workspace_scope.view<Real1D>("calc_sfc_varnce::ustar",
                                              ngrdcol); // Surface friction velocity, u* [m/s]
  Real1D zeta  = workspace_scope.view<Real1D>(
      "calc_sfc_varnce::zeta",
      ngrdcol); // Dimensionless height z_const/Lngth, where z_const = 1 m.  [-]
  Real1D wp2_splat_sfc_correction =
      workspace_scope.view<Real1D>("calc_sfc_varnce::wp2_splat_sfc_correction",
                                   ngrdcol); // Reduction in wp2_sfc due to splatting [m^2/s^2]
  int i{-1};
  // Capture component Views and zero-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm          = gr.k_lb_zm;
  const auto gr_zm               = gr.zm;
  const auto gr_grid_dir_indx    = gr.grid_dir_indx;
  const auto gr_k_lb_zt          = gr.k_lb_zt;
  const auto sclr_idx_iisclr_rt  = sclr_idx.iisclr_rt;
  const auto sclr_idx_iisclr_thl = sclr_idx.iisclr_thl;

  // Reflect surface varnce changes in budget
  if (stats.l_sample) {
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D thlp2_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::thlp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::thlp2_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column              = (cell % ngrdcol);
            const int level               = (cell / ngrdcol);
            thlp2_sf_begin(column, level) = thlp2(column, level) / dt;
          });
      stats_begin_budget("thlp2_sf", thlp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtp2_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::rtp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtp2_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column             = (cell % ngrdcol);
            const int level              = (cell / ngrdcol);
            rtp2_sf_begin(column, level) = rtp2(column, level) / dt;
          });
      stats_begin_budget("rtp2_sf", rtp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtpthlp_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::rtpthlp_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtpthlp_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column                = (cell % ngrdcol);
            const int level                 = (cell / ngrdcol);
            rtpthlp_sf_begin(column, level) = rtpthlp(column, level) / dt;
          });
      stats_begin_budget("rtpthlp_sf", rtpthlp_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D up2_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::up2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::up2_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            up2_sf_begin(column, level) = up2(column, level) / dt;
          });
      stats_begin_budget("up2_sf", up2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D vp2_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::vp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::vp2_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            vp2_sf_begin(column, level) = vp2(column, level) / dt;
          });
      stats_begin_budget("vp2_sf", vp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D wp2_sf_begin =
          workspace_scope.view<Real2D>("calc_sfc_varnce::wp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::wp2_sf_begin", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            wp2_sf_begin(column, level) = wp2(column, level) / dt;
          });
      stats_begin_budget("wp2_sf", wp2_sf_begin, stats);
    }
  }
  Kokkos::parallel_for(
      "calc_sfc_varnce::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int k{-1};
        // Find thickness of layer near surface with positive heat flux.
        // This is used when l_vary_convect_depth=.true. in order to determine wp2.
        if ((wpthlp(i, gr_k_lb_zm) <= zero)) {
          depth_pos_wpthlp(i) = one;
          // When sfc heat flux is negative, set depth to 1 m.
        } else {
          // When sfc heat flux is positive, march up sounding until wpthlp 1st becomes negative.
          k = gr_k_lb_zm;
          while (((wpthlp(i, k) > zero) && ((gr_zm(i, k) - sfc_elevation(i)) < real_t(1000.)))) {
            k = (gr_grid_dir_indx + k);
          }
          depth_pos_wpthlp(i) =
              (one > (gr_zm(i, k) - sfc_elevation(i)) ? one : (gr_zm(i, k) - sfc_elevation(i)));
        }
      });
  // a_const used to be set here; now it is a tunable parameter
  // if ( .not. l_vary_convect_depth ) then
  // a_const = 1.8_core_rknd
  // else
  // a_const = 0.6_core_rknd
  // end if
  if (l_andre_1978) {
    // Calculate <u>^2 and <v>^2.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_1", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i   = column;
          um_sfc_sqd(i) = (um(i, gr_k_lb_zt) * um(i, gr_k_lb_zt));
          vm_sfc_sqd(i) = (vm(i, gr_k_lb_zt) * vm(i, gr_k_lb_zt));
        });
    // Calculate surface friction velocity, u*.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_2", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          ustar(i) = (Kokkos::pow(((upwp_sfc(i) * upwp_sfc(i)) + (vpwp_sfc(i) * vpwp_sfc(i))),
                                  one_fourth) > ufmin
                          ? Kokkos::pow(((upwp_sfc(i) * upwp_sfc(i)) + (vpwp_sfc(i) * vpwp_sfc(i))),
                                        one_fourth)
                          : ufmin);
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_3", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          real_t Lngth{};
          if ((Kokkos::abs(wpthlp(i, gr_k_lb_zm)) > eps)) {
            // Find Monin-Obukhov Length (Andre et al., 1978, p. 1866).
            Lngth = (-(((ustar(i) * ustar(i)) * ustar(i))) /
                     (((real_t(0.35) * (one / T0)) * grav) * wpthlp(i, gr_k_lb_zm)));
            // Find the value of dimensionless height zeta
            // (Andre et al., 1978, p. 1866).
            zeta(i) = (z_const / Lngth);
          } else {
            // wpthlp = 0
            // The value of Monin-Obukhov length is +inf when ustar < 0 and -inf
            // when ustar > 0.  Either way, zeta = 0.
            zeta(i) = zero;
          }
          // wpthlp /= 0
        });
    // Andre et al, 1978, Eq. 29.
    // Notes:  1) "reduce_coef" is a reduction coefficient intended to make
    // the values of rtp2, thlp2, and rtpthlp smaller at the
    // surface.
    // 2) With the reduction coefficient having a value of 0.2, the
    // surface correlations of both w & rt and w & thl have a value
    // of about 0.845.  These correlations are greater if zeta < 0.
    // The correlations have a value greater than 1 if
    // zeta <= -0.212.
    // 3) The surface correlation of rt & thl is 1.
    // Brian Griffin; February 2, 2008.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_4", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          if ((zeta(i) < zero)) {
            thlp2(i, gr_k_lb_zm) =
                (((reduce_coef *
                   ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i)))) *
                  four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
            rtp2(i, gr_k_lb_zm) =
                (((reduce_coef * ((wprtp_sfc(i) * wprtp_sfc(i)) / (ustar(i) * ustar(i)))) * four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
            rtpthlp(i, gr_k_lb_zm) =
                (((reduce_coef * ((wprtp_sfc(i) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i)))) *
                  four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
            wp2(i, gr_k_lb_zm) = ((ustar(i) * ustar(i)) *
                                  (real_t(1.75) + (two * Kokkos::pow(-(zeta(i)), two_thirds))));
          } else {
            thlp2(i, gr_k_lb_zm) =
                ((reduce_coef * four) *
                 ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i))));
            rtp2(i, gr_k_lb_zm) =
                ((reduce_coef * four) * ((wprtp_sfc(i) * wprtp_sfc(i)) / (ustar(i) * ustar(i))));
            rtpthlp(i, gr_k_lb_zm) =
                ((reduce_coef * four) *
                 ((wprtp_sfc(i) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i))));
            wp2(i, gr_k_lb_zm) = (real_t(1.75) * (ustar(i) * ustar(i)));
          }
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_5", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          thlp2(i, gr_k_lb_zm) =
              ((thl_tol * thl_tol) > thlp2(i, gr_k_lb_zm) ? (thl_tol * thl_tol)
                                                          : thlp2(i, gr_k_lb_zm));
          rtp2(i, gr_k_lb_zm) =
              ((rt_tol * rt_tol) > rtp2(i, gr_k_lb_zm) ? (rt_tol * rt_tol) : rtp2(i, gr_k_lb_zm));
        });
    // Calculate wstar following Andre et al., 1978, p. 1866.
    // w* = ( ( 1 / T0 ) * g * <w'thl'>|_sfc * z_i )^(1/3);
    // where z_i is the height of the mixed layer.  The value of CLUBB's
    // upward component of mixing length, Lscale_up, at the surface will be
    // used as z_i.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_6", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          real_t wstar{};
          wstar = Kokkos::pow(
              ((((one / T0) * grav) * wpthlp(i, gr_k_lb_zm)) * Lscale_up(i, gr_k_lb_zt)),
              one_third);
          // Andre et al., 1978, Eq. 29.
          // Andre et al. (1978) defines horizontal wind surface variances in terms
          // of orientation with the mean surface wind.  The vector u_s is the wind
          // vector oriented with the mean surface wind.  The vector v_s is the wind
          // vector oriented perpendicular to the mean surface wind.  Thus, <u_s> is
          // equal to the mean surface wind (both in speed and direction), and <v_s>
          // is 0.  Equation 29 gives the formula for the variance of u_s, which is
          // <u_s'^2> (usp2_sfc in the code), and the formula for the variance of
          // v_s, which is <v_s'^2> (vsp2_sfc in the code).
          if ((wpthlp(i, gr_k_lb_zm) > zero)) {
            usp2_sfc(i) = ((four * (ustar(i) * ustar(i))) + (real_t(0.3) * (wstar * wstar)));
            vsp2_sfc(i) =
                ((real_t(1.75) * (ustar(i) * ustar(i))) + (real_t(0.3) * (wstar * wstar)));
          } else {
            usp2_sfc(i) = (four * (ustar(i) * ustar(i)));
            vsp2_sfc(i) = (real_t(1.75) * (ustar(i) * ustar(i)));
          }
        });
    // Add effect of vertical compression of eddies on horizontal gustiness.
    // First, ensure that wp2 does not make the correlation
    // of (w,thl) or (w,rt) outside (-1,1).
    // Perhaps in the future we should also ensure that the correlations
    // of (w,u) and (w,v) are not outside (-1,1).
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_7", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          min_wp2_sfc_val(i) =
              ((w_tol_sqd > ((wprtp_sfc(i) * wprtp_sfc(i)) /
                             (rtp2(i, gr_k_lb_zm) *
                              (max_mag_correlation_flux * max_mag_correlation_flux)))
                    ? w_tol_sqd
                    : ((wprtp_sfc(i) * wprtp_sfc(i)) /
                       (rtp2(i, gr_k_lb_zm) *
                        (max_mag_correlation_flux * max_mag_correlation_flux)))) >
                       ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) /
                        (thlp2(i, gr_k_lb_zm) *
                         (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? (w_tol_sqd > ((wprtp_sfc(i) * wprtp_sfc(i)) /
                                   (rtp2(i, gr_k_lb_zm) *
                                    (max_mag_correlation_flux * max_mag_correlation_flux)))
                          ? w_tol_sqd
                          : ((wprtp_sfc(i) * wprtp_sfc(i)) /
                             (rtp2(i, gr_k_lb_zm) *
                              (max_mag_correlation_flux * max_mag_correlation_flux))))
                   : ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) /
                      (thlp2(i, gr_k_lb_zm) *
                       (max_mag_correlation_flux * max_mag_correlation_flux))));
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_8", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          if (((wp2(i, gr_k_lb_zm) - ((tau_zm(i, gr_k_lb_zm) * lhs_splat_wp2(i, gr_k_lb_zm)) *
                                      wp2(i, gr_k_lb_zm))) < min_wp2_sfc_val(i))) {
            // if ( wp2(i,gr%k_lb_zm) + tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm) &
            // < min_wp2_sfc_val(i) ) then
            // splatting correction drives wp2 to overly small value
            wp2_splat_sfc_correction(i) = (-(wp2(i, gr_k_lb_zm)) + min_wp2_sfc_val(i));
            wp2(i, gr_k_lb_zm)          = min_wp2_sfc_val(i);
          } else {
            wp2_splat_sfc_correction(i) =
                ((tau_zm(i, gr_k_lb_zm) * lhs_splat_wp2(i, gr_k_lb_zm)) * wp2(i, gr_k_lb_zm));
            // wp2_splat_sfc_correction(i) = tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm)
            wp2(i, gr_k_lb_zm) = (wp2(i, gr_k_lb_zm) + wp2_splat_sfc_correction(i));
          }
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_9", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          usp2_sfc(i) = (usp2_sfc(i) - (real_t(0.5) * wp2_splat_sfc_correction(i)));
          vsp2_sfc(i) = (vsp2_sfc(i) - (real_t(0.5) * wp2_splat_sfc_correction(i)));
        });
    // Variance of u, <u'^2>, at the surface can be found from <u_s'^2>,
    // <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
    // <u'^2>|_sfc = <u_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
    // + <v_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
    // where <u>^2 + <v>^2 /= 0.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_10", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          up2(i, gr_k_lb_zm) =
              ((usp2_sfc(i) *
                (um_sfc_sqd(i) /
                 ((um_sfc_sqd(i) + vm_sfc_sqd(i)) > eps ? (um_sfc_sqd(i) + vm_sfc_sqd(i)) : eps))) +
               (vsp2_sfc(i) *
                (vm_sfc_sqd(i) /
                 ((um_sfc_sqd(i) + vm_sfc_sqd(i)) > eps ? (um_sfc_sqd(i) + vm_sfc_sqd(i)) : eps))));
        });
    // Variance of v, <v'^2>, at the surface can be found from <u_s'^2>,
    // <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
    // <v'^2>|_sfc = <v_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
    // + <u_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
    // where <u>^2 + <v>^2 /= 0.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_11", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          vp2(i, gr_k_lb_zm) =
              ((vsp2_sfc(i) *
                (um_sfc_sqd(i) /
                 ((um_sfc_sqd(i) + vm_sfc_sqd(i)) > eps ? (um_sfc_sqd(i) + vm_sfc_sqd(i)) : eps))) +
               (usp2_sfc(i) *
                (vm_sfc_sqd(i) /
                 ((um_sfc_sqd(i) + vm_sfc_sqd(i)) > eps ? (um_sfc_sqd(i) + vm_sfc_sqd(i)) : eps))));
        });
    // Passive scalars
    if ((sclr_dim > 0)) {
      {
        const int level_count = sclr_dim;
        Kokkos::parallel_for(
            "calc_sfc_varnce::cells_12", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
              const int i    = (cell % ngrdcol);
              const int sclr = (cell / ngrdcol);
              // Notes:  1) "reduce_coef" is a reduction coefficient intended to
              // make the values of sclrprtp, sclrpthlp, and sclrp2
              // smaller at the surface.
              // 2) With the reduction coefficient having a value of 0.2,
              // the surface correlation of w & sclr has a value of
              // about 0.845.  The correlation is greater if zeta < 0.
              // The correlation has a value greater than 1 if
              // zeta <= -0.212.
              // 3) The surface correlations of both rt & sclr and
              // thl & sclr are 1.
              // Brian Griffin; February 2, 2008.
              if ((zeta(i) < zero)) {
                sclrprtp(i, gr_k_lb_zm, sclr) =
                    (((reduce_coef *
                       ((wpsclrp_sfc(i, sclr) * wprtp_sfc(i)) / (ustar(i) * ustar(i)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
                sclrpthlp(i, gr_k_lb_zm, sclr) =
                    (((reduce_coef *
                       ((wpsclrp_sfc(i, sclr) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
                sclrp2(i, gr_k_lb_zm, sclr) =
                    (((reduce_coef *
                       ((wpsclrp_sfc(i, sclr) * wpsclrp_sfc(i, sclr)) / (ustar(i) * ustar(i)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i))), -(two_thirds)));
              } else {
                sclrprtp(i, gr_k_lb_zm, sclr) =
                    ((reduce_coef * four) *
                     ((wpsclrp_sfc(i, sclr) * wprtp_sfc(i)) / (ustar(i) * ustar(i))));
                sclrpthlp(i, gr_k_lb_zm, sclr) =
                    ((reduce_coef * four) *
                     ((wpsclrp_sfc(i, sclr) * wpthlp(i, gr_k_lb_zm)) / (ustar(i) * ustar(i))));
                sclrp2(i, gr_k_lb_zm, sclr) =
                    ((reduce_coef * four) *
                     ((wpsclrp_sfc(i, sclr) * wpsclrp_sfc(i, sclr)) / (ustar(i) * ustar(i))));
              }
            });
      }
      // sclr = 1, sclr_dim
    }
  } else {
    // Previous code.
    // Compute ustar^2
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_13", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          real_t ustar2{};
          real_t wstar{};
          ustar2 = Kokkos::sqrt(((upwp_sfc(i) * upwp_sfc(i)) + (vpwp_sfc(i) * vpwp_sfc(i))));
          // Compute wstar following Andre et al., 1976
          if ((wpthlp(i, gr_k_lb_zm) > zero)) {
            if (!l_vary_convect_depth) {
              wstar =
                  Kokkos::pow(((((one / T0) * grav) * wpthlp(i, gr_k_lb_zm)) * z_const), one_third);
            } else {
              wstar = Kokkos::pow((((((one / T0) * grav) * wpthlp(i, gr_k_lb_zm)) * real_t(0.2)) *
                                   depth_pos_wpthlp(i)),
                                  one_third);
            }
          } else {
            wstar = zero;
          }
          // Surface friction velocity following Andre et al. 1978
          if (!l_vary_convect_depth) {
            uf(i) = Kokkos::sqrt((ustar2 + ((real_t(0.3) * wstar) * wstar)));
          } else {
            uf(i) = Kokkos::sqrt((ustar2 + (wstar * wstar)));
          }
          uf(i) = (ufmin > uf(i) ? ufmin : uf(i));
        });
    // Compute estimate for surface second order moments
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_14", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i        = column;
          wp2(i, gr_k_lb_zm) = (a_const(i) * (uf(i) * uf(i)));
          up2(i, gr_k_lb_zm) = ((up2_sfc_coef(i) * a_const(i)) * (uf(i) * uf(i)));
          // From Andre, et al. 1978
          vp2(i, gr_k_lb_zm) = ((up2_sfc_coef(i) * a_const(i)) * (uf(i) * uf(i)));
          // "  "
        });
    // Notes:  1) With "a" having a value of 1.8, the surface correlations of
    // both w & rt and w & thl have a value of about 0.878.
    // 2) The surface correlation of rt & thl is 0.5.
    // Brian Griffin; February 2, 2008.
    if (!l_vary_convect_depth) {
      Kokkos::parallel_for(
          "calc_sfc_varnce::columns_15", ngrdcol, KOKKOS_LAMBDA(const int column) {
            const int i = column;
            // Round each Native power base before squaring. NVHPC otherwise
            // rewrites repeated C++ quotients as flux^2/uf^2 (BFB-015).
            const real_t wpthlp_over_uf = wpthlp(i, gr_k_lb_zm) / uf(i);
            const real_t wprtp_over_uf  = wprtp_sfc(i) / uf(i);
            thlp2(i, gr_k_lb_zm) = (real_t(0.4) * a_const(i)) * (wpthlp_over_uf * wpthlp_over_uf);
            rtp2(i, gr_k_lb_zm)  = (real_t(0.4) * a_const(i)) * (wprtp_over_uf * wprtp_over_uf);
            rtpthlp(i, gr_k_lb_zm) =
                (((real_t(0.2) * a_const(i)) * (wpthlp(i, gr_k_lb_zm) / uf(i))) *
                 (wprtp_sfc(i) / uf(i)));
          });
    } else {
      Kokkos::parallel_for(
          "calc_sfc_varnce::columns_16", ngrdcol, KOKKOS_LAMBDA(const int column) {
            const int i = column;
            // Preserve the same Native quotient-power boundary (BFB-015).
            const real_t wpthlp_over_uf = wpthlp(i, gr_k_lb_zm) / uf(i);
            const real_t wprtp_over_uf  = wprtp_sfc(i) / uf(i);
            thlp2(i, gr_k_lb_zm) =
                (wpthlp_over_uf * wpthlp_over_uf) /
                ((max_mag_correlation_flux * max_mag_correlation_flux) * a_const(i));
            rtp2(i, gr_k_lb_zm) =
                (wprtp_over_uf * wprtp_over_uf) /
                ((max_mag_correlation_flux * max_mag_correlation_flux) * a_const(i));
            rtpthlp(i, gr_k_lb_zm) = (max_mag_correlation_flux *
                                      Kokkos::sqrt((thlp2(i, gr_k_lb_zm) * rtp2(i, gr_k_lb_zm))));
          });
    }
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_17", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          thlp2(i, gr_k_lb_zm) =
              ((thl_tol * thl_tol) > thlp2(i, gr_k_lb_zm) ? (thl_tol * thl_tol)
                                                          : thlp2(i, gr_k_lb_zm));
          rtp2(i, gr_k_lb_zm) =
              ((rt_tol * rt_tol) > rtp2(i, gr_k_lb_zm) ? (rt_tol * rt_tol) : rtp2(i, gr_k_lb_zm));
        });
    // Add effect of vertical compression of eddies on horizontal gustiness.
    // First, ensure that wp2 does not make the correlation
    // of (w,thl) or (w,rt) outside (-1,1).
    // Perhaps in the future we should also ensure that the correlations
    // of (w,u) and (w,v) are not outside (-1,1).
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_18", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          min_wp2_sfc_val(i) =
              ((w_tol_sqd > ((wprtp_sfc(i) * wprtp_sfc(i)) /
                             (rtp2(i, gr_k_lb_zm) *
                              (max_mag_correlation_flux * max_mag_correlation_flux)))
                    ? w_tol_sqd
                    : ((wprtp_sfc(i) * wprtp_sfc(i)) /
                       (rtp2(i, gr_k_lb_zm) *
                        (max_mag_correlation_flux * max_mag_correlation_flux)))) >
                       ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) /
                        (thlp2(i, gr_k_lb_zm) *
                         (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? (w_tol_sqd > ((wprtp_sfc(i) * wprtp_sfc(i)) /
                                   (rtp2(i, gr_k_lb_zm) *
                                    (max_mag_correlation_flux * max_mag_correlation_flux)))
                          ? w_tol_sqd
                          : ((wprtp_sfc(i) * wprtp_sfc(i)) /
                             (rtp2(i, gr_k_lb_zm) *
                              (max_mag_correlation_flux * max_mag_correlation_flux))))
                   : ((wpthlp(i, gr_k_lb_zm) * wpthlp(i, gr_k_lb_zm)) /
                      (thlp2(i, gr_k_lb_zm) *
                       (max_mag_correlation_flux * max_mag_correlation_flux))));
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_19", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          if (((wp2(i, gr_k_lb_zm) - ((tau_zm(i, gr_k_lb_zm) * lhs_splat_wp2(i, gr_k_lb_zm)) *
                                      wp2(i, gr_k_lb_zm))) < min_wp2_sfc_val(i))) {
            // if ( wp2(i,gr%k_lb_zm) + tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm) &
            // < min_wp2_sfc_val(i) ) then
            // splatting correction drives wp2 to overly small values
            wp2_splat_sfc_correction(i) = (-(wp2(i, gr_k_lb_zm)) + min_wp2_sfc_val(i));
            wp2(i, gr_k_lb_zm)          = min_wp2_sfc_val(i);
          } else {
            wp2_splat_sfc_correction(i) =
                ((tau_zm(i, gr_k_lb_zm) * lhs_splat_wp2(i, gr_k_lb_zm)) * wp2(i, gr_k_lb_zm));
            // wp2_splat_sfc_correction(i) = tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm)
            wp2(i, gr_k_lb_zm) = (wp2(i, gr_k_lb_zm) + wp2_splat_sfc_correction(i));
          }
          up2(i, gr_k_lb_zm) = (up2(i, gr_k_lb_zm) - (real_t(0.5) * wp2_splat_sfc_correction(i)));
          vp2(i, gr_k_lb_zm) = (vp2(i, gr_k_lb_zm) - (real_t(0.5) * wp2_splat_sfc_correction(i)));
        });
    // Passive scalars
    if ((sclr_dim > 0)) {
      {
        const int level_count = sclr_dim;
        Kokkos::parallel_for(
            "calc_sfc_varnce::cells_20", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
              const int i    = (cell % ngrdcol);
              const int sclr = (cell / ngrdcol);
              // Vince Larson changed coeffs to make correlations between [-1,1].
              // 31 Jan 2008
              // sclrprtp(i) &
              // = a * (wprtp_sfc / uf) * (wpsclrp(i) / uf)
              // sclrpthlp(i) &
              // = a * (wpthlp / uf) * (wpsclrp(i) / uf)
              // sclrp2(i) &
              // = sclr_var_coef * a * ( wpsclrp(i) / uf )**2
              // Notes:  1) With "a" having a value of 1.8 and "sclr_var_coef"
              // having a value of 0.4, the surface correlation of
              // w & sclr has a value of about 0.878.
              // 2) With "sclr_var_coef" having a value of 0.4, the
              // surface correlations of both rt & sclr and
              // thl & sclr are 0.5.
              // Brian Griffin; February 2, 2008.
              // We use the following if..then's to make sclr_rt and sclr_thl
              // close to the actual thlp2/rtp2 at the surface.
              // -dschanen 25 Sep 08
              if ((sclr == sclr_idx_iisclr_rt)) {
                // If we are trying to emulate rt with the scalar, then we
                // use the variance coefficient from above
                sclrprtp(i, gr_k_lb_zm, sclr) =
                    (((real_t(0.4) * a_const(i)) * (wprtp_sfc(i) / uf(i))) *
                     (wpsclrp_sfc(i, sclr) / uf(i)));
              } else {
                sclrprtp(i, gr_k_lb_zm, sclr) =
                    (((real_t(0.2) * a_const(i)) * (wprtp_sfc(i) / uf(i))) *
                     (wpsclrp_sfc(i, sclr) / uf(i)));
              }
              if ((sclr == sclr_idx_iisclr_thl)) {
                // As above, but for thetal
                sclrpthlp(i, gr_k_lb_zm, sclr) =
                    (((real_t(0.4) * a_const(i)) * (wpthlp(i, gr_k_lb_zm) / uf(i))) *
                     (wpsclrp_sfc(i, sclr) / uf(i)));
              } else {
                sclrpthlp(i, gr_k_lb_zm, sclr) =
                    (((real_t(0.2) * a_const(i)) * (wpthlp(i, gr_k_lb_zm) / uf(i))) *
                     (wpsclrp_sfc(i, sclr) / uf(i)));
              }
              // Preserve the rounded Native power base for passive scalars too (BFB-015).
              const real_t wpsclrp_over_uf = wpsclrp_sfc(i, sclr) / uf(i);
              sclrp2(i, gr_k_lb_zm, sclr) =
                  (sclr_var_coef * a_const(i)) * (wpsclrp_over_uf * wpsclrp_over_uf);
              // End Vince Larson's change.
            });
      }
      // 1,...sclr_dim
    }
    // sclr_dim > 0
  }
  // l_andre_1978
  // Clip wp2 at wp2_max, same as in advance_wp2_wp3
  Kokkos::parallel_for(
      "calc_sfc_varnce::columns_21", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i        = column;
        wp2(i, gr_k_lb_zm) = (wp2(i, gr_k_lb_zm) < wp2_max ? wp2(i, gr_k_lb_zm) : wp2_max);
      });
  Kokkos::parallel_for(
      "calc_sfc_varnce::columns_22", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        if ((Kokkos::abs((gr_zm(i, gr_k_lb_zm) - sfc_elevation(i))) >
             ((Kokkos::abs((gr_zm(i, gr_k_lb_zm) + sfc_elevation(i))) * eps) / 2))) {
          // Variances for cases where the lowest level is not at the surface.
          // Eliminate surface effects on lowest level variances.
          wp2(i, gr_k_lb_zm)     = w_tol_sqd;
          up2(i, gr_k_lb_zm)     = w_tol_sqd;
          vp2(i, gr_k_lb_zm)     = w_tol_sqd;
          thlp2(i, gr_k_lb_zm)   = (thl_tol * thl_tol);
          rtp2(i, gr_k_lb_zm)    = (rt_tol * rt_tol);
          rtpthlp(i, gr_k_lb_zm) = real_t(0.0);
        }
        // gr%zm(1,gr%k_lb_zm) == sfc_elevation
      });
  if ((sclr_dim > 0)) {
    {
      const int level_count = sclr_dim;
      Kokkos::parallel_for(
          "calc_sfc_varnce::cells_23", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int sclr = (cell / ngrdcol);
            if ((Kokkos::abs((gr_zm(i, gr_k_lb_zm) - sfc_elevation(i))) >
                 ((Kokkos::abs((gr_zm(i, gr_k_lb_zm) + sfc_elevation(i))) * eps) / 2))) {
              // Variances for cases where the lowest level is not at the surface.
              // Eliminate surface effects on lowest level variances.
              sclrp2(i, gr_k_lb_zm, sclr)    = real_t(0.0);
              sclrprtp(i, gr_k_lb_zm, sclr)  = real_t(0.0);
              sclrpthlp(i, gr_k_lb_zm, sclr) = real_t(0.0);
            }
          });
    }
  }
  if (stats.l_sample) {
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D thlp2_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::thlp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::thlp2_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column              = (cell % ngrdcol);
            const int level               = (cell / ngrdcol);
            thlp2_sf_final(column, level) = thlp2(column, level) / dt;
          });
      stats_finalize_budget("thlp2_sf", thlp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtp2_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::rtp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtp2_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column             = (cell % ngrdcol);
            const int level              = (cell / ngrdcol);
            rtp2_sf_final(column, level) = rtp2(column, level) / dt;
          });
      stats_finalize_budget("rtp2_sf", rtp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtpthlp_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::rtpthlp_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtpthlp_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column                = (cell % ngrdcol);
            const int level                 = (cell / ngrdcol);
            rtpthlp_sf_final(column, level) = rtpthlp(column, level) / dt;
          });
      stats_finalize_budget("rtpthlp_sf", rtpthlp_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D up2_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::up2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::up2_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            up2_sf_final(column, level) = up2(column, level) / dt;
          });
      stats_finalize_budget("up2_sf", up2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D vp2_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::vp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::vp2_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            vp2_sf_final(column, level) = vp2(column, level) / dt;
          });
      stats_finalize_budget("vp2_sf", vp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D wp2_sf_final =
          workspace_scope.view<Real2D>("calc_sfc_varnce::wp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::wp2_sf_final", (ngrdcol * nzm), KOKKOS_LAMBDA(const int cell) {
            const int column            = (cell % ngrdcol);
            const int level             = (cell / ngrdcol);
            wp2_sf_final(column, level) = wp2(column, level) / dt;
          });
      stats_finalize_budget("wp2_sf", wp2_sf_final, stats);
    }
  }
  if (clubb_at_least_debug_level_api(2)) {
    // Native diagnostics run on host. Stage each checked field once here;
    // preserve the Native per-column checker and its all-column fatal semantics.
    const auto wp2_debug_host       = to_host_space(workspace, wp2);
    const auto up2_debug_host       = to_host_space(workspace, up2);
    const auto vp2_debug_host       = to_host_space(workspace, vp2);
    const auto thlp2_debug_host     = to_host_space(workspace, thlp2);
    const auto rtp2_debug_host      = to_host_space(workspace, rtp2);
    const auto rtpthlp_debug_host   = to_host_space(workspace, rtpthlp);
    const auto sclrp2_debug_host    = to_host_space(workspace, sclrp2);
    const auto sclrprtp_debug_host  = to_host_space(workspace, sclrprtp);
    const auto sclrpthlp_debug_host = to_host_space(workspace, sclrpthlp);
    for (i = 0; i < ngrdcol; i += 1) {
      sfc_varnce_check(workspace, sclr_dim, wp2_debug_host(i, gr_k_lb_zm),
                       up2_debug_host(i, gr_k_lb_zm), vp2_debug_host(i, gr_k_lb_zm),
                       thlp2_debug_host(i, gr_k_lb_zm), rtp2_debug_host(i, gr_k_lb_zm),
                       rtpthlp_debug_host(i, gr_k_lb_zm),
                       Kokkos::subview(sclrp2_debug_host, i, gr_k_lb_zm, Kokkos::ALL),
                       Kokkos::subview(sclrprtp_debug_host, i, gr_k_lb_zm, Kokkos::ALL),
                       Kokkos::subview(sclrpthlp_debug_host, i, gr_k_lb_zm, Kokkos::ALL), err_info);
      int column_status;
      Kokkos::deep_copy(column_status, Kokkos::subview(err_info.err_code, i));
      if ((column_status == clubb_fatal_error)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          const auto header_host = to_host_space(workspace, err_info.err_header);
          diagnostic.write(header_host.data() + (header_host.extent(0) * i), header_host.extent(0));
          std::cerr << diagnostic.str() << '\n';
        }
      }
    }
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Error in calc_sfc_varnce" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Intent(in)" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "upwp_sfc = " << ' ';
        {
          const auto field_host = to_host_space(workspace, upwp_sfc);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "vpwp_sfc = " << ' ';
        {
          const auto field_host = to_host_space(workspace, vpwp_sfc);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wpthlp = " << ' ';
        {
          const auto field_host = to_host_space(workspace, wpthlp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wprtp_sfc = " << ' ';
        {
          const auto field_host = to_host_space(workspace, wprtp_sfc);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      if ((sclr_dim > 0)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "wpsclrp = " << ' ';
          {
            const auto field_host = to_host_space(workspace, wpsclrp_sfc);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          std::cerr << diagnostic.str() << '\n';
        }
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Intent(out)" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2 = " << ' ';
        {
          const auto field_host = to_host_space(workspace, wp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "up2 = " << ' ';
        {
          const auto field_host = to_host_space(workspace, up2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "vp2 = " << ' ';
        {
          const auto field_host = to_host_space(workspace, vp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "thlp2 = " << ' ';
        {
          const auto field_host = to_host_space(workspace, thlp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "rtp2 = " << ' ';
        {
          const auto field_host = to_host_space(workspace, rtp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "rtpthlp = " << ' ';
        {
          const auto field_host = to_host_space(workspace, rtpthlp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      if ((sclr_dim > 0)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "sclrp2 = " << ' ';
          {
            const auto field_host = to_host_space(workspace, sclrp2);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "sclrprtp = " << ' ';
          {
            const auto field_host = to_host_space(workspace, sclrprtp);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "sclrpthlp = " << ' ';
          {
            const auto field_host = to_host_space(workspace, sclrpthlp);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          std::cerr << diagnostic.str() << '\n';
        }
      }
    }
    // err_code == clubb_fatal_error
  }
  // clubb_at_least_debug_level_api ( 2 )! Update surface stats
  return;
}

} // namespace clubb_kokkos
