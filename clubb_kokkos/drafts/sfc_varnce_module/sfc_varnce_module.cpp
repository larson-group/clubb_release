// Surface velocity/scalar second moments and their Native budget updates.
// The public contract is in sfc_varnce_module.hpp; formulas and branches mirror
// src/CLUBB_core/sfc_varnce_module.F90, including its dormant Andre1978 scheme.
#include "sfc_varnce_module.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace clubb_kokkos {

void calc_sfc_varnce(int nzm, int nzt, int ngrdcol, int sclr_dim, const sclr_idx_type& sclr_idx,
                     const grid& gr, real_t dt, Real1D sfc_elevation, Real1D upwp_sfc,
                     Real1D vpwp_sfc, Real2D wpthlp, Real1D wprtp_sfc, Real2D um, Real2D vm,
                     Real2D Lscale_up, Real2D wpsclrp_sfc, Real2D lhs_splat_wp2, Real2D tau_zm,
                     bool l_vary_convect_depth, real_t T0, Real1D up2_sfc_coef, Real1D a_const,
                     f2k::stats_type& stats, Real2D wp2, Real2D up2, Real2D vp2, Real2D thlp2,
                     Real2D rtp2, Real2D rtpthlp, Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp,
                     err_info_type& err_info) {
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
  // wp2_splat,    & ! Tendency of <w'^2> due to splatting of eddies at zm(1) [m^2/s^3]
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

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one, including scalar-emulation indices.
  constexpr bool l_andre_1978 = false;
  constexpr real_t z_const = one;               // Defined height of 1 meter                [m]
  constexpr real_t ufmin = real_t(0.01);        // Minimum allowable value of u*   [m/s]
  constexpr real_t sclr_var_coef = real_t(0.4); // This value is made up! - Vince Larson 12 Jul 2005
  constexpr real_t reduce_coef = real_t(0.2);
  Real1D uf("calc_sfc_varnce::uf", ngrdcol); // Surface friction vel., u*, in older version   [m/s]
  Real1D depth_pos_wpthlp("calc_sfc_varnce::depth_pos_wpthlp",
                          ngrdcol); // Thickness of the layer near the surface with wpthlp > 0 [m]
  Real1D min_wp2_sfc_val("calc_sfc_varnce::min_wp2_sfc_val",
                         ngrdcol); // Minimum value of wp2_sfc that guarantees  [m^2/s^2]
  Real1D um_sfc_sqd("calc_sfc_varnce::um_sfc_sqd", ngrdcol); // Surface value of <u>^2 [m^2/s^2]
  Real1D vm_sfc_sqd("calc_sfc_varnce::vm_sfc_sqd", ngrdcol); // Surface value of <v>^2 [m^2/s^2]
  Real1D usp2_sfc("calc_sfc_varnce::usp2_sfc",
                  ngrdcol); // u_s (vector oriented w/ mean sfc. wind) variance [m^2/s^2]
  Real1D vsp2_sfc("calc_sfc_varnce::vsp2_sfc",
                  ngrdcol); // v_s (vector perpen. to mean sfc. wind) variance  [m^2/s^2]
  Real1D ustar("calc_sfc_varnce::ustar", ngrdcol); // Surface friction velocity, u* [m/s]
  Real1D zeta("calc_sfc_varnce::zeta",
              ngrdcol); // Dimensionless height z_const/Lngth, where z_const = 1 m.  [-]
  Real1D wp2_splat_sfc_correction("calc_sfc_varnce::wp2_splat_sfc_correction",
                                  ngrdcol); // Reduction in wp2_sfc due to splatting [m^2/s^2]
  int i{};
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_zm = gr.zm;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto sclr_idx_iisclr_rt = sclr_idx.iisclr_rt;
  const auto sclr_idx_iisclr_thl = sclr_idx.iisclr_thl;

  // Reflect surface varnce changes in budget
  if (stats.flags().l_sample) {
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D thlp2_sf_begin("calc_sfc_varnce::thlp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::thlp2_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            thlp2_sf_begin(column, level) = thlp2(column, level) / dt;
          });
      f2k::stats_begin_budget("thlp2_sf", thlp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtp2_sf_begin("calc_sfc_varnce::rtp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtp2_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            rtp2_sf_begin(column, level) = rtp2(column, level) / dt;
          });
      f2k::stats_begin_budget("rtp2_sf", rtp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtpthlp_sf_begin("calc_sfc_varnce::rtpthlp_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtpthlp_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            rtpthlp_sf_begin(column, level) = rtpthlp(column, level) / dt;
          });
      f2k::stats_begin_budget("rtpthlp_sf", rtpthlp_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D up2_sf_begin("calc_sfc_varnce::up2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::up2_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            up2_sf_begin(column, level) = up2(column, level) / dt;
          });
      f2k::stats_begin_budget("up2_sf", up2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D vp2_sf_begin("calc_sfc_varnce::vp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::vp2_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            vp2_sf_begin(column, level) = vp2(column, level) / dt;
          });
      f2k::stats_begin_budget("vp2_sf", vp2_sf_begin, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D wp2_sf_begin("calc_sfc_varnce::wp2_sf_begin", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::wp2_sf_begin", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            wp2_sf_begin(column, level) = wp2(column, level) / dt;
          });
      f2k::stats_begin_budget("wp2_sf", wp2_sf_begin, stats);
    }
  }
  Kokkos::parallel_for(
      "calc_sfc_varnce::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column + 1;
        int k{};
        // Find thickness of layer near surface with positive heat flux.
        // This is used when l_vary_convect_depth=.true. in order to determine wp2.
        if ((wpthlp(i - 1, gr_k_lb_zm - 1) <= zero)) {
          depth_pos_wpthlp(i - 1) = one;
          // When sfc heat flux is negative, set depth to 1 m.
        } else {
          // When sfc heat flux is positive, march up sounding until wpthlp 1st becomes negative.
          k = gr_k_lb_zm;
          while (((wpthlp(i - 1, k - 1) > zero) &&
                  ((gr_zm(i - 1, k - 1) - sfc_elevation(i - 1)) < real_t(1000.)))) {
            k = (k + gr_grid_dir_indx);
          }
          depth_pos_wpthlp(i - 1) = (one > (gr_zm(i - 1, k - 1) - sfc_elevation(i - 1))
                                         ? one
                                         : (gr_zm(i - 1, k - 1) - sfc_elevation(i - 1)));
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
          const int i = column + 1;
          um_sfc_sqd(i - 1) = (um(i - 1, gr_k_lb_zt - 1) * um(i - 1, gr_k_lb_zt - 1));
          vm_sfc_sqd(i - 1) = (vm(i - 1, gr_k_lb_zt - 1) * vm(i - 1, gr_k_lb_zt - 1));
        });
    // Calculate surface friction velocity, u*.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_2", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          ustar(i - 1) = (Kokkos::pow(((upwp_sfc(i - 1) * upwp_sfc(i - 1)) +
                                       (vpwp_sfc(i - 1) * vpwp_sfc(i - 1))),
                                      one_fourth) > ufmin
                              ? Kokkos::pow(((upwp_sfc(i - 1) * upwp_sfc(i - 1)) +
                                             (vpwp_sfc(i - 1) * vpwp_sfc(i - 1))),
                                            one_fourth)
                              : ufmin);
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_3", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          real_t Lngth{};
          if ((Kokkos::abs(wpthlp(i - 1, gr_k_lb_zm - 1)) > eps)) {
            // Find Monin-Obukhov Length (Andre et al., 1978, p. 1866).
            Lngth = (-(((ustar(i - 1) * ustar(i - 1)) * ustar(i - 1))) /
                     (((real_t(0.35) * (one / T0)) * grav) * wpthlp(i - 1, gr_k_lb_zm - 1)));
            // Find the value of dimensionless height zeta
            // (Andre et al., 1978, p. 1866).
            zeta(i - 1) = (z_const / Lngth);
          } else {
            // wpthlp = 0
            // The value of Monin-Obukhov length is +inf when ustar < 0 and -inf
            // when ustar > 0.  Either way, zeta = 0.
            zeta(i - 1) = zero;
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
          const int i = column + 1;
          if ((zeta(i - 1) < zero)) {
            thlp2(i - 1, gr_k_lb_zm - 1) =
                (((reduce_coef * ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                                  (ustar(i - 1) * ustar(i - 1)))) *
                  four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
            rtp2(i - 1, gr_k_lb_zm - 1) =
                (((reduce_coef *
                   ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) / (ustar(i - 1) * ustar(i - 1)))) *
                  four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
            rtpthlp(i - 1, gr_k_lb_zm - 1) =
                (((reduce_coef * ((wprtp_sfc(i - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                                  (ustar(i - 1) * ustar(i - 1)))) *
                  four) *
                 Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
            wp2(i - 1, gr_k_lb_zm - 1) =
                ((ustar(i - 1) * ustar(i - 1)) *
                 (real_t(1.75) + (two * Kokkos::pow(-(zeta(i - 1)), two_thirds))));
          } else {
            thlp2(i - 1, gr_k_lb_zm - 1) =
                ((reduce_coef * four) *
                 ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                  (ustar(i - 1) * ustar(i - 1))));
            rtp2(i - 1, gr_k_lb_zm - 1) =
                ((reduce_coef * four) *
                 ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) / (ustar(i - 1) * ustar(i - 1))));
            rtpthlp(i - 1, gr_k_lb_zm - 1) =
                ((reduce_coef * four) * ((wprtp_sfc(i - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                                         (ustar(i - 1) * ustar(i - 1))));
            wp2(i - 1, gr_k_lb_zm - 1) = (real_t(1.75) * (ustar(i - 1) * ustar(i - 1)));
          }
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_5", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          thlp2(i - 1, gr_k_lb_zm - 1) =
              ((thl_tol * thl_tol) > thlp2(i - 1, gr_k_lb_zm - 1) ? (thl_tol * thl_tol)
                                                                  : thlp2(i - 1, gr_k_lb_zm - 1));
          rtp2(i - 1, gr_k_lb_zm - 1) =
              ((rt_tol * rt_tol) > rtp2(i - 1, gr_k_lb_zm - 1) ? (rt_tol * rt_tol)
                                                               : rtp2(i - 1, gr_k_lb_zm - 1));
        });
    // Calculate wstar following Andre et al., 1978, p. 1866.
    // w* = ( ( 1 / T0 ) * g * <w'thl'>|_sfc * z_i )^(1/3);
    // where z_i is the height of the mixed layer.  The value of CLUBB's
    // upward component of mixing length, Lscale_up, at the surface will be
    // used as z_i.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_6", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          real_t wstar{};
          wstar = Kokkos::pow(((((one / T0) * grav) * wpthlp(i - 1, gr_k_lb_zm - 1)) *
                               Lscale_up(i - 1, gr_k_lb_zt - 1)),
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
          if ((wpthlp(i - 1, gr_k_lb_zm - 1) > zero)) {
            usp2_sfc(i - 1) =
                ((four * (ustar(i - 1) * ustar(i - 1))) + (real_t(0.3) * (wstar * wstar)));
            vsp2_sfc(i - 1) =
                ((real_t(1.75) * (ustar(i - 1) * ustar(i - 1))) + (real_t(0.3) * (wstar * wstar)));
          } else {
            usp2_sfc(i - 1) = (four * (ustar(i - 1) * ustar(i - 1)));
            vsp2_sfc(i - 1) = (real_t(1.75) * (ustar(i - 1) * ustar(i - 1)));
          }
        });
    // Add effect of vertical compression of eddies on horizontal gustiness.
    // First, ensure that wp2 does not make the correlation
    // of (w,thl) or (w,rt) outside (-1,1).
    // Perhaps in the future we should also ensure that the correlations
    // of (w,u) and (w,v) are not outside (-1,1).
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_7", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          min_wp2_sfc_val(i - 1) =
              ((w_tol_sqd > ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                             (rtp2(i - 1, gr_k_lb_zm - 1) *
                              (max_mag_correlation_flux * max_mag_correlation_flux)))
                    ? w_tol_sqd
                    : ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                       (rtp2(i - 1, gr_k_lb_zm - 1) *
                        (max_mag_correlation_flux * max_mag_correlation_flux)))) >
                       ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                        (thlp2(i - 1, gr_k_lb_zm - 1) *
                         (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? (w_tol_sqd > ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                                   (rtp2(i - 1, gr_k_lb_zm - 1) *
                                    (max_mag_correlation_flux * max_mag_correlation_flux)))
                          ? w_tol_sqd
                          : ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                             (rtp2(i - 1, gr_k_lb_zm - 1) *
                              (max_mag_correlation_flux * max_mag_correlation_flux))))
                   : ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                      (thlp2(i - 1, gr_k_lb_zm - 1) *
                       (max_mag_correlation_flux * max_mag_correlation_flux))));
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_8", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          if (((wp2(i - 1, gr_k_lb_zm - 1) -
                ((tau_zm(i - 1, gr_k_lb_zm - 1) * lhs_splat_wp2(i - 1, gr_k_lb_zm - 1)) *
                 wp2(i - 1, gr_k_lb_zm - 1))) < min_wp2_sfc_val(i - 1))) {
            // if ( wp2(i,gr%k_lb_zm) + tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm) &
            // < min_wp2_sfc_val(i) ) then
            // splatting correction drives wp2 to overly small value
            wp2_splat_sfc_correction(i - 1) =
                (-(wp2(i - 1, gr_k_lb_zm - 1)) + min_wp2_sfc_val(i - 1));
            wp2(i - 1, gr_k_lb_zm - 1) = min_wp2_sfc_val(i - 1);
          } else {
            wp2_splat_sfc_correction(i - 1) =
                ((tau_zm(i - 1, gr_k_lb_zm - 1) * lhs_splat_wp2(i - 1, gr_k_lb_zm - 1)) *
                 wp2(i - 1, gr_k_lb_zm - 1));
            // wp2_splat_sfc_correction(i) = tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm)
            wp2(i - 1, gr_k_lb_zm - 1) =
                (wp2(i - 1, gr_k_lb_zm - 1) + wp2_splat_sfc_correction(i - 1));
          }
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_9", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          usp2_sfc(i - 1) = (usp2_sfc(i - 1) - (real_t(0.5) * wp2_splat_sfc_correction(i - 1)));
          vsp2_sfc(i - 1) = (vsp2_sfc(i - 1) - (real_t(0.5) * wp2_splat_sfc_correction(i - 1)));
        });
    // Variance of u, <u'^2>, at the surface can be found from <u_s'^2>,
    // <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
    // <u'^2>|_sfc = <u_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
    // + <v_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
    // where <u>^2 + <v>^2 /= 0.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_10", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          up2(i - 1, gr_k_lb_zm - 1) =
              ((usp2_sfc(i - 1) *
                (um_sfc_sqd(i - 1) / ((um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1)) > eps
                                          ? (um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1))
                                          : eps))) +
               (vsp2_sfc(i - 1) *
                (vm_sfc_sqd(i - 1) / ((um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1)) > eps
                                          ? (um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1))
                                          : eps))));
        });
    // Variance of v, <v'^2>, at the surface can be found from <u_s'^2>,
    // <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
    // <v'^2>|_sfc = <v_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
    // + <u_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
    // where <u>^2 + <v>^2 /= 0.
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_11", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          vp2(i - 1, gr_k_lb_zm - 1) =
              ((vsp2_sfc(i - 1) *
                (um_sfc_sqd(i - 1) / ((um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1)) > eps
                                          ? (um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1))
                                          : eps))) +
               (usp2_sfc(i - 1) *
                (vm_sfc_sqd(i - 1) / ((um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1)) > eps
                                          ? (um_sfc_sqd(i - 1) + vm_sfc_sqd(i - 1))
                                          : eps))));
        });
    // Passive scalars
    if ((sclr_dim > 0)) {
      {
        const int level_count = sclr_dim;
        Kokkos::parallel_for(
            "calc_sfc_varnce::cells_12", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
              const int i = cell % ngrdcol + 1;
              const int sclr = cell / ngrdcol + 1;
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
              if ((zeta(i - 1) < zero)) {
                sclrprtp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((reduce_coef * ((wpsclrp_sfc(i - 1, sclr - 1) * wprtp_sfc(i - 1)) /
                                      (ustar(i - 1) * ustar(i - 1)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
                sclrpthlp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((reduce_coef *
                       ((wpsclrp_sfc(i - 1, sclr - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                        (ustar(i - 1) * ustar(i - 1)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
                sclrp2(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((reduce_coef *
                       ((wpsclrp_sfc(i - 1, sclr - 1) * wpsclrp_sfc(i - 1, sclr - 1)) /
                        (ustar(i - 1) * ustar(i - 1)))) *
                      four) *
                     Kokkos::pow((one - (real_t(8.3) * zeta(i - 1))), -(two_thirds)));
              } else {
                sclrprtp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    ((reduce_coef * four) * ((wpsclrp_sfc(i - 1, sclr - 1) * wprtp_sfc(i - 1)) /
                                             (ustar(i - 1) * ustar(i - 1))));
                sclrpthlp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    ((reduce_coef * four) *
                     ((wpsclrp_sfc(i - 1, sclr - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                      (ustar(i - 1) * ustar(i - 1))));
                sclrp2(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    ((reduce_coef * four) *
                     ((wpsclrp_sfc(i - 1, sclr - 1) * wpsclrp_sfc(i - 1, sclr - 1)) /
                      (ustar(i - 1) * ustar(i - 1))));
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
          const int i = column + 1;
          real_t ustar2{};
          real_t wstar{};
          ustar2 = Kokkos::sqrt(
              ((upwp_sfc(i - 1) * upwp_sfc(i - 1)) + (vpwp_sfc(i - 1) * vpwp_sfc(i - 1))));
          // Compute wstar following Andre et al., 1976
          if ((wpthlp(i - 1, gr_k_lb_zm - 1) > zero)) {
            if (!(l_vary_convect_depth)) {
              wstar = Kokkos::pow(((((one / T0) * grav) * wpthlp(i - 1, gr_k_lb_zm - 1)) * z_const),
                                  one_third);
            } else {
              wstar = Kokkos::pow(
                  (((((one / T0) * grav) * wpthlp(i - 1, gr_k_lb_zm - 1)) * real_t(0.2)) *
                   depth_pos_wpthlp(i - 1)),
                  one_third);
            }
          } else {
            wstar = zero;
          }
          // Surface friction velocity following Andre et al. 1978
          if (!(l_vary_convect_depth)) {
            uf(i - 1) = Kokkos::sqrt((ustar2 + ((real_t(0.3) * wstar) * wstar)));
          } else {
            uf(i - 1) = Kokkos::sqrt((ustar2 + (wstar * wstar)));
          }
          uf(i - 1) = (ufmin > uf(i - 1) ? ufmin : uf(i - 1));
        });
    // Compute estimate for surface second order moments
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_14", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          wp2(i - 1, gr_k_lb_zm - 1) = (a_const(i - 1) * (uf(i - 1) * uf(i - 1)));
          up2(i - 1, gr_k_lb_zm - 1) =
              ((up2_sfc_coef(i - 1) * a_const(i - 1)) * (uf(i - 1) * uf(i - 1)));
          // From Andre, et al. 1978
          vp2(i - 1, gr_k_lb_zm - 1) =
              ((up2_sfc_coef(i - 1) * a_const(i - 1)) * (uf(i - 1) * uf(i - 1)));
          // "  "
        });
    // Notes:  1) With "a" having a value of 1.8, the surface correlations of
    // both w & rt and w & thl have a value of about 0.878.
    // 2) The surface correlation of rt & thl is 0.5.
    // Brian Griffin; February 2, 2008.
    if (!(l_vary_convect_depth)) {
      Kokkos::parallel_for(
          "calc_sfc_varnce::columns_15", ngrdcol, KOKKOS_LAMBDA(const int column) {
            const int i = column + 1;
            thlp2(i - 1, gr_k_lb_zm - 1) =
                ((real_t(0.4) * a_const(i - 1)) * ((wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1)) *
                                                   (wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1))));
            rtp2(i - 1, gr_k_lb_zm - 1) =
                ((real_t(0.4) * a_const(i - 1)) *
                 ((wprtp_sfc(i - 1) / uf(i - 1)) * (wprtp_sfc(i - 1) / uf(i - 1))));
            rtpthlp(i - 1, gr_k_lb_zm - 1) =
                (((real_t(0.2) * a_const(i - 1)) * (wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1))) *
                 (wprtp_sfc(i - 1) / uf(i - 1)));
          });
    } else {
      Kokkos::parallel_for(
          "calc_sfc_varnce::columns_16", ngrdcol, KOKKOS_LAMBDA(const int column) {
            const int i = column + 1;
            thlp2(i - 1, gr_k_lb_zm - 1) =
                (((wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1)) *
                  (wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1))) /
                 ((max_mag_correlation_flux * max_mag_correlation_flux) * a_const(i - 1)));
            rtp2(i - 1, gr_k_lb_zm - 1) =
                (((wprtp_sfc(i - 1) / uf(i - 1)) * (wprtp_sfc(i - 1) / uf(i - 1))) /
                 ((max_mag_correlation_flux * max_mag_correlation_flux) * a_const(i - 1)));
            rtpthlp(i - 1, gr_k_lb_zm - 1) =
                (max_mag_correlation_flux *
                 Kokkos::sqrt((thlp2(i - 1, gr_k_lb_zm - 1) * rtp2(i - 1, gr_k_lb_zm - 1))));
          });
    }
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_17", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          thlp2(i - 1, gr_k_lb_zm - 1) =
              ((thl_tol * thl_tol) > thlp2(i - 1, gr_k_lb_zm - 1) ? (thl_tol * thl_tol)
                                                                  : thlp2(i - 1, gr_k_lb_zm - 1));
          rtp2(i - 1, gr_k_lb_zm - 1) =
              ((rt_tol * rt_tol) > rtp2(i - 1, gr_k_lb_zm - 1) ? (rt_tol * rt_tol)
                                                               : rtp2(i - 1, gr_k_lb_zm - 1));
        });
    // Add effect of vertical compression of eddies on horizontal gustiness.
    // First, ensure that wp2 does not make the correlation
    // of (w,thl) or (w,rt) outside (-1,1).
    // Perhaps in the future we should also ensure that the correlations
    // of (w,u) and (w,v) are not outside (-1,1).
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_18", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          min_wp2_sfc_val(i - 1) =
              ((w_tol_sqd > ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                             (rtp2(i - 1, gr_k_lb_zm - 1) *
                              (max_mag_correlation_flux * max_mag_correlation_flux)))
                    ? w_tol_sqd
                    : ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                       (rtp2(i - 1, gr_k_lb_zm - 1) *
                        (max_mag_correlation_flux * max_mag_correlation_flux)))) >
                       ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                        (thlp2(i - 1, gr_k_lb_zm - 1) *
                         (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? (w_tol_sqd > ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                                   (rtp2(i - 1, gr_k_lb_zm - 1) *
                                    (max_mag_correlation_flux * max_mag_correlation_flux)))
                          ? w_tol_sqd
                          : ((wprtp_sfc(i - 1) * wprtp_sfc(i - 1)) /
                             (rtp2(i - 1, gr_k_lb_zm - 1) *
                              (max_mag_correlation_flux * max_mag_correlation_flux))))
                   : ((wpthlp(i - 1, gr_k_lb_zm - 1) * wpthlp(i - 1, gr_k_lb_zm - 1)) /
                      (thlp2(i - 1, gr_k_lb_zm - 1) *
                       (max_mag_correlation_flux * max_mag_correlation_flux))));
        });
    Kokkos::parallel_for(
        "calc_sfc_varnce::columns_19", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column + 1;
          if (((wp2(i - 1, gr_k_lb_zm - 1) -
                ((tau_zm(i - 1, gr_k_lb_zm - 1) * lhs_splat_wp2(i - 1, gr_k_lb_zm - 1)) *
                 wp2(i - 1, gr_k_lb_zm - 1))) < min_wp2_sfc_val(i - 1))) {
            // if ( wp2(i,gr%k_lb_zm) + tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm) &
            // < min_wp2_sfc_val(i) ) then
            // splatting correction drives wp2 to overly small values
            wp2_splat_sfc_correction(i - 1) =
                (-(wp2(i - 1, gr_k_lb_zm - 1)) + min_wp2_sfc_val(i - 1));
            wp2(i - 1, gr_k_lb_zm - 1) = min_wp2_sfc_val(i - 1);
          } else {
            wp2_splat_sfc_correction(i - 1) =
                ((tau_zm(i - 1, gr_k_lb_zm - 1) * lhs_splat_wp2(i - 1, gr_k_lb_zm - 1)) *
                 wp2(i - 1, gr_k_lb_zm - 1));
            // wp2_splat_sfc_correction(i) = tau_zm(i,gr%k_lb_zm) * wp2_splat(i,gr%k_lb_zm)
            wp2(i - 1, gr_k_lb_zm - 1) =
                (wp2(i - 1, gr_k_lb_zm - 1) + wp2_splat_sfc_correction(i - 1));
          }
          up2(i - 1, gr_k_lb_zm - 1) =
              (up2(i - 1, gr_k_lb_zm - 1) - (real_t(0.5) * wp2_splat_sfc_correction(i - 1)));
          vp2(i - 1, gr_k_lb_zm - 1) =
              (vp2(i - 1, gr_k_lb_zm - 1) - (real_t(0.5) * wp2_splat_sfc_correction(i - 1)));
        });
    // Passive scalars
    if ((sclr_dim > 0)) {
      {
        const int level_count = sclr_dim;
        Kokkos::parallel_for(
            "calc_sfc_varnce::cells_20", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
              const int i = cell % ngrdcol + 1;
              const int sclr = cell / ngrdcol + 1;
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
                sclrprtp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((real_t(0.4) * a_const(i - 1)) * (wprtp_sfc(i - 1) / uf(i - 1))) *
                     (wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1)));
              } else {
                sclrprtp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((real_t(0.2) * a_const(i - 1)) * (wprtp_sfc(i - 1) / uf(i - 1))) *
                     (wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1)));
              }
              if ((sclr == sclr_idx_iisclr_thl)) {
                // As above, but for thetal
                sclrpthlp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((real_t(0.4) * a_const(i - 1)) *
                      (wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1))) *
                     (wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1)));
              } else {
                sclrpthlp(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                    (((real_t(0.2) * a_const(i - 1)) *
                      (wpthlp(i - 1, gr_k_lb_zm - 1) / uf(i - 1))) *
                     (wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1)));
              }
              sclrp2(i - 1, gr_k_lb_zm - 1, sclr - 1) =
                  ((sclr_var_coef * a_const(i - 1)) * ((wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1)) *
                                                       (wpsclrp_sfc(i - 1, sclr - 1) / uf(i - 1))));
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
        const int i = column + 1;
        wp2(i - 1, gr_k_lb_zm - 1) =
            (wp2(i - 1, gr_k_lb_zm - 1) < wp2_max ? wp2(i - 1, gr_k_lb_zm - 1) : wp2_max);
      });
  Kokkos::parallel_for(
      "calc_sfc_varnce::columns_22", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column + 1;
        if ((Kokkos::abs((gr_zm(i - 1, gr_k_lb_zm - 1) - sfc_elevation(i - 1))) >
             ((Kokkos::abs((gr_zm(i - 1, gr_k_lb_zm - 1) + sfc_elevation(i - 1))) * eps) / 2))) {
          // Variances for cases where the lowest level is not at the surface.
          // Eliminate surface effects on lowest level variances.
          wp2(i - 1, gr_k_lb_zm - 1) = w_tol_sqd;
          up2(i - 1, gr_k_lb_zm - 1) = w_tol_sqd;
          vp2(i - 1, gr_k_lb_zm - 1) = w_tol_sqd;
          thlp2(i - 1, gr_k_lb_zm - 1) = (thl_tol * thl_tol);
          rtp2(i - 1, gr_k_lb_zm - 1) = (rt_tol * rt_tol);
          rtpthlp(i - 1, gr_k_lb_zm - 1) = real_t(0.0);
        }
        // gr%zm(1,gr%k_lb_zm) == sfc_elevation
      });
  if ((sclr_dim > 0)) {
    {
      const int level_count = sclr_dim;
      Kokkos::parallel_for(
          "calc_sfc_varnce::cells_23", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int sclr = cell / ngrdcol + 1;
            if ((Kokkos::abs((gr_zm(i - 1, gr_k_lb_zm - 1) - sfc_elevation(i - 1))) >
                 ((Kokkos::abs((gr_zm(i - 1, gr_k_lb_zm - 1) + sfc_elevation(i - 1))) * eps) /
                  2))) {
              // Variances for cases where the lowest level is not at the surface.
              // Eliminate surface effects on lowest level variances.
              sclrp2(i - 1, gr_k_lb_zm - 1, sclr - 1) = real_t(0.0);
              sclrprtp(i - 1, gr_k_lb_zm - 1, sclr - 1) = real_t(0.0);
              sclrpthlp(i - 1, gr_k_lb_zm - 1, sclr - 1) = real_t(0.0);
            }
          });
    }
  }
  if (stats.flags().l_sample) {
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D thlp2_sf_final("calc_sfc_varnce::thlp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::thlp2_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            thlp2_sf_final(column, level) = thlp2(column, level) / dt;
          });
      f2k::stats_finalize_budget("thlp2_sf", thlp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtp2_sf_final("calc_sfc_varnce::rtp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtp2_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            rtp2_sf_final(column, level) = rtp2(column, level) / dt;
          });
      f2k::stats_finalize_budget("rtp2_sf", rtp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D rtpthlp_sf_final("calc_sfc_varnce::rtpthlp_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::rtpthlp_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            rtpthlp_sf_final(column, level) = rtpthlp(column, level) / dt;
          });
      f2k::stats_finalize_budget("rtpthlp_sf", rtpthlp_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D up2_sf_final("calc_sfc_varnce::up2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::up2_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            up2_sf_final(column, level) = up2(column, level) / dt;
          });
      f2k::stats_finalize_budget("up2_sf", up2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D vp2_sf_final("calc_sfc_varnce::vp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::vp2_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            vp2_sf_final(column, level) = vp2(column, level) / dt;
          });
      f2k::stats_finalize_budget("vp2_sf", vp2_sf_final, stats);
    }
    {
      // Materialize the Native array/dt expression only at this sampled stats call.
      Real2D wp2_sf_final("calc_sfc_varnce::wp2_sf_final", ngrdcol, nzm);
      Kokkos::parallel_for(
          "calc_sfc_varnce::wp2_sf_final", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int level = cell / ngrdcol;
            wp2_sf_final(column, level) = wp2(column, level) / dt;
          });
      f2k::stats_finalize_budget("wp2_sf", wp2_sf_final, stats);
    }
  }
  if (f2k::clubb_at_least_debug_level_api(2)) {
    // Native diagnostics run on host. Stage each checked field once here;
    // preserve the Native per-column checker and its all-column fatal semantics.
    const auto wp2_debug_host = to_host_space(wp2);
    const auto up2_debug_host = to_host_space(up2);
    const auto vp2_debug_host = to_host_space(vp2);
    const auto thlp2_debug_host = to_host_space(thlp2);
    const auto rtp2_debug_host = to_host_space(rtp2);
    const auto rtpthlp_debug_host = to_host_space(rtpthlp);
    const auto sclrp2_debug_host = to_host_space(sclrp2);
    const auto sclrprtp_debug_host = to_host_space(sclrprtp);
    const auto sclrpthlp_debug_host = to_host_space(sclrpthlp);
    for (i = 1; i <= (ngrdcol); i += 1) {
      // TODO: Keep the public Native checker until its diagnostic owner is ported.
      f2k::sfc_varnce_check(
          sclr_dim, wp2_debug_host(i - 1, gr_k_lb_zm - 1), up2_debug_host(i - 1, gr_k_lb_zm - 1),
          vp2_debug_host(i - 1, gr_k_lb_zm - 1), thlp2_debug_host(i - 1, gr_k_lb_zm - 1),
          rtp2_debug_host(i - 1, gr_k_lb_zm - 1), rtpthlp_debug_host(i - 1, gr_k_lb_zm - 1),
          Kokkos::subview(sclrp2_debug_host, i - 1, gr_k_lb_zm - 1, Kokkos::ALL),
          Kokkos::subview(sclrprtp_debug_host, i - 1, gr_k_lb_zm - 1, Kokkos::ALL),
          Kokkos::subview(sclrpthlp_debug_host, i - 1, gr_k_lb_zm - 1, Kokkos::ALL), err_info);
      if ((to_host_space(err_info.err_code)(i - 1) == clubb_fatal_error)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          const auto header_host = to_host_space(err_info.err_header);
          diagnostic.write(header_host.data() + ((i)-1) * header_host.extent(0),
                           header_host.extent(0));
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
          const auto field_host = to_host_space(upwp_sfc);
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
          const auto field_host = to_host_space(vpwp_sfc);
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
          const auto field_host = to_host_space(wpthlp);
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
          const auto field_host = to_host_space(wprtp_sfc);
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
            const auto field_host = to_host_space(wpsclrp_sfc);
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
          const auto field_host = to_host_space(wp2);
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
          const auto field_host = to_host_space(up2);
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
          const auto field_host = to_host_space(vp2);
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
          const auto field_host = to_host_space(thlp2);
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
          const auto field_host = to_host_space(rtp2);
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
          const auto field_host = to_host_space(rtpthlp);
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
            const auto field_host = to_host_space(sclrp2);
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
            const auto field_host = to_host_space(sclrprtp);
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
            const auto field_host = to_host_space(sclrpthlp);
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
