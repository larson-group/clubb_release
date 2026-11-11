// Compute surface exchange for the GABLS3 nighttime LES case.
//
// gabls3_night_sfclyr interpolates theta_l [K] and water [kg/kg] at roughness
// height, evaluates the SAM land-flux formulation, and either interpolates or
// diagnoses momentum fluxes. Inputs include time [s], winds [m/s], theta_l [K],
// water [kg/kg] and lowest height [m]. Outputs are momentum[(m/s)^2], heat [K m/s],
// water[(kg/kg) m/s] and ustar [m/s], one value per column. The appended case
// context replaces Native module storage. If l_t_dependent is false, outputs
// remain unchanged. Private similarity functions retain Native order/formulas.
// References: http://www4.ncsu.edu/~sbasu5/GABLS3/; SAM 6.7.5 landflx;
// Joshua Fasching's January 2009 GABLS3 night heat/moisture modifications.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

KOKKOS_INLINE_FUNCTION static void landflx(real_t th, real_t ts, real_t qh, real_t qs,
                                           real_t uh, real_t vh, real_t h, real_t z0,
                                           real_t& shf, real_t& lhf, real_t& vel, real_t& ustar);

//-----------------------------------------------------------------------
inline void gabls3_night_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time,
                                Real1D um_sfc, Real1D vm_sfc, Real1D thlm_sfc, Real1D rtm_sfc,
                                Real1D lowest_level, Real1D upwp_sfc, Real1D vpwp_sfc,
                                Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar,
                                const time_dependent_input_state& input) {
  // Description:
  // Compute heat, moisture and momentum exchange for GABLS3 nighttime using
  // supplied roughness-height conditions and optional supplied momentum fluxes.
  auto workspace_scope = workspace.scope();
  constexpr real_t z0 = real_t(0.15); // Roughness length [m].
  if (input.l_t_dependent) {
    real_t time_frac;
    int before_time, after_time;
    // Identify the neighboring times and their interpolation fraction.
    time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
                input.time_sfc_given, before_time, after_time, time_frac);
    const auto surface_theta = to_host_space(workspace, input.thlm_sfc_given);
    const real_t ts =
        linear_interp_factor(time_frac, surface_theta(after_time), surface_theta(before_time));
    const auto surface_water = to_host_space(workspace, input.rtm_sfc_given);
    const real_t qs =
        linear_interp_factor(time_frac, surface_water(after_time), surface_water(before_time));

    const bool l_input_xpwp_sfc = input.l_input_xpwp_sfc;
    // Compute heat and moisture fluxes.
    using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
    Kokkos::parallel_for(
        "gabls3_night_sfclyr::land", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          real_t ubar;
          landflx(thlm_sfc(i), ts, rtm_sfc(i), qs, um_sfc(i), vm_sfc(i), lowest_level(i), z0,
                  wpthlp_sfc(i), wprtp_sfc(i), ubar, ustar(i));
          if (!l_input_xpwp_sfc) {
            upwp_sfc(i) = -um_sfc(i) * (ustar(i) * ustar(i)) / ubar;
            vpwp_sfc(i) = -vm_sfc(i) * (ustar(i) * ustar(i)) / ubar;
          }
        });
    if (input.l_input_xpwp_sfc) {
      // Feed in the supplied momentum fluxes.
      const auto u_flux = to_host_space(workspace, input.upwp_sfc_given);
      const real_t upwp_sfc_interp =
          linear_interp_factor(time_frac, u_flux(after_time), u_flux(before_time));
      const auto v_flux = to_host_space(workspace, input.vpwp_sfc_given);
      const real_t vpwp_sfc_interp =
          linear_interp_factor(time_frac, v_flux(after_time), v_flux(before_time));
      Kokkos::parallel_for(
          "gabls3_night_sfclyr::input_momentum", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
            upwp_sfc(i) = upwp_sfc_interp;
            vpwp_sfc(i) = vpwp_sfc_interp;
          });
    }
  }
}

//-----------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static real_t psi_h(real_t x, real_t xlmo) {
  // Description:
  // Return the heat similarity correction at height x [m] for Obukhov length
  // xlmo [m], as used in the GABLS3 night flux denominator.
  return (-real_t(5.0) * x) / xlmo;
}

//-----------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static real_t gm1(real_t x) {
  // Description:
  // Evaluate the unstable momentum similarity function at dimensionless x.
  return Kokkos::pow(real_t(1.0) - real_t(15.0) * x, real_t(0.25));
}

//-----------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static real_t gh1(real_t x) {
  // Description:
  // Evaluate the unstable heat similarity function at dimensionless x.
  return Kokkos::sqrt(real_t(1.0) - real_t(9.0) * x) / real_t(0.74);
}

//-----------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static real_t fm1(real_t x) {
  // Description:
  // Integrate the momentum similarity function using the Native log/atan form.
  const real_t pii = Kokkos::acos(-real_t(1.0)) / real_t(2.0);
  return real_t(2.0) * Kokkos::log((real_t(1.0) + x) / real_t(2.0)) +
         Kokkos::log((real_t(1.0) + x * x) / real_t(2.0)) - real_t(2.0) * Kokkos::atan(x) + pii;
}

//-----------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static real_t fh1(real_t x) {
  // Description:
  // Integrate the heat similarity function at the supplied dimensionless x.
  return real_t(2.0) * Kokkos::log((real_t(1.0) + real_t(0.74) * x) / real_t(2.0));
}

//------------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION static void landflx(real_t th, real_t ts, real_t qh, real_t qs,
                                           real_t uh, real_t vh, real_t h, real_t z0,
                                           real_t& shf, real_t& lhf, real_t& vel, real_t& ustar) {
  // Description:
  // SAM 6.7.5 landflx: derive kinematic sensible/latent heat fluxes from theta [K],
  // vapor [kg/kg], wind [m/s] at height h [m], and conditions at roughness z0[m].
  // Return shf [K m/s], lhf[(kg/kg) m/s], wind speed and friction velocity [m/s].
  const real_t zody = Kokkos::log(h / z0);
  vel               = Kokkos::sqrt(Kokkos::max(real_t(0.5), uh * uh + vh * vh));
  const real_t r    = real_t(9.81) / ts *
                      (th * (real_t(1.0) + ep1 * qh) - ts * (real_t(1.0) + ep1 * qs)) * h /
                      (vel * vel);
  int iter          = 0;
  real_t xsi, xsi1, xsi2, xm, xh, fm, fh;
  if (r < real_t(0.0)) {
    xsi  = real_t(0.0);
    iter = iter + 1;
    xm   = gm1(xsi);
    xh   = gh1(xsi);
    fm   = zody - fm1(xm);
    fh   = real_t(0.74) * (zody - fh1(xh));
    xsi1 = r / fh * (fm * fm);
    xsi  = xsi1;

    xsi  = -Kokkos::abs(xsi);
    iter = iter + 1;
    xm   = gm1(xsi);
    xh   = gh1(xsi);
    fm   = zody - fm1(xm);
    fh   = real_t(0.74) * (zody - fh1(xh));
    xsi1 = r / fh * (fm * fm);
    xsi  = xsi1;

    xsi  = -Kokkos::abs(xsi);
    iter = iter + 1;
    xm   = gm1(xsi);
    xh   = gh1(xsi);
    fm   = zody - fm1(xm);
    fh   = real_t(0.74) * (zody - fh1(xh));
    xsi1 = r / fh * (fm * fm);
    xsi  = xsi1;
  } else {
    const real_t a = real_t(4.8) * real_t(4.8) * r - real_t(1.00) * real_t(6.35);
    const real_t b = (real_t(2.0) * r * real_t(4.8) - real_t(1.00)) * zody;
    const real_t c = r * (zody * zody);
    const real_t d = Kokkos::sqrt(b * b - real_t(4.0) * a * c);
    xsi1           = (-b + d) / a / real_t(2.0);
    xsi2           = (-b - d) / a / real_t(2.0);
    xsi            = Kokkos::max(xsi1, xsi2);
    fm             = zody + real_t(4.8) * xsi;
    fh             = real_t(1.00) * (zody + real_t(7.8) * xsi);
  }
  (void)iter;

  vel = Kokkos::sqrt(uh * uh + vh * vh);
  // GABLS3 night modification: use the case specification for both fluxes,
  // rather than the disabled SAM 0.4^2/fm/fh formulation (J. Fasching, Jan 2009).
  ustar = real_t(0.4) / fm * vel;
  if (xsi >= real_t(0.0))
    xsi = Kokkos::max(real_t(1.e-5), xsi);
  else
    xsi = Kokkos::min(-real_t(1.e-5), xsi);
  const real_t xlmo = h / xsi;
  shf               = (real_t(0.4) * ustar * (ts - th)) /
                      (Kokkos::log(h / real_t(0.25)) - psi_h(h, xlmo) + psi_h(real_t(0.25), xlmo));
  lhf               = (real_t(0.4) * ustar * (qs - qh)) /
                      (Kokkos::log(h / real_t(0.25)) - psi_h(h, xlmo) + psi_h(real_t(0.25), xlmo));
}

} // namespace clubb_kokkos
