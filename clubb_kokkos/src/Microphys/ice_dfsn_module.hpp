// Microphysics calculations translated from the same-named Native source.
//
// Public coordinators take caller-owned column-first LayoutLeft Views and retain
// Native call, statistics and recurrence order. Scalar helpers execute on the
// Kokkos backend; scratch allocations live only for the current call.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/microphys_diagnostics.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/index_mapping.hpp"
#include "CLUBB_core/T_in_K_module.hpp"
#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/fill_holes.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/matrix_solver_wrapper.hpp"
#include "CLUBB_core/diffusion.hpp"
#include "CLUBB_core/mean_adv.hpp"
#include "Microphys/parameters_microphys.hpp"
namespace clubb_kokkos {
using Int1D = Kokkos::View<int*, Kokkos::LayoutLeft>;
using Int3D = Kokkos::View<int***, Kokkos::LayoutLeft>;
inline void ice_dfsn(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                     Real2D thlm, Real2D rcm, Real2D exner, Real2D p_in_Pa, Real2D rho,
                     int saturation_formula, stats_type& stats, Real2D rcm_icedfsn,
                     Real2D thlm_icedfsn);
KOKKOS_INLINE_FUNCTION real_t Diff_denom(real_t T_in_K, real_t p_in_Pa, real_t e_i);
inline void ice_dfsn(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                     Real2D thlm, Real2D rcm, Real2D exner, Real2D p_in_Pa, Real2D rho,
                     int saturation_formula, stats_type& stats, Real2D rcm_icedfsn,
                     Real2D thlm_icedfsn) {
  // Description:
  //   This subroutine is based on a COAMPS subroutine (nov11_icedfs)
  //   written by Adam Smith and Vince Larson to calculate the
  //   depletion of cloud water by the diffusional growth of ice.
  //   The calculation is applied to all model columns.
  //
  //---------------Brian's comment--------------------------------------!
  // This code does not use actual microphysics.  Diffusional growth of !
  // ice is supposed to be the growth of ice due to diffusion of water  !
  // vapor.  Liquid water is not involved in diffusional growth.        !
  // However, in mixed phase clouds (both ice and liquid water), most   !
  // of the water vapor condenses onto the liquid droplets due to the   !
  // fact that they have so much more available surface area.  This     !
  // brings the amount of water vapor in the atmosphere to the          !
  // saturation level with respect to liquid water.  However, since the !
  // saturation vapor pressure with respect to ice is less than the     !
  // saturation vapor pressure with respect to liquid water, a          !
  // saturated atmosphere with respect to liquid water is still         !
  // supersaturated with respect to ice.  As a result, ice still grows  !
  // due to diffusion.  When this happens, the environmental vapor      !
  // pressure drops to the point of saturation with respect to ice.     !
  // This leaves the atmosphere subsaturated with respect to liquid     !
  // water.  As a result, some of the liquid water evaporates until     !
  // the atmosphere becomes saturated with respect to liquid water      !
  // again.  The process then repeats itself.  As a result, the ice     !
  // essentially grows at the expense of the liquid water.  This is     !
  // why the diffusional growth of ice is being deducted from liquid    !
  // water in this subroutine.
  //-------------------------------------------------------------------------------
  // References:
  //   Section 4.2 of Larson et al. (2006), "What determines altocumulus
  //     dissipation time?", J. Geophys. Res., Vol. 111, D19207.
  //
  //   Mitchell, D. L. (1996), "Use of mass- and area- ...", J. Atmos. Sci.
  //     Vol. 53, 1710--1723.
  //
  //   Rogers and Yau (1989), "A Short Course in Cloud Physics", 3rd. Ed.
  //
  //   Fleishauer et al. (2002), "Observed microphysical structure of
  //     midlevel, mixed-phase clouds", J. Atmos. Sci., Vol. 59,
  //     pp. 1779--1804.
  //-------------------------------------------------------------------------------
  // Constant Parameters
  // Number of ice crystals per unit volume of air    [m^{-3}]
  // Vince Larson avgd legs 2 and 7 (Fleishauer et al)  21 Jan 2005
  //---------------------- Input variables ----------------------
  //------------------------ Input/Output Variables ------------------------
  //--------------------------- Output Variables ---------------------------
  //--------------------------- Local Variables ---------------------------
  // diam(:) ! Icedfs diameter; Michael Falk, 1 Nov 2006
  // dqc_dt_icedfsn(:) ! Icedfs change in liquid; Michael Falk, 1 Nov 2006
  // u_T_cm(:)        ! Icedfs fallspeed (cm/s); Michael Falk, 1 Nov 2006
  auto workspace_scope = workspace.scope();
  const bool l_sample = stats.l_sample;

  constexpr real_t N_i = real_t(2000.);
  Real2D T_in_K        = workspace_scope.view<Real2D>(
      "ice_dfsn::T_in_K", ngrdcol, gr.nzt); // Absolute temperature                        [K]
  Real2D mass_ice_cryst =
      workspace_scope.view<Real2D>("ice_dfsn::mass_ice_cryst", ngrdcol,
                                   gr.nzt); // Mass of a single ice crystal                [kg]
  Real2D r_s = workspace_scope.view<Real2D>(
      "ice_dfsn::r_s", ngrdcol, gr.nzt); // Saturation mixing ratio over vapor          [kg kg^{-1}]
  Real2D e_s = workspace_scope.view<Real2D>(
      "ice_dfsn::e_s", ngrdcol, gr.nzt); // Saturation vapor pressure over liquid       [Pa]
  Real2D e_i = workspace_scope.view<Real2D>(
      "ice_dfsn::e_i", ngrdcol, gr.nzt); // Saturation vapor pressure over ice          [Pa]
  Real2D S_i = workspace_scope.view<Real2D>(
      "ice_dfsn::S_i", ngrdcol, gr.nzt); // Ratio of saturation w.r.t. liquid to that for ice []
  Real2D Denom = workspace_scope.view<Real2D>(
      "ice_dfsn::Denom", ngrdcol,
      gr.nzt); // Denominator of diffusional growth equation  [m s kg^{-1}]
  Real2D dmass_ice_cryst = workspace_scope.view<Real2D>(
      "ice_dfsn::dmass_ice_cryst", ngrdcol,
      gr.nzt); // Change in ice mass over vertical grid box   [kg m^{-1}]
  Real2D diam = workspace_scope.view<Real2D>(
      "ice_dfsn::diam", ngrdcol, gr.nzt); // Diameter of ice crystal                     [m]
  Real2D u_T_cm = workspace_scope.view<Real2D>(
      "ice_dfsn::u_T_cm", ngrdcol,
      gr.nzt);           // Fallspeed of ice crystal in cm/s            [cm s^{-1}]
  real_t a_coef;         // Pre-factor for mass-diameter relationship, Mitchell (1996) [kg]
  real_t b_expn;         // Exponential for mass-diameter relationship, Mitchell (1996) []
  real_t k_u_coef;       // Pre-factor for fallspeed-diameter formula                  [m s^{-1}]
  real_t q_expn;         // Exponential of density in fallspeed-diameter formula       []
  real_t n_expn;         // Exponential of diameter in fallspeed-diameter formula      []
  const auto gr_nzt       = gr.nzt;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  //---------------------- Begin Code ----------------------

  // Determine absolute temperature
  thlm2T_in_K_api(gr.nzt, ngrdcol, thlm, exner, rcm, T_in_K);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //                                                                     !
  // Coefficients for mass-diameter relationship, Mitchell (1996)        !
  // mass = a (diameter/(1 meter))^b,  [a] = kg, [b] = []                !
  //                                                                     !
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  a_coef = real_t(2.05e-3);
  b_expn = real_t(1.8);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //                                                                     !
  //  Coefficients for mass-diameter relationship, Kajikawa (1989)       !
  //  mass = a (diam/(1m))^b,  [a] = kg, [b] = []                        !
  //                                                                     !
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //       a = 2.50e-4_core_rknd
  //       b = 1.4_core_rknd

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //                                                                     !
  //  Coefficients for fallspeed-diameter relationship, Mitchell (1996)  !
  //  u_T = k_u rho^{-q} (diameter/(1 meter))^n,                         !
  //       [k_u] = m/s, [q] = [], [n] = [], [rho] = kg m^{-3}            !
  //                                                                     !
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  k_u_coef = real_t(55.);
  q_expn   = real_t(0.17);
  n_expn   = real_t(0.70);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //                                                                     !
  //  Coefficients for fallspeed-diameter relationship, Kajikawa (1989)  !
  //  u_T = k_u rho^{-q} (diam/(1m))^n,  [k_u] = m/s, [q] = [], [n] = [] !
  //       [rho] = kg m^{-3}                                             !
  //                                                                     !
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //       k_u = 0.438_core_rknd
  //       q = 0.0_core_rknd
  //       n = 0.0742_core_rknd

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //                                                                     !
  //  Initialize ice particle mass                                       !
  //                                                                     !
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Kokkos::parallel_for(
      "ice_dfsn::loop1", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Vertical-level loop index [-]
        const int i = column;
        for (k = 0; k < gr_nzt; k += 1) {
          mass_ice_cryst((i), (k)) = real_t(1.0e-11);
        }
      });

  Kokkos::parallel_for(
      "ice_dfsn::loop2", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Vertical-level loop index [-]
        const int i = column;
        for (k = gr_nzt - 1; k >= 0; k += (-1)) {

          // Check whether we're in cloud and below freezing.
          // Note:  A value of 1.0E-5 kg/kg is used as a threshold value
          // for rcm because the CLUBB model shows a small amount of liquid
          // water all the way to the model top, which messes with the
          // ice diffusion calculations.
          if (((rcm((i), (k)) >= real_t(1.0E-5)) && (T_in_K((i), (k)) < T_freeze_K))) {

            // Find saturation mixing ratio over vapor [kg kg^{-1}]
            r_s((i), (k)) =
                sat_mixrat_liq_api(p_in_Pa((i), (k)), T_in_K((i), (k)), saturation_formula);

            // Saturation vapor pressure over liquid in Pa
            e_s((i), (k)) = ((r_s((i), (k)) * p_in_Pa((i), (k))) / (ep + r_s((i), (k))));

            // Saturation vapor pressure over ice in Pa, Eq. 2.15 Rogers and Yau
            e_i((i), (k)) =
                (e_s((i), (k)) / Kokkos::exp(((Lf / (Rv * T_freeze_K)) *
                                              ((T_freeze_K / T_in_K((i), (k))) - real_t(1.0)))));

            // Saturation ratio in a liquid-saturated cloud, p. 158 Rogers and Yau
            //---------------Brian's comment--------------------------------------!
            // The actual formula is:  Si = e/ei = (e/esat)*(esat/ei) = S*(esat/ei)     !
            // It is assumed that any supersaturation forms liquid water and that !
            // the atmosphere is then saturated with respect to liquid water.     !
            // Therefore, S = 1.0, allowing Si = esat/ei.                           !
            //--------------------------------------------------------------------!
            S_i((i), (k)) = (e_s((i), (k)) / e_i((i), (k)));

            // Denominator of diffusional growth equation, 9.4 of Rogers and Yau
            Denom((i), (k)) = Diff_denom(T_in_K((i), (k)), p_in_Pa((i), (k)), e_i((i), (k)));

            // Change in mass of a single ice crystal, m,
            // as it falls a distance gr%invrs_dzt(1,:) in meters

            //---------------Brian's comment--------------------------------------!
            // dm/dt = 4*pi*C*(Si-1)/Denom; Rogers and Yau, Eq. 9.4.              !
            // For plate-type ice crystals, C = 2r/pi (Rogers and Yau, p. 159).   !
            // Since 2r = D, C = D/pi, and the equation becomes:                  !
            // dm/dt = 4*D*(Si-1)/Denom.                                          !
            // The mass-diameter relationship for an ice crystal is:              !
            // D = (m/a)^(1/b); Rogers and Yau, Eq. 9.7.  This means:             !
            // dm/dt = [(m/a)^(1/b)]*4*(Si-1)/Denom;                              !
            // Dividing by rho yields the change in mixing ratio over time        !
            // for an individual crystal.  Multiplying that by the ice crystal    !
            // concentration yields the overall change in mixing ratio over time. !
            //--------------------------------------------------------------------!
            rcm_icedfsn((i), (k)) =
                (((-(N_i / rho((i), (k)))) *
                  ((real_t(4.) * (S_i((i), (k)) - real_t(1.))) / Denom((i), (k)))) *
                 Kokkos::pow((mass_ice_cryst((i), (k)) / a_coef), (real_t(1.) / b_expn)));

            // Ensure that liquid is not over-depleted
            if (((rcm((i), (k)) + (rcm_icedfsn((i), (k)) * dt)) < real_t(0.0))) {
              rcm_icedfsn((i), (k)) = ((-rcm((i), (k))) / dt);
            }

            //---------------Brian's comment-----------------------------------!
            // dm = (dm/dt)*(dt/dz)*dz                                         !
            // dm = (dm/dt)*(1/u_T)*dz                                         !
            //-----------------------------------------------------------------!
            dmass_ice_cryst((i), (k)) =
                ((((((real_t(4.) * (S_i((i), (k)) - real_t(1.))) / Denom((i), (k))) *
                    (real_t(1.0) / k_u_coef)) *
                   Kokkos::pow(rho((i), (k)), q_expn)) *
                  Kokkos::pow((mass_ice_cryst((i), (k)) / a_coef),
                              ((real_t(1.0) - n_expn) / b_expn))) *
                 (real_t(1.0) / gr_invrs_dzm((i), (k - 1))));

            if ((k > 0)) {
              mass_ice_cryst((i), (k - 1)) = (mass_ice_cryst((i), (k)) + dmass_ice_cryst((i), (k)));
            }

            // Diameter of ice crystal in meters.
            diam((i), (k)) =
                Kokkos::pow((mass_ice_cryst((i), (k)) / a_coef), (real_t(1.) / b_expn));

            // Fallspeed of ice crystal in cm/s.
            u_T_cm((i), (k)) =
                (((cm_per_m * k_u_coef) *
                  Kokkos::pow((mass_ice_cryst((i), (k)) / a_coef), (n_expn / b_expn))) *
                 Kokkos::pow(rho((i), (k)), (-q_expn)));

          } else {
            // There's no liquid and/or ice present; assume no ice growth

            if ((k > 0)) {
              mass_ice_cryst((i), (k - 1)) = mass_ice_cryst((i), (k));
            }
            rcm_icedfsn((i), (k)) = real_t(0.0);
            // Set zero to remind that we don't grow ice
            diam((i), (k)) = real_t(0.0);
            // Set zero to remind that we don't grow ice
            u_T_cm((i), (k)) = real_t(0.0);
          }

          // k = gr%nzt, 1, -1
        }
      });

  if (stats.l_sample) {
    stats_update(std::string("rcm_icedfs"), rcm_icedfsn, stats);
    stats_update(std::string("diam"), diam, stats);
    stats_update(std::string("mass_ice_cryst"), mass_ice_cryst, stats);
    stats_update(std::string("u_T_cm"), u_T_cm, stats);
  }

  // Determine time tendency of liquid potential temperature
  Kokkos::parallel_for(
      "ice_dfsn::loop3", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Vertical-level loop index [-]
        const int i = column;
        for (k = 0; k < gr_nzt; k += 1) {
          thlm_icedfsn((i), (k)) = ((-(Lv / (Cp * exner((i), (k))))) * rcm_icedfsn((i), (k)));
        }
      });

  return;
}
KOKKOS_INLINE_FUNCTION real_t Diff_denom(real_t T_in_K, real_t p_in_Pa, real_t e_i) {
  // Description:
  //   Compute denominator of diffusional growth equation
  //
  // References:
  //   Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
  //
  //-----------------------------------------------------------------------------
  // Reference:  Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
  // Constant Parameters
  //   real, parameter :: Ls = 2.834e6

  real_t Diff_denom; // Denominator of diffusional growth equation  [m s kg^{-1}]
  real_t Ka;         // Ka in J./(m.*sec.*K)
  real_t Dv;
  real_t Fk;
  real_t Fd;
  real_t Celsius;

  // ---- Begin Code ----

  Celsius = (T_in_K - T_freeze_K);

  // Ka in cal./(cm.*sec.*C)
  Ka = ((real_t(5.69) + (real_t(0.017) * Celsius)) * real_t(0.00001));
  // Ka in J./(m.*sec.*K)
  Ka = ((real_t(4.1868) * real_t(100.0)) * Ka);

  Dv = ((real_t(0.221) * Kokkos::pow((T_in_K / T_freeze_K), real_t(1.94))) *
        (real_t(101325.0) / p_in_Pa));
  // Dv in (cm.^2)/sec.  ! .221 is correct.
  // Dv in (m.^2)/sec.
  Dv = (Dv / real_t(10000.0));

  Fk = ((((Ls / (Rv * T_in_K)) - real_t(1.0)) * Ls) / (Ka * T_in_K));
  Fd = ((Rv * T_in_K) / (Dv * e_i));

  Diff_denom = (Fk + Fd);

  return Diff_denom;
}
// Host/strided callers stage through the existing array owner; device Views alias.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
inline void ice_dfsn(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                     Array0 thlm, Array1 rcm, Array2 exner, Array3 p_in_Pa, Array4 rho,
                     int saturation_formula, stats_type& stats, Array5 rcm_icedfsn,
                     Array6 thlm_icedfsn) {
  // Native argument order is unchanged; only memory-space staging occurs here.
  auto workspace_scope        = workspace.scope();
  auto thlm_execution         = to_execution_space(workspace, thlm, true);
  auto rcm_execution          = to_execution_space(workspace, rcm, true);
  auto exner_execution        = to_execution_space(workspace, exner, true);
  auto p_in_Pa_execution      = to_execution_space(workspace, p_in_Pa, true);
  auto rho_execution          = to_execution_space(workspace, rho, true);
  auto rcm_icedfsn_execution  = to_execution_space(workspace, rcm_icedfsn, false);
  auto thlm_icedfsn_execution = to_execution_space(workspace, thlm_icedfsn, false);
  ice_dfsn(workspace, gr, ngrdcol, dt, thlm_execution, rcm_execution, exner_execution,
           p_in_Pa_execution, rho_execution, saturation_formula, stats, rcm_icedfsn_execution,
           thlm_icedfsn_execution);
  from_execution_space(workspace, rcm_icedfsn, rcm_icedfsn_execution);
  from_execution_space(workspace, thlm_icedfsn, thlm_icedfsn_execution);
}
} // namespace clubb_kokkos
