// Sounding coordinate/thermodynamic interpretation; see input_interpret.hpp.
#include "scratch_workspace.hpp"
#include "input_interpret.hpp"
#include "input_names.hpp"
#include "hydrostatic_module.hpp"
#include "CLUBB_core/calc_pressure.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/saturation.hpp"
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

//=============================================================================
void read_z_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                    const std::vector<one_dim_read_var>& retVars, real_t p_sfc,
                    real_t zm_init,                                       // In
                    int saturation_formula,                               // In
                    Real1D z, Real1D p_in_Pa, std::string& altitude_type) // Out
{
  // Description:
  // Searches for the variable specified by either 'z[m]' or 'Press[Pa]' in the
  // collection of retVars.  If the subroutine finds the variable indicated by
  // 'z[m]', then it returns it.  If the subroutine finds 'Press[Pa]' then it
  // converts it to values of altitude in meters.  If it does not find either
  // or finds both the program using this subroutine will exit gracefully with
  // a warning message.

  //-----------------------------------------------------------------------

  // p_sfc[Pa],zm_init[m]; z[m],p_in_Pa[Pa]. Scratch is(column1,sounding level).
  auto workspace_scope  = workspace.scope();
  const bool l_z        = get_target_index(nvar, z_name, retVars) >= 0;
  const bool l_pressure = get_target_index(nvar, pressure_name, retVars) >= 0;
  // Native leaves outputs undefined if both coordinate headers exist. Retain
  // this exact source guard; do not invent precedence or a new input repair.
  if (int(l_z) + int(l_pressure) <= 1) {
    if (l_z) {
      // The input sounding is given in terms of altitude.
      altitude_type = z_name;
      // Obtain the value of altitude at each sounding level.
      read_x_profile(nvar, nsize, altitude_type, retVars, z);
      // Set pressure at the sounding levels to the fill value.
      Kokkos::deep_copy(p_in_Pa, real_t(-999.9));
    } else if (l_pressure) {
      // The input sounding is given in terms of pressure.
      altitude_type = pressure_name;
      // Obtain total pressure at each pressure sounding level.
      read_x_profile(nvar, nsize, altitude_type, retVars, p_in_Pa);
      const int nlevels = retVars[0].values.extent(0);
      const real_t zero_threshold = clubb_kokkos::zero_threshold;
      Real2D exner                = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                 "read_z_profile::exner", 1, nsize);
      Real2D thvm                 = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                                 "read_z_profile::thvm", 1, nsize);
      Real2D rcm = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "read_z_profile::rcm",
                                                1, nsize);
      Real2D theta = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                  "read_z_profile::theta", 1, nsize);
      Real2D rtm  = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "read_z_profile::rtm",
                                                 1, nsize);
      Real2D thlm = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                 "read_z_profile::thlm", 1, nsize);
      std::string temperature_type;
      // Obtain theta,theta_l or temperature at each pressure sounding level.
      auto theta_profile = Kokkos::subview(theta, 0, Kokkos::ALL);
      read_theta_profile(workspace, nvar, nsize, retVars, temperature_type, theta_profile);
      // Obtain total water mixing ratio at each pressure sounding level.
      auto rtm_profile = Kokkos::subview(rtm, 0, Kokkos::ALL);
      read_x_profile(nvar, nsize, rt_name, retVars, rtm_profile);
      // Calculate Exner from pressure.
      Kokkos::parallel_for(
          "read_z_profile::exner", nlevels,
          KOKKOS_LAMBDA(int k) { exner(0, k) = Kokkos::pow(p_in_Pa(k) / p0, kappa); });
      if (temperature_type == temperature_name) {
        // theta actually contains temperature[K] here. Determine initial cloud
        // water from temperature/rtm, then convert temperature to theta and
        // theta_l=theta-[Lv/(Cp*Exner)]*rcm.
        Kokkos::parallel_for(
            "read_z_profile::temperature", nlevels, KOKKOS_LAMBDA(int k) {
              rcm(0, k) = Kokkos::max(
                  rtm(0, k) - sat_mixrat_liq_api(p_in_Pa(k), theta(0, k), saturation_formula),
                  zero_threshold);
              theta(0, k) = theta(0, k) / exner(0, k);
              thlm(0, k) = theta(0, k) - Lv / (Cp * exner(0, k)) * rcm(0, k);
            });
      } else if (temperature_type == theta_name) {
        // theta does contain potential temperature. Determine initial cloud
        // water from theta,Exner and rtm, then recover theta_l.
        Kokkos::parallel_for(
            "read_z_profile::theta", nlevels, KOKKOS_LAMBDA(int k) {
              rcm(0, k) =
                  Kokkos::max(rtm(0, k) - sat_mixrat_liq_api(p_in_Pa(k), theta(0, k) * exner(0, k),
                                                             saturation_formula),
                              zero_threshold);
              thlm(0, k) = theta(0, k) - Lv / (Cp * exner(0, k)) * rcm(0, k);
            });
      } else if (temperature_type == thetal_name) {
        // theta contains theta_l. Unsaturated theta=theta_l; saturated profiles
        // require the Native iterative cloudwater adjustment. Recover theta
        // afterward:theta=theta_l+[Lv/(Cp*Exner)]*rcm.
        Kokkos::parallel_for(
            "read_z_profile::thetal", nlevels, KOKKOS_LAMBDA(int k) {
              thlm(0, k) = theta(0, k);
              rcm(0, k) =
                  rcm_sat_adj(thlm(0, k), rtm(0, k), p_in_Pa(k), exner(0, k), saturation_formula);
              theta(0, k) = thlm(0, k) + Lv / (Cp * exner(0, k)) * rcm(0, k);
            });
      } else {
        std::cerr << "Invalid temperature_type: " << temperature_type << '\n';
        throw std::runtime_error("");
      }
      // Now, the variable "theta" contains potential temperature, theta.
      // Compute initial theta_v based on theta, total water mixing ratio, and
      // cloud water mixing ratio.
      //
      // theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
      //
      // which can be rearranged as:
      //
      // theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
      //
      // The exact form of the equation for theta_v (including water vapor and
      // cloud water) is:
      //
      // theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
      //
      // which can be rearranged as:
      //
      // theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
      //                         * { r_v / ( 1 + r_v + r_c ) }
      //                     - { r_c / ( 1 + r_v + r_c ) } ].
      //
      // This version is written with r_t = r_v + r_c, such that:
      //
      // theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
      //                         * { ( r_t - r_c ) / ( 1 + r_t ) }
      //                     - { r_c / ( 1 + r_t ) } ].
      //
      // To use theta_l instead of theta, simply substitute the following for
      // theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
      //
      // The CLUBB code uses a linearized version of the above equation (in
      // order to calculate thv' terms -- such as w'thv', etc.) for theta_v
      // throughout the model code, such that:
      //
      // theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
      //                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
      //
      // where thv_ds is used as a reference value to approximate theta_l.
      // The source's array expression materializes one reference-temperature
      // profile; no scientific expression is hidden behind an adapter.
      Real2D thv_ds_zt = workspace_scope.view<Real2D>("read_z_profile::thv_ds_zt", 1, nlevels);
      Kokkos::parallel_for(
          "read_z_profile::reference", nlevels, KOKKOS_LAMBDA(int k) {
            thv_ds_zt(0, k) = theta(0, k) * Kokkos::pow(one + ep2 * (rtm(0, k) - rcm(0, k)), kappa);
          });
      calculate_thvm(nlevels, 1, thlm, rtm, rcm, exner, // In
                     thv_ds_zt,                         // In
                     thvm);                             // Out
      // Find pressure-sounding altitudes. The Native actuals are leading
      // contiguous prefixes despite their declared rank2 scratch storage.
      auto thvm_prefix = Kokkos::subview(thvm, 0, Kokkos::pair<int, int>(0, nlevels));
      auto exner_prefix = Kokkos::subview(exner, 0, Kokkos::pair<int, int>(0, nlevels));
      auto z_prefix = Kokkos::subview(z, Kokkos::pair<int, int>(0, nlevels));
      inverse_hydrostatic(workspace, p_sfc, zm_init, nlevels, thvm_prefix, exner_prefix, // In
                          z_prefix);                                                     // Out
    } else {
      throw std::runtime_error("Could not read sounding vertical-coordinate variable");
    }
  }
}

//=============================================================================
void read_theta_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                        const std::vector<one_dim_read_var>& retVars, // In
                        std::string& temperature_type, Real1D theta)  // Out
{
  (void)workspace;
  //
  //  Description: Searches for the variable specified by either 'thetal[K]' or 'theta[K]' in the
  //  collection of retVars. If the function finds the variable then it returns
  //  it. If it does not the program using this function will exit gracefully
  //  with a warning message.
  //
  //-----------------------------------------------------------------------------------------------

  const bool l_theta       = get_target_index(nvar, theta_name, retVars) >= 0;
  const bool l_thetal      = get_target_index(nvar, thetal_name, retVars) >= 0;
  const bool l_temperature = get_target_index(nvar, temperature_name, retVars) >= 0;
  // Retain Native's ambiguous-header no-write boundary, as in read_z_profile.
  if (int(l_theta) + int(l_thetal) + int(l_temperature) <= 1) {
    if (l_theta)
      temperature_type = theta_name;
    else if (l_thetal)
      temperature_type = thetal_name;
    else if (l_temperature)
      temperature_type = temperature_name;
    else
      throw std::runtime_error("Could not read theta compatable variable");
    read_x_profile(nvar, nsize, temperature_type, retVars, theta);
  }
}

//=============================================================================
void read_subs_profile(ScratchWorkspace& workspace, int nvar, int nsize,
                       const std::vector<one_dim_read_var>& retVars, // In
                       std::string& subs_type, Real1D subs)          // Out
{
  (void)workspace;
  //
  //  Description: Searches for the variable specified by either 'w[m\s]' or 'omega[Pa\s]' in the
  //  collection of retVars. If the function finds the variable then it returns
  //  it. If it does not the program using this function will exit gracefully
  //  with a warning message.
  //
  //-----------------------------------------------------------------------------------------------

  const bool l_w     = get_target_index(nvar, wm_name, retVars) >= 0;
  const bool l_omega = get_target_index(nvar, omega_name, retVars) >= 0;
  // Retain Native's ambiguous-header no-write boundary.
  if (int(l_w) + int(l_omega) <= 1) {
    if (l_w)
      subs_type = wm_name;
    else if (l_omega)
      subs_type = omega_name;
    else
      throw std::runtime_error("Could not read vertical velocity compatable variable");
    read_x_profile(nvar, nsize, subs_type, retVars, subs);
  }
}

} // namespace clubb_kokkos
