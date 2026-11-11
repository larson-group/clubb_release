// Read and apply case-owned time-dependent input; see time_dependent_input.hpp.
#include "scratch_workspace.hpp"
#include "time_dependent_input.hpp"
#include "sfc_flux.hpp"
#include "Input_fields/input_names.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/remapping_module.hpp"
#include "readiopdata_module.hpp"
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

static void initialize_t_dependent_sfc(ScratchWorkspace& workspace, int iunit, bool l_readiopdata,
                                       Real2D rho, const std::string& input_file,
                                       time_dependent_input_state& input, readiopdata_state& iop);
static void initialize_t_dependent_forcings(ScratchWorkspace& workspace, int iunit,
                                            const std::string& input_file, int nz, Real1D grid,
                                            Real1D p_in_Pa, bool l_add_dycore_grid,
                                            int grid_adapt_in_time_method,
                                            time_dependent_input_state& input);
static void finalize_t_dependent_forcings(time_dependent_input_state& input);
static void finalize_t_dependent_sfc(time_dependent_input_state& input);
template <class Grid>
static void read_to_grid(ScratchWorkspace& workspace, int ntwo_dim_vars, int dim_size,
                         int other_dim_size, int nz, Grid grid,
                         const std::vector<two_dim_read_var>& two_dim_vars,
                         const std::string& target_name, Real2D var);
static void apply_time_dependent_forcings_from_array(
    ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt, int sclr_dim, int edsclr_dim,
    const sclr_idx_type& sclr_idx, bool l_do_expldiff_rtm_thlm, const grid& gr, Real2D rtm,
    Real2D rho, Real2D exner, Real3D forcings_array, Real2D thlm_f, Real2D rtm_f, Real2D um_ref,
    Real2D vm_ref, Real2D um_f, Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg,
    Real3D sclrm_forcing, Real3D edsclrm_forcing, const time_dependent_input_state& input);

//================================================================================================
void initialize_t_dependent_input(ScratchWorkspace& workspace, // InOut
                                  int iunit, const std::string& runtype, bool l_readiopdata,
                                  Real2D rho,                                               // In
                                  int nz, Real1D grid, Real1D p_in_Pa,                      // In
                                  bool l_add_dycore_grid, int grid_adapt_in_time_method,    // In
                                  time_dependent_input_state& input, readiopdata_state& iop // InOut
) {
  // Description:
  // Read time-dependent case information into its owning context. The two
  // appended contexts replace Native module storage, including the IOP reader.
  auto workspace_scope                = workspace.scope();
  constexpr const char* input_path    = "../input/case_data/";
  constexpr const char* forcings_path = "_forcings.in";
  constexpr const char* sfc_path      = "_sfc.in";

  // Preserve Native TRIM(runtype) in each file path; leading blanks remain.
  const auto last_case_character = runtype.find_last_not_of(' ');
  const std::string case_name    = last_case_character == std::string::npos
                                       ? std::string{}
                                       : runtype.substr(0, (last_case_character + 1));

  if (input.l_ignore_forcings)
    input.forcing_times = {};
  if (!input.l_ignore_forcings) {
    initialize_t_dependent_forcings(workspace, iunit, input_path + case_name + forcings_path, nz,
                                    grid, p_in_Pa, l_add_dycore_grid, grid_adapt_in_time_method,
                                    input);
  }
  if (!input.l_sfc_already_initialized) {
    initialize_t_dependent_sfc(workspace, iunit, l_readiopdata, rho,
                               input_path + case_name + sfc_path, input, iop);
  }
}

//================================================================================================
void finalize_t_dependent_input(time_dependent_input_state& input) {
  // Description:
  // Free the time-dependent forcing and surface data in Native lifecycle order.
  if (!input.l_ignore_forcings)
    finalize_t_dependent_forcings(input);
  input.forcing_times = {};
  finalize_t_dependent_sfc(input);
}

//================================================================================================
static void initialize_t_dependent_sfc(ScratchWorkspace& workspace, int iunit, bool l_readiopdata,
                                       Real2D rho, const std::string& input_file,
                                       time_dependent_input_state& input, readiopdata_state& iop) {
  // Description:
  // Read surface input varying in one dimension. Convert NetCDF IOP heat
  // fluxes with the first density, or retain each named formatted input series.
  auto workspace_scope = workspace.scope();
  if (l_readiopdata) {
    Real1D time_execution;
    readiopdata_sfc(workspace, &input.T_sfc_given, nullptr, &input.sens_ht_given,
                    &input.latent_ht_given, &time_execution, iop);
    input.time_sfc_given   = to_host_space(time_execution);
    const int ntimes       = static_cast<int>(input.time_sfc_given.extent(0));
    input.wpthlp_sfc_given =
        Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "wpthlp_sfc_given"), ntimes);
    input.wpqtp_sfc_given =
        Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "wpqtp_sfc_given"), ntimes);
    const auto heat        = input.sens_ht_given;
    const auto water       = input.latent_ht_given;
    const auto heat_flux   = input.wpthlp_sfc_given;
    const auto water_flux  = input.wpqtp_sfc_given;
    using Policy           = Kokkos::RangePolicy<typename Real1D::execution_space>;
    Kokkos::parallel_for(
        "initialize_t_dependent_sfc::iop", Policy(0, ntimes), KOKKOS_LAMBDA(int i) {
          heat_flux(i)  = convert_sens_ht_to_km_s(heat(i), rho(0, 0));
          water_flux(i) = convert_latent_ht_to_m_s(water(i), rho(0, 0));
        });
  } else {
    const int nforcings = count_columns(iunit, input_file);
    std::vector<one_dim_read_var> retVars(nforcings);

    // Read the surface.in file, then fill -999.9 blanks by linear interpolation.
    read_one_dim_file(iunit, nforcings, input_file, retVars);
    fill_blanks_one_dim_vars(nforcings, retVars);
    const int dim_size = static_cast<int>(retVars[0].values.extent(0));

    // Store the data read from the file in each surface series, if present.
    if (get_target_index(nforcings, time_name, retVars) >= 0) {
      input.time_sfc_given = decltype(input.time_sfc_given)(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "time_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, time_name, retVars, input.time_sfc_given, input_file);
    }
    if (get_target_index(nforcings, latent_ht_name, retVars) >= 0) {
      input.latent_ht_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "latent_ht_given"), dim_size);
      read_x_profile(nforcings, dim_size, latent_ht_name, retVars, input.latent_ht_given,
                     input_file);
    }
    if (get_target_index(nforcings, sens_ht_name, retVars) >= 0) {
      input.sens_ht_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "sens_ht_given"), dim_size);
      read_x_profile(nforcings, dim_size, sens_ht_name, retVars, input.sens_ht_given, input_file);
    }
    if (get_target_index(nforcings, thetal_name, retVars) >= 0) {
      input.thlm_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "thlm_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, thetal_name, retVars, input.thlm_sfc_given, input_file);
    }
    if (get_target_index(nforcings, rt_name, retVars) >= 0) {
      input.rtm_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rtm_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, rt_name, retVars, input.rtm_sfc_given, input_file);
    }
    if (get_target_index(nforcings, CO2_umol_name, retVars) >= 0) {
      input.CO2_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "CO2_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, CO2_umol_name, retVars, input.CO2_sfc_given, input_file);
    }
    if (get_target_index(nforcings, upwp_sfc_name, retVars) >= 0) {
      input.upwp_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "upwp_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, upwp_sfc_name, retVars, input.upwp_sfc_given, input_file);
    }
    if (get_target_index(nforcings, vpwp_sfc_name, retVars) >= 0) {
      input.vpwp_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "vpwp_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, vpwp_sfc_name, retVars, input.vpwp_sfc_given, input_file);
    }
    if (get_target_index(nforcings, T_sfc_name, retVars) >= 0) {
      input.T_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "T_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, T_sfc_name, retVars, input.T_sfc_given, input_file);
    }
    if (get_target_index(nforcings, wpthlp_sfc_name, retVars) >= 0) {
      input.wpthlp_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "wpthlp_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, wpthlp_sfc_name, retVars, input.wpthlp_sfc_given,
                     input_file);
    }
    if (get_target_index(nforcings, wpqtp_sfc_name, retVars) >= 0) {
      input.wpqtp_sfc_given =
          Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "wpqtp_sfc_given"), dim_size);
      read_x_profile(nforcings, dim_size, wpqtp_sfc_name, retVars, input.wpqtp_sfc_given,
                     input_file);
    }

    // Deallocate the call-local input records.
    deallocate_one_dim_vars(nforcings, retVars);
  }
  input.l_sfc_already_initialized = true;
}

//================================================================================================
static void initialize_t_dependent_forcings(ScratchWorkspace& workspace, int iunit,
                                            const std::string& input_file, int nz, Real1D grid,
                                            Real1D p_in_Pa, bool l_add_dycore_grid,
                                            int grid_adapt_in_time_method,
                                            time_dependent_input_state& input) {
  // Description:
  // Read and interpolate a two-dimensional forcing table to the current grid.
  // Retain raw input only when later adaptive-grid initialization needs it.
  auto workspace_scope                  = workspace.scope();
  auto& t_dependent_forcing_data        = input.t_dependent_forcing_data;
  auto& t_dependent_forcing_data_f_grid = input.t_dependent_forcing_data_f_grid;
  auto& dimension_var                   = input.dimension_var;
  if (!t_dependent_forcing_data_f_grid[0].values.data()) {
    // Read only if the raw forcing data are not already stored in this case.
    read_two_dim_file(iunit, nforcings, input_file, t_dependent_forcing_data_f_grid, dimension_var);
  }
  const int n_f_grid_z = static_cast<int>(t_dependent_forcing_data_f_grid[0].values.extent(0));
  const int n_f_grid_t = static_cast<int>(dimension_var.values.extent(0));

  // Whole -999.9 profiles remain blank after the interpolation of other blanks.
  fill_blanks_two_dim_vars(nforcings, dimension_var, t_dependent_forcing_data_f_grid);
  input.forcing_times = to_host_space(dimension_var.values);
  for (int i = 0; i < nforcings; ++i) {
    if (!t_dependent_forcing_data[i].values.data()) {
      t_dependent_forcing_data[i].values =
          Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "t_dependent_forcing_data"), nz,
                 n_f_grid_t);
    }
  }
  const auto target_grid = t_dependent_forcing_data[0].values;
  using Policy           = Kokkos::RangePolicy<typename Real2D::execution_space>;
  if (t_dependent_forcing_data_f_grid[0].name == z_name) {
    t_dependent_forcing_data[0].name = z_name;
    Kokkos::parallel_for(
        "initialize_t_dependent_forcings::height", Policy(0, (n_f_grid_t * nz)),
        KOKKOS_LAMBDA(int index) { target_grid((index % nz), (index / nz)) = grid((index % nz)); });
  } else if (t_dependent_forcing_data_f_grid[0].name == pressure_name) {
    t_dependent_forcing_data[0].name = pressure_name;
    Kokkos::parallel_for(
        "initialize_t_dependent_forcings::pressure", Policy(0, (n_f_grid_t * nz)),
        KOKKOS_LAMBDA(int index) {
          target_grid((index % nz), (index / nz)) = -p_in_Pa((index % nz));
        });
    const auto raw_grid = t_dependent_forcing_data_f_grid[0].values;
    Kokkos::parallel_for(
        "initialize_t_dependent_forcings::negative_raw_pressure",
        Policy(0, (n_f_grid_t * n_f_grid_z)), KOKKOS_LAMBDA(int index) {
          const int k    = (index % n_f_grid_z);
          const int t    = (index / n_f_grid_z);
          raw_grid(k, t) = -raw_grid(k, t);
        });
  } else {
    throw std::runtime_error("Incompatible grid type in first element of t_dependent_forcings.");
  }

  // Interpolate the time-dependent input data to the appropriate grid.
  const auto grid_profile = Kokkos::subview(target_grid, Kokkos::ALL(), 0);
  for (int i = 1; i < nforcings; ++i) {
    t_dependent_forcing_data[i].name = t_dependent_forcing_data_f_grid[i].name;
    read_to_grid(workspace, nforcings, n_f_grid_z, n_f_grid_t, nz, grid_profile,
                 t_dependent_forcing_data_f_grid, t_dependent_forcing_data[i].name,
                 t_dependent_forcing_data[i].values);
  }
  if (l_add_dycore_grid || grid_adapt_in_time_method == no_grid_adaptation) {
    // Raw input is needed only for adaptation without simulated dycore forcing.
    for (int i = 0; i < nforcings; ++i)
      t_dependent_forcing_data_f_grid[i].values = Real2D{};
  }
}

//================================================================================================
static void finalize_t_dependent_forcings(time_dependent_input_state& input) {
  // Description:
  // Release raw forcing data retained for adaptation, current-grid profiles,
  // and the other-dimension record, preserving the remaining metadata/flags.
  for (int i = 0; i < nforcings; ++i)
    input.t_dependent_forcing_data_f_grid[i].values = Real2D{};
  for (int i = 0; i < nforcings; ++i)
    input.t_dependent_forcing_data[i].values = Real2D{};
  input.dimension_var.values = Real1D{};
  input.forcing_times        = {};
}

//================================================================================================
static void finalize_t_dependent_sfc(time_dependent_input_state& input) {
  // Description:
  // Clear all surface series allocated during input initialization. Native's
  // l_sfc_already_initialized flag is reset by the owning driver's next init.
  input.time_sfc_given   = {};
  input.latent_ht_given  = Real1D{};
  input.sens_ht_given    = Real1D{};
  input.thlm_sfc_given   = Real1D{};
  input.rtm_sfc_given    = Real1D{};
  input.CO2_sfc_given    = Real1D{};
  input.upwp_sfc_given   = Real1D{};
  input.vpwp_sfc_given   = Real1D{};
  input.T_sfc_given      = Real1D{};
  input.wpthlp_sfc_given = Real1D{};
  input.wpqtp_sfc_given  = Real1D{};
}

//================================================================================================
template <class Grid>
static void read_to_grid(ScratchWorkspace& workspace, int ntwo_dim_vars, int dim_size,
                         int other_dim_size, int nz, Grid grid,
                         const std::vector<two_dim_read_var>& two_dim_vars,
                         const std::string& target_name, Real2D var) {
  // Description:
  // Translate a forcing variable from its supplied grid to the model grid by
  // applying the permanent Native-compatible vertical interpolation per time.
  auto workspace_scope = workspace.scope();
  auto temp_var = workspace.view<Real2D>(Kokkos::WithoutInitializing, "read_to_grid::temp_var",
                                         dim_size, other_dim_size);
  read_x_table(ntwo_dim_vars, dim_size, other_dim_size, target_name, two_dim_vars, temp_var);
  for (int i = 0; i < other_dim_size; ++i) {
    zlinterp_fnc(nz, dim_size, grid, Kokkos::subview(two_dim_vars[0].values, Kokkos::ALL(), i),
                 Kokkos::subview(temp_var, Kokkos::ALL(), i),
                 Kokkos::subview(var, Kokkos::ALL(), i));
  }
}

//================================================================================================
static void apply_time_dependent_forcings_from_array(
    ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt, int sclr_dim, int edsclr_dim,
    const sclr_idx_type& sclr_idx, bool l_do_expldiff_rtm_thlm, const grid& gr, Real2D rtm,
    Real2D rho, Real2D exner, Real3D forcings_array, Real2D thlm_f, Real2D rtm_f, Real2D um_ref,
    Real2D vm_ref, Real2D um_f, Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg,
    Real3D sclrm_forcing, Real3D edsclrm_forcing, const time_dependent_input_state& input) {
  // Description:
  // Parse one set of forcing profiles per column and update the corresponding
  // model tendencies/references/winds. Any dummy value suppresses that whole
  // variable profile, preserving legitimate forcing from other sources.
  auto workspace_scope = workspace.scope();
  (void)sclr_dim;
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Parse the stored values into CLUBB-compatible forcing data.
  for (int n = 1; n < nforcings; ++n) {
    const auto temp_array = Kokkos::subview(forcings_array, Kokkos::ALL(), n, Kokkos::ALL());

    // A dummy profile must not overwrite information from another source.
    int dummy_count = 0;
    Kokkos::parallel_reduce(
        "apply_time_dependent_forcings_from_array::dummy", Policy(0, (ngrdcol * nzt)),
        KOKKOS_LAMBDA(int index, int& count) {
          const real_t value = temp_array((index % ngrdcol), (index / ngrdcol));
          if (Kokkos::abs(value - (-real_t(999.9))) <
              Kokkos::abs(value + (-real_t(999.9))) / 2 * eps)
            ++count;
        },
        dummy_count);
    if (dummy_count != 0)
      continue;

    const auto& name = input.t_dependent_forcing_data[n].name;
    if (name == temperature_f_name || name == theta_f_name || name == thetal_f_name) {
      if (name == temperature_f_name) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::temperature", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i  = (index % ngrdcol);
              const int k  = (index / ngrdcol);
              thlm_f(i, k) = temp_array(i, k) / exner(i, k);
            });
      } else if (name == theta_f_name) {
        // Native applies no additional conversion for potential temperature.
        Kokkos::deep_copy(thlm_f, temp_array);
      } else if (name == thetal_f_name) {
        Kokkos::deep_copy(thlm_f, temp_array);
      }
      if (sclr_idx.iisclr_thl > -1) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::passive_thl", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i                              = (index % ngrdcol);
              const int k                              = (index / ngrdcol);
              sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_f(i, k);
            });
      }
      if (sclr_idx.iiedsclr_thl > -1) {
        // With explicit diffusion, advance_windm_edsclrm supplies a copy of
        // thlm already forced in advance_xm_wpxp: do not force the copy twice.
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::eddy_thl", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i = (index % ngrdcol);
              const int k = (index / ngrdcol);
              if (l_do_expldiff_rtm_thlm && edsclr_dim > 1)
                edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = zero;
              else
                edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_f(i, k);
            });
      }
    } else if (name == rt_f_name || name == sp_humidity_f_name) {
      if (name == sp_humidity_f_name) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::specific_humidity",
            Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
              const int i         = (index % ngrdcol);
              const int k         = (index / ngrdcol);
              const real_t factor = real_t(1.0) + rtm(i, k);
              rtm_f(i, k)         = temp_array(i, k) * (factor * factor);
            });
      } else if (name == rt_f_name) {
        Kokkos::deep_copy(rtm_f, temp_array);
      }
      if (sclr_idx.iisclr_rt > -1) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::passive_rt", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i                             = (index % ngrdcol);
              const int k                             = (index / ngrdcol);
              sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_f(i, k);
            });
      }
      if (sclr_idx.iiedsclr_rt > -1) {
        // Copy of rtm with explicit diffusion: see the theta_l case above.
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::eddy_rt", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i = (index % ngrdcol);
              const int k = (index / ngrdcol);
              if (l_do_expldiff_rtm_thlm && edsclr_dim > 1)
                edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = zero;
              else
                edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_f(i, k);
            });
      }
    } else if (name == um_ref_name) {
      Kokkos::deep_copy(um_ref, temp_array);
    } else if (name == vm_ref_name) {
      Kokkos::deep_copy(vm_ref, temp_array);
    } else if (name == um_f_name) {
      Kokkos::deep_copy(um_f, temp_array);
    } else if (name == vm_f_name) {
      Kokkos::deep_copy(vm_f, temp_array);
    } else if (name == wm_name || name == omega_name || name == omega_mb_hr_name) {
      if (name == wm_name) {
        Kokkos::deep_copy(wm_zt, temp_array);
      } else if (name == omega_name) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::omega", Policy(0, (ngrdcol * nzt)),
            KOKKOS_LAMBDA(int index) {
              const int i = (index % ngrdcol);
              const int k = (index / ngrdcol);
              wm_zt(i, k) = -temp_array(i, k) / (grav * rho(i, k));
            });
      } else if (name == omega_mb_hr_name) {
        Kokkos::parallel_for(
            "apply_time_dependent_forcings_from_array::omega_mb_hr", Policy(0, ngrdcol * nzt),
            KOKKOS_LAMBDA(int index) {
              const int i = index % ngrdcol;
              const int k = index / ngrdcol;
              const real_t pressure_tendency = temp_array(i, k) * pascal_per_mb / sec_per_hr;
              wm_zt(i, k)                    = -pressure_tendency / (grav * rho(i, k));
            });
      }
      zt2zm_api(nzm, nzt, ngrdcol, gr, wm_zt, wm_zm);
    } else if (name == ug_name) {
      Kokkos::deep_copy(ug, temp_array);
    } else if (name == vg_name) {
      Kokkos::deep_copy(vg, temp_array);
    } else {
      std::cerr << "Incompatable forcing type: " << name << '\n';
      throw std::runtime_error("");
    }
  }
}

//================================================================================================
void apply_time_dependent_forcings(ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt,
                                   int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                                   bool l_do_expldiff_rtm_thlm, const grid& gr, double time,
                                   Real2D rtm, Real2D rho, Real2D exner, Real2D thlm_f,
                                   Real2D rtm_f, Real2D um_ref, Real2D vm_ref, Real2D um_f,
                                   Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg,
                                   Real3D sclrm_forcing, Real3D edsclrm_forcing,
                                   const time_dependent_input_state& input) {
  // Description:
  // Interpolate stored forcing information in time and apply it in the format
  // used by CLUBB. Each column receives the same supplied physics-grid profile.
  auto workspace_scope = workspace.scope();
  // Slot zero is grid metadata and unread here; all consumed forcing cells are filled below.
  Real3D forcings_array = workspace_scope.view<Real3D>(
      Kokkos::WithoutInitializing, "apply_time_dependent_forcings::forcings_array", ngrdcol,
      nforcings, nzt);
  real_t time_frac = -one;
  int before_time, after_time;
  time_select(time, static_cast<int>(input.forcing_times.extent(0)), input.forcing_times,
              before_time, after_time, time_frac);

  // The input is on one grid, so every column gets the same profile.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  for (int n = 1; n < nforcings; ++n) {
    const auto data = input.t_dependent_forcing_data[n].values;
    Kokkos::parallel_for(
        "apply_time_dependent_forcings::interpolate", Policy(0, ngrdcol * nzt),
        KOKKOS_LAMBDA(int index) {
          const int i = index % ngrdcol;
          const int k = index / ngrdcol;
          forcings_array(i, n, k) =
              linear_interp_factor(time_frac, data(k, after_time), data(k, before_time));
        });
  }

  apply_time_dependent_forcings_from_array(
      workspace, ngrdcol, nzm, nzt, sclr_dim, edsclr_dim, sclr_idx, l_do_expldiff_rtm_thlm, gr, rtm,
      rho, exner, forcings_array, thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, wm_zt, wm_zm, ug, vg,
      sclrm_forcing, edsclrm_forcing, input);
}

//================================================================================================
void apply_time_dependent_forcings_from_dycore(
    ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt, int sclr_dim, int edsclr_dim,
    const sclr_idx_type& sclr_idx, bool l_do_expldiff_rtm_thlm, const grid& gr,
    const grid& gr_dycore, double time, Real2D rtm, Real2D rho, Real2D exner, int grid_remap_method,
    int total_idx_rho_lin_spline, Real2D rho_lin_spline_vals, Real2D rho_lin_spline_levels,
    Real1D p_sfc, Real2D thlm_f, Real2D rtm_f, Real2D um_ref, Real2D vm_ref, Real2D um_f,
    Real2D vm_f, Real2D wm_zt, Real2D wm_zm, Real2D ug, Real2D vg, Real3D sclrm_forcing,
    Real3D edsclrm_forcing, const time_dependent_input_state& input) {
  // Description:
  // Interpolate time-dependent forcing stored on the dycore grid, remap each
  // profile to the current physics grid, then apply the corresponding fields.
  auto workspace_scope     = workspace.scope();
  const int dycore_nzt = gr_dycore.nzt;
  Real2D temp_array_dycore = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "apply_time_dependent_forcings_from_dycore::temp_array", ngrdcol,
      dycore_nzt);
  // Slot zero is grid metadata and unread here; all consumed forcing cells are filled below.
  Real3D forcings_array = workspace_scope.view<Real3D>(
      Kokkos::WithoutInitializing, "apply_time_dependent_forcings_from_dycore::forcings_array",
      ngrdcol, nforcings, nzt);
  real_t time_frac         = -one;
  const int iv             = 1;
  const bool l_zt_variable = true;
  int before_time, after_time;
  time_select(time, static_cast<int>(input.forcing_times.extent(0)), input.forcing_times,
              before_time, after_time, time_frac);

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  for (int n = 1; n < nforcings; ++n) {
    // Each column starts with the same dycore profile and is remapped to its
    // own physics grid, preserving the Native per-variable remap call order.
    const auto data = input.t_dependent_forcing_data[n].values;
    Kokkos::parallel_for(
        "apply_time_dependent_forcings_from_dycore::interpolate", Policy(0, ngrdcol * dycore_nzt),
        KOKKOS_LAMBDA(int index) {
          const int i = index % ngrdcol;
          const int k = index / ngrdcol;
          temp_array_dycore(i, k) =
              linear_interp_factor(time_frac, data(k, after_time), data(k, before_time));
        });
    remap_vals_to_target(workspace, ngrdcol, gr_dycore, gr, dycore_nzt, temp_array_dycore, gr.nzt,
                         total_idx_rho_lin_spline, rho_lin_spline_vals, rho_lin_spline_levels, iv,
                         p_sfc, grid_remap_method, l_zt_variable,
                         Kokkos::subview(forcings_array, Kokkos::ALL(), n, Kokkos::ALL()));
  }

  apply_time_dependent_forcings_from_array(
      workspace, ngrdcol, nzm, nzt, sclr_dim, edsclr_dim, sclr_idx, l_do_expldiff_rtm_thlm, gr, rtm,
      rho, exner, forcings_array, thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, wm_zt, wm_zm, ug, vg,
      sclrm_forcing, edsclrm_forcing, input);
}

//================================================================================================
void time_select(
    double time, int nvar,                                                         // In
    Kokkos::View<const real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> time_array, // In
    int& before_time, int& after_time, real_t& time_frac                           // Out
) {
  // Description:
  // Select zero-based bracketing slots and the fraction of time [s] between
  // them. Preserve Native endpoint selection, full sortedness check and fatal
  // out-of-range diagnostics. All supplied time metadata resides on the host.
  // time/time_array are in seconds; time_frac is dimensionless.
  before_time            = -999;
  after_time             = -999;
  const real_t core_time = static_cast<real_t>(time);
  if (core_time < time_array(0)) {
    throw std::runtime_error("In subroutine time_select:\n"
                             "Selected time is before the first (begin) time\n"
                             "at which data are available.  Cannot interpolate.");
  } else if (core_time > time_array(nvar - 1)) {
    throw std::runtime_error("In subroutine time_select:\n"
                             "Selected time is after the last (end) time\n"
                             "at which data are available.  Cannot interpolate.");
  } else {
    for (int k = 0; k < nvar - 1; ++k) {
      if (time_array(k) >= time_array(k + 1)) {
        throw std::runtime_error("In subroutine time_select:\n"
                                 "times are not sorted. Check (case)_sfc.in\n"
                                 "and (case)_forcings.in, located in input/case_data.");
      }
      if (before_time == -999 && core_time >= time_array(k) && core_time <= time_array(k + 1)) {
        before_time = k;
        after_time  = k + 1;
      }
    }
  }

  // Compute the position of time between the selected times as a fraction.
  time_frac =
      (core_time - time_array(before_time)) / (time_array(after_time) - time_array(before_time));
}

} // namespace clubb_kokkos
