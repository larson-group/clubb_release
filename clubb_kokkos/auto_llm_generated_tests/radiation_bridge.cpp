// Verify that C++ radiation retains heating between radiation calls.
// advance_clubb_radiation updates radht only on radiation timesteps. With stats
// and soil/vegetation disabled, an intervening call must preserve every value.
#include "scratch_workspace_owner.hpp"
#include "Radiation/radiation_module.hpp"
#include "CLUBB_core/error_code.hpp"
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
  // Exercise the ported routine with l_rad_itime=false and a nonzero,
  // column-dependent heating rate; no radiation solver needs initialization.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    grid gr;
    gr.nzt             = 3;
    gr.nzm             = 4;
    int ngrdcol        = 2;
    int hydromet_dim   = 0;
    int pdf_dim        = 0;
    int lh_num_samples = 0;
    bool l_rad_itime   = false;
    real_t dt_main     = real_t(30);
    int day{};
    int month{};
    int year{};
    real_t lat_vals{};
    real_t lon_vals{};
    double time_current{};
    double time_initial{};
    Kokkos::View<real_t**, Kokkos::LayoutLeft> rho("rho", ngrdcol, gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> rho_zm("rho_zm", ngrdcol, gr.nzm);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> p_in_Pa("p_in_Pa", ngrdcol, gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> exner("exner", ngrdcol, gr.nzt);
    Kokkos::View<real_t*, Kokkos::LayoutLeft> wpthlp_sfc("wpthlp_sfc", ngrdcol);
    Kokkos::View<real_t*, Kokkos::LayoutLeft> wprtp_sfc("wprtp_sfc", ngrdcol);
    Kokkos::View<real_t*, Kokkos::LayoutLeft> p_sfc("p_sfc", ngrdcol);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> cloud_frac("cloud_frac", ngrdcol, gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> ice_supersat_frac("ice_supersat_frac", ngrdcol,
                                                                 gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> thlm("thlm", ngrdcol, gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> rtm("rtm", ngrdcol, gr.nzt);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> rcm("rcm", ngrdcol, gr.nzt);
    Kokkos::View<real_t****, Kokkos::LayoutLeft> X_nl_all_levs("X_nl_all_levs", ngrdcol,
                                                               lh_num_samples, gr.nzt, pdf_dim);
    Kokkos::View<real_t***, Kokkos::LayoutLeft> lh_rt_clipped("lh_rt_clipped", ngrdcol,
                                                              lh_num_samples, gr.nzt);
    Kokkos::View<real_t***, Kokkos::LayoutLeft> lh_thl_clipped("lh_thl_clipped", ngrdcol,
                                                               lh_num_samples, gr.nzt);
    Kokkos::View<real_t***, Kokkos::LayoutLeft> lh_rc_clipped("lh_rc_clipped", ngrdcol,
                                                              lh_num_samples, gr.nzt);
    Kokkos::View<real_t***, Kokkos::LayoutLeft> lh_sample_point_weights(
        "lh_sample_point_weights", ngrdcol, lh_num_samples, gr.nzt);
    Kokkos::View<real_t***, Kokkos::LayoutLeft> hydromet("hydromet", ngrdcol, gr.nzt, hydromet_dim);
    hm_metadata_type hm_metadata{};
    clubb_kokkos::stats_type stats{};
    err_info_type err_info{};
    Kokkos::View<real_t*, Kokkos::LayoutLeft> deep_soil_T_in_K("deep_soil_T_in_K", ngrdcol);
    Kokkos::View<real_t*, Kokkos::LayoutLeft> sfc_soil_T_in_K("sfc_soil_T_in_K", ngrdcol);
    Kokkos::View<real_t*, Kokkos::LayoutLeft> veg_T_in_K("veg_T_in_K", ngrdcol);
    Kokkos::View<real_t**, Kokkos::LayoutLeft> radht("radht", ngrdcol, gr.nzt);

    parameters_radiation_state parameters_radiation;
    radiation_variables_state radiation_variables;
    simple_rad_state simple_rad;
    soil_vegetation_state soil_vegetation;
    initialize_soil_veg(ngrdcol, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
                        soil_vegetation);
    init_default_err_info_api(workspace, ngrdcol, err_info);
    auto before = Kokkos::create_mirror_view(radht);
    for (int k = 0; k < gr.nzt; ++k)
      for (int i = 0; i < ngrdcol; ++i)
        before(i, k) = real_t(11 + i + 10 * k);
    Kokkos::deep_copy(radht, before);

    advance_clubb_radiation(workspace, gr, ngrdcol, hydromet_dim, pdf_dim, lh_num_samples, // In
                            l_rad_itime, dt_main, day, month, year,                        // In
                            lat_vals, lon_vals, time_current, time_initial,                // In
                            rho, rho_zm, p_in_Pa, exner,                                   // In
                            wpthlp_sfc, wprtp_sfc, p_sfc,                                  // In
                            cloud_frac, ice_supersat_frac,                                 // In
                            thlm, rtm, rcm,                                                // In
                            X_nl_all_levs, lh_rt_clipped, lh_thl_clipped,
                            lh_rc_clipped,                                  // In
                            lh_sample_point_weights, hydromet, hm_metadata, // In
                            stats, err_info,                                // InOut
                            deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,  // InOut
                            radht,                                          // InOut
                            parameters_radiation, radiation_variables, simple_rad,
                            soil_vegetation); // Native module state

    const auto after = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, radht);
    for (int k = 0; k < gr.nzt; ++k)
      for (int i = 0; i < ngrdcol; ++i)
        if (after(i, k) != before(i, k))
          throw std::runtime_error("C++ radiation discarded the retained heating rate");
    if (any_equal(err_info.err_code, clubb_fatal_error))
      throw std::runtime_error("C++ radiation reported a fatal error");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
