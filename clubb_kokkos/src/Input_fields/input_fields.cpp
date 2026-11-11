// Read restart fields and interpret external profile files for one CLUBB case.
//
// This mirrors all eight Native inputfields bodies. Public routines keep their
// Native order and one-column profile ranks; the appended inputfields_state
// replaces actual module settings. Raw namelist/NetCDF operations alone use f2k.
#include "scratch_workspace.hpp"
#include "input_fields.hpp"
#include "stat_file_utils.hpp"
#include "extrapolation.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "driver_io.hpp"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
namespace clubb_kokkos {
namespace {
constexpr int coamps_input_type = 1, clubb_input_type = 2;
constexpr int sam_input_type = 3, rams_input_type = 4;
constexpr int clubb_zt = 1, clubb_zm = 2, clubb_sfc = 3;
constexpr int coamps_sm = 0, coamps_sw = 1, sam_file = 0, rams_file = 0;
constexpr int num_sam_inputfields = 70, num_coamps = 61, num_rams_inputfields = 4;
// Native private metadata is call-local. The View aliases the destination field.
struct input_field {
  bool l_input_var;
  std::string input_name, clubb_name; // Native fixed50.
  std::string clubb_grid_type;        // Native fixed3: zt or zm.
  int input_file_index;               // Zero-based stat_files slot.
  Real1D clubb_var;
  real_t adjustment;
};
} // namespace
static void get_input_variables_interp(
    ScratchWorkspace& workspace,                                                        // InOut
    const grid& gr, int nvars, const std::vector<input_field>& variables, int timestep, // In
    std::vector<int>& k_lowest_zt, std::vector<int>& k_highest_zt,                      // Out
    std::vector<int>& k_lowest_zm, std::vector<int>& k_highest_zm, bool& l_error,       // Out
    const inputfields_state& inputfields); // Native module In
template <class LESProfile>
static void inputfields_interp_and_adjust(ScratchWorkspace& workspace, const grid& gr,
                                          const std::string& clubb_grid_type,
                                          const stat_file& fread_var, // In
                                          real_t adjustment, LESProfile LES_tmp, int k_lowest_input,
                                          int k_highest_input, // In
                                          Real1D clubb_var);   // InOut

void set_filenames(const std::string& file_prefix, // In
                   inputfields_state& inputfields) // Native module InOut
{
  // Set the NetCDF file names used by inputfields and restart runs.
  // Native fixed120 assignment truncation remains at each file-name boundary.
  inputfields.stat_files.clear();
  inputfields.stats_input_type = -1;
  std::string input_type       = inputfields.input_type;
  input_type.erase(input_type.find_last_not_of(' ') + 1);
  std::string prefix = file_prefix;
  prefix.erase(prefix.find_last_not_of(' ') + 1);
  if (input_type == "coamps_les" || input_type == "COAMPS_LES" || input_type == "coamps") {
    inputfields.stat_files.resize(2);
    inputfields.stat_files[coamps_sm] = (prefix + "_coamps_sm.nc").substr(0, 120);
    inputfields.stat_files[coamps_sw] = (prefix + "_coamps_sw.nc").substr(0, 120);
    inputfields.stats_input_type      = coamps_input_type;
  } else if (input_type == "clubb_scm" || input_type == "CLUBB_SCM" || input_type == "clubb") {
    inputfields.stat_files.resize(3);
    inputfields.stat_files[(clubb_zt - 1)]  = (prefix + "_zt.nc").substr(0, 120);
    inputfields.stat_files[(clubb_zm - 1)]  = (prefix + "_zm.nc").substr(0, 120);
    inputfields.stat_files[(clubb_sfc - 1)] = (prefix + "_sfc.nc").substr(0, 120);
    inputfields.stats_input_type            = clubb_input_type;
    const bool l_stats_netcdf = std::filesystem::exists(inputfields.stat_files[(clubb_zt - 1)]);
    if (!l_stats_netcdf) {
      std::cerr << "inputfields: Expected netCDF stats file not found: " << prefix << "_stats.nc\n";
      inputfields.stat_files[(clubb_zt - 1)]  = (prefix + "_stats.nc").substr(0, 120);
      inputfields.stat_files[(clubb_zm - 1)]  = (prefix + "_stats.nc").substr(0, 120);
      inputfields.stat_files[(clubb_sfc - 1)] = (prefix + "_stats.nc").substr(0, 120);
    }
  } else if (input_type == "SAM" || input_type == "sam") {
    inputfields.stat_files.resize(1);
    // SAM uses one NetCDF output file; file_prefix already includes its suffix.
    inputfields.stat_files[sam_file] = prefix.substr(0, 120);
    inputfields.stats_input_type     = sam_input_type;
  } else if (input_type == "rams_les" || input_type == "RAMS_LES" || input_type == "rams") {
    inputfields.stat_files.resize(1);
    // RAMS uses one output file.
    inputfields.stat_files[rams_file] = (prefix + "_rams.nc").substr(0, 120);
    inputfields.stats_input_type      = rams_input_type;
  } else {
    std::cerr << "Don't know how to handle input_type = " << inputfields.input_type << '\n';
    throw std::runtime_error("STOP in set_filenames");
  }
}

void stat_fields_reader(ScratchWorkspace& workspace, // InOut
                        const grid& gr, int timestep, int hydromet_dim,
                        const hm_metadata_type& hm_metadata,                         // In
                        const std::string& microphys_scheme, bool l_predict_Nc,      // In
                        Real1D um, Real1D upwp, Real1D vm, Real1D vpwp,              // InOut
                        Real1D up2, Real1D vp2, Real1D rtm,                          // InOut
                        Real1D wprtp, Real1D thlm, Real1D wpthlp,                    // InOut
                        Real1D rtp2, Real1D rtp3,                                    // InOut
                        Real1D thlp2, Real1D thlp3, Real1D rtpthlp,                  // InOut
                        Real1D wp2, Real1D wp3,                                      // InOut
                        Real1D p_in_Pa, Real1D exner, Real1D rcm, Real1D cloud_frac, // InOut
                        Real1D wpthvp, Real1D wp2thvp, Real1D wp2up, Real1D rtpthvp,
                        Real1D thlpthvp,                                               // InOut
                        Real1D wp2rtp, Real1D wp2thlp, Real1D uprcp, Real1D vprcp,     // InOut
                        Real1D rc_coef_zm, Real1D wp4, Real1D wpup2,                   // InOut
                        Real1D wpvp2, Real1D wp2up2,                                   // InOut
                        Real1D wp2vp2, Real1D ice_supersat_frac,                       // InOut
                        Real1D wm_zt, Real1D rho, Real1D rho_zm, Real1D rho_ds_zm,     // InOut
                        Real1D rho_ds_zt, Real1D thv_ds_zm, Real1D thv_ds_zt,          // InOut
                        Real1D thlm_forcing, Real1D rtm_forcing, Real1D wprtp_forcing, // InOut
                        Real1D wpthlp_forcing, Real1D rtp2_forcing,                    // InOut
                        Real1D thlp2_forcing, Real1D rtpthlp_forcing,                  // InOut
                        Real2D hydromet, Real2D hydrometp2, Real2D wphydrometp,        // InOut
                        Real1D Ncm, Real1D Nccnm, Real1D thvm, Real1D em,              // InOut
                        Real1D tau_zm, Real1D tau_zt,                                  // InOut
                        Real1D Kh_zt, Real1D Kh_zm, Real1D ug, Real1D vg,              // InOut
                        Real1D thlprcp,                                                // InOut
                        Real1D sigma_sqd_w, Real1D sigma_sqd_w_zt, Real1D radht,       // InOut
                        real_t& deep_soil_T_in_K, real_t& sfc_soil_T_in_K,
                        real_t& veg_T_in_K,                                      // InOut
                        pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm, // InOut
                        const inputfields_state& inputfields)                    // Native module In
{
  // Description:
  //   Read variables for the model from statistical data.
  // The Native API is a single-column profile boundary; PDF writes retain row1.
  // um: eastward grid-mean wind component (thermo. levs.)  [m/s]
  // vm: northward grid-mean wind component (thermo. levs.) [m/s]
  // rtm: total water mixing ratio, r_t (thermo. levels)     [kg/kg]
  // thlm: liq. water pot. temp., th_l (thermo. levels)       [K]
  // rtp3: r_t'^3 (thermodynamic levels)                      [(kg/kg)^3]
  // thlp3: th_l'^3 (thermodynamic levels)                     [K^3]
  // wp3: w'^3 (thermodynamic levels)                        [m^3/s^3]
  // p_in_Pa: Air pressure (thermodynamic levels)                [Pa]
  // exner: Exner function (thermodynamic levels)              [-]
  // rcm: cloud water mixing ratio, r_c (thermo. levels)     [kg/kg]
  // cloud_frac: cloud fraction (thermodynamic levels)              [-]
  // wp2thvp: < w'^2 th_v' > (thermodynamic levels)              [m^2/s^2 K]
  // wp2up: < w'^2 u' > (thermodynamic levels)                 [m^3/s^3]
  // wp2rtp: w'^2 rt' (thermodynamic levels)                    [m^2/s^2 kg/kg]
  // wp2thlp: w'^2 thl' (thermodynamic levels)                   [m^2/s^2 K]
  // wpup2: w'u'^2 (thermodynamic levels)                      [m^3/s^3]
  // wpvp2: w'v'^2 (thermodynamic levels)                      [m^3/s^3]
  // ice_supersat_frac: ice cloud fraction (thermo. levels)                [-]
  // upwp: u'w' (momentum levels)                         [m^2/s^2]
  // vpwp: v'w' (momentum levels)                         [m^2/s^2]
  // up2: u'^2 (momentum levels)                         [m^2/s^2]
  // vp2: v'^2 (momentum levels)                         [m^2/s^2]
  // wprtp: w' r_t' (momentum levels)                      [kg/kg m/s]
  // wpthlp: w'th_l' (momentum levels)                      [(m/s) K]
  // rtp2: r_t'^2 (momentum levels)                       [(kg/kg)^2]
  // thlp2: th_l'^2 (momentum levels)                      [K^2]
  // rtpthlp: r_t'th_l' (momentum levels)                    [(kg/kg) K]
  // wp2: w'^2 (momentum levels)                         [m^2/s^2]
  // wpthvp: < w' th_v' > (momentum levels)                 [kg/kg K]
  // rtpthvp: < r_t' th_v' > (momentum levels)               [kg/kg K]
  // thlpthvp: < th_l' th_v' > (momentum levels)              [K^2]
  // uprcp: < u' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
  // vprcp: < v' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
  // rc_coef_zm: Coef of X'r_c' in Eq. (34) (m-levs.)           [K/(kg/kg)]
  // wp4: w'^4 (momentum levels)                         [m^4/s^4]
  // wp2up2: w'^2 u'^2 (momentum levels)                    [m^4/s^4]
  // wp2vp2: w'^2 v'^2 (momentum levels)                    [m^4/s^4]
  // wm_zt: vertical mean wind component on thermo. levels  [m/s]
  // rho: Air density on thermodynamic levels             [kg/m^3]
  // rho_ds_zt: Dry, static density on thermo. levels           [kg/m^3]
  // thv_ds_zt: Dry, base-state theta_v on thermo levels        [K]
  // rho_zm: Air density on momentum levels                  [kg/m^3]
  // rho_ds_zm: Dry, static density on momentum levels          [kg/m^3]
  // thv_ds_zm: Dry, base-state theta_v on momentum levels      [K]
  // thlm_forcing: liquid potential temp. forcing (thermo. levels) [K/s]
  // rtm_forcing: total water forcing (thermo. levels)      [(kg/kg)/s]
  // wprtp_forcing: total water turbulent flux forcing (m-levs) [m*K/s^2]
  // wpthlp_forcing: liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
  // rtp2_forcing: total water variance forcing (m-levs)   [(kg/kg)^2/s]
  // thlp2_forcing: liq pot temp variance forcing (m-levs)  [K^2/s]
  // rtpthlp_forcing: <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]
  // hydromet: Array of hydrometeors                [hm units]
  // hydrometp2: Variance of a hydrometeor (m-levs.)  [<hm units>^2]
  // wphydrometp: Covariance of w and a hydrometeor    [(m/s) <hm units>]
  // Ncm: Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
  // Nccnm: Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
  // thvm: Virtual potential temperature                        [K]
  // tau_zt: Eddy dissipation time scale on thermodynamic levels  [s]
  // Kh_zt: Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
  // ug: u geostrophic wind                                   [m/s]
  // vg: v geostrophic wind                                   [m/s]
  // sigma_sqd_w_zt: PDF width parameter interpolated to t-levs.          [-]
  // radht: SW + LW heating rate                                 [K/s]
  // thlprcp: thl'rc'                                              [K kg/kg]
  // sigma_sqd_w: PDF width parameter (momentum levels)                [-]
  // em: Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
  // tau_zm: Eddy dissipation time scale on momentum levels       [s]
  // Kh_zm: Eddy diffusivity coefficient on momentum levels      [m^2/s]
  // pdf_params: PDF parameters (thermodynamic levels)    [units vary]
  // pdf_params_zm: PDF parameters on momentum levels        [units vary]
  auto workspace_scope = workspace.scope();
  bool l_read_error, l_fatal_error;
  int k, iirr, iiNr, iirs, iiri, iirg, iiNi, iiNg, iiNs;
  Real1D tmp1 =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "stat_fields_reader::tmp1", gr.nzm);
  Real1D temp_Nrm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Nrm", gr.nzt);
  Real1D temp_Ncm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Ncm", gr.nzt);
  Real1D temp_rgm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_rgm", gr.nzt);
  Real1D temp_rim      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_rim", gr.nzt);
  Real1D temp_Ngm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Ngm", gr.nzt);
  Real1D temp_Nsm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Nsm", gr.nzt);
  Real1D temp_rrm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_rrm", gr.nzt);
  Real1D temp_rsm      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_rsm", gr.nzt);
  Real1D temp_Nim      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Nim", gr.nzt);
  Real1D up2_zt        = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::up2_zt", gr.nzt);
  Real1D vp2_zt        = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::vp2_zt", gr.nzt);
  Real1D temp_tke      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_tke", gr.nzm);
  Real1D temp_wpup_sgs = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_wpup_sgs", gr.nzm);
  Real1D temp_wpvp_sgs = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_wpvp_sgs", gr.nzm);
  Real1D temp_rrp2     = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_rrp2", gr.nzm);
  Real1D temp_Nrp2     = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_Nrp2", gr.nzm);
  Real1D temp_wprrp    = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_wprrp", gr.nzm);
  Real1D temp_wpNrp    = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "stat_fields_reader::temp_wpNrp", gr.nzm);
  Real1D zero_height   = workspace_scope.view<Real1D>("stat_fields_reader::zero_height", 1);
  const auto gr_zt = gr.zt, gr_zm = gr.zm;
  std::vector<int> k_lowest_zt(inputfields.stat_files.size()),
      k_highest_zt(inputfields.stat_files.size()), k_lowest_zm(inputfields.stat_files.size()),
      k_highest_zm(inputfields.stat_files.size());
  std::vector<input_field> SAM_variables, coamps_variables, RAMS_variables;
  iirr = hm_metadata.iirr;
  iiNr = hm_metadata.iiNr;
  iirs = hm_metadata.iirs;
  iiri = hm_metadata.iiri;
  iirg = hm_metadata.iirg;
  iiNi = hm_metadata.iiNi;
  iiNg = hm_metadata.iiNg;
  iiNs = hm_metadata.iiNs;

  switch (inputfields.stats_input_type) {

  //-------------------------------------
  // CLUBB stats data
  //-------------------------------------
  case clubb_input_type: {

    //  Thermo grid - zt file

    // Initialize l_fatal_error for case clubb_input_type
    l_fatal_error = false;

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_um, inputfields.stat_files[(clubb_zt - 1)], "um", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), um, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_vm, inputfields.stat_files[(clubb_zt - 1)], "vm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), vm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rtm, inputfields.stat_files[(clubb_zt - 1)], "rtm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), rtm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_thlm, inputfields.stat_files[(clubb_zt - 1)], "thlm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), thlm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wp3, inputfields.stat_files[(clubb_zt - 1)], "wp3", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wp3, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_tau_zt,
                                    inputfields.stat_files[(clubb_zt - 1)], "tau_zt", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), tau_zt,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rrm, inputfields.stat_files[(clubb_zt - 1)], "rrm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_rrm) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iirr),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }
    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rsm, inputfields.stat_files[(clubb_zt - 1)], "rsm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_rsm) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iirs),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }
    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rim, inputfields.stat_files[(clubb_zt - 1)], "rim", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_rim) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iiri),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }
    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rgm, inputfields.stat_files[(clubb_zt - 1)], "rgm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_rgm) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iirg),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }
    l_fatal_error = (l_fatal_error || l_read_error);

    //--------------------------------------------------------
    // Added variables for clubb_restart
    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_p, inputfields.stat_files[(clubb_zt - 1)], "p_in_Pa", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), p_in_Pa, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_exner, inputfields.stat_files[(clubb_zt - 1)], "exner",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), exner, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_ug, inputfields.stat_files[(clubb_zt - 1)], "ug", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), ug, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_vg, inputfields.stat_files[(clubb_zt - 1)], "vg", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), vg, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rcm, inputfields.stat_files[(clubb_zt - 1)], "rcm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), rcm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wm_zt, inputfields.stat_files[(clubb_zt - 1)], "wm_zt",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wm_zt, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rho, inputfields.stat_files[(clubb_zt - 1)], "rho", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), rho, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rho_ds_zt,
                                    inputfields.stat_files[(clubb_zt - 1)], "rho_ds_zt", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), rho_ds_zt,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thv_ds_zt,
                                    inputfields.stat_files[(clubb_zt - 1)], "thv_ds_zt", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), thv_ds_zt,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Kh_zt, inputfields.stat_files[(clubb_zt - 1)], "Kh_zt",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), Kh_zt, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_thvm, inputfields.stat_files[(clubb_zt - 1)], "thvm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), thvm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thlm_forcing,
                                    inputfields.stat_files[(clubb_zt - 1)], "thlm_forcing", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), thlm_forcing,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rtm_forcing,
                                    inputfields.stat_files[(clubb_zt - 1)], "rtm_forcing", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), rtm_forcing,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Ncm, inputfields.stat_files[(clubb_zt - 1)], "Ncm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_Ncm) {
      Kokkos::deep_copy(Kokkos::subview(Ncm, Kokkos::pair<int, int>(0, gr.nzt)),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Nccnm, inputfields.stat_files[(clubb_zt - 1)], "Nccnm",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), Nccnm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Nim, inputfields.stat_files[(clubb_zt - 1)], "Nim", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_Nim) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iiNi),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_cloud_frac,
                                    inputfields.stat_files[(clubb_zt - 1)], "cloud_frac", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), cloud_frac,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Nrm, inputfields.stat_files[(clubb_zt - 1)], "Nrm", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)), l_read_error, inputfields);
    if (inputfields.l_input_Nrm) {
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::pair<int, int>(0, gr.nzt), iiNr),
                        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, gr.nzt)));
    }

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_sigma_sqd_w_zt,
                                    inputfields.stat_files[(clubb_zt - 1)], "sigma_sqd_w_zt",
                                    gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
                                    sigma_sqd_w_zt, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wp2thvp,
                                    inputfields.stat_files[(clubb_zt - 1)], "wp2thvp", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wp2thvp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wp2up, inputfields.stat_files[(clubb_zt - 1)], "wp2up",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wp2up, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wp2rtp,
                                    inputfields.stat_files[(clubb_zt - 1)], "wp2rtp", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wp2rtp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wp2thlp,
                                    inputfields.stat_files[(clubb_zt - 1)], "wp2thlp", gr.nzt,
                                    timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wp2thlp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wpup2, inputfields.stat_files[(clubb_zt - 1)], "wpup2",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wpup2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wpvp2, inputfields.stat_files[(clubb_zt - 1)], "wpvp2",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), wpvp2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_iss_frac,
                                    inputfields.stat_files[(clubb_zt - 1)], "ice_supersat_frac",
                                    gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
                                    ice_supersat_frac, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_radht, inputfields.stat_files[(clubb_zt - 1)], "radht",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL), radht, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    // PDF Parameters (needed for CLUBB restarts)
    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_w_1, inputfields.stat_files[(clubb_zt - 1)], "w_1", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.w_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_w_2, inputfields.stat_files[(clubb_zt - 1)], "w_2", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.w_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_w_1, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_w_1", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_w_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_w_2, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_w_2", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_w_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rt_1, inputfields.stat_files[(clubb_zt - 1)], "rt_1", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.rt_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rt_2, inputfields.stat_files[(clubb_zt - 1)], "rt_2", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.rt_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_rt_1, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_rt_1", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_rt_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_rt_2, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_rt_2", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_rt_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_thl_1, inputfields.stat_files[(clubb_zt - 1)], "thl_1",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.thl_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_thl_2, inputfields.stat_files[(clubb_zt - 1)], "thl_2",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.thl_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_thl_1, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_thl_1", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_thl_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_thl_2, inputfields.stat_files[(clubb_zt - 1)],
        "varnce_thl_2", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.varnce_thl_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_mixt_frac, inputfields.stat_files[(clubb_zt - 1)],
        "mixt_frac", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.mixt_frac, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_chi_1, inputfields.stat_files[(clubb_zt - 1)], "chi_1",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.chi_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_chi_2, inputfields.stat_files[(clubb_zt - 1)], "chi_2",
        gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.chi_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_stdev_chi_1, inputfields.stat_files[(clubb_zt - 1)],
        "stdev_chi_1", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.stdev_chi_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_stdev_chi_2, inputfields.stat_files[(clubb_zt - 1)],
        "stdev_chi_2", gr.nzt, timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.stdev_chi_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rc_1, inputfields.stat_files[(clubb_zt - 1)], "rc_1", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.rc_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rc_2, inputfields.stat_files[(clubb_zt - 1)], "rc_2", gr.nzt,
        timestep, Kokkos::subview(gr.zt, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params.rc_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    //--------------------------------------------------------

    // Read in the zm file

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wp2, inputfields.stat_files[(clubb_zm - 1)], "wp2", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wp2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wprtp, inputfields.stat_files[(clubb_zm - 1)], "wprtp",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wprtp, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wpthlp,
                                    inputfields.stat_files[(clubb_zm - 1)], "wpthlp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wpthlp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wpthvp,
                                    inputfields.stat_files[(clubb_zm - 1)], "wpthvp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wpthvp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rtpthvp,
                                    inputfields.stat_files[(clubb_zm - 1)], "rtpthvp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rtpthvp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thlpthvp,
                                    inputfields.stat_files[(clubb_zm - 1)], "thlpthvp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), thlpthvp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_rtp2, inputfields.stat_files[(clubb_zm - 1)], "rtp2", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rtp2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_thlp2, inputfields.stat_files[(clubb_zm - 1)], "thlp2",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), thlp2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rtpthlp,
                                    inputfields.stat_files[(clubb_zm - 1)], "rtpthlp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rtpthlp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_upwp, inputfields.stat_files[(clubb_zm - 1)], "upwp", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), upwp, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_vpwp, inputfields.stat_files[(clubb_zm - 1)], "vpwp", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), vpwp, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    //-----------------------------------------------------------
    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_em, inputfields.stat_files[(clubb_zm - 1)], "em", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), em, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rho_zm,
                                    inputfields.stat_files[(clubb_zm - 1)], "rho_zm", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rho_zm,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rho_ds_zm,
                                    inputfields.stat_files[(clubb_zm - 1)], "rho_ds_zm", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rho_ds_zm,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thv_ds_zm,
                                    inputfields.stat_files[(clubb_zm - 1)], "thv_ds_zm", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), thv_ds_zm,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_Kh_zm, inputfields.stat_files[(clubb_zm - 1)], "Kh_zm",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), Kh_zm, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_tau_zm,
                                    inputfields.stat_files[(clubb_zm - 1)], "tau_zm", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), tau_zm,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_up2, inputfields.stat_files[(clubb_zm - 1)], "up2", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), up2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_vp2, inputfields.stat_files[(clubb_zm - 1)], "vp2", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), vp2, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_wp4, inputfields.stat_files[(clubb_zm - 1)], "wp4", gr.nzm,
        timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wp4, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_uprcp, inputfields.stat_files[(clubb_zm - 1)], "uprcp",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), uprcp, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_vprcp, inputfields.stat_files[(clubb_zm - 1)], "vprcp",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), vprcp, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wp2up2,
                                    inputfields.stat_files[(clubb_zm - 1)], "wp2up2", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wp2up2,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wp2vp2,
                                    inputfields.stat_files[(clubb_zm - 1)], "wp2vp2", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wp2vp2,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_sigma_sqd_w,
                                    inputfields.stat_files[(clubb_zm - 1)], "sigma_sqd_w", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), sigma_sqd_w,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wprtp_forcing,
                                    inputfields.stat_files[(clubb_zm - 1)], "wprtp_forcing", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), wprtp_forcing,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_wpthlp_forcing,
                                    inputfields.stat_files[(clubb_zm - 1)], "wpthlp_forcing",
                                    gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
                                    wpthlp_forcing, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rtp2_forcing,
                                    inputfields.stat_files[(clubb_zm - 1)], "rtp2_forcing", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rtp2_forcing,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thlp2_forcing,
                                    inputfields.stat_files[(clubb_zm - 1)], "thlp2_forcing", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), thlp2_forcing,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rtpthlp_forcing,
                                    inputfields.stat_files[(clubb_zm - 1)], "rtpthlp_forcing",
                                    gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
                                    rtpthlp_forcing, l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_thlprcp,
                                    inputfields.stat_files[(clubb_zm - 1)], "thlprcp", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), thlprcp,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(workspace, inputfields.l_input_rc_coef_zm,
                                    inputfields.stat_files[(clubb_zm - 1)], "rc_coef_zm", gr.nzm,
                                    timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL), rc_coef_zm,
                                    l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    // PDF Parameters (needed for CLUBB restarts)
    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_w_1_zm, inputfields.stat_files[(clubb_zm - 1)], "w_1_zm",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params_zm.w_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_w_2_zm, inputfields.stat_files[(clubb_zm - 1)], "w_2_zm",
        gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params_zm.w_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_w_1_zm, inputfields.stat_files[(clubb_zm - 1)],
        "varnce_w_1_zm", gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params_zm.varnce_w_1, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_varnce_w_2_zm, inputfields.stat_files[(clubb_zm - 1)],
        "varnce_w_2_zm", gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params_zm.varnce_w_2, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_mixt_frac_zm, inputfields.stat_files[(clubb_zm - 1)],
        "mixt_frac_zm", gr.nzm, timestep, Kokkos::subview(gr.zm, 0, Kokkos::ALL),
        Kokkos::subview(pdf_params_zm.mixt_frac, 0, Kokkos::ALL), l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);

    //-----------------------------------------------------------

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_veg_T_in_K, inputfields.stat_files[(clubb_sfc - 1)],
        "veg_T_in_K", 1, timestep, zero_height, Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)),
        l_read_error, inputfields);

    l_fatal_error = (l_fatal_error || l_read_error);
    veg_T_in_K    = to_host_space(Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)))(0);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_deep_soil_T_in_K, inputfields.stat_files[(clubb_sfc - 1)],
        "deep_soil_T_in_K", 1, timestep, zero_height,
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)), l_read_error, inputfields);

    l_fatal_error    = (l_fatal_error || l_read_error);
    deep_soil_T_in_K = to_host_space(Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)))(0);

    get_clubb_variable_interpolated(
        workspace, inputfields.l_input_sfc_soil_T_in_K, inputfields.stat_files[(clubb_sfc - 1)],
        "sfc_soil_T_in_K", 1, timestep, zero_height,
        Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)), l_read_error, inputfields);

    l_fatal_error   = (l_fatal_error || l_read_error);
    sfc_soil_T_in_K = to_host_space(Kokkos::subview(tmp1, Kokkos::pair<int, int>(0, 1)))(0);

    if (l_fatal_error) {
      throw std::runtime_error("oops, get_grads_var failed in stat_fields_reader");
    }

    //---------------------------------------
    // COAMPS LES stats data
    //--------------------------------------

    break;
  }
  case coamps_input_type: {
    l_fatal_error = false;

    coamps_variables.resize(num_coamps);

    k = 0;

    coamps_variables[k].l_input_var      = inputfields.l_input_um;
    coamps_variables[k].input_name       = "um";
    coamps_variables[k].clubb_var        = um;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_vm;
    coamps_variables[k].input_name       = "vm";
    coamps_variables[k].clubb_var        = vm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rtm;
    coamps_variables[k].input_name       = "qtm";
    coamps_variables[k].clubb_var        = rtm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_thlm;
    coamps_variables[k].input_name       = "thlm";
    coamps_variables[k].clubb_var        = thlm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    // We obtain wp2 from stats_sw

    coamps_variables[k].l_input_var      = inputfields.l_input_wp3;
    coamps_variables[k].input_name       = "wp3";
    coamps_variables[k].clubb_var        = wp3;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_wprtp;
    coamps_variables[k].input_name       = "wpqtp";
    coamps_variables[k].clubb_var        = wprtp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_wpthlp;
    coamps_variables[k].input_name       = "wpthlp";
    coamps_variables[k].clubb_var        = wpthlp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rtp2;
    coamps_variables[k].input_name       = "qtp2";
    coamps_variables[k].clubb_var        = rtp2;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_thlp2;
    coamps_variables[k].input_name       = "thlp2";
    coamps_variables[k].clubb_var        = thlp2;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rtpthlp;
    coamps_variables[k].input_name       = "qtpthlp";
    coamps_variables[k].clubb_var        = rtpthlp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_ug;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "ug";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_vg;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "vg";

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rcm;
    coamps_variables[k].input_name       = "qcm";
    coamps_variables[k].clubb_var        = rcm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_wm_zt;
    coamps_variables[k].input_name       = "wlsm";
    coamps_variables[k].clubb_var        = wm_zt;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_exner;
    coamps_variables[k].input_name       = "ex0";
    coamps_variables[k].clubb_var        = exner;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    // Note that em (clubb) = em (coamps) + tke (coamps), so just store tke
    // in a temporary variable for now and add them together later.
    coamps_variables[k].l_input_var      = inputfields.l_input_em;
    coamps_variables[k].input_name       = "em";
    coamps_variables[k].clubb_var        = em;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    Kokkos::deep_copy(temp_tke, real_t(0.0));
    // Initialize temp to 0.0

    coamps_variables[k].l_input_var      = inputfields.l_input_em;
    coamps_variables[k].input_name       = "tke";
    coamps_variables[k].clubb_var        = temp_tke;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_p;
    coamps_variables[k].input_name       = "pm";
    coamps_variables[k].clubb_var        = p_in_Pa;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rho;
    coamps_variables[k].input_name       = "dn0";
    coamps_variables[k].clubb_var        = rho;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    // rho_zm is in stats_sw

    coamps_variables[k].l_input_var      = inputfields.l_input_Kh_zt;
    coamps_variables[k].input_name       = "kh";
    coamps_variables[k].clubb_var        = Kh_zt;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_Kh_zm;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "Kh_zm";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_tau_zt;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "tau_zt";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_tau_zm;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "tau_zm";

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_wpthvp;
    coamps_variables[k].input_name       = "wpthvp";
    coamps_variables[k].clubb_var        = wpthvp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_thl_1;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "thl_1";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_thl_2;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "thl_2";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_mixt_frac;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "mixt_frac";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_chi_1;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "chi_1";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_chi_2;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "chi_2";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_stdev_chi_1;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "stdev_chi_1";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_stdev_chi_2;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "stdev_chi_2";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_rc_1;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "rc_1";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_rc_2;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "rc_2";

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_thvm;
    coamps_variables[k].input_name       = "thvm";
    coamps_variables[k].clubb_var        = thvm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    Kokkos::deep_copy(temp_rrm, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_rrm) {
      if ((iirr < 0)) {
        std::cerr << (("Rain water mixing ratio cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_rrm;
        coamps_variables[k].input_name       = "qrm";
        coamps_variables[k].clubb_var        = temp_rrm;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    Kokkos::deep_copy(temp_Nrm, real_t(0.0));
    // Initialize to 0.0

    if (inputfields.l_input_Nrm) {
      if ((iiNr < 0)) {
        std::cerr << (("Rain droplet number conc. cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_Nrm;
        coamps_variables[k].input_name       = "nrm";
        coamps_variables[k].clubb_var        = temp_Nrm;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    Kokkos::deep_copy(temp_Ncm, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_Ncm) {
      if ((!l_predict_Nc)) {
        std::cerr << (("Cloud droplet number conc. cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_Ncm;
        coamps_variables[k].input_name       = "ncm";
        coamps_variables[k].clubb_var        = temp_Ncm;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    Kokkos::deep_copy(temp_rsm, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_rsm) {
      if ((iirs < 0)) {
        std::cerr << (("Snow mixing ratio cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_rsm;
        coamps_variables[k].input_name       = "qsm";
        coamps_variables[k].clubb_var        = temp_rsm;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    Kokkos::deep_copy(temp_rim, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_rim) {
      if ((iiri < 0)) {
        std::cerr << (("Ice mixing ratio cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_rim;
        coamps_variables[k].input_name       = "qim";
        coamps_variables[k].clubb_var        = temp_rim;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    Kokkos::deep_copy(temp_rgm, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_rgm) {
      if ((iirg < 0)) {
        std::cerr << (("Graupel mixing ratio cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_rgm;
        coamps_variables[k].input_name       = "qgm";
        coamps_variables[k].clubb_var        = temp_rgm;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    coamps_variables[k].l_input_var      = inputfields.l_input_Nccnm;
    coamps_variables[k].input_name       = "nccnm";
    coamps_variables[k].clubb_var        = Nccnm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zt";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    Kokkos::deep_copy(temp_Nim, real_t(0.0));
    // initialize to 0.0

    if (inputfields.l_input_Nim) {
      if ((iiNi < 0)) {
        std::cerr << (("Ice number conc. cannot be input with"
                       " microphys_scheme = ") +
                      microphys_scheme)
                  << "\n";
        l_fatal_error = true;
      } else {
        coamps_variables[k].l_input_var      = inputfields.l_input_Nim;
        coamps_variables[k].input_name       = "nim";
        coamps_variables[k].clubb_var        = temp_Nim;
        coamps_variables[k].adjustment       = real_t(1.0);
        coamps_variables[k].clubb_grid_type  = "zt";
        coamps_variables[k].input_file_index = coamps_sm;

        k = (k + 1);
      }
    }

    coamps_variables[k].l_input_var = inputfields.l_input_thlm_forcing;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "thlm_forcing";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_rtm_forcing;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "rtm_forcing";

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_up2;
    coamps_variables[k].input_name       = "up2";
    coamps_variables[k].clubb_var        = up2;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_vp2;
    coamps_variables[k].input_name       = "vp2";
    coamps_variables[k].clubb_var        = vp2;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sm;

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_sigma_sqd_w;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "sigma_sqd_w";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_sigma_sqd_w_zt;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "sigma_sqd_w_zt";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_cloud_frac;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "cloud_frac";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_veg_T_in_K;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "veg_T_in_K";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_deep_soil_T_in_K;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "deep_soil_T_in_K";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_sfc_soil_T_in_K;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "sfc_soil_T_in_K";

    k = (k + 1);

    coamps_variables[k].l_input_var = inputfields.l_input_radht;
    coamps_variables[k].input_name  = "none";
    coamps_variables[k].clubb_name  = "radht";

    k = (k + 1);

    // Note:  wpup_sgs and wpvp_sgs must be added to make the u'w' and v'w' terms
    //        as they are in CLUBB.
    coamps_variables[k].l_input_var      = inputfields.l_input_upwp;
    coamps_variables[k].input_name       = "wpup";
    coamps_variables[k].clubb_var        = upwp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    k = (k + 1);

    Kokkos::deep_copy(temp_wpup_sgs, real_t(0.0));
    // initialize to 0.0

    coamps_variables[k].l_input_var      = inputfields.l_input_upwp;
    coamps_variables[k].input_name       = "wpup_sgs";
    coamps_variables[k].clubb_var        = temp_wpup_sgs;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_vpwp;
    coamps_variables[k].input_name       = "wpvp";
    coamps_variables[k].clubb_var        = vpwp;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    k = (k + 1);

    Kokkos::deep_copy(temp_wpvp_sgs, real_t(0.0));
    // initialize to 0.0

    coamps_variables[k].l_input_var      = inputfields.l_input_vpwp;
    coamps_variables[k].input_name       = "wpvp_sgs";
    coamps_variables[k].clubb_var        = temp_wpvp_sgs;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_wp2;
    coamps_variables[k].input_name       = "wp2";
    coamps_variables[k].clubb_var        = wp2;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    k = (k + 1);

    coamps_variables[k].l_input_var      = inputfields.l_input_rho_zm;
    coamps_variables[k].input_name       = "dn0";
    coamps_variables[k].clubb_var        = rho_zm;
    coamps_variables[k].adjustment       = real_t(1.0);
    coamps_variables[k].clubb_grid_type  = "zm";
    coamps_variables[k].input_file_index = coamps_sw;

    // Initialize l_read_error for case ( "les" )
    l_read_error = false;

    if (l_fatal_error) {
      throw std::runtime_error("oops, bad input for inputfields");
    }

    get_input_variables_interp(workspace, gr, ((k + 1)), coamps_variables, timestep, k_lowest_zt,
                               k_highest_zt, k_lowest_zm, k_highest_zm, l_fatal_error, inputfields);

    if (inputfields.l_input_rrm) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iirr),
          Kokkos::subview(temp_rrm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_Nrm) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iiNr),
          Kokkos::subview(temp_Nrm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_Ncm) {
      Kokkos::deep_copy(
          Kokkos::subview(
              Ncm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1))),
          Kokkos::subview(temp_Ncm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_rsm) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iirs),
          Kokkos::subview(temp_rsm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_rim) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iiri),
          Kokkos::subview(temp_rim, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_rgm) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iirg),
          Kokkos::subview(temp_rgm, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_Nim) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet,
              Kokkos::pair<int, int>(k_lowest_zt[coamps_sm], (k_highest_zt[coamps_sm] + 1)), iiNi),
          Kokkos::subview(temp_Nim, Kokkos::pair<int, int>(k_lowest_zt[coamps_sm],
                                                           (k_highest_zt[coamps_sm] + 1))));
    }

    if (inputfields.l_input_wprtp) {
      // When this is a standard scenario, where CLUBB momentum level 2 is the
      // first momentum level above the lowest LES level, use the values of
      // wprtp at momentum levels 3 and 2 to find the value at momentum level 1
      // through the use of a linear extension.  It should be pointed out that
      // the boundary flux is usually solved in LES or CLUBB via a subroutine
      // like surface_varnce.
      if ((k_lowest_zm[coamps_sm] == 1)) {
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_1", 1, KOKKOS_LAMBDA(int) {
              wprtp(0) =
                  lin_ext_zm_bottom(wprtp(2), wprtp(1), gr_zm(0, 2), gr_zm(0, 1), gr_zm(0, 0));
            });
      }
    }

    if (inputfields.l_input_wpthlp) {
      // When this is a standard scenario, where CLUBB momentum level 2 is the
      // first momentum level above the lowest LES level, use the values of
      // wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
      // through the use of a linear extension.  It should be pointed out that
      // the boundary flux is usually solved in LES or CLUBB via a subroutine
      // like surface_varnce.
      if ((k_lowest_zm[coamps_sm] == 1)) {
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_2", 1, KOKKOS_LAMBDA(int) {
              wpthlp(0) =
                  lin_ext_zm_bottom(wpthlp(2), wpthlp(1), gr_zm(0, 2), gr_zm(0, 1), gr_zm(0, 0));
            });
      }
    }

    if (inputfields.l_input_rtp2) {
      // When this is a standard scenario, where CLUBB momentum level 2 is the
      // first momentum level above the lowest LES level, set the value of rtp2
      // at momentum level 1 to the value at momentum level 2.
      // Using a linear extension here resulted in negatives.
      if ((k_lowest_zm[coamps_sm] == 1)) {
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_3", 1, KOKKOS_LAMBDA(int) { rtp2(0) = rtp2(1); });
      }

      Kokkos::parallel_for(
          "stat_fields_reader::clip_4", rtp2.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((rtp2(level_index) < (rt_tol * rt_tol))) {
              rtp2(level_index) = (rt_tol * rt_tol);
            }
          });
    }

    if (inputfields.l_input_thlp2) {
      // When this is a standard scenario, where CLUBB momentum level 2 is the
      // first momentum level above the lowest LES level, set the value of thlp2
      // at momentum level 1 to the value at momentum level 2.
      // Using a linear extension here resulted in negatives.
      if ((k_lowest_zm[coamps_sm] == 1)) {
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_5", 1, KOKKOS_LAMBDA(int) { thlp2(0) = thlp2(1); });
      }

      Kokkos::parallel_for(
          "stat_fields_reader::clip_6", thlp2.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((thlp2(level_index) < (thl_tol * thl_tol))) {
              thlp2(level_index) = (thl_tol * thl_tol);
            }
          });
    }

    if (inputfields.l_input_rtpthlp) {
      // When this is a standard scenario, where CLUBB momentum level 2 is the
      // first momentum level above the lowest LES level, use the values of
      // rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
      // through the use of a linear extension.  It should be pointed out that
      // the boundary flux is usually solved in LES or CLUBB via a subroutine
      // like surface_varnce.
      if ((k_lowest_zm[coamps_sm] == 1)) {
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_7", 1, KOKKOS_LAMBDA(int) {
              rtpthlp(0) =
                  lin_ext_zm_bottom(rtpthlp(2), rtpthlp(1), gr_zm(0, 2), gr_zm(0, 1), gr_zm(0, 0));
            });
      }
    }

    if (inputfields.l_input_em) {
      const int range_low_8 = k_lowest_zm[coamps_sm], range_high_8 = k_highest_zm[coamps_sm];
      Kokkos::parallel_for(
          "stat_fields_reader::cells_8", (range_high_8 - range_low_8 + 1),
          KOKKOS_LAMBDA(int level_index) {
            const int k = (level_index + range_low_8);
            em(k)       = (em(k) + temp_tke(k));
          });

      Kokkos::parallel_for(
          "stat_fields_reader::clip_9", em.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((em(level_index) < em_min)) {
              em(level_index) = em_min;
            }
          });
    }

    if (inputfields.l_input_up2) {
      // Clip up2 to be no smaller than w_tol_sqd
      Kokkos::parallel_for(
          "stat_fields_reader::clip_10", up2.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((up2(level_index) < w_tol_sqd)) {
              up2(level_index) = w_tol_sqd;
            }
          });
    }

    if (inputfields.l_input_vp2) {
      Kokkos::parallel_for(
          "stat_fields_reader::clip_11", vp2.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((vp2(level_index) < w_tol_sqd)) {
              vp2(level_index) = w_tol_sqd;
            }
          });
    }

    if (inputfields.l_input_upwp) {
      const int range_low_12 = k_lowest_zm[coamps_sw], range_high_12 = k_highest_zm[coamps_sw];
      Kokkos::parallel_for(
          "stat_fields_reader::cells_12", (range_high_12 - range_low_12 + 1),
          KOKKOS_LAMBDA(int level_index) {
            const int k = (level_index + range_low_12);
            upwp(k)     = temp_wpup_sgs(k);
          });
    }

    if (inputfields.l_input_vpwp) {
      const int range_low_13 = k_lowest_zm[coamps_sw], range_high_13 = k_highest_zm[coamps_sw];
      Kokkos::parallel_for(
          "stat_fields_reader::cells_13", (range_high_13 - range_low_13 + 1),
          KOKKOS_LAMBDA(int level_index) {
            const int k = (level_index + range_low_13);
            vpwp(k)     = temp_wpvp_sgs(k);
          });
    }

    if (inputfields.l_input_wp2) {
      Kokkos::parallel_for(
          "stat_fields_reader::clip_14", wp2.extent(0), KOKKOS_LAMBDA(int level_index) {
            if ((wp2(level_index) < w_tol_sqd)) {
              wp2(level_index) = w_tol_sqd;
            }
          });
    }

    coamps_variables.clear();

    //-------------------------------------------
    // SAM stats data
    //-------------------------------------------
    break;
  }
  case sam_input_type: {

    SAM_variables.resize(num_sam_inputfields);

    l_fatal_error = false;

    k = 0;

    SAM_variables[k].l_input_var      = inputfields.l_input_um;
    SAM_variables[k].input_name       = "U";
    SAM_variables[k].clubb_var        = um;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_vm;
    SAM_variables[k].input_name       = "V";
    SAM_variables[k].clubb_var        = vm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rtm;
    SAM_variables[k].input_name  = "RTM";
    SAM_variables[k].clubb_var   = rtm;
    SAM_variables[k].adjustment  = real_t(1.0);
    //1.0e-3_core_rknd
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_thlm;
    SAM_variables[k].input_name       = "THLM";
    SAM_variables[k].clubb_var        = thlm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_wp2;
    SAM_variables[k].input_name       = "WP2";
    SAM_variables[k].clubb_var        = wp2;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_rho;
    SAM_variables[k].input_name       = "RHO";
    SAM_variables[k].clubb_var        = rho;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    // Note that this needs to be adjusted by 1/(RHO * Lv)
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_wprtp;
    SAM_variables[k].input_name       = "WPRTP";
    SAM_variables[k].clubb_var        = wprtp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    // Note that this needs to be adjusted by 1/(RHO * Cp)
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_wpthlp;
    SAM_variables[k].input_name       = "WPTHLP";
    SAM_variables[k].clubb_var        = wpthlp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_wp3;
    SAM_variables[k].input_name       = "W3";
    SAM_variables[k].clubb_var        = wp3;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rtp2;
    SAM_variables[k].input_name  = "RTP2";
    SAM_variables[k].clubb_var   = rtp2;
    SAM_variables[k].adjustment  = real_t(1.0);
    //1.0e-6_core_rknd
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_rtp3;
    SAM_variables[k].input_name       = "QTO3";
    SAM_variables[k].clubb_var        = rtp3;
    SAM_variables[k].adjustment       = real_t(1.0e-9);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_thlp2;
    SAM_variables[k].input_name       = "THLP2";
    SAM_variables[k].clubb_var        = thlp2;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_thlp3;
    SAM_variables[k].input_name       = "THEL3";
    SAM_variables[k].clubb_var        = thlp3;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_upwp;
    SAM_variables[k].input_name       = "UW";
    SAM_variables[k].clubb_var        = upwp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_vpwp;
    SAM_variables[k].input_name       = "VW";
    SAM_variables[k].clubb_var        = vpwp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_ug;
    SAM_variables[k].input_name       = "UOBS";
    SAM_variables[k].clubb_var        = ug;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_vg;
    SAM_variables[k].input_name       = "VOBS";
    SAM_variables[k].clubb_var        = vg;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_rcm;
    SAM_variables[k].input_name       = "QCL";
    SAM_variables[k].clubb_var        = rcm;
    SAM_variables[k].adjustment       = real_t(1.0e-3);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_wm_zt;
    SAM_variables[k].input_name       = "WOBS";
    SAM_variables[k].clubb_var        = wm_zt;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_exner;
    SAM_variables[k].clubb_name  = "exner";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_em;
    SAM_variables[k].input_name       = "TKE";
    SAM_variables[k].clubb_var        = em;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_p;
    SAM_variables[k].input_name       = "PRES";
    SAM_variables[k].clubb_var        = p_in_Pa;
    SAM_variables[k].adjustment       = pascal_per_mb;
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rho_zm;
    SAM_variables[k].clubb_name  = "rho_zm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rho_ds_zm;
    SAM_variables[k].clubb_name  = "rho_ds_zm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rho_ds_zt;
    SAM_variables[k].clubb_name  = "rho_ds_zt";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_thv_ds_zm;
    SAM_variables[k].clubb_name  = "thv_ds_zm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_thv_ds_zt;
    SAM_variables[k].clubb_name  = "thv_ds_zt";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_Kh_zt;
    SAM_variables[k].input_name       = "TKH";
    SAM_variables[k].clubb_var        = Kh_zt;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_Kh_zm;
    SAM_variables[k].clubb_name  = "Kh_zm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_tau_zm;
    SAM_variables[k].clubb_name  = "tau_zm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_tau_zt;
    SAM_variables[k].clubb_name  = "tau_zt";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_thl_1;
    SAM_variables[k].clubb_name  = "thl_1";
    SAM_variables[k].input_name  = "none";

    k                            = (k + 1);
    SAM_variables[k].l_input_var = inputfields.l_input_thl_2;
    SAM_variables[k].clubb_name  = "thl_2";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_mixt_frac;
    SAM_variables[k].clubb_name  = "mixt_frac";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_chi_1;
    SAM_variables[k].clubb_name  = "chi_1";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_chi_2;
    SAM_variables[k].clubb_name  = "chi_2";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_stdev_chi_1;
    SAM_variables[k].clubb_name  = "stdev_chi_1";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_stdev_chi_2;
    SAM_variables[k].clubb_name  = "stdev_chi_2";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rc_1;
    SAM_variables[k].clubb_name  = "rc_1";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rc_2;
    SAM_variables[k].clubb_name  = "rc_2";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_thvm;
    SAM_variables[k].input_name       = "THETAV";
    SAM_variables[k].clubb_var        = thvm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_rrm, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_rrm;
    SAM_variables[k].input_name       = "QPL";
    SAM_variables[k].clubb_var        = temp_rrm;
    SAM_variables[k].adjustment       = (real_t(1.0) / g_per_kg);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_Nrm, real_t(0.0));
    // initialize to 0.0

    // Note that this needs to be adjusted by 1e6/RHO
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var = inputfields.l_input_Nrm;
    SAM_variables[k].input_name  = "NR";
    // SAM Morrison microphysics
    //SAM_variables(k)%input_name = "CONP" ! SAM KK microphysics
    SAM_variables[k].clubb_var        = temp_Nrm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);
    Kokkos::deep_copy(temp_Ncm, real_t(0.0));
    // initialize to 0.0

    // Note that this needs to be adjusted by 1e6/RHO
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_Ncm;
    SAM_variables[k].input_name       = "NC";
    SAM_variables[k].clubb_var        = temp_Ncm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_rsm, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_rsm;
    SAM_variables[k].input_name       = "QS";
    SAM_variables[k].clubb_var        = temp_rsm;
    SAM_variables[k].adjustment       = (real_t(1.0) / g_per_kg);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_rim, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_rim;
    SAM_variables[k].input_name       = "QI";
    SAM_variables[k].clubb_var        = temp_rim;
    SAM_variables[k].adjustment       = (real_t(1.0) / g_per_kg);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_rgm, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_rgm;
    SAM_variables[k].input_name       = "QG";
    SAM_variables[k].clubb_var        = temp_rgm;
    SAM_variables[k].adjustment       = (real_t(1.0) / g_per_kg);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_Nccnm;
    SAM_variables[k].clubb_name  = "Nccnm";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    Kokkos::deep_copy(temp_Nim, real_t(0.0));
    // initialize to 0.0

    // Note that this needs to be adjusted by 1e6/RHO
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_Nim;
    SAM_variables[k].input_name       = "NI";
    SAM_variables[k].clubb_var        = temp_Nim;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_Nsm, real_t(0.0));
    // initialize to 0.0

    // Note that this needs to be adjusted by 1e6/RHO
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_Nsm;
    SAM_variables[k].input_name       = "NS";
    SAM_variables[k].clubb_var        = temp_Nsm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_Ngm, real_t(0.0));
    // initialize to 0.0

    // Note that this needs to be adjusted by 1e6/RHO
    // This will need to be adjusted outside of get_sam_variable_interpolated
    SAM_variables[k].l_input_var      = inputfields.l_input_Ngm;
    SAM_variables[k].input_name       = "NG";
    SAM_variables[k].clubb_var        = temp_Ngm;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_rrp2, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var = inputfields.l_input_rrp2;
    SAM_variables[k].input_name  = "rrp2";
    // SAM Morrison microphysics
    //SAM_variables(k)%input_name = "RRP2" ! SAM KK microphysics
    SAM_variables[k].clubb_var        = temp_rrp2;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_Nrp2, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var = inputfields.l_input_Nrp2;
    SAM_variables[k].input_name  = "Nrp2";
    // SAM Morrison microphysics
    //SAM_variables(k)%input_name = "NRP2" ! SAM KK microphysics
    SAM_variables[k].clubb_var        = temp_Nrp2;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_wprrp, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_wprrp;
    SAM_variables[k].input_name       = "WPRRP";
    SAM_variables[k].clubb_var        = temp_wprrp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    Kokkos::deep_copy(temp_wpNrp, real_t(0.0));
    // initialize to 0.0

    SAM_variables[k].l_input_var      = inputfields.l_input_wpNrp;
    SAM_variables[k].input_name       = "WPNRP";
    SAM_variables[k].clubb_var        = temp_wpNrp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_thlm_forcing;
    SAM_variables[k].clubb_name  = "thlm_forcing";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_rtm_forcing;
    SAM_variables[k].clubb_name  = "rtm_forcing";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    // Note that U2 is on the thermodynamic grid in SAM while up2 is
    // on the momentum grid in CLUBB. This will have to be interpolated.
    SAM_variables[k].l_input_var      = inputfields.l_input_up2;
    SAM_variables[k].input_name       = "U2";
    SAM_variables[k].clubb_var        = up2_zt;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    // Note that V2 is on the thermodynamic grid in SAM while vp2 is
    // on the momentum grid in CLUBB. This will have to be interpolated.
    SAM_variables[k].l_input_var      = inputfields.l_input_vp2;
    SAM_variables[k].input_name       = "V2";
    SAM_variables[k].clubb_var        = vp2_zt;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_sigma_sqd_w;
    SAM_variables[k].clubb_name  = "sigma_sqd_w";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_cloud_frac;
    SAM_variables[k].input_name       = "CLD";
    SAM_variables[k].clubb_var        = cloud_frac;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_sigma_sqd_w_zt;
    SAM_variables[k].clubb_name  = "sigma_sqd_w_zt";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_veg_T_in_K;
    SAM_variables[k].clubb_name  = "veg_T_in_K";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_deep_soil_T_in_K;
    SAM_variables[k].clubb_name  = "deep_soil_T_in_K";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var = inputfields.l_input_sfc_soil_T_in_K;
    SAM_variables[k].clubb_name  = "sfc_soil_T_in_K";
    SAM_variables[k].input_name  = "none";

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_radht;
    SAM_variables[k].clubb_name       = "radht";
    SAM_variables[k].input_name       = "RADQR";
    SAM_variables[k].clubb_var        = radht;
    SAM_variables[k].adjustment       = (real_t(1.0) / sec_per_day);
    SAM_variables[k].clubb_grid_type  = "zt";
    SAM_variables[k].input_file_index = sam_file;

    k = (k + 1);

    SAM_variables[k].l_input_var      = inputfields.l_input_rtpthlp;
    SAM_variables[k].clubb_name       = "rtpthlp";
    SAM_variables[k].input_name       = "RTPTHLP";
    SAM_variables[k].clubb_var        = rtpthlp;
    SAM_variables[k].adjustment       = real_t(1.0);
    SAM_variables[k].clubb_grid_type  = "zm";
    SAM_variables[k].input_file_index = sam_file;

    get_input_variables_interp(workspace, gr, ((k + 1)), SAM_variables, timestep, k_lowest_zt,
                               k_highest_zt, k_lowest_zm, k_highest_zm, l_fatal_error, inputfields);

    if (l_fatal_error) {
      throw std::runtime_error("Error reading SAM input fields");
    }

    // Perform non-constant adjustments
    const bool flag_l_input_Ncm = inputfields.l_input_Ncm;
    const bool flag_l_input_Ngm = inputfields.l_input_Ngm;
    const bool flag_l_input_Nim = inputfields.l_input_Nim;
    const bool flag_l_input_Nrm = inputfields.l_input_Nrm;
    const bool flag_l_input_Nsm = inputfields.l_input_Nsm;
    const int range_low_15 = k_lowest_zt[sam_file], range_high_15 = k_highest_zt[sam_file];
    Kokkos::parallel_for(
        "stat_fields_reader::cells_15", (range_high_15 - range_low_15 + 1),
        KOKKOS_LAMBDA(int level_index) {
          const int k = (level_index + range_low_15);
          if (flag_l_input_Nrm) {
            // Nrm = NR * (1.0e6 / RHO)
            temp_Nrm(k) = (temp_Nrm(k) * (real_t(1.0e6) / rho(k)));
          }

          if (flag_l_input_Ncm) {
            // Ncm = NC * (1.0e6 / RHO)
            temp_Ncm(k) = (temp_Ncm(k) * (real_t(1.0e6) / rho(k)));
          }

          if (flag_l_input_Nim) {
            // Nim = NI * (1.0e6 / RHO)
            temp_Nim(k) = (temp_Nim(k) * (real_t(1.0e7) / rho(k)));
          }

          if (flag_l_input_Nsm) {
            // Nsm = NI * (1.0e6 / RHO)
            temp_Nsm(k) = (temp_Nsm(k) * (real_t(1.0e7) / rho(k)));
          }

          if (flag_l_input_Ngm) {
            // Ngm = NG * (1.0e6 / RHO)
            temp_Ngm(k) = (temp_Ngm(k) * (real_t(1.0e7) / rho(k)));
          }
        });

    for (k = k_lowest_zm[sam_file]; k <= k_highest_zm[sam_file]; ++k) {
      //if(l_input_wprtp) then
      // wprtp = QTFLUX / (RHO*Lv)
      //wprtp(k) = wprtp(k) / (rho(k)*Lv) weberjk, should not need to convert with new fields
      //end if

      //if(l_input_wpthlp) then
      // wpthlp = TLFLUX / (RHO*Cp)
      //wpthlp(k) = wpthlp(k) / (rho(k)*Cp) weberjk, should not need to convert with new fields
      //end if

      // Use linear interpolation to convert up2 and vp2 from the thermodynamic
      // grid to the momentum grid. See CLUBB Equations 170 and 171.
      if (inputfields.l_input_up2) {
        const real_t scalar_value = zt2zm_api(gr, up2_zt, k);
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_16", 1, KOKKOS_LAMBDA(int) { up2(k) = scalar_value; });
      }

      if (inputfields.l_input_vp2) {
        const real_t scalar_value = zt2zm_api(gr, vp2_zt, k);
        Kokkos::parallel_for(
            "stat_fields_reader::scalar_17", 1, KOKKOS_LAMBDA(int) { vp2(k) = scalar_value; });
      }
    }

    // Add hydrometeor variables.
    if ((inputfields.l_input_Nrm && (iiNr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iiNr),
          Kokkos::subview(temp_Nrm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_Ncm && l_predict_Nc)) {
      Kokkos::deep_copy(
          Kokkos::subview(
              Ncm, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1))),
          Kokkos::subview(temp_Ncm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_Nim && (iiNi > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iiNi),
          Kokkos::subview(temp_Nim, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_rrm && (iirr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iirr),
          Kokkos::subview(temp_rrm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_Ngm && (iiNg > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iiNg),
          Kokkos::subview(temp_Ngm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_rgm && (iirg > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iirg),
          Kokkos::subview(temp_rgm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_rim && (iiri > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iiri),
          Kokkos::subview(temp_rim, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_Nsm && (iiNs > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iiNs),
          Kokkos::subview(temp_Nsm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_rsm && (iirs > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydromet, Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)),
              iirs),
          Kokkos::subview(temp_rsm, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                           (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_rrp2 && (iirr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydrometp2,
              Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)), iirr),
          Kokkos::subview(temp_rrp2, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                            (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_Nrp2 && (iiNr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              hydrometp2,
              Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)), iiNr),
          Kokkos::subview(temp_Nrp2, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                            (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_wprrp && (iirr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              wphydrometp,
              Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)), iirr),
          Kokkos::subview(temp_wprrp, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                             (k_highest_zt[sam_file] + 1))));
    }

    if ((inputfields.l_input_wpNrp && (iiNr > -1))) {
      Kokkos::deep_copy(
          Kokkos::subview(
              wphydrometp,
              Kokkos::pair<int, int>(k_lowest_zt[sam_file], (k_highest_zt[sam_file] + 1)), iiNr),
          Kokkos::subview(temp_wpNrp, Kokkos::pair<int, int>(k_lowest_zt[sam_file],
                                                             (k_highest_zt[sam_file] + 1))));
    }

    if (l_fatal_error) {
      throw std::runtime_error("Failed to read inputfields for SAM.");
    }

    SAM_variables.clear();

    //---------------------------------------------------
    // RAMS LES stats data
    //---------------------------------------------------
    break;
  }
  case rams_input_type: {
    RAMS_variables.resize(num_rams_inputfields);

    l_fatal_error = false;

    k = 0;

    RAMS_variables[k].l_input_var      = inputfields.l_input_wp3;
    RAMS_variables[k].input_name       = "wp3";
    RAMS_variables[k].clubb_var        = wp3;
    RAMS_variables[k].adjustment       = real_t(1.0);
    RAMS_variables[k].clubb_grid_type  = "zt";
    RAMS_variables[k].input_file_index = rams_file;

    k = (k + 1);

    RAMS_variables[k].l_input_var      = inputfields.l_input_wp2;
    RAMS_variables[k].input_name       = "wp2";
    RAMS_variables[k].clubb_var        = wp2;
    RAMS_variables[k].adjustment       = real_t(1.0);
    RAMS_variables[k].clubb_grid_type  = "zm";
    RAMS_variables[k].input_file_index = rams_file;

    k = (k + 1);

    RAMS_variables[k].l_input_var = inputfields.l_input_up2;
    RAMS_variables[k].input_name  = "up2";

    RAMS_variables[k].clubb_var        = up2;
    RAMS_variables[k].adjustment       = real_t(1.0);
    RAMS_variables[k].clubb_grid_type  = "zm";
    RAMS_variables[k].input_file_index = rams_file;

    k = (k + 1);

    RAMS_variables[k].l_input_var      = inputfields.l_input_vp2;
    RAMS_variables[k].input_name       = "vp2";
    RAMS_variables[k].clubb_var        = vp2;
    RAMS_variables[k].adjustment       = real_t(1.0);
    RAMS_variables[k].clubb_grid_type  = "zm";
    RAMS_variables[k].input_file_index = rams_file;

    get_input_variables_interp(workspace, gr, ((k + 1)), RAMS_variables, timestep, k_lowest_zt,
                               k_highest_zt, k_lowest_zm, k_highest_zm, l_fatal_error, inputfields);

    if (l_fatal_error) {
      throw std::runtime_error("Error reading RAMS LES input fields");
    }

    RAMS_variables.clear();

    break;
  }
  }

  // Clipping on the variance of u, v and w
  if (inputfields.l_input_wp2) {
    Kokkos::parallel_for(
        "stat_fields_reader::clip_18", wp2.extent(0), KOKKOS_LAMBDA(int level_index) {
          if ((wp2(level_index) < w_tol_sqd)) {
            wp2(level_index) = w_tol_sqd;
          }
        });
  }

  if (inputfields.l_input_up2) {
    Kokkos::parallel_for(
        "stat_fields_reader::clip_19", up2.extent(0), KOKKOS_LAMBDA(int level_index) {
          if ((up2(level_index) < w_tol_sqd)) {
            up2(level_index) = w_tol_sqd;
          }
        });
  }

  if (inputfields.l_input_vp2) {
    Kokkos::parallel_for(
        "stat_fields_reader::clip_20", vp2.extent(0), KOKKOS_LAMBDA(int level_index) {
          if ((vp2(level_index) < w_tol_sqd)) {
            vp2(level_index) = w_tol_sqd;
          }
        });
  }
}

void compute_timestep(ScratchWorkspace& workspace, const std::string& filename,
                      bool l_restart,                       // In
                      double time,                          // In
                      int& nearest_timestep,                // Out
                      const inputfields_state& inputfields) // Native module In
{
  // Given a time [s], determine the closest NetCDF output record.
  // Keep Native NINT/no+1 record convention and restart minute guards; JAX's
  // coordinate-nearest-record adaptation is deliberately not applied here.
  stat_file fread_var;
  bool l_error = false;
  if (!l_netcdf_file(filename)) {
    std::cerr << "compute_timestep: Expected netCDF file, got: " << filename << '\n';
    throw std::runtime_error("Fatal error: Non-netCDF input.");
  }
  if (inputfields.stats_input_type == sam_input_type)
    open_netcdf_read(workspace, "U", filename, fread_var, l_error, inputfields);
  else
    open_netcdf_read(workspace, "thlm", filename, fread_var, l_error, inputfields);
  if (l_error) {
    std::cerr << "Error reading file " << filename << '\n';
    throw std::runtime_error("Fatal error");
  }
  // Eric Raut June2013; Joshua Fasching March2008.
  const real_t delta_time = real_t(time - fread_var.time);
  if (l_restart) {
    std::cout << "Initial time of netCDF reference file [seconds since midnight]: "
              << fread_var.time << '\n';
    std::cout << "Model restart time [s]: " << time << '\n';
    std::cout << "Elapsed time between initial time of ref file and restart time [s]: "
              << delta_time << '\n';
    std::cout << "netCDF file output time interval [s]: " << fread_var.dtwrite << '\n';
    if (Kokkos::fmod(delta_time, fread_var.dtwrite) > real_t(1e-8) ||
        Kokkos::fmod(delta_time, fread_var.dtwrite) < -real_t(1e-8)) {
      std::cerr << "Error: Elapsed time is not a multiple of the reference netCDF output time "
                   "interval.\n";
      std::cerr << "Elapsed time [s] = " << delta_time << '\n';
      std::cerr << "netCDF output time interval = " << fread_var.dtwrite << '\n';
      throw std::runtime_error("STOP in compute_timestep");
    }
    if (Kokkos::fmod(delta_time, sec_per_min) > real_t(1e-8) ||
        Kokkos::fmod(delta_time, sec_per_min) < -real_t(1e-8)) {
      std::cerr << "Error: Elapsed time is not a multiple of one minute.\n";
      std::cerr << "Elapsed time [s] = " << delta_time << '\n';
      throw std::runtime_error("STOP in compute_timestep");
    }
  }
  // NINT rounds halfway away from zero; C++ round retains that rule.
  nearest_timestep = int(Kokkos::round(delta_time / fread_var.dtwrite));
  if (l_restart) {
    std::cout << "Elapsed time between initial time of ref file and restart time "
                 "rounded to nearest minute: "
              << nearest_timestep << '\n';
    // Joshua Fasching March2008: keep the source diagnostic expression.
    std::cout << "Nearest netCDF output time iteration [ ]: "
              << int(Kokkos::round(real_t(nearest_timestep) / fread_var.dtwrite / sec_per_min)) - 1
              << '\n';
  }
  close_netcdf_read(fread_var);
}

void get_clubb_variable_interpolated(ScratchWorkspace& workspace, bool l_input_var,
                                     const std::string& filename,
                                     const std::string& varname, // In
                                     int vardim, int timestep,
                                     ConstInputProfile clubb_heights,      // In
                                     InputProfile variable_interpolated,   // InOut
                                     bool& l_read_error,                   // Out
                                     const inputfields_state& inputfields) // Native module In
{
  // Obtain a CLUBB variable profile from NetCDF and interpolate if needed.
  constexpr int npower = 1; // Not applied.
  if (l_input_var) {
    real_t first_height;
    Kokkos::deep_copy(first_height, Kokkos::subview(clubb_heights, 0));
    const bool l_spec_bound_cond = first_height < real_t(0.0);
    stat_file_average(workspace, filename, vardim, timestep, timestep, clubb_heights, varname,
                      npower, l_spec_bound_cond, l_read_error, inputfields, variable_interpolated);
  } else {
    l_read_error = false;
  }
}

static void
get_input_variables_interp(ScratchWorkspace& workspace, const grid& gr, int nvars,
                           const std::vector<input_field>& variables, int timestep,       // In
                           std::vector<int>& k_lowest_zt, std::vector<int>& k_highest_zt, // Out
                           std::vector<int>& k_lowest_zm, std::vector<int>& k_highest_zm,
                           bool& l_error,                        // Out
                           const inputfields_state& inputfields) // Native module In
{
  // Obtain COAMPS/SAM/RAMS profiles and interpolate to the CLUBB grid.
  // Metadata remains host-owned; scientific array data stays in Kokkos Views.
  auto workspace_scope = workspace.scope();
  l_error              = false;
  std::vector<stat_file> fread_vars(inputfields.stat_files.size());
  bool l_internal_error;
  for (int file_index = 0; file_index < int(inputfields.stat_files.size()); ++file_index) {
    open_netcdf_read(workspace, "U", inputfields.stat_files[file_index], // In
                     fread_vars[file_index], l_internal_error,           // Out
                     inputfields);                                       // Native module In
    l_error = l_error || l_internal_error;
  }
  if (l_error) {
    std::cerr << "A fatal error occured while reading the input file.\n";
    throw std::runtime_error("Fatal error in input_fields.get_input_variables_interp");
  }
  int min_ia = 9998, max_iz = -2;
  const auto zt = Kokkos::subview(gr.zt, 0, Kokkos::ALL);
  const auto zm = Kokkos::subview(gr.zm, 0, Kokkos::ALL);
  for (int file_index = 0; file_index < int(inputfields.stat_files.size()); ++file_index) {
    CLUBB_levels_within_LES_domain(workspace, fread_vars[file_index], zt, // In
                                   k_lowest_zt[file_index],
                                   k_highest_zt[file_index]);             // Out
    CLUBB_levels_within_LES_domain(workspace, fread_vars[file_index], zm, // In
                                   k_lowest_zm[file_index],
                                   k_highest_zm[file_index]); // Out
    min_ia = std::min(min_ia, fread_vars[file_index].ia);
    max_iz = std::max(max_iz, fread_vars[file_index].iz);
  }
  // Preserve Native (file, vertical) order and explicit min_ia lower bound.
  Real2D LES_tmp = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                "get_input_variables_interp::LES_tmp",
                                                inputfields.stat_files.size(), max_iz - min_ia + 1);
  for (int i = 0; i < nvars; ++i) {
    const input_field current_var = variables[i];
    if (current_var.input_name == "none" && current_var.l_input_var) {
      std::cerr << "The variable " << current_var.clubb_name << " is not setup for input_type "
                << inputfields.input_type << '\n';
      l_error = true;
      continue;
    }
    if (current_var.l_input_var) {
      const int file_index        = current_var.input_file_index;
      const auto& fread_var       = fread_vars[file_index];
      const bool l_convert_to_MKS = false;
      auto profile =
          Kokkos::subview(LES_tmp, file_index,
                          Kokkos::pair<int, int>(fread_var.ia - min_ia, fread_var.iz - min_ia + 1));
      get_netcdf_var(workspace, fread_var, current_var.input_name, timestep, l_convert_to_MKS,
                     profile, l_internal_error);
      l_error = l_error || l_internal_error;
      // Native passes the whole row to a dummy declared ia:iz. Retain sequence
      // association by taking its first iz-ia+1 entries, not shifting the row.
      const auto LES_row = Kokkos::subview(
          LES_tmp, file_index, Kokkos::pair<int, int>(0, fread_var.iz - fread_var.ia + 1));
      if (current_var.clubb_grid_type == "zt")
        inputfields_interp_and_adjust(workspace, gr, current_var.clubb_grid_type, // In
                                      fread_var, current_var.adjustment, LES_row, // In
                                      k_lowest_zt[file_index],
                                      k_highest_zt[file_index], // In
                                      current_var.clubb_var);   // InOut
      else if (current_var.clubb_grid_type == "zm")
        inputfields_interp_and_adjust(workspace, gr, current_var.clubb_grid_type, // In
                                      fread_var, current_var.adjustment, LES_row, // In
                                      k_lowest_zm[file_index],
                                      k_highest_zm[file_index], // In
                                      current_var.clubb_var);   // InOut
    }
  }
  for (int file_index = 0; file_index < int(inputfields.stat_files.size()); ++file_index)
    close_netcdf_read(fread_vars[file_index]);
  LES_tmp = Real2D{};
}

template <class LESProfile>
static void inputfields_interp_and_adjust(ScratchWorkspace& workspace, const grid& gr,
                                          const std::string& clubb_grid_type,
                                          const stat_file& fread_var, // In
                                          real_t adjustment, LESProfile LES_tmp, int k_lowest_input,
                                          int k_highest_input, // In
                                          Real1D clubb_var)    // InOut
{
  // Interpolate an input profile and apply its constant unit adjustment.
  // Batched bracket outputs stay on device; only the Native per-level scan is serial.
  auto workspace_scope = workspace.scope();
  ConstInputProfile grid_var;
  if (clubb_grid_type == "zt")
    grid_var = Kokkos::subview(gr.zt, 0, Kokkos::ALL);
  else if (clubb_grid_type == "zm")
    grid_var = Kokkos::subview(gr.zm, 0, Kokkos::ALL);
  // Native leaves grid_var unassociated for any other grid type.
  const int nlevels = k_highest_input - k_lowest_input + 1;
  auto exact  = workspace_scope.view<InputIndices>(Kokkos::WithoutInitializing,
                                                   "inputfields_interp_and_adjust::exact", nlevels);
  auto lower  = workspace_scope.view<InputIndices>(Kokkos::WithoutInitializing,
                                                   "inputfields_interp_and_adjust::lower", nlevels);
  auto upper  = workspace_scope.view<InputIndices>(Kokkos::WithoutInitializing,
                                                   "inputfields_interp_and_adjust::upper", nlevels);
  auto linear = workspace_scope.view<Bool1D>(Kokkos::WithoutInitializing,
                                             "inputfields_interp_and_adjust::linear", nlevels);
  LES_grid_to_CLUBB_grid(fread_var, grid_var, k_lowest_input, k_highest_input, exact, lower, upper,
                         linear);
  const auto z = fread_var.z;
  const int ia = fread_var.ia;
  Kokkos::parallel_for(
      "inputfields_interp_and_adjust::apply", nlevels, KOKKOS_LAMBDA(int index) {
        const int k = index + k_lowest_input;
        if (linear(index))
          clubb_var(k) =
              lin_interpolate_two_points(grid_var(k), z(upper(index)), z(lower(index)),
                                         LES_tmp(upper(index) - ia), LES_tmp(lower(index) - ia));
        else
          clubb_var(k) = LES_tmp(exact(index) - ia);
        clubb_var(k) = clubb_var(k) * adjustment;
      });
}

void inputfields_init(int iunit, const std::string& namelist_filename, int day, int month,
                      int year,                       // In
                      inputfields_state& inputfields) // Native module InOut
{
  // Initialize inputfields settings and read the Native /setfields/ namelist.
  // Defaults/state updates are permanent C++; only raw formatted parsing uses
  // the existing transient descriptor owner. Flags omitted by Native's reset
  // list retain their case-owned values, rather than gaining a new reset here.
  inputfields.clubb_day                = day;
  inputfields.clubb_month              = month;
  inputfields.clubb_year               = year;
  std::string datafile                 = ""; // Native fixed120 local input prefix.
  inputfields.input_type               = "clubb";
  inputfields.l_input_um               = false;
  inputfields.l_input_vm               = false;
  inputfields.l_input_rtm              = false;
  inputfields.l_input_thlm             = false;
  inputfields.l_input_wp2              = false;
  inputfields.l_input_wprtp            = false;
  inputfields.l_input_wpthlp           = false;
  inputfields.l_input_wp3              = false;
  inputfields.l_input_rtp2             = false;
  inputfields.l_input_rtp3             = false;
  inputfields.l_input_thlp2            = false;
  inputfields.l_input_thlp3            = false;
  inputfields.l_input_rtpthlp          = false;
  inputfields.l_input_upwp             = false;
  inputfields.l_input_vpwp             = false;
  inputfields.l_input_ug               = false;
  inputfields.l_input_vg               = false;
  inputfields.l_input_rcm              = false;
  inputfields.l_input_wm_zt            = false;
  inputfields.l_input_exner            = false;
  inputfields.l_input_em               = false;
  inputfields.l_input_p                = false;
  inputfields.l_input_rho              = false;
  inputfields.l_input_Lscale           = false;
  inputfields.l_input_Lscale_up        = false;
  inputfields.l_input_Lscale_down      = false;
  inputfields.l_input_Kh_zt            = false;
  inputfields.l_input_Kh_zm            = false;
  inputfields.l_input_tau_zm           = false;
  inputfields.l_input_tau_zt           = false;
  inputfields.l_input_wpthvp           = false;
  inputfields.l_input_thl_1            = false;
  inputfields.l_input_thl_2            = false;
  inputfields.l_input_mixt_frac        = false;
  inputfields.l_input_chi_1            = false;
  inputfields.l_input_chi_2            = false;
  inputfields.l_input_stdev_chi_1      = false;
  inputfields.l_input_stdev_chi_2      = false;
  inputfields.l_input_rc_1             = false;
  inputfields.l_input_rc_2             = false;
  inputfields.l_input_thvm             = false;
  inputfields.l_input_rrm              = false;
  inputfields.l_input_Nrm              = false;
  inputfields.l_input_rsm              = false;
  inputfields.l_input_rim              = false;
  inputfields.l_input_rgm              = false;
  inputfields.l_input_rrp2             = false;
  inputfields.l_input_Nrp2             = false;
  inputfields.l_input_wprrp            = false;
  inputfields.l_input_wpNrp            = false;
  inputfields.l_input_thlm_forcing     = false;
  inputfields.l_input_rtm_forcing      = false;
  inputfields.l_input_up2              = false;
  inputfields.l_input_vp2              = false;
  inputfields.l_input_sigma_sqd_w      = false;
  inputfields.l_input_Ncm              = false;
  inputfields.l_input_Nccnm            = false;
  inputfields.l_input_Nim              = false;
  inputfields.l_input_Nsm              = false;
  inputfields.l_input_Ngm              = false;
  inputfields.l_input_cloud_frac       = false;
  inputfields.l_input_sigma_sqd_w_zt   = false;
  inputfields.l_input_veg_T_in_K       = false;
  inputfields.l_input_deep_soil_T_in_K = false;
  inputfields.l_input_sfc_soil_T_in_K  = false;
  inputfields.l_input_wprtp_forcing    = false;
  inputfields.l_input_wpthlp_forcing   = false;
  inputfields.l_input_rtp2_forcing     = false;
  inputfields.l_input_thlp2_forcing    = false;
  inputfields.l_input_rtpthlp_forcing  = false;
  inputfields.l_input_thlprcp          = false;
  inputfields.l_input_rcm_mc           = false;
  inputfields.l_input_rvm_mc           = false;
  inputfields.l_input_thlm_mc          = false;
  inputfields.l_input_wprtp_mc         = false;
  inputfields.l_input_wpthlp_mc        = false;
  inputfields.l_input_rtp2_mc          = false;
  inputfields.l_input_thlp2_mc         = false;
  inputfields.l_input_rtpthlp_mc       = false;
  std::cout << namelist_filename << "\n";
  read_setfields_nl(
      iunit, namelist_filename, datafile, inputfields.input_type, inputfields.l_input_um,
      inputfields.l_input_vm, inputfields.l_input_rtm, inputfields.l_input_thlm,
      inputfields.l_input_wp2, inputfields.l_input_wprtp, inputfields.l_input_wpthlp,
      inputfields.l_input_wp3, inputfields.l_input_rtp2, inputfields.l_input_rtp3,
      inputfields.l_input_thlp2, inputfields.l_input_thlp3, inputfields.l_input_rtpthlp,
      inputfields.l_input_upwp, inputfields.l_input_vpwp, inputfields.l_input_ug,
      inputfields.l_input_vg, inputfields.l_input_rcm, inputfields.l_input_wm_zt,
      inputfields.l_input_exner, inputfields.l_input_em, inputfields.l_input_p,
      inputfields.l_input_rho, inputfields.l_input_rho_zm, inputfields.l_input_Lscale,
      inputfields.l_input_Lscale_up, inputfields.l_input_Lscale_down, inputfields.l_input_Kh_zt,
      inputfields.l_input_Kh_zm, inputfields.l_input_tau_zm, inputfields.l_input_tau_zt,
      inputfields.l_input_wpthvp, inputfields.l_input_radht, inputfields.l_input_thl_1,
      inputfields.l_input_thl_2, inputfields.l_input_mixt_frac, inputfields.l_input_chi_1,
      inputfields.l_input_chi_2, inputfields.l_input_stdev_chi_1, inputfields.l_input_stdev_chi_2,
      inputfields.l_input_rc_1, inputfields.l_input_rc_2, inputfields.l_input_thvm,
      inputfields.l_input_rrm, inputfields.l_input_Nrm, inputfields.l_input_rsm,
      inputfields.l_input_rim, inputfields.l_input_rgm, inputfields.l_input_rrp2,
      inputfields.l_input_Nrp2, inputfields.l_input_wprrp, inputfields.l_input_wpNrp,
      inputfields.l_input_thlm_forcing, inputfields.l_input_rtm_forcing, inputfields.l_input_up2,
      inputfields.l_input_vp2, inputfields.l_input_sigma_sqd_w, inputfields.l_input_Ncm,
      inputfields.l_input_Nccnm, inputfields.l_input_Nim, inputfields.l_input_Ngm,
      inputfields.l_input_Nsm, inputfields.l_input_cloud_frac, inputfields.l_input_sigma_sqd_w_zt,
      inputfields.l_input_veg_T_in_K, inputfields.l_input_deep_soil_T_in_K,
      inputfields.l_input_sfc_soil_T_in_K, inputfields.l_input_wprtp_forcing,
      inputfields.l_input_wpthlp_forcing, inputfields.l_input_rtp2_forcing,
      inputfields.l_input_thlp2_forcing, inputfields.l_input_rtpthlp_forcing,
      inputfields.l_input_thlprcp, inputfields.l_input_rcm_mc, inputfields.l_input_rvm_mc,
      inputfields.l_input_thlm_mc, inputfields.l_input_wprtp_mc, inputfields.l_input_wpthlp_mc,
      inputfields.l_input_rtp2_mc, inputfields.l_input_thlp2_mc, inputfields.l_input_rtpthlp_mc);
  set_filenames(datafile, inputfields);
}

void cleanup_input_fields(inputfields_state& inputfields) // Native module InOut
{
  // Release the input-file names allocated during setup, before CLUBB ends.
  // Native DEALLOCATE requires this allocation to exist.
  if (inputfields.stat_files.empty())
    throw std::runtime_error("DEALLOCATE of unallocated stat_files in cleanup_input_fields");
  inputfields.stat_files.clear();
}

} // namespace clubb_kokkos
