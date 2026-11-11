#pragma once

#include "scratch_workspace.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"
#include "array_types.hpp"
#include "err_info_type_module.hpp"
#include <array>
#include <optional>

namespace clubb_kokkos {

struct hmp2_ip_on_hmm2_ip_ratios_type {
  // Prescribed parameters for hydrometeor values of <hm|_ip'^2> / <hm|_ip>^2,
  // where  <hm|_ip'^2>  is the in-precip. variance of the hydrometeor and
  //        <hm|_ip>     is the in-precip. mean of the hydrometeor
  //
  // These values are dependent on the horizontal grid spacing of the run, and are calculated
  // using a slope and intercept corresponding to each hydrometer.
  real_t rr = real_t(1.0); // For rain water mixing ratio [-]
  real_t Nr = real_t(1.0); // For rain drop concentration [-]
  real_t ri = real_t(1.0); // For ice mixing ratio        [-]
  real_t Ni = real_t(1.0); // For ice concentration       [-]
  real_t rs = real_t(1.0); // For snow mixing ratio       [-]
  real_t Ns = real_t(1.0); // For snow concentration      [-]
  real_t rg = real_t(1.0); // For graupel mixing ratio    [-]
  real_t Ng = real_t(1.0); // For graupel concentration   [-]
};

// These slopes and intercepts below are used to calculate the hmp2_ip_on_hmm2_ip_ratios_type
// values that are defined above. This functionality is described by equations 8, 10, & 11
// in ``Parameterization of the Spatial Variability of Rain for Large-Scale Models and
// Remote Sensing`` Lebo, et al, October 2015
// https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-15-0066.1
// see clubb:ticket:830 for more detail
//
// hmp2_ip_on_hmm2_ip(iirr) = hmp2_ip_on_hmm2_ip_intrcpt%rr + &
//                            hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )
//
// In Lebo et al. the suggested values were
//     slope = 2.12e-5 [1/m]
//     intercept = 0.54 [-]
//
// In CLUBB standalone, these parameters can be set based on the value for a
// given case in the CASE_model.in file.

struct hmp2_ip_on_hmm2_ip_slope_type {
  real_t rr = real_t(2.12e-5); // For rain water mixing ratio [1/m]
  real_t Nr = real_t(2.12e-5); // For rain drop concentration [1/m]
  real_t ri = real_t(2.12e-5); // For ice mixing ratio        [1/m]
  real_t Ni = real_t(2.12e-5); // For ice concentration       [1/m]
  real_t rs = real_t(2.12e-5); // For snow mixing ratio       [1/m]
  real_t Ns = real_t(2.12e-5); // For snow concentration      [1/m]
  real_t rg = real_t(2.12e-5); // For graupel mixing ratio    [1/m]
  real_t Ng = real_t(2.12e-5); // For graupel concentration   [1/m]
};

struct hmp2_ip_on_hmm2_ip_intrcpt_type {
  real_t rr = real_t(0.54); // For rain water mixing ratio [-]
  real_t Nr = real_t(0.54); // For rain drop concentration [-]
  real_t ri = real_t(0.54); // For ice mixing ratio        [-]
  real_t Ni = real_t(0.54); // For ice concentration       [-]
  real_t rs = real_t(0.54); // For snow mixing ratio       [-]
  real_t Ns = real_t(0.54); // For snow concentration      [-]
  real_t rg = real_t(0.54); // For graupel mixing ratio    [-]
  real_t Ng = real_t(0.54); // For graupel concentration   [-]
};

// BEGIN GENERATED HM METADATA TYPE
// Mirrors src/CLUBB_core/corr_varnce_module.F90; arrays use Fortran-order LayoutLeft.
struct hm_metadata_type {
  int iirr = -1;
  int iirs = -1;
  int iiri = -1;
  int iirg = -1;
  int iiNr = -1;
  int iiNs = -1;
  int iiNi = -1;
  int iiNg = -1;
  Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft> l_frozen_hm;
  Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft> l_mix_rat_hm;
  Kokkos::View<char**, Kokkos::LayoutLeft> hydromet_list;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> hydromet_tol;
  int iiPDF_chi = -2;
  int iiPDF_eta = -2;
  int iiPDF_w   = -2;
  int iiPDF_rr  = -2;
  int iiPDF_rs  = -2;
  int iiPDF_ri  = -2;
  int iiPDF_rg  = -2;
  int iiPDF_Nr  = -2;
  int iiPDF_Ns  = -2;
  int iiPDF_Ni  = -2;
  int iiPDF_Ng  = -2;
  int iiPDF_Ncn = -2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> hmp2_ip_on_hmm2_ip;
  real_t Ncnp2_on_Ncnm2 = 1.0;
};

// END GENERATED HM METADATA TYPE

// Repository standalone/JAX code only reads these default tables. They are
// immutable C++ constants; external mutation of the public Native arrays is
// not a supported C++ API. Fortran-flat order retains its RESHAPE semantics.
inline constexpr int d_var_total = 12; // Size of the default correlation arrays.

inline constexpr std::array<real_t, d_var_total * d_var_total> corr_array_n_cloud_def = {
    real_t(1.),   real_t(-.6),  real_t(.09),   real_t(.09),
    real_t(.788), real_t(.675), real_t(.240),  real_t(.222),
    real_t(.240), real_t(.222), real_t(.240),  real_t(.222), // chi Native column
    real_t(0.),   real_t(1.),   real_t(.027),  real_t(.027),
    real_t(.114), real_t(.115), real_t(-.029), real_t(.093),
    real_t(.022), real_t(.013), real_t(0.),    real_t(0.), // eta Native column
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.34),
    real_t(.315), real_t(.270), real_t(.120),  real_t(.167),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.), // w Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(0.),   real_t(0.),   real_t(.464),  real_t(.320),
    real_t(.168), real_t(.232), real_t(0.),    real_t(0.), // Ncn Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.821), real_t(0.),    real_t(0.),
    real_t(.173), real_t(.164), real_t(.319),  real_t(.308), // rr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.152),  real_t(.143),
    real_t(0.),   real_t(0.),   real_t(.285),  real_t(.273), // Nr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758),
    real_t(.585), real_t(.571), real_t(.379),  real_t(.363), // ri Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(.571), real_t(.550), real_t(.363),  real_t(.345), // Ni Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.758), real_t(.485),  real_t(.470), // rs Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.470),  real_t(.450), // Ns Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758), // rg Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.), // Ng Native column
};

inline constexpr std::array<real_t, d_var_total * d_var_total> corr_array_n_below_def = {
    real_t(1.),   real_t(.3),   real_t(.09),   real_t(.09),
    real_t(.788), real_t(.675), real_t(.240),  real_t(.222),
    real_t(.240), real_t(.222), real_t(.240),  real_t(.222), // chi Native column
    real_t(0.),   real_t(1.),   real_t(.027),  real_t(.027),
    real_t(.114), real_t(.115), real_t(-.029), real_t(.093),
    real_t(.022), real_t(.013), real_t(0.),    real_t(0.), // eta Native column
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.34),
    real_t(.315), real_t(.270), real_t(.120),  real_t(.167),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.), // w Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(0.),   real_t(0.),   real_t(.464),  real_t(.320),
    real_t(.168), real_t(.232), real_t(0.),    real_t(0.), // Ncn Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.821), real_t(0.),    real_t(0.),
    real_t(.173), real_t(.164), real_t(.319),  real_t(.308), // rr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.152),  real_t(.143),
    real_t(0.),   real_t(0.),   real_t(.285),  real_t(.273), // Nr Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758),
    real_t(.585), real_t(.571), real_t(.379),  real_t(.363), // ri Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.),
    real_t(.571), real_t(.550), real_t(.363),  real_t(.345), // Ni Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(1.),   real_t(.758), real_t(.485),  real_t(.470), // rs Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(1.),   real_t(.470),  real_t(.450), // Ns Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(1.),    real_t(.758), // rg Native column
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(0.),
    real_t(0.),   real_t(0.),   real_t(0.),    real_t(1.), // Ng Native column
};

int get_corr_var_index(const std::string& var_name, const hm_metadata_type& hm_metadata);

void init_pdf_hydromet_arrays_api(
    real_t host_dx, real_t host_dy, int hydromet_dim, // In
    int iirr, int iiNr, int iiri, int iiNi,           // In
    int iirs, int iiNs, int iirg, int iiNg,           // In
    real_t Ncnp2_on_Ncnm2,                            // In
    hm_metadata_type& hm_metadata, int& pdf_dim,      // Out
    std::optional<hmp2_ip_on_hmm2_ip_slope_type> hmp2_ip_on_hmm2_ip_slope_in =
        std::nullopt, // Optional In
    std::optional<hmp2_ip_on_hmm2_ip_intrcpt_type> hmp2_ip_on_hmm2_ip_intrcpt_in =
        std::nullopt); // Optional In

void setup_corr_varnce_array_api(
    ScratchWorkspace& workspace,                                      // InOut
    int pdf_dim, const hm_metadata_type& hm_metadata,                 // In
    bool l_fix_w_chi_eta_correlations,                                // In
    Real2D corr_array_n_cloud, Real2D corr_array_n_below,             // Out
    std::optional<Real2D> corr_array_n_cloud_default = std::nullopt,  // Optional In
    std::optional<Real2D> corr_array_n_below_default = std::nullopt); // Optional In

void assert_corr_symmetric(
    ScratchWorkspace& workspace,
    Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> corr_array_n, // In
    int pdf_dim,                                                                // In
    err_info_type& err_info);                                                   // InOut
void print_corr_matrix(int pdf_dim,                                             // In
                       Real2D corr_array_n);                                    // In

// Existing shared array ownership adapts host/strided caller sections once.
template <class Cloud, class Below, class DefaultCloud = Real2D, class DefaultBelow = Real2D>
void setup_corr_varnce_array_api(
    ScratchWorkspace& workspace,                                           // InOut
    int pdf_dim, const hm_metadata_type& hm_metadata,                      // In
    bool l_fix_w_chi_eta_correlations,                                     // In
    Cloud corr_array_n_cloud, Below corr_array_n_below,                    // Out
    std::optional<DefaultCloud> corr_array_n_cloud_default = std::nullopt, // Optional In
    std::optional<DefaultBelow> corr_array_n_below_default = std::nullopt) // Optional In
{
  // Keep caller storage while staging only the existing public array boundary.
  auto workspace_scope = workspace.scope();
  auto cloud_exec      = to_execution_space(workspace, corr_array_n_cloud, false);
  auto below_exec      = to_execution_space(workspace, corr_array_n_below, false);
  std::optional<Real2D> cloud_default_exec, below_default_exec;
  if (corr_array_n_cloud_default.has_value() && corr_array_n_below_default.has_value()) {
    cloud_default_exec = to_execution_space(workspace, *corr_array_n_cloud_default);
    below_default_exec = to_execution_space(workspace, *corr_array_n_below_default);
  }
  setup_corr_varnce_array_api(workspace, pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations,
                              cloud_exec, below_exec, cloud_default_exec, below_default_exec);
  from_execution_space(workspace, corr_array_n_cloud, cloud_exec);
  from_execution_space(workspace, corr_array_n_below, below_exec);
}
template <class Matrix>
void assert_corr_symmetric(ScratchWorkspace& workspace, Matrix corr_array_n, int pdf_dim,
                           err_info_type& err_info) {
  // Preserve an independent snapshot for the ordered host correlation check.
  auto workspace_scope = workspace.scope();
  assert_corr_symmetric(workspace, to_host_space(workspace, corr_array_n), pdf_dim, err_info);
}
template <class Matrix> void print_corr_matrix(int pdf_dim, Matrix corr_array_n) {
  // Preserve the existing public array boundary and Native matrix orientation.
  print_corr_matrix(pdf_dim, to_execution_space(corr_array_n));
}

} // namespace clubb_kokkos
