// Initialize hydrometeor PDF metadata and normal-space correlation matrices.
//
// get_corr_var_index maps case-sensitive names to Native one-based indices;
// init_pdf_hydromet_arrays_api owns per-call metadata initialization and optional
// variance parameters. setup_corr_varnce_array_api fills/mirrors caller matrices;
// assert_corr_symmetric reports the first fatal error and print_corr_matrix
// displays the matrix. Small metadata/matrix work stages locally on the host.
// The generated hm_metadata_type remains the sole metadata owner. Declaration
// fragments and generator-preservation requirements are in INTEGRATION.md.
#include "corr_varnce_module.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "err_info_type_module.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp" // Pending public lower matrix contract.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>

namespace clubb_kokkos {
namespace {
using HostMatrix = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
//-----------------------------------------------------------------------------
int def_corr_idx(int iiPDF_x, const hm_metadata_type& hm_metadata) {
  // Description:
  //   Map from a iiPDF index to the corresponding index in the default
  //   correlation arrays.
  //-----------------------------------------------------------------------------

  // Indices represent the order in the default correlation arrays:
  // chi (old s), eta (old t), w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng.
  constexpr int ii_chi = 1;
  constexpr int ii_eta = 2;
  constexpr int ii_w = 3;
  constexpr int ii_Ncn = 4;
  constexpr int ii_rr = 5;
  constexpr int ii_Nr = 6;
  constexpr int ii_ri = 7;
  constexpr int ii_Ni = 8;
  constexpr int ii_rs = 9;
  constexpr int ii_Ns = 10;
  constexpr int ii_rg = 11;
  constexpr int ii_Ng = 12;
  int ii_def_corr = -1;
  if (iiPDF_x == hm_metadata.iiPDF_chi) {
    ii_def_corr = ii_chi;
  } else if (iiPDF_x == hm_metadata.iiPDF_eta) {
    ii_def_corr = ii_eta;
  } else if (iiPDF_x == hm_metadata.iiPDF_w) {
    ii_def_corr = ii_w;
  } else if (iiPDF_x == hm_metadata.iiPDF_Ncn) {
    ii_def_corr = ii_Ncn;
  } else if (iiPDF_x == hm_metadata.iiPDF_rr) {
    ii_def_corr = ii_rr;
  } else if (iiPDF_x == hm_metadata.iiPDF_Nr) {
    ii_def_corr = ii_Nr;
  } else if (iiPDF_x == hm_metadata.iiPDF_ri) {
    ii_def_corr = ii_ri;
  } else if (iiPDF_x == hm_metadata.iiPDF_Ni) {
    ii_def_corr = ii_Ni;
  } else if (iiPDF_x == hm_metadata.iiPDF_rs) {
    ii_def_corr = ii_rs;
  } else if (iiPDF_x == hm_metadata.iiPDF_Ns) {
    ii_def_corr = ii_Ns;
  } else if (iiPDF_x == hm_metadata.iiPDF_rg) {
    ii_def_corr = ii_rg;
  } else if (iiPDF_x == hm_metadata.iiPDF_Ng) {
    ii_def_corr = ii_Ng;
  }
  return ii_def_corr;
}

//-----------------------------------------------------------------------------
void set_corr_arrays_to_default(int pdf_dim, const hm_metadata_type& hm_metadata,             // In
                                HostMatrix corr_array_n_cloud, HostMatrix corr_array_n_below) // Out
{
  // Description:
  //   If there are no corr_array.in files for the current case, default
  //   correlations are used.
  //-----------------------------------------------------------------------------

  Kokkos::deep_copy(corr_array_n_cloud, zero);
  Kokkos::deep_copy(corr_array_n_below, zero);
  for (int i = 1; i <= pdf_dim; ++i) {
    corr_array_n_cloud(i - 1, i - 1) = one;
    corr_array_n_below(i - 1, i - 1) = one;
  }

  // Immutable Native tables need no global Kokkos allocation. Keep their const
  // unmanaged matrix aliases local to this metadata call; first index is fast.
  using DefaultMatrix = Kokkos::View<const real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  const DefaultMatrix corr_array_n_cloud_def_view(corr_array_n_cloud_def.data(), d_var_total,
                                                  d_var_total);
  const DefaultMatrix corr_array_n_below_def_view(corr_array_n_below_def.data(), d_var_total,
                                                  d_var_total);
  for (int i = 1; i <= pdf_dim - 1; ++i) {
    for (int j = i + 1; j <= pdf_dim; ++j) {
      const int idx_i = def_corr_idx(i, hm_metadata);
      const int idx_j = def_corr_idx(j, hm_metadata);
      if (idx_i > idx_j) {
        corr_array_n_cloud(j - 1, i - 1) = corr_array_n_cloud_def_view(idx_i - 1, idx_j - 1);
        corr_array_n_below(j - 1, i - 1) = corr_array_n_below_def_view(idx_i - 1, idx_j - 1);
      } else {
        corr_array_n_cloud(j - 1, i - 1) = corr_array_n_cloud_def_view(idx_j - 1, idx_i - 1);
        corr_array_n_below(j - 1, i - 1) = corr_array_n_below_def_view(idx_j - 1, idx_i - 1);
      }
    }
  }
}
} // namespace

//--------------------------------------------------------------------------
int get_corr_var_index(const std::string& var_name, const hm_metadata_type& hm_metadata) {
  // Description:
  //   Returns the index for a variable based on its name.
  //--------------------------------------------------------------------------

  // var_name is case sensitive; the result is Native one-based, or -1 if unknown.
  int i = -1;
  // Native TRIM removes trailing blanks only; leading blanks/case are retained.
  std::string trimmed = var_name;
  while (!trimmed.empty() && trimmed.back() == ' ')
    trimmed.pop_back();
  if (trimmed == "chi") {
    i = hm_metadata.iiPDF_chi;
  } else if (trimmed == "eta") {
    i = hm_metadata.iiPDF_eta;
  } else if (trimmed == "w") {
    i = hm_metadata.iiPDF_w;
  } else if (trimmed == "Ncn") {
    i = hm_metadata.iiPDF_Ncn;
  } else if (trimmed == "rr") {
    i = hm_metadata.iiPDF_rr;
  } else if (trimmed == "Nr") {
    i = hm_metadata.iiPDF_Nr;
  } else if (trimmed == "ri") {
    i = hm_metadata.iiPDF_ri;
  } else if (trimmed == "Ni") {
    i = hm_metadata.iiPDF_Ni;
  } else if (trimmed == "rs") {
    i = hm_metadata.iiPDF_rs;
  } else if (trimmed == "Ns") {
    i = hm_metadata.iiPDF_Ns;
  } else if (trimmed == "rg") {
    i = hm_metadata.iiPDF_rg;
  } else if (trimmed == "Ng") {
    i = hm_metadata.iiPDF_Ng;
  }
  return i;
}

//================================================================================================
// subroutine init_pdf_hydromet_arrays_api
//
// DESCRIPTION:
//     This subroutine intializes the hydromet arrays(iirr, iiNr, etc.) to the values specified by
//     the input arguements, this determines which hyrometeors are to be used by the microphysics
//     scheme. It also sets up the corresponding pdf and hydromet arrays, and calculates the
//     subgrid variance ratio for each hydrometeor.
//
// OPTIONAL FUNCTIONALITY:
//     The subgrid variance ratio for each hydrometeor is calculated based on the grid spacing
//     defined by the host model. The calculation is a linear equation defined by a slope and
//     intercept, each of which may or may not be passed in to this subroutine. If the slope
//     and/or intercept are not passed in through the arguement list the default values, which
//     are set in the corresponding type definitions, will be used. Otherwise the values
//     specified by the aruements will be used.
//
// NOTES:
//     'hmp2_ip_on_hmm2_ip_slope_in' is of type 'hmp2_ip_on_hmm2_ip_slope_type' and
//     'hmp2_ip_on_hmm2_ip_intrcpt_in' is of type 'hmp2_ip_on_hmm2_ip_intrcpt_type', both of which
//     are deinfed in corr_varnce_module.F90, and made public through this API.
//
//     If full control over the hydrometeor variance ratios is desired, pass in slopes that are
//     initialized to 0.0, this causes the ratios to no longer depend on the grid spacing. Then
//     pass in the intercepts set to the values of the desired ratios.
//
// ARGUEMENTS:
//     host_dx (real) - Horizontal grid spacings
//     host_dy (real)
//
//     hydromet_dim (integer) - Number of enabled hydrometeors
//
//         Each of these is an index value corresponding to a hydrometeor,
//         used to index the hydrometeor arrays. Each index has to be unqiue
//         for each different hyrometeor that is enabled. Setting one of these
//         indices to -1 disables that hydrometeor
//     iirr (integer) - Index of rain water mixing ratio
//     iiNr (integer) - Index of rain drop concentration
//     iiri (integer) - Index of ice mixing ratio
//     iiNi (integer) - Index of ice crystal concentration
//     iirs (integer) - Index of snow mixing ratio
//     iiNs (integer) - Index of snow flake concentration
//     iirg (integer) - Index of graupel mixing ratio
//     iiNg (integer) - Index of graupel concentration
//
//     hmp2_ip_on_hmm2_ip_slope_in (hmp2_ip_on_hmm2_ip_slope_type) - Custom slope values
//     hmp2_ip_on_hmm2_ip_intrcpt_in (hmp2_ip_on_hmm2_ip_intrcpt_type) - Custom intercept values
//
//================================================================================================
void init_pdf_hydromet_arrays_api(
    real_t host_dx, real_t host_dy, int hydromet_dim,                             // In
    int iirr, int iiNr, int iiri, int iiNi,                                       // In
    int iirs, int iiNs, int iirg, int iiNg,                                       // In
    real_t Ncnp2_on_Ncnm2,                                                        // In
    hm_metadata_type& hm_metadata, int& pdf_dim,                                  // Out
    std::optional<hmp2_ip_on_hmm2_ip_slope_type> hmp2_ip_on_hmm2_ip_slope_in,     // Optional In
    std::optional<hmp2_ip_on_hmm2_ip_intrcpt_type> hmp2_ip_on_hmm2_ip_intrcpt_in) // Optional In
{
  // Description:
  //
  // Initialization for the hydromet arrays. How the arrays are initialized is
  // determined by how the hydromet indicies (iirr, iiNr, etc.) have been set,
  // which are determined by the microphysics scheme.
  //-------------------------------------------------------------------------------

  // host_dx/host_dy: horizontal grid spacing [m]. hydromet_dim counts enabled
  // hydrometeor fields, including mixing ratios and concentrations. Species
  // indices are unique Native one-based slots; negative values disable fields.
  // Optional slopes [1/m] and intercepts [-] override the Native type defaults.
  // Native INTENT(OUT) deallocates old array components and restores declared
  // scalar defaults before initialization; reuse the generated type itself.
  hm_metadata = {};
  int pdf_count;
  hmp2_ip_on_hmm2_ip_slope_type hmp2_ip_on_hmm2_ip_slope;
  hmp2_ip_on_hmm2_ip_intrcpt_type hmp2_ip_on_hmm2_ip_intrcpt;

  // Save indices
  hm_metadata.iirr = iirr;
  hm_metadata.iiri = iiri;
  hm_metadata.iirs = iirs;
  hm_metadata.iirg = iirg;
  hm_metadata.iiNr = iiNr;
  hm_metadata.iiNi = iiNi;
  hm_metadata.iiNs = iiNs;
  hm_metadata.iiNg = iiNg;

  // Overwrite default Ncnp2_on_Ncnm2
  hm_metadata.Ncnp2_on_Ncnm2 = Ncnp2_on_Ncnm2;

  //-----------------------------------------------------------
  // Calculate the subgrid variances of the hydrometeors
  //-----------------------------------------------------------

  // If slope and intercept are present in call, then overwrite default values
  if (hmp2_ip_on_hmm2_ip_slope_in.has_value()) {
    hmp2_ip_on_hmm2_ip_slope = *hmp2_ip_on_hmm2_ip_slope_in;
  }

  if (hmp2_ip_on_hmm2_ip_intrcpt_in.has_value()) {
    hmp2_ip_on_hmm2_ip_intrcpt = *hmp2_ip_on_hmm2_ip_intrcpt_in;
  }

  hm_metadata.hmp2_ip_on_hmm2_ip = decltype(hm_metadata.hmp2_ip_on_hmm2_ip)(
      "init_pdf_hydromet_arrays_api::hmp2_ip_on_hmm2_ip", hydromet_dim);
  auto hmp2_ip_on_hmm2_ip_host = Kokkos::create_mirror_view(hm_metadata.hmp2_ip_on_hmm2_ip);
  Kokkos::deep_copy(hmp2_ip_on_hmm2_ip_host, real_t(0.0));

  if (iirr > 0) {
    hmp2_ip_on_hmm2_ip_host(iirr - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.rr + hmp2_ip_on_hmm2_ip_slope.rr * std::max(host_dx, host_dy);
  }

  if (iirs > 0) {
    hmp2_ip_on_hmm2_ip_host(iirs - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.rs + hmp2_ip_on_hmm2_ip_slope.rs * std::max(host_dx, host_dy);
  }

  if (iiri > 0) {
    hmp2_ip_on_hmm2_ip_host(iiri - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.ri + hmp2_ip_on_hmm2_ip_slope.ri * std::max(host_dx, host_dy);
  }

  if (iirg > 0) {
    hmp2_ip_on_hmm2_ip_host(iirg - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.rg + hmp2_ip_on_hmm2_ip_slope.rg * std::max(host_dx, host_dy);
  }

  if (iiNr > 0) {
    hmp2_ip_on_hmm2_ip_host(iiNr - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.Nr + hmp2_ip_on_hmm2_ip_slope.Nr * std::max(host_dx, host_dy);
  }

  if (iiNs > 0) {
    hmp2_ip_on_hmm2_ip_host(iiNs - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.Ns + hmp2_ip_on_hmm2_ip_slope.Ns * std::max(host_dx, host_dy);
  }

  if (iiNi > 0) {
    hmp2_ip_on_hmm2_ip_host(iiNi - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.Ni + hmp2_ip_on_hmm2_ip_slope.Ni * std::max(host_dx, host_dy);
  }

  if (iiNg > 0) {
    hmp2_ip_on_hmm2_ip_host(iiNg - 1) =
        hmp2_ip_on_hmm2_ip_intrcpt.Ng + hmp2_ip_on_hmm2_ip_slope.Ng * std::max(host_dx, host_dy);
  }

  //-----------------------------------------------------------
  // Set up predictive precipitating hydrometeor arrays.
  //-----------------------------------------------------------
  hm_metadata.hydromet_list = decltype(hm_metadata.hydromet_list)(
      "init_pdf_hydromet_arrays_api::hydromet_list", 10, hydromet_dim);
  auto hydromet_list_host = Kokkos::create_mirror_view(hm_metadata.hydromet_list);
  hm_metadata.hydromet_tol = decltype(hm_metadata.hydromet_tol)(
      "init_pdf_hydromet_arrays_api::hydromet_tol", hydromet_dim);
  auto hydromet_tol_host = Kokkos::create_mirror_view(hm_metadata.hydromet_tol);
  hm_metadata.l_mix_rat_hm = decltype(hm_metadata.l_mix_rat_hm)(
      "init_pdf_hydromet_arrays_api::l_mix_rat_hm", hydromet_dim);
  auto l_mix_rat_hm_host = Kokkos::create_mirror_view(hm_metadata.l_mix_rat_hm);
  hm_metadata.l_frozen_hm =
      decltype(hm_metadata.l_frozen_hm)("init_pdf_hydromet_arrays_api::l_frozen_hm", hydromet_dim);
  auto l_frozen_hm_host = Kokkos::create_mirror_view(hm_metadata.l_frozen_hm);

  if (iirr > 0) {
    // The microphysics scheme predicts rain water mixing ratio, rr.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iirr - 1) = "rrm       "[character];
    }
    l_mix_rat_hm_host(iirr - 1) = true;
    l_frozen_hm_host(iirr - 1) = false;
    hydromet_tol_host(iirr - 1) = rr_tol;
  }

  if (iiri > 0) {
    // The microphysics scheme predicts ice mixing ratio, ri.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iiri - 1) = "rim       "[character];
    }
    l_mix_rat_hm_host(iiri - 1) = true;
    l_frozen_hm_host(iiri - 1) = true;
    hydromet_tol_host(iiri - 1) = ri_tol;
  }

  if (iirs > 0) {
    // The microphysics scheme predicts snow mixing ratio, rs.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iirs - 1) = "rsm       "[character];
    }
    l_mix_rat_hm_host(iirs - 1) = true;
    l_frozen_hm_host(iirs - 1) = true;
    hydromet_tol_host(iirs - 1) = rs_tol;
  }

  if (iirg > 0) {
    // The microphysics scheme predicts graupel mixing ratio, rg.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iirg - 1) = "rgm       "[character];
    }
    l_mix_rat_hm_host(iirg - 1) = true;
    l_frozen_hm_host(iirg - 1) = true;
    hydromet_tol_host(iirg - 1) = rg_tol;
  }

  if (iiNr > 0) {
    // The microphysics scheme predicts rain drop concentration, Nr.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iiNr - 1) = "Nrm       "[character];
    }
    l_frozen_hm_host(iiNr - 1) = false;
    l_mix_rat_hm_host(iiNr - 1) = false;
    hydromet_tol_host(iiNr - 1) = Nr_tol;
  }

  if (iiNi > 0) {
    // The microphysics scheme predicts ice concentration, Ni.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iiNi - 1) = "Nim       "[character];
    }
    l_mix_rat_hm_host(iiNi - 1) = false;
    l_frozen_hm_host(iiNi - 1) = true;
    hydromet_tol_host(iiNi - 1) = Ni_tol;
  }

  if (iiNs > 0) {
    // The microphysics scheme predicts snowflake concentration, Ns.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iiNs - 1) = "Nsm       "[character];
    }
    l_mix_rat_hm_host(iiNs - 1) = false;
    l_frozen_hm_host(iiNs - 1) = true;
    hydromet_tol_host(iiNs - 1) = Ns_tol;
  }

  if (iiNg > 0) {
    // The microphysics scheme predicts graupel concentration, Ng.
    for (int character = 0; character < 10; ++character) {
      hydromet_list_host(character, iiNg - 1) = "Ngm       "[character];
    }
    l_mix_rat_hm_host(iiNg - 1) = false;
    l_frozen_hm_host(iiNg - 1) = true;
    hydromet_tol_host(iiNg - 1) = Ng_tol;
  }

  //-----------------------------------------------------------
  // Set up the PDF indices
  //-----------------------------------------------------------

  hm_metadata.iiPDF_chi = 1; // Extended liquid water mixing ratio, chi
  hm_metadata.iiPDF_eta = 2; // 'eta' orthogonal to 'chi'
  hm_metadata.iiPDF_w = 3;   // vertical velocity
  hm_metadata.iiPDF_Ncn = 4; // Simplified cloud nuclei concentration or extended Nc.

  pdf_count = hm_metadata.iiPDF_Ncn;

  // Loop over hydrometeors.
  // Hydrometeor indices in the PDF arrays should be in the same order as
  // found in the hydrometeor arrays.
  if (hydromet_dim > 0) {

    for (int i = 1; i <= hydromet_dim; ++i) {

      if (i == iirr) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_rr = pdf_count;
      }

      if (i == iiNr) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_Nr = pdf_count;
      }

      if (i == iiri) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_ri = pdf_count;
      }

      if (i == iiNi) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_Ni = pdf_count;
      }

      if (i == iirs) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_rs = pdf_count;
      }

      if (i == iiNs) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_Ns = pdf_count;
      }

      if (i == iirg) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_rg = pdf_count;
      }

      if (i == iiNg) {
        pdf_count = pdf_count + 1;
        hm_metadata.iiPDF_Ng = pdf_count;
      }

    } // i = 1, hydromet_dim, 1

  } // hydromet_dim > 0

  pdf_dim = pdf_count;

  // Restore the generated View components once at this metadata boundary.
  Kokkos::deep_copy(hm_metadata.hmp2_ip_on_hmm2_ip, hmp2_ip_on_hmm2_ip_host);
  Kokkos::deep_copy(hm_metadata.hydromet_list, hydromet_list_host);
  Kokkos::deep_copy(hm_metadata.hydromet_tol, hydromet_tol_host);
  Kokkos::deep_copy(hm_metadata.l_mix_rat_hm, l_mix_rat_hm_host);
  Kokkos::deep_copy(hm_metadata.l_frozen_hm, l_frozen_hm_host);
}

//===============================================================================
void setup_corr_varnce_array_api(int pdf_dim, const hm_metadata_type& hm_metadata,     // In
                                 bool l_fix_w_chi_eta_correlations,                    // In
                                 Real2D corr_array_n_cloud, Real2D corr_array_n_below, // Out
                                 std::optional<Real2D> corr_array_n_cloud_default, // Optional In
                                 std::optional<Real2D> corr_array_n_below_default) // Optional In
{
  // Description:
  //   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
  //   elements to be correlations between various variables.

  // References:
  //   None.
  //-------------------------------------------------------------------------------

  // l_fix_w_chi_eta_correlations selects fixed Mellor chi/eta correlations.
  // Both optional normal-space matrices must be present to override defaults.
  // Metadata matrices are small host-side setup values; keep only call-local
  // mirrors and restore the caller's allocations after setup/diagnostics.
  HostMatrix corr_array_n_cloud_host("setup_corr_varnce_array_api::cloud", pdf_dim, pdf_dim);
  HostMatrix corr_array_n_below_host("setup_corr_varnce_array_api::below", pdf_dim, pdf_dim);
  if (corr_array_n_cloud_default.has_value() && corr_array_n_below_default.has_value()) {
    Kokkos::deep_copy(corr_array_n_cloud_host, *corr_array_n_cloud_default);
    Kokkos::deep_copy(corr_array_n_below_host, *corr_array_n_below_default);
  } else {
    set_corr_arrays_to_default(pdf_dim, hm_metadata,                              // In
                               corr_array_n_cloud_host, corr_array_n_below_host); // Out
  }

  // Mirror the correlation matrices.
  // TODO: Replace this explicit public lower boundary after its owner is ported.
  f2k::mirror_lower_triangular_matrix(pdf_dim,                  // In
                                      corr_array_n_cloud_host); // InOut
  f2k::mirror_lower_triangular_matrix(pdf_dim,                  // In
                                      corr_array_n_below_host); // InOut

  // Sanity check to avoid confusing non-convergence results.
  // TODO: Replace this host query after the Native error-code owner is ported.
  if (f2k::clubb_at_least_debug_level_api(2)) {
    if (!l_fix_w_chi_eta_correlations && hm_metadata.iiPDF_Ncn > 0) {
      bool l_warning = false;
      for (int i = 1; i <= pdf_dim; ++i) {
        if (((std::abs(corr_array_n_cloud_host(i - 1, hm_metadata.iiPDF_Ncn - 1)) > eps) ||
             (std::abs(corr_array_n_below_host(i - 1, hm_metadata.iiPDF_Ncn - 1)) > eps)) &&
            i != hm_metadata.iiPDF_Ncn) {
          l_warning = true;
        }
      }
      if (l_warning) {
        std::cerr << "Warning: the specified correlations for chi"
                     " (old s) and Ncn are non-zero.\n";
        std::cerr << "The latin hypercube code will not converge to"
                     " the analytic solution using these settings.\n";
      }
    }
  }
  from_execution_space(corr_array_n_cloud, corr_array_n_cloud_host);
  from_execution_space(corr_array_n_below, corr_array_n_below_host);
}

//-----------------------------------------------------------------------------
void assert_corr_symmetric(Real2D corr_array_n,     // In
                           int pdf_dim,             // In
                           err_info_type& err_info) // InOut
{
  // Description:
  //   Asserts that corr_matrix(i,j) == corr_matrix(j,i) for all indeces
  //   in the correlation array. On the first failure, marks all columns fatal and returns.
  // References:
  //   None
  //---------------------------------------------------------------------------

  // pdf_dim is the number of variables. err_info receives the fatal status
  // and supplies per-column/global diagnostic headers. Matrices are dimensionless.
  constexpr real_t tol = real_t(1.0e-6); // Tolerance for real-precision testing.
  const auto corr_array_n_host = to_host_space(corr_array_n);

  // Do the check in Native column/row order, with asymmetry before diagonal.
  for (int n_col = 1; n_col <= pdf_dim; ++n_col) {
    for (int n_row = 1; n_row <= pdf_dim; ++n_row) {
      // These errors occur inside an individual grid column.
      if (std::abs(corr_array_n_host(n_col - 1, n_row - 1) -
                   corr_array_n_host(n_row - 1, n_col - 1)) > tol) {
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        const auto err_header_host = to_host_space(err_info.err_header);
        std::ostringstream message;
        for (std::size_t column = 0; column < err_header_host.extent(1); ++column) {
          for (std::size_t character = 0; character < err_header_host.extent(0); ++character)
            message << err_header_host(character, column);
          message << ' ';
        }
        message << '\n';
        message << "in subroutine assert_corr_symmetric: "
                   "Correlation array is non symmetric.\n";
        // Native whole-array output traverses the first index fastest.
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        for (int column = 0; column < pdf_dim; ++column)
          for (int row = 0; row < pdf_dim; ++row)
            message << corr_array_n_host(row, column) << ' ';
        message << '\n';
        std::cerr << message.str();
        return;
      }

      if (n_col == n_row && std::abs(corr_array_n_host(n_col - 1, n_row - 1) - one) > eps) {
        Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
        std::ostringstream message;
        message << err_info.err_header_global << '\n';
        message << "in subroutine assert_corr_symmetric: "
                   "Correlation array is formatted incorrectly.\n";
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        for (int column = 0; column < pdf_dim; ++column)
          for (int row = 0; row < pdf_dim; ++row)
            message << corr_array_n_host(row, column) << ' ';
        message << '\n';
        std::cerr << message.str();
        return;
      }
    }
  }
}

//-----------------------------------------------------------------------------
void print_corr_matrix(int pdf_dim,         // In
                       Real2D corr_array_n) // In
{
  // Description:
  //   Prints the correlation matrix to the console.
  // References:
  //   None
  //---------------------------------------------------------------------------

  // pdf_dim is the number of variables; corr_array_n is dimensionless.
  const auto corr_array_n_host = to_host_space(corr_array_n);
  int current_character_index = 0;
  std::string current_line(72, ' ');
  for (int n = 1; n <= pdf_dim; ++n) {
    for (int m = 1; m <= pdf_dim; ++m) {
      // Preserve Native F5.2 in LEN10, six-character strides and LEN72 lines.
      // C++ streams require explicit width-overflow stars and blank padding.
      std::ostringstream value;
      value << std::fixed << std::setprecision(2) << std::setw(5)
            << corr_array_n_host(m - 1, n - 1);
      std::string str_array_value = value.str();
      if (str_array_value.size() > 5)
        str_array_value = "*****";
      str_array_value.resize(10, ' ');
      current_line =
          (current_line.substr(0, current_character_index) + str_array_value).substr(0, 72);
      current_line.resize(72, ' ');
      current_character_index = current_character_index + 6;
    }
    std::cout << current_line << '\n';
    current_line.assign(72, ' ');
    current_character_index = 0;
  }
}
//-----------------------------------------------------------------------------
} // namespace clubb_kokkos
