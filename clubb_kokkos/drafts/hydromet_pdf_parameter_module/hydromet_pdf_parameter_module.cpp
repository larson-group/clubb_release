// Initialize hydrometeor PDF parameters and precipitation fractions.
//
// Public routines:
//   init_hydromet_pdf_params: zero fixed hydrometeor parameters or a whole grid collection.
//   init_precip_fracs_api: allocate and zero three precipitation-fraction fields.
//   zero_precip_fracs_api: reset those fields while retaining their allocations.
//
// Existing generated types own all storage. Grid fields use (column,level)
// LayoutLeft Views. Precipitation initialization replaces the caller's prior
// fields, whereas reset routines preserve the existing case allocations.
#include "hydromet_pdf_parameter_module.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <stdexcept>

namespace clubb_kokkos {

// =============================================================================
void init_hydromet_pdf_params(
    hydromet_pdf_parameter& hydromet_pdf_params // InOut
) {
  // Description:
  // Initialize the elements of hydromet_pdf_params.
  // References:
  //-----------------------------------------------------------------------
  // hydromet_pdf_params: Hydrometeor PDF parameters [units vary].

  // Native fixed-size array fields always exist. Allocate their C++ View
  // storage once; preserve matching caller allocations when resetting values.
  constexpr int max_hydromet_dim = hydromet_pdf_parameter_collection::max_hydromet_dim;

  if (!hydromet_pdf_params.hm_1.is_allocated()) {
    hydromet_pdf_params.hm_1 =
        decltype(hydromet_pdf_params.hm_1)("hydromet_pdf_params::hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field hm_1");

  if (!hydromet_pdf_params.hm_2.is_allocated()) {
    hydromet_pdf_params.hm_2 =
        decltype(hydromet_pdf_params.hm_2)("hydromet_pdf_params::hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field hm_2");

  if (!hydromet_pdf_params.mu_hm_1.is_allocated()) {
    hydromet_pdf_params.mu_hm_1 =
        decltype(hydromet_pdf_params.mu_hm_1)("hydromet_pdf_params::mu_hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.mu_hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field mu_hm_1");

  if (!hydromet_pdf_params.mu_hm_2.is_allocated()) {
    hydromet_pdf_params.mu_hm_2 =
        decltype(hydromet_pdf_params.mu_hm_2)("hydromet_pdf_params::mu_hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.mu_hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field mu_hm_2");

  if (!hydromet_pdf_params.sigma_hm_1.is_allocated()) {
    hydromet_pdf_params.sigma_hm_1 = decltype(hydromet_pdf_params.sigma_hm_1)(
        "hydromet_pdf_params::sigma_hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.sigma_hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field sigma_hm_1");

  if (!hydromet_pdf_params.sigma_hm_2.is_allocated()) {
    hydromet_pdf_params.sigma_hm_2 = decltype(hydromet_pdf_params.sigma_hm_2)(
        "hydromet_pdf_params::sigma_hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.sigma_hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field sigma_hm_2");

  if (!hydromet_pdf_params.corr_w_hm_1.is_allocated()) {
    hydromet_pdf_params.corr_w_hm_1 = decltype(hydromet_pdf_params.corr_w_hm_1)(
        "hydromet_pdf_params::corr_w_hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_w_hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_w_hm_1");

  if (!hydromet_pdf_params.corr_w_hm_2.is_allocated()) {
    hydromet_pdf_params.corr_w_hm_2 = decltype(hydromet_pdf_params.corr_w_hm_2)(
        "hydromet_pdf_params::corr_w_hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_w_hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_w_hm_2");

  if (!hydromet_pdf_params.corr_chi_hm_1.is_allocated()) {
    hydromet_pdf_params.corr_chi_hm_1 = decltype(hydromet_pdf_params.corr_chi_hm_1)(
        "hydromet_pdf_params::corr_chi_hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_chi_hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_chi_hm_1");

  if (!hydromet_pdf_params.corr_chi_hm_2.is_allocated()) {
    hydromet_pdf_params.corr_chi_hm_2 = decltype(hydromet_pdf_params.corr_chi_hm_2)(
        "hydromet_pdf_params::corr_chi_hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_chi_hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_chi_hm_2");

  if (!hydromet_pdf_params.corr_eta_hm_1.is_allocated()) {
    hydromet_pdf_params.corr_eta_hm_1 = decltype(hydromet_pdf_params.corr_eta_hm_1)(
        "hydromet_pdf_params::corr_eta_hm_1", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_eta_hm_1.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_eta_hm_1");

  if (!hydromet_pdf_params.corr_eta_hm_2.is_allocated()) {
    hydromet_pdf_params.corr_eta_hm_2 = decltype(hydromet_pdf_params.corr_eta_hm_2)(
        "hydromet_pdf_params::corr_eta_hm_2", max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_eta_hm_2.extent_int(0) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_eta_hm_2");

  if (!hydromet_pdf_params.corr_hmx_hmy_1.is_allocated()) {
    hydromet_pdf_params.corr_hmx_hmy_1 = decltype(hydromet_pdf_params.corr_hmx_hmy_1)(
        "hydromet_pdf_params::corr_hmx_hmy_1", max_hydromet_dim, max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_hmx_hmy_1.extent_int(0) != max_hydromet_dim ||
      hydromet_pdf_params.corr_hmx_hmy_1.extent_int(1) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_hmx_hmy_1");

  if (!hydromet_pdf_params.corr_hmx_hmy_2.is_allocated()) {
    hydromet_pdf_params.corr_hmx_hmy_2 = decltype(hydromet_pdf_params.corr_hmx_hmy_2)(
        "hydromet_pdf_params::corr_hmx_hmy_2", max_hydromet_dim, max_hydromet_dim);
  }

  if (hydromet_pdf_params.corr_hmx_hmy_2.extent_int(0) != max_hydromet_dim ||
      hydromet_pdf_params.corr_hmx_hmy_2.extent_int(1) != max_hydromet_dim)
    throw std::invalid_argument("init_hydromet_pdf_params: fixed field corr_hmx_hmy_2");

  // Initialize hydromet_pdf_params in Native field order.
  Kokkos::deep_copy(hydromet_pdf_params.hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.mu_hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.mu_hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.sigma_hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.sigma_hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_w_hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_w_hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_chi_hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_chi_hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_eta_hm_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_eta_hm_2, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_hmx_hmy_1, zero);
  Kokkos::deep_copy(hydromet_pdf_params.corr_hmx_hmy_2, zero);
  hydromet_pdf_params.mu_Ncn_1 = zero;
  hydromet_pdf_params.mu_Ncn_2 = zero;
  hydromet_pdf_params.sigma_Ncn_1 = zero;
  hydromet_pdf_params.sigma_Ncn_2 = zero;
}

// =============================================================================
void init_precip_fracs_api(
    int nzt, int ngrdcol, // In
    precipitation_fractions& precip_fracs // Out
) {
  // Description:
  // Zero the precipitation-fraction fields in precip_fracs.
  // References:
  //-----------------------------------------------------------------------
  // nzt: Number of thermodynamic grid levels [-].
  // ngrdcol: Number of grid columns [-].
  // precip_fracs: Precipitation fractions [-].

  // Native IntentOut releases the caller object before allocation.
  precip_fracs = {};

  // Allocate precip frac arrays.
  precip_fracs.precip_frac =
      decltype(precip_fracs.precip_frac)("precip_fracs::precip_frac", ngrdcol, nzt);
  precip_fracs.precip_frac_1 =
      decltype(precip_fracs.precip_frac_1)("precip_fracs::precip_frac_1", ngrdcol, nzt);
  precip_fracs.precip_frac_2 =
      decltype(precip_fracs.precip_frac_2)("precip_fracs::precip_frac_2", ngrdcol, nzt);

  // Set metadata values.
  precip_fracs.ngrdcol = ngrdcol;
  precip_fracs.nzt = nzt;

  // Zero precip_fracs.
  zero_precip_fracs_api(
      precip_fracs // InOut
  );
}

// =============================================================================
void zero_precip_fracs_api(
    precipitation_fractions& precip_fracs // InOut
) {
  // Description:
  // Zero the precipitation-fraction fields in precip_fracs.
  // References:
  //-----------------------------------------------------------------------
  // precip_fracs: Precipitation fractions [-].

  const auto precip_fracs_precip_frac = precip_fracs.precip_frac;
  const auto precip_fracs_precip_frac_1 = precip_fracs.precip_frac_1;
  const auto precip_fracs_precip_frac_2 = precip_fracs.precip_frac_2;
  const int ngrdcol = precip_fracs.ngrdcol;
  const int nzt = precip_fracs.nzt;

  // Initialize precip_fracs.
  Kokkos::parallel_for(
      "zero_precip_fracs_api", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        precip_fracs_precip_frac(i, k) = zero;
        precip_fracs_precip_frac_1(i, k) = zero;
        precip_fracs_precip_frac_2(i, k) = zero;
      });
}

} // namespace clubb_kokkos
