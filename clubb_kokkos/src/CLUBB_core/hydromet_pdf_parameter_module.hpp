#pragma once

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"

#include "array_types.hpp"
#include <optional>
#include <type_traits>

namespace clubb_kokkos {

// BEGIN GENERATED HYDROMET PDF TYPES
// Mirrors src/CLUBB_core/hydromet_pdf_parameter_module.F90; arrays use Fortran-order LayoutLeft.
struct hydromet_pdf_parameter {
  Kokkos::View<real_t*, Kokkos::LayoutLeft> hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> hm_2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> mu_hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> mu_hm_2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> sigma_hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> sigma_hm_2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_w_hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_w_hm_2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_chi_hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_chi_hm_2;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_eta_hm_1;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> corr_eta_hm_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_hmx_hmy_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_hmx_hmy_2;
  real_t mu_Ncn_1{};
  real_t mu_Ncn_2{};
  real_t sigma_Ncn_1{};
  real_t sigma_Ncn_2{};
};

// Device storage for Fortran hydromet_pdf_params(column, level).
// Each component is one View so adjacent columns stay contiguous.
struct hydromet_pdf_parameter_collection {
  static constexpr int max_hydromet_dim = 8;

  int ngrdcol{};
  int nzt{};
  Kokkos::View<real_t***, Kokkos::LayoutLeft> hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> hm_2;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> mu_hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> mu_hm_2;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sigma_hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sigma_hm_2;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_w_hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_w_hm_2;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_chi_hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_chi_hm_2;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_eta_hm_1;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> corr_eta_hm_2;
  Kokkos::View<real_t****, Kokkos::LayoutLeft> corr_hmx_hmy_1;
  Kokkos::View<real_t****, Kokkos::LayoutLeft> corr_hmx_hmy_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> mu_Ncn_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> mu_Ncn_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sigma_Ncn_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sigma_Ncn_2;

  void allocate(int columns, int levels) {
    // Allocate one View per Fortran component with columns contiguous.
    ngrdcol = columns;
    nzt = levels;
    hm_1 = decltype(hm_1)(
        "hydromet_hm_1", columns, levels, max_hydromet_dim);
    hm_2 = decltype(hm_2)(
        "hydromet_hm_2", columns, levels, max_hydromet_dim);
    mu_hm_1 = decltype(mu_hm_1)(
        "hydromet_mu_hm_1", columns, levels, max_hydromet_dim);
    mu_hm_2 = decltype(mu_hm_2)(
        "hydromet_mu_hm_2", columns, levels, max_hydromet_dim);
    sigma_hm_1 = decltype(sigma_hm_1)(
        "hydromet_sigma_hm_1", columns, levels, max_hydromet_dim);
    sigma_hm_2 = decltype(sigma_hm_2)(
        "hydromet_sigma_hm_2", columns, levels, max_hydromet_dim);
    corr_w_hm_1 = decltype(corr_w_hm_1)(
        "hydromet_corr_w_hm_1", columns, levels, max_hydromet_dim);
    corr_w_hm_2 = decltype(corr_w_hm_2)(
        "hydromet_corr_w_hm_2", columns, levels, max_hydromet_dim);
    corr_chi_hm_1 = decltype(corr_chi_hm_1)(
        "hydromet_corr_chi_hm_1", columns, levels, max_hydromet_dim);
    corr_chi_hm_2 = decltype(corr_chi_hm_2)(
        "hydromet_corr_chi_hm_2", columns, levels, max_hydromet_dim);
    corr_eta_hm_1 = decltype(corr_eta_hm_1)(
        "hydromet_corr_eta_hm_1", columns, levels, max_hydromet_dim);
    corr_eta_hm_2 = decltype(corr_eta_hm_2)(
        "hydromet_corr_eta_hm_2", columns, levels, max_hydromet_dim);
    corr_hmx_hmy_1 = decltype(corr_hmx_hmy_1)(
        "hydromet_corr_hmx_hmy_1", columns, levels, max_hydromet_dim, max_hydromet_dim);
    corr_hmx_hmy_2 = decltype(corr_hmx_hmy_2)(
        "hydromet_corr_hmx_hmy_2", columns, levels, max_hydromet_dim, max_hydromet_dim);
    mu_Ncn_1 = decltype(mu_Ncn_1)(
        "hydromet_mu_Ncn_1", columns, levels);
    mu_Ncn_2 = decltype(mu_Ncn_2)(
        "hydromet_mu_Ncn_2", columns, levels);
    sigma_Ncn_1 = decltype(sigma_Ncn_1)(
        "hydromet_sigma_Ncn_1", columns, levels);
    sigma_Ncn_2 = decltype(sigma_Ncn_2)(
        "hydromet_sigma_Ncn_2", columns, levels);
  }

  void zero() {
    // Reset every component before the next logical runtime batch.
    Kokkos::deep_copy(hm_1, real_t{});
    Kokkos::deep_copy(hm_2, real_t{});
    Kokkos::deep_copy(mu_hm_1, real_t{});
    Kokkos::deep_copy(mu_hm_2, real_t{});
    Kokkos::deep_copy(sigma_hm_1, real_t{});
    Kokkos::deep_copy(sigma_hm_2, real_t{});
    Kokkos::deep_copy(corr_w_hm_1, real_t{});
    Kokkos::deep_copy(corr_w_hm_2, real_t{});
    Kokkos::deep_copy(corr_chi_hm_1, real_t{});
    Kokkos::deep_copy(corr_chi_hm_2, real_t{});
    Kokkos::deep_copy(corr_eta_hm_1, real_t{});
    Kokkos::deep_copy(corr_eta_hm_2, real_t{});
    Kokkos::deep_copy(corr_hmx_hmy_1, real_t{});
    Kokkos::deep_copy(corr_hmx_hmy_2, real_t{});
    Kokkos::deep_copy(mu_Ncn_1, real_t{});
    Kokkos::deep_copy(mu_Ncn_2, real_t{});
    Kokkos::deep_copy(sigma_Ncn_1, real_t{});
    Kokkos::deep_copy(sigma_Ncn_2, real_t{});
  }
};

// Mirrors src/CLUBB_core/hydromet_pdf_parameter_module.F90; arrays use Fortran-order LayoutLeft.
struct precipitation_fractions {
  int ngrdcol{};
  int nzt{};
  Kokkos::View<real_t**, Kokkos::LayoutLeft> precip_frac;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> precip_frac_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> precip_frac_2;
};

// END GENERATED HYDROMET PDF TYPES

void init_hydromet_pdf_params(
    hydromet_pdf_parameter& hydromet_pdf_params // InOut
);
void init_precip_fracs_api(
    int nzt, int ngrdcol, // In
    precipitation_fractions& precip_fracs // Out
);
void zero_precip_fracs_api(
    precipitation_fractions& precip_fracs // InOut
);
// Native repeated per-cell initialization through the existing SoA owner.
inline void init_hydromet_pdf_params(
    hydromet_pdf_parameter_collection& hydromet_pdf_params // InOut
) {
  // Reset the existing case-owned collection; its schema owns all 18 fields.
  hydromet_pdf_params.zero();
}


}  // namespace clubb_kokkos
