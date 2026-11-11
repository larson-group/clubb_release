#pragma once

#include "scratch_workspace.hpp"

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"

#include "array_types.hpp"
#include <optional>
#include <type_traits>

namespace clubb_kokkos {

// BEGIN GENERATED PDF PARAMETER TYPES
// Mirrors src/CLUBB_core/pdf_parameter_module.F90; arrays use Fortran-order LayoutLeft.
struct pdf_parameter {
  int ngrdcol{};
  int nz{};
  Kokkos::View<real_t**, Kokkos::LayoutLeft> w_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> w_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_w_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_w_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rt_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rt_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_rt_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_rt_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_thl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> varnce_thl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_rt_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_rt_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_thl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_thl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_rt_thl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_rt_thl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> alpha_thl;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> alpha_rt;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> crt_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> crt_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> cthl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> cthl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> chi_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> chi_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stdev_chi_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stdev_chi_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stdev_eta_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stdev_eta_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> covar_chi_eta_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> covar_chi_eta_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_chi_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_chi_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_eta_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_w_eta_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_chi_eta_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> corr_chi_eta_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rsatl_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rsatl_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rc_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rc_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> cloud_frac_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> cloud_frac_2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> mixt_frac;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ice_supersat_frac_1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> ice_supersat_frac_2;
};

// Mirrors src/CLUBB_core/pdf_parameter_module.F90; arrays use Fortran-order LayoutLeft.
struct implicit_coefs_terms {
  int ngrdcol{};
  int nz{};
  int sclr_dim{};
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp4_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2rtp_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2rtp_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2thlp_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2thlp_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2up_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2up_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wp2vp_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wp2vp_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtp2_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtp2_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpthlp2_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpthlp2_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtpthlp_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtpthlp_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpup2_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpup2_explicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpvp2_implicit;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpvp2_explicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> coef_wp2sclrp_implicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> term_wp2sclrp_explicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> coef_wpsclrp2_implicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> term_wpsclrp2_explicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> coef_wprtpsclrp_implicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> term_wprtpsclrp_explicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> coef_wpthlpsclrp_implicit;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> term_wpthlpsclrp_explicit;
};

// END GENERATED PDF PARAMETER TYPES

inline constexpr int num_pdf_params = 47;      // Native packed field count.
void init_pdf_params(int nz, int ngrdcol,      // In
                     pdf_parameter& pdf_params // Out
);
void zero_pdf_params_api(pdf_parameter& pdf_params // InOut
);
void init_pdf_implicit_coefs_terms_api(int nz, int ngrdcol, int sclr_dim,             // In
                                       implicit_coefs_terms& pdf_implicit_coefs_terms // Out
);
void zero_pdf_implicit_coefs_terms_api(implicit_coefs_terms& pdf_implicit_coefs_terms // InOut
);
void print_pdf_params(ScratchWorkspace& workspace,                 // InOut
                      const pdf_parameter& pdf_params, int ngrdcol // In
);
void pack_pdf_params_api(const pdf_parameter& pdf_params, int nz,      // In
                         Real2D r_param_array,                         // Out
                         std::optional<int> k_start_in = std::nullopt, // Optional In
                         std::optional<int> k_end_in   = std::nullopt  // Optional In
);
void unpack_pdf_params_api(Real2D r_param_array, int nz,                 // In
                           pdf_parameter& pdf_params,                    // InOut
                           std::optional<int> k_start_in = std::nullopt, // Optional In
                           std::optional<int> k_end_in   = std::nullopt  // Optional In
);
template <class Buffer>
inline void pack_pdf_params_api(const pdf_parameter& pdf_params, int nz,      // In
                                Buffer r_param_array,                         // Out
                                std::optional<int> k_start_in = std::nullopt, // Optional In
                                std::optional<int> k_end_in   = std::nullopt  // Optional In
) {
  // Preserve host/strided packet ownership through the existing array owner.
  const auto buffer_exec = to_execution_space(r_param_array, false);
  pack_pdf_params_api(pdf_params, nz,      // In
                      buffer_exec,         // Out
                      k_start_in, k_end_in // Optional In
  );
  from_execution_space(r_param_array, buffer_exec);
}
template <class Buffer>
inline void unpack_pdf_params_api(Buffer r_param_array, int nz,                 // In
                                  pdf_parameter& pdf_params,                    // InOut
                                  std::optional<int> k_start_in = std::nullopt, // Optional In
                                  std::optional<int> k_end_in   = std::nullopt  // Optional In
) {
  // Preserve host/strided packet ownership through the existing array owner.
  const auto buffer_exec = to_execution_space(r_param_array);
  unpack_pdf_params_api(buffer_exec, nz,     // In
                        pdf_params,          // InOut
                        k_start_in, k_end_in // Optional In
  );
}

} // namespace clubb_kokkos
