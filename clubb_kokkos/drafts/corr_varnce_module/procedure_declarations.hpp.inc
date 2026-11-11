// Merge inside clubb_kokkos after the generated metadata type and scalar types.
// Header includes required: array_types.hpp, err_info_type_module.hpp, <optional>.
// Metadata indices remain Native one-based; outputs retain caller allocations.
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
    int pdf_dim, const hm_metadata_type& hm_metadata,                 // In
    bool l_fix_w_chi_eta_correlations,                                // In
    Real2D corr_array_n_cloud, Real2D corr_array_n_below,             // Out
    std::optional<Real2D> corr_array_n_cloud_default = std::nullopt,  // Optional In
    std::optional<Real2D> corr_array_n_below_default = std::nullopt); // Optional In

void assert_corr_symmetric(Real2D corr_array_n,      // In
                           int pdf_dim,              // In
                           err_info_type& err_info); // InOut
void print_corr_matrix(int pdf_dim,                  // In
                       Real2D corr_array_n);         // In

// Existing shared array ownership adapts host/strided caller sections once.
template <class Cloud, class Below, class DefaultCloud = Real2D, class DefaultBelow = Real2D>
void setup_corr_varnce_array_api(
    int pdf_dim, const hm_metadata_type& hm_metadata,                      // In
    bool l_fix_w_chi_eta_correlations,                                     // In
    Cloud corr_array_n_cloud, Below corr_array_n_below,                    // Out
    std::optional<DefaultCloud> corr_array_n_cloud_default = std::nullopt, // Optional In
    std::optional<DefaultBelow> corr_array_n_below_default = std::nullopt) // Optional In
{
  // Keep caller storage while staging only the existing public array boundary.
  auto cloud_exec = to_execution_space(corr_array_n_cloud, false);
  auto below_exec = to_execution_space(corr_array_n_below, false);
  std::optional<Real2D> cloud_default_exec, below_default_exec;
  if (corr_array_n_cloud_default.has_value() && corr_array_n_below_default.has_value()) {
    cloud_default_exec = to_execution_space(*corr_array_n_cloud_default);
    below_default_exec = to_execution_space(*corr_array_n_below_default);
  }
  setup_corr_varnce_array_api(pdf_dim, hm_metadata, l_fix_w_chi_eta_correlations, cloud_exec,
                              below_exec, cloud_default_exec, below_default_exec);
  from_execution_space(corr_array_n_cloud, cloud_exec);
  from_execution_space(corr_array_n_below, below_exec);
}
template <class Matrix>
void assert_corr_symmetric(Matrix corr_array_n, int pdf_dim, err_info_type& err_info) {
  // Reuse the existing owner for host/strided input; error Views remain shared.
  assert_corr_symmetric(to_execution_space(corr_array_n), pdf_dim, err_info);
}
template <class Matrix> void print_corr_matrix(int pdf_dim, Matrix corr_array_n) {
  // Preserve the existing public array boundary and Native matrix orientation.
  print_corr_matrix(pdf_dim, to_execution_space(corr_array_n));
}
