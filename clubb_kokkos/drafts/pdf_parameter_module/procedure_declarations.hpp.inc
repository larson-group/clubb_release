// Merge inside clubb_kokkos after the generated types.
inline constexpr int num_pdf_params = 47; // Native packed field count.
void init_pdf_params(
    int nz, int ngrdcol, // In
    pdf_parameter& pdf_params // Out
);
void zero_pdf_params_api(
    pdf_parameter& pdf_params // InOut
);
void init_pdf_implicit_coefs_terms_api(
    int nz, int ngrdcol, int sclr_dim, // In
    implicit_coefs_terms& pdf_implicit_coefs_terms // Out
);
void zero_pdf_implicit_coefs_terms_api(
    implicit_coefs_terms& pdf_implicit_coefs_terms // InOut
);
void print_pdf_params(
    const pdf_parameter& pdf_params, int ngrdcol // In
);
void pack_pdf_params_api(
    const pdf_parameter& pdf_params, int nz, // In
    Real2D r_param_array, // Out
    std::optional<int> k_start_in = std::nullopt, // Optional In
    std::optional<int> k_end_in = std::nullopt // Optional In
);
void unpack_pdf_params_api(
    Real2D r_param_array, int nz, // In
    pdf_parameter& pdf_params, // InOut
    std::optional<int> k_start_in = std::nullopt, // Optional In
    std::optional<int> k_end_in = std::nullopt // Optional In
);
template <class Buffer>
inline void pack_pdf_params_api(
    const pdf_parameter& pdf_params, int nz, // In
    Buffer r_param_array, // Out
    std::optional<int> k_start_in = std::nullopt, // Optional In
    std::optional<int> k_end_in = std::nullopt // Optional In
) {
  // Preserve host/strided packet ownership through the existing array owner.
  const auto buffer_exec = to_execution_space(r_param_array, false);
  pack_pdf_params_api(
      pdf_params, nz, // In
      buffer_exec, // Out
      k_start_in, k_end_in // Optional In
  );
  from_execution_space(r_param_array, buffer_exec);
}
template <class Buffer>
inline void unpack_pdf_params_api(
    Buffer r_param_array, int nz, // In
    pdf_parameter& pdf_params, // InOut
    std::optional<int> k_start_in = std::nullopt, // Optional In
    std::optional<int> k_end_in = std::nullopt // Optional In
) {
  // Preserve host/strided packet ownership through the existing array owner.
  const auto buffer_exec = to_execution_space(r_param_array);
  unpack_pdf_params_api(
      buffer_exec, nz, // In
      pdf_params, // InOut
      k_start_in, k_end_in // Optional In
  );
}
