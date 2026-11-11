// Load prescribed normal-space correlations in the Native metadata ordering.
#pragma once
#include "input_reader.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/matrix_operations.hpp"
#include <iostream>
#include <stdexcept>
namespace clubb_kokkos {

//=============================================================================
template <class Matrix>
void read_correlation_matrix(ScratchWorkspace& workspace,                      // InOut
                             int iunit, const std::string& input_file,         // In
                             int pdf_dim, const hm_metadata_type& hm_metadata, // In
                             Matrix corr_array_n)                              // InOut
{
  // Description:
  //   Reads a correlation variance array from a file and stores it in an array.
  //-----------------------------------------------------------------------------

  auto workspace_scope = workspace.scope();
  const int nCols      = count_columns(iunit, input_file);
  // Allocate all arrays based on pdf_dim.
  std::vector<one_dim_read_var> retVars(nCols);
  // Initializing tozero assumes unspecified correlations arezero.
  Kokkos::deep_copy(corr_array_n, real_t(0.0));
  // Set main diagonal toone.
  Kokkos::parallel_for(
      "read_correlation_matrix::diagonal",
      Kokkos::RangePolicy<typename Matrix::execution_space>(0, pdf_dim),
      KOKKOS_LAMBDA(int i) { corr_array_n(i, i) = real_t(1.0); });
  // Read values from the specified file.
  read_one_dim_file(iunit, nCols, input_file, // In
                    retVars);                 // Out
  if (static_cast<int>(retVars[0].values.extent(0)) != nCols) {
    std::cerr << "Correlation matrix must have an equal number of rows and cols in file "
              << input_file << '\n';
    throw std::runtime_error("Bad data in correlation file.");
  }
  // Skip the first coordinate/diagonal column, as in the Native file reader.
  for (int i = 1; i < nCols; ++i) {
    const int var_index1 = ((get_corr_var_index(retVars[i].name, hm_metadata) + 1)) - 1;
    if (var_index1 > -2) {
      // Metadata lookup stays onhost; stage this single inputprofile once for
      // its scalar coefficients. Matrixupdates use the canonical C++owner.
      auto profile_scope = workspace.scope();
      const auto values  = to_host_space(workspace, retVars[i].values);
      for (int j = 0; j < i; ++j) {
        const int var_index2 = ((get_corr_var_index(retVars[j].name, hm_metadata) + 1)) - 1;
        if (var_index2 > -2) {
          set_lower_triangular_matrix(pdf_dim, var_index1, var_index2, values(j), // In
                                      corr_array_n);                              // InOut
        }
      }
    }
  }
  deallocate_one_dim_vars(nCols,    // In
                          retVars); // InOut
}
} // namespace clubb_kokkos
