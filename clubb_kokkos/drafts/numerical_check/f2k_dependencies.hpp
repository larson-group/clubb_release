// Declaration-only pending public Native lower calls for this isolated draft.
// Add implementations through the existing transient owner, then remove this
// staging header. No private bindings or global metadata accessors are added.
#pragma once
#include "f2k/procedure_bridge.hpp"
namespace f2k {
// Native error_code.F90:69–94; bare no-argument subroutine declaration.
// The current parser requires parentheses and cannot read this exact signature.
void initialize_error_headers();
// Native parameters_tunable.F90:check_parameters_api; all arguments required.
template <class Params>
void check_parameters_api(int ngrdcol, const Params& clubb_params,
                          clubb_kokkos::err_info_type& err_info);
} // namespace f2k
