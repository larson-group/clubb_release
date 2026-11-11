// Declaration-only Native lower boundaries for this isolated draft.
// Root must add missing definitions through the existing f2k generator owner,
// then remove this staging header. These declarations are not bridge implementations.
#pragma once
#include "f2k/procedure_bridge.hpp"
namespace f2k {
// src/CLUBB_core/stats_netcdf.F90:stats_update_budget_2d; all arguments required.
template <class Array0>
void stats_update_budget(const std::string& name, const Array0& values, stats_type& stats);
} // namespace f2k
