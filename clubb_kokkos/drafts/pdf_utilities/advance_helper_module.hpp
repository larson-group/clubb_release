// Smoothed scalar maximum from src/CLUBB_core/advance_helper_module.F90.
//
// This draft contains only smooth_max_scalars and the scalar smooth_max generic
// needed by pdf_utilities::smooth_corr_quotient. Other native routines remain
// outside this partial header. Root owns merging it into the permanent module
// in native procedure order. This draft has not been compiled or tested.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t smooth_max_scalars(real_t input_var1, real_t input_var2,
                                                 real_t smth_coef) {
  // Compute a smoothed maximum of two scalar values, returning the same units
  // as the inputs. smth_coef controls the smoothing and should have a similar
  // magnitude to input_var1 and input_var2. Reference: CLUBB ticket 965.
  real_t output_var; // Smoothed maximum [input units]

  //----------------------------- Begin Code -----------------------------
  output_var =
      one_half *
      ((input_var1 + input_var2) +
       Kokkos::sqrt((input_var1 - input_var2) * (input_var1 - input_var2) + smth_coef * smth_coef));
  return output_var;
}

// C++ overload dispatch for the native scalar smooth_max generic.
KOKKOS_INLINE_FUNCTION real_t smooth_max(real_t input_var1, real_t input_var2, real_t smth_coef) {
  // Select the native scalar smoothed-maximum formula.
  return smooth_max_scalars(input_var1, input_var2, smth_coef);
}

} // namespace clubb_kokkos
