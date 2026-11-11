// Extend bottom staggered-grid values linearly from the two levels above.
#pragma once
#include "CLUBB_core/array_types.hpp"
namespace clubb_kokkos {

//=============================================================================
KOKKOS_INLINE_FUNCTION real_t lin_ext_zm_bottom(real_t var_zmp2, real_t var_zmp1,          // In
                                                real_t zmp2, real_t zmp1, real_t stats_zm) // In
{
  // Description:
  //   This function computes the value of a momentum-level variable at a bottom
  //   grid level by using a linear extension of the values of the variable at
  //   the two levels immediately above the level where the result value is
  //   needed.

  // References:
  //   None
  //-----------------------------------------------------------------------

  // var hascallerunits; allaltitudes aremeters.
  const real_t var_zm = ((var_zmp2 - var_zmp1) / (zmp2 - zmp1)) * (stats_zm - zmp1) + var_zmp1;
  return var_zm;
}

//=============================================================================
KOKKOS_INLINE_FUNCTION real_t lin_ext_zt_bottom(real_t var_ztp2, real_t var_ztp1,          // In
                                                real_t ztp2, real_t ztp1, real_t stats_zt) // In
{
  // Description:
  //   This function computes the value of a thermodynamic-level variable at a
  //   bottom grid level by using a linear extension of the values of the
  //   variable at the two levels immediately above the level where the result
  //   value is needed.
  //
  // References:
  //   None
  //-----------------------------------------------------------------------

  // var hascallerunits; allaltitudes aremeters.
  const real_t var_zt = ((var_ztp2 - var_ztp1) / (ztp2 - ztp1)) * (stats_zt - ztp1) + var_ztp1;
  return var_zt;
}
} // namespace clubb_kokkos
