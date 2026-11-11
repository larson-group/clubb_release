#pragma once

namespace clubb_kokkos {
#if CLUBB_REAL_TYPE == 4
using real_t = float;
#elif CLUBB_REAL_TYPE == 8
using real_t = double;
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
}  // namespace clubb_kokkos
