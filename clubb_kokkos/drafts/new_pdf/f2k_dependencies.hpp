// Pending public complex-root adapter; declaration only, with no fallback.
// Root owns promotion into src/f2k or replacement by the calc_roots C++ port.
#pragma once
#include "array_types.hpp"
#include <complex>
namespace f2k {
using HostComplex2D =
    Kokkos::View<std::complex<clubb_kokkos::real_t>**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using HostReal1D = Kokkos::View<clubb_kokkos::real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
HostComplex2D cubic_solve(int nz, HostReal1D a_coef, HostReal1D b_coef, HostReal1D c_coef,
                          HostReal1D d_coef);
} // namespace f2k
