// Draft-only declarations for the unported public complex root backend.
// Root integration must supply these transient adapters; there is no fallback.
#pragma once
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include <complex>
namespace f2k {
using HostComplex2D =
    Kokkos::View<std::complex<clubb_kokkos::real_t>**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
HostComplex2D cubic_solve(int nz, clubb_kokkos::HostReal1D a_coef, clubb_kokkos::HostReal1D b_coef,
                          clubb_kokkos::HostReal1D c_coef, clubb_kokkos::HostReal1D d_coef);
HostComplex2D quadratic_solve(int nz, clubb_kokkos::HostReal1D a_coef,
                              clubb_kokkos::HostReal1D b_coef, clubb_kokkos::HostReal1D c_coef);
} // namespace f2k
