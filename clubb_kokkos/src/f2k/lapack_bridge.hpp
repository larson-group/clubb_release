// BFB-only C++/Fortran ABI adaptation for LAPACK PDF matrix routines.
// Host buffers are contiguous in LAPACK layout. Keep character truncation,
// scalar references and all raw C declarations inside the transient bridge.
#pragma once

#include "clubb_real.hpp"

namespace f2k {
extern "C" {
void f2k_library_potrf_c(char, int, clubb_kokkos::real_t*, int, int*);
void f2k_library_poequ_c(int, const clubb_kokkos::real_t*, int, clubb_kokkos::real_t*, clubb_kokkos::real_t*, clubb_kokkos::real_t*, int*);
void f2k_library_laqsy_c(char, int, clubb_kokkos::real_t*, int, const clubb_kokkos::real_t*, clubb_kokkos::real_t, clubb_kokkos::real_t, char*);
void f2k_library_syev_c(char, char, int, clubb_kokkos::real_t*, int, clubb_kokkos::real_t*, clubb_kokkos::real_t*, int, int*);
}
} // namespace f2k

namespace f2k {
using clubb_kokkos::real_t;

inline void lapack_potrf(const char* uplo, int n, real_t* a, int lda, int& info) {
  // Forward caller-owned host buffers to the precision-selected library.
  f2k_library_potrf_c(uplo[0], n, a, lda, &info);
}

inline void lapack_poequ(int n, const real_t* a, int lda, real_t* s, real_t& scond,
                         real_t& amax, int& info) {
  // Return the library's diagonal scaling and status without retained state.
  f2k_library_poequ_c(n, a, lda, s, &scond, &amax, &info);
}

inline void lapack_laqsy(const char* uplo, int n, real_t* a, int lda, const real_t* s,
                         real_t scond, real_t amax, char& equed) {
  // Preserve LAPACK's in-place equilibration and single-character result.
  f2k_library_laqsy_c(uplo[0], n, a, lda, s, scond, amax, &equed);
}

inline void lapack_syev(const char* jobz, const char* uplo, int n, real_t* a, int lda,
                        real_t* w, real_t* work, int lwork, int& info) {
  // The caller owns the symmetric matrix, eigenvalues and workspace.
  f2k_library_syev_c(jobz[0], uplo[0], n, a, lda, w, work, lwork, &info);
}

} // namespace f2k
