// BFB-only C++ entry points to LAPACK PDF matrix factorization and diagnostics.
//
// Cholesky/scaling/eigenvalue entries accept contiguous column-major host
// buffers in core real precision. The f2k shim owns the Fortran library ABI.
#pragma once

#include "clubb_real.hpp"
#include "f2k/lapack_bridge.hpp"

namespace clubb_kokkos {

inline void lapack_potrf(const char* uplo, int n, real_t* a, int lda, int& info) {
  // Factor the selected triangle of a caller-owned contiguous host matrix.
  f2k::lapack_potrf(uplo, n, a, lda, info);
}

inline void lapack_poequ(int n, const real_t* a, int lda, real_t* s, real_t& scond,
                         real_t& amax, int& info) {
  // Compute the Native diagonal scaling factors and condition ratio.
  f2k::lapack_poequ(n, a, lda, s, scond, amax, info);
}

inline void lapack_laqsy(const char* uplo, int n, real_t* a, int lda, const real_t* s,
                         real_t scond, real_t amax, char& equed) {
  // Equilibrate the symmetric matrix in place when the library selects scaling.
  f2k::lapack_laqsy(uplo, n, a, lda, s, scond, amax, equed);
}

inline void lapack_syev(const char* jobz, const char* uplo, int n, real_t* a, int lda,
                        real_t* w, real_t* work, int lwork, int& info) {
  // Compute the symmetric eigenvalues used by Native matrix diagnostics.
  f2k::lapack_syev(jobz, uplo, n, a, lda, w, work, lwork, info);
}

} // namespace clubb_kokkos
