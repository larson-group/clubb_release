// Missing public typed entries in the existing lapack_interfaces.hpp owner.
// Integration must extend f2k/lapack_bridge.hpp/.F90 raw-library shims; no
// generated CLUBB procedure bridge or raw Fortran symbol belongs in physics.
#pragma once
#include "clubb_real.hpp"
namespace clubb_kokkos {
void lapack_potrf(const char* uplo, int n, real_t* a, int lda, int& info);
void lapack_poequ(int n, const real_t* a, int lda, real_t* s, real_t& scond, real_t& amax,
                  int& info);
void lapack_laqsy(const char* uplo, int n, real_t* a, int lda, const real_t* s, real_t scond,
                  real_t amax, char& equed);
void lapack_syev(const char* jobz, const char* uplo, int n, real_t* a, int lda, real_t* w,
                 real_t* work, int lwork, int& info);
} // namespace clubb_kokkos
