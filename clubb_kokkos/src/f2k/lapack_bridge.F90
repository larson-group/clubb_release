! C bindings for LAPACK PDF matrix factorization and diagnostics.
! These routines only adapt the library ABI; they do not call CLUBB's solver
! wrappers or retain any state. C++ owns packing, workspaces and error handling.
module f2k_lapack_bridge
  use, intrinsic :: iso_c_binding, only: c_int, c_char
  use clubb_precision, only: core_rknd
  implicit none
  private
  public :: f2k_library_potrf_c, f2k_library_poequ_c
  public :: f2k_library_laqsy_c, f2k_library_syev_c
contains

  subroutine f2k_library_potrf_c(uplo, n, a, lda, info) bind(C)
    ! Description: Call SPOTRF/DPOTRF on the caller's selected matrix triangle.
    character(kind=c_char), value, intent(in) :: uplo ! Selected triangle
    integer(c_int), value, intent(in) :: n, lda ! Matrix order and leading dimension
    real(core_rknd), intent(inout) :: a(lda,n) ! Contiguous host matrix
    integer(c_int), intent(out) :: info ! External library status
#if CLUBB_REAL_TYPE == 8
    call dpotrf(uplo, n, a, lda, info)
#elif CLUBB_REAL_TYPE == 4
    call spotrf(uplo, n, a, lda, info)
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
  end subroutine f2k_library_potrf_c

  subroutine f2k_library_poequ_c(n, a, lda, s, scond, amax, info) bind(C)
    ! Description: Call SPOEQU/DPOEQU and return diagonal equilibration data.
    integer(c_int), value, intent(in) :: n, lda ! Matrix order and leading dimension
    real(core_rknd), intent(in) :: a(lda,n) ! Contiguous host matrix
    real(core_rknd), intent(out) :: s(n), scond, amax ! Scaling, ratio and maximum
    integer(c_int), intent(out) :: info ! External library status
#if CLUBB_REAL_TYPE == 8
    call dpoequ(n, a, lda, s, scond, amax, info)
#elif CLUBB_REAL_TYPE == 4
    call spoequ(n, a, lda, s, scond, amax, info)
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
  end subroutine f2k_library_poequ_c

  subroutine f2k_library_laqsy_c(uplo, n, a, lda, s, scond, amax, equed) bind(C)
    ! Description: Call SLAQSY/DLAQSY with Native equilibration selectors.
    character(kind=c_char), value, intent(in) :: uplo ! Selected triangle
    integer(c_int), value, intent(in) :: n, lda ! Matrix order and leading dimension
    real(core_rknd), intent(inout) :: a(lda,n) ! Contiguous host matrix
    real(core_rknd), intent(in) :: s(n) ! Diagonal scaling
    real(core_rknd), value, intent(in) :: scond, amax ! Scaling ratio and maximum
    character(kind=c_char), intent(out) :: equed ! Whether scaling was applied
#if CLUBB_REAL_TYPE == 8
    call dlaqsy(uplo, n, a, lda, s, scond, amax, equed)
#elif CLUBB_REAL_TYPE == 4
    call slaqsy(uplo, n, a, lda, s, scond, amax, equed)
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
  end subroutine f2k_library_laqsy_c

  subroutine f2k_library_syev_c(jobz, uplo, n, a, lda, w, work, lwork, info) bind(C)
    ! Description: Call SSYEV/DSYEV for Native correlation-matrix diagnostics.
    character(kind=c_char), value, intent(in) :: jobz, uplo ! Eigenvector/triangle modes
    integer(c_int), value, intent(in) :: n, lda, lwork ! Matrix/workspace dimensions
    real(core_rknd), intent(inout) :: a(lda,n), work(lwork) ! Matrix and library scratch
    real(core_rknd), intent(out) :: w(n) ! Eigenvalues
    integer(c_int), intent(out) :: info ! External library status
#if CLUBB_REAL_TYPE == 8
    call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
#elif CLUBB_REAL_TYPE == 4
    call ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
#else
#error Unsupported CLUBB_REAL_TYPE
#endif
  end subroutine f2k_library_syev_c
end module f2k_lapack_bridge
