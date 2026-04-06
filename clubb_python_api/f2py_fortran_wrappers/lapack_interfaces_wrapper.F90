subroutine f2py_lapack_gtsv(n, nrhs, dl, d, du, b, soln, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_gtsv

  implicit none

  integer, intent(in) :: n, nrhs
  real(core_rknd), dimension(n-1), intent(in) :: dl, du
  real(core_rknd), dimension(n), intent(in) :: d
  real(core_rknd), dimension(n, nrhs), intent(in) :: b
  real(core_rknd), dimension(n, nrhs), intent(out) :: soln
  integer, intent(out) :: info

  real(core_rknd), dimension(n-1) :: dl_local, du_local
  real(core_rknd), dimension(n) :: d_local

  dl_local = dl
  d_local = d
  du_local = du
  soln = b

  call lapack_gtsv(n, nrhs, dl_local, d_local, du_local, soln, n, info)

end subroutine f2py_lapack_gtsv

subroutine f2py_lapack_gtsvx(n, nrhs, dl, d, du, b, soln, rcond, ferr, berr, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_gtsvx

  implicit none

  integer, intent(in) :: n, nrhs
  real(core_rknd), dimension(n-1), intent(in) :: dl, du
  real(core_rknd), dimension(n), intent(in) :: d
  real(core_rknd), dimension(n, nrhs), intent(in) :: b
  real(core_rknd), dimension(n, nrhs), intent(out) :: soln
  real(core_rknd), intent(out) :: rcond
  real(core_rknd), dimension(nrhs), intent(out) :: ferr, berr
  integer, intent(out) :: info

  real(core_rknd), dimension(n-1) :: dl_local, du_local, dlf, duf
  real(core_rknd), dimension(n) :: d_local, df
  real(core_rknd), dimension(max(1, n-2)) :: du2
  integer, dimension(n) :: ipiv, iwork
  real(core_rknd), dimension(3*n) :: work

  dl_local = dl
  d_local = d
  du_local = du

  call lapack_gtsvx("Not Factored", "No Transpose lhs", n, nrhs, &
                    dl_local, d_local, du_local, dlf, df, duf, du2, ipiv, &
                    b, n, soln, n, rcond, ferr, berr, work, iwork, info)

end subroutine f2py_lapack_gtsvx

subroutine f2py_lapack_gbsv(n, kl, ku, nrhs, ab, b, soln, ipiv, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_gbsv

  implicit none

  integer, intent(in) :: n, kl, ku, nrhs
  real(core_rknd), dimension(2*kl+ku+1, n), intent(in) :: ab
  real(core_rknd), dimension(n, nrhs), intent(in) :: b
  real(core_rknd), dimension(n, nrhs), intent(out) :: soln
  integer, dimension(n), intent(out) :: ipiv
  integer, intent(out) :: info

  real(core_rknd), dimension(2*kl+ku+1, n) :: ab_local

  ab_local = ab
  soln = b

  call lapack_gbsv(n, kl, ku, nrhs, ab_local, 2*kl+ku+1, ipiv, soln, n, info)

end subroutine f2py_lapack_gbsv

subroutine f2py_lapack_gbsvx(n, kl, ku, nrhs, ab, b, soln, rcond, ferr, berr, equed, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_gbsvx

  implicit none

  integer, intent(in) :: n, kl, ku, nrhs
  real(core_rknd), dimension(kl+ku+1, n), intent(in) :: ab
  real(core_rknd), dimension(n, nrhs), intent(in) :: b
  real(core_rknd), dimension(n, nrhs), intent(out) :: soln
  real(core_rknd), intent(out) :: rcond
  real(core_rknd), dimension(nrhs), intent(out) :: ferr, berr
  character(len=1), intent(out) :: equed
  integer, intent(out) :: info

  real(core_rknd), dimension(kl+ku+1, n) :: ab_local
  real(core_rknd), dimension(2*kl+ku+1, n) :: afb
  integer, dimension(n) :: ipiv, iwork
  real(core_rknd), dimension(n) :: r, c
  real(core_rknd), dimension(3*n) :: work

  ab_local = ab

  call lapack_gbsvx("Equilibrate lhs", "No Transpose lhs", n, kl, ku, nrhs, &
                    ab_local, kl+ku+1, afb, 2*kl+ku+1, ipiv, equed, r, c, &
                    b, n, soln, n, rcond, ferr, berr, work, iwork, info)

end subroutine f2py_lapack_gbsvx

subroutine f2py_lapack_potrf(uplo, n, a, a_factor, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_potrf

  implicit none

  character(len=*), intent(in) :: uplo
  integer, intent(in) :: n
  real(core_rknd), dimension(n, n), intent(in) :: a
  real(core_rknd), dimension(n, n), intent(out) :: a_factor
  integer, intent(out) :: info

  a_factor = a
  call lapack_potrf(uplo, n, a_factor, n, info)

end subroutine f2py_lapack_potrf

subroutine f2py_lapack_poequ(n, a, s, scond, amax, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_poequ

  implicit none

  integer, intent(in) :: n
  real(core_rknd), dimension(n, n), intent(in) :: a
  real(core_rknd), dimension(n), intent(out) :: s
  real(core_rknd), intent(out) :: scond, amax
  integer, intent(out) :: info

  call lapack_poequ(n, a, n, s, scond, amax, info)

end subroutine f2py_lapack_poequ

subroutine f2py_lapack_laqsy(uplo, n, a, s, scond, amax, a_equed, equed)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_laqsy

  implicit none

  character(len=*), intent(in) :: uplo
  integer, intent(in) :: n
  real(core_rknd), dimension(n, n), intent(in) :: a
  real(core_rknd), dimension(n), intent(in) :: s
  real(core_rknd), intent(in) :: scond, amax
  real(core_rknd), dimension(n, n), intent(out) :: a_equed
  character(len=1), intent(out) :: equed

  a_equed = a
  call lapack_laqsy(uplo, n, a_equed, n, s, scond, amax, equed)

end subroutine f2py_lapack_laqsy

subroutine f2py_lapack_syev(jobz, uplo, n, a, a_out, w, info)

  use clubb_precision, only: core_rknd
  use lapack_interfaces, only: lapack_syev

  implicit none

  character(len=*), intent(in) :: jobz, uplo
  integer, intent(in) :: n
  real(core_rknd), dimension(n, n), intent(in) :: a
  real(core_rknd), dimension(n, n), intent(out) :: a_out
  real(core_rknd), dimension(n), intent(out) :: w
  integer, intent(out) :: info

  integer :: lwork
  real(core_rknd), dimension(max(1, 3*n-1)) :: work

  a_out = a
  lwork = max(1, 3*n-1)
  call lapack_syev(jobz, uplo, n, a_out, n, w, work, lwork, info)

end subroutine f2py_lapack_syev
