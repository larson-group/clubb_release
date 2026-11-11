! calc_roots_wrapper.F90 — wrappers for module calc_roots

subroutine f2py_cubic_solve(nz, a_coef, b_coef, c_coef, d_coef, roots_real, roots_imag)

  use clubb_precision, only: core_rknd
  use calc_roots, only: cubic_solve

  implicit none

  integer, intent(in) :: nz
  real(kind=core_rknd), dimension(nz), intent(in) :: a_coef, b_coef, c_coef, d_coef
  real(kind=core_rknd), dimension(nz, 3), intent(out) :: roots_real, roots_imag

  complex(kind=core_rknd), dimension(nz, 3) :: roots

  roots = cubic_solve(nz, a_coef, b_coef, c_coef, d_coef)
  roots_real = real(roots, kind=core_rknd)
  roots_imag = aimag(roots)

end subroutine f2py_cubic_solve


subroutine f2py_quadratic_solve(nz, a_coef, b_coef, c_coef, roots_real, roots_imag)

  use clubb_precision, only: core_rknd
  use calc_roots, only: quadratic_solve

  implicit none

  integer, intent(in) :: nz
  real(kind=core_rknd), dimension(nz), intent(in) :: a_coef, b_coef, c_coef
  real(kind=core_rknd), dimension(nz, 2), intent(out) :: roots_real, roots_imag

  complex(kind=core_rknd), dimension(nz, 2) :: roots

  roots = quadratic_solve(nz, a_coef, b_coef, c_coef)
  roots_real = real(roots, kind=core_rknd)
  roots_imag = aimag(roots)

end subroutine f2py_quadratic_solve
