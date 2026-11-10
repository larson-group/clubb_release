! calc_roots_wrapper.F90 — wrappers for module calc_roots

subroutine f2py_cubic_solve(nz, a_coef, b_coef, c_coef, d_coef, roots_real, roots_imag)

  ! Description:
  !   Split the elemental polynomial roots into real and imaginary F2PY arrays.

  use clubb_precision, only: core_rknd
  use calc_roots, only: cubic_solve

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of coefficient tuples [-]
  ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0    [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: a_coef
  ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0    [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: b_coef
  ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0      [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: c_coef
  ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0             [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: d_coef

  !--------------------------- Output Variables ---------------------------
  real(kind=core_rknd), dimension(nz, 3), intent(out) :: roots_real ! Real parts of the roots [-]
  ! Imaginary parts of the roots [-]
  real(kind=core_rknd), dimension(nz, 3), intent(out) :: roots_imag

  !--------------------------- Local Variables ---------------------------
  complex(kind=core_rknd), dimension(nz, 3) :: roots ! Complex roots before splitting for F2PY [-]

  !--------------------------- Begin Code ---------------------------

  call cubic_solve( a_coef, b_coef, c_coef, d_coef, & ! In
                    roots(:,1), roots(:,2), & ! Out
                    roots(:,3) ) ! Out
  roots_real = real(roots, kind=core_rknd)
  roots_imag = aimag(roots)

end subroutine f2py_cubic_solve

subroutine f2py_quadratic_solve(nz, a_coef, b_coef, c_coef, roots_real, roots_imag)

  ! Description:
  !   Split the elemental polynomial roots into real and imaginary F2PY arrays.

  use clubb_precision, only: core_rknd
  use calc_roots, only: quadratic_solve

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: nz ! Number of coefficient tuples [-]
  ! Coefficient a (of x^2) in a*x^2 + b*x + c = 0    [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: a_coef
  ! Coefficient b (of x) in a*x^2 + b*x + c = 0      [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: b_coef
  ! Coefficient c in a*x^2 + b*x + c = 0             [-]
  real(kind=core_rknd), dimension(nz), intent(in) :: c_coef

  !--------------------------- Output Variables ---------------------------
  real(kind=core_rknd), dimension(nz, 2), intent(out) :: roots_real ! Real parts of the roots [-]
  ! Imaginary parts of the roots [-]
  real(kind=core_rknd), dimension(nz, 2), intent(out) :: roots_imag

  !--------------------------- Local Variables ---------------------------
  complex(kind=core_rknd), dimension(nz, 2) :: roots ! Complex roots before splitting for F2PY [-]

  !--------------------------- Begin Code ---------------------------

  call quadratic_solve(a_coef, b_coef, c_coef, & ! In
    roots(:,1), roots(:,2)) ! Out
  roots_real = real(roots, kind=core_rknd)
  roots_imag = aimag(roots)

end subroutine f2py_quadratic_solve
