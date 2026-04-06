subroutine f2py_cholesky_factor(ndim, a_input, a_scaling, a_cholesky, l_scaled)

  use clubb_precision, only: core_rknd
  use matrix_operations, only: cholesky_factor

  implicit none

  integer, intent(in) :: ndim
  real(core_rknd), dimension(ndim, ndim), intent(in) :: a_input
  real(core_rknd), dimension(ndim), intent(out) :: a_scaling
  real(core_rknd), dimension(ndim, ndim), intent(out) :: a_cholesky
  logical, intent(out) :: l_scaled

  call cholesky_factor(ndim, a_input, a_scaling, a_cholesky, l_scaled)

end subroutine f2py_cholesky_factor

subroutine f2py_mirror_lower_triangular_matrix(nvars, matrix)

  use clubb_precision, only: core_rknd
  use matrix_operations, only: mirror_lower_triangular_matrix

  implicit none

  integer, intent(in) :: nvars
  real(core_rknd), dimension(nvars, nvars), intent(inout) :: matrix

  call mirror_lower_triangular_matrix(nvars, matrix)

end subroutine f2py_mirror_lower_triangular_matrix
