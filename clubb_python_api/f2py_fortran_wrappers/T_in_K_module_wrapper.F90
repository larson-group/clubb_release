! T_in_K_module_wrapper.F90 — F2PY wrapper for T_in_K_module.F90
!
! F2PY cannot return arrays from Fortran functions, so this thin wrapper
! converts the function call into a subroutine call with an intent(out) arg.
!
module T_in_K_wrapper

  use clubb_precision, only: core_rknd

  implicit none
  private
  public :: f2py_thlm2T_in_K_1D

contains

  subroutine f2py_thlm2T_in_K_1D(nz, thlm, exner, rcm, T_in_K) &
    bind(C, name="f2py_thlm2t_in_k_1d_")

    use T_in_K_module, only: thlm2T_in_K_api

    integer, intent(in) :: nz
    real(core_rknd), dimension(nz), intent(in) :: thlm, exner, rcm
    real(core_rknd), dimension(nz), intent(out) :: T_in_K

    T_in_K = thlm2T_in_K_api(nz, thlm, exner, rcm)

  end subroutine f2py_thlm2T_in_K_1D

end module T_in_K_wrapper
