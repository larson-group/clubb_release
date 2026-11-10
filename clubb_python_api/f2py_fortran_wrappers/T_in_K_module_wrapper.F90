! T_in_K_module_wrapper.F90 — F2PY wrapper for T_in_K_module.F90
!
! F2PY cannot return arrays from Fortran functions, so this thin wrapper
! converts the function call into a subroutine call with an intent(out) arg.
!
module T_in_K_wrapper

  use clubb_precision, only: core_rknd

  implicit none
  private
  public :: f2py_thlm2T_in_K

contains

  subroutine f2py_thlm2T_in_K( nz, ngrdcol, thlm, exner, rcm, & ! In
                               T_in_K ) bind(C, name="f2py_thlm2t_in_k_") ! Out

    ! Description:
    !   Call thlm2T_in_K with complete column batches through the F2PY bridge.

    use T_in_K_module, only: thlm2T_in_K_api

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: nz ! Number of vertical levels [-]
    integer, intent(in) :: ngrdcol ! Number of grid columns [-]
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: thlm ! Liquid potential temperature  [K]
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: exner ! Exner function                [-]
    real(core_rknd), dimension(ngrdcol, nz), intent(in) :: rcm ! Liquid water mixing ratio     [kg/kg]

    !--------------------------- Output Variables ---------------------------
    real(core_rknd), dimension(ngrdcol, nz), intent(out) :: T_in_K ! Result temperature [K]

    !--------------------------- Begin Code ---------------------------

    T_in_K = thlm2T_in_K_api(nz, ngrdcol, thlm, exner, rcm)

  end subroutine f2py_thlm2T_in_K

end module T_in_K_wrapper
