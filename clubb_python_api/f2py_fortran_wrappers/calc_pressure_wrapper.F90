! calc_pressure_wrapper.F90 — F2PY wrapper for calc_pressure.F90
!
! Named after the source file it wraps. Imports stored_grid from
! derived_type_storage and passes it where the original routine expects
! type(grid).
!
module calc_pressure_wrapper

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid

  implicit none
  private
  public :: f2py_init_pressure

contains

  subroutine f2py_init_pressure(ngrdcol, nzt, nzm, thvm, p_sfc, &
                                p_in_Pa, exner, p_in_Pa_zm, exner_zm) &
    bind(C, name="f2py_init_pressure_")

    use calc_pressure, only: init_pressure

    integer, intent(in) :: ngrdcol, nzt, nzm
    real(core_rknd), dimension(ngrdcol, nzt), intent(in)  :: thvm
    real(core_rknd), dimension(ngrdcol),      intent(in)  :: p_sfc
    real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: p_in_Pa, exner
    real(core_rknd), dimension(ngrdcol, nzm), intent(out) :: p_in_Pa_zm, exner_zm

    ! The ONLY thing this wrapper does is pass stored_grid where the
    ! original routine expects type(grid). Everything else passes through.
    call init_pressure(ngrdcol, stored_grid, thvm, p_sfc, &
                       p_in_Pa, exner, p_in_Pa_zm, exner_zm)

  end subroutine f2py_init_pressure

end module calc_pressure_wrapper
subroutine f2py_calculate_thvm(nzt, ngrdcol, &
    thlm, rtm, rcm, exner, thv_ds_zt, thvm)

  use clubb_precision, only: core_rknd
  use calc_pressure, only: calculate_thvm

  implicit none

  integer, intent(in) :: nzt, ngrdcol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    thlm, rtm, rcm, exner, thv_ds_zt
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: thvm

  call calculate_thvm(nzt, ngrdcol, thlm, rtm, rcm, exner, thv_ds_zt, thvm)

end subroutine f2py_calculate_thvm

