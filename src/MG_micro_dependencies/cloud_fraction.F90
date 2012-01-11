! $Id$
module cloud_fraction
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: cldfrc_getparams

  contains

  subroutine cldfrc_getparams (rhminl_out, rhminh_out, premit_out)
  !-------------------------------------------------------------------------------------------------
  ! Return cldfrc tuning parameters.
  ! This is just a dummy subroutine to allow the code to compile. It should never be called.
  !-------------------------------------------------------------------------------------------------

    real(8),          intent(out), optional :: rhminl_out
    real(8),          intent(out), optional :: rhminh_out
    real(8),          intent(out), optional :: premit_out

    stop "ERROR: cldfrc_getparams is a dummy routine and should not be called"

    if ( present(rhminl_out) )      rhminl_out = -9999.99
    if ( present(rhminh_out) )      rhminh_out = -9999.99
    if ( present(premit_out) )      premit_out = -9999.99

  end subroutine cldfrc_getparams

end module cloud_fraction
