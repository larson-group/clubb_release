!----------------------------------------------------------------------
! $Id$
module jun25

! Description:
!   Contains subroutines for the Clex9 Nov 02 case.

! References:
!   See below.
!----------------------------------------------------------------------

  implicit none

  public :: jun25_altocu_read_t_dependent

  ! Note: bottom point not at surface, so there is no sfc
  ! subroutine

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine jun25_altocu_read_t_dependent( time, &
                                           T_sfc, SH, LH )

! Description:
!   This subroutine reads in the values from the _surface.in file for
!   this case.
!
! References:
!   None
!--------------------------------------------------------------------------

    use stats_precision, only: &
      time_precision ! Variable(s)

    use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                    SH_given, LH_given, &
                                    time_select                   ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    implicit none

    ! Input variables
    real(kind=time_precision), intent(in) :: &
      time             ! Current time          [s]

    ! Output variables
    real, intent(out) :: &
      T_sfc, &   ! surface temperature [K]
      SH,    &   ! sensible heat flux [W/m^2]
      LH         ! latent heat flux [W/m^2]

    ! Local variables
    real :: &
      time_frac ! time fraction used for interpolation

    integer :: &
      before_time, after_time  ! time indexes used for interpolation


    ! ---- Begin Code ----

    ! interpolate T_sfc from time_dependent_input

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                      before_time, after_time, time_frac )

    T_sfc = factor_interp( time_frac, T_sfc_given(after_time), &
                                      T_sfc_given(before_time) )
  
    SH = factor_interp( time_frac, SH_given(after_time), &
                                   SH_given(before_time) )

    LH = factor_interp( time_frac, LH_given(after_time), &
                                   LH_given(before_time) )

    return

  end subroutine jun25_altocu_read_t_dependent 


end module jun25
