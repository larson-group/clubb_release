!$Id: Skw_func.F90 2871 2008-09-06 13:36:17Z griffinb $
!-------------------------------------------------------------------------------
module Skw_module

  implicit none

  private ! Default Scope

  public :: Skw_func

  contains

!-------------------------------------------------------------------------------
  elemental function Skw_func( wp2, wp3 )  &
    result( Skw )

! Description:
!   Calculate the skewness of w, Skw.

! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only:  &
      wtol_sqd,  &! Constant for w_{tol}^2, i.e. threshold for vertical velocity
      Skw_max_mag ! Max magnitude of skewness

    implicit none

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    logical, parameter ::  & 
      l_clipping_kluge = .false.

    ! Input Variables
    real, intent(in) :: & 
      wp2,  & ! w'^2    [m^2/s^2]
      wp3     ! w'^3    [m^3/s^3]

    ! Output Variable
    real :: & 
      Skw     ! Result Skw [-]

    ! ---- Begin Code ----

    Skw = wp3 / ( max( wp2, wtol_sqd ) )**1.5

    ! This is no longer need since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code
    if ( l_clipping_kluge ) then
      Skw = min( max( Skw, -Skw_max_mag ), Skw_max_mag )
    end if

    return
  end function Skw_func
!-----------------------------------------------------------------------

end module Skw_module
