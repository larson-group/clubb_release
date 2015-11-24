! $Id$
module code_timer_module

! Description:
!   This module contains a diagnostic timer utility that can be used
!   to time a piece of code.

  implicit none

  private ! Set default scope

  ! A large integer kind is used for timing to avoid
  ! overflow.
  integer, parameter, public :: &
    timer_kind = selected_int_kind(16)

  ! A timer!!
  type timer_t
    integer(kind=timer_kind) :: time_elapsed        ! Time elapsed [msec]
    integer :: secstart, msstart                    ! Timer starting times
  end type timer_t

  public :: timer_t, timer_start, timer_stop

  contains

  !-----------------------------------------------------------------------
  subroutine timer_start( timer )

  ! Description:
  !   Starts the timer

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variables
    type(timer_t), intent(inout) :: timer

    ! Local Variables
    integer, dimension(8) :: values
  !-----------------------------------------------------------------------
    !----- Begin Code -----
    call date_and_time(VALUES=values)

    timer%secstart = values(7)
    timer%msstart  = values(8)

    return
  end subroutine timer_start
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine timer_stop( timer )

  ! Description:
  !   Stops the timer

  ! References:
  !   None
  !-----------------------------------------------------------------------
    implicit none

    ! Input/Output Variables
    type(timer_t), intent(inout) :: timer

    ! Local Variables
    integer :: values(8), secend, msend, secdif, msdif

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    call date_and_time(VALUES=values)

    secend = values(7)
    msend  = values(8)

    secdif = secend - timer%secstart
    if (secdif < 0) secdif = secdif + 60

    msdif = msend - timer%msstart
    timer%time_elapsed = timer%time_elapsed + int( secdif*1000 + msdif, kind=timer_kind )

    return
  end subroutine timer_stop
  !-----------------------------------------------------------------------

end module code_timer_module
