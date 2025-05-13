!-------------------------------------------------------------------------------
! $Id$
!-------------------------------------------------------------------------------

module error_code

! Description:
!   Since f90/95 lacks enumeration, we're stuck numbering each
!   error code by hand like this.

!   We are "enumerating" error codes to be used with CLUBB. Adding
!   additional codes is as simple adding an additional integer
!   parameter. The error codes are ranked by severity, the higher
!   number being more servere. When two errors occur, assign the
!   most servere to the output.

!   This code also handles subroutines related to debug_level. See
!   the 'set_clubb_debug_level_api' description for more detail.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    private ! Default Scope

    public :: &
        clubb_at_least_debug_level_api, &
        set_clubb_debug_level_api, &
        initialize_error_headers

    private :: clubb_debug_level

    ! Model-Wide Debug Level
    integer, save :: clubb_debug_level = 0
    !$acc declare create(clubb_debug_level)

    ! Error Code Values
    integer, parameter, public :: & 
        clubb_no_error                 = 0,  &
        clubb_generalized_grd_test_err = 50, &
        clubb_fatal_error              = 99

    character(len=35), public :: err_header

    !$omp threadprivate(err_header)

    contains
!-------------------------------------------------------------------------------
! Description:
!   Checks to see if clubb has been set to a specified debug level
!-------------------------------------------------------------------------------
    logical function clubb_at_least_debug_level_api( level )

        implicit none

        !$acc routine seq

        ! Input variable
        integer, intent(in) :: level   ! The debug level being checked against the current setting

        ! ---- Begin Code ----

        clubb_at_least_debug_level_api = ( level <= clubb_debug_level )

        return

    end function clubb_at_least_debug_level_api


    subroutine initialize_error_headers

        implicit none

#ifdef _OPENMP
        integer :: omp_get_thread_num
        write(err_header,'(A7,I7,A20)') "Thread ", omp_get_thread_num(), " -- CLUBB -- ERROR: "
#else
#ifndef CLUBB_CAM
        ! This code cannot be used for CAM because
        ! it causes issues when tested with the
        ! NAG compiler.
        integer :: getpid
        write(err_header,'(A7,I7,A20)') "Process ", getpid(), " -- CLUBB -- ERROR: "
#else
        write(err_header,'(A20)') " -- CLUBB -- ERROR: "
#endif /* CLUBB_CAM */
#endif

    end subroutine initialize_error_headers


!-------------------------------------------------------------------------------
!  Description:
!    Accessor for clubb_debug_level
!
!   0 => Print no debug messages to the screen
!   1 => Print lightweight debug messages, e.g. print statements
!   2 => Print debug messages that require extra testing,
!        e.g. checks for NaNs and spurious negative values.
!  References:
!    None
!-------------------------------------------------------------------------------
    subroutine set_clubb_debug_level_api( level )

        implicit none

        ! Input variable
        integer, intent(in) :: level ! The debug level being checked against the current setting

        ! ---- Begin Code ----

        clubb_debug_level = max(level,0)
        !$acc update device(clubb_debug_level)


        return
        end subroutine set_clubb_debug_level_api

    end module error_code
