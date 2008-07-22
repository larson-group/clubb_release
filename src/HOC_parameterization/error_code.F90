!-----------------------------------------------------------------------
! $Id: error_code.F90,v 1.1 2008-07-22 16:04:23 faschinj Exp $
!-----------------------------------------------------------------------

        module error_code

!       Description:
!       Since f90/95 lacks enumeration, we're stuck numbering each
!       parameter by hand like this.

!       We are "enumerating" error codes to be used with CLUBB. Adding
!       additional codes is as simple adding an additional integer
!       parameter. The error codes are ranked by severity, the higher
!       number being more servere. When two errors occur, assign the
!       most servere to the output.

!       This code also handles subroutines related to debug_level. See
!       the 'set_clubb_debug_level' description for more detail.        
!-----------------------------------------------------------------------

        implicit none

        private ! Default Scope

        public :: & 
        clubb_no_error,  & 
        clubb_var_less_than_zero, & 
        clubb_var_equals_NaN,  & 
        clubb_singular_matrix, & 
        clubb_bad_lapack_arg, & 
        clubb_rtm_level_not_found, & 
        clubb_var_out_of_bounds, & 
        reportError,  & 
        fatal_error, & 
        lapack_error,     & 
        clubb_at_debug_level,  & 
        set_clubb_debug_level, & 
        clubb_debug

        private :: clubb_debug_level
        
        ! Model-Wide Debug Level
        integer :: clubb_debug_level   =  0

        ! Error Code Values
        integer, parameter :: & 
        clubb_no_error                 =  0,  & 
        clubb_var_less_than_zero       =  1, & 
        clubb_var_equals_NaN           =  2,  & 
        clubb_singular_matrix          =  3, & 
        clubb_bad_lapack_arg           =  4, & 
        clubb_rtm_level_not_found      =  5, & 
        clubb_var_out_of_bounds        =  6 
     
        contains

!-----------------------------------------------
        subroutine reportError( err_code )
!
!       Description: Reports meaning of error code to console.
!
!-----------------------------------------------        

        use constants, only: & 
            fstderr ! Variable(s)

        implicit none
        
        ! Input Variable
        integer, intent(in) :: err_code ! Error Code being examined

        select case ( err_code )

        case ( clubb_no_error )
                write(fstderr,*) "No errors reported."

        case ( clubb_var_less_than_zero )
                write(fstderr,*) "Variable in CLUBB is less than zero."

        case ( clubb_singular_matrix )
                write(fstderr,*) "Singular Matrix in CLUBB."

        case ( clubb_var_equals_NaN )
                write(fstderr,*) "Variable in CLUBB is NaN."

        case ( clubb_bad_lapack_arg )
                write(fstderr,*)  & 
                    "Argument used in LAPACK procedure is invalid."
                
        case ( clubb_rtm_level_not_found )
                write(fstderr,*) "rtm level not found"
                
        case ( clubb_var_out_of_bounds )
                write(fstderr,*) "Input variable is out of bounds."

        case default
                write(fstderr,*) "Unknown error: ", err_code

        end select

        end subroutine reportError
!---------------------------------------------------------------------
        logical function lapack_error( err_code )
!
!       Description: Checks to see if the err_code is equal to one
!       caused by an error encountered using lapack
!---------------------------------------------------------------------        
        implicit none
        
        ! Input variable
        integer,intent(in) :: err_code ! Error Code being examined
        
        lapack_error = (err_code == clubb_singular_matrix .or. & 
            err_code == clubb_bad_lapack_arg ) 
        
        end function lapack_error

!---------------------------------------------------------------------       
        logical function fatal_error( err_code )
!
!       Description: Checks to see if the err_code is one that usually
!       causes an exit in other parts of CLUBB.
!---------------------------------------------------------------------        
        implicit none

        ! Input Variable
        integer, intent(in) :: err_code ! Error Code being examined

        fatal_error = ( err_code == clubb_singular_matrix     .or. & 
                      err_code == clubb_bad_lapack_arg      .or. & 
                      err_code == clubb_var_equals_NaN      .or. & 
                      err_code == clubb_rtm_level_not_found .or. & 
                      err_code == clubb_var_out_of_bounds )


        end function fatal_error

!------------------------------------------------------------------	
        logical function clubb_at_debug_level( level )
!       
!       Description:
!       Checks to see if clubb has been set to a specified debug level
!------------------------------------------------------------------
        implicit none

        ! Input variable
        integer, intent(in) :: level   ! The debug level being checked against the current setting

        clubb_at_debug_level = ( level <= clubb_debug_level )

        end function clubb_at_debug_level

!----------------------------------------------------------------------
        subroutine set_clubb_debug_level( level )
!
!       Description:
!       Accessor for clubb_debug_level
!
!        0 => Print no debug messages to the screen
!        1 => Print lightweight debug messages, e.g. print statements
!        2 => Print debug messages that require extra testing,
!                e.g. checks for NaNs and spurious negative values.

!----------------------------------------------------------------------
        implicit none
       
        ! Input variable
        integer, intent(in) :: level ! The debug level being checked against the current setting 
        clubb_debug_level = level

        end subroutine set_clubb_debug_level

!----------------------------------------------------------------------
        subroutine clubb_debug( level, str )
!
!       Description:
!       Prints a message to file unit fstderr if the level is greater
!       than or equal to the current debug level.        
!----------------------------------------------------------------------
        use constants, only: & 
            fstderr ! Variable(s)

        implicit none

        ! Input Variables        
        character*(*), intent(in) :: str ! The message being reported

        integer, intent(in)      :: level ! The debug level being checked against the current setting

        if (level <= clubb_debug_level) then
                write(*,*) str
        endif
        
        end subroutine clubb_debug

        end module error_code
!-----------------------------------------------------------------------
