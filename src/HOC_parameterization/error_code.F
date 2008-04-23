!-----------------------------------------------------------------------
! $Id: error_code.F,v 1.3 2008-04-23 17:55:07 faschinj Exp $
!-----------------------------------------------------------------------
        module error_code

!       Description:
!       Since f90/95 lacks enumeration, we're stuck numbering each
!       parameter by hand like this.

!       We are "enumerating" error codes to be used with CLUBB. Adding
!       additional codes is as simple adding an additional integer
!       parameter. The error codes are ranked by severity, the higher
!       number being more servere. When two errors occur, assign the
!       most servere to the output.`
!-----------------------------------------------------------------------

        implicit none

        private

        public ::
     .  CLUBB_no_error, 
     .  CLUBB_var_less_than_zero,
     .  CLUBB_var_equals_NaN, 
     .  CLUBB_singular_matrix,
     .  CLUBB_bad_lapack_arg,
     .  CLUBB_rtm_level_not_found,
     .  CLUBB_var_out_of_bounds,
     .  reportError, 
     .  any_error,
     .  lapack_error,    
     .  clubb_at_debug_level, 
     .  set_clubb_debug_level    

        integer :: CLUBB_debug_level   =  0

        integer, parameter ::
     .  CLUBB_no_error                 =  0, 
     .  CLUBB_var_less_than_zero       =  1,
     .  CLUBB_var_equals_NaN           =  2, 
     .  CLUBB_singular_matrix          =  3,
     .  CLUBB_bad_lapack_arg           =  4,
     .  CLUBB_rtm_level_not_found      =  5,
     .  CLUBB_var_out_of_bounds        =  6 
     
        contains

!-----------------------------------------------
        subroutine reportError( err_code )
!
!       Description: Reports meaning of error code to console.
!
!-----------------------------------------------        

        use constants

        implicit none
        
        integer, intent(in) :: err_code

        select case ( err_code )

        case ( CLUBB_no_error )
                write(fstderr,*) "No errors reported."

        case ( CLUBB_var_less_than_zero )
                write(fstderr,*) "Variable in CLUBB is less than zero."

        case ( CLUBB_singular_matrix )
                write(fstderr,*) "Singular Matrix in CLUBB."

        case ( CLUBB_var_equals_NaN )
                write(fstderr,*) "Variable in CLUBB is NaN."

        case ( CLUBB_bad_lapack_arg )
                write(fstderr,*) 
     .              "Argument used in LAPACK procedure is invalid."
                
        case ( CLUBB_rtm_level_not_found )
                write(fstderr,*) "rtm level not found"
                
        case ( CLUBB_var_out_of_bounds )
                write(fstderr,*) "Input variable is out of bounds."

        case default
                write(fstderr,*) "Unknown error: ", err_code

        end select

        end subroutine reportError

        logical function lapack_error( err_code )
        
        implicit none
        
        integer,intent(in) :: err_code
        
        lapack_error = (err_code == CLUBB_singular_matrix .or.
     .      err_code == CLUBB_bad_lapack_arg ) 
        
        end function

!---------------------------------------------------------------------       
        logical function any_error( err_code )
!---------------------------------------------------------------------        
        implicit none

        integer, intent(in) :: err_code

        any_error = ( err_code == CLUBB_singular_matrix     .or.
     .                err_code == CLUBB_bad_lapack_arg      .or.
     .                err_code == CLUBB_var_equals_NaN      .or.
     .                err_code == CLUBB_rtm_level_not_found .or.
     .                err_code == CLUBB_var_out_of_bounds )


        end function

!------------------------------------------------------------------	
        logical function clubb_at_debug_level( level )
!       
!       Description:
!       Checks to see if CLUBB has been set to a specified debug level
!
!------------------------------------------------------------------
        implicit none

        ! Input variable
        integer, intent(in) :: level   ! The debug level being checked against the current setting

        clubb_at_debug_level = ( level <= CLUBB_debug_level )

        end function clubb_at_debug_level

!----------------------------------------------------------------------
        subroutine set_clubb_debug_level( level )
!
!       Description:
!       Accessor for CLUBB_debug_level
!
!----------------------------------------------------------------------
        implicit none
       
        integer, intent(in) :: level 
        CLUBB_debug_level = level

        end subroutine set_clubb_debug_level

!----------------------------------------------------------------------
        subroutine clubb_debug( level, str )
!
!       Description:
!       Prints a message to file unit fstderr if the level is greater
!       than or equal to the current debug level.        
!----------------------------------------------------------------------

        
        use constants, only:fstderr
        implicit none
        
        character(*), intent(in) :: str
        integer, intent(in)      :: level

        if (level <= CLUBB_debug_level) then
                write(fstderr,*) str
        endif
        end subroutine

        
       
        end module error_code
!-----------------------------------------------------------------------
