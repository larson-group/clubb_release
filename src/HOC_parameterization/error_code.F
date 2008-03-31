!-----------------------------------------------------------------------
! $Id: error_code.F,v 1.1 2008-03-31 17:40:16 faschinj Exp $
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
     .  reportError   
        
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
        subroutine reportError( code )
!
!       Description: Reports meaning of error code to console.
!
!-----------------------------------------------        

        use constants

        implicit none
        
        integer, intent(in) :: code

        select case ( code )

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
                write(fstderr,*) "Unknown error: ", code

        end select

        end subroutine reportError
        
        end module error_code
!-----------------------------------------------------------------------
