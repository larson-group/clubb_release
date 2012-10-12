! $Id$
module abortutils
!
! Dummy module for importing variables into morrison_gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: endrun
  
  contains
  
!================================================================================================
  subroutine endrun (msg)
    !
    !  Description: The original subroutine aborts the model for abnormal termination.
    !               This dummy subroutine simply prints out the error message.
    !
    !---------------------------------------------------------------------------------

      character(len=*), intent(in), optional :: msg    ! string to be printed
      
      ! This subroutine should never be called. It is only here to allow MG to compile correctly.
      write(*,*) msg
        
      return

  end subroutine endrun
  
end module abortutils
