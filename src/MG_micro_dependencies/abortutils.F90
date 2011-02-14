! $Id$
module abortutils
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: endrun
  
  contains
  
!================================================================================================
  subroutine endrun (msg)
    !
    !  Description: The original subroutine aborts the model for abnormal termination.
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------

      character(len=*), intent(in), optional :: msg    ! string to be printed
      
      ! This subroutine should never be called. It is only here to allow MG to compile correctly.
      write(*,*) "WARNING: endrun from abortutils dummy file called. This subroutine should" &
        // " never be called. It is only present to allow MG to compile correctly."
        
      return

  end subroutine endrun
  
end module abortutils
