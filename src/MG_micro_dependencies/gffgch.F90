! $Id$
subroutine gffgch(t       ,es      ,itype   )
    !
    !  Description: The original subroutine computes saturation vapor pressure over
    !               water and/or ice using Goff & Gratch (1946) relationships.
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------
    implicit none
    
    real, intent(in) :: t ! Temperature
    integer, intent(inout) :: itype ! Flag for ice phase and associated transistion
    real, intent(out) :: es ! Saturation vapor pressure
    
    ! Avoid compiler warnings
    es = 0
    
    ! This subroutine should never be called. It is only here to allow MG to compile correctly.
    write(*,*) "WARNING: gffgch dummy file called. This subroutine should" &
      // " never be called. It is only present to allow MG to compile correctly."
    
    return
    
end subroutine gffgch
