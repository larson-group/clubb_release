! $Id$
subroutine gffgch(t       ,es      ,itype   )
    !
    !  Description: The original subroutine computes saturation vapor pressure over
    !               water and/or ice using Goff & Gratch (1946) relationships.
    !               In our case, we don't directly use this, so this subroutine does
    !               nothing.
    !
    !---------------------------------------------------------------------------------
    implicit none
    
    real, intent(in) :: t ! Temperature
    integer, intent(inout) :: itype ! Flag for ice phase and associated transistion
    real, intent(out) :: es ! Saturation vapor pressure
    
    ! Avoid compiler warnings
    if(.false.) then
      es = 0
    end if
    
    return
    
end subroutine gffgch
