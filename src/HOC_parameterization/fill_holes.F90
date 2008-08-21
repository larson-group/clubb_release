!-----------------------------------------------------------------------
! $Id$

module fill_holes

implicit none

public :: fill_holes_driver

private :: fill_holes_multiplicative, vertical_avg

private ! Set Default Scope

contains

!--------------------------------------------------------------------------
subroutine fill_holes_driver( num_pts, threshold, field_grid, &
                              field )

!       Description:
!       This subroutine clips values of 'field' that are below 'threshold' as 
!       much as possible (i.e. "fills holes"), but conserves the total 
!       integrated mass of 'field'.  This prevents clipping from acting as a 
!       spurious source.
!
!       Mass is conserved by reducing the clipped field everywhere by a constant
!       multiplicative coefficient.
!
!       This subroutine does not guarantee that the clipped field will exceed 
!       threshold everywhere; blunt clipping is needed for that.
!
!       This subroutine doesn't account for variation in air density with 
!       altitude, and therefore won't properly conserve the vertical integral of
!       field if we upgrade to theanelastic equation!!!

!       References:
!       ``Numerical Methods for Wave Equations in Geophysical Fluid
!       Dynamics", Durran (1999), p. 292.
!-----------------------------------------------------------------------

use grid_class, only: & 
   gr ! Variable

implicit none

! Input variables
integer, intent(in) :: & 
  num_pts  ! The number of points on either side of the hole;
         ! Mass is drawn from these points to fill the hole.  []  

real, intent(in) :: & 
  threshold  ! A threshold (e.g. wtol*wtol) below which field must not
           ! fall                           [Units vary; same as field]

character(len=2), intent(in) :: & 
  field_grid ! The grid of the field, either zt or zm

! Input/Output variable
real, dimension(gr%nnzp), intent(inout) :: & 
  field  ! The field (e.g. wp2) that contains holes
       !                                    [Units same as threshold]

! Local Variables
integer :: & 
  k,          & ! Loop index for absolute grid level              []
  begin_idx,  & ! Lower grid level of hole-filling range          []
  end_idx    ! Upper grid level of hole-filling rang           []

!-----------------------------------------------------------------------

! Check whether any holes exist in the entire profile.
! The lowest level (k=1) should not be included, as the hole-filling 
! scheme should not alter the set value of 'field' at the surface.
if ( any( field( 2:gr%nnzp ) < threshold ) ) then

   ! Make one pass up the profile, filling holes as much as we can 
   ! using nearby mass.
   ! The lowest level (k=1) should not be included in the loop, as
   ! the hole-filling scheme should not alter the set value of 'field' 
   ! at the surface.
   do k = 2+num_pts, gr%nnzp-num_pts, 1
   
      begin_idx = k - num_pts
      end_idx   = k + num_pts 

      if ( any( field( begin_idx:end_idx ) < threshold ) ) then

         call fill_holes_multiplicative( begin_idx, end_idx,  & 
                                threshold, field_grid,  & 
                                field( begin_idx:end_idx ) )

      endif
    
   enddo 

   ! Fill holes globally, to maximize the chance that all holes are filled
   ! The lowest level (k=1) should not be included, as the hole-filling 
   ! scheme should not alter the set value of 'field' at the surface.
   if ( any( field( 2:gr%nnzp ) < threshold ) ) then

      call fill_holes_multiplicative( 2, gr%nnzp,  & 
                                threshold, field_grid,  & 
                                field( 2:gr%nnzp ) )

   endif  

endif  ! End overall check for existence of holes      
    
return

end subroutine fill_holes_driver

!-----------------------------------------------------------------------
subroutine fill_holes_multiplicative & 
          ( begin_idx, end_idx, threshold, field_grid, field )

!       Description: 
!       This subroutine clips values of 'field' that are below 'threshold' as 
!       much as possible (i.e. "fills holes"), but conserves the total 
!       integrated mass of 'field'.  This prevents clipping from acting as a 
!       spurious source.
!
!       Mass is conserved by reducing the clipped field everywhere by a constant
!       multiplicative coefficient.
!
!       This subroutine does not guarantee that the clipped field will exceed 
!       threshold everywhere; blunt clipping is needed for that.
!
!       This subroutine doesn't account for variation in air density with 
!       altitude, and therefore won't properly conserve the vertical integral 
!       of field if we upgrade to the anelastic equation!!!

!       References:
!       ``Numerical Methods for Wave Equations in Geophysical Fluid
!       Dynamics", Durran (1999), p. 292.  
!-----------------------------------------------------------------------

implicit none

! Input variables
integer, intent(in) :: & 
  begin_idx,  & ! The beginning index (e.g. k=2) of the range of hole-filling 
  end_idx    ! The end index (e.g. k=gr%nnzp) of the range of hole-filling

real, intent(in) :: & 
  threshold  ! A threshold (e.g. wtol*wtol) below which field must not 
           ! fall                           [Units vary; same as field]

character(len=2), intent(in) :: & 
  field_grid ! The grid of the field, either zt or zm

! Input/Output variable
real, dimension(end_idx-begin_idx+1), intent(inout) ::  & 
  field  ! The field (e.g. wp2) that contains holes 
       !                                    [Units same as threshold]

! Local Variables
real, dimension(end_idx-begin_idx+1)  ::  & 
  field_clipped  ! The raw field (e.g. wp2) that contains no holes 
               !                          [Units same as threshold]

real ::  & 
  field_avg,  & ! Vertical average of field [Units of field]
  field_clipped_avg,   & ! Vertical average of clipped field [Units of field]
  mass_fraction  ! Coefficient that multiplies clipped field 
               ! in order to conserve mass.                      []

!-----------------------------------------------------------------------

! Compute the field's vertical average, which we must conserve.
field_avg = vertical_avg( begin_idx, end_idx,  & 
                          field_grid, field )

! Clip small or negative values from field.
if ( field_avg >= threshold ) then  
   ! We know we can fill in holes completely
   field_clipped = max( threshold, field )
else 
   ! We can only fill in holes partly;
   ! to do so, we remove all mass above threshold.
   field_clipped = min( threshold, field )
end if

! Compute the clipped field's vertical integral.
! clipped_total_mass >= original_total_mass
field_clipped_avg = vertical_avg( begin_idx, end_idx,  & 
                                  field_grid, field_clipped )

! If the difference between the field_clipped_avg and the 
! threshold is so small that it falls within numerical 
! round-off, return to the parent subroutine without altering
! the field in order to avoid divide-by-zero error.
!        if ( abs(field_clipped_avg - threshold)  &
!              < threshold*epsilon(threshold) ) then
if ( abs(field_clipped_avg - threshold) == 0.0 ) then
   return
endif

! Compute coefficient that makes the clipped field have the same
! mass as the original field.
! We should always have mass_fraction > 0.
mass_fraction = ( field_avg - threshold ) / & 
                      ( field_clipped_avg - threshold )

! Output normalized, filled field
field = mass_fraction * ( field_clipped - threshold )  & 
             + threshold


return
end subroutine fill_holes_multiplicative

!-----------------------------------------------------------------------
function vertical_avg( begin_idx, end_idx, field_grid, field )

!       Description:
!       Compute the vertical average of a field.  
!       Does not account for air density that varies with altitude!!!

!       References:
!-----------------------------------------------------------------------

use grid_class, only: & 
   gr ! Variable

use error_code, only: & 
   clubb_debug  ! Subroutine

implicit none

! Input variables
integer, intent(in) :: & 
  begin_idx,  & ! The beginning index (e.g. 2) of the range of integration of field 
  end_idx    ! The end index (e.g. gr%nnzp) of the range of integration of field


character(len=2), intent(in) :: & 
  field_grid ! The grid of the field, either zt or zm

real, dimension(end_idx-begin_idx+1), intent(in) ::  & 
  field  ! The field (e.g. wp2) that we want to average  [Units vary]
       ! The field points need to be arranged from lowest to highest in altitude

! Output variable
real :: & 
  vertical_avg  ! Vertical average of field    [Units of field]

! Local variables
real :: & 
  vertical_integral,  & ! Vertical integral of height [(length) * (units of field)]
  height             ! Altitude range over which we integrate  [m]

integer ::  & 
  k,      & ! Loop index for absolute grid level
  k_rel  ! Loop index for grid level relative to starting idx of field

!-------------------------------------------------------------------------------

!  Assertion checks: that begin_idx <= gr%nnzp - 1 
!                    that end_idx   >= 2
!                    that begin_idx <= end_idx

! Initialize variables
vertical_integral = 0.0
height = 0.0

! If field is on the zt grid
if ( field_grid == "zt" ) then

   k_rel = 1
   do k = max(2,begin_idx), end_idx, 1

      vertical_integral = vertical_integral  & 
                               +  field(k_rel) / gr%dzt(k)
      height = height + 1.0 / gr%dzt(k)
      k_rel = k_rel + 1

   enddo

elseif ( field_grid == "zm" )  then

   k_rel = 1
   do k = max(2,begin_idx), min(gr%nnzp-1,end_idx), 1

      vertical_integral = vertical_integral  & 
                                 +  field(k_rel) / gr%dzm(k)
      height = height + 1.0 / gr%dzm(k)
      k_rel = k_rel + 1

   enddo

else

  call clubb_debug( 0,  & 
     "Neither zt nor zm grid is specified in vert_integrate" ) 
  vertical_integral = -9.0e20
  height = 1.0

endif

vertical_avg = vertical_integral / height


return
end function vertical_avg

!-----------------------------------------------------------------------

end module fill_holes
