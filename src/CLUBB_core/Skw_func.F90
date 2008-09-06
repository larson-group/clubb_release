!$Id$
module Skw

implicit none

private ! Default Scope

public :: Skw_func

contains

!-----------------------------------------------------------------------
function Skw_func( wp2, wp3 )  &
result( Skw )

! Description:
! Calculate the skewness of w, Skw.

! References:
!-----------------------------------------------------------------------

use constants, only:  &
    wtol_sqd ! Variable(s)

implicit none

! External
intrinsic :: min, max

! Parameter Constants
logical, parameter ::  & 
  l_clipping_kluge = .false.

! Input Variables
real, intent(in) :: & 
  wp2,  & ! w'^2    [m^2/s^2]
  wp3     ! w'^3    [m^3/s^3]

! Output Variable
real :: & 
  Skw     ! Result Skw [-]

Skw = wp3 / ( max( wp2, wtol_sqd ) )**1.5

if ( l_clipping_kluge ) then
  Skw = min( max( Skw, -4.5 ), 4.5)
endif

return
end function Skw_func
!-----------------------------------------------------------------------

end module Skw
