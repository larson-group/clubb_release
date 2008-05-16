         module Skw
         
         implicit none
         
         private

         public :: Skw_func
         
         contains

!-----------------------------------------------------------------------
         function Skw_func( wp2, wp3, wtol ) 
     .   result( Skw )

!        Description:
!        Calculate Skw

!        References:
!-----------------------------------------------------------------------
         implicit none

         ! External
         intrinsic :: min, max

         ! Parameter Constants
         logical, parameter :: 
     .   clipping_kluge = .false.

         ! Input 
         real, intent(in) ::
     .   wp2, ! w'^2    [m^2/s^2]
     .   wp3  ! w'^3    [m^3/s^3]
         real, intent(in) ::
     .   wtol ! w tol.  [m/s]

         real ::
     .   Skw ! Result Skw [-]

         Skw = wp3 / ( max( wp2, wtol**2 ) )**1.5

         if ( clipping_kluge ) then
           Skw = min( max( Skw, -4.5 ), 4.5)
         end if

         return
         end function Skw_func
!-----------------------------------------------------------------------

         end module Skw
