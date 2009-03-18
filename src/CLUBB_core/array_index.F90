!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module array_index

! Description:
!   Contains indices to variables in larger arrays.
!   Note that the 'ii' is necessary because 'i' is used in
!   statistics to track locations in the zt/zm/sfc derived types.

! References:
!   None
!-----------------------------------------------------------------------
  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrainm, iirsnowm, iiricem, iirgraupelm ! [kg/kg]

  ! Microphysics number concentration
  integer, public :: &
    iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm ! [#/kg]

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iiCO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedCO2 ! "    "

  private ! Default Scope

end module array_index
!-----------------------------------------------------------------------
