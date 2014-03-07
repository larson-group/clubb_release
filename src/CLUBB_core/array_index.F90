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
    iirrainm,    & ! Rain water mixing ratio  [kg/kg]
    iirsnowm,    & ! Snow mixing ratio        [kg/kg]
    iiricem,     & ! Ice mixing ratio         [kg/kg]
    iirgraupelm    ! Graupel mixing ratio     [kg/kg]
!$omp threadprivate(iirrainm, iirsnowm, iiricem, iirgraupelm)

  ! Microphysics concentrations
  integer, public :: &
    iiNrm,       & ! Rain drop concentration                       [num/kg]
    iiNsnowm,    & ! Snow concentration                            [num/kg]
    iiNim,       & ! Ice concentration                             [num/kg]
    iiNgraupelm    ! Graupel concentration                         [num/kg]
!$omp threadprivate(iiNrm, iiNsnowm, iiNim, iiNgraupelm)

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
!$omp threadprivate(iisclr_rt, iisclr_thl, iisclr_CO2, &
!$omp   iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2)

  ! Logical fields
  logical, dimension(:), allocatable, public :: &
    l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
    l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio
!$omp threadprivate(l_frozen_hm, l_mix_rat_hm)

  private ! Default Scope

end module array_index
!-----------------------------------------------------------------------
