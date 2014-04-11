!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module array_index

  ! Description:
  ! Contains indices to variables in larger arrays.
  ! Note that the 'ii' is necessary because 'i' is used in
  ! statistics to track locations in the zt/zm/sfc derived types.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrainm,    & ! Hydrometeor array index for rain water mixing ratio, rr
    iirsnowm,    & ! Hydrometeor array index for snow mixing ratio, rs
    iiricem,     & ! Hydrometeor array index for ice mixing ratio, ri
    iirgraupelm    ! Hydrometeor array index for graupel mixing ratio, rg
!$omp threadprivate(iirrainm, iirsnowm, iiricem, iirgraupelm)

  ! Microphysics concentrations
  integer, public :: &
    iiNrm,       & ! Hydrometeor array index for rain drop concentration, Nr
    iiNsnowm,    & ! Hydrometeor array index for snow concentration, Ns
    iiNim,       & ! Hydrometeor array index for ice concentration, Ni
    iiNgraupelm    ! Hydrometeor array index for graupel concentration, Ng
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

!===============================================================================

end module array_index
