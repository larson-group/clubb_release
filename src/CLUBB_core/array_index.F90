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

  use clubb_precision, only: &
    core_rknd      ! Precision

  implicit none

  private ! Default Scope

  type sclr_idx_type

    ! Scalar quantities
    integer, public :: & 
      iisclr_rt,    & ! [kg/kg]/[K]/[1e6 mol/mol]
      iisclr_thl,   & ! [kg/kg]/[K]/[1e6 mol/mol]
      iisclr_CO2,   & ! [kg/kg]/[K]/[1e6 mol/mol]
      iiedsclr_rt,  & ! [kg/kg]/[K]/[1e6 mol/mol]
      iiedsclr_thl, & ! [kg/kg]/[K]/[1e6 mol/mol]
      iiedsclr_CO2    ! [kg/kg]/[K]/[1e6 mol/mol]

  end type sclr_idx_type
!$omp declare mapper (sclr_idx_type::x) map ( &
!$omp  x%iiedsclr_co2 &
!$omp )

  public :: sclr_idx_type

!===============================================================================

end module array_index


