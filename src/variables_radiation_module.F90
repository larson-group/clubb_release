!---------------------------------------------------------------
! $Id$
!===============================================================================
module variables_radiation_module

!   This module contains definitions of all radiation arrays
!   used in the single column model, as well as subroutines to
!   allocate, deallocate, and initialize them.
!---------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd, & ! Variable(s)
    dp

  implicit none

  public :: &
    setup_radiation_variables, &
    cleanup_radiation_variables

  private ! Set Default Scoping

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    radht_LW, & ! LW heating rate   [K/s]
    radht_SW, & ! SW heating rate   [K/s]
    Frad_SW,  & ! SW radiative flux [W/m^2]
    Frad_LW     ! LW radiative flux [W/m^2]

!$omp threadprivate(radht_LW, radht_SW, Frad_SW, Frad_LW)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    T_in_K,   & ! Temperature        [K]
    rcil,     & ! Ice mixing ratio   [kg/kg]
    o3l         ! Ozone mixing ratio [kg/kg]

!$omp threadprivate(T_in_K, rcil, o3l)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    rsm_2d,& ! Two-dimensional copies of the input parameters
    rcm_in_cloud_2d, &
    cloud_frac_2d, &
    ice_supersat_frac_2d

!$omp threadprivate(rsm_2d, rcm_in_cloud_2d, cloud_frac_2d, &
!$omp   ice_supersat_frac_2d)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    radht_SW_2d, & ! SW Radiative heating rate  [W/m^2]
    radht_LW_2d    ! LW Radiative heating rate  [W/m^2]

!$omp threadprivate(radht_SW_2d, radht_LW_2d)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    p_in_mb, &   ! Pressure in millibars        [mb]
    sp_humidity  ! Specific humidity            [kg/kg]

!$omp threadprivate(p_in_mb, sp_humidity)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    Frad_uLW, & ! LW upwelling flux         [W/m^2]
    Frad_dLW, & ! LW downwelling flux       [W/m^2]
    Frad_uSW, & ! SW upwelling flux         [W/m^2]
    Frad_dSW    ! SW downwelling flux       [W/m^2]

!$omp threadprivate(Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW)

  real(kind = dp), public, dimension(:,:), allocatable :: &
     fdswcl, &  !Downward clear-sky SW flux                 (W/m^-2).
     fuswcl, &  !Upward clear-sky SW flux                   (W/m^-2).
     fdlwcl, &  !Downward clear-sky LW flux                 (W/m^-2).
     fulwcl     !Upward clear-sky LW flux                   (W/m^-2).

!$omp threadprivate(fdswcl, fuswcl, fdlwcl, fulwcl)

  ! Constant parameters
  integer, public, parameter :: &
    nlen = 1, &   ! Length of the total domain
    slen = 1      ! Length of the sub domain

  contains

  !---------------------------------------------------------------------
  subroutine setup_radiation_variables( nzmax, lin_int_buffer, &
                                        extended_atmos_range_size )
  ! Description:
  !   Allocates and initializes prognostic scalar and array variables
  !   for the CLUBB model code.
  !---------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzmax, & ! Number of grid levels [-]
      lin_int_buffer,& ! Number of interpolated levels between the computational
                       ! grid and the extended atmosphere [-]
      extended_atmos_range_size ! The number of levels in the extended atmosphere [-]

    ! Local Variables

    integer :: rad_zt_dim, rad_zm_dim ! Dimensions of the radiation grid

    !----------------------------BEGIN CODE-------------------------------

    rad_zt_dim = (nzmax-1)+lin_int_buffer+extended_atmos_range_size
    rad_zm_dim = (nzmax-1)+lin_int_buffer+extended_atmos_range_size+1


    ! --- Allocation ---

    allocate( radht_SW(1:nzmax) )
    allocate( radht_LW(1:nzmax) )
    allocate( Frad_SW(1:nzmax) )
    allocate( Frad_LW(1:nzmax) )

    allocate( T_in_K(nlen, rad_zt_dim ) )
    allocate( rcil(nlen, rad_zt_dim ) )
    allocate( o3l(nlen, rad_zt_dim ) )

    allocate( rsm_2d(nlen, rad_zt_dim ) )
    allocate( rcm_in_cloud_2d(nlen, rad_zt_dim ) )
    allocate( cloud_frac_2d(nlen, rad_zt_dim ) )
    allocate( ice_supersat_frac_2d(nlen, rad_zt_dim ) )

    allocate( p_in_mb(nlen, rad_zt_dim) )
    allocate( sp_humidity(nlen, rad_zt_dim) )

    allocate( radht_SW_2d(nlen, rad_zt_dim ) )
    allocate( radht_LW_2d(nlen, rad_zt_dim ) )

    allocate( Frad_uLW(nlen, rad_zm_dim ) )
    allocate( Frad_dLW(nlen, rad_zm_dim ) )
    allocate( Frad_uSW(nlen, rad_zm_dim ) )
    allocate( Frad_dSW(nlen, rad_zm_dim ) )

    allocate( fdswcl(slen, rad_zm_dim ) )
    allocate( fuswcl(slen, rad_zm_dim ) )
    allocate( fdlwcl(slen, rad_zm_dim ) )
    allocate( fulwcl(slen, rad_zm_dim ) )

    ! --- Initialization ---

    ! CLUBB zt
    radht_SW = 0.0_core_rknd
    radht_LW = 0.0_core_rknd
    Frad_SW = 0.0_core_rknd
    Frad_LW = 0.0_core_rknd

    ! CLUBB zt + extended levels
    T_in_K = 0.0_dp
    rcil = 0.0_dp
    o3l = 0.0_dp
    rsm_2d = 0.0_dp
    rcm_in_cloud_2d = 0.0_dp
    cloud_frac_2d = 0.0_dp
    ice_supersat_frac_2d = 0.0_dp
    radht_SW_2d = 0.0_dp
    radht_LW_2d = 0.0_dp
    p_in_mb = 0.0_dp
    sp_humidity = 0.0_dp

    ! CLUBB zm + extended levels
    Frad_uLW = 0.0_dp
    Frad_dLW = 0.0_dp
    Frad_uSW = 0.0_dp
    Frad_dSW = 0.0_dp
    fdswcl = 0.0_dp
    fuswcl = 0.0_dp
    fdlwcl = 0.0_dp
    fulwcl = 0.0_dp

    return
  end subroutine setup_radiation_variables

  !---------------------------------------------------------------------
  subroutine cleanup_radiation_variables( )
  
  ! Description:
  !   Subroutine to deallocate variables defined in module global
  !---------------------------------------------------------------------
 
    implicit none

    ! --- Deallocate ---

    if (allocated(radht_SW)) deallocate( radht_SW )
    if (allocated(radht_LW)) deallocate( radht_LW )
    if (allocated(Frad_SW)) deallocate( Frad_SW )
    if (allocated(Frad_LW)) deallocate( Frad_LW )

    if (allocated(T_in_K)) deallocate( T_in_K )
    if (allocated(rcil)) deallocate( rcil )
    if (allocated(o3l)) deallocate( o3l )

    if (allocated(rsm_2d)) deallocate( rsm_2d )
    if (allocated(rcm_in_cloud_2d)) deallocate( rcm_in_cloud_2d )
    if (allocated(cloud_frac_2d)) deallocate( cloud_frac_2d )
    if (allocated(ice_supersat_frac_2d)) deallocate( ice_supersat_frac_2d )

    if (allocated(p_in_mb)) deallocate( p_in_mb )
    if (allocated(sp_humidity)) deallocate( sp_humidity )

    if (allocated(radht_SW_2d)) deallocate( radht_SW_2d )
    if (allocated(radht_LW_2d)) deallocate( radht_LW_2d )

    if (allocated(Frad_uLW)) deallocate( Frad_uLW )
    if (allocated(Frad_dLW)) deallocate( Frad_dLW )
    if (allocated(Frad_uSW)) deallocate( Frad_uSW )
    if (allocated(Frad_dSW)) deallocate( Frad_dSW )

    if (allocated(fdswcl)) deallocate( fdswcl )
    if (allocated(fuswcl)) deallocate( fuswcl )
    if (allocated(fdlwcl)) deallocate( fdlwcl )
    if (allocated(fulwcl)) deallocate( fulwcl )

    return
  end subroutine cleanup_radiation_variables
    
end module variables_radiation_module
