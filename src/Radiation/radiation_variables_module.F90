module radiation_variables_module

  use clubb_precision, only: &
    core_rknd, &
    dp

  use extended_atmosphere_module, only: &
    determine_extended_atmos_bounds

  implicit none

  public :: &
    setup_radiation_variables, &
    reset_radiation_variables, &
    cleanup_radiation_variables, &
    setup_bugsrad_variables, &
    reset_bugsrad_variables, &
    cleanup_bugsrad_variables

  private

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    Frad, &
    Frad_SW_up, &
    Frad_LW_up, &
    Frad_SW_down, &
    Frad_LW_down

!$omp threadprivate(Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    radht_LW, &
    radht_SW, &
    Frad_SW, &
    Frad_LW

!$omp threadprivate(radht_LW, radht_SW, Frad_SW, Frad_LW)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    fdswcl, &
    fuswcl, &
    fdlwcl, &
    fulwcl

!$omp threadprivate(fdswcl, fuswcl, fdlwcl, fulwcl)

  integer, public :: &
    extended_atmos_bottom_level = 0, &
    extended_atmos_top_level = 0, &
    extended_atmos_range_size = 0, &
    lin_int_buffer = 0

!$omp threadprivate(extended_atmos_bottom_level, extended_atmos_top_level, &
!$omp              extended_atmos_range_size, lin_int_buffer)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    T_in_K,   &
    rcil,     &
    o3l

!$omp threadprivate(T_in_K, rcil, o3l)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    rsm_rad, &
    rcm_in_cloud_rad, &
    cloud_frac_rad, &
    ice_supersat_frac_rad

!$omp threadprivate(rsm_rad, rcm_in_cloud_rad, cloud_frac_rad, ice_supersat_frac_rad)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    radht_SW_rad, &
    radht_LW_rad

!$omp threadprivate(radht_SW_rad, radht_LW_rad)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    p_in_mb, &
    sp_humidity

!$omp threadprivate(p_in_mb, sp_humidity)

  real( kind = dp ), public, dimension(:,:), allocatable :: &
    Frad_uLW, &
    Frad_dLW, &
    Frad_uSW, &
    Frad_dSW

!$omp threadprivate(Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW)

contains

  subroutine setup_radiation_variables( ngrdcol, nzm, nzt )

    implicit none

    integer, intent(in) :: &
      ngrdcol, &
      nzm, &
      nzt

    integer :: &
      i, &
      k

    allocate( Frad(ngrdcol, nzm) )
    allocate( Frad_SW_up(ngrdcol, nzm) )
    allocate( Frad_LW_up(ngrdcol, nzm) )
    allocate( Frad_SW_down(ngrdcol, nzm) )
    allocate( Frad_LW_down(ngrdcol, nzm) )

    allocate( radht_LW(ngrdcol, nzt) )
    allocate( radht_SW(ngrdcol, nzt) )
    allocate( Frad_SW(ngrdcol, nzm) )
    allocate( Frad_LW(ngrdcol, nzm) )

    !$acc enter data create( Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down, &
    !$acc                    radht_LW, radht_SW, Frad_SW, Frad_LW )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        Frad(i,k) = 0._core_rknd
        Frad_SW_up(i,k) = 0._core_rknd
        Frad_LW_up(i,k) = 0._core_rknd
        Frad_SW_down(i,k) = 0._core_rknd
        Frad_LW_down(i,k) = 0._core_rknd
        Frad_SW(i,k) = 0._core_rknd
        Frad_LW(i,k) = 0._core_rknd
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        radht_LW(i,k) = 0._core_rknd
        radht_SW(i,k) = 0._core_rknd
      end do
    end do
    !$acc end parallel loop

  end subroutine setup_radiation_variables

  subroutine reset_radiation_variables()

    implicit none

    integer :: i, k

    if ( .not. allocated( Frad ) ) return

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, size( Frad, 2 )
      do i = 1, size( Frad, 1 )
        Frad(i,k) = 0._core_rknd
        Frad_SW_up(i,k) = 0._core_rknd
        Frad_LW_up(i,k) = 0._core_rknd
        Frad_SW_down(i,k) = 0._core_rknd
        Frad_LW_down(i,k) = 0._core_rknd
        Frad_SW(i,k) = 0._core_rknd
        Frad_LW(i,k) = 0._core_rknd
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, size( radht_LW, 2 )
      do i = 1, size( radht_LW, 1 )
        radht_LW(i,k) = 0._core_rknd
        radht_SW(i,k) = 0._core_rknd
      end do
    end do
    !$acc end parallel loop

  end subroutine reset_radiation_variables

  subroutine setup_bugsrad_variables( &
             ngrdcol, nzmax, zt_grid, zm_grid, zm_grid_spacing, p_in_Pa_zm )

    implicit none

    integer, intent(in) :: &
      ngrdcol, &
      nzmax

    real( kind = core_rknd ), dimension(nzmax-1), intent(in) :: &
      zt_grid, &
      zm_grid_spacing

    real( kind = core_rknd ), dimension(nzmax), intent(in) :: &
      zm_grid, &
      p_in_Pa_zm

    integer :: &
      extended_atmos_bottom_level_local, &
      extended_atmos_top_level_local, &
      extended_atmos_range_size_local, &
      lin_int_buffer_local, &
      rad_zt_dim, &
      rad_zm_dim, &
      i, &
      k

    call determine_extended_atmos_bounds( nzmax, zt_grid, zm_grid, zm_grid_spacing, p_in_Pa_zm, &
                                          extended_atmos_bottom_level_local, &
                                          extended_atmos_top_level_local, &
                                          extended_atmos_range_size_local, &
                                          lin_int_buffer_local )

    extended_atmos_bottom_level = extended_atmos_bottom_level_local
    extended_atmos_top_level = extended_atmos_top_level_local
    extended_atmos_range_size = extended_atmos_range_size_local
    lin_int_buffer = lin_int_buffer_local

    rad_zt_dim = (nzmax-1) + lin_int_buffer + extended_atmos_range_size
    rad_zm_dim = (nzmax-1) + lin_int_buffer + extended_atmos_range_size + 1

    allocate( T_in_K(ngrdcol, rad_zt_dim ) )
    allocate( rcil(ngrdcol, rad_zt_dim ) )
    allocate( o3l(ngrdcol, rad_zt_dim ) )

    allocate( rsm_rad(ngrdcol, rad_zt_dim ) )
    allocate( rcm_in_cloud_rad(ngrdcol, rad_zt_dim ) )
    allocate( cloud_frac_rad(ngrdcol, rad_zt_dim ) )
    allocate( ice_supersat_frac_rad(ngrdcol, rad_zt_dim ) )

    allocate( p_in_mb(ngrdcol, rad_zt_dim) )
    allocate( sp_humidity(ngrdcol, rad_zt_dim) )

    allocate( radht_SW_rad(ngrdcol, rad_zt_dim ) )
    allocate( radht_LW_rad(ngrdcol, rad_zt_dim ) )

    allocate( Frad_uLW(ngrdcol, rad_zm_dim ) )
    allocate( Frad_dLW(ngrdcol, rad_zm_dim ) )
    allocate( Frad_uSW(ngrdcol, rad_zm_dim ) )
    allocate( Frad_dSW(ngrdcol, rad_zm_dim ) )

    allocate( fdswcl(ngrdcol, rad_zm_dim) )
    allocate( fuswcl(ngrdcol, rad_zm_dim) )
    allocate( fdlwcl(ngrdcol, rad_zm_dim) )
    allocate( fulwcl(ngrdcol, rad_zm_dim) )

    do k = 1, rad_zt_dim
      do i = 1, ngrdcol
        T_in_K(i,k) = 0._dp
        rcil(i,k) = 0._dp
        o3l(i,k) = 0._dp
        rsm_rad(i,k) = 0._dp
        rcm_in_cloud_rad(i,k) = 0._dp
        cloud_frac_rad(i,k) = 0._dp
        ice_supersat_frac_rad(i,k) = 0._dp
        radht_SW_rad(i,k) = 0._dp
        radht_LW_rad(i,k) = 0._dp
        p_in_mb(i,k) = 0._dp
        sp_humidity(i,k) = 0._dp
      end do
    end do
    do k = 1, rad_zm_dim
      do i = 1, ngrdcol
        Frad_uLW(i,k) = 0._dp
        Frad_dLW(i,k) = 0._dp
        Frad_uSW(i,k) = 0._dp
        Frad_dSW(i,k) = 0._dp
        fdswcl(i,k) = 0._dp
        fuswcl(i,k) = 0._dp
        fdlwcl(i,k) = 0._dp
        fulwcl(i,k) = 0._dp
      end do
    end do

  end subroutine setup_bugsrad_variables

  subroutine reset_bugsrad_variables()

    implicit none

    integer :: &
      i, &
      k

    if ( .not. allocated( T_in_K ) ) return

    do k = 1, size( T_in_K, 2 )
      do i = 1, size( T_in_K, 1 )
        T_in_K(i,k) = 0._dp
        rcil(i,k) = 0._dp
        o3l(i,k) = 0._dp
        rsm_rad(i,k) = 0._dp
        rcm_in_cloud_rad(i,k) = 0._dp
        cloud_frac_rad(i,k) = 0._dp
        ice_supersat_frac_rad(i,k) = 0._dp
        p_in_mb(i,k) = 0._dp
        sp_humidity(i,k) = 0._dp
        radht_SW_rad(i,k) = 0._dp
        radht_LW_rad(i,k) = 0._dp
      end do
    end do

    do k = 1, size( Frad_uLW, 2 )
      do i = 1, size( Frad_uLW, 1 )
        Frad_uLW(i,k) = 0._dp
        Frad_dLW(i,k) = 0._dp
        Frad_uSW(i,k) = 0._dp
        Frad_dSW(i,k) = 0._dp
        fdswcl(i,k) = 0._dp
        fuswcl(i,k) = 0._dp
        fdlwcl(i,k) = 0._dp
        fulwcl(i,k) = 0._dp
      end do
    end do

  end subroutine reset_bugsrad_variables

  subroutine cleanup_radiation_variables()

    implicit none

    if ( allocated( Frad ) ) then
      !$acc exit data delete( Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down, &
      !$acc                       radht_LW, radht_SW, Frad_SW, Frad_LW )
      deallocate( Frad )
      deallocate( Frad_SW_up )
      deallocate( Frad_LW_up )
      deallocate( Frad_SW_down )
      deallocate( Frad_LW_down )
      deallocate( radht_LW )
      deallocate( radht_SW )
      deallocate( Frad_SW )
      deallocate( Frad_LW )
    end if

  end subroutine cleanup_radiation_variables

  subroutine cleanup_bugsrad_variables()

    implicit none

    extended_atmos_bottom_level = 0
    extended_atmos_top_level = 0
    extended_atmos_range_size = 0
    lin_int_buffer = 0

    if ( allocated( T_in_K ) ) then
      deallocate( T_in_K )
      deallocate( rcil )
      deallocate( o3l )
      deallocate( rsm_rad )
      deallocate( rcm_in_cloud_rad )
      deallocate( cloud_frac_rad )
      deallocate( ice_supersat_frac_rad )
      deallocate( p_in_mb )
      deallocate( sp_humidity )
      deallocate( radht_SW_rad )
      deallocate( radht_LW_rad )
      deallocate( Frad_uLW )
      deallocate( Frad_dLW )
      deallocate( Frad_uSW )
      deallocate( Frad_dSW )
      deallocate( fdswcl )
      deallocate( fuswcl )
      deallocate( fdlwcl )
      deallocate( fulwcl )
    end if

  end subroutine cleanup_bugsrad_variables

end module radiation_variables_module
