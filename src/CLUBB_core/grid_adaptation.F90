!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_adaptation_module


  use clubb_precision, only: &
      core_rknd ! Variable(s)

  use grid_class, only: &
      grid ! Type

  use constants_clubb, only: &
      one, fstdout, fstderr ! Constants

  use error_code, only: &
      clubb_at_least_debug_level

  implicit none

  public :: setup_gr_dycore, &
            remapping_matrix_zm_values, remapping_matrix_zt_values, &
            remap_vals_to_target, &
            normalize_grid_density, &
            adapt_grid, &
            calc_grid_dens, &
            clean_up_grid_adaptation_module

  private :: check_remap_conservation, check_remap_consistency, check_remap_monotonicity, &
             remapping_matrix, remap_vals_to_target_helper, &
             check_remap_consistency_w_vals, &
             check_mass_conservation, &
             check_vertical_integral_conservation, &
             check_remapped_val_for_monotonicity, &
             calc_mass_over_grid_intervals, vertical_integral_conserve_mass, &
             calc_integral, create_fixed_min_gr_dens_func, &
             normalize_min_grid_density, normalize_grid_density_helper, &
             create_grid_from_normalized_grid_density_func, &
             create_grid_from_grid_density_func, check_grid

  private ! Default Scoping

  real( kind = core_rknd ), parameter ::  &
    tol = 1.0e-6_core_rknd ! tolerance to check if to real numbers are equal

  !real( kind = core_rknd ), parameter :: &
  !  lambda = 0.3

  integer :: &
    fixed_min_gr_dens_idx ! number of elements in the vertical axis of fixed_min_gr_dens

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    fixed_min_gr_dens_z, &  ! z coordinates of the piecewise linear function [m]
    fixed_min_gr_dens       ! density values given at the z coordinates [# levs/m]

  integer :: &
    gr_dens_old_idx_global ! number of levels in gr_dens_old_global and gr_dens_old_z_global
  
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    gr_dens_old_global ! grid density from adaptation step before [# levs/m]

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    gr_dens_old_z_global ! grid density altitudes from adaptation step before [m]

  integer :: &
    Lscale_counter, &
    brunt_counter, &
    chi_counter

  real( kind = core_rknd ), dimension(:), allocatable :: &
    cumulative_Lscale ! cumulative Lscale sum [m]

  real( kind = core_rknd ), dimension(:), allocatable :: &
    cumulative_brunt ! cumulative brunt term sum [m]

  real( kind = core_rknd ), dimension(:), allocatable :: &
    cumulative_chi ! cumulative chi term sum [m]


  ! TODO remove those three variables
  integer :: &
    max_history_Lscale = 1000

  integer :: &
    ind_Lscale_history

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    Lscale_history ! Lscale values of the last iterations [m]





  real(core_rknd), parameter ::  D0_0                    =  0.0_core_rknd
 real(core_rknd), parameter ::  D1EM14                  =  1.0e-14_core_rknd
 real(core_rknd), parameter ::  D0_125                  =  0.125_core_rknd
 real(core_rknd), parameter ::  D0_1875                 =  0.1875_core_rknd
 real(core_rknd), parameter ::  D0_25                   =  0.25_core_rknd
 real(core_rknd), parameter ::  D0_5                    =  0.5_core_rknd
 real(core_rknd), parameter ::  D1_0                    =  1.0_core_rknd
 real(core_rknd), parameter ::  D1_5                    =  1.5_core_rknd
 real(core_rknd), parameter ::  D2_0                    =  2.0_core_rknd
 real(core_rknd), parameter ::  D3_0                    =  3.0_core_rknd
 real(core_rknd), parameter ::  D4_0                    =  4.0_core_rknd
 real(core_rknd), parameter ::  D5_0                    =  5.0_core_rknd
 real(core_rknd), parameter ::  D8_0                    =  8.0_core_rknd
 real(core_rknd), parameter ::  D12_0                   = 12.0_core_rknd






  contains

  !=============================================================================

  subroutine setup_gr_dycore( iunit, ngrdcol, grid_sfc, grid_top, gr )

    use grid_class, only: &
      read_grid_heights, &
      setup_grid

    use error_code, only: &
      clubb_fatal_error
    
    implicit none

    !--------------------- Input Variables ---------------------
    integer :: &
      iunit, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
        grid_sfc, &    ! grids surface height for the dycore grid [m]
        grid_top       ! grids highest level for the dycore grid [m]

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- Local Variables ---------------------
    integer :: i
    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer :: &
      grid_type, &
      err_code, &
      nzmax, &
      nlevel, &
      level_min

    character(len=32) :: & 
      zm_grid_fname,  & ! Path and filename of file for momentum level altitudes
      zt_grid_fname     ! Path and filename of file for thermodynamic level altitudes

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
        sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical :: l_implemented

    real( kind = core_rknd ), dimension(ngrdcol,73) :: & ! since dycore grid has 73 levels in total
      momentum_heights      ! Momentum level altitudes (file input)      [m]
    
    real( kind = core_rknd ), dimension(ngrdcol,73-1) :: &! since dycore grid has 73 levels in total
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      deltaz ! vertical grid spacing [m]

    !--------------------- Begin Code ---------------------

    nzmax = 73 ! since dycore grid has in total 73 (zm) levels
    grid_type = 3
    zm_grid_fname = '../input/grid/dycore.grd'
    zt_grid_fname = ''

    do i = 1, ngrdcol
      call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                              zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                              iunit, &                            ! Intent(in)
                              momentum_heights(i,:), &            ! Intent(out)
                              thermodynamic_heights(i,:), &       ! Intent(out)
                              err_code )                          ! Intent(inout)
      sfc_elevation(i) = 0.0_core_rknd
      deltaz(i) = 0.0_core_rknd
    end do

    if ( err_code == clubb_fatal_error ) then
      write(fstderr, *) "Error in read_grid_heights"
      error stop
    end if

    nlevel = nzmax
    do while( momentum_heights(1,nlevel) > grid_top(1) .and. nlevel > 1 )

      nlevel = nlevel - 1

    end do

    do i = 1, ngrdcol
      momentum_heights(i,nlevel) = grid_top(i)
    end do

    level_min = nlevel
    do while( momentum_heights(1,level_min) > grid_sfc(1) .and. level_min > 1 )

      level_min = level_min - 1

    end do

    do i = 1, ngrdcol
      momentum_heights(i,level_min) = grid_sfc(i)
    end do


    l_implemented = .false.
    
    call setup_grid( nlevel, ngrdcol, sfc_elevation, l_implemented, &   ! intent(in)
                     grid_type, deltaz, grid_sfc, grid_top, &           ! intent(in)
                     momentum_heights, thermodynamic_heights, &         ! intent(in)
                     gr, err_code )                                     ! intent(inout)

    if ( err_code == clubb_fatal_error ) then
      error stop "Error in CLUBB calling setup_grid"
    end if

  end subroutine setup_gr_dycore


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!                                REMAPPING
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine check_remap_conservation( R_ij, dim1, dim2, &
                                       levels_source, levels_target, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels )
    ! checks conservation with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (9)
    ! Defintion of conservation by Ullrich and Taylor: A remapping operator R is conservative if
    ! the global mass of any field is maintained across the remapping operation

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    real( kind = core_rknd ), dimension(dim2+1), intent(in) :: &
      levels_source ! the height of the levels in the source grid [m]

    real( kind = core_rknd ), dimension(dim1+1), intent(in) :: &
      levels_target ! the height of the levels in the target grid [m]

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: weighted_col_sum

    real( kind = core_rknd ), dimension(dim1) :: &
        J_t ! local weights on the target grid (mass over every grid cell) []

    real( kind = core_rknd ), dimension(dim2) :: &
        J_s ! local weights on the source grid (mass over every grid cell) []

    logical :: conservative

    !--------------------- Begin Code ---------------------
    conservative = .true.

    J_t = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         dim1+1, levels_target )
    
    J_s = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         dim2+1, levels_source )

    j = 1
    do while (conservative .and. j <= dim2)
      weighted_col_sum = 0
      do i = 1, dim1
        weighted_col_sum = weighted_col_sum + R_ij(i,j) * J_t(i)
      end do
      if (abs(weighted_col_sum - J_s(j)) > tol) then
        conservative = .false.
      end if
      j = j + 1
    end do

    if (.not. conservative) then
      write(fstderr,*) 'WARNING! The remapping operator is not conservative'
      error stop 'Operator should be conservative, something went wrong...'
    end if

  end subroutine check_remap_conservation

  subroutine check_remap_consistency( R_ij, dim1, dim2 )
    ! checks consistency with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (12)
    ! Defintion of consistency by Ullrich and Taylor: A remapping operator R is consistent if the 
    ! constant field is maintained across the remapping operation
    ! So for example, if we have a discrete function of only ones and apply the remapping operator, 
    ! then the result should also be constant one

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: row_sum

    logical :: l_consistent

    !--------------------- Begin Code ---------------------
    l_consistent = .true.

    i = 1
    do while (l_consistent .and. i <= dim1)
      row_sum = 0
      do j = 1, dim2
        row_sum = row_sum + R_ij(i,j)
      end do
      if (abs(row_sum - 1) > tol) then
        l_consistent = .false.
      end if
      i = i + 1
    end do

    if (.not. l_consistent) then
      write(fstderr,*) 'WARNING! The remapping operator is not consistent'
      error stop 'Operator should be consistent, something went wrong...'
    end if

  end subroutine check_remap_consistency

  subroutine check_remap_monotonicity( R_ij, dim1, dim2 )
    ! checks monotonicity with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (14)
    ! Defintion of monotonicity by Ullrich and Taylor: A remapping operator R is monotone if the 
    ! remapping operation cannot introduce additional global extrema

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: i, j

    logical :: l_monotone

    !--------------------- Begin Code ---------------------
    l_monotone = .true.

    i = 1
    do while (l_monotone .and. i <= dim1)
      do j = 1, dim2
        if (R_ij(i,j) < 0) then
          l_monotone = .false.
        end if
      end do
      i = i + 1
    end do

    if (.not. l_monotone) then
      write(fstderr,*) 'WARNING! The remapping operator is not monotone'
      error stop 'Operator should be monotone, something went wrong...'
    end if

  end subroutine check_remap_monotonicity

  function remapping_matrix( nlevel_source, nlevel_target, &
                             levels_source, levels_target, &
                             total_idx_rho_lin_spline, &
                             rho_lin_spline_vals, &
                             rho_lin_spline_levels )
    ! implements the remapping scheme proposed by Ullrich et al. in 
    !'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part II' 
    ! in formula (30), with the addition that omega incorporates in our case also rho, so it is not 
    ! the geometric region, but the mass

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nlevel_source, & ! number of levels in the target grid []
      nlevel_target    ! number of levels in the source grid []

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid [m]

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid [m]

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(nlevel_target-1,nlevel_source-1) :: &
      remapping_matrix ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: k, j

    real( kind = core_rknd ) :: omega_ov, omega_ov_upper, omega_ov_lower

    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      omega_ts ! densities of all intervals in the target grid [kg]

    real( kind = core_rknd ), dimension(2) :: &
      omega_ov_zm ! helper variable to hold density of the overlap region

    !--------------------- Begin Code ---------------------

    omega_ts = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              nlevel_target, levels_target )
    do k = 1, (nlevel_target-1)
        do j = 1, (nlevel_source-1)
            omega_ov_upper = min(levels_target(k+1), levels_source(j+1))
            omega_ov_lower = max(levels_target(k), levels_source(j))
            if (omega_ov_upper > omega_ov_lower) then
                omega_ov_zm(1) = omega_ov_lower
                omega_ov_zm(2) = omega_ov_upper
                omega_ov = sum( calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                                               rho_lin_spline_vals, &
                                                               rho_lin_spline_levels, &
                                                               2, omega_ov_zm ) )
            else
                omega_ov = 0
            end if
            remapping_matrix(k,j) = omega_ov/omega_ts(k)
        end do
    end do

    if ( clubb_at_least_debug_level( 2 ) ) then
      call check_mass_conservation( nlevel_source, nlevel_target, &
                                    levels_source, levels_target, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels )
      call check_remap_conservation( remapping_matrix, (nlevel_target-1), (nlevel_source-1), &
                                     levels_source, levels_target, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels)
      call check_remap_consistency( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
      call check_remap_monotonicity( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
    end if

  end function remapping_matrix

  subroutine remapping_matrix_zm_values( ngrdcol, &
                                         gr_source, gr_target, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         levels_source, levels_target, &
                                         R_ij )

    implicit none
    !--------------------- Input Variables ---------------------
    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      ngrdcol, &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_target%nzm,gr_source%nzm), intent(out) :: &
      R_ij ! matrix to apply to input values for remapping (R_ij) []

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzt+2), intent(out) :: &
      levels_source ! levels on source grid [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr_target%nzt+2), intent(out) :: &
      levels_target ! levels on target grid [m]

    !--------------------- Local Variables ---------------------
    integer :: i, j

    integer :: &
      nlevels_source, & ! number of levels in source grid []
      nlevels_target    ! number of levels in target grid []

    !--------------------- Begin Code ---------------------
    ! since we have the values given on the zm levels, we dont have surrounding zt levels for the
    ! first and last level, so for the first value (at the bottom) we have the cell given by
    ! [zm(1),zt(1)] and for the last value we have the cell given by [zt(n),zm(n+1)] if we have
    ! n value given (so zt(n) is the last zt level and zm(n+1) is the last zm level),
    ! so we basically assume that the value given on the first and last zm
    ! level aren't on the zm level itself, but we shift them in our artifically cotructed interval
    ! to the middle, but this doesn't result in problems regarding conservation, since the method
    ! assumes that the value is constant over the whole cell, and so does the vertical integral
    ! (as defined in the paper by Ullrich etc)

    nlevels_source = gr_source%nzt+2
    nlevels_target = gr_target%nzt+2

    do i = 1, ngrdcol

      levels_source(i,1) = gr_source%zm(i,1)
      do j=1, gr_source%nzt
          levels_source(i,j+1) = gr_source%zt(i,j)
      end do

      levels_source(i,nlevels_source) = gr_source%zm(i,gr_source%nzm)
      levels_target(i,1) = gr_target%zm(i,1)

      do j=1, gr_target%nzt
          levels_target(i,j+1) = gr_target%zt(i,j)
      end do

      levels_target(i,nlevels_target) = gr_target%zm(i,gr_target%nzm)

      R_ij(i,:,:) = remapping_matrix( nlevels_source, nlevels_target, &
                                      levels_source(i,:), levels_target(i,:), &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals(i,:), &
                                      rho_lin_spline_levels(i,:) )
    
    end do

  end subroutine remapping_matrix_zm_values

  subroutine remapping_matrix_zt_values( ngrdcol, &
                                         gr_source, gr_target, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         R_ij )

    implicit none
    !--------------------- Input Variables ---------------------
    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      ngrdcol, &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_target%nzt,gr_source%nzt), intent(out) :: &
      R_ij ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      R_ij(i,:,:) = remapping_matrix( gr_source%nzm, gr_target%nzm, &
                                      gr_source%zm(i,:), gr_target%zm(i,:), &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals(i,:), &
                                      rho_lin_spline_levels(i,:) )
    end do

  end subroutine remapping_matrix_zt_values

  function remap_vals_to_target_helper( ngrdcol, &
                                        nz_source, nz_target, &
                                        source_values, &
                                        R_ij )
    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      nz_source, & ! number of levels in the target grid []
      nz_target    ! number of levels in the source grid []

    real( kind = core_rknd ), dimension(ngrdcol,nz_source), intent(in) :: &
      source_values  ! given values on the source grid that should be remapped 
                     ! to the target grid

    real( kind = core_rknd ), dimension(ngrdcol,nz_target,nz_source), intent(in) :: &
      R_ij

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz_target) :: &
      remap_vals_to_target_helper ! remapped values with the dimension of the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k, j

    !--------------------- Begin Code ---------------------

    do i = 1, ngrdcol
      ! matrix vector multiplication
      do k = 1, (nz_target)
        remap_vals_to_target_helper(i,k) = 0
        do j = 1, (nz_source)
          remap_vals_to_target_helper(i,k) = remap_vals_to_target_helper(i,k) &
                                             + R_ij(i,k,j)*source_values(i,j)
        end do
      end do
    end do

  end function remap_vals_to_target_helper

  function remap_vals_to_target( ngrdcol, &
                                 nlevel_source, nlevel_target, &
                                 levels_source, levels_target, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 source_values, &
                                 iv, &
                                 R_ij, p_sfc )

    use constants_clubb, only: &
      grav
    
    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      nlevel_source, & ! number of levels in the target grid []
      nlevel_target    ! number of levels in the source grid []

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nlevel_source) :: &
      levels_source           ! altitudes for the source grid [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nlevel_target) :: &
      levels_target           ! altitudes for the target grid [m]

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_source-1), intent(in) :: &
      source_values  ! given values on the source grid that should be remapped 
                     ! to the target grid
    
    integer, intent(in) :: &
      iv ! -1: Winds
         !  0: positive definite scalars
         !  1: other variables

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_target-1,nlevel_source-1), intent(in) :: &
      R_ij

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nlevel_target-1) :: &
      remap_vals_to_target_flipped, &
      remap_vals_to_target ! remapped values with the dimension of the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_source-1) :: &
      mass_on_source_cells

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_target-1) :: &
      mass_on_target_cells

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_source) :: &
      pressure_levels_source

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_target) :: &
      pressure_levels_target

    real( kind = core_rknd ), dimension(ngrdcol,nlevel_source-1) :: &
      source_values_flipped

    integer :: &
      kord ! order: should be 4 for vertical remapping

    !--------------------- Begin Code ---------------------
    !p_sfc = 1029.e2
    kord = 4
    do i = 1, ngrdcol
      mass_on_source_cells(i,:) = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 nlevel_source, levels_source )
      mass_on_target_cells(i,:) = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 nlevel_target, levels_target )

      ! TODO check if grid need to start at 0, if p_sfc is given for 0 and we set lowest level to p_sfc
      pressure_levels_source(i,nlevel_source) = p_sfc(i)
      do k = 2, nlevel_source
        pressure_levels_source(i,nlevel_source-k+1) = pressure_levels_source(i,nlevel_source-k+2) - mass_on_source_cells(i,k-1)*grav
      end do

      ! TODO check if grid need to start at 0, if p_sfc is given for 0 and we set lowest level to p_sfc
      pressure_levels_target(i,nlevel_target) = p_sfc(i)
      do k = 2, nlevel_target
        pressure_levels_target(i,nlevel_target-k+1) = pressure_levels_target(i,nlevel_target-k+2) - mass_on_target_cells(i,k-1)*grav
      end do

      do k = 1, nlevel_source-1
        source_values_flipped(i,k) = source_values(i,nlevel_source-k)
      end do

      call map1_ppm( nlevel_source-1,   pressure_levels_source(i,:),   source_values_flipped(i,:),  &
                   nlevel_target-1,   pressure_levels_target(i,:),   remap_vals_to_target_flipped(i,:), &
                         0, 0, 1, 1, 1,                         &
                         1, 1, 1, iv, kord) ! TODO use named variables

      do k = 1, nlevel_target-1
        remap_vals_to_target(i,k) = remap_vals_to_target_flipped(i,nlevel_target-k)
      end do
    end do

    
    
    !call map1_ppm( nlevel_source-1,   levels_source,   source_values,  &
    !               nlevel_target-1,   levels_target,   remap_vals_to_target, &
    !                     0, 0, 1, 1, 1,                         &
    !                     1, 1, 1, 0, 3)

    ! TODO flip source_values and after remapping remap_vals_to_target
    !do i = 1, ngrdcol
    !  do k = 1, nlevel_source-1
    !    source_values_flipped(i,k) = source_values(i,nlevel_source-k)
    !  end do
    !end do

    !call map1_ppm( nlevel_source-1,   pressure_levels_source,   source_values_flipped,  &
    !               nlevel_target-1,   pressure_levels_target,   remap_vals_to_target_flipped, &
    !                     0, 0, 1, 1, 1,                         &
    !                     1, 1, 1, 0, 3)
    !
    !do i = 1, ngrdcol
    !  do k = 1, nlevel_target-1
    !    remap_vals_to_target(i,k) = remap_vals_to_target_flipped(i,nlevel_target-k)
    !  end do
    !end do

    !write(*,*) 'source values: ', source_values
    !write(*,*) 'target values: ', remap_vals_to_target
    !remap_vals_to_target = remap_vals_to_target_helper( ngrdcol, &
    !                                                    nlevel_source-1, &
    !                                                    nlevel_target-1, &
    !                                                    source_values, &
    !                                                    R_ij )

    if ( clubb_at_least_debug_level( 2 ) .and. .true. ) then ! TODO build check in again
      do i = 1, ngrdcol
        ! check conservation
        call check_vertical_integral_conservation( total_idx_rho_lin_spline, &                ! In
                                                   rho_lin_spline_vals(i,:), &                ! In
                                                   rho_lin_spline_levels(i,:), &              ! In
                                                   nlevel_source, nlevel_target, &            ! In
                                                   levels_source(i,:), levels_target(i,:), &  ! In
                                                   source_values(i,:), &                      ! In
                                                   remap_vals_to_target(i,:) )                ! In

        ! check consistency
        ! TODO adjust for ppm remapping
        call check_remap_consistency_w_vals( nlevel_source-1, nlevel_target-1, &            ! In
                                             R_ij(i,:,:) )                                  ! In

        ! check monotonicity
        ! check monotonicity if iv is 1
        ! ppm is monotone on the inner levels (from 3 to n-2) if iv=1 and kord=4 or higher but it also seems to be in most cases monotone for iv=0 and kord=3, the relevant variable is lmt in kmppm, lmt=0 is not always monotone, lmt=1 is always monotone for the inner levels, but unfortunately i havent found a configuration where all levels, including the boundarties are monotone
        if ( iv == 1 .and. kord >= 4 ) then
          call check_remapped_val_for_monotonicity( 1, &                                   ! In
                                                    nlevel_source-1, nlevel_target-1, &    ! In
                                                    source_values(i,:), &                  ! In
                                                    remap_vals_to_target(i,:) )            ! In
        end if
      end do
    end if

  end function remap_vals_to_target

  subroutine check_remap_consistency_w_vals( n_source_vals, n_target_vals, &
                                             R_ij )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      n_source_vals, & ! number of values on the target grid []
      n_target_vals    ! number of values on the source grid []

    real( kind = core_rknd ), dimension(n_target_vals,n_source_vals), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(n_source_vals) :: &
      ones

    real( kind = core_rknd ), dimension(1,n_target_vals) :: &
      output_ones

    integer :: k

    logical :: l_consistent

    !--------------------- Begin Code ---------------------
    do k = 1, (n_source_vals)
      ones(k) = 1
    end do

    output_ones = remap_vals_to_target_helper( 1, &
                                               n_source_vals, n_target_vals, &
                                               ones, &
                                               R_ij )

    l_consistent = .true.
    k = 1
    do while (l_consistent .and. k <= n_target_vals)
      if (abs( output_ones(1,k) - 1 ) > tol) then
        l_consistent = .false.
      end if
      k = k + 1
    end do

    if (.not. l_consistent) then
      write(fstderr,*) 'WARNING! The remap_vals_to_target_helper function is not consistent'
      write(fstderr,*) 'The output should be all ones, but was: ', output_ones(1,:)
    end if

  end subroutine check_remap_consistency_w_vals

  subroutine check_mass_conservation( nlevel_source, nlevel_target, &
                                      levels_source, levels_target, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nlevel_source, & ! number of levels in the target grid []
      nlevel_target    ! number of levels in the source grid []

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid [m]

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid [m]

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      sum_source_mass, &
      sum_target_mass

    real( kind = core_rknd ), dimension(nlevel_source-1) :: source_mass

    real( kind = core_rknd ), dimension(nlevel_target-1) :: target_mass

    !--------------------- Begin Code ---------------------

    ! check mass for grid intervals of values at zt levels
    source_mass = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                 rho_lin_spline_levels, &
                                                 nlevel_source, levels_source )

    target_mass = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                 rho_lin_spline_levels, &
                                                 nlevel_target, levels_target )

    sum_source_mass = sum( source_mass )
    sum_target_mass = sum( target_mass )

    if (abs(sum_source_mass - sum_target_mass) > tol) then

      write(fstderr,*) "WARNING! Mass for different grids was not the same."
      write(fstderr,*) "Mass was ", sum_target_mass, " instead of ", sum_source_mass
      error stop 'Mass should be conserved, something went wrong...'

    end if

  end subroutine check_mass_conservation

  subroutine check_vertical_integral_conservation( total_idx_rho_lin_spline, &
                                                   rho_lin_spline_vals, &
                                                   rho_lin_spline_levels, &
                                                   nlevel_source, nlevel_target, &
                                                   levels_source, levels_target, &
                                                   field_source, field_target )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    integer, intent(in) :: &
      nlevel_source, &       ! number of levels of the source grid []
      nlevel_target          ! number of levels of the target grid []

    real( kind = core_rknd ), intent(in), dimension(nlevel_source) :: &
      levels_source           ! altitudes for the source grid [m]

    real( kind = core_rknd ), intent(in), dimension(nlevel_target) :: &
      levels_target           ! altitudes for the target grid [m]

    real( kind = core_rknd ), intent(in), dimension(nlevel_source-1) :: &
      field_source            ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(nlevel_target-1) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      integral_source, &
      integral_target

    !--------------------- Begin Code ---------------------
    integral_source = vertical_integral_conserve_mass( total_idx_rho_lin_spline, &
                                                       rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       nlevel_source, levels_source, &
                                                       field_source )
 
    integral_target = vertical_integral_conserve_mass( total_idx_rho_lin_spline, &
                                                       rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       nlevel_target, levels_target, &
                                                       field_target )

    ! set tolerance to 1.0e-3 since there are some extremely large variables
    ! that are getting written to file
    if (abs(integral_target - integral_source) > 1.0e-3) then

      write(fstderr,*) "WARNING! Integral for field was not conserved."
      write(fstderr,*) "Integral was ", integral_target, " instead of ", integral_source

    end if

  end subroutine check_vertical_integral_conservation

  subroutine check_remapped_val_for_monotonicity( ngrdcol, &
                                                  nz_source, nz_target, &
                                                  field_source, field_target )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol, nz_source, nz_target

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, nz_source) :: &
      field_source           ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, nz_target) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k

    real( kind = core_rknd ) :: &
      source_maxval, & ! maximum value of the source fields
      source_minval    ! minimum value of the source fields

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      source_maxval = maxval(field_source(i,:))
      source_minval = minval(field_source(i,:))

      ! TODO only restrict the monotonicity test to the inner levels if we use ppm
      do k = 3, nz_target-2
        ! check if remapped target field is outside the region the global extrema
        ! of the source values and if deviation is greater than the tolerance
        if ( ( field_target(i,k) < source_minval &
               .and. abs(field_target(i,k) - source_minval ) > tol ) &
             .or. &
             ( field_target(i,k) > source_maxval &
               .and. abs(field_target(i,k) - source_maxval ) > tol ) ) then
          write(fstderr,*) 'WARNING! The remapped values are not monotone.'
          write(fstderr,*) 'The extrema on the source field were: ', source_minval, ' and ', &
                            source_maxval, ', but the remapped value is: ', field_target(i,k)
        end if

      end do

    end do
  end subroutine check_remapped_val_for_monotonicity

  function calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                          lin_spline_rho_vals, &
                                          lin_spline_rho_levels, &
                                          nlevel, grid_levels &
                                        ) result( mass_per_interval )

    ! Description:
    ! Calculate the mass over every interval (grid_levels(i) to grid_levels(i+1))
    ! of the grid levels (grid_levels).
    ! This function assumes a linear spline for rho, defined by the values of 
    ! rho (lin_spline_rho_vals) given at the altitudes (lin_spline_rho_levels).
    ! So with the assumption that the linear spline is the exact rho function, this function 
    ! computes the exact mass over each interval, such that any target grid with the same lowest 
    ! and highest level would give the same total mass, if the linear spline is the same.

    ! Returns an array with the dimension nlevel-1 where at each index the mass of that 
    ! interval is stored, starting at the bottom level.
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use interpolation, only: &
        lin_interpolate_two_points ! Procedure

    use grid_class, only:  &
        grid    ! Type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      total_idx_lin_spline, &  ! The total numer of indices of the spline points []
      nlevel                   ! The total numer of indices of the grid levels []

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values [m]
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    real( kind = core_rknd ), dimension(nlevel), intent(in) :: &
      grid_levels ! levels of grid over which the masses shoud be calculated [m]

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      val_below, & ! The rho value of the level or spline connection point below [m]
      level_below, & ! The altitude of the level or spline connection point below [m]
      mass_over_interval, & ! The total mass over one interval (zm(i) to zm(i+1)) in the grid [kg]
      upper_zm_level_rho

    real( kind = core_rknd ), dimension(nlevel-1) :: &
      mass_per_interval ! Array to store the mass of each interval
                        ! (grid_levels(i) to grid_levels(i+1)) [kg]

    integer :: i, j

    logical :: level_found

    !--------------------- Begin Code ---------------------
    ! check if highest value of lin spline is higher or equal to highest level of target grid
    if ( (grid_levels(nlevel) &
          - lin_spline_rho_levels(total_idx_lin_spline) ) > tol ) then
      write(fstderr,*) 'Wrong input: highest level of lin_spline_rho_levels must', &
                       ' be higher or equal to highest zm level of grid_levels'
      write(fstderr,*) 'Cannot compute the mass to get the remapping operator if the', & 
                       'density function is not defined over the whole grid region'
      error stop 'Cannot compute the remapping operator'
    end if

    j = 1
    level_found = .false.
    ! find first level and set initial val_below and level_below
    do while ( (.not. level_found) .and. (j <= total_idx_lin_spline) )
      if ( abs( grid_levels(1) - lin_spline_rho_levels(j) ) < tol ) then
        val_below = lin_spline_rho_vals(j)
        level_below = grid_levels(1)
        level_found = .true.
      else if ( grid_levels(1) < lin_spline_rho_levels(j) ) then ! j > 1
        if (j <= 1) then
          write(fstderr,*) 'Wrong input: lowest level of lin_spline_rho_levels must', &
                           ' be lower or equal to lowest zm level of grid_levels'
          write(fstderr,*) 'Cannot compute the mass to get the remapping operator if the', & 
                           'density function is not defined over the whole grid region'
          error stop 'Cannot compute the remapping operator'
        else
          val_below = lin_interpolate_two_points( grid_levels(1), &
                                                  lin_spline_rho_levels(j), &
                                                  lin_spline_rho_levels(j-1), &
                                                  lin_spline_rho_vals(j), &
                                                  lin_spline_rho_vals(j-1) )
          level_below = grid_levels(1)
          level_found = .true.
          j = j - 1
        end if
      end if
      j = j + 1
    end do
    ! find all spline connection points between the level_below and the next grid level
    ! in grid_levels and add up their masses, finally the rho value on the next grid level is
    ! interpolated and the last part of mass added to the total mass of that interval
    do i = 1, nlevel-1
      mass_over_interval = 0.0_core_rknd
      do while ( (lin_spline_rho_levels(j) < grid_levels(i+1)) &
                 .and. (j < total_idx_lin_spline) )
        mass_over_interval = mass_over_interval &
                             + (lin_spline_rho_levels(j) - level_below) &
                              *( val_below + lin_spline_rho_vals(j) )/2
        level_below = lin_spline_rho_levels(j)
        val_below = lin_spline_rho_vals(j)
        j = j + 1
      end do
      if ( abs( grid_levels(i+1) - lin_spline_rho_levels(j) ) < tol ) then 
        ! have approx equal altitude, so just copy value
        upper_zm_level_rho = lin_spline_rho_vals(j)
      else ! lin spline at j is above grid_levels(i+1) in grid_levels
        ! lin interpolation between j lin spline node and lin spine node j-1
        upper_zm_level_rho = lin_interpolate_two_points( grid_levels(i+1), &
                                                         lin_spline_rho_levels(j), &
                                                         lin_spline_rho_levels(j-1), &
                                                         lin_spline_rho_vals(j), &
                                                         lin_spline_rho_vals(j-1) )
      end if
      mass_over_interval = mass_over_interval &
                           + (grid_levels(i+1) - level_below) &
                            *( val_below + upper_zm_level_rho )/2

      mass_per_interval(i) = mass_over_interval

      val_below = upper_zm_level_rho
      level_below = grid_levels(i+1)
    end do
    return
  end function calc_mass_over_grid_intervals

  function vertical_integral_conserve_mass( total_idx_lin_spline, &
                                            lin_spline_rho_vals, &
                                            lin_spline_rho_levels, &
                                            nlevel, grid_levels, &
                                            field )

    ! Description:
    ! Computes the vertical integral. lin_spline_rho_vals and lin_spline_rho_levels must be
    ! of size total_idx_lin_spline and should start at the same index.
    ! field must be of size nlevel-1
    ! For more conditions to those parameters read description of calc_mass_over_grid_intervals

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only:  &
        grid    ! Type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      total_idx_lin_spline  ! The total numer of indices of the spline points []

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values [m]
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    integer, intent(in) :: & 
      nlevel  ! The total numer of indices for grid_levels []

    real( kind = core_rknd ), dimension(nlevel-1), intent(in) ::  &
      field      ! The field to be vertically averaged   [Units vary]

    real( kind = core_rknd ), dimension(nlevel), intent(in) ::  &
      grid_levels ! [m]

    ! Note:  The field and grid_levels points need to be arranged from
    !        lowest to highest in altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel-1) :: &
      mass_per_interval ! Array to store the mass of each interval of the target_grid [kg]

    real( kind = core_rknd ) :: &
      vertical_integral_conserve_mass ! [kg]

    !--------------------- Begin Code ---------------------
    ! Initializing vertical_integral_conserve_mass to avoid a compiler warning.
    vertical_integral_conserve_mass = 0.0_core_rknd

    ! Compute the integral.
    ! Multiply the field at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    mass_per_interval = calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                                       lin_spline_rho_vals, &
                                                       lin_spline_rho_levels, &
                                                       nlevel, grid_levels )
    vertical_integral_conserve_mass = sum( field * mass_per_interval )

    return
  end function vertical_integral_conserve_mass


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!                    GRID CONSTRUCTION AND ADAPTATION
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  function calc_integral( g_idx, g_x, g_y )
    ! Description:
    ! Computes the exact integral, assuming the given points build a piecewise linear function.
    ! g_x and g_y_idx must both have g_idx elements

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      g_idx  ! The total numer of indices of g_x and g_y []

    real( kind = core_rknd ), dimension(g_idx), intent(in) ::  &
      g_x,  &    ! the x coordinates of the connection points of the piecewise linear function [m]
      g_y        ! the y coordinates of the connection points of the piecewise linear function

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: &
      calc_integral ! The integral

    !--------------------- Begin Code ---------------------
    ! Initializing calc_integral to avoid a compiler warning.
    calc_integral = 0.0_core_rknd

    ! calculate each linear segments exact integral and add them up
    do i = 1, (g_idx-1)
        calc_integral = calc_integral + (g_y(i) + g_y(i+1))/2 * (g_x(i+1) - g_x(i))
    end do

    return

  end function calc_integral

  subroutine create_fixed_min_gr_dens_func( iunit, ngrdcol, &
                                            grid_sfc, grid_top )

    ! Description:
    ! Creates and allocates the fixed minimum grid density function from the given dycore grid.
    ! Takes the dycore inverse grid distance between zm levels.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      iunit, &
      ngrdcol

    real( kind = core_rknd ), intent(in) :: &
      grid_sfc, & ! altitude of the grids surface
      grid_top    ! altitude of the grids top

    !--------------------- Local Variables ---------------------
    integer :: i, j

    type( grid ) :: gr_dycore

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      grid_sfc_arr, &
      grid_top_arr

    !--------------------- Begin Code ---------------------
    if ( (.not. allocated( fixed_min_gr_dens_z )) .or. (.not. allocated( fixed_min_gr_dens )) ) then
      do i = 1, ngrdcol
        grid_sfc_arr(i) = grid_sfc
        grid_top_arr(i) = grid_top
      end do
      ! Initialize the dycore grid to use as a profile for minimum grid density function
      call setup_gr_dycore( iunit, ngrdcol, &              ! Intent(in)
                            grid_sfc_arr, grid_top_arr, &  ! Intent(in)
                            gr_dycore )                    ! Intent(out)
      fixed_min_gr_dens_idx = gr_dycore%nzt+2
    end if

    if ( .not. allocated( fixed_min_gr_dens_z ) ) then
      allocate( fixed_min_gr_dens_z(ngrdcol,gr_dycore%nzt+2) )
      do i = 1, ngrdcol
        fixed_min_gr_dens_z(i,1) = grid_sfc
        do j = 1, gr_dycore%nzt
          fixed_min_gr_dens_z(i,j+1) = gr_dycore%zt(i,j)
        end do
        fixed_min_gr_dens_z(i,gr_dycore%nzt+2) = grid_top
      end do
    end if

    if ( .not. allocated( fixed_min_gr_dens ) ) then
      allocate( fixed_min_gr_dens(ngrdcol,gr_dycore%nzt+2) )
      do i = 1, ngrdcol
        ! use linear extrapolation to calculate points on the outer zm levels
        fixed_min_gr_dens(i,1) = gr_dycore%invrs_dzt(i,1) &
                                 + (fixed_min_gr_dens_z(i,1) - gr_dycore%zt(i,1)) &
                                   /(gr_dycore%zt(i,2) - gr_dycore%zt(i,1)) &
                                    *(gr_dycore%invrs_dzt(i,2) - gr_dycore%invrs_dzt(i,1))

        if ( clubb_at_least_debug_level( 2 ) ) then
          if ( fixed_min_gr_dens(i,1) <= 0 ) then
            error stop 'Initial minimum grid density function needs to be positive.'
          end if
        end if
        
        do j = 1, gr_dycore%nzt
          fixed_min_gr_dens(i,j+1) = gr_dycore%invrs_dzt(i,j)
          if ( clubb_at_least_debug_level( 2 ) ) then
            if ( fixed_min_gr_dens(i,j+1) <= 0 ) then
              error stop 'Initial minimum grid density function needs to be positive.'
            end if
          end if
        end do

        ! use linear extrapolation to calculate points on the outer zm levels
        fixed_min_gr_dens(i,gr_dycore%nzt+2) = gr_dycore%invrs_dzt(i,gr_dycore%nzt-1) &
                                               + (fixed_min_gr_dens_z(i,gr_dycore%nzt+2) &
                                                  - gr_dycore%zt(i,gr_dycore%nzt-1)) &
                                                 /(gr_dycore%zt(i,gr_dycore%nzt) &
                                                   - gr_dycore%zt(i,gr_dycore%nzt-1)) &
                                                  *(gr_dycore%invrs_dzt(i,gr_dycore%nzt) &
                                                    - gr_dycore%invrs_dzt(i,gr_dycore%nzt-1))

        if ( clubb_at_least_debug_level( 2 ) ) then
          if ( fixed_min_gr_dens(i,gr_dycore%nzt+2) <= 0 ) then
            error stop 'Initial minimum grid density function needs to be positive.'
          end if
        end if

      end do
    end if
    
  end subroutine create_fixed_min_gr_dens_func

  subroutine normalize_min_grid_density( ngrdcol, &
                                         min_gr_dens_idx, &
                                         min_gr_dens_z, min_gr_dens, &
                                         lambda, num_levels, &
                                         min_gr_dens_norm_z, min_gr_dens_norm )

    ! Description:
    ! Takes the initial minimum grid density function and normalizes it.
    ! Normalized means in this case, that the integral over
    ! the grid region is between 0 (excluded) and num_levels-1 (included). What value
    ! should be taken inbetween is controlled with lambda.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      min_gr_dens_idx  ! number of levels in min_gr_dens_z

    real( kind = core_rknd ), dimension(ngrdcol,min_gr_dens_idx), intent(in) :: &
      min_gr_dens_z, &  ! z coordinates of the piecewise linear minimum grid density function [m]
      min_gr_dens       ! density values given at the z coordintas [# levs/m]

    real( kind = core_rknd ), intent(in) :: &
      lambda      ! number between 0 (excluded) and 1 (included) used as an adjustable factor

    integer, intent(in) :: &
      num_levels  ! number of wanted grid levels for the new grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,min_gr_dens_idx), intent(out) :: &
      min_gr_dens_norm_z, & ! the normalized minimal grid density z coordinates [m]
      min_gr_dens_norm      ! the normalized minimal grid density values given at
                            ! the z coordinates [# levs/m]

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: &
      integral, &
      norm_factor

    !--------------------- Begin Code ---------------------
    ! check if lambda is in (0,1]
    if ( lambda <= 0.0_core_rknd .or. lambda > one ) then
      error stop 'lambda needs to be between 0 (excluded) and 1 (included)'
    end if

    do i = 1, ngrdcol
      ! calculate factor to normalize minimum grid density function
      integral = calc_integral( min_gr_dens_idx, min_gr_dens_z(i,:), min_gr_dens(i,:) )
      norm_factor = lambda*(num_levels - 1)/integral

      ! normalize minimum grid density function with norm_factor and copy z coordinates
      do j = 1, (min_gr_dens_idx)
        min_gr_dens_norm_z(i,j) = min_gr_dens_z(i,j)
        min_gr_dens_norm(i,j) = norm_factor*min_gr_dens(i,j)
      end do
    end do

    if ( clubb_at_least_debug_level( 2 ) ) then

      ! check if integral is lambda*(n-1) as wanted
      integral = calc_integral(min_gr_dens_idx, min_gr_dens_norm_z, min_gr_dens_norm)
      if ( abs(integral - (lambda*(num_levels-1))) > tol ) then
        write(fstderr,*) 'Warning! Integral in normalize_min_grid_density', &
                         ' should be something different.'
        error stop 'Something went wrong, integral is different than it should be.'
      end if

      ! check if normalized function is always >0
      do i = 1, ngrdcol
        do j = 1, min_gr_dens_idx
          if ( min_gr_dens_norm(i,j) <= 0.0_core_rknd ) then
            error stop 'Minimum grid density function in needs to be positive.'
          end if
        end do
      end do

    end if
    
  end subroutine normalize_min_grid_density

  subroutine normalize_grid_density_helper( ngrdcol, &
                                            gr_dens_idx, &
                                            gr_dens_z, gr_dens, &
                                            min_gr_dens_norm, &
                                            lambda, num_levels, &
                                            gr_dens_norm_z, gr_dens_norm )

    ! Description:
    ! Takes the linear piecewise grid density function (gr_dens_z, gr_dens) each of size
    ! gr_dens_idx and normalizes this density, such that the integral is the
    ! number of desired grid levels (num_levels)-1 and the minimum is above the minimal
    ! density function (gr_dens_z,min_gr_dens_norm)

    ! Note: The minimum grid density function needs to be normalized before with
    !       normalize_min_grid_density using it here.

    ! Note: The initial grid density function needs to be non-negative.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use interpolation, only: &
      zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      gr_dens_idx, &           ! number of elements in gr_dens_z, gr_dens, min_gr_dens_norm_z
                               ! and min_gr_dens_norm []
      num_levels               ! desired number of grid levels []

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) :: &
      gr_dens_z, & ! the grid density z coordinates [m]
      gr_dens      ! the grid density given at the z coordinates [# levs/meter]
      ! gr_dens_z and gr_dens should be ordered from the bottom to the top

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) :: &
      min_gr_dens_norm  ! the minimum grid density given at the z coordinates of gr_dens_z
                        ! [# levs/meter]
      ! min_gr_dens_norm should be ordered from the bottom to the top
      ! Note: min_gr_dens_norm should be given on the gr_dens_z altitudes

    real( kind = core_rknd ), intent(in) :: &
      lambda   ! a factor for how close you want to get to an equidistant grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(out) :: &
      gr_dens_norm_z, & ! the normalized grid density z coordinates [m]
      gr_dens_norm      ! the normalized grid density values given at the z coordinates [# levs/m]

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: &
      integral_gr_dens, &
      integral_min_gr_dens, &
      integral_test, &
      norm_factor, &
      shift, &
      grid_sfc, &
      grid_top

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      ! copy z coordinates and shift initial function down, such that lowest point is zero
      do j = 1, gr_dens_idx
        gr_dens_norm_z(i,j) = gr_dens_z(i,j)
        ! if grid density profile is too noise consider removing this shift and just
        ! copy value instead
        gr_dens_norm(i,j) = gr_dens(i,j) !- minval(gr_dens(i,:))
      end do

      ! calculate integrals and check
      integral_gr_dens = calc_integral( gr_dens_idx, gr_dens_norm_z(i,:), gr_dens_norm(i,:) )
      integral_min_gr_dens = calc_integral( gr_dens_idx, &
                                            gr_dens_norm_z(i,:), min_gr_dens_norm(i,:) )

      if ( abs(integral_min_gr_dens - (lambda*(num_levels-1))) > tol ) then
        write(fstderr,*) 'Warning! The minimum grid density function has not the correct integral.'
        error stop 'Normalize the minimum grid density function before using it.'
      end if

      if ( integral_gr_dens < 1.0e-8 ) then
        ! this happens if initial grid density function was constant
        ! in that case we get the normalized minimum grid density function shifted,
        ! such that the inetgral condition is fulfilled

        ! calculate shift value
        grid_sfc = gr_dens_norm_z(1,1)
        grid_top = gr_dens_norm_z(1,gr_dens_idx)
        shift = (num_levels-1 - integral_min_gr_dens)/(grid_top - grid_sfc)

        ! normalize function
        do j = 1, gr_dens_idx
          gr_dens_norm(i,j) = min_gr_dens_norm(i,j) + shift
        end do

      else
        ! calculate factor for normalization
        norm_factor = (num_levels-1 - integral_min_gr_dens)/integral_gr_dens

        ! normalize function
        do j = 1, gr_dens_idx
          gr_dens_norm(i,j) = norm_factor*gr_dens_norm(i,j) + min_gr_dens_norm(i,j)
        end do

      end if

      ! run checks
      if ( clubb_at_least_debug_level( 2 ) ) then

        ! check if integral is (n-1) as wanted
        integral_test = calc_integral( gr_dens_idx, gr_dens_norm_z(i,:), gr_dens_norm(i,:) )
        if ( abs(integral_test - (num_levels-1)) > tol ) then
          write(fstderr,*) 'Warning! Integral in normalize_grid_density_helper', &
                           ' should be something different.'
          error stop 'Something went wrong, integral is different than it should be.'
        end if

        ! check if normalized function is always >=minimum_density
        do j = 1, gr_dens_idx
          if ( gr_dens_norm(i,j) < min_gr_dens_norm(i,j) .and. &
               (min_gr_dens_norm(i,j) - gr_dens_norm(i,j)) > tol ) then
            write(fstderr,*) 'Normalized function was below minimum grid density function.' 
            write(fstderr,*) 'Normalized minimum grid density function is ', &
                              min_gr_dens_norm(i,j), &
                             ' while normalized grid density function is ', gr_dens_norm(i,j)
            error stop 'Normalized function should be above or equal to the minimum density.'
          end if
        end do
      
      end if

    end do

  end subroutine normalize_grid_density_helper

  function decide_if_grid_adapt_helper( ngrdcol, &
                                        gr_dens_old_idx, &
                                        gr_dens_old_z, &
                                        gr_dens_old, &
                                        gr_dens_new_idx, &
                                        gr_dens_new_z, &
                                        gr_dens_new, &
                                        threshold ) &
                                      result (l_adapt_grid)
    ! Description:
    ! Checks how different the two density functions are and returns depending on the threshold
    ! if the grid should be adapted or not.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use interpolation, only: &
        zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      gr_dens_old_idx, &  ! The total numer of indices of gr_dens_old_z
                          ! and gr_dens_old []
      gr_dens_new_idx     ! The total numer of indices of gr_dens_new_z
                          ! and gr_dens_new []

    real( kind = core_rknd ), dimension(ngrdcol, gr_dens_old_idx), intent(in) ::  &
      gr_dens_old_z,  &  ! the z coordinates of the connection points of the
                         ! piecewise linear grid density function from the adaptation
                         ! step before [m]
      gr_dens_old        ! the density values at the given z coordinates of the connection
                         ! points of the piecewise linear grid density
                         ! function from the adaptation step before [# levs/meter]

    real( kind = core_rknd ), dimension(ngrdcol, gr_dens_new_idx), intent(in) ::  &
      gr_dens_new_z,  &  ! the z coordinates of the connection points of the
                         ! piecewise linear grid density function from the current adaptation
                         ! step [m]
      gr_dens_new        ! the density values at the given z coordinates of the connection
                         ! points of the piecewise linear grid density
                         ! function from the current adaptation step [# levs/meter]

    real( kind = core_rknd ), intent(in) ::  &
      threshold   ! threshold which decides if grid should be adapted or not []

    !--------------------- Output Variable ---------------------
    logical, dimension(ngrdcol) :: &
      l_adapt_grid ! true or false, whether grid should be adapted or not

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(gr_dens_old_idx + gr_dens_new_idx) ::  &
      common_grid

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      gr_dens_old_interp

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      gr_dens_new_interp

    integer :: i, j, k, l

    logical :: l_adapt_grid_tmp

    real( kind = core_rknd ) :: &
      mean_gr_dens_old_interp, &
      mean_gr_dens_new_interp

    real( kind = core_rknd ) :: &
      sum

    real( kind = core_rknd ) :: &
      max_dist

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      l_adapt_grid_tmp = .false.

      ! create a common grid for both function
      j = 0
      k = 1
      l = 1
      do while ( k < gr_dens_old_idx .and. l < gr_dens_new_idx )
        j = j + 1
        if ( gr_dens_old_z(i,k) < gr_dens_new_z(i,l) - tol ) then
          common_grid(j) = gr_dens_old_z(i,k)
          k = k + 1
        else if ( gr_dens_new_z(i,l) < gr_dens_old_z(i,k) - tol ) then
          common_grid(j) = gr_dens_new_z(i,l)
          l = l + 1
        else
          ! both have the a level with the same altitude within the tolerance
          common_grid(j) = ( gr_dens_old_z(i,k) + gr_dens_new_z(i,l) )/2
          k = k + 1
          l = l + 1
        end if
      end do

      ! now the actual number of levels in the common grid is j which
      ! might be smaller than gr_dens_old_idx + gr_dens_new_idx
      allocate( gr_dens_old_interp(j) )
      allocate( gr_dens_new_interp(j) )

      ! interpolate both density functions to the common grid
      gr_dens_old_interp = zlinterp_fnc( j, gr_dens_old_idx, &
                                         common_grid, gr_dens_old_z, &
                                         gr_dens_old )
      gr_dens_new_interp = zlinterp_fnc( j, gr_dens_new_idx, &
                                         common_grid, gr_dens_new_z, &
                                         gr_dens_new )

      ! take average absolute squared distance as measure
      !sum = 0
      !do k = 1, j
      !  sum = sum + ( gr_dens_old_interp(k) - gr_dens_new_interp(k) )**2
      !end do

      ! other option for grid adaptation trigger
      ! TODO remove again if it doesnt work
      !sum = 0
      !do k = 1, j
      !  sum = sum + ( gr_dens_old_interp(k) - gr_dens_new_interp(k) )**2/gr_dens_old_interp(k)
      !end do

      ! other option for grid adaptation trigger
      ! TODO remove again if it doesnt work
      sum = 0
      do k = 1, j
        sum = sum + ( gr_dens_old_interp(k) - gr_dens_new_interp(k) )**2/gr_dens_old_interp(k)**2
      end do
      
      
      write(*,*) 'difference: ', sum/j
      write(*,*) 'threshold: ', threshold
      if ( sum/j > threshold ) then
        l_adapt_grid_tmp = .true.
      else
        write(*,*) '-----------------not adapted'
      end if

      !! take maximum absolute distance as measure
      !max_dist = abs( gr_dens_old_interp(1) - gr_dens_new_interp(1) )
      !do k = 2, j
      !  if ( abs( gr_dens_old_interp(k) - gr_dens_new_interp(k) ) > max_dist ) then
      !    max_dist = abs( gr_dens_old_interp(k) - gr_dens_new_interp(k) )
      !  end if
      !end do
!
      !if ( max_dist > threshold ) then
      !  l_adapt_grid_tmp = .true.
      !end if

      !! take average weighted squared distance as measure
      !sum = 0
      !do k = 1, j
      !  sum = sum + ( (gr_dens_old_interp(k) - gr_dens_new_interp(k)) &
      !                * 1/((gr_dens_old_interp(k) + gr_dens_new_interp(k)/2)) )**2
      !end do
      !write(*,*) 'difference: ', sum/j
      !write(*,*) 'threshold: ', threshold
      !if ( sum/j > threshold ) then
      !  l_adapt_grid_tmp = .true.
      !end if
!


      deallocate( gr_dens_old_interp )
      deallocate( gr_dens_new_interp )

      l_adapt_grid(i) = l_adapt_grid_tmp
    end do

    return

  end function decide_if_grid_adapt_helper

  function decide_if_grid_adapt( ngrdcol, &
                                 gr_dens_new_idx, &
                                 gr_dens_new_z, &
                                 gr_dens_new, &
                                 threshold ) &
                               result (l_adapt_grid)
    ! Description:
    ! Checks how different the two density functions are and returns depending on the threshold
    ! if the grid should be adapted or not.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use interpolation, only: &
        zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      gr_dens_new_idx     ! The total numer of indices of gr_dens_new_z
                          ! and gr_dens_new []

    real( kind = core_rknd ), dimension(ngrdcol, gr_dens_new_idx), intent(in) ::  &
      gr_dens_new_z,  &  ! the z coordinates of the connection points of the
                         ! piecewise linear grid density function from the current adaptation
                         ! step [m]
      gr_dens_new        ! the density values at the given z coordinates of the connection
                         ! points of the piecewise linear grid density
                         ! function from the current adaptation step [# levs/meter]

    real( kind = core_rknd ), intent(in) ::  &
      threshold   ! threshold which decides if grid should be adapted or not []

    !--------------------- Output Variable ---------------------
    logical, dimension(ngrdcol) :: &
      l_adapt_grid ! true or false, whether grid should be adapted or not

    !--------------------- Begin Code ---------------------

    l_adapt_grid = decide_if_grid_adapt_helper( ngrdcol, &
                                                gr_dens_old_idx_global, &
                                                gr_dens_old_z_global, &
                                                gr_dens_old_global, &
                                                gr_dens_new_idx, &
                                                gr_dens_new_z, &
                                                gr_dens_new, &
                                                threshold )    

    return

  end function decide_if_grid_adapt

  function create_grid_from_normalized_grid_density_func( num_levels, &
                                                          gr_dens_norm_idx, &
                                                          gr_dens_norm_z, &
                                                          gr_dens_norm ) &
                                                        result (grid_heights)
    ! Description:
    ! Creates the grid from the normalized grid density function following the equi-distribution
    ! principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      gr_dens_norm_idx, &  ! The total numer of indices of gr_dens_norm_z and gr_dens_norm []
      num_levels           ! number of levels the new grid should have []

    real( kind = core_rknd ), dimension(gr_dens_norm_idx), intent(in) ::  &
      gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized piecewise
                          ! linear grid density function [m]
      gr_dens_norm        ! the density values at the given z coordinates of the connection points
                          ! of the normalized piecewise linear grid density function [# levs/meter]
      ! The grid density function needs to be normalized with the number of grid levels num_levels,
      ! before the function can be applied here

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    !--------------------- Local Variables ---------------------
    integer :: i, &
               prev_x_ind  ! the index of the last used x coordinate []

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid [m]
      grid_top, &               ! the highest point of the grid  [m]
      area_up_to_prev_x, &      ! the integral of the function up to the last used x coordinate
                                ! [# levs]
      new_x_area_increment, &   ! integral from grid density function on [prev_x_ind, prev_x_ind+1]
                                ! [# levs]
      desired_area_up_to_ith_level, & ! the desired area up to the ith grid level, should be exactly
                                      ! i, since the grid density function is normalized so that the
                                      ! integral of the whole function is num_levels-1 and since we
                                      ! follow the equi-distribution approach [# levs]
      A, grid_level, slope, intercept, p_pq_formula, q_pq_formula

    logical :: grid_level_found

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_norm_z(1)
    grid_top = gr_dens_norm_z(gr_dens_norm_idx)

    grid_heights(1) = grid_sfc

    ! Initializing calc_integral to avoid a compiler warning.
    prev_x_ind = 1
    area_up_to_prev_x = 0.0_core_rknd

    do i = 2, num_levels-1
        grid_level_found = .false.
        do while (.not. grid_level_found)
            new_x_area_increment = (gr_dens_norm(prev_x_ind) + gr_dens_norm(prev_x_ind+1))/2 &
                                   * (gr_dens_norm_z(prev_x_ind+1) - gr_dens_norm_z(prev_x_ind))
            desired_area_up_to_ith_level = i-1
            if (area_up_to_prev_x + new_x_area_increment >= desired_area_up_to_ith_level) then
                A = desired_area_up_to_ith_level - area_up_to_prev_x
                ! check if gr_dens_norm(prev_x_ind+1) == gr_dens_norm(prev_x_ind)
                if (abs(gr_dens_norm(prev_x_ind+1) - gr_dens_norm(prev_x_ind)) < tol) then
                    ! gr_dens_norm(prev_x_ind+1) == gr_dens_norm(prev_x_ind)
                    grid_level = A/gr_dens_norm(prev_x_ind) + gr_dens_norm_z(prev_x_ind)
                else
                    ! calculate slope and intercept of the linear function
                    slope = (gr_dens_norm(prev_x_ind+1) - gr_dens_norm(prev_x_ind)) &
                            /(gr_dens_norm_z(prev_x_ind+1) - gr_dens_norm_z(prev_x_ind))
                    intercept = gr_dens_norm(prev_x_ind) - gr_dens_norm_z(prev_x_ind)*slope
                    ! calculate p and q of second-order polynomial to solve for x
                    ! solve x^2+px+q=0 with x=-p/2(+/-)sqrt((p/2)^2-q)
                    p_pq_formula = 2*intercept/slope
                    q_pq_formula = -1*gr_dens_norm_z(prev_x_ind)**2 &
                                   - 2/slope*intercept*gr_dens_norm_z(prev_x_ind) - 2/slope*A
                    ! since we have second-order polynomial, we have two roots, so we need to
                    ! identify the right solution - use pq formula to get solutions
                    grid_level = -1*p_pq_formula/2 + sqrt((p_pq_formula/2)**2 - q_pq_formula)

                    ! TODO maybe set grid_level directly to gr_dens_norm_z if it is the same with tol
                    if ( ( grid_level < gr_dens_norm_z(prev_x_ind) &
                           .and. (gr_dens_norm_z(prev_x_ind) - grid_level) > tol ) &
                         .or. &
                         ( grid_level > gr_dens_norm_z(prev_x_ind+1) &
                           .and. (grid_level - gr_dens_norm_z(prev_x_ind+1)) > tol ) & 
                         .or. (grid_level < grid_heights(i-1) ) ) then
                        ! first solution was not the one we were looking for
                        grid_level = -1*p_pq_formula/2 - sqrt((p_pq_formula/2)**2 - q_pq_formula)
                        
                        if ( ( grid_level < gr_dens_norm_z(prev_x_ind) &
                               .and. (gr_dens_norm_z(prev_x_ind) - grid_level) > tol ) &
                             .or. &
                             ( grid_level > gr_dens_norm_z(prev_x_ind+1) &
                               .and. (grid_level - gr_dens_norm_z(prev_x_ind+1)) > tol ) & 
                             .or. (grid_level < grid_heights(i-1) ) ) then
                            write(fstderr,*) "None of the two solutions works. ", &
                                             "Something went wrong."
                            error stop 'Something went wrong with generating the grid'
                        end if
                    end if
                end if

                grid_heights(i) = grid_level
                grid_level_found = .true.
            else
                area_up_to_prev_x = area_up_to_prev_x + new_x_area_increment
                prev_x_ind = prev_x_ind + 1
            end if
        end do
    end do

    grid_heights(num_levels) = grid_top

    return

  end function create_grid_from_normalized_grid_density_func

  subroutine normalize_grid_density( ngrdcol, &
                                     iunit, itime, &
                                     gr_dens_idx, &
                                     gr_dens_z, gr_dens, &
                                     lambda, &
                                     num_levels, &
                                     norm_min_grid_dens, &
                                     norm_grid_dens )
    ! Description:
    ! Creates the grid from some unnormalized  minimum density function and the unnormalized
    ! grid density function following the equidistribution principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use interpolation, only: &
      zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      iunit, &
      itime, &
      gr_dens_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have []

    real( kind = core_rknd ) :: &
      lambda   ! a factor for defining how close you want to get to an equidistant grid

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) ::  &
      gr_dens_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function [m]
      gr_dens         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function
                      ! [# levs/meter]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(out) :: &
      norm_min_grid_dens, & ! the density at the given z coordinates of the connection points
                            ! of the normalized piecewise linear grid density function
                            ! [# levs/meter]
      norm_grid_dens        ! the density at the given z coordinates of the connection points of the
                            ! normalized piecewise linear grid density function [# levs/meter]

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx, i, j

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens        ! the density at the given z coordinates of the connection points
                         ! of the normalized piecewise linear grid density function
                         ! [# levs/meter]

    ! TODO remove this variable
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens_norm_z      ! the z coordinates of the connection points of the normalized
                              ! piecewise linear grid density function [m]

    ! TODO remove this variable
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      gr_dens_norm_z      ! the z coordinates of the connection points of the normalized piecewise
                          ! linear grid density function [m]

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_z(1,1)
    grid_top = gr_dens_z(1,gr_dens_idx)

    call create_fixed_min_gr_dens_func( iunit+1, ngrdcol, &
                                        grid_sfc, grid_top )

    ! set the minimum grid density profile to be the linear piecewise function of the
    ! original minimum grid density function evaluated at the current grid levels
    ! this way the function can be executed with all checks and works out exactly
    do i = 1, ngrdcol
      min_gr_dens(i,:) = zlinterp_fnc( gr_dens_idx, fixed_min_gr_dens_idx, &
                                       gr_dens_z(i,:), fixed_min_gr_dens_z(i,:), &
                                       fixed_min_gr_dens(i,:) )
    end do

    ! TODO remove min_gr_dens_norm_z as output in normalize_min_grid_density
    ! normalize the minimum grid density function
    call normalize_min_grid_density( ngrdcol, &                        ! Intent(in)
                                     gr_dens_idx, &                    ! Intent(in)
                                     gr_dens_z, min_gr_dens, &         ! Intent(in)
                                     lambda, num_levels, &             ! Intent(in)
                                     min_gr_dens_norm_z, &             ! Intent(out)
                                     norm_min_grid_dens )              ! Intent(out)

    ! TODO remove gr_dens_norm_z as output in normalize_grid_density_helper
    ! normalize the grid density function with the normalized minimum grid density function
    call normalize_grid_density_helper( ngrdcol, &                               ! Intent(in)
                                        gr_dens_idx, &                           ! Intent(in)
                                        gr_dens_z, gr_dens, &                    ! Intent(in)
                                        norm_min_grid_dens, &                    ! Intent(in)
                                        lambda, num_levels, &                    ! Intent(in)
                                        gr_dens_norm_z, norm_grid_dens )         ! Intent(out)

    write(iunit, *) 'gr_dens_z', itime, gr_dens_norm_z
    write(iunit, *) 'gr_dens', itime, norm_grid_dens
    write(iunit, *) 'min_gr_dens_z', itime, min_gr_dens_norm_z
    write(iunit, *) 'min_gr_dens', itime, norm_min_grid_dens

    return

  end subroutine normalize_grid_density

  subroutine create_grid_from_normalized_grid_density( ngrdcol, &
                                                       gr_dens_norm_idx, &
                                                       gr_dens_norm_z, gr_dens_norm, &
                                                       min_gr_dens_norm_z, min_gr_dens_norm, &
                                                       num_levels, &
                                                       threshold, &
                                                       grid_heights, &
                                                       l_adapt_grid )
    ! Description:
    ! Creates the grid from some unnormalized  minimum density function and the unnormalized
    ! grid density function following the equidistribution principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use interpolation, only: &
      zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      gr_dens_norm_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have []

    real( kind = core_rknd ), intent(in) :: &
      threshold ! threshold which decides if grid should be adapted or not

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(in) ::  &
      min_gr_dens_norm_z,  & 
      min_gr_dens_norm,  & 
      gr_dens_norm_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function [m]
      gr_dens_norm         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function
                      ! [# levs/meter]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_levels), intent(out) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    logical, dimension(ngrdcol), intent(out) :: &
      l_adapt_grid ! whether the grid should be adapted or not

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx, i, j

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code --------------------- ! TODO fofo
    grid_sfc = gr_dens_norm_z(1,1)
    grid_top = gr_dens_norm_z(1,gr_dens_norm_idx)

    ! decide if grid should be adapted or not
    if ( .not. allocated( gr_dens_old_z_global ) .and. .not. allocated( gr_dens_old_global ) ) then
      error stop 'gr_dens_old_z_global and gr_dens_old_global were not allocated...'
    end if
    l_adapt_grid = decide_if_grid_adapt_helper( ngrdcol, &
                                                gr_dens_old_idx_global, &
                                                gr_dens_old_z_global, &
                                                gr_dens_old_global, &
                                                gr_dens_norm_idx, &
                                                gr_dens_norm_z, &
                                                gr_dens_norm, &
                                                threshold )

    ! create grid from normalized grid density function
    do i = 1, ngrdcol
      if ( l_adapt_grid(i) ) then
        ! store new normalized grid density as old grid density for next adaptation iteration
        gr_dens_old_z_global(i,:) = gr_dens_norm_z(i,:)
        gr_dens_old_global(i,:) = gr_dens_norm(i,:)
        ! only create grid if it should actually be adapted
        grid_heights(i,:) = create_grid_from_normalized_grid_density_func( num_levels, &
                                                                           gr_dens_norm_idx, &
                                                                           gr_dens_norm_z(i,:), &
                                                                           gr_dens_norm(i,:) )
      end if
    end do

    if ( clubb_at_least_debug_level( 2 ) ) then

        grid_heights_idx = size(grid_heights, 2)

        do i = 1, ngrdcol
          if ( l_adapt_grid(i) ) then
            ! only check if grid was actually adapted
            call check_grid( grid_heights_idx, grid_heights(i,:), &            ! Intent(in)
                             num_levels, &                                     ! Intent(in)
                             gr_dens_norm_idx, &                                    ! Intent(in)
                             min_gr_dens_norm_z(i,:), min_gr_dens_norm(i,:), & ! Intent(in)
                             grid_sfc, grid_top )                              ! Intent(in)
          end if
        end do

    end if 

    return

  end subroutine create_grid_from_normalized_grid_density

  subroutine create_grid_from_grid_density_func( ngrdcol, &
                                                 iunit, itime, &
                                                 gr_dens_idx, &
                                                 gr_dens_z, gr_dens, &
                                                 lambda, &
                                                 num_levels, &
                                                 threshold, &
                                                 grid_heights, &
                                                 l_adapt_grid )
    ! Description:
    ! Creates the grid from some unnormalized  minimum density function and the unnormalized
    ! grid density function following the equidistribution principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use interpolation, only: &
      zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      iunit, &
      itime, &
      gr_dens_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have []

    real( kind = core_rknd ) :: &
      lambda, & ! a factor for defining how close you want to get to an equidistant grid
      threshold ! threshold which decides if grid should be adapted or not

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) ::  &
      gr_dens_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function [m]
      gr_dens         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function
                      ! [# levs/meter]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_levels), intent(out) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    logical, dimension(ngrdcol), intent(out) :: &
      l_adapt_grid ! whether the grid should be adapted or not

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx, i, j

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens_z,  &  ! the z coordinates of the connection points of the normalized
                         ! piecewise linear grid density function [m]
      min_gr_dens        ! the density at the given z coordinates of the connection points
                         ! of the normalized piecewise linear grid density function
                         ! [# levs/meter]

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized
                              ! piecewise linear grid density function [m]
      min_gr_dens_norm        ! the density at the given z coordinates of the connection points
                              ! of the normalized piecewise linear grid density function
                              ! [# levs/meter]

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized piecewise
                          ! linear grid density function [m]
      gr_dens_norm        ! the density at the given z coordinates of the connection points of the
                          ! normalized piecewise linear grid density function [# levs/meter]

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_z(1,1)
    grid_top = gr_dens_z(1,gr_dens_idx)

    call create_fixed_min_gr_dens_func( iunit+1, ngrdcol, &
                                        grid_sfc, grid_top )

    ! set the minimum grid density profile to be the linear piecewise function of the
    ! original minimum grid density function evaluated at the current grid levels
    ! this way the function can be executed with all checks and works out exactly
    do i = 1, ngrdcol
      do j = 1, gr_dens_idx
        min_gr_dens_z(i,j) = gr_dens_z(i,j)
      end do
      min_gr_dens(i,:) = zlinterp_fnc( gr_dens_idx, fixed_min_gr_dens_idx, &
                                       gr_dens_z(i,:), fixed_min_gr_dens_z(i,:), &
                                       fixed_min_gr_dens(i,:) )
    end do

    ! normalize the minimum grid density function
    call normalize_min_grid_density( ngrdcol, &                        ! Intent(in)
                                     gr_dens_idx, &                    ! Intent(in)
                                     min_gr_dens_z, min_gr_dens, &     ! Intent(in)
                                     lambda, num_levels, &             ! Intent(in)
                                     min_gr_dens_norm_z, &             ! Intent(out)
                                     min_gr_dens_norm )                ! Intent(out)

    ! normalize the grid density function with the normalized minimum grid density function
    call normalize_grid_density_helper( ngrdcol, &                               ! Intent(in)
                                        gr_dens_idx, &                           ! Intent(in)
                                        gr_dens_z, gr_dens, &                    ! Intent(in)
                                        min_gr_dens_norm, &                      ! Intent(in)
                                        lambda, num_levels, &                    ! Intent(in)
                                        gr_dens_norm_z, gr_dens_norm )           ! Intent(out)

    ! decide if grid should be adapted or not
    if ( .not. allocated( gr_dens_old_z_global ) .and. .not. allocated( gr_dens_old_global ) ) then
      error stop 'gr_dens_old_z_global and gr_dens_old_global were not allocated...'
    end if
    l_adapt_grid = decide_if_grid_adapt_helper( ngrdcol, &
                                                gr_dens_old_idx_global, &
                                                gr_dens_old_z_global, &
                                                gr_dens_old_global, &
                                                gr_dens_idx, &
                                                gr_dens_norm_z, &
                                                gr_dens_norm, &
                                                threshold )

    write(iunit, *) 'gr_dens_z', itime, gr_dens_norm_z
    write(iunit, *) 'gr_dens', itime, gr_dens_norm
    write(iunit, *) 'min_gr_dens_z', itime, min_gr_dens_norm_z
    write(iunit, *) 'min_gr_dens', itime, min_gr_dens_norm

    ! create grid from normalized grid density function
    do i = 1, ngrdcol
      if ( l_adapt_grid(i) ) then
        ! store new normalized grid density as old grid density for next adaptation iteration
        gr_dens_old_z_global(i,:) = gr_dens_norm_z(i,:)
        gr_dens_old_global(i,:) = gr_dens_norm(i,:)
        ! only create grid if it should actually be adapted
        grid_heights(i,:) = create_grid_from_normalized_grid_density_func( num_levels, &
                                                                           gr_dens_idx, &
                                                                           gr_dens_norm_z(i,:), &
                                                                           gr_dens_norm(i,:) )
      end if
    end do

    if ( clubb_at_least_debug_level( 2 ) ) then

        grid_heights_idx = size(grid_heights, 2)

        do i = 1, ngrdcol
          if ( l_adapt_grid(i) ) then
            ! only check if grid was actually adapted
            call check_grid( grid_heights_idx, grid_heights(i,:), &            ! Intent(in)
                             num_levels, &                                     ! Intent(in)
                             gr_dens_idx, &                                    ! Intent(in)
                             min_gr_dens_norm_z(i,:), min_gr_dens_norm(i,:), & ! Intent(in)
                             grid_sfc, grid_top )                              ! Intent(in)
          end if
        end do

    end if 

    return

  end subroutine create_grid_from_grid_density_func

  subroutine check_grid( grid_heights_idx, grid_heights, &
                         desired_num_levels, &
                         desired_min_dens_idx, &
                         desired_min_dens_z, desired_min_dens, &
                         desired_grid_sfc, desired_grid_top )

    ! Description:
    ! Checks if the grid is valid.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use interpolation, only: &
      zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      grid_heights_idx, &      ! number of elements in grid_heights []
      desired_min_dens_idx, &  ! number of elements in desired_min_dens_z and desired_min_dens []
      desired_num_levels       ! desired number of grid levels []
    ! Note: if everything worked fine, those two integers should be the same

    real( kind = core_rknd ), dimension(grid_heights_idx), intent(in) :: &
      grid_heights ! the grid heights ordered from bottom to top [m]

    real( kind = core_rknd ), dimension(desired_min_dens_idx), intent(in) :: &
      desired_min_dens_z, & ! the z coordinates of the piecewise linear desired minimum grid
                            ! density function [m]
      desired_min_dens      ! the density values on the given z coordinates [# levs/m]

    real( kind = core_rknd ), intent(in) :: &
      desired_grid_sfc, &   ! desired surface height of the grid [m]
      desired_grid_top      ! desired top height of the grid [m]

    !--------------------- Local Variables ---------------------
    integer :: j, k

    real( kind = core_rknd ) :: min_desired_min_dens

    real( kind = core_rknd ), dimension(desired_num_levels) :: desired_min_dens_on_new_gr

    !--------------------- Begin Code ---------------------
    ! check if first grid element is grid_sfc
    if (abs(grid_heights(1) - desired_grid_sfc) > tol) then
        write(fstderr,*) "WARNING! The first grid level is not correct. It should be the grids", &
                         " surface: ", desired_grid_sfc, " and not ", grid_heights(1)
        error stop 'Something went wrong with generating the grid'
    end if

    do j = 1, (grid_heights_idx-1)
        ! check if grid heights are getting bigger
        if (grid_heights(j) >= grid_heights(j+1)) then
            write(fstderr,*) "WARNING! The grid levels are not in the right order. ", &
                             "They should get bigger with the index."
            error stop 'Something went wrong with generating the grid'
        end if

        ! check if grid fulfills min_dens
        ! find minimal min_dens in grid cell
        ! check if minimum is in cell or on the boundaries
        ! first check values in interval
        min_desired_min_dens = maxval(desired_min_dens_z)
        do k = 1, desired_min_dens_idx
          if ( (desired_min_dens_z(k) > grid_heights(j)) &
               .and. (desired_min_dens_z(k) < grid_heights(j+1)) ) then
            if ( desired_min_dens(k) < min_desired_min_dens ) then
              min_desired_min_dens = desired_min_dens(k)
            end if
          end if
        end do

        ! then if maximum is on boundary
        desired_min_dens_on_new_gr = zlinterp_fnc( desired_num_levels, desired_min_dens_idx, &
                                                   grid_heights, desired_min_dens_z, &
                                                   desired_min_dens )
        if ( desired_min_dens_on_new_gr(j) < min_desired_min_dens ) then
          min_desired_min_dens = desired_min_dens_on_new_gr(j)
        end if
        if ( desired_min_dens_on_new_gr(j+1) < min_desired_min_dens ) then
          min_desired_min_dens = desired_min_dens_on_new_gr(j+1)
        end if

        if ((grid_heights(j+1) - grid_heights(j)) > one/(min_desired_min_dens) + tol) then
            write(fstderr,*) "WARNING! The grid has a distance=", &
                             (grid_heights(j+1) - grid_heights(j)), &
                             " which is bigger than the max_dist=", one/(min_desired_min_dens), &
                             "that is defined by min_dens."
            error stop 'Something went wrong with generating the grid'
        end if
    end do

    ! check if grid has the correct number of levels
    if (grid_heights_idx /= desired_num_levels) then
        write(fstderr,*) "WARNING! The grid has not the right number of levels. There are ", &
                         grid_heights_idx, " level, but there should be ", desired_num_levels, &
                         " level."
        error stop 'Something went wrong with generating the grid'
    end if

    ! check if last grid level is grid_top
    if (abs(grid_heights(desired_num_levels) - desired_grid_top) > tol) then
        write(fstderr,*) "WARNING! The last grid level is not correct. It should be the grids", &
                         " top: ", desired_grid_top, " and not ", grid_heights(desired_num_levels)
        error stop 'Something went wrong with generating the grid'
    end if

  end subroutine check_grid

  subroutine adapt_grid( iunit, itime, ngrdcol, gr_dens_norm_idx, &
                         gr_dens_norm_z, gr_dens_norm, &
                         min_gr_dens_norm_z, min_gr_dens_norm, &
                         sfc_elevation, l_implemented, &
                         hydromet_dim, sclr_dim, edsclr_dim, &
                         thvm, idx_thvm, p_sfc, &
                         grid_remap_method, &
                         gr, &
                         thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
                         sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
                         wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
                         rtpthlp_forcing, wm_zm, wm_zt, &
                         rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
                         p_in_Pa, rho_zm, rho, exner, &
                         rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                         invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                         rfrzm, wphydrometp, &
                         wp2hmp, rtphmp_zt, thlphmp_zt, &
                         um, vm, upwp, vpwp, up2, vp2, up3, vp3, &
                         thlm, rtm, wprtp, wpthlp, &
                         wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
                         sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &
                         wpsclrp, edsclrm, &
                         rcm, cloud_frac, &
                         wpthvp, wp2thvp, rtpthvp, thlpthvp, &
                         sclrpthvp, &
                         wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &
                         wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &
                         um_pert, vm_pert, upwp_pert, vpwp_pert, &
                         Kh_zm, Kh_zt, &
                         thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &
                         cloudy_updraft_frac, cloudy_downdraft_frac, &
                         rcm_in_layer, cloud_cover, invrs_tau_zm, &
                         Lscale )
    ! Description:
    ! Adapts the grid based on the density function and interpolates all values to the new grid.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)
    
    use grid_class, only: &
      setup_grid

    use calc_pressure, only: &
      init_pressure    ! Procedure(s)

    use model_flags, only: &
      cons_ullrich_remap

    use error_code, only: &
      clubb_fatal_error

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      iunit, &
      itime, &
      ngrdcol, &
      hydromet_dim, &
      sclr_dim, &
      edsclr_dim, &
      gr_dens_norm_idx, & ! total numer of indices of density_func_z and density_func_dens []
      idx_thvm   ! numer of indices of thvm

    real( kind = core_rknd ), dimension(gr_dens_norm_idx), intent(in) ::  &
      min_gr_dens_norm_z, min_gr_dens_norm, &
      gr_dens_norm_z,  &   ! the height values of the connection points of the piecewise linear
                           ! grid density function/profile [m]
      gr_dens_norm    ! the density values of the connection points of the piecewise linear
                           ! grid density function/profile [# levs/m]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      sfc_elevation, p_sfc

    logical, intent(in) :: l_implemented

    real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
      thvm

    integer, intent(in) :: &
      grid_remap_method ! specifies what remapping method should be used

    !--------------------- Output Variable ---------------------

    !--------------------- In/Out Variable ---------------------
    type( grid ), intent(inout) :: gr

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm, hydromet_dim), intent(inout) :: &
      wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt, hydromet_dim), intent(inout) :: &
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and a hydrometeor [K <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) :: &
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
        edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Lscale     ! Length scale         [m]

    !--------------------- Local Variables ---------------------
    integer :: &
        num_levels, &                ! number of levels the new grid should have []
        total_idx_rho_lin_spline, &  ! total number of indices of the rho density
                                     ! piecewise linear function []
        grid_type, &
        err_code, &
        i, k

    real( kind = core_rknd ) ::  &
      equi_dens   ! density of the equidistant grid [# levs/meter]

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      deltaz

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) ::  &
      new_gr_zm      ! zm levels for the adapted grid based on the density function [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) ::  &
      thermodynamic_heights_placeholder  ! placeholder for the zt levels, but not needed

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      rho_lin_spline_vals, &
      rho_lin_spline_levels
    
    type( grid ) :: new_gr

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzt) ::  &
      thvm_new_grid

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      p_in_Pa_zm, exner_zm ! just placeholder variables; are only needed for subroutine call

    real( kind = core_rknd ) ::  &
      lambda, & ! a factor between 0 and 1 defining how close you want to get to an equidistant grid
      threshold ! some threshold to decide wether grid should be adapted or not

    real( kind = core_rknd ), dimension(gr%nzm,gr%nzm) :: &
      R_ij_zm

    real( kind = core_rknd ), dimension(gr%nzt,gr%nzt) :: &
      R_ij_zt

    real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_source_zm_vals

    real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_target_zm_vals

    logical, dimension(ngrdcol) :: l_adapt_grid

    real( kind = core_rknd ), dimension(1,10) ::  &
      input_test, output_test

    !--------------------- Begin Code ---------------------
    ! Setup the new grid
    num_levels = gr%nzm
    equi_dens = (num_levels-1)/(gr%zm(1,gr%nzm) - gr%zm(1,1))
    ! lalala
    !threshold = 0.0
    !threshold = 4.0e-6
    
    !threshold = 5.0e-6   !!
    !threshold = 6.0e-6
    !threshold = 7.0e-6
    
    !threshold = 8.0e-6   !!!
    
    !threshold = 9.0e-6
    !threshold = 1.0e-5
    !threshold = 2.0e-5

    !lambda = 0.00000000000000000001
    !lambda = 0.3

    !threshold = 1.0e-3
    threshold = 5.0e-4 !!!!!
    threshold = 4.0e-4
    threshold = 6.0e-4 !!!!
    threshold = 7.0e-4 !!!!!!!!!!
    threshold = 6.5e-4 !!!!!!!!!!!!!!
    
    ! trigger threshold
    threshold = 5.0e-3   !!!!!!!!!!!!!!!!!!!!!
    threshold = 6.0e-3
    !threshold = 7.0e-3

    ! Allocate and set gr_dens_old_global if not already allocated
    if ( .not. allocated( gr_dens_old_global ) ) then
      gr_dens_old_idx_global = gr%nzm
      allocate( gr_dens_old_global(ngrdcol,gr%nzm) )
      gr_dens_old_global = gr%invrs_dzm
    end if

    if ( .not. allocated( gr_dens_old_z_global ) ) then
      allocate( gr_dens_old_z_global(ngrdcol,gr%nzm) )
      gr_dens_old_z_global = gr%zm
    end if

    !call create_grid_from_grid_density_func( ngrdcol, &
    !                                         iunit, itime, &
    !                                         total_idx_density_func, &
    !                                         density_func_z, density_func_dens, &
    !                                         lambda, &
    !                                         num_levels, &
    !                                         threshold, &
    !                                         new_gr_zm, &
    !                                         l_adapt_grid )

    call create_grid_from_normalized_grid_density( ngrdcol, &
                                                   gr_dens_norm_idx, &
                                                   gr_dens_norm_z, gr_dens_norm, &
                                                   min_gr_dens_norm_z, min_gr_dens_norm, &
                                                   num_levels, &
                                                   threshold, &
                                                   new_gr_zm, &
                                                   l_adapt_grid )
    

    do i = 1, ngrdcol
      deltaz(i) = 0.0_core_rknd
    end do
    grid_type = 3

    ! TODOlala adjust for ngrdcol > 1
    if ( l_adapt_grid(1) ) then
      call setup_grid( num_levels, ngrdcol, sfc_elevation, l_implemented, &   ! intent(in)
                       grid_type, deltaz, gr%zm(:,1), gr%zm(:,num_levels), &  ! intent(in)
                       new_gr_zm, thermodynamic_heights_placeholder, &        ! intent(in)
                       new_gr, err_code )                                     ! intent(inout)

      if ( err_code == clubb_fatal_error ) then
        error stop "Error in CLUBB calling setup_grid"
      end if

      Lscale_counter = 0
      brunt_counter = 0
      chi_counter = 0
      do k = 1, gr%nzm
        cumulative_Lscale(k) = 0.0_core_rknd
        cumulative_brunt(k) = 0.0_core_rknd
        cumulative_chi(k) = 0.0_core_rknd
      end do

      ! Set the density values to use for interpolation for mass calculation
      total_idx_rho_lin_spline = gr%nzm
      rho_lin_spline_vals = rho_ds_zm
      rho_lin_spline_levels = gr%zm

      call remap_all_clubb_core_vals( ngrdcol, total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, rho_lin_spline_levels, &
                                      l_implemented, &
                                      hydromet_dim, sclr_dim, edsclr_dim, &
                                      thvm, idx_thvm, p_sfc, &
                                      grid_remap_method, &
                                      gr, new_gr, &
                                      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
                                      sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
                                      wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
                                      rtpthlp_forcing, wm_zm, wm_zt, &
                                      rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
                                      p_in_Pa, rho_zm, rho, exner, &
                                      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                                      invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                                      rfrzm, wphydrometp, &
                                      wp2hmp, rtphmp_zt, thlphmp_zt, &
                                      um, vm, upwp, vpwp, up2, vp2, up3, vp3, &
                                      thlm, rtm, wprtp, wpthlp, &
                                      wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
                                      sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &
                                      wpsclrp, edsclrm, &
                                      rcm, cloud_frac, &
                                      wpthvp, wp2thvp, rtpthvp, thlpthvp, &
                                      sclrpthvp, &
                                      wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &
                                      wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &
                                      um_pert, vm_pert, upwp_pert, vpwp_pert, &
                                      Kh_zm, Kh_zt, &
                                      thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &
                                      cloudy_updraft_frac, cloudy_downdraft_frac, &
                                      rcm_in_layer, cloud_cover, invrs_tau_zm, &
                                      Lscale )

      
      gr = new_gr
    end if

  end subroutine adapt_grid

  !subroutine adapt_grid_old( iunit, itime, ngrdcol, total_idx_density_func, &
  !                       density_func_z, density_func_dens, &
  !                       sfc_elevation, l_implemented, &
  !                       hydromet_dim, sclr_dim, edsclr_dim, &
  !                       thvm, idx_thvm, p_sfc, &
  !                       grid_remap_method, &
  !                       gr, &
  !                       thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
  !                       sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
  !                       wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
  !                       rtpthlp_forcing, wm_zm, wm_zt, &
  !                       rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
  !                       p_in_Pa, rho_zm, rho, exner, &
  !                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
  !                       invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
  !                       rfrzm, wphydrometp, &
  !                       wp2hmp, rtphmp_zt, thlphmp_zt, &
  !                       um, vm, upwp, vpwp, up2, vp2, up3, vp3, &
  !                       thlm, rtm, wprtp, wpthlp, &
  !                       wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
  !                       sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &
  !                       wpsclrp, edsclrm, &
  !                       rcm, cloud_frac, &
  !                       wpthvp, wp2thvp, rtpthvp, thlpthvp, &
  !                       sclrpthvp, &
  !                       wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &
  !                       wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &
  !                       um_pert, vm_pert, upwp_pert, vpwp_pert, &
  !                       Kh_zm, Kh_zt, &
  !                       thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &
  !                       cloudy_updraft_frac, cloudy_downdraft_frac, &
  !                       rcm_in_layer, cloud_cover, invrs_tau_zm, &
  !                       Lscale )
  !  ! Description:
  !  ! Adapts the grid based on the density function and interpolates all values to the new grid.
  !  !-----------------------------------------------------------------------
!
  !  use clubb_precision, only: &
  !    core_rknd ! Variable(s)
  !  
  !  use grid_class, only: &
  !    setup_grid
!
  !  use calc_pressure, only: &
  !    init_pressure    ! Procedure(s)
!
  !  use model_flags, only: &
  !    cons_ullrich_remap
!
  !  use error_code, only: &
  !    clubb_fatal_error
!
  !  implicit none
!
  !  !--------------------- Input Variables ---------------------
  !  integer, intent(in) :: &
  !    iunit, &
  !    itime, &
  !    ngrdcol, &
  !    hydromet_dim, &
  !    sclr_dim, &
  !    edsclr_dim, &
  !    total_idx_density_func, & ! total numer of indices of density_func_z and density_func_dens []
  !    idx_thvm   ! numer of indices of thvm
!
  !  real( kind = core_rknd ), dimension(total_idx_density_func), intent(in) ::  &
  !    density_func_z,  &   ! the height values of the connection points of the piecewise linear
  !                         ! grid density function/profile [m]
  !    density_func_dens    ! the density values of the connection points of the piecewise linear
  !                         ! grid density function/profile [# levs/m]
!
  !  real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
  !    sfc_elevation, p_sfc
!
  !  logical, intent(in) :: l_implemented
!
  !  real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
  !    thvm
!
  !  integer, intent(in) :: &
  !    grid_remap_method ! specifies what remapping method should be used
!
  !  !--------------------- Output Variable ---------------------
!
  !  !--------------------- In/Out Variable ---------------------
  !  type( grid ), intent(inout) :: gr
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
  !    rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
  !    um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
  !    vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
  !    wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
  !    rho,             & ! Air density on thermodynamic levels       [kg/m^3]
  !    rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
  !    invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
  !    thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
  !    rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
  !    wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
  !    wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
  !    rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
  !    thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
  !    rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
  !    wm_zm,           & ! w mean wind component on momentum levels  [m/s]
  !    rho_zm,          & ! Air density on momentum levels            [kg/m^3]
  !    rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
  !    invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
  !    thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzm, hydromet_dim), intent(inout) :: &
  !    wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzt, hydromet_dim), intent(inout) :: &
  !    wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
  !    rtphmp_zt,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
  !    thlphmp_zt     ! Covariance of thl and a hydrometeor [K <hm units>]
!
  !  ! Passive scalar variables
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
  !    sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]
!
  !  ! Eddy passive scalar variables
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
  !    edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]
!
  !  ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) ::  &
  !    rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
  !    thlm_ref, & ! Initial liquid water potential temperature   [K]
  !    um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
  !    vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
  !    ug,       & ! u geostrophic wind                           [m/s]
  !    vg          ! v geostrophic wind                           [m/s]
!
  !  ! These are prognostic or are planned to be in the future
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    um,      & ! u mean wind component (thermodynamic levels)   [m/s]
  !    vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
  !    up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
  !    vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
  !    rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
  !    thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
  !    rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
  !    thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
  !    wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
  !    upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
  !    vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
  !    up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
  !    vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
  !    wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
  !    wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
  !    rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
  !    thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
  !    rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
  !    wp2        ! w'^2 (momentum levels)                         [m^2/s^2]
!
  !  ! Passive scalar variables
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
  !    sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
  !    sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
  !    wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
  !    sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
  !    sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
  !    sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
  !    exner      ! Exner function (thermodynamic levels)     [-]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
  !    cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
  !    wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
  !    wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
  !    rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
  !    thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
  !    sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
  !    wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
  !    wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
  !    wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
  !    ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
  !    uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
  !    vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
  !    rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
  !    wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
  !    wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
  !    wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]
!
  !  ! Variables used to track perturbed version of winds.
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) :: &
  !    um_pert,   & ! perturbed <u>       [m/s]
  !    vm_pert      ! perturbed <v>       [m/s]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) :: &
  !    upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
  !    vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
  !      edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
  !    cloud_cover     ! cloud cover                                     [-]
!
  !  ! Variables that need to be output for use in host models
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
  !    wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
  !    invrs_tau_zm             ! One divided by tau on zm levels       [1/s]
!
  !  real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
  !    w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
  !    w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
  !    cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
  !    cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
  !    Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
  !    Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
  !    thlprcp    ! thl'rc'              [K kg/kg]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
  !    Lscale     ! Length scale         [m]
!
  !  !--------------------- Local Variables ---------------------
  !  integer :: &
  !      num_levels, &                ! number of levels the new grid should have []
  !      total_idx_rho_lin_spline, &  ! total number of indices of the rho density
  !                                   ! piecewise linear function []
  !      grid_type, &
  !      err_code, &
  !      i, k
!
  !  real( kind = core_rknd ) ::  &
  !    equi_dens   ! density of the equidistant grid [# levs/meter]
!
  !  real( kind = core_rknd ), dimension(ngrdcol) ::  &
  !    deltaz
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) ::  &
  !    new_gr_zm      ! zm levels for the adapted grid based on the density function [m]
!
  !  real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) ::  &
  !    thermodynamic_heights_placeholder  ! placeholder for the zt levels, but not needed
!
  !  real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
  !    rho_lin_spline_vals, &
  !    rho_lin_spline_levels
  !  
  !  type( grid ) :: new_gr
!
  !  real( kind = core_rknd ), dimension(ngrdcol, gr%nzt) ::  &
  !    thvm_new_grid
!
  !  real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
  !    p_in_Pa_zm, exner_zm ! just placeholder variables; are only needed for subroutine call
!
  !  real( kind = core_rknd ) ::  &
  !    lambda, & ! a factor between 0 and 1 defining how close you want to get to an equidistant grid
  !    threshold ! some threshold to decide wether grid should be adapted or not
!
  !  real( kind = core_rknd ), dimension(gr%nzm,gr%nzm) :: &
  !    R_ij_zm
!
  !  real( kind = core_rknd ), dimension(gr%nzt,gr%nzt) :: &
  !    R_ij_zt
!
  !  real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_source_zm_vals
!
  !  real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_target_zm_vals
!
  !  logical, dimension(ngrdcol) :: l_adapt_grid
!
  !  real( kind = core_rknd ), dimension(1,10) ::  &
  !    input_test, output_test
!
  !  integer :: &
  !    iv_wind = -1, &
  !    iv_pos_def = 0, &
  !    iv_other = 1
!
  !  !--------------------- Begin Code ---------------------
  !  ! Setup the new grid
  !  num_levels = gr%nzm
  !  equi_dens = (num_levels-1)/(gr%zm(1,gr%nzm) - gr%zm(1,1))
  !  ! lalala
  !  !threshold = 3.0e-6
  !  threshold = 4.0e-6
  !  !threshold = 5.0e-6
  !  lambda = 0.00000000000000000001
  !  lambda = 0.3
!
  !  ! Allocate and set gr_dens_old_global if not already allocated
  !  if ( .not. allocated( gr_dens_old_global ) ) then
  !    gr_dens_old_idx_global = gr%nzm
  !    allocate( gr_dens_old_global(ngrdcol,gr%nzm) )
  !    gr_dens_old_global = gr%invrs_dzm
  !  end if
!
  !  if ( .not. allocated( gr_dens_old_z_global ) ) then
  !    allocate( gr_dens_old_z_global(ngrdcol,gr%nzm) )
  !    gr_dens_old_z_global = gr%zm
  !  end if
!
  !  call create_grid_from_grid_density_func( ngrdcol, &
  !                                           iunit, itime, &
  !                                           total_idx_density_func, &
  !                                           density_func_z, density_func_dens, &
  !                                           lambda, &
  !                                           num_levels, &
  !                                           threshold, &
  !                                           new_gr_zm, &
  !                                           l_adapt_grid )
  !  
!
  !  do i = 1, ngrdcol
  !    deltaz(i) = 0.0_core_rknd
  !  end do
  !  grid_type = 3
!
  !  ! TODOlala adjust for ngrdcol > 1
  !  if ( l_adapt_grid(1) ) then
  !    call setup_grid( num_levels, ngrdcol, sfc_elevation, l_implemented, &   ! intent(in)
  !                     grid_type, deltaz, gr%zm(:,1), gr%zm(:,num_levels), &  ! intent(in)
  !                     new_gr_zm, thermodynamic_heights_placeholder, &        ! intent(in)
  !                     new_gr, err_code )                                     ! intent(inout)
!
  !    if ( err_code == clubb_fatal_error ) then
  !      error stop "Error in CLUBB calling setup_grid"
  !    end if
!
  !    Lscale_counter = 0
  !    do k = 1, gr%nzm
  !      cumulative_Lscale(k) = 0.0_core_rknd
  !    end do
!
  !    ! Set the density values to use for interpolation for mass calculation
  !    total_idx_rho_lin_spline = gr%nzm
  !    rho_lin_spline_vals = rho_ds_zm
  !    rho_lin_spline_levels = gr%zm
!
  !    !write(*,*) 'input ppm: ', thlm_forcing
!
!
!
  !    !call map1_ppm( num_levels,   gr%zm(1,:),   thlm_forcing,  num_levels,   new_gr%zm(1,:),   output_test,                &
  !    !                   0, 0, 1, 1, 1,                         &
  !    !                   1, 1, 1, 0, 3)
!
  !    !call map1_ppm( 10,   \[0, 1, 3, 4, 5, 8, 10, 20, 30, 50\],   \[0.5, 0.8, 0.9, 1.2, 1.5, 1.7, 2.0, 1.0, 0.8, 1.0\],  10,   \[0, 0.5, 1, 1.5, 3, 5, 8, 17, 25, 50\],   output_test,                &
  !    !                   0, 0, 1, 1, 1,                         &
  !    !                   1, 1, 1, 0, 3)
!
  !    !write(*,*) 'output ppm: ', output_test
!
  !    call remapping_matrix_zm_values( ngrdcol, &                                      ! Intent(in)
  !                                     gr, new_gr, &                                   ! Intent(in)
  !                                     total_idx_rho_lin_spline, &                     ! Intent(in)
  !                                     rho_lin_spline_vals, &                          ! Intent(in)
  !                                     rho_lin_spline_levels, &                        ! Intent(in)
  !                                     levels_source_zm_vals, levels_target_zm_vals, & ! Intent(out)
  !                                     R_ij_zm )                                       ! Intent(out)
!
  !    call remapping_matrix_zt_values( ngrdcol, &                     ! Intent(in)
  !                                     gr, new_gr, &                  ! Intent(in)
  !                                     total_idx_rho_lin_spline, &    ! Intent(in)
  !                                     rho_lin_spline_vals, &         ! Intent(in)
  !                                     rho_lin_spline_levels, &       ! Intent(in)
  !                                     R_ij_zt )                      ! Intent(out)
!
  !    if ( grid_remap_method == cons_ullrich_remap ) then
!
  !      ! Remap all zt values
  !      thlm_forcing = remap_vals_to_target( ngrdcol, &
  !                                           gr%nzm, new_gr%nzm, &
  !                                           gr%zm, new_gr%zm, &
  !                                           total_idx_rho_lin_spline, &
  !                                           rho_lin_spline_vals, &
  !                                           rho_lin_spline_levels, &
  !                                           thlm_forcing, &
  !                                           iv_other, R_ij_zt, p_sfc )
  !      rtm_forcing = remap_vals_to_target( ngrdcol, &
  !                                          gr%nzm, new_gr%nzm, &
  !                                          gr%zm, new_gr%zm, &
  !                                          total_idx_rho_lin_spline, &
  !                                          rho_lin_spline_vals, &
  !                                          rho_lin_spline_levels, &
  !                                          rtm_forcing, &
  !                                          iv_other, R_ij_zt, p_sfc )
  !      um_forcing = remap_vals_to_target( ngrdcol, &
  !                                         gr%nzm, new_gr%nzm, &
  !                                         gr%zm, new_gr%zm, &
  !                                         total_idx_rho_lin_spline, &
  !                                         rho_lin_spline_vals, &
  !                                         rho_lin_spline_levels, &
  !                                         um_forcing, &
  !                                         iv_other, R_ij_zt, p_sfc )
  !      vm_forcing = remap_vals_to_target( ngrdcol, &
  !                                         gr%nzm, new_gr%nzm, &
  !                                         gr%zm, new_gr%zm, &
  !                                         total_idx_rho_lin_spline, &
  !                                         rho_lin_spline_vals, &
  !                                         rho_lin_spline_levels, &
  !                                         vm_forcing, &
  !                                         iv_other, R_ij_zt, p_sfc )
  !      wm_zt = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    wm_zt, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      rho = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  rho, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      rho_ds_zt = remap_vals_to_target( ngrdcol, &
  !                                        gr%nzm, new_gr%nzm, &
  !                                        gr%zm, new_gr%zm, &
  !                                        total_idx_rho_lin_spline, &
  !                                        rho_lin_spline_vals, &
  !                                        rho_lin_spline_levels, &
  !                                        rho_ds_zt, &
  !                                        iv_other, R_ij_zt, p_sfc )
  !      invrs_rho_ds_zt = remap_vals_to_target( ngrdcol, &
  !                                              gr%nzm, new_gr%nzm, &
  !                                              gr%zm, new_gr%zm, &
  !                                              total_idx_rho_lin_spline, &
  !                                              rho_lin_spline_vals, &
  !                                              rho_lin_spline_levels, &
  !                                              invrs_rho_ds_zt, &
  !                                              iv_other, R_ij_zt, p_sfc )
  !      thv_ds_zt = remap_vals_to_target( ngrdcol, &
  !                                        gr%nzm, new_gr%nzm, &
  !                                        gr%zm, new_gr%zm, &
  !                                        total_idx_rho_lin_spline, &
  !                                        rho_lin_spline_vals, &
  !                                        rho_lin_spline_levels, &
  !                                        thv_ds_zt, &
  !                                        iv_other, R_ij_zt, p_sfc )
  !      rfrzm = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    rfrzm, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      do i = 1, hydromet_dim
  !          wp2hmp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                gr%nzm, new_gr%nzm, &
  !                                                gr%zm, new_gr%zm, &
  !                                                total_idx_rho_lin_spline, &
  !                                                rho_lin_spline_vals, &
  !                                                rho_lin_spline_levels, &
  !                                                wp2hmp(:,:,i), &
  !                                                iv_other, R_ij_zt, p_sfc )
  !          rtphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                   gr%nzm, new_gr%nzm, &
  !                                                   gr%zm, new_gr%zm, &
  !                                                   total_idx_rho_lin_spline, &
  !                                                   rho_lin_spline_vals, &
  !                                                   rho_lin_spline_levels, &
  !                                                   rtphmp_zt(:,:,i), &
  !                                                   iv_other, R_ij_zt, p_sfc )
  !          thlphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                    gr%nzm, new_gr%nzm, &
  !                                                    gr%zm, new_gr%zm, &
  !                                                    total_idx_rho_lin_spline, &
  !                                                    rho_lin_spline_vals, &
  !                                                    rho_lin_spline_levels, &
  !                                                    thlphmp_zt(:,:,i), &
  !                                                    iv_other, R_ij_zt, p_sfc )
  !      end do
  !      do i = 1, sclr_dim
  !          sclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                       gr%nzm, new_gr%nzm, &
  !                                                       gr%zm, new_gr%zm, &
  !                                                       total_idx_rho_lin_spline, &
  !                                                       rho_lin_spline_vals, &
  !                                                       rho_lin_spline_levels, &
  !                                                       sclrm_forcing(:,:,i), &
  !                                                       iv_other, R_ij_zt, p_sfc )
  !          sclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                               gr%nzm, new_gr%nzm, &
  !                                               gr%zm, new_gr%zm, &
  !                                               total_idx_rho_lin_spline, &
  !                                               rho_lin_spline_vals, &
  !                                               rho_lin_spline_levels, &
  !                                               sclrm(:,:,i), &
  !                                               iv_other, R_ij_zt, p_sfc )
  !          sclrp3(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                gr%nzm, new_gr%nzm, &
  !                                                gr%zm, new_gr%zm, &
  !                                                total_idx_rho_lin_spline, &
  !                                                rho_lin_spline_vals, &
  !                                                rho_lin_spline_levels, &
  !                                                sclrp3(:,:,i), &
  !                                                iv_other, R_ij_zt, p_sfc )
  !      end do
  !      do i = 1, edsclr_dim
  !          edsclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                         gr%nzm, new_gr%nzm, &
  !                                                         gr%zm, new_gr%zm, &
  !                                                         total_idx_rho_lin_spline, &
  !                                                         rho_lin_spline_vals, &
  !                                                         rho_lin_spline_levels, &
  !                                                         edsclrm_forcing(:,:,i), &
  !                                                         iv_other, R_ij_zt, p_sfc )
  !          edsclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                 gr%nzm, new_gr%nzm, &
  !                                                 gr%zm, new_gr%zm, &
  !                                                 total_idx_rho_lin_spline, &
  !                                                 rho_lin_spline_vals, &
  !                                                 rho_lin_spline_levels, &
  !                                                 edsclrm(:,:,i), &
  !                                                 iv_other, R_ij_zt, p_sfc )
  !      end do
  !      rtm_ref = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzm, new_gr%nzm, &
  !                                      gr%zm, new_gr%zm, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      rtm_ref, &
  !                                      iv_other, R_ij_zt, p_sfc )
  !      thlm_ref = remap_vals_to_target( ngrdcol, &
  !                                       gr%nzm, new_gr%nzm, &
  !                                       gr%zm, new_gr%zm, &
  !                                       total_idx_rho_lin_spline, &
  !                                       rho_lin_spline_vals, &
  !                                       rho_lin_spline_levels, &
  !                                       thlm_ref, &
  !                                       iv_other, R_ij_zt, p_sfc )
  !      um_ref = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzm, new_gr%nzm, &
  !                                     gr%zm, new_gr%zm, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     um_ref, &
  !                                     iv_other, R_ij_zt, p_sfc )
  !      vm_ref = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzm, new_gr%nzm, &
  !                                     gr%zm, new_gr%zm, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     vm_ref, &
  !                                     iv_other, R_ij_zt, p_sfc )
  !      ug = remap_vals_to_target( ngrdcol, &
  !                                 gr%nzm, new_gr%nzm, &
  !                                 gr%zm, new_gr%zm, &
  !                                 total_idx_rho_lin_spline, &
  !                                 rho_lin_spline_vals, &
  !                                 rho_lin_spline_levels, &
  !                                 ug, &
  !                                 iv_other, R_ij_zt, p_sfc )
  !      vg = remap_vals_to_target( ngrdcol, &
  !                                 gr%nzm, new_gr%nzm, &
  !                                 gr%zm, new_gr%zm, &
  !                                 total_idx_rho_lin_spline, &
  !                                 rho_lin_spline_vals, &
  !                                 rho_lin_spline_levels, &
  !                                 vg, &
  !                                 iv_other, R_ij_zt, p_sfc )
  !      um = remap_vals_to_target( ngrdcol, &
  !                                 gr%nzm, new_gr%nzm, &
  !                                 gr%zm, new_gr%zm, &
  !                                 total_idx_rho_lin_spline, &
  !                                 rho_lin_spline_vals, &
  !                                 rho_lin_spline_levels, &
  !                                 um, &
  !                                 iv_other, R_ij_zt, p_sfc )
  !      vm = remap_vals_to_target( ngrdcol, &
  !                                 gr%nzm, new_gr%nzm, &
  !                                 gr%zm, new_gr%zm, &
  !                                 total_idx_rho_lin_spline, &
  !                                 rho_lin_spline_vals, &
  !                                 rho_lin_spline_levels, &
  !                                 vm, &
  !                                 iv_other, R_ij_zt, p_sfc )
  !      up3 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  up3, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      vp3 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  vp3, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      rtm = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  rtm, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      thlm = remap_vals_to_target( ngrdcol, &
  !                                   gr%nzm, new_gr%nzm, &
  !                                   gr%zm, new_gr%zm, &
  !                                   total_idx_rho_lin_spline, &
  !                                   rho_lin_spline_vals, &
  !                                   rho_lin_spline_levels, &
  !                                   thlm, &
  !                                   iv_other, R_ij_zt, p_sfc )
  !      rtp3 = remap_vals_to_target( ngrdcol, &
  !                                   gr%nzm, new_gr%nzm, &
  !                                   gr%zm, new_gr%zm, &
  !                                   total_idx_rho_lin_spline, &
  !                                   rho_lin_spline_vals, &
  !                                   rho_lin_spline_levels, &
  !                                   rtp3, &
  !                                   iv_other, R_ij_zt, p_sfc )
  !      thlp3 = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    thlp3, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      wp3 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  wp3, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      ! remap thvm values to new grid, to calculate new p_in_Pa
  !      thvm_new_grid = remap_vals_to_target( ngrdcol, &
  !                                            gr%nzm, new_gr%nzm, &
  !                                            gr%zm, new_gr%zm, &
  !                                            total_idx_rho_lin_spline, &
  !                                            rho_lin_spline_vals, &
  !                                            rho_lin_spline_levels, &
  !                                            thvm, &
  !                                            iv_other, R_ij_zt, p_sfc )
  !      ! calculate p_in_Pa instead of remapping directly since it can run into problems if for
  !      ! example the two highest levels have the same pressure value, which could be happening with
  !      ! the remapping, also calculate exner accordingly
  !      call init_pressure( ngrdcol, new_gr, thvm_new_grid, p_sfc, &
  !                          p_in_Pa, exner, p_in_Pa_zm, exner_zm )
  !      rcm = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzm, new_gr%nzm, &
  !                                  gr%zm, new_gr%zm, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  rcm, &
  !                                  iv_other, R_ij_zt, p_sfc )
  !      cloud_frac = remap_vals_to_target( ngrdcol, &
  !                                         gr%nzm, new_gr%nzm, &
  !                                         gr%zm, new_gr%zm, &
  !                                         total_idx_rho_lin_spline, &
  !                                         rho_lin_spline_vals, &
  !                                         rho_lin_spline_levels, &
  !                                         cloud_frac, &
  !                                         iv_other, R_ij_zt, p_sfc )
  !      wp2thvp = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzm, new_gr%nzm, &
  !                                      gr%zm, new_gr%zm, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      wp2thvp, &
  !                                      iv_other, R_ij_zt, p_sfc )
  !      wp2rtp = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzm, new_gr%nzm, &
  !                                     gr%zm, new_gr%zm, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     wp2rtp, &
  !                                     iv_other, R_ij_zt, p_sfc )
  !      wp2thlp = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzm, new_gr%nzm, &
  !                                      gr%zm, new_gr%zm, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      wp2thlp, &
  !                                      iv_other, R_ij_zt, p_sfc )
  !      wpup2 = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    wpup2, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      wpvp2 = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    wpvp2, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      ice_supersat_frac = remap_vals_to_target( ngrdcol, &
  !                                                gr%nzm, new_gr%nzm, &
  !                                                gr%zm, new_gr%zm, &
  !                                                total_idx_rho_lin_spline, &
  !                                                rho_lin_spline_vals, &
  !                                                rho_lin_spline_levels, &
  !                                                ice_supersat_frac, &
  !                                                iv_other, R_ij_zt, p_sfc )
  !      um_pert = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzm, new_gr%nzm, &
  !                                      gr%zm, new_gr%zm, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      um_pert, &
  !                                      iv_other, R_ij_zt, p_sfc )
  !      vm_pert = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzm, new_gr%nzm, &
  !                                      gr%zm, new_gr%zm, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      vm_pert, &
  !                                      iv_other, R_ij_zt, p_sfc )
  !      rcm_in_layer = remap_vals_to_target( ngrdcol, &
  !                                           gr%nzm, new_gr%nzm, &
  !                                           gr%zm, new_gr%zm, &
  !                                           total_idx_rho_lin_spline, &
  !                                           rho_lin_spline_vals, &
  !                                           rho_lin_spline_levels, &
  !                                           rcm_in_layer, &
  !                                           iv_other, R_ij_zt, p_sfc )
  !      cloud_cover = remap_vals_to_target( ngrdcol, &
  !                                          gr%nzm, new_gr%nzm, &
  !                                          gr%zm, new_gr%zm, &
  !                                          total_idx_rho_lin_spline, &
  !                                          rho_lin_spline_vals, &
  !                                          rho_lin_spline_levels, &
  !                                          cloud_cover, &
  !                                          iv_other, R_ij_zt, p_sfc )
  !      w_up_in_cloud = remap_vals_to_target( ngrdcol, &
  !                                            gr%nzm, new_gr%nzm, &
  !                                            gr%zm, new_gr%zm, &
  !                                            total_idx_rho_lin_spline, &
  !                                            rho_lin_spline_vals, &
  !                                            rho_lin_spline_levels, &
  !                                            w_up_in_cloud, &
  !                                            iv_other, R_ij_zt, p_sfc )
  !      w_down_in_cloud = remap_vals_to_target( ngrdcol, &
  !                                              gr%nzm, new_gr%nzm, &
  !                                              gr%zm, new_gr%zm, &
  !                                              total_idx_rho_lin_spline, &
  !                                              rho_lin_spline_vals, &
  !                                              rho_lin_spline_levels, &
  !                                              w_down_in_cloud, &
  !                                              iv_other, R_ij_zt, p_sfc )
  !      cloudy_updraft_frac = remap_vals_to_target( ngrdcol, &
  !                                                  gr%nzm, new_gr%nzm, &
  !                                                  gr%zm, new_gr%zm, &
  !                                                  total_idx_rho_lin_spline, &
  !                                                  rho_lin_spline_vals, &
  !                                                  rho_lin_spline_levels, &
  !                                                  cloudy_updraft_frac, &
  !                                                  iv_other, R_ij_zt, p_sfc )
  !      cloudy_downdraft_frac = remap_vals_to_target( ngrdcol, &
  !                                                    gr%nzm, new_gr%nzm, &
  !                                                    gr%zm, new_gr%zm, &
  !                                                    total_idx_rho_lin_spline, &
  !                                                    rho_lin_spline_vals, &
  !                                                    rho_lin_spline_levels, &
  !                                                    cloudy_downdraft_frac, &
  !                                                    iv_other, R_ij_zt, p_sfc )
  !      Kh_zt = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzm, new_gr%nzm, &
  !                                    gr%zm, new_gr%zm, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    Kh_zt, &
  !                                    iv_other, R_ij_zt, p_sfc )
  !      Lscale = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzm, new_gr%nzm, &
  !                                     gr%zm, new_gr%zm, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     Lscale, &
  !                                     iv_other, R_ij_zt, p_sfc )
!
!
  !      ! Remap all zm values
  !      wprtp_forcing = remap_vals_to_target( ngrdcol, &
  !                                            gr%nzt+2, new_gr%nzt+2, &
  !                                            levels_source_zm_vals, &
  !                                            levels_target_zm_vals, &
  !                                            total_idx_rho_lin_spline, &
  !                                            rho_lin_spline_vals, &
  !                                            rho_lin_spline_levels, &
  !                                            wprtp_forcing, &
  !                                            iv_other, R_ij_zm, p_sfc )
  !      wpthlp_forcing = remap_vals_to_target( ngrdcol, &
  !                                             gr%nzt+2, new_gr%nzt+2, &
  !                                             levels_source_zm_vals, &
  !                                             levels_target_zm_vals, &
  !                                             total_idx_rho_lin_spline, &
  !                                             rho_lin_spline_vals, &
  !                                             rho_lin_spline_levels, &
  !                                             wpthlp_forcing, &
  !                                             iv_other, R_ij_zm, p_sfc )
  !      rtp2_forcing = remap_vals_to_target( ngrdcol, &
  !                                           gr%nzt+2, new_gr%nzt+2, &
  !                                           levels_source_zm_vals, &
  !                                           levels_target_zm_vals, &
  !                                           total_idx_rho_lin_spline, &
  !                                           rho_lin_spline_vals, &
  !                                           rho_lin_spline_levels, &
  !                                           rtp2_forcing, &
  !                                           iv_other, R_ij_zm, p_sfc )
  !      thlp2_forcing = remap_vals_to_target( ngrdcol, &
  !                                            gr%nzt+2, new_gr%nzt+2, &
  !                                            levels_source_zm_vals, &
  !                                            levels_target_zm_vals, &
  !                                            total_idx_rho_lin_spline, &
  !                                            rho_lin_spline_vals, &
  !                                            rho_lin_spline_levels, &
  !                                            thlp2_forcing, &
  !                                            iv_other, R_ij_zm, p_sfc )
  !      rtpthlp_forcing = remap_vals_to_target( ngrdcol, &
  !                                              gr%nzt+2, new_gr%nzt+2, &
  !                                              levels_source_zm_vals, &
  !                                              levels_target_zm_vals, &
  !                                              total_idx_rho_lin_spline, &
  !                                              rho_lin_spline_vals, &
  !                                              rho_lin_spline_levels, &
  !                                              rtpthlp_forcing, &
  !                                              iv_other, R_ij_zm, p_sfc )
  !      wm_zm = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzt+2, new_gr%nzt+2, &
  !                                    levels_source_zm_vals, &
  !                                    levels_target_zm_vals, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    wm_zm, &
  !                                    iv_other, R_ij_zm, p_sfc )
  !      rho_zm = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzt+2, new_gr%nzt+2, &
  !                                     levels_source_zm_vals, &
  !                                     levels_target_zm_vals, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     rho_zm, &
  !                                     iv_other, R_ij_zm, p_sfc )
  !      rho_ds_zm = remap_vals_to_target( ngrdcol, &
  !                                        gr%nzt+2, new_gr%nzt+2, &
  !                                        levels_source_zm_vals, &
  !                                        levels_target_zm_vals, &
  !                                        total_idx_rho_lin_spline, &
  !                                        rho_lin_spline_vals, &
  !                                        rho_lin_spline_levels, &
  !                                        rho_ds_zm, &
  !                                        iv_other, R_ij_zm, p_sfc )
  !      invrs_rho_ds_zm = remap_vals_to_target( ngrdcol, &
  !                                              gr%nzt+2, new_gr%nzt+2, &
  !                                              levels_source_zm_vals, &
  !                                              levels_target_zm_vals, &
  !                                              total_idx_rho_lin_spline, &
  !                                              rho_lin_spline_vals, &
  !                                              rho_lin_spline_levels, &
  !                                              invrs_rho_ds_zm, &
  !                                              iv_other, R_ij_zm, p_sfc )
  !      thv_ds_zm = remap_vals_to_target( ngrdcol, &
  !                                        gr%nzt+2, new_gr%nzt+2, &
  !                                        levels_source_zm_vals, &
  !                                        levels_target_zm_vals, &
  !                                        total_idx_rho_lin_spline, &
  !                                        rho_lin_spline_vals, &
  !                                        rho_lin_spline_levels, &
  !                                        thv_ds_zm, &
  !                                        iv_other, R_ij_zm, p_sfc )
  !      do i = 1, hydromet_dim
  !          wphydrometp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                     gr%nzt+2, new_gr%nzt+2, &
  !                                                     levels_source_zm_vals, &
  !                                                     levels_target_zm_vals, &
  !                                                     total_idx_rho_lin_spline, &
  !                                                     rho_lin_spline_vals, &
  !                                                     rho_lin_spline_levels, &
  !                                                     wphydrometp(:,:,i), &
  !                                                     iv_other, R_ij_zm, p_sfc )
  !      end do
  !      upwp = remap_vals_to_target( ngrdcol, &
  !                                   gr%nzt+2, new_gr%nzt+2, &
  !                                   levels_source_zm_vals, &
  !                                   levels_target_zm_vals, &
  !                                   total_idx_rho_lin_spline, &
  !                                   rho_lin_spline_vals, &
  !                                   rho_lin_spline_levels, &
  !                                   upwp, &
  !                                   iv_other, R_ij_zm, p_sfc )
  !      vpwp = remap_vals_to_target( ngrdcol, &
  !                                   gr%nzt+2, new_gr%nzt+2, &
  !                                   levels_source_zm_vals, &
  !                                   levels_target_zm_vals, &
  !                                   total_idx_rho_lin_spline, &
  !                                   rho_lin_spline_vals, &
  !                                   rho_lin_spline_levels, &
  !                                   vpwp, &
  !                                   iv_other, R_ij_zm, p_sfc )
  !      up2 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzt+2, new_gr%nzt+2, &
  !                                  levels_source_zm_vals, &
  !                                  levels_target_zm_vals, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  up2, &
  !                                  iv_other, R_ij_zm, p_sfc )
  !      vp2 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzt+2, new_gr%nzt+2, &
  !                                  levels_source_zm_vals, &
  !                                  levels_target_zm_vals, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  vp2, &
  !                                  iv_other, R_ij_zm, p_sfc )
  !      wprtp = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzt+2, new_gr%nzt+2, &
  !                                    levels_source_zm_vals, &
  !                                    levels_target_zm_vals, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    wprtp, &
  !                                    iv_other, R_ij_zm, p_sfc )
  !      wpthlp = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzt+2, new_gr%nzt+2, &
  !                                     levels_source_zm_vals, &
  !                                     levels_target_zm_vals, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     wpthlp, &
  !                                     iv_other, R_ij_zm, p_sfc )
  !      rtp2 = remap_vals_to_target( ngrdcol, &
  !                                   gr%nzt+2, new_gr%nzt+2, &
  !                                   levels_source_zm_vals, &
  !                                   levels_target_zm_vals, &
  !                                   total_idx_rho_lin_spline, &
  !                                   rho_lin_spline_vals, &
  !                                   rho_lin_spline_levels, &
  !                                   rtp2, &
  !                                   iv_other, R_ij_zm, p_sfc )
  !      thlp2 = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzt+2, new_gr%nzt+2, &
  !                                    levels_source_zm_vals, &
  !                                    levels_target_zm_vals, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    thlp2, &
  !                                    iv_other, R_ij_zm, p_sfc )
  !      rtpthlp = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzt+2, new_gr%nzt+2, &
  !                                      levels_source_zm_vals, &
  !                                      levels_target_zm_vals, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      rtpthlp, &
  !                                      iv_other, R_ij_zm, p_sfc )
  !      wp2 = remap_vals_to_target( ngrdcol, &
  !                                  gr%nzt+2, new_gr%nzt+2, &
  !                                  levels_source_zm_vals, &
  !                                  levels_target_zm_vals, &
  !                                  total_idx_rho_lin_spline, &
  !                                  rho_lin_spline_vals, &
  !                                  rho_lin_spline_levels, &
  !                                  wp2, &
  !                                  iv_other, R_ij_zm, p_sfc )
  !      do i = 1, sclr_dim
  !          wpsclrp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                 gr%nzt+2, new_gr%nzt+2, &
  !                                                 levels_source_zm_vals, &
  !                                                 levels_target_zm_vals, &
  !                                                 total_idx_rho_lin_spline, &
  !                                                 rho_lin_spline_vals, &
  !                                                 rho_lin_spline_levels, &
  !                                                 wpsclrp(:,:,i), &
  !                                                 iv_other, R_ij_zm, p_sfc )
  !          sclrp2(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                gr%nzt+2, new_gr%nzt+2, &
  !                                                levels_source_zm_vals, &
  !                                                levels_target_zm_vals, &
  !                                                total_idx_rho_lin_spline, &
  !                                                rho_lin_spline_vals, &
  !                                                rho_lin_spline_levels, &
  !                                                sclrp2(:,:,i), &
  !                                                iv_other, R_ij_zm, p_sfc )
  !          sclrprtp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                  gr%nzt+2, new_gr%nzt+2, &
  !                                                  levels_source_zm_vals, &
  !                                                  levels_target_zm_vals, &
  !                                                  total_idx_rho_lin_spline, &
  !                                                  rho_lin_spline_vals, &
  !                                                  rho_lin_spline_levels, &
  !                                                  sclrprtp(:,:,i), &
  !                                                  iv_other, R_ij_zm, p_sfc )
  !          sclrpthlp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                   gr%nzt+2, new_gr%nzt+2, &
  !                                                   levels_source_zm_vals, &
  !                                                   levels_target_zm_vals, &
  !                                                   total_idx_rho_lin_spline, &
  !                                                   rho_lin_spline_vals, &
  !                                                   rho_lin_spline_levels, &
  !                                                   sclrpthlp(:,:,i), &
  !                                                   iv_other, R_ij_zm, p_sfc )
  !          sclrpthvp(:,:,i) = remap_vals_to_target( ngrdcol, &
  !                                                   gr%nzt+2, new_gr%nzt+2, &
  !                                                   levels_source_zm_vals, &
  !                                                   levels_target_zm_vals, &
  !                                                   total_idx_rho_lin_spline, &
  !                                                   rho_lin_spline_vals, &
  !                                                   rho_lin_spline_levels, &
  !                                                   sclrpthvp(:,:,i), &
  !                                                   iv_other, R_ij_zm, p_sfc )
  !      end do
  !      wpthvp = remap_vals_to_target( ngrdcol, &
  !                                     gr%nzt+2, new_gr%nzt+2, &
  !                                     levels_source_zm_vals, &
  !                                     levels_target_zm_vals, &
  !                                     total_idx_rho_lin_spline, &
  !                                     rho_lin_spline_vals, &
  !                                     rho_lin_spline_levels, &
  !                                     wpthvp, &
  !                                     iv_other, R_ij_zm, p_sfc )
  !      rtpthvp = remap_vals_to_target( ngrdcol, &
  !                                      gr%nzt+2, new_gr%nzt+2, &
  !                                      levels_source_zm_vals, &
  !                                      levels_target_zm_vals, &
  !                                      total_idx_rho_lin_spline, &
  !                                      rho_lin_spline_vals, &
  !                                      rho_lin_spline_levels, &
  !                                      rtpthvp, &
  !                                      iv_other, R_ij_zm, p_sfc )
  !      thlpthvp = remap_vals_to_target( ngrdcol, &
  !                                       gr%nzt+2, new_gr%nzt+2, &
  !                                       levels_source_zm_vals, &
  !                                       levels_target_zm_vals, &
  !                                       total_idx_rho_lin_spline, &
  !                                       rho_lin_spline_vals, &
  !                                       rho_lin_spline_levels, &
  !                                       thlpthvp, &
  !                                       iv_other, R_ij_zm, p_sfc )
  !      uprcp = remap_vals_to_target( ngrdcol, &
  !                                    gr%nzt+2, new_gr%nzt+2, &
  !                                    levels_source_zm_vals, &
  !                                    levels_target_zm_vals, &
  !                                    total_idx_rho_lin_spline, &
  !                                    rho_lin_spline_vals, &
  !                                    rho_lin_spline_levels, &
  !                                    uprcp, &
  !                                    iv_other, R_ij_zm, p_sfc )
  !      vprcp = remap_vals_to_target( ngrdcol, &
   !                                   gr%nzt+2, new_gr%nzt+2, &
   !                                   levels_source_zm_vals, &
   !                                   levels_target_zm_vals, &
   !                                   total_idx_rho_lin_spline, &
   !                                   rho_lin_spline_vals, &
   !                                   rho_lin_spline_levels, &
   !                                   vprcp, &
   !                                   iv_other, R_ij_zm, p_sfc )
   !     rc_coef_zm = remap_vals_to_target( ngrdcol, &
   !                                        gr%nzt+2, new_gr%nzt+2, &
   !                                        levels_source_zm_vals, &
   !                                        levels_target_zm_vals, &
   !                                        total_idx_rho_lin_spline, &
   !                                        rho_lin_spline_vals, &
   !                                        rho_lin_spline_levels, &
   !                                        rc_coef_zm, &
   !                                        iv_other, R_ij_zm, p_sfc )
   !     wp4 = remap_vals_to_target( ngrdcol, &
   !                                 gr%nzt+2, new_gr%nzt+2, &
   !                                 levels_source_zm_vals, &
   !                                 levels_target_zm_vals, &
   !                                 total_idx_rho_lin_spline, &
   !                                 rho_lin_spline_vals, &
   !                                 rho_lin_spline_levels, &
   !                                 wp4, &
   !                                 iv_other, R_ij_zm, p_sfc )
   !     wp2up2 = remap_vals_to_target( ngrdcol, &
   !                                    gr%nzt+2, new_gr%nzt+2, &
   !                                    levels_source_zm_vals, &
   !                                    levels_target_zm_vals, &
   !                                    total_idx_rho_lin_spline, &
   !                                    rho_lin_spline_vals, &
   !                                    rho_lin_spline_levels, &
   !                                    wp2up2, &
   !                                    iv_other, R_ij_zm, p_sfc )
   !     wp2vp2 = remap_vals_to_target( ngrdcol, &
   !                                    gr%nzt+2, new_gr%nzt+2, &
   !                                    levels_source_zm_vals, &
   !                                    levels_target_zm_vals, &
   !                                    total_idx_rho_lin_spline, &
   !                                    rho_lin_spline_vals, &
   !                                    rho_lin_spline_levels, &
   !                                    wp2vp2, &
   !                                    iv_other, R_ij_zm, p_sfc )
   !     upwp_pert = remap_vals_to_target( ngrdcol, &
   !                                       gr%nzt+2, new_gr%nzt+2, &
   !                                       levels_source_zm_vals, &
   !                                       levels_target_zm_vals, &
   !                                       total_idx_rho_lin_spline, &
   !                                       rho_lin_spline_vals, &
   !                                       rho_lin_spline_levels, &
   !                                       upwp_pert, &
   !                                       iv_other, R_ij_zm, p_sfc )
   !     vpwp_pert = remap_vals_to_target( ngrdcol, &
   !                                       gr%nzt+2, new_gr%nzt+2, &
   !                                       levels_source_zm_vals, &
   !                                       levels_target_zm_vals, &
   !                                       total_idx_rho_lin_spline, &
   !                                       rho_lin_spline_vals, &
   !                                       rho_lin_spline_levels, &
   !                                       vpwp_pert, &
   !                                       iv_other, R_ij_zm, p_sfc )
   !     wprcp = remap_vals_to_target( ngrdcol, &
   !                                   gr%nzt+2, new_gr%nzt+2, &
   !                                   levels_source_zm_vals, &
   !                                   levels_target_zm_vals, &
   !                                   total_idx_rho_lin_spline, &
   !                                   rho_lin_spline_vals, &
   !                                   rho_lin_spline_levels, &
   !                                   wprcp, &
   !                                   iv_other, R_ij_zm, p_sfc )
   !     invrs_tau_zm = remap_vals_to_target( ngrdcol, &
   !                                          gr%nzt+2, new_gr%nzt+2, &
   !                                          levels_source_zm_vals, &
   !                                          levels_target_zm_vals, &
   !                                          total_idx_rho_lin_spline, &
   !                                          rho_lin_spline_vals, &
   !                                          rho_lin_spline_levels, &
   !                                          invrs_tau_zm, &
   !                                          iv_other, R_ij_zm, p_sfc )
   !     Kh_zm = remap_vals_to_target( ngrdcol, &
   !                                   gr%nzt+2, new_gr%nzt+2, &
   !                                   levels_source_zm_vals, &
   !                                   levels_target_zm_vals, &
   !                                   total_idx_rho_lin_spline, &
   !                                   rho_lin_spline_vals, &
   !                                   rho_lin_spline_levels, &
   !                                   Kh_zm, &
   !                                   iv_other, R_ij_zm, p_sfc )
   !     thlprcp = remap_vals_to_target( ngrdcol, &
   !                                     gr%nzt+2, new_gr%nzt+2, &
   !                                     levels_source_zm_vals, &
   !                                     levels_target_zm_vals, &
   !                                     total_idx_rho_lin_spline, &
   !                                     rho_lin_spline_vals, &
   !                                     rho_lin_spline_levels, &
   !                                     thlprcp, &
   !                                     iv_other, R_ij_zm, p_sfc )
   !   else
   !     write(fstderr,*) 'There is no method implemented for grid_remap_method=', &
   !                      grid_remap_method, '. Try another integer value.'
   !     error stop 'Invalid value for grid_remap_method.'
   !   end if
!
   !   gr = new_gr
   ! end if
!
  !end subroutine adapt_grid_old

  subroutine remap_all_clubb_core_vals( ngrdcol, total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, rho_lin_spline_levels, &
                                        l_implemented, &
                                        hydromet_dim, sclr_dim, edsclr_dim, &
                                        thvm, idx_thvm, p_sfc, &
                                        grid_remap_method, &
                                        gr, new_gr, &
                                        thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
                                        sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
                                        wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
                                        rtpthlp_forcing, wm_zm, wm_zt, &
                                        rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
                                        p_in_Pa, rho_zm, rho, exner, &
                                        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                                        invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                                        rfrzm, wphydrometp, &
                                        wp2hmp, rtphmp_zt, thlphmp_zt, &
                                        um, vm, upwp, vpwp, up2, vp2, up3, vp3, &
                                        thlm, rtm, wprtp, wpthlp, &
                                        wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
                                        sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &
                                        wpsclrp, edsclrm, &
                                        rcm, cloud_frac, &
                                        wpthvp, wp2thvp, rtpthvp, thlpthvp, &
                                        sclrpthvp, &
                                        wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &
                                        wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &
                                        um_pert, vm_pert, upwp_pert, vpwp_pert, &
                                        Kh_zm, Kh_zt, &
                                        thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &
                                        cloudy_updraft_frac, cloudy_downdraft_frac, &
                                        rcm_in_layer, cloud_cover, invrs_tau_zm, &
                                        Lscale )
    ! Description:
    ! Adapts the grid based on the density function and interpolates all values to the new grid.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)
    
    use grid_class, only: &
      setup_grid

    use calc_pressure, only: &
      init_pressure    ! Procedure(s)

    use model_flags, only: &
      cons_ullrich_remap

    use error_code, only: &
      clubb_fatal_error

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      hydromet_dim, &
      sclr_dim, &
      edsclr_dim, &
      total_idx_rho_lin_spline, &  ! total number of indices of the rho density
                                   ! piecewise linear function []
      idx_thvm   ! numer of indices of thvm

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) ::  &
      rho_lin_spline_vals, &
      rho_lin_spline_levels
    
    ! TODO rename to gr_source and gr_target?
    type( grid ), intent(in) :: &
      gr, &
      new_gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      p_sfc

    logical, intent(in) :: l_implemented

    real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
      thvm

    integer, intent(in) :: &
      grid_remap_method ! specifies what remapping method should be used

    !--------------------- Output Variable ---------------------

    !--------------------- In/Out Variable ---------------------
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm, hydromet_dim), intent(inout) :: &
      wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt, hydromet_dim), intent(inout) :: &
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and a hydrometeor [K <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) :: &
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt,edsclr_dim) :: &
        edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Lscale     ! Length scale         [m]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzt) ::  &
      thvm_new_grid

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      p_in_Pa_zm, exner_zm ! just placeholder variables; are only needed for subroutine call

    real( kind = core_rknd ), dimension(gr%nzm,gr%nzm) :: &
      R_ij_zm

    real( kind = core_rknd ), dimension(gr%nzt,gr%nzt) :: &
      R_ij_zt

    real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_source_zm_vals

    real( kind=core_rknd ), dimension(ngrdcol,gr%nzt+2) :: levels_target_zm_vals

    integer :: &
      iv_wind = -1, &
      iv_pos_def = 0, &
      iv_other = 1

    !--------------------- Begin Code ---------------------

      call remapping_matrix_zm_values( ngrdcol, &                                      ! Intent(in)
                                       gr, new_gr, &                                   ! Intent(in)
                                       total_idx_rho_lin_spline, &                     ! Intent(in)
                                       rho_lin_spline_vals, &                          ! Intent(in)
                                       rho_lin_spline_levels, &                        ! Intent(in)
                                       levels_source_zm_vals, levels_target_zm_vals, & ! Intent(out)
                                       R_ij_zm )                                       ! Intent(out)

      call remapping_matrix_zt_values( ngrdcol, &                     ! Intent(in)
                                       gr, new_gr, &                  ! Intent(in)
                                       total_idx_rho_lin_spline, &    ! Intent(in)
                                       rho_lin_spline_vals, &         ! Intent(in)
                                       rho_lin_spline_levels, &       ! Intent(in)
                                       R_ij_zt )                      ! Intent(out)

      ! integrate grid_remap_method into remap_vals_to_target function
      if ( grid_remap_method == cons_ullrich_remap ) then

        ! Remap all zt values
        thlm_forcing = remap_vals_to_target( ngrdcol, &
                                             gr%nzm, new_gr%nzm, &
                                             gr%zm, new_gr%zm, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             thlm_forcing, &
                                             iv_other, R_ij_zt, p_sfc )
        rtm_forcing = remap_vals_to_target( ngrdcol, &
                                            gr%nzm, new_gr%nzm, &
                                            gr%zm, new_gr%zm, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            rtm_forcing, &
                                            iv_other, R_ij_zt, p_sfc )
        um_forcing = remap_vals_to_target( ngrdcol, &
                                           gr%nzm, new_gr%nzm, &
                                           gr%zm, new_gr%zm, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           um_forcing, &
                                           iv_other, R_ij_zt, p_sfc )
        vm_forcing = remap_vals_to_target( ngrdcol, &
                                           gr%nzm, new_gr%nzm, &
                                           gr%zm, new_gr%zm, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           vm_forcing, &
                                           iv_other, R_ij_zt, p_sfc )
        wm_zt = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wm_zt, &
                                      iv_wind, R_ij_zt, p_sfc )
        rho = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    rho, &
                                    iv_pos_def, R_ij_zt, p_sfc )
        rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                          gr%nzm, new_gr%nzm, &
                                          gr%zm, new_gr%zm, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          rho_ds_zt, &
                                          iv_pos_def, R_ij_zt, p_sfc )
        invrs_rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                                gr%nzm, new_gr%nzm, &
                                                gr%zm, new_gr%zm, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                invrs_rho_ds_zt, &
                                                iv_pos_def, R_ij_zt, p_sfc )
        thv_ds_zt = remap_vals_to_target( ngrdcol, &
                                          gr%nzm, new_gr%nzm, &
                                          gr%zm, new_gr%zm, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          thv_ds_zt, &
                                          iv_pos_def, R_ij_zt, p_sfc )
        rfrzm = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      rfrzm, &
                                      iv_pos_def, R_ij_zt, p_sfc )
        do i = 1, hydromet_dim
            wp2hmp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                  gr%nzm, new_gr%nzm, &
                                                  gr%zm, new_gr%zm, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  wp2hmp(:,:,i), &
                                                  iv_other, R_ij_zt, p_sfc )
            rtphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                     gr%nzm, new_gr%nzm, &
                                                     gr%zm, new_gr%zm, &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals, &
                                                     rho_lin_spline_levels, &
                                                     rtphmp_zt(:,:,i), &
                                                     iv_other, R_ij_zt, p_sfc )
            thlphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                      gr%nzm, new_gr%nzm, &
                                                      gr%zm, new_gr%zm, &
                                                      total_idx_rho_lin_spline, &
                                                      rho_lin_spline_vals, &
                                                      rho_lin_spline_levels, &
                                                      thlphmp_zt(:,:,i), &
                                                      iv_other, R_ij_zt, p_sfc )
        end do
        do i = 1, sclr_dim
            sclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                         gr%nzm, new_gr%nzm, &
                                                         gr%zm, new_gr%zm, &
                                                         total_idx_rho_lin_spline, &
                                                         rho_lin_spline_vals, &
                                                         rho_lin_spline_levels, &
                                                         sclrm_forcing(:,:,i), &
                                                         iv_other, R_ij_zt, p_sfc )
            sclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                 gr%nzm, new_gr%nzm, &
                                                 gr%zm, new_gr%zm, &
                                                 total_idx_rho_lin_spline, &
                                                 rho_lin_spline_vals, &
                                                 rho_lin_spline_levels, &
                                                 sclrm(:,:,i), &
                                                 iv_other, R_ij_zt, p_sfc )
            sclrp3(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                  gr%nzm, new_gr%nzm, &
                                                  gr%zm, new_gr%zm, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  sclrp3(:,:,i), &
                                                  iv_other, R_ij_zt, p_sfc )
        end do
        do i = 1, edsclr_dim
            edsclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                           gr%nzm, new_gr%nzm, &
                                                           gr%zm, new_gr%zm, &
                                                           total_idx_rho_lin_spline, &
                                                           rho_lin_spline_vals, &
                                                           rho_lin_spline_levels, &
                                                           edsclrm_forcing(:,:,i), &
                                                           iv_other, R_ij_zt, p_sfc )
            edsclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                   gr%nzm, new_gr%nzm, &
                                                   gr%zm, new_gr%zm, &
                                                   total_idx_rho_lin_spline, &
                                                   rho_lin_spline_vals, &
                                                   rho_lin_spline_levels, &
                                                   edsclrm(:,:,i), &
                                                   iv_other, R_ij_zt, p_sfc )
        end do
        rtm_ref = remap_vals_to_target( ngrdcol, &
                                        gr%nzm, new_gr%nzm, &
                                        gr%zm, new_gr%zm, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        rtm_ref, &
                                        iv_pos_def, R_ij_zt, p_sfc )
        thlm_ref = remap_vals_to_target( ngrdcol, &
                                         gr%nzm, new_gr%nzm, &
                                         gr%zm, new_gr%zm, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         thlm_ref, &
                                         iv_pos_def, R_ij_zt, p_sfc )
        um_ref = remap_vals_to_target( ngrdcol, &
                                       gr%nzm, new_gr%nzm, &
                                       gr%zm, new_gr%zm, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       um_ref, &
                                       iv_wind, R_ij_zt, p_sfc )
        vm_ref = remap_vals_to_target( ngrdcol, &
                                       gr%nzm, new_gr%nzm, &
                                       gr%zm, new_gr%zm, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       vm_ref, &
                                       iv_wind, R_ij_zt, p_sfc )
        ug = remap_vals_to_target( ngrdcol, &
                                   gr%nzm, new_gr%nzm, &
                                   gr%zm, new_gr%zm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   ug, &
                                   iv_wind, R_ij_zt, p_sfc )
        vg = remap_vals_to_target( ngrdcol, &
                                   gr%nzm, new_gr%nzm, &
                                   gr%zm, new_gr%zm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   vg, &
                                   iv_wind, R_ij_zt, p_sfc )
        um = remap_vals_to_target( ngrdcol, &
                                   gr%nzm, new_gr%nzm, &
                                   gr%zm, new_gr%zm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   um, &
                                   iv_wind, R_ij_zt, p_sfc )
        vm = remap_vals_to_target( ngrdcol, &
                                   gr%nzm, new_gr%nzm, &
                                   gr%zm, new_gr%zm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   vm, &
                                   iv_wind, R_ij_zt, p_sfc )
        up3 = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    up3, &
                                    iv_other, R_ij_zt, p_sfc )
        vp3 = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    vp3, &
                                    iv_other, R_ij_zt, p_sfc )
        rtm = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    rtm, &
                                    iv_pos_def, R_ij_zt, p_sfc )
        thlm = remap_vals_to_target( ngrdcol, &
                                     gr%nzm, new_gr%nzm, &
                                     gr%zm, new_gr%zm, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     thlm, &
                                     iv_pos_def, R_ij_zt, p_sfc )
        rtp3 = remap_vals_to_target( ngrdcol, &
                                     gr%nzm, new_gr%nzm, &
                                     gr%zm, new_gr%zm, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     rtp3, &
                                     iv_other, R_ij_zt, p_sfc )
        thlp3 = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      thlp3, &
                                      iv_other, R_ij_zt, p_sfc )
        wp3 = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    wp3, &
                                    iv_other, R_ij_zt, p_sfc )
        ! remap thvm values to new grid, to calculate new p_in_Pa
        thvm_new_grid = remap_vals_to_target( ngrdcol, &
                                              gr%nzm, new_gr%nzm, &
                                              gr%zm, new_gr%zm, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              thvm, &
                                              iv_pos_def, R_ij_zt, p_sfc )
        ! calculate p_in_Pa instead of remapping directly since it can run into problems if for
        ! example the two highest levels have the same pressure value, which could be happening with
        ! the remapping, also calculate exner accordingly
        call init_pressure( ngrdcol, new_gr, thvm_new_grid, p_sfc, &
                            p_in_Pa, exner, p_in_Pa_zm, exner_zm )
        rcm = remap_vals_to_target( ngrdcol, &
                                    gr%nzm, new_gr%nzm, &
                                    gr%zm, new_gr%zm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    rcm, &
                                    iv_pos_def, R_ij_zt, p_sfc )
        cloud_frac = remap_vals_to_target( ngrdcol, &
                                           gr%nzm, new_gr%nzm, &
                                           gr%zm, new_gr%zm, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           cloud_frac, &
                                           iv_pos_def, R_ij_zt, p_sfc )
        wp2thvp = remap_vals_to_target( ngrdcol, &
                                        gr%nzm, new_gr%nzm, &
                                        gr%zm, new_gr%zm, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        wp2thvp, &
                                        iv_other, R_ij_zt, p_sfc )
        wp2rtp = remap_vals_to_target( ngrdcol, &
                                       gr%nzm, new_gr%nzm, &
                                       gr%zm, new_gr%zm, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wp2rtp, &
                                       iv_other, R_ij_zt, p_sfc )
        wp2thlp = remap_vals_to_target( ngrdcol, &
                                        gr%nzm, new_gr%nzm, &
                                        gr%zm, new_gr%zm, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        wp2thlp, &
                                        iv_other, R_ij_zt, p_sfc )
        wpup2 = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wpup2, &
                                      iv_other, R_ij_zt, p_sfc )
        wpvp2 = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wpvp2, &
                                      iv_other, R_ij_zt, p_sfc )
        ice_supersat_frac = remap_vals_to_target( ngrdcol, &
                                                  gr%nzm, new_gr%nzm, &
                                                  gr%zm, new_gr%zm, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  ice_supersat_frac, &
                                                  iv_pos_def, R_ij_zt, p_sfc )
        um_pert = remap_vals_to_target( ngrdcol, &
                                        gr%nzm, new_gr%nzm, &
                                        gr%zm, new_gr%zm, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        um_pert, &
                                        iv_wind, R_ij_zt, p_sfc )
        vm_pert = remap_vals_to_target( ngrdcol, &
                                        gr%nzm, new_gr%nzm, &
                                        gr%zm, new_gr%zm, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        vm_pert, &
                                        iv_wind, R_ij_zt, p_sfc )
        rcm_in_layer = remap_vals_to_target( ngrdcol, &
                                             gr%nzm, new_gr%nzm, &
                                             gr%zm, new_gr%zm, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             rcm_in_layer, &
                                             iv_pos_def, R_ij_zt, p_sfc )
        cloud_cover = remap_vals_to_target( ngrdcol, &
                                            gr%nzm, new_gr%nzm, &
                                            gr%zm, new_gr%zm, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            cloud_cover, &
                                            iv_pos_def, R_ij_zt, p_sfc )
        w_up_in_cloud = remap_vals_to_target( ngrdcol, &
                                              gr%nzm, new_gr%nzm, &
                                              gr%zm, new_gr%zm, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              w_up_in_cloud, &
                                              iv_other, R_ij_zt, p_sfc )
        w_down_in_cloud = remap_vals_to_target( ngrdcol, &
                                                gr%nzm, new_gr%nzm, &
                                                gr%zm, new_gr%zm, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                w_down_in_cloud, &
                                                iv_other, R_ij_zt, p_sfc )
        cloudy_updraft_frac = remap_vals_to_target( ngrdcol, &
                                                    gr%nzm, new_gr%nzm, &
                                                    gr%zm, new_gr%zm, &
                                                    total_idx_rho_lin_spline, &
                                                    rho_lin_spline_vals, &
                                                    rho_lin_spline_levels, &
                                                    cloudy_updraft_frac, &
                                                    iv_pos_def, R_ij_zt, p_sfc )
        cloudy_downdraft_frac = remap_vals_to_target( ngrdcol, &
                                                      gr%nzm, new_gr%nzm, &
                                                      gr%zm, new_gr%zm, &
                                                      total_idx_rho_lin_spline, &
                                                      rho_lin_spline_vals, &
                                                      rho_lin_spline_levels, &
                                                      cloudy_downdraft_frac, &
                                                      iv_pos_def, R_ij_zt, p_sfc )
        Kh_zt = remap_vals_to_target( ngrdcol, &
                                      gr%nzm, new_gr%nzm, &
                                      gr%zm, new_gr%zm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      Kh_zt, &
                                      iv_pos_def, R_ij_zt, p_sfc )
        Lscale = remap_vals_to_target( ngrdcol, &
                                       gr%nzm, new_gr%nzm, &
                                       gr%zm, new_gr%zm, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       Lscale, &
                                       iv_pos_def, R_ij_zt, p_sfc )


        ! Remap all zm values
        wprtp_forcing = remap_vals_to_target( ngrdcol, &
                                              gr%nzt+2, new_gr%nzt+2, &
                                              levels_source_zm_vals, &
                                              levels_target_zm_vals, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              wprtp_forcing, &
                                              iv_other, R_ij_zm, p_sfc )
        wpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                               gr%nzt+2, new_gr%nzt+2, &
                                               levels_source_zm_vals, &
                                               levels_target_zm_vals, &
                                               total_idx_rho_lin_spline, &
                                               rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               wpthlp_forcing, &
                                               iv_other, R_ij_zm, p_sfc )
        rtp2_forcing = remap_vals_to_target( ngrdcol, &
                                             gr%nzt+2, new_gr%nzt+2, &
                                             levels_source_zm_vals, &
                                             levels_target_zm_vals, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             rtp2_forcing, &
                                             iv_other, R_ij_zm, p_sfc )
        thlp2_forcing = remap_vals_to_target( ngrdcol, &
                                              gr%nzt+2, new_gr%nzt+2, &
                                              levels_source_zm_vals, &
                                              levels_target_zm_vals, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              thlp2_forcing, &
                                              iv_other, R_ij_zm, p_sfc )
        rtpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                                gr%nzt+2, new_gr%nzt+2, &
                                                levels_source_zm_vals, &
                                                levels_target_zm_vals, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                rtpthlp_forcing, &
                                                iv_other, R_ij_zm, p_sfc )
        wm_zm = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wm_zm, &
                                      iv_wind, R_ij_zm, p_sfc )
        rho_zm = remap_vals_to_target( ngrdcol, &
                                       gr%nzt+2, new_gr%nzt+2, &
                                       levels_source_zm_vals, &
                                       levels_target_zm_vals, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       rho_zm, &
                                       iv_pos_def, R_ij_zm, p_sfc )
        rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                          gr%nzt+2, new_gr%nzt+2, &
                                          levels_source_zm_vals, &
                                          levels_target_zm_vals, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          rho_ds_zm, &
                                          iv_pos_def, R_ij_zm, p_sfc )
        invrs_rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                                gr%nzt+2, new_gr%nzt+2, &
                                                levels_source_zm_vals, &
                                                levels_target_zm_vals, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                invrs_rho_ds_zm, &
                                                iv_pos_def, R_ij_zm, p_sfc )
        thv_ds_zm = remap_vals_to_target( ngrdcol, &
                                          gr%nzt+2, new_gr%nzt+2, &
                                          levels_source_zm_vals, &
                                          levels_target_zm_vals, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          thv_ds_zm, &
                                          iv_pos_def, R_ij_zm, p_sfc )
        do i = 1, hydromet_dim
            wphydrometp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                       gr%nzt+2, new_gr%nzt+2, &
                                                       levels_source_zm_vals, &
                                                       levels_target_zm_vals, &
                                                       total_idx_rho_lin_spline, &
                                                       rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       wphydrometp(:,:,i), &
                                                       iv_other, R_ij_zm, p_sfc )
        end do
        upwp = remap_vals_to_target( ngrdcol, &
                                     gr%nzt+2, new_gr%nzt+2, &
                                     levels_source_zm_vals, &
                                     levels_target_zm_vals, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     upwp, &
                                     iv_other, R_ij_zm, p_sfc )
        vpwp = remap_vals_to_target( ngrdcol, &
                                     gr%nzt+2, new_gr%nzt+2, &
                                     levels_source_zm_vals, &
                                     levels_target_zm_vals, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     vpwp, &
                                     iv_other, R_ij_zm, p_sfc )
        up2 = remap_vals_to_target( ngrdcol, &
                                    gr%nzt+2, new_gr%nzt+2, &
                                    levels_source_zm_vals, &
                                    levels_target_zm_vals, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    up2, &
                                    iv_pos_def, R_ij_zm, p_sfc )
        vp2 = remap_vals_to_target( ngrdcol, &
                                    gr%nzt+2, new_gr%nzt+2, &
                                    levels_source_zm_vals, &
                                    levels_target_zm_vals, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    vp2, &
                                    iv_pos_def, R_ij_zm, p_sfc )
        wprtp = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wprtp, &
                                      iv_other, R_ij_zm, p_sfc )
        wpthlp = remap_vals_to_target( ngrdcol, &
                                       gr%nzt+2, new_gr%nzt+2, &
                                       levels_source_zm_vals, &
                                       levels_target_zm_vals, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wpthlp, &
                                       iv_other, R_ij_zm, p_sfc )
        rtp2 = remap_vals_to_target( ngrdcol, &
                                     gr%nzt+2, new_gr%nzt+2, &
                                     levels_source_zm_vals, &
                                     levels_target_zm_vals, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     rtp2, &
                                     iv_pos_def, R_ij_zm, p_sfc )
        thlp2 = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      thlp2, &
                                      iv_pos_def, R_ij_zm, p_sfc )
        rtpthlp = remap_vals_to_target( ngrdcol, &
                                        gr%nzt+2, new_gr%nzt+2, &
                                        levels_source_zm_vals, &
                                        levels_target_zm_vals, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        rtpthlp, &
                                        iv_other, R_ij_zm, p_sfc )
        wp2 = remap_vals_to_target( ngrdcol, &
                                    gr%nzt+2, new_gr%nzt+2, &
                                    levels_source_zm_vals, &
                                    levels_target_zm_vals, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    wp2, &
                                    iv_pos_def, R_ij_zm, p_sfc )
        do i = 1, sclr_dim
            wpsclrp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                   gr%nzt+2, new_gr%nzt+2, &
                                                   levels_source_zm_vals, &
                                                   levels_target_zm_vals, &
                                                   total_idx_rho_lin_spline, &
                                                   rho_lin_spline_vals, &
                                                   rho_lin_spline_levels, &
                                                   wpsclrp(:,:,i), &
                                                   iv_other, R_ij_zm, p_sfc )
            sclrp2(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                  gr%nzt+2, new_gr%nzt+2, &
                                                  levels_source_zm_vals, &
                                                  levels_target_zm_vals, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  sclrp2(:,:,i), &
                                                  iv_other, R_ij_zm, p_sfc )
            sclrprtp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                    gr%nzt+2, new_gr%nzt+2, &
                                                    levels_source_zm_vals, &
                                                    levels_target_zm_vals, &
                                                    total_idx_rho_lin_spline, &
                                                    rho_lin_spline_vals, &
                                                    rho_lin_spline_levels, &
                                                    sclrprtp(:,:,i), &
                                                    iv_other, R_ij_zm, p_sfc )
            sclrpthlp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                     gr%nzt+2, new_gr%nzt+2, &
                                                     levels_source_zm_vals, &
                                                     levels_target_zm_vals, &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals, &
                                                     rho_lin_spline_levels, &
                                                     sclrpthlp(:,:,i), &
                                                     iv_other, R_ij_zm, p_sfc )
            sclrpthvp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                     gr%nzt+2, new_gr%nzt+2, &
                                                     levels_source_zm_vals, &
                                                     levels_target_zm_vals, &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals, &
                                                     rho_lin_spline_levels, &
                                                     sclrpthvp(:,:,i), &
                                                     iv_other, R_ij_zm, p_sfc )
        end do
        wpthvp = remap_vals_to_target( ngrdcol, &
                                       gr%nzt+2, new_gr%nzt+2, &
                                       levels_source_zm_vals, &
                                       levels_target_zm_vals, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wpthvp, &
                                       iv_other, R_ij_zm, p_sfc )
        rtpthvp = remap_vals_to_target( ngrdcol, &
                                        gr%nzt+2, new_gr%nzt+2, &
                                        levels_source_zm_vals, &
                                        levels_target_zm_vals, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        rtpthvp, &
                                        iv_other, R_ij_zm, p_sfc )
        thlpthvp = remap_vals_to_target( ngrdcol, &
                                         gr%nzt+2, new_gr%nzt+2, &
                                         levels_source_zm_vals, &
                                         levels_target_zm_vals, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         thlpthvp, &
                                         iv_other, R_ij_zm, p_sfc )
        uprcp = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      uprcp, &
                                      iv_other, R_ij_zm, p_sfc )
        vprcp = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      vprcp, &
                                      iv_other, R_ij_zm, p_sfc )
        rc_coef_zm = remap_vals_to_target( ngrdcol, &
                                           gr%nzt+2, new_gr%nzt+2, &
                                           levels_source_zm_vals, &
                                           levels_target_zm_vals, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           rc_coef_zm, &
                                           iv_pos_def, R_ij_zm, p_sfc )
        wp4 = remap_vals_to_target( ngrdcol, &
                                    gr%nzt+2, new_gr%nzt+2, &
                                    levels_source_zm_vals, &
                                    levels_target_zm_vals, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    wp4, &
                                    iv_pos_def, R_ij_zm, p_sfc )
        wp2up2 = remap_vals_to_target( ngrdcol, &
                                       gr%nzt+2, new_gr%nzt+2, &
                                       levels_source_zm_vals, &
                                       levels_target_zm_vals, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wp2up2, &
                                       iv_pos_def, R_ij_zm, p_sfc )
        wp2vp2 = remap_vals_to_target( ngrdcol, &
                                       gr%nzt+2, new_gr%nzt+2, &
                                       levels_source_zm_vals, &
                                       levels_target_zm_vals, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wp2vp2, &
                                       iv_pos_def, R_ij_zm, p_sfc )
        upwp_pert = remap_vals_to_target( ngrdcol, &
                                          gr%nzt+2, new_gr%nzt+2, &
                                          levels_source_zm_vals, &
                                          levels_target_zm_vals, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          upwp_pert, &
                                          iv_other, R_ij_zm, p_sfc )
        vpwp_pert = remap_vals_to_target( ngrdcol, &
                                          gr%nzt+2, new_gr%nzt+2, &
                                          levels_source_zm_vals, &
                                          levels_target_zm_vals, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          vpwp_pert, &
                                          iv_other, R_ij_zm, p_sfc )
        wprcp = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      wprcp, &
                                      iv_other, R_ij_zm, p_sfc )
        invrs_tau_zm = remap_vals_to_target( ngrdcol, &
                                             gr%nzt+2, new_gr%nzt+2, &
                                             levels_source_zm_vals, &
                                             levels_target_zm_vals, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             invrs_tau_zm, &
                                             iv_pos_def, R_ij_zm, p_sfc )
        Kh_zm = remap_vals_to_target( ngrdcol, &
                                      gr%nzt+2, new_gr%nzt+2, &
                                      levels_source_zm_vals, &
                                      levels_target_zm_vals, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      Kh_zm, &
                                      iv_pos_def, R_ij_zm, p_sfc )
        thlprcp = remap_vals_to_target( ngrdcol, &
                                        gr%nzt+2, new_gr%nzt+2, &
                                        levels_source_zm_vals, &
                                        levels_target_zm_vals, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        thlprcp, &
                                        iv_other, R_ij_zm, p_sfc )
      else
        write(fstderr,*) 'There is no method implemented for grid_remap_method=', &
                         grid_remap_method, '. Try another integer value.'
        error stop 'Invalid value for grid_remap_method.'
      end if

  end subroutine remap_all_clubb_core_vals

  subroutine calc_grid_dens( ngrdcol, gr, &
                             um, vm, &
                             Lscale_zt, wp2, &
                             pdf_params, &
                             thlm, exner, rtm, &
                             rcm, p_in_Pa, thvm, &
                             ice_supersat_frac, &
                             bv_efold, &
                             clubb_config_flags, &
                             gr_dens_z, gr_dens, &
                             alt_term, Lscale_term, &
                             Lscale_term_time_avg, &
                             chi_term, &
                             chi_term_time_avg, &
                             brunt_term, &
                             brunt_term_time_avg )
    ! Description:
    ! Calculate the necessary variables and then construct grid density from them.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use grid_class, only: &
        zt2zm, & ! Procedures
        ddzt

    use model_flags, only: &
        clubb_config_flags_type  ! Type

    use advance_helper_module, only: &
        calc_brunt_vaisala_freq_sqd ! Procedure

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol
        
    type( grid ), intent(in) :: &
      gr
    
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      Lscale_zt,         &  ! Length scale   [m]
      um,                & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,                & ! northward grid-mean wind component (thermo. levs.) [m/s]
      thvm,              & ! Virtual potential temperature                        [K]
      exner,             & ! Exner function (thermodynamic levels)       [-]
      rtm,               & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,              & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rcm,               & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      p_in_Pa,           & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      ice_supersat_frac    ! w'^3 (momentum levels)                    [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm), intent(in) ::  &
      wp2                  ! w'^2 on thermo. grid [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      bv_efold             ! Control parameter for inverse e-folding of cloud fraction
                           ! in the mixed Brunt Vaisala frequency  [-]

    type(pdf_parameter), intent(in) :: & 
      pdf_params           ! PDF parameters

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags   ! Derived type holding all configurable CLUBB flags

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(gr%nzm), intent(out) ::  &
      gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
                       ! grid density function [m]
      gr_dens  ! the values of the connection points of the piecewise linear
               ! grid density function, given on the z values of gr_dens_z [# levs/meter]

    ! TODO add descriptions and units
    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      alt_term, &
      Lscale_term, &
      Lscale_term_time_avg, &
      chi_term_time_avg, &
      brunt_term_time_avg, &
      chi_term, &
      brunt_term

    !--------------------- Local Variables ---------------------
    integer :: k, i

    real( kind = core_rknd ), dimension(gr%nzt) :: &
      chi_zt   ! The variable 's' in Mellor (1977)    [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      Lscale,                       & ! Length scale                                 [m]
      ddzt_um,                      & ! Vertical derivative of um                    [s^-1]
      ddzt_vm,                      & ! Vertical derivative of vm                    [s^-1]
      ddzt_umvm_sqd,                & ! Squared vertical norm of derivative of
                                      ! mean horizontal wind speed                   [s^-2]
      brunt_vaisala_freq_sqd,       & ! Buoyancy frequency squared, N^2              [s^-2]
      brunt_vaisala_freq_sqd_zt,    & ! Buoyancy frequency squared (on zt grid), N^2 [s^-2]
      brunt_vaisala_freq_sqd_mixed, & ! A mixture of dry and moist N^2               [s^-2]
      brunt_vaisala_freq_sqd_dry,   & ! dry N^2                                      [s^-2]
      brunt_vaisala_freq_sqd_moist, & ! moist N^2                                    [s^-2]
      brunt_vaisala_freq_sqd_smth     ! Mix between dry and moist N^2 that is
                                      ! smoothed in the vertical                     [s^-2]

    real( kind = core_rknd ), dimension(gr%nzm) :: &
      chi   ! The variable 's' in Mellor (1977)    [kg/kg]

    !--------------------- Begin Code ---------------------

    ! Calculate variables needed for the refinement criterions:
    ! Lscale, chi, ddzt_umvm_sqd, brunt_vaisala_freq_sqd

    ! Calculate Lscale
    Lscale = zt2zm( gr%nzm, gr%nzt, ngrdcol, gr, Lscale_zt )

    ! Calculate ddzt_umvm_sqd
    ddzt_um = ddzt( gr%nzm, gr%nzt, ngrdcol, gr, um )
    ddzt_vm = ddzt( gr%nzm, gr%nzt, ngrdcol, gr, vm )
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        ddzt_umvm_sqd(i,k) = ddzt_um(i,k)**2 + ddzt_vm(i,k)**2
      end do
    end do

    ! Calculate chi
    chi_zt(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
                  + ( one - pdf_params%mixt_frac(1,:) ) * pdf_params%chi_2(1,:)
    chi = zt2zm( gr, chi_zt )

    ! Calculate brunt_vaisala_freq_sqd
    call calc_brunt_vaisala_freq_sqd( gr%nzm, gr%nzt, ngrdcol, gr, thlm,                  & ! In
                                      exner, rtm, rcm, p_in_Pa, thvm,                     & ! In
                                      ice_supersat_frac,                                  & ! In
                                      clubb_config_flags%saturation_formula,              & ! In
                                      clubb_config_flags%l_brunt_vaisala_freq_moist,      & ! In
                                      clubb_config_flags%l_use_thvm_in_bv_freq,           & ! In
                                      clubb_config_flags%l_modify_limiters_for_cnvg_test, & ! In
                                      bv_efold,                                           & ! In
                                      brunt_vaisala_freq_sqd,                             & ! Out
                                      brunt_vaisala_freq_sqd_mixed,                       & ! Out
                                      brunt_vaisala_freq_sqd_dry,                         & ! Out
                                      brunt_vaisala_freq_sqd_moist,                       & ! Out
                                      brunt_vaisala_freq_sqd_smth )                         ! Out

    ! Calculate the grid density from Lscale, chi, ddzt_umvm_sqd and brunt_vaisala_freq_sqd
    call calc_grid_dens_helper( ngrdcol, gr%nzm, gr%zm, &
                                Lscale, wp2, &
                                chi, &
                                brunt_vaisala_freq_sqd, &
                                ddzt_umvm_sqd, &
                                gr_dens_z, gr_dens, &
                                alt_term, Lscale_term, &
                                Lscale_term_time_avg, &
                                chi_term, &
                                chi_term_time_avg, &
                                brunt_term, &
                                brunt_term_time_avg )

  end subroutine calc_grid_dens

  subroutine calc_grid_dens_helper( ngrdcol, nzm, zm, &
                                    Lscale, wp2, &
                                    chi, &
                                    brunt_vaisala_freq_sqd, &
                                    ddzt_umvm_sqd, &
                                    gr_dens_z, gr_dens, &
                                    alt_term, Lscale_term, &
                                    Lscale_term_time_avg, &
                                    chi_term, &
                                    chi_term_time_avg, &
                                    brunt_term, &
                                    brunt_term_time_avg )
    ! Description:
    ! Calculates the non-normalized grid density from Lscale, chi,
    ! ddzt_umvm_sqd and brunt_vaisala_freq_sqd.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, & 
      nzm

    real( kind = core_rknd ), dimension(ngrdcol, nzm), intent(in) ::  &
      zm,                     & ! levels at which the values are given         [m]
      Lscale,                 & ! Length scale                                 [m]
      wp2,                    & ! w'^2 on thermo. grid                         [m^2/s^2]
      brunt_vaisala_freq_sqd, & ! Buoyancy frequency squared, N^2              [s^-2]
      ddzt_umvm_sqd             ! Squared vertical norm of derivative of
                                ! mean horizontal wind speed                   [s^-2]

    real( kind = core_rknd ), dimension(nzm), intent(in) ::  &
      chi ! The variable 's' in Mellor (1977) given on zm levels    [kg/kg]

    !--------------------- Output Variable ---------------------
    ! TODO mazbe remove gr_dens_z, since it should be the same as the zm that was handed in?
    real( kind = core_rknd ), dimension(nzm), intent(out) ::  &
      gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
                       ! grid density function [m]
      gr_dens  ! the values of the connection points of the piecewise linear
               ! grid density function, given on the z values of gr_dens_z [# levs/meter]

    ! TODO add description and units
    real( kind = core_rknd ), dimension(nzm), intent(out) :: &
      alt_term, &
      Lscale_term, &
      Lscale_term_time_avg, &
      chi_term_time_avg, &
      brunt_term_time_avg, &
      chi_term, &
      brunt_term

    !--------------------- Local Variables ---------------------
    integer :: k

    real( kind = core_rknd ) :: &
      wp2_threshold   ! threshold for wp2 whether Lscale is used or not [m^2/s^2]

    real( kind = core_rknd ), dimension(nzm) :: &
      wp2_term

    !--------------------- Begin Code ---------------------
    ! Prepare cumulative_Lscale to use for time averaged Lscale
    if ( .not. allocated( cumulative_Lscale ) ) then
      allocate( cumulative_Lscale(nzm) )
      Lscale_counter = 0
      do k = 1, nzm
        cumulative_Lscale(k) = 0.0_core_rknd
      end do
    end if
    Lscale_counter = Lscale_counter + 1

    if ( .not. allocated( cumulative_brunt ) ) then
      allocate( cumulative_brunt(nzm) )
      brunt_counter = 0
      do k = 1, nzm
        cumulative_brunt(k) = 0.0_core_rknd
      end do
    end if
    brunt_counter = brunt_counter + 1

    if ( .not. allocated( cumulative_chi ) ) then
      allocate( cumulative_chi(nzm) )
      chi_counter = 0
      do k = 1, nzm
        cumulative_chi(k) = 0.0_core_rknd
      end do
    end if
    chi_counter = chi_counter + 1


    wp2_threshold = 0.001

    ! Calculate each refinement criterion term
    do k = 1, nzm
      gr_dens_z(k) = zm(1,k)

      ! Calculate alt_term
      alt_term(k) = 1/(gr_dens_z(k) + 100) ! 100

      ! Calculate Lscale_term
      if ( wp2(1,k) > wp2_threshold ) then
        Lscale_term(k) = 1/(Lscale(1,k)+10.0) ! 0.1/ ...3000
        brunt_term(k) = maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))]) &
                      / ( (maxval([0.0_core_rknd,(ddzt_umvm_sqd(1,k))]) + 1.0e-5) &
                          * (gr_dens_z(k) + 100) ) ! the 20 was 100 before...
      else
        Lscale_term(k) = 0.0_core_rknd
        brunt_term(k) = 0.0_core_rknd
      end if

      ! Calculate chi_term
      chi_term(k) = exp(1000 * chi(k))/(gr_dens_z(k) + 1000.0) ! + 1000, was 100 before

      ! Calculate brunt_term using also ddzt_umvm_sqd
      ! TODO maybe split into more variables to make more readable?
      !brunt_term(k) = maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))]) &
      !                / ( (maxval([0.0_core_rknd,(ddzt_umvm_sqd(1,k))]) + 1.0e-5) &
      !                    * (gr_dens_z(k) + 20) ) ! the 20 was 100 before...

      ! Calculate time averaged Lscale_term, by dividing the cumulative sum by
      ! the number of elements in that sum
      ! (if the grid is adapted this cumulative sum and the counter are set back to 0)
      cumulative_Lscale(k) = cumulative_Lscale(k) + Lscale_term(k)
      cumulative_chi(k) = cumulative_chi(k) + chi_term(k)
      cumulative_brunt(k) = cumulative_brunt(k) + brunt_term(k)

      Lscale_term_time_avg(k) = cumulative_Lscale(k) / Lscale_counter
      chi_term_time_avg(k) = cumulative_chi(k) / chi_counter
      brunt_term_time_avg(k) = cumulative_brunt(k) / brunt_counter
    
    end do
   
    ! Weigh each term to get the final density as sum of the individual terms
    do k = 1, nzm
        ! ref crit
        gr_dens(k) = 0.0 * Lscale_term_time_avg(k) & ! 15
                     + 3.0 * alt_term(k) & ! 1.0
                     + 100.0 * chi_term_time_avg(k) &
                     + 7.5 * brunt_term_time_avg(k) ! 7.5
    end do

    ! Smooth the grid density
    do k = 1, 3
      gr_dens = moving_average( nzm, &
                                gr_dens_z, gr_dens )
    end do

  end subroutine calc_grid_dens_helper

  ! TODO remove function
  !subroutine calc_grid_dens_func( ngrdcol, nzm, zm, &
  !                                gr, &
  !                                gr_dycore, &
  !                                Lscale, wp2, &
  !                                ddzt_umvm_sqd, &
  !                                grid_sfc, grid_top, &
  !                                num_levels, &
  !                                pdf_params, &
  !                                brunt_vaisala_freq_sqd, &
  !                                gr_dens_z, gr_dens, &
  !                                inv_alt_term, lscale_term, &
  !                                lscale_term_time_avg, &
  !                                chi_term, brunt_vaisala_term )
  !  ! Description:
  !  ! Creates an unnormalized grid density function from Lscale
!
  !  ! References:
  !  ! None
  !  !-----------------------------------------------------------------------
!
  !  use clubb_precision, only: &
  !      core_rknd ! Variable(s)
!
  !  use pdf_parameter_module, only:  &
  !      pdf_parameter  ! Type
!
  !  use grid_class, only: &
  !      zt2zm
!
  !  implicit none
!
  !  !--------------------- Input Variables ---------------------
  !  type( grid ) :: &
  !    gr, &
  !    gr_dycore
!
  !  integer, intent(in) :: & 
  !    ngrdcol, & 
  !    nzm, &
  !    num_levels ! the desired number of levels
!
  !  real( kind = core_rknd ), intent(in) ::  &
  !    grid_sfc, &  ! the grids surface; so the first level in the grid
  !                 ! density function has this height [m]
  !    grid_top     ! the grids top; so the last level in the grid
  !                 ! density function has this height [m]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzm), intent(in) ::  &
  !    zm, &      ! levels at which the values are given [m]
  !    Lscale, &  ! Length scale   [m]
  !    ddzt_umvm_sqd, &
  !    wp2     ! w'^2 on thermo. grid [m^2/s^2]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzm), intent(in) ::  &
  !    brunt_vaisala_freq_sqd
!
  !  type(pdf_parameter), intent(in) :: & 
  !    pdf_params     ! PDF parameters
!
  !  !--------------------- Output Variable ---------------------
  !  real( kind = core_rknd ), dimension(nzm), intent(out) ::  &
  !    gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
  !                     ! grid density function [m]
  !    gr_dens  ! the values of the connection points of the piecewise linear
  !             ! grid density function, given on the z values of gr_dens_z [# levs/meter]
!
  !  real( kind = core_rknd ), dimension(nzm), intent(out) :: &
  !    inv_alt_term, &
  !    chi_term, &
  !    brunt_vaisala_term, &
  !    lscale_term, &
  !    lscale_term_time_avg
!
  !  !--------------------- Local Variables ---------------------
  !  real( kind = core_rknd ), dimension(nzm) :: &
  !    lscale_term_dycore, &
  !    chi_zm   ! The variable 's' in Mellor (1977) given on zm levels    [kg/kg]
!
  !  integer :: k, l
!
  !  real( kind = core_rknd ) :: &
  !    threshold, &
  !    inv_alt_norm_factor, &
  !    chi_norm_factor, &
  !    brunt_vaisala_norm_factor
!
  !  real( kind = core_rknd ), dimension(nzm-1) :: &
  !    chi   ! The variable 's' in Mellor (1977)    [kg/kg]
!
  !  real( kind = core_rknd ), dimension(nzm) :: &
  !    Lscale_history_sum
!
  !  !--------------------- Begin Code ---------------------
  !  if ( .not. allocated( cumulative_Lscale ) ) then
  !    allocate( cumulative_Lscale(nzm) )
  !    Lscale_counter = 0
  !    do k = 1, nzm
  !      cumulative_Lscale(k) = 0.0_core_rknd
  !    end do
  !  end if
!
  !  if ( .not. allocated( Lscale_history ) ) then
  !    allocate( Lscale_history(max_history_Lscale,nzm) )
  !    ind_Lscale_history = 1
  !    do l = 1, max_history_Lscale
  !      do k = 1, nzm
  !        Lscale_history(l,k) = 0.0_core_rknd
  !      end do
  !    end do
  !  end if
!
  !  chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
  !                + ( one - pdf_params%mixt_frac(1,:) ) * pdf_params%chi_2(1,:)
!
  !  chi_zm = zt2zm( gr, chi )
!
  !  threshold = 0.001
  !  Lscale_counter = Lscale_counter + 1
!
  !  ! set each refinement criterion term
  !  do k = 1, nzm
  !    gr_dens_z(k) = zm(1,k)
  !    inv_alt_term(k) = 1/(gr_dens_z(k) + 100)
  !    if ( wp2(1,k) > threshold ) then
  !      lscale_term(k) = 1/(Lscale(1,k)+10.0) ! 0.1/ ...3000
  !    else
  !      lscale_term(k) = 0.0_core_rknd
  !    end if
!
  !    cumulative_Lscale(k) = cumulative_Lscale(k) + lscale_term(k)
  !    lscale_term_time_avg(k) = cumulative_Lscale(k) / Lscale_counter
!
  !    ! TODO remap to common grid before writing to array
  !    !Lscale_history(mod(ind_Lscale_history,max_history_Lscale)+1, k) = lscale_term(k)
  !    !Lscale_history_sum = sum(Lscale_history,dim=1)
  !    !lscale_term_time_avg(k) = Lscale_history_sum(k) / max_history_Lscale
!
  !    chi_term(k) = exp(1000 * chi_zm(k))/(gr_dens_z(k) + 1000.0) ! + 1000, was 100 before
  !    brunt_vaisala_term(k) = maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))]) &
  !                            / ( (maxval([0.0_core_rknd,(ddzt_umvm_sqd(1,k))]) + 1.0e-5) &
  !                                * (gr_dens_z(k) + 20) ) ! the 20 was 100 before...
!
  !  end do
!
  !  ! build time average for lscale_term
  !  ! TODO instead of map1_ppm use remap_vals_to_target function
  !  
  !  !call map1_ppm( nzm-1,   zm,   lscale_term,  &
  !  !               gr_dycore%nzm-1,   gr_dycore%zm,   lscale_term_dycore, &
  !  !                     0, 0, 1, 1, 1,                         &
  !  !                     1, 1, 1, 0, 3)
  !  !do k = 1, nzm
  !  !  Lscale_history(mod(ind_Lscale_history,max_history_Lscale)+1, k) = lscale_term_dycore(k)
  !  !end do
!!
  !  !ind_Lscale_history = ind_Lscale_history + 1
  !  !!Lscale_history_sum = sum(Lscale_history,dim=1)
  !  !call map1_ppm( gr_dycore%nzm-1,   gr_dycore%zm,   sum(Lscale_history,dim=1),  &
  !  !               nzm-1,   zm,   Lscale_history_sum, &
  !  !                     0, 0, 1, 1, 1,                         &
  !  !                     1, 1, 1, 0, 3)
!!
  !  !do k = 1, nzm
  !  !  lscale_term_time_avg(k) = Lscale_history_sum(k) / max_history_Lscale
  !  !end do
!
  !  ! calculate noramlization factors, so each refinement criterion has same magnitude
  !  !inv_alt_norm_factor = 1/calc_integral( nzm, zm(1,:), inv_alt_term )
  !  !chi_norm_factor = 1/calc_integral( nzm, zm(1,:), chi_term )
  !  !brunt_vaisala_norm_factor = 1/calc_integral( nzm, zm(1,:), brunt_vaisala_term )
!
  !  do k = 1, nzm
  !  
  !      !gr_dens(k) = 0.2 * inv_alt_norm_factor * inv_alt_term(k) &              ! 0.2
  !      !             + 0.5 * chi_norm_factor * chi_term(k) &                    ! 0.5
  !      !             + 0.3 * brunt_vaisala_norm_factor * brunt_vaisala_term(k)  ! 0.3
!
  !      ! lololo
  !      !gr_dens(k) = 5.0 * lscale_term(k) &
  !      !             + 3.0 * inv_alt_term(k) & 
  !      !             + 40.0 * chi_term(k) &
  !      !             + 1.0 * brunt_vaisala_term(k)
!
  !      !do l = 1, 1
  !      !  lscale_term = moving_average( nzm, &
  !      !                                gr_dens_z, lscale_term )
  !      !end do
!
  !      ! !!!!!!! those coefficients work good !!!!
  !      !gr_dens(k) = 15.0 * lscale_term_time_avg(k) &
  !      !             + 1.0 * inv_alt_term(k) & 
  !      !             + 100.0 * chi_term(k) &
  !      !             + 0.0 * brunt_vaisala_term(k)
!
  !      gr_dens(k) = 15.0 * lscale_term_time_avg(k) &
  !                   + 1.0 * inv_alt_term(k) &
  !                   + 100.0 * chi_term(k) &
  !                   + 0.0 * brunt_vaisala_term(k)
!
  !      !gr_dens(k) = lscale_term_time_avg(k)
!
  !      ! cond0: only adapt every 80 iteration -> reduces noise and shows some good results, but cond4 looks better
  !      ! cond1: brunt_vaisala/2 and ddzt_umvm*2 adapt every 80th iteration -> results were worse than before
  !      ! cond2: ddzt_umvm*2 adapt every 80th iteration -> results were worse than before
  !      ! cond3: gr_dens(k) = gr_dens(k) + 0.5/(gr_dens_z(k)+1.0) adapt every 80th iteration -> results were a bit worse than before but noise was reduced for lwp, but more noise in wp2
  !      ! cond4: gr_dens(k) = gr_dens(k) + 0.05/(gr_dens_z(k)+1.0) adapt every 120th iteration -> !!!
  !  end do
!
  !  if (gr_dens_z(1) > grid_sfc) then
  !      gr_dens_z(1) = grid_sfc
  !  end if
!
  !  if (gr_dens_z(nzm) < grid_top) then
  !      gr_dens_z(nzm) = grid_top
  !  end if
!
  !  do k = 1, 3
  !    gr_dens = moving_average( nzm, &
  !                              gr_dens_z, gr_dens )
  !  end do
!
  !end subroutine calc_grid_dens_func
!
  !! gabls2 case v1
  !!subroutine calc_grid_dens_func( ngrdcol, nzm, zm, &
  !                                gr, &
  !                                Lscale, wp2, &
  !                                ddzt_umvm_sqd, &
  !                                grid_sfc, grid_top, &
  !                                num_levels, &
  !                                pdf_params, &
  !                                brunt_vaisala_freq_sqd, &
  !                                gr_dens_z, gr_dens )
  !  ! Description:
  !  ! Creates an unnormalized grid density function from Lscale
!
  !  ! References:
  !  ! None
  !  !-----------------------------------------------------------------------
!
  !  use clubb_precision, only: &
  !      core_rknd ! Variable(s)
!
  !  use pdf_parameter_module, only:  &
  !      pdf_parameter  ! Type
!
  !  use grid_class, only: &
  !      zt2zm
!
  !  implicit none
!
  !  !--------------------- Input Variables ---------------------
  !  type( grid ) :: &
  !    gr
!
  !  integer, intent(in) :: & 
  !    ngrdcol, & 
  !    nzm, &
  !    num_levels ! the desired number of levels
!
  !  real( kind = core_rknd ), intent(in) ::  &
  !    grid_sfc, &  ! the grids surface; so the first level in the grid
  !                 ! density function has this height [m]
  !    grid_top     ! the grids top; so the last level in the grid
  !                 ! density function has this height [m]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzm), intent(in) ::  &
  !    zm, &      ! levels at which the values are given [m]
  !    Lscale, &  ! Length scale   [m]
  !    ddzt_umvm_sqd, &
  !    wp2     ! w'^2 on thermo. grid [m^2/s^2]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzm), intent(in) ::  &
  !    brunt_vaisala_freq_sqd
!
  !  type(pdf_parameter), intent(in) :: & 
  !    pdf_params     ! PDF parameters
!
  !  !--------------------- Output Variable ---------------------
  !  real( kind = core_rknd ), dimension(nzm), intent(out) ::  &
  !    gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
  !                     ! grid density function [m]
  !    gr_dens  ! the values of the connection points of the piecewise linear
  !             ! grid density function, given on the z values of gr_dens_z [# levs/meter]
!
  !  !--------------------- Local Variables ---------------------
  !  integer :: k
!
  !  real( kind = core_rknd ) :: threshold
!
  !  real( kind = core_rknd ), dimension(nzm-1) :: &
  !    chi   ! The variable 's' in Mellor (1977)    [kg/kg]
!
  !  real( kind = core_rknd ), dimension(nzm) :: &
  !    chi_zm   ! The variable 's' in Mellor (1977) given on zm levels    [kg/kg]
!
  !  !--------------------- Begin Code ---------------------
  !  chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
  !                + ( one - pdf_params%mixt_frac(1,:) ) * pdf_params%chi_2(1,:)
!
  !  chi_zm = zt2zm( gr, chi )
!
  !  threshold = 0.001
  !  !threshold = 0.1
  !  do k = 1, nzm
  !      gr_dens_z(k) = zm(1,k)
  !      if ( wp2(1,k) > threshold ) then
  !        !gr_dens(k) = 1.0/(Lscale(1,k)+10.0) + 3.0/(gr_dens_z(k)+20.0)
  !        !gr_dens(k) = 1.0/(Lscale(1,k)+10.0)! + 10.0/(gr_dens_z(k)+1.0) ! +20 before ! if we just use this condition, then at least some of the resuts in some of the timeframes show better results (e.g. wp3 in 0-500 or 1000-2000)
  !        gr_dens(k) = 1.0/(Lscale(1,k)+10.0)
  !      
  !      else
  !        gr_dens(k) = 0.0_core_rknd*(num_levels - 1)/(grid_top - grid_sfc) ! 0.1 before
  !      end if
  !      gr_dens(k) = gr_dens(k) + 2.0/1000.0 * exp(100 * chi_zm(k))
!
  !      gr_dens(k) = gr_dens(k) &
  !                   + 1.0e2_core_rknd*maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))])
  !      gr_dens(k) = gr_dens(k) + 4.0/(gr_dens_z(k)+1.0)
  !      gr_dens(k) = gr_dens(k) &
  !                   + 2.0e1_core_rknd*maxval([0.0_core_rknd,(ddzt_umvm_sqd(1,k))])
  !      !gr_dens(k) = gr_dens(k) + 4.0/(gr_dens_z(k)+1.0)
!
!
!
  !      !gr_dens(k) = gr_dens(k) &
  !      !             + 2.0e5_core_rknd*maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))**2])
  !  end do
!
  !  !do k = 1, nzm
  !  !    gr_dens_z(k) = zm(1,k)
  !  !    gr_dens(k) = 1.0/10.0 * exp(100 * chi_zm(k))
  !  !end do
!
  !  !threshold = 0.001
  !  !do i = 1, nzm
  !  !    gr_dens_z(i) = zm(1,i)f
  !  !    if ( wp2(1,i) > threshold ) then
  !  !      gr_dens(i) = 1/Lscale(1,i)
  !  !    else
  !  !      gr_dens(i) = 0.75*(num_levels - 1)/(grid_top - grid_sfc)
  !  !    end if
  !  !end do
!
  !  ! TODO find better way to ensure grid density function goes to zm(1) and zm(n)
  !  ! TODO can be removed if I just use values on zm levels (just use zt2zm...)
  !  if (gr_dens_z(1) > grid_sfc) then
  !      gr_dens_z(1) = grid_sfc
  !  end if
!
  !  if (gr_dens_z(nzm) < grid_top) then
  !      gr_dens_z(nzm) = grid_top
  !  end if
!
  !  do k = 1, 3
  !  gr_dens = moving_average( nzm, &
  !                            gr_dens_z, gr_dens )
  !  end do
  !      
  !  
!
  !end subroutine calc_grid_dens_func

  ! optimized for astex case
  !subroutine calc_grid_dens_func_old( ngrdcol, nzt, zt, &
  !                                    Lscale, wp2_zt, &
  !                                    grid_sfc, grid_top, &
  !                                    num_levels, &
  !                                    pdf_params, &
  !                                    brunt_vaisala_freq_sqd_zt, &
  !                                    gr_dens_z, gr_dens )
  !  ! Description:
  !  ! Creates an unnormalized grid density function from Lscale
!
  !  ! References:
  !  ! None
  !  !-----------------------------------------------------------------------
!
  !  use clubb_precision, only: &
  !      core_rknd ! Variable(s)
!
  !  use pdf_parameter_module, only:  &
  !      pdf_parameter  ! Type
!
  !  implicit none
!
  !  !--------------------- Input Variables ---------------------
  !  integer, intent(in) :: & 
  !    ngrdcol, & 
  !    nzt, &
  !    num_levels ! the desired number of levels
!
  !  real( kind = core_rknd ), intent(in) ::  &
  !    grid_sfc, &  ! the grids surface; so the first level in the grid
  !                 ! density function has this height [m]
  !    grid_top     ! the grids top; so the last level in the grid
  !                 ! density function has this height [m]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzt), intent(in) ::  &
  !    zt, &      ! levels at which the values are given [m]
  !    Lscale, &  ! Length scale   [m]
  !    wp2_zt     ! w'^2 on thermo. grid [m^2/s^2]
!
  !  real( kind = core_rknd ), dimension(ngrdcol, nzt), intent(in) ::  &
  !    brunt_vaisala_freq_sqd_zt
!
  !  type(pdf_parameter), intent(in) :: & 
  !    pdf_params     ! PDF parameters
!
  !  !--------------------- Output Variable ---------------------
  !  real( kind = core_rknd ), dimension(nzt), intent(out) ::  &
  !    gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
  !                     ! grid density function [m]
  !    gr_dens  ! the values of the connection points of the piecewise linear
  !             ! grid density function, given on the z values of gr_dens_z [# levs/meter]
!
  !  !--------------------- Local Variables ---------------------
  !  integer :: k
!
  !  real( kind = core_rknd ) :: threshold
!
  !  real( kind = core_rknd ), dimension(nzt) :: &
  !    chi   ! The variable 's' in Mellor (1977)    [kg/kg]
!
  !  !--------------------- Begin Code ---------------------
  !  !lambda = 0.01
  !  !lambda = 0.5
  !  chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
  !                + ( one - pdf_params%mixt_frac(1,:) ) * pdf_params%chi_2(1,:)
!
  !  threshold = 0.001
  !  do k = 1, nzt
  !      gr_dens_z(k) = zt(1,k)
  !      if ( wp2_zt(1,k) > threshold ) then
  !        !gr_dens(k) = 1.0/(Lscale(1,k)+10.0) + 3.0/(gr_dens_z(k)+20.0)
  !        gr_dens(k) = 1.0/(Lscale(1,k)+10.0) + 4.0/(gr_dens_z(k)+1.0) ! +20 before
  !      else
  !        gr_dens(k) = 0.0_core_rknd*(num_levels - 1)/(grid_top - grid_sfc) ! 0.1 before
  !      end if
  !      gr_dens(k) = gr_dens(k) + 2.0/1000.0 * exp(100 * chi(k))
  !      !gr_dens(k) = gr_dens(k) &
  !      !             + 2.0e5_core_rknd*maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd_zt(1,k))**2])
  !      gr_dens(k) = gr_dens(k) &
  !                   + 1.0e2_core_rknd*1.0/3.0e-4*maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd_zt(1,k))**2])
  !  end do
!
  !  !do k = 1, nzt
  !  !    gr_dens_z(k) = zt(1,k)
  !  !    gr_dens(k) = 1.0/10.0 * exp(100 * chi(k))
  !  !end do
!
  !  !threshold = 0.001
  !  !do i = 1, nzt
  !  !    gr_dens_z(i) = zt(1,i)
  !  !    if ( wp2_zt(1,i) > threshold ) then
  !  !      gr_dens(i) = 1/Lscale(1,i)
  !  !    else
  !  !      gr_dens(i) = 0.75*(num_levels - 1)/(grid_top - grid_sfc)
  !  !    end if
  !  !end do
!
  !  ! TODO find better way to ensure grid density function goes to zm(1) and zm(n)
  !  if (gr_dens_z(1) > grid_sfc) then
  !      gr_dens_z(1) = grid_sfc
  !  end if
!
  !  if (gr_dens_z(nzt) < grid_top) then
  !      gr_dens_z(nzt) = grid_top
  !  end if
!
  !  gr_dens = moving_average( nzt, &
  !                            gr_dens_z, gr_dens )
!
  !end subroutine calc_grid_dens_func_old

  subroutine clean_up_grid_adaptation_module()

    implicit none

    if ( allocated( fixed_min_gr_dens ) ) then
      deallocate( fixed_min_gr_dens )
    end if

    if ( allocated( fixed_min_gr_dens_z ) ) then
      deallocate( fixed_min_gr_dens_z )
    end if

    if ( allocated( gr_dens_old_global ) ) then
      deallocate( gr_dens_old_global )
    end if

    if ( allocated( gr_dens_old_z_global ) ) then
      deallocate( gr_dens_old_z_global )
    end if

    if ( allocated( cumulative_Lscale ) ) then
      deallocate( cumulative_Lscale )
    end if

    if ( allocated( cumulative_brunt ) ) then
      deallocate( cumulative_brunt )
    end if

    if ( allocated( cumulative_chi ) ) then
      deallocate( cumulative_chi )
    end if

  end subroutine

  function moving_average( nz, &
                           gr_dens_z, gr_dens )

    ! Description:
    ! Applies a moving average filter to the given function

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      gr_dens_z, & ! the z value coordinates of the connection points of the piecewise linear
                   ! grid density function [m]
      gr_dens  ! the values of the connection points of the piecewise linear
               ! grid density function, given on the z values of gr_dens_z [# levs/meter]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(nz) ::  &
      moving_average ! the smoothed values of the connection points of the piecewise linear
                     ! grid density function, given on the z values of gr_dens_z [# levs/meter]

    !--------------------- Local Variables ---------------------
    integer :: &
      window_size, &
      start_index_full_window, &
      end_index_full_window, &
      k, j

    real( kind = core_rknd ) :: &
      avg, &
      sum_weights, &
      weight_for_middle_point ! a weight for the point in the middle

    real( kind = core_rknd ), dimension(:), allocatable :: &
      weights ! the weights for one filter, changes for each point

    !--------------------- Begin Code ---------------------
    window_size = 3 ! needs to be uneven, so on both sides equally many points are
                    ! taking into consideration

    allocate( weights(window_size) )

    weight_for_middle_point = 0.5 ! needs to be betweeen [0,1]

    
    start_index_full_window = (window_size-1)/2 + 1 ! is how many points on each side are
                                                    ! taking into consideration

    end_index_full_window = nz - (window_size-1)/2  ! is how many points on each side are
                                                    ! taking into consideration

    ! just copy values where no full window is available
    do k = 1, start_index_full_window - 1
      moving_average(k) = gr_dens(k)
    end do

    do k = end_index_full_window+1, nz
      moving_average(k) = gr_dens(k)
    end do

    do k = start_index_full_window, end_index_full_window
      
      do j = k-(window_size-1)/2, k+(window_size-1)/2

        weights(j-(k-(window_size-1)/2)+1) = 0.0_core_rknd
        if ( j /= k ) then
          ! so we are not on the middle point of the current filter
          weights(j-(k-(window_size-1)/2)+1) = 1/abs( gr_dens_z(k) - gr_dens_z(j) )
        end if

      end do

      ! normalize weights
      sum_weights = sum( weights )
      do j = k-(window_size-1)/2, k+(window_size-1)/2
        weights(j-(k-(window_size-1)/2)+1) = (1-weight_for_middle_point)/sum_weights &
                                             * weights(j-(k-(window_size-1)/2)+1)
      end do

      weights((window_size-1)/2 + 1) = weight_for_middle_point
      
      if ( (sum(weights) - one) > tol ) then
        error stop 'Sum of weights in moving_average must be 1.0.'
      end if

      ! apply weighted moving average filter
      avg = 0.0_core_rknd
      do j = k-(window_size-1)/2, k+(window_size-1)/2
        avg = avg + weights(j-(k-(window_size-1)/2)+1) * gr_dens(j)
      end do
      moving_average(k) = avg
    end do

    deallocate( weights )

  end function

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm(dm, a4, itot, lmt)

! !USES:
 implicit none

! !INPUT PARAMETERS:
      real(core_rknd), intent(in)::     dm(*)     ! ??????
      integer, intent(in) ::     itot      ! Total Longitudes
      integer, intent(in) ::     lmt       ! 0: Standard PPM constraint
                                           ! 1: Improved full monotonicity constraint (Lin)
                                           ! 2: Positive definite constraint
                                           ! 3: do nothing (return immediately)

! !INPUT/OUTPUT PARAMETERS:
      real(core_rknd), intent(inout) :: a4(4,*)   ! ???????
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
!    Writes a standard set of data to the history buffer. 
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real(core_rknd)       r12
      parameter (r12 = D1_0/D12_0)

      real(core_rknd) qmp
      integer i
      real(core_rknd) da1, da2, a6da
      real(core_rknd) fmin

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == D0_0) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = D0_0
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = D3_0*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = D3_0*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = D2_0*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = D3_0*( D2_0*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+D0_25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < D0_0 ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = D0_0
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = D3_0*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = D3_0*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

      return
!EOC
 end subroutine kmppm
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

! !USES:
   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real(core_rknd), intent(in) ::  dp(i1:i2,km)       ! grid size
      real(core_rknd), intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real(core_rknd), intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real(core_rknd), intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real(core_rknd), intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real(core_rknd), intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real(core_rknd) alfa(i1:i2,km)
      real(core_rknd)    f(i1:i2,km)
      real(core_rknd)  rat(i1:i2,km)
      real(core_rknd)  dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) = (rat(i,k+1) - rat(i,k))                             &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<D0_0 .and. df2(i,k)/=D0_0) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(D0_0, min(D0_5, -D0_1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = D0_0
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (D1_0-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

      return
!EOC
 end subroutine steepz
!----------------------------------------------------------------------- 

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! !USES:
 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real (core_rknd), intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real (core_rknd), intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
!   02.04.23    Sawyer     Incorporated minor algorithmic change to 
!                          maintain CAM zero diffs (see comments inline)
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real(core_rknd)  dc(i1:i2,km)
      real(core_rknd)  h2(i1:i2,km)
      real(core_rknd) delq(i1:i2,km)
      real(core_rknd) df2(i1:i2,km)
      real(core_rknd) d4(i1:i2,km)

! local scalars:
      real(core_rknd) fac
      real(core_rknd) a1, a2, c1, c2, c3, d1, d2
      real(core_rknd) qmax, qmin, cmax, cmin
      real(core_rknd) qm, dq, tmp

      integer i, k, km1, lmt
      real(core_rknd) qmp, pmp
      real(core_rknd) lac
      integer it

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
            c1  = (delp(i,k-1)+D0_5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+D0_5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
            dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
      a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
      a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + D2_0/(d4(i,k-1)+d4(i,k+1)) *    &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = D2_0*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = D4_0*(a4(2,i,3)-qm-d2*dq) / ( d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)) )
      c3 = dq - D0_5*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
      a4(2,i,2) = qm - D0_25*c1*d1*d2*(d2+D3_0*d1)
      a4(2,i,1) = d1*(D2_0*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

      if( iv == 0 ) then
         do i=i1,i2
!
! WS: 02.04.23  Algorithmic difference with FVGCM.  FVGCM does this:
!
!!!            a4(2,i,1) = a4(1,i,1)
!!!            a4(3,i,1) = a4(1,i,1)
!
!     CAM does this:
!
            a4(2,i,1) = max(D0_0,a4(2,i,1))
            a4(2,i,2) = max(D0_0,a4(2,i,2))
         enddo
      elseif ( iv == -1 ) then
! Winds:
        if( km > 32 ) then
          do i=i1,i2
! More dampping: top layer as the sponge
             a4(2,i,1) = a4(1,i,1)
             a4(3,i,1) = a4(1,i,1)
          enddo
        else
          do i=i1,i2
             if( a4(1,i,1)*a4(2,i,1) <=  D0_0 ) then
                 a4(2,i,1) = D0_0
             else
                 a4(2,i,1) = sign(min(abs(a4(1,i,1)),    &
                                      abs(a4(2,i,1))),   &
                                          a4(1,i,1)  )
            endif
          enddo
        endif
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = D2_0*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)))
         c3 = dq - D2_0*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+D3_0*d1)
         a4(3,i,km) = d1*(D8_0*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

! Enforce constraint at the surface

      if ( iv == 0 ) then
! Positive definite scalars:
           do i=i1,i2
              a4(3,i,km) = max(D0_0, a4(3,i,km))
           enddo
      elseif ( iv == -1 ) then
! Winds:
           do i=i1,i2
              if( a4(1,i,km)*a4(3,i,km) <=  D0_0 ) then
                  a4(3,i,km) = D0_0
              else
                  a4(3,i,km) = sign( min(abs(a4(1,i,km)),   &
                                         abs(a4(3,i,km))),  &
                                             a4(1,i,km)  )
              endif
           enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = D2_0*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+D0_5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = D1_5           ! original quasi-monotone
      else
         fac = D0_125         ! full monotone
      endif

      do k=3, km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + D2_0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + D0_5*delq(i,k-1)
!
         pmp   = D2_0*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - D2_0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - D0_5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=3, km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = D3_0*(D2_0*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
!EOC
 end subroutine ppm2m
!-----------------------------------------------------------------------

  !----------------------------------------------------------------------- 
  !BOP
  ! !ROUTINE:  map1_ppm --- Piecewise parabolic mapping, variant 1
  !
  ! !INTERFACE:
  subroutine map1_ppm( km, pe1, q1, kn, pe2, q2, &
                       ng_s, ng_n, itot, i1, i2, &
                       j, jfirst, jlast, iv, kord )

        implicit none

  ! !INPUT PARAMETERS:
        integer, intent(in) :: i1                ! Starting longitude
        integer, intent(in) :: i2                ! Finishing longitude
        integer, intent(in) :: itot              ! Total latitudes
        integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
        integer, intent(in) :: kord              ! Method order
        integer, intent(in) :: j                 ! Current latitude
        integer, intent(in) :: jfirst            ! Starting latitude
        integer, intent(in) :: jlast             ! Finishing latitude
        integer, intent(in) :: ng_s              ! Ghosted latitudes south
        integer, intent(in) :: ng_n              ! Ghosted latitudes north
        integer, intent(in) :: km                ! Original vertical dimension
        integer, intent(in) :: kn                ! Target vertical dimension

        real( core_rknd ), intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                                 ! (from model top to bottom surface)
                                                 ! in the original vertical coordinate
        real( core_rknd ), intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                                 ! (from model top to bottom surface)
                                                 ! in the new vertical coordinate
        real( core_rknd ), intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km) ! Field input

  ! !INPUT/OUTPUT PARAMETERS:
        real( core_rknd ), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

  ! !DESCRIPTION:
  !
  !     Perform piecewise parabolic method on a given latitude    
  ! IV = 0: constituents
  ! pe1: pressure at layer edges (from model top to bottom surface)
  !      in the original vertical coordinate
  ! pe2: pressure at layer edges (from model top to bottom surface)
  !      in the new vertical coordinate
  !
  ! !REVISION HISTORY: 
  !    00.04.24   Lin       Last modification
  !    01.03.26   Sawyer    Added ProTeX documentation
  !    02.04.04   Sawyer    incorporated latest FVGCM version
  !    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
  !    03.07.22   Parks     Cleaned main loop, removed gotos
  !    05.05.25   Sawyer    Merged CAM and GEOS5 versions
  !
  !EOP
  !-----------------------------------------------------------------------
  !BOC
  !
  ! !LOCAL VARIABLES:
        real(core_rknd)       r3, r23
        parameter (r3 = D1_0/D3_0, r23 = D2_0/D3_0)
        real(core_rknd)   dp1(i1:i2,km)
        real(core_rknd)  q4(4,i1:i2,km)

        integer i, k, kk, kl, k0(i1:i2,0:kn+1), k0found
        real(core_rknd)    pl, pr, qsum, qsumk(i1:i2,kn), delp, esl

        do k=1,km
           do i=i1,i2
              dp1(i,k) = pe1(i,k+1) - pe1(i,k)
              q4(1,i,k) = q1(i,j,k)
           enddo
        enddo

  ! Mapping
  ! Compute vertical subgrid distribution
        call ppm2m( q4, dp1, km, i1, i2, iv, kord )

  ! For each pe2(i,k), determine lowest pe1 interval = smallest k0 (= k0(i,k))
  !   such that pe1(i,k0) <= pe2(i,k) <= pe1(i,k0+1)
  !   Note that pe2(i,1)==pe1(i,1) and pe2(i,kn+1)==pe1(i,kn+1)
  !   Note also that pe1, pe2 are assumed to be monotonically increasing
  !#if defined( UNICOSMP ) || defined ( NEC_SX )
  !      do kk = km, 1, -1
  !         do k = 1, kn+1
  !!dir$ prefervector
  !            do i = i1, i2
  !               if (pe2(i,k) <= pe1(i,kk+1)) then
  !                  k0(i,k) = kk
  !                  write(*,*) 'kk: ', kk
  !               endif
  !            enddo
  !         enddo
  !      enddo
  !#else
        do i = i1, i2
           k0(i,0) = 1
           do k = 1, kn+1
              k0found = -1
              do kk = k0(i,k-1), km
                 if (pe2(i,k) <= pe1(i,kk+1)) then
                    k0(i,k) = kk
                    k0found = kk
                    exit
                 endif
              enddo
              if (k0found .lt. 0) then
                 write(fstderr,*) 'mapz error - k0found i j k (kk,pe1,pe2) = ',   &
                    k0found, i, j, k, (kk,pe1(i,kk),pe2(i,kk),kk=1,km+1)
                 !call endrun('MAPZ_MODULE')
                 return
              endif
           enddo
        enddo
  !#endif

  ! Interpolate
        do k = 1, kn

  ! Prepare contribution between pe1(i,ko(i,k)+1) and pe1(i,k0(i,k+1))
           qsumk(:,k) = D0_0
           do i = i1, i2
              do kl = k0(i,k)+1, k0(i,k+1)-1
                 qsumk(i,k) = qsumk(i,k) + dp1(i,kl)*q4(1,i,kl)
              enddo
           enddo

           do i = i1, i2
              kk = k0(i,k)
  ! Consider contribution between pe1(i,kk) and pe2(i,k)
              pl = (pe2(i,k)-pe1(i,kk)) / dp1(i,kk)
  ! Check to see if pe2(i,k+1) and pe2(i,k) are in same pe1 interval
              if (k0(i,k+1) == k0(i,k)) then
                 pr = (pe2(i,k+1)-pe1(i,kk)) / dp1(i,kk)
                 q2(i,j,k) = q4(2,i,kk) + D0_5*(q4(4,i,kk)+q4(3,i,kk)-q4(2,i,kk))  &
                    *(pr+pl) - q4(4,i,kk)*r3*(pr*(pr+pl)+pl**2)
              else
  ! Consider contribution between pe2(i,k) and pe1(i,kk+1)
                 qsum = (pe1(i,kk+1)-pe2(i,k))*(q4(2,i,kk)+D0_5*(q4(4,i,kk)+       &
                    q4(3,i,kk)-q4(2,i,kk))*(D1_0+pl)-q4(4,i,kk)*                    &
                    (r3*(D1_0+pl*(D1_0+pl))))
  ! Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
                 qsum = qsum + qsumk(i,k)
  ! Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
                 kl = k0(i,k+1)
                 delp = pe2(i,k+1)-pe1(i,kl)
                 esl = delp / dp1(i,kl)
                 qsum = qsum + delp*(q4(2,i,kl)+D0_5*esl*                          &
                    (q4(3,i,kl)-q4(2,i,kl)+q4(4,i,kl)*(D1_0-r23*esl)))
                 q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
              endif
           enddo
        enddo

        return
  !EOC
  end subroutine map1_ppm
!----------------------------------------------------------------------- 







!===============================================================================

end module grid_adaptation_module