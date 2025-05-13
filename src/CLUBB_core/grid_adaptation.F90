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
      clubb_at_least_debug_level_api

  implicit none

  real( kind = core_rknd ) ::  &
        tol = 1.0e-6_core_rknd ! tolerance to check if to real numbers are equal

  public :: setup_simple_gr_dycore, &
            check_remap_conservation, check_remap_consistency, check_remap_monotonicity, &
            remapping_matrix, remap_vals_to_target, &
            remap_zt_values, remap_zm_values, remap_forcings, &
            check_remap_for_consistency, check_remap_for_consistency_zt_zm, &
            check_mass_conservation, check_mass_conservation_zt_zm, &
            check_vertical_integral_conservation, &
            check_vertical_integral_conservation_all_zt_values, &
            check_vertical_integral_conservation_all_zm_values, &
            calc_mass_over_grid_intervals, vertical_integral_conserve_mass, &
            adapt_grid, calc_grid_dens_func

  private :: calc_integral, normalize_grid_density, &
             create_grid_from_normalized_grid_density_func, &
             check_grid, create_grid_from_grid_density_func

  private ! Default Scoping

  contains

  !=============================================================================
  subroutine setup_simple_gr_dycore( nlevels, grid_sfc, grid_top, gr, err_info )

    use constants_clubb, only : &
        fstderr     ! Fortran file unit I/O constant

    use clubb_api_module, only: &
      setup_grid_api

    use error_code, only: &
      clubb_at_least_debug_level_api,  & ! Procedure
      clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: nlevels  ! number of levels the dycore grid should have []

    real( kind = core_rknd ) :: &
        grid_sfc, &    ! grids surface height for the dycore grid [m]
        grid_top       ! grids highest level for the dycore grid [m]

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- In/Output Variables ---------------------
    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Local Variables ---------------------
    integer :: i
    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer :: grid_type

    real( kind = core_rknd ) ::  &
        sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical :: l_implemented

    real( kind = core_rknd ) ::  & 
        deltaz,   & ! Vertical grid spacing                  [m]
        zm_init,  & ! Initial grid altitude (momentum level) [m]
        zm_top      ! Maximum grid altitude (momentum level) [m]

    real( kind = core_rknd ), dimension(nlevels) :: & 
      momentum_heights      ! Momentum level altitudes (file input)      [m]
    
    real( kind = core_rknd ), dimension(nlevels-1) :: &
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

    !--------------------- Begin Code ---------------------

    sfc_elevation = 0.0
    l_implemented = .false.
    grid_type = 3
    deltaz = (grid_top-grid_sfc)/(nlevels-1)
    zm_init = grid_sfc
    zm_top = grid_top+1 ! make sure the highest level in momentum_heights is included in the grid

    do i = 1, nlevels-1
        momentum_heights(i) = grid_sfc + (i-1)*deltaz
    end do

    momentum_heights(nlevels) = grid_top

    call setup_grid_api( nlevels, sfc_elevation, l_implemented, &       ! intent(in)
                         .true., grid_type, deltaz, zm_init, zm_top, &  ! intent(in)
                         momentum_heights, thermodynamic_heights, &     ! intent(in)
                         gr, err_info )                                 ! intent(inout)

    if ( clubb_at_least_debug_level_api(0) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid_api in setup_simple_gr_dycore"
      end if
    end if

  end subroutine setup_simple_gr_dycore

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

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel_target-1,nlevel_source-1) :: &
      remapping_matrix ! matrix to apply to input values for remapping (R_ij) []

    integer :: i, j

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
    do i = 1, (nlevel_target-1)
        do j = 1, (nlevel_source-1)
            omega_ov_upper = min(levels_target(i+1), levels_source(j+1))
            omega_ov_lower = max(levels_target(i), levels_source(j))
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
            remapping_matrix(i,j) = omega_ov/omega_ts(i)
        end do
    end do

    if ( clubb_at_least_debug_level_api( 2 ) ) then
      call check_remap_conservation( remapping_matrix, (nlevel_target-1), (nlevel_source-1), &
                                     levels_source, levels_target, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels)
      call check_remap_consistency( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
      call check_remap_monotonicity( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
    end if

  end function remapping_matrix

  function remap_vals_to_target( nlevel_source, nlevel_target, &
                                 levels_source, levels_target, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 gr_source_values )

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

    real( kind = core_rknd ), dimension(nlevel_source-1) :: &
      gr_source_values  ! given values on the gr_source grid that should be interpolated 
                        ! to the gr_target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      remap_vals_to_target ! interpolated values with the dimension of gr_target

    integer :: i, j

    real( kind = core_rknd ), dimension(nlevel_target-1,nlevel_source-1) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Begin Code ---------------------

    R_ij = remapping_matrix( nlevel_source, nlevel_target, &
                             levels_source, levels_target, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels )
    ! matrix vector multiplication
    do i = 1, (nlevel_target-1)
      remap_vals_to_target(i) = 0
      do j = 1, (nlevel_source-1)
        remap_vals_to_target(i) = remap_vals_to_target(i) + R_ij(i,j)*gr_source_values(j)
      end do
    end do

  end function remap_vals_to_target

  function remap_zt_values( ngrdcol, gr_source, gr_target, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            values_source )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzt), intent(in) :: &
      values_source ! values given on the source grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_target%nzt) :: &
      remap_zt_values ! interpolated values on target grid

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    do i=1, ngrdcol
        remap_zt_values(i,:) = remap_vals_to_target( gr_source%nzm, gr_target%nzm, &
                                                     gr_source%zm(i,:), gr_target%zm(i,:), &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals(i,:), &
                                                     rho_lin_spline_levels(i,:), &
                                                     values_source(i,:) )
    end do
  end function remap_zt_values

  function remap_zm_values( ngrdcol, gr_source, gr_target, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            values_source )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzm), intent(in) :: &
      values_source ! values given on the source grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_target%nzm) :: &
      remap_zm_values ! interpolated values on target grid

    !--------------------- Local Variables ---------------------
    integer :: i, j

    integer :: &
      nlevels_source ! number of levels in source grid []

    integer :: &
      nlevels_target ! number of levels in target grid []

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! levels on source grid [m]

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! levels on target grid [m]

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
    
    do i=1, ngrdcol

        levels_source(1) = gr_source%zm(i,1)
        do j=1, gr_source%nzt
            levels_source(j+1) = gr_source%zt(i,j)
        end do
        levels_source(gr_source%nzt+2) = gr_source%zm(i,gr_source%nzm)

        levels_target(1) = gr_target%zm(i,1)
        do j=1, gr_target%nzt
            levels_target(j+1) = gr_target%zt(i,j)
        end do
        levels_target(gr_target%nzt+2) = gr_target%zm(i,gr_target%nzm)

        remap_zm_values(i,:) = remap_vals_to_target( nlevels_source, nlevels_target, &
                                                     levels_source, levels_target, &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals(i,:), &
                                                     rho_lin_spline_levels(i,:), &
                                                     values_source(i,:) )
    end do
  end function remap_zm_values

  subroutine remap_forcings( ngrdcol, gr_dycore, gr_clubb, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             thlm_forcing_dycore, rtm_forcing_dycore, &
                             thlm_forcing, rtm_forcing )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_dycore, gr_clubb

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_dycore%nzt), intent(in) :: &
      thlm_forcing_dycore, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing_dycore     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_clubb%nzt), intent(out) :: &
      thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Begin Code ---------------------
    thlm_forcing = remap_zt_values( ngrdcol, gr_dycore, gr_clubb, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    thlm_forcing_dycore )
    rtm_forcing = remap_zt_values( ngrdcol, gr_dycore, gr_clubb, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   rtm_forcing_dycore )
  end subroutine remap_forcings

  subroutine check_remap_for_consistency( nlevel_source, nlevel_target, &
                                          levels_source, levels_target, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
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
    real( kind = core_rknd ), dimension(nlevel_source-1) :: &
      ones

    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      output_ones

    integer :: i, j

    logical :: l_consistent

    !--------------------- Begin Code ---------------------
    do i = 1, (nlevel_source-1)
      ones(i) = 1
    end do

    output_ones = remap_vals_to_target( nlevel_source, nlevel_target, &
                                        levels_source, levels_target, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        ones )

    l_consistent = .true.
    j = 1
    do while (l_consistent .and. j <= nlevel_target-1)
      if (abs( output_ones(j) - 1 ) > tol) then
        l_consistent = .false.
      end if
      j = j + 1
    end do

    if (.not. l_consistent) then
      write(fstderr,*) 'WARNING! The remap_vals_to_target function is not l_consistent'
      error stop 'Function should be l_consistent, something went wrong...'
    end if

  end subroutine check_remap_for_consistency

  subroutine check_remap_for_consistency_zt_zm( ngrdcol, gr_source, gr_target, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type(grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! the height of the levels in the source grid [m]

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! the height of the levels in the target grid [m]

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      ! check for values given on zt levels
      call check_remap_for_consistency( gr_source%nzm, gr_target%nzm, &
                                        gr_source%zm(i,:), gr_target%zm(i,:), &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals(i,:), &
                                        rho_lin_spline_levels(i,:) )

      ! check for values given on zm levels
      
      ! create grid levels for those grids for values on the zm levels
      levels_source(1) = gr_source%zm(i,1)
      do j=1, gr_source%nzt
          levels_source(j+1) = gr_source%zt(i,j)
      end do
      levels_source(gr_source%nzt+2) = gr_source%zm(i,gr_source%nzm)

      levels_target(1) = gr_target%zm(i,1)
      do j=1, gr_target%nzt
          levels_target(j+1) = gr_target%zt(i,j)
      end do
      levels_target(gr_target%nzt+2) = gr_target%zm(i,gr_target%nzm)

      call check_remap_for_consistency( gr_source%nzt+2, gr_target%nzt+2, &
                                        levels_source, levels_target, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals(i,:), &
                                        rho_lin_spline_levels(i,:) )
      
    end do
  end subroutine check_remap_for_consistency_zt_zm

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

  subroutine check_mass_conservation_zt_zm( ngrdcol, gr_source, gr_target, &
                                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                            rho_lin_spline_levels )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type( grid ), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: levels_source ! [m]

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: levels_target ! [m]

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      ! check mass for grids that build cells for values on the zt levels
      call check_mass_conservation( gr_source%nzm, gr_target%nzm, &
                                    gr_source%zm(i,:), gr_target%zm(i,:), &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                    rho_lin_spline_levels(i,:) )

      ! check mass for grids that build cells for values on the zm levels
      levels_source(1) = gr_source%zm(i,1)
      do j=1, gr_source%nzt
          levels_source(j+1) = gr_source%zt(i,j)
      end do
      levels_source(gr_source%nzt+2) = gr_source%zm(i,gr_source%nzm)

      levels_target(1) = gr_target%zm(i,1)
      do j=1, gr_target%nzt
          levels_target(j+1) = gr_target%zt(i,j)
      end do
      levels_target(gr_target%nzt+2) = gr_target%zm(i,gr_target%nzm)

      call check_mass_conservation( gr_source%nzt+2, gr_target%nzt+2, &
                                    levels_source, levels_target, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                    rho_lin_spline_levels(i,:) )
    end do

  end subroutine check_mass_conservation_zt_zm

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

    if (abs(integral_target - integral_source) > tol) then

      write(fstderr,*) "WARNING! Integral for field was not conserved."
      write(fstderr,*) "Integral was ", integral_target, " instead of ", integral_source
      error stop 'Integral should be conserved, something went wrong...'

    end if

  end subroutine check_vertical_integral_conservation

  subroutine check_vertical_integral_conservation_all_zt_values( ngrdcol, &
                                                                 gr_source, gr_target, &
                                                                 total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 field_source, &
                                                                 field_target )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type(grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_source%nzt) :: &
      field_source           ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_target%nzt) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      call check_vertical_integral_conservation( total_idx_rho_lin_spline, &            ! intent(in)
                                                 rho_lin_spline_vals(i,:), &            ! intent(in)
                                                 rho_lin_spline_levels(i,:), &          ! intent(in)
                                                 gr_source%nzm, gr_target%nzm, &        ! intent(in)
                                                 gr_source%zm(i,:), gr_target%zm(i,:), &! intent(in)
                                                 field_source(i,:), field_target(i,:) ) ! intent(in)
    end do
  end subroutine check_vertical_integral_conservation_all_zt_values

  subroutine check_vertical_integral_conservation_all_zm_values( ngrdcol, &
                                                                 gr_source, gr_target, &
                                                                 total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 field_source, &
                                                                 field_target )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type(grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_source%nzm) :: &
      field_source           ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_target%nzm) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! heights of the levels of the source grid [m]

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! heights of the levels of the target grid [m]

    !--------------------- Begin Code ---------------------
    do i = 1, ngrdcol
      ! build the grid for values given on the zm levels
      levels_source(1) = gr_source%zm(i,1)
      do j=1, gr_source%nzt
          levels_source(j+1) = gr_source%zt(i,j)
      end do
      levels_source(gr_source%nzt+2) = gr_source%zm(i,gr_source%nzm)

      levels_target(1) = gr_target%zm(i,1)
      do j=1, gr_target%nzt
          levels_target(j+1) = gr_target%zt(i,j)
      end do
      levels_target(gr_target%nzt+2) = gr_target%zm(i,gr_target%nzm)

      call check_vertical_integral_conservation( total_idx_rho_lin_spline, &            ! intent(in)
                                                 rho_lin_spline_vals(i,:), &            ! intent(in)
                                                 rho_lin_spline_levels(i,:), &          ! intent(in)
                                                 gr_source%nzt+2, gr_target%nzt+2, &    ! intent(in)
                                                 levels_source, levels_target, &        ! intent(in)
                                                 field_source(i,:), field_target(i,:) ) ! intent(in)
    end do
  end subroutine check_vertical_integral_conservation_all_zm_values

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

  subroutine normalize_grid_density( gr_dens_idx, gr_dens_z, gr_dens, &
                                     min_dens, num_levels, &
                                     gr_dens_norm_z, gr_dens_norm )

    ! Description:
    ! Takes the linear piecewise grid density function (gr_dens_z, gr_dens) each of size
    ! gr_dens_idx and normalizes this density, such that the integral is the
    ! number of desired grid levels (num_levels) -1 and the minimum is above the minimal
    ! density (min_dens)

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      gr_dens_idx, &   ! number of elements in gr_dens_z and gr_dens []
      num_levels       ! desired number of grid levels []

    real( kind = core_rknd ), dimension(gr_dens_idx), intent(in) :: &
      gr_dens_z, & ! the grid density z coordinates [m]
      gr_dens      ! the grid density given at the z coordinates [# levs/meter]
      ! gr_dens_z and gr_dens should be ordered from the bottom to the top

    real( kind = core_rknd ), intent(in) :: &
      min_dens   ! the desired minimum grid density [# levs/meter]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(gr_dens_idx), intent(out) :: &
      gr_dens_norm_z, & ! the normalized grid density z coordinates [m]
      gr_dens_norm      ! the normalized grid density values given at the z coordinates [# levs/m]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid [m]
      grid_top, &               ! the highest point of the grid [m]
      delta_s, &                ! delta by which the function first gets shifted,
                                ! so the minimum of the function is above min_dens [m]
      integral_before_shift, &  ! the integral of the shifted density function before normalizing it
      equi_dens, &              ! the grid density for an equi-distant grid [# levs/m]
      factor, &                 ! factor to multiply shifted function with to end up with an
                                ! integral of num_levels-1
      h, &                      ! absolut distance from lowest to highest point of the grid [m]
      diff, zero

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_z(1)
    grid_top = gr_dens_z(gr_dens_idx)
    diff = (min_dens - minval(gr_dens))
    zero = 0.0_core_rknd
    delta_s = maxval([diff, zero])
    h = grid_top - grid_sfc
    equi_dens = (num_levels-1)/h

    ! first shift function such that minimum is above min_dens, then shift, so that min_dens is zero
    do i = 1, gr_dens_idx
        gr_dens_norm_z(i) = gr_dens_z(i)
        gr_dens_norm(i) = gr_dens(i) + delta_s - min_dens
    end do

    integral_before_shift = calc_integral(gr_dens_idx, gr_dens_norm_z, gr_dens_norm)

    ! if the integral is zero, then all gr_dens_norm were the same and smaller or equal to min_dens
    ! so in that case we have an equi-distant grid
    if (integral_before_shift < tol) then
        ! set to equi-distant grid
        do i = 1, gr_dens_idx
            gr_dens_norm(i) = equi_dens
        end do
    else
        factor = (num_levels-1-h*min_dens)/integral_before_shift
        do i = 1, gr_dens_idx
            gr_dens_norm(i) = factor*(gr_dens_norm(i)) + min_dens
        end do
    end if

    if ( clubb_at_least_debug_level_api( 2 ) ) then
      ! check if minimum of function is actually bigger or equal to min_dens
      if ((minval(gr_dens_norm) + tol) < min_dens) then
        write(fstderr,*) "WARNING! The minimum of the normalized function is below min_dens."
        error stop 'Something went wrong with calculating the normalized function...'
      end if
      ! check if the integral of the normalized function is actually num_levels-1
      if (abs(calc_integral(gr_dens_idx, gr_dens_norm_z, gr_dens_norm) - (num_levels-1)) > tol) then
        write(fstderr,*) "WARNING! The integral of the normalized function is not num_levels-1."
        error stop 'Something went wrong with calculating the normalized function...'
      end if
    end if

  end subroutine normalize_grid_density

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

                    if ((grid_level < gr_dens_norm_z(prev_x_ind)) &
                        .or. (grid_level > gr_dens_norm_z(prev_x_ind+1)) & 
                        .or. (grid_level < grid_heights(i-1))) then
                        ! first solution was not the one we were looking for
                        grid_level = -1*p_pq_formula/2 - sqrt((p_pq_formula/2)**2 - q_pq_formula)
                        if ((grid_level < gr_dens_norm_z(prev_x_ind)) &
                            .or. (grid_level > gr_dens_norm_z(prev_x_ind+1)) & 
                            .or. (grid_level < grid_heights(i-1))) then
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
  
  subroutine check_grid( grid_heights_idx, grid_heights, &
                         desired_num_levels, desired_min_dens, &
                         desired_grid_sfc, desired_grid_top )

    ! Description:
    ! Checks if the grid is valid.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      grid_heights_idx, &   ! number of elements in grid_heights []
      desired_num_levels    ! desired number of grid levels []
    ! Note: if everything worked fine, those two integers should be the same

    real( kind = core_rknd ), dimension(grid_heights_idx), intent(in) :: &
      grid_heights ! the grid heights ordered from bottom to top [m]

    real( kind = core_rknd ), intent(in) :: &
      desired_min_dens, &   ! the desired minimum grid density [# levs/m]
      desired_grid_sfc, &   ! desired surface height of the grid [m]
      desired_grid_top      ! desired top height of the grid [m]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: &
      max_dist  ! the maximum distance between two grid levels, depending on desired_min_dens [m]

    !--------------------- Begin Code ---------------------
    ! check if first grid element is grid_sfc
    if (abs(grid_heights(1) - desired_grid_sfc) > tol) then
        write(fstderr,*) "WARNING! The first grid level is not correct. It should be the grids", &
                         " surface: ", desired_grid_sfc, " and not ", grid_heights(1)
        error stop 'Something went wrong with generating the grid'
    end if

    max_dist = one/desired_min_dens

    do i = 1, (grid_heights_idx-1)
        ! check if grid heights are getting bigger
        if (grid_heights(i) >= grid_heights(i+1)) then
            write(fstderr,*) "WARNING! The grid levels are not in the right order. ", &
                             "They should get bigger with the index."
            error stop 'Something went wrong with generating the grid'
        end if
        ! check if grid fulfills min_dens
        if ((grid_heights(i+1) - grid_heights(i)) > max_dist + tol) then
            write(fstderr,*) "WARNING! The grid has a distance=", &
                             (grid_heights(i+1) - grid_heights(i)), &
                             " which is bigger than the max_dist=", max_dist, &
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

  function create_grid_from_grid_density_func( gr_dens_idx, &
                                               gr_dens_z, gr_dens, &
                                               min_dens, &
                                               num_levels ) &
                                             result (grid_heights)
    ! Description:
    ! Creates the grid from the unnormalized grid density function following the equi-distribution
    ! principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      gr_dens_idx, &  ! total numer of indices of gr_dens_z and gr_dens []
      num_levels      ! number of levels the new grid should have []

    real( kind = core_rknd ), intent(in) ::  &
      min_dens  ! the minimum density the new grid should have [# levs/meter]

    real( kind = core_rknd ), dimension(gr_dens_idx), intent(in) ::  &
      gr_dens_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function [m]
      gr_dens         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function
                      ! [# levs/meter]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx

    real( kind = core_rknd ), dimension(gr_dens_idx) ::  &
      gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized piecewise
                          ! linear grid density function [m]
      gr_dens_norm        ! the density at the given z coordinates of the connection points of the
                          ! normalized piecewise linear grid density function [# levs/meter]

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_z(1)
    grid_top = gr_dens_z(gr_dens_idx)

    call normalize_grid_density( gr_dens_idx, gr_dens_z, gr_dens, &
                                 min_dens, num_levels, &
                                 gr_dens_norm_z, gr_dens_norm )

    grid_heights = create_grid_from_normalized_grid_density_func( num_levels, &
                                                                  gr_dens_idx, &
                                                                  gr_dens_norm_z, &
                                                                  gr_dens_norm )

    ! TODO how to read may change if we incorporate ngrdcol!,
    ! TODO since size does not work then, but should be size/ngrdcol?,
    ! TODO since all columns should have same number of grid levels?
    !grid_heights_idx = shape(grid_heights)
    grid_heights_idx = size(grid_heights) ! counts the scalar values

    if ( clubb_at_least_debug_level_api( 2 ) ) then
        call check_grid( grid_heights_idx, grid_heights, &
                         num_levels, min_dens, &
                         grid_sfc, grid_top )
    end if 

    return

  end function create_grid_from_grid_density_func


  subroutine adapt_grid( ngrdcol, total_idx_density_func, &
                         density_func_z, density_func_dens, &
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
                         Lscale, err_info )
    ! Description:
    ! Adapts the grid based on the density function and interpolates all values to the new grid.
    !-----------------------------------------------------------------------

    use constants_clubb, only : &
        fstderr     ! Fortran file unit I/O constant

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error

    use clubb_api_module, only: &
      setup_grid_api

    use calc_pressure, only: &
        init_pressure    ! Procedure(s)

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      hydromet_dim, &
      sclr_dim, &
      edsclr_dim, &
      total_idx_density_func, & ! total numer of indices of density_func_z and density_func_dens []
      idx_thvm   ! numer of indices of thvm

    real( kind = core_rknd ), dimension(total_idx_density_func), intent(in) ::  &
      density_func_z,  &   ! the height values of the connection points of the piecewise linear
                           ! grid density function/profile [m]
      density_func_dens    ! the density values of the connection points of the piecewise linear
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

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Local Variables ---------------------
    integer :: &
        num_levels, &                ! number of levels the new grid should have []
        total_idx_rho_lin_spline, &  ! total number of indices of the rho density
                                     ! piecewise linear function []
        grid_type, &
        i

    real( kind = core_rknd ) ::  &
      equi_dens,  &    ! density of the equidistant grid [# levs/meter]
      min_dens,   &    ! minimum boundary of the grids density [# levs/meter]
      deltaz

    real( kind = core_rknd ), dimension(gr%nzm) ::  &
      new_gr_zm      ! zm levels for the adapted grid based on the density function [m]

    real( kind = core_rknd ), dimension(gr%nzt) ::  &
      thermodynamic_heights_placeholder  ! placeholder for the zt levels, but not needed

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      rho_lin_spline_vals, rho_lin_spline_levels
    
    type( grid ) :: new_gr

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzt) ::  &
      thvm_new_grid

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      p_in_Pa_zm, exner_zm ! just placeholder variables; are only needed for subroutine call

    !--------------------- Begin Code ---------------------
    ! Setup the new grid
    num_levels = gr%nzm
    equi_dens = (num_levels-1)/(gr%zm(1,gr%nzm) - gr%zm(1,1))
    min_dens = 0.5*equi_dens ! TODO find better way to determine a good min_dens value
    new_gr_zm = create_grid_from_grid_density_func( total_idx_density_func, &
                                                    density_func_z, density_func_dens, &
                                                    min_dens, &
                                                    num_levels )
    deltaz = 0.0
    grid_type = 3
    call setup_grid_api( num_levels, sfc_elevation(1), l_implemented, &         ! intent(in)
                         .true., grid_type, deltaz, gr%zm(1,1), gr%zm(1,num_levels), & ! intent(in)
                         new_gr_zm, thermodynamic_heights_placeholder, &        ! intent(in)
                         new_gr, err_info )                                     ! intent(inout)

    if ( clubb_at_least_debug_level_api(0) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid_api in adapt_grid"
        return
      end if
    end if

    ! Set the density values to use for interpolation for mass calculation
    total_idx_rho_lin_spline = gr%nzm
    rho_lin_spline_vals = rho_ds_zm
    rho_lin_spline_levels = gr%zm

    if ( grid_remap_method == 1 ) then
      ! Remap all zt values
      thlm_forcing = remap_zt_values( ngrdcol, gr, new_gr, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      thlm_forcing )
      rtm_forcing = remap_zt_values( ngrdcol, gr, new_gr, &
                                     total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     rtm_forcing )
      um_forcing = remap_zt_values( ngrdcol, gr, new_gr, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    um_forcing )
      vm_forcing = remap_zt_values( ngrdcol, gr, new_gr, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    vm_forcing )
      wm_zt = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wm_zt )
      rho = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             rho )
      rho_ds_zt = remap_zt_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   rho_ds_zt )
      invrs_rho_ds_zt = remap_zt_values( ngrdcol, gr, new_gr, &
                                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         invrs_rho_ds_zt )
      thv_ds_zt = remap_zt_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   thv_ds_zt )
      rfrzm = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               rfrzm )
      do i = 1, hydromet_dim
          wp2hmp(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                           total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           wp2hmp(:,:,i) )
          rtphmp_zt(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              rtphmp_zt(:,:,i) )
          thlphmp_zt(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               thlphmp_zt(:,:,i) )
      end do
      do i = 1, sclr_dim
          sclrm_forcing(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  sclrm_forcing(:,:,i) )
          sclrm(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          sclrm(:,:,i) )
          sclrp3(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           sclrp3(:,:,i) )
      end do
      do i = 1, edsclr_dim
          edsclrm_forcing(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                                    total_idx_rho_lin_spline, &
                                                    rho_lin_spline_vals, &
                                                    rho_lin_spline_levels, &
                                                    edsclrm_forcing(:,:,i) )
          edsclrm(:,:,i) = remap_zt_values( ngrdcol, gr, new_gr, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            edsclrm(:,:,i) )
      end do
      rtm_ref = remap_zt_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 rtm_ref )
      thlm_ref = remap_zt_values( ngrdcol, gr, new_gr, &
                                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  thlm_ref )
      um_ref = remap_zt_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                um_ref )
      vm_ref = remap_zt_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                vm_ref )
      ug = remap_zt_values( ngrdcol, gr, new_gr, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            ug )
      vg = remap_zt_values( ngrdcol, gr, new_gr, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            vg )
      um = remap_zt_values( ngrdcol, gr, new_gr, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            um )
      vm = remap_zt_values( ngrdcol, gr, new_gr, &
                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                            rho_lin_spline_levels, &
                            vm )
      up3 = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             up3 )
      vp3 = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             vp3 )
      rtm = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             rtm )
      thlm = remap_zt_values( ngrdcol, gr, new_gr, &
                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              thlm )
      rtp3 = remap_zt_values( ngrdcol, gr, new_gr, &
                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              rtp3 )
      thlp3 = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               thlp3 )
      wp3 = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             wp3 )
      ! remap thvm values to new grid, to calculate new p_in_Pa
      thvm_new_grid = remap_zt_values( ngrdcol, gr, new_gr, &
                                       total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       thvm )
      ! calculate p_in_Pa instead of remapping directly since it can run into problems if
      ! for example the two highest levels have the same pressure value, which could be happening
      ! with the remapping, also calculate exner accordingly
      call init_pressure( ngrdcol, new_gr, thvm_new_grid, p_sfc, &
                          p_in_Pa, exner, p_in_Pa_zm, exner_zm )
      rcm = remap_zt_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             rcm )
      cloud_frac = remap_zt_values( ngrdcol, gr, new_gr, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    cloud_frac )
      wp2thvp = remap_zt_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 wp2thvp )
      wp2rtp = remap_zt_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                wp2rtp )
      wp2thlp = remap_zt_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 wp2thlp )
      wpup2 = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wpup2 )
      wpvp2 = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wpvp2 )
      ice_supersat_frac = remap_zt_values( ngrdcol, gr, new_gr, &
                                           total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           ice_supersat_frac )
      um_pert = remap_zt_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 um_pert )
      vm_pert = remap_zt_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 vm_pert )
      rcm_in_layer = remap_zt_values( ngrdcol, gr, new_gr, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      rcm_in_layer )
      cloud_cover = remap_zt_values( ngrdcol, gr, new_gr, &
                                     total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     cloud_cover )
      w_up_in_cloud = remap_zt_values( ngrdcol, gr, new_gr, &
                                       total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       w_up_in_cloud )
      w_down_in_cloud = remap_zt_values( ngrdcol, gr, new_gr, &
                                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         w_down_in_cloud )
      cloudy_updraft_frac = remap_zt_values( ngrdcol, gr, new_gr, &
                                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             cloudy_updraft_frac )
      cloudy_downdraft_frac = remap_zt_values( ngrdcol, gr, new_gr, &
                                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               cloudy_downdraft_frac )
      Kh_zt = remap_zt_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               Kh_zt )
      Lscale = remap_zt_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                Lscale )

    
      ! Remap all zm values
      wprtp_forcing = remap_zm_values( ngrdcol, gr, new_gr, &
                                       total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       wprtp_forcing )
      wpthlp_forcing = remap_zm_values( ngrdcol, gr, new_gr, &
                                        total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        wpthlp_forcing )
      rtp2_forcing = remap_zm_values( ngrdcol, gr, new_gr, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      rtp2_forcing )
      thlp2_forcing = remap_zm_values( ngrdcol, gr, new_gr, &
                                       total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       thlp2_forcing )
      rtpthlp_forcing = remap_zm_values( ngrdcol, gr, new_gr, &
                                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         rtpthlp_forcing )
      wm_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wm_zm )
      rho_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                rho_zm )
      rho_ds_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   rho_ds_zm )
      invrs_rho_ds_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         invrs_rho_ds_zm )
      thv_ds_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   thv_ds_zm )
      do i = 1, hydromet_dim
          wphydrometp(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                wphydrometp(:,:,i) )
      end do
      upwp = remap_zm_values( ngrdcol, gr, new_gr, &
                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              upwp )
      vpwp = remap_zm_values( ngrdcol, gr, new_gr, &
                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              vpwp )
      up2 = remap_zm_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             up2 )
      vp2 = remap_zm_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             vp2 )
      wprtp = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wprtp )
      wpthlp = remap_zm_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                wpthlp )
      rtp2 = remap_zm_values( ngrdcol, gr, new_gr, &
                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                              rho_lin_spline_levels, &
                              rtp2 )
      thlp2 = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               thlp2 )
      rtpthlp = remap_zm_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 rtpthlp )
      wp2 = remap_zm_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             wp2 )
      do i = 1, sclr_dim
          wpsclrp(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                            total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            wpsclrp(:,:,i) )
          sclrp2(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                           total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           sclrp2(:,:,i) )
          sclrprtp(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             sclrprtp(:,:,i) )
          sclrpthlp(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              sclrpthlp(:,:,i) )
          sclrpthvp(:,:,i) = remap_zm_values( ngrdcol, gr, new_gr, &
                                              total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              sclrpthvp(:,:,i) )
      end do
      wpthvp = remap_zm_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                wpthvp )
      rtpthvp = remap_zm_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 rtpthvp )
      thlpthvp = remap_zm_values( ngrdcol, gr, new_gr, &
                                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  thlpthvp )
      uprcp = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               uprcp )
      vprcp = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               vprcp )
      rc_coef_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                    total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    rc_coef_zm )
      wp4 = remap_zm_values( ngrdcol, gr, new_gr, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels, &
                             wp4 )
      wp2up2 = remap_zm_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                wp2up2 )
      wp2vp2 = remap_zm_values( ngrdcol, gr, new_gr, &
                                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                wp2vp2 )
      upwp_pert = remap_zm_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   upwp_pert )
      vpwp_pert = remap_zm_values( ngrdcol, gr, new_gr, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   vpwp_pert )
      wprcp = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               wprcp )
      invrs_tau_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      invrs_tau_zm )
      Kh_zm = remap_zm_values( ngrdcol, gr, new_gr, &
                               total_idx_rho_lin_spline, rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               Kh_zm )
      thlprcp = remap_zm_values( ngrdcol, gr, new_gr, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 thlprcp )
    else
      write(fstderr,*) 'There is no method implemented for grid_remap_method=', &
                       grid_remap_method, '. Try another integer value.'
      error stop 'Invalid value for grid_remap_method.'
    end if
    
    gr = new_gr
  end subroutine adapt_grid

  subroutine calc_grid_dens_func( ngrdcol, nzt, zt, &
                                  Lscale, wp2_zt, &
                                  grid_sfc, grid_top, &
                                  gr_dens_z, gr_dens )
    ! Description:
    ! Creates an unnormalized grid density function from Lscale

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, & 
      nzt

    real( kind = core_rknd ), intent(in) ::  &
      grid_sfc, &  ! the grids surface; so the first level in the grid
                   ! density function has this height [m]
      grid_top     ! the grids top; so the last level in the grid
                   ! density function has this height [m]

    real( kind = core_rknd ), dimension(ngrdcol, nzt), intent(in) ::  &
      zt, &      ! levels at which the values are given [m]
      Lscale, &  ! Length scale   [m]
      wp2_zt     ! w'^2 on thermo. grid [m^2/s^2]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(nzt), intent(out) ::  &
      gr_dens_z,  &    ! the z value coordinates of the connection points of the piecewise linear
                       ! grid density function [m]
      gr_dens  ! the values of the connection points of the piecewise linear
               ! grid density function, given on the z values of gr_dens_z [# levs/meter]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: threshold

    !--------------------- Begin Code ---------------------
    threshold = 0.001
    do i = 1, nzt
        gr_dens_z(i) = zt(1,i)
        if ( wp2_zt(1,i) > threshold ) then
          gr_dens(i) = 1/Lscale(1,i)
        else
          gr_dens(i) = threshold
        end if
    end do

    ! TODO find better way to ensure grid density function goes to zm(1) and zm(n)
    if (gr_dens_z(1) > grid_sfc) then
        gr_dens_z(1) = grid_sfc
    end if

    if (gr_dens_z(nzt) < grid_top) then
        gr_dens_z(nzt) = grid_top
    end if

  end subroutine calc_grid_dens_func
!===============================================================================

end module grid_adaptation_module
