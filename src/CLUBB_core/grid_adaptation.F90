!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_adaptation_module


  use clubb_precision, only: &
      core_rknd ! Variable(s)

  use grid_class, only: &
      grid ! Type

  use constants_clubb, only: &
      fstdout, fstderr ! Constants

  use error_code, only: &
      clubb_at_least_debug_level

  implicit none

  public :: setup_simple_gr_dycore, lin_interpolate, &
            check_conservation, check_consistent, check_monotone, &
            remapping_matrix, remap, &
            interpolate_values_zt, interpolate_values_zm, interpolate_forcings, &
            check_remap_for_consistency, check_remap_for_consistency_all, &
            check_mass_conservation, check_mass_conservation_all, &
            check_vertical_integral_conservation, &
            check_vertical_integral_conservation_all_zt_values, &
            check_vertical_integral_conservation_all_zm_values, &
            calc_mass_over_grid_intervals, vertical_integral_conserve_mass

  private :: setup_gr

  private ! Default Scoping  


  contains

  !=============================================================================
  subroutine setup_gr( nzmax, grid_type, &
                       zm_grid_fname, zt_grid_fname, &
                       file_unit, sfc_elevation, &
                       l_implemented, deltaz, &
                       zm_init, zm_top, &
                       gr )

    use grid_class, only: &
      read_grid_heights

    use clubb_api_module, only: &
      setup_grid_api

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: nzmax

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    character(len=*), intent(in) :: & 
      zm_grid_fname,  & ! Path and filename of file for momentum level altitudes
      zt_grid_fname     ! Path and filename of file for thermodynamic level altitudes

    integer, intent(in) :: &
      file_unit   ! Unit number for zt_grid_fname & zm_grid_fname (based on the OpenMP thread)

    real( kind = core_rknd ), intent(in) ::  &
        sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical, intent(in) :: l_implemented

    real( kind = core_rknd ), intent(in) ::  & 
        deltaz,   & ! Vertical grid spacing                  [m]
        zm_init,  & ! Initial grid altitude (momentum level) [m]
        zm_top      ! Maximum grid altitude (momentum level) [m]

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- Local Variables ---------------------
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), dimension(nzmax) :: & 
      momentum_heights      ! Momentum level altitudes (file input)      [m]
    real( kind = core_rknd ), dimension(nzmax-1) :: &
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

    !--------------------- Begin Code ---------------------
    call read_grid_heights( nzmax, grid_type,  &             ! intent(in)
                            zm_grid_fname, zt_grid_fname, &  ! intent(in)
                            file_unit, &                     ! intent(in)
                            momentum_heights, &              ! intent(out)
                            thermodynamic_heights )          ! intent(out)

    call setup_grid_api( nzmax, sfc_elevation, l_implemented, &     ! intent(in)
                         grid_type, deltaz, zm_init, zm_top, &      ! intent(in)
                         momentum_heights, thermodynamic_heights, & ! intent(in)
                         gr )                                       ! intent(inout)

  end subroutine setup_gr

  subroutine setup_simple_gr_dycore( gr )

    implicit none

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- Local Variables ---------------------
    integer :: nzmax

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer :: grid_type

    character(len=0) :: & 
      zm_grid_fname,  & ! Path and filename of file for momentum level altitudes
      zt_grid_fname     ! Path and filename of file for thermodynamic level altitudes

    integer :: &
      file_unit   ! Unit number for zt_grid_fname & zm_grid_fname (based on the OpenMP thread)

    real( kind = core_rknd ) ::  &
        sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical :: l_implemented

    real( kind = core_rknd ) ::  & 
        deltaz,   & ! Vertical grid spacing                  [m]
        zm_init,  & ! Initial grid altitude (momentum level) [m]
        zm_top      ! Maximum grid altitude (momentum level) [m]

    nzmax = 31
    grid_type = 1
    zm_grid_fname = ''
    zt_grid_fname = ''
    file_unit = 1
    sfc_elevation = 0.0
    l_implemented = .false.
    deltaz = 132.0
    zm_init = 0.0
    zm_top = 3960.0

    call setup_gr( nzmax, grid_type, &
                   zm_grid_fname, zt_grid_fname, &
                   file_unit, sfc_elevation, &
                   l_implemented, deltaz, &
                   zm_init, zm_top, &
                   gr )

  end subroutine setup_simple_gr_dycore

  function lin_interpolate( size_interpolate, size_current, &
                            interpolate_altitudes, current_altitudes, &
                            current_values )
    ! interpolates from the values (gr_current_values) on the current grid (gr_current) to the
    ! gr_interpolate grid

    use interpolation, only: &
      lin_interpolate_two_points

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
        size_interpolate, &
        size_current

    real( kind = core_rknd ), dimension(size_interpolate) :: &
      interpolate_altitudes ! altitudes where the values should be interpolated to

    real( kind = core_rknd ), dimension(size_current) :: &
      current_altitudes ! altitudes where the given values are currently on

    real( kind = core_rknd ), dimension(size_current) :: &
      current_values ! given values on the current_altitudes altitudes that should be interpolated 
                     ! to the interpolate_altitudes altitudes

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(size_interpolate) :: &
      lin_interpolate ! interpolated values with the dimension of gr_interpolate
    integer :: i, k
    logical :: calc_done
    real( kind = core_rknd ) :: epsilon

    !--------------------- Begin Code ---------------------

    epsilon = 0.01

    do i = 1, size_interpolate
      calc_done = .false.
      ! find nearest zt level of gr_current for zt(i) of gr_interpolate grid
      ! if there is no gr_current zt level above or below the gr_interpolate grid level, 
      ! then the gr_current value of the closest level is taken
      k = 1
      do while ((.not. calc_done) .and. (k <= size_current))
        if (abs(interpolate_altitudes(i) - current_altitudes(k)) < epsilon) then
          lin_interpolate(i) = current_values(k)
          calc_done = .true.
        else if (interpolate_altitudes(i) < current_altitudes(k)) then
          if (k > 1) then
            lin_interpolate(i) = lin_interpolate_two_points( interpolate_altitudes(i), &
                                                             current_altitudes(k), &
                                                             current_altitudes(k-1), &
                                                             current_values(k), &
                                                             current_values(k-1) )
          else ! k=1
            lin_interpolate(i) = current_values(1)
          endif
          calc_done = .true.
        endif
        k = k + 1
      end do
      if ((.not. calc_done) .and. (k > size_current)) then
        lin_interpolate(i) = current_values(size_current)
        calc_done = .true.
      endif
    end do

  end function lin_interpolate

  subroutine check_conservation( R_ij, dim1, dim2, &
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
      R_ij  ! matrix to apply to input values for remapping (R_ij)

    real( kind = core_rknd ), dimension(dim2+1), intent(in) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(dim1+1), intent(in) :: &
      levels_target ! the height of the levels in the target grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: epsilon, weighted_col_sum

    real( kind = core_rknd ), dimension(dim1) :: &
        Jt ! local weights on the target grid (mass over every grid cell)

    real( kind = core_rknd ), dimension(dim2) :: &
        Js ! local weights on the source grid (mass over every grid cell)

    logical :: conservative

    !--------------------- Begin Code ---------------------
    conservative = .true.
    epsilon = 0.000001

    Jt = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        dim1+1, levels_target )
    
    Js = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        dim2+1, levels_source )

    j = 1
    do while (conservative .and. j <= dim2)
      weighted_col_sum = 0
      do i = 1, dim1
        weighted_col_sum = weighted_col_sum + R_ij(i,j) * Jt(i)
      end do
      if (abs(weighted_col_sum - Js(j)) > epsilon) then
        conservative = .false.
      end if
      j = j + 1
    end do

    if (.not. conservative) then
      write(*,*) 'remapping operator is not conservative'
    end if

  end subroutine check_conservation

  subroutine check_consistent( R_ij, dim1, dim2 )
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
      R_ij  ! matrix to apply to input values for remapping (R_ij)

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: row_sum, epsilon

    logical :: consistent

    !--------------------- Begin Code ---------------------
    consistent = .true.
    epsilon = 0.000001

    i = 1
    do while (consistent .and. i <= dim1)
      row_sum = 0
      do j = 1, dim2
        row_sum = row_sum + R_ij(i,j)
      end do
      if (abs(row_sum - 1) > epsilon) then
        consistent = .false.
      end if
      i = i + 1
    end do

    if (.not. consistent) then
      write(*,*) 'remapping operator is not consistent'
    end if

  end subroutine check_consistent

  subroutine check_monotone( R_ij, dim1, dim2 )
    ! checks monotonicity with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (14)
    ! Defintion of monotonicity by Ullrich and Taylor: A remapping operator R is monotone if the 
    ! remapping operation cannot introduce additional global extrema

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij)

    !--------------------- Local Variables ---------------------
    integer :: i, j

    logical :: monotone

    !--------------------- Begin Code ---------------------
    monotone = .true.

    i = 1
    do while (monotone .and. i <= dim1)
      do j = 1, dim2
        if (R_ij(i,j) < 0) then
          monotone = .false.
        end if
      end do
      i = i + 1
    end do

    if (.not. monotone) then
      write(*,*) 'remapping operator is not monotone'
    end if

  end subroutine check_monotone

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
      nlevel_source, & ! number of levels in the target grid
      nlevel_target    ! number of levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel_target-1,nlevel_source-1) :: &
      remapping_matrix ! matrix to apply to input values for remapping (R_ij)

    integer :: i, j

    real( kind = core_rknd ) :: omega_ov, omega_ov_upper, omega_ov_lower

    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      omega_ts ! densities of all intervals in the target grid

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

    if ( clubb_at_least_debug_level( 2 ) ) then
      call check_conservation( remapping_matrix, (nlevel_target-1), (nlevel_source-1), &
                               levels_source, levels_target, &
                               total_idx_rho_lin_spline, &
                               rho_lin_spline_vals, &
                               rho_lin_spline_levels)
      call check_consistent( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
      call check_monotone( remapping_matrix, (nlevel_target-1), (nlevel_source-1) )
    end if

  end function remapping_matrix

  function remap( nlevel_source, nlevel_target, &
                  levels_source, levels_target, &
                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                  rho_lin_spline_levels, &
                  gr_source_values )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nlevel_source, & ! number of levels in the target grid
      nlevel_target    ! number of levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(nlevel_source-1) :: &
      gr_source_values  ! given values on the gr_source grid that should be interpolated 
                        ! to the gr_target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      remap ! interpolated values with the dimension of gr_target

    integer :: i, j

    real( kind = core_rknd ), dimension(nlevel_target-1,nlevel_source-1) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij)

    !--------------------- Begin Code ---------------------

    R_ij = remapping_matrix( nlevel_source, nlevel_target, &
                             levels_source, levels_target, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels )
    ! matrix vector multiplication
    do i = 1, (nlevel_target-1)
      remap(i) = 0
      do j = 1, (nlevel_source-1)
        remap(i) = remap(i) + R_ij(i,j)*gr_source_values(j)
      end do
    end do

  end function remap

  function interpolate_values_zt( ngrdcol, gr_source, gr_target, &
                                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  values_source )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzt), intent(in) :: &
      values_source ! values given on the source grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_target%nzt) :: &
      interpolate_values_zt ! interpolated values on target grid

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    do i=1, ngrdcol
        interpolate_values_zt(i,:) = remap( gr_source%nzm, gr_target%nzm, &
                                            gr_source%zm(i,:), gr_target%zm(i,:), &
                                            total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                            rho_lin_spline_levels(i,:), &
                                            values_source(i,:) )
    end do
  end function interpolate_values_zt

  function interpolate_values_zm( ngrdcol, gr_source, gr_target, &
                                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  values_source )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzm), intent(in) :: &
      values_source ! values given on the source grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_target%nzm) :: &
      interpolate_values_zm ! interpolated values on target grid

    !--------------------- Local Variables ---------------------
    integer :: i, j

    integer :: &
      nlevels_source ! number of levels in source grid

    integer :: &
      nlevels_target ! number of levels in target grid

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! levels on source grid

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! levels on target grid

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

        interpolate_values_zm(i,:) = remap( nlevels_source, nlevels_target, &
                                            levels_source, levels_target, &
                                            total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                            rho_lin_spline_levels(i,:), &
                                            values_source(i,:) )
    end do
  end function interpolate_values_zm

  subroutine interpolate_forcings( ngrdcol, gr_dycore, gr_clubb, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   thlm_forcing_dycore, rtm_forcing_dycore, &
                                   thlm_forcing, rtm_forcing )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type (grid), intent(in) :: gr_dycore, gr_clubb

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(ngrdcol, gr_dycore%nzt), intent(in) :: &
      thlm_forcing_dycore, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing_dycore     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol, gr_clubb%nzt), intent(out) :: &
      thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Begin Code ---------------------
    thlm_forcing = interpolate_values_zt( ngrdcol, gr_dycore, gr_clubb, &
                                          total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          thlm_forcing_dycore )
    rtm_forcing = interpolate_values_zt( ngrdcol, gr_dycore, gr_clubb, &
                                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         rtm_forcing_dycore )
  end subroutine interpolate_forcings

  subroutine check_remap_for_consistency( nlevel_source, nlevel_target, &
                                          levels_source, levels_target, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nlevel_source, & ! number of levels in the target grid
      nlevel_target    ! number of levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel_source-1) :: &
      ones

    real( kind = core_rknd ), dimension(nlevel_target-1) :: &
      output_ones

    integer :: i, j

    real( kind = core_rknd ) :: epsilon

    logical :: consistent

    !--------------------- Begin Code ---------------------
    do i = 1, (nlevel_source-1)
      ones(i) = 1
    end do

    output_ones = remap( nlevel_source, nlevel_target, &
                         levels_source, levels_target, &
                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                         rho_lin_spline_levels, &
                         ones )

    consistent = .true.
    epsilon = 0.000001
    j = 1
    do while (consistent .and. j <= nlevel_target-1)
      if (abs( output_ones(j) - 1 ) > epsilon) then
        consistent = .false.
      end if
      j = j + 1
    end do

    if (.not. consistent) then
      write(*,*) 'remap function is not consistent'
    end if

  end subroutine check_remap_for_consistency

  subroutine check_remap_for_consistency_all( ngrdcol, gr_source, gr_target, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels )

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type(grid), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! the height of the levels in the target grid

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
  end subroutine check_remap_for_consistency_all

  subroutine check_mass_conservation( nlevel_source, nlevel_target, &
                                      levels_source, levels_target, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nlevel_source, & ! number of levels in the target grid
      nlevel_target    ! number of levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_source), intent(in) :: &
      levels_source ! the height of the levels in the source grid

    real( kind = core_rknd ), dimension(nlevel_target), intent(in) :: &
      levels_target ! the height of the levels in the target grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, &
      sum_source_mass, &
      sum_target_mass

    real( kind = core_rknd ), dimension(nlevel_source-1) :: source_mass

    real( kind = core_rknd ), dimension(nlevel_target-1) :: target_mass

    !--------------------- Begin Code ---------------------
    tol = 0.000001

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

      write(fstderr,*) "Mass for different grids was not the same."
      write(fstderr,*) "Mass was ", sum_target_mass, " instead of ", sum_source_mass

    end if

  end subroutine check_mass_conservation

  subroutine check_mass_conservation_all( ngrdcol, gr_source, gr_target, &
                                          total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                          rho_lin_spline_levels )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol

    type( grid ), intent(in) :: gr_source, gr_target

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: levels_source

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: levels_target

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

  end subroutine check_mass_conservation_all

  subroutine check_vertical_integral_conservation( total_idx_rho_lin_spline, &
                                                   rho_lin_spline_vals, &
                                                   rho_lin_spline_levels, &
                                                   nlevel_source, nlevel_target, &
                                                   levels_source, levels_target, &
                                                   field_source, field_target )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    integer, intent(in) :: &
      nlevel_source, &       ! number of levels of the source grid
      nlevel_target          ! number of levels of the target grid

    real( kind = core_rknd ), intent(in), dimension(nlevel_source) :: &
      levels_source           ! altitudes for the source grid

    real( kind = core_rknd ), intent(in), dimension(nlevel_target) :: &
      levels_target           ! altitudes for the target grid

    real( kind = core_rknd ), intent(in), dimension(nlevel_source-1) :: &
      field_source            ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(nlevel_target-1) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, &
      integral_source, &
      integral_target

    !--------------------- Begin Code ---------------------
    tol = 0.000001
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

      write(fstderr,*) "Integral for field was not conserved."
      write(fstderr,*) "Integral was ", integral_target, " instead of ", integral_source

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
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
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
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_source%nzm) :: &
      field_source           ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(ngrdcol, gr_target%nzm) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ), dimension(gr_source%nzt+2) :: &
      levels_source ! heights of the levels of the source grid

    real( kind = core_rknd ), dimension(gr_target%nzt+2) :: &
      levels_target ! heights of the levels of the target grid

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
      total_idx_lin_spline, &  ! The total numer of indices of the spline points
      nlevel                   ! The total numer of indices of the grid levels

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    real( kind = core_rknd ), dimension(nlevel), intent(in) :: &
      grid_levels ! levels of grid over which the masses shoud be calculated

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, & ! The difference tolerance for two real numbers to be considered the same
      val_below, & ! The rho value of the level or spline connection point below
      level_below, & ! The altitude of the level or spline connection point below
      mass_over_interval, & ! The total mass over one interval (zm(i) to zm(i+1)) in the grid
      upper_zm_level_rho

    real( kind = core_rknd ), dimension(nlevel-1) :: &
      mass_per_interval ! Array to store the mass of each interval
                        ! (grid_levels(i) to grid_levels(i+1))

    integer :: i, j

    logical :: level_found

    !--------------------- Begin Code ---------------------
    tol = 0.000000001
    ! check if highest value of lin spline is higher or equal to highest level of target grid
    if ( (grid_levels(nlevel) &
          - lin_spline_rho_levels(total_idx_lin_spline) ) > tol ) then
      write(*,*) 'Wrong input: highest level of lin_spline_rho_levels must', &
                 ' be higher or equal to highest zm level of grid_levels'
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
          write(*,*) 'Wrong input: lowest level of lin_spline_rho_levels must', &
                     ' be lower or equal to lowest zm level of grid_levels'
          call exit(1)
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
      total_idx_lin_spline  ! The total numer of indices of the spline points

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    integer, intent(in) :: & 
      nlevel  ! The total numer of indices for grid_levels

    real( kind = core_rknd ), dimension(nlevel-1), intent(in) ::  &
      field      ! The field to be vertically averaged   [Units vary]

    real( kind = core_rknd ), dimension(nlevel), intent(in) ::  &
      grid_levels

    ! Note:  The field and grid_levels points need to be arranged from
    !        lowest to highest in altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(nlevel-1) :: &
      mass_per_interval ! Array to store the mass of each interval of the target_grid

    real( kind = core_rknd ) :: &
      vertical_integral_conserve_mass

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
      g_idx  ! The total numer of indices of g_x and g_y

    real( kind = core_rknd ), dimension(g_idx), intent(in) ::  &
      g_x,  &    ! the x coordinates of the connection points of the piecewise linear function
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

  subroutine normalize_density( gd_idx, gd_x, gd_y, &
                                min_dens, num_levels, &
                                gd_norm_x, gd_norm_y )

    ! Description:
    ! Takes the linear piecewise grid density function (gd_x, gd_y) each of size gd_idx and
    ! normalizes this density, such that the integral is the
    ! number of desired grid levels (num_levels) -1 and the minimum is above the minimal
    ! density (min_dens)

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      gd_idx, &   ! number of elements in gd_x and gd_y
      num_levels  ! desired number of grid levels

    real( kind = core_rknd ), dimension(gd_idx), intent(in) :: &
      gd_x, & ! the grid density x coordinates
      gd_y    ! the grid density y coordinates
      ! gd_x and gd_y should be ordered from the bottom to the top

    real( kind = core_rknd ), intent(in) :: &
      min_dens   ! the desired minimum grid density

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(gd_idx), intent(out) :: &
      gd_norm_x, & ! the normalized grid density x coordinates
      gd_norm_y    ! the normalized grid density y coordinates

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid
      grid_top, &               ! the highest point of the grid
      delta_s, &                ! delta by which the function first gets shifted,
                                ! so the minimum of the function is above min_dens
      integral_before_shift, &  ! the integral of the shifted density function before normalizing it
      equi_dens, &              ! the grid density for an equi-distant grid
      factor, &                 ! factor to multiply shifted function with to end up with an
                                ! integral of num_levels-1
      h, &                      ! absolut distance from lowest to highest point of the grid
      diff, zero

    !--------------------- Begin Code ---------------------
    grid_sfc = gd_x(1)
    grid_top = gd_x(gd_idx)
    diff = (min_dens - minval(gd_y))
    zero = 0.0_core_rknd
    delta_s = maxval([diff, zero])
    h = grid_top - grid_sfc
    equi_dens = (num_levels-1)/h

    ! first shift function such that minimum is above min_dens, then shift, so that min_dens is zero
    do i = 1, gd_idx
        gd_norm_x(i) = gd_x(i)
        gd_norm_y(i) = gd_y(i) + delta_s - min_dens
    end do

    integral_before_shift = calc_integral(gd_idx, gd_norm_x, gd_norm_y)

    ! if the integral is zero, then all gd_norm_y were the same and smaller or equal to min_dens
    ! so in that case we have an equi-distant grid
    if (integral_before_shift < 1.0e-6_core_rknd) then
        ! set to equi-distant grid
        do i = 1, gd_idx
            gd_norm_y(i) = equi_dens
        end do
    else
        factor = (num_levels-1-h*min_dens)/integral_before_shift
        do i = 1, gd_idx
            gd_norm_y(i) = factor*(gd_norm_y(i)) + min_dens
        end do
    end if

    if ( clubb_at_least_debug_level( 2 ) ) then
      ! check if minimum of function is actually bigger or equal to min_dens
      if (minval(gd_norm_y) < min_dens) then
        write(fstderr,*) "Warning! The minimum of the normalized function is below min_dens."
      end if
      ! check if the integral of the normalized function is actually num_levels-1
      if (abs(calc_integral(gd_idx, gd_norm_x, gd_norm_y) - (num_levels-1)) > 0.000001) then
        write(fstderr,*) "Warning! The integral of the normalized function is not num_levels-1."
      end if
    end if

  end subroutine normalize_density

  function create_grid_from_normalized_grid_density_func( num_levels, &
                                                          gd_norm_idx, &
                                                          gd_norm_x, &
                                                          gd_norm_y ) &
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
      gd_norm_idx, &  ! The total numer of indices of gd_norm_x and gd_norm_y
      num_levels      ! number of levels the new grid should have

    real( kind = core_rknd ), dimension(gd_norm_idx), intent(in) ::  &
      gd_norm_x,  &    ! the x coordinates of the connection points of the normalized piecewise
                       ! linear grid density function
      gd_norm_y        ! the y coordinates of the connection points of the normalized piecewise
                       ! linear grid density function
      ! The grid density function needs to be normalized with the number of grid levels num_levels,
      ! before the function can be applied here

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top

    !--------------------- Local Variables ---------------------
    integer :: i, &
               prev_x_ind  ! the index of the last used x coordinate

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid
      grid_top, &               ! the highest point of the grid
      area_up_to_prev_x, &      ! the integral of the function up to the last used x coordinate
      new_x_area_increment, &   ! integral from grid density function on [prev_x_ind, prev_x_ind+1]
      desired_area_up_to_ith_level, & ! the desired area up to the ith grid level, should be exactly
                                      ! i, since the grid density function is normalized so that the
                                      ! integral of the whole function is num_levels-1 and since we
                                      ! follow the equi-distribution approach
      A, grid_level, m, b, p, q

    logical :: grid_level_found

    !--------------------- Begin Code ---------------------
    grid_sfc = gd_norm_x(1)
    grid_top = gd_norm_x(gd_norm_idx)

    grid_heights(1) = grid_sfc

    ! Initializing calc_integral to avoid a compiler warning.
    prev_x_ind = 1
    area_up_to_prev_x = 0.0_core_rknd

    do i = 2, num_levels-1
        grid_level_found = .false.
        do while (.not. grid_level_found)
            new_x_area_increment = (gd_norm_y(prev_x_ind) + gd_norm_y(prev_x_ind+1))/2 &
                                   * (gd_norm_x(prev_x_ind+1) - gd_norm_x(prev_x_ind))
            desired_area_up_to_ith_level = i-1
            if (area_up_to_prev_x + new_x_area_increment >= desired_area_up_to_ith_level) then
                A = desired_area_up_to_ith_level - area_up_to_prev_x
                ! check if gd_norm_y(prev_x_ind+1) == gd_norm_y(prev_x_ind)
                if (abs(gd_norm_y(prev_x_ind+1) - gd_norm_y(prev_x_ind)) < 1.0e-6_core_rknd) then
                    ! gd_norm_y(prev_x_ind+1) == gd_norm_y(prev_x_ind)
                    grid_level = A/gd_norm_y(prev_x_ind) + gd_norm_x(prev_x_ind)
                else
                    ! calculate slope and absolute value of linear function
                    m = (gd_norm_y(prev_x_ind+1) - gd_norm_y(prev_x_ind)) &
                        /(gd_norm_x(prev_x_ind+1) - gd_norm_x(prev_x_ind))
                    b = gd_norm_y(prev_x_ind) - gd_norm_x(prev_x_ind)*m
                    ! calculate p and q of second-order polynomial to solve for x
                    p = 2*b/m
                    q = -1*gd_norm_x(prev_x_ind)**2 - 2/m*b*gd_norm_x(prev_x_ind) - 2/m*A
                    ! since we have second-order polynomial, we have two roots, so we need to
                    ! identify the right solution - use pq formula to get solutions
                    grid_level = -1*p/2 + sqrt((p/2)**2 - q)

                    if ((grid_level < gd_norm_x(prev_x_ind)) &
                        .or. (grid_level > gd_norm_x(prev_x_ind+1)) & 
                        .or. (grid_level < grid_heights(i-1))) then
                        ! first solution was not the one we were looking for
                        grid_level = -1*p/2 - sqrt((p/2)**2 - q)
                        if ((grid_level < gd_norm_x(prev_x_ind)) &
                            .or. (grid_level > gd_norm_x(prev_x_ind+1)) & 
                            .or. (grid_level < grid_heights(i-1))) then
                            write(fstderr,*) "None of the two solutions works. ", &
                                             "Something went wrong."
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
      grid_heights_idx, &   ! number of elements in grid_heights
      desired_num_levels    ! desired number of grid levels
    ! Note: if everything worked fine, those two integers should be the same

    real( kind = core_rknd ), dimension(grid_heights_idx), intent(in) :: &
      grid_heights ! the grid heights ordered from bottom to top

    real( kind = core_rknd ), intent(in) :: &
      desired_min_dens, &   ! the desired minimum grid density
      desired_grid_sfc, &   ! desired surface height of the grid
      desired_grid_top      ! desired top height of the grid

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ) :: &
      max_dist  ! the maximum distance between two grid levels, depending on desired_min_dens

    !--------------------- Begin Code ---------------------
    ! check if first grid element is grid_sfc
    if (abs(grid_heights(1) - desired_grid_sfc) > 0.000001) then
        write(fstderr,*) "Warning! The first grid level is not correct. It should be the grids", &
                         " surface: ", desired_grid_sfc, " and not ", grid_heights(1)
    end if

    max_dist = 1/desired_min_dens

    do i = 1, (grid_heights_idx-1)
        ! check if grid heights are getting bigger
        if (grid_heights(i) >= grid_heights(i+1)) then
            write(fstderr,*) "Warning! The grid levels are not in the right order. ", &
                             "They should get bigger with the index."
        end if
        ! check if grid fulfills min_dense
        if ((grid_heights(i+1) - grid_heights(i)) > max_dist) then
            write(fstderr,*) "Warning! The grid has a distance, ", &
                             "that is bigger than defined by min_dens."
        end if
    end do

    ! check if grid has the correct number of levels
    if (grid_heights_idx /= desired_num_levels) then
        write(fstderr,*) "Warning! The grid has not the right number of levels. There are ", &
                         grid_heights_idx, " level, but there should be ", desired_num_levels, &
                         " level."
    end if

    ! check if last grid level is grid_top
    if (abs(grid_heights(desired_num_levels) - desired_grid_top) > 0.000001) then
        write(fstderr,*) "Warning! The last grid level is not correct. It should be the grids", &
                         " top: ", desired_grid_top, " and not ", grid_heights(desired_num_levels)
    end if

  end subroutine check_grid

  function create_grid_from_grid_density_func( gd_idx, &
                                               gd_x, gd_y, &
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
      gd_idx, &       ! total numer of indices of gd_x and gd_y
      num_levels      ! number of levels the new grid should have

    real( kind = core_rknd ), intent(in) ::  &
      min_dens  ! the minimum density the new grid should have

    real( kind = core_rknd ), dimension(gd_idx), intent(in) ::  &
      gd_x,  &    ! the x coordinates of the connection points of the piecewise linear
                  ! grid density function
      gd_y        ! the y coordinates of the connection points of the piecewise linear
                  ! grid density function

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx

    real( kind = core_rknd ), dimension(gd_idx) ::  &
      gd_norm_x,  &    ! the x coordinates of the connection points of the normalized piecewise
                       ! linear grid density function
      gd_norm_y        ! the y coordinates of the connection points of the normalized piecewise
                       ! linear grid density function

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface
      grid_top        ! height of the top of the grid

    !--------------------- Begin Code ---------------------
    grid_sfc = gd_x(1)
    grid_top = gd_x(gd_idx)

    call normalize_density( gd_idx, gd_x, gd_y, &
                            min_dens, num_levels, &
                            gd_norm_x, gd_norm_y )

    grid_heights = create_grid_from_normalized_grid_density_func( num_levels, &
                                                                  gd_idx, &
                                                                  gd_norm_x, &
                                                                  gd_norm_y )

    ! TODO how to read may change if we incorporate ngrdcol!,
    ! TODO since size does not work then, but should be size/ngrdcol?,
    ! TODO since all columns should have same number of grid levels?
    !grid_heights_idx = shape(grid_heights)
    grid_heights_idx = size(grid_heights) ! counts the scalar values

    if ( clubb_at_least_debug_level( 2 ) ) then
        call check_grid( grid_heights_idx, grid_heights, &
                         num_levels, min_dens, &
                         grid_sfc, grid_top )
    end if 

    return

  end function create_grid_from_grid_density_func

!===============================================================================

end module grid_adaptation_module