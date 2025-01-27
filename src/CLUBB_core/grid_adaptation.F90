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

  public :: setup_simple_gr_dycore, lin_interpolate, check_consistent, check_monotone, &
            remapping_matrix, remap, interpolate_forcings, check_mass_conservation, &
            check_conservation, calc_mass_over_grid_intervals, vertical_integral_conserve_mass, &
            check_remap_for_consistency

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

  function remapping_matrix( gr_target, gr_source, &
                             total_idx_rho_lin_spline, &
                             rho_lin_spline_vals, &
                             rho_lin_spline_levels )
    ! implements the remapping scheme proposed by Ullrich et al. in 
    !'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part II' 
    ! in formula (30), with the addition that omega incorporates in our case also rho, so it is not 
    ! the geometric region, but the mass

    implicit none
    !--------------------- Input Variables ---------------------
    type (grid) :: &
      gr_target, & ! grid where the values should be interpolated to
      gr_source    ! grid where the given values are currently on

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(gr_target%nzt,gr_source%nzt) :: &
      remapping_matrix ! matrix to apply to input values for remapping (R_ij)

    integer :: i, j

    real( kind = core_rknd ) :: omega_ov, omega_ov_upper, omega_ov_lower

    real( kind = core_rknd ), dimension(gr_target%nzt) :: &
      omega_ts ! densities of all intervals in the target grid

    real( kind = core_rknd ), dimension(2) :: &
      omega_ov_zm ! helper variable to hold density of the overlap region

    !--------------------- Begin Code ---------------------

    omega_ts = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              gr_target%nzm, gr_target%zm )
    do i = 1, gr_target%nzt
        do j = 1, gr_source%nzt
            omega_ov_upper = min(gr_target%zt(1,i)+0.5*gr_target%dzt(1,i), gr_source%zt(1,j)+0.5*gr_source%dzt(1,j))
            omega_ov_lower = max(gr_target%zt(1,i)-0.5*gr_target%dzt(1,i), gr_source%zt(1,j)-0.5*gr_source%dzt(1,j))
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
      call check_consistent( remapping_matrix, gr_target%nzt, gr_source%nzt )
      call check_monotone( remapping_matrix, gr_target%nzt, gr_source%nzt )
    end if

  end function remapping_matrix

  function remap( gr_target, gr_source, &
                  total_idx_rho_lin_spline, rho_lin_spline_vals, &
                  rho_lin_spline_levels, &
                  gr_source_values )

    implicit none
    !--------------------- Input Variables ---------------------
    type (grid) :: &
      gr_target, & ! grid where the values should be interpolated to
      gr_source    ! grid where the given values are currently on

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    real( kind = core_rknd ), dimension(gr_source%nzt) :: &
      gr_source_values  ! given values on the gr_source grid that should be interpolated 
                        ! to the gr_target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(gr_target%nzt) :: &
      remap ! interpolated values with the dimension of gr_target

    integer :: i, j

    real( kind = core_rknd ), dimension(gr_target%nzt,gr_source%nzt) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij)

    !--------------------- Begin Code ---------------------
    !$acc enter data create( remap, i, j, R_ij )

    R_ij = remapping_matrix( gr_target, gr_source, &
                             total_idx_rho_lin_spline, rho_lin_spline_vals, &
                             rho_lin_spline_levels )
    ! matrix vector multiplication
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, gr_target%nzt 
      remap(i) = 0
      do j = 1, gr_source%nzt
        remap(i) = remap(i) + R_ij(i,j)*gr_source_values(j)
      end do
    end do
    !$acc end parallel loop

  end function remap

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

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    do i=1, ngrdcol
        thlm_forcing(i,:) = remap( gr_clubb, gr_dycore, &
                                   total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                   rho_lin_spline_levels(i,:), &
                                   thlm_forcing_dycore(i,:) )
        rtm_forcing(i,:) = remap( gr_clubb, gr_dycore, &
                                  total_idx_rho_lin_spline, rho_lin_spline_vals(i,:), &
                                  rho_lin_spline_levels(i,:), &
                                  rtm_forcing_dycore(i,:) )
    end do
  end subroutine interpolate_forcings

  subroutine check_remap_for_consistency( gr_dycore, gr_clubb, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels )

    implicit none
    !--------------------- Input Variables ---------------------
    type (grid), intent(in) :: gr_dycore, gr_clubb

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(gr_dycore%nzt) :: &
      ones

    real( kind = core_rknd ), dimension(gr_clubb%nzt) :: &
      output_ones

    integer :: i, j

    real( kind = core_rknd ) :: epsilon

    logical :: consistent

    !--------------------- Begin Code ---------------------
    do i = 1, gr_dycore%nzt
      ones(i) = 1
    end do

    output_ones = remap( gr_clubb, gr_dycore, &
                         total_idx_rho_lin_spline, rho_lin_spline_vals, &
                         rho_lin_spline_levels, &
                         ones )

    consistent = .true.
    epsilon = 0.000001
    j = 1
    do while (consistent .and. j <= gr_clubb%nzt)
      if (abs( output_ones(j) - 1 ) > epsilon) then
        consistent = .false.
      end if
      j = j + 1
    end do

    if (.not. consistent) then
      write(*,*) 'remap function is not consistent'
    end if

  end subroutine check_remap_for_consistency

  subroutine check_mass_conservation( gr, gr_dycore, &
                                      total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                      rho_lin_spline_levels )

    implicit none

    !--------------------- Input Variables ---------------------
    type( grid ), intent(in) :: gr, gr_dycore

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, &
      sum_mass, &
      sum_dycore_mass

    real( kind = core_rknd ), dimension(gr%nzm-1) :: mass

    real( kind = core_rknd ), dimension(gr_dycore%nzm-1) :: dycore_mass

    !--------------------- Begin Code ---------------------
    tol = 0.000001

    mass = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          gr%nzm, gr%zm )

    dycore_mass = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                 rho_lin_spline_levels, &
                                                 gr_dycore%nzm, gr_dycore%zm )

    sum_mass = sum( mass )
    sum_dycore_mass = sum( dycore_mass )

    if (abs(sum_mass - sum_dycore_mass) > tol) then

      write(fstderr,*) "Mass for different grids was not conserved."
      write(fstderr,*) "Mass was ", sum_mass, " instead of ", sum_dycore_mass

    end if

  end subroutine check_mass_conservation

  subroutine check_conservation( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 total_field_zm_idx, total_field_zm_idx_dycore, &
                                 field_zm, field_zm_dycore, &
                                 field, field_dycore )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    integer, intent(in) :: &
      total_field_zm_idx, &       ! number of zm levels of the CLUBB grid
      total_field_zm_idx_dycore   ! number of zm levels of the dycore grid

    real( kind = core_rknd ), intent(in), dimension(total_field_zm_idx) :: &
      field_zm                  ! altitudes for the given field values on the CLUBB grid

    real( kind = core_rknd ), intent(in), dimension(total_field_zm_idx_dycore) :: &
      field_zm_dycore           ! altitudes for the given field values on the dycore grid

    real( kind = core_rknd ), intent(in), dimension(total_field_zm_idx-1) :: &
      field                  ! values for a variable on the CLUBB grid

    real( kind = core_rknd ), intent(in), dimension(total_field_zm_idx_dycore-1) :: &
      field_dycore           ! values for a variable on the dycore grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, &
      integral, &
      dycore_integral

    !--------------------- Begin Code ---------------------
    tol = 0.000001

    integral = vertical_integral_conserve_mass( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                total_field_zm_idx, field_zm, field )
    dycore_integral = vertical_integral_conserve_mass( total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       total_field_zm_idx_dycore, field_zm_dycore, &
                                                       field_dycore )

    if (abs(integral - dycore_integral) > tol) then

      write(fstderr,*) "Integral for field was not conserved."
      write(fstderr,*) "Integral was ", integral, " instead of ", dycore_integral

    end if

  end subroutine check_conservation

  function calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                          lin_spline_rho_vals, &
                                          lin_spline_rho_levels, &
                                          target_total_idx, target_zm &
                                        ) result( mass_per_interval )

    ! Description:
    ! Calculate the mass over every interval (zm(i) to zm(i+1)) of the zm levels (target_zm).
    ! This function assumes a linear spline for rho, defined by the values of 
    ! rho (lin_spline_rho_vals) given at the altitudes (lin_spline_rho_levels).
    ! So with the assumption that the linear spline is the exact rho function, this function 
    ! computes the exact mass over each interval, such that any target grid with the same lowest 
    ! and highest level would give the same total mass, if the linear spline is the same.

    ! Returns an array with the dimension target_total_idx-1 where at each index the mass of that 
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
      target_total_idx         ! The total numer of indices of the zm levels in the target grid

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    real( kind = core_rknd ), dimension(target_total_idx), intent(in) :: &
      target_zm ! zm levels of target grid over which the masses shoud be calculated

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      tol, & ! The difference tolerance for two real numbers to be considered the same
      val_below, & ! The rho value of the level or spline connection point below
      level_below, & ! The altitude of the level or spline connection point below
      mass_over_interval, & ! The total mass over one interval (zm(i) to zm(i+1)) in the grid
      upper_zm_level_rho

    real( kind = core_rknd ), dimension(target_total_idx-1) :: &
      mass_per_interval ! Array to store the mass of each interval (target_zm(i) to target_zm(i+1))

    integer :: i, j

    logical :: level_found

    !--------------------- Begin Code ---------------------
    tol = 0.000000001
    ! check if highest value of lin spline is higher or equal to highest level of target grid
    if ( (target_zm(target_total_idx) &
          - lin_spline_rho_levels(total_idx_lin_spline) ) > tol ) then
      write(*,*) 'Wrong input: highest level of lin_spline_rho_levels must', &
                 ' be higher or equal to highest zm level of target_zm'
    end if

    j = 1
    level_found = .false.
    ! find first level and set initial val_below and level_below
    do while ( (.not. level_found) .and. (j <= total_idx_lin_spline) )
      if ( abs( target_zm(1) - lin_spline_rho_levels(j) ) < tol ) then
        val_below = lin_spline_rho_vals(j)
        level_below = target_zm(1)
        level_found = .true.
      else if ( target_zm(1) < lin_spline_rho_levels(j) ) then ! j > 1
        if (j <= 1) then
          write(*,*) 'Wrong input: lowest level of lin_spline_rho_levels must', &
                     ' be lower or equal to lowest zm level of target_zm'
          call exit(1)
        else
          val_below = lin_interpolate_two_points( target_zm(1), &
                                                  lin_spline_rho_levels(j), &
                                                  lin_spline_rho_levels(j-1), &
                                                  lin_spline_rho_vals(j), &
                                                  lin_spline_rho_vals(j-1) )
          level_below = target_zm(1)
          level_found = .true.
          j = j - 1
        end if
      end if
      j = j + 1
    end do
    ! find all spline connection points between the level_below and the next zm level
    ! in target_zm and add up their masses, finally the rho value on the next zm level is
    ! interpolated and the last part of mass added to the total mass of that interval
    do i = 1, target_total_idx-1
      mass_over_interval = 0.0_core_rknd
      do while ( (lin_spline_rho_levels(j) < target_zm(i+1)) &
                 .and. (j < total_idx_lin_spline) )
        mass_over_interval = mass_over_interval &
                             + (lin_spline_rho_levels(j) - level_below) &
                              *( val_below + lin_spline_rho_vals(j) )/2
        level_below = lin_spline_rho_levels(j)
        val_below = lin_spline_rho_vals(j)
        j = j + 1
      end do
      if ( abs( target_zm(i+1) - lin_spline_rho_levels(j) ) < tol ) then 
        ! have approx equal altitude, so just copy value
        upper_zm_level_rho = lin_spline_rho_vals(j)
      else ! lin spline at j is above zm(i+1) in target_zm
        ! lin interpolation between j lin spline node and lin spine node j-1
        upper_zm_level_rho = lin_interpolate_two_points( target_zm(i+1), &
                                                         lin_spline_rho_levels(j), &
                                                         lin_spline_rho_levels(j-1), &
                                                         lin_spline_rho_vals(j), &
                                                         lin_spline_rho_vals(j-1) )
      end if
      mass_over_interval = mass_over_interval &
                           + (target_zm(i+1) - level_below) &
                            *( val_below + upper_zm_level_rho )/2

      mass_per_interval(i) = mass_over_interval

      val_below = upper_zm_level_rho
      level_below = target_zm(i+1)
    end do
    return
  end function calc_mass_over_grid_intervals

  function vertical_integral_conserve_mass( total_idx_lin_spline, lin_spline_rho_vals, &
                                            lin_spline_rho_levels, &
                                            total_field_zm_idx, field_zm, field )

    ! Description:
    ! Computes the vertical integral. lin_spline_rho_vals and lin_spline_rho_levels must be
    ! of size total_idx_lin_spline and should start at the same index.
    ! field must be of size total_field_zm_idx-1
    ! For more conditions to thos parameters read description of calc_mass_over_grid_intervals

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
      total_field_zm_idx  ! The total numer of indices for field_zm within the range of averaging

    real( kind = core_rknd ), dimension(total_field_zm_idx-1), intent(in) ::  &
      field      ! The field to be vertically averaged   [Units vary]

    real( kind = core_rknd ), dimension(total_field_zm_idx), intent(in) ::  &
      field_zm

    ! Note:  The field and field_zm points need to be arranged from
    !        lowest to highest in altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(total_field_zm_idx-1) :: &
      mass_per_interval ! Array to store the mass of each interval of the target_grid

    real( kind = core_rknd ) :: &
      vertical_integral_conserve_mass

    !--------------------- Begin Code ---------------------

    ! Initializing vertical_integral_conserve_mass to avoid a compiler warning.
    vertical_integral_conserve_mass = 0.0_core_rknd

    ! Compute the integral.
    ! Multiply the field at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    ! Note:  The values of the field and rho_ds are passed into this function
    !        so that field(1) and rho_ds(1) are actually the field and rho_ds
    !        at level k_start.
    mass_per_interval = calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                                       lin_spline_rho_vals, &
                                                       lin_spline_rho_levels, &
                                                       total_field_zm_idx, field_zm )
    vertical_integral_conserve_mass = sum( field * mass_per_interval )

    return
  end function vertical_integral_conserve_mass

!===============================================================================

end module grid_adaptation_module