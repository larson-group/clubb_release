!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module remapping_module

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  use constants_clubb, only: &
      fstderr ! Constants

  implicit none

  public :: remap_vals_to_target

  private :: calc_mass_over_grid_intervals, matrix_vector_mult, vertical_integral_conserve_mass, &
             check_conservation, check_monotonicity, check_consistency, &
             remap_vals_to_target_helper, &
             check_remap_matrix_conservation, check_remap_matrix_consistency, &
             check_remap_matrix_monotonicity, remapping_matrix, &
             kmppm, steepz, ppm2m, map1_ppm, remap_vals_ppm

  private

  real( kind = core_rknd ), parameter ::  &
    tol = 1.0e-6_core_rknd ! tolerance to check if two real numbers are equal

  contains

  !=============================================================================

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !
  !
  !                GENERAL REMAPPING FUNCTIONS AND SUBROUTINES
  !
  !
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  function calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                          lin_spline_rho_vals, &
                                          lin_spline_rho_levels, &
                                          grid_levels_idx, grid_levels &
                                        ) result( mass_per_interval )

    ! Description:
    ! Calculate the mass over every interval (grid_levels(k) to grid_levels(k+1))
    ! of the grid levels (grid_levels).
    ! This function assumes a linear spline for rho, defined by the values of 
    ! rho (lin_spline_rho_vals) given at the altitudes (lin_spline_rho_levels).
    ! So with the assumption that the linear spline is the exact rho function, this function 
    ! computes the exact mass over each interval, such that any target grid with the same lowest 
    ! and highest level would give the same total mass, if the linear spline is the same.

    ! Returns an array with the dimension grid_levels_idx-1 where at each index the mass of that 
    ! interval is stored, starting at the bottom level.
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use interpolation, only: &
        lin_interpolate_two_points ! Procedure

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      total_idx_lin_spline, &  ! The total numer of indices of the spline points []
      grid_levels_idx          ! The total numer of indices of the grid levels   []

    real( kind = core_rknd ), dimension(total_idx_lin_spline), intent(in) ::  &
      lin_spline_rho_vals,  &    ! Dry, static density                   [kg/m^3]
      lin_spline_rho_levels      ! The altitudes of the given rho values [m]
    ! Note:  The lin_spline_rho_levels and lin_spline_rho_vals need to be arranged from
    !        lowest to highest in altitude

    real( kind = core_rknd ), dimension(grid_levels_idx), intent(in) :: &
      grid_levels ! levels of grid over which the masses shoud be calculated [m]

    !--------------------- Output Variable ----------------------
    real( kind = core_rknd ), dimension(grid_levels_idx-1) :: &
      mass_per_interval ! Array to store the mass of each interval
                        ! (grid_levels(k) to grid_levels(k+1)) [kg]

    !--------------------- Local Variables ----------------------
    real( kind = core_rknd ) :: &
      val_below, &          ! The rho value of the level or spline connection point below  [kg/m^3]
      level_below, &        ! The altitude of the level or spline connection point below   [m]
      mass_over_interval, & ! The total mass over one grid cell                            [kg]
      upper_zm_level_rho    ! The rho value of the upper level of the grid cell            [m]

    integer :: k, l

    logical :: l_level_found

    !--------------------- Begin Code ---------------------
    ! check if highest value of lin spline is higher or equal to highest level of target grid
    if ( (grid_levels(grid_levels_idx) &
          - lin_spline_rho_levels(total_idx_lin_spline) ) > tol ) then
      write(fstderr,*) 'Wrong input: highest level of lin_spline_rho_levels must', &
                       ' be higher or equal to highest zm level of grid_levels'
      write(fstderr,*) 'Cannot compute the mass to get the remapping operator if the', & 
                       'density function is not defined over the whole grid region'
      error stop 'Cannot compute the remapping operator'
    end if

    l = 1
    l_level_found = .false.
    ! find first level and set initial val_below and level_below
    do while ( (.not. l_level_found) .and. (l <= total_idx_lin_spline) )
      if ( abs( grid_levels(1) - lin_spline_rho_levels(l) ) < tol ) then
        val_below = lin_spline_rho_vals(l)
        level_below = grid_levels(1)
        l_level_found = .true.
      else if ( grid_levels(1) < lin_spline_rho_levels(l) ) then ! l > 1
        if (l <= 1) then
          write(fstderr,*) 'Wrong input: lowest level of lin_spline_rho_levels must', &
                           ' be lower or equal to lowest zm level of grid_levels'
          write(fstderr,*) 'Cannot compute the mass to get the remapping operator if the', & 
                           'density function is not defined over the whole grid region'
          error stop 'Cannot compute the remapping operator'
        else
          val_below = lin_interpolate_two_points( grid_levels(1), &
                                                  lin_spline_rho_levels(l), &
                                                  lin_spline_rho_levels(l-1), &
                                                  lin_spline_rho_vals(l), &
                                                  lin_spline_rho_vals(l-1) )
          level_below = grid_levels(1)
          l_level_found = .true.
          l = l - 1
        end if
      end if
      l = l + 1
    end do
    ! find all spline connection points between the level_below and the next grid level
    ! in grid_levels and add up their masses, finally the rho value on the next grid level is
    ! interpolated and the last part of mass added to the total mass of that interval
    do k = 1, grid_levels_idx-1
      mass_over_interval = 0.0_core_rknd
      do while ( (lin_spline_rho_levels(l) < grid_levels(k+1)) &
                 .and. (l < total_idx_lin_spline) )
        mass_over_interval = mass_over_interval &
                             + (lin_spline_rho_levels(l) - level_below) &
                              *( val_below + lin_spline_rho_vals(l) )/2
        level_below = lin_spline_rho_levels(l)
        val_below = lin_spline_rho_vals(l)
        l = l + 1
      end do
      if ( abs( grid_levels(k+1) - lin_spline_rho_levels(l) ) < tol ) then 
        ! have approx equal altitude, so just copy value
        upper_zm_level_rho = lin_spline_rho_vals(l)
      else ! lin spline at l is above grid_levels(k+1) in grid_levels
        ! lin interpolation between l lin spline node and lin spine node l-1
        upper_zm_level_rho = lin_interpolate_two_points( grid_levels(k+1), &
                                                         lin_spline_rho_levels(l), &
                                                         lin_spline_rho_levels(l-1), &
                                                         lin_spline_rho_vals(l), &
                                                         lin_spline_rho_vals(l-1) )
      end if
      mass_over_interval = mass_over_interval &
                           + (grid_levels(k+1) - level_below) &
                            *( val_below + upper_zm_level_rho )/2

      mass_per_interval(k) = mass_over_interval

      val_below = upper_zm_level_rho
      level_below = grid_levels(k+1)
    end do

  end function calc_mass_over_grid_intervals

  function matrix_vector_mult( ngrdcol, &
                               dim_input, dim_output, &
                               x_vectors, &
                               A_matrices ) result ( y_vectors )

    ! Description:
    ! Calculates the matrix-vector product y = Ax, where A is a matrix with dim_output x dim_input
    ! and y,x are vectors with x having dim_input many entries and y having dim_output many entries.
    ! The result y is calculated for every ngrdcol

    implicit none
    
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      dim_input, & ! dimension of the input vector []
      dim_output   ! dimension of the output vector []

    real( kind = core_rknd ), dimension(ngrdcol,dim_input), intent(in) :: &
      x_vectors  ! values of the x vectors (ngrdcol many vectors)

    real( kind = core_rknd ), dimension(ngrdcol,dim_output,dim_input), intent(in) :: &
      A_matrices ! values of the matrices (ngrdcol many matrices)

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,dim_output) :: &
      y_vectors ! values of the y vectors (ngrdcol many vectors)

    !--------------------- Local Variables ---------------------
    integer :: i, k, j

    !--------------------- Begin Code ---------------------

    do i = 1, ngrdcol
      ! matrix vector multiplication
      do k = 1, dim_output
        y_vectors(i,k) = 0
        do j = 1, dim_input
          y_vectors(i,k) = y_vectors(i,k) + A_matrices(i,k,j)*x_vectors(i,j)
        end do
      end do
    end do

  end function matrix_vector_mult

  function vertical_integral_conserve_mass( total_idx_lin_spline, &
                                            lin_spline_rho_vals, &
                                            lin_spline_rho_levels, &
                                            grid_levels_idx, grid_levels, &
                                            field )

    ! Description:
    ! Computes the vertical integral. lin_spline_rho_vals and lin_spline_rho_levels must be
    ! of size total_idx_lin_spline and should start at the same index.
    ! field must be of size grid_levels_idx-1
    ! For more conditions to those parameters read description of calc_mass_over_grid_intervals

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero ! Constant

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
      grid_levels_idx  ! The total numer of indices for grid_levels []

    real( kind = core_rknd ), dimension(grid_levels_idx-1), intent(in) ::  &
      field      ! The field to be vertically averaged   [Units vary]

    real( kind = core_rknd ), dimension(grid_levels_idx), intent(in) ::  &
      grid_levels ! [m]

    ! Note:  The field and grid_levels points need to be arranged from
    !        lowest to highest in altitude

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(grid_levels_idx-1) :: &
      mass_per_interval ! Array to store the mass of each interval of the target_grid [kg]

    real( kind = core_rknd ) :: &
      vertical_integral_conserve_mass ! [kg]

    !--------------------- Begin Code ---------------------
    ! Initializing vertical_integral_conserve_mass to avoid a compiler warning.
    vertical_integral_conserve_mass = zero

    ! Compute the integral.
    ! Multiply the field at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    mass_per_interval = calc_mass_over_grid_intervals( total_idx_lin_spline, &
                                                       lin_spline_rho_vals, &
                                                       lin_spline_rho_levels, &
                                                       grid_levels_idx, grid_levels )
    vertical_integral_conserve_mass = sum( field * mass_per_interval )

  end function vertical_integral_conserve_mass

  subroutine check_conservation( total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 levels_source_idx, levels_target_idx, &
                                 levels_source, levels_target, &
                                 field_source, field_target )

    ! Description:
    ! Check if the vertical integral is conserved during remapping. It is checked, if the vertical
    ! integral over the remapped values is the same as the vertical integral over the source
    ! values. The vertical integral is exactly calculated, assuming the rho values build a
    ! piecewise linear function. 

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude

    integer, intent(in) :: &
      levels_source_idx, &       ! number of levels of the source grid []
      levels_target_idx          ! number of levels of the target grid []

    real( kind = core_rknd ), intent(in), dimension(levels_source_idx) :: &
      levels_source           ! altitudes for the source grid [m]

    real( kind = core_rknd ), intent(in), dimension(levels_target_idx) :: &
      levels_target           ! altitudes for the target grid [m]

    real( kind = core_rknd ), intent(in), dimension(levels_source_idx-1) :: &
      field_source            ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(levels_target_idx-1) :: &
      field_target           ! values for a variable on the target grid

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ) :: &
      integral_source, &  ! vertical integral over the source values
      integral_target, &  ! vertical integral over the target values
      tol_cons,        &  ! fixed tolerance for conservation
      err_percentage      ! the percentage the vertical integral of the remapped values is
                          ! allowed to deviate from the original vertical inetgral due to
                          ! numerical inaccuracies (with [0,1])

    !--------------------- Begin Code ---------------------
    integral_source = vertical_integral_conserve_mass( total_idx_rho_lin_spline, &
                                                       rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       levels_source_idx, levels_source, &
                                                       field_source )
 
    integral_target = vertical_integral_conserve_mass( total_idx_rho_lin_spline, &
                                                       rho_lin_spline_vals, &
                                                       rho_lin_spline_levels, &
                                                       levels_target_idx, levels_target, &
                                                       field_target )

    tol_cons = 1.e-3_core_rknd
    err_percentage = 1.e-7_core_rknd
    
    ! not use fixed global tolerance, since there are some extremely large vertical integrals
    if ( abs( integral_target - integral_source ) > abs( integral_source*err_percentage ) .and. &
         abs( integral_target - integral_source ) > tol_cons ) then

      write(fstderr,*) "WARNING! Integral for field was not conserved."
      write(fstderr,*) "Integral was ", integral_target, " instead of ", integral_source
    end if

  end subroutine check_conservation

  subroutine check_monotonicity( field_source_idx, field_target_idx, &
                                 field_source, field_target, &
                                 grid_remap_method )

    ! Description:
    ! Check if remapped values are globally monotone. So, check if all remapped values on the
    ! target grid are within the range of the minimum and maximum of the source values.

    use model_flags, only: &
        ppm_remap

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      field_source_idx, & ! dimension of field_source []
      field_target_idx    ! dimension of field_target []

    real( kind = core_rknd ), intent(in), dimension(field_source_idx) :: &
      field_source  ! values for a variable on the source grid

    real( kind = core_rknd ), intent(in), dimension(field_target_idx) :: &
      field_target  ! remapped values for a variable on the target grid

    integer, intent(in) :: &
      grid_remap_method ! integer specifying which remapping method was used

    !--------------------- Local Variables ---------------------
    integer :: k, &
               start_index, & ! start index to check target values for monotonicity
               end_index      ! end index to check target values for monotonicity

    real( kind = core_rknd ) :: &
      source_maxval, & ! maximum value of the source field values
      source_minval    ! minimum value of the source field values

    !--------------------- Begin Code ---------------------
    start_index = 1
    end_index = field_target_idx

    if ( grid_remap_method == ppm_remap ) then
      ! Restrict the monotonicity test to the inner levels if PPM is used
      start_index = 4
      end_index = field_target_idx-3
    end if

    source_maxval = maxval( field_source )
    source_minval = minval( field_source )

    do k = start_index, end_index
      ! check if remapped target field is outside the region the global extrema
      ! of the source values and if deviation is greater than the tolerance
      if ( ( field_target(k) < source_minval &
             .and. abs(field_target(k) - source_minval ) > tol ) &
           .or. &
           ( field_target(k) > source_maxval &
             .and. abs(field_target(k) - source_maxval ) > tol ) ) then
        write(fstderr,*) 'WARNING! The remapped values are not monotone.'
        write(fstderr,*) 'The extrema on the source field were: ', source_minval, ' and ', &
                          source_maxval, ', but the remapped value is: ', field_target(k)
      end if
    end do

  end subroutine check_monotonicity

  subroutine check_consistency( dim_input, dim_output, &
                                R_ij )

    ! Description:
    ! Check if matrix used to remap values remaps values in a consistent way. So, check if
    ! the matrix is used to remap a vector of ones, the remapped values are also ones.

    use constants_clubb, only: &
        one

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      dim_input, & ! dimension of the input values on the source grid []
      dim_output   ! dimension of the reampped values on the target grid []

    real( kind = core_rknd ), dimension(dim_output,dim_input), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(dim_input) :: &
      ones ! array with ones

    real( kind = core_rknd ), dimension(1,dim_output) :: &
      output_ones ! remapped array of ones

    integer :: k

    logical :: l_consistent

    !--------------------- Begin Code ---------------------
    do k = 1, (dim_input)
      ones(k) = one
    end do

    output_ones = matrix_vector_mult( 1, &
                                      dim_input, dim_output, &
                                      ones, &
                                      R_ij )

    l_consistent = .true.
    k = 1
    do while (l_consistent .and. k <= dim_output)
      if (abs( output_ones(1,k) - one ) > tol) then
        l_consistent = .false.
      end if
      k = k + 1
    end do

    if (.not. l_consistent) then
      write(fstderr,*) 'WARNING! The matrix_vector_mult function is not consistent'
      write(fstderr,*) 'The output should be all ones, but was: ', output_ones(1,:)
    end if

  end subroutine check_consistency

  function remap_vals_to_target_helper( ngrdcol, &
                                        levels_source_idx, levels_target_idx, &
                                        levels_source, levels_target, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        source_values, &
                                        iv, p_sfc, &
                                        grid_remap_method, &
                                        R_ij ) result( remap_vals_to_target )

    use constants_clubb, only: &
        grav ! Constant

    use model_flags, only: &
        cons_ullrich_remap, &
        ppm_remap ! Constants

    use error_code, only: &
        clubb_at_least_debug_level
    
    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      levels_source_idx, & ! number of levels in the target grid []
      levels_target_idx    ! number of levels in the source grid []

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,levels_source_idx) :: &
      levels_source           ! altitudes for the source grid [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,levels_target_idx) :: &
      levels_target           ! altitudes for the target grid [m]

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! dimension of the linear spline definition arrays for rho []

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes  [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude
    !       rho_lin_spline_vals, rho_lin_spline_levels build an array of (x,y) points that are
    !       used to construct a piecewise linear function for rho, such that the exact
    !       integral can be calculated

    real( kind = core_rknd ), dimension(ngrdcol,levels_source_idx-1), intent(in) :: &
      source_values  ! given values on the source grid that should be remapped 
                     ! to the target grid
    
    integer, intent(in) :: &
      iv ! setting for E3SM's PPM method, that sets what type of value is remapped,
         ! the following are the three valid settings:
         ! iv =-1: winds
         ! iv = 0: positive definite scalars
         ! iv = 1: others

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc ! prescribed pressure at the surface [Pa]

    integer, intent(in) :: &
      grid_remap_method ! specifies which remapping method should be used
                        ! options: 1. Ullrich linear remapping scheme (first order)
                        !          2. Piecewise Parabolic Method (PPM) (high order)

    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx-1,levels_source_idx-1), &
    intent(in), optional :: &
      R_ij      ! remapping matrix for Ullrich linear remapping scheme
                ! if remapping matrix was already calculated before it can be used again,
                ! to save computation time, otherwise it will be calculated

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx-1) :: &
      remap_vals_to_target ! remapped values with the dimension of the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k

    real( kind = core_rknd ), dimension(ngrdcol,levels_source_idx-1) :: &
      mass_on_source_cells ! mass over every source grid cell [kg]

    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx-1) :: &
      mass_on_target_cells ! mass over every target grid cell [kg]

    real( kind = core_rknd ), dimension(ngrdcol,levels_source_idx) :: &
      pressure_levels_source ! pressure levels of the source grid [Pa]

    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx) :: &
      pressure_levels_target ! pressure levels of the target grid [Pa]

    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx-1,levels_source_idx-1) :: &
      R_ij_internal   ! this variable is used for the remapping matrix in this function to also
                      ! set it if optional values was not present

    !--------------------- Begin Code ---------------------

    do i = 1, ngrdcol

      ! Calculate pressure levels from altitudes and densities
      mass_on_source_cells(i,:) = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 levels_source_idx, levels_source )
      mass_on_target_cells(i,:) = calc_mass_over_grid_intervals( total_idx_rho_lin_spline, &
                                                                 rho_lin_spline_vals, &
                                                                 rho_lin_spline_levels, &
                                                                 levels_target_idx, levels_target )
      pressure_levels_source(i,1) = p_sfc(i)
      do k = 2, levels_source_idx
        pressure_levels_source(i,k) = pressure_levels_source(i,k-1) &
                                      - mass_on_source_cells(i,k-1)*grav
      end do
      pressure_levels_target(i,1) = p_sfc(i)
      do k = 2, levels_target_idx
        pressure_levels_target(i,k) = pressure_levels_target(i,k-1) &
                                      - mass_on_target_cells(i,k-1)*grav
      end do

      ! Construct remapping matrix for cons_ullrich_remap method, if not
      ! already given as an optional parameter
      if ( grid_remap_method == cons_ullrich_remap .and. .not. present( R_ij ) ) then
        R_ij_internal(i,:,:) = remapping_matrix( levels_source_idx, levels_target_idx, &
                                                 pressure_levels_source(i,:), &
                                                 pressure_levels_target(i,:) )

      end if

    end do

    ! Remap values with the selected remapping method
    if ( grid_remap_method == cons_ullrich_remap ) then

      if ( present( R_ij ) ) then
        remap_vals_to_target = matrix_vector_mult( ngrdcol, &
                                                   levels_source_idx-1, &
                                                   levels_target_idx-1, &
                                                   source_values, &
                                                   R_ij )
      else
        remap_vals_to_target = matrix_vector_mult( ngrdcol, &
                                                   levels_source_idx-1, &
                                                   levels_target_idx-1, &
                                                   source_values, &
                                                   R_ij_internal )
      end if

    else if ( grid_remap_method == ppm_remap ) then

      remap_vals_to_target = remap_vals_ppm( ngrdcol, &
                                             levels_source_idx, levels_target_idx, &
                                             pressure_levels_source, pressure_levels_target, &
                                             source_values, iv )

    else 

        write(fstderr,*) 'There is no remapping method implemented for grid_remap_method=', &
                         grid_remap_method
        error stop 'Please adjust the configurable model flags.'
    
    end if

    ! Check remapped values
    if ( clubb_at_least_debug_level( 2 ) ) then
      do i = 1, ngrdcol

        call check_conservation( total_idx_rho_lin_spline, &                ! In
                                 rho_lin_spline_vals(i,:), &                ! In
                                 rho_lin_spline_levels(i,:), &              ! In
                                 levels_source_idx, levels_target_idx, &    ! In
                                 levels_source(i,:), levels_target(i,:), &  ! In
                                 source_values(i,:), &                      ! In
                                 remap_vals_to_target(i,:) )                ! In

        call check_monotonicity( levels_source_idx-1, levels_target_idx-1, & ! In
                                 source_values(i,:), &                       ! In
                                 remap_vals_to_target(i,:), &                ! In
                                 grid_remap_method )                         ! In

        if ( grid_remap_method == cons_ullrich_remap ) then
          if ( present( R_ij ) ) then
            call check_consistency( levels_source_idx-1, levels_target_idx-1, &  ! In
                                    R_ij(i,:,:) )                                ! In
          else
            call check_consistency( levels_source_idx-1, levels_target_idx-1, &  ! In
                                    R_ij_internal(i,:,:) )                       ! In
          end if
        end if
                                 
      end do
    end if

  end function remap_vals_to_target_helper

  function remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 source_values_idx, &
                                 source_values, &
                                 target_values_idx, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )

    use grid_class, only: &
        grid ! Type 
    
    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      source_values_idx, & ! dimension of source_values []
      target_values_idx    ! dimension of remapped source_values on target grid []

    type( grid ) :: &
      gr_source, & ! source grid object, where source_values are given
      gr_target    ! target grid object, where source_values are remapped to

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays []

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes [kg/m^3]
      rho_lin_spline_levels  ! altitudes for the given rho values [m]
    ! Note: both these arrays need to be sorted from low to high altitude
    !       rho_lin_spline_vals, rho_lin_spline_levels build an array of (x,y) points that are
    !       used to construct a piecewise linear function for rho, such that the exact
    !       integral can be calculated

    real( kind = core_rknd ), dimension(ngrdcol,source_values_idx), intent(in) :: &
      source_values  ! given values on the source grid that should be remapped 
                     ! to the target grid
    
    integer, intent(in) :: &
      iv ! setting for E3SM's PPM method, that sets what type of value is remapped,
         ! the following are the three valid settings:
         ! iv =-1: winds
         ! iv = 0: positive definite scalars
         ! iv = 1: others

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc ! prescribed pressure at the surface [Pa]

    integer, intent(in) :: &
      grid_remap_method ! specifies which remapping method should be used
                        ! options: 1. Ullrich linear remapping scheme (first order)
                        !          2. Piecewise Parabolic Method (PPM) (high order)

    logical, intent(in) :: &
      l_zt_variable ! if .true., the variable that was input is given on the zt levels of
                    ! the grid, if it is .false., it is given on the zm levels

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,target_values_idx) :: &
      remap_vals_to_target ! remapped values with the dimension of the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k

    real( kind = core_rknd ), dimension(ngrdcol,source_values_idx+1) :: &
      levels_source           ! altitudes for the source grid [m]

    real( kind = core_rknd ), dimension(ngrdcol,target_values_idx+1) :: &
      levels_target           ! altitudes for the target grid [m]

    !--------------------- Begin Code ---------------------

    if ( l_zt_variable ) then

      ! if variable is given on the zt levels, just take the surrounding zm levels to build
      ! the corresponding grid cell
      remap_vals_to_target = remap_vals_to_target_helper( ngrdcol, &
                                                          gr_source%nzm, gr_target%nzm, &
                                                          gr_source%zm, gr_target%zm, &
                                                          total_idx_rho_lin_spline, &
                                                          rho_lin_spline_vals, &
                                                          rho_lin_spline_levels, &
                                                          source_values, &
                                                          iv, p_sfc, &
                                                          grid_remap_method )

    else if ( .not. l_zt_variable ) then

      ! construct surrounding grid levels for values given on the zm levels
      ! for inner levels take the surrounding zt levels, for the outer levels,
      ! take the zt level and the zm level
      do i = 1, ngrdcol
        levels_source(i,1) = gr_source%zm(i,1)
        levels_target(i,1) = gr_target%zm(i,1)
        do k = 1, gr_source%nzt
          levels_source(i,k+1) = gr_source%zt(i,k)
        end do
        do k = 1, gr_target%nzt
          levels_target(i,k+1) = gr_target%zt(i,k)
        end do
        levels_source(i,gr_source%nzt+2) = gr_source%zm(i,gr_source%nzm)
        levels_target(i,gr_target%nzt+2) = gr_target%zm(i,gr_target%nzm)
      end do

      remap_vals_to_target = remap_vals_to_target_helper( ngrdcol, &
                                                          gr_source%nzt+2, gr_target%nzt+2, &
                                                          levels_source, levels_target, &
                                                          total_idx_rho_lin_spline, &
                                                          rho_lin_spline_vals, &
                                                          rho_lin_spline_levels, &
                                                          source_values, &
                                                          iv, p_sfc, &
                                                          grid_remap_method )

    end if

  end function remap_vals_to_target

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !
  !                LINEAR REMAPPING METHOD BY ULLRICH ET AL.
  !        formula (30) in 'Arbitrary-Order Conservative and Consistent
  !                         Remapping and a Theory of Linear Maps: Part II'
  !
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  subroutine check_remap_matrix_conservation( R_ij, dim1, dim2, &
                                              levels_source, levels_target )
    ! checks conservation with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (9)
    ! Defintion of conservation by Ullrich and Taylor: A remapping operator R is conservative if
    ! the global mass of any field is maintained across the remapping operation

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2 ! dimensions of the remapping matrix []

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    real( kind = core_rknd ), dimension(dim2+1), intent(in) :: &
      levels_source ! the pressure levels in the source grid [kg/(ms^2)]

    real( kind = core_rknd ), dimension(dim1+1), intent(in) :: &
      levels_target ! the pressure levels in the target grid [kg/(ms^2)]

    !--------------------- Local Variables ---------------------
    integer :: k, l

    real( kind = core_rknd ) :: &
      weighted_col_sum

    logical :: l_conservative

    !--------------------- Begin Code ---------------------
    l_conservative = .true.

    l = 1
    do while (l_conservative .and. l <= dim2)
      weighted_col_sum = 0
      do k = 1, dim1
        weighted_col_sum = weighted_col_sum + R_ij(k,l) * abs(levels_target(k+1) - levels_target(k))
      end do
      if (abs(weighted_col_sum - abs(levels_source(l+1) - levels_source(l))) > tol) then
        l_conservative = .false.
      end if
      l = l + 1
    end do

    if (.not. l_conservative) then
      write(fstderr,*) 'WARNING! The remapping operator is not conservative'
      error stop 'Operator should be conservative, something went wrong...'
    end if

  end subroutine check_remap_matrix_conservation

  subroutine check_remap_matrix_consistency( R_ij, dim1, dim2 )
    ! checks consistency with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (12)
    ! Defintion of consistency by Ullrich and Taylor: A remapping operator R is consistent if the 
    ! constant field is maintained across the remapping operation
    ! So for example, if we have a discrete function of only ones and apply the remapping operator, 
    ! then the result should also be constant one

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2 ! dimension of the remapping matrix []

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: k, l

    real( kind = core_rknd ) :: row_sum

    logical :: l_consistent

    !--------------------- Begin Code ---------------------
    l_consistent = .true.

    k = 1
    do while (l_consistent .and. k <= dim1)
      row_sum = 0
      do l = 1, dim2
        row_sum = row_sum + R_ij(k,l)
      end do
      if (abs(row_sum - 1) > tol) then
        l_consistent = .false.
      end if
      k = k + 1
    end do

    if (.not. l_consistent) then
      write(fstderr,*) 'WARNING! The remapping operator is not consistent'
      error stop 'Operator should be consistent, something went wrong...'
    end if

  end subroutine check_remap_matrix_consistency

  subroutine check_remap_matrix_monotonicity( R_ij, dim1, dim2 )
    ! checks monotonicity with condition defined by Ullrich and Taylor in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part I' 
    ! in formula (14)
    ! Defintion of monotonicity by Ullrich and Taylor: A remapping operator R is monotone if the 
    ! remapping operation cannot introduce additional global extrema

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: dim1, dim2 ! dimensions of the remapping matrix []

    real( kind = core_rknd ), dimension(dim1,dim2), intent(in) :: &
      R_ij  ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: k, l

    logical :: l_monotone

    !--------------------- Begin Code ---------------------
    l_monotone = .true.

    k = 1
    do while (l_monotone .and. k <= dim1)
      do l = 1, dim2
        if (R_ij(k,l) < 0) then
          l_monotone = .false.
        end if
      end do
      k = k + 1
    end do

    if (.not. l_monotone) then
      write(fstderr,*) 'WARNING! The remapping operator is not monotone'
      error stop 'Operator should be monotone, something went wrong...'
    end if

  end subroutine check_remap_matrix_monotonicity

  function remapping_matrix( levels_source_idx, levels_target_idx, &
                             levels_source, levels_target )
    ! implements the remapping scheme proposed by Ullrich et al. in 
    ! 'Arbitrary-Order Conservative and Consistent Remapping and a Theory of Linear Maps: Part II' 
    ! in formula (30)

    use error_code, only: &
        clubb_at_least_debug_level

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      levels_source_idx, & ! number of levels in the target grid []
      levels_target_idx    ! number of levels in the source grid []

    real( kind = core_rknd ), dimension(levels_source_idx), intent(in) :: &
      levels_source ! the pressure levels in the source grid [kg/(ms^2)]

    real( kind = core_rknd ), dimension(levels_target_idx), intent(in) :: &
      levels_target ! the pressure levels in the target grid [kg/(ms^2)]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(levels_target_idx-1,levels_source_idx-1) :: &
      remapping_matrix ! matrix to apply to input values for remapping (R_ij) []

    !--------------------- Local Variables ---------------------
    integer :: k, l

    real( kind = core_rknd ) :: omega_ov, omega_ov_upper, omega_ov_lower

    !--------------------- Begin Code ---------------------

    do k = 1, (levels_target_idx-1)
        do l = 1, (levels_source_idx-1)
            omega_ov_upper = max(levels_target(k+1), levels_source(l+1))
            omega_ov_lower = min(levels_target(k), levels_source(l))
            if (omega_ov_upper < omega_ov_lower) then
                omega_ov = abs(omega_ov_upper - omega_ov_lower)
            else
                omega_ov = 0
            end if
            remapping_matrix(k,l) = omega_ov/abs(levels_target(k+1) - levels_target(k))
        end do
    end do

    if ( clubb_at_least_debug_level( 2 ) ) then
      call check_remap_matrix_conservation( remapping_matrix, &
                                            (levels_target_idx-1), (levels_source_idx-1), &
                                            levels_source, levels_target )
      call check_remap_matrix_consistency( remapping_matrix, &
                                           (levels_target_idx-1), (levels_source_idx-1) )
      call check_remap_matrix_monotonicity( remapping_matrix, &
                                            (levels_target_idx-1), (levels_source_idx-1) )
    end if

  end function remapping_matrix

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !
  !                PIECEWISE PARABOLIC METHOD BY ULLRICH AND WOODWARD
  !        in 'The Piecewise Parabolic Method (PPM) for gas-dynamical simulations'
  !
  ! The implementation of the PPM was obtained from the Energy Exascale Earth
  ! System Model project, sponsored by the U.S.Department of Energy, Office of
  ! Science, Office of Biological and Environmental Research Earth Systems Model
  ! Development Program area of Earth and Environmental System Modeling.
  !
  ! https://github.com/E3SM-Project/E3SM/blob/master/components/eam/src/dynamics/fv/mapz_module.F90
  !
  ! The subroutines kmppm, steepz, ppm2m and map1_ppm were obtained
  ! from the E3SM code.
  !
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  subroutine kmppm( dm, a4, itot, lmt )
    ! the monotonicity constrained defined in the PPM paper

    use constants_clubb, only: &
        zero, &
        one_fourth, &
        one, &
        two, &
        three, &
        twelve ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    real(core_rknd), intent(in) :: &
        dm(*)
    integer, intent(in) :: &
        itot      ! Total Longitudes
    integer, intent(in) :: &
        lmt     ! 0: Standard PPM constraint
                ! 1: Improved full monotonicity constraint (Lin)
                ! 2: Positive definite constraint
                ! 3: do nothing (return immediately)

    !--------------------- Input/Output Variables ---------------------
    real(core_rknd), intent(inout) :: &
        a4(4,*)   ! AA <-- a4(1,i)
                  ! AL <-- a4(2,i)
                  ! AR <-- a4(3,i)
                  ! A6 <-- a4(4,i)

    !--------------------- Local Variables ---------------------
    real(core_rknd) :: r12
    parameter (r12 = one/twelve)
    real(core_rknd) :: qmp
    integer :: i
    real(core_rknd) :: da1, da2, a6da
    real(core_rknd) :: fmin

    !--------------------- Begin Code ---------------------
    ! Developer: S.-J. Lin, NASA-GSFC

    if ( lmt == 3 ) return

    if( lmt == 0 ) then
        ! Standard PPM constraint
        do i=1,itot
            if(dm(i) == zero) then
                a4(2,i) = a4(1,i)
                a4(3,i) = a4(1,i)
                a4(4,i) = zero
            else
                da1  = a4(3,i) - a4(2,i)
                da2  = da1**2
                a6da = a4(4,i)*da1
                if(a6da < -da2) then
                    a4(4,i) = three*(a4(2,i)-a4(1,i))
                    a4(3,i) = a4(2,i) - a4(4,i)
                elseif(a6da > da2) then
                    a4(4,i) = three*(a4(3,i)-a4(1,i))
                    a4(2,i) = a4(3,i) - a4(4,i)
                endif
            endif
        enddo

    elseif ( lmt == 1 ) then

        ! Improved full monotonicity constraint (Lin)
        ! Note: no need to provide first guess of A6 <-- a4(4,i)
        do i=1, itot
            qmp = two*dm(i)
            a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
            a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
            a4(4,i) = three*( two*a4(1,i) - (a4(2,i)+a4(3,i)) )
        enddo

    elseif (lmt == 2) then

        ! Positive definite constraint
        do i=1,itot
            if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
                fmin = a4(1,i)+one_fourth*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
                if( fmin < zero ) then
                    if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
                        a4(3,i) = a4(1,i)
                        a4(2,i) = a4(1,i)
                        a4(4,i) = zero
                    elseif(a4(3,i) > a4(2,i)) then
                        a4(4,i) = three*(a4(2,i)-a4(1,i))
                        a4(3,i) = a4(2,i) - a4(4,i)
                    else
                        a4(4,i) = three*(a4(3,i)-a4(1,i))
                        a4(2,i) = a4(3,i) - a4(4,i)
                    endif
                endif
            endif
        enddo

    endif

    return
  end subroutine kmppm

  subroutine steepz( i1, i2, km, a4, df2, dm, dq, dp, d4 )
    ! the discontinuity adjustment defined in the PPM paper

    use constants_clubb, only: &
        zero, &
        three_sixteenth, &
        one_half, &
        one ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: km                   ! Total levels
    integer, intent(in) :: i1                   ! Starting longitude
    integer, intent(in) :: i2                   ! Finishing longitude
    real(core_rknd), intent(in) ::  dp(i1:i2,km)       ! grid size
    real(core_rknd), intent(in) ::  dq(i1:i2,km)       ! backward diff of q
    real(core_rknd), intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
    real(core_rknd), intent(in) :: df2(i1:i2,km)       ! first guess mismatch
    real(core_rknd), intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

    !--------------------- Input/Output Variables ---------------------
    real(core_rknd), intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

    !--------------------- Local Variables ---------------------
    integer :: i, k
    real(core_rknd) :: alfa(i1:i2,km)
    real(core_rknd) :: f(i1:i2,km)
    real(core_rknd) :: rat(i1:i2,km)
    real(core_rknd) :: dg2

    !--------------------- Begin Code ---------------------

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
            if(f(i,k+1)*f(i,k-1)<zero .and. df2(i,k)/=zero) then
                dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                      + d4(i,k)*d4(i,k+1) )
                alfa(i,k) = max(zero, min(one_half, -three_sixteenth*dg2/df2(i,k))) 
            else
                alfa(i,k) = zero
            endif
        enddo
    enddo

    do k=4,km-2
        do i=i1,i2
            a4(2,i,k) = (one-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
        enddo
    enddo

    return
  end subroutine steepz

  subroutine ppm2m( a4, delp, km, i1, i2, iv, kord )

    use constants_clubb, only: &
        zero, &
        one_eighth, &
        one_fourth, &
        one_half, &
        three_halves, &
        two, &
        three, &
        four, &
        five, &
        eight ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in):: iv    ! iv =-1: winds
                                ! iv = 0: positive definite scalars
                                ! iv = 1: others
    integer, intent(in):: i1   ! Starting longitude
    integer, intent(in):: i2   ! Finishing longitude
    integer, intent(in):: km   ! vertical dimension
    integer, intent(in):: kord    ! Order (or more accurately method no.):
    real (core_rknd), intent(in):: delp(i1:i2,km)     ! layer pressure thickness

    !--------------------- Input/Output Variables ---------------------
    real (core_rknd), intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

    !--------------------- Local Variables ---------------------
    real(core_rknd) :: dc(i1:i2,km)
    real(core_rknd) :: h2(i1:i2,km)
    real(core_rknd) :: delq(i1:i2,km)
    real(core_rknd) :: df2(i1:i2,km)
    real(core_rknd) :: d4(i1:i2,km)

    real(core_rknd) :: fac
    real(core_rknd) :: a1, a2, c1, c2, c3, d1, d2
    real(core_rknd) :: qmax, qmin, cmax, cmin
    real(core_rknd) :: qm, dq, tmp

    integer :: i, k, km1, lmt
    real(core_rknd) :: qmp, pmp
    real(core_rknd) :: lac
    integer :: it

    !--------------------- Begin Code ---------------------

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
            c1  = (delp(i,k-1)+one_half*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+one_half*delp(i,k))/d4(i,k)
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
            a4(2,i,k) = a4(1,i,k-1) + c1 + two/(d4(i,k-1)+d4(i,k+1)) *    &
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
        dq = two*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
        c1 = four*(a4(2,i,3)-qm-d2*dq) / ( d2*(two*d2*d2+d1*(d2+three*d1)) )
        c3 = dq - one_half*c1*(d2*(five*d1+d2)-three*d1**2)
        a4(2,i,2) = qm - one_fourth*c1*d1*d2*(d2+three*d1)
        a4(2,i,1) = d1*(two*c1*d1**2-c3) + a4(2,i,2)
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
            a4(2,i,1) = max(zero,a4(2,i,1))
            a4(2,i,2) = max(zero,a4(2,i,2))
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
                if( a4(1,i,1)*a4(2,i,1) <=  zero ) then
                    a4(2,i,1) = zero
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
        dq = two*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
        c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(two*d2*d2+d1*(d2+three*d1)))
        c3 = dq - two*c1*(d2*(five*d1+d2)-three*d1**2)
        a4(2,i,km) = qm - c1*d1*d2*(d2+three*d1)
        a4(3,i,km) = d1*(eight*c1*d1**2-c3) + a4(2,i,km)
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
            a4(3,i,km) = max(zero, a4(3,i,km))
        enddo
    elseif ( iv == -1 ) then
        ! Winds:
        do i=i1,i2
            if( a4(1,i,km)*a4(3,i,km) <=  zero ) then
                a4(3,i,km) = zero
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
            a4(4,i,k) = three*(two*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
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
                !           h2(i,k) = two*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
                !    &               / ( delp(i,k)+one_half*(delp(i,k-1)+delp(i,k+1)) )
                !    &               * delp(i,k)**2
                ! Method#3
                h2(i,k) = dc(i,k+1) - dc(i,k-1)
            enddo
        enddo

        if( kord == 7 ) then
            fac = three_halves           ! original quasi-monotone
        else
            fac = one_eighth         ! full monotone
        endif

        do k=3, km-2
            do i=i1,i2
                ! Right edges
                !        qmp   = a4(1,i,k) + two*delq(i,k-1)
                !        lac   = a4(1,i,k) + fac*h2(i,k-1) + one_half*delq(i,k-1)
                !
                pmp   = two*dc(i,k)
                qmp   = a4(1,i,k) + pmp
                lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
                qmin  = min(a4(1,i,k), qmp, lac)
                qmax  = max(a4(1,i,k), qmp, lac)
                a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
                ! Left  edges
                !        qmp   = a4(1,i,k) - two*delq(i,k)
                !        lac   = a4(1,i,k) + fac*h2(i,k+1) - one_half*delq(i,k)
                !
                qmp   = a4(1,i,k) - pmp
                lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
                qmin  = min(a4(1,i,k), qmp, lac)
                qmax  = max(a4(1,i,k), qmp, lac)
                a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
                ! Recompute A6
                a4(4,i,k) = three*(two*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
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
                    a4(4,i,k) = three*(two*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo
            endif
            call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
        enddo
    endif

    do k=km1,km
        do i=i1,i2
            a4(4,i,k) = three*(two*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
        call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
    enddo

    return
  end subroutine ppm2m

  subroutine map1_ppm( km, pe1, q1, kn, pe2, q2, &
                       ng_s, ng_n, itot, i1, i2, &
                       j, jfirst, jlast, iv, kord )

    use constants_clubb, only: &
        zero, &
        one_half, &
        one, &
        two, &
        three ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
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

    !--------------------- Input/Output Variables ---------------------
    real( core_rknd ), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

    ! !DESCRIPTION:
    !
    !     Perform piecewise parabolic method on a given latitude    
    ! IV = 0: constituents
    ! pe1: pressure at layer edges (from model top to bottom surface)
    !      in the original vertical coordinate
    ! pe2: pressure at layer edges (from model top to bottom surface)
    !      in the new vertical coordinate
  
    !--------------------- Local Variables ---------------------
    real(core_rknd) :: r3, r23
    parameter (r3 = one/three, r23 = two/three)
    real(core_rknd) :: dp1(i1:i2,km)
    real(core_rknd) :: q4(4,i1:i2,km)

    integer :: i, k, kk, kl, k0(i1:i2,0:kn+1), k0found
    real(core_rknd) :: pl, pr, qsum, qsumk(i1:i2,kn), delp, esl

    !--------------------- Begin Code ---------------------
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
            if (k0found < 0) then
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
        qsumk(:,k) = zero
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
                q2(i,j,k) = q4(2,i,kk) + one_half*(q4(4,i,kk)+q4(3,i,kk)-q4(2,i,kk))  &
                   *(pr+pl) - q4(4,i,kk)*r3*(pr*(pr+pl)+pl**2)
            else
                ! Consider contribution between pe2(i,k) and pe1(i,kk+1)
                qsum = (pe1(i,kk+1)-pe2(i,k))*(q4(2,i,kk)+one_half*(q4(4,i,kk)+       &
                   q4(3,i,kk)-q4(2,i,kk))*(one+pl)-q4(4,i,kk)*                    &
                   (r3*(one+pl*(one+pl))))
                ! Next consider contribution between pe1(i,kk+1) and pe1(i,k0(i,k+1))
                qsum = qsum + qsumk(i,k)
                ! Now consider contribution between pe1(i,k0(i,k+1)) and pe2(i,k+1)
                kl = k0(i,k+1)
                delp = pe2(i,k+1)-pe1(i,kl)
                esl = delp / dp1(i,kl)
                qsum = qsum + delp*(q4(2,i,kl)+one_half*esl*                          &
                   (q4(3,i,kl)-q4(2,i,kl)+q4(4,i,kl)*(one-r23*esl)))
                q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
            endif
        enddo
    enddo
  end subroutine map1_ppm

  function remap_vals_ppm( ngrdcol, &
                           levels_source_idx, levels_target_idx, &
                           levels_source, levels_target, &
                           source_values, iv )

    ! Description:
    ! This is a wrapper function for E3SM's subroutine implementing the
    ! Piecewise Parabolic Method (PPM) that takes in the pressure levels from surface to top
    !-----------------------------------------------------------------------

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      levels_source_idx, & ! number of levels in the target grid []
      levels_target_idx    ! number of levels in the source grid []

    real( kind = core_rknd ), dimension(ngrdcol,levels_source_idx), intent(in) :: &
      levels_source ! the pressure levels of the source grid [kg/(ms^2)]

    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx), intent(in) :: &
      levels_target ! the pressure levels of the target grid [kg/(ms^2)]

    real( kind = core_rknd ), dimension(ngrdcol,levels_source_idx-1), intent(in) :: &
      source_values  ! given values on the source grid that should be remapped 
                     ! to the target grid

    integer, intent(in) :: &
      iv  ! setting for E3SM's PPM method, that sets what type of value is remapped,
          ! the following are the three valid settings:
          ! iv =-1: winds
          ! iv = 0: positive definite scalars
          ! iv = 1: others

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,levels_target_idx-1) :: &
      remap_vals_ppm ! remapped values with the dimension following the target grid

    !--------------------- Local Variables ---------------------
    integer :: i, k

    real( kind = core_rknd ), dimension(levels_source_idx) :: &
      source_flipped ! levels_source flipped [kg/(ms^2)]

    real( kind = core_rknd ), dimension(levels_target_idx) :: &
      target_flipped ! levels_target flipped [kg/(ms^2)]

    real( kind = core_rknd ), dimension(levels_source_idx-1) :: &
      source_values_flipped  ! source_values flipped

    real( kind = core_rknd ), dimension(levels_target_idx-1) :: &
      target_values_flipped  ! remapped values on the target grid

    ! configurations for the map1_ppm method
    integer :: &
      ghost_lat_south, &  ! Ghosted latitudes south
      ghost_lat_north, &  ! Ghosted latitudes north
      total_lat, &        ! Total latitudes
      start_long, &       ! Starting longitude
      finish_long, &      ! Finishing longitude
      current_lat, &      ! Current latitude
      start_lat, &        ! Starting latitude
      finish_lat, &       ! Finishing latitude
      kord                ! order: should be 4 for vertical remapping

    !--------------------- Begin Code ---------------------

    ghost_lat_south = 0
    ghost_lat_north = 0
    total_lat = 1
    start_long = 1
    finish_long = 1
    current_lat = 1
    start_lat = 1
    finish_lat = 1
    kord = 4

    do i = 1, ngrdcol

      ! the pressure levels and values are flipped, since map1_ppm takes in the grids
      ! from top to surface instead of surface to top, like CLUBB does
      do k = 1, levels_source_idx
          source_flipped(k) = levels_source(i,levels_source_idx - k + 1)
      end do
      do k = 1, levels_target_idx
          target_flipped(k) = levels_target(i,levels_target_idx - k + 1)
      end do
      do k = 1, levels_source_idx-1
          source_values_flipped(k) = source_values(i,levels_source_idx - k)
      end do

      ! remap values using E3SM's PPM implementation
      call map1_ppm( levels_source_idx-1, source_flipped, source_values_flipped, &
                     levels_target_idx-1, target_flipped, target_values_flipped, &
                     ghost_lat_south, ghost_lat_north, total_lat, start_long, finish_long, &
                     current_lat, start_lat, finish_lat, iv, kord )

      ! flip output values back so they correspond to the surface to top ordering in CLUBB's grid
      do k = 1, levels_target_idx-1
          remap_vals_ppm(i,k) = target_values_flipped(levels_target_idx - k)
      end do

    end do

  end function remap_vals_ppm

end module remapping_module