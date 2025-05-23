!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_adaptation_module


  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: setup_gr_dycore, &
            normalize_grid_density, &
            adapt_grid, &
            calc_grid_dens, &
            clean_up_grid_adaptation_module

  private :: setup_gr_min, calc_integral, create_fixed_min_gr_dens_func, &
             normalize_min_grid_density, normalize_grid_density_helper, &
             decide_if_grid_adapt, apply_equidistribution_principle, &
             create_grid_from_normalized_grid_density, check_grid, &
             calc_grid_dens_helper, moving_average, remap_all_clubb_core_vals

  private ! Default Scoping

  real( kind = core_rknd ), parameter ::  &
    tol = 1.0e-6_core_rknd ! tolerance to check if to real numbers are equal

  ! Save prescribed minimum grid density profile, this is then interpolated every timestep
  ! to the current physics grid
  integer :: &
    fixed_min_gr_dens_idx ! number of elements in the vertical axis of fixed_min_gr_dens

  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    fixed_min_gr_dens_z, &  ! z coordinates of the piecewise linear minimum
                            ! grid density function                                       [m]
    fixed_min_gr_dens       ! minimum grid density values given at the z coordinates      [1/m]

  ! Save the normalized grid density from the last grid adaptation, to check if new grid would
  ! be significantly different (grid adaptation trigger)
  integer :: &
    gr_dens_old_idx_global ! number of levels in gr_dens_old_global and gr_dens_old_z_global
  
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
    gr_dens_old_z_global, &  ! grid density altitudes from adaptation step before [m]
    gr_dens_old_global       ! grid density from adaptation step before           [1/m]

  ! Save the calculated values for the refinement criterion terms as a cumulative sum starting
  ! from the last grid adaptation, to take the time average over that region, to have a
  ! smoother evolution in time
  integer :: &
    Lscale_counter, &          ! the number of values cumulated in cumulative_Lscale
    richardson_num_counter, &  ! the number of values cumulated in cumulative_richardson_num
    chi_counter                ! the number of values cumulated in cumulative_chi

  real( kind = core_rknd ), dimension(:), allocatable :: &
    cumulative_Lscale, &          ! cumulative sum of the Lscale term in the
                                  ! refinement cirterion                                  [1/m]
    cumulative_richardson_num, &  ! cumulative sum of the richardson number term in the
                                  ! refinement cirterion                                  [1/m]
    cumulative_chi                ! cumulative sum of the chi term in the
                                  ! refinement cirterion                                  [1/m]

  contains

  !=============================================================================

  ! TODO find better way to set these up, since both do the same thing, just with different zm
  ! file paths, and maybe some solution where parameters are not hardcoded, so maybe something
  ! where the parameters are read from model.in file and dont use hardcoded 73 in dimension of
  ! heights vectors
  subroutine setup_gr_dycore( iunit, ngrdcol, grid_sfc, grid_top, gr )

    ! Description:
    ! Thid subroutine is used to set up the dycore grid for the grid adaptation

    use grid_class, only: &
        read_grid_heights, &
        setup_grid, & ! Procedures
        grid          ! Type

    use error_code, only: &
        clubb_at_least_debug_level_api, &
        clubb_fatal_error

    use constants_clubb, only: &
        fstderr ! Constants
    
    use err_info_type_module, only: &
        err_info_type        ! Type

    implicit none

    !--------------------- Input Variables ---------------------
    integer :: &
      iunit, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
        grid_sfc, &    ! grids surface height for the dycore grid [m]
        grid_top       ! grids highest level for the dycore grid  [m]

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- In/Output Variables ---------------------

    !--------------------- Local Variables ---------------------

    integer :: &
      i, &
      nlevel, &
      level_min

    real( kind = core_rknd ), dimension(ngrdcol,73) :: & ! since dycore grid has 73 levels in total
      momentum_heights      ! Momentum level altitudes (file input)      [m]
    
    real( kind = core_rknd ), dimension(ngrdcol,73-1) :: &! since dycore grid has 73 levels in total
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      deltaz ! vertical grid spacing [m]

    ! Configuration for read_grid_heights and setup_grid
    integer, parameter :: &
      grid_type = 3, &
      nzmax = 73

    logical, parameter :: &
      l_implemented = .false.

    character(len=100), parameter :: & 
      zm_grid_fname  = '../input/grid/dycore.grd',       & ! Path and filename of file for
                                                           ! momentum level altitudes
      zt_grid_fname = ''                                   ! Path and filename of file for
                                                           ! thermodynamic level altitudes

    type(err_info_type) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Begin Code ---------------------

    do i = 1, ngrdcol
      call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                              zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                              iunit, &                            ! Intent(in)
                              momentum_heights(i,:), &            ! Intent(out)
                              thermodynamic_heights(i,:), &       ! Intent(out)
                              err_info )                          ! Intent(inout)
      sfc_elevation(i) = 0.0_core_rknd
      deltaz(i) = 0.0_core_rknd
    end do

    if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling read_grid_heights in setup_gr_dycore"
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
    
    call setup_grid( nlevel, ngrdcol, sfc_elevation, l_implemented, &           ! intent(in)
                     .true., grid_type, deltaz, grid_sfc, grid_top, &           ! intent(in)
                     momentum_heights, thermodynamic_heights, &                 ! intent(in)
                     gr, err_info )                                             ! intent(inout)

    if ( clubb_at_least_debug_level_api(0) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid in setup_gr_dycore"
      end if
    end if

  end subroutine setup_gr_dycore

  subroutine setup_gr_min( iunit, ngrdcol, grid_sfc, grid_top, gr )

    ! Description:
    ! Thid subroutine is used to get the minimum grid density profile from a file

    use grid_class, only: &
        read_grid_heights, &
        setup_grid, & ! Procedures
        grid          ! Type

    use error_code, only: &
        clubb_fatal_error

    use constants_clubb, only: &
        fstderr ! Constants

    use err_info_type_module, only: &
        err_info_type        ! Type
    
    implicit none

    !--------------------- Input Variables ---------------------
    integer :: &
      iunit, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
        grid_sfc, &    ! grids surface height for the dycore grid [m]
        grid_top       ! grids highest level for the dycore grid  [m]

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- Local Variables ---------------------

    integer :: &
      i, &
      nlevel, &
      level_min

    real( kind = core_rknd ), dimension(ngrdcol,73) :: & ! since dycore grid has 73 levels in total
      momentum_heights      ! Momentum level altitudes (file input)      [m]
    
    real( kind = core_rknd ), dimension(ngrdcol,73-1) :: &! since dycore grid has 73 levels in total
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      deltaz ! vertical grid spacing [m]

    ! Configuration for read_grid_heights and setup_grid
    integer, parameter :: &
      grid_type = 3, &
      nzmax = 73

    logical, parameter :: &
      l_implemented = .false.

    character(len=100), parameter :: & 
      zm_grid_fname  = '../input/grid/gr_min.grd',  & ! Path and filename of file for
                                                      ! momentum level altitudes
      zt_grid_fname = ''                              ! Path and filename of file for
                                                      ! thermodynamic level altitudes

    type(err_info_type) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Begin Code ---------------------

    do i = 1, ngrdcol
      call read_grid_heights( nzmax, grid_type, &                 ! Intent(in)
                              zm_grid_fname, zt_grid_fname, &     ! Intent(in)
                              iunit, &                            ! Intent(in)
                              momentum_heights(i,:), &            ! Intent(out)
                              thermodynamic_heights(i,:), &       ! Intent(out)
                              err_info )                          ! Intent(inout)
      sfc_elevation(i) = 0.0_core_rknd
      deltaz(i) = 0.0_core_rknd
    end do

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "Fatal error calling read_grid_heights in setup_gr_min"
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
    
    call setup_grid( nlevel, ngrdcol, sfc_elevation, l_implemented, &           ! intent(in)
                     .true., grid_type, deltaz, grid_sfc, grid_top, &           ! intent(in)
                     momentum_heights, thermodynamic_heights, &                 ! intent(in)
                     gr, err_info )                                             ! intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "Fatal error calling setup_grid in setup_gr_min"
    end if

  end subroutine setup_gr_min

  function calc_integral( g_idx, g_x, g_y )
    ! Description:
    ! Computes the exact integral, assuming the given points build a piecewise linear function.
    ! g_x and g_y must both have g_idx elements

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
      zero

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      g_idx  ! The total numer of indices of g_x and g_y []

    real( kind = core_rknd ), dimension(g_idx), intent(in) ::  &
      g_x,  &    ! the x coordinates of the connection points of the piecewise linear function [m]
      g_y        ! the y coordinates of the connection points of the piecewise linear function

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ) :: &
      calc_integral ! The integral

    !--------------------- Local Variables ---------------------
    integer :: i

    !--------------------- Begin Code ---------------------
    ! Initializing calc_integral to avoid a compiler warning.
    calc_integral = zero

    ! Calculate each linear segments exact integral and add them
    do i = 1, (g_idx-1)
        calc_integral = calc_integral + (g_y(i) + g_y(i+1))/2 * (g_x(i+1) - g_x(i))
    end do

  end function calc_integral

  subroutine create_fixed_min_gr_dens_func( iunit, ngrdcol, &
                                            grid_sfc, grid_top )

    ! Description:
    ! Creates and allocates the fixed minimum grid density function from read in minimum grid.
    ! Takes the inverse grid distance between zm levels.

    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use error_code, only: &
        clubb_at_least_debug_level_api

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      iunit, &
      ngrdcol

    real( kind = core_rknd ), intent(in) :: &
      grid_sfc, & ! altitude of the grids surface [m]
      grid_top    ! altitude of the grids top     [m]

    !--------------------- Local Variables ---------------------
    integer :: i, j

    type( grid ) :: gr_min ! The grid that is read in from a file and that is used for the
                           ! prescribed non-normalized minimum grid density

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      grid_sfc_arr, & ! an array with grid_sfc
      grid_top_arr    ! an array with grid_top

    !--------------------- Begin Code ---------------------
    if ( (.not. allocated( fixed_min_gr_dens_z )) .or. (.not. allocated( fixed_min_gr_dens )) ) then
      do i = 1, ngrdcol
        grid_sfc_arr(i) = grid_sfc
        grid_top_arr(i) = grid_top
      end do

      ! Initialize the minimum grid to use as a profile for the prescribed minimum grid density
      call setup_gr_min( iunit, ngrdcol, &              ! Intent(in)
                         grid_sfc_arr, grid_top_arr, &  ! Intent(in)
                         gr_min )                       ! Intent(out)
      fixed_min_gr_dens_idx = gr_min%nzt+2
    end if

    if ( .not. allocated( fixed_min_gr_dens_z ) ) then
      allocate( fixed_min_gr_dens_z(ngrdcol,gr_min%nzt+2) )
      do i = 1, ngrdcol
        fixed_min_gr_dens_z(i,1) = grid_sfc
        do j = 1, gr_min%nzt
          fixed_min_gr_dens_z(i,j+1) = gr_min%zt(i,j)
        end do
        fixed_min_gr_dens_z(i,gr_min%nzt+2) = grid_top
      end do
    end if

    ! TODO maybe just use invrs_dzm instead since it has already nzm levels,
    ! but this will probably change results slightly, but then fixed_min_gr_dens_idx and
    ! fixed_min_gr_dens_z would also change
    if ( .not. allocated( fixed_min_gr_dens ) ) then
      allocate( fixed_min_gr_dens(ngrdcol,gr_min%nzt+2) )
      do i = 1, ngrdcol
        ! Use linear extrapolation to calculate points on the outer zm levels
        fixed_min_gr_dens(i,1) = gr_min%invrs_dzt(i,1) &
                                 + (fixed_min_gr_dens_z(i,1) - gr_min%zt(i,1)) &
                                   /(gr_min%zt(i,2) - gr_min%zt(i,1)) &
                                    *(gr_min%invrs_dzt(i,2) - gr_min%invrs_dzt(i,1))

        ! Check if extrapolated value is non-negative
        if ( clubb_at_least_debug_level_api( 2 ) ) then
          if ( fixed_min_gr_dens(i,1) <= 0 ) then
            error stop 'Initial minimum grid density function needs to be positive.'
          end if
        end if
        
        do j = 1, gr_min%nzt
          fixed_min_gr_dens(i,j+1) = gr_min%invrs_dzt(i,j)

          ! Check if value is non-negative
          if ( clubb_at_least_debug_level_api( 2 ) ) then
            if ( fixed_min_gr_dens(i,j+1) <= 0 ) then
              error stop 'Initial minimum grid density function needs to be positive.'
            end if
          end if
        end do

        ! Use linear extrapolation to calculate points on the outer zm levels
        fixed_min_gr_dens(i,gr_min%nzt+2) = gr_min%invrs_dzt(i,gr_min%nzt-1) &
                                               + (fixed_min_gr_dens_z(i,gr_min%nzt+2) &
                                                  - gr_min%zt(i,gr_min%nzt-1)) &
                                                 /(gr_min%zt(i,gr_min%nzt) &
                                                   - gr_min%zt(i,gr_min%nzt-1)) &
                                                  *(gr_min%invrs_dzt(i,gr_min%nzt) &
                                                    - gr_min%invrs_dzt(i,gr_min%nzt-1))

        ! Check if extrapolated value is non-negative
        if ( clubb_at_least_debug_level_api( 2 ) ) then
          if ( fixed_min_gr_dens(i,gr_min%nzt+2) <= 0 ) then
            error stop 'Initial minimum grid density function needs to be positive.'
          end if
        end if

      end do
    end if
    
  end subroutine create_fixed_min_gr_dens_func

  function normalize_min_grid_density( ngrdcol, &
                                       min_gr_dens_idx, &
                                       min_gr_dens_z, min_gr_dens, &
                                       lambda, num_levels ) result( min_gr_dens_norm )

    ! Description:
    ! Takes the prescribed non-normalized minimum grid density and normalizes it.
    ! Normalized means in this case, that the integral over
    ! the grid region is between 0 (excluded) and num_levels-1 (included). What value
    ! should be taken inbetween is controlled with lambda.

    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero, one, fstderr ! Constants

    use error_code, only: &
        clubb_at_least_debug_level_api

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      min_gr_dens_idx  ! number of levels in min_gr_dens_z

    ! These two arrays define the points of the piecewise linear prescribed non-normalized
    ! minimum grid density
    real( kind = core_rknd ), dimension(ngrdcol,min_gr_dens_idx), intent(in) :: &
      min_gr_dens_z, &  ! z coordinates of the piecewise linear minimum grid density [m]
      min_gr_dens       ! density values given at the z coordintas                   [1/m]

    real( kind = core_rknd ), intent(in) :: &
      lambda      ! number between 0 (excluded) and 1 (included) used as an adjustable factor;
                  ! lambda controls the integral of the normalized minimum grid density;
                  ! after the normalization the integral of the minimum grid density
                  ! will be lambda*(num_levels-1)

    integer, intent(in) :: &
      num_levels  ! number of wanted grid levels for the new grid

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,min_gr_dens_idx) :: &
      min_gr_dens_norm      ! the normalized minimal grid density values given at
                            ! the z coordinates of min_gr_dens_z                     [1/m]

    !--------------------- Local Variables ---------------------
    integer :: i, j

    real( kind = core_rknd ) :: &
      integral, &
      norm_factor

    !--------------------- Begin Code ---------------------
    ! check if lambda is in (0,1]
    if ( lambda <= zero .or. lambda > one ) then
      error stop 'lambda needs to be between 0 (excluded) and 1 (included)'
    end if

    do i = 1, ngrdcol
      ! calculate factor to normalize minimum grid density function
      integral = calc_integral( min_gr_dens_idx, min_gr_dens_z(i,:), min_gr_dens(i,:) )
      norm_factor = lambda*(num_levels - 1)/integral

      ! normalize minimum grid density function with norm_factor and copy z coordinates
      do j = 1, (min_gr_dens_idx)
        min_gr_dens_norm(i,j) = norm_factor*min_gr_dens(i,j)
      end do
    end do

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      ! check if integral is lambda*(n-1) as wanted
      integral = calc_integral( min_gr_dens_idx, min_gr_dens_z, min_gr_dens_norm )
      if ( abs(integral - (lambda*(num_levels-1))) > tol ) then
        write(fstderr,*) 'Warning! Integral in normalize_min_grid_density', &
                         ' should be', lambda*(num_levels-1), ' but is ', integral
        error stop 'Something went wrong, integral is different than it should be.'
      end if

      ! check if normalized function is always >0
      do i = 1, ngrdcol
        do j = 1, min_gr_dens_idx
          if ( min_gr_dens_norm(i,j) <= zero ) then
            error stop 'Minimum grid density function in needs to be positive.'
          end if
        end do
      end do

    end if
    
  end function normalize_min_grid_density

  function normalize_grid_density_helper( ngrdcol, &
                                            gr_dens_idx, &
                                            gr_dens_z, gr_dens, &
                                            min_gr_dens_norm, &
                                            lambda, num_levels ) result( gr_dens_norm )

    ! Description:
    ! Takes the linear piecewise grid density function (gr_dens_z, gr_dens) each of size
    ! gr_dens_idx and normalizes this density, such that the integral is the
    ! number of desired grid levels (num_levels)-1 and the minimum is above the minimal
    ! density function (gr_dens_z,min_gr_dens_norm)

    ! Note: The minimum grid density function needs to be normalized with
    !       normalize_min_grid_density before using it here.

    ! Note: The initial grid density function needs to be non-negative.

    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constants

    use error_code, only: &
        clubb_at_least_debug_level_api

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      gr_dens_idx, &           ! number of elements in gr_dens_z, gr_dens, min_gr_dens_norm_z
                               ! and min_gr_dens_norm
      num_levels               ! desired number of grid levels

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) :: &
      gr_dens_z, & ! the prescribed non-normalized grid density z coordinates              [m]
      gr_dens      ! the prescribed non-normalized grid density given at the z coordinates [1/meter]
      ! gr_dens_z and gr_dens should be ordered from the bottom to the top

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) :: &
      min_gr_dens_norm  ! the minimum grid density given at the z coordinates of gr_dens_z [1/m]
      ! min_gr_dens_norm should be ordered from the bottom to the top
      ! Note: min_gr_dens_norm should be given on the gr_dens_z altitudes

    real( kind = core_rknd ), intent(in) :: &
      lambda   ! a factor to adjust the integral of the minimum grid density

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) :: &
      gr_dens_norm      ! the normalized grid density values given at the z coordinates
                        ! of gr_dens_z [1/m]

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

      ! calculate integrals and check
      integral_gr_dens = calc_integral( gr_dens_idx, gr_dens_z(i,:), gr_dens(i,:) )
      integral_min_gr_dens = calc_integral( gr_dens_idx, &
                                            gr_dens_z(i,:), min_gr_dens_norm(i,:) )

      if ( abs(integral_min_gr_dens - (lambda*(num_levels-1))) > tol ) then
        write(fstderr,*) 'Warning! The minimum grid density function has not the correct integral.'
        error stop 'Normalize the minimum grid density function before using it.'
      end if

      if ( integral_gr_dens < 1e-8_core_rknd ) then
        ! this happens if the prescribed non-normalized grid density was zero at all altitudes

        ! calculate shift value
        grid_sfc = gr_dens_z(1,1)
        grid_top = gr_dens_z(1,gr_dens_idx)
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
          gr_dens_norm(i,j) = norm_factor*gr_dens(i,j) + min_gr_dens_norm(i,j)
        end do

      end if

      ! run checks
      if ( clubb_at_least_debug_level_api( 2 ) ) then

        ! check if integral is (n-1) as wanted
        integral_test = calc_integral( gr_dens_idx, gr_dens_z(i,:), gr_dens_norm(i,:) )
        if ( abs(integral_test - (num_levels-1)) > tol ) then
          write(fstderr,*) 'Warning! Integral in normalize_grid_density_helper', &
                           ' should be something different.'
          error stop 'Something went wrong, integral is different than it should be.'
        end if

        ! check if normalized function is always >=min_gr_dens_norm
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

  end function normalize_grid_density_helper

  subroutine normalize_grid_density( ngrdcol, iunit, &
                                     gr_dens_idx, &
                                     gr_dens_z, gr_dens, &
                                     lambda, &
                                     num_levels, &
                                     norm_min_grid_dens, &
                                     norm_grid_dens )
    ! Description:
    ! Normalizes the prescribed grid density.

    !-----------------------------------------------------------------------

    use interpolation, only: &
        zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      iunit, &
      gr_dens_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have       []

    real( kind = core_rknd ) :: &
      lambda   ! a factor to adjust the integral of the minimum grid density

    ! These two arrays build the piecewise linear prescribed non-normalized grid density
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) ::  &
      gr_dens_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function                                                  [m]
      gr_dens         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function        [1/m]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(out) :: &
      norm_min_grid_dens, & ! the density at the given z coordinates of the connection points
                            ! of the normalized piecewise linear grid density function         [1/m]
      norm_grid_dens        ! the density at the given z coordinates of the connection points
                            ! of the normalized piecewise linear grid density function         [1/m]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens        ! the density at the given z coordinates of the connection points
                         ! of the normalized piecewise linear grid density function          [1/m]

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface   [m]
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

    norm_min_grid_dens = normalize_min_grid_density( ngrdcol, &
                                                     gr_dens_idx, &
                                                     gr_dens_z, min_gr_dens, &
                                                     lambda, num_levels )

    norm_grid_dens = normalize_grid_density_helper( ngrdcol, &
                                                    gr_dens_idx, &
                                                    gr_dens_z, gr_dens, &
                                                    norm_min_grid_dens, &
                                                    lambda, num_levels )

  end subroutine normalize_grid_density
  
  function decide_if_grid_adapt( ngrdcol, &
                                 gr_dens_old_idx, &
                                 gr_dens_old_z, &
                                 gr_dens_old, &
                                 gr_dens_new_idx, &
                                 gr_dens_new_z, &
                                 gr_dens_new, &
                                 threshold ) result( l_adapt_grid )

    ! Description:
    ! Checks how different the two noramlized grid densities are and returns depending
    ! on the threshold if the grid should be adapted or not.

    !-----------------------------------------------------------------------

    use interpolation, only: &
        zlinterp_fnc

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      gr_dens_old_idx, &  ! The number of indices of gr_dens_old_z
                          ! and gr_dens_old []
      gr_dens_new_idx     ! The number of indices of gr_dens_new_z
                          ! and gr_dens_new []

    ! These two arrays build the piecewise linear normalized grid density from
    ! the last grid adaptation
    real( kind = core_rknd ), dimension(ngrdcol, gr_dens_old_idx), intent(in) ::  &
      gr_dens_old_z,  &  ! the z coordinates of the connection points of the
                         ! piecewise linear grid density function from the adaptation
                         ! step before                                                       [m]
      gr_dens_old        ! the density values at the given z coordinates of the connection
                         ! points of the piecewise linear grid density
                         ! function from the adaptation step before                          [1/m]

    ! These two arrays build the piecewise linear normalized grid density from
    ! the current iteration
    real( kind = core_rknd ), dimension(ngrdcol, gr_dens_new_idx), intent(in) ::  &
      gr_dens_new_z,  &  ! the z coordinates of the connection points of the
                         ! piecewise linear grid density function from the current
                         ! adaptation step                                                  [m]
      gr_dens_new        ! the density values at the given z coordinates of the
                         ! connection points of the piecewise linear grid density
                         ! function from the current adaptation step                        [1/m]

    real( kind = core_rknd ), intent(in) ::  &
      threshold   ! threshold which decides if grid should be adapted or not []

    !--------------------- Output Variable ---------------------
    logical, dimension(ngrdcol) :: &
      l_adapt_grid ! true or false, whether grid should be adapted or not

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(gr_dens_old_idx + gr_dens_new_idx) ::  &
      common_grid ! the grid tha results from combining the grid the old grid density was
                  ! given on and the grid the current grid density is given on              [m]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      gr_dens_old_interp ! old normalized grid density interpolated to the common grid  [1/m]

    real( kind = core_rknd ), dimension(:), allocatable ::  &
      gr_dens_new_interp ! new normalized grid density interpolated to the common grid  [1/m]

    integer :: i, j, k, l

    logical :: l_adapt_grid_tmp

    real( kind = core_rknd ) :: &
      sum

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

      ! calculate normalized mean squared error sum
      sum = 0
      do k = 1, j
        sum = sum + ( gr_dens_old_interp(k) - gr_dens_new_interp(k) )**2/gr_dens_old_interp(k)**2
      end do
      
      ! check if mean value is bigger than threshold, if so, adapt
      if ( sum/j > threshold ) then
        l_adapt_grid_tmp = .true.
      end if

      deallocate( gr_dens_old_interp )
      deallocate( gr_dens_new_interp )

      l_adapt_grid(i) = l_adapt_grid_tmp
    end do

  end function decide_if_grid_adapt

  function apply_equidistribution_principle( num_levels, &
                                             gr_dens_norm_idx, &
                                             gr_dens_norm_z, &
                                             gr_dens_norm ) result( grid_heights )
    ! Description:
    ! Creates the grid from the normalized grid density function following the equidistribution
    ! principle.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      gr_dens_norm_idx, &  ! The total numer of indices of gr_dens_norm_z and gr_dens_norm []
      num_levels           ! number of levels the new grid should have []

    real( kind = core_rknd ), dimension(gr_dens_norm_idx), intent(in) ::  &
      gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized
                          ! piecewise linear grid density function                          [m]
      gr_dens_norm        ! the density values at the given z coordinates of the
                          ! density function connection points of the normalized
                          ! piecewise linear grid                                           [1/m]
      ! The grid density function needs to be normalized with the number of grid levels num_levels,
      ! before the function can be applied here

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    !--------------------- Local Variables ---------------------
    integer :: &
      i, &
      prev_x_ind  ! the index of the last used x coordinate []

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid              [m]
      grid_top, &               ! the highest point of the grid               [m]
      area_up_to_prev_x, &      ! the integral of the function up
                                ! to the last used x coordinate               [# levs]
      new_x_area_increment, &   ! integral from grid density
                                ! function on [prev_x_ind, prev_x_ind+1]      [# levs]
      desired_area_up_to_ith_level, & ! the desired area up to the ith grid level, should be exactly
                                      ! i, since the grid density function is normalized so that the
                                      ! integral of the whole function is num_levels-1 and since we
                                      ! follow the equi-distribution approach [# levs]
      A, grid_level, slope, intercept, p_pq_formula, q_pq_formula

    logical :: l_grid_level_found

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_norm_z(1)
    grid_top = gr_dens_norm_z(gr_dens_norm_idx)

    grid_heights(1) = grid_sfc

    ! Initializing calc_integral to avoid a compiler warning.
    prev_x_ind = 1
    area_up_to_prev_x = 0.0_core_rknd

    do i = 2, num_levels-1
        l_grid_level_found = .false.
        do while (.not. l_grid_level_found)
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
                l_grid_level_found = .true.
            else
                area_up_to_prev_x = area_up_to_prev_x + new_x_area_increment
                prev_x_ind = prev_x_ind + 1
            end if
        end do
    end do

    grid_heights(num_levels) = grid_top

  end function apply_equidistribution_principle

  subroutine create_grid_from_normalized_grid_density( ngrdcol, &
                                                       gr_dens_norm_idx, &
                                                       gr_dens_norm_z, gr_dens_norm, &
                                                       min_gr_dens_norm_z, min_gr_dens_norm, &
                                                       num_levels, &
                                                       threshold, &
                                                       grid_heights, &
                                                       l_adapt_grid )
    ! Description:
    ! Creates the grid from the non-normalized minimum grid density and the non-normalized
    ! grid density following the equidistribution principle.

    !-----------------------------------------------------------------------

    use error_code, only: &
        clubb_at_least_debug_level_api

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      gr_dens_norm_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have []

    real( kind = core_rknd ), intent(in) :: &
      threshold ! threshold which decides if grid should be adapted or not

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(in) ::  &
      ! The two arrays min_gr_dens_norm_z and min_gr_dens_norm build the piecewise linear
      ! normalized minimum grid density
      min_gr_dens_norm_z,  & ! altitudes [m]
      min_gr_dens_norm,  &   ! densities [1/m]
      ! The two arrays gr_dens_norm_z and gr_dens_norm build the piecewise linear
      ! normalized grid density
      gr_dens_norm_z,  &     ! altitudes [m]
      gr_dens_norm           ! desnities [1/m]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_levels), intent(out) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    logical, dimension(ngrdcol), intent(out) :: &
      l_adapt_grid ! whether the grid should be adapted or not

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx, i

    real( kind = core_rknd ) ::  &
      grid_sfc,  &    ! height of the grids surface   [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_norm_z(1,1)
    grid_top = gr_dens_norm_z(1,gr_dens_norm_idx)

    ! decide if grid should be adapted or not
    if ( .not. allocated( gr_dens_old_z_global ) .and. .not. allocated( gr_dens_old_global ) ) then
      error stop 'gr_dens_old_z_global and gr_dens_old_global were not allocated...'
    end if
    l_adapt_grid = decide_if_grid_adapt( ngrdcol, &
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
        grid_heights(i,:) = apply_equidistribution_principle( num_levels, &
                                                                           gr_dens_norm_idx, &
                                                                           gr_dens_norm_z(i,:), &
                                                                           gr_dens_norm(i,:) )
      end if
    end do

    if ( clubb_at_least_debug_level_api( 2 ) ) then

        grid_heights_idx = size(grid_heights, 2)

        do i = 1, ngrdcol
          if ( l_adapt_grid(i) ) then
            ! only check if grid was actually adapted
            call check_grid( grid_heights_idx, grid_heights(i,:), &            ! Intent(in)
                             num_levels, &                                     ! Intent(in)
                             gr_dens_norm_idx, &                               ! Intent(in)
                             min_gr_dens_norm_z(i,:), min_gr_dens_norm(i,:), & ! Intent(in)
                             grid_sfc, grid_top )                              ! Intent(in)
          end if
        end do

    end if

  end subroutine create_grid_from_normalized_grid_density

  subroutine check_grid( grid_heights_idx, grid_heights, &
                         desired_num_levels, &
                         desired_min_dens_idx, &
                         desired_min_dens_z, desired_min_dens, &
                         desired_grid_sfc, desired_grid_top )

    ! Description:
    ! Checks if the grid is valid.

    !-----------------------------------------------------------------------

    use interpolation, only: &
        zlinterp_fnc ! Procedure

    use constants_clubb, only: &
        one, fstderr ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      grid_heights_idx, &      ! number of elements in grid_heights                            []
      desired_min_dens_idx, &  ! number of elements in desired_min_dens_z and desired_min_dens []
      desired_num_levels       ! desired number of grid levels                                 []
    ! Note: if everything worked fine, those two integers should be the same

    real( kind = core_rknd ), dimension(grid_heights_idx), intent(in) :: &
      grid_heights ! the grid heights ordered from bottom to top [m]

    real( kind = core_rknd ), dimension(desired_min_dens_idx), intent(in) :: &
      desired_min_dens_z, & ! the z coordinates of the piecewise linear desired
                            ! minimum grid density function                         [m]
      desired_min_dens      ! the density values on the given z coordinates         [1/m]

    real( kind = core_rknd ), intent(in) :: &
      desired_grid_sfc, &   ! desired surface height of the grid [m]
      desired_grid_top      ! desired top height of the grid     [m]

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

  ! TODO add support for ngrdcol
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
                             richardson_num_term, &
                             richardson_num_term_time_avg )
    ! Description:
    ! Calculate the necessary variables and then construct the non-normalized
    ! grid density from them.

    !-----------------------------------------------------------------------

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use grid_class, only: &
        grid, &  ! Type
        zt2zm_api, & ! Procedures
        ddzt

    use model_flags, only: &
        clubb_config_flags_type  ! Type

    use advance_helper_module, only: &
        calc_brunt_vaisala_freq_sqd ! Procedure

    use constants_clubb, only: &
        one

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol
        
    type( grid ), intent(in) :: &
      gr ! grid object on which the values are given
    
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      Lscale_zt,         & ! Length scale given on the zt levels                [m]
      um,                & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,                & ! northward grid-mean wind component (thermo. levs.) [m/s]
      thvm,              & ! Virtual potential temperature                      [K]
      exner,             & ! Exner function (thermodynamic levels)              [-]
      rtm,               & ! total water mixing ratio, r_t (thermo. levels)     [kg/kg]
      thlm,              & ! liq. water pot. temp., th_l (thermo. levels)       [K]
      rcm,               & ! cloud water mixing ratio, r_c (thermo. levels)     [kg/kg]
      p_in_Pa,           & ! pressure on the thermodynamic levels               [Pa]
      ice_supersat_frac    ! ice cloud fraction                                 [-]

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
    ! These two arrays build the piecewise linear grid density
    real( kind = core_rknd ), dimension(gr%nzm), intent(out) ::  &
      gr_dens_z,  &    ! altitudes [m]
      gr_dens          ! densities [1/m]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      alt_term, &                  ! altitude term of the refinement criterion               [1/m]
      Lscale_term, &               ! Lscale term of the refinement criterion                 [1/m]
      chi_term, &                  ! chi term of the refinement criterion                    [1/m]
      richardson_num_term, &       ! richardson number term of the refinement criterion      [1/m]
      Lscale_term_time_avg, &      ! time averaged Lscale term of the refinement criterion   [1/m]
      chi_term_time_avg, &         ! time averaged chi term of the refinement criterion      [1/m]
      richardson_num_term_time_avg ! time averaged richardson number term of the refinement
                                   ! criterion                                               [1/m]

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

    ! Interpolate Lscale to momentum levels (zm)
    Lscale = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, Lscale_zt )

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
    chi = zt2zm_api( gr, chi_zt )

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
    call calc_grid_dens_helper( ngrdcol, gr%nzm, gr%zm,        & ! In
                                Lscale, wp2,                   & ! In
                                chi,                           & ! In
                                brunt_vaisala_freq_sqd,        & ! In
                                ddzt_umvm_sqd,                 & ! In
                                gr_dens_z, gr_dens,            & ! Out
                                alt_term, Lscale_term,         & ! Out
                                Lscale_term_time_avg,          & ! Out
                                chi_term,                      & ! Out
                                chi_term_time_avg,             & ! Out
                                richardson_num_term,           & ! Out
                                richardson_num_term_time_avg )   ! Out

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
                                    richardson_num_term, &
                                    richardson_num_term_time_avg )
    ! Description:
    ! Calculates the non-normalized grid density from Lscale, chi,
    ! ddzt_umvm_sqd and brunt_vaisala_freq_sqd.

    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero

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
    ! These two arrays build the piecewise linear non-normalized grid density
    real( kind = core_rknd ), dimension(nzm), intent(out) ::  &
      gr_dens_z,  &    ! altitudes [m]
      gr_dens          ! densities [1/m]

    real( kind = core_rknd ), dimension(nzm), intent(out) :: &
      alt_term, &                  ! altitude term of the refinement criterion               [1/m]
      Lscale_term, &               ! Lscale term of the refinement criterion                 [1/m]
      chi_term, &                  ! chi term of the refinement criterion                    [1/m]
      richardson_num_term, &       ! richardson number term of the refinement criterion      [1/m]
      Lscale_term_time_avg, &      ! time averaged Lscale term of the refinement criterion   [1/m]
      chi_term_time_avg, &         ! time averaged chi term of the refinement criterion      [1/m]
      richardson_num_term_time_avg ! time averaged richardson number term of the refinement
                                   ! criterion                                               [1/m]

    !--------------------- Local Variables ---------------------
    integer :: k

    real( kind = core_rknd ), parameter :: &
      wp2_threshold = 0.001 ! threshold for wp2 that decides whether Lscale is used or not [m^2/s^2]

    !--------------------- Begin Code ---------------------
    ! Prepare cumulative_Lscale to use for time averaged Lscale term
    if ( .not. allocated( cumulative_Lscale ) ) then
      allocate( cumulative_Lscale(nzm) )
      Lscale_counter = 0
      do k = 1, nzm
        cumulative_Lscale(k) = zero
      end do
    end if
    Lscale_counter = Lscale_counter + 1

    ! Prepare cumulative_richardson_num to use for time averaged richardson number term
    if ( .not. allocated( cumulative_richardson_num ) ) then
      allocate( cumulative_richardson_num(nzm) )
      richardson_num_counter = 0
      do k = 1, nzm
        cumulative_richardson_num(k) = zero
      end do
    end if
    richardson_num_counter = richardson_num_counter + 1

    ! Prepare cumulative_chi to use for time averaged chi term
    if ( .not. allocated( cumulative_chi ) ) then
      allocate( cumulative_chi(nzm) )
      chi_counter = 0
      do k = 1, nzm
        cumulative_chi(k) = zero
      end do
    end if
    chi_counter = chi_counter + 1

    ! Calculate each refinement criterion term
    do k = 1, nzm
      gr_dens_z(k) = zm(1,k)

      ! Calculate alt_term
      alt_term(k) = 1/(gr_dens_z(k) + 100)

      ! Calculate Lscale_term and richardson number term
      if ( wp2(1,k) > wp2_threshold ) then
        Lscale_term(k) = 1/(Lscale(1,k)+10.0)
        richardson_num_term(k) = maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(1,k))]) &
                                 / ( (maxval([0.0_core_rknd,(ddzt_umvm_sqd(1,k))]) + 1.0e-5) &
                                     * (gr_dens_z(k) + 100) )
      else
        Lscale_term(k) = zero
        richardson_num_term(k) = zero
      end if

      ! Calculate chi_term
      chi_term(k) = exp(1000 * chi(k))/(gr_dens_z(k) + 1000.0)

      ! add calculated values to the cumulative sums
      cumulative_Lscale(k) = cumulative_Lscale(k) + Lscale_term(k)
      cumulative_chi(k) = cumulative_chi(k) + chi_term(k)
      cumulative_richardson_num(k) = cumulative_richardson_num(k) + richardson_num_term(k)

      ! Calculate time averaged Lscale_term, chi_term and richardson_num_term,
      ! by dividing the cumulative sum by the number of elements in that sum
      ! (if the grid is adapted this cumulative sum and the counter are set back to 0)
      Lscale_term_time_avg(k) = cumulative_Lscale(k) / Lscale_counter
      chi_term_time_avg(k) = cumulative_chi(k) / chi_counter
      richardson_num_term_time_avg(k) = cumulative_richardson_num(k) / richardson_num_counter
    
    end do
   
    ! Weigh each term to get the final density as sum of the individual terms
    do k = 1, nzm
        ! refinement criterion
        gr_dens(k) =   0.0   * Lscale_term_time_avg(k) &
                     + 3.0   * alt_term(k) &
                     + 100.0 * chi_term_time_avg(k) &
                     + 7.5   * richardson_num_term_time_avg(k)
    end do

    ! Smooth the grid density
    do k = 1, 3
      gr_dens = moving_average( nzm, &
                                gr_dens_z, gr_dens )
    end do

  end subroutine calc_grid_dens_helper

  function moving_average( nz, &
                           gr_dens_z, gr_dens )

    ! Description:
    ! Applies a weighted moving average filter to the given function

    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero, one_half, one

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      nz

    ! These two arrays define the piecewise linear grid density function
    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      gr_dens_z, & ! altitudes [m]
      gr_dens      ! densities [1/m]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(nz) ::  &
      moving_average ! the smoothed densitiy values given at the altitudes gr_dens_z [1/m]

    !--------------------- Local Variables ---------------------
    integer :: &
      window_size, &              ! size of the moving average filter
      start_index_full_window, &  ! used to store the index where the first
                                  ! moving average can be applied
      end_index_full_window, &    ! used to store the index where the last
                                  ! moving average can be applied
      k, j

    real( kind = core_rknd ) :: &
      avg, &                   ! moving average for every indexs
      sum_weights, &           ! sum of the moving average filter weights
      weight_for_middle_point  ! prescribed weight for the point in the middle

    real( kind = core_rknd ), dimension(:), allocatable :: &
      weights ! the weights for one filter, changes for each point

    !--------------------- Begin Code ---------------------
    window_size = 3 ! needs to be uneven, so on both sides equally many points are
                    ! taking into consideration

    allocate( weights(window_size) )

    weight_for_middle_point = one_half ! needs to be betweeen [0,1]

    
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

        weights(j-(k-(window_size-1)/2)+1) = zero
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
      avg = zero
      do j = k-(window_size-1)/2, k+(window_size-1)/2
        avg = avg + weights(j-(k-(window_size-1)/2)+1) * gr_dens(j)
      end do
      moving_average(k) = avg
    end do

    deallocate( weights )

  end function moving_average

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

    if ( allocated( cumulative_richardson_num ) ) then
      deallocate( cumulative_richardson_num )
    end if

    if ( allocated( cumulative_chi ) ) then
      deallocate( cumulative_chi )
    end if

  end subroutine clean_up_grid_adaptation_module

  subroutine adapt_grid( ngrdcol, gr_dens_norm_idx, &
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
                         Lscale, err_info )

    ! Description:
    ! Adapts the grid based on the density function and remaps all values to the new grid.
    !-----------------------------------------------------------------------
    
    use grid_class, only: &
        setup_grid ! Procedure

    use error_code, only: &
        clubb_fatal_error
    
    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        zero

    use err_info_type_module, only: &
        err_info_type        ! Type

    use constants_clubb, only: &
        fstderr ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      hydromet_dim, &
      sclr_dim, &
      edsclr_dim, &
      gr_dens_norm_idx, & ! total numer of indices of density_func_z and density_func_dens []
      idx_thvm   ! numer of indices of thvm

    real( kind = core_rknd ), dimension(gr_dens_norm_idx), intent(in) ::  &
      ! The two arrays min_gr_dens_norm_z and min_gr_dens_norm build the piecewise
      ! linear normalized minimum grid density.
      min_gr_dens_norm_z, & ! altitudes [m]
      min_gr_dens_norm, &   ! densities [1/m]
      ! The two arrays gr_dens_norm_z and gr_dens_norm build the piecewise
      ! linear normalized grid density.
      gr_dens_norm_z,  & ! altitudes [m]
      gr_dens_norm       ! densities [1/m]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      sfc_elevation, &   ! Elevation of ground level    [m AMSL]
      p_sfc              ! Pressure at the surface      [Pa]

    logical, intent(in) :: &
      l_implemented ! Parameter for srtting up the grid

    real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
      thvm  ! Virtual potential temperature             [K]

    integer, intent(in) :: &
      grid_remap_method ! specifies what remapping method should be used

    !--------------------- Output Variable ---------------------

    !--------------------- In/Out Variable ---------------------
    type( grid ), intent(inout) :: gr ! the current grid, all values are given on

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
      num_levels, &                ! number of levels the new grid should have     []
      total_idx_rho_lin_spline, &  ! total number of indices of the rho density
                                   ! piecewise linear function                     []
      grid_type, &                 ! the grid type for setting up the new grid     []
      i, k

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      deltaz ! grid spacings for setting up the new grid  [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) ::  &
      new_gr_zm  ! zm levels for the adapted grid based on the density function [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) ::  &
      thermodynamic_heights_placeholder  ! placeholder for the zt levels, but not needed

    ! These two arrays build the piecewise linear function of the rho values
    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      rho_lin_spline_vals, &  ! rho values [kg/m^3]
      rho_lin_spline_levels   ! altitudes  [m]
    
    type( grid ) :: new_gr ! the new grid object

    real( kind = core_rknd ) ::  &
      threshold ! some threshold to decide wether grid should be adapted or not

    logical, dimension(ngrdcol) :: l_adapt_grid

    !--------------------- Begin Code ---------------------
    num_levels = gr%nzm
    ! trigger threshold
    threshold = 6.0e-3

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

    call create_grid_from_normalized_grid_density( ngrdcol, &
                                                   gr_dens_norm_idx, &
                                                   gr_dens_norm_z, gr_dens_norm, &
                                                   min_gr_dens_norm_z, min_gr_dens_norm, &
                                                   num_levels, &
                                                   threshold, &
                                                   new_gr_zm, &
                                                   l_adapt_grid )
    

    do i = 1, ngrdcol
      deltaz(i) = zero
    end do
    grid_type = 3

    ! TODO adjust for ngrdcol > 1
    if ( l_adapt_grid(1) ) then
      call setup_grid( num_levels, ngrdcol, sfc_elevation, l_implemented, &          ! intent(in)
                       .true., grid_type, deltaz, gr%zm(:,1), gr%zm(:,num_levels), & ! intent(in)
                       new_gr_zm, thermodynamic_heights_placeholder, &               ! intent(in)
                       new_gr, err_info )                                            ! intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "Fatal error calling setup_grid in adapt_grid"
    end if

      Lscale_counter = 0
      richardson_num_counter = 0
      chi_counter = 0
      do k = 1, gr%nzm
        cumulative_Lscale(k) = zero
        cumulative_richardson_num(k) = zero
        cumulative_chi(k) = zero
      end do

      ! Set the density values to use for interpolation for mass calculation
      total_idx_rho_lin_spline = gr%nzm
      rho_lin_spline_vals = rho_ds_zm
      rho_lin_spline_levels = gr%zm

      call remap_all_clubb_core_vals( ngrdcol, total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, rho_lin_spline_levels, &
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

  subroutine remap_all_clubb_core_vals( ngrdcol, total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, rho_lin_spline_levels, &
                                        hydromet_dim, sclr_dim, edsclr_dim, &
                                        thvm, idx_thvm, p_sfc, &
                                        grid_remap_method, &
                                        gr_source, gr_target, &
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
    ! Remaps all values from the old to the new grid.
    !-----------------------------------------------------------------------

    use calc_pressure, only: &
        init_pressure  ! Procedure

    use remapping_module, only: &
        remap_vals_to_target ! Procedure

    use grid_class, only: &
        grid ! Type

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

    ! These two arrays build the piecewise linear function of the rho values
    real( kind = core_rknd ), dimension(ngrdcol, total_idx_rho_lin_spline), intent(in) ::  &
      rho_lin_spline_vals, &  ! rho values [kg/m^3]
      rho_lin_spline_levels   ! altitudes  [m]
    
    type( grid ), intent(in) :: &
      gr_source, & ! grid object values are currently given on
      gr_target    ! grid object where values should be remapped to

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      p_sfc ! Pressure at surface level [Pa]

    real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
      thvm  ! Virtual potential temperature      [K]

    integer, intent(in) :: &
      grid_remap_method ! specifies what remapping method should be used

    !--------------------- Output Variable ---------------------

    !--------------------- In/Out Variable ---------------------
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
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

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) ::  &
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

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzm, hydromet_dim), intent(inout) :: &
      wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzt, hydromet_dim), intent(inout) :: &
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and a hydrometeor [K <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzt), intent(inout) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) ::  &
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
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) ::  &
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) ::  &
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) :: &
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt,edsclr_dim) :: &
        edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr_source%nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzt), intent(inout) :: &
      Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzm), intent(inout) :: &
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzm), intent(inout) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr_source%nzt), intent(inout) :: &
      Lscale     ! Length scale         [m]

    !--------------------- Local Variables ---------------------
    integer :: i

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzt) ::  &
      thvm_new_grid

    real( kind = core_rknd ), dimension(ngrdcol, gr_source%nzm) ::  &
      p_in_Pa_zm, exner_zm ! just placeholder variables; are only needed for subroutine call

    integer :: &
      iv_wind = -1, &
      iv_pos_def = 0, &
      iv_other = 1

    logical :: l_zt_variable

    !--------------------- Begin Code ---------------------

    ! Remap all zt values
    l_zt_variable = .true.
    thlm_forcing = remap_vals_to_target( ngrdcol, &
                                         gr_source, gr_target, &
                                         gr_source%nzt, &
                                         thlm_forcing, &
                                         gr_target%nzt, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         iv_other, p_sfc, &
                                         grid_remap_method, &
                                         l_zt_variable )
    rtm_forcing = remap_vals_to_target( ngrdcol, &
                                        gr_source, gr_target, &
                                        gr_source%nzt, &
                                        rtm_forcing, &
                                        gr_target%nzt, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        iv_other, p_sfc, &
                                        grid_remap_method, &
                                        l_zt_variable )
    um_forcing = remap_vals_to_target( ngrdcol, &
                                       gr_source, gr_target, &
                                       gr_source%nzt, &
                                       um_forcing, &
                                       gr_target%nzt, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       iv_other, p_sfc, &
                                       grid_remap_method, &
                                       l_zt_variable )
    vm_forcing = remap_vals_to_target( ngrdcol, &
                                       gr_source, gr_target, &
                                       gr_source%nzt, &
                                       vm_forcing, &
                                       gr_target%nzt, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       iv_other, p_sfc, &
                                       grid_remap_method, &
                                       l_zt_variable )
    wm_zt = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  wm_zt, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_wind, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    rho = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                rho, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzt, &
                                      rho_ds_zt, &
                                      gr_target%nzt, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_pos_def, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    invrs_rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzt, &
                                            invrs_rho_ds_zt, &
                                            gr_target%nzt, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_pos_def, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
    thv_ds_zt = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzt, &
                                      thv_ds_zt, &
                                      gr_target%nzt, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_pos_def, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    rfrzm = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  rfrzm, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_pos_def, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    do i = 1, hydromet_dim
      wp2hmp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzt, &
                                            wp2hmp(:,:,i), &
                                            gr_target%nzt, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_other, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
      rtphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                               gr_source, gr_target, &
                                               gr_source%nzt, &
                                               rtphmp_zt(:,:,i), &
                                               gr_target%nzt, &
                                               total_idx_rho_lin_spline, &
                                               rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               iv_other, p_sfc, &
                                               grid_remap_method, &
                                               l_zt_variable )
      thlphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                gr_source, gr_target, &
                                                gr_source%nzt, &
                                                thlphmp_zt(:,:,i), &
                                                gr_target%nzt, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                iv_other, p_sfc, &
                                                grid_remap_method, &
                                                l_zt_variable )
    end do
    do i = 1, sclr_dim
      sclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                   gr_source, gr_target, &
                                                   gr_source%nzt, &
                                                   sclrm_forcing(:,:,i), &
                                                   gr_target%nzt, &
                                                   total_idx_rho_lin_spline, &
                                                   rho_lin_spline_vals, &
                                                   rho_lin_spline_levels, &
                                                   iv_other, p_sfc, &
                                                   grid_remap_method, &
                                                   l_zt_variable )
      sclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                           gr_source, gr_target, &
                                           gr_source%nzt, &
                                           sclrm(:,:,i), &
                                           gr_target%nzt, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           iv_other, p_sfc, &
                                           grid_remap_method, &
                                           l_zt_variable )
      sclrp3(:,:,i) = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzt, &
                                            sclrp3(:,:,i), &
                                            gr_target%nzt, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_other, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
    end do
    do i = 1, edsclr_dim
      edsclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                     gr_source, gr_target, &
                                                     gr_source%nzt, &
                                                     edsclrm_forcing(:,:,i), &
                                                     gr_target%nzt, &
                                                     total_idx_rho_lin_spline, &
                                                     rho_lin_spline_vals, &
                                                     rho_lin_spline_levels, &
                                                     iv_other, p_sfc, &
                                                     grid_remap_method, &
                                                     l_zt_variable )
      edsclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                             gr_source, gr_target, &
                                             gr_source%nzt, &
                                             edsclrm(:,:,i), &
                                             gr_target%nzt, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             iv_other, p_sfc, &
                                             grid_remap_method, &
                                             l_zt_variable )
    end do
    rtm_ref = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzt, &
                                    rtm_ref, &
                                    gr_target%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    thlm_ref = remap_vals_to_target( ngrdcol, &
                                     gr_source, gr_target, &
                                     gr_source%nzt, &
                                     thlm_ref, &
                                     gr_target%nzt, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     iv_pos_def, p_sfc, &
                                     grid_remap_method, &
                                     l_zt_variable )
    um_ref = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzt, &
                                   um_ref, &
                                   gr_target%nzt, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_wind, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    vm_ref = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzt, &
                                   vm_ref, &
                                   gr_target%nzt, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_wind, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    ug = remap_vals_to_target( ngrdcol, &
                               gr_source, gr_target, &
                               gr_source%nzt, &
                               ug, &
                               gr_target%nzt, &
                               total_idx_rho_lin_spline, &
                               rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               iv_wind, p_sfc, &
                               grid_remap_method, &
                               l_zt_variable )
    vg = remap_vals_to_target( ngrdcol, &
                               gr_source, gr_target, &
                               gr_source%nzt, &
                               vg, &
                               gr_target%nzt, &
                               total_idx_rho_lin_spline, &
                               rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               iv_wind, p_sfc, &
                               grid_remap_method, &
                               l_zt_variable )
    um = remap_vals_to_target( ngrdcol, &
                               gr_source, gr_target, &
                               gr_source%nzt, &
                               um, &
                               gr_target%nzt, &
                               total_idx_rho_lin_spline, &
                               rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               iv_wind, p_sfc, &
                               grid_remap_method, &
                               l_zt_variable )
    vm = remap_vals_to_target( ngrdcol, &
                               gr_source, gr_target, &
                               gr_source%nzt, &
                               vm, &
                               gr_target%nzt, &
                               total_idx_rho_lin_spline, &
                               rho_lin_spline_vals, &
                               rho_lin_spline_levels, &
                               iv_wind, p_sfc, &
                               grid_remap_method, &
                               l_zt_variable )
    up3 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                up3, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_other, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    vp3 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                vp3, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_other, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    rtm = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                rtm, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    thlm = remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 gr_source%nzt, &
                                 thlm, &
                                 gr_target%nzt, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv_pos_def, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )
    rtp3 = remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 gr_source%nzt, &
                                 rtp3, &
                                 gr_target%nzt, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv_other, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )
    thlp3 = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  thlp3, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    wp3 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                wp3, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_other, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    ! remap thvm values to new grid, to calculate new p_in_Pa
    thvm_new_grid = remap_vals_to_target( ngrdcol, &
                                          gr_source, gr_target, &
                                          gr_source%nzt, &
                                          thvm, &
                                          gr_target%nzt, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          iv_pos_def, p_sfc, &
                                          grid_remap_method, &
                                          l_zt_variable )
    ! calculate p_in_Pa instead of remapping directly since it can run into problems if for
    ! example the two highest levels have the same pressure value, which could be happening with
    ! the remapping, also calculate exner accordingly
    call init_pressure( ngrdcol, gr_target, thvm_new_grid, p_sfc, &
                        p_in_Pa, exner, p_in_Pa_zm, exner_zm )
    rcm = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzt, &
                                rcm, &
                                gr_target%nzt, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    cloud_frac = remap_vals_to_target( ngrdcol, &
                                       gr_source, gr_target, &
                                       gr_source%nzt, &
                                       cloud_frac, &
                                       gr_target%nzt, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       iv_pos_def, p_sfc, &
                                       grid_remap_method, &
                                       l_zt_variable )
    wp2thvp = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzt, &
                                    wp2thvp, &
                                    gr_target%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    wp2rtp = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzt, &
                                   wp2rtp, &
                                   gr_target%nzt, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_other, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    wp2thlp = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzt, &
                                    wp2thlp, &
                                    gr_target%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    wpup2 = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  wpup2, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    wpvp2 = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  wpvp2, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    ice_supersat_frac = remap_vals_to_target( ngrdcol, &
                                              gr_source, gr_target, &
                                              gr_source%nzt, &
                                              ice_supersat_frac, &
                                              gr_target%nzt, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              iv_pos_def, p_sfc, &
                                              grid_remap_method, &
                                              l_zt_variable )
    um_pert = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzt, &
                                    um_pert, &
                                    gr_target%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    vm_pert = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzt, &
                                    vm_pert, &
                                    gr_target%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    rcm_in_layer = remap_vals_to_target( ngrdcol, &
                                         gr_source, gr_target, &
                                         gr_source%nzt, &
                                         rcm_in_layer, &
                                         gr_target%nzt, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         iv_pos_def, p_sfc, &
                                         grid_remap_method, &
                                         l_zt_variable )
    cloud_cover = remap_vals_to_target( ngrdcol, &
                                        gr_source, gr_target, &
                                        gr_source%nzt, &
                                        cloud_cover, &
                                        gr_target%nzt, &
                                        total_idx_rho_lin_spline, &
                                        rho_lin_spline_vals, &
                                        rho_lin_spline_levels, &
                                        iv_pos_def, p_sfc, &
                                        grid_remap_method, &
                                        l_zt_variable )
    w_up_in_cloud = remap_vals_to_target( ngrdcol, &
                                          gr_source, gr_target, &
                                          gr_source%nzt, &
                                          w_up_in_cloud, &
                                          gr_target%nzt, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          iv_other, p_sfc, &
                                          grid_remap_method, &
                                          l_zt_variable )
    w_down_in_cloud = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzt, &
                                            w_down_in_cloud, &
                                            gr_target%nzt, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_other, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
    cloudy_updraft_frac = remap_vals_to_target( ngrdcol, &
                                                gr_source, gr_target, &
                                                gr_source%nzt, &
                                                cloudy_updraft_frac, &
                                                gr_target%nzt, &
                                                total_idx_rho_lin_spline, &
                                                rho_lin_spline_vals, &
                                                rho_lin_spline_levels, &
                                                iv_pos_def, p_sfc, &
                                                grid_remap_method, &
                                                l_zt_variable )
    cloudy_downdraft_frac = remap_vals_to_target( ngrdcol, &
                                                  gr_source, gr_target, &
                                                  gr_source%nzt, &
                                                  cloudy_downdraft_frac, &
                                                  gr_target%nzt, &
                                                  total_idx_rho_lin_spline, &
                                                  rho_lin_spline_vals, &
                                                  rho_lin_spline_levels, &
                                                  iv_pos_def, p_sfc, &
                                                  grid_remap_method, &
                                                  l_zt_variable )
    Kh_zt = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzt, &
                                  Kh_zt, &
                                  gr_target%nzt, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_pos_def, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    Lscale = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzt, &
                                   Lscale, &
                                   gr_target%nzt, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_pos_def, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )

    ! Remap all zm values
    l_zt_variable = .false.
    wprtp_forcing = remap_vals_to_target( ngrdcol, &
                                          gr_source, gr_target, &
                                          gr_source%nzm, &
                                          wprtp_forcing, &
                                          gr_target%nzm, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          iv_other, p_sfc, &
                                          grid_remap_method, &
                                          l_zt_variable )
    wpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                           gr_source, gr_target, &
                                           gr_source%nzm, &
                                           wpthlp_forcing, &
                                           gr_target%nzm, &
                                           total_idx_rho_lin_spline, &
                                           rho_lin_spline_vals, &
                                           rho_lin_spline_levels, &
                                           iv_other, p_sfc, &
                                           grid_remap_method, &
                                           l_zt_variable )
    rtp2_forcing = remap_vals_to_target( ngrdcol, &
                                         gr_source, gr_target, &
                                         gr_source%nzm, &
                                         rtp2_forcing, &
                                         gr_target%nzm, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         iv_other, p_sfc, &
                                         grid_remap_method, &
                                         l_zt_variable )
    thlp2_forcing = remap_vals_to_target( ngrdcol, &
                                          gr_source, gr_target, &
                                          gr_source%nzm, &
                                          thlp2_forcing, &
                                          gr_target%nzm, &
                                          total_idx_rho_lin_spline, &
                                          rho_lin_spline_vals, &
                                          rho_lin_spline_levels, &
                                          iv_other, p_sfc, &
                                          grid_remap_method, &
                                          l_zt_variable )
    rtpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzm, &
                                            rtpthlp_forcing, &
                                            gr_target%nzm, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_other, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
    wm_zm = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  wm_zm, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_wind, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    rho_zm = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzm, &
                                   rho_zm, &
                                   gr_target%nzm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_pos_def, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzm, &
                                      rho_ds_zm, &
                                      gr_target%nzm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_pos_def, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    invrs_rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzm, &
                                            invrs_rho_ds_zm, &
                                            gr_target%nzm, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_pos_def, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
    thv_ds_zm = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzm, &
                                      thv_ds_zm, &
                                      gr_target%nzm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_pos_def, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    do i = 1, hydromet_dim
      wphydrometp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                                 gr_source, gr_target, &
                                                 gr_source%nzm, &
                                                 wphydrometp(:,:,i), &
                                                 gr_target%nzm, &
                                                 total_idx_rho_lin_spline, &
                                                 rho_lin_spline_vals, &
                                                 rho_lin_spline_levels, &
                                                 iv_other, p_sfc, &
                                                 grid_remap_method, &
                                                 l_zt_variable )
    end do
    upwp = remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 gr_source%nzm, &
                                 upwp, &
                                 gr_target%nzm, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv_other, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )
    vpwp = remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 gr_source%nzm, &
                                 vpwp, &
                                 gr_target%nzm, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv_other, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )
    up2 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzm, &
                                up2, &
                                gr_target%nzm, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    vp2 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzm, &
                                vp2, &
                                gr_target%nzm, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    wprtp = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  wprtp, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    wpthlp = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzm, &
                                   wpthlp, &
                                   gr_target%nzm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_other, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    rtp2 = remap_vals_to_target( ngrdcol, &
                                 gr_source, gr_target, &
                                 gr_source%nzm, &
                                 rtp2, &
                                 gr_target%nzm, &
                                 total_idx_rho_lin_spline, &
                                 rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 iv_pos_def, p_sfc, &
                                 grid_remap_method, &
                                 l_zt_variable )
    thlp2 = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  thlp2, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_pos_def, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    rtpthlp = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzm, &
                                    rtpthlp, &
                                    gr_target%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    wp2 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzm, &
                                wp2, &
                                gr_target%nzm, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    do i = 1, sclr_dim
      wpsclrp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                             gr_source, gr_target, &
                                             gr_source%nzm, &
                                             wpsclrp(:,:,i), &
                                             gr_target%nzm, &
                                             total_idx_rho_lin_spline, &
                                             rho_lin_spline_vals, &
                                             rho_lin_spline_levels, &
                                             iv_other, p_sfc, &
                                             grid_remap_method, &
                                             l_zt_variable )
      sclrp2(:,:,i) = remap_vals_to_target( ngrdcol, &
                                            gr_source, gr_target, &
                                            gr_source%nzm, &
                                            sclrp2(:,:,i), &
                                            gr_target%nzm, &
                                            total_idx_rho_lin_spline, &
                                            rho_lin_spline_vals, &
                                            rho_lin_spline_levels, &
                                            iv_other, p_sfc, &
                                            grid_remap_method, &
                                            l_zt_variable )
      sclrprtp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                              gr_source, gr_target, &
                                              gr_source%nzm, &
                                              sclrprtp(:,:,i), &
                                              gr_target%nzm, &
                                              total_idx_rho_lin_spline, &
                                              rho_lin_spline_vals, &
                                              rho_lin_spline_levels, &
                                              iv_other, p_sfc, &
                                              grid_remap_method, &
                                              l_zt_variable )
      sclrpthlp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                               gr_source, gr_target, &
                                               gr_source%nzm, &
                                               sclrpthlp(:,:,i), &
                                               gr_target%nzm, &
                                               total_idx_rho_lin_spline, &
                                               rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               iv_other, p_sfc, &
                                               grid_remap_method, &
                                               l_zt_variable )
      sclrpthvp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                               gr_source, gr_target, &
                                               gr_source%nzm, &
                                               sclrpthvp(:,:,i), &
                                               gr_target%nzm, &
                                               total_idx_rho_lin_spline, &
                                               rho_lin_spline_vals, &
                                               rho_lin_spline_levels, &
                                               iv_other, p_sfc, &
                                               grid_remap_method, &
                                               l_zt_variable )
    end do
    wpthvp = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzm, &
                                   wpthvp, &
                                   gr_target%nzm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_other, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    rtpthvp = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzm, &
                                    rtpthvp, &
                                    gr_target%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
    thlpthvp = remap_vals_to_target( ngrdcol, &
                                     gr_source, gr_target, &
                                     gr_source%nzm, &
                                     thlpthvp, &
                                     gr_target%nzm, &
                                     total_idx_rho_lin_spline, &
                                     rho_lin_spline_vals, &
                                     rho_lin_spline_levels, &
                                     iv_other, p_sfc, &
                                     grid_remap_method, &
                                     l_zt_variable )
    uprcp = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  uprcp, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    vprcp = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  vprcp, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    rc_coef_zm = remap_vals_to_target( ngrdcol, &
                                       gr_source, gr_target, &
                                       gr_source%nzm, &
                                       rc_coef_zm, &
                                       gr_target%nzm, &
                                       total_idx_rho_lin_spline, &
                                       rho_lin_spline_vals, &
                                       rho_lin_spline_levels, &
                                       iv_pos_def, p_sfc, &
                                       grid_remap_method, &
                                       l_zt_variable )
    wp4 = remap_vals_to_target( ngrdcol, &
                                gr_source, gr_target, &
                                gr_source%nzm, &
                                wp4, &
                                gr_target%nzm, &
                                total_idx_rho_lin_spline, &
                                rho_lin_spline_vals, &
                                rho_lin_spline_levels, &
                                iv_pos_def, p_sfc, &
                                grid_remap_method, &
                                l_zt_variable )
    wp2up2 = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzm, &
                                   wp2up2, &
                                   gr_target%nzm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_pos_def, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    wp2vp2 = remap_vals_to_target( ngrdcol, &
                                   gr_source, gr_target, &
                                   gr_source%nzm, &
                                   wp2vp2, &
                                   gr_target%nzm, &
                                   total_idx_rho_lin_spline, &
                                   rho_lin_spline_vals, &
                                   rho_lin_spline_levels, &
                                   iv_pos_def, p_sfc, &
                                   grid_remap_method, &
                                   l_zt_variable )
    upwp_pert = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzm, &
                                      upwp_pert, &
                                      gr_target%nzm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_other, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    vpwp_pert = remap_vals_to_target( ngrdcol, &
                                      gr_source, gr_target, &
                                      gr_source%nzm, &
                                      vpwp_pert, &
                                      gr_target%nzm, &
                                      total_idx_rho_lin_spline, &
                                      rho_lin_spline_vals, &
                                      rho_lin_spline_levels, &
                                      iv_other, p_sfc, &
                                      grid_remap_method, &
                                      l_zt_variable )
    wprcp = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  wprcp, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_other, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    invrs_tau_zm = remap_vals_to_target( ngrdcol, &
                                         gr_source, gr_target, &
                                         gr_source%nzm, &
                                         invrs_tau_zm, &
                                         gr_target%nzm, &
                                         total_idx_rho_lin_spline, &
                                         rho_lin_spline_vals, &
                                         rho_lin_spline_levels, &
                                         iv_pos_def, p_sfc, &
                                         grid_remap_method, &
                                         l_zt_variable )
    Kh_zm = remap_vals_to_target( ngrdcol, &
                                  gr_source, gr_target, &
                                  gr_source%nzm, &
                                  Kh_zm, &
                                  gr_target%nzm, &
                                  total_idx_rho_lin_spline, &
                                  rho_lin_spline_vals, &
                                  rho_lin_spline_levels, &
                                  iv_pos_def, p_sfc, &
                                  grid_remap_method, &
                                  l_zt_variable )
    thlprcp = remap_vals_to_target( ngrdcol, &
                                    gr_source, gr_target, &
                                    gr_source%nzm, &
                                    thlprcp, &
                                    gr_target%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

  end subroutine remap_all_clubb_core_vals

!===============================================================================

end module grid_adaptation_module
