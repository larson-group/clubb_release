!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_adaptation_module


  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: setup_gr_fixed_min, &
            setup_gr_dycore, &
            normalize_grid_density, &
            adapt_grid, &
            calc_grid_dens, &                   ! Procedures
            lscale_term_weight, &
            alt_term_weight, &
            chi_term_weight, &
            richardson_num_term_weight          ! Parameters

  private :: calc_integral, &
             normalize_min_grid_density, normalize_grid_density_helper, &
             decide_if_grid_adapt, apply_equidistribution_principle, &
             create_grid_from_normalized_grid_density, check_grid, &
             calc_grid_dens_helper, moving_average, remap_all_clubb_core_vals

  private ! Default Scoping

  real( kind = core_rknd ), parameter ::  &
    tol = 1.0e-6_core_rknd ! tolerance to check if to real numbers are equal

  ! weights for each term in the refinement criterion
  real( kind = core_rknd ), parameter ::  &
    lscale_term_weight = 0.0, &
    alt_term_weight = 3.0, &
    chi_term_weight = 100.0, &
    richardson_num_term_weight = 7.5

  contains

  !=============================================================================

  ! TODO find better way to set these up, since both do the same thing, just with different zm
  ! file paths, and maybe some solution where parameters are not hardcoded, so maybe something
  ! where the parameters are read from model.in file and dont use hardcoded 73 in dimension of
  ! heights vectors
  subroutine setup_gr_dycore( iunit, ngrdcol, grid_sfc, grid_top, l_ascending_grid, gr, err_info )

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
    integer, intent(in) :: &
      iunit, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      grid_sfc, &    ! grids surface height for the dycore grid [m]
      grid_top       ! grids highest level for the dycore grid  [m]

    logical, intent(in) :: &
      l_ascending_grid

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- In/Output Variables ---------------------

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Local Variables ---------------------

    integer :: &
      i, k ! Column and vertical loop indices [-]
    integer, dimension(ngrdcol) :: &
      nlevel, & ! Last input level, then the retained level count [-]
      level_min ! First retained input level [-]

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
    do i = 1, ngrdcol
      do while( momentum_heights(i,nlevel(i)) > grid_top(i) .and. nlevel(i) > 1 )
        nlevel(i) = nlevel(i) - 1
      end do
    end do

    do i = 1, ngrdcol
      momentum_heights(i,nlevel(i)) = grid_top(i)
    end do

    level_min = nlevel
    do i = 1, ngrdcol
      do while( momentum_heights(i,level_min(i)) > grid_sfc(i) .and. level_min(i) > 1 )
        level_min(i) = level_min(i) - 1
      end do
    end do

    do i = 1, ngrdcol
      momentum_heights(i,level_min(i)) = grid_sfc(i)
    end do
    
    nlevel = nlevel - level_min + 1
    if ( any(nlevel /= nlevel(1)) ) then
      write(fstderr,*) "Fixed grids require the same number of levels in every column"
      err_info%err_code = clubb_fatal_error
      return
    end if

    do k = 1, nlevel(1)
      do i = 1, ngrdcol
        momentum_heights(i,k) = momentum_heights(i,level_min(i)+k-1)
      end do
    end do

    call setup_grid( nlevel(1), ngrdcol, sfc_elevation, l_implemented, &           ! In
                     l_ascending_grid, grid_type, deltaz, grid_sfc, grid_top, & ! In
                     momentum_heights, thermodynamic_heights, &                 ! In
                     gr, err_info )                                             ! InOut

    if ( clubb_at_least_debug_level_api(0) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid in setup_gr_dycore"
      end if
    end if

  end subroutine setup_gr_dycore

  subroutine setup_gr_fixed_min( iunit, ngrdcol, &
                                 grid_first, grid_last, &
                                 l_ascending_grid, &
                                 gr, err_info )

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
    integer, intent(in) :: &
      iunit, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      grid_first, &    ! grids first momentum height for the fixed minimum grid [m]
      grid_last        ! grids last momentum level for the fixed minimum grid  [m]

    logical, intent(in) :: &
      l_ascending_grid

    !--------------------- Output Variables ---------------------
    type (grid), intent(out) :: gr

    !--------------------- Input/Output Variables ---------------------
    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    !--------------------- Local Variables ---------------------

    integer :: &
      i, k ! Column and vertical loop indices [-]
    integer, dimension(ngrdcol) :: &
      nlevel, & ! Last input level, then the retained level count [-]
      level_min ! First retained input level [-]

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

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      grid_sfc, &    ! grids surface height for the dycore grid [m]
      grid_top       ! grids highest level for the dycore grid  [m]

    ! these are the grid_sfc/grid_top variable -/+ a small value, to ensure that the
    ! setup_grid subroutine always includes the exact grid_sfc and grid_top values
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      grid_sfc_outer, &   !  [m]
      grid_top_outer      !  [m]

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
      write(fstderr, *) "Fatal error calling read_grid_heights in setup_gr_fixed_min"
    end if

    if ( l_ascending_grid ) then
      do i = 1, ngrdcol
        grid_sfc(i) = grid_first(i)
        grid_top(i) = grid_last(i)
      end do
    else
      do i = 1, ngrdcol
        grid_sfc(i) = grid_last(i)
        grid_top(i) = grid_first(i)
      end do
    end if

    nlevel = nzmax
    do i = 1, ngrdcol
      do while( momentum_heights(i,nlevel(i)) > grid_top(i) .and. nlevel(i) > 1 )
        nlevel(i) = nlevel(i) - 1
      end do
    end do

    if ( any(momentum_heights(:,nzmax) < grid_top) .or. any(momentum_heights(:,1) > grid_sfc) ) then
      write(fstderr, *) "Error in setup_gr_fixed_min in grid_adaptation_module."
      write(fstderr, *) "The fixed minimum grid density momentum heights"
      write(fstderr, *) "need to cover at least the physics-grid domain."
    end if

    do i = 1, ngrdcol
      momentum_heights(i,nlevel(i)) = grid_top(i)
    end do

    level_min = nlevel
    do i = 1, ngrdcol
      do while( momentum_heights(i,level_min(i)) > grid_sfc(i) .and. level_min(i) > 1 )
        level_min(i) = level_min(i) - 1
      end do
    end do

    do i = 1, ngrdcol
      momentum_heights(i,level_min(i)) = grid_sfc(i)
    end do

    nlevel = nlevel - level_min + 1
    if ( any(nlevel /= nlevel(1)) ) then
      write(fstderr,*) "Fixed grids require the same number of levels in every column"
      err_info%err_code = clubb_fatal_error
      return
    end if

    do k = 1, nlevel(1)
      do i = 1, ngrdcol
        momentum_heights(i,k) = momentum_heights(i,level_min(i)+k-1)
      end do
    end do

    do i = 1, ngrdcol
      grid_sfc_outer(i) = grid_sfc(i) - 0.001
      grid_top_outer(i) = grid_top(i) + 0.001
    end do
    
    call setup_grid( nlevel(1), ngrdcol, sfc_elevation, l_implemented, &           ! In
                     l_ascending_grid, grid_type, deltaz, &                     ! In
                     grid_sfc_outer, grid_top_outer, &                          ! In
                     momentum_heights, thermodynamic_heights, &                 ! In
                     gr, err_info )                                             ! InOut

    if ( any(err_info%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info%err_header_global
      write(fstderr, *) "Fatal error calling setup_grid in setup_gr_fixed_min"
    end if

  end subroutine setup_gr_fixed_min

  function calc_integral( ngrdcol, g_idx, g_x, g_y )
    ! Description:
    ! Computes the exact integral, assuming the given points build a piecewise linear function.
    ! Each column of g_x and g_y has g_idx elements.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
      zero

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol ! Number of grid columns
    integer, intent(in) :: & 
      g_idx  ! The total numer of indices of g_x and g_y []

    real( kind = core_rknd ), dimension(ngrdcol,g_idx), intent(in) ::  &
      g_x,  &    ! the x coordinates of the connection points of the piecewise linear function [m]
      g_y        ! the y coordinates of the connection points of the piecewise linear function

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      calc_integral ! The integral

    !--------------------- Local Variables ---------------------
    integer :: i, k

    !--------------------- Begin Code ---------------------
    ! Initializing calc_integral to avoid a compiler warning.
    calc_integral = zero

    ! Calculate each linear segments exact integral and add them
    do k = 1, g_idx-1
      do i = 1, ngrdcol
        calc_integral(i) = calc_integral(i) + (g_y(i,k) + g_y(i,k+1))/2 * (g_x(i,k+1) - g_x(i,k))
      end do
    end do

  end function calc_integral

  function normalize_min_grid_density( ngrdcol, min_gr_dens_idx, min_gr_dens_z, & ! In
                                       min_gr_dens, & ! In
                                       lambda, & ! In
                                       num_levels ) result( min_gr_dens_norm ) ! In

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

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      integral, &
      norm_factor

    !--------------------- Begin Code ---------------------
    ! check if lambda is in (0,1]
    if ( lambda <= zero .or. lambda > one ) then
      error stop 'lambda needs to be between 0 (excluded) and 1 (included)'
    end if

    ! Normalize every column without changing its vertical summation order.
    integral = calc_integral( ngrdcol, min_gr_dens_idx, min_gr_dens_z, min_gr_dens )
    norm_factor = lambda*(num_levels - 1)/integral
    do j = 1, min_gr_dens_idx
      do i = 1, ngrdcol
        min_gr_dens_norm(i,j) = norm_factor(i)*min_gr_dens(i,j)
      end do
    end do

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      ! check if integral is lambda*(n-1) as wanted
      integral = calc_integral( ngrdcol, min_gr_dens_idx, min_gr_dens_z, min_gr_dens_norm )
      if ( any( abs(integral - (lambda*(num_levels-1))) > tol ) ) then
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

  function normalize_grid_density_helper( ngrdcol, gr_dens_idx, gr_dens_z, gr_dens, & ! In
    min_gr_dens_norm, & ! In
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

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      integral_gr_dens, &
      integral_min_gr_dens, &
      integral_test, &
      norm_factor, & ! Density scaling factor in each column [-]
      shift          ! Uniform density added in each column [1/m]

    real( kind = core_rknd ) :: &
      grid_sfc, &
      grid_top

    !--------------------- Begin Code ---------------------
    integral_gr_dens = calc_integral( ngrdcol, gr_dens_idx, gr_dens_z, gr_dens )
    integral_min_gr_dens = calc_integral( ngrdcol, gr_dens_idx, gr_dens_z, min_gr_dens_norm )

    do i = 1, ngrdcol
      if ( abs(integral_min_gr_dens(i) - (lambda*(num_levels-1))) > tol ) then
        write(fstderr,*) 'Warning! The minimum grid density function has not the correct integral.'
        error stop 'Normalize the minimum grid density function before using it.'
      end if

      if ( integral_gr_dens(i) < 1e-8_core_rknd ) then
        ! this happens if the prescribed non-normalized grid density was zero at all altitudes

        ! calculate shift value
        grid_sfc = gr_dens_z(i,1)
        grid_top = gr_dens_z(i,gr_dens_idx)
        shift(i) = (num_levels-1 - integral_min_gr_dens(i))/(grid_top - grid_sfc)

      else
        ! calculate factor for normalization
        norm_factor(i) = (num_levels-1 - integral_min_gr_dens(i))/integral_gr_dens(i)

      end if

    end do ! i

    do j = 1, gr_dens_idx
      do i = 1, ngrdcol
        if ( integral_gr_dens(i) < 1e-8_core_rknd ) then
          gr_dens_norm(i,j) = min_gr_dens_norm(i,j) + shift(i)
        else
          gr_dens_norm(i,j) = norm_factor(i)*gr_dens(i,j) + min_gr_dens_norm(i,j)
        end if
      end do
    end do

    if ( clubb_at_least_debug_level_api( 2 ) ) then
      do j = 1, gr_dens_idx
        do i = 1, ngrdcol
          if ( gr_dens_norm(i,j) < min_gr_dens_norm(i,j) .and. &
            (min_gr_dens_norm(i,j) - gr_dens_norm(i,j)) > tol ) then
            write(fstderr,*) 'Normalized function was below minimum grid density function.'
            write(fstderr,*) 'Normalized minimum grid density function is ', &
              min_gr_dens_norm(i,j), &
              ' while normalized grid density function is ', gr_dens_norm(i,j)
            error stop 'Normalized function should be above or equal to the minimum density.'
          end if
        end do
      end do
    end if

    if ( clubb_at_least_debug_level_api( 2 ) ) then
      integral_test = calc_integral( ngrdcol, gr_dens_idx, gr_dens_z, gr_dens_norm )
      if ( any( abs(integral_test - (num_levels-1)) > tol ) ) then
        error stop 'Something went wrong, integral is different than it should be.'
      end if
    end if

  end function normalize_grid_density_helper

  subroutine normalize_grid_density( ngrdcol, gr_dens_idx, num_levels, gr_dens_z, gr_dens, & ! In
                                     lambda, & ! In
                                     gr_fixed_min, & ! In
                                     norm_min_grid_dens, norm_grid_dens ) ! Out
    ! Description:
    ! Normalizes the prescribed grid density.

    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use interpolation, only: &
        zlinterp

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: & 
      ngrdcol, &
      gr_dens_idx, &     ! total numer of indices of gr_dens_z and gr_dens []
      num_levels         ! number of levels the new grid should have       []

    real( kind = core_rknd ), intent(in) :: &
      lambda   ! a factor to adjust the integral of the minimum grid density

    ! These two arrays build the piecewise linear prescribed non-normalized grid density
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(in) ::  &
      gr_dens_z,  &   ! the z coordinates of the connection points of the piecewise linear
                      ! grid density function                                                  [m]
      gr_dens         ! the values of the density function at the given z coordinates of the
                      ! connection points of the piecewise linear grid density function        [1/m]

    type( grid ), intent(in) :: &
      gr_fixed_min

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx), intent(out) :: &
      norm_min_grid_dens, & ! the density at the given z coordinates of the connection points
                            ! of the normalized piecewise linear grid density function         [1/m]
      norm_grid_dens        ! the density at the given z coordinates of the connection points
                            ! of the normalized piecewise linear grid density function         [1/m]

    !--------------------- Local Variables ---------------------

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_idx) ::  &
      min_gr_dens        ! the density at the given z coordinates of the connection points
                         ! of the normalized piecewise linear grid density function          [1/m]

    !--------------------- Begin Code ---------------------

    ! set the minimum grid density profile to be the linear piecewise function of the
    ! original minimum grid density function evaluated at the current grid levels
    ! this way the function can be executed with all checks and works out exactly
    min_gr_dens = zlinterp( ngrdcol, gr_dens_idx, & ! In
                            gr_fixed_min%nzm, gr_dens_z, & ! In
                            gr_fixed_min%zm, gr_fixed_min%invrs_dzm ) ! In

    norm_min_grid_dens = normalize_min_grid_density( ngrdcol, gr_dens_idx, gr_dens_z, & ! In
                                                     min_gr_dens, & ! In
                                                     lambda, & ! In
                                                     num_levels ) ! In

    norm_grid_dens = normalize_grid_density_helper( ngrdcol, gr_dens_idx, gr_dens_z, gr_dens, & ! In
                                                    norm_min_grid_dens, & ! In
                                                    lambda, num_levels ) ! In

  end subroutine normalize_grid_density
  
  function decide_if_grid_adapt( ngrdcol, gr_dens_old_idx, gr_dens_old_z, & ! In
                                 gr_dens_old, & ! In
                                 gr_dens_new_idx, gr_dens_new_z, gr_dens_new, & ! In
                                 threshold ) result( l_adapt_grid ) ! In

    ! Description:
    ! Checks how different the two noramlized grid densities are and returns depending
    ! on the threshold if the grid should be adapted or not.

    !-----------------------------------------------------------------------

    use interpolation, only: &
        zlinterp

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
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_old_idx + gr_dens_new_idx) ::  &
      common_grid ! the grid tha results from combining the grid the old grid density was
                  ! given on and the grid the current grid density is given on              [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_old_idx + gr_dens_new_idx) ::  &
      gr_dens_old_interp ! old normalized grid density interpolated to the common grid  [1/m]

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_old_idx + gr_dens_new_idx) ::  &
      gr_dens_new_interp ! new normalized grid density interpolated to the common grid  [1/m]

    integer :: i, j, k ! Column, merged-level and reduction indices [-]
    integer :: old_cursor, new_cursor ! Next level in each input grid [-]
    integer :: common_count(ngrdcol) ! Active levels of each merged grid

    real( kind = core_rknd ) :: &
      sum(ngrdcol) ! Sum of normalized squared density differences in each column [-]

    !--------------------- Begin Code ---------------------
    common_count = 0

    ! Interior heights differ, so merge each column's two grids independently.
    do i = 1, ngrdcol
      old_cursor = 1
      new_cursor = 1
      do while ( old_cursor < gr_dens_old_idx .and. new_cursor < gr_dens_new_idx )
        j = common_count(i) + 1

        if ( gr_dens_old_z(i,old_cursor) < gr_dens_new_z(i,new_cursor) - tol ) then
          common_grid(i,j) = gr_dens_old_z(i,old_cursor)
          old_cursor = old_cursor + 1
        else if ( gr_dens_new_z(i,new_cursor) < gr_dens_old_z(i,old_cursor) - tol ) then
          common_grid(i,j) = gr_dens_new_z(i,new_cursor)
          new_cursor = new_cursor + 1
        else
          common_grid(i,j) = ( gr_dens_old_z(i,old_cursor) + gr_dens_new_z(i,new_cursor) )/2
          old_cursor = old_cursor + 1
          new_cursor = new_cursor + 1
        end if

        common_count(i) = j
      end do
    end do

    gr_dens_old_interp = zlinterp( ngrdcol, gr_dens_old_idx + gr_dens_new_idx, & ! In
                                   gr_dens_old_idx, common_grid, & ! In
                                   gr_dens_old_z, gr_dens_old, & ! In
                                   common_count ) ! Optional In
    gr_dens_new_interp = zlinterp( ngrdcol, gr_dens_old_idx + gr_dens_new_idx, & ! In
                                   gr_dens_new_idx, common_grid, & ! In
                                   gr_dens_new_z, gr_dens_new, & ! In
                                   common_count ) ! Optional In

    ! Retain the merged-level reduction order in each column.
    sum = 0
    do k = 1, gr_dens_old_idx + gr_dens_new_idx
      do i = 1, ngrdcol
        if ( k > common_count(i) ) cycle
        sum(i) = sum(i) + ( gr_dens_old_interp(i,k) - gr_dens_new_interp(i,k) &
          )**2/gr_dens_old_interp(i,k)**2
      end do
    end do

    do i = 1, ngrdcol
      l_adapt_grid(i) = sum(i)/common_count(i) > threshold
    end do

  end function decide_if_grid_adapt

  function apply_equidistribution_principle( ngrdcol, num_levels, gr_dens_norm_idx, & ! In
    gr_dens_norm_z, & ! In
                                             gr_dens_norm, l_adapt_grid ) result( grid_heights )
    ! Description:
    ! Redistribute interior levels independently in each selected column.
    ! The level count and domain endpoints remain fixed.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol ! Number of grid columns
    integer, intent(in) :: & 
      gr_dens_norm_idx, &  ! The total numer of indices of gr_dens_norm_z and gr_dens_norm []
      num_levels           ! number of levels the new grid should have []

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(in) ::  &
      gr_dens_norm_z,  &  ! the z coordinates of the connection points of the normalized
                          ! piecewise linear grid density function                          [m]
      gr_dens_norm        ! the density values at the given z coordinates of the
                          ! density function connection points of the normalized
                          ! piecewise linear grid                                           [1/m]
      ! The grid density function needs to be normalized with the number of grid levels num_levels,
      ! before the function can be applied here

    logical, intent(in) :: l_adapt_grid(ngrdcol) ! Columns selected for adaptation

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_levels) :: &
      grid_heights ! the heights of the newly created grid, from bottom to top [m]

    !--------------------- Local Variables ---------------------
    integer :: &
      k ! Vertical level index [-]

    integer :: prev_x_ind ! Last density coordinate used in the current column [-]

    real( kind = core_rknd ) :: &
      grid_sfc, &               ! the surface height of the grid              [m]
      grid_top, &               ! the highest point of the grid               [m]
      new_x_area_increment, &   ! integral from grid density
                                ! function on [prev_x_ind, prev_x_ind+1]      [# levs]
      desired_area_up_to_ith_level, & ! the desired area up to the ith grid level, should be exactly
                                      ! i, since the grid density function is normalized so that the
                                      ! integral of the whole function is num_levels-1 and since we
                                      ! follow the equi-distribution approach [# levs]
      A, grid_level, slope, intercept, p_pq_formula, q_pq_formula

    integer :: i ! Column index

    logical :: l_grid_level_found ! Whether the current interior level is located [-]

    real( kind = core_rknd ) :: &
      area_up_to_prev_x ! Integral up to the last density coordinate [# levs]

    !--------------------- Begin Code ---------------------
    grid_heights = 0.0_core_rknd
    do i = 1, ngrdcol
      if ( .not. l_adapt_grid(i) ) cycle
      grid_sfc = gr_dens_norm_z(i,1)
      grid_top = gr_dens_norm_z(i,gr_dens_norm_idx)

      grid_heights(i,1) = grid_sfc

      ! Start the density search at the surface of each active column.
      prev_x_ind = 1
      area_up_to_prev_x = 0.0_core_rknd

      grid_heights(i,num_levels) = grid_top

      do k = 2, num_levels-1
        l_grid_level_found = .false.
        do while ( .not. l_grid_level_found )

          new_x_area_increment = (gr_dens_norm(i,prev_x_ind) + &
            gr_dens_norm(i,prev_x_ind+1))/2 &
            * (gr_dens_norm_z(i,prev_x_ind+1) - gr_dens_norm_z(i,prev_x_ind))
          desired_area_up_to_ith_level = k-1
          if (area_up_to_prev_x + new_x_area_increment >= desired_area_up_to_ith_level) then
            A = desired_area_up_to_ith_level - area_up_to_prev_x
            ! check if gr_dens_norm(i,prev_x_ind+1) == gr_dens_norm(i,prev_x_ind)
            if (abs(gr_dens_norm(i,prev_x_ind+1) - gr_dens_norm(i,prev_x_ind)) < tol) then
              ! gr_dens_norm(i,prev_x_ind+1) == gr_dens_norm(i,prev_x_ind)
              grid_level = A/gr_dens_norm(i,prev_x_ind) + gr_dens_norm_z(i,prev_x_ind)
            else
              ! calculate slope and intercept of the linear function
              slope = (gr_dens_norm(i,prev_x_ind+1) - gr_dens_norm(i,prev_x_ind)) &
                /(gr_dens_norm_z(i,prev_x_ind+1) - gr_dens_norm_z(i,prev_x_ind))
              intercept = gr_dens_norm(i,prev_x_ind) - gr_dens_norm_z(i,prev_x_ind)*slope
              ! calculate p and q of second-order polynomial to solve for x
              ! solve x^2+px+q=0 with x=-p/2(+/-)sqrt((p/2)^2-q)
              p_pq_formula = 2*intercept/slope
              q_pq_formula = -1*gr_dens_norm_z(i,prev_x_ind)**2 &
                - 2/slope*intercept*gr_dens_norm_z(i,prev_x_ind) - 2/slope*A
              ! since we have second-order polynomial, we have two roots, so we need to
              ! identify the right solution - use pq formula to get solutions
              grid_level = -1*p_pq_formula/2 + sqrt((p_pq_formula/2)**2 - q_pq_formula)

              if ( ( grid_level < gr_dens_norm_z(i,prev_x_ind) &
                .and. (gr_dens_norm_z(i,prev_x_ind) - grid_level) > tol ) &
                .or. &
                ( grid_level > gr_dens_norm_z(i,prev_x_ind+1) &
                .and. (grid_level - gr_dens_norm_z(i,prev_x_ind+1)) > tol ) &
                .or. (grid_level < grid_heights(i,k-1) ) ) then
                ! first solution was not the one we were looking for
                grid_level = -1*p_pq_formula/2 - sqrt((p_pq_formula/2)**2 - q_pq_formula)

                if ( ( grid_level < gr_dens_norm_z(i,prev_x_ind) &
                  .and. (gr_dens_norm_z(i,prev_x_ind) - grid_level) > tol ) &
                  .or. &
                  ( grid_level > gr_dens_norm_z(i,prev_x_ind+1) &
                  .and. (grid_level - gr_dens_norm_z(i,prev_x_ind+1)) > tol ) &
                  .or. (grid_level < grid_heights(i,k-1) ) ) then
                  write(fstderr,*) "None of the two solutions works. ", &
                    "Something went wrong."
                  error stop 'Something went wrong with generating the grid'
                end if
              end if
            end if

            grid_heights(i,k) = grid_level
            l_grid_level_found = .true.
          else
            area_up_to_prev_x = area_up_to_prev_x + new_x_area_increment
            prev_x_ind = prev_x_ind + 1
          end if
        end do
      end do

    end do

  end function apply_equidistribution_principle

  subroutine create_grid_from_normalized_grid_density( ngrdcol, gr_dens_norm_idx, num_levels, & ! In
                                                       gr_dens_norm_z, gr_dens_norm, & ! In
                                                       min_gr_dens_norm_z, min_gr_dens_norm, & ! In
                                                       threshold, & ! In
                                                       gr_dens_old_z, gr_dens_old, & ! InOut
                                                       grid_heights, l_adapt_grid ) ! Out
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

    !--------------------- Input/Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(inout) ::  &
      ! The two arrays gr_dens_old_z and gr_dens_old build the piecewise linear
      ! normalized grid density, from the last time the grid was adapted
      gr_dens_old_z,  & ! altitudes [m]
      gr_dens_old       ! densities [1/m]

    !--------------------- Local Variables ---------------------
    integer :: grid_heights_idx, i, k ! Grid size, column and level indices [-]

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      grid_sfc,  &    ! height of the grids surface   [m]
      grid_top        ! height of the top of the grid [m]

    !--------------------- Begin Code ---------------------
    grid_sfc = gr_dens_norm_z(:,1)
    grid_top = gr_dens_norm_z(:,gr_dens_norm_idx)

    ! decide if grid should be adapted or not
    l_adapt_grid = decide_if_grid_adapt( ngrdcol, gr_dens_norm_idx, gr_dens_old_z, & ! In
                                         gr_dens_old, & ! In
                                         gr_dens_norm_idx, gr_dens_norm_z, gr_dens_norm, & ! In
                                         threshold ) ! In

    ! create grid from normalized grid density function
    do k = 1, gr_dens_norm_idx
      do i = 1, ngrdcol
        if ( l_adapt_grid(i) ) then
          ! Save the density used by the next adaptation decision.
          gr_dens_old_z(i,k) = gr_dens_norm_z(i,k)
          gr_dens_old(i,k) = gr_dens_norm(i,k)
        end if
      end do
    end do
    grid_heights = apply_equidistribution_principle( ngrdcol, num_levels, gr_dens_norm_idx, & ! In
                                                     gr_dens_norm_z, & ! In
                                                     gr_dens_norm, l_adapt_grid ) ! In

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      grid_heights_idx = size(grid_heights, 2)

      call check_grid( ngrdcol, grid_heights_idx, num_levels, gr_dens_norm_idx, & ! In
                       l_adapt_grid, grid_heights, min_gr_dens_norm_z, min_gr_dens_norm, & ! In
                       grid_sfc, grid_top ) ! In

    end if

  end subroutine create_grid_from_normalized_grid_density

  subroutine check_grid( ngrdcol, grid_heights_idx, desired_num_levels, desired_min_dens_idx, & ! In
                         l_adapt_grid, grid_heights, desired_min_dens_z, desired_min_dens, & ! In
                         desired_grid_sfc, desired_grid_top ) ! In

    ! Description:
    ! Checks if the grid is valid.

    !-----------------------------------------------------------------------

    use interpolation, only: &
        zlinterp ! Procedure

    use constants_clubb, only: &
        one, fstderr ! Constants

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol ! Number of grid columns
    integer, intent(in) :: &
      grid_heights_idx, &      ! number of elements in grid_heights                            []
      desired_min_dens_idx, &  ! number of elements in desired_min_dens_z and desired_min_dens []
      desired_num_levels       ! desired number of grid levels                                 []
    ! Note: if everything worked fine, those two integers should be the same

    real( kind = core_rknd ), dimension(ngrdcol,grid_heights_idx), intent(in) :: &
      grid_heights ! the grid heights ordered from bottom to top [m]

    real( kind = core_rknd ), dimension(ngrdcol,desired_min_dens_idx), intent(in) :: &
      desired_min_dens_z, & ! the z coordinates of the piecewise linear desired
                            ! minimum grid density function                         [m]
      desired_min_dens      ! the density values on the given z coordinates         [1/m]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      desired_grid_sfc, &   ! desired surface height of the grid [m]
      desired_grid_top      ! desired top height of the grid     [m]

    logical, intent(in) :: l_adapt_grid(ngrdcol) ! Columns whose grid was rebuilt

    !--------------------- Local Variables ---------------------

    integer :: i, j, k

    real( kind = core_rknd ) :: min_desired_min_dens

    real( kind = core_rknd ), dimension(ngrdcol,desired_num_levels) :: desired_min_dens_on_new_gr

    !--------------------- Begin Code ---------------------
    desired_min_dens_on_new_gr = zlinterp( ngrdcol, desired_num_levels, & ! In
                                           desired_min_dens_idx, grid_heights, & ! In
                                           desired_min_dens_z, desired_min_dens ) ! In
    do i = 1, ngrdcol
      if ( .not. l_adapt_grid(i) ) cycle
      ! check if first grid element is grid_sfc
      if (abs(grid_heights(i,1) - desired_grid_sfc(i)) > tol) then
        write(fstderr,*) "WARNING! The first grid level is not correct. It should be the grids", &
        " surface: ", desired_grid_sfc(i), " and not ", grid_heights(i,1)
        error stop 'Something went wrong with generating the grid'
      end if

      do j = 1, (grid_heights_idx-1)
        ! check if grid heights are getting bigger
        if (grid_heights(i,j) >= grid_heights(i,j+1)) then
          write(fstderr,*) "WARNING! The grid levels are not in the right order. ", &
          "They should get bigger with the index."
          error stop 'Something went wrong with generating the grid'
        end if

        ! check if grid fulfills min_dens
        ! find minimal min_dens in grid cell
        ! check if minimum is in cell or on the boundaries
        ! first check values in interval
        min_desired_min_dens = maxval(desired_min_dens_z(i,:))
        do k = 1, desired_min_dens_idx
          if ( (desired_min_dens_z(i,k) > grid_heights(i,j)) &
            .and. (desired_min_dens_z(i,k) < grid_heights(i,j+1)) ) then
            if ( desired_min_dens(i,k) < min_desired_min_dens ) then
              min_desired_min_dens = desired_min_dens(i,k)
            end if
          end if
        end do

        ! then if maximum is on boundary
        if ( desired_min_dens_on_new_gr(i,j) < min_desired_min_dens ) then
          min_desired_min_dens = desired_min_dens_on_new_gr(i,j)
        end if
        if ( desired_min_dens_on_new_gr(i,j+1) < min_desired_min_dens ) then
          min_desired_min_dens = desired_min_dens_on_new_gr(i,j+1)
        end if

        if ((grid_heights(i,j+1) - grid_heights(i,j)) > one/(min_desired_min_dens) + tol) then
          write(fstderr,*) "WARNING! The grid has a distance=", &
          (grid_heights(i,j+1) - grid_heights(i,j)), &
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
      if (abs(grid_heights(i,desired_num_levels) - desired_grid_top(i)) > tol) then
        write(fstderr,*) "WARNING! The last grid level is not correct. It should be the grids", &
        " top: ", desired_grid_top(i), " and not ", grid_heights(i,desired_num_levels)
        error stop 'Something went wrong with generating the grid'
      end if

    end do

  end subroutine check_grid

  subroutine calc_grid_dens( ngrdcol, gr, & ! In
                             um, vm, Lscale_zt, wp2, pdf_params, thlm, exner, rtm, rcm, & ! In
                             p_in_Pa, thvm, ice_supersat_frac, bv_efold, T0, & ! In
                             clubb_config_flags, & ! In
                             Lscale_counter, chi_counter, richardson_num_counter, & ! InOut
                             cumulative_Lscale, cumulative_chi, cumulative_richardson_num, & ! InOut
                             gr_dens_z, gr_dens, alt_term, Lscale_term, Lscale_term_time_avg, & ! Out
                             chi_term, chi_term_time_avg, richardson_num_term, & ! Out
                             richardson_num_term_time_avg ) ! Out
    ! Description:
    ! Calculate the necessary variables and then construct the non-normalized
    ! grid density independently in each column.

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

    real( kind = core_rknd ), intent(in) :: &
      T0                     ! Reference temperature (usually 300)  [K]

    type(pdf_parameter), intent(in) :: & 
      pdf_params           ! PDF parameters

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags   ! Derived type holding all configurable CLUBB flags

    !--------------------- Input/Output Variable ---------------------
    ! Save the calculated values for the refinement criterion terms as a cumulative sum starting
    ! from the last grid adaptation, to take the time average over that region, to have a
    ! smoother evolution in time
    integer, dimension(ngrdcol), intent(inout) :: &
      Lscale_counter, &          ! Steps accumulated for the mixing-length term [-]
      chi_counter, &             ! Steps accumulated for the chi term [-]
      richardson_num_counter     ! Steps accumulated for the Richardson-number term [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      cumulative_Lscale, &          ! Sum of mixing-length refinement terms [1/m]
      cumulative_chi, &             ! Sum of chi refinement terms [1/m]
      cumulative_richardson_num     ! Sum of Richardson-number refinement terms [1/m]

    !--------------------- Output Variable ---------------------
    ! These two arrays build the piecewise linear grid density
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) ::  &
      gr_dens_z,  &    ! altitudes [m]
      gr_dens          ! densities [1/m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) :: &
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

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      chi_zt   ! The variable 's' in Mellor (1977)    [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol, gr%nzm) ::  &
      Lscale,                       & ! Length scale                                 [m]
      ddzt_um,                      & ! Vertical derivative of um                    [s^-1]
      ddzt_vm,                      & ! Vertical derivative of vm                    [s^-1]
      ddzt_umvm_sqd,                & ! Squared vertical norm of derivative of
                                      ! mean horizontal wind speed                   [s^-2]
      brunt_vaisala_freq_sqd,       & ! Buoyancy frequency squared, N^2              [s^-2]
      brunt_vaisala_freq_sqd_mixed, & ! A mixture of dry and moist N^2               [s^-2]
      brunt_vaisala_freq_sqd_smth     ! Mix between dry and moist N^2 that is
                                      ! smoothed in the vertical                     [s^-2]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
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
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        chi_zt(i,k) = pdf_params%mixt_frac(i,k) * pdf_params%chi_1(i,k) &
          + ( one - pdf_params%mixt_frac(i,k) ) * pdf_params%chi_2(i,k)
      end do
    end do
    chi = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, chi_zt )

    ! Calculate brunt_vaisala_freq_sqd
    call calc_brunt_vaisala_freq_sqd( gr%nzm, gr%nzt, ngrdcol, gr, thlm,                  & ! In
      exner, rtm, rcm, p_in_Pa, thvm,                     & ! In
      ice_supersat_frac,                                  & ! In
      clubb_config_flags%saturation_formula,              & ! In
      clubb_config_flags%l_brunt_vaisala_freq_moist,      & ! In
      clubb_config_flags%l_use_thvm_in_bv_freq,           & ! In
      clubb_config_flags%l_modify_limiters_for_cnvg_test, & ! In
      bv_efold, T0,                                       & ! In
      brunt_vaisala_freq_sqd,                             & ! Out
      brunt_vaisala_freq_sqd_mixed,                       & ! Out
      brunt_vaisala_freq_sqd_smth )                         ! Out

    ! Calculate the grid density from Lscale, chi, ddzt_umvm_sqd and brunt_vaisala_freq_sqd
    call calc_grid_dens_helper( ngrdcol, gr%nzm, gr%zm, Lscale, wp2, chi, & ! In
                                brunt_vaisala_freq_sqd, & ! In
                                ddzt_umvm_sqd, & ! In
                                Lscale_counter, chi_counter, richardson_num_counter, & ! InOut
                                cumulative_Lscale, cumulative_chi, & ! InOut
                                cumulative_richardson_num, & ! InOut
                                gr_dens_z, gr_dens, alt_term, Lscale_term, & ! Out
                                Lscale_term_time_avg, & ! Out
                                chi_term, chi_term_time_avg, richardson_num_term, & ! Out
                                richardson_num_term_time_avg ) ! Out

  end subroutine calc_grid_dens

  subroutine calc_grid_dens_helper( ngrdcol, nzm, zm, Lscale, wp2, chi, & ! In
                                    brunt_vaisala_freq_sqd, & ! In
                                    ddzt_umvm_sqd, & ! In
                                    Lscale_counter, chi_counter, richardson_num_counter, & ! InOut
                                    cumulative_Lscale, cumulative_chi, & ! InOut
                                    cumulative_richardson_num, & ! InOut
                                    gr_dens_z, gr_dens, alt_term, Lscale_term, & ! Out
                                    Lscale_term_time_avg, & ! Out
                                    chi_term, chi_term_time_avg, richardson_num_term, & ! Out
                                    richardson_num_term_time_avg ) ! Out
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

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      chi ! The variable 's' in Mellor (1977) given on zm levels    [kg/kg]

    !--------------------- Input/Output Variable ---------------------
    ! Save the calculated values for the refinement criterion terms as a cumulative sum starting
    ! from the last grid adaptation, to take the time average over that region, to have a
    ! smoother evolution in time
    integer, dimension(ngrdcol), intent(inout) :: &
      Lscale_counter, &          ! Steps accumulated for the mixing-length term [-]
      chi_counter, &             ! Steps accumulated for the chi term [-]
      richardson_num_counter     ! Steps accumulated for the Richardson-number term [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      cumulative_Lscale, &          ! Sum of mixing-length refinement terms [1/m]
      cumulative_chi, &             ! Sum of chi refinement terms [1/m]
      cumulative_richardson_num     ! Sum of Richardson-number refinement terms [1/m]

    !--------------------- Output Variable ---------------------
    ! These two arrays build the piecewise linear non-normalized grid density
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      gr_dens_z,  &    ! altitudes [m]
      gr_dens          ! densities [1/m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
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

    real( kind = core_rknd ), parameter :: &
      wp2_threshold = 0.001 ! threshold for wp2 that decides whether Lscale is used or not [m^2/s^2]

    !--------------------- Begin Code ---------------------

    Lscale_counter = Lscale_counter + 1
    richardson_num_counter = richardson_num_counter + 1
    chi_counter = chi_counter + 1

    ! Calculate each refinement criterion term
    do k = 1, nzm
      do i = 1, ngrdcol
        gr_dens_z(i,k) = zm(i,k)

        ! Calculate alt_term
        alt_term(i,k) = 1/(gr_dens_z(i,k) + 100)

        ! Calculate Lscale_term and richardson number term
        if ( wp2(i,k) > wp2_threshold ) then
          Lscale_term(i,k) = 1/(Lscale(i,k)+10.0)
          richardson_num_term(i,k) = maxval([0.0_core_rknd,(brunt_vaisala_freq_sqd(i,k))]) &
            / ( (maxval([0.0_core_rknd,(ddzt_umvm_sqd(i,k))]) + 1.0e-5) &
            * (gr_dens_z(i,k) + 100) )
        else
          Lscale_term(i,k) = zero
          richardson_num_term(i,k) = zero
        end if

        ! Calculate chi_term
        chi_term(i,k) = exp(1000 * chi(i,k))/(gr_dens_z(i,k) + 1000.0)

        ! add calculated values to the cumulative sums
        cumulative_Lscale(i,k) = cumulative_Lscale(i,k) + Lscale_term(i,k)
        cumulative_chi(i,k) = cumulative_chi(i,k) + chi_term(i,k)
        cumulative_richardson_num(i,k) = cumulative_richardson_num(i,k) + richardson_num_term(i,k)

        ! Calculate time averaged Lscale_term, chi_term and richardson_num_term,
        ! by dividing the cumulative sum by the number of elements in that sum
        ! (if the grid is adapted this cumulative sum and the counter are set back to 0)
        Lscale_term_time_avg(i,k) = cumulative_Lscale(i,k) / Lscale_counter(i)
        chi_term_time_avg(i,k) = cumulative_chi(i,k) / chi_counter(i)
        richardson_num_term_time_avg(i,k) &
          = cumulative_richardson_num(i,k) / richardson_num_counter(i)

      end do
    end do

    ! Weigh each term to get the final density as sum of the individual terms
    do k = 1, nzm
      do i = 1, ngrdcol
        ! refinement criterion
        gr_dens(i,k) =   lscale_term_weight           * Lscale_term_time_avg(i,k) &
          + alt_term_weight              * alt_term(i,k) &
          + chi_term_weight              * chi_term_time_avg(i,k) &
          + richardson_num_term_weight   * richardson_num_term_time_avg(i,k)
      end do
    end do

    ! Smooth the grid density
    do k = 1, 3
      gr_dens = moving_average( ngrdcol, nzm, gr_dens_z, gr_dens ) ! In
    end do

  end subroutine calc_grid_dens_helper

  function moving_average( ngrdcol, nz, gr_dens_z, gr_dens )

    ! Description:
    ! Applies a weighted moving average filter to the given function

    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero, one_half, one

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol ! Number of grid columns
    integer, intent(in) :: & 
      nz

    ! These two arrays define the piecewise linear grid density function
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      gr_dens_z, & ! altitudes [m]
      gr_dens      ! densities [1/m]

    !--------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      moving_average ! the smoothed densitiy values given at the altitudes gr_dens_z [1/m]

    !--------------------- Local Variables ---------------------
    integer :: &
      window_size, &              ! size of the moving average filter
      start_index_full_window, &  ! used to store the index where the first
                                  ! moving average can be applied
      end_index_full_window, &    ! used to store the index where the last
                                  ! moving average can be applied
      i, k, j

    real( kind = core_rknd ) :: &
      weight_for_middle_point  ! prescribed weight for the point in the middle

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      weights ! Filter weights for each column [-]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      avg, &         ! Weighted density sum in each column [1/m]
      sum_weights    ! Sum of unnormalized filter weights in each column [1/m]

    !--------------------- Begin Code ---------------------
    window_size = 3 ! needs to be uneven, so on both sides equally many points are
    ! taking into consideration

    allocate( weights(ngrdcol,window_size) )

    weight_for_middle_point = one_half ! needs to be betweeen [0,1]


    start_index_full_window = (window_size-1)/2 + 1 ! is how many points on each side are
    ! taking into consideration

    end_index_full_window = nz - (window_size-1)/2  ! is how many points on each side are
    ! taking into consideration

    ! Endpoints do not have a full filter window.
    do k = 1, start_index_full_window - 1
      do i = 1, ngrdcol
        moving_average(i,k) = gr_dens(i,k)
      end do
    end do

    do k = end_index_full_window+1, nz
      do i = 1, ngrdcol
        moving_average(i,k) = gr_dens(i,k)
      end do
    end do

    do k = start_index_full_window, end_index_full_window
      do j = k-(window_size-1)/2, k+(window_size-1)/2
        do i = 1, ngrdcol
          weights(i,j-(k-(window_size-1)/2)+1) = zero
          if ( j /= k ) then
            weights(i,j-(k-(window_size-1)/2)+1) = 1/abs( gr_dens_z(i,k) - gr_dens_z(i,j) )
          end if
        end do
      end do

      ! Retain the original reduction order within each filter window.
      do i = 1, ngrdcol
        sum_weights(i) = sum( weights(i,:) )
      end do

      do j = k-(window_size-1)/2, k+(window_size-1)/2
        do i = 1, ngrdcol
          weights(i,j-(k-(window_size-1)/2)+1) = (1-weight_for_middle_point)/sum_weights(i) &
            * weights(i,j-(k-(window_size-1)/2)+1)
        end do
      end do

      do i = 1, ngrdcol
        weights(i,(window_size-1)/2 + 1) = weight_for_middle_point
        if ( (sum(weights(i,:)) - one) > tol ) then
          error stop 'Sum of weights in moving_average must be 1.0.'
        end if
      end do

      avg = zero
      do j = k-(window_size-1)/2, k+(window_size-1)/2
        do i = 1, ngrdcol
          avg(i) = avg(i) + weights(i,j-(k-(window_size-1)/2)+1) * gr_dens(i,j)
        end do
      end do

      do i = 1, ngrdcol
        moving_average(i,k) = avg(i)
      end do
    end do

    deallocate( weights )

  end function moving_average

  subroutine adapt_grid( ngrdcol, gr_dens_norm_idx, gr_dens_norm_z, gr_dens_norm, & ! In
                         min_gr_dens_norm_z, & ! In
                         min_gr_dens_norm, sfc_elevation, l_implemented, hydromet_dim, & ! In
                         sclr_dim, & ! In
                         edsclr_dim, thvm, idx_thvm, p_sfc, & ! In
                         grid_remap_method, & ! In
                         gr, gr_dens_old_z, gr_dens_old, Lscale_counter, chi_counter, & ! InOut
                         richardson_num_counter, & ! InOut
                         cumulative_Lscale, cumulative_chi, & ! In
                         cumulative_richardson_num, & ! In
                         thlm_forcing, rtm_forcing, um_forcing, & ! InOut
                         vm_forcing, sclrm_forcing, edsclrm_forcing, wprtp_forcing, & ! InOut
                         wpthlp_forcing, & ! InOut
                         rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wm_zm, wm_zt, rtm_ref, & ! InOut
                         thlm_ref, um_ref, vm_ref, ug, vg, p_in_Pa, rho_zm, rho, exner, & ! InOut
                         rho_ds_zm, & ! InOut
                         rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! InOut
                         rfrzm, & ! InOut
                         wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, um, vm, upwp, vpwp, up2, & ! InOut
                         vp2, & ! InOut
                         up3, vp3, thlm, rtm, wprtp, wpthlp, wp2, wp3, rtp2, rtp3, thlp2, & ! InOut
                         thlp3, & ! InOut
                         rtpthlp, sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, wpsclrp, & ! InOut
                         edsclrm, & ! InOut
                         rcm, & ! InOut
                         cloud_frac, wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, sclrpthvp, & ! InOut
                         wp2rtp, & ! InOut
                         wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, wpup2, wpvp2, wp2up2, & ! InOut
                         wp2vp2, & ! InOut
                         ice_supersat_frac, um_pert, vm_pert, upwp_pert, vpwp_pert, Kh_zm, & ! InOut
                         Kh_zt, & ! InOut
                         thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, & ! InOut
                         cloudy_downdraft_frac, rcm_in_layer, cloud_cover, invrs_tau_zm, & ! InOut
                         Lscale, & ! InOut
                         err_info ) ! InOut

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

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(in) ::  &
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
      l_implemented ! True if CLUBB is being implemented and run in a host model

    real( kind = core_rknd ), dimension(ngrdcol, idx_thvm), intent(in) ::  &
      thvm  ! Virtual potential temperature             [K]

    integer, intent(in) :: &
      grid_remap_method ! specifies what remapping method should be used

    !--------------------- Output Variable ---------------------

    !--------------------- In/Out Variable ---------------------
    type( grid ), intent(inout) :: gr ! the current grid, all values are given on

    real( kind = core_rknd ), dimension(ngrdcol,gr_dens_norm_idx), intent(inout) ::  &
      ! The two arrays gr_dens_old_z and gr_dens_old build the piecewise linear
      ! normalized grid density, from the last time the grid was adapted
      gr_dens_old_z,  & ! altitudes [m]
      gr_dens_old       ! densities [1/m]

    ! Save the calculated values for the refinement criterion terms as a cumulative sum starting
    ! from the last grid adaptation, to take the time average over that region, to have a
    ! smoother evolution in time
    integer, dimension(ngrdcol), intent(inout) :: &
      Lscale_counter, &          ! the number of values cumulated in cumulative_Lscale
      richardson_num_counter, &  ! the number of values cumulated in cumulative_richardson_num
      chi_counter                ! the number of values cumulated in cumulative_chi
  
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      cumulative_Lscale, &          ! cumulative sum of the Lscale term in the
                                    ! refinement cirterion                                  [1/m]
      cumulative_richardson_num, &  ! cumulative sum of the richardson number term in the
                                    ! refinement cirterion                                  [1/m]
      cumulative_chi                ! cumulative sum of the chi term in the
                                    ! refinement cirterion                                  [1/m]

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
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      wp2up         ! < w'^2 u' > (thermodynamic levels)             [m^3/s^3]

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
      i, j, k                    ! Column, interpolation-weight and level indices [-]

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

    call create_grid_from_normalized_grid_density( ngrdcol, gr_dens_norm_idx, num_levels, & ! In
                                                   gr_dens_norm_z, gr_dens_norm, & ! In
                                                   min_gr_dens_norm_z, min_gr_dens_norm, & ! In
                                                   threshold, & ! In
                                                   gr_dens_old_z, gr_dens_old, & ! InOut
                                                   new_gr_zm, l_adapt_grid ) ! Out
    

    do i = 1, ngrdcol
      deltaz(i) = zero
    end do
    grid_type = 3

    if ( any(l_adapt_grid) ) then

      ! Unchanged columns retain their grid while selected columns are rebuilt.
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          if ( .not. l_adapt_grid(i) ) new_gr_zm(i,k) = gr%zm(i,k)
        end do
      end do
      call setup_grid( num_levels, ngrdcol, sfc_elevation, l_implemented, &          ! In
                       .true., grid_type, deltaz, gr%zm(:,1), gr%zm(:,num_levels), & ! In
                       new_gr_zm, thermodynamic_heights_placeholder, &               ! In
                       new_gr, err_info )                                            ! InOut

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling setup_grid in adapt_grid"
        return
      end if

      ! A retained thermodynamic grid need not use momentum-level midpoints.
      ! Preserve its geometry and interpolation weights along with its fields.
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          if ( l_adapt_grid(i) ) cycle
          new_gr%zt(i,k) = gr%zt(i,k)
          new_gr%dzt(i,k) = gr%dzt(i,k)
          new_gr%invrs_dzt(i,k) = gr%invrs_dzt(i,k)
        end do
      end do

      do k = 1, gr%nzm
        do i = 1, ngrdcol
          if ( l_adapt_grid(i) ) cycle
          new_gr%zm(i,k) = gr%zm(i,k)
          new_gr%dzm(i,k) = gr%dzm(i,k)
          new_gr%invrs_dzm(i,k) = gr%invrs_dzm(i,k)
        end do
      end do

      do j = lbound(gr%weights_zm2zt,3), ubound(gr%weights_zm2zt,3)
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( .not. l_adapt_grid(i) ) new_gr%weights_zm2zt(i,k,j) = gr%weights_zm2zt(i,k,j)
          end do
        end do
      end do

      do j = lbound(gr%weights_zt2zm,3), ubound(gr%weights_zt2zm,3)
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( .not. l_adapt_grid(i) ) new_gr%weights_zt2zm(i,k,j) = gr%weights_zt2zm(i,k,j)
          end do
        end do
      end do

      where ( l_adapt_grid )
        Lscale_counter = 0
        richardson_num_counter = 0
        chi_counter = 0
      end where

      do k = 1, gr%nzm
        do i = 1, ngrdcol
          if ( .not. l_adapt_grid(i) ) cycle
          cumulative_Lscale(i,k) = zero
          cumulative_richardson_num(i,k) = zero
          cumulative_chi(i,k) = zero
        end do
      end do

      ! Set the density values to use for interpolation for mass calculation
      total_idx_rho_lin_spline = gr%nzm
      rho_lin_spline_vals = rho_ds_zm
      rho_lin_spline_levels = gr%zm

      call remap_all_clubb_core_vals( ngrdcol, l_adapt_grid, total_idx_rho_lin_spline, &
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
                                      wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
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

  subroutine remap_all_clubb_core_vals( ngrdcol, l_adapt_grid, total_idx_rho_lin_spline, &
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
                                        wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &
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

    logical, intent(in) :: l_adapt_grid(ngrdcol) ! Columns selected for remapping [-]

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
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      wp2up         ! < w'^2 u' > (thermodynamic levels)             [m^3/s^3]

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
      thvm_new_grid, & ! Virtual potential temperature on the target grid [K]
      p_in_Pa_new, &   ! Recomputed pressure in adapting columns [Pa]
      exner_new        ! Recomputed Exner function in adapting columns [-]

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
    thlm_forcing = remap_vals_to_target( ngrdcol, & ! In
                                         gr_source, gr_target, & ! In
                                         gr_source%nzt, & ! In
                                         thlm_forcing, & ! In
                                         gr_target%nzt, & ! In
                                         total_idx_rho_lin_spline, & ! In
                                         rho_lin_spline_vals, & ! In
                                         rho_lin_spline_levels, & ! In
                                         iv_other, p_sfc, & ! In
                                         grid_remap_method, & ! In
                                         l_zt_variable, & ! In
                                         l_adapt_grid ) ! Optional In
    rtm_forcing = remap_vals_to_target( ngrdcol, & ! In
                                        gr_source, gr_target, & ! In
                                        gr_source%nzt, & ! In
                                        rtm_forcing, & ! In
                                        gr_target%nzt, & ! In
                                        total_idx_rho_lin_spline, & ! In
                                        rho_lin_spline_vals, & ! In
                                        rho_lin_spline_levels, & ! In
                                        iv_other, p_sfc, & ! In
                                        grid_remap_method, & ! In
                                        l_zt_variable, & ! In
                                        l_adapt_grid ) ! Optional In
    um_forcing = remap_vals_to_target( ngrdcol, & ! In
                                       gr_source, gr_target, & ! In
                                       gr_source%nzt, & ! In
                                       um_forcing, & ! In
                                       gr_target%nzt, & ! In
                                       total_idx_rho_lin_spline, & ! In
                                       rho_lin_spline_vals, & ! In
                                       rho_lin_spline_levels, & ! In
                                       iv_other, p_sfc, & ! In
                                       grid_remap_method, & ! In
                                       l_zt_variable, & ! In
                                       l_adapt_grid ) ! Optional In
    vm_forcing = remap_vals_to_target( ngrdcol, & ! In
                                       gr_source, gr_target, & ! In
                                       gr_source%nzt, & ! In
                                       vm_forcing, & ! In
                                       gr_target%nzt, & ! In
                                       total_idx_rho_lin_spline, & ! In
                                       rho_lin_spline_vals, & ! In
                                       rho_lin_spline_levels, & ! In
                                       iv_other, p_sfc, & ! In
                                       grid_remap_method, & ! In
                                       l_zt_variable, & ! In
                                       l_adapt_grid ) ! Optional In
    wm_zt = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  wm_zt, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_wind, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    rho = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                rho, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    rho_ds_zt = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzt, & ! In
                                      rho_ds_zt, & ! In
                                      gr_target%nzt, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_pos_def, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    invrs_rho_ds_zt = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzt, & ! In
                                            invrs_rho_ds_zt, & ! In
                                            gr_target%nzt, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_pos_def, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid ) ! Optional In
    thv_ds_zt = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzt, & ! In
                                      thv_ds_zt, & ! In
                                      gr_target%nzt, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_pos_def, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    rfrzm = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  rfrzm, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_pos_def, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    do i = 1, hydromet_dim
      wp2hmp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzt, & ! In
                                            wp2hmp(:,:,i), & ! In
                                            gr_target%nzt, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_other, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid )
      rtphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                               gr_source, gr_target, & ! In
                                               gr_source%nzt, & ! In
                                               rtphmp_zt(:,:,i), & ! In
                                               gr_target%nzt, & ! In
                                               total_idx_rho_lin_spline, & ! In
                                               rho_lin_spline_vals, & ! In
                                               rho_lin_spline_levels, & ! In
                                               iv_other, p_sfc, & ! In
                                               grid_remap_method, & ! In
                                               l_zt_variable, & ! In
                                               l_adapt_grid )
      thlphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                                gr_source, gr_target, & ! In
                                                gr_source%nzt, & ! In
                                                thlphmp_zt(:,:,i), & ! In
                                                gr_target%nzt, & ! In
                                                total_idx_rho_lin_spline, & ! In
                                                rho_lin_spline_vals, & ! In
                                                rho_lin_spline_levels, & ! In
                                                iv_other, p_sfc, & ! In
                                                grid_remap_method, & ! In
                                                l_zt_variable, & ! In
                                                l_adapt_grid )
    end do
    do i = 1, sclr_dim
      sclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                                   gr_source, gr_target, & ! In
                                                   gr_source%nzt, & ! In
                                                   sclrm_forcing(:,:,i), & ! In
                                                   gr_target%nzt, & ! In
                                                   total_idx_rho_lin_spline, & ! In
                                                   rho_lin_spline_vals, & ! In
                                                   rho_lin_spline_levels, & ! In
                                                   iv_other, p_sfc, & ! In
                                                   grid_remap_method, & ! In
                                                   l_zt_variable, & ! In
                                                   l_adapt_grid )
      sclrm(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                           gr_source, gr_target, & ! In
                                           gr_source%nzt, & ! In
                                           sclrm(:,:,i), & ! In
                                           gr_target%nzt, & ! In
                                           total_idx_rho_lin_spline, & ! In
                                           rho_lin_spline_vals, & ! In
                                           rho_lin_spline_levels, & ! In
                                           iv_other, p_sfc, & ! In
                                           grid_remap_method, & ! In
                                           l_zt_variable, & ! In
                                           l_adapt_grid )
      sclrp3(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzt, & ! In
                                            sclrp3(:,:,i), & ! In
                                            gr_target%nzt, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_other, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid )
    end do
    do i = 1, edsclr_dim
      edsclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                                     gr_source, gr_target, & ! In
                                                     gr_source%nzt, & ! In
                                                     edsclrm_forcing(:,:,i), & ! In
                                                     gr_target%nzt, & ! In
                                                     total_idx_rho_lin_spline, & ! In
                                                     rho_lin_spline_vals, & ! In
                                                     rho_lin_spline_levels, & ! In
                                                     iv_other, p_sfc, & ! In
                                                     grid_remap_method, & ! In
                                                     l_zt_variable, & ! In
                                                     l_adapt_grid )
      edsclrm(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                             gr_source, gr_target, & ! In
                                             gr_source%nzt, & ! In
                                             edsclrm(:,:,i), & ! In
                                             gr_target%nzt, & ! In
                                             total_idx_rho_lin_spline, & ! In
                                             rho_lin_spline_vals, & ! In
                                             rho_lin_spline_levels, & ! In
                                             iv_other, p_sfc, & ! In
                                             grid_remap_method, & ! In
                                             l_zt_variable, & ! In
                                             l_adapt_grid )
    end do
    rtm_ref = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzt, & ! In
                                    rtm_ref, & ! In
                                    gr_target%nzt, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_pos_def, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    thlm_ref = remap_vals_to_target( ngrdcol, & ! In
                                     gr_source, gr_target, & ! In
                                     gr_source%nzt, & ! In
                                     thlm_ref, & ! In
                                     gr_target%nzt, & ! In
                                     total_idx_rho_lin_spline, & ! In
                                     rho_lin_spline_vals, & ! In
                                     rho_lin_spline_levels, & ! In
                                     iv_pos_def, p_sfc, & ! In
                                     grid_remap_method, & ! In
                                     l_zt_variable, & ! In
                                     l_adapt_grid ) ! Optional In
    um_ref = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzt, & ! In
                                   um_ref, & ! In
                                   gr_target%nzt, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_wind, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    vm_ref = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzt, & ! In
                                   vm_ref, & ! In
                                   gr_target%nzt, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_wind, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    ug = remap_vals_to_target( ngrdcol, & ! In
                               gr_source, gr_target, & ! In
                               gr_source%nzt, & ! In
                               ug, & ! In
                               gr_target%nzt, & ! In
                               total_idx_rho_lin_spline, & ! In
                               rho_lin_spline_vals, & ! In
                               rho_lin_spline_levels, & ! In
                               iv_wind, p_sfc, & ! In
                               grid_remap_method, & ! In
                               l_zt_variable, & ! In
                               l_adapt_grid ) ! Optional In
    vg = remap_vals_to_target( ngrdcol, & ! In
                               gr_source, gr_target, & ! In
                               gr_source%nzt, & ! In
                               vg, & ! In
                               gr_target%nzt, & ! In
                               total_idx_rho_lin_spline, & ! In
                               rho_lin_spline_vals, & ! In
                               rho_lin_spline_levels, & ! In
                               iv_wind, p_sfc, & ! In
                               grid_remap_method, & ! In
                               l_zt_variable, & ! In
                               l_adapt_grid ) ! Optional In
    um = remap_vals_to_target( ngrdcol, & ! In
                               gr_source, gr_target, & ! In
                               gr_source%nzt, & ! In
                               um, & ! In
                               gr_target%nzt, & ! In
                               total_idx_rho_lin_spline, & ! In
                               rho_lin_spline_vals, & ! In
                               rho_lin_spline_levels, & ! In
                               iv_wind, p_sfc, & ! In
                               grid_remap_method, & ! In
                               l_zt_variable, & ! In
                               l_adapt_grid ) ! Optional In
    vm = remap_vals_to_target( ngrdcol, & ! In
                               gr_source, gr_target, & ! In
                               gr_source%nzt, & ! In
                               vm, & ! In
                               gr_target%nzt, & ! In
                               total_idx_rho_lin_spline, & ! In
                               rho_lin_spline_vals, & ! In
                               rho_lin_spline_levels, & ! In
                               iv_wind, p_sfc, & ! In
                               grid_remap_method, & ! In
                               l_zt_variable, & ! In
                               l_adapt_grid ) ! Optional In
    up3 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                up3, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_other, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    vp3 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                vp3, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_other, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    rtm = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                rtm, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    thlm = remap_vals_to_target( ngrdcol, & ! In
                                 gr_source, gr_target, & ! In
                                 gr_source%nzt, & ! In
                                 thlm, & ! In
                                 gr_target%nzt, & ! In
                                 total_idx_rho_lin_spline, & ! In
                                 rho_lin_spline_vals, & ! In
                                 rho_lin_spline_levels, & ! In
                                 iv_pos_def, p_sfc, & ! In
                                 grid_remap_method, & ! In
                                 l_zt_variable, & ! In
                                 l_adapt_grid ) ! Optional In
    rtp3 = remap_vals_to_target( ngrdcol, & ! In
                                 gr_source, gr_target, & ! In
                                 gr_source%nzt, & ! In
                                 rtp3, & ! In
                                 gr_target%nzt, & ! In
                                 total_idx_rho_lin_spline, & ! In
                                 rho_lin_spline_vals, & ! In
                                 rho_lin_spline_levels, & ! In
                                 iv_other, p_sfc, & ! In
                                 grid_remap_method, & ! In
                                 l_zt_variable, & ! In
                                 l_adapt_grid ) ! Optional In
    thlp3 = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  thlp3, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    wp3 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                wp3, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_other, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    ! remap thvm values to new grid, to calculate new p_in_Pa
    thvm_new_grid = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          gr_source%nzt, & ! In
                                          thvm, & ! In
                                          gr_target%nzt, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv_pos_def, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable, & ! In
                                          l_adapt_grid ) ! Optional In
    ! calculate p_in_Pa instead of remapping directly since it can run into problems if for
    ! example the two highest levels have the same pressure value, which could be happening with
    ! the remapping, also calculate exner accordingly
    call init_pressure( ngrdcol, gr_target, thvm_new_grid, p_sfc, &
                        p_in_Pa_new, exner_new, p_in_Pa_zm, exner_zm )

    where ( spread(l_adapt_grid, 2, gr_target%nzt) )
      p_in_Pa = p_in_Pa_new
      exner = exner_new
    end where
    rcm = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzt, & ! In
                                rcm, & ! In
                                gr_target%nzt, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    cloud_frac = remap_vals_to_target( ngrdcol, & ! In
                                       gr_source, gr_target, & ! In
                                       gr_source%nzt, & ! In
                                       cloud_frac, & ! In
                                       gr_target%nzt, & ! In
                                       total_idx_rho_lin_spline, & ! In
                                       rho_lin_spline_vals, & ! In
                                       rho_lin_spline_levels, & ! In
                                       iv_pos_def, p_sfc, & ! In
                                       grid_remap_method, & ! In
                                       l_zt_variable, & ! In
                                       l_adapt_grid ) ! Optional In
    wp2thvp = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzt, & ! In
                                    wp2thvp, & ! In
                                    gr_target%nzt, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_other, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    wp2up = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  wp2up, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    wp2rtp = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzt, & ! In
                                   wp2rtp, & ! In
                                   gr_target%nzt, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_other, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    wp2thlp = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzt, & ! In
                                    wp2thlp, & ! In
                                    gr_target%nzt, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_other, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    wpup2 = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  wpup2, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    wpvp2 = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  wpvp2, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    ice_supersat_frac = remap_vals_to_target( ngrdcol, & ! In
                                              gr_source, gr_target, & ! In
                                              gr_source%nzt, & ! In
                                              ice_supersat_frac, & ! In
                                              gr_target%nzt, & ! In
                                              total_idx_rho_lin_spline, & ! In
                                              rho_lin_spline_vals, & ! In
                                              rho_lin_spline_levels, & ! In
                                              iv_pos_def, p_sfc, & ! In
                                              grid_remap_method, & ! In
                                              l_zt_variable, & ! In
                                              l_adapt_grid ) ! Optional In
    um_pert = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzt, & ! In
                                    um_pert, & ! In
                                    gr_target%nzt, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_wind, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    vm_pert = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzt, & ! In
                                    vm_pert, & ! In
                                    gr_target%nzt, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_wind, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    rcm_in_layer = remap_vals_to_target( ngrdcol, & ! In
                                         gr_source, gr_target, & ! In
                                         gr_source%nzt, & ! In
                                         rcm_in_layer, & ! In
                                         gr_target%nzt, & ! In
                                         total_idx_rho_lin_spline, & ! In
                                         rho_lin_spline_vals, & ! In
                                         rho_lin_spline_levels, & ! In
                                         iv_pos_def, p_sfc, & ! In
                                         grid_remap_method, & ! In
                                         l_zt_variable, & ! In
                                         l_adapt_grid ) ! Optional In
    cloud_cover = remap_vals_to_target( ngrdcol, & ! In
                                        gr_source, gr_target, & ! In
                                        gr_source%nzt, & ! In
                                        cloud_cover, & ! In
                                        gr_target%nzt, & ! In
                                        total_idx_rho_lin_spline, & ! In
                                        rho_lin_spline_vals, & ! In
                                        rho_lin_spline_levels, & ! In
                                        iv_pos_def, p_sfc, & ! In
                                        grid_remap_method, & ! In
                                        l_zt_variable, & ! In
                                        l_adapt_grid ) ! Optional In
    w_up_in_cloud = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          gr_source%nzt, & ! In
                                          w_up_in_cloud, & ! In
                                          gr_target%nzt, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv_other, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable, & ! In
                                          l_adapt_grid ) ! Optional In
    w_down_in_cloud = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzt, & ! In
                                            w_down_in_cloud, & ! In
                                            gr_target%nzt, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_other, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid ) ! Optional In
    cloudy_updraft_frac = remap_vals_to_target( ngrdcol, & ! In
                                                gr_source, gr_target, & ! In
                                                gr_source%nzt, & ! In
                                                cloudy_updraft_frac, & ! In
                                                gr_target%nzt, & ! In
                                                total_idx_rho_lin_spline, & ! In
                                                rho_lin_spline_vals, & ! In
                                                rho_lin_spline_levels, & ! In
                                                iv_pos_def, p_sfc, & ! In
                                                grid_remap_method, & ! In
                                                l_zt_variable, & ! In
                                                l_adapt_grid ) ! Optional In
    cloudy_downdraft_frac = remap_vals_to_target( ngrdcol, & ! In
                                                  gr_source, gr_target, & ! In
                                                  gr_source%nzt, & ! In
                                                  cloudy_downdraft_frac, & ! In
                                                  gr_target%nzt, & ! In
                                                  total_idx_rho_lin_spline, & ! In
                                                  rho_lin_spline_vals, & ! In
                                                  rho_lin_spline_levels, & ! In
                                                  iv_pos_def, p_sfc, & ! In
                                                  grid_remap_method, & ! In
                                                  l_zt_variable, & ! In
                                                  l_adapt_grid ) ! Optional In
    Kh_zt = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzt, & ! In
                                  Kh_zt, & ! In
                                  gr_target%nzt, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_pos_def, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    Lscale = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzt, & ! In
                                   Lscale, & ! In
                                   gr_target%nzt, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_pos_def, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In

    ! Remap all zm values
    l_zt_variable = .false.
    wprtp_forcing = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          gr_source%nzm, & ! In
                                          wprtp_forcing, & ! In
                                          gr_target%nzm, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv_other, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable, & ! In
                                          l_adapt_grid ) ! Optional In
    wpthlp_forcing = remap_vals_to_target( ngrdcol, & ! In
                                           gr_source, gr_target, & ! In
                                           gr_source%nzm, & ! In
                                           wpthlp_forcing, & ! In
                                           gr_target%nzm, & ! In
                                           total_idx_rho_lin_spline, & ! In
                                           rho_lin_spline_vals, & ! In
                                           rho_lin_spline_levels, & ! In
                                           iv_other, p_sfc, & ! In
                                           grid_remap_method, & ! In
                                           l_zt_variable, & ! In
                                           l_adapt_grid ) ! Optional In
    rtp2_forcing = remap_vals_to_target( ngrdcol, & ! In
                                         gr_source, gr_target, & ! In
                                         gr_source%nzm, & ! In
                                         rtp2_forcing, & ! In
                                         gr_target%nzm, & ! In
                                         total_idx_rho_lin_spline, & ! In
                                         rho_lin_spline_vals, & ! In
                                         rho_lin_spline_levels, & ! In
                                         iv_other, p_sfc, & ! In
                                         grid_remap_method, & ! In
                                         l_zt_variable, & ! In
                                         l_adapt_grid ) ! Optional In
    thlp2_forcing = remap_vals_to_target( ngrdcol, & ! In
                                          gr_source, gr_target, & ! In
                                          gr_source%nzm, & ! In
                                          thlp2_forcing, & ! In
                                          gr_target%nzm, & ! In
                                          total_idx_rho_lin_spline, & ! In
                                          rho_lin_spline_vals, & ! In
                                          rho_lin_spline_levels, & ! In
                                          iv_other, p_sfc, & ! In
                                          grid_remap_method, & ! In
                                          l_zt_variable, & ! In
                                          l_adapt_grid ) ! Optional In
    rtpthlp_forcing = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzm, & ! In
                                            rtpthlp_forcing, & ! In
                                            gr_target%nzm, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_other, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid ) ! Optional In
    wm_zm = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  wm_zm, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_wind, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    rho_zm = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzm, & ! In
                                   rho_zm, & ! In
                                   gr_target%nzm, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_pos_def, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    rho_ds_zm = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzm, & ! In
                                      rho_ds_zm, & ! In
                                      gr_target%nzm, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_pos_def, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    invrs_rho_ds_zm = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzm, & ! In
                                            invrs_rho_ds_zm, & ! In
                                            gr_target%nzm, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_pos_def, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid ) ! Optional In
    thv_ds_zm = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzm, & ! In
                                      thv_ds_zm, & ! In
                                      gr_target%nzm, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_pos_def, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    do i = 1, hydromet_dim
      wphydrometp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                                 gr_source, gr_target, & ! In
                                                 gr_source%nzm, & ! In
                                                 wphydrometp(:,:,i), & ! In
                                                 gr_target%nzm, & ! In
                                                 total_idx_rho_lin_spline, & ! In
                                                 rho_lin_spline_vals, & ! In
                                                 rho_lin_spline_levels, & ! In
                                                 iv_other, p_sfc, & ! In
                                                 grid_remap_method, & ! In
                                                 l_zt_variable, & ! In
                                                 l_adapt_grid )
    end do
    upwp = remap_vals_to_target( ngrdcol, & ! In
                                 gr_source, gr_target, & ! In
                                 gr_source%nzm, & ! In
                                 upwp, & ! In
                                 gr_target%nzm, & ! In
                                 total_idx_rho_lin_spline, & ! In
                                 rho_lin_spline_vals, & ! In
                                 rho_lin_spline_levels, & ! In
                                 iv_other, p_sfc, & ! In
                                 grid_remap_method, & ! In
                                 l_zt_variable, & ! In
                                 l_adapt_grid ) ! Optional In
    vpwp = remap_vals_to_target( ngrdcol, & ! In
                                 gr_source, gr_target, & ! In
                                 gr_source%nzm, & ! In
                                 vpwp, & ! In
                                 gr_target%nzm, & ! In
                                 total_idx_rho_lin_spline, & ! In
                                 rho_lin_spline_vals, & ! In
                                 rho_lin_spline_levels, & ! In
                                 iv_other, p_sfc, & ! In
                                 grid_remap_method, & ! In
                                 l_zt_variable, & ! In
                                 l_adapt_grid ) ! Optional In
    up2 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzm, & ! In
                                up2, & ! In
                                gr_target%nzm, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    vp2 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzm, & ! In
                                vp2, & ! In
                                gr_target%nzm, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    wprtp = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  wprtp, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    wpthlp = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzm, & ! In
                                   wpthlp, & ! In
                                   gr_target%nzm, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_other, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    rtp2 = remap_vals_to_target( ngrdcol, & ! In
                                 gr_source, gr_target, & ! In
                                 gr_source%nzm, & ! In
                                 rtp2, & ! In
                                 gr_target%nzm, & ! In
                                 total_idx_rho_lin_spline, & ! In
                                 rho_lin_spline_vals, & ! In
                                 rho_lin_spline_levels, & ! In
                                 iv_pos_def, p_sfc, & ! In
                                 grid_remap_method, & ! In
                                 l_zt_variable, & ! In
                                 l_adapt_grid ) ! Optional In
    thlp2 = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  thlp2, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_pos_def, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    rtpthlp = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzm, & ! In
                                    rtpthlp, & ! In
                                    gr_target%nzm, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_other, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    wp2 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzm, & ! In
                                wp2, & ! In
                                gr_target%nzm, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    do i = 1, sclr_dim
      wpsclrp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                             gr_source, gr_target, & ! In
                                             gr_source%nzm, & ! In
                                             wpsclrp(:,:,i), & ! In
                                             gr_target%nzm, & ! In
                                             total_idx_rho_lin_spline, & ! In
                                             rho_lin_spline_vals, & ! In
                                             rho_lin_spline_levels, & ! In
                                             iv_other, p_sfc, & ! In
                                             grid_remap_method, & ! In
                                             l_zt_variable, & ! In
                                             l_adapt_grid )
      sclrp2(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                            gr_source, gr_target, & ! In
                                            gr_source%nzm, & ! In
                                            sclrp2(:,:,i), & ! In
                                            gr_target%nzm, & ! In
                                            total_idx_rho_lin_spline, & ! In
                                            rho_lin_spline_vals, & ! In
                                            rho_lin_spline_levels, & ! In
                                            iv_other, p_sfc, & ! In
                                            grid_remap_method, & ! In
                                            l_zt_variable, & ! In
                                            l_adapt_grid )
      sclrprtp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                              gr_source, gr_target, & ! In
                                              gr_source%nzm, & ! In
                                              sclrprtp(:,:,i), & ! In
                                              gr_target%nzm, & ! In
                                              total_idx_rho_lin_spline, & ! In
                                              rho_lin_spline_vals, & ! In
                                              rho_lin_spline_levels, & ! In
                                              iv_other, p_sfc, & ! In
                                              grid_remap_method, & ! In
                                              l_zt_variable, & ! In
                                              l_adapt_grid )
      sclrpthlp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                               gr_source, gr_target, & ! In
                                               gr_source%nzm, & ! In
                                               sclrpthlp(:,:,i), & ! In
                                               gr_target%nzm, & ! In
                                               total_idx_rho_lin_spline, & ! In
                                               rho_lin_spline_vals, & ! In
                                               rho_lin_spline_levels, & ! In
                                               iv_other, p_sfc, & ! In
                                               grid_remap_method, & ! In
                                               l_zt_variable, & ! In
                                               l_adapt_grid )
      sclrpthvp(:,:,i) = remap_vals_to_target( ngrdcol, & ! In
                                               gr_source, gr_target, & ! In
                                               gr_source%nzm, & ! In
                                               sclrpthvp(:,:,i), & ! In
                                               gr_target%nzm, & ! In
                                               total_idx_rho_lin_spline, & ! In
                                               rho_lin_spline_vals, & ! In
                                               rho_lin_spline_levels, & ! In
                                               iv_other, p_sfc, & ! In
                                               grid_remap_method, & ! In
                                               l_zt_variable, & ! In
                                               l_adapt_grid )
    end do
    wpthvp = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzm, & ! In
                                   wpthvp, & ! In
                                   gr_target%nzm, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_other, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    rtpthvp = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzm, & ! In
                                    rtpthvp, & ! In
                                    gr_target%nzm, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_other, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In
    thlpthvp = remap_vals_to_target( ngrdcol, & ! In
                                     gr_source, gr_target, & ! In
                                     gr_source%nzm, & ! In
                                     thlpthvp, & ! In
                                     gr_target%nzm, & ! In
                                     total_idx_rho_lin_spline, & ! In
                                     rho_lin_spline_vals, & ! In
                                     rho_lin_spline_levels, & ! In
                                     iv_other, p_sfc, & ! In
                                     grid_remap_method, & ! In
                                     l_zt_variable, & ! In
                                     l_adapt_grid ) ! Optional In
    uprcp = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  uprcp, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    vprcp = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  vprcp, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    rc_coef_zm = remap_vals_to_target( ngrdcol, & ! In
                                       gr_source, gr_target, & ! In
                                       gr_source%nzm, & ! In
                                       rc_coef_zm, & ! In
                                       gr_target%nzm, & ! In
                                       total_idx_rho_lin_spline, & ! In
                                       rho_lin_spline_vals, & ! In
                                       rho_lin_spline_levels, & ! In
                                       iv_pos_def, p_sfc, & ! In
                                       grid_remap_method, & ! In
                                       l_zt_variable, & ! In
                                       l_adapt_grid ) ! Optional In
    wp4 = remap_vals_to_target( ngrdcol, & ! In
                                gr_source, gr_target, & ! In
                                gr_source%nzm, & ! In
                                wp4, & ! In
                                gr_target%nzm, & ! In
                                total_idx_rho_lin_spline, & ! In
                                rho_lin_spline_vals, & ! In
                                rho_lin_spline_levels, & ! In
                                iv_pos_def, p_sfc, & ! In
                                grid_remap_method, & ! In
                                l_zt_variable, & ! In
                                l_adapt_grid ) ! Optional In
    wp2up2 = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzm, & ! In
                                   wp2up2, & ! In
                                   gr_target%nzm, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_pos_def, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    wp2vp2 = remap_vals_to_target( ngrdcol, & ! In
                                   gr_source, gr_target, & ! In
                                   gr_source%nzm, & ! In
                                   wp2vp2, & ! In
                                   gr_target%nzm, & ! In
                                   total_idx_rho_lin_spline, & ! In
                                   rho_lin_spline_vals, & ! In
                                   rho_lin_spline_levels, & ! In
                                   iv_pos_def, p_sfc, & ! In
                                   grid_remap_method, & ! In
                                   l_zt_variable, & ! In
                                   l_adapt_grid ) ! Optional In
    upwp_pert = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzm, & ! In
                                      upwp_pert, & ! In
                                      gr_target%nzm, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_other, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    vpwp_pert = remap_vals_to_target( ngrdcol, & ! In
                                      gr_source, gr_target, & ! In
                                      gr_source%nzm, & ! In
                                      vpwp_pert, & ! In
                                      gr_target%nzm, & ! In
                                      total_idx_rho_lin_spline, & ! In
                                      rho_lin_spline_vals, & ! In
                                      rho_lin_spline_levels, & ! In
                                      iv_other, p_sfc, & ! In
                                      grid_remap_method, & ! In
                                      l_zt_variable, & ! In
                                      l_adapt_grid ) ! Optional In
    wprcp = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  wprcp, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_other, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    invrs_tau_zm = remap_vals_to_target( ngrdcol, & ! In
                                         gr_source, gr_target, & ! In
                                         gr_source%nzm, & ! In
                                         invrs_tau_zm, & ! In
                                         gr_target%nzm, & ! In
                                         total_idx_rho_lin_spline, & ! In
                                         rho_lin_spline_vals, & ! In
                                         rho_lin_spline_levels, & ! In
                                         iv_pos_def, p_sfc, & ! In
                                         grid_remap_method, & ! In
                                         l_zt_variable, & ! In
                                         l_adapt_grid ) ! Optional In
    Kh_zm = remap_vals_to_target( ngrdcol, & ! In
                                  gr_source, gr_target, & ! In
                                  gr_source%nzm, & ! In
                                  Kh_zm, & ! In
                                  gr_target%nzm, & ! In
                                  total_idx_rho_lin_spline, & ! In
                                  rho_lin_spline_vals, & ! In
                                  rho_lin_spline_levels, & ! In
                                  iv_pos_def, p_sfc, & ! In
                                  grid_remap_method, & ! In
                                  l_zt_variable, & ! In
                                  l_adapt_grid ) ! Optional In
    thlprcp = remap_vals_to_target( ngrdcol, & ! In
                                    gr_source, gr_target, & ! In
                                    gr_source%nzm, & ! In
                                    thlprcp, & ! In
                                    gr_target%nzm, & ! In
                                    total_idx_rho_lin_spline, & ! In
                                    rho_lin_spline_vals, & ! In
                                    rho_lin_spline_levels, & ! In
                                    iv_other, p_sfc, & ! In
                                    grid_remap_method, & ! In
                                    l_zt_variable, & ! In
                                    l_adapt_grid ) ! Optional In

  end subroutine remap_all_clubb_core_vals

!===============================================================================

end module grid_adaptation_module
