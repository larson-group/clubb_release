!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_adaptation_module


  use clubb_precision, only: &
      core_rknd ! Variable(s)

  use grid_class, only: &
      grid ! Type

  use constants_clubb, only: &
      one, fstderr ! Constants

  use error_code, only: &
      clubb_at_least_debug_level

  implicit none

  public :: setup_gr_dycore, &
            normalize_grid_density, &
            adapt_grid, &
            calc_grid_dens, &
            clean_up_grid_adaptation_module

  !private :: check_remap_conservation, check_remap_consistency, check_remap_monotonicity, &
  !           remapping_matrix, matrix_vector_mult, &
  !           check_remap_consistency_w_vals, &
  !           check_mass_conservation

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

    character(len=100) :: & 
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
    !nzmax = 31 ! since dycore grid has in total 73 (zm) levels
    grid_type = 3
    !zm_grid_fname = '../input/grid/dycore.grd'
    !zm_grid_fname = '../input/grid/31_level_zm_grid.grd'
    zm_grid_fname = '../input/grid/test_dycore.grd'
    !zm_grid_fname = '../input/grid/dycore_e3sm_coarsened_6_0.grd'
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

  subroutine setup_min_gr( iunit, ngrdcol, grid_sfc, grid_top, gr )

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

    character(len=100) :: & 
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
    !nzmax = 31 ! since dycore grid has in total 73 (zm) levels
    grid_type = 3
    !zm_grid_fname = '../input/grid/dycore.grd'
    !zm_grid_fname = '../input/grid/31_level_zm_grid.grd'
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
    !write(*,*) 'nlevel: ', nlevel
    !write(*,*) 'momentum_heights: ', momentum_heights
    !write(*,*) 'grid_sfc: ', grid_sfc
    !write(*,*) 'grid_top: ', grid_top
    call setup_grid( nlevel, ngrdcol, sfc_elevation, l_implemented, &   ! intent(in)
                     grid_type, deltaz, grid_sfc, grid_top, &           ! intent(in)
                     momentum_heights, thermodynamic_heights, &         ! intent(in)
                     gr, err_code )                                     ! intent(inout)

    if ( err_code == clubb_fatal_error ) then
      error stop "Error in CLUBB calling setup_grid"
    end if

  end subroutine setup_min_gr

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
      !call setup_gr_dycore( iunit, ngrdcol, &              ! Intent(in)
      !                      grid_sfc_arr, grid_top_arr, &  ! Intent(in)
      !                      gr_dycore )                    ! Intent(out)
      call setup_min_gr( iunit, ngrdcol, &              ! Intent(in)
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

    ! TODO maybe just use invrs_dzm instead since it has already nzm levels
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
      
      
      !write(*,*) 'difference: ', sum/j
      !write(*,*) 'threshold: ', threshold
      if ( sum/j > threshold ) then
        l_adapt_grid_tmp = .true.
      else
        !write(*,*) '-----------------not adapted'
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

    !write(iunit, *) 'gr_dens_z', itime, gr_dens_norm_z
    !write(iunit, *) 'gr_dens', itime, norm_grid_dens
    !write(iunit, *) 'min_gr_dens_z', itime, min_gr_dens_norm_z
    !write(iunit, *) 'min_gr_dens', itime, norm_min_grid_dens

    !call sleep(5)

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

    !write(iunit, *) 'gr_dens_z', itime, gr_dens_norm_z
    !write(iunit, *) 'gr_dens', itime, gr_dens_norm
    !write(iunit, *) 'min_gr_dens_z', itime, min_gr_dens_norm_z
    !write(iunit, *) 'min_gr_dens', itime, min_gr_dens_norm

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
                         time_remap, &
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
      cons_ullrich_remap, &
      ppm_remap

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
    real( kind = core_rknd ), intent(inout) :: &
      time_remap
    
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

    real( kind = core_rknd ) :: &
      time_stop, time_start    ! help variables to measure the time [s]

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

      call cpu_time(time_start)

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

      call cpu_time(time_stop)
      time_remap = time_remap + time_stop - time_start

      
      gr = new_gr
    end if

  end subroutine adapt_grid

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
    
    use grid_class, only: &
        setup_grid

    use calc_pressure, only: &
        init_pressure    ! Procedure(s)

    use model_flags, only: &
        cons_ullrich_remap, &
        ppm_remap

    use error_code, only: &
        clubb_fatal_error

    use remapping_module, only: &
        remap_vals_to_target

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

    logical :: l_zt_variable

    !--------------------- Begin Code ---------------------

        ! Remap all zt values
        l_zt_variable = .true.
        thlm_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thlm_forcing, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
         rtm_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rtm_forcing, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         um_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    um_forcing, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         vm_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vm_forcing, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wm_zt = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wm_zt, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         rho = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rho, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rho_ds_zt, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         invrs_rho_ds_zt = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    invrs_rho_ds_zt, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         thv_ds_zt = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thv_ds_zt, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
         rfrzm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rfrzm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        do i = 1, hydromet_dim
            
             wp2hmp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wp2hmp(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
             rtphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rtphmp_zt(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
             thlphmp_zt(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thlphmp_zt(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        end do
        do i = 1, sclr_dim
            
             sclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    sclrm_forcing(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
             sclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    sclrm(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
             sclrp3(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    sclrp3(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        end do
        do i = 1, edsclr_dim
            
             edsclrm_forcing(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    edsclrm_forcing(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
             edsclrm(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    edsclrm(:,:,i), &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        end do
        
         rtm_ref = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rtm_ref, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
         thlm_ref = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thlm_ref, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         um_ref = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    um_ref, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         vm_ref = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vm_ref, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         ug = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    ug, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         vg = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vg, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         um = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    um, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         vm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         up3 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    up3, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         vp3 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vp3, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         rtm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rtm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         thlm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thlm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         rtp3 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rtp3, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         thlp3 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thlp3, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wp3 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wp3, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        ! remap thvm values to new grid, to calculate new p_in_Pa
         thvm_new_grid = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    thvm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        ! calculate p_in_Pa instead of remapping directly since it can run into problems if for
        ! example the two highest levels have the same pressure value, which could be happening with
        ! the remapping, also calculate exner accordingly
        call init_pressure( ngrdcol, new_gr, thvm_new_grid, p_sfc, &
                            p_in_Pa, exner, p_in_Pa_zm, exner_zm )
        
         rcm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rcm, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         cloud_frac = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    cloud_frac, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wp2thvp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wp2thvp, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wp2rtp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wp2rtp, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wp2thlp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wp2thlp, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wpup2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wpup2, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         wpvp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    wpvp2, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

         ice_supersat_frac = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    ice_supersat_frac, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

         um_pert = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    um_pert, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

         vm_pert = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    vm_pert, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         rcm_in_layer = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    rcm_in_layer, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         cloud_cover = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    cloud_cover, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         w_up_in_cloud = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    w_up_in_cloud, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         w_down_in_cloud = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    w_down_in_cloud, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         cloudy_updraft_frac = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    cloudy_updraft_frac, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         cloudy_downdraft_frac = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    cloudy_downdraft_frac, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
         Kh_zt = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    Kh_zt, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

         Lscale = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzt, &
                                    Lscale, &
                                    new_gr%nzt, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )


        ! Remap all zm values
        l_zt_variable = .false.
        wprtp_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wprtp_forcing, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wpthlp_forcing, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rtp2_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rtp2_forcing, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        thlp2_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    thlp2_forcing, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rtpthlp_forcing = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rtpthlp_forcing, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wm_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wm_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_wind, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rho_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rho_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rho_ds_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        invrs_rho_ds_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    invrs_rho_ds_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        thv_ds_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    thv_ds_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        do i = 1, hydromet_dim
            
            wphydrometp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wphydrometp(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        end do
        
        upwp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    upwp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        vpwp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    vpwp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        up2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    up2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        vp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    vp2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wprtp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wprtp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wpthlp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wpthlp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rtp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rtp2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        thlp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    thlp2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rtpthlp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rtpthlp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wp2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        do i = 1, sclr_dim
            
            wpsclrp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wpsclrp(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
            sclrp2(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    sclrp2(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
            sclrprtp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    sclrprtp(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
            sclrpthlp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    sclrpthlp(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
            
            sclrpthvp(:,:,i) = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    sclrpthvp(:,:,i), &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        end do
        
        wpthvp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wpthvp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rtpthvp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rtpthvp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        thlpthvp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    thlpthvp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        uprcp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    uprcp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        vprcp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    vprcp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        rc_coef_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    rc_coef_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wp4 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wp4, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wp2up2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wp2up2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wp2vp2 = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wp2vp2, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        upwp_pert = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    upwp_pert, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        vpwp_pert = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    vpwp_pert, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        wprcp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    wprcp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        invrs_tau_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    invrs_tau_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        Kh_zm = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    Kh_zm, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_pos_def, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )
        
        thlprcp = remap_vals_to_target( ngrdcol, &
                                    gr, new_gr, &
                                    gr%nzm, &
                                    thlprcp, &
                                    new_gr%nzm, &
                                    total_idx_rho_lin_spline, &
                                    rho_lin_spline_vals, &
                                    rho_lin_spline_levels, &
                                    iv_other, p_sfc, &
                                    grid_remap_method, &
                                    l_zt_variable )

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


!===============================================================================

end module grid_adaptation_module