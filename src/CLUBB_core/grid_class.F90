!------------------------------------------------------------------------
! $Id$
!=============================================================================== 
module grid_class

  ! Description:
  !
  ! Definition of a grid class and associated functions
  !
  ! The grid specification is as follows:
  !
  !     +                ================== zm(nz) =========GP=========
  !     |
  !     |
  ! 1/dzt(nz) +          ------------------ zt(nz) ---------GP---------
  !     |        |
  !     |        |
  !     + 1/dzm(nz-1)    ================== zm(nz-1) ==================
  !              |
  !              |
  !              +       ------------------ zt(nz-1) ------------------
  !
  !                                           .
  !                                           .
  !                                           .
  !                                           .
  !
  !                      ================== zm(k+1) ===================
  !
  !
  !              +       ------------------ zt(k+1) -------------------
  !              |
  !              |
  !     +    1/dzm(k)    ================== zm(k) =====================
  !     |        |
  !     |        |
  ! 1/dzt(k)     +       ------------------ zt(k) ---------------------
  !     |
  !     |
  !     +                ================== zm(k-1) ===================
  !
  !
  !                      ------------------ zt(k-1) -------------------
  !
  !                                           .
  !                                           .
  !                                           .
  !                                           .
  !
  !     +                ================== zm(2) =====================
  !     |
  !     |
  ! 1/dzt(2)     +       ------------------ zt(2) ---------------------
  !     |        |
  !     |        |
  !     +    1/dzm(1)    ================== zm(1) ============GP=======  zm_init
  !              |       //////////////////////////////////////////////  surface
  !              |
  !              +       ------------------ zt(1) ------------GP-------
  !
  !
  ! The variable zm(k) stands for the momentum level altitude at momentum
  ! level k; the variable zt(k) stands for the thermodynamic level altitude at
  ! thermodynamic level k; the variable invrs_dzt(k) is the inverse distance
  ! between momentum levels (over a central thermodynamic level k); and the
  ! variable invrs_dzm(k) is the inverse distance between thermodynamic levels
  ! (over a central momentum level k).  Please note that in the above diagram,
  ! "invrs_dzt" is denoted "dzt", and "invrs_dzm" is denoted "dzm", such that
  ! 1/dzt is the distance between successive momentum levels k-1 and k (over a
  ! central thermodynamic level k), and 1/dzm is the distance between successive
  ! thermodynamic levels k and k+1 (over a central momentum level k).
  !
  ! The grid setup is compatible with a stretched (unevely-spaced) grid.  Thus,
  ! the distance between successive grid levels may not always be constant.
  !
  ! The following diagram is an example of a stretched grid that is defined on
  ! momentum levels.  The thermodynamic levels are placed exactly halfway
  ! between the momentum levels.  However, the momentum levels do not fall
  ! halfway between the thermodynamic levels.
  !
  !        =============== zm(k+1) ===============
  !
  !
  !
  !        --------------- zt(k+1) ---------------
  !
  !
  !
  !        ===============  zm(k)  ===============
  !
  !        ---------------  zt(k)  ---------------
  !
  !        =============== zm(k-1) ===============
  !
  ! The following diagram is an example of a stretched grid that is defined on
  ! thermodynamic levels.  The momentum levels are placed exactly halfway
  ! between the thermodynamic levels.  However, the thermodynamic levels do not
  ! fall halfway between the momentum levels.
  !
  !        --------------- zt(k+1) ---------------
  !
  !
  !
  !        ===============  zm(k)  ===============
  !
  !
  !
  !        ---------------  zt(k)  ---------------
  !
  !        =============== zm(k-1) ===============
  !
  !        --------------- zt(k-1) ---------------
  !
  ! NOTE:  Any future code written for use in the CLUBB parameterization should
  !        use interpolation formulas consistent with a stretched grid.  The
  !        simplest way to do so is to call the appropriate interpolation
  !        function from this module.  Interpolations should *not* be handled in
  !        the form of:  ( var_zm(k) + var_zm(k-1) ) / 2; *nor* in the form of:
  !        0.5_core_rknd*( var_zt(k+1) + var_zt(k) ).  Rather, all explicit interpolations
  !        should call zt2zm or zm2zt; while interpolations for a variable being
  !        solved for implicitly in the code should use gr%weights_zt2zm (which
  !        refers to interp_weights_zt2zm_imp), or gr%weights_zm2zt (which
  !        refers to interp_weights_zm2zt_imp).
  !
  ! Momentum level 1 is placed at altitude zm_init, which is usually at the
  ! surface.  However, in general, zm_init can be at any altitude defined by the
  ! user.
  !
  ! GP indicates ghost points. Variables located at those levels are not
  ! prognosed, but only used for boundary conditions.
  !
  ! Chris Golaz, 7/17/99
  ! modified 9/10/99
  ! schemena, modified 6/11/2014 - Restructered code to add cubic/linear flag

  !  References:
  !
  !  https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_grid
  !
  !  Section 3c, p. 3548 /Numerical discretization/ of:
  !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  !     Method and Model Description'' Golaz, et al. (2002)
  !     JAS, Vol. 59, pp. 3540--3551.
  !-----------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: grid, zt2zm, zm2zt, zt2zm2zt, zm2zt2zm, & 
            ddzm, ddzt, & 
            setup_grid, cleanup_grid, setup_grid_heights, &
            read_grid_heights, flip

  private :: t_above, t_below, m_above, m_below

  private ! Default Scoping

  ! Constant parameters
  integer, parameter :: & 
    t_above = 1, & ! Upper thermodynamic level index (gr%weights_zt2zm).
    t_below = 2, & ! Lower thermodynamic level index (gr%weights_zt2zm).
    m_above = 1, & ! Upper momentum level index (gr%weights_zm2zt).
    m_below = 2    ! Lower momentum level index (gr%weights_zm2zt).
  
  type grid

    integer :: nz ! Number of points in the grid
    !   Note: Fortran 90/95 prevents an allocatable array from appearing
    !   within a derived type.  However, Fortran 2003 does not!!!!!!!!
    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      zm, & ! Momentum grid
      zt    ! Thermo grid
      
    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      invrs_dzm, & ! The inverse spacing between thermodynamic grid
                   ! levels; centered over momentum grid levels.
      invrs_dzt    ! The inverse spacing between momentum grid levels;
                   ! centered over thermodynamic grid levels.

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      dzm, &  ! Spacing between thermodynamic grid levels; centered over
              ! momentum grid levels
      dzt     ! Spcaing between momentum grid levels; centered over
              ! thermodynamic grid levels

    ! These weights are normally used in situations
    ! where a momentum level variable is being
    ! solved for implicitly in an equation and
    ! needs to be interpolated to the thermodynamic grid levels.
    real( kind = core_rknd ), allocatable, dimension(:,:,:) :: weights_zm2zt, & 
    ! These weights are normally used in situations where a
    ! thermodynamic level variable is being solved for implicitly in an equation
    ! and needs to be interpolated to the momentum grid levels.
                                     weights_zt2zm

  end type grid

  !   The grid is defined here so that it is common throughout the module.
  !   The implication is that only one grid can be defined !


!   Modification for using CLUBB in a host model (i.e. one grid per column)


  ! Interfaces provided for function overloading

  interface zt2zm
    ! For l_cubic_interp = .true.
    ! This version uses cublic spline interpolation of Stefen (1990).
    !
    ! For l_cubic_interp = .false.
    ! This performs a linear extension at the highest grid level and therefore
    ! does not guarantee, for positive definite quantities (e.g. wp2), that the
    ! extended point is indeed positive definite.  Positive definiteness can be
    ! ensured with a max statement.
    ! In the future, we could add a flag (lposdef) and, when needed, apply the
    ! max statement directly within interpolated_azm and interpolated_azmk.
    module procedure redirect_interpolated_azm_k    ! Works over a single vertical level
    module procedure redirect_interpolated_azm_1D   ! Works over all vertical levels 
    module procedure redirect_interpolated_azm_2D   ! Works over all vertical levels and columns
  end interface

  interface zm2zt
    ! For l_cubic_interp = .true.
    ! This version uses cublic spline interpolation of Stefen (1990).
    !
    ! For l_cubic_interp = .false.
    ! This performs a linear extension at the lowest grid level and therefore
    ! does not guarantee, for positive definite quantities (e.g. wp2), that the
    ! extended point is indeed positive definite.  Positive definiteness can be
    ! ensured with a max statement.
    ! In the future, we could add a flag (lposdef) and, when needed, apply the
    ! max statement directly within interpolated_azt and interpolated_aztk.
    module procedure redirect_interpolated_azt_k    ! Works over a single vertical level
    module procedure redirect_interpolated_azt_1D   ! Works over all vertical levels 
    module procedure redirect_interpolated_azt_2D   ! Works over all vertical levels and columns
  end interface

  ! Vertical derivative functions
  interface ddzm
    module procedure gradzm_1D    ! Works over all vertical levels 
    module procedure gradzm_2D    ! Works over all vertical levels and columns
  end interface

  interface ddzt
    module procedure gradzt_1D    ! Works over all vertical levels 
    module procedure gradzt_2D    ! Works over all vertical levels and columns
  end interface

  contains

  !=============================================================================
  subroutine setup_grid( nzmax, ngrdcol, sfc_elevation, l_implemented,  &
                         grid_type, deltaz, zm_init, zm_top,            &
                         momentum_heights, thermodynamic_heights,       &
                         gr )

    ! Description:
    !   Grid Constructor
    !
    !   This subroutine sets up the CLUBB vertical grid.
    !
    ! References:
    !   ``Equations for CLUBB'',  Sec. 8,  Grid Configuration.
    !-----------------------------------------------------------------------

    use constants_clubb, only:  & 
        fstderr ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        err_code, &                     ! Error indicator
        clubb_fatal_error, &            ! Constant
        err_header                      ! String

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      NWARNING = 250 ! Issue a warning if nzmax exceeds this number.

    ! Input Variables
    integer, intent(in) ::  & 
      nzmax, &  ! Number of vertical levels in grid      [#]
      ngrdcol

    type(grid), target, intent(inout) :: gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) ::  & 
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    ! Local Variables
    integer :: &
      begin_height, &  ! Lower bound for *_heights arrays [-]
      end_height       ! Upper bound for *_heights arrays [-]

    integer :: ierr, & ! Allocation stat
               i, k    ! Loop index


    ! ---- Begin Code ----

    ! Define the grid size

    if ( nzmax > NWARNING .and. clubb_at_least_debug_level( 1 ) ) then
      write(fstderr,*) "Warning:  running with vertical grid "// & 
                       "which is larger than", NWARNING, "levels."
      write(fstderr,*) "This may take a lot of CPU time and memory."
    end if

    gr%nz = nzmax

    ! Default bounds
    begin_height = 1
    end_height   = gr%nz

    !---------------------------------------------------
    if ( .not. l_implemented ) then
      
      ! Since l_implemented=.false., we have to calculate the number of vertical levels
      ! depending on the grid type. Currently, clubb can only have multiple columns
      ! if implemented in a host model, but if in the future we have a standalone 
      ! version with multiple columns, we may end up with columns that have a
      ! different number of vertical levels. CLUBB cannot handle this however, so
      ! we use only the first column to determine the number of vertical levels.

      if ( grid_type == 1 ) then

        ! Determine the number of grid points given the spacing
        ! to fit within the bounds without going over.
        gr%nz = floor( ( zm_top(1) - zm_init(1) + deltaz(1) ) / deltaz(1) )

      else if( grid_type == 2 ) then! Thermo

        ! Find begin_height (lower bound)
        k = gr%nz

        do while( thermodynamic_heights(1,k) >= zm_init(1) .and. k > 1 )

          k = k - 1

        end do

        if( thermodynamic_heights(1,k) >= zm_init(1) ) then

          write(fstderr,*) err_header, "Stretched zt grid cannot fulfill zm_init requirement"
          err_code = clubb_fatal_error
          return

        else

          begin_height = k

        end if

        ! Find end_height (upper bound)
        k = gr%nz

        do while( thermodynamic_heights(1,k) > zm_top(1) .and. k > 1 )

          k = k - 1

        end do

        if( zm_top(1) < thermodynamic_heights(1,k) ) then

          write(fstderr,*) err_header, "Stretched zt grid cannot fulfill zm_top requirement"
          err_code = clubb_fatal_error
          return

        else

          end_height = k

          gr%nz = size( thermodynamic_heights(1,begin_height:end_height ) )

        end if

      else if( grid_type == 3 ) then ! Momentum

        ! Find begin_height (lower bound)
        k = 1

        do while( momentum_heights(1,k) < zm_init(1) .and. k < gr%nz )

          k = k + 1

        end do

        if( momentum_heights(1,k) < zm_init(1) ) then

          write(fstderr,*) err_header, "Stretched zm grid cannot fulfill zm_init requirement"
          err_code = clubb_fatal_error
          return

        else

          begin_height = k

        end if

        ! Find end_height (lower bound)
        k = gr%nz

        do while( momentum_heights(1,k) > zm_top(1) .and. k > 1 )

          k = k - 1

        end do

        if( momentum_heights(1,k) > zm_top(1) ) then

          write(fstderr,*) err_header, "Stretched zm grid cannot fulfill zm_top requirement"
          err_code = clubb_fatal_error
          return

        else

          end_height = k

          gr%nz = size( momentum_heights(1,begin_height:end_height ) )

        end if

      endif ! grid_type

    endif ! .not. l_implemented

    !---------------------------------------------------

    ! Allocate memory for the grid levels
    allocate( gr%zm(ngrdcol,gr%nz), gr%zt(ngrdcol,gr%nz), & 
              gr%dzm(ngrdcol,gr%nz), gr%dzt(ngrdcol,gr%nz), &
              gr%invrs_dzm(ngrdcol,gr%nz), gr%invrs_dzt(ngrdcol,gr%nz),  & 
              gr%weights_zm2zt(ngrdcol,gr%nz,m_above:m_below), & 
              gr%weights_zt2zm(ngrdcol,gr%nz,t_above:t_below), & 
              stat=ierr )

    if ( ierr /= 0 ) then
      write(fstderr,*) err_header, "In setup_grid: allocation of grid variables failed."
      err_code = clubb_fatal_error
      return
    end if

    ! Set the values for the derived types used for heights, derivatives, and
    ! interpolation from the momentum/thermodynamic grid
    call setup_grid_heights( &
                  gr%nz, ngrdcol, & ! intent(in)
                  l_implemented, grid_type,  & ! intent(in)
                  deltaz, zm_init,  & ! intent(in)
                  momentum_heights(:,begin_height:end_height),  & ! intent(in) 
                  thermodynamic_heights(:,begin_height:end_height), & ! intent(in)
                  gr ) ! intent(inout)

    do i = 1, ngrdcol
      if ( sfc_elevation(i) > gr%zm(i,1) ) then
        write(fstderr,*) "The altitude of the lowest momentum level, "        &
                         // "gr%zm(1,1), must be at or above the altitude of "  &
                         // "the surface, sfc_elevation.  The lowest model "  &
                         // "momentum level cannot be below the surface."
        write(fstderr,*) "Altitude of lowest momentum level =", gr%zm(i,1)
        write(fstderr,*) "Altitude of the surface =", sfc_elevation(i)
        err_code = clubb_fatal_error
        return
      endif
    end do

    return

  end subroutine setup_grid

  !=============================================================================
  subroutine cleanup_grid( gr )

    ! Description:
    !   De-allocates the memory for the grid
    !
    ! References:
    !   None
    !------------------------------------------------------------------------------
    use constants_clubb, only: &
        fstderr  ! Constant(s)

    implicit none
  
    type(grid), target, intent(inout) :: gr

    ! Local Variable(s)
    integer :: ierr

    ! ----- Begin Code -----

    ! Allocate memory for grid levels
    deallocate( gr%zm, gr%zt, & 
                gr%dzm, gr%dzt, &
                gr%invrs_dzm, gr%invrs_dzt,  & 
                gr%weights_zm2zt, gr%weights_zt2zm, & 
                stat=ierr )

    if ( ierr /= 0 ) then
      write(fstderr,*) "Grid deallocation failed."
    end if


    return

  end subroutine cleanup_grid

  !=============================================================================
  subroutine setup_grid_heights( nz, ngrdcol, &
                                 l_implemented, grid_type,  & 
                                 deltaz, zm_init, momentum_heights,  & 
                                 thermodynamic_heights, &
                                 gr )

    ! Description:
    !   Sets the heights and interpolation weights of the column.
    !   This is seperated from setup_grid for those host models that have heights
    !   that vary with time.
    ! References:
    !   None
    !------------------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use error_code, only: &
        err_code, &                     ! Error indicator
        clubb_fatal_error, &            ! Constant
        err_header                      ! String

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(inout) :: gr

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init     ! Initial grid altitude (momentum level) [m]


    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) ::  & 
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    integer :: k, i

    ! ---- Begin Code ----

    if ( .not. l_implemented ) then


      if ( grid_type == 1 ) then

        ! Evenly-spaced grid.
        ! Momentum level altitudes are defined based on the grid starting
        ! altitude, zm_init, the constant grid-spacing, deltaz, and the number
        ! of grid levels, gr%nz.

        ! Define momentum level altitudes. The first momentum level is at
        ! altitude zm_init.
        do k = 1, nz, 1
          do i = 1, ngrdcol
            gr%zm(i,k) = zm_init(i) + real( k-1, kind = core_rknd ) * deltaz(i)
          end do
        end do

        ! Define thermodynamic level altitudes.  Thermodynamic level altitudes
        ! are located at the central altitude levels, exactly halfway between
        ! momentum level altitudes.  The lowermost thermodynamic level is
        ! found by taking 1/2 the altitude difference between the bottom two
        ! momentum levels and subtracting that value from the bottom momentum
        ! level.  The first thermodynamic level is below zm_init.
        do i = 1, ngrdcol
          gr%zt(i,1) = zm_init(i) - ( 0.5_core_rknd * deltaz(i) )
        end do
        
        do k = 2, nz, 1
          do i = 1, ngrdcol
            gr%zt(i,k) = 0.5_core_rknd * ( gr%zm(i,k) + gr%zm(i,k-1) )
          end do
        enddo


      elseif ( grid_type == 2 ) then

        ! Stretched (unevenly-spaced) grid:  stretched thermodynamic levels.
        ! Thermodynamic levels are defined according to a stretched grid that
        ! is entered through the use of an input file.  This is similar to a
        ! SAM-style stretched grid.

        ! Define thermodynamic level altitudes.
        do k = 1, nz, 1
          do i = 1, ngrdcol
            gr%zt(i,k) = thermodynamic_heights(i,k)
          end do
        enddo

        ! Define momentum level altitudes.  Momentum level altitudes are
        ! located at the central altitude levels, exactly halfway between
        ! thermodynamic level altitudes.  The uppermost momentum level
        ! altitude is found by taking 1/2 the altitude difference between the
        ! top two thermodynamic levels and adding that value to the top
        ! thermodynamic level.
        do k = 1, nz-1, 1
          do i = 1, ngrdcol
            gr%zm(i,k) = 0.5_core_rknd * ( gr%zt(i,k+1) + gr%zt(i,k) )
          end do
        enddo
        
        do i = 1, ngrdcol
          gr%zm(i,nz) = gr%zt(i,nz) +  & 
               0.5_core_rknd * ( gr%zt(i,nz) - gr%zt(i,nz-1) )
        end do

      elseif ( grid_type == 3 ) then

        ! Stretched (unevenly-spaced) grid:  stretched momentum levels.
        ! Momentum levels are defined according to a stretched grid that is
        ! entered through the use of an input file.  This is similar to a
        ! WRF-style stretched grid.

        ! Define momentum level altitudes.
        do k = 1, nz, 1
          do i = 1, ngrdcol
            gr%zm(i,k) = momentum_heights(i,k)
          end do
        enddo

        ! Define thermodynamic level altitudes.  Thermodynamic level altitudes
        ! are located at the central altitude levels, exactly halfway between
        ! momentum level altitudes.  The lowermost thermodynamic level
        ! altitude is found by taking 1/2 the altitude difference between the
        ! bottom two momentum levels and subtracting that value from the
        ! bottom momentum level.
        do i = 1, ngrdcol
          gr%zt(i,1) = gr%zm(i,1) - 0.5_core_rknd * ( gr%zm(i,2) - gr%zm(i,1) )
        end do
        
        do k = 2, nz, 1
          do i = 1, ngrdcol
            gr%zt(i,k) = 0.5_core_rknd * ( gr%zm(i,k) + gr%zm(i,k-1) )
          end do
        enddo


      else

        ! Invalid grid type.
        write(fstderr,*) err_header, "Invalid grid type: ", grid_type, & 
                         ".  Valid options are 1, 2, or 3."
        err_code = clubb_fatal_error
        return


      endif


    else

      ! The CLUBB parameterization is implemented in a host model.
      ! Use the host model's momentum level altitudes and thermodynamic level
      ! altitudes to set up the CLUBB grid.

      ! Momentum level altitudes from host model.
      do k = 1, nz, 1
        do i = 1, ngrdcol
          gr%zm(i,k) = momentum_heights(i,k)
        end do
      enddo

      ! Thermodynamic level altitudes from host model after possible grid-index
      ! adjustment for CLUBB interface.
      do k = 1, nz, 1
        do i = 1, ngrdcol
          gr%zt(i,k) = thermodynamic_heights(i,k)
        end do
      enddo


    endif ! not l_implemented


    ! Define dzm, the spacing between thermodynamic grid levels; centered over
    ! momentum grid levels
    do k=1,nz-1
      do i = 1, ngrdcol
        gr%dzm(i,k) = gr%zt(i,k+1) - gr%zt(i,k)
      end do
    enddo
    
    do i = 1, ngrdcol
      gr%dzm(i,nz) = gr%dzm(i,nz-1)
    end do

    ! Define dzt, the spacing between momentum grid levels; centered over
    ! thermodynamic grid levels
    do k=2,nz
      do i = 1, ngrdcol
        gr%dzt(i,k) = gr%zm(i,k) - gr%zm(i,k-1)
      end do
    enddo
    
    do i = 1, ngrdcol
      gr%dzt(i,1) = gr%dzt(i,2)
    end do

    ! Define invrs_dzm, which is the inverse spacing between thermodynamic grid
    ! levels; centered over momentum grid levels.
    do k=1,nz-1
      do i = 1, ngrdcol
        gr%invrs_dzm(i,k) = 1._core_rknd / ( gr%zt(i,k+1) - gr%zt(i,k) )
      end do
    enddo
    
    do i = 1, ngrdcol
      gr%invrs_dzm(i,nz) = gr%invrs_dzm(i,nz-1)
    end do


    ! Define invrs_dzt, which is the inverse spacing between momentum grid
    ! levels; centered over thermodynamic grid levels.
    do k=2,nz
      do i = 1, ngrdcol
        gr%invrs_dzt(i,k) = 1._core_rknd / ( gr%zm(i,k) - gr%zm(i,k-1) )
      end do
    enddo
    
    do i = 1, ngrdcol
      gr%invrs_dzt(i,1) = gr%invrs_dzt(i,2)
    end do


    ! Interpolation Weights: zm grid to zt grid.
    ! The grid index (k) is the index of the level on the thermodynamic (zt)
    ! grid.  The result is the weights of the upper and lower momentum levels
    ! (that sandwich the thermodynamic level) applied to that thermodynamic
    ! level.  These weights are normally used in situations where a momentum
    ! level variable is being solved for implicitly in an equation, and the
    ! aforementioned variable needs to be interpolated from three successive
    ! momentum levels (the central momentum level, as well as one momentum level
    ! above and below the central momentum level) to the intermediate
    ! thermodynamic grid levels that sandwich the central momentum level.
    ! For more information, see the comments in function interpolated_aztk_imp.
    call calc_zm2zt_weights( nz, ngrdcol, & ! In
                             gr )           ! InOut


    ! Interpolation Weights: zt grid to zm grid.
    ! The grid index (k) is the index of the level on the momentum (zm) grid.
    ! The result is the weights of the upper and lower thermodynamic levels
    ! (that sandwich the momentum level) applied to that momentum level.  These
    ! weights are normally used in situations where a thermodynamic level
    ! variable is being solved for implicitly in an equation, and the
    ! aforementioned variable needs to be interpolated from three successive
    ! thermodynamic levels (the central thermodynamic level, as well as one
    ! thermodynamic level above and below the central thermodynamic level) to
    ! the intermediate momentum grid levels that sandwich the central
    ! thermodynamic level.
    ! For more information, see the comments in function interpolated_azmk_imp.
    call calc_zt2zm_weights( nz, ngrdcol, & ! In
                             gr )           ! InOut

    return
  end subroutine setup_grid_heights

  !=============================================================================
  subroutine read_grid_heights( nzmax, grid_type,  & 
                                zm_grid_fname, zt_grid_fname, & 
                                file_unit, &
                                momentum_heights, & 
                                thermodynamic_heights )

    ! Description:
    ! This subroutine is used foremost in cases where the grid_type corresponds
    ! with the stretched (unevenly-spaced) grid options (either grid_type = 2 or
    ! grid_type = 3).  This subroutine reads in the values of the stretched grid
    ! altitude levels for either the thermodynamic level grid or the momentum
    ! level grid.  This subroutine also handles basic error checking for all
    ! three grid types.
    !------------------------------------------------------------------------

    use constants_clubb, only:  & 
        fstderr ! Variable(s)

    use file_functions, only:  & 
        file_read_1d ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        err_code, &                     ! Error indicator
        clubb_fatal_error, &            ! Constant
        err_header                      ! String

    implicit none

    ! Input Variables.

    ! Declared number of vertical levels.
    integer, intent(in) :: & 
      nzmax

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: & 
      grid_type

    character(len=*), intent(in) :: & 
      zm_grid_fname,  & ! Path and filename of file for momentum level altitudes
      zt_grid_fname     ! Path and filename of file for thermodynamic level altitudes

    integer, intent(in) :: &
      file_unit   ! Unit number for zt_grid_fname & zm_grid_fname (based on the OpenMP thread)

    ! Output Variables.

    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), dimension(nzmax), intent(out) :: & 
      momentum_heights,   & ! Momentum level altitudes (file input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

    ! Local Variables.

    integer :: &
      zt_level_count,  & ! Number of altitudes found in zt_grid_fname
      zm_level_count     ! Number of altitudes found in zm_grid_fname

    integer :: input_status   ! Status of file being read:
    ! > 0 ==> error reading file.
    ! = 0 ==> no error and more file to be read.
    ! < 0 ==> end of file indicator.

    ! Generic variable for storing file data while counting the number
    ! of file entries.
    real( kind = core_rknd ) :: generic_input_item

    integer :: k   ! Loop index

    ! ---- Begin Code ----

    ! Declare the momentum level altitude array and the thermodynamic level
    ! altitude array to be 0 until overwritten.
    momentum_heights(1:nzmax) = 0.0_core_rknd
    thermodynamic_heights(1:nzmax) = 0.0_core_rknd

    ! Avoid uninitialized memory
    generic_input_item = 0.0_core_rknd


    if ( grid_type == 1 ) then

      ! Evenly-spaced grid.
      ! Grid level altitudes are based on a constant distance between them and
      ! a starting point for the bottom of the grid.

      ! As a way of error checking, make sure that there isn't any file entry
      ! for either momentum level altitudes or thermodynamic level altitudes.
      if ( zm_grid_fname /= '' ) then
        write(fstderr,*) err_header, & 
           "An evenly-spaced grid has been selected. " & 
           // " Please reset zm_grid_fname to ''."
        err_code = clubb_fatal_error
        return
      endif
      if ( zt_grid_fname /= '' ) then
        write(fstderr,*) err_header, & 
           "An evenly-spaced grid has been selected. " & 
           // " Please reset zt_grid_fname to ''."
        err_code = clubb_fatal_error
        return
      endif


    elseif ( grid_type == 2 ) then

      ! Stretched (unevenly-spaced) grid:  stretched thermodynamic levels.
      ! Thermodynamic levels are defined according to a stretched grid that is
      ! entered through the use of an input file.  Momentum levels are set
      ! halfway between thermodynamic levels.  This is similar to a SAM-style
      ! stretched grid.

      ! As a way of error checking, make sure that there isn't any file entry
      ! for momentum level altitudes.
      if ( zm_grid_fname /= '' ) then
        write(fstderr,*) err_header, & 
           "Thermodynamic level altitudes have been selected " & 
           // "for use in a stretched (unevenly-spaced) grid. " & 
           // " Please reset zm_grid_fname to ''."
        err_code = clubb_fatal_error
        return
      endif

!$omp critical
      ! Open the file zt_grid_fname.
      open( unit=file_unit, file=zt_grid_fname,  &
            status='old', action='read' )

      ! Find the number of thermodynamic level altitudes listed
      ! in file zt_grid_fname.
      zt_level_count = 0
      do
        read( unit=file_unit, fmt=*, iostat=input_status )  & 
           generic_input_item
        if ( input_status < 0 ) exit   ! end of file indicator
        if ( input_status > 0 ) then
          write(fstderr,*) err_header, &  ! error reading input
             "Error reading thermodynamic level input file."
          err_code = clubb_fatal_error
          exit
        end if
        zt_level_count = zt_level_count + 1
      enddo

      ! Close the file zt_grid_fname.
      close( unit=file_unit )
!$omp end critical
      
      if ( err_code == clubb_fatal_error ) return

      ! Check that the number of thermodynamic grid altitudes in the input file
      ! matches the declared number of CLUBB grid levels (nzmax).
      if ( zt_level_count /= nzmax ) then
        write(fstderr,*)  & 
           "The number of thermodynamic grid altitudes " & 
           // "listed in file " // trim(zt_grid_fname)  & 
           // " does not match the number of CLUBB grid " & 
           // "levels specified in the model.in file."
        write(fstderr,*) & 
           "Number of thermodynamic grid altitudes listed:  ", & 
           zt_level_count
        write(fstderr,*) & 
           "Number of CLUBB grid levels specified:  ", nzmax
        err_code = clubb_fatal_error
        return
      endif

      ! Read the thermodynamic level altitudes from zt_grid_fname.
      call file_read_1d( file_unit, zt_grid_fname, nzmax, 1,  & ! intent(in)
                         thermodynamic_heights ) ! intent(out)

      ! Check that each thermodynamic level altitude increases
      ! in height as the thermodynamic level grid index increases.
      do k = 2, nzmax, 1
        if ( thermodynamic_heights(k)  & 
             <= thermodynamic_heights(k-1) ) then
          write(fstderr,*)  & 
             "The declared thermodynamic level grid " & 
             // "altitudes are not increasing in height " & 
             // "as grid level index increases."
          write(fstderr,*) & 
             "Grid index:  ", k-1, ";", & 
             "  Thermodynamic level altitude:  ", & 
             thermodynamic_heights(k-1)
          write(fstderr,*) & 
             "Grid index:  ", k, ";", & 
             "  Thermodynamic level altitude:  ", & 
             thermodynamic_heights(k)
          err_code = clubb_fatal_error
          return
        endif
      enddo


    elseif ( grid_type == 3 ) then

      ! Stretched (unevenly-spaced) grid:  stretched momentum levels.
      ! Momentum levels are defined according to a stretched grid that is
      ! entered through the use of an input file.  Thermodynamic levels are set
      ! halfway between momentum levels.  This is similar to a WRF-style
      ! stretched grid.

      ! As a way of error checking, make sure that there isn't any file entry
      ! for thermodynamic level altitudes.
      if ( zt_grid_fname /= '' ) then
        write(fstderr,*) & 
           "Momentum level altitudes have been selected " & 
           // "for use in a stretched (unevenly-spaced) grid. " & 
           // " Please reset zt_grid_fname to ''."
        err_code = clubb_fatal_error
        return
      endif

      ! Open the file zm_grid_fname.
      open( unit=file_unit, file=zm_grid_fname,  & 
            status='old', action='read' )

      ! Find the number of momentum level altitudes
      ! listed in file zm_grid_fname.
      zm_level_count = 0
      do
        read( unit=file_unit, fmt=*, iostat=input_status ) & 
           generic_input_item
        if ( input_status < 0 ) exit   ! end of file indicator
        if ( input_status > 0 ) then

          write(fstderr,*) err_header, &! error reading input
                           "Error reading momentum level input file."
          err_code = clubb_fatal_error
          return
        end if
        zm_level_count = zm_level_count + 1
      enddo

      ! Close the file zm_grid_fname.
      close( unit=file_unit )

      ! Check that the number of momentum grid altitudes in the input file
      ! matches the declared number of CLUBB grid levels (nzmax).
      if ( zm_level_count /= nzmax ) then
        write(fstderr,*) & 
           "The number of momentum grid altitudes " & 
           // "listed in file " // trim(zm_grid_fname) & 
           // " does not match the number of CLUBB grid " & 
           // "levels specified in the model.in file."
        write(fstderr,*) & 
           "Number of momentum grid altitudes listed:  ", & 
           zm_level_count
        write(fstderr,*) & 
           "Number of CLUBB grid levels specified:  ", nzmax
        err_code = clubb_fatal_error
        return
      endif

      ! Read the momentum level altitudes from zm_grid_fname.
      call file_read_1d( file_unit, zm_grid_fname, nzmax, 1, & ! intent(in) 
                         momentum_heights ) ! intent(out)

      ! Check that each momentum level altitude increases in height as the
      ! momentum level grid index increases.
      do k = 2, nzmax, 1
        if ( momentum_heights(k)  & 
             <= momentum_heights(k-1) ) then
          write(fstderr,*)  & 
             "The declared momentum level grid " & 
             // "altitudes are not increasing in height " & 
             // "as grid level index increases."
          write(fstderr,*) & 
             "Grid index:  ", k-1, ";", & 
             "  Momentum level altitude:  ", & 
             momentum_heights(k-1)
          write(fstderr,*) & 
             "Grid index:  ", k, ";", & 
             "  Momentum level altitude:  ", & 
             momentum_heights(k)
          err_code = clubb_fatal_error
          return
        endif
      enddo


    endif


    ! The purpose of this if statement is to avoid a compiler warning.
    if ( generic_input_item > 0.0_core_rknd ) then
      ! Do nothing
    endif
    ! Joshua Fasching  June 2008

    return

  end subroutine read_grid_heights
  
  !=============================================================================
  ! Wrapped in interface zt2zm
  function redirect_interpolated_azm_k( gr, azt, k )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none

    !---------------------------- Input Variables ----------------------------
    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]
      
    integer, intent(in) :: &
      k

    !---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ) :: &
      redirect_interpolated_azm_k    ! Variable when interp. to momentum levels

    !---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(gr%nz) :: &
      redirect_interpolated_azm_col    ! Variable when interp. to momentum levels
      
    !---------------------------- Begin Code ----------------------------

    redirect_interpolated_azm_col = redirect_interpolated_azm_1D(gr, azt)

    redirect_interpolated_azm_k = redirect_interpolated_azm_col(k)

    return
  end function redirect_interpolated_azm_k
  
  !=============================================================================
  ! Wrapped in interface zt2zm
  function redirect_interpolated_azm_1D( gr, azt )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none

    !---------------------------- Input Variables ----------------------------
    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    !---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ), dimension(gr%nz) :: &
      redirect_interpolated_azm_1D    ! Variable when interp. to momentum levels

    !---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      azt_col    ! Variable on thermodynamic grid levels    [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      redirect_interpolated_azm_col    ! Variable when interp. to momentum levels

    !---------------------------- Begin Code ----------------------------

    azt_col(1,:) = azt

    redirect_interpolated_azm_col = redirect_interpolated_azm_2D(gr%nz, 1, gr, azt_col)

    redirect_interpolated_azm_1D = redirect_interpolated_azm_col(1,:)

    return
  end function redirect_interpolated_azm_1D

  !=============================================================================
  ! Wrapped in interface zt2zm
  function redirect_interpolated_azm_2D( nz, ngrdcol, gr, azt )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none
    
    !---------------------------- Input Variables ----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    !---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      redirect_interpolated_azm_2D    ! Variable when interp. to momentum levels
      
    !---------------------------- Begin Code ----------------------------

    ! Sanity Check
    if (l_quintic_poly_interp) then
      if (.not. l_cubic_interp) then
        write (fstdout, *) "Error: Model flag l_quintic_poly_interp should not be true if "&
                         //"l_cubic_interp is false."
        error stop
      end if
    end if

    ! Redirect
    if ( l_cubic_interp .and. nz >= 3 ) then
      redirect_interpolated_azm_2D = cubic_interpolated_azm_2D( nz, ngrdcol, gr, azt )
    else
      call linear_interpolated_azm_2D( nz, ngrdcol, gr, azt, &
                                       redirect_interpolated_azm_2D )
    end if

    return
  end function redirect_interpolated_azm_2D
  
  !=============================================================================
  ! Wrapped in interface zm2zt
  function redirect_interpolated_azt_k( gr, azm, k )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none

    type (grid), target, intent(in) :: gr

    !---------------------------- Input Variables ----------------------------
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]
      
    integer, intent(in) :: &
      k   ! Vertical level

    !---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ) :: &
      redirect_interpolated_azt_k    ! Variable when interp. to momentum levels

    !---------------------------- Local Variables ----------------------------

    real( kind = core_rknd ), dimension(gr%nz) :: &
      redirect_interpolated_azt_col    ! Variable when interp. to momentum levels
    
    !---------------------------- Begin Code ----------------------------

    redirect_interpolated_azt_col = redirect_interpolated_azt_1D(gr, azm)

    redirect_interpolated_azt_k = redirect_interpolated_azt_col(k)

    return
  end function redirect_interpolated_azt_k
  

  !=============================================================================
  ! Wrapped in interface zm2zt
  function redirect_interpolated_azt_1D( gr, azm )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none

    type (grid), target, intent(in) :: gr

    !---------------------------- Input Variables ----------------------------
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    !---------------------------- Return Variable ----------------------------
    real( kind = core_rknd ), dimension(gr%nz) :: &
      redirect_interpolated_azt_1D    ! Variable when interp. to momentum levels

    !---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      azm_col    ! Variable on thermodynamic grid levels    [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      redirect_interpolated_azt_col    ! Variable when interp. to momentum levels
    
    !---------------------------- Begin Code ----------------------------

    azm_col(1,:) = azm

    redirect_interpolated_azt_col = redirect_interpolated_azt_2D(gr%nz, 1, gr, azm_col)

    redirect_interpolated_azt_1D = redirect_interpolated_azt_col(1,:)

    return
  end function redirect_interpolated_azt_1D

  !=============================================================================
  ! Wrapped in interface zm2zt
  function redirect_interpolated_azt_2D( nz, ngrdcol, gr, azm )

    ! Description:
    ! Calls the appropriate corresponding function based on l_cubic_temp
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_cubic_interp, & ! Variable(s)
        l_quintic_poly_interp

    use constants_clubb, only: &
        fstdout ! Variable

    implicit none
    
    !---------------------------- Input Variables ----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    !---------------------------- Return Variable----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      redirect_interpolated_azt_2D    ! Variable when interp. to momentum levels
    
    !---------------------------- Begin Code ----------------------------

    ! Sanity Check
    if (l_quintic_poly_interp) then
      if (.not. l_cubic_interp) then
        write (fstdout, *) "Error: Model flag l_quintic_poly_interp should not be true if "&
                         //"l_cubic_interp is false."
        error stop
      end if
    end if

    ! Redirect
    if (l_cubic_interp .and. nz >= 3 ) then
      redirect_interpolated_azt_2D = cubic_interpolated_azt_2D( nz, ngrdcol, gr, azm )
    else
      call linear_interpolated_azt_2D( nz, ngrdcol, gr, azm, &
                                       redirect_interpolated_azt_2D )
    end if

    return
  end function redirect_interpolated_azt_2D
  
  !=============================================================================
  subroutine linear_interpolated_azm_2D( nz, ngrdcol, gr, azt, &
                                              linear_interpolated_azm )

    ! Description:
    ! Function to interpolate a variable located on the thermodynamic grid
    ! levels (azt) to the momentum grid levels (azm).  This function inputs the
    ! entire azt array and outputs the results as an azm array.  The
    ! formulation used is compatible with a stretched (unevenly-spaced) grid.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use interpolation, only: &
        linear_interp_factor  ! Procedure(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    ! ------------------------------ Input Variable ------------------------------
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! ------------------------------ Return Variable ------------------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz) :: &
      linear_interpolated_azm    ! Variable when interp. to momentum levels

    ! ------------------------------ Local Variable ------------------------------
    integer :: i, k  ! Grid level loop index

    ! ------------------------------ Begin Code ------------------------------

    !$acc data copyin( azt, gr, gr%weights_zt2zm, gr%zt, gr%zm ) &
    !$acc      copyout( linear_interpolated_azm )

    ! Interpolate the value of a thermodynamic-level variable to the central
    ! momentum level, k, between two successive thermodynamic levels using
    ! linear interpolation.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz-1
      do i = 1, ngrdcol
        linear_interpolated_azm(i,k) = gr%weights_zt2zm(i,k,1) &
                                          * ( azt(i,k+1) - azt(i,k) ) + azt(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Set the value of the thermodynamic-level variable, azt, at the uppermost
    ! level of the model, which is a momentum level.  The name of the variable
    ! when interpolated/extended to momentum levels is azm.
    ! Use a linear extension based on the values of azt at levels gr%nz and
    ! gr%nz-1 to find the value of azm at level gr%nz (the uppermost level
    ! in the model).
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      linear_interpolated_azm(i,nz) &
        = ( ( azt(i,nz) - azt(i,nz-1) ) / ( gr%zt(i,nz) - gr%zt(i,nz-1) ) ) & 
          * ( gr%zm(i,nz) - gr%zt(i,nz) ) + azt(i,nz)
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine linear_interpolated_azm_2D

  !=============================================================================
  function zt2zm2zt( nz, ngrdcol, gr, azt )

    ! Description:
    !    Function to interpolate a variable located on the thermodynamic grid
    !    levels (azt) to the momentum grid levels (azm), then interpolate back
    !    to thermodynamic grid levels (azt).
    !
    ! Note:
    !   This is intended for smoothing variables.
    !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none
    
    ! ------------------------------ Input Variable ------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! ------------------------------ Return Variable ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zt2zm2zt    ! Variable when interp. to momentum levels

    ! ------------------------------ Local Variable ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      azt_zm

    ! ------------------------------ Begin Code ------------------------------

    !$acc data copyin( azt ) &
    !$acc     copyout( zt2zm2zt ) &
    !$acc      create( azt_zm )

    ! Interpolate azt to momentum levels 
    azt_zm = zt2zm( nz, ngrdcol, gr, azt )

    ! Interpolate back to termodynamic levels
    zt2zm2zt = zm2zt( nz, ngrdcol, gr, azt_zm )

    !$acc end data

    return 

  end function zt2zm2zt
  
  !=============================================================================
  function zm2zt2zm( nz, ngrdcol, gr, azm )

    ! Description:
    !    Function to interpolate a variable located on the momentum grid 
    !    levels(azm) to thermodynamic grid levels (azt), then interpolate 
    !    back to momentum grid levels (azm).
    !
    ! Note:
    !   This is intended for smoothing variables.
    !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none
    
    ! ------------------------------ Input Variable ------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! ------------------------------ Return Variable ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zm2zt2zm    ! Variable when interp. to momentum levels

    ! ------------------------------ Local Variable ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      azm_zt

    ! ------------------------------ Begin Code ------------------------------

    !$acc data copyin( azm ) &
    !$acc     copyout( zm2zt2zm ) &
    !$acc      create( azm_zt )

    ! Interpolate azt to termodynamic levels 
    azm_zt = zm2zt( nz, ngrdcol, gr, azm )

    ! Interpolate back to momentum levels
    zm2zt2zm = zt2zm( nz, ngrdcol, gr, azm_zt )

    !$acc end data

    return 

  end function zm2zt2zm

  !=============================================================================
  function cubic_interpolated_azm_2D( nz, ngrdcol, gr, azt )

    ! Description:
    !   Function to interpolate a variable located on the momentum grid
    !   levels (azt) to the thermodynamic grid levels (azm).  This function outputs the
    !   value of azt at a all grid levels using Steffen's monotonic cubic
    !   interpolation implemented by Tak Yamaguchi.
    ! 
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use interpolation, only: &
        mono_cubic_interp  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azt

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      cubic_interpolated_azm_2D

    ! Local Variable(s)
    integer :: &
      k, i
      
    integer :: km1, k00, kp1, kp2

    ! ---- Begin Code ----

    do k = 1, nz 
      do i = 1, ngrdcol

        ! k levels are based on Tak's find_indices subroutine -dschanen 24 Oct 2011
        if ( k == nz-1 ) then
          km1 = nz-2
          kp1 = nz
          kp2 = nz
          k00 = nz-1
        else if ( k == nz ) then ! Extrapolation
          km1 = nz
          kp1 = nz
          kp2 = nz
          k00 = nz-1
        else if ( k == 1 ) then
          km1 = 1
          kp1 = 2
          kp2 = 3
          k00 = 1
        else
          km1 = k-1
          kp1 = k+1
          kp2 = k+2
          k00 = k
        end if
        
        ! Do the actual interpolation.
        ! Use a cubic monotonic spline interpolation.
        cubic_interpolated_azm_2D(i,k) = &
          mono_cubic_interp( gr%zm(i,k), km1, k00, kp1, kp2, &
                             gr%zt(i,km1), gr%zt(i,k00), gr%zt(i,kp1), gr%zt(i,kp2), &
                             azt(i,km1), azt(i,k00), azt(i,kp1), azt(i,kp2) )
      end do
    end do

    return

  end function cubic_interpolated_azm_2D

  !=============================================================================
  subroutine calc_zt2zm_weights( nz, ngrdcol, &
                                      gr ) 

    ! Description:
    ! Function used to help in an interpolation of a variable (var_zt) located
    ! on the thermodynamic grid levels (azt) to the momentum grid levels (azm).
    ! This function computes a weighting factor for both the upper thermodynamic
    ! level (k+1) and the lower thermodynamic level (k) applied to the central
    ! momentum level (k).  For the uppermost momentum grid level (k=gr%nz), a
    ! weighting factor for both the thermodynamic level at gr%nz and the
    ! thermodynamic level at gr%nz-1 are calculated based on the use of a
    ! linear extension.  This function outputs the weighting factors at a single
    ! momentum grid level (k).  The formulation used is compatible with a
    ! stretched (unevenly-spaced) grid.  The weights are defined as follows:
    !
    ! ---var_zt(k+1)------------------------------------------- t(k+1)
    !                       azt_weight(t_above) = factor
    ! ===========var_zt(interp)================================ m(k)
    !                       azt_weight(t_below) = 1 - factor
    ! ---var_zt(k)--------------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! For all levels k < gr%nz:
    !
    ! The formula for a linear interpolation is given by:
    !
    ! var_zt( interp to zm(k) )
    ! = [ ( var_zt(k+1) - var_zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !     * ( zm(k) - zt(k) ) + var_zt(k);
    !
    ! which can be rewritten as:
    !
    ! var_zt( interp to zm(k) )
    ! = [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !     * ( var_zt(k+1) - var_zt(k) ) + var_zt(k).
    !
    ! Furthermore, the formula can be rewritten as:
    !
    ! var_zt( interp to zm(k) )
    ! = factor * var_zt(k+1) + ( 1 - factor ) * var_zt(k);
    !
    ! where:
    !
    ! factor = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ).
    !
    ! One of the important uses of this function is in situations where the
    ! variable to be interpolated is being treated IMPLICITLY in an equation.
    ! Usually, the variable to be interpolated is involved in a derivative (such
    ! as d(var_zt)/dz in the diagram below).  For the term of the equation
    ! containing the derivative, grid weights are needed for two interpolations,
    ! rather than just one interpolation.  Thus, four grid weights (labeled
    ! A(k), B(k), C(k), and D(k) in the diagram below) are needed.
    !
    ! ---var_zt(k+1)------------------------------------------- t(k+1)
    !                                       A(k)
    ! ===========var_zt(interp)================================ m(k)
    !                                       B(k) = 1 - A(k)
    ! ---var_zt(k)-----------d(var_zt)/dz---------------------- t(k)
    !                                       C(k)
    ! ===========var_zt(interp)================================ m(k-1)
    !                                       D(k) = 1 - C(k)
    ! ---var_zt(k-1)------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! The grid weights, indexed around the central thermodynamic level (k), are
    ! defined as follows:
    !
    ! A(k) = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) );
    !
    ! which is the same as "factor" for the interpolation to momentum
    ! level (k).  In the code, this interpolation is referenced as
    ! gr%weights_zt2zm(t_above,mk), which can be read as "grid weight in a zt2zm
    ! interpolation of the thermodynamic level above momentum level (k) (applied
    ! to momentum level (k))".
    !
    ! B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !      = 1 - A(k);
    !
    ! which is the same as "1 - factor" for the interpolation to momentum
    ! level (k).  In the code, this interpolation is referenced as
    ! gr%weights_zt2zm(t_below,mk), which can be read as "grid weight in a zt2zm
    ! interpolation of the thermodynamic level below momentum level (k) (applied
    ! to momentum level (k))".
    !
    ! C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) );
    !
    ! which is the same as "factor" for the interpolation to momentum
    ! level (k-1).  In the code, this interpolation is referenced as
    ! gr%weights_zt2zm(t_above,mkm1), which can be read as "grid weight in a
    ! zt2zm interpolation of the thermodynamic level above momentum level (k-1)
    ! (applied to momentum level (k-1))".
    !
    ! D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !      = 1 - C(k);
    !
    ! which is the same as "1 - factor" for the interpolation to momentum
    ! level (k-1).  In the code, this interpolation is referenced as
    ! gr%weights_zt2zm(t_below,mkm1), which can be read as "grid weight in a
    ! zt2zm interpolation of the thermodynamic level below momentum level (k-1)
    ! (applied to momentum level (k-1))".
    !
    ! Additionally, as long as the central thermodynamic level (k) in the above
    ! scenario is not the uppermost thermodynamic level or the lowermost
    ! thermodynamic level (k /= gr%nz and k /= 1), the four weighting factors
    ! have the following relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    !
    !
    ! Special condition for uppermost grid level, k = gr%nz:
    !
    ! The uppermost momentum grid level is above the uppermost thermodynamic
    ! grid level.  Thus, a linear extension is used at this level.
    !
    ! For level k = gr%nz:
    !
    ! The formula for a linear extension is given by:
    !
    ! var_zt( extend to zm(k) )
    ! = [ ( var_zt(k) - var_zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !     * ( zm(k) - zt(k-1) ) + var_zt(k-1);
    !
    ! which can be rewritten as:
    !
    ! var_zt( extend to zm(k) )
    ! = [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !     * ( var_zt(k) - var_zt(k-1) ) + var_zt(k-1).
    !
    ! Furthermore, the formula can be rewritten as:
    !
    ! var_zt( extend to zm(k) )
    ! = factor * var_zt(k) + ( 1 - factor ) * var_zt(k-1);
    !
    ! where:
    !
    ! factor = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ).
    !
    ! Due to the fact that a linear extension is being used, the value of factor
    ! will be greater than 1.  The weight of thermodynamic level k = gr%nz on
    ! momentum level k = gr%nz equals the value of factor.  The weight of
    ! thermodynamic level k = gr%nz-1 on momentum level k = gr%nz equals
    ! 1 - factor, which is less than 0.  However, the sum of the two weights
    ! equals 1.
    !
    !
    ! Brian Griffin; September 12, 2008.
    !
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      t_above = 1,  & ! Upper thermodynamic level.
      t_below = 2     ! Lower thermodynamic level.

    ! Input Variable
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    ! Input/Output Variable
    type (grid), target, intent(inout) :: gr

    ! Local Variables
    real( kind = core_rknd ) :: factor

    integer :: i, k

    do k = 1, nz-1
      do i = 1, ngrdcol
        
        ! The top model level (gr%nz) is formulated differently because the top
        ! momentum level is above the top thermodynamic level.  A linear
        ! extension is required, rather than linear interpolation.
        ! Note:  Variable "factor" will be greater than 1 in this situation.
        gr%weights_zt2zm(i,k,t_above) = ( gr%zm(i,k) - gr%zt(i,k) ) &
                                          / ( gr%zt(i,k+1) - gr%zt(i,k) )

        ! Weight of lower thermodynamic level on momentum level.
        gr%weights_zt2zm(i,k,t_below) = one - gr%weights_zt2zm(i,k,t_above)
        
      end do
    end do
    
    ! At most levels, the momentum level is found in-between two
    ! thermodynamic levels.  Linear interpolation is used.
    do i = 1, ngrdcol
      gr%weights_zt2zm(i,nz,t_above) = ( gr%zm(i,nz) - gr%zt(i,nz-1) ) &
                                        / ( gr%zt(i,nz) - gr%zt(i,nz-1) )
                                        
      gr%weights_zt2zm(i,nz,t_below) = one - gr%weights_zt2zm(i,nz,t_above)
    end do

    return

  end subroutine calc_zt2zm_weights
  
  !=============================================================================
  subroutine linear_interpolated_azt_2D( nz, ngrdcol, gr, azm, &
                                              linear_interpolated_azt )

    ! Description:
    ! Function to interpolate a variable located on the momentum grid levels
    ! (azm) to the thermodynamic grid levels (azt).  This function inputs the
    ! entire azm array and outputs the results as an azt array.  The formulation
    ! used is compatible with a stretched (unevenly-spaced) grid.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use interpolation, only: &
        linear_interp_factor  ! Procedure(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    ! ------------------------------ Input Variable ------------------------------
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! ------------------------------ Output Variable ------------------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz) :: &
      linear_interpolated_azt    ! Variable when interp. to thermodynamic levels

    ! ------------------------------ Local Variable ------------------------------
    integer :: i, k  ! Grid level loop index

    ! ------------------------------ Begin Code ------------------------------

    !$acc data copyin( azm, gr, gr%weights_zm2zt, gr%zt, gr%zm ) &
    !$acc      copyout( linear_interpolated_azt )

    ! Set the value of the momentum-level variable, azm, at the lowermost level
    ! of the model (below the model lower boundary), which is a thermodynamic
    ! level.  The name of the variable when interpolated/extended to
    ! thermodynamic levels is azt.
    ! Use a linear extension based on the values of azm at levels 1 and 2 to
    ! find the value of azt at level 1 (the lowermost level in the model).
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      linear_interpolated_azt(i,1) &
        = ( ( azm(i,2) - azm(i,1) ) / ( gr%zm(i,2) - gr%zm(i,1) ) ) & 
          * ( gr%zt(i,1) - gr%zm(i,1) ) + azm(i,1)
    end do
    !$acc end parallel loop

    ! Interpolate the value of a momentum-level variable to the central
    ! thermodynamic level, k, between two successive momentum levels using
    ! linear interpolation.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 2, nz
      do i = 1, ngrdcol
        linear_interpolated_azt(i,k) = gr%weights_zm2zt(i,k,1) &
                                       * ( azm(i,k) - azm(i,k-1) ) + azm(i,k-1)
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine linear_interpolated_azt_2D
  
  !=============================================================================
  function cubic_interpolated_azt_2D( nz, ngrdcol, gr, azm )

    ! Description:
    !   Function to interpolate a variable located on the momentum grid
    !   levels (azm) to the thermodynamic grid levels (azt).  This function outputs the
    !   value of azt at a all grid levels using Steffen's monotonic cubic
    !   interpolation implemented by Tak Yamaguchi.
    ! 
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use interpolation, only: &
        mono_cubic_interp  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      cubic_interpolated_azt_2D

    ! Local Variable(s)
    integer :: &
      k, i
      
    integer :: km1, k00, kp1, kp2

    ! ---- Begin Code ----

    do k = 1, nz 
      do i = 1, ngrdcol

        ! k levels are based on Tak's find_indices subroutine -dschanen 24 Oct 2011
        if ( k == nz ) then
          km1 = nz-2
          kp1 = nz
          kp2 = nz
          k00 = nz-1
        else if ( k == 2 ) then
          km1 = 1
          kp1 = 2
          kp2 = 3
          k00 = 1
        else if ( k == 1 ) then ! Extrapolation for the ghost point
          km1 = nz
          k00 = 1
          kp1 = 2
          kp2 = 3
        else
          km1 = k-2
          kp1 = k
          kp2 = k+1
          k00 = k-1
        end if
        
        ! Do the actual interpolation.
        ! Use a cubic monotonic spline interpolation.
        cubic_interpolated_azt_2D(i,k) = &
          mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, &
                             gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                             azm(i,km1), azm(i,k00), azm(i,kp1), azm(i,kp2) )
      end do
    end do

    return

  end function cubic_interpolated_azt_2D

  !=============================================================================
  subroutine calc_zm2zt_weights( nz, ngrdcol, &
                                      gr )

    ! Description:
    ! Function used to help in an interpolation of a variable (var_zm) located
    ! on the momentum grid levels (azm) to the thermodynamic grid levels (azt).
    ! This function computes a weighting factor for both the upper momentum
    ! level (k) and the lower momentum level (k-1) applied to the central
    ! thermodynamic level (k).  For the lowermost thermodynamic grid
    ! level (k=1), a weighting factor for both the momentum level at 1 and the
    ! momentum level at 2 are calculated based on the use of a linear extension.
    ! This function outputs the weighting factors at a single thermodynamic grid
    ! level (k).   The formulation used is compatible with a stretched
    ! (unevenly-spaced) grid.  The weights are defined as follows:
    !
    ! ===var_zm(k)============================================= m(k)
    !                       azm_weight(m_above) = factor
    ! -----------var_zm(interp)-------------------------------- t(k)
    !                       azm_weight(m_below) = 1 - factor
    ! ===var_zm(k-1)=========================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! For all levels k > 1:
    !
    ! The formula for a linear interpolation is given by:
    !
    ! var_zm( interp to zt(k) )
    ! = [ ( var_zm(k) - var_zm(k-1) ) / ( zm(k) - zm(k-1) ) ]
    !     * ( zt(k) - zm(k-1) ) + var_zm(k-1);
    !
    ! which can be rewritten as:
    !
    ! var_zm( interp to zt(k) )
    ! = [ ( zt(k) - zm(k-1) ) / ( zm(k) - zm(k-1) ) ]
    !     * ( var_zm(k) - var_zm(k-1) ) + var_zm(k-1).
    !
    ! Furthermore, the formula can be rewritten as:
    !
    ! var_zm( interp to zt(k) )
    ! = factor * var_zm(k) + ( 1 - factor ) * var_zm(k-1);
    !
    ! where:
    !
    ! factor = ( zt(k) - zm(k-1) ) / ( zm(k) - zm(k-1) ).
    !
    ! One of the important uses of this function is in situations where the
    ! variable to be interpolated is being treated IMPLICITLY in an equation.
    ! Usually, the variable to be interpolated is involved in a derivative (such
    ! as d(var_zm)/dz in the diagram below).  For the term of the equation
    ! containing the derivative, grid weights are needed for two interpolations,
    ! rather than just one interpolation.  Thus, four grid weights (labeled
    ! A(k), B(k), C(k), and D(k) in the diagram below) are needed.
    !
    ! ===var_zm(k+1)=========================================== m(k+1)
    !                                       A(k)
    ! -----------var_zm(interp)-------------------------------- t(k+1)
    !                                       B(k) = 1 - A(k)
    ! ===var_zm(k)===========d(var_zm)/dz====================== m(k)
    !                                       C(k)
    ! -----------var_zm(interp)-------------------------------- t(k)
    !                                       D(k) = 1 - C(k)
    ! ===var_zm(k-1)=========================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! The grid weights, indexed around the central momentum level (k), are
    ! defined as follows:
    !
    ! A(k) = ( zt(k+1) - zm(k) ) / ( zm(k+1) - zm(k) );
    !
    ! which is the same as "factor" for the interpolation to thermodynamic
    ! level (k+1).  In the code, this interpolation is referenced as
    ! gr%weights_zm2zt(m_above,tkp1), which can be read as "grid weight in a
    ! zm2zt interpolation of the momentum level above thermodynamic
    ! level (k+1) (applied to thermodynamic level (k+1))".
    !
    ! B(k) = 1 - [ ( zt(k+1) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
    !      = 1 - A(k);
    !
    ! which is the same as "1 - factor" for the interpolation to thermodynamic
    ! level (k+1).  In the code, this interpolation is referenced as
    ! gr%weights_zm2zt(m_below,tkp1), which can be read as "grid weight in a
    ! zm2zt interpolation of the momentum level below thermodynamic
    ! level (k+1) (applied to thermodynamic level (k+1))".
    !
    ! C(k) = ( zt(k) - zm(k-1) ) / ( zm(k) - zm(k-1) );
    !
    ! which is the same as "factor" for the interpolation to thermodynamic
    ! level (k).  In the code, this interpolation is referenced as
    ! gr%weights_zm2zt(m_above,tk), which can be read as "grid weight in a zm2zt
    ! interpolation of the momentum level above thermodynamic level (k) (applied
    ! to thermodynamic level (k))".
    !
    ! D(k) = 1 - [ ( zt(k) - zm(k-1) ) / ( zm(k) - zm(k-1) ) ]
    !      = 1 - C(k);
    !
    ! which is the same as "1 - factor" for the interpolation to thermodynamic
    ! level (k).  In the code, this interpolation is referenced as
    ! gr%weights_zm2zt(m_below,tk), which can be read as "grid weight in a zm2zt
    ! interpolation of the momentum level below thermodynamic level (k) (applied
    ! to thermodynamic level (k))".
    !
    ! Additionally, as long as the central momentum level (k) in the above
    ! scenario is not the lowermost momentum level or the uppermost momentum
    ! level (k /= 1 and k /= gr%nz), the four weighting factors have the
    ! following relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    !
    !
    ! Special condition for lowermost grid level, k = 1:
    !
    ! The lowermost thermodynamic grid level is below the lowermost momentum
    ! grid level.  Thus, a linear extension is used at this level.  It should
    ! be noted that the thermodynamic level k = 1 is considered to be below the
    ! model lower boundary, which is defined to be at momentum level k = 1.
    ! Thus, the values of most variables at thermodynamic level k = 1 are not
    ! often needed or referenced.
    !
    ! For level k = 1:
    !
    ! The formula for a linear extension is given by:
    !
    ! var_zm( extend to zt(k) )
    ! = [ ( var_zm(k+1) - var_zm(k) ) / ( zm(k+1) - zm(k) ) ]
    !     * ( zt(k) - zm(k) ) + var_zm(k);
    !
    ! which can be rewritten as:
    !
    ! var_zm( extend to zt(k) )
    ! = [ ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ) ]
    !     * ( var_zm(k+1) - var_zm(k) ) + var_zm(k).
    !
    ! Furthermore, the formula can be rewritten as:
    !
    ! var_zm( extend to zt(k) )
    ! = factor * var_zm(k+1) + ( 1 - factor ) * var_zm(k);
    !
    ! where:
    !
    ! factor = ( zt(k) - zm(k) ) / ( zm(k+1) - zm(k) ).
    !
    ! Due to the fact that a linear extension is being used, the value of factor
    ! will be less than 0.  The weight of the upper momentum level, which is
    ! momentum level k = 2, on thermodynamic level k = 1 equals the value of
    ! factor.  The weight of the lower momentum level, which is momentum level
    ! k = 1, on thermodynamic level k = 1 equals 1 - factor, which is greater
    ! than 1.  However, the sum of the weights equals 1.
    !
    !
    ! Brian Griffin; September 12, 2008.
    !
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      m_above = 1,  & ! Upper momentum level.
      m_below = 2     ! Lower momentum level.

    ! Input Variable
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    ! Input/Output Variable
    type (grid), target, intent(inout) :: gr

    integer :: i, k

    
    ! At most levels, the momentum level is found in-between two
    ! thermodynamic levels.  Linear interpolation is used.
    do i = 1, ngrdcol
      gr%weights_zm2zt(i,1,m_above) = ( gr%zt(i,1) - gr%zm(i,1) ) &
                                        / ( gr%zm(i,2) - gr%zm(i,1) )
                                        
      gr%weights_zm2zt(i,1,m_below) = one - gr%weights_zm2zt(i,1,m_above)
    end do

    do k = 2, nz
      do i = 1, ngrdcol
        
        ! The top model level (gr%nz) is formulated differently because the top
        ! momentum level is above the top thermodynamic level.  A linear
        ! extension is required, rather than linear interpolation.
        ! Note:  Variable "factor" will be greater than 1 in this situation.
        gr%weights_zm2zt(i,k,m_above) = ( gr%zt(i,k) - gr%zm(i,k-1) ) &
                                         / ( gr%zm(i,k) - gr%zm(i,k-1) )

        ! Weight of lower thermodynamic level on momentum level.
        gr%weights_zm2zt(i,k,m_below) = one - gr%weights_zm2zt(i,k,m_above)
        
      end do
    end do
    
    return

  end subroutine calc_zm2zt_weights
  
  !=============================================================================
  ! Wrapped in interface ddzm
  function gradzm_2D( nz, ngrdcol, gr, azm )

    ! Description:
    !  2D version of gradzm
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      gradzm_2D    ! Vertical derivative of azm    [units vary / m]

    ! Local Variable
    integer :: i, k  ! Grid level loop index
    
    !$acc data copyin( gr, gr%invrs_dzt, azm ) &
    !$acc     copyout( gradzm_2D )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      gradzm_2D(i,1) = ( azm(i,2) - azm(i,1) ) * gr%invrs_dzt(i,2)
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 2, nz
      do i = 1, ngrdcol
        gradzm_2D(i,k) = ( azm(i,k) - azm(i,k-1) ) * gr%invrs_dzt(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function gradzm_2D
  
  !=============================================================================
  ! Wrapped in interface ddzm
  function gradzm_1D( gr, azm )

    ! Description:
    !  2D version of gradzm
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      gradzm_1D    ! Vertical derivative of azm    [units vary / m]
      
    ! Local Variables
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      azm_col    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      gradzm_1D_col    ! Vertical derivative of azm    [units vary / m]

    azm_col(1,:) = azm
    
    gradzm_1D_col = gradzm_2D(gr%nz, 1, gr, azm_col)
      
    gradzm_1D = gradzm_1D_col(1,:)

    return

  end function gradzm_1D
  
  !=============================================================================
  ! Wrapped in interface ddzt
  function gradzt_2D( nz, ngrdcol, gr, azt )

    ! Description:
    !  2D version of gradzt
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      gradzt_2D    ! Vertical derivative of azt    [units vary / m]

    ! Local Variable
    integer :: i, k  ! Grid level loop index

    !$acc data copyin( gr, gr%invrs_dzm, azt ) &
    !$acc     copyout( gradzt_2D )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      gradzt_2D(i,nz) = ( azt(i,nz) - azt(i,nz-1) ) * gr%invrs_dzm(i,nz-1)
    end do
    !$acc end parallel loop
    
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz-1
      do i = 1, ngrdcol
        gradzt_2D(i,k) = ( azt(i,k+1) - azt(i,k) ) * gr%invrs_dzm(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function gradzt_2D
  
  !=============================================================================
  ! Wrapped in interface ddzt
  function gradzt_1D( gr, azt )

    ! Description:
    !  2D version of gradzt
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      gradzt_1D    ! Vertical derivative of azt    [units vary / m]

    ! Local Variables

    ! Input Variable
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      azt_col    ! Variable on thermodynamic grid levels    [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      gradzt_1D_col    ! Vertical derivative of azt    [units vary / m]

    azt_col(1,:) = azt
    
    gradzt_1D_col = gradzt_2D(gr%nz, 1, gr, azt_col)
    
    gradzt_1D = gradzt_1D_col(1,:)

    return

  end function gradzt_1D

  !=============================================================================
  function flip( x, xdim )

    ! Description:
    !   Flips a single dimension array (i.e. a vector), so the first element
    !   becomes the last and vice versa for the whole column.  This is a
    !   necessary part of the code because BUGSrad and CLUBB store altitudes in
    !   reverse order.
    !
    ! References:
    !   None
    !-------------------------------------------------------------------------
    
    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    integer, intent(in) :: xdim

    real(kind = dp), dimension(xdim), intent(in) :: x

    ! Output Variables
    real(kind = dp), dimension(xdim) :: flip

    ! Local Variables
    real(kind = dp), dimension(xdim) :: tmp

    integer :: indx


    ! Get rid of an annoying compiler warning.
    indx = 1
    indx = indx

    forall ( indx = 1 : xdim )
       tmp(indx) = x((xdim+1) - (indx))
    end forall

    flip = tmp


    return

  end function flip

!===============================================================================

end module grid_class
