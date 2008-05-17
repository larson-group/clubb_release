!------------------------------------------------------------------------
! $Id: grid.F,v 1.18 2008-05-17 21:33:53 griffinb Exp $
!===============================================================================
        module grid_class

! Description:
!
! Definition of a grid class and associated functions
!
! The grid specification is as follows:
!
!     +                ================== zm(nnzp) =========GP=======
!     |
!     |
! 1/dzt(nnzp)  +       ------------------ zt(nnzp) ---------GP-------
!     |        |
!     |        |
!     +  1/dzm(nnzp-1) ================== zm(nnzp-1) ================
!              |
!              |
!              +       ------------------ zt(nnzp-1) ----------------
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
! level k; the variable zt(k) stands for the thermodynamic level altitude
! at thermodynamic level k; the variable dzt(k) is the inverse distance 
! between momentum levels (over a central thermodynamic level k); and
! the variable dzm(k) is the inverse distance between thermodynamic 
! levels (over a central momentum level k).
!
! The grid setup is compatible with a stretched (unevely-spaced) grid.
! Thus, the distance between successive grid levels may not always be
! constant.
!
! NOTE:  Any future code written for use in the CLUBB parameterization
!        should use interpolation formulas consistent with a stretched
!        grid.  The simplest way to do so is to call the appropriate
!        interpolation function from this module.  Interpolations 
!        should not be handled in the form of:
!        ( var_zm(k) + var_zm(k-1) ) / 2; nor in the form of:
!        0.5*( var_zt(k+1) + var_zt(k) ).  Rather, all explicit 
!        interpolations should call zt2zm or zm2zt; while interpolations
!        for a variable being solved for implicitly in the code should 
!        use gr%weights_zt2zm (which refers to interp_weights_zt2zm_imp),
!        or gr%weights_zm2zt (which refers to interp_weights_zm2zt_imp).
!
! Momentum level 1 is placed at altitude zm_init, which is usually at
! the surface.  However, zm_init can be at any level defined by the 
! user.
!
! GP indicates ghost points. Variables located at those levels are
! not prognosed, but only used for boundary conditions.
!
! Chris Golaz, 7/17/99
! modified 9/10/99

!  References:

!  Section 3c, p. 3548 /Numerical discretization/ of:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!     Method and Model Description'' Golaz, et al. (2002)
!     JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

           implicit none

           public :: gr, grid, zt2zm, interp_weights_zt2zm_imp, zm2zt,
     .               interp_weights_zm2zt_imp, ddzm, ddzt,
     .               gridsetup, read_grid_heights

           private :: interpolated_azm, interpolated_azmk,
     .                interpolated_azmk_imp, interpolated_azt,
     .                interpolated_aztk, interpolated_aztk_imp,
     .                gradzm, gradzt

           private ! Default Scoping

           type grid
              integer :: nnzp
!   Note: Fortran 90/95 prevent an allocatable array from appearing
!   within a derived type.  However, a pointer can be used in the same
!   manner as an allocatable array, as we have done here (the grid 
!   pointers are always allocated rather than assigned and nullified 
!   like real pointers).
              real, pointer, dimension(:) :: zm, zt
              real, pointer, dimension(:) :: dzm, dzt
              real, pointer, dimension(:,:) :: weights_zm2zt,
     .                                         weights_zt2zm
           end type grid

!   The grid is defined here so that it is common throughout the module. 
!   The implication is that only one grid can be defined !

           type (grid) gr  

!   Modification for using HOC in a host model (i.e. one grid per column)
!$omp   threadprivate(gr)

           ! Interfaces provided for function overloading

           ! Interpolation functions
           interface zt2zm
              module procedure interpolated_azm, interpolated_azmk
           end interface

           interface interp_weights_zt2zm_imp
              module procedure interpolated_azmk_imp
           end interface

           interface zm2zt
              module procedure interpolated_azt, interpolated_aztk
           end interface

           interface interp_weights_zm2zt_imp
              module procedure interpolated_aztk_imp
           end interface

           ! Vertical derivative functions
           interface ddzm
              module procedure gradzm
           end interface

           interface ddzt
              module procedure gradzt
           end interface

           contains

!===============================================================================
           subroutine gridsetup( nnzp, implemented, grid_type, 
     .                           deltaz, zm_init, momentum_heights, 
     .                           thermodynamic_heights )

!          Description:
!          Grid Constructor
!
!          This subroutine sets up the CLUBB vertical grid.
!
!-----------------------------------------------------------------------

           use constants, only: fstderr

           implicit none

           ! Constant parameters
           ! Issue a warning if nnzp exceeds this number.
           integer, parameter ::
     .     NWARNING = 250  

           integer, parameter ::
     .     t_above = 1,   ! Upper thermodynamic level index (gr%weights_zt2zm).
     .     t_below = 2,   ! Lower thermodynamic level index (gr%weights_zt2zm).
     .     m_above = 1,   ! Upper momentum level index (gr%weights_zm2zt).
     .     m_below = 2    ! Lower momentum level index (gr%weights_zm2zt).

           ! Input Variables
           integer, intent(in) :: 
     .     nnzp     ! Number of vertical levels in grid      [#]

           ! Flag to see if CLUBB is running on it's own, 
           ! or if it's implemented as part of a host model.
           logical, intent(in) :: implemented

           ! If CLUBB is running on it's own, this option determines
           ! if it is using:
           ! 1) an evenly-spaced grid,
           ! 2) a stretched (unevenly-spaced) grid on the 
           !    thermodynamic grid levels, or
           ! 3) a stretched (unevenly-spaced) grid on the
           !    momentum grid levels.
           integer, intent(in) :: grid_type
         
           ! If the CLUBB model is running by itself, and is using an
           ! evenly-spaced grid (grid_type = 1), it needs the vertical 
           ! grid spacing and momentum-level starting altitude as input.
           real, intent(in) :: 
     .     deltaz,  ! Vertical grid spacing                  [m]
     .     zm_init  ! Initial grid altitude (momentum level) [m]

           ! If the CLUBB parameterization is implemented in a host model,
           ! it needs to use the host model's momentum level altitudes 
           ! and thermodynamic level altitudes.
           ! If the CLUBB model is running by itself, but is using a 
           ! stretched grid for thermodynamic levels (grid_type = 2),
           ! it needs to use the thermodynamic level altitudes as input.
           ! If the CLUBB model is running by itself, but is using a 
           ! stretched grid for momentum levels (grid_type = 3),
           ! it needs to use the momentum level altitudes as input.
           real, intent(in), dimension(nnzp) :: 
     .     momentum_heights,     ! Momentum level altitudes (input)      [m]
     .     thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

           ! Local Variables
           integer :: k, ierr ! Loop index and allocation stat


           ! Define the grid size

           if ( nnzp > NWARNING ) then
             write(fstderr,*) "Warning: running with vertical grid "//
     .         "which is larger than ", NWARNING
           endif

           gr%nnzp = nnzp

           allocate( gr%zm(1:nnzp), gr%zt(1:nnzp),
     .               gr%dzm(1:nnzp), gr%dzt(1:nnzp), 
     .               gr%weights_zm2zt(m_above:m_below,1:nnzp),
     .               gr%weights_zt2zm(t_above:t_below,1:nnzp),
     .               stat=ierr )

           if ( ierr /= 0 ) stop "Grid allocation failed."


           if ( .not. implemented ) then


              if ( grid_type == 1 ) then

                 ! Evenly-spaced grid.
                 ! Momentum level altitudes are defined based on the grid 
                 ! starting altitude, zm_init, the constant grid-spacing, 
                 ! deltaz, and the number of grid levels, gr%nnzp.

                 ! Define momentum level altitudes. The first momentum 
                 ! level is at altitude zm_init.
                 do k = 1, gr%nnzp, 1
                    gr%zm(k) = zm_init + (k-1) * deltaz
                 enddo

                 ! Define thermodynamic level altitudes.  Thermodynamic 
                 ! level altitudes are located at the central altitude 
                 ! levels, exactly halfway between momentum level 
                 ! altitudes.  The lowermost thermodynamic level is found
                 ! by taking 1/2 the altitude difference between the
                 ! bottom two momentum levels and subtracting that value
                 ! from the bottom momentum level.  The first 
                 ! thermodynamic level is below zm_init.
                 gr%zt(1) = zm_init - ( 0.5 * deltaz )
                 do k = 2, gr%nnzp, 1
                    gr%zt(k) = 0.5 * ( gr%zm(k) + gr%zm(k-1) )
                 enddo


              elseif ( grid_type == 2 ) then

                 ! Stretched (unevenly-spaced) grid:  stretched 
                 ! thermodynamic levels.  Thermodynamic levels are 
                 ! defined according to a stretched grid that is entered 
                 ! through the use an input file.  This is similar to a 
                 ! SAM-style stretched grid.

                 ! Define thermodynamic level altitudes.
                 do k = 1, gr%nnzp, 1
                    gr%zt(k) = thermodynamic_heights(k)
                 enddo

                 ! Define momentum level altitudes.  Momentum level 
                 ! altitudes are located at the central altitude levels,
                 ! exactly halfway between thermodynamic level altitudes.
                 ! The uppermost momentum level altitude is found by 
                 ! taking 1/2 the altitude difference between the top 
                 ! two thermodynamic levels and adding that value to the
                 ! top thermodynamic level.
                 do k = 1, gr%nnzp-1, 1
                    gr%zm(k) = 0.5 * ( gr%zt(k+1) + gr%zt(k) )
                 enddo
                 gr%zm(gr%nnzp) = gr%zt(gr%nnzp) + 
     .                0.5 * ( gr%zt(gr%nnzp) - gr%zt(gr%nnzp-1) )


              elseif ( grid_type == 3 ) then

                 ! Stretched (unevenly-spaced) grid:  stretched 
                 ! momentum levels.  Momentum levels are defined 
                 ! according to a stretched grid that is entered through
                 ! the use an input file.  This is similar to a WRF-style
                 ! stretched grid.

                 ! Define momentum level altitudes.
                 do k = 1, gr%nnzp, 1
                    gr%zm(k) = momentum_heights(k)
                 enddo

                 ! Define thermodynamic level altitudes.  Thermodynamic
                 ! level altitudes are located at the central altitude
                 ! levels, exactly halfway between momentum level 
                 ! altitudes.  The lowermost thermodynamic level is found
                 ! by taking 1/2 the altitude difference between the 
                 ! bottom two momentum levels and subtracting that value
                 ! from the bottom momentum level.
                 gr%zt(1) = gr%zm(1) - 0.5 * ( gr%zm(2) - gr%zm(1) )
                 do k = 2, gr%nnzp, 1
                    gr%zt(k) = 0.5 * ( gr%zm(k) + gr%zm(k-1) )
                 enddo


              else

                 ! Invalid grid type.
                 write(fstderr,*) "Invalid grid type: ", grid_type,
     .                            ".  Valid options are 1, 2, or 3."
                 stop


              endif


           else

              ! The CLUBB parameterization is implemented in a host 
              ! model.  Use the host model's momentum level altitudes
              ! and thermodynamic level altitudes to set up the CLUBB
              ! grid.

              ! Momentum level altitudes from host model.
              do k = 1, gr%nnzp, 1
                 gr%zm(k) = momentum_heights(k)
              enddo

              ! Thermodynamic level altitudes from host model after 
              ! possible grid-index adjustment for CLUBB interface.
              do k = 1, gr%nnzp, 1
                 gr%zt(k) = thermodynamic_heights(k)
              enddo


           endif


           ! Define dzm, which is the inverse spacing between 
           ! thermodynamic grid levels.
           do k=1,gr%nnzp-1
              gr%dzm(k) = 1. / ( gr%zt(k+1) - gr%zt(k) )
           enddo
           gr%dzm(gr%nnzp) = gr%dzm(gr%nnzp-1)


           ! Define dzt, which is the inverse spacing between 
           ! momentum grid levels.
           do k=2,gr%nnzp
              gr%dzt(k) = 1. / ( gr%zm(k) - gr%zm(k-1) )
           enddo
           gr%dzt(1) = gr%dzt(2)


           ! Interpolation Weights: zm grid to zt grid.
           ! The grid index (k) is the index of the level on the 
           ! thermodynamic (zt) grid.
           ! The result is the weights of the upper and lower momentum
           ! levels on the thermodynamic level.
           ! These weights are normally used in situations where a 
           ! momentum level variable is being solved for implicitly in 
           ! an equation and needs to be interpolated to the 
           ! thermodynamic grid levels.
           do k = 1, gr%nnzp, 1
              gr%weights_zm2zt(m_above:m_below,k) 
     .               = interp_weights_zm2zt_imp( k )
           enddo


           ! Interpolation Weights: zt grid to zm grid.
           ! The grid index (k) is the index of the level on the 
           ! momentum (zm) grid.
           ! The result is the weights of the upper and lower 
           ! thermodynamic levels on the momentum level.
           ! These weights are normally used in situations where a 
           ! thermodynamic level variable is being solved for implicitly
           ! in an equation and needs to be interpolated to the 
           ! momentum grid levels.
           do k = 1, gr%nnzp, 1
              gr%weights_zt2zm(t_above:t_below,k) 
     .               = interp_weights_zt2zm_imp( k )
           enddo


           return

           end subroutine gridsetup

!===============================================================================
           subroutine read_grid_heights( nzmax, grid_type, 
     .                                   zm_grid_fname, zt_grid_fname,
     .                                   momentum_heights,
     .                                   thermodynamic_heights )

!          Description:
!          This subroutine is used foremost in cases where the grid_type
!          corresponds with the stretched (unevenly-spaced) grid options
!          (either grid_type = 2 or grid_type = 3).  This subroutine 
!          reads in the values of the stretched grid altitude levels for
!          either the thermodynamic level grid or the momentum level 
!          grid.  This subroutine also handles basic error checking for
!          all three grid types.
!------------------------------------------------------------------------

           use constants, only: fstderr
           use file_functions, only: file_read_1d

           implicit none

           ! Input Variables.

           ! Declared number of vertical levels.
           integer, intent(in) ::
     .     nzmax

           ! If CLUBB is running on it's own, this option determines
           ! if it is using:
           ! 1) an evenly-spaced grid,
           ! 2) a stretched (unevenly-spaced) grid on the
           !    thermodynamic grid levels, or
           ! 3) a stretched (unevenly-spaced) grid on the
           !    momentum grid levels.
           integer, intent(in) ::
     .     grid_type

           character(len=50), intent(in) ::
     .     zm_grid_fname, ! Path and filename of file for momentum level altitudes
     .     zt_grid_fname  ! Path and filename of file for thermodynamic level altitudes

           ! Output Variables.

           ! If the CLUBB model is running by itself, but is using a
           ! stretched grid for thermodynamic levels (grid_type = 2),
           ! it needs to use the thermodynamic level altitudes as input.
           ! If the CLUBB model is running by itself, but is using a
           ! stretched grid for momentum levels (grid_type = 3),
           ! it needs to use the momentum level altitudes as input.
           real, dimension(nzmax), intent(out) ::
     .     momentum_heights,     ! Momentum level altitudes (file input)      [m]
     .     thermodynamic_heights ! Thermodynamic level altitudes (file input) [m]

           ! Local Variables.

           integer, parameter :: zt_unit = 25, ! Unit number for zt_grid_fname
     .                           zm_unit = 35  ! Unit number for zm_grid_fname

           integer :: zt_level_count, ! Number of altitudes found in zt_grid_fname
     .                zm_level_count  ! Number of altitudes found in zm_grid_fname

           integer :: input_status   ! Status of file being read: 
                                     ! > 0 ==> error reading file.
                                     ! = 0 ==> no error and more file to be read.
                                     ! < 0 ==> end of file indicator.

           ! Generic variable for storing file data while counting the number
           ! of file entries.
           real :: generic_input_item

           integer :: k   ! Loop index


           ! Declare the momentum level altitude array and the 
           ! thermodynamic level altitude array to be 0 until overwritten.
           momentum_heights(1:nzmax) = 0.0
           thermodynamic_heights(1:nzmax) = 0.0


           if ( grid_type == 1 ) then

              ! Evenly-spaced grid.
              ! Grid level altitudes are based on a constant distance 
              ! between them and a starting point for the bottom of the
              ! grid.

              ! As a way of error checking, make sure that there isn't any
              ! file entry for either momentum level altitudes or
              ! thermodynamic level altitudes.
              if ( zm_grid_fname /= '' ) then
                 write(fstderr,*)
     .              "An evenly-spaced grid has been selected. "
     .              // " Please reset zm_grid_fname to ''."
                 stop
              endif
              if ( zt_grid_fname /= '' ) then
                 write(fstderr,*)
     .              "An evenly-spaced grid has been selected. "
     .              // " Please reset zt_grid_fname to ''."
                 stop
              endif


           elseif ( grid_type == 2 ) then

              ! Stretched (unevenly-spaced) grid:  stretched
              ! thermodynamic levels.  Thermodynamic levels are
              ! defined according to a stretched grid that is entered
              ! through the use an input file.  This is similar to a
              ! SAM-style stretched grid.

              ! As a way of error checking, make sure that there isn't any
              ! file entry for momentum level altitudes.
              if ( zm_grid_fname /= '' ) then
                 write(fstderr,*)
     .              "Thermodynamic level altitudes have been selected "
     .              // "for use in a stretched (unevenly-spaced) grid. "
     .              // " Please reset zm_grid_fname to ''."
                 stop
              endif

              ! Open the file zt_grid_fname.
              open( unit=zt_unit, file=zt_grid_fname, 
     .              status='old', action='read' )

              ! Find the number of thermodynamic level altitudes 
              ! listed in file zt_grid_fname.
              zt_level_count = 0
              do
                 read( unit=zt_unit, fmt=*, iostat=input_status ) 
     .              generic_input_item
                 if ( input_status < 0 ) exit   ! end of file indicator
                 if ( input_status > 0 ) stop   ! error reading input
     .                "Error reading thermodynamic level input file."
                 zt_level_count = zt_level_count + 1
              enddo

              ! Close the file zt_grid_fname.
              close( unit=zt_unit )

              ! Check that the number of thermodynamic grid altitudes
              ! in the input file matches the declared number of CLUBB
              ! grid levels (nzmax).
              if ( zt_level_count /= nzmax ) then
                 write(fstderr,*) 
     .              "The number of thermodynamic grid altitudes "
     .              // "listed in file " // trim(zt_grid_fname) 
     .              // " does not match the number of CLUBB grid "
     .              // "levels specified in the model.in file."
                 write(fstderr,*)
     .              "Number of thermodynamic grid altitudes listed:  ",
     .              zt_level_count
                 write(fstderr,*)
     .              "Number of CLUBB grid levels specified:  ", nzmax
                 stop
              endif

              ! Read the thermodynamic level altitudes from zt_grid_fname.
              call file_read_1d( zt_unit, zt_grid_fname, nzmax, 1, 
     .                           thermodynamic_heights )

              ! Check that each thermodynamic level altitude increases
              ! in height as the thermodynamic level grid index increases.
              do k = 2, nzmax, 1
                 if ( thermodynamic_heights(k) 
     .                <= thermodynamic_heights(k-1) ) then
                    write(fstderr,*) 
     .                 "The declared thermodynamic level grid "
     .                 // "altitudes are not increasing in height "
     .                 // "as grid level index increases."
                    write(fstderr,*)
     .                 "Grid index:  ", k-1, ";",
     .                 "  Thermodynamic level altitude:  ",
     .                 thermodynamic_heights(k-1)
                    write(fstderr,*)
     .                 "Grid index:  ", k, ";",
     .                 "  Thermodynamic level altitude:  ",
     .                 thermodynamic_heights(k)
                    stop
                 endif
              enddo


           elseif ( grid_type == 3 ) then

              ! Stretched (unevenly-spaced) grid:  stretched
              ! momentum levels.  Momentum levels are defined
              ! according to a stretched grid that is entered through
              ! the use an input file.  This is similar to a WRF-style
              ! stretched grid.

              ! As a way of error checking, make sure that there isn't any
              ! file entry for thermodynamic level altitudes.
              if ( zt_grid_fname /= '' ) then
                 write(fstderr,*)
     .              "Momentum level altitudes have been selected "
     .              // "for use in a stretched (unevenly-spaced) grid. "
     .              // " Please reset zt_grid_fname to ''."
                 stop
              endif

              ! Open the file zm_grid_fname.
              open( unit=zm_unit, file=zm_grid_fname, 
     .              status='old', action='read' )

              ! Find the number of momentum level altitudes
              ! listed in file zm_grid_fname.
              zm_level_count = 0
              do
                 read( unit=zm_unit, fmt=*, iostat=input_status )
     .              generic_input_item
                 if ( input_status < 0 ) exit   ! end of file indicator
                 if ( input_status > 0 ) stop   ! error reading input
     .                "Error reading momentum level input file."
                 zm_level_count = zm_level_count + 1
              enddo

              ! Close the file zm_grid_fname.
              close( unit=zm_unit )

              ! Check that the number of momentum grid altitudes
              ! in the input file matches the declared number of CLUBB
              ! grid levels (nzmax).
              if ( zm_level_count /= nzmax ) then
                 write(fstderr,*)
     .              "The number of momentum grid altitudes "
     .              // "listed in file " // trim(zm_grid_fname)
     .              // " does not match the number of CLUBB grid "
     .              // "levels specified in the model.in file."
                 write(fstderr,*)
     .              "Number of momentum grid altitudes listed:  ",
     .              zm_level_count
                 write(fstderr,*)
     .              "Number of CLUBB grid levels specified:  ", nzmax
                 stop
              endif

              ! Read the momentum level altitudes from zm_grid_fname.
              call file_read_1d( zm_unit, zm_grid_fname, nzmax, 1,
     .                           momentum_heights )

              ! Check that each momentum level altitude increases
              ! in height as the momentum level grid index increases.
              do k = 2, nzmax, 1
                 if ( momentum_heights(k) 
     .                <= momentum_heights(k-1) ) then
                    write(fstderr,*) 
     .                 "The declared momentum level grid "
     .                 // "altitudes are not increasing in height "
     .                 // "as grid level index increases."
                    write(fstderr,*)
     .                 "Grid index:  ", k-1, ";",
     .                 "  Momentum level altitude:  ",
     .                 momentum_heights(k-1)
                    write(fstderr,*)
     .                 "Grid index:  ", k, ";",
     .                 "  Momentum level altitude:  ",
     .                 momentum_heights(k)
                    stop
                 endif
              enddo


           endif


           return

           end subroutine read_grid_heights

!===============================================================================
           pure function interpolated_azm( azt )

!          Description:
!          Function to interpolate a variable located on the 
!          thermodynamic grid levels (azt) to the momentum grid 
!          levels (azm).  This function inputs the entire azt array and
!          outputs the results as an azm array.  The formulation used is
!          compatible with a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azt

           ! Return Variable
           real, dimension(gr%nnzp) :: interpolated_azm
        
           ! Local Variable
           integer :: k

           ! Do the actual interpolation.
           ! Use linear interpolation.
           do k = 1, gr%nnzp-1, 1
              interpolated_azm(k) = 
     .             ( ( azt(k+1)-azt(k) ) / ( gr%zt(k+1)-gr%zt(k) ) )
     .              * ( gr%zm(k)-gr%zt(k) ) + azt(k)
           enddo
!           ! Set the value of azm at level gr%nnzp (the uppermost level
!           ! in the model) to the value of azt at level gr%nnzp.
!           interpolated_azm(gr%nnzp) = azt(gr%nnzp)
           ! Use a linear extension based on the values of azt at levels
           ! gr%nnzp and gr%nnzp-1 to find the value of azm at level 
           ! gr%nnzp (the uppermost level in the model).
           interpolated_azm(gr%nnzp) = 
     .          ( ( azt(gr%nnzp)-azt(gr%nnzp-1) )
     .           / ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) ) )
     .           * ( gr%zm(gr%nnzp)-gr%zt(gr%nnzp) ) + azt(gr%nnzp)

           return

           end function interpolated_azm

!===============================================================================
           pure function interpolated_azmk( azt, k )

!          Description:
!          Function to interpolate a variable located on the 
!          thermodynamic grid levels (azt) to the momentum grid 
!          levels (azm).  This function outputs the value of azm at a 
!          single grid level (k) after interpolating using values of azt
!          at two grid levels.  The formulation used is compatible with 
!          a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variables
           real, intent(in), dimension(gr%nnzp) :: azt

           integer, intent(in) :: k

           ! Return Variable
           real :: interpolated_azmk

           ! Do the actual interpolation.
           ! Use linear interpolation.
           if ( k /= gr%nnzp ) then
              interpolated_azmk = 
     .             ( ( azt(k+1)-azt(k) ) / ( gr%zt(k+1)-gr%zt(k) ) )
     .              * ( gr%zm(k)-gr%zt(k) ) + azt(k)
           else
!              ! Set the value of azm at level gr%nnzp (the uppermost 
!              ! level in the model) to the value of azt at level gr%nnzp.
!              interpolated_azmk = azt(gr%nnzp)
              ! Use a linear extension based on the values of azt at 
              ! levels gr%nnzp and gr%nnzp-1 to find the value of azm at
              ! level gr%nnzp (the uppermost level in the model).
              interpolated_azmk = 
     .             ( ( azt(gr%nnzp)-azt(gr%nnzp-1) )
     .              / ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) ) )
     .              * ( gr%zm(gr%nnzp)-gr%zt(gr%nnzp) ) + azt(gr%nnzp)
           endif

           return

           end function interpolated_azmk
        
!===============================================================================
           pure function interpolated_azmk_imp( m_lev )
     .     result( azt_weight )

!          Description:
!          Function used to help in an interpolation of a variable 
!          located on the thermodynamic grid levels (azt) to the 
!          momentum grid levels (azm).   This function computes a 
!          weighting factor for both the upper thermodynamic level (k+1)
!          and the lower thermodynamic level (k) on the central momentum
!          level (k).  For the uppermost momentum grid level (k=gr%nnzp),
!          a weighting factor for both the thermodynamic level at gr%nnzp
!          and the thermodynamic level at gr%nnzp-1 are computed based on
!          the use of a linear extension.  This function outputs the 
!          weighting factors at a single grid level (k).  This function 
!          is used in situations where the variable to be interpolated is
!          being treated implicitly in an equation.  The formulation used
!          is compatible with a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Constant parameters
           integer, parameter ::
     .     t_above = 1,   ! Upper thermodynamic level.
     .     t_below = 2    ! Lower thermodynamic level.

           ! Input
           integer, intent(in) :: m_lev  ! Momentum level index

           ! Output
           real, dimension(2) :: azt_weight  ! Weights of the
                                             ! thermodynamic levels.

           ! Local Variables
           real :: factor
           integer :: k

           ! Compute the weighting factors at momentum level k.
           k = m_lev

           if ( k /= gr%nnzp ) then
              ! At most levels, the momentum level is found
              ! in-between two thermodynamic levels.  Linear
              ! interpolation is used.
              factor = 
     .             ( gr%zm(k)-gr%zt(k) ) / ( gr%zt(k+1)-gr%zt(k) )
           else
              ! The top model level (gr%nnzp) is formulated
              ! differently because the top momentum level is
              ! above the top thermodynamic level.  A linear
              ! extension is required, rather than linear
              ! interpolation.
              factor = 
     .             ( gr%zm(gr%nnzp)-gr%zt(gr%nnzp-1) ) 
     .             / ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
           endif

           ! Weight of upper thermodynamic level on momentum level.
           azt_weight(t_above) = factor
           ! Weight of lower thermodynamic level on momentum level.
           azt_weight(t_below) = 1.0 - factor

           return

           end function interpolated_azmk_imp

!===============================================================================
           pure function interpolated_azt( azm )

!          Description:
!          Function to interpolate a variable located on the
!          momentum grid levels (azm) to the thermodynamic grid
!          levels (azt).  This function inputs the entire azm array and
!          outputs the results as an azt array.  The formulation used is
!          compatible with a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azm

           ! Output Variable
           real, dimension(gr%nnzp) :: interpolated_azt
        
           ! Local Variable
           integer :: k  ! Index

           ! Do actual interpolation.
           ! Use linear interpolation.
           do k = gr%nnzp, 2, -1
              interpolated_azt(k) =
     .             ( ( azm(k)-azm(k-1) ) / ( gr%zm(k)-gr%zm(k-1) ) )
     .              * ( gr%zt(k)-gr%zm(k-1) ) + azm(k-1)
           enddo
!           ! Set the value of azt at level 1 (the lowermost level in the 
!           ! model) to the value of azm at level 1.
!           interpolated_azt(1) = azm(1)
           ! Use a linear extension based on the values of azm at levels
           ! 1 and 2 to find the value of azt at level 1 (the lowermost 
           ! level in the model).
           interpolated_azt(1) =
     .          ( ( azm(2)-azm(1) ) / ( gr%zm(2)-gr%zm(1) ) )
     .           * ( gr%zt(1)-gr%zm(1) ) + azm(1)

           return

           end function interpolated_azt

!===============================================================================
           pure function interpolated_aztk( azm, k )

!          Description:
!          Function to interpolate a variable located on the
!          momentum grid levels (azm) to the thermodynamic grid
!          levels (azt).  This function outputs the value of azt at a
!          single grid level (k) after interpolating using values of azm
!          at two grid levels.  The formulation used is compatible with 
!          a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variables
           real, intent(in), dimension(gr%nnzp) :: azm

           integer, intent(in) :: k

           ! Return Variables
           real :: interpolated_aztk
        
           ! Do actual interpolation.
           ! Use linear interpolation.
           if ( k /= 1 ) then
              interpolated_aztk =
     .             ( ( azm(k)-azm(k-1) ) / ( gr%zm(k)-gr%zm(k-1) ) )
     .              * ( gr%zt(k)-gr%zm(k-1) ) + azm(k-1)
           else
!              ! Set the value of azt at level 1 (the lowermost level in 
!              ! the model) to the value of azm at level 1.
!              interpolated_aztk = azm(1)
              ! Use a linear extension based on the values of azm at 
              ! levels 1 and 2 to find the value of azt at level 1 (the 
              ! lowermost level in the model).
              interpolated_aztk =
     .             ( ( azm(2)-azm(1) ) / ( gr%zm(2)-gr%zm(1) ) )
     .              * ( gr%zt(1)-gr%zm(1) ) + azm(1)
           endif

           return

           end function interpolated_aztk

!===============================================================================
           pure function interpolated_aztk_imp( t_lev )
     .     result( azm_weight )

!          Description:
!          Function used to help in an interpolation of a variable
!          located on the momentum grid levels (azm) to the 
!          thermodynamic grid levels (azt).   This function computes a 
!          weighting factor for both the upper momentum level (k) and 
!          the lower momentum level (k-1) on the central thermodynamic 
!          level (k).  For the lowermost thermodynamic grid level (k=1), 
!          a weighting factor for both the momentum level at 1 and the 
!          momentum level at 2 are computed based on the use of a linear 
!          extension.  This function outputs the weighting factors at a 
!          single grid level (k).  This function is used in situations 
!          where the variable to be interpolated is being treated 
!          implicitly in an equation.  The formulation used is compatible
!          with a stretched (unevenly-spaced) grid.
!-----------------------------------------------------------------------

           implicit none

           ! Constant parameters
           integer, parameter ::
     .     m_above = 1,   ! Upper momentum level.
     .     m_below = 2    ! Lower momentum level.

           ! Input
           integer, intent(in) :: t_lev  ! Thermodynamic level index

           ! Output
           real, dimension(2) :: azm_weight  ! Weights of the
                                             ! momentum levels.

           ! Local Variables
           real :: factor
           integer :: k

           ! Compute the weighting factors at thermodynamic level k.
           k = t_lev

           if ( k /= 1 ) then
              ! At most levels, the thermodynamic level is
              ! found in-between two momentum levels.  Linear
              ! interpolation is used.
              factor = 
     .             ( gr%zt(k)-gr%zm(k-1) ) / ( gr%zm(k)-gr%zm(k-1) )
           else
              ! The bottom model level (1) is formulated
              ! differently because the bottom thermodynamic level
              ! is below the bottom momentum level.  A linear
              ! extension is required, rather than linear
              ! interpolation.
              factor = 
     .             ( gr%zt(1)-gr%zm(1) ) / ( gr%zm(2)-gr%zm(1) )
           endif

           ! Weight of upper momentum level on thermodynamic level.
           azm_weight(m_above) = factor
           ! Weight of lower momentum level on thermodynamic level.
           azm_weight(m_below) = 1.0 - factor

           return

           end function interpolated_aztk_imp

!===============================================================================
           pure function gradzm( azm )

!          Description:
!          Function to compute the vertical derivative of a variable (azm)
!          located on the momentum grid.  The results are returned in an 
!          array defined on the thermodynamic grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azm

           ! Return Variable
           real, dimension(gr%nnzp) :: gradzm

           ! Local Variable
           integer :: k

           ! Compute vertical derivatives.
           do k = gr%nnzp, 2, -1
              ! Take derivative of momentum-level variable azm
              ! over the central thermodynamic level (k).
              gradzm(k) = ( azm(k) - azm(k-1) ) * gr%dzt(k)
           enddo
!           ! Thermodynamic level 1 is located below momentum level 1,
!           ! so there is not enough information to calculate the 
!           ! derivative over thermodynamic level 1.  Thus, the value of 
!           ! the derivative at thermodynamic level 1 is set equal to 0. 
!           ! This formulation is consistent with setting the value of 
!           ! the variable azm below the model grid to the value of the
!           ! variable azm at the lowest grid level.
!           gradzm(1) = 0.
           ! Thermodynamic level 1 is located below momentum level 1,
           ! so there is not enough information to calculate the 
           ! derivative over thermodynamic level 1.  Thus, the value of 
           ! the derivative at thermodynamic level 1 is set equal to the 
           ! value of the derivative at thermodynamic level 2.  This 
           ! formulation is consistent with using a linear extension to 
           ! find the values of the variable azm below the model grid.
           gradzm(1) = gradzm(2)

           return

           end function gradzm

!===============================================================================
           pure function gradzt( azt )

!          Description:
!          Function to compute the vertical derivative of a variable (azt)
!          located on the thermodynamic grid.  The results are returned in
!          an array defined on the momentum grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azt

           ! Output Variable
           real, dimension(gr%nnzp) :: gradzt

           ! Local Variable
           integer :: k

           ! Compute vertical derivative.
           do k = 1, gr%nnzp-1, 1
              ! Take derivative of thermodynamic-level variable azt
              ! over the central momentum level (k).
              gradzt(k) = ( azt(k+1) - azt(k) ) * gr%dzm(k)
           enddo
!           ! Momentum level gr%nnzp is located above thermodynamic level
!           ! gr%nnzp, so there is not enough information to calculate the
!           ! derivative over momentum level gr%nnzp.  Thus, the value of
!           ! the derivative at momentum level gr%nnzp is set equal to 0.
!           ! This formulation is consistent with setting the value of
!           ! the variable azt above the model grid to the value of the
!           ! variable azt at the highest grid level.
!           gradzt(gr%nnzp) = 0.
           ! Momentum level gr%nnzp is located above thermodynamic level
           ! gr%nnzp, so there is not enough information to calculate the
           ! derivative over momentum level gr%nnzp.  Thus, the value of
           ! the derivative at momentum level gr%nnzp is set equal to the
           ! value of the derivative at momentum level gr%nnzp-1.  This
           ! formulation is consistent with using a linear extension to
           ! find the values of the variable azt above the model grid.
           gradzt(gr%nnzp) = gradzt(gr%nnzp-1)

           return

           end function gradzt

!===============================================================================

        end module grid_class
