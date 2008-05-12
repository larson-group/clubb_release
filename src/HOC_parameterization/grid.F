!------------------------------------------------------------------------
! $Id: grid.F,v 1.13 2008-05-12 18:29:26 dschanen Exp $

        module grid_class

! Description:
!
! Definition of a grid class and associated functions
!
! The grid specification is as follows:
!
!   +              ================== zm(nnzp) =========GP=======
!   |
!   |
! dzt(nnzp)  +     ------------------ zt(nnzp) ---------GP-------
!   |        |
!   |        |
!   +  dzm(nnzp-1) ================== zm(nnzp-1) ================
!            |
!            |
!            +     ------------------ zt(nnzp-1) ----------------
! 
!                                       .
!                                       .
!                                       .
!                                       .
! 
!   +              ================== zm(2) =====================
!   |
!   |
! dzt(2)     +     ------------------ zt(2) ---------------------
!   |        |
!   |        |
!   +      dzm(1)  ================== zm(1) ============GP=======  surface
!            |     //////////////////////////////////////////////
!            |
!            +     ------------------ zt(1) ------------GP-------
!
!
! GP indicates ghost points. Variables located at those levels are
! not prognosed, but only used for boundary conditions.
!
! Chris Golaz, 7/17/99
! modified 9/10/99
!
! modified 2/28/01 to get rid of UNDEFINED.
! Experience has shown that the 'UNDEFINED' idea was not shuch a good idea
! after all.

!  References:

!  Section 3c, p. 3548 /Numerical discretization/ of:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!     Method and Model Description'' Golaz, et al. (2002)
!     JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

           implicit none

           public :: gr, grid, zt2zm, interp_weights_zt2zm_imp, zm2zt,
     .               interp_weights_zm2zt_imp, ddzm, ddzt,
     .               gridsetup

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

!-----------------------------------------------------------------------
           subroutine gridsetup( nnzp, deltaz, zm_init, implemented,
     .                           host_momentum, host_thermodynamic )

!          Description:
!          Grid Constructor
!
!          This subroutine sets up a uniformly spaced vertical grid 
!          similar to the RAMS grid.
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
     .     nnzp ! Number of vertical levels in grid     [#]

           real, intent(in) :: 
     .     deltaz,  ! Vertical grid spacing             [m]
     .     zm_init  ! Initial grid point                [m]

           ! Flag to see if HOC is running on it's own or it it's part
           ! of a larger model.
           logical, intent(in) :: implemented

           ! If it's part of a larger model, it needs to use that models
           ! grid heights.
           real, intent(in), dimension(nnzp) :: 
     .     host_momentum,     ! Host momentum grid              [m]
     .     host_thermodynamic ! Host thermodynamic grid         [m]

           ! Local Variables
           integer :: k, ierr ! Loop index and allocation stat


           ! Define the grid size

           if ( nnzp > NWARNING ) then
             write(fstderr,*) "Warning: running with vertical grid "//
     .         "which is larger than ", NWARNING
           end if

           gr%nnzp = nnzp

           allocate( gr%zm(1:nnzp), gr%zt(1:nnzp),
     .               gr%dzm(1:nnzp), gr%dzt(1:nnzp), 
     .               gr%weights_zm2zt(m_above:m_below,1:nnzp),
     .               gr%weights_zt2zm(t_above:t_below,1:nnzp),
     .               stat=ierr )

           if ( ierr /= 0 ) stop "Grid allocation failed."

           IF ( .NOT. implemented ) THEN

              ! Original code

              ! Define momentum grid points. The first momentum level is 
              ! at the surface
              DO k = 1, gr%nnzp, 1
                 gr%zm(k) = zm_init + (k-1) * deltaz
              END DO

             ! Define thermodynamic grid. For the time being, themodynamic 
             ! grid points are located in the center point between momentum 
             ! points. The first thermodynmic point is below the surface.
              gr%zt(1) = zm_init - ( 0.5 * deltaz )
              DO k = 2, gr%nnzp, 1
                 gr%zt(k) = 0.5 * ( gr%zm(k) + gr%zm(k-1) )
              END DO

           ELSE

              ! Momentum levels from host model.
              DO k = 1, gr%nnzp, 1
                 gr%zm(k) = host_momentum(k)
              END DO

              ! Thermodynamic levels from host model after possible 
              ! grid-index adjustment for HOC interface.
              DO k = 1, gr%nnzp, 1
                 gr%zt(k) = host_thermodynamic(k)
              END DO

           END IF

           ! Define dzm. dzm is the inverse spacing between 
           ! thermodynamic points

           do k=1,gr%nnzp-1
              gr%dzm(k) = 1. / ( gr%zt(k+1) - gr%zt(k) )
           enddo
           gr%dzm(gr%nnzp) = gr%dzm(gr%nnzp-1)

           ! Define dzt. dzt is the inverse spacing between momentum points

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

!-----------------------------------------------------------------------
           pure function interpolated_azm( azt )

!          Description:
!          Function to interpolate a variable located on the 
!          thermodynamic grid (azt) to a variable on the momentum 
!          grid (azm)

!          Notes:
!          Modified by Brian for a non-uniform grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azt

           ! Return Variable
           real, dimension(gr%nnzp) :: interpolated_azm
        
           ! Local Variable
           integer :: k

           ! Do the actual interpolation 

           do k=1, gr%nnzp-1
!              interpolated_azm(k) = 0.5 * ( azt(k) + azt(k+1) )
              interpolated_azm(k) = 
     .                ( ( azt(k+1)-azt(k) )/( gr%zt(k+1)-gr%zt(k) ) )
     .               * ( gr%zm(k)-gr%zt(k) ) + azt(k)
           end do
!           interpolated_azm(gr%nnzp) = UNDEFINED
!           interpolated_azm(gr%nnzp) = azt(gr%nnzp)
           interpolated_azm(gr%nnzp) = 
     .                ( ( azt(gr%nnzp)-azt(gr%nnzp-1) )
     .                 /( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) ) )
     .               * ( gr%zm(gr%nnzp)-gr%zt(gr%nnzp) ) + azt(gr%nnzp)

           return
           end function interpolated_azm

!-----------------------------------------------------------------------
           pure function interpolated_azmk( azt, k )

!          Description:
!          Function to interpolate a variable located on the
!          thermodynamic grid (azt) to a variable on the momentum
!          grid (azm). The function interpolates only for level k

!          Notes:
!          Modified by Brian for a non-uniform grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azt

           integer, intent(in) :: k

           ! Return Variable
           real :: interpolated_azmk

           ! Do the actual interpolation

           if ( k /= gr%nnzp ) then
!              interpolated_azmk = 0.5 * ( azt(k) + azt(k+1) )
              interpolated_azmk = 
     .                ( ( azt(k+1)-azt(k) )/( gr%zt(k+1)-gr%zt(k) ) )
     .               * ( gr%zm(k)-gr%zt(k) ) + azt(k)
           else
!              interpolated_azmk = UNDEFINED
!              interpolated_azmk = azt(gr%nnzp)
              interpolated_azmk = 
     .                ( ( azt(gr%nnzp)-azt(gr%nnzp-1) )
     .                 /( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) ) )
     .               * ( gr%zm(gr%nnzp)-gr%zt(gr%nnzp) ) + azt(gr%nnzp)
           end if

           return
           end function interpolated_azmk
        
!-----------------------------------------------------------------------
           pure function interpolated_azmk_imp( m_lev )
     .     result( azt_weight )

!          Description:
!          Function used to help in an interpolation of a variable 
!          located on the thermodynamic grid (azt) to a variable on the 
!          momentum grid (azm).   This function is used in 
!          situations where the variable to be interpolated is being 
!          treated implicitly in an equation.  This function computes a
!          weighting factor.  This function interpolates only for 
!          level k.

!          Notes:
!          Created by Brian for a non-uniform grid.
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

           k = m_lev

           ! Do actual interpolation

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
     .             /( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
           endif

           ! Weight of upper thermodynamic level on momentum level.
           azt_weight(t_above) = factor
           ! Weight of lower thermodynamic level on momentum level.
           azt_weight(t_below) = 1.0 - factor

           return
           end function interpolated_azmk_imp

!-----------------------------------------------------------------------
           pure function interpolated_azt( azm )

!          Description:
!          Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt)

!          Notes:
!          Modified by Brian for a non-uniform grid.
!-----------------------------------------------------------------------

           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azm

           ! Output Variable
           real, dimension(gr%nnzp) :: interpolated_azt
        
           ! Local Variable
           integer :: k  ! Index

           ! Do actual interpolation 

           do k=gr%nnzp,2,-1
!              interpolated_azt(k) = 0.5 * ( azm(k-1) + azm(k) )
              interpolated_azt(k) =
     .                ( ( azm(k)-azm(k-1) )/( gr%zm(k)-gr%zm(k-1) ) )
     .               * ( gr%zt(k)-gr%zm(k-1) ) + azm(k-1)
           end do

!           interpolated_azt(1) = UNDEFINED
!           interpolated_azt(1) = azm(1)
           interpolated_azt(1) =
     .                ( ( azm(2)-azm(1) )/( gr%zm(2)-gr%zm(1) ) )
     .               * ( gr%zt(1)-gr%zm(1) ) + azm(1)

           return
           end function interpolated_azt

!-----------------------------------------------------------------------
           pure function interpolated_aztk( azm, k )

!          Description:
!          Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt). This function interpolates only for level k.

!          Notes:
!          Modified by Brian for a non-uniform grid.
!-----------------------------------------------------------------------
           implicit none

           ! Input Variables
           real, intent(in), dimension(gr%nnzp) :: azm

           integer, intent(in) :: k

           ! Return Variables
           real :: interpolated_aztk
        
           ! Do actual interpolation 

           if ( k /= 1 ) then
!              interpolated_aztk = 0.5 * ( azm(k-1) + azm(k) )
              interpolated_aztk =
     .                ( ( azm(k)-azm(k-1) )/( gr%zm(k)-gr%zm(k-1) ) )
     .               * ( gr%zt(k)-gr%zm(k-1) ) + azm(k-1)
           else
!              interpolated_aztk = UNDEFINED
!              interpolated_aztk = azm(1)
              interpolated_aztk =
     .                ( ( azm(2)-azm(1) )/( gr%zm(2)-gr%zm(1) ) )
     .               * ( gr%zt(1)-gr%zm(1) ) + azm(1)
           end if

           return
           end function interpolated_aztk

!-----------------------------------------------------------------------
           pure function interpolated_aztk_imp( t_lev )
     .     result( azm_weight )

!          Description:
!          Function used to help in an interpolation of a variable 
!          located on the momentum grid (azm) to a variable on the 
!          thermodynamic grid (azt).   This function is used in 
!          situations where the variable to be interpolated is being 
!          treated implicitly in an equation.  This function computes a
!          weighting factor.  This function interpolates only for 
!          level k.

!          Notes:
!          Created by Brian for a non-uniform grid.
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

           k = t_lev

           ! Do actual interpolation

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
     .             ( gr%zt(1)-gr%zm(1) )/( gr%zm(2)-gr%zm(1) )
           endif

           ! Weight of upper momentum level on thermodynamic level.
           azm_weight(m_above) = factor
           ! Weight of lower momentum level on thermodynamic level.
           azm_weight(m_below) = 1.0 - factor

           return
           end function interpolated_aztk_imp

!-----------------------------------------------------------------------
           pure function gradzm( azm )

!          Description:
!          Function to compute the vertical derivative of a variable (azm)
!          located on the momentum grid. The results is returned in an 
!          array defined on the thermodynamic grid.
!-----------------------------------------------------------------------
           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azm

           ! Return Variable
           real, dimension(gr%nnzp) :: gradzm

           ! Local Variable
           integer :: k

           ! Compute vertical derivative

           do k=gr%nnzp, 2, -1
              gradzm(k) = ( azm(k) - azm(k-1) ) * gr%dzt(k)
           end do

!           gradzm(1) = UNDEFINED
!           gradzm(1) = 0.
           gradzm(1) = gradzm(2)

           return
           end function gradzm

!-----------------------------------------------------------------------
           pure function gradzt( azt )

!          Description:
!          Function to compute the vertical derivative of a variable (azt)
!          located on the thermodynamic grid. The results is returned in
!          an array defined on the momentum grid.
!-----------------------------------------------------------------------
           implicit none

           ! Input Variable
           real, intent(in), dimension(gr%nnzp) :: azt

           ! Output Variable
           real, dimension(gr%nnzp) :: gradzt

           ! Local Variable
           integer :: k

           ! Compute vertical derivative

           do k=1,gr%nnzp-1
              gradzt(k) = ( azt(k+1) - azt(k) ) * gr%dzm(k)
           end do
!           gradzt(gr%nnzp) = UNDEFINED
!           gradzt(gr%nnzp) = 0.
           gradzt(gr%nnzp) = gradzt(gr%nnzp-1)

           return
           end function gradzt

        end module grid_class

!------------------------------------------------------------------------
