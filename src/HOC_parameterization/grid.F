!------------------------------------------------------------------------
! $Id: grid.F,v 1.3 2006-07-18 01:36:41 ajsmith4 Exp $

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
!-----------------------------------------------------------------------
        module grid_class

           implicit none

           type grid
              integer :: nnzp
!   Note: Fortran 90/95 prevent an allocatable array from appearing
!   within a derived type.  However, a pointer can be used in the same
!   manner as an allocatable array, as we have done here (the grid 
!   pointers are always allocated rather than assigned and nullified 
!   like real pointers).
              real, pointer, dimension(:) :: zm, zt
              real, pointer, dimension(:) :: dzm, dzt
           end type grid

           type (grid) gr       ! The grid is defined here so that it
                                ! is common throughout the module
                                ! The implication is that only one grid
                                ! can be defined !

           interface zt2zm
              module procedure interpolated_azm, interpolated_azmk
           end interface

           interface zm2zt
              module procedure interpolated_azt, interpolated_aztk
           end interface

           interface ddzm
              module procedure gradzm
           end interface

           interface ddzt
              module procedure gradzt
           end interface

           contains

!-----------------------------------------------------------------------
!           Grid Constructor
!
!           This subroutine sets up a uniformly spaced vertical grid similar
!           to RAMS grid.
!
!          Input parameters:
!
!             integer nnzp        number of vertical levels in grid
!             real deltaz        vertical grid spacing
!-----------------------------------------------------------------------
           subroutine gridsetup( nnzp, deltaz, implemented,
     .                           host_momentum, host_thermodynamic )

           implicit none

!          Parameter

c           integer, parameter :: NWARNING = 200

c Modification by Adam Smith, 27 June 2006
c We want to increase the maximum vertical grid size
           integer, parameter :: NWARNING = 250

!          Input

           integer, intent(in) :: nnzp
           real, intent(in)    :: deltaz

           ! Flag to see if HOC is running on it's own or it it's part
           ! of a larger model.
           logical, intent(in) :: implemented
           ! If it's part of a larger model, it needs to use that models
           ! grid heights.
           real, dimension(1:nnzp), intent(in) :: host_momentum
           real, dimension(1:nnzp), intent(in) :: host_thermodynamic

!          Internal

           integer i, ierr

!          Define grid size

           if ( nnzp > NWARNING ) then
             write(*,*) "Warning: running with vertical grid "//
     .                  "which is larger than ", NWARNING
           endif

           gr%nnzp = nnzp
           allocate( gr%zm(1:nnzp), gr%zt(1:nnzp),
     .               gr%dzm(1:nnzp), gr%dzt(1:nnzp), stat=ierr )
           if ( ierr /= 0 ) stop "grid allocation failure"

           IF ( .NOT. implemented ) THEN

!             Original code
!             Define momentum grid points. The first momentum level is 
!             at the surface
              DO i = 1, gr%nnzp, 1
                 gr%zm(i) = (i-1) * deltaz
              ENDDO

!             Define thermodynamic grid. For the time being, themodynamic 
!             grid points are located in the center point between momentum 
!             points. The first thermodynmic point is below the surface.
              gr%zt(1) = -0.5 * deltaz
              DO i = 2, gr%nnzp, 1
                 gr%zt(i) = 0.5 * ( gr%zm(i) + gr%zm(i-1) )
              ENDDO

           ELSE

!             Momentum levels from host model.
              DO i = 1, gr%nnzp, 1
                 gr%zm(i) = host_momentum(i)
              ENDDO

!             Thermodynamic levels from host model after possible 
!             grid-index adjustment for HOC interface.
              DO i = 1, gr%nnzp, 1
                 gr%zt(i) = host_thermodynamic(i)
              ENDDO

           ENDIF

!          Define dzm. dzm is the inverse spacing between thermodynamic 
!          points

           do i=1,gr%nnzp-1
              gr%dzm(i) = 1. / ( gr%zt(i+1) - gr%zt(i) )
           end do
           gr%dzm(gr%nnzp) = gr%dzm(gr%nnzp-1)

!          Define dzt. dzt is the inverse spacing between momentum points

           do i=2,gr%nnzp
              gr%dzt(i) = 1. / ( gr%zm(i) - gr%zm(i-1) )
           end do
           gr%dzt(1) = gr%dzt(2)

           return
           end subroutine gridsetup

!-----------------------------------------------------------------------
!           Function to interpolate a variable located on the 
!          thermodynamic grid (azt) to a variable on the momentum 
!          grid (azm)
!-----------------------------------------------------------------------

           function interpolated_azm( azt )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azt

!          Output

           real, dimension(1:gr%nnzp) :: interpolated_azm
        
!          Internal

           integer i

!          Do actual interpolation 

           do i=1, gr%nnzp-1
!              interpolated_azm(i) = 0.5 * ( azt(i) + azt(i+1) )
              interpolated_azm(i) = 
     .                ( ( azt(i+1)-azt(i) )/( gr%zt(i+1)-gr%zt(i) ) )
     .               * ( gr%zm(i)-gr%zt(i) ) + azt(i)
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
!           Function to interpolate a variable located on the
!           thermodynamic grid (azt) to a variable on the momentum
!           grid (azm). The function interpolates only for level k
!-----------------------------------------------------------------------

           function interpolated_azmk( azt, k )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azt
           integer, intent(in)                    :: k

!           Output

           real interpolated_azmk

!           Do actual interpolation

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
!           Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt)
!-----------------------------------------------------------------------

           function interpolated_azt( azm )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azm

!          Output

           real, dimension(1:gr%nnzp) :: interpolated_azt
        
!          Internal

           integer i

!          Do actual interpolation 

           do i=gr%nnzp,2,-1
!              interpolated_azt(i) = 0.5 * ( azm(i-1) + azm(i) )
              interpolated_azt(i) =
     .                ( ( azm(i)-azm(i-1) )/( gr%zm(i)-gr%zm(i-1) ) )
     .               * ( gr%zt(i)-gr%zm(i-1) ) + azm(i-1)
           end do
!           interpolated_azt(1) = UNDEFINED
!           interpolated_azt(1) = azm(1)
           interpolated_azt(1) =
     .                ( ( azm(2)-azm(1) )/( gr%zm(2)-gr%zm(1) ) )
     .               * ( gr%zt(1)-gr%zm(1) ) + azm(1)

           return
           end function interpolated_azt

!-----------------------------------------------------------------------
!           Function to interpolate a variable located on the 
!          momentum grid (azm) to a variable on the thermodynamic
!          grid (azt). This function interpolates only for level k.
!-----------------------------------------------------------------------

           function interpolated_aztk( azm, k )
           implicit none

!          Input

           real, dimension(1:gr%nnzp), intent(in) :: azm
           integer, intent(in)                    :: k

!          Output

           real interpolated_aztk
        
!          Do actual interpolation 

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
!          Function to compute the vertical derivative of a variable (azm)
!          located on the momentum grid. The results is returned in an 
!          array defined on the thermodynamic grid.
!-----------------------------------------------------------------------

           function gradzm( azm )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azm

!          Output

           real, dimension(1:gr%nnzp) :: gradzm

!           Internal

           integer i

!           Compute vertical derivative

           do i=gr%nnzp, 2, -1
              gradzm(i) = ( azm(i) - azm(i-1) ) * gr%dzt(i)
           end do
!           gradzm(1) = UNDEFINED
!           gradzm(1) = 0.
           gradzm(1) = gradzm(2)

           return
           end function gradzm

!-----------------------------------------------------------------------
!           Function to compute the vertical derivative of a variable (azt)
!          located on the thermodynamic grid. The results is returned in
!          an array defined on the momentum grid.
!-----------------------------------------------------------------------

           function gradzt( azt )
           implicit none

!           Input

           real, dimension(1:gr%nnzp), intent(in) :: azt

!          Output

           real, dimension(1:gr%nnzp) :: gradzt

!           Internal

           integer i

!          Compute vertical derivative

           do i=1,gr%nnzp-1
              gradzt(i) = ( azt(i+1) - azt(i) ) * gr%dzm(i)
           end do
!           gradzt(gr%nnzp) = UNDEFINED
!           gradzt(gr%nnzp) = 0.
           gradzt(gr%nnzp) = gradzt(gr%nnzp-1)

           return
           end function gradzt

        end module grid_class

!------------------------------------------------------------------------
