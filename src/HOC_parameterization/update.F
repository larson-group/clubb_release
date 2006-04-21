!------------------------------------------------------------------------
! $Id: update.F,v 1.1 2006-04-21 21:50:56 dschanen Exp $

! SUBROUTINE update_ed_var( )
! Description:
! Computes wind or other variable eddy diffusivity variable
!------------------------------------------------------------------------
        subroutine update_ed_var( x, xt, xw0, K, dt, xw, isValid )

        use grid_class

        implicit none

!       Input/Output array

        real, dimension(1:gr%nnzp), intent(inout) :: x  ! Prognostic array on
                                                        ! thermodynamic grid

!       Input arrays

        real, dimension(1:gr%nnzp), intent(in) :: xt    ! x tendency

        real, intent(in)                       :: xw0   ! sfc flux

        real, dimension(1:gr%nnzp), intent(in) :: K     ! Diffusion coefficient
                                                        ! on momentum grid

        real, intent(in) :: dt                          ! Timestep
        
!       Output array

        real, dimension(1:gr%nnzp), intent(out) :: xw   ! momentum flux

!       Internal

        real, allocatable :: a(:), b(:), c(:), r(:)

        real atmp, ctmp

        integer i

!       Miscellaneous
        logical, intent(out) :: isValid

!       Allocate memory

        allocate( a(1:gr%nnzp) ) ! sub diagonal
        allocate( b(1:gr%nnzp) ) ! main diagonal
        allocate( c(1:gr%nnzp) ) ! super diagonal
        allocate( r(1:gr%nnzp) ) ! right hand side

!       Prepare tridiagonal system

!       zt(1) is below ground, we don't have to worry about it

        a(1) = 0.
        b(1) = 1./dt
        c(1) = 0.
        r(1) = 1./dt

!       zt(2) is the first active model layer. We need to impose the surface
!       momentum flux xw0.

        atmp = 0.
        ctmp = -0.5 * K(2) * gr%dzm(2) * gr%dzt(2)

        a(2) = atmp
        c(2) = ctmp
        b(2) = - ctmp + 1./dt
        r(2) = ( ctmp + 1./dt ) * x(2)
     .         - ctmp * x(3)
     .         + xt(2)
     .         + xw0 * gr%dzt(2)

!       loop for levels 3 to gr%nnzp

        do i=3,gr%nnzp-1

           atmp = -0.5 * K(i-1) * gr%dzt(i) * gr%dzm(i-1)
           ctmp = -0.5 * K(i) * gr%dzt(i) * gr%dzm(i)

           a(i) = atmp
           c(i) = ctmp
           b(i) = - atmp - ctmp + 1./dt
           r(i) = - atmp * x(i-1)
     .            + ( atmp + ctmp + 1./dt ) * x(i)
     .            - ctmp * x(i+1)
     .            + xt(i)

        end do

!       level gr%nnzp. We impose zero flux from model top

        atmp = -0.5 * K(gr%nnzp-1) * gr%dzm(gr%nnzp) * gr%dzt(gr%nnzp-1)
        ctmp = 0.

        a(gr%nnzp) = atmp
        c(gr%nnzp) = ctmp
        b(gr%nnzp) = - atmp + 1./dt
        r(gr%nnzp) = - atmp * x(gr%nnzp-1)
     .               + ( atmp + 1./dt ) * x(gr%nnzp)
!    Caused problems with DYCOMS II RF02
!     .               + xt(gr%nnzp)

!    Attempted to compensate for the DYCOMS problem using the code below.
!    Doesn't actually seem to make a difference
!       atmp = -0.5 * K(gr%nnzp-1) * gr%dzt(gr%nnzp) * gr%dzm(gr%nnzp-1)

!       a(gr%nnzp) = atmp
!       b(gr%nnzp) = - c(gr%nnzp-1) + 1./dt
!       c(gr%nnzp) = UNDEFINED
!       r(gr%nnzp) = - atmp * x(gr%nnzp-1)
!    .               + ( atmp + 1./dt ) * x(gr%nnzp)
!    .               + xt(gr%nnzp)



!       Store momentum flux (explicit component)

        xw(1) = xw0
        do i=2,gr%nnzp-1
           xw(i) = -0.5 * K(i) * gr%dzm(i) * ( x(i+1) - x(i) ) 
        end do
        xw(gr%nnzp) = 0.

!       Solve tridiagonal system
        call tridag( a, b, c, r, x, gr%nnzp, isValid)
        if (.not. isValid ) then
           print *, "update: tridag failed"
           deallocate( a, b, c, r )
           return
        endif

!       Second part of momentum (implicit component)

        do i=2,gr%nnzp-1
           xw(i) = xw(i) - 0.5 * K(i) * gr%dzm(i) * ( x(i+1) - x(i) ) 
        end do

!       Deallocate memory

        deallocate( a, b, c, r )

        return
        end subroutine update_ed_var
!-----------------------------------------------------------------------
