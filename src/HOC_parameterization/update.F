!------------------------------------------------------------------------
! $Id: update.F,v 1.2 2006-09-05 22:39:04 dschanen Exp $

        subroutine update_ed_var( solve_type, x, xt, xw0, Khm, 
     .                            dt, xw, isValid )
!       Description:
!       Computes wind or other eddy diffusivity variable.

!       References:
!       Eqn. 8 & 9 on p. 3545 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!------------------------------------------------------------------------

        use grid_class
        use lapack_wrap, only: tridag_solve

        implicit none

        ! Input Variables
        character(len=*), intent(in) ::
     .  solve_type ! Desc. of what is being solved for

        real, intent(in) :: 
     .  dt,  ! Timestep                                 [s]
     .  xw0  ! sfc flux                                 [units vary]

        real, dimension(1:gr%nnzp), intent(in) :: 
     .  xt, ! x tendency                                [units vary]
     .  Khm ! Diffusion coefficient on momentum grid    [m^2/s]

        ! Input/Output Variables
        real, dimension(1:gr%nnzp), intent(inout) :: 
     .  x  ! Prognostic array on the thermodynamic grid [units vary]

        
        ! Output Variables
        real, dimension(1:gr%nnzp), intent(out) :: 
     .  xw   ! Momentum flux                            [units vary]

        logical, intent(out) :: isValid ! false when matrix is singular

        ! Local Variables
        real,  dimension(:), allocatable :: a, b, c, r

        real :: atmp, ctmp

        integer :: k ! Index


        ! Allocate memory for the tridiagonal matrix
        allocate( a(1:gr%nnzp) ) ! sub diagonal
        allocate( b(1:gr%nnzp) ) ! main diagonal
        allocate( c(1:gr%nnzp) ) ! super diagonal
        allocate( r(1:gr%nnzp) ) ! right hand side

        ! Prepare tridiagonal system

        ! zt(1) is below ground, we don't have to worry about it

        a(1) = 0.
        b(1) = 1./dt
        c(1) = 0.
        r(1) = 1./dt

        ! zt(2) is the first active model layer. We need to impose the 
        ! surface momentum flux xw0.

        atmp = 0.
        ctmp = -0.5 * Khm(2) * gr%dzm(2) * gr%dzt(2)

        a(2) = atmp
        c(2) = ctmp
        b(2) = - ctmp + 1./dt
        r(2) = ( ctmp + 1./dt ) * x(2)
     .         - ctmp * x(3)
     .         + xt(2)
     .         + xw0 * gr%dzt(2)

        ! Loop from level 3 to gr%nnzp

        do k=3, gr%nnzp-1, 1

           atmp = -0.5 * Khm(k-1) * gr%dzt(k) * gr%dzm(k-1)
           ctmp = -0.5 * Khm(k) * gr%dzt(k) * gr%dzm(k)

           a(k) = atmp
           c(k) = ctmp
           b(k) = - atmp - ctmp + 1./dt
           r(k) = - atmp * x(k-1)
     .            + ( atmp + ctmp + 1./dt ) * x(k)
     .            - ctmp * x(k+1)
     .            + xt(k)

        end do

        ! Level gr%nnzp. We impose zero flux from model top

        atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzm(gr%nnzp) 
     .              * gr%dzt(gr%nnzp-1)
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
!       atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzt(gr%nnzp) * gr%dzm(gr%nnzp-1)

!       a(gr%nnzp) = atmp
!       b(gr%nnzp) = - c(gr%nnzp-1) + 1./dt
!       c(gr%nnzp) = UNDEFINED
!       r(gr%nnzp) = - atmp * x(gr%nnzp-1)
!    .               + ( atmp + 1./dt ) * x(gr%nnzp)
!    .               + xt(gr%nnzp)

        ! Store momentum flux (explicit component)
        xw(1) = xw0
        do k=2,gr%nnzp-1
          xw(k) = -0.5 * Khm(k) * gr%dzm(k) * ( x(k+1) - x(k) ) 
        end do
        xw(gr%nnzp) = 0.

        ! Solve tridiagonal system
        call tridag_solve( solve_type, gr%nnzp, 1, c, b, a, 
     .                     r, x, isValid )
        if (.not. isValid ) then
          deallocate( a, b, c, r )
          return
        end if

        ! Second part of momentum (implicit component)
        do k=2, gr%nnzp-1, 1
          xw(k) = xw(k) - 0.5 * Khm(k) * gr%dzm(k) * ( x(k+1) - x(k) ) 
        end do

        ! Deallocate memory
        deallocate( a, b, c, r )

        return
        end subroutine update_ed_var
!-----------------------------------------------------------------------
