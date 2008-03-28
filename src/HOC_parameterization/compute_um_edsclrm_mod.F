!------------------------------------------------------------------------
! $Id: compute_um_edsclrm_mod.F,v 1.1 2008-03-28 14:51:14 griffinb Exp $
!------------------------------------------------------------------------
        module compute_um_edsclrm_mod

        public :: compute_um_edsclrm

        contains

!------------------------------------------------------------------------
        subroutine compute_um_edsclrm( solve_type, xpwp_sfc, xm_tndcy, 
     .                                 Khm, dt, xm, xpwp, isValid )
!       Description:
!         Prognoses a horizontal wind component or other 
!         eddy diffusivity variable.
!         Diagnoses the turbulent flux.
!         Uses a Crank-Nicholson time-stepping algorithm.

!       References:
!       Eqn. 8 & 9 on p. 3545 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!------------------------------------------------------------------------

        use grid_class
        use lapack_wrap, only: tridag_solve
#ifdef STATS
        use stats_hoc
#endif
        use constants, only: fstderr
        
        implicit none

        ! Input Variables
        character(len=*), intent(in) ::
     .  solve_type ! Desc. of what is being solved for

        real, intent(in) :: 
     .  dt,          ! Timestep                             [s]
     .  xpwp_sfc     ! sfc flux                             [units vary]

        real, dimension(gr%nnzp), intent(in) :: 
     .  xm_tndcy, ! x tendency                              [units vary]
     .  Khm       ! Diffusion coefficient on momentum grid  [m^2/s]

        ! Input/Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  xm ! Prognostic array on the thermodynamic grid     [units vary]
        
        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  xpwp   ! Momentum flux                              [units vary]

        logical, intent(out) ::
     .  isValid ! false when matrix is singular

        ! Local Variables
        real,  dimension(gr%nnzp) :: a, b, c, rhs

        real :: atmp, ctmp

        integer :: k, kp1, km1 ! Indices

#ifdef STATS
        integer :: ixm_ta, ixm_bt

        select case ( trim( solve_type ) )
        case ( "um" )
          ixm_bt = ium_bt
          ixm_ta = ium_ta
        case ( "vm" )
          ixm_bt = ivm_bt
          ixm_ta = ivm_ta
        case default  ! Eddy scalars
          ixm_bt = 0
          ixm_ta = 0
        end select
#endif

#ifdef STATS
        if ( lstats_samp ) then
          ! xm total time tendency ( 1st calculation)
          if ( ixm_bt > 0 ) then
            zt%x(:,ixm_bt) = zt%x(:,ixm_bt) - xm / dt
          end if
        end if
#endif /*STATS*/

        ! Prepare tridiagonal system

        ! zt(1) is below ground, we don't have to worry about it

        a(1)   = 0.
        b(1)   = 1./dt
        c(1)   = 0.
        rhs(1) = 1./dt

        ! zt(2) is the first active model layer. We need to impose the 
        ! surface momentum flux xpwp_sfc.

        atmp = 0.
        ctmp = -0.5 * Khm(2) * gr%dzm(2) * gr%dzt(2)

        a(2)   = atmp
        c(2)   = ctmp
        b(2)   = - ctmp + 1./dt
        rhs(2) = ( ctmp + 1./dt ) * xm(2)
     .           - ctmp * xm(3)
     .           + xm_tndcy(2)
     .           + xpwp_sfc * gr%dzt(2)

#ifdef STATS
           if ( lstats_samp .and. ixm_ta > 0 ) then
             ztscr01(1) = 0.0
             ztscr02(1) = 0.0
             ztscr03(1) = 0.0

             ztscr01(2) = -atmp
             ztscr02(2) = ctmp
             ztscr03(2) = -ctmp
           end if
#endif
        ! Loop from level 3 to gr%nnzp

        do k=3, gr%nnzp-1, 1

           atmp = -0.5 * Khm(k-1) * gr%dzt(k) * gr%dzm(k-1)
           ctmp = -0.5 * Khm(k) * gr%dzt(k) * gr%dzm(k)

           a(k)   = atmp
           c(k)   = ctmp
           b(k)   = - atmp - ctmp + 1./dt
           rhs(k) = - atmp * xm(k-1)
     .              + ( atmp + ctmp + 1./dt ) * xm(k)
     .              - ctmp * xm(k+1)
     .              + xm_tndcy(k)
#ifdef STATS
           if ( lstats_samp .and. ixm_ta > 0 ) then
             ztscr01(k) = -atmp
             ztscr02(k) =  atmp + ctmp
             ztscr03(k) = -ctmp
           end if
#endif
        end do

        ! Level gr%nnzp. We impose zero flux from model top

! Vince Larson simplified the upper BC.  It still imposes
!    zero flux. 7 Jul 2007 
!        atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzm(gr%nnzp) 
!     .              * gr%dzt(gr%nnzp-1)
!        ctmp = 0.

!        a(gr%nnzp)   = atmp
!        c(gr%nnzp)   = ctmp
!        b(gr%nnzp)   = - atmp + 1./dt
!        rhs(gr%nnzp) = - atmp * xm(gr%nnzp-1)
!     .               + ( atmp + 1./dt ) * xm(gr%nnzp)

        a(gr%nnzp)   = -1.
        c(gr%nnzp)   =  0.
        b(gr%nnzp)   =  1.
        rhs(gr%nnzp) =  0.
! End of Vince Larson's change


#ifdef STATS
           if ( lstats_samp .and. ixm_ta > 0 ) then
             ztscr01(gr%nnzp) = -atmp
             ztscr02(gr%nnzp) =  atmp
             ztscr03(gr%nnzp) = -ctmp
           end if
#endif
!    Caused problems with DYCOMS II RF02
!     .               + xm_tndcy(gr%nnzp)

!    Attempted to compensate for the DYCOMS problem using the code below.
!    Doesn't actually seem to make a difference
!       atmp = -0.5 * Khm(gr%nnzp-1) * gr%dzt(gr%nnzp) * gr%dzm(gr%nnzp-1)

!       a(gr%nnzp)   = atmp
!       b(gr%nnzp)   = - c(gr%nnzp-1) + 1./dt
!       c(gr%nnzp)   = UNDEFINED
!       rhs(gr%nnzp) = - atmp * x(gr%nnzp-1)
!    .                 + ( atmp + 1./dt ) * x(gr%nnzp)
!    .                 + xm_tndcy(gr%nnzp)

        ! Store momentum flux (explicit component)
        xpwp(1) = xpwp_sfc
        do k=2,gr%nnzp-1
          xpwp(k) = -0.5 * Khm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) ) 
        end do
        xpwp(gr%nnzp) = 0.

#ifdef STATS
           ! Turbulent transport (explicit component)
           if ( lstats_samp .and. ixm_ta > 0 ) then
             do k=1,gr%nnzp,1
               km1 = max( k-1, 1 )
               kp1 = min( k+1, gr%nnzp )

               zt%x(k,ixm_ta) = zt%x(k,ixm_ta)
     .         + ztscr01(k) * xm(km1)
     .         + ztscr02(k) * xm(k)
     .         + ztscr03(k) * xm(kp1)
             end do
           end if
#endif

        ! Solve tridiagonal system
        call tridag_solve( solve_type, gr%nnzp, 1, c, b, a, 
     .                     rhs, xm, isValid )
        if (.not. isValid ) then
          return
        end if

        ! Second part of momentum (implicit component)
        do k=2, gr%nnzp-1, 1
          xpwp(k) = xpwp(k) 
     .               - 0.5 * Khm(k) * gr%dzm(k) * ( xm(k+1) - xm(k) )
        end do

#ifdef STATS
        if ( lstats_samp ) then
          do k = 1, gr%nnzp, 1
            km1 = max( k-1, 1 )
            kp1 = min( k+1, gr%nnzp )

            ! x time tendency (2nd calculation)
            if ( ixm_bt > 0 ) then
              zt%x(k,ixm_bt) = zt%x(k,ixm_bt) + xm(k) / dt
              zt%n(k,ixm_bt) = zt%n(k,ixm_bt) + 1
            end if
            ! x turbulent transport (implicit component)
            if ( ixm_ta > 0 ) then
              zt%x(k,ixm_ta) = zt%x(k,ixm_ta)
     .        + ztscr01(k) * xm(km1)
     .        + ztscr02(k) * xm(k)
     .        + ztscr03(k) * xm(kp1)
              zt%n(k,ixm_ta) = zt%n(k,ixm_ta) + 1
             end if
           end do

         end if
#endif

!       Error Report
!       Joshua Fasching February 2008
        if ( .not. isValid ) then
                
           write(fstderr,*) "Error in compute_um_edsclrm"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "dt = ", dt
           write(fstderr,*) "xpwp_sfc = ", xpwp_sfc
           write(fstderr,*) "xm_tndcy = ", xm_tndcy
           write(fstderr,*) "Khm = ", Khm
           
           write(fstderr,*) "Intent(inout)"
           
           write(fstderr,*) "xm = ", xm
           
           write(fstderr,*) "Intent(out)"
           
           write(fstderr,*) "xpwp = ", xpwp
           
        endif

        return

        end subroutine compute_um_edsclrm
!-----------------------------------------------------------------------

        end module compute_um_edsclrm_mod
