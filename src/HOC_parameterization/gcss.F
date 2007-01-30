#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: gcss.F,v 1.28 2007-01-30 17:50:00 dschanen Exp $ 

        module gcss

!       Description:
!       Old file for the GCSS forcings, etc.  Add new cases to the GCSS
!       directory

!       References:
!       <http://www.gewex.org/gcss.html>
!----------------------------------------------------------------------
        use model_flags

        implicit none

        private ! Default to private

        public ::
     .  rico_tndcy,        ! Michael Falk, 13 Dec 2006
     .  rico_sfclyr        ! Michael Falk, 13 Dec 2006

        contains

!----------------------------------------------------------------------
        subroutine rico_tndcy
     .  ( dt, exner, wmt, wmm, theta_tndcy,
     .    rt_tndcy, radht,
     .    sclr_tndcy )

!        Description:
!          Subroutine to apply case-specific forcings to RICO case
!          (Michael Falk, 13 Dec 2006).
!
!        References:
!-----------------------------------------------------------------------

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif /*STATS*/
        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  dt            ! Current length of timestep      [s]

        real, dimension(gr%nnzp), intent(in) ::
     .  exner        ! Exner function                         [-]

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  wmt,         ! Large-scale vertical motion on t grid   [m/s]
     .  wmm,         ! Large-scale vertical motion on m grid   [m/s]
     .  theta_tndcy, ! Large-scale thlm tendency               [K/s]
     .  rt_tndcy,    ! Large-scale rtm tendency                [kg/kg/s]
     .  radht        ! dT/dt, then d Theta/dt, due to rad.     [K/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar LS tendency            [units/s]

        ! Local Variables, general
        integer :: i, k ! Loop indices
        real    :: t_tendency ! Temperature (not potential temperature) tendency

        ! Compute vertical motion
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            wmt(k) = -(0.005 / 2100) * gr%zt(k)
          else if (gr%zt(k) < 3000) then
            wmt(k) = -0.005 + (0.005 / (3000 - 2100))
     .               * (gr%zt(k) - 2100)
          else
            wmt(k) = 0.
          end if
        end do
        wmm = zt2zm( wmt )

        ! Boundary conditions on vertical motion.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        

        ! Compute large-scale horizontal temperature advection
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            t_tendency = - 0.8 / 86400  ! Units [K/s]
          else if (gr%zt(k) < 3000) then
            t_tendency = - 0.8 / 86400
     .                   + (0.8 / 86400) / (3000 - 2100)
     .                   * (gr%zt(k) - 2100)  ! Units [K/s]
          else
            t_tendency = 0.  ! Units [K/s]
          end if
          ! Convert to units of [K/s] but potential T instead of T
!          theta_tndcy(k) = (t_tendency * ((psfc/p(k)) ** (Rd/Cp)))
          theta_tndcy(k) = (t_tendency / exner(k))
        end do


        ! Compute large-scale horizontal moisture advection [g/kg/s]
        do k=1,gr%nnzp
          if (gr%zt(k) < 550) then
            rt_tndcy(k) = - 1.0 / 86400  ! Units [g/kg/s]
          else if (gr%zt(k) < 1500) then
            rt_tndcy(k) = - 1.0 / 86400
     .                    + ((0.0 + 1.0) / 86400)
     .                    / (1500 - 550)
     .                   * (gr%zt(k) - 550)  ! Units [g/kg/s]
          else
            rt_tndcy(k) = 0.  ! Units [K/s]
          end if
          rt_tndcy(k) = rt_tndcy(k) / 1000.  ! Converts [g/kg/s] to [kg/kg/s]
        end do


        ! Compute radiative cooling
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            radht(k) = - (2 / 86400)
     .              + ((-1.3 + 2) / 86400) / 2100 * gr%zt(k)  ! Units [K/s]
          else if (gr%zt(k) < 3000) then
            radht(k) = - (1.3 / 86400)
     .              + (1.3 / 86400) / (3000 - 2100) * (gr%zt(k) - 2100) ! Units [K/s]
          else
            radht(k) = 0.
          end if
!          radht(k) = (radht(k) * ((psfc/p(k)) ** (Rd/Cp)))
          radht(k) = (radht(k) / exner(k))
          theta_tndcy(k) = theta_tndcy(k) + radht(k)
        end do


        ! Set scalar tendencies
        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(:,SCLR_THETA) = theta_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine rico_tndcy

!----------------------------------------------------------------------
        subroutine rico_sfclyr( up, vp, thp, rtp,
     .                          lowestlevel, sst, psfc,
     .                          uw, vw, wtp, wrp,
     .                          sclrm_sfc, wpsclrp_sfc )

!----------------------------------------------------------------------
!        Description:
!          Surface forcing subroutine for MPACE case.  Written July-
!          November 2006 by Michael Falk.
!
!        References:
!          MPACE specification, arm.gov
!-----------------------------------------------------------------------

        use constants
        implicit none

!       External
        real, external :: rsat

!       Constants
        real, parameter :: ubmin = 0.25    ! I don't know where this comes from.  All of the other cases use it, though.
        real, parameter :: ustar = 0.3     ! Defined by ATEX specification
        real, parameter :: C_10  = 0.0013  ! Defined by ATEX specification
        real, parameter :: z0    = 0.00015 ! Roughness length, defined by ATEX specification

!       Internal variables
        real ubar,Cz

!       Input variables
        real, intent(in) :: up,vp,thp,rtp,lowestlevel,sst,psfc

!       Input variables (optional)
        real, optional, intent(in) :: sclrm_sfc(sclr_dim)

!       Output variables
        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)
        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

! Compute heat and moisture fluxes

        ubar = max(ubmin, sqrt(up*up + vp*vp))
        Cz   = C_10 * ((log(10/z0))/(log(lowestlevel/z0))) * 
     .         ((log(10/z0))/(log(lowestlevel/z0)))              ! Modification in case lowest model level isn't at 10 m, from ATEX specification
        wtp  = -Cz * ubar * ( thp - sst * (p0/psfc)**kappa )     ! K * m/s
        wrp  = -Cz * ubar * ( rtp - rsat(psfc,sst) )             ! kg/kg * m/s

! Optional passive scalar
        if ( present( wpsclrp_sfc ) .and. present( sclrm_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = -C_10 * ubar * 
     .           ( sclrm_sfc(SCLR_THETA) - sst * (p0/psfc)**kappa )

          wpsclrp_sfc(SCLR_RT) = -C_10 * ubar * 
     .                     ( sclrm_sfc(SCLR_RT) - rsat(psfc,sst) )
        end if

! Compute momentum fluxes

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        ! Compute scalar fluxes
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine rico_sfclyr

      end module gcss
!------------------------------------------------------------------------
