#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: cobra.F,v 1.1 2007-01-30 17:47:48 dschanen Exp $
        module cobra
!       Description:
!       Contains subroutines for the COBRA CO2 case.
!----------------------------------------------------------------------

        implicit none

        private

        public ::
     .  cobra_tndcy,
     .  cobra_sfclyr

        contains

!----------------------------------------------------------------------
        subroutine cobra_tndcy( time, wmt, wmm, 
     .                          thlm_forcing, rtm_forcing,
     .                          sclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for COBRA CO2 case

!       References:
!       None
!----------------------------------------------------------------------
        use grid_class
        use constants, only: sclr_dim

        implicit none

        ! Input Variables
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid                 [m/s]
     .  wmm,         ! w wind on momentum grid                      [m/s]
     .  thlm_forcing,! Liquid water potential temperature tendency  [K/s]
     .  rtm_forcing  ! Total water mixing ratio tendency            [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing ! Passive scalar tendency [units vary]

        ! Local Variables
        integer :: k

        ! Large-scale subsidence

        DO k = 2, gr%nnzp, 1
          wmt(k) = -5.0E-6 * gr%zt(k)
        END DO

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       No large-scale water tendency or cooling

        rtm_forcing  = 0.0
        thlm_forcing = 0.0

        if ( present( sclrm_forcing ) ) then
!          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = 0.0
          sclrm_forcing(1:gr%nnzp,SCLR_RT) 
     .    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine cobra_tndcy

!-----------------------------------------------------------------------
        subroutine cobra_sfclyr( time, z, dn0, thp, up, vp,
     .                           uw, vw, wtp, wrp, ustar,
     .                           wpsclrp_sfc, wpedsclrp_sfc )

!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to the format used for the GCSS ARM 
!       case.
!
!       The data has been altered so it can be used for the COBRA CO2 
!       case.
!-----------------------------------------------------------------------

        use constants, only: Cp, Lv, grav, sclr_dim

        implicit none

        ! External
        real, external :: diag_ustar

        ! Parameter Constants

        real, parameter ::
     .  ubmin = 0.25

        ! Input variables

        real, intent(in) :: time
        real, intent(in) :: z
        real, intent(in) :: dn0
        real, intent(in) :: thp
        real, intent(in) :: up,vp

        real, intent(inout) :: ustar

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)
        real, optional, intent(out) :: wpedsclrp_sfc(sclr_dim)

!       Internal variables

        integer i, j
        real usfc, vsfc, ubar
        real true_time
        real heat_flx, moisture_flx
        real heat_flx2, moisture_flx2
        real bflx

        real CO2_flx, CO2_flx2
        real, parameter :: M_da = 0.02897  ! Molecular weight of dry air.

!       COBRA roughness height
!        real, parameter :: z0 = 0.035  ! ARM momentum roughness height
        real, parameter :: z0 = 1.75   ! momentum roughness height

! Compute heat and moisture fluxes from ARM data in (W/m2)

        true_time = time
        call cobra_sfcflx( true_time, heat_flx, moisture_flx
#ifdef SCALARS
     .                    ,CO2_flx
#endif /* SCALARS */
     .                   )

! Compute momentum fluxes

!       Convert heat_flx and moisture_flx to natural units
        heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
        moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

#ifdef SCALARS
!       Convert CO2 surface flux to natural units.
!       The CO2 flux has been given in units of:  umol/(m^2 s).
!       umol stands for micromoles.  The CO2 concentration in
!       this code is in units of ppmv, which is also the molar
!       mixing ratio times 10^6.
!       The units are:  10^6 * [ mol (CO2) / mol (dry air) ].
!       w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
!       where M is the molecular weight of dry air.
        CO2_flx2 = CO2_flx * ( M_da / dn0 )
#endif /* SCALARS */

!       Heat flux in units of (m2/s3) (needed by diag_ustar)
        bflx = grav/thp * heat_flx2

!       Sfc winds
        usfc = up
        vsfc = vp
        ubar = max( ubmin, sqrt( usfc ** 2 + vsfc ** 2 ) )

!       Compute ustar
        ustar = diag_ustar( z, bflx, ubar, z0 )

!       Assign fluxes

        uw  = -usfc/ubar * ustar * ustar
        vw  = -vsfc/ubar * ustar * ustar
        wtp = heat_flx2
        wrp = moisture_flx2

!        if (present( wpsclrp_sfc ))
!     .                 wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc ))
     .                 wpsclrp_sfc(SCLR_THETA) = CO2_flx2
        if (present( wpsclrp_sfc ))
     .                 wpsclrp_sfc(SCLR_RT) = wrp

!        if (present( wpedsclrp_sfc ))
!     .                 wpedsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpedsclrp_sfc ))
     .                 wpedsclrp_sfc(SCLR_THETA) = CO2_flx2
        if (present( wpedsclrp_sfc ))
     .                 wpedsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine cobra_sfclyr

!-----------------------------------------------------------------------
! This subroutine computes surface heat and moisture for a specific time
! according to the format used for the GCSS ARM case.
! Fluxes returned are in (W/m2).
!
! The data has been altered so it can be used for the COBRA CO2 case.

        subroutine cobra_sfcflx( time, heat_flx, moisture_flx,
     .                           CO2_flx )

        implicit none

!       Input variable
        real, intent(in) :: time       ! time in seconds

!       Output variables
        real, intent(out) :: heat_flx, moisture_flx

!       Output variables (optional)
        real, optional, intent(out) :: CO2_flx

!       Internal variables
        integer i1,i2
        real a

        integer ntimes
        parameter ( ntimes = 6 )
        real times(ntimes), H(ntimes), LE(ntimes)
#ifdef SCALARS
        real CO2(ntimes)
#endif /* SCALARS */
        ! Best fit line
!        data times / 57600.0, 64800.0, 73800.0, 75600.0 /
        ! Every hour (3 half-hour point weighted average)
        data times / 57600.0, 61200.0, 64800.0, 68400.0,
     .               72000.0, 75600.0 /

!       H and LE specifications
        ! Every hour (3 half-hour point weighted average)
        data H   / 415.3, 388.3, 404.0, 326.8, 
     .             236.9, 220.1 /
        data LE  / 121.2, 119.1, 122.9, 102.7,
     .             104.8,  87.1 /

#ifdef SCALARS
!       CO2 specifications
        ! Every hour (3 half-hour point weighted average)
        data CO2 / -12.4, -11.3, -10.4,  -8.0,
     .              -7.2,  -7.1 /
#endif /* SCALARS */

        if ( time <= times(1) ) then
           heat_flx     = H(1)
           moisture_flx = LE(1)
        else if ( time >= times(ntimes) ) then
           heat_flx     = H(ntimes)
           moisture_flx = LE(ntimes)
        else
           i1 = 1
           do while ( i1 <= ntimes-1 )
              i2 = i1 + 1
              if ( time >= times(i1) .and. time < times(i2) ) then
                 a            = (time-times(i1))/(times(i2)-times(i1))
                 heat_flx     = ( 1. - a ) * H(i1) + a * H(i2)
                 moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
#ifdef SCALARS
                 CO2_flx      = ( 1. - a ) * CO2(i1) + a * CO2(i2)
#endif /* SCALARS */
                 i1           = ntimes
              end if
              i1 = i2
           end do
        end if

        return
        end subroutine cobra_sfcflx

        end module cobra
