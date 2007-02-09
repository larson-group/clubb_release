#define SCLR_CO2 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: cobra.F,v 1.2 2007-02-09 22:10:11 dschanen Exp $
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
          sclrm_forcing(1:gr%nnzp,SCLR_CO2) = 0.0
          sclrm_forcing(1:gr%nnzp,SCLR_RT) 
     .    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine cobra_tndcy

!-----------------------------------------------------------------------
        subroutine cobra_sfclyr( time, z, dn0, thlm_sfc, um_sfc, vm_sfc,
     .                           upwp_sfc, vpwp_sfc, 
     .                           wpthlp_sfc, wprtp_sfc, ustar,
     .                           wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to the format used for the GCSS ARM 
!       case.

!       Notes:
!       The data has been altered so it can be used for the COBRA CO2 
!       case.
!-----------------------------------------------------------------------

        use constants, only: Cp, Lv, grav, sclr_dim

        implicit none

        ! External
        real, external :: diag_ustar

        intrinsic :: sqrt, max, present

        ! Parameter Constants
        real, parameter ::
     .  ubmin = 0.25,
     .  M_da  = 0.02897  ! Molecular weight of dry air.

        ! Input variables
        real, intent(in) :: 
     .  time,     ! Current time                [s]
     .  z,        ! Elevation at zt=2           [m]
     .  dn0,      ! Air density at surface      [kg/m^3]
     .  thlm_sfc, ! Theta_l at zt(2)            [K]
     .  um_sfc,   ! u wind at zt(2)             [m/s]
     .  vm_sfc    ! v wind at zt(2)             [m/s]

        real, intent(inout) :: ustar

        ! Output variables
        real, intent(out) :: 
     .  upwp_sfc,   ! u'w' at surface           [m^2/s^2]
     .  vpwp_sfc,   ! v'w' at surface           [m^2/s^2]
     .  wpthlp_sfc, ! w'theta_l' surface flux   [(m K)/s]
     .  wprtp_sfc   ! w'rt' surface flux        [(m kg)/(kg s)]

        ! Output variables (optional)
        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! w'sclr' surface flux          [units m/s]
     .  wpedsclrp_sfc ! w' edsclr' surface flux       [units m/s]

        ! Local variables
        integer :: i, j

        real :: 
     .  usfc, vsfc, ubar,
     .  true_time,
     .  heat_flx, moisture_flx,
     .  heat_flx2, moisture_flx2,
     .  bflx

        real :: CO2_flx, CO2_flx2

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

        ! Convert heat_flx and moisture_flx to natural units
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

        ! Heat flux in units of (m2/s3) (needed by diag_ustar)
        bflx = grav/thlm_sfc * heat_flx2

        ! Surface winds
        usfc = um_sfc
        vsfc = vm_sfc
        ubar = max( ubmin, sqrt( usfc ** 2 + vsfc ** 2 ) )

        ! Compute ustar
        ustar = diag_ustar( z, bflx, ubar, z0 )

        ! Assign fluxes
        upwp_sfc   = -usfc/ubar * ustar * ustar
        vpwp_sfc   = -vsfc/ubar * ustar * ustar
        wpthlp_sfc = heat_flx2
        wprtp_sfc  = moisture_flx2

        if ( present( wpsclrp_sfc ) ) then
!         wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_CO2) = CO2_flx2
          wpsclrp_sfc(SCLR_RT)  = wprtp_sfc
        end if

        if ( present( wpedsclrp_sfc ) ) then
!         wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_CO2) = CO2_flx2
          wpedsclrp_sfc(SCLR_RT)  = wprtp_sfc
        end if

        return
        end subroutine cobra_sfclyr

!-----------------------------------------------------------------------
        subroutine cobra_sfcflx( time, heat_flx, moisture_flx,
     .                           CO2_flx )
!       Description:
!       This subroutine computes surface heat and moisture for a specific
!       time according to the format used for the GCSS ARM case.
!       Heat and moisture fluxes are returned in (W/m2).

!       Notes:
!       The data has been altered so it can be used for the COBRA CO2 case.
!-----------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: present

        ! Constant parameters
        integer, parameter :: ntimes = 6

        real, parameter, dimension(ntimes) ::
        ! Best fit line
!    .   times = (/ 57600.0, 64800.0, 73800.0, 75600.0 /)
        ! Every hour (3 half-hour point weighted average)
     .  times = (/ 57600.0, 61200.0, 64800.0, 
     .             68400.0, 72000.0, 75600.0 /),
        ! H and LE specifications
        ! Every hour (3 half-hour point weighted average)
     .  H   = (/ 415.3, 388.3, 404.0, 
     .           326.8, 236.9, 220.1 /),
     .  LE  = (/ 121.2, 119.1, 122.9, 
     .           102.7, 104.8,  87.1 /),
        ! CO2 specifications
        ! Every hour (3 half-hour point weighted average)
     .  CO2 = (/ -12.4, -11.3, -10.4,
     .            -8.0,  -7.2,  -7.1 /)

        ! Input variable
        real, intent(in) :: time ! Current time [s]

        ! Output variables
        real, intent(out) ::
     .  heat_flx,       ! Heat flux             [W/m^2]
     .  moisture_flx    ! Moisture flux         [W/m^2]

        ! Output variables (optional)
        real, optional, intent(out) :: 
     .  CO2_flx  ! Carbon dioxide flux  [umol/(m^2 s)]

        ! Local variables
        integer :: i1, i2

        real :: a

        if ( time <= times(1) ) then
           heat_flx     = H(1)
           moisture_flx = LE(1)

        else if ( time >= times(ntimes) ) then
           heat_flx     = H(ntimes)
           moisture_flx = LE(ntimes)

        else  ! time > times(1) and time < times(ntimes)
          i1 = 1

          do while ( i1 <= ntimes-1 )
            i2 = i1 + 1

            if ( time >= times(i1) .and. time < times(i2) ) then
              a = (time-times(i1))/(times(i2)-times(i1))

              heat_flx     = ( 1. - a ) * H(i1) + a * H(i2)
              moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)

              if ( present( CO2_flx ) ) then
                CO2_flx = ( 1. - a ) * CO2(i1) + a * CO2(i2)
              end if ! present

              i1           = ntimes
            end if ! time >= times(i1) & time < times(i2)

            i1 = i2
          end do ! while i1 <= ntimes-1

        end if ! else 

        return
        end subroutine cobra_sfcflx

        end module cobra
