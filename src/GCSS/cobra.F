!----------------------------------------------------------------------
! $Id: cobra.F,v 1.14 2008-07-08 04:44:10 dschanen Exp $
        module cobra
!       Description:
!       Contains subroutines for the COBRA CO2 case.
!----------------------------------------------------------------------

        implicit none

        private ! Default Scope

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

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use constants, only: fstderr ! Variable(s)

        use parameters, only: sclr_dim ! Variable(s)

        use stats_precision, only: time_precision ! Variable(s)

        use array_index, only: 
     .      iisclr_thl, iisclr_rt, iiCO2 ! Variable(s)

        implicit none

        ! Input Variables
        real(kind=time_precision), intent(in) :: 
     .  time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid                 [m/s]
     .  wmm,         ! w wind on momentum grid                      [m/s]
     .  thlm_forcing,! Liquid water potential temperature tendency  [K/s]
     .  rtm_forcing  ! Total water mixing ratio tendency            [kg/kg/s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing ! Passive scalar tendency [units vary]

        ! Local Variables
        integer :: k

        ! Large-scale subsidence

        DO k = 2, gr%nnzp, 1

!  Vince Larson tapered the subsidence to zero to avoid 
!     constant thlm and instability aloft. 
!  6 Mar 2007
!          wmt(k) = -5.0E-6 * gr%zt(k)

           if ( gr%zt(k) >= 0. .and. gr%zt(k) < 3000. ) then
              wmt(k) = -5.0E-6 * gr%zt(k)
           else if ( gr%zt(k) >= 3000. ) then
              wmt(k)
     .          = - 0.0150
     .            + 0.0150 * ( gr%zt(k) - 3000. ) / ( 4000. - 3000. )
           else
              write(fstderr,*) "cobra_tndcy:"
     .           //" error in subsidence profile."
              write(fstderr,*) 'Altitude gr%zt = ', gr%zt(k)
              stop
           end if
! End Vince Larson's change

        END DO

        ! Boundary condition.
        wmt(1) = 0.0        ! Below surface

        ! Interpolation
        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        ! No large-scale water tendency or cooling
        rtm_forcing  = 0.0
        thlm_forcing = 0.0

        ! Setup passive scalars, if they're enabled
        if ( iiCO2 > 0 ) sclrm_forcing(:,iiCO2) = 0.0
        if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
        if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt) = rtm_forcing

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

        use constants, only: Cp, Lv, grav ! Variable(s)

        use parameters, only: sclr_dim ! Variable(s)

        use stats_precision, only: time_precision ! Variable(s)

        use diag_ustar_mod, only: diag_ustar ! Variable(s)

        use array_index, only: iisclr_rt, iisclr_thl, iiCO2
        
        implicit none

        intrinsic :: sqrt, max

        ! Parameter Constants
        real, parameter ::
     .  ubmin = 0.25,
     .  M_da  = 0.02897  ! Molecular weight of dry air.

        ! Input variables
        real(kind=time_precision), intent(in) :: 
     .  time      ! Current time                [s]

        real, intent(in) ::
     .  z,        ! Elevation at zt=2           [m]
     .  dn0,      ! Air density at surface      [kg/m^3]
     .  thlm_sfc, ! Theta_l at zt(2)            [K]
     .  um_sfc,   ! u wind at zt(2)             [m/s]
     .  vm_sfc    ! v wind at zt(2)             [m/s]

        ! Output variables
        real, intent(out) :: 
     .  upwp_sfc,   ! u'w' at surface           [m^2/s^2]
     .  vpwp_sfc,   ! v'w' at surface           [m^2/s^2]
     .  wpthlp_sfc, ! w'theta_l' surface flux   [(m K)/s]
     .  wprtp_sfc,  ! w'rt' surface flux        [(m kg)/(kg s)]
     .  ustar       ! surface friction velocity [m/s]

        ! Output variables
        real, intent(out), dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! w'sclr' surface flux          [units m/s]
     .  wpedsclrp_sfc ! w' edsclr' surface flux       [units m/s]

        ! Local variables
        real :: 
     .  usfc, vsfc, ubar,
     .  true_time,
     .  heat_flx, moisture_flx,
     .  heat_flx2, moisture_flx2,
     .  bflx

        real :: CO2_flx

        real :: CO2_flx2


!       COBRA roughness height
!        real, parameter :: z0 = 0.035  ! ARM momentum roughness height
        real, parameter :: z0 = 1.75   ! momentum roughness height

! Compute heat and moisture fluxes from ARM data in (W/m2)

        CO2_flx2 = 0.0 ! Default initialization 

        true_time = time
        call cobra_sfcflx( true_time, heat_flx, moisture_flx,
     .                    CO2_flx
     .                   )

        ! Compute momentum fluxes

        ! Convert heat_flx and moisture_flx to natural units
        heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
        moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

!       Convert CO2 surface flux to natural units.
!       The CO2 flux has been given in units of:  umol/(m^2 s).
!       umol stands for micromoles.  The CO2 concentration in
!       this code is in units of ppmv, which is also the molar
!       mixing ratio times 10^6.
!       The units are:  10^6 * [ mol (CO2) / mol (dry air) ].
!       w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
!       where M is the molecular weight of dry air.
        CO2_flx2 = CO2_flx * ( M_da / dn0 )

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

        if ( iiCO2 > 0 ) wpsclrp_sfc(iiCO2) = CO2_flx2
        if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
        if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

        if ( iiCO2 > 0 ) wpedsclrp_sfc(iiCO2) = CO2_flx2
        if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
        if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

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

        ! Constant parameters
        integer, parameter :: ntimes = 49

        real, parameter, dimension(ntimes) ::
     .  times = (/  57600.0,  61200.0,  64800.0,  68400.0,  ! Hrs 0, 1, 2, 3
     .              72000.0,  75600.0,  79200.0,  82800.0,  ! Hrs 4, 5, 6, 7
     .              86400.0,  90000.0,  93600.0,  97200.0,  ! Hrs 8, 9, 10, 11
     .             100800.0, 104400.0, 108000.0, 111600.0,  ! Hrs 12, 13, 14, 15
     .             115200.0, 118800.0, 122400.0, 126000.0,  ! Hrs 16, 17, 18, 19
     .             129600.0, 133200.0, 136800.0, 140400.0,  ! Hrs 20, 21, 22, 23
     .             144000.0, 147600.0, 151200.0, 154800.0,  ! Hrs 24, 25, 26, 27
     .             158400.0, 162000.0, 165600.0, 169200.0,  ! Hrs 28, 29, 30, 31
     .             172800.0, 176400.0, 180000.0, 183600.0,  ! Hrs 32, 33, 34, 35
     .             187200.0, 190800.0, 194400.0, 198000.0,  ! Hrs 36, 37, 38, 39
     .             201600.0, 205200.0, 208800.0, 212400.0,  ! Hrs 40, 41, 42, 43
     .             216000.0, 219600.0, 223200.0, 226800.0,  ! Hrs 44, 45, 46, 47
     .             230400.0 /),                             ! Hr  48

        ! H = sensible heat flux at surface   [W m^{-2}]
        ! LE = latent heat flux at surface    [W m^{-2}]
        ! Every hour (3 half-hour point weighted average)
        ! Note:  For H and LE, the "***" symbol next to certain rows indicates 
        !        that the values were missing from the data listed in those 
        !        rows.  The values listed here were interpolated from what was 
        !        provided in order to follow the general trend.

! Vince Larson changed fluxes to diurnal composite fluxes over Jun 2004
!    in order to avoid large, negative surface sensible heat fluxes at night.
!    29 Mar 2007

! Fluxes for two days following COBRA aircraft measurements at 16 GMT on 11 Jun 2004
!     .  H   = (/ 415.3, 388.3, 404.0, 326.8, ! Hrs 0, 1, 2, 3
!     .           236.9, 220.1, 128.9,  31.4, ! Hrs 4, 5, 6, 7
!     .           -40.1, -64.4, -51.9, -38.1, ! Hrs 8, 9, 10, 11
!     .           -50.2, -45.7, -66.3, -78.1, ! Hrs 12, 13, 14, 15
!     .           -61.7, -58.8, -30.8,  46.3, ! Hrs 16, 17, 18, 19
!     .           148.2, 258.4, 362.4, 349.6, ! Hrs 20, 21, 22, 23
!     .           359.3, 429.4, 418.9, 364.9, ! Hrs 24, 25, 26, 27
!     .           242.8, 192.0,  83.1,  18.0, ! Hrs 28, 29, 30, 31
!     .             5.0,  -5.0,  -8.0, -11.0, ! Hrs 32, 33, 34, 35 ***
!     .           -10.0,  -9.0,  -6.0,  -3.0, ! Hrs 36, 37, 38, 39 ***
!     .             0.0,   5.0,  15.0,  30.0, ! Hrs 40, 41, 42, 43 ***
!     .           154.5, 268.7, 263.7, 265.8, ! Hrs 44, 45, 46, 47
!     .           361.5 /),                   ! Hr  48
!     .  LE  = (/ 121.2, 119.1, 122.9, 102.7, ! Hrs 0, 1, 2, 3
!     .           104.8,  87.1,  80.8,  77.9, ! Hrs 4, 5, 6, 7
!     .            23.6,  -1.3,  -3.8,   4.0, ! Hrs 8, 9, 10, 11
!     .            -0.1,  -4.1,   3.1,  -0.8, ! Hrs 12, 13, 14, 15
!     .            -6.0,   1.9,   3.7,  15.7, ! Hrs 16, 17, 18, 19
!     .            49.3,  67.8,  88.4,  84.7, ! Hrs 20, 21, 22, 23
!     .           103.8, 137.0, 149.4, 138.2, ! Hrs 24, 25, 26, 27
!     .           108.1,  91.4,  69.6,  54.3, ! Hrs 28, 29, 30, 31
!     .            20.0,   8.0,   2.0,  -1.0, ! Hrs 32, 33, 34, 35 ***
!     .            -3.0,  -4.0,  -5.0,  -3.0, ! Hrs 36, 37, 38, 39 ***
!     .            -1.0,   3.0,   9.0,  18.0, ! Hrs 40, 41, 42, 43 ***
!     .            37.3,  94.8, 120.6, 132.8, ! Hrs 44, 45, 46, 47
!     .           120.4 /),                   ! Hr  48
        ! CO2 specifications
        ! Every hour (3 half-hour point weighted average)
!     .  CO2 = (/ -12.4, -11.3, -10.4,  -8.0, ! Hrs 0, 1, 2, 3
!     .            -7.2,  -7.1,  -5.3,  -2.6, ! Hrs 4, 5, 6, 7
!     .             1.7,   5.2,   2.7,   3.8, ! Hrs 8, 9, 10, 11
!     .             4.4,   4.6,   4.4,   3.8, ! Hrs 12, 13, 14, 15
!     .             2.8,   2.5,   0.5,  -4.0, ! Hrs 16, 17, 18, 19
!     .            -7.4, -10.6, -13.2, -11.2, ! Hrs 20, 21, 22, 23
!     .           -12.1, -13.9, -13.3, -11.3, ! Hrs 24, 25, 26, 27
!     .            -8.4,  -5.9,  -3.9,  -2.5, ! Hrs 28, 29, 30, 31
!     .             0.4,   4.0,   4.4,   4.4, ! Hrs 32, 33, 34, 35
!     .             4.4,   4.4,   4.5,   4.5, ! Hrs 36, 37, 38, 39
!     .             4.5,   3.4,  -0.6,  -4.8, ! Hrs 40, 41, 42, 43
!     .            -7.3, -19.4, -16.4, -13.1, ! Hrs 44, 45, 46, 47
!     .           -11.2 /)                    ! Hr  48

! Initial point corresponds to 16 GMT = noon EDT.
! Each point is a one-hour average.
     .  H   = (/ 257.5, 274.0, 234.5, 195.0, ! Hrs 0, 1, 2, 3
     .           165.8, 117.0,  61.7,   4.9, ! Hrs 4, 5, 6, 7
     .           -24.1, -26.6, -25.3, -30.4, ! Hrs 8, 9, 10, 11
     .           -32.3, -31.4, -34.5, -16.8, ! Hrs 12, 13, 14, 15
     .           -31.6, -33.4, -13.0,  38.4, ! Hrs 16, 17, 18, 19
     .           111.1, 156.4, 205.8, 218.6, ! Hrs 20, 21, 22, 23
     .           257.5, 274.0, 234.5, 195.0, ! Hrs 24, 25, 26, 27
     .           165.8, 117.0,  61.7,   4.9, ! Hrs 28, 29, 30, 31
     .           -24.1, -26.6, -25.3, -30.4, ! Hrs 32, 33, 34, 35 ***
     .           -32.3, -31.4, -34.5, -16.8, ! Hrs 36, 37, 38, 39 ***
     .           -31.6, -33.4, -13.0,  38.4, ! Hrs 40, 41, 42, 43 ***
     .           111.1, 156.4, 205.8, 218.6, ! Hrs 44, 45, 46, 47
     .           257.5 /),                   ! Hr  48
     .  LE  = (/ 116.8, 117.8, 116.7, 105.5, ! Hrs 0, 1, 2, 3
     .            92.8,  79.8,  58.4,  32.4, ! Hrs 4, 5, 6, 7
     .             8.0,   3.5,   2.5,   5.1, ! Hrs 8, 9, 10, 11
     .             1.1,   1.6,   0.5,  -0.5, ! Hrs 12, 13, 14, 15
     .            -1.4,   1.4,   0.2,  11.0, ! Hrs 16, 17, 18, 19
     .            32.1,  48.2,  74.9,  96.0, ! Hrs 20, 21, 22, 23
     .           116.8, 117.8, 116.7, 105.5, ! Hrs 24, 25, 26, 27
     .            92.8,  79.8,  58.4,  32.4, ! Hrs 28, 29, 30, 31
     .             8.0,   3.5,   2.5,   5.1, ! Hrs 32, 33, 34, 35 ***
     .             1.1,   1.6,   0.5,  -0.5, ! Hrs 36, 37, 38, 39 ***
     .            -1.4,   1.4,   0.2,  11.0, ! Hrs 40, 41, 42, 43 ***
     .            32.1,  48.2,  74.9,  96.0, ! Hrs 44, 45, 46, 47
     .           116.8 /),                   ! Hr  48
        ! CO2 = CO2 hourly averaged surface flux      [umol/(m^2 s)]
     .  CO2 = (/ -12.0, -12.1, -10.8,  -9.3, ! Hrs 0, 1, 2, 3
     .            -8.6,  -6.3,  -3.8,  -0.9, ! Hrs 4, 5, 6, 7
     .             2.4,   3.2,   4.2,   4.5, ! Hrs 8, 9, 10, 11
     .             4.5,   4.2,   4.3,   4.5, ! Hrs 12, 13, 14, 15
     .             4.5,   3.8,   1.6,  -2.3, ! Hrs 16, 17, 18, 19
     .            -6.5, -10.2, -10.4, -11.5, ! Hrs 20, 21, 22, 23
     .           -12.0, -12.1, -10.8,  -9.3, ! Hrs 24, 25, 26, 27
     .            -8.6,  -6.3,  -3.8,  -0.9, ! Hrs 28, 29, 30, 31
     .             2.4,   3.2,   4.2,   4.5, ! Hrs 32, 33, 34, 35
     .             4.5,   4.2,   4.3,   4.5, ! Hrs 36, 37, 38, 39
     .             4.5,   3.8,   1.6,  -2.3, ! Hrs 40, 41, 42, 43
     .            -6.5, -10.2, -10.4, -11.5, ! Hrs 44, 45, 46, 47
     .           -12.0 /)                    ! Hr  48
! End of Vince Larson's change.

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

              CO2_flx = ( 1. - a ) * CO2(i1) + a * CO2(i2)

              i1           = ntimes
            end if ! time >= times(i1) & time < times(i2)

            i1 = i2
          end do ! while i1 <= ntimes-1

        end if ! else 

        return
        end subroutine cobra_sfcflx

        end module cobra
