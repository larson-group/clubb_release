#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: arm_97.F,v 1.1 2007-07-19 19:06:44 dschanen Exp $
        module arm_97

!       Description:
!       Contains subroutines for the July 27 1997 ARM case.
!----------------------------------------------------------------------

        public :: arm_97_tndcy, arm_97_sfclyr

        private

        contains

!----------------------------------------------------------------------
        subroutine arm_97_tndcy( time, wmt, wmm, 
     .                           thlm_forcing, rtm_forcing,
     .                           sclrm_forcing )
!       Description:
!       Subroutine to set thetal and total water tendencies for ARM 97 case

!       References:
!       None
!----------------------------------------------------------------------

        use grid_class
        use constants, only: sclr_dim
        use zlinterp, only: zlinterp_fnc

        implicit none

        ! External
        intrinsic :: present

        ! Constant Parameters
        integer, parameter :: ntimes = 2

        real, dimension(ntimes), parameter ::
     .  times = (/ 0., 10800./)   ! Time index          [s]

        real, parameter, dimension(18) :: 
     .  z1 =
     . (/ 59.210,   524.304,  1008.614,  1515.234,  2047.661,  2609.434,
     .  3203.960,  3835.436,  4509.983,  5235.930,  6023.204,  6882.909,
     .  7828.958,  8881.194, 10071.930, 11458.434, 13147.343, 15385.763
     .  /),
     .  thl_ls1 =
     .  (/ -0.11645E-04, 0.12399E-04, 0.43716E-04, 
     .      0.71247E-04, 0.74449E-04, 0.59121E-04, 
     .      0.35906E-04, 0.12283E-04,-0.12166E-04,
     .     -0.29332E-04,-0.33809E-04,-0.25980E-04,
     .     -0.20418E-04,-0.64505E-05, 0.31542E-04,
     .      0.54968E-04, 0.25066E-04, 0.17456E-05 /),
     .  rt_ls1 = 
     .  (/ -0.99861E-08,-0.12141E-07,-0.27897E-07,
     .     -0.51029E-07,-0.60910E-07,-0.55165E-07,
     .     -0.43801E-07,-0.38822E-07,-0.30447E-07,
     .     -0.18448E-07,-0.59411E-08,-0.10767E-08,
     .     -0.62347E-09,-0.20151E-08,-0.12447E-08,
     .     -0.40283E-09,-0.62638E-10,-0.11792E-10 /)

        real, parameter, dimension(18) :: 
     .  z2 =
     .  (/   55.768,  520.329, 1005.165, 1512.854, 2046.478, 2609.350,
     .     3204.743, 3836.836, 4511.717, 5237.574, 6024.259, 6883.026,
     .     7828.236, 8880.256,10071.867,11460.698,13153.004,15395.084/),
     .  thl_ls2 =
     .  (/ 0.25600E-08, 0.23874E-04, 0.39514E-04,
     .     0.45791E-04, 0.46159E-04, 0.41388E-04,
     .     0.33821E-04, 0.30719E-04, 0.34638E-04,
     .     0.32793E-04, 0.24066E-04, 0.22069E-04,
     .     0.22446E-04, 0.40732E-04, 0.68543E-04,
     .     0.70464E-04, 0.30210E-04, 0.38782E-05/),
     .  rt_ls2 =
     .  (/ 0.38979E-08, 0.13248E-07, 0.13695E-07,
     .    -0.82704E-08,-0.29288E-07,-0.41506E-07,
     .    -0.43908E-07,-0.43992E-07,-0.38316E-07,
     .    -0.25748E-07,-0.13755E-07,-0.75518E-08,
     .    -0.41940E-08,-0.29164E-08,-0.10718E-08,
     .    -0.27600E-09,-0.77380E-11,-0.88493E-11/)


        ! Input Variables
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  thlm_forcing, ! Liquid water potential temperature tendency     [K/s]
     .  wmm,          ! Vertical velocity on moment. grid               [m/s]
     .  wmt,          ! Vertical velocity on thermo. grid               [m/s]
     .  rtm_forcing   ! Total water mixing ratio tendency               [kg/kg/s]

        ! Output (optional) Variables
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) :: 
     .  sclrm_forcing ! Passive scalar tendency [units vary]

        real :: a
        integer :: i1, i2

!-----------------------------------------------------------------------

        ! Thetal forcing is equal to the LS tendency given here and the
        ! interactive calculation done in BUGSrad

!       if ( time == times(1) ) then
!         a = 0.0
!       else if ( time >= times(ntimes) ) then
!         a = 1.0
!       else
!         i1 = 1
!         do while ( i1 <= ntimes-1 )
!           i2 = i1 + 1
!           if ( time >= times(i1) .and. time < times(i2) ) then
!             a            = (time-times(i1))/(times(i2)-times(i1))
!             heat_flx     = ( 1. - a ) * SE(i1) + a * SE(i2)
!             moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
!             i1           = ntimes
!           end if
!             i1 = i2
!          end do
!        end if ! time <= times(1)

        i1 = 1
        i2 = 2

        a = (time-times(i1))/(times(i2)-times(i1))

        ! Interpolate LS thetal tendency to the HOC grid
        thlm_forcing(:) = 
     .  (1.-a)*zlinterp_fnc( gr%nnzp, 18, gr%zt(:), z1(:), thl_ls1(:) )
     .  + a*zlinterp_fnc( gr%nnzp, 18, gr%zt(:), z2(:), thl_ls2(:) )

        ! Interpolate LS rt tendency to the HOC grid
        rtm_forcing(:) = 
     .  (1.-a)*zlinterp_fnc( gr%nnzp, 18, gr%zt(:), z1(:), rt_ls1(:) )
     .  + a*zlinterp_fnc( gr%nnzp, 18, gr%zt(:), z2(:), rt_ls2(:) )

        ! Zero imposed subsidence (always)
        wmm(:) = 0.
        wmt(:) = 0.

        if ( present( sclrm_forcing ) ) then 
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine arm_97_tndcy
!----------------------------------------------------------------------
        subroutine arm_97_sfclyr
     .             ( time, z, rho0,
     .               thlm_sfc, um_sfc, vm_sfc, 
     .               upwp_sfc, vpwp_sfc,
     .               wpthlp_sfc, wprtp_sfc,
     .               wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

        use constants, only: Cp, Lv, grav, sclr_dim

        implicit none

        ! External
        real, external :: diag_ustar

        intrinsic :: max, sqrt, present

         ! Constant Parameters
        integer, parameter :: ntimes = 2

        real, parameter :: 
     .  ubmin = 0.25,   ! Minimum value for ubar 
     .  z0    = 0.035   ! ARM Cu mom. roughness height

        real, dimension(ntimes), parameter ::
     .  times = (/ 0., 10800./)   ! Time index          [s]

        real, dimension(ntimes), parameter :: 
     .  SE = (/ 51.482, 13.389/), ! Sensible heat flux  [W/m^2]
     .  LE = (/143.459, 63.825/)  ! Evaporation         [W/m^2]

        ! Input Variables
        real, intent(in) :: 
     .  time,     ! Current time        [s] 
     .  z,        ! Height at zt=2      [s] 
     .  rho0,     ! Density at zm=1     [kg/m^3] 
     .  um_sfc,   ! um at (2)           [m/s]
     .  vm_sfc,   ! vm at (2)           [m/s]
     .  thlm_sfc  ! thlm at (2)         [m/s]

        ! Output variables
        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc    ! w'r_t'(1) at (1) [(m kg)/(s kg)]

        ! Output variables (optional)
        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! Passive scalar surface flux      [units m/s] 
     .  wpedsclrp_sfc ! Passive eddy-scalar surface flux [units m/s]

        ! Local variables
        real :: ubar, ustar, bflx, heat_flx, moisture_flx, a
        integer :: i1, i2

        if ( time <= times(1) ) then
          heat_flx     = SE(1)
          moisture_flx = LE(1)
        else if ( time >= times(ntimes) ) then
          heat_flx     = SE(ntimes)
          moisture_flx = LE(ntimes)
        else
          i1 = 1
          do while ( i1 <= ntimes-1 )
            i2 = i1 + 1
            if ( time >= times(i1) .and. time < times(i2) ) then
              a            = (time-times(i1))/(times(i2)-times(i1))
              heat_flx     = ( 1. - a ) * SE(i1) + a * SE(i2)
              moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
              i1           = ntimes
            end if
              i1 = i2
           end do
         end if ! time <= times(1)

        ! Convert W/m^2 into w'thl' w'rt' units
        wpthlp_sfc = heat_flx / ( Cp * rho0 )     ! (K m/s)
        wprtp_sfc  = moisture_flx / ( Lv * rho0 ) ! (kg m/ kg s)

        ! Let passive scalars be equal to rt and theta_l for now 
        if ( present( wpsclrp_sfc ) ) then 
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        if ( present( wpsclrp_sfc ) ) then 
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        ! Compute momentum fluxes using ARM formulae

        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        bflx = grav/thlm_sfc * wpthlp_sfc

        ! Compute ustar
        ustar = diag_ustar( z, bflx, ubar, z0 )

        upwp_sfc = -um_sfc * ustar**2 / ubar
        vpwp_sfc = -vm_sfc * ustar**2 / ubar

        return
        end subroutine arm_97_sfclyr

        end module arm_97
