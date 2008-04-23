#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: arm_97.F,v 1.15 2008-04-23 18:09:57 faschinj Exp $
        module arm_97

!       Description:
!       Contains subroutines for the July 26-30 1997 ARM IOP A case.
!----------------------------------------------------------------------

        public :: arm_97_tndcy, arm_97_sfclyr, arm_97_init

        private
        ! Constant Parameters
        integer, parameter :: ntimes = 169, nz = 18,
     .   per_line = 5
        ! File path for the forcing files
        character(*), parameter ::
     .   file_path = '../model/arm_97_forcings/'
        
        real, dimension(ntimes) :: times       ! Time from day0      [s]
        real, dimension(nz, ntimes) :: z       ! Height              [m]
        real, dimension(nz, ntimes) :: thl_ls  ! Potential Temperature
                                               ! Tendency            [K/s]
        real, dimension(nz, ntimes) :: rt_ls   ! Water Vapor Advective
                                               ! Tendency            [Kg/Kg/s]
        real, dimension(nz, ntimes) :: um_obs  ! Obs. wind u         [m/s]
        real, dimension(nz, ntimes) :: vm_obs  ! Obs. wind v         [m/s]
        real, dimension(ntimes) :: SE          ! Sensible heat flux  [W/m^2]
        real, dimension(ntimes) :: LE          ! Evaporation         [W/m^2] 


        contains

!----------------------------------------------------------------------
        subroutine arm_97_tndcy( time, wmt, wmm, 
     .                           thlm_forcing, rtm_forcing,
!              Joshua Fasching Oct 27 2007 
     .                           um_hoc_grid, vm_hoc_grid,
     .                           sclrm_forcing )
!       Description:
!       Subroutine to set thetal and total water tendencies for ARM 97 case

!       References:
!       None
!----------------------------------------------------------------------

        use grid_class
        use constants, only: sclr_dim
        use interpolation, only: zlinterp_fnc
        use stats_prec, only: time_prec

        implicit none

        ! External
        intrinsic :: present


        ! Input Variables
        real(kind=time_prec), intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  thlm_forcing, ! Liquid water potential temperature tendency  [K/s]
     .  wmm,          ! Vertical velocity on moment. grid            [m/s]
     .  wmt,          ! Vertical velocity on thermo. grid            [m/s]
     .  rtm_forcing,  ! Total water mixing ratio tendency            [kg/kg/s]
     .  um_hoc_grid,  ! Observed wind, for nudging                   [m/s]
     .  vm_hoc_grid   ! Observed wind, for nudging                   [m/s]

        ! Output (optional) Variables
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) :: 
     .  sclrm_forcing ! Passive scalar tendency [units vary]

        real :: time_frac
        integer :: i1, i2

        real, dimension(nz) :: thlm_t_interp,
     .   rtm_t_interp, um_obs_t_interp, vm_obs_t_interp


!-----------------------------------------------------------------------

        ! Thetal forcing is equal to the LS tendency given here and the
        ! interactive calculation done in BUGSrad

        ! Compute time_frac, where 0 < time_frac < 1.
        !     time_frac measures the fractional time interval
        !     between the current time (time), and the previous
        !     time at which forcing data are available (times(i1))
        !     and the next time (times(i2)).
        if ( time == times(1) ) then
          time_frac = 0.0
          i1 = 1
          i2 = 2
        else if ( time >= times(ntimes) ) then
          time_frac = 1.0
          i1 = ntimes-1
          i2 = ntimes
        else
          i1 = 1
          do while ( i1 <= ntimes-1 )
            i2 = i1 + 1
            if ( time >= times(i1) .and. time < times(i2) ) then
              time_frac = (time-times(i1))/(times(i2)-times(i1))
              exit
            end if
              i1 = i2
          end do
        end if ! time <= times(1)
        
        ! Interpolate LS thetal tendency to the HOC grid
        ! Time
        thlm_t_interp = (1.-time_frac) * thl_ls(:,i1) + time_frac * 
     .     thl_ls(:,i2)
        ! Vertical
        thlm_forcing(1:gr%nnzp) = zlinterp_fnc
     .           ( gr%nnzp, nz, gr%zt(:), z(:,i1), thlm_t_interp )

        ! Interpolate LS rt tendency to the HOC grid
        ! Time
        rtm_t_interp = (1.-time_frac) * rt_ls(:,i1) + time_frac * 
     .     rt_ls(:,i2)
        ! Vertical
        rtm_forcing(1:gr%nnzp) = zlinterp_fnc
     .           ( gr%nnzp, nz, gr%zt(:), z(:,i1), rtm_t_interp )
        ! Modified by Joshua Fasching (October 2007)

        ! Interpolate um observed to the HOC grid
        ! Time
        um_obs_t_interp = (1.-time_frac) * um_obs(:,i1) + time_frac * 
     .     um_obs(:,i2)
        ! Vertical
        um_hoc_grid(1:gr%nnzp) = zlinterp_fnc
     .           ( gr%nnzp, nz, gr%zt(:), z(:,i1), um_obs_t_interp )
        ! Added by Joshua Fasching (October 27 2007)

        ! Interpolate vm observed to the HOC grid
        ! Time
        vm_obs_t_interp = (1.-time_frac) * vm_obs(:,i1) + time_frac * 
     .     vm_obs(:,i2)
        ! Vertical
        vm_hoc_grid(1:gr%nnzp) = zlinterp_fnc
     .           ( gr%nnzp, nz, gr%zt(:), z(:,i1), vm_obs_t_interp )
        ! Added by Joshua Fasching (October 27 2007)
        um_hoc_grid (1) = um_hoc_grid(2)
        vm_hoc_grid (1) = vm_hoc_grid(2)

        ! The specification calls for zero imposed subsidence (always)

        wmm(:) = 0.
        wmt(:) = 0.

        if ( present( sclrm_forcing ) ) then 
          ! Test scalars with thetal and rt if desired
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
     .               wpthlp_sfc, wprtp_sfc, ustar,
     .               wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

        use constants, only: Cp, Lv, grav, sclr_dim
        use stats_hoc, only: iustar, sfc, lstats_samp
        use stats_prec, only: time_prec

        implicit none

        ! External
        real, external :: diag_ustar

        intrinsic :: max, sqrt, present


        real, parameter :: 
     .  ubmin = 0.25,   ! Minimum value for ubar 
     .  z0    = 0.035   ! ARM Cu mom. roughness height


        ! Input Variables
        real(time_prec), intent(in) :: 
     .  time      ! Current time        [s] 

        real, intent(in) :: 
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
     .  wprtp_sfc,   ! w'r_t'(1) at (1) [(m kg)/(s kg)]
     .  ustar        ! surface friction velocity [m/s]

        ! Output variables (optional)
        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! Passive scalar surface flux      [units m/s] 
     .  wpedsclrp_sfc ! Passive eddy-scalar surface flux [units m/s]

        ! Local variables
!        real :: ubar, ustar, bflx, heat_flx, moisture_flx, time_frac
        real :: ubar, bflx, heat_flx, moisture_flx, time_frac
        integer :: i1, i2
!----------------------------------------------------------------------

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
              time_frac            = (time-times(i1))/(times(i2)
     .           - times(i1))
              heat_flx     = ( 1. - time_frac ) * SE(i1) + 
     .           time_frac * SE(i2)
              moisture_flx = ( 1. - time_frac ) * LE(i1) + 
     .           time_frac * LE(i2)
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

        ! Compute momentum fluxes using ARM Cu formulae

        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        bflx = grav/thlm_sfc * wpthlp_sfc

        ! Compute ustar
        ustar = diag_ustar( z, bflx, ubar, z0 )

#ifdef STATS
! Save ustar for statistics
        if ( lstats_samp .and. iustar > 0 ) then
          sfc%x(1,iustar) = sfc%x(1,iustar) + ustar
          sfc%n(1,iustar) = sfc%n(1,iustar) + 1
        end if
#endif
        upwp_sfc = -um_sfc * ustar**2 / ubar
        vpwp_sfc = -vm_sfc * ustar**2 / ubar

        return
        end subroutine arm_97_sfclyr
!----------------------------------------------------------------------
        subroutine arm_97_init()
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------------

           use file_functions, only: file_read_1d, file_read_2d
           implicit none


           call file_read_1d(10,
     .      file_path//'arm_97_times.dat',
     .      ntimes, per_line, times)

           call file_read_2d(10,
     .      file_path//'arm_97_heights.dat',
     .      nz, ntimes, per_line, z)

           call file_read_1d(10,
     .      file_path//'arm_97_LE.dat',
     .      ntimes, per_line, LE)

           call file_read_1d(10,
     .      file_path//'arm_97_SE.dat',
     .      ntimes, per_line, SE)

           call file_read_2d(10,
     .      file_path//'arm_97_dTdt.dat',
     .      nz, ntimes, per_line, thl_ls)

           call file_read_2d(10,
     .       file_path//'arm_97_dqdt.dat',
     .       nz, ntimes, per_line, rt_ls)

           call file_read_2d(10,
     .       file_path//'arm_97_um_obs.dat',
     .       nz, ntimes, per_line, um_obs)

           call file_read_2d(10,
     .       file_path//'arm_97_vm_obs.dat',
     .       nz, ntimes, per_line, vm_obs)

        end subroutine arm_97_init
!----------------------------------------------------------------------
        end module arm_97

