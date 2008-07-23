!$Id: gabls2.F90,v 1.3 2008-07-23 17:38:08 faschinj Exp $
!----------------------------------------------------------------------
        module gabls2

!       Description:
!       Contains subroutines for the GABLS2 case.
!----------------------------------------------------------------------

        implicit none

        public :: gabls2_tndcy, gabls2_sfclyr

        private ! Default Scope

        contains

!----------------------------------------------------------------------
        subroutine gabls2_tndcy & 
        ( time, time_initial, & 
          rhot, rcm, kk_rain, wmt, wmm, & 
          thlm_forcing, rtm_forcing, radht, Ncm, & 
          sclrm_forcing )

!        Description:
!          Subroutine to apply case-specific forcings to GABLS2 case
!          (Michael Falk, 29 Dec 2006).
!
!        References:
!          http://people.su.se/~gsven/gabls/
!-----------------------------------------------------------------------

        use parameters, only: sclr_dim ! Variable(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use saturation, only: sat_mixrat_liq ! Procedure(s)

        use stats_precision, only: time_precision ! Variable(s)

        use array_index, only: iisclr_rt, iisclr_thl

        implicit none

        ! Input Variables
        real(kind=time_precision), intent(in) :: & 
          time,           & ! Current length of timestep      [s]
          time_initial   ! Current length of timestep      [s]

        real, dimension(gr%nnzp), intent(in) :: & 
          rhot,           & ! Air density on t-grid           [kg/m^3]
          rcm            ! Cloud water mixing ratio        [kg/kg]

        logical, intent(in) :: & 
          kk_rain       !  Logical variable-- are we using KK scheme?

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: & 
          wmt,          & ! Large-scale vertical motion on t grid   [m/s]
          wmm,          & ! Large-scale vertical motion on m grid   [m/s]
          thlm_forcing, & ! Large-scale thlm tendency               [K/s]
          rtm_forcing,  & ! Large-scale rtm tendency                [kg/kg/s]
          radht,        & ! dT/dt, then d Theta/dt, due to rad.     [K/s]
          Ncm          ! Number of cloud droplets                [#/kg]

        ! Output Variables (optional)
        real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
          sclrm_forcing ! Passive scalar LS tendency            [units/s]

        ! Local Variables, general
        integer :: k ! Loop index

        ! Compute vertical motion
        if (time > (time_initial + 93600.)) then ! That is, after 26 hours of model time; 
                                                 ! per GABLS2 specification
          do k=1,gr%nnzp
            if (gr%zt(k) <= 1000) then
              wmt(k) = -0.005 * (gr%zt(k) / 1000)
            else
              wmt(k) = -0.005
            end if
          end do
        else
          do k=1,gr%nnzp
            wmt(k) = 0.
          end do
        end if

        wmm = zt2zm( wmt )

        ! Boundary conditions on vertical motion.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top


        ! Compute large-scale horizontal temperature advection
        do k=1,gr%nnzp
          thlm_forcing(k) = 0.
          radht(k) = 0.
        end do


        ! Compute large-scale horizontal moisture advection [g/kg/s]
        do k=1,gr%nnzp
          rtm_forcing(k) = 0.
        end do

        ! Compute initial cloud droplet concentration
        if ( kk_rain ) then
          do k=1,gr%nnzp
            if ( rcm(k) >= 1.e-6 ) then
              ! Ncm is in units of kg^-1.  If the coefficient is in m^-3, then
              ! it needs to be divided by rhot in order to get units of kg^-1.
              ! Brian.  Sept. 8, 2007.
!              Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6
!     .                 * rhot(k)
              Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6 & 
                       / rhot(k)
            else
              Ncm(k) = 0.
            end if
          end do
        end if

        ! Test scalars with thetal and rt if desired
        if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
        if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

        return
        end subroutine gabls2_tndcy
!----------------------------------------------------------------------




!----------------------------------------------------------------------
        subroutine gabls2_sfclyr & 
                   ( time, time_initial, lowest_level, psfc, & 
                     um, vm, thlm, rtm, & 
                     upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc,  & 
                     ustar, & 
                     wpsclrp_sfc, wpedsclrp_sfc )

!----------------------------------------------------------------------
!        Description:
!          Surface forcing subroutine for GABLS2 case.  Written
!          29 December 2006 by Michael Falk.
!
!        References:
!          http://people.su.se/~gsven/gabls/
!-----------------------------------------------------------------------

          use constants, only: Cp, Rd, p0, kappa, grav ! Variable(s)

          use parameters, only: sclr_dim ! Variable(s)

          use saturation, only: sat_mixrat_liq ! Procedure(s)

          use stats_precision, only: time_precision ! Variable(s)

          use diag_ustar_mod, only: diag_ustar ! Variable(s)
          
          use array_index, only: iisclr_rt, iisclr_thl

          implicit none

          ! Local constants
          real, parameter ::     & 
            ubmin   = 0.25,      & ! Minimum value for ubar.
            C_10    = 0.0013,    & ! Drag coefficient, defined by ATEX specification
            z0      = 0.03         ! Roughness length, defined by GABLS2 specification

          ! Input variables
          real(kind=time_precision), intent(in) :: & 
            time,                & ! Time elapsed since 0.0 s          [s]
            time_initial           ! Initial time of model integration [s]

          real, intent(in) ::    & 
            psfc,                & ! Surface pressure [Pa]
            lowest_level,        & ! Height of lowest above-ground gridpoint [m]
            um,                  & ! u at the lowest above-ground model level.  [m/s]
            vm,                  & ! v at the lowest above-ground model level.  [m/s]
            thlm,                & ! theta-l at the lowest above-ground model level. 
                                   ! (theta = theta-l because there's no liquid in this case)  [K]
            rtm                    ! rt at the lowest above-ground model level.  [kg/kg]
          ! Output variables
          real, intent(out) :: & 
            upwp_sfc,     & ! The turbulent upward flux of u-momentum         [(m/s)^2]
            vpwp_sfc,     & ! The turbulent upward flux of v-momentum         [(m/s)^2]
            wpthlp_sfc,   & ! The turbulent upward flux of theta-l            [K m/s]
            wprtp_sfc,    & ! The turbulent upward flux of rtm (total water)  [kg/kg m/s]
            ustar           ! surface friction velocity                       [m/s]

          ! Output variables (optional)
          real, optional, intent(out), dimension(sclr_dim) :: & 
            wpsclrp_sfc,    & ! The upward flux of the scalars       [units m/s]
            wpedsclrp_sfc     ! The upward flux of the eddy-scalars  [units m/s]

          ! Local variables
          real :: & 
            ubar,                & ! Root (u^2 + v^2), per ATEX and RICO spec.
!     .      ustar,              ! Friction velocity, computed from diag_ustar.
            Cz,                  & ! C_10 scaled to the height of the lowest 
                                   ! model level. (Per ATEX spec)
            time_in_hours,       & ! time in hours from 00 local on first day of experiment 
                                   ! (experiment starts at 14)
            sst,                 & ! Sea surface temperature [K].
            sstheta,             & ! Sea surface potential temperature [K].
            bflx                   ! Needed for diag_ustar; equal to wpthlp_sfc * (g/theta)

          ! Define variable values
          ubar = max( ubmin, sqrt( um*um + vm*vm ) )
          Cz   = C_10 * ((log( 10/z0 ))/(log( lowest_level/z0 ))) * & 
                 ((log( 10/z0 ))/(log( lowest_level/z0 ))) ! Modification in case 
                                                           ! lowest model level isn't at 10 m,
                                                           ! from ATEX specification
          time_in_hours = real((time - time_initial) / 3600. + 14.) ! at initial time, 
                                                                    ! time_in_hours = 14 
                                                                    ! (14 local; 19 UTC)


          ! Compute sea surface temperature
          if (time_in_hours <= 17.4) then
            sst = -10 - (25*cos(time_in_hours*0.22 + 0.2)) ! SST in celsius per GABLS2 spec
          else if (time_in_hours <= 30.0) then
            sst = (-0.54 * time_in_hours) + 15.2
          else if (time_in_hours <= 41.9) then
            sst = -7 - (25*cos(time_in_hours*0.21 + 1.8))
          else if (time_in_hours <= 53.3) then
            sst = (-0.37 * time_in_hours) + 18.0
          else if (time_in_hours <= 65.6) then
            sst = -4 - (25*cos(time_in_hours*0.22 + 2.5))
          else
            sst = 4.4
          end if

          sst     = sst + 273.15
          sstheta = sst * ((p0 / psfc)**(Rd/Cp))

          ! Compute heat and moisture fluxes
          wpthlp_sfc  = -Cz * ubar * ( thlm - sst * (p0/psfc)**kappa ) ! [K * m/s]
          wprtp_sfc  = -Cz * ubar *  & 
                        ( rtm - sat_mixrat_liq(psfc,sst) ) * 0.025   ! [kg/kg * m/s]
                                                                     ! 2.5% factor from 
                                                                     ! GABLS2 specification
          ! Compute momentum fluxes
          bflx  = wpthlp_sfc * grav / sstheta
          ustar = diag_ustar(lowest_level,bflx,ubar,z0)
          upwp_sfc    = -um * ustar*ustar / ubar                       ! [(m/s)^2]
          vpwp_sfc    = -vm * ustar*ustar / ubar                       ! [(m/s)^2]

          ! Let passive scalars be equal to rt and theta_l for now
          if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
          if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

          if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
          if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

        return
        end subroutine gabls2_sfclyr

        end module gabls2
