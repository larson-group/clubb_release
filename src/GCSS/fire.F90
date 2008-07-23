!----------------------------------------------------------------------
! $Id: fire.F90,v 1.2 2008-07-23 17:38:08 faschinj Exp $
        module fire

!       Description:
!       Contains subroutines for the GCSS FIRE case.
!----------------------------------------------------------------------

        implicit none

        public :: fire_tndcy, sfc_momentum_fluxes, sfc_thermo_fluxes

        private ! Default Scope

        contains

!----------------------------------------------------------------------
        subroutine fire_tndcy & 
                   ( rhot, rcm, exner,  & 
                     wmt, wmm, Frad, radht,  & 
                     thlm_forcing, rtm_forcing, & 
                     sclrm_forcing )
!       Description:
!       Subroutine to large-scale subsidence for FIRE case. Calls
!       cloud_rad for computing radiation

!       References:
!       None
!----------------------------------------------------------------------

        use parameters, only: sclr_dim ! Variable(s)

        use model_flags, only: lbugsrad  ! Variable(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use atex_cloud_rad, only: cloud_rad ! Procedure(s)

        use stats_precision, only: time_precision ! Variable(s)

        use array_index, only: iisclr_rt, iisclr_thl

#ifdef STATS
        use stats_type, only: stat_update_var ! Procedure(s)

        use stats_variables, only: zt, iradht_LW, lstats_samp ! Variable(s)
#endif

        implicit none

        ! Input Variables
        real, intent(in), dimension(gr%nnzp) :: & 
        rhot,  & ! Density                         [kg/m^3]
        rcm,   & ! Liquid water mixing ratio       [kg/kg]
        exner ! Exner function                  [-]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: & 
        wmt,          & ! w wind on thermodynamic grid     [m/s]
        wmm,          & ! w wind on momentum grid          [m/s]
        Frad,         & ! Radiative flux                   [W/m^2]
        radht,        & ! Radiative heating rate           [K/s]
        thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
        rtm_forcing  ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
        sclrm_forcing ! Passive scalar tendency [units vary]

!       Internal variables

        integer :: k

!       Large-scale subsidence

        do k = 2, gr%nnzp

           if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
              wmt(k) = - 5.e-6 * gr%zt(k)
           end if

        end do

        ! Boundary condition.
        wmt(1) = 0.0        ! Below surface

        wmm = zt2zm( wmt )

        ! Boundary condition.
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        ! Radiative theta-l tendency is computed interactively elsewhere

        thlm_forcing = 0.0

        ! Large scale advective moisture tendency

        rtm_forcing = 0.0

        ! Use cloud_rad to compute radiation
        if ( .not. lbugsrad ) then
          call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing )
        end if

#ifdef STATS
        if ( .not. lbugsrad .and. lstats_samp ) then
          call stat_update_var( iradht_LW, radht, zt )
        end if
#endif

        ! Test scalars with thetal and rt if desired
        if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
        if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

        return
        end subroutine fire_tndcy

!------------------------------------------------------------------------
        subroutine sfc_momentum_fluxes( u, v, upwp_sfc, vpwp_sfc,  & 
                                        ustar )

!       Description:
!       This subroutine computes surface momentum fluxes using aerodynamic
!       formulas.

!       References:
!       None
!------------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: sqrt

        ! Constant parameter
!        real, intent(out) :: 
!     .  ustar = 0.3

        ! Input variables
        real, intent(in) ::  & 
        u,  & ! u wind first level above ground    [m/s]
        v  ! v wind first level above ground    [m/s]

        ! Output Variables
        real, intent(out) ::  & 
        upwp_sfc, & ! sfc u momentum flux (m^2/s^2)
        vpwp_sfc, & ! sfc v momentum flux (m^2/s^2)
        ustar    ! surface friction velocity [m/s]

        ! Local Variables
        real :: M ! total wind speed above ground

        ! Declare the value of ustar
        ustar = 0.3

        ! Computes fluxes

        M = sqrt( u*u + v*v )
        upwp_sfc = - ustar*ustar * u / M
        vpwp_sfc = - ustar*ustar * v / M

        return
        end subroutine sfc_momentum_fluxes

!------------------------------------------------------------------------
        subroutine sfc_thermo_fluxes( u, v, Tsfc, psfc, thlair, rtair, & 
                                      wpthlp_sfc, wprtp_sfc, & 
                                      wpsclrp_sfc, & 
                                      wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of heat and moisture 
!       using aerodynamic formulas.

!       References:
!       None
!------------------------------------------------------------------------

        use constants, only: kappa, p0 ! Variable(s)

        use parameters, only: sclr_dim ! Variable(s)

        use saturation, only: sat_mixrat_liq ! Procedure(s)

        use array_index, only: iisclr_rt, iisclr_thl

        implicit none
        
        ! External
        intrinsic :: present, sqrt
        
        ! Parameter
        real, parameter :: C = 1.3e-3

        ! Input Variables
        real, intent(in) ::  & 
        u,       & ! u wind                        [m/s]
        v,       & ! u wind                        [m/s]
        Tsfc,    & ! Surface temperature           [K]
        psfc,    & ! Surface pressure              [Pa]
        thlair,  & ! theta_l at first model layer  [K]
        rtair   ! rt at first model layer       [kg/kg]


        ! Output Variables
        real, intent(out) ::  & 
        wpthlp_sfc, & ! surface thetal flux        [K m/s]
        wprtp_sfc  ! surface moisture flux      [kg/kg m/s]

        ! Output Variables (optional) 
        real, optional, intent(out), dimension(sclr_dim) ::  & 
        wpsclrp_sfc,    & ! scalar surface flux            [units m/s]
        wpedsclrp_sfc  ! eddy-scalar surface flux       [units m/s]

        ! Local Variables
        real :: M  ! Total wind speed above ground

        ! Compute fluxes
        M = sqrt( u*u + v*v )
        wpthlp_sfc = -C * M * ( thlair - Tsfc * (psfc/p0)**kappa )
        wprtp_sfc  = -C * M * ( rtair - sat_mixrat_liq( psfc, Tsfc ) )

        ! Let passive scalars be equal to rt and theta_l for now
        if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
        if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

        if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
        if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

        return
        end subroutine sfc_thermo_fluxes

        end module fire
