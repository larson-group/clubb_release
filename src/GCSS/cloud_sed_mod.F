!----------------------------------------------------------------------
!$Id: cloud_sed_mod.F,v 1.5 2008-06-24 18:13:57 faschinj Exp $
        module cloud_sed_mod

        implicit none

        public :: cloud_drop_sed

        private ! Default Scope

        contains
!----------------------------------------------------------------------
        subroutine cloud_drop_sed( rcm, Ncm, rhom, rhot, exner,
     .                             rtm_mc, thlm_mc )

!       Description:
!       Account for cloud droplet sedimentation in cases like DYCOMS II RF 02

!       References:
!       None
!----------------------------------------------------------------------

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm, ddzm ! Procedure(s)

        use constants, only: rho_lw, pi, Cp, Lv ! Variable(s)

#ifdef STATS
        use stats_type, only: stat_update_var ! Procedure(s)

        use stats_variables, only: ised_rcm, iFcsed, zt, zm, lstats_samp ! Variable(s)
#endif

        implicit none

        ! Constant parameters
        ! sigma_g is the Geometric standard deviation of cloud droplets falling
        ! in a Stokes regime.
        real, parameter :: sigma_g = 1.5  ! The DYCOMS2 RF02 intercomparison
                                          ! specified value is 1.5.
!        real, parameter :: sigma_g = 1.2  ! Bjorn Stevens suggests a value of
!                                          ! 1.2 is more appropriate for RF02.

        ! External
        intrinsic :: exp, log
        
        ! Input Variables
        real, intent(in), dimension(gr%nnzp) ::
     .  rcm,     ! Liquid water mixing ratio.   [kg/kg]
     .  rhom,    ! Density on moment. grid      [kg/m^3]
     .  rhot,    ! Density on thermo. grid      [kg/m^3]
     .  exner,   ! Exner function               [-]
     .  Ncm      ! Cloud droplet number conc.   [num/kg]

        ! Input/Output Variables
        real, intent(inout), dimension(gr%nnzp) :: 
     .  rtm_mc, ! r_t change due to microphysics     [kg/kg)/s] 
     .  thlm_mc ! thlm change due to microphysics    [K/s] 

        ! Local Variables
        ! Addition for DYCOMS_2
        real, dimension(gr%nnzp) :: 
     .  Fcsed,    ! Cloud water sedimentation flux       [kg/(m^2 s)]
     .  sed_rcm   ! d(rcm)/dt due to cloud sedimentation [kg/(m^2 s)]

        integer :: k


        !=====================================================================
        ! NOTE:  ADDITION OF RAIN EFFECTS AND CLOUD WATER SEDIMENTATION ON
        !        RTM AND THLM.
        !=====================================================================
        !
        ! Equations:  rtm = rvm + rcm
        !            thlm = thm - ( Lv / (Cp*exner) ) * rcm
        !
        ! When water condenses, latent heat is given off and theta (thm)
        ! increases by a factor of ( Lv / (Cp*exner) ) * rcm(condensed).
        ! The opposite effect occurs with evaporation.
        !
        !=====================================================================
        !||     Effect     |  rvm  |  rcm  |  rtm  |  rrm  |  thm  |  thlm  ||
        !|===================================================================|
        !|| Sedimentation  |       |       |       |       |       |        ||
        !|| Effects of     | stays | incr. | incr. | stays | stays | decr.  ||
        !|| Cloud Water.   | same  |       |       | same  | same  |        ||
        !|| sed_rcm > 0    |       |       |       |       |       |        ||
        !|===================================================================|
        !|| Evaporation    |       |       |       |       |       |        ||
        !|| of rain to     | incr. | stays | incr. | decr. | decr. | decr.  ||
        !|| water vapor.   |       | same  |       |       |       |        ||
        !|| cond_rrm < 0   |       |       |       |       |       |        ||
        !|-------------------------------------------------------------------|
        !|| Autoconversion |       |       |       |       |       |        ||
        !|| of cloud water | stays | decr. | decr. | incr. | stays | incr.  ||
        !|| to rain water. | same  |       |       |       | same  |        ||
        !|| auto_rrm > 0   |       |       |       |       |       |        ||
        !|-------------------------------------------------------------------|
        !|| Accretion of   |       |       |       |       |       |        ||
        !|| cloud water by | stays | decr. | decr. | incr. | stays | incr.  ||
        !|| rain water.    | same  |       |       |       | same  |        ||
        !|| accr_rrm > 0   |       |       |       |       |       |        ||
        !=====================================================================
        !
        ! Note: In HOC, cond_rrm will always be either negative or zero.
        !
        ! Overall effects of rain and cloud water sedimentation:
        !
        ! (drtm/dt)t  = (drtm/dt)0 
        !                       + sed_rcm - cond_rrm - auto_rrm - accr_rrm
        ! (dthlm/dt)t = (dthlm/dt)0  -  ( Lv / (Cp*exner) ) 
        !                       * ( sed_rcm - cond_rrm - auto_rrm - accr_rrm )
        !
        ! Note by Brian Griffin.
        !=====================================================================

        ! Code addition by Brian for cloud water sedimentation.
        !
        ! Sedimentation flux of cloud droplets should be treated by assuming
        ! a log-normal size distribution of droplets falling in a Stoke's
        ! regime, in which the sedimentation flux (Fcsed) is given by:
        !
        ! Sedimentation Flux = constant
        !                     *[(3/(4*pi*rho_lw*NcV))^(2/3)]
        !                     *[(rho*rc)^(5/3)]
        !                     *EXP[5*((LOG(sigma_g))^2)]
        ! constant = 1.19 x 10^8 (m^-1 s^-1)
        ! sigma_g:  geometric standard deviation
        !
        ! When written for a mass-dependent cloud droplet concentration, Nc:
        !
        ! Sedimentation Flux = constant
        !                     *[(3/(4*pi*rho_lw*Nc*rho))^(2/3)]
        !                     *[(rho*rc)^(5/3)]
        !                     *EXP[5*((LOG(sigma_g))^2)]
        !
        ! According to the above equation, sedimentation flux
        ! is defined positive downwards.  Therefore, 
        !
        ! (drc/dt)Fcsed = (1.0/rho) * d(Fcsed)/dz


        ! Define cloud water sedimentation flux on momentum levels.

        DO k = 2, gr%nnzp-1, 1

          IF ( zt2zm(rcm,k) > 0.0 .AND. zt2zm( Ncm, k ) > 0.0 ) THEN
            Fcsed(k) = 1.19E8
     .                  * (   (  3.0 / ( 4.0*pi*rho_lw
     .                                   *zt2zm( Ncm, k )*rhom(k) )  )
     .                     **(2.0/3.0)   )
     .                  * ( ( rhom(k)*zt2zm( rcm, k ) )**(5.0/3.0) )
     .                  * EXP( 5.0*( (LOG( sigma_g ))**(2) ) )
          ELSE
            Fcsed(k) = 0.0
          END IF

        END DO ! k=2..gr%nnzp-1

        ! Boundary conditions.
        Fcsed(1)       = 0.0
        Fcsed(gr%nnzp) = 0.0

        ! Find drc/dt due to cloud water sedimentation flux.
        ! This value is defined on thermodynamic levels.
        ! Fcsed units:  kg (liquid) / [ m^2 * s ]
        ! Multiply by Lv for units of W / m^2.
        ! sed_rcm units:  [ kg (liquid) / kg (air) ] / s
        sed_rcm = (1.0/rhot) * ddzm( Fcsed )

#ifdef STATS
        if ( lstats_samp ) then
          call stat_update_var( ised_rcm, sed_rcm, zt )

          call stat_update_var( iFcsed, Fcsed, zm ) 
        end if
#endif


        ! + thlm/rtm_microphysics -- cloud water sedimentation.
        ! Code addition by Brian for cloud water sedimentation.

        rtm_mc  = rtm_mc + sed_rcm
        thlm_mc = thlm_mc - ( Lv / (Cp*exner) ) * sed_rcm

        return
        end subroutine cloud_drop_sed

        end module cloud_sed_mod
