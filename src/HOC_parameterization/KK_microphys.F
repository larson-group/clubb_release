!----------------------------------------------------------------------
! $Id: KK_microphys.F,v 1.17 2008-07-21 21:47:31 faschinj Exp $

      module rain_equations

!     Description:
!     Functions and subroutines for the Khairoutdinov & Kogan
!     rain parameterization.

!     References:
!     None
!----------------------------------------------------------------------

      implicit none

      public :: kk_microphys
      private :: mean_volume_radius, cond_evap_rrm, cond_evap_Nrm
      private :: autoconv_rrm, autoconv_Nrm, accretion_rrm
      PRIVATE :: G_T_p
!      PRIVATE :: PDF_TRIVAR_2G_LN_LN, PDF_BIVAR_2G_LN, PDF_BIVAR_LN_LN
      PRIVATE :: Dv_fnc

      PRIVATE ! Set default scope
    
      ! Statistical rain parameters.
      ! Old parameterization.
!      REAL, PARAMETER:: rrp2_rrm2 = 1.2
!      REAL, PARAMETER:: Nrp2_Nrm2 = 1.0
!      REAL, PARAMETER:: Ncp2_Ncm2 = 0.07
!      REAL, PARAMETER:: corr_rrNr_LL =  0.85
!      REAL, PARAMETER:: corr_srr_NL =  0.10
!      REAL, PARAMETER:: corr_sNr_NL = -0.15
!      REAL, PARAMETER:: corr_sNc_NL =  0.45
      ! Parameters for in-cloud (from SAM RF02 DO).
      REAL, PARAMETER:: rrp2_rrm2_cloud = 0.766
      REAL, PARAMETER:: Nrp2_Nrm2_cloud = 0.429
      REAL, PARAMETER:: Ncp2_Ncm2_cloud = 0.003
      REAL, PARAMETER:: corr_rrNr_LL_cloud = 0.786
      REAL, PARAMETER:: corr_srr_NL_cloud = 0.242
      REAL, PARAMETER:: corr_sNr_NL_cloud = 0.285
      REAL, PARAMETER:: corr_sNc_NL_cloud = 0.433
      ! Parameters for below-cloud (from SAM RF02 DO).
      REAL, PARAMETER:: rrp2_rrm2_below = 8.97
      REAL, PARAMETER:: Nrp2_Nrm2_below = 12.03
      REAL, PARAMETER:: Ncp2_Ncm2_below = 0.00 ! Not applicable below cloud.
      REAL, PARAMETER:: corr_rrNr_LL_below = 0.886
      REAL, PARAMETER:: corr_srr_NL_below = 0.056
      REAL, PARAMETER:: corr_sNr_NL_below = 0.015
      REAL, PARAMETER:: corr_sNc_NL_below = 0.00 ! Not applicable below cloud.
      REAL, PARAMETER:: C_evap = 0.86    ! Khairoutdinov and Kogan (2000) 
                                         ! ratio of drizzle drop mean 
                                         ! geometric radius to drizzle drop
                                         ! mean volume radius. Khairoutdinov
                                         ! and Kogan (2000); p. 233.
      !REAL, PARAMETER:: C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
      !REAL, PARAMETER:: C_evap = 0.55    ! KK 2000, Marshall-Palmer (1948) value.
      ! Vince Larson set r_0=28mum to agree with COAMPS-LES formula. 15 April 2005
      REAL, PARAMETER:: r_0 = 25.0e-6   ! Assumed radius of all new drops; m.
                                        ! Value specified in KK (2000); p. 235.
      !REAL, PARAMETER:: r_0 = 28.0e-6   ! Assumed radius of all new drops; m.
      !                                  ! Value that COAMPS LES has in it.
      !REAL, PARAMETER:: r_0 = 30.0e-6   ! Assumed radius of all new drops; m.
      !                                  ! Khairoutdinov said it was okay!
      ! End Vince Larson's change.

!!!!!!! COMMENTS FOR RICO MARINE CUMULUS CASE
!!!!!!! In order to correctly run the RICO case, the parameters above
!!!!!!! need to be changed.  You should comment these in and 
!!!!!!! comment the above copies out.  These are the constants needed to 
!!!!!!! run the RICO case, as I did in my case rico0206, submitted:
!      REAL, PARAMETER:: rrp2_rrm2_cloud = 30.
!      REAL, PARAMETER:: Nrp2_Nrm2_cloud = 10.
!      REAL, PARAMETER:: Ncp2_Ncm2_cloud = 60.
!      REAL, PARAMETER:: corr_rrNr_LL_cloud = 0.8
!      REAL, PARAMETER:: corr_srr_NL_cloud = 0.25
!      REAL, PARAMETER:: corr_sNr_NL_cloud = 0.3
!      REAL, PARAMETER:: corr_sNc_NL_cloud = 0.24

!      REAL, PARAMETER:: rrp2_rrm2_below = 20.
!      REAL, PARAMETER:: Nrp2_Nrm2_below = 2.
!      REAL, PARAMETER:: Ncp2_Ncm2_below = 0.
!      REAL, PARAMETER:: corr_rrNr_LL_below = 0.78
!      REAL, PARAMETER:: corr_srr_NL_below = 0.
!      REAL, PARAMETER:: corr_sNr_NL_below = 0.
!      REAL, PARAMETER:: corr_sNc_NL_below = 0.

!      REAL, PARAMETER:: C_evap = 0.5
!      REAL, PARAMETER:: r_0 = 100.0e-6
!!!!!!! END RICO COMMENTS, Michael Falk, 27-28 November 2007

      CONTAINS

        subroutine kk_microphys
     .             ( thlm, T_in_K, p, exner, rhot, 
     .               thl1, thl2, a, rc1, rc2, s1,
     .               s2, ss1, ss2, rtm, rcm, Ncm, rrm, Nrm, 
     .               lsample, AKm_est, AKm,
     .               rrm_mc_tndcy, Nrm_mc_tndcy, 
     .               hm_rt_tndcy, hm_thl_tndcy,
     .               Vrr, VNr )
!       Description:
!       This uses the code from Brian Griffin's old subroutine rain
!       and computes microphysical tendencies and fall speeds
!       per se rather than adding in the effects of diffusion, advection, 
!       and sedimentation.  This is for the purpose of generalizing the
!       code to use different microphysical schemes with the same 
!       diffusion, advection, and sedimentation code.

!       References:
!       ``A New Cloud Physics Parameterization in a Large-Eddy Simulation
!         Model of Marine Stratocumulus''  Khairoutdinov and Kogan. (2000)
!       Monthly Weather Review, Volume 128, Issue 1 pp. 229–-243
!----------------------------------------------------------------------

        use model_flags, only:
     .      LH_on ! Variable(s)

        use grid_class, only:
     .      gr, ! Variable(s)
     .      zt2zm ! Procedure(s)

        use constants, only:
     .      rc_tol, ! Variable(s)
     .      ep,
     .      Cp,
     .      Lv,
     .      kappa,
     .      p0,
     .      Rd,
     .      Rv

        use saturation, only: 
     .      sat_mixrat_liq ! Procedure(s)
        
#ifdef STATS
        use stats_type, only:
     .      stat_update_var, stat_update_var_pt ! Procedure(s)

        use stats_variables, only:
     .      zt,  ! Variable(s)
     .      imean_vol_rad_rain,
     .      irrm_cond,
     .      irrm_auto,
     .      irrm_accr,
     .      iNrm_cond,
     .      iNrm_auto,
     .      lstats_samp
#endif

        implicit none


        ! Input
        real, intent(in), dimension(gr%nnzp) ::
     .  thlm,      ! Theta-l                            [K]
     .  T_in_K,    ! Temperature                        [K]
     .  p,         ! Pressure                           [Pa]
     .  exner,     ! Exner function                     [-]
     .  rhot,      ! Density on thermo. grid            [kg/m^3]
     .  thl1, thl2,! PDF parameters thl1 &thl2          [K]
     .  a,         ! PDF parameter a                    [-]
     .  s1, s2,    ! PDF parameters s1 & s2             [kg/kg]
     .  ss1, ss2,  ! PDF parameters ss1 & ss2           [kg/kg]
     .  rc1, rc2,  ! PDF parameters rc1 & rc2           [kg/kg]
     .  rtm,       ! Total water mixing ratio           [kg/kg]
     .  rcm,       ! Cloud water mixing ratio           [kg/kg]
     .  Ncm,       ! Cloud droplet number conc.         [number/kg]
     .  rrm,       ! Rain water mixing ratio            [kg/kg]
     .  Nrm        ! Rain drop number conc.             [number/kg]

        ! Latin hypercube variables - Vince Larson 22 May 2005
        real, intent(in), dimension(gr%nnzp) :: 
     .  Akm,    ! Latin hypercube estimate of Kessler autoconversion   
     .  AKm_est ! Analytic calculation of Kessler autoconversion

        logical, intent(in) ::
     .  lsample ! Whether to sample stats for this call

        ! Output
        real, intent(out), dimension(gr%nnzp) :: 
     .  Vrr,         ! Mean sedimentation velocity of rrm       [m/s]    
     .  VNr,         ! Mean sedimentation velocity of Nrm       [m/s]
     .  rrm_mc_tndcy,! Rain water microphysical tendency        [(kg/kg)/s]
     .  Nrm_mc_tndcy,! Rain drop number conc. micro. tend.      [(num/kg)/s]
     .  hm_rt_tndcy, ! Contributions to total water from micro. [(kg/kg)/s]
     .  hm_thl_tndcy ! Contributions to theta_l from micro.     [K/s]

        ! Local variables
        real, dimension(gr%nnzp) :: 
     .  rrm_cond,    ! Change in rrm due to condensation        [(kg/kg)/s]
     .  rrm_auto,    ! Change in rrm due to autoconversion      [(kg/kg)/s]
     .  rrm_accr,    ! Change in rrm due to accretion           [(kg/kg)/s]
     .  Nrm_cond,    ! Change in Nrm due to condensation        [(num/kg)/s]
     .  Nrm_auto,    ! Change in Nrm due to autoconversion      [(num/kg)/s]
     .  mean_vol_rad,! Mean volume radius                       [m]
!     .  rvm,         ! Water vapor mixing ratio                 [kg/kg]
!     .  esat,          ! Saturation vapor pressure                [Pa]
     .  Supsat,      ! Supersaturation                          [-]
     .  rsat          ! Saturation mixing ratio                  [kg/kg]

        real, dimension(gr%nnzp) :: 
     .  rrp2_rrm2,   ! rrp2/rrm^2               []
     .  Nrp2_Nrm2,   ! Nrp2/Nrm^2               []
     .  Ncp2_Ncm2,   ! Ncp2/Ncm^2               []
     .  corr_rrNr_LL,! Correlation of rr and Nr []
     .  corr_srr_NL, ! Correlation of s and rr  []
     .  corr_sNr_NL, ! Correlation of s and Nr  []
     .  corr_sNc_NL  ! Correlation of s and Nc  []

        real :: 
!     .  e,      ! Vapor pressure                [Pa]
     .  Beta_T  ! Beta_T                        [kg/kg]

        ! Array indices
        integer :: k

! DESCRIPTION:
!
! This subroutine computes microphysical contributions of 
! the new rain water mixing ratio, rr.
! It computes it according to the equation:
!
! drr/dt(explicit) = (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! It also computes Vrr, which is applied to the LHS of the tridiagonal
! system in subroutine timestep_microphys like so:
!
! drr/dt(implicit) = -w*(drr/dz) + Vrr*(drr/dz) + (K+nu_r)*(d2rr/dz2)
!
! where Vrr is the mean sedimentation velocity of the rain drops.  It
! is determined by Vrr = 0.012*rvr - 0.2, with Vrr in m/s and rvr, the
! rain drop mean volume radius, in um.  In turn the rain drop mean
! volume radius is given by the following equation:
!
! rvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^(1/3) * NrV^(-1/3);
!
! where NrV is the rain drop concentration (num/m^3).  rvr is given in
! meters.
!
! It is important to note that a one-sided upwind advection scheme is used
! for sedimentation velocity, with Vrr being from the level above the level
! which is being dealt with.  This will increase the diffusivity of the
! model.
!
! For the other terms:
!
! (drr/dt)cond = 
!      3*Cevap*G(T,p)*[(4*PI*rho_lw)/(3*rho)]^(2/3) * rr^(1/3) * NrV^(2/3) * S
!
! G(T,p) = 1/(Fk + Fd); where
!
! Fk = [Lv/(Rv*T) - 1]*(Lv*rho_lw)/(Ka*T)
! Fd = (rho_lw*Rv*T)/(Dv*esat(T))
!
! Ka is the Coefficient of Thermal Conductivity of Air, and Dv is the
! Coefficient of Diffusion of Water Vapor in Air.
!
! (drr/dt)auto = 1350 * rc^2.47 * NcV^-1.79;
!
! with cloud drop concentration in (num/cm^3) for this particular calculation.
!
! (drr/dt)accr = 67*(rc*rr)^1.15
!
! The problem is solved by moving around some terms to result in the following
! equation:
!
! drr/dt + w*(drr/dz) - Vrr*(drr/dz) - nu_r*(d2rr/dz2)
!                      = (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! The equation is then broken down as following:
! (n is the current timestep; n+1 is the next timestep; 
! k is the vertical level).
!
! rr(n+1,k)/dt
! + w(n,k)*(rr(n+1,k+1)-rr(n+1,k-1))/(del_zm(n+1,k)-del_zm(n+1,k-1))
! - Vr(n,k+1)*(rr(n+1,k+1)-rr(n+1,k))/(del_zm(n+z,k))
! + nu_r*(1/del_zt(n+1,k)) 
!      * [  (rr(n+1,k+1)-rr(n+1,k))/del_zm(n+1,k) 
!         - (rr(n+1,k)-rr(n+1,k-1))/del_zm(n+1,k-1)   ]
! = rr(n,k)/dt + (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! The equation is then broken down into a matrix:
!
! LHS*rr = RHS
!
! A tridiagonal matrix is used, with 
! subdiagonal aa that deals with rr(n+1,k-1) terms,
! main diagonal bb that deals with rr(n+1,k) terms,
! superdiagonal cc that deals with rr(n+1,k+1) terms, and
! dd, which is the right-hand side matrix with rr(n,k) terms.
!
! IMPORTANT NOTES
!
! The equations for mean volume radius (and therefore sedimentation
! velocity), condensation/evaporation, autoconversion, and accretion
! stated above are all equations from Khairoutdinov and Kogan (2000).
! These are called the "local" equations, because they show the local
! value at a grid point with high resolution.  These equations would
! be used in a LES model.  HOC is a one-dimensional model that uses
! a Probability Density Function (PDF) to determine variance of values
! in the horizontal directions.  When put into any 3-dimensional
! model, HOC can use the PDF to show the subgrid variability of many
! values.  In order to determine rainfall due to subgrid variability,
! one needs to put Khairoutdinov and Kogan equations into a PDF
! form.  Here is a brief description.
!
! a) Mean Volume Radius (and therefore sedimentation velocity) uses a
!    bivariate lognormal distribution because the factors that make it
!    up (rr and Nr) are both distributed lognormally.
!
! b) Condensation/evaporation uses a single normal-lognormal-lognormal
!    PDF due to the fact that supersaturation, S, (one value that makes
!    it up) follows a truncated double Gaussian, rr (another factor that
!    makes it up) follows a single lognormal distribution, and Nr (the
!    last factor that makes it up) also follows a single lognormal
!    distribution.
!
! c) Autoconversion uses a single normal-lognormal PDF due to the fact
!    that rc (one factor that makes it up) follows a truncated double
!    Gaussian and Nc (the other factor that makes it up) follows a
!    single lognormal distribution.
!
! d) Accretion uses a single normal-lognormal PDF due to the fact
!    that rc (one factor that makes it up) follows a truncated double
!    Gaussian and rr (the other factor that makes it up) follows a
!    single lognormal distribution.
!
! The above PDF-ed Khairoutdinov and Kogan equations are used for the
! "non-local" formula, which is the one used in HOC.
!
! This subroutine also solves for rain drop number concentration (num/m^3)
! using the exact same formula as noted in this description.  The left-hand
! side of the equation is exactly the same as for rr, except for the fact
! that VNr is used instead of Vrr.  VNr also uses the mean volume radius
! for it's calculation.  The equation is VNr = 0.007*rvr - 0.1, with VNr
! in m/s and rvr in micrometers.  The right-hand side of the Nr equation
! only contains terms for time tendency, condensation/evaporation, and
! autoconversion.
!
! (dNrV/dt)cond = (NrV/rr) * (drr/dt)cond
!
! (dNrV/dt)auto = (drr/dt)auto / ( [(4*PI*rho_lw)/(3*rho)] * r_o^3 )
!
! where r_o is the assumed radius of all newly-formed rain droplets (m).
! It is given a value of 25 um.
!
! This information is plugged into the tridiagonal matrix in order to
! solve for NrV at the next timestep.
!
! Description and notes written by Brian Griffin.
!
! Changes for kk_microphys -dschanen 14 Feb 2007

        ! Set up the values of the statistical correlations and variances.
        ! Since we currently do not have enough variables to compute the 
        ! correlations and variances directly, we have obtained these values 
        ! by analyzing LES runs of certain cases.  We have divided those 
        ! results into an inside-cloud average and an outside-cloud (or 
        ! below-cloud) average.  This coding leaves the software architecture 
        ! in place in case we ever have the variables in place to compute 
        ! these values directly.  It also allows us to use separate 
        ! inside-cloud and outside-cloud parameter values.  
        ! Brian Griffin; February 3, 2007.
        !
        ! Set the value of the parameters based on whether the altitude
        ! is above or below cloud  base.
        ! Determine whether there is cloud at any given vertical level.  In 
        ! order for a vertical level to have cloud, the amount of cloud water
        ! (rcm) must be greater than or equal to the tolerance level (rc_tol).
        ! If there is cloud at a given vertical level, then the ###_cloud value
        ! is used.  Otherwise, the ###_below value is used.
        DO k = 1, gr%nnzp, 1
           IF ( rcm(k) >= rc_tol ) THEN
              rrp2_rrm2(k)    = rrp2_rrm2_cloud
              Nrp2_Nrm2(k)    = Nrp2_Nrm2_cloud
              Ncp2_Ncm2(k)    = Ncp2_Ncm2_cloud
              corr_rrNr_LL(k) = corr_rrNr_LL_cloud
              corr_srr_NL(k)  = corr_srr_NL_cloud
              corr_sNr_NL(k)  = corr_sNr_NL_cloud
              corr_sNc_NL(k)  = corr_sNc_NL_cloud
           ELSE
              rrp2_rrm2(k)    = rrp2_rrm2_below
              Nrp2_Nrm2(k)    = Nrp2_Nrm2_below
              Ncp2_Ncm2(k)    = Ncp2_Ncm2_below
              corr_rrNr_LL(k) = corr_rrNr_LL_below
              corr_srr_NL(k)  = corr_srr_NL_below
              corr_sNr_NL(k)  = corr_sNr_NL_below
              corr_sNc_NL(k)  = corr_sNc_NL_below
           ENDIF
        ENDDO

        ! Find the drop mean volume radius.  It is calculated using
        ! the rain water mixing ratio, the rain droplet concentration,
        ! and the air density.  These values are taken from the previous
        ! timestep.  It is located on thermodynamic levels.
        do k = 1, gr%nnzp, 1

           mean_vol_rad(k) 
     .     = mean_volume_radius
     .              ( rrm(k), Nrm(k), rhot(k), rrp2_rrm2(k), 
     .                           Nrp2_Nrm2(k), corr_rrNr_LL(k) )

        end do

#ifdef STATS
        ! Save mean volume radius for stats purposes
        ! Note: added lsample for latin hypercube sampling -dschanen
        if ( lsample .and. lstats_samp ) then
          call stat_update_var(imean_vol_rad_rain, mean_vol_rad, zt )
        end if

#endif /*STATS*/

        ! The sedimentation velocity is found from the drop mean volume radius.
        ! It is located on the momentum levels.  This is due to the fact that
        ! momentum levels are where other vertical velocities are stored, such
        ! as the vertical component of wind velocity (w).
        !
        ! Khairoutdinov and Kogan Sedimentation Velocity Calculation.
        ! sedimentation velocity of rain drops (in m/s).
        ! mean volume radius converted to um.
        ! Note 1: positive sedimentation velocity means downwards.
        ! Note 2: Vrr(m/s) = 0.012*rvr(um) - 0.2
        !
        ! Note 3:
        ! Changed sedimentation to be negative and consistent with
        ! the COAMPS microphysics.
        ! -dschanen 5 Dec 2006
        !
        ! If rvr(um) is too small, the equation will result in an
        ! upwards (positive) sedimentation velocity.  The limiter
        ! is put in there for that reason.
        !
        ! Note 4:
        ! The original Khairoutdinov and Kogan prognostic equations for rain 
        ! water mixing ratio (eq. 8) and rain drop concentration (eq. 9)
        ! defined sedimentation velocity as positive downwards and read:
        !
        ! d(rr)/dt = -w (drr/dz) + V_rr (drr/dz) + ...;
        !
        ! However, to be consistent with COAMPS microphysical schemes, we 
        ! have defined sedimentation velocity to be positive upwards (or 
        ! negative downwards), which alters the prognostic equation to read:
        !
        ! d(rr)/dt = -w (drr/dz) - V_rr (drr/dz) + ...;
        !
        ! This enters into computation in microphysics.F.  Brian Griffin.

        do k = 1, gr%nnzp-1, 1

           ! rrm sedimentation velocity.
           Vrr(k) = 0.012 
     .             * ( 1.0e6 * zt2zm(mean_vol_rad,k)  )
     .             - 0.2
           ! Negative meaning a downward velocity now -dschanen 5 Dec 2006
           Vrr(k) = -max( Vrr(k), 0.0 )

           ! Nrm sedimentation velocity.
           VNr(k) = 0.007 
     .             * ( 1.0e6 * zt2zm(mean_vol_rad,k)  )
     .             - 0.1

           ! Negative meaning a downward velocity now -dschanen 5 Dec 2006
           VNr(k) = -max( VNr(k), 0.0 )

        end do ! 1..gr%nnzp

        ! The flux of rain water through the model top is 0.
        ! Vrr and VNr are set to 0 at the highest model level.
        Vrr(gr%nnzp) = 0.0
        VNr(gr%nnzp) = 0.0

        ! Find values for other variables.
        do k = 2, gr%nnzp, 1

           ! Find the important conditions

!           rvm(k) = rtm(k) - rcm(k)   ! rvm is water vapor mixing ratio.

           ! Find current vapor pressure.
!           e = (p(k)*rvm(k))/(ep + rvm(k))

           ! Saturation mixing ratio
           rsat(k) = sat_mixrat_liq( p(k), T_in_K(k) )

           ! Find saturation vapor pressure.
!           esat(k) = (p(k)*rsat(k))/(ep + rsat(k))

        end do ! 2..gr%nnzp

        ! Set the boundary conditions
!        rvm(1)          = 0.0
        rsat(1)          = 0.0
!        esat(1)           = 0.0
        Supsat(1)       = 0.0
        Supsat(gr%nnzp) = 0.0

        do k = 2, gr%nnzp-1, 1

           ! Compute supersaturation via s1, s2. 
           !     Larson et al 2002, JAS, Vol 59, p 3534.
           ! This allows a more direct comparison of local, nonlocal formulas.
           Beta_T = (Rd/Rv) * ( Lv/(Rd*T_in_K(k)) ) 
     .              * ( Lv/(Cp*T_in_K(k)) )

           Supsat(k) = ( a(k)*s1(k) + (1-a(k))*s2(k) )
     .                 *( ( 1.0 + Beta_T*rsat(k) ) / rsat(k) )

           ! Now find the elements that make up the right-hand side of the
           ! equation, dd, for rain water mixing ratio, rrm.

           rrm_cond(k) 
     .     = cond_evap_rrm
     .       ( rrm(k), Nrm(k),
     .         s1(k), ss1(k), s2(k), ss2(k),
     .         thl1(k), thl2(k), rc1(k), rc2(k), a(k),
     .         p(k), exner(k), rhot(k), T_in_K(k), Supsat(k), 
     .         rrp2_rrm2(k), Nrp2_Nrm2(k), corr_srr_NL(k),
     .         corr_sNr_NL(k), corr_rrNr_LL(k) )


        ! Vince Larson added option to call LH sampled Kessler autoconversion.
        ! 22 May 2005
!           rrm_auto(k) = autoconv_rrm( rcm(k), Ncm(k), rhot(k),
!     .                        a, s1, s2, ss1, ss2, rc1, rc2 )
           if ( LH_on ) then

!              rrm_auto(k) = AKm_est(k)
              rrm_auto(k) = AKm(k)

           else

              rrm_auto(k) 
     .        = autoconv_rrm( rcm(k), Ncm(k), s1(k), ss1(k), 
     .                        s2(k), ss2(k), a(k), rhot(k),
     .                        Ncp2_Ncm2(k), corr_sNc_NL(k) )

           end if ! LH_on
        ! End Vince Larson's addition

           rrm_accr(k) 
     .     = accretion_rrm( rcm(k), rrm(k), s1(k), ss1(k),
     .                      s2(k), ss2(k), a(k), 
     .                      rrp2_rrm2(k), corr_srr_NL(k) )

           ! Now find the elements that make up the right-hand side of the
           ! equation, rr, for rain drop number concentration, Nrm.

           Nrm_cond(k) = cond_evap_Nrm( rrm_cond(k), Nrm(k), rrm(k) )

           Nrm_auto(k) = autoconv_Nrm( rrm_auto(k) )

#ifdef STATS
           if ( lsample .and. lstats_samp ) then
             ! Explicit contributions to rrm.
             call stat_update_var_pt( irrm_cond, k, rrm_cond(k), zt )

             call stat_update_var_pt( irrm_auto, k, rrm_auto(k), zt )

             call stat_update_var_pt( irrm_accr, k, rrm_accr(k), zt )

             ! Explicit contributions to Nrm.
             call stat_update_var_pt( iNrm_cond, k, Nrm_cond(k), zt )

             call stat_update_var_pt( iNrm_auto, k, Nrm_auto(k), zt )

           end if ! lstats_samp and lsample
#endif /*STATS*/

          rrm_mc_tndcy(k) = rrm_cond(k) + rrm_auto(k) + rrm_accr(k)

          Nrm_mc_tndcy(k) = Nrm_cond(k) + Nrm_auto(k)

          ! Explicit contributions to thlm and rtm from the microphysics
          hm_rt_tndcy(k)  = -rrm_mc_tndcy(k)
          hm_thl_tndcy(k) = ( Lv / ( Cp*exner(k) ) ) * rrm_mc_tndcy(k)
        end do ! k=2..gr%nnzp-1


        ! Boundary conditions

        ! Explicit contributions to rrm and Nrm from microphysics
        ! auto_rrm, cond_rrm, and accr_rrm are not set at nz=1 or nnzp
        ! There is no condensation/evaporation, autoconversion, or accretion
        ! at level 1, which is below the ground surface.  Brian.
        rrm_mc_tndcy(1) = 0.0
        Nrm_mc_tndcy(1) = 0.0

        rrm_mc_tndcy(gr%nnzp) = 0.0
        Nrm_mc_tndcy(gr%nnzp) = 0.0

        ! Contributions to theta_l and rt.  See further comments 
        ! about this in the subroutine mixing_solve.
        hm_rt_tndcy(1)  = 0.0
        hm_thl_tndcy(1) = 0.0

        hm_rt_tndcy(gr%nnzp)  = 0.0
        hm_thl_tndcy(gr%nnzp) = 0.0

        return
        end subroutine kk_microphys

!===============================================================================
        !
        ! FUNCTION mean_volume_radius
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! This function calculates the mean volume radius of a rain drop from
        ! Khairoutdinov and Kogan (2000) equation 3. That equation is:
        !
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^(1/3) * NrV^(-1/3)
        !
        ! This turns into:
        !
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^alpha * NrV^beta
        !
        ! where alpha = 1/3 and beta = -1/3
        !
        ! When applied to a Probability Density Function (PDF), it becomes:
        !
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * Y1^alpha * Y2^beta
        !
        ! Y1 is used for rr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        ! Y2 is used for NrV.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Notes:  For our model, Nrm is in concentration per mass instead of 
        !         concentration per volume.  NrV is divided by air density 
        !         (rho) in order to achieve that value.  That also explains why
        !         rho is removed from the denominator of the factor 
        !         [(4*PI*rho_lw)/(3*rho)]^(2/3).  The factor becomes 
        !         [(4*PI*rho_lw)/(3)]^(2/3).
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION mean_volume_radius( rrm, Nrm, rhot, rrp2_rrm2, 
     .                               Nrp2_Nrm2, corr_rrNr_LL )

        USE constants, only:
     .      Nr_tol, ! Variable(s)
     .      rr_tol,
     .      rho_lw,
     .      pi
        
        USE model_flags, only:
     .      local_kk ! Variable(s)

        implicit none

        ! Input variables.
        REAL, INTENT(IN):: rrm          ! Grid-box average rrm     [kg kg^-1]
!        REAL, INTENT(IN):: rrp2         ! Grid-box rr variance     [kg^2 kg^-2]
        REAL, INTENT(IN):: Nrm          ! Grid-box average Nrm     [kg^-1]
!        REAL, INTENT(IN):: Nrp2         ! Grid-box Nr variance     [kg^-2]
        REAL, INTENT(IN):: rhot         ! Grid-box average density (t-level)
                                        !                          [kg m^-3]
        REAL, INTENT(IN):: rrp2_rrm2    ! rrp2/rrm^2               []
        REAL, INTENT(IN):: Nrp2_Nrm2    ! Nrp2/Nrm^2               []
        REAL, INTENT(IN):: corr_rrNr_LL ! Correlation of rr and Nr []

        ! Output variables.
        REAL:: mean_volume_radius       ! [m]

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of rr
        REAL:: beta_exp   ! Exponent of Nr

        ! Original terms.
        REAL:: mu_rr       ! Grid-box average of rr                [kg kg^-1]
        REAL:: sigma_rr    ! Grid-box standard deviation of rr     [kg kg^-1]
        REAL:: mu_Nr       ! Grid-box average of Nr                [kg^-1]
        REAL:: sigma_Nr    ! Grid-box standard deviation of Nr     [kg^-1]
        REAL:: corr_rrNr   ! Correlation of rr and Nr              []

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              mean_volume_radius = 
     .           (  ( (4.0/3.0)*pi*rho_lw )**(-1.0/3.0)  )
     .          * rrm**(1.0/3.0) * Nrm**(-1.0/3.0)

           ELSE

              ! If either rrm or Nrm are 0.
              mean_volume_radius = 0.0

           ENDIF

        ELSEIF ( .NOT. local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              ! Exponents on rr and Nr, respectively.
              alpha_exp = (1.0/3.0)
              beta_exp  = -(1.0/3.0)

              ! rr is distributed Lognormally.
              mu_rr = rrm
!              sigma_rr = SQRT(rrp2)
              sigma_rr = rrm * SQRT(rrp2_rrm2)

              ! Nr is distributed Lognormally.
              mu_Nr = Nrm
!              sigma_Nr = SQRT(Nrp2)
              sigma_Nr = Nrm * SQRT(Nrp2_Nrm2)

              ! Correlations.
              corr_rrNr = corr_rrNr_LL


              mean_volume_radius =
     .           (  ( (4.0/3.0)*pi*rho_lw )**(-1.0/3.0)  )
     .           * PDF_BIVAR_LN_LN ( mu_rr, mu_Nr, sigma_rr, sigma_Nr,
     .                               corr_rrNr, alpha_exp, beta_exp )

           ELSE

              ! If either rrm or Nrm are 0.
              mean_volume_radius = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION mean_volume_radius

!===============================================================================
        !
        ! FUNCTION cond_evap_rrm
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! This function calculates the evaporation of rain water based on
        ! Khairoutdinov and Kogan (2000) equation 22. That equation is:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * rr^(1/3) * NrV^(2/3) * S
        !
        ! where S = (e/esat) - 1
        !
        ! since S = [ ( 1 + Beta_T * rs ) / rs ] * s
        !
        ! (where s is the extended liquid water specific humidity),
        !
        ! the equation becomes:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * rr^(1/3) * NrV^(2/3) * [ ( 1 + Beta_T * rs ) / rs ] * s
        !
        ! This turns into:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * [ ( 1 + Beta_T * rs ) / rs ]
        !      * s^alpha * rr^beta * NrV^gamma
        !
        ! where alpha = 1.0, beta = 1/3, and gamma = 2/3
        !
        ! When applied to a Probability Density Function (PDF), it becomes:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * [ ( 1 + Beta_T * rs ) / rs ]
        !      * Y1^alpha * Y2^beta * Y3^gamma
        !
        ! Y1 is used for s.  It is distributed as a double Gaussian with a
        ! domain from negative infinity to zero.  It is truncated at zero
        ! because we are only interested in cases of evaporation, not
        ! condensation.
        ! Y2 is used for rr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        ! Y3 is used for NrV.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Notes:  For our model, Nrm is in concentration per mass instead of 
        !         concentration per volume.  NrV is divided by air density 
        !         (rho) in order to achieve that value.  That also explains why
        !         rho is removed from the denominator of the factor 
        !         [(4*PI*rho_lw)/(3*rho)]^(2/3).  The factor becomes 
        !         [(4*PI*rho_lw)/(3)]^(2/3).
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION cond_evap_rrm( rrm, Nrm,
     .                          s1, ss1, s2, ss2,
     .                          thl1, thl2, rc1, rc2, a,
     .                          press, exner, rhot, T_in_K, Supsat, 
     .                          rrp2_rrm2, Nrp2_Nrm2, corr_srr_NL,
     .                          corr_sNr_NL, corr_rrNr_LL  )

        USE constants, only:
     .      Nr_tol, ! Variable(s)
     .      rr_tol,
     .      rho_lw,
     .      kappa,
     .      p0,
     .      pi,
     .      Cp,
     .      Lv,
     .      Rd,
     .      Rv
        
        USE model_flags, only:
     .      local_kk ! Variable(s)

        USE saturation, only: 
     .      sat_mixrat_liq ! Variable(s)

        implicit none

        ! Input variables.
        REAL, INTENT(IN):: rrm          ! Grid-box average rrm     [kg kg^-1]
!        REAL, INTENT(IN):: rrp2         ! Grid-box rr variance     [kg^2 kg^-2]
        REAL, INTENT(IN):: Nrm          ! Grid-box average Nrm     [kg^-1]
!        REAL, INTENT(IN):: Nrp2         ! Grid-box Nr variance     [kg^-2]
        REAL, INTENT(IN):: s1           ! Plume 1 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss1          ! Plume 1 sigma s1 (not sigma^2 s1)
                                        !                          [kg kg^-1]
        REAL, INTENT(IN):: s2           ! Plume 2 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss2          ! Plume 2 sigma s2 (not sigma^2 s2)
                                        !                          [kg kg^-1]
        REAL, INTENT(IN):: thl1         ! Plume 1 average theta-l  [K]
        REAL, INTENT(IN):: thl2         ! Plume 2 average theta-l  [K]
        REAL, INTENT(IN):: rc1          ! Plume 1 average rc       [kg kg^-1]
        REAL, INTENT(IN):: rc2          ! Plume 2 average rc       [kg kg^-1]
        REAL, INTENT(IN):: a            ! Relative weight of each individual
                                        ! Gaussian "plume."        []
        REAL, INTENT(IN):: press        ! Grid-box average pressure [Pa]
        REAL, INTENT(IN):: exner        ! Grid-box average exner function [-]
        REAL, INTENT(IN):: rhot         ! Grid-box average density (t-level)
                                        !                          [kg m^-3]
        REAL, INTENT(IN):: T_in_K         ! Grid-box average Temperature [K]
        REAL, INTENT(IN):: Supsat       ! Grid-box average Supersaturation []
        REAL, INTENT(IN):: rrp2_rrm2    ! rrp2/rrm^2               []
        REAL, INTENT(IN):: Nrp2_Nrm2    ! Nrp2/Nrm^2               []
        REAL, INTENT(IN):: corr_srr_NL  ! Correlation of s and rr  []
        REAL, INTENT(IN):: corr_sNr_NL  ! Correlation of s and Nr  []
        REAL, INTENT(IN):: corr_rrNr_LL ! Correlation of rr and Nr []

        ! Output variables.
        REAL:: cond_evap_rrm  ! [kg kg^-1 s^-1]

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s
        REAL:: beta_exp   ! Exponent of rr
        REAL:: gamma_exp  ! Exponent of Nr

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s               [kg kg^-1]
        REAL:: sigma_s1    ! Plume 1 standard deviation of s    [kg kg^-1]
        REAL:: mu_s2       ! Plume 2 average of s               [kg kg^-1]
        REAL:: sigma_s2    ! Plume 2 standard deviation of s    [kg kg^-1]
        REAL:: mu_rr       ! Grid-box average of rr             [kg kg^-1]
        REAL:: sigma_rr    ! Grid-box standard deviation of rr  [kg kg^-1]
        REAL:: mu_Nr       ! Grid-box average of Nr             [kg^-1]
        REAL:: sigma_Nr    ! Grid-box standard deviation of Nr  [kg^-1]
        REAL:: corr_srr    ! Correlation of s and rr            []
        REAL:: corr_sNr    ! Correlation of s and Nr            []
        REAL:: corr_rrNr   ! Correlation of rr and Nr           []

        REAL:: Tl_1, Tl_2, T_1, T_2, rsl_1, rsl_2, Beta_T1, Beta_T2
        REAL:: plume_1_constants, plume_2_constants

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              IF ( Supsat < 0.0 ) THEN
              ! The air is not saturated, therefore evaporation can take place.

                 cond_evap_rrm =
     .            3.0 * C_evap * G_T_p( T_in_K, press )
     .                * ( ( (4.0/3.0)*pi*rho_lw )**(2.0/3.0) )
     .                * (rrm**(1.0/3.0)) * (Nrm**(2.0/3.0)) * Supsat

              ELSE
              ! The air is saturated, so there is no evaporation.
              ! In the HOC code, any condensation is added into cloud water
              ! instead of rain water.

                 cond_evap_rrm = 0.0

              ENDIF

           ELSE

              ! If either rrm or Nrm are 0.
              cond_evap_rrm = 0.0

           ENDIF

        ELSEIF ( .NOT. local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              !!! Define plume constants for each plume

              ! Gaussian "plume" 1
              Tl_1 = thl1 * exner !((press/p0)**kappa)

              T_1 = Tl_1 + (Lv/Cp)*rc1

              rsl_1 = sat_mixrat_liq(press, Tl_1)

              Beta_T1 = (Rd/Rv) * ( Lv/(Rd*Tl_1) ) 
     .                          * ( Lv/(Cp*Tl_1) )

              plume_1_constants = 
     .           3.0 * C_evap * G_T_p( T_1, press )
     .               * ( ( (4.0/3.0)*pi*rho_lw )**(2.0/3.0) )
     .               * ( (1.0 + Beta_T1*rsl_1) / rsl_1 )

              ! Gaussian "plume" 2
              Tl_2 = thl2 * exner !((press/p0)**kappa)

              T_2 = Tl_2 + (Lv/Cp)*rc2

              rsl_2 = sat_mixrat_liq(press, Tl_2)

              Beta_T2 = (Rd/Rv) * ( Lv/(Rd*Tl_2) ) 
     .                          * ( Lv/(Cp*Tl_2) )

              plume_2_constants = 
     .           3.0 * C_evap * G_T_p( T_2, press )
     .               * ( ( (4.0/3.0)*pi*rho_lw )**(2.0/3.0) )
     .               * ( (1.0 + Beta_T2*rsl_2) / rsl_2 )

              ! Exponents on s, rr, and Nr, respectively.
              alpha_exp = 1.0
              beta_exp  = (1.0/3.0)
              gamma_exp = (2.0/3.0)

              ! "s" is a truncated double Gaussian.
              mu_s1 = s1
              sigma_s1 = ss1
              mu_s2 = s2
              sigma_s2 = ss2

              ! rr is distributed Lognormally.
              mu_rr = rrm
!              sigma_rr = SQRT(rrp2)
              sigma_rr = rrm * SQRT(rrp2_rrm2)

              ! Nr is distributed Lognormally.
              mu_Nr = Nrm
!              sigma_Nr = SQRT(Nrp2)
              sigma_Nr = Nrm * SQRT(Nrp2_Nrm2)

              ! Correlations.
              corr_srr = corr_srr_NL
              corr_sNr = corr_sNr_NL
              corr_rrNr = corr_rrNr_LL


              cond_evap_rrm = 
     .         ( a ) * plume_1_constants
     .         * PDF_TRIVAR_2G_LN_LN ( mu_s1, mu_rr, mu_Nr,
     .                                 sigma_s1, sigma_rr, sigma_Nr,
     .                                 corr_srr, corr_sNr, corr_rrNr,
     .                                 alpha_exp, beta_exp, gamma_exp )
     .        +
     .         (1-a) * plume_2_constants
     .         * PDF_TRIVAR_2G_LN_LN ( mu_s2, mu_rr, mu_Nr,
     .                                 sigma_s2, sigma_rr, sigma_Nr,
     .                                 corr_srr, corr_sNr, corr_rrNr,
     .                                 alpha_exp, beta_exp, gamma_exp )

           ELSE

              ! If either rrm or Nrm are 0.
              cond_evap_rrm = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION cond_evap_rrm

!===============================================================================

        FUNCTION cond_evap_Nrm( cond_rrm, Nrm, rrm )

        IMPLICIT NONE

        REAL, INTENT(IN):: rrm      ! [kg kg^-1]
        REAL, INTENT(IN):: Nrm      ! [kg^-1]
        REAL, INTENT(IN):: cond_rrm ! [kg kg^-1 s^-1]
        REAL:: cond_evap_Nrm        ! [kg^-1 s^-1]

        IF (rrm > 0.0 .AND. Nrm > 0.0) THEN

           cond_evap_Nrm = ( Nrm / rrm ) * cond_rrm

        ELSE

           cond_evap_Nrm = 0.0

        ENDIF

        RETURN
        END FUNCTION cond_evap_Nrm

!===============================================================================
        !
        ! FUNCTION autoconv_rrm
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! This function calculates the autoconversion of cloud water to
        ! rain water based on Khairoutdinov and Kogan (2000) equation 29.
        ! That equation is:
        !
        ! (drr/dt)auto = 1350 * rc^2.47 * NcV^-1.79
        !
        ! where NcV is in concentration per cm^3.  Since there are 10^-6
        ! cubic meters in one cubic centimenter, the equation becomes:
        !
        ! (drr/dt)auto = 1350 * rc^2.47 * ((10^-6)NcV)^-1.79
        !
        ! simplified:
        !
        ! (drr/dt)auto = 1350 * (10^-6)^-1.79 * rc^2.47 * NcV^-1.79
        !
        ! This turns into:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rc^2.47 * NcV^-1.79
        !
        ! with NcV in concentration per m^3.
        !
        ! For our model, concentration per mass, Nc, is used instead of 
        ! concentration per volume, NcV.  We end up with:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rc^2.47 * (rho*Nc)^-1.79
        !
        ! simplified:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rho^-1.79 * rc^2.47 * Nc^-1.79
        !
        ! with Nc in concentration per kg.
        !
        ! This equation becomes:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rho^-1.79 * rc^alpha * Nc^beta
        !
        ! where alpha = 2.47 and beta = -1.79
        !
        ! s (extended liquid water specific humidity) is equal to rc anywhere
        ! that rc is greater than 0.
        ! 
        ! (drr/dt)auto = 1350 * 10^10.74 * rho^-1.79 * s^alpha * Nc^beta
        !
        ! When applied to a Probability Density Function (PDF), it becomes:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rho^-1.79 * Y1^alpha * Y2^beta
        !
        ! Y1 is used for s.  It is distributed as a double Gaussian with a
        ! domain from zero to infinity.  It is truncated at zero because we 
        ! are only interested in cases where there is cloud water.
        ! Y2 is used for Nc.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION autoconv_rrm( rcm, Ncm, s1, ss1, 
     .                         s2, ss2, a, rhot,
     .                         Ncp2_Ncm2, corr_sNc_NL )

        USE constants, only:
     .      Nc_tol, ! Variable(s)
     .      rc_tol
        USE model_flags, only:
     .      local_kk ! Variable(s)

        implicit none

        ! Input variables.
        REAL, INTENT(IN):: rcm         ! Grid-box average rcm     [kg kg^-1]
!        REAL, INTENT(IN):: rcp2        ! Grid-box rc variance     [kg^2 kg^-2]
        REAL, INTENT(IN):: Ncm         ! Grid-box average Ncm     [kg^-1]
!        REAL, INTENT(IN):: Ncp2        ! Grid-box Nc variance     [kg^-2]
        REAL, INTENT(IN):: s1          ! Plume 1 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss1         ! Plume 1 sigma s1 (not sigma^2 s1)
                                       !                          [kg kg^-1]
        REAL, INTENT(IN):: s2          ! Plume 2 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss2         ! Plume 2 sigma s2 (not sigma^2 s2)
                                       !                          [kg kg^-1]
        REAL, INTENT(IN):: a           ! Relative weight of each individual
                                       ! Gaussian "plume."        []
        REAL, INTENT(IN):: rhot        ! Grid-box average density (t-level)
                                       !                          [kg m^-3]
        REAL, INTENT(IN):: Ncp2_Ncm2   ! Ncp2/Ncm^2               []
        REAL, INTENT(IN):: corr_sNc_NL ! Correlation of s and Nc  []

        ! Output variables.
        REAL:: autoconv_rrm  ! [kg kg^-1 s^-1]

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s1
        REAL:: beta_exp   ! Exponent of Nc

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s1              [kg kg^-1]
        REAL:: sigma_s1    ! Plume 1 standard deviation of s1   [kg kg^-1]
        REAL:: mu_s2       ! Plume 2 average of s2              [kg kg^-1]
        REAL:: sigma_s2    ! Plume 2 standard deviation of s2   [kg kg^-1]
        REAL:: mu_Nc       ! Grid-box average of Nc             [kg^-1]
        REAL:: sigma_Nc    ! Grid-box standard deviation of Nc  [kg^-1]
        REAL:: corr_sNc    ! Correlation of s and Nc            []

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rcm > rc_tol .AND. Ncm > Nc_tol ) THEN

              autoconv_rrm = 7.4188E13 * rhot**(-1.79) 
     .                          * rcm**2.47 * Ncm**(-1.79)

           ELSE

              ! If either rcm or Ncm are 0.
              autoconv_rrm = 0.0

           ENDIF

        ELSEIF ( .NOT. local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rcm > rc_tol .AND. Ncm > Nc_tol ) THEN

              ! Exponents on s and Nc, respectively.
              alpha_exp = 2.47
              beta_exp  = -1.79

              ! "s" is a truncated double Gaussian.
              mu_s1 = s1
              sigma_s1 = ss1
              mu_s2 = s2
              sigma_s2 = ss2

              ! Nc is distributed Lognormally.
              mu_Nc = Ncm
!              sigma_Nc = SQRT(Ncp2)
              sigma_Nc = Ncm * SQRT(Ncp2_Ncm2)

              ! Correlations.
              corr_sNc = corr_sNc_NL


              autoconv_rrm = 7.4188E13 * rhot**(-1.79) * (
     .             ( a ) 
     .           * PDF_BIVAR_2G_LN ( mu_s1, mu_Nc, sigma_s1, sigma_Nc,
     .                               corr_sNc, alpha_exp, beta_exp )
     .         +   (1-a)
     .           * PDF_BIVAR_2G_LN ( mu_s2, mu_Nc, sigma_s2, sigma_Nc,
     .                               corr_sNc, alpha_exp, beta_exp )
     .                                                   )

           ELSE

              ! If either rcm or Ncm are 0.
              autoconv_rrm = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION autoconv_rrm

!===============================================================================

        FUNCTION autoconv_Nrm( auto_rrm )

        USE constants, only:
     .      rho_lw, ! Variable(s)
     .      pi

        IMPLICIT NONE

        REAL, INTENT(IN):: auto_rrm  ! [kg kg^-1 s^-1]
        REAL:: autoconv_Nrm          ! [kg^-1 s^-1]

        autoconv_Nrm = auto_rrm /
     .                 ( ( (4.0/3.0)*pi*rho_lw ) * (r_0**3) )

        RETURN
        END FUNCTION autoconv_Nrm

!===============================================================================
        !
        ! FUNCTION accretion_rrm
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! This function calculates the accretion of rain water based on
        ! Khairoutdinov and Kogan (2000) equation 33. That equation is:
        !
        ! (drr/dt)accr = 67 * (rc*rr)^1.15
        !
        ! the equation becomes:
        !
        ! (drr/dt)accr = 67 * rc^1.15 * rr^1.15
        !
        ! This turns into:
        !
        ! (drr/dt)accr = 67 * rc^alpha * rr^beta
        !
        ! where alpha = 1.15 and beta = 1.15
        !
        ! s (extended liquid water specific humidity) is equal to rc anywhere
        ! that rc is greater than 0.
        ! 
        ! (drr/dt)accr = 67 * s^alpha * rr^beta
        !
        ! When applied to a Probability Density Function (PDF), it becomes:
        !
        ! (drr/dt)accr = 67 * Y1^alpha * Y2^beta
        !
        ! Y1 is used for s.  It is distributed as a double Gaussian with a
        ! domain from zero to infinity.  It is truncated at zero because we 
        ! are only interested in cases where there is cloud water.
        ! Y2 is used for rr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Notes:  For our code using the local version, we use KK(2000)
        !         equation 33 exactly, as listed above.  We do not mess around
        !         with changing anything.
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION accretion_rrm( rcm, rrm, s1, ss1,
     .                          s2, ss2, a, 
     .                          rrp2_rrm2, corr_srr_NL )

        USE constants, only:
     .      rr_tol, ! Variable(s)
     .      rc_tol

        USE model_flags, only:
     .      local_kk ! Variable(s)

        implicit none

        ! Input variables.
        REAL, INTENT(IN):: rcm         ! Grid-box average rcm     [kg kg^-1]
!        REAL, INTENT(IN):: rcp2        ! Grid-box rc variance     [kg^2 kg^-2]
        REAL, INTENT(IN):: rrm         ! Grid-box average rrm     [kg kg^-1]
!        REAL, INTENT(IN):: rrp2        ! Grid-box rr variance     [kg^2 kg^-2]
        REAL, INTENT(IN):: s1          ! Plume 1 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss1         ! Plume 1 sigma s1 (not sigma^2 s1)
                                       !                          [kg kg^-1]
        REAL, INTENT(IN):: s2          ! Plume 2 average s        [kg kg^-1]
        REAL, INTENT(IN):: ss2         ! Plume 2 sigma s2 (not sigma^2 s2)
                                       !                          [kg kg^-1]
        REAL, INTENT(IN):: a           ! Relative weight of each individual
                                       ! Gaussian "plume."        []
        REAL, INTENT(IN):: rrp2_rrm2   ! rrp2/rrm^2               []
        REAL, INTENT(IN):: corr_srr_NL ! Correlation of s and rr  []

        ! Output variables.
        REAL:: accretion_rrm  ! [kg kg^-1 s^-1]

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s
        REAL:: beta_exp   ! Exponent of rr

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s               [kg kg^-1]
        REAL:: sigma_s1    ! Plume 1 standard deviation of s    [kg kg^-1]
        REAL:: mu_s2       ! Plume 2 average of s               [kg kg^-1]
        REAL:: sigma_s2    ! Plume 2 standard deviation of s    [kg kg^-1]
        REAL:: mu_rr       ! Grid-box average of rr             [kg kg^-1]
        REAL:: sigma_rr    ! Grid-box standard deviation of rr  [kg kg^-1]
        REAL:: corr_srr    ! Correlation of s and rr            []

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rcm > rc_tol .AND. rrm > rr_tol ) THEN

              accretion_rrm = 67.0 * (rcm*rrm)**1.15

           ELSE

              ! If either rcm or rrm are 0.
              accretion_rrm = 0.0

           ENDIF

        ELSEIF ( .NOT. local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rcm > rc_tol .AND. rrm > rr_tol ) THEN

              ! Exponents on s and rr, respectively.
              alpha_exp = 1.15
              beta_exp  = 1.15

              ! "s" is a truncated double Gaussian.
              mu_s1 = s1
              sigma_s1 = ss1
              mu_s2 = s2
              sigma_s2 = ss2

              ! rr is distributed Lognormally.
              mu_rr = rrm
!              sigma_rr = SQRT(rrp2)
              sigma_rr = rrm * SQRT(rrp2_rrm2)

              ! Correlations.
              corr_srr = corr_srr_NL


              accretion_rrm = 67.0 * (
     .             ( a ) 
     .           * PDF_BIVAR_2G_LN ( mu_s1, mu_rr, sigma_s1, sigma_rr,
     .                               corr_srr, alpha_exp, beta_exp )
     .         +   (1-a)
     .           * PDF_BIVAR_2G_LN ( mu_s2, mu_rr, sigma_s2, sigma_rr,
     .                               corr_srr, alpha_exp, beta_exp )
     .                               )

           ELSE

              ! If either rcm or rrm are 0.
              accretion_rrm = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION accretion_rrm

!===============================================================================

        FUNCTION G_T_p( T_in_K, press )

        USE constants, only:
     .      ep, ! Variable(s)
     .      rho_lw,
     .      Lv,
     .      Rv
        USE saturation, only:
     .      sat_mixrat_liq ! Procedure(s)

        
        IMPLICIT NONE


        ! Here we compute G(T,p) as in KK (17) 
        !      and Eq. (7.17) of Rogers and Yau 1989



        ! Input      
        REAL, INTENT(IN):: T_in_K   ! [K]
        REAL, INTENT(IN):: press  ! [Pa]

        ! Ouput
        REAL:: G_T_p  ! [m^2 s^-1]

        ! Internal
        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: esat, rs, Celsius

        Celsius = T_in_K - 273.16

        Ka = (5.69 + 0.017*Celsius)*0.00001  ! Ka in cal./(cm.*sec.*C)
        Ka = 4.1868*100.0*Ka  ! Ka in J./(m.*sec.*K)

        Dv = 0.221*((T_in_K/273.16)**1.94)*(101325.0/press)  
                                ! Dv in (cm.^2)/sec.  ! .221 is correct.
        Dv = Dv/10000.0  ! Dv in (m.^2)/sec.

        rs = sat_mixrat_liq(press, T_in_K)
        esat = (press*rs)/(ep + rs)

        Fk = (Lv/(Rv*T_in_K) - 1.0)*(Lv*rho_lw)/(Ka*T_in_K)
        Fd = (rho_lw*Rv*T_in_K)/(Dv*esat)

        G_T_p = 1.0/(Fk + Fd)

        RETURN

        END FUNCTION G_T_p

!===============================================================================
        !
        ! FUNCTION PDF_TRIVAR_2G_LN_LN
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! TRIVARIATE DOUBLE GAUSSIAN -- LOG-NORMAL -- LOG-NORMAL PDF
        !
        ! This function helps solve for the expression:
        !
        ! overbar{ Y1^alpha Y2^beta Y3^gamma }; where:
        !
        ! Y1 is distributed as a Double Gaussian with a domain from negative
        ! infinity to zero.
        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
        ! Y3 is distributed as Log-normal with a domain from zero to infinity.
        !
        ! Y1, Y2, and Y3 are then transformed into X1, X2, and X3 respectively.
        !
        ! X1 is distributed as a Double Gaussian with a domain from negative
        ! infinity to zero.  X1 = Y1
        ! X2 is distributed as Normal (Gaussian) with a domain from negative
        ! infinity to infinity.  X2 = LN( Y2 )
        ! X3 is distributed as Normal (Gaussian) with a domain from negative
        ! infinity to infinity.  X3 = LN( Y3 )
        !
        ! The resulting equation that needs to be solved for is:
        !
        ! overbar{ Y1^alpha Y2^beta Y3^gamma } =
        ! INT(-inf:0) INT(-inf:inf) INT(-inf:inf)
        !             X1^alpha exp( beta*X2 + gamma*X3 ) P(X1,X2,X3) dX3 dX2 dX1
        !
        ! since X1 is a double Gaussian, P(X1,X2,X3) =  ( a ) P_1(X1,X2,X3) 
        !                                             + (1-a) P_2(X1,X2,X3)
        !
        ! where "a" is a constant and is the relative weight of each individual
        ! Gaussian "plume".
        !
        ! P_1(X1,X2,X3) and P_2(X1,X2,X3) are simply the equation for a
        ! Trivariate Gaussian Distribution.  The only difference between the
        ! two is the dependence on each individual plume for the plume average
        ! and plume variance for the X1 variable.  There is also a dependence
        ! on each individual plume for the within-the-plume-correlation 
        ! (a.k.a. intra-Gaussian correlation) of X1 with any other variable.
        !
        ! The Probability Density Function for a Trivariate Normal Distribution:
        !
        ! P_i(X1,X2,X3) = 
        !
        ! 1 /
        !     {
        !       (2*PI)^(3/2) * sigmaX1i * sigmaX2 * sigmaX3 
        !       * SQRT[ 1 - ( corrX1iX2^2 + corrX1iX3^2 + corrX2X3^2 ) 
        !                 + 2*corrX1iX2*corrX1iX3*corrX2X3 ]
        !     }
        ! * EXP{ -(1/2) * phi }
        !
        ! where,
        !
        ! phi =
        !
        ! 1 /
        !     [
        !       1 - ( corrX1iX2^2 + corrX1iX3^2 + corrX2X3^2 ) 
        !         + 2*corrX1iX2*corrX1iX3*corrX2X3
        !     ]
        ! * {
        !       [ (1-corrX2X3^2) / sigmaX1i^2 ] * ( X1i - muX1i )^2
        !     + [ (1-corrX1iX3^2) / sigmaX2^2 ] * ( X2 - muX2 )^2
        !     + [ (1-corrX1iX2^2) / sigmaX3^2 ] * ( X3 - muX3 )^2
        !     + [ 2*(corrX1iX3*corrX2X3-corrX1iX2) / (sigmaX1i*sigmaX2) ]
        !       * ( X1i - muX1i ) * ( X2 - muX2 )
        !     + [ 2*(corrX1iX2*corrX2X3-corrX1iX3) / (sigmaX1i*sigmaX3) ]
        !       * ( X1i - muX1i ) * ( X3 - muX3 )
        !     + [ 2*(corrX1iX2*corrX1iX3-corrX2X3) / (sigmaX2*sigmaX3) ]
        !       * ( X2 - muX2 ) * ( X3 - muX3 )
        !   }
        !
        ! Definitions:
        !
        ! muX1i:  The plume average of the X1 (or Y1) variable.
        ! muX2:   The average of the X2 (or ln Y2) variable.
        ! muX3:   The average of the X3 (or ln Y3) variable.
        !
        ! sigmaX1i:  The plume standard deviation (square-root of variance)
        !            of the X1 (or Y1) variable.
        ! sigmaX2:   The standard deviation (square-root of variance) of the
        !            X2 (or ln Y2) variable.
        ! sigmaX3:   The standard deviation (square-root of variance) of the 
        !            X3 (or ln Y3) variable.
        !
        ! corrX1iX2:  Intra-Gaussian correlation between X1 and X2 (or between
        !             Y1 and ln Y2).
        ! corrX1iX3:  Intra-Gaussian correlation between X1 and X3 (or between
        !             Y1 and ln Y3).
        ! corrX2X3:   Intra-Gaussian correlation between X2 and X3 (or between
        !             ln Y2 and ln Y3).
        !
        ! After careful (and long) triple-integration of the above equation,
        ! one gets the result for each plume of the trivariate normal
        ! -- log-normal -- log-normal equation.
        !
        ! The result for each plume is (Equation #1):
        !
        !   [ 1 / SQRT(2*PI) ] * ( -sigmaX1i )^alpha 
        ! * EXP{ muX2*beta + muX3*gamma }
        ! * EXP{ (1/2) * [  (1-corrX1iX2^2)*sigmaX2^2*beta^2 
        !                 + (1-corrX1iX3^2)*sigmaX3^2*gamma^2 
        !                 + 2*(corrX2X3-corrX1iX2*corrX1iX3)
        !                    *sigmaX2*beta*sigmaX3*gamma
        !                ]
        !      }
        ! * EXP{  (1/4) * [  ( muX1i / sigmaX1i ) 
        !                  + corrX1iX2*sigmaX2*beta
        !                  + corrX1iX3*sigmaX3*gamma
        !                 ]^2
        !       - ( muX1i / sigmaX1i ) *
        !                 [  ( muX1i / sigmaX1i ) 
        !                  + corrX1iX2*sigmaX2*beta
        !                  + corrX1iX3*sigmaX3*gamma
        !                 ]
        !       + (1/2) * ( muX1i^2 / sigmaX1i^2 )
        !      }
        ! * GAMMA_FNC { alpha + 1 }
        ! * PARAB_CYL_FNC_[-(alpha+1)] {  ( muX1i / sigmaX1i ) 
        !                               + corrX1iX2*sigmaX2*beta
        !                               + corrX1iX3*sigmaX3*gamma
        !                              }
        !
        ! The sum of the weighted results for each plume give us the answer for:
        !
        ! overbar{ Y1^alpha Y2^beta Y3^gamma }
        !
        ! The purpose of this function is to give results for each plume.
        !
        ! Brian Griffin.  November 4, 2006.
        !
        !-----------------------------------------------------------------------
        !
        ! NOTE -- special case if sigmaX1i = 0
        ! ------------------------------------
        !
        ! The equation above is only a valid answer if sigmaX1i is not equal
        ! to 0.  The only variable that can be found in the denominator of any
        ! term or factor of the equation is sigmaX1i.  This is also the only
        ! case where the equation above is not a valid answer.
        !
        ! When sigmaX1i = 0, the Gaussian for X1i becomes a "spike" at muX1i
        ! with no variance.  This is also known as a Delta function.  The
        ! domain of X1i is from negative infinity to zero.  If sigmaX1i = 0
        ! and muX1i >= 0, then the entire Gaussian lies outside the domain,
        ! and the result is 0 for that individual plume.
        !
        ! However, if sigmaX1i = 0 and muX1i < 0, then the entire Gaussian
        ! lies within the domain.  A close approximation for the parabolic
        ! cylinder function is made.  This approximation is only valid in
        ! cases where the input into the parabolic cylinder function is
        ! extremely large in magnitude and very much greater in magnitude than
        ! the order of the parabolic cylinder function -- as it is in this case
        ! where sigmaX1i = 0 means that the magnitude of the input is infinite
        ! and very much greater than the magnitude of the order of the function
        ! which comes to | - (alpha + 1) | = alpha + 1.
        !
        ! The result for each plume is (Equation #2):
        !
        !   ( muX1i )^alpha 
        ! * EXP{ muX2*beta + muX3*gamma }
        ! * EXP{ (1/2) * [  (1-corrX1iX2^2)*sigmaX2^2*beta^2 
        !                 + (1-corrX1iX3^2)*sigmaX3^2*gamma^2 
        !                 + 2*(corrX2X3-corrX1iX2*corrX1iX3)
        !                    *sigmaX2*beta*sigmaX3*gamma
        !                ]
        !      }
        !
        ! One also needs to note an issue with numerical error.  If sigmaX1i
        ! is very small, but not quite 0, one would still be able to compute
        ! the answer mathematically using the equation for sigmaX1i /= 0.
        ! However, ( muX1i / sigmaX1i ) becomes very large.  If it becomes
        ! large enough, the result of the parabolic cylinder function becomes
        ! greater than the computer can represent.  The greatest value that
        ! can be represented by a DOUBLE PRECISION variable is of the order
        ! of 10^308.  Therefore, a sufficiently small sigmaX1i must be treated
        ! the same as if sigmaX1i = 0.  A test is needed to determine if the
        ! computer can represent the number or not.  Secondly, there is a case
        ! where sigmaX1i^2 is in the denominator.  While sigmaX1i may not equal
        ! 0, it may be small enough so that sigmaX1i^2 can not be represented
        ! by the computer, thereby returning a value of 0 for sigmaX1i^2.
        ! A tolerance value for sigmaX1i is used so that sigmaX1i^2 does not
        ! return a value of 0 when sigmaX1i /= 0.
        !
        ! In Summary:
        ! 
        ! When sigmaX1i /= 0 (or is not sufficiently small):  Equation #1
        ! When sigmaX1i = 0 (or is sufficiently small)
        !                                     and muX1i < 0:  Equation #2
        ! When sigmaX1i = 0 (or is sufficiently small)
        !                                    and muX1i >= 0:  0
        !
        ! In Equation #1, the factor ( -sigmaX1i )^alpha restricts the value
        ! of alpha so that the factor does not become one of an i-term.
        ! Ex: alpha = (1/2) is not allowed.  alpha = 1 is fine.
        !
        ! Brian Griffin.  November 10, 2006.
        !
        !-----------------------------------------------------------------------
        !
        ! NOTE #2 -- variable notation
        ! ----------------------------
        !
        ! In this particular case, Y1 stands for "s," Y2 stands for "rr," and
        ! Y3 stands for "Nr."  X1 stands for a Gaussianized "s," although s
        ! already is a Gaussian, so nothing is changed.  X2 stands for a
        ! Gaussianized "rr" and X3 stands for a Gaussianized "Nr."
        !
        ! Furthermore, mu stands for an average, sigma for a standard deviation,
        ! corr for a correlation. The letter "i" (by either X1, Y1, or "s")
        ! stands for the individual Gaussian plume (1 or 2).
        !
        ! While the equations written above are in a more generalized form, the
        ! actual code written below is written in terms of the actual names
        ! of the individual variables in order to decrease confusion and to be
        ! comparable with the notes by Griffin and Larson.
        !
        ! Brian Griffin.  November 18, 2006.
        !
        !-----------------------------------------------------------------------


        FUNCTION PDF_TRIVAR_2G_LN_LN ( mu_si, mu_rr, mu_Nr,
     .                                 sigma_si, sigma_rr, sigma_Nr,
     .                                 corr_sirr, corr_siNr, corr_rrNr,
     .                                 alpha_exp, beta_exp, gamma_exp )

        USE constants, only:
     .      pi ! Variable(s)
        
        USE parabolic, ONLY: 
     .      gamma ! Variable(s)

        IMPLICIT NONE

        ! Input variables.
        REAL, INTENT(IN):: mu_si     ! Plume (i) average of s
        REAL, INTENT(IN):: mu_rr     ! Average of rr
        REAL, INTENT(IN):: mu_Nr     ! Average of Nr
        REAL, INTENT(IN):: sigma_si  ! Plume (i) standard deviation of s
        REAL, INTENT(IN):: sigma_rr  ! Standard deviation of rr
        REAL, INTENT(IN):: sigma_Nr  ! Standard deviation of Nr
        REAL, INTENT(IN):: corr_sirr ! Intra-Gaussian correlation between s, rr
        REAL, INTENT(IN):: corr_siNr ! Intra-Gaussian correlation between s, Nr
        REAL, INTENT(IN):: corr_rrNr ! Intra-Gaussian correlation between rr, Nr
        REAL, INTENT(IN):: alpha_exp ! Exponent associated with "s" variable.
        REAL, INTENT(IN):: beta_exp  ! Exponent associated with rr variable.
        REAL, INTENT(IN):: gamma_exp ! Exponent associated with Nr variable.
        
        ! Output variable.
        REAL:: PDF_TRIVAR_2G_LN_LN

        ! Local variables.

        ! Variables Gaussianized and Converted.
        REAL:: mu_siG      ! Plume (i) average of s(G)
        REAL:: mu_rrG      ! Average of rr(G)
        REAL:: mu_NrG      ! Average of Nr(G)
        REAL:: sigma_siG   ! Plume (i) standard deviation of s(G)
        REAL:: sigma_rrG   ! Standard deviation of rr(G)
        REAL:: sigma_NrG   ! Standard deviation of Nr(G)
        REAL:: corr_sirrG  ! Intra-Gaussian correlation between s(G), rr(G)
        REAL:: corr_siNrG  ! Intra-Gaussian correlation between s(G), Nr(G)
        REAL:: corr_rrNrG  ! Intra-Gaussian correlation between rr(G), Nr(G)

        DOUBLE PRECISION:: gamma_fnc_input
        DOUBLE PRECISION:: parab_cyl_fnc_ord
        DOUBLE PRECISION:: scci
        DOUBLE PRECISION:: test
        LOGICAL:: use_eq_1
!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308.0
!       Found that above is not valid on most compilers -dschanen
        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308
        REAL, PARAMETER:: sigs_tol = 10.0**(-18)

        !----- Section #1 ------------------------------------------------------

        !!!!!!!!!! Convert all means, standard deviations, and correlations
        !!!!!!!!!! to Gaussian terms.

        ! si is a truncated double Gaussian.
        ! Therefore, si is still a truncated double Gaussian.
        ! It does not change.
        mu_siG = mu_si
        sigma_siG = sigma_si

        ! rr is a Lognormal.  It is converted to rr(G), which is a Gaussian.
        ! rr(G) = ln rr
        mu_rrG = LOG(  mu_rr * (  1.0 
     .                          + ( (sigma_rr**2.0) / (mu_rr**2.0) ) 
     .                         )**(-1.0/2.0)
     .              )
        sigma_rrG = SQRT( LOG(  1.0 
     .                        + ( (sigma_rr**2.0) / (mu_rr**2.0) )
     .                       )
     .                  )

        ! Nr is a Lognormal.  It is converted to Nr(G), which is a Gaussian.
        ! Nr(G) = ln Nr
        mu_NrG = LOG(  mu_Nr * (  1.0 
     .                          + ( (sigma_Nr**2.0) / (mu_Nr**2.0) ) 
     .                         )**(-1.0/2.0)
     .              )
        sigma_NrG = SQRT( LOG(  1.0
     .                        + ( (sigma_Nr**2.0) / (mu_Nr**2.0) )
     .                       )
     .                  )

        !!! Intra-Gaussian correlations.

        ! corr_sirr is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_sirrG, which is a correlation between
        ! two Gaussians.
        corr_sirrG = corr_sirr
     .                   * SQRT( EXP(sigma_rrG**2.0) - 1.0 ) / sigma_rrG
        ! corr_siNr is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_siNrG, which is a correlation between
        ! two Gaussians.
        corr_siNrG = corr_siNr
     .                   * SQRT( EXP(sigma_NrG**2.0) - 1.0 ) / sigma_NrG
        ! corr_rrNr is a correlation between two Log-normals.  It must be 
        ! converted to corr_rrNrG, which is a correlation between two Gaussians.
        corr_rrNrG = LOG( 1.0 + corr_rrNr
     .                             * SQRT( EXP(sigma_rrG**2.0) - 1.0 )
     .                             * SQRT( EXP(sigma_NrG**2.0) - 1.0 )
     .                  ) / ( sigma_rrG * sigma_NrG )

        !----- Section #2 ------------------------------------------------------

        !!!!!!!!!! Use the appropriate equation and solve.

        ! Initialize logical to false.
        use_eq_1 = .false.

        ! Input to the gamma function.
        gamma_fnc_input = alpha_exp + 1.0

        ! Order of the parabolic cylinder function.
        parab_cyl_fnc_ord = - ( alpha_exp + 1.0 )

        ! If sigma_siG = 0, then equation #1 cannot be used.
        ! A tolerance value (a small number) is used instead of zero in order
        ! to prevent numerical errors.
        IF ( sigma_siG < sigs_tol ) THEN
           use_eq_1 = .false.
        ELSE
           ! Input to the parabolic cylinder function.
           scci =  (mu_siG/sigma_siG)
     .           + corr_sirrG*sigma_rrG*beta_exp 
     .           + corr_siNrG*sigma_NrG*gamma_exp

           ! Test to see whether the value of sigma_siG is sufficiently small
           ! enough to cause the parabolic cylinder function to produce a
           ! result to large to be represented numerically by the computer.
           test = Dv_fnc( parab_cyl_fnc_ord, scci )

           IF ( test >= 0.0d0 .AND. test < limit ) THEN
              use_eq_1 = .true.
           ELSE
              use_eq_1 = .false.
           ENDIF

        ENDIF


        IF ( use_eq_1 ) THEN

           ! Case where sigma_siG is not equal to 0 (and sigma_siG is not
           ! sufficiently small enough to cause the parabolic cylinder function
           ! to produce a result greater than what the computer can represent 
           ! numerically).

           PDF_TRIVAR_2G_LN_LN = real(
     .       (1.0/SQRT(2.0*pi)) * (-sigma_siG)**alpha_exp 
     .     * EXP( mu_rrG*beta_exp + mu_NrG*gamma_exp )
     .     * EXP( (1.0/2.0) * (
     .                            ( 1.0 - corr_sirrG**2.0 )
     .                           *(sigma_rrG**2.0)*(beta_exp**2.0)
     .                         +  ( 1.0 - corr_siNrG**2.0 )
     .                           *(sigma_NrG**2.0)*(gamma_exp**2.0)
     .                         + 2.0
     .                           *( corr_rrNrG - corr_sirrG*corr_siNrG )
     .                           *sigma_rrG*beta_exp*sigma_NrG*gamma_exp
     .                        )
     .          )
     .     * EXP(   (1.0/4.0) * ( scci )**2.0
     .            - (mu_siG/sigma_siG) * ( scci )
     .            + (1.0/2.0) * ( (mu_siG**2.0) / (sigma_siG**2.0) )
     .          )
     .     * GAMMA( gamma_fnc_input )
     .     * Dv_fnc( parab_cyl_fnc_ord, scci ))

        ELSEIF ( .not. use_eq_1 .AND. mu_siG < 0.0 ) THEN

           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is less than 0.

           PDF_TRIVAR_2G_LN_LN =
     .       (mu_siG)**alpha_exp 
     .     * EXP( mu_rrG*beta_exp + mu_NrG*gamma_exp )
     .     * EXP( (1.0/2.0) * (
     .                            ( 1.0 - corr_sirrG**2.0 )
     .                           *(sigma_rrG**2.0)*(beta_exp**2.0)
     .                         +  ( 1.0 - corr_siNrG**2.0 )
     .                           *(sigma_NrG**2.0)*(gamma_exp**2.0)
     .                         + 2.0
     .                           *( corr_rrNrG - corr_sirrG*corr_siNrG )
     .                           *sigma_rrG*beta_exp*sigma_NrG*gamma_exp
     .                        )
     .          )
         ELSE
!        ELSEIF ( .not. use_eq_1 .AND. mu_siG >=  0.0 ) THEN

           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is greater than or equal to 0.

           PDF_TRIVAR_2G_LN_LN = 0.0

        ENDIF

        RETURN

        END FUNCTION PDF_TRIVAR_2G_LN_LN

!===============================================================================
        !
        ! FUNCTION PDF_BIVAR_2G_LN
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! BIVARIATE DOUBLE GAUSSIAN -- LOG-NORMAL PDF
        !
        ! This function helps solve for the expression:
        !
        ! overbar{ Y1^alpha Y2^beta }; where:
        !
        ! Y1 is distributed as a Double Gaussian with a domain from zero to
        ! infinity.
        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
        !
        ! Y1 and Y2 are then transformed into X1 and X2 respectively.
        !
        ! X1 is distributed as a Double Gaussian with a domain from zero to
        ! infinity.  X1 = Y1
        ! X2 is distributed as Normal (Gaussian) with a domain from negative
        ! infinity to infinity.  X2 = LN( Y2 )
        !
        ! The resulting equation that needs to be solved for is:
        !
        ! overbar{ Y1^alpha Y2^beta } =
        ! INT(0:inf) INT(-inf:inf)
        !             X1^alpha exp( beta*X2 ) P(X1,X2) dX2 dX1
        !
        ! since X1 is a double Gaussian, P(X1,X2) =  ( a ) P_1(X1,X2) 
        !                                          + (1-a) P_2(X1,X2)
        !
        ! where "a" is a constant and is the relative weight of each individual
        ! Gaussian "plume".
        !
        ! P_1(X1,X2) and P_2(X1,X2) are simply the equation for a
        ! Bivariate Gaussian Distribution.  The only difference between the
        ! two is the dependence on each individual plume for the plume average
        ! and plume variance for the X1 variable.  There is also a dependence
        ! on each individual plume for the within-the-plume-correlation
        ! (a.k.a. intra-Gaussian correlation) of X1 with any other variable.
        !
        ! The Probability Density Function for a Bivariate Normal Distribution:
        !
        ! P_i(X1,X2) = 
        !
        ! 1 / { 2 * PI * sigmaX1i * sigmaX2 * SQRT[ 1 - corrX1iX2^2 ] }
        ! * EXP{ -(1/2) * phi }
        !
        ! where,
        !
        ! phi =
        !
        ! 1 / [ 1 - corrX1iX2^2 ]
        ! * {
        !       [ 1 / sigmaX1i^2 ] * ( X1i - muX1i )^2
        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
        !     + [ 2*corrX1iX2 / (sigmaX1i*sigmaX2) ]
        !       * ( X1i - muX1i ) * ( X2 - muX2 )
        !   }
        !
        ! Definitions:
        !
        ! muX1i:  The plume average of the X1 (or Y1) variable.
        ! muX2:   The average of the X2 (or ln Y2) variable.
        !
        ! sigmaX1i:  The plume standard deviation (square-root of variance)
        !            of the X1 (or Y1) variable.
        ! sigmaX2:   The standard deviation (square-root of variance) of the
        !            X2 (or ln Y2) variable.
        !
        ! corrX1iX2:  Intra-Gaussian correlation between X1 and X2 (or between
        !             Y1 and ln Y2).
        !
        ! After careful (and long) double-integration of the above equation,
        ! one gets the result for each plume of the bivariate normal
        ! -- log-normal equation.
        !
        ! The result for each plume is (Equation #1):
        !
        !   [ 1 / SQRT(2*PI) ] * ( sigmaX1i )^alpha 
        ! * EXP{    muX2*beta + (1/2)*sigmaX2^2*beta^2
        !        - (1/4) * [ ( muX1i / sigmaX1i ) + corrX1iX2*sigmaX2*beta ]^2
        !      }
        ! * GAMMA_FNC { alpha + 1 }
        ! * PARAB_CYL_FNC_[-(alpha+1)] { - [  ( muX1i / sigmaX1i )
        !                                   + corrX1iX2*sigmaX2*beta ]
        !                              }
        !
        ! The sum of the weighted results for each plume give us the answer for:
        !
        ! overbar{ Y1^alpha Y2^beta }
        !
        ! The purpose of this function is to give results for each plume.
        !
        ! Brian Griffin.  November 10, 2006.
        !
        !-----------------------------------------------------------------------
        !
        ! NOTE -- special case if sigmaX1i = 0
        ! ------------------------------------
        !
        ! The equation above is only a valid answer if sigmaX1i is not equal
        ! to 0.  The only variable that can be found in the denominator of any
        ! term or factor of the equation is sigmaX1i.  This is also the only
        ! case where the equation above is not a valid answer.
        !
        ! When sigmaX1i = 0, the Gaussian for X1i becomes a "spike" at muX1i
        ! with no variance.  This is also known as a Delta function.  The
        ! domain of X1i is from zero to infinity.  If sigmaX1i = 0 and 
        ! muX1i <= 0, then the entire Gaussian lies outside the domain,
        ! and the result is 0 for that individual plume.
        !
        ! However, if sigmaX1i = 0 and muX1i > 0, then the entire Gaussian
        ! lies within the domain.  A close approximation for the parabolic
        ! cylinder function is made.  This approximation is only valid in
        ! cases where the input into the parabolic cylinder function is
        ! extremely large in magnitude and very much greater in magnitude than
        ! the order of the parabolic cylinder function -- as it is in this case
        ! where sigmaX1i = 0 means that the magnitude of the input is infinite
        ! and very much greater than the magnitude of the order of the function
        ! which comes to | - (alpha + 1) | = alpha + 1.
        !
        ! The result for each plume is (Equation #2):
        !
        ! ( muX1i )^alpha * EXP{ muX2*beta + (1/2)*sigmaX2^2*beta^2 }
        !
        ! One also needs to note an issue with numerical error.  If sigmaX1i
        ! is very small, but not quite 0, one would still be able to compute
        ! the answer mathematically using the equation for sigmaX1i /= 0.
        ! However, ( muX1i / sigmaX1i ) becomes very large.  If it becomes
        ! large enough, the result of the parabolic cylinder function becomes
        ! greater than the computer can represent.  The greatest value that
        ! can be represented by a DOUBLE PRECISION variable is of the order
        ! of 10^308.  Therefore, a sufficiently small sigmaX1i must be treated
        ! the same as if sigmaX1i = 0.  A test is needed to determine if the
        ! computer can represent the number or not.  Secondly, there is a case
        ! where sigmaX1i^2 is in the denominator.  While sigmaX1i may not equal
        ! 0, it may be small enough so that sigmaX1i^2 can not be represented
        ! by the computer, thereby returning a value of 0 for sigmaX1i^2.
        ! A tolerance value for sigmaX1i is used so that sigmaX1i^2 does not
        ! return a value of 0 when sigmaX1i /= 0.
        !
        ! In Summary:
        ! 
        ! When sigmaX1i /= 0 (or is not sufficiently small):  Equation #1
        ! When sigmaX1i = 0 (or is sufficiently small)
        !                                     and muX1i > 0:  Equation #2
        ! When sigmaX1i = 0 (or is sufficiently small)
        !                                    and muX1i <= 0:  0
        !
        ! Brian Griffin.  November 10, 2006.
        !
        !-----------------------------------------------------------------------
        !
        ! NOTE #2 -- variable notation
        ! ----------------------------
        !
        ! In this particular case, Y1 stands for "s" and Y2 stands for "xx,"
        ! which can be either "Nc" (as in the case of autoconversion) or "rr"
        ! (as in the case of accretion).  X1 stands for a Gaussianized "s,"
        ! although s already is a Gaussian, so nothing is changed.  X2 stands
        ! for a Gaussianized "xx," for both "Nc" and "rr" are lognormal and
        ! need to be Gaussianized.
        !
        ! Furthermore, mu stands for an average, sigma for a standard deviation,
        ! corr for a correlation. The letter "i" (by either X1, Y1, or "s")
        ! stands for the individual Gaussian plume (1 or 2).
        !
        ! While the equations written above are in a more generalized form, the
        ! actual code written below is written in terms of the actual names
        ! of the individual variables in order to decrease confusion and to be
        ! comparable with the notes by Griffin and Larson.
        !
        ! Brian Griffin.  November 18, 2006.
        !
        !-----------------------------------------------------------------------


        FUNCTION PDF_BIVAR_2G_LN ( mu_si, mu_xx, sigma_si, sigma_xx,
     .                             corr_sixx, alpha_exp, beta_exp )

        USE constants, only:
     .      pi ! Variable(s)
        USE parabolic, ONLY: 
     .      gamma ! Variable(s)

        IMPLICIT NONE

        ! Input variables.
        REAL, INTENT(IN):: mu_si     ! Plume (i) average of s
        REAL, INTENT(IN):: mu_xx     ! Average of xx
        REAL, INTENT(IN):: sigma_si  ! Plume (i) standard deviation of s
        REAL, INTENT(IN):: sigma_xx  ! Standard deviation of Y2
        REAL, INTENT(IN):: corr_sixx ! Intra-Gaussian correlation between s, xx
        REAL, INTENT(IN):: alpha_exp ! Exponent associated with "s" variable.
        REAL, INTENT(IN):: beta_exp  ! Exponent associated with xx variable.
        
        ! Output variable.
        REAL:: PDF_BIVAR_2G_LN

        ! Local variables.

        ! Variables Gaussianized and Converted.
        REAL:: mu_siG      ! Plume (i) average of s(G)
        REAL:: mu_xxG      ! Average of xx(G)
        REAL:: sigma_siG   ! Plume (i) standard deviation of s(G)
        REAL:: sigma_xxG   ! Standard deviation of xx(G)
        REAL:: corr_sixxG  ! Intra-Gaussian correlation between s(G), xx(G)

        DOUBLE PRECISION:: gamma_fnc_input
        DOUBLE PRECISION:: parab_cyl_fnc_ord
        DOUBLE PRECISION:: sci
        DOUBLE PRECISION:: test
        LOGICAL:: use_eq_1
        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308
        REAL, PARAMETER:: sigs_tol = 10.0**(-18)

        sci = 0.0 ! Arbitrary default value
        ! Joshua Fasching June 2008

        !----- Section #1 ------------------------------------------------------

        !!!!!!!!!! Convert all means, standard deviations, and correlations
        !!!!!!!!!! to Gaussian terms.

        ! si is a truncated double Gaussian.
        ! Therefore, si(G) is still a truncated double Gaussian.
        ! It does not change.
        mu_siG = mu_si
        sigma_siG = sigma_si

        ! xx is a Lognormal.  It is converted to xx(G), which is a Gaussian.
        ! xx(G) = ln xx
        mu_xxG = LOG(  mu_xx * (  1.0 
     .                          + ( (sigma_xx**2.0) / (mu_xx**2.0) ) 
     .                         )**(-1.0/2.0)
     .              )
        sigma_xxG = SQRT( LOG(  1.0 
     .                        + ( (sigma_xx**2.0) / (mu_xx**2.0) )
     .                       )
     .                  )

        !!! Intra-Gaussian correlation.

        ! corr_sixx is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_sixxG, which is a correlation between
        ! two Gaussians.
        corr_sixxG = corr_sixx
     .                   * SQRT( EXP(sigma_xxG**2.0) - 1.0 ) / sigma_xxG

        !----- Section #2 ------------------------------------------------------

        !!!!!!!!!! Use the appropriate equation and solve.

        ! Initialize logical to false.
        use_eq_1 = .false.

        ! Input to the gamma function.
        gamma_fnc_input = alpha_exp + 1.0

        ! Order of the parabolic cylinder function.
        parab_cyl_fnc_ord = - ( alpha_exp + 1.0 )

        ! If sigma_siG = 0, then equation #1 cannot be used.
        ! A tolerance value (a small number) is used instead of zero in order
        ! to prevent numerical errors.
        IF (sigma_siG < sigs_tol ) THEN
           use_eq_1 = .false.
        ELSE
           ! Input to the parabolic cylinder function (without minus sign).
           sci =  (mu_siG/sigma_siG)
     .          + corr_sixxG*sigma_xxG*beta_exp 

           ! Test to see whether the value of sigma_siG is sufficiently small
           ! enough to cause the parabolic cylinder function to produce a
           ! result to large to be represented numerically by the computer.
           test = Dv_fnc( parab_cyl_fnc_ord, -sci )

           IF ( test >= 0.0d0 .AND. test < limit ) THEN
              use_eq_1 = .true.
           ELSE
              use_eq_1 = .false.
           ENDIF

        ENDIF


        IF ( use_eq_1 ) THEN

           ! Case where sigma_siG is not equal to 0 (and sigma_siG is not
           ! sufficiently small enough to cause the parabolic cylinder function
           ! to produce a result greater than what the computer can represent 
           ! numerically).

           PDF_BIVAR_2G_LN = real(
     .       (1.0/SQRT(2.0*pi)) * (sigma_siG)**alpha_exp 
     .     * EXP(   mu_xxG*beta_exp
     .            + (1.0/2.0)*(sigma_xxG**2.0)*(beta_exp**2.0)
     .            - (1.0/4.0)*( sci )**2.0
     .          )
     .     * GAMMA( gamma_fnc_input )
     .     * Dv_fnc( parab_cyl_fnc_ord, -sci ))

        ELSEIF ( .not. use_eq_1 .AND. mu_siG > 0.0 ) THEN

           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is greater than 0.

           PDF_BIVAR_2G_LN =
     .       (mu_siG)**alpha_exp 
     .     * EXP(   mu_xxG*beta_exp
     .            + (1.0/2.0)*(sigma_xxG**2.0)*(beta_exp**2.0) )

!        ELSEIF ( .not. use_eq_1 .AND. mu_siG <=  0.0 ) THEN
         ELSE
           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is less than or equal to 0.

           PDF_BIVAR_2G_LN = 0.0

        ENDIF

        RETURN

        END FUNCTION PDF_BIVAR_2G_LN

!===============================================================================
        !
        ! FUNCTION PDF_BIVAR_LN_LN
        !-----------------------------------------------------------------------
        !
        ! DESCRIPTION
        !------------
        !
        ! BIVARIATE LOG-NORMAL -- LOG-NORMAL PDF
        !
        ! This function helps solve for the expression:
        !
        ! overbar{ Y1^alpha Y2^beta }; where:
        !
        ! Y1 is distributed as Log-normal with a domain from zero to infinity.
        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
        !
        ! Y1 and Y2 are then transformed into X1 and X2 respectively.
        !
        ! X1 is distributed as Normal (Gaussian) with a domain from negative
        ! infinity to infinity.  X1 = LN( Y1 )
        ! X2 is distributed as Normal (Gaussian) with a domain from negative
        ! infinity to infinity.  X2 = LN( Y2 )
        !
        ! The resulting equation that needs to be solved for is:
        !
        ! overbar{ Y1^alpha Y2^beta } =
        ! INT(-inf:inf) INT(-inf:inf)
        !             exp( alpha*X1 + beta*X2 ) P(X1,X2) dX2 dX1
        !
        ! P(X1,X2) is simply the equation for a Bivariate Gaussian Distribution.
        !
        ! The Probability Density Function for a Bivariate Normal Distribution:
        !
        ! P(X1,X2) = 
        !
        ! 1 / { 2 * PI * sigmaX1 * sigmaX2 * SQRT[ 1 - corrX1X2^2 ] }
        ! * EXP{ -(1/2) * phi }
        !
        ! where,
        !
        ! phi =
        !
        ! 1 / [ 1 - corrX1X2^2 ]
        ! * {
        !       [ 1 / sigmaX1^2 ] * ( X1 - muX1 )^2
        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
        !     + [ 2*corrX1X2 / (sigmaX1*sigmaX2) ]
        !       * ( X1 - muX1 ) * ( X2 - muX2 )
        !   }
        !
        ! Definitions:
        !
        ! muX1:  The average of the X1 (or ln Y1) variable.
        ! muX2:  The average of the X2 (or ln Y2) variable.
        !
        ! sigmaX1:  The standard deviation (square-root of variance) of the
        !           X1 (or ln Y1) variable.
        ! sigmaX2:  The standard deviation (square-root of variance) of the
        !           X2 (or ln Y2) variable.
        !
        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
        !            ln Y1 and ln Y2).
        !
        ! After careful double-integration of the above equation, one gets the
        ! result for each plume of the bivariate log-normal -- log-normal
        ! equation.
        !
        ! The result is:
        !
        ! * EXP{    muX1*alpha + muX2*beta
        !        + (1/2)*sigmaX1^2*alpha^2 + (1/2)*sigmaX2^2*beta^2
        !        + corrX1X2*sigmaX1*alpha*sigmaX2*beta
        !      }
        !
        ! The above equation gives us the answer for:
        !
        ! overbar{ Y1^alpha Y2^beta }
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        !
        ! NOTE #2 -- variable notation
        ! ----------------------------
        !
        ! In this particular case, Y1 stands for "rr" and Y2 stands for "Nr."
        ! X1 stands for a Gaussianized "rr" and X2 stands for a
        ! Gaussianized "Nr."
        !
        ! Furthermore, mu stands for an average, sigma for a standard deviation,
        ! corr for a correlation.
        !
        ! While the equations written above are in a more generalized form, the
        ! actual code written below is written in terms of the actual names
        ! of the individual variables in order to decrease confusion and to be
        ! comparable with the notes by Griffin and Larson.
        !
        ! Brian Griffin.  November 18, 2006.
        !
        !-----------------------------------------------------------------------


        FUNCTION PDF_BIVAR_LN_LN ( mu_rr, mu_Nr, sigma_rr, sigma_Nr,
     .                             corr_rrNr, alpha_exp, beta_exp   )

        USE constants, only: 
     .      pi ! Variable(s)

        IMPLICIT NONE

        ! Input variables.
        REAL, INTENT(IN):: mu_rr     ! Average of rr
        REAL, INTENT(IN):: mu_Nr     ! Average of Nr
        REAL, INTENT(IN):: sigma_rr  ! Standard deviation of rr
        REAL, INTENT(IN):: sigma_Nr  ! Standard deviation of Nr
        REAL, INTENT(IN):: corr_rrNr ! Intra-Gaussian correlation between rr, Nr
        REAL, INTENT(IN):: alpha_exp ! Exponent associated with rr variable.
        REAL, INTENT(IN):: beta_exp  ! Exponent associated with Nr variable.
        
        ! Output variable.
        REAL:: PDF_BIVAR_LN_LN

        ! Local variables.

        ! Variables Gaussianized and Converted.
        REAL:: mu_rrG      ! Average of rr(G)
        REAL:: mu_NrG      ! Average of Nr(G)
        REAL:: sigma_rrG   ! Standard deviation of rr(G)
        REAL:: sigma_NrG   ! Standard deviation of Nr(G)
        REAL:: corr_rrNrG  ! Intra-Gaussian correlation between rr(G), X2(G)

        !----- Section #1 ------------------------------------------------------

        !!!!!!!!!! Convert all means, standard deviations, and correlations
        !!!!!!!!!! to Gaussian terms.

        ! rr is a Lognormal.  It is converted to rr(G), which is a Gaussian.
        ! rr(G) = ln rr
        mu_rrG = LOG(  mu_rr * (  1.0 
     .                          + ( (sigma_rr**2.0) / (mu_rr**2.0) ) 
     .                         )**(-1.0/2.0)
     .              )
        sigma_rrG = SQRT( LOG(  1.0 
     .                        + ( (sigma_rr**2.0) / (mu_rr**2.0) )
     .                       )
     .                  )

        ! Nr is a Lognormal.  It is converted to Nr(G), which is a Gaussian.
        ! Nr(G) = ln Nr
        mu_NrG = LOG(  mu_Nr * (  1.0 
     .                          + ( (sigma_Nr**2.0) / (mu_Nr**2.0) ) 
     .                         )**(-1.0/2.0)
     .              )
        sigma_NrG = SQRT( LOG(  1.0
     .                        + ( (sigma_Nr**2.0) / (mu_Nr**2.0) )
     .                       )
     .                  )

        !!! Intra-Gaussian correlation.
        ! corr_rrNr is a correlation between two Log-normals.  It must be 
        ! converted to corr_rrNrG, which is a correlation between two Gaussians.
        corr_rrNrG = LOG( 1.0 + corr_rrNr
     .                             * SQRT( EXP(sigma_rrG**2.0) - 1.0 )
     .                             * SQRT( EXP(sigma_NrG**2.0) - 1.0 )
     .                  ) / ( sigma_rrG * sigma_NrG )

        !----- Section #2 ------------------------------------------------------

        !!!!!!!!!! Solve.

        PDF_BIVAR_LN_LN =
     .     EXP(   mu_rrG*alpha_exp + mu_NrG*beta_exp
     .          + (1.0/2.0)*(sigma_rrG**2.0)*(alpha_exp**2.0)
     .          + (1.0/2.0)*(sigma_NrG**2.0)*(beta_exp**2.0)
     .          + corr_rrNrG*sigma_rrG*alpha_exp*sigma_NrG*beta_exp
     .        )

        RETURN

        END FUNCTION PDF_BIVAR_LN_LN

!!===============================================================================
!        !
!        ! FUNCTION PDF_TRIVAR_2G_LN_LN
!        !-----------------------------------------------------------------------
!        !
!        ! DESCRIPTION
!        !------------
!        !
!        ! TRIVARIATE DOUBLE GAUSSIAN -- LOG-NORMAL -- LOG-NORMAL PDF
!        !
!        ! This function helps solve for the expression:
!        !
!        ! overbar{ Y1^alpha Y2^beta Y3^gamma }; where:
!        !
!        ! Y1 is distributed as a Double Gaussian with a domain from negative
!        ! infinity to zero.
!        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
!        ! Y3 is distributed as Log-normal with a domain from zero to infinity.
!        !
!        ! Y1, Y2, and Y3 are then transformed into X1, X2, and X3 respectively.
!        !
!        ! X1 is distributed as a Double Gaussian with a domain from negative
!        ! infinity to zero.  X1 = Y1
!        ! X2 is distributed as Normal (Gaussian) with a domain from negative
!        ! infinity to infinity.  X2 = LN( Y2 )
!        ! X3 is distributed as Normal (Gaussian) with a domain from negative
!        ! infinity to infinity.  X3 = LN( Y3 )
!        !
!        ! The resulting equation that needs to be solved for is:
!        !
!        ! overbar{ Y1^alpha Y2^beta Y3^gamma } =
!        ! INT(-inf:0) INT(-inf:inf) INT(-inf:inf)
!        !             X1^alpha exp( beta*X2 + gamma*X3 ) P(X1,X2,X3) dX3 dX2 dX1
!        !
!        ! since X1 is a double Gaussian, P(X1,X2,X3) =  ( a ) P_1(X1,X2,X3) 
!        !                                             + (1-a) P_2(X1,X2,X3)
!        !
!        ! where "a" is a constant and is the relative weight of each individual
!        ! Gaussian "plume".
!        !
!        ! P_1(X1,X2,X3) and P_2(X1,X2,X3) are simply the equation for a
!        ! Trivariate Gaussian Distribution.  The only difference between the
!        ! two is the dependence on each individual plume for the plume average
!        ! and plume variance for the X1 variable.  There is also a dependence
!        ! on each individual plume for the within-the-plume-correlation 
!        ! (a.k.a. intra-Gaussian correlation) of X1 with any other variable.
!        !
!        ! The Probability Density Function for a Trivariate Normal Distribution:
!        !
!        ! P_i(X1,X2,X3) = 
!        !
!        ! 1 /
!        !     {
!        !       (2*PI)^(3/2) * sigmaX1i * sigmaX2 * sigmaX3 
!        !       * SQRT[ 1 - ( corrX1iX2^2 + corrX1iX3^2 + corrX2X3^2 ) 
!        !                 + 2*corrX1iX2*corrX1iX3*corrX2X3 ]
!        !     }
!        ! * EXP{ -(1/2) * phi }
!        !
!        ! where,
!        !
!        ! phi =
!        !
!        ! 1 /
!        !     [
!        !       1 - ( corrX1iX2^2 + corrX1iX3^2 + corrX2X3^2 ) 
!        !         + 2*corrX1iX2*corrX1iX3*corrX2X3
!        !     ]
!        ! * {
!        !       [ (1-corrX2X3^2) / sigmaX1i^2 ] * ( X1i - muX1i )^2
!        !     + [ (1-corrX1iX3^2) / sigmaX2^2 ] * ( X2 - muX2 )^2
!        !     + [ (1-corrX1iX2^2) / sigmaX3^2 ] * ( X3 - muX3 )^2
!        !     + [ 2*(corrX1iX3*corrX2X3-corrX1iX2) / (sigmaX1i*sigmaX2) ]
!        !       * ( X1i - muX1i ) * ( X2 - muX2 )
!        !     + [ 2*(corrX1iX2*corrX2X3-corrX1iX3) / (sigmaX1i*sigmaX3) ]
!        !       * ( X1i - muX1i ) * ( X3 - muX3 )
!        !     + [ 2*(corrX1iX2*corrX1iX3-corrX2X3) / (sigmaX2*sigmaX3) ]
!        !       * ( X2 - muX2 ) * ( X3 - muX3 )
!        !   }
!        !
!        ! Definitions:
!        !
!        ! muX1i:  The plume average of the X1 (or Y1) variable.
!        ! muX2:   The average of the X2 (or ln Y2) variable.
!        ! muX3:   The average of the X3 (or ln Y3) variable.
!        !
!        ! sigmaX1i:  The plume standard deviation (square-root of variance)
!        !            of the X1 (or Y1) variable.
!        ! sigmaX2:   The standard deviation (square-root of variance) of the
!        !            X2 (or ln Y2) variable.
!        ! sigmaX3:   The standard deviation (square-root of variance) of the 
!        !            X3 (or ln Y3) variable.
!        !
!        ! corrX1iX2:  Intra-Gaussian correlation between X1 and X2 (or between
!        !             Y1 and ln Y2).
!        ! corrX1iX3:  Intra-Gaussian correlation between X1 and X3 (or between
!        !             Y1 and ln Y3).
!        ! corrX2X3:   Intra-Gaussian correlation between X2 and X3 (or between
!        !             ln Y2 and ln Y3).
!        !
!        ! After careful (and long) triple-integration of the above equation,
!        ! one gets the result for each plume of the trivariate normal
!        ! -- log-normal -- log-normal equation.
!        !
!        ! The result for each plume is (Equation #1):
!        !
!        !   [ 1 / SQRT(2*PI) ] * ( -sigmaX1i )^alpha 
!        ! * EXP{ muX2*beta + muX3*gamma }
!        ! * EXP{ (1/2) * [  (1-corrX1iX2^2)*sigmaX2^2*beta^2 
!        !                 + (1-corrX1iX3^2)*sigmaX3^2*gamma^2 
!        !                 + 2*(corrX2X3-corrX1iX2*corrX1iX3)
!        !                    *sigmaX2*beta*sigmaX3*gamma
!        !                ]
!        !      }
!        ! * EXP{  (1/4) * [  ( muX1i / sigmaX1i ) 
!        !                  + corrX1iX2*sigmaX2*beta
!        !                  + corrX1iX3*sigmaX3*gamma
!        !                 ]^2
!        !       - ( muX1i / sigmaX1i ) *
!        !                 [  ( muX1i / sigmaX1i ) 
!        !                  + corrX1iX2*sigmaX2*beta
!        !                  + corrX1iX3*sigmaX3*gamma
!        !                 ]
!        !       + (1/2) * ( muX1i^2 / sigmaX1i^2 )
!        !      }
!        ! * GAMMA_FNC { alpha + 1 }
!        ! * PARAB_CYL_FNC_[-(alpha+1)] {  ( muX1i / sigmaX1i ) 
!        !                               + corrX1iX2*sigmaX2*beta
!        !                               + corrX1iX3*sigmaX3*gamma
!        !                              }
!        !
!        ! The sum of the weighted results for each plume give us the answer for:
!        !
!        ! overbar{ Y1^alpha Y2^beta Y3^gamma }
!        !
!        ! The purpose of this function is to give results for each plume.
!        !
!        ! Brian Griffin.  November 4, 2006.
!        !
!        !-----------------------------------------------------------------------
!        !
!        ! NOTE -- special case if sigmaX1i = 0
!        ! ------------------------------------
!        !
!        ! The equation above is only a valid answer if sigmaX1i is not equal
!        ! to 0.  The only variable that can be found in the denominator of any
!        ! term or factor of the equation is sigmaX1i.  This is also the only
!        ! case where the equation above is not a valid answer.
!        !
!        ! When sigmaX1i = 0, the Gaussian for X1i becomes a "spike" at muX1i
!        ! with no variance.  This is also known as a Delta function.  The
!        ! domain of X1i is from negative infinity to zero.  If sigmaX1i = 0
!        ! and muX1i >= 0, then the entire Gaussian lies outside the domain,
!        ! and the result is 0 for that individual plume.
!        !
!        ! However, if sigmaX1i = 0 and muX1i < 0, then the entire Gaussian
!        ! lies within the domain.  A close approximation for the parabolic
!        ! cylinder function is made.  This approximation is only valid in
!        ! cases where the input into the parabolic cylinder function is
!        ! extremely large in magnitude and very much greater in magnitude than
!        ! the order of the parabolic cylinder function -- as it is in this case
!        ! where sigmaX1i = 0 means that the magnitude of the input is infinite
!        ! and very much greater than the magnitude of the order of the function
!        ! which comes to | - (alpha + 1) | = alpha + 1.
!        !
!        ! The result for each plume is (Equation #2):
!        !
!        !   ( muX1i )^alpha 
!        ! * EXP{ muX2*beta + muX3*gamma }
!        ! * EXP{ (1/2) * [  (1-corrX1iX2^2)*sigmaX2^2*beta^2 
!        !                 + (1-corrX1iX3^2)*sigmaX3^2*gamma^2 
!        !                 + 2*(corrX2X3-corrX1iX2*corrX1iX3)
!        !                    *sigmaX2*beta*sigmaX3*gamma
!        !                ]
!        !      }
!        !
!        ! One also needs to note an issue with numerical error.  If sigmaX1i
!        ! is very small, but not quite 0, one would still be able to compute
!        ! the answer mathematically using the equation for sigmaX1i /= 0.
!        ! However, ( muX1i / sigmaX1i ) becomes very large.  If it becomes
!        ! large enough, the result of the parabolic cylinder function becomes
!        ! greater than the computer can represent.  The greatest value that
!        ! can be represented by a DOUBLE PRECISION variable is of the order
!        ! of 10^308.  Therefore, a sufficiently small sigmaX1i must be treated
!        ! the same as if sigmaX1i = 0.  A test is needed to determine if the
!        ! computer can represent the number or not.  Secondly, there is a case
!        ! where sigmaX1i^2 is in the denominator.  While sigmaX1i may not equal
!        ! 0, it may be small enough so that sigmaX1i^2 can not be represented
!        ! by the computer, thereby returning a value of 0 for sigmaX1i^2.
!        ! A tolerance value for sigmaX1i is used so that sigmaX1i^2 does not
!        ! return a value of 0 when sigmaX1i /= 0.
!        !
!        ! In Summary:
!        ! 
!        ! When sigmaX1i /= 0 (or is not sufficiently small):  Equation #1
!        ! When sigmaX1i = 0 (or is sufficiently small)
!        !                                     and muX1i < 0:  Equation #2
!        ! When sigmaX1i = 0 (or is sufficiently small)
!        !                                    and muX1i >= 0:  0
!        !
!        ! In Equation #1, the factor ( -sigmaX1i )^alpha restricts the value
!        ! of alpha so that the factor does not become one of an i-term.
!        ! Ex: alpha = (1/2) is not allowed.  alpha = 1 is fine.
!        !
!        ! Brian Griffin.  November 10, 2006.
!        !
!        !-----------------------------------------------------------------------
!
!
!        FUNCTION PDF_TRIVAR_2G_LN_LN ( muY1i, muY2, muY3,
!     .                                 sigmaY1i, sigmaY2, sigmaY3,
!     .                                 corrY1iY2, corrY1iY3, corrY2Y3,
!     .                                 alpha_exp, beta_exp, gamma_exp )
!
!        USE constants
!        USE polpak_gamma, ONLY: gamma
!
!        IMPLICIT NONE
!
!        ! Input variables.
!        REAL, INTENT(IN):: muY1i     ! Plume average of Y1
!        REAL, INTENT(IN):: muY2      ! Average of Y2
!        REAL, INTENT(IN):: muY3      ! Average of Y3
!        REAL, INTENT(IN):: sigmaY1i  ! Plume standard deviation of Y1
!        REAL, INTENT(IN):: sigmaY2   ! Standard deviation of Y2
!        REAL, INTENT(IN):: sigmaY3   ! Standard deviation of Y3
!        REAL, INTENT(IN):: corrY1iY2 ! Intra-Gaussian correlation between Y1,Y2
!        REAL, INTENT(IN):: corrY1iY3 ! Intra-Gaussian correlation between Y1,Y3
!        REAL, INTENT(IN):: corrY2Y3  ! Intra-Gaussian correlation between Y2,Y3
!        REAL, INTENT(IN):: alpha_exp ! Exponent associated with Y1 variable.
!        REAL, INTENT(IN):: beta_exp  ! Exponent associated with Y2 variable.
!        REAL, INTENT(IN):: gamma_exp ! Exponent associated with Y3 variable.
!        
!        ! Output variable.
!        REAL:: PDF_TRIVAR_2G_LN_LN
!
!        ! Local variables.
!
!        ! "Y" variables Gaussianized and converted to "X" variables.
!        REAL:: muX1i     ! Plume average of X1
!        REAL:: muX2      ! Average of X2
!        REAL:: muX3      ! Average of X3
!        REAL:: sigmaX1i  ! Plume standard deviation of X1
!        REAL:: sigmaX2   ! Standard deviation of X2
!        REAL:: sigmaX3   ! Standard deviation of X3
!        REAL:: corrX1iX2 ! Intra-Gaussian correlation between X1,X2
!        REAL:: corrX1iX3 ! Intra-Gaussian correlation between X1,X3
!        REAL:: corrX2X3  ! Intra-Gaussian correlation between X2,X3
!
!        DOUBLE PRECISION:: gamma_fnc_input
!        DOUBLE PRECISION:: parab_cyl_fnc_ord
!        DOUBLE PRECISION:: parab_cyl_fnc_input
!        DOUBLE PRECISION:: test
!        LOGICAL:: use_eq_1
!!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308.0
!!       Found that above is not valid on most compilers -dschanen
!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308
!        REAL, PARAMETER:: sigX1_tol = 10.0**(-18)
!
!        !----- Section #1 ------------------------------------------------------
!
!        !!!!!!!!!! Convert all means, standard deviations, and correlations
!        !!!!!!!!!! to Gaussian terms.
!
!        ! Y1i is a truncated double Gaussian.
!        ! Therefore, X1i is still a truncated double Gaussian.
!        ! It does not change.
!        muX1i = muY1i
!        sigmaX1i = sigmaY1i
!
!        ! Y2 is a Lognormal.  It is converted to X2, which is a Gaussian.
!        ! X2 = ln Y2
!        muX2 = LOG(  muY2 * (  1.0 
!     .                       + ( (sigmaY2**2.0) / (muY2**2.0) ) 
!     .                      )**(-1.0/2.0)
!     .            )
!        sigmaX2 = SQRT( LOG(  1.0 + ( (sigmaY2**2.0) / (muY2**2.0) ) )
!     .                )
!
!        ! Y3 is a Lognormal.  It is converted to X3, which is a Gaussian.
!        ! X3 = ln Y3
!        muX3 = LOG(  muY3 * (  1.0 
!     .                       + ( (sigmaY3**2.0) / (muY3**2.0) ) 
!     .                      )**(-1.0/2.0)
!     .            )
!        sigmaX3 = SQRT( LOG(  1.0 + ( (sigmaY3**2.0) / (muY3**2.0) ) )
!     .                )
!
!        !!! Intra-Gaussian correlations.
!
!        ! corrY1iY2 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1iX2, which is a correlation between
!        ! two Gaussians.
!        corrX1iX2 = corrY1iY2 
!     .                   * SQRT( EXP(sigmaX2**2.0) - 1.0 ) / sigmaX2
!        ! corrY1iY3 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1iX3, which is a correlation between
!        ! two Gaussians.
!        corrX1iX3 = corrY1iY3
!     .                   * SQRT( EXP(sigmaX3**2.0) - 1.0 ) / sigmaX3
!        ! corrY2Y3 is a correlation between two Log-normals.  It must be 
!        ! converted to corrX2X3, which is a correlation between two Gaussians.
!        corrX2X3 = LOG( 1.0 + corrY2Y3
!     .                             * SQRT( EXP(sigmaX2**2.0) - 1.0 )
!     .                             * SQRT( EXP(sigmaX3**2.0) - 1.0 )
!     .                ) / ( sigmaX2 * sigmaX3 )
!
!        !----- Section #2 ------------------------------------------------------
!
!        !!!!!!!!!! Use the appropriate equation and solve.
!
!        ! Initialize logical to false.
!        use_eq_1 = .false.
!
!        ! Input to the gamma function.
!        gamma_fnc_input = alpha_exp + 1.0
!
!        ! Order of the parabolic cylinder function.
!        parab_cyl_fnc_ord = - ( alpha_exp + 1.0 )
!
!        ! If sigmaX1i = 0, then equation #1 cannot be used.
!        ! A tolerance value (a small number) is used instead of zero in order
!        ! to prevent numerical errors.
!        IF ( sigmaX1i < sigX1_tol ) THEN
!           use_eq_1 = .false.
!        ELSE
!           ! Input to the parabolic cylinder function.
!           parab_cyl_fnc_input =  (muX1i/sigmaX1i)
!     .                          + corrX1iX2*sigmaX2*beta_exp 
!     .                          + corrX1iX3*sigmaX3*gamma_exp
!
!           ! Test to see whether the value of sigmaX1i is sufficiently small
!           ! enough to cause the parabolic cylinder function to produce a
!           ! result to large to be represented numerically by the computer.
!           test = Dv_fnc( parab_cyl_fnc_ord, parab_cyl_fnc_input )
!
!           IF ( test >= 0.0d0 .AND. test < limit ) THEN
!              use_eq_1 = .true.
!           ELSE
!              use_eq_1 = .false.
!           ENDIF
!
!        ENDIF
!
!
!        IF ( use_eq_1 ) THEN
!
!           ! Case where sigmaX1i is not equal to 0 (and sigmaX1i is not
!           ! sufficiently small enough to cause the parabolic cylinder function
!           ! to produce a result greater than what the computer can represent 
!           ! numerically).
!
!           PDF_TRIVAR_2G_LN_LN =
!     .       (1.0/SQRT(2.0*pi)) * (-sigmaX1i)**alpha_exp 
!     .     * EXP( muX2*beta_exp + muX3*gamma_exp )
!     .     * EXP( (1.0/2.0) * (
!     .                            ( 1.0 - corrX1iX2**2.0 )
!     .                           *(sigmaX2**2.0)*(beta_exp**2.0)
!     .                         +  ( 1.0 - corrX1iX3**2.0 )
!     .                           *(sigmaX3**2.0)*(gamma_exp**2.0)
!     .                         + 2.0*( corrX2X3 - corrX1iX2*corrX1iX3 )
!     .                              *sigmaX2*beta_exp*sigmaX3*gamma_exp
!     .                        )
!     .          )
!     .     * EXP(   (1.0/4.0) * (
!     .                             (muX1i/sigmaX1i)
!     .                           + corrX1iX2*sigmaX2*beta_exp 
!     .                           + corrX1iX3*sigmaX3*gamma_exp
!     .                          )**2.0
!     .            - (muX1i/sigmaX1i) * (
!     .                                    (muX1i/sigmaX1i)
!     .                                  + corrX1iX2*sigmaX2*beta_exp 
!     .                                  + corrX1iX3*sigmaX3*gamma_exp
!     .                                 )
!     .            + (1.0/2.0) * ( (muX1i**2.0) / (sigmaX1i**2.0) )
!     .          )
!     .     * GAMMA( gamma_fnc_input )
!     .     * Dv_fnc( parab_cyl_fnc_ord, parab_cyl_fnc_input )
!
!        ELSEIF ( .not. use_eq_1 .AND. muX1i < 0.0 ) THEN
!
!           ! Case where sigmaX1i is equal to 0 (or sigmaX1i is sufficiently
!           ! small enough to cause the parabolic cylinder function to produce a
!           ! result greater than what the computer can represent numerically)
!           ! and muX1i is less than 0.
!
!           PDF_TRIVAR_2G_LN_LN =
!     .       (muX1i)**alpha_exp 
!     .     * EXP( muX2*beta_exp + muX3*gamma_exp )
!     .     * EXP( (1.0/2.0) * (
!     .                            ( 1.0 - corrX1iX2**2.0 )
!     .                           *(sigmaX2**2.0)*(beta_exp**2.0)
!     .                         +  ( 1.0 - corrX1iX3**2.0 )
!     .                           *(sigmaX3**2.0)*(gamma_exp**2.0)
!     .                         + 2.0*( corrX2X3 - corrX1iX2*corrX1iX3 )
!     .                              *sigmaX2*beta_exp*sigmaX3*gamma_exp
!     .                        )
!     .          )
!
!        ELSEIF ( .not. use_eq_1 .AND. muX1i >= 0.0 ) THEN
!           ! Case where sigmaX1i is equal to 0 (or sigmaX1i is sufficiently
!           ! small enough to cause the parabolic cylinder function to produce a
!           ! result greater than what the computer can represent numerically)
!           ! and muX1i is greater than or equal to 0.
!
!           PDF_TRIVAR_2G_LN_LN = 0.0
!
!        ENDIF
!
!        RETURN
!
!        END FUNCTION PDF_TRIVAR_2G_LN_LN
!
!!===============================================================================
!        !
!        ! FUNCTION PDF_BIVAR_2G_LN
!        !-----------------------------------------------------------------------
!        !
!        ! DESCRIPTION
!        !------------
!        !
!        ! BIVARIATE DOUBLE GAUSSIAN -- LOG-NORMAL PDF
!        !
!        ! This function helps solve for the expression:
!        !
!        ! overbar{ Y1^alpha Y2^beta }; where:
!        !
!        ! Y1 is distributed as a Double Gaussian with a domain from zero to
!        ! infinity.
!        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
!        !
!        ! Y1 and Y2 are then transformed into X1 and X2 respectively.
!        !
!        ! X1 is distributed as a Double Gaussian with a domain from zero to
!        ! infinity.  X1 = Y1
!        ! X2 is distributed as Normal (Gaussian) with a domain from negative
!        ! infinity to infinity.  X2 = LN( Y2 )
!        !
!        ! The resulting equation that needs to be solved for is:
!        !
!        ! overbar{ Y1^alpha Y2^beta } =
!        ! INT(0:inf) INT(-inf:inf)
!        !             X1^alpha exp( beta*X2 ) P(X1,X2) dX2 dX1
!        !
!        ! since X1 is a double Gaussian, P(X1,X2) =  ( a ) P_1(X1,X2) 
!        !                                          + (1-a) P_2(X1,X2)
!        !
!        ! where "a" is a constant and is the relative weight of each individual
!        ! Gaussian "plume".
!        !
!        ! P_1(X1,X2) and P_2(X1,X2) are simply the equation for a
!        ! Bivariate Gaussian Distribution.  The only difference between the
!        ! two is the dependence on each individual plume for the plume average
!        ! and plume variance for the X1 variable.  There is also a dependence
!        ! on each individual plume for the within-the-plume-correlation
!        ! (a.k.a. intra-Gaussian correlation) of X1 with any other variable.
!        !
!        ! The Probability Density Function for a Bivariate Normal Distribution:
!        !
!        ! P_i(X1,X2) =
!        !
!        ! 1 / { 2 * PI * sigmaX1i * sigmaX2 * SQRT[ 1 - corrX1iX2^2 ] }
!        ! * EXP{ -(1/2) * phi }
!        !
!        ! where,
!        !
!        ! phi =
!        !
!        ! 1 / [ 1 - corrX1iX2^2 ]
!        ! * {
!        !       [ 1 / sigmaX1i^2 ] * ( X1i - muX1i )^2
!        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
!        !     + [ 2*corrX1iX2 / (sigmaX1i*sigmaX2) ]
!        !       * ( X1i - muX1i ) * ( X2 - muX2 )
!        !   }
!        !
!        ! Definitions:
!        !
!        ! muX1i:  The plume average of the X1 (or Y1) variable.
!        ! muX2:   The average of the X2 (or ln Y2) variable.
!        !
!        ! sigmaX1i:  The plume standard deviation (square-root of variance)
!        !            of the X1 (or Y1) variable.
!        ! sigmaX2:   The standard deviation (square-root of variance) of the
!        !            X2 (or ln Y2) variable.
!        !
!        ! corrX1iX2:  Intra-Gaussian correlation between X1 and X2 (or between
!        !             Y1 and ln Y2).
!        !
!        ! After careful (and long) double-integration of the above equation,
!        ! one gets the result for each plume of the bivariate normal
!        ! -- log-normal equation.
!        !
!        ! The result for each plume is (Equation #1):
!        !
!        !   [ 1 / SQRT(2*PI) ] * ( sigmaX1i )^alpha
!        ! * EXP{    muX2*beta + (1/2)*sigmaX2^2*beta^2
!        !        - (1/4) * [ ( muX1i / sigmaX1i ) + corrX1iX2*sigmaX2*beta ]^2
!        !      }
!        ! * GAMMA_FNC { alpha + 1 }
!        ! * PARAB_CYL_FNC_[-(alpha+1)] { - [  ( muX1i / sigmaX1i )
!        !                                   + corrX1iX2*sigmaX2*beta ]
!        !                              }
!        !
!        ! The sum of the weighted results for each plume give us the answer for:
!        !
!        ! overbar{ Y1^alpha Y2^beta }
!        !
!        ! The purpose of this function is to give results for each plume.
!        !
!        ! Brian Griffin.  November 10, 2006.
!        !
!        !-----------------------------------------------------------------------
!        !
!        ! NOTE -- special case if sigmaX1i = 0
!        ! ------------------------------------
!        !
!        ! The equation above is only a valid answer if sigmaX1i is not equal
!        ! to 0.  The only variable that can be found in the denominator of any
!        ! term or factor of the equation is sigmaX1i.  This is also the only
!        ! case where the equation above is not a valid answer.
!        !
!        ! When sigmaX1i = 0, the Gaussian for X1i becomes a "spike" at muX1i
!        ! with no variance.  This is also known as a Delta function.  The
!        ! domain of X1i is from zero to infinity.  If sigmaX1i = 0 and 
!        ! muX1i <= 0, then the entire Gaussian lies outside the domain,
!        ! and the result is 0 for that individual plume.
!        !
!        ! However, if sigmaX1i = 0 and muX1i > 0, then the entire Gaussian
!        ! lies within the domain.  A close approximation for the parabolic
!        ! cylinder function is made.  This approximation is only valid in
!        ! cases where the input into the parabolic cylinder function is
!        ! extremely large in magnitude and very much greater in magnitude than
!        ! the order of the parabolic cylinder function -- as it is in this case
!        ! where sigmaX1i = 0 means that the magnitude of the input is infinite
!        ! and very much greater than the magnitude of the order of the function
!        ! which comes to | - (alpha + 1) | = alpha + 1.
!        !
!        ! The result for each plume is (Equation #2):
!        !
!        ! ( muX1i )^alpha * EXP{ muX2*beta + (1/2)*sigmaX2^2*beta^2 }
!        !
!        ! One also needs to note an issue with numerical error.  If sigmaX1i
!        ! is very small, but not quite 0, one would still be able to compute
!        ! the answer mathematically using the equation for sigmaX1i /= 0.
!        ! However, ( muX1i / sigmaX1i ) becomes very large.  If it becomes
!        ! large enough, the result of the parabolic cylinder function becomes
!        ! greater than the computer can represent.  The greatest value that
!        ! can be represented by a DOUBLE PRECISION variable is of the order
!        ! of 10^308.  Therefore, a sufficiently small sigmaX1i must be treated
!        ! the same as if sigmaX1i = 0.  A test is needed to determine if the
!        ! computer can represent the number or not.  Secondly, there is a case
!        ! where sigmaX1i^2 is in the denominator.  While sigmaX1i may not equal
!        ! 0, it may be small enough so that sigmaX1i^2 can not be represented
!        ! by the computer, thereby returning a value of 0 for sigmaX1i^2.
!        ! A tolerance value for sigmaX1i is used so that sigmaX1i^2 does not
!        ! return a value of 0 when sigmaX1i /= 0.
!        !
!        ! In Summary:
!        ! 
!        ! When sigmaX1i /= 0 (or is not sufficiently small):  Equation #1
!        ! When sigmaX1i = 0 (or is sufficiently small)
!        !                                     and muX1i > 0:  Equation #2
!        ! When sigmaX1i = 0 (or is sufficiently small)
!        !                                    and muX1i <= 0:  0
!        !
!        ! Brian Griffin.  November 10, 2006.
!        !
!        !-----------------------------------------------------------------------
!
!
!        FUNCTION PDF_BIVAR_2G_LN ( muY1i, muY2, sigmaY1i, sigmaY2,
!     .                             corrY1iY2, alpha_exp, beta_exp )
!
!        USE constants
!        USE polpak_gamma, ONLY: gamma
!
!        IMPLICIT NONE
!
!        ! Input variables.
!        REAL, INTENT(IN):: muY1i     ! Plume average of Y1
!        REAL, INTENT(IN):: muY2      ! Average of Y2
!        REAL, INTENT(IN):: sigmaY1i  ! Plume standard deviation of Y1
!        REAL, INTENT(IN):: sigmaY2   ! Standard deviation of Y2
!        REAL, INTENT(IN):: corrY1iY2 ! Intra-Gaussian correlation between Y1,Y2
!        REAL, INTENT(IN):: alpha_exp ! Exponent associated with Y1 variable.
!        REAL, INTENT(IN):: beta_exp  ! Exponent associated with Y2 variable.
!        
!        ! Output variable.
!        REAL:: PDF_BIVAR_2G_LN
!
!        ! Local variables.
!
!        ! "Y" variables Gaussianized and converted to "X" variables.
!        REAL:: muX1i     ! Plume average of X1
!        REAL:: muX2      ! Average of X2
!        REAL:: sigmaX1i  ! Plume standard deviation of X1
!        REAL:: sigmaX2   ! Standard deviation of X2
!        REAL:: corrX1iX2 ! Intra-Gaussian correlation between X1,X2
!
!        DOUBLE PRECISION:: gamma_fnc_input
!        DOUBLE PRECISION:: parab_cyl_fnc_ord
!        DOUBLE PRECISION:: parab_cyl_fnc_input
!        DOUBLE PRECISION:: test
!        LOGICAL:: use_eq_1
!!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308.0
!!       Found that above is not valid on most compilers -dschanen
!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308
!        REAL, PARAMETER:: sigX1_tol = 10.0**(-18)
!
!        !----- Section #1 ------------------------------------------------------
!
!        !!!!!!!!!! Convert all means, standard deviations, and correlations
!        !!!!!!!!!! to Gaussian terms.
!
!        ! Y1i is a truncated double Gaussian.
!        ! Therefore, X1i is still a truncated double Gaussian.
!        ! It does not change.
!        muX1i = muY1i
!        sigmaX1i = sigmaY1i
!
!        ! Y2 is a Lognormal.  It is converted to X2, which is a Gaussian.
!        ! X2 = ln Y2
!        muX2 = LOG(  muY2 * (  1.0 
!     .                       + ( (sigmaY2**2.0) / (muY2**2.0) ) 
!     .                      )**(-1.0/2.0)
!     .            )
!        sigmaX2 = SQRT( LOG(  1.0 + ( (sigmaY2**2.0) / (muY2**2.0) ) )
!     .                )
!
!        !!! Intra-Gaussian correlation.
!
!        ! corrY1iY2 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1iX2, which is a correlation between
!        ! two Gaussians.
!        corrX1iX2 = corrY1iY2 
!     .                   * SQRT( EXP(sigmaX2**2.0) - 1.0 ) / sigmaX2
!
!        !----- Section #2 ------------------------------------------------------
!
!        !!!!!!!!!! Use the appropriate equation and solve.
!
!        ! Initialize logical to false.
!        use_eq_1 = .false.
!
!        ! Input to the gamma function.
!        gamma_fnc_input = alpha_exp + 1.0
!
!        ! Order of the parabolic cylinder function.
!        parab_cyl_fnc_ord = - ( alpha_exp + 1.0 )
!
!        ! If sigmaX1i = 0, then equation #1 cannot be used.
!        ! A tolerance value (a small number) is used instead of zero in order
!        ! to prevent numerical errors.
!        IF ( sigmaX1i < sigX1_tol ) THEN
!           use_eq_1 = .false.
!        ELSE
!           ! Input to the parabolic cylinder function.
!           parab_cyl_fnc_input = - (  (muX1i/sigmaX1i)
!     .                              + corrX1iX2*sigmaX2*beta_exp )
!
!           ! Test to see whether the value of sigmaX1i is sufficiently small
!           ! enough to cause the parabolic cylinder function to produce a
!           ! result to large to be represented numerically by the computer.
!           test = Dv_fnc( parab_cyl_fnc_ord, parab_cyl_fnc_input )
!
!           IF ( test >= 0.0d0 .AND. test < limit ) THEN
!              use_eq_1 = .true.
!           ELSE
!              use_eq_1 = .false.
!           ENDIF
!
!        ENDIF
!
!
!        IF ( use_eq_1 ) THEN
!
!           ! Case where sigmaX1i is not equal to 0 (and sigmaX1i is not
!           ! sufficiently small enough to cause the parabolic cylinder function
!           ! to produce a result greater than what the computer can represent 
!           ! numerically).
!
!           PDF_BIVAR_2G_LN =
!     .       (1.0/SQRT(2.0*pi)) * (sigmaX1i)**alpha_exp 
!     .     * EXP(   muX2*beta_exp
!     .            + (1.0/2.0)*(sigmaX2**2.0)*(beta_exp**2.0)
!     .            - (1.0/4.0)*(  (muX1i/sigmaX1i)
!     .                         + corrX1iX2*sigmaX2*beta_exp )**2.0
!     .          )
!     .     * GAMMA( gamma_fnc_input )
!     .     * Dv_fnc( parab_cyl_fnc_ord, parab_cyl_fnc_input )
!
!        ELSEIF ( .not. use_eq_1 .AND. muX1i > 0.0 ) THEN
!
!           ! Case where sigmaX1i is equal to 0 (or sigmaX1i is sufficiently
!           ! small enough to cause the parabolic cylinder function to produce a
!           ! result greater than what the computer can represent numerically)
!           ! and muX1i is greater than 0.
!
!           PDF_BIVAR_2G_LN =
!     .       (muX1i)**alpha_exp 
!     .     * EXP(   muX2*beta_exp
!     .            + (1.0/2.0)*(sigmaX2**2.0)*(beta_exp**2.0) )
!
!        ELSEIF ( .not. use_eq_1 .AND. muX1i <=  0.0 ) THEN
!           ! Case where sigmaX1i is equal to 0 (or sigmaX1i is sufficiently
!           ! small enough to cause the parabolic cylinder function to produce a
!           ! result greater than what the computer can represent numerically)
!           ! and muX1i is less than or equal to 0.
!
!           PDF_BIVAR_2G_LN = 0.0
!
!        ENDIF
!
!        RETURN
!
!        END FUNCTION PDF_BIVAR_2G_LN
!
!!===============================================================================
!        !
!        ! FUNCTION PDF_BIVAR_LN_LN
!        !-----------------------------------------------------------------------
!        !
!        ! DESCRIPTION
!        !------------
!        !
!        ! BIVARIATE LOG-NORMAL -- LOG-NORMAL PDF
!        !
!        ! This function helps solve for the expression:
!        !
!        ! overbar{ Y1^alpha Y2^beta }; where:
!        !
!        ! Y1 is distributed as Log-normal with a domain from zero to infinity.
!        ! Y2 is distributed as Log-normal with a domain from zero to infinity.
!        !
!        ! Y1 and Y2 are then transformed into X1 and X2 respectively.
!        !
!        ! X1 is distributed as Normal (Gaussian) with a domain from negative
!        ! infinity to infinity.  X1 = LN( Y1 )
!        ! X2 is distributed as Normal (Gaussian) with a domain from negative
!        ! infinity to infinity.  X2 = LN( Y2 )
!        !
!        ! The resulting equation that needs to be solved for is:
!        !
!        ! overbar{ Y1^alpha Y2^beta } =
!        ! INT(-inf:inf) INT(-inf:inf)
!        !             exp( alpha*X1 + beta*X2 ) P(X1,X2) dX2 dX1
!        !
!        ! P(X1,X2) is simply the equation for a Bivariate Gaussian Distribution.
!        !
!        ! The Probability Density Function for a Bivariate Normal Distribution:
!        !
!        ! P(X1,X2) = 
!        !
!        ! 1 / { 2 * PI * sigmaX1 * sigmaX2 * SQRT[ 1 - corrX1X2^2 ] }
!        ! * EXP{ -(1/2) * phi }
!        !
!        ! where,
!        !
!        ! phi =
!        !
!        ! 1 / [ 1 - corrX1X2^2 ]
!        ! * {
!        !       [ 1 / sigmaX1^2 ] * ( X1 - muX1 )^2
!        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
!        !     + [ 2*corrX1X2 / (sigmaX1*sigmaX2) ]
!        !       * ( X1 - muX1 ) * ( X2 - muX2 )
!        !   }
!        !
!        ! Definitions:
!        !
!        ! muX1:  The average of the X1 (or ln Y1) variable.
!        ! muX2:  The average of the X2 (or ln Y2) variable.
!        !
!        ! sigmaX1:  The standard deviation (square-root of variance) of the
!        !           X1 (or ln Y1) variable.
!        ! sigmaX2:  The standard deviation (square-root of variance) of the
!        !           X2 (or ln Y2) variable.
!        !
!        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
!        !            ln Y1 and ln Y2).
!        !
!        ! After careful double-integration of the above equation, one gets the
!        ! result for each plume of the bivariate log-normal -- log-normal
!        ! equation.
!        !
!        ! The result is:
!        !
!        ! * EXP{    muX1*alpha + muX2*beta
!        !        + (1/2)*sigmaX1^2*alpha^2 + (1/2)*sigmaX2^2*beta^2
!        !        + corrX1X2*sigmaX1*alpha*sigmaX2*beta
!        !      }
!        !
!        ! The above equation gives us the answer for:
!        !
!        ! overbar{ Y1^alpha Y2^beta }
!        !
!        ! Brian Griffin.  November 14, 2006.
!        !
!        !-----------------------------------------------------------------------
!
!
!        FUNCTION PDF_BIVAR_LN_LN ( muY1, muY2, sigmaY1, sigmaY2,
!     .                             corrY1Y2, alpha_exp, beta_exp )
!
!        USE constants
!
!        IMPLICIT NONE
!
!        ! Input variables.
!        REAL, INTENT(IN):: muY1      ! Average of Y1
!        REAL, INTENT(IN):: muY2      ! Average of Y2
!        REAL, INTENT(IN):: sigmaY1   ! Standard deviation of Y1
!        REAL, INTENT(IN):: sigmaY2   ! Standard deviation of Y2
!        REAL, INTENT(IN):: corrY1Y2  ! Intra-Gaussian correlation between Y1,Y2
!        REAL, INTENT(IN):: alpha_exp ! Exponent associated with Y1 variable.
!        REAL, INTENT(IN):: beta_exp  ! Exponent associated with Y2 variable.
!        
!        ! Output variable.
!        REAL:: PDF_BIVAR_LN_LN
!
!        ! Local variables.
!
!        ! "Y" variables Gaussianized and converted to "X" variables.
!        REAL:: muX1      ! Average of X1
!        REAL:: muX2      ! Average of X2
!        REAL:: sigmaX1   ! Standard deviation of X1
!        REAL:: sigmaX2   ! Standard deviation of X2
!        REAL:: corrX1X2  ! Intra-Gaussian correlation between X1,X2
!
!        !----- Section #1 ------------------------------------------------------
!
!        !!!!!!!!!! Convert all means, standard deviations, and correlations
!        !!!!!!!!!! to Gaussian terms.
!
!        ! Y1 is a Lognormal.  It is converted to X1, which is a Gaussian.
!        ! X1 = ln Y1
!        muX1 = LOG(  muY1 * (  1.0 
!     .                       + ( (sigmaY1**2.0) / (muY1**2.0) ) 
!     .                      )**(-1.0/2.0)
!     .            )
!        sigmaX1 = SQRT( LOG(  1.0 + ( (sigmaY1**2.0) / (muY1**2.0) ) )
!     .                )
!
!        ! Y2 is a Lognormal.  It is converted to X2, which is a Gaussian.
!        ! X2 = ln Y2
!        muX2 = LOG(  muY2 * (  1.0 
!     .                       + ( (sigmaY2**2.0) / (muY2**2.0) ) 
!     .                      )**(-1.0/2.0)
!     .            )
!        sigmaX2 = SQRT( LOG(  1.0 + ( (sigmaY2**2.0) / (muY2**2.0) ) )
!     .                )
!
!        !!! Intra-Gaussian correlation.
!        ! corrY1Y2 is a correlation between two Log-normals.  It must be 
!        ! converted to corrX1X2, which is a correlation between two Gaussians.
!        corrX1X2 = LOG( 1.0 + corrY1Y2
!     .                             * SQRT( EXP(sigmaX1**2.0) - 1.0 )
!     .                             * SQRT( EXP(sigmaX2**2.0) - 1.0 )
!     .                ) / ( sigmaX1 * sigmaX2 )
!
!        !----- Section #2 ------------------------------------------------------
!
!        !!!!!!!!!! Solve.
!
!        PDF_BIVAR_LN_LN =
!     .     EXP(   muX1*alpha_exp + muX2*beta_exp
!     .          + (1.0/2.0)*(sigmaX1**2.0)*(alpha_exp**2.0)
!     .          + (1.0/2.0)*(sigmaX2**2.0)*(beta_exp**2.0)
!     .          + corrX1X2*sigmaX1*alpha_exp*sigmaX2*beta_exp
!     .        )
!
!        RETURN
!
!        END FUNCTION PDF_BIVAR_LN_LN
!
!!===============================================================================

!===============================================================================
        double precision function Dv_fnc( order, argument )

!       Description:
!       Compute the parabolic cylinder function in terms of Dv
!       using an Algorithm from ACM TOMS.  Replaces the more expensive
!       D_fnc function used previously.

!       References:
!       Algorithm 850, collected algorithms from ACM.
!       ACM Transactions on Mathematical Software,
!       Vol. 32, No. 1, March 2006 pp. 102--112
!===============================================================================
        use parabolic, only: 
     .      gamma, ! Variable(s) 
     .      parab

        use constants, only: 
     .      pi_dp ! Variable(s)

        implicit none

        ! External
        intrinsic :: sin

        ! Parameter constants
        integer, parameter ::
     .  scaling = 0 ! 0 = Unscaled functions, 1 = scaled functions

        double precision, parameter :: limit = 10.0d0**308

        ! Input Variables
        double precision, intent(in) ::
     .  order,   ! Order 'a' of Dv(a,x)         [-]
     .  argument ! Argument 'x' of Dv(a,x)      [-]

        ! Local Variables
        double precision, dimension(2) ::
     .  uaxx,  ! U(a,x), U'(a,x)                [-]
     .  vaxx   ! V(a,x), V'(a,x)                [-]
        ! Where a is the order and x is the argument

        integer :: ierr ! Error condition

        if ( argument <= 0.0d0 ) then
          call parab( -order-0.5, -argument, scaling, uaxx, vaxx, ierr )
          Dv_fnc = vaxx(1) / ( (1.0d0/pi_dp) * gamma( -order ) )
     .           - sin( pi_dp * ( -order ) ) * uaxx(1)
        else
          call parab( -order-0.5, argument, scaling, uaxx, vaxx, ierr )
          Dv_fnc = uaxx(1)
        end if

        ! Handle the overflow condition
        if ( ierr /= 0 ) then 
          Dv_fnc = limit
        end if

        return
        end function Dv_fnc

       END MODULE rain_equations
