!----------------------------------------------------------------------
! $Id: rain.F,v 1.21 2006-12-20 18:53:02 griffinb Exp $
!----------------------------------------------------------------------

      MODULE rain_equations

      IMPLICIT NONE

      ! Statistical rain parameters.
!      ! Old parameterization.
!      REAL, PARAMETER:: rrp2_rrm2 = 1.2
!      REAL, PARAMETER:: Nrp2_Nrm2 = 1.0
!      REAL, PARAMETER:: Ncp2_Ncm2 = 0.07
!      REAL, PARAMETER:: corr_rrNr_LL =  0.85
!      REAL, PARAMETER:: corr_rcrr_NL =  0.10
!      REAL, PARAMETER:: corr_rcNr_NL = -0.15
!      REAL, PARAMETER:: corr_rcNc_NL =  0.45
!      ! Parameters for case with cloud water sedimentation.
!      REAL, PARAMETER:: rrp2_rrm2 = 1.7123
!      REAL, PARAMETER:: Nrp2_Nrm2 = 1.0618
!      REAL, PARAMETER:: Ncp2_Ncm2 = 0.2964
!      REAL, PARAMETER:: corr_rrNr_LL = 0.872580
!      REAL, PARAMETER:: corr_rcrr_NL = 0.393905
!      REAL, PARAMETER:: corr_rcNr_NL = 0.402948
!      REAL, PARAMETER:: corr_rcNc_NL = 0.670482
!      ! Parameters for case without cloud water sedimentation.
!      REAL, PARAMETER:: rrp2_rrm2 = 1.4584
!      REAL, PARAMETER:: Nrp2_Nrm2 = 1.0285
!      REAL, PARAMETER:: Ncp2_Ncm2 = 0.1034
!      REAL, PARAMETER:: corr_rrNr_LL = 0.861597
!      REAL, PARAMETER:: corr_rcrr_NL = 0.252875
!      REAL, PARAMETER:: corr_rcNr_NL = 0.196055
!      REAL, PARAMETER:: corr_rcNc_NL = 0.616638
      ! Rounded average of the above two sets (to nearest 0.05).
      REAL, PARAMETER:: rrp2_rrm2 = 1.60
      REAL, PARAMETER:: Nrp2_Nrm2 = 1.05
      REAL, PARAMETER:: Ncp2_Ncm2 = 0.20
      REAL, PARAMETER:: corr_rrNr_LL = 0.85
      REAL, PARAMETER:: corr_rcrr_NL = 0.30
      REAL, PARAMETER:: corr_rcNr_NL = 0.30
      REAL, PARAMETER:: corr_rcNc_NL = 0.65
      !REAL, PARAMETER:: C_evap = 0.86    ! Khairoutdinov and Kogan (2000) 
      !                                   ! ratio of drizzle drop mean 
      !                                   ! geometric radius to drizzle drop
      !                                   ! mean volume radius. Khairoutdinov
      !                                   ! and Kogan (2000); p. 233.
      !REAL, PARAMETER:: C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
      REAL, PARAMETER:: C_evap = 0.55    ! KK 2000, Marshall-Palmer (1948) value.
      ! Vince Larson set r_0=28mum to agree with COAMPS-LES formula. 15 April 2005
      REAL, PARAMETER:: r_0 = 25.0e-6   ! Assumed radius of all new drops; m.
                                        ! Value specified in KK (2000); p. 235.
      !REAL, PARAMETER:: r_0 = 28.0e-6   ! Assumed radius of all new drops; m.
                                         ! Value that COAMPS LES has in it.
      !REAL, PARAMETER:: r_0 = 30.0e-6   ! Assumed radius of all new drops; m.
                                         ! Khairoutdinov said it was okay!
      ! End Vince Larson's change.

      PUBLIC :: rain, kk_microphys
      PRIVATE :: mean_volume_radius, cond_evap_rrm, cond_evap_Nrm
      PRIVATE :: autoconv_rrm, autoconv_Nrm, accretion_rrm
      PRIVATE :: G_T_p
!      PRIVATE :: PDF_TRIVAR_2G_LN_LN, PDF_BIVAR_2G_LN, PDF_BIVAR_LN_LN
      PRIVATE :: Dv_fnc

      PRIVATE ! Default

      CONTAINS
!----------------------------------------------------------------------
! This subroutine advances drizzle (rrm and Nrm) one timestep.
! It uses implicit discretization and formulas by 
!     Khairoutdinov and Kogan (2000)
!----------------------------------------------------------------------
        SUBROUTINE rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .                   rtm, rcm, Ncm, rrm, Nrm, cond_rrm,
     .                   auto_rrm, accr_rrm, Vrr, mean_vol_rad,
     .                   AKm_est, AKm, isValid )

        USE model_flags
        USE grid_class
        USE constants
        USE lapack_wrap, ONLY: tridag_solve
#ifdef STATS
        USE statistics
#endif /*STATS*/

        implicit none

        ! External functions

        real rsat
        external rsat

        ! Input
        REAL, INTENT(IN):: dt
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: thlm   ! theta-l
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: p      ! pressure
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner  ! exner function
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot   ! density
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: wmt    ! vertical velocity
        REAL, DIMENSION(1:gr%nnzp, 1:26), INTENT(IN):: pdf_parms
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rtm    ! total water mixing ratio
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm    ! cloud water mixing ratio
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: Ncm    ! cloud droplet number conc.
                                                        !   number/m^3

        ! Input/Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: rrm ! rain water mixing ratio
        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: Nrm ! rain droplet number conc.

        ! Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: cond_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: auto_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: accr_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Vrr
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: mean_vol_rad

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: rvm          ! water vapor mixing ratio; kg/kg
        REAL, DIMENSION(1:gr%nnzp):: Temp         ! Temperature; K
        REAL, DIMENSION(1:gr%nnzp):: es           ! Saturation vapor pressure; Pa
        REAL, DIMENSION(1:gr%nnzp):: Supsat       ! Supersaturation
        REAL, DIMENSION(1:gr%nnzp):: VNr          ! Number conc. sedimentation 
                                                  !   velocity; m/s 
       
        REAL, DIMENSION(1:gr%nnzp):: aa, bb, cc, dd  ! Tridiagonals
        REAL, DIMENSION(1:gr%nnzp):: oo, pp, qq, rr

        REAL, DIMENSION(1:gr%nnzp):: rs           ! Saturation mixing ratio; kg/kg

        REAL, DIMENSION(1:gr%nnzp):: cond_Nrm
        REAL, DIMENSION(1:gr%nnzp):: auto_Nrm


        REAL:: e             ! Vapor pressure; Pa
        REAL:: thm           ! Potential temperature; K


        REAL:: a, thl1, thl2, rc1, rc2, s1, s2, ss1, ss2

        INTEGER:: i, k
        INTEGER:: km1, kp1

        REAL :: Beta_T


        ! Function Calls
!        REAL:: cond_evap_rrm, autoconv_rrm, accretion_rrm
!        REAL:: cond_evap_Nrm, autoconv_Nrm

! Whether this is valid
        LOGICAL, INTENT(INOUT)                    :: isValid

! Latin hypercube estimate of Kessler autoconversion   Vince Larson 22 May 2005
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN) :: AKm_est
! Analytic calculation of Kessler autoconversion       Vince Larson 22 May 2005
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN) :: AKm



! DESCRIPTION
!
! This subroutine computes the new rain water mixing ratio, rr.
! It computes it according to the equation:
!
! drr/dt = -w*(drr/dz) + Vrr*(drr/dz) + nu_r*(d2rr/dz2)
!                      + (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! where Vrr is the mean sedimentation velocity of the rain drops.  It
! is determined by Vrr = 0.012*rvr - 0.2, with Vrr in m/s and rvr, the
! rain drop mean volume radius, in um.  In turn the rain drop mean
! volume radius is given by the following equation:
!
! rvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^(1/3) * Nr^(-1/3);
!
! where Nr is the rain drop concentration (num/m^3).  rvr is given in
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
!      3*Cevap*G(T,p)*[(4*PI*rho_lw)/(3*rho)]^(2/3) * rr^(1/3) * Nr^(2/3) * S
!
! G(T,p) = 1/(Fk + Fd); where
!
! Fk = [Lv/(Rv*T) - 1]*(Lv*rho_lw)/(Ka*T)
! Fd = (rho_lw*Rv*T)/(Dv*es(T))
!
! Ka is the Coefficient of Thermal Conductivity of Air, and Dv is the
! Coefficient of Diffusion of Water Vapor in Air.
!
! (drr/dt)auto = 1350 * rc^2.47 * Nc ^ -1.79;
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
! (dNr/dt)cond = (Nr/rr) * (drr/dt)cond
!
! (dNr/dt)auto = (drr/dt)auto / ( [(4*PI*rho_lw)/(3*rho)] * r_o^3 )
!
! where r_o is the assumed radius of all newly-formed rain droplets (m).
! It is given a value of 25 um.
!
! This information is plugged into the tridiagonal matrix in order to
! solve for Nr at the next timestep.
!
! Description and notes written by Brian Griffin.

#ifdef STATS

        if (lstats_samp) then

         if ( irrm_bt > 0 ) then
           zt%x(:,irrm_bt) = zt%x(:,irrm_bt) - rrm/dt
         endif

         if ( irrm_bt > 0 ) then
           zt%x(:,iNrm_bt) = zt%x(:,iNrm_bt) - Nrm/dt
         endif

        endif

#endif /*STATS*/

        ! Find the drop mean volume radius.  It is calculated using
        ! the rain water mixing ratio, the rain droplet concentration,
        ! and the air density.  These values are taken from the previous
        ! timestep.  It is located on thermodynamic levels.
        DO k = 1, gr%nnzp, 1

           mean_vol_rad(k) = mean_volume_radius( rrm(k), Nrm(k),
     .                                           rhot(k)         )

        ENDDO


        ! The sedimentation velocity is found from the drop mean volume radius.
        ! It is located on the momentum levels.  This is due to the fact that
        ! momentum levels are where other vertical velocities are stored, such
        ! as the vertical component of wind velocity (w).

        ! Khairoutdinov and Kogan Sedimentation Velocity Calculation.
        ! sedimentation velocity of rain drops (in m/s).
        ! mean volume radius converted to um.
        ! Note 1: positive sedimentation velocity means downwards.
        ! Note 2: Vrr(m/s) = 0.012*rvr(um) - 0.2
        ! If rvr(um) is too small, the equation will result in an
        ! upwards (negative) sedimentation velocity.  The limiter
        ! is put in there for that reason.

        DO k = 1, gr%nnzp, 1

           kp1 = MIN(k+1,gr%nnzp)

           ! rrm sedimentation velocity.
           Vrr(k) = 0.012 
     .             * ( 1000000.0 * (1.0/2.0) 
     .             * ( mean_vol_rad(kp1) + mean_vol_rad(k) ) )
     .             - 0.2

           Vrr(k) = MAX( Vrr(k), 0.0 )

           ! Nrm sedimentation velocity.
           VNr(k) = 0.007 
     .             * ( 1000000.0 * (1.0/2.0) 
     .             * ( mean_vol_rad(kp1) + mean_vol_rad(k) ) )
     .             - 0.1

           VNr(k) = MAX( VNr(k), 0.0 )

        ENDDO


        ! Find values for other variables.
        DO k = 2, gr%nnzp, 1

           ! Find the important conditions

           rvm(k) = rtm(k) - rcm(k)   ! rvm is water vapor mixing ratio.

           ! Find current vapor pressure.
           e = (p(k)*rvm(k))/(ep + rvm(k))

           ! Get theta from theta-l
           thm = thlm(k) + (Lv/(Cp*exner(k)))*rcm(k)

           ! Saturation mixing ratio
           rs(k) = rsat(p(k),thm*exner(k))
           ! Find temperature and saturation vapor pressure.
           Temp(k) = thm/((p0/p(k))**kappa)
           es(k) = (p(k)*rs(k))/(ep + rs(k))

        ENDDO

        ! Set the boundary conditions
        rvm(1)          = 0.
        rs(1)           = 0.
        Temp(1)         = 0.
        es(1)           = 0.
        Supsat(1)       = 0.
        Supsat(gr%nnzp) = 0.

        ! Boundary Conditions.
        ! These have been chosen so that the column totals in the LHS matrix
        ! of the diffusional term are equal to zero.  This leads to 
        ! "conservative differencing" of the diffusion term for equally spaced grids.     
        ! The ground is impermeable to diffusion of rrm, Nrm.   Brian.

        ! rain water mixing ratio (rrm)
        ! lower boundary
        aa(1) = 0.0                                                    ! aa(1) term

        bb(1) = 
     .         + (1.0/dt)                                              ! bb(1) time tendency
     .         - ( - nu_r * gr%dzt(2) * gr%dzm(1) )                    ! minus aa(2) diffusion

        cc(1) = 
     .         + wmt(1) * (1.0/2.0) * gr%dzt(1)                        ! cc(1) advection
     .         - (1.0/4.0) * ( Vrr(1) + Vrr(1) ) * gr%dzt(1)           ! cc(1) sedimentation
     .         - ( + nu_r * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) ) )    ! minus bb(2) diffusion
     .         - ( - nu_r * gr%dzt(3) * gr%dzm(2) )                    ! minus aa(3) diffusion

        ! upper boundary
        aa(gr%nnzp) =
     .         - (   nu_r*gr%dzt(gr%nnzp-1)
     .              *(gr%dzm(gr%nnzp-1) + gr%dzm(gr%nnzp-2)) )  ! bb(gr%nnzp-1) diffusion
     .         - ( - nu_r*gr%dzt(gr%nnzp-2)*gr%dzm(gr%nnzp-2) ) ! cc(gr%nnzp-2) diffusion
     .         - wmt(gr%nnzp) 
     .          / (1.0/gr%dzm(gr%nnzp) + 1.0/gr%dzm(gr%nnzp-1)) ! aa(gr%nnzp) advection


        bb(gr%nnzp) =
     .         - ( - nu_r*gr%dzt(gr%nnzp-1)*gr%dzm(gr%nnzp-1) ) ! cc(gr%nnzp-1) diffusion
     .         + (1.0/dt)                                       ! bb(gr%nnzp) time tendency
                                                                ! bb(gr%nnzp) sedimentation = 0.

        cc(gr%nnzp) = 0.0                                       ! cc(gr%nnzp) term

        ! Implicit contributions to rrm.
#ifdef STATS

        if (lstats_samp) then

         if ( irrm_ma > 0 ) then
           ztscr01(1) = 0.0
           ztscr02(1) = - wmt(1) * (1.0/2.0) * gr%dzt(1)
         endif

         if ( irrm_sd > 0 ) then
           ztscr03(1) = 0.0
           ztscr04(1) = + (1.0/4.0) * ( Vrr(1) + Vrr(1) ) * gr%dzt(1)
         endif

         if ( irrm_dff > 0 ) then
           ztscr05(1) = 0.0
           ztscr06(1) = - nu_r * gr%dzt(2) * gr%dzm(1)
           ztscr07(1) = + nu_r * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) )
     .                  - nu_r * gr%dzt(3) * gr%dzm(2)
         endif

        endif

#endif /*STATS*/

        ! rain drop concentration (Nrm)
        ! lower boundary
        oo(1) = 0.0                                                    ! oo(1) term

        pp(1) = 
     .         + (1.0/dt)                                              ! pp(1) time tendency
     .         - ( - nu_r * gr%dzt(2) * gr%dzm(1) )                    ! minus oo(2) diffusion

        qq(1) = 
     .         + wmt(1) * (1.0/2.0) * gr%dzt(1)                        ! qq(1) advection
     .         - (1.0/4.0) * ( VNr(1) + VNr(1) ) * gr%dzt(1)           ! qq(1) sedimentation
     .         - ( + nu_r * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) ) )    ! minus pp(2) diffusion
     .         - ( - nu_r * gr%dzt(3) * gr%dzm(2) )                    ! minus oo(3) diffusion

        ! upper boundary
        oo(gr%nnzp) =
     .         - (   nu_r*gr%dzt(gr%nnzp-1)
     .              *(gr%dzm(gr%nnzp-1) + gr%dzm(gr%nnzp-2)) )  ! pp(gr%nnzp-1) diffusion
     .         - ( - nu_r*gr%dzt(gr%nnzp-2)*gr%dzm(gr%nnzp-2) ) ! qq(gr%nnzp-2) diffusion
     .         - wmt(gr%nnzp) 
     .          / (1.0/gr%dzm(gr%nnzp) + 1.0/gr%dzm(gr%nnzp-1)) ! oo(gr%nnzp) advection


        pp(gr%nnzp) =
     .         - ( - nu_r*gr%dzt(gr%nnzp-1)*gr%dzm(gr%nnzp-1) ) ! qq(gr%nnzp-1) diffusion
     .         + (1.0/dt)                                       ! pp(gr%nnzp) time tendency
                                                                ! pp(gr%nnzp) sedimentation = 0.

        qq(gr%nnzp) = 0.0                                       ! qq(gr%nnzp) term

        ! Implicit contributions to Nrm.
#ifdef STATS

        if (lstats_samp) then

         if ( iNrm_ma > 0 ) then
           ztscr08(1) = 0.0
           ztscr09(1) = - wmt(1) * (1.0/2.0) * gr%dzt(1)
         endif

         if ( iNrm_sd > 0 ) then
           ztscr10(1) = 0.0
           ztscr11(1) = + (1.0/4.0) * ( VNr(1) + VNr(1) ) * gr%dzt(1)
         endif

         if ( iNrm_dff > 0 ) then
           ztscr12(1) = 0.0
           ztscr13(1) = - nu_r * gr%dzt(2) * gr%dzm(1)
           ztscr14(1) = + nu_r * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) )
     .                  - nu_r * gr%dzt(3) * gr%dzm(2)
         endif

        endif

#endif /*STATS*/



        ! Main Loop
        DO k = 2, gr%nnzp-1, 1

           km1 = MAX(k-1,1)
           kp1 = MIN(k+1,gr%nnzp)

           ! LHS terms for rain water mixing ratio, rrm.

           ! Find subdiagonal aa
           aa(k) = 
     .            - wmt(k) * (1.0/2.0) * gr%dzt(k)                  ! advection
     .            + (1.0/4.0) * ( Vrr(k) + Vrr(km1) ) * gr%dzt(k)   ! sedimentation
     .            - nu_r * gr%dzt(k) * gr%dzm(km1)                  ! diffusion
     
           ! Find main diagonal bb
           bb(k) = 
     .            + (1.0/dt)                                         ! time tendency
     .            + nu_r * gr%dzt(k) * ( gr%dzm(k) + gr%dzm(km1) )   ! diffusion

           ! Find superdiagonal cc
           cc(k) = 
     .            + wmt(k) * (1.0/2.0) * gr%dzt(k)                  ! advection
     .            - (1.0/4.0) * ( Vrr(k) + Vrr(km1) ) * gr%dzt(k)   ! sedimentation
     .            - nu_r * gr%dzt(k) * gr%dzm(k)                    ! diffusion


        ! Implicit contributions to rrm.
#ifdef STATS

        if (lstats_samp) then

         if ( irrm_ma > 0 ) then
           ztscr01(k) = + wmt(k) * (1.0/2.0) * gr%dzt(k)
           ztscr02(k) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
         endif

         if ( irrm_sd > 0 ) then
           ztscr03(k) = - (1.0/4.0) * ( Vrr(k) + Vrr(km1) ) * gr%dzt(k)
           ztscr04(k) = + (1.0/4.0) * ( Vrr(k) + Vrr(km1) ) * gr%dzt(k)
         endif

         if ( irrm_dff > 0 ) then
           ztscr05(k) = + nu_r * gr%dzt(k) * gr%dzm(km1)
           ztscr06(k) = - nu_r * gr%dzt(k) * ( gr%dzm(k) + gr%dzm(km1) )
           ztscr07(k) = + nu_r * gr%dzt(k) * gr%dzm(k)
         endif

        endif

#endif /*STATS*/


           ! LHS terms for rain droplet number concentration, Nrm.

           ! Find subdiagonal oo
           oo(k) = 
     .            - wmt(k) * (1.0/2.0) * gr%dzt(k)                  ! advection
     .            + (1.0/4.0) * ( VNr(k) + VNr(km1) ) * gr%dzt(k)   ! sedimentation
     .            - nu_r * gr%dzt(k) * gr%dzm(km1)                  ! diffusion

           ! Find main diagonal pp
           pp(k) = 
     .            + (1.0/dt)                                         ! time tendency
     .            + nu_r * gr%dzt(k) * ( gr%dzm(k) + gr%dzm(km1) )   ! diffusion

           ! Find superdiagonal qq
           qq(k) = 
     .            + wmt(k) * (1.0/2.0) * gr%dzt(k)                  ! advection
     .            - (1.0/4.0) * ( VNr(k) + VNr(km1) ) * gr%dzt(k)   ! sedimentation
     .            - nu_r * gr%dzt(k) * gr%dzm(k)                    ! diffusion


        ! Implicit contributions to Nrm.
#ifdef STATS

        if (lstats_samp) then

         if ( iNrm_ma > 0 ) then
           ztscr08(k) = + wmt(k) * (1.0/2.0) * gr%dzt(k)
           ztscr09(k) = - wmt(k) * (1.0/2.0) * gr%dzt(k)
         endif

         if ( iNrm_sd > 0 ) then
           ztscr10(k) = - (1.0/4.0) * ( VNr(k) + VNr(km1) ) * gr%dzt(k)
           ztscr11(k) = + (1.0/4.0) * ( VNr(k) + VNr(km1) ) * gr%dzt(k)
         endif

         if ( iNrm_dff > 0 ) then
           ztscr12(k) = + nu_r * gr%dzt(k) * gr%dzm(km1)
           ztscr13(k) = - nu_r * gr%dzt(k) * ( gr%dzm(k) + gr%dzm(km1) )
           ztscr14(k) = + nu_r * gr%dzt(k) * gr%dzm(k)
         endif

        endif

#endif /*STATS*/


           ! PDF parameters from pdf_closure_new for use in calculation
           ! of condensation/evaporation, autoconversion, accretion.

           a    = pdf_parms(k,13)
           thl1 = pdf_parms(k,9)
           thl2 = pdf_parms(k,10)
           rc1  = pdf_parms(k,14)
           rc2  = pdf_parms(k,15)
           s1   = pdf_parms(k,20)
           s2   = pdf_parms(k,21)
           ss1  = pdf_parms(k,22)
           ss2  = pdf_parms(k,23)

           ! Compute supersaturation via s1, s2. 
           !     Larson et al 2002, JAS, Vol 59, p 3534.
           ! This allows a more direct comparison of local, nonlocal formulas.
           Beta_T = (Rd/Rv) * ( Lv/(Rd*Temp(k)) ) * ( Lv/(Cp*Temp(k)) )

           Supsat(k) = (a*s1+(1-a)*s2)*((1.0 + Beta_T*rs(k))/rs(k))

           ! Now find the elements that make up the right-hand side of the
           ! equation, dd, for rain water mixing ratio, rrm.

           cond_rrm(k) = cond_evap_rrm( rrm(k), Nrm(k),
     .                     s1, ss1, s2, ss2, thl1, thl2, rc1, rc2, a,
     .                     p(k), rhot(k), Temp(k), Supsat(k)  )

! Vince Larson added option to call LH sampled Kessler autoconversion.
! 22 May 2005
!           auto_rrm(k) = autoconv_rrm( rcm(k), Ncm(k), rhot(k),
!     .                        a, s1, s2, ss1, ss2, rc1, rc2 )
           if ( LH_on ) then

!              auto_rrm(k) = AKm_est(k)
              auto_rrm(k) = AKm(k)

           else

              auto_rrm(k) = autoconv_rrm( rcm(k), Ncm(k),
     .                           s1, ss1, s2, ss2, a, rhot(k) )

           endif
! End Vince Larson's addition

            accr_rrm(k) = accretion_rrm( rcm(k), rrm(k),
     .                                   s1, ss1, s2, ss2, a )

        ! Explicit contributions to rrm.
#ifdef STATS

        if (lstats_samp) then

         if ( irrm_cond > 0 ) then
           zt%x(k,irrm_cond) = zt%x(k,irrm_cond) + cond_rrm(k)
           zt%n(k,irrm_cond) = zt%n(k,irrm_cond) + 1
         endif

         if ( irrm_auto > 0 ) then
           zt%x(k,irrm_auto) = zt%x(k,irrm_auto) + auto_rrm(k)
           zt%n(k,irrm_auto) = zt%n(k,irrm_auto) + 1
         endif

         if ( irrm_accr > 0 ) then
           zt%x(k,irrm_accr) = zt%x(k,irrm_accr) + accr_rrm(k)
           zt%n(k,irrm_accr) = zt%n(k,irrm_accr) + 1
         endif

        endif

#endif /*STATS*/


           ! Now find the elements that make up the right-hand side of the
           ! equation, rr, for rain droplet number concentration, Nrm.

           cond_Nrm(k) = cond_evap_Nrm( cond_rrm(k), Nrm(k), rrm(k) )

           auto_Nrm(k) = autoconv_Nrm( auto_rrm(k), rhot(k) )


        ! Explicit contributions to Nrm.
#ifdef STATS

        if (lstats_samp) then

         if ( iNrm_cond > 0 ) then
           zt%x(k,iNrm_cond) = zt%x(k,iNrm_cond) + cond_Nrm(k)
           zt%n(k,iNrm_cond) = zt%n(k,iNrm_cond) + 1
         endif

         if ( iNrm_auto > 0 ) then
           zt%x(k,iNrm_auto) = zt%x(k,iNrm_auto) + auto_Nrm(k)
           zt%n(k,iNrm_auto) = zt%n(k,iNrm_auto) + 1
         endif

        endif

#endif /*STATS*/


        ENDDO ! k=2..gr%nnzp-1



        ! Enter the entire dd array, which is the RHS of the equation.
        dd(2:gr%nnzp-1) = (1.0/dt) * rrm(2:gr%nnzp-1)
     .                  + cond_rrm(2:gr%nnzp-1) 
     .                  + auto_rrm(2:gr%nnzp-1) + accr_rrm(2:gr%nnzp-1)
        ! auto_rrm, cond_rrm, and accr_rrm are not set at nz=1 or nnzp
        ! There is no condensation/evaporation, autoconversion, or accretion
        ! at level 1, which is below the ground surface.  Brian.
        dd(1)       = (1.0/dt) * rrm(1)
        dd(gr%nnzp) = (1.0/dt) * rrm(gr%nnzp) + cond_rrm(gr%nnzp-1)
     .              + auto_rrm(gr%nnzp-1) + accr_rrm(gr%nnzp-1)

        ! Enter the entire rr array, which is the RHS of the equation.
        rr(2:gr%nnzp-1) = (1.0/dt) * Nrm(2:gr%nnzp-1) 
     .                  + cond_Nrm(2:gr%nnzp-1) + auto_Nrm(2:gr%nnzp-1)
        ! auto_rrm and cond_Nrm are not set at nz=1 or nnzp
        ! There is no condensation/evaporation or autoconversion
        ! at level 1, which is below the ground surface.  Brian.
        rr(1)       = (1.0/dt) * Nrm(1)
        rr(gr%nnzp) = (1.0/dt) * Nrm(gr%nnzp) + cond_Nrm(gr%nnzp-1)
     .              + auto_Nrm(gr%nnzp-1)


!       The Boundary condition setup is now located above the main loop.
!       Brian.  December 15, 2005.


!       solve tridiagonal systems

        !call tridag(aa, bb, cc, dd, rrm, gr%nnzp, isValid)
        call tridag_solve( "rrm", gr%nnzp, 1, cc, bb, aa, dd, 
     .                      rrm, isValid)

        if (.not. isValid ) return

        !call tridag(oo, pp, qq, rr, Nrm, gr%nnzp, isValid)
        call tridag_solve( "Nrm", gr%nnzp, 1, qq, pp, oo, rr, 
     .                      Nrm, isValid)

        if (.not. isValid ) return


        DO k = 1, gr%nnzp, 1

           km1 = MAX(k-1,1)
           kp1 = MIN(k+1,gr%nnzp)

#ifdef STATS

        if (lstats_samp) then

        ! Finalize implicit contributions to rrm.
         if ( irrm_ma > 0 ) then
           zt%x(k,irrm_ma) = zt%x(k,irrm_ma)
     .       + ztscr01(k) * rrm(km1)
     .       + ztscr02(k) * rrm(kp1)
           zt%n(k,irrm_ma) = zt%n(k,irrm_ma) + 1
         endif

         if ( irrm_sd > 0 ) then
           zt%x(k,irrm_sd) = zt%x(k,irrm_sd)
     .       + ztscr03(k) * rrm(km1)
     .       + ztscr04(k) * rrm(kp1)
           zt%n(k,irrm_sd) = zt%n(k,irrm_sd) + 1
         endif

         if ( irrm_dff > 0 ) then
           zt%x(k,irrm_dff) = zt%x(k,irrm_dff)
     .       + ztscr05(k) * rrm(km1)
     .       + ztscr06(k) * rrm(k)
     .       + ztscr07(k) * rrm(kp1)
           zt%n(k,irrm_dff) = zt%n(k,irrm_dff) + 1
         endif


        ! Finalize implicit contributions to Nrm.
         if ( iNrm_ma > 0 ) then
           zt%x(k,iNrm_ma) = zt%x(k,iNrm_ma)
     .       + ztscr08(k) * Nrm(km1)
     .       + ztscr09(k) * Nrm(kp1)
           zt%n(k,iNrm_ma) = zt%n(k,iNrm_ma) + 1
         endif

         if ( iNrm_sd > 0 ) then
           zt%x(k,iNrm_sd) = zt%x(k,iNrm_sd)
     .       + ztscr10(k) * Nrm(km1)
     .       + ztscr11(k) * Nrm(kp1)
           zt%n(k,iNrm_sd) = zt%n(k,iNrm_sd) + 1
         endif

         if ( iNrm_dff > 0 ) then
           zt%x(k,iNrm_dff) = zt%x(k,iNrm_dff)
     .       + ztscr12(k) * Nrm(km1)
     .       + ztscr13(k) * Nrm(k)
     .       + ztscr14(k) * Nrm(kp1)
           zt%n(k,iNrm_dff) = zt%n(k,iNrm_dff) + 1
         endif

        endif

#endif /*STATS*/

           ! Value of rain water mixing ratio cannot fall below 0.
           IF (rrm(k) < 0.0) THEN
              rrm(k) = 0.0
           ENDIF

           ! Value of rain droplet number concentration cannot fall below 0.
           IF (Nrm(k) < 0.0) THEN
              Nrm(k) = 0.0
           ENDIF

        ENDDO  ! k = 1, gr%nnzp, 1


!       boundary conditions on results

        ! rrm
        rrm(1) = rrm(2)
        rrm(gr%nnzp) = rrm(gr%nnzp-1)

        ! Nrm
        Nrm(1) = Nrm(2)
        Nrm(gr%nnzp) = Nrm(gr%nnzp-1)


#ifdef STATS

        if (lstats_samp) then

         if ( irrm_bt > 0 ) then
           zt%x(:,irrm_bt) = zt%x(:,irrm_bt) + rrm/dt
           zt%n(:,irrm_bt) = zt%n(:,irrm_bt) + 1
         endif

         if ( irrm_bt > 0 ) then
           zt%x(:,iNrm_bt) = zt%x(:,iNrm_bt) + Nrm/dt
           zt%n(:,iNrm_bt) = zt%n(:,iNrm_bt) + 1
         endif

        endif

#endif /*STATS*/


        RETURN
        END SUBROUTINE rain

!----------------------------------------------------------------------
        subroutine kk_microphys
     .             ( thlm, p, exner, rhot, thl1, thl2, a, rc1, rc2, s1,
     .               s2, ss1, ss2, rtm, rcm, Ncm, rrm, Nrm, 
     .               lsample, AKm_est, AKm,
     .               rrm_mc_tndcy, Nrm_mc_tndcy, 
     .               hm_rt_tndcy, hm_thl_tndcy,
     .               Vrr, VNr )
!       Description:
!       This code takes Brian Griffin's code from subroutine rain()
!       (above) and computes microphysical tendencies and fall speeds
!       per se rather than adding in the effects of diffusion, advection, 
!       and sedimentation.  This is for the purpose of testing different
!       microphysical schemes with the same diff., advect., and sed.
!----------------------------------------------------------------------

        use model_flags
        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif

        implicit none

        ! External functions
        real, external :: rsat

        ! Input
        real, intent(in), dimension(1:gr%nnzp) ::
     .  thlm,      ! theta-l                            [K]
     .  p,         ! pressure                           [Pa]
     .  exner,     ! exner function                     [-]
     .  rhot,      ! density on thermo. grid            [kg/m^3]
     .  thl1, thl2,! PDF parameters thl1 &thl2          [K]
     .  a,         ! PDF parameter a                    [-]
     .  s1, s2,    ! PDF parameters s1 & s2             [kg/kg]
     .  ss1, ss2,  ! PDF parameters ss1 & ss2           [kg/kg]
     .  rc1, rc2,  ! PDF parameters rc1 & rc2           [kg/kg]
     .  rtm,       ! Total water mixing ratio           [kg/kg]
     .  rcm,       ! Cloud water mixing ratio           [kg/kg]
     .  Ncm,       ! Cloud droplet number conc.         [number/m^3]
     .  rrm,       ! rain water mixing ratio            [kg/kg]
     .  Nrm        ! rain droplet number conc.          [number/m^3]

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
     .  Nrm_mc_tndcy,! Rain droplet number conc. micro. tend.   [(num/m^3)/s]
     .  hm_rt_tndcy, ! Contributions to total water from micro. [(kg/kg)/s]
     .  hm_thl_tndcy ! Contributions to theta_l from micro.     [K/s]

        ! Local variables
        real, dimension(gr%nnzp):: 
     .  rrm_cond,    ! Change in rrm due to condensation        [(kg/kg)/s]
     .  rrm_auto,    ! Change in rrm due to autoconversion      [(kg/kg)/s]
     .  rrm_accr,    ! Change in rrm due to accretion           [(kg/kg)/s]
     .  Nrm_cond,    ! Change in Nrm due to condensation        [(num/m^3)/s]
     .  Nrm_auto,    ! Change in Nrm due to autoconversion      [(num/m^3)/s]
     .  mean_vol_rad,! Mean volume radius                       [m]
     .  rvm,         ! Water vapor mixing ratio                 [kg/kg]
     .  Temp,        ! Temperature                              [K]
     .  es,          ! Saturation vapor pressure                [Pa]
     .  Supsat,      ! Supersaturation                          [-]
     .  rs           ! Saturation mixing ratio                  [kg/kg]

        real :: 
     .  e,      ! Vapor pressure                [Pa]
     .  thm,    ! Potential temperature         [K]
     .  Beta_T  ! Beta_T                        [kg/kg]

        ! Array indices
        integer :: k, km1, kp1, i


        ! Find the drop mean volume radius.  It is calculated using
        ! the rain water mixing ratio, the rain droplet concentration,
        ! and the air density.  These values are taken from the previous
        ! timestep.  It is located on thermodynamic levels.
        do k = 1, gr%nnzp, 1

           mean_vol_rad(k) 
     .     = mean_volume_radius( rrm(k), Nrm(k), rhot(k) )

        end do

#ifdef STATS
        ! Save mean volume radius for stats purposes
        ! Note: added lsample for latin hypercube sampling -dschanen

        if ( lsample .and. lstats_samp 
     .       .and. imean_vol_rad_rain > 0 ) then ! Brian

          zt%x(:,imean_vol_rad_rain) 
     .    = zt%x(:,imean_vol_rad_rain) + mean_vol_rad
          zt%n(:,imean_vol_rad_rain) = zt%n(:,imean_vol_rad_rain) + 1

        end if

#endif /*STATS*/

        ! The sedimentation velocity is found from the drop mean volume radius.
        ! It is located on the momentum levels.  This is due to the fact that
        ! momentum levels are where other vertical velocities are stored, such
        ! as the vertical component of wind velocity (w).

        ! Khairoutdinov and Kogan Sedimentation Velocity Calculation.
        ! sedimentation velocity of rain drops (in m/s).
        ! mean volume radius converted to um.
        ! Note 1: positive sedimentation velocity means downwards.
        ! Note 2: Vrr(m/s) = 0.012*rvr(um) - 0.2

        ! Note 3:
        ! Changed sedimentation to be negative and consistent with
        ! the COAMPS microphysics.
        ! -dschanen 5 Dec 2006

        ! If rvr(um) is too small, the equation will result in an
        ! upwards (positive) sedimentation velocity.  The limiter
        ! is put in there for that reason.

        do k = 1, gr%nnzp, 1

           kp1 = min( k+1, gr%nnzp )


           ! rrm sedimentation velocity.
           Vrr(k) = 0.012 
     .             * ( 1000000.0 * (1.0/2.0) 
     .             * ( mean_vol_rad(kp1) + mean_vol_rad(k) ) )
     .             - 0.2
           ! Negative meaning a downward velocity now -dschanen 5 Dec 2006
           Vrr(k) = -max( Vrr(k), 0.0 )

           ! Nrm sedimentation velocity.
           VNr(k) = 0.007 
     .             * ( 1000000.0 * (1.0/2.0) 
     .             * ( mean_vol_rad(kp1) + mean_vol_rad(k) ) )
     .             - 0.1

           ! Negative meaning a downward velocity now -dschanen 5 Dec 2006
           VNr(k) = -max( VNr(k), 0.0 )

        end do ! 1..gr%nnzp


        ! Find values for other variables.
        do k = 2, gr%nnzp, 1

           ! Find the important conditions

           rvm(k) = rtm(k) - rcm(k)   ! rvm is water vapor mixing ratio.

           ! Find current vapor pressure.
           e = (p(k)*rvm(k))/(ep + rvm(k))

           ! Get theta from theta-l
           thm = thlm(k) + (Lv/(Cp*exner(k)))*rcm(k)

           ! Saturation mixing ratio
           rs(k) = rsat( p(k),thm*exner(k) )

           ! Find temperature and saturation vapor pressure.
           Temp(k) = thm/((p0/p(k))**kappa)
           es(k) = (p(k)*rs(k))/(ep + rs(k))

        end do ! 2..gr%nnzp

        ! Set the boundary conditions
        rvm(1)          = 0.0
        rs(1)           = 0.0
        Temp(1)         = 0.0
        es(1)           = 0.0
        Supsat(1)       = 0.0
        Supsat(gr%nnzp) = 0.0

        do k = 2, gr%nnzp-1, 1

           ! Compute supersaturation via s1, s2. 
           !     Larson et al 2002, JAS, Vol 59, p 3534.
           ! This allows a more direct comparison of local, nonlocal formulas.
           Beta_T = (Rd/Rv) * ( Lv/(Rd*Temp(k)) ) * ( Lv/(Cp*Temp(k)) )

           Supsat(k) = ( a(k)*s1(k) + (1-a(k))*s2(k) )
     .                 *( ( 1.0 + Beta_T*rs(k) ) / rs(k) )

           ! Now find the elements that make up the right-hand side of the
           ! equation, dd, for rain water mixing ratio, rrm.

           rrm_cond(k) 
     .     = cond_evap_rrm( rrm(k), Nrm(k),
     .                      s1(k), ss1(k), s2(k), ss2(k), 
     .                      thl1(k), thl2(k), rc1(k), rc2(k), a(k),
     .                      p(k), rhot(k), Temp(k), Supsat(k)  )


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
     .                        s2(k), ss2(k), a(k), rhot(k) )

           end if ! LH_on
        ! End Vince Larson's addition

           rrm_accr(k) 
     .     = accretion_rrm( rcm(k), rrm(k), s1(k), ss1(k), 
     .                      s2(k), ss2(k), a(k) )

           ! Now find the elements that make up the right-hand side of the
           ! equation, rr, for rain droplet number concentration, Nrm.

           Nrm_cond(k) = cond_evap_Nrm( rrm_cond(k), Nrm(k), rrm(k) )

           Nrm_auto(k) = autoconv_Nrm( rrm_auto(k), rhot(k) )

#ifdef STATS
           if ( lsample .and. lstats_samp ) then
              ! Explicit contributions to rrm.
              if ( irrm_cond > 0 ) then
                zt%x(k,irrm_cond) = zt%x(k,irrm_cond) + rrm_cond(k)
                zt%n(k,irrm_cond) = zt%n(k,irrm_cond) + 1
              end if

              if ( irrm_auto > 0 ) then
                zt%x(k,irrm_auto) = zt%x(k,irrm_auto) + rrm_auto(k)
                zt%n(k,irrm_auto) = zt%n(k,irrm_auto) + 1
              end if

              if ( irrm_accr > 0 ) then
                zt%x(k,irrm_accr) = zt%x(k,irrm_accr) + rrm_accr(k)
                zt%n(k,irrm_accr) = zt%n(k,irrm_accr) + 1
              end if

             ! Explicit contributions to Nrm.
             if ( iNrm_cond > 0 ) then
               zt%x(k,iNrm_cond) = zt%x(k,iNrm_cond) + Nrm_cond(k)
               zt%n(k,iNrm_cond) = zt%n(k,iNrm_cond) + 1
             end if

             if ( iNrm_auto > 0 ) then
               zt%x(k,iNrm_auto) = zt%x(k,iNrm_auto) + Nrm_auto(k)
               zt%n(k,iNrm_auto) = zt%n(k,iNrm_auto) + 1
             end if

           end if ! lstats_samp and lsample
#endif /*STATS*/

          ! Explicit contributions to rrm and Nrm from microphysics
          rrm_mc_tndcy(k) = rrm_cond(k) + rrm_auto(k) + rrm_accr(k)
          Nrm_mc_tndcy(k) = Nrm_cond(k) + Nrm_auto(k)

          ! Explicit contributions to thlm and rtm from the microphysics
          hm_rt_tndcy(k)  = -rrm_mc_tndcy(k)
          hm_thl_tndcy(k) = ( Lv / ( Cp*exner(k) ) ) * rrm_mc_tndcy(k)
        end do ! k=2..gr%nnzp-1

        ! Boundary conditions
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
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^(1/3) * Nr^(-1/3)
        !
        ! This turns into:
        !
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^alpha * Nr^beta
        !
        ! where alpha = 1/3 and beta = -1/3
        !
        ! When applied to a Probability Density Function (PDF), it becomes:
        !
        ! mvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * Y1^alpha * Y2^beta
        !
        ! Y1 is used for rr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        ! Y2 is used for Nr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Notes:  For our PDF-version (non-local), we turn Nrm into
        !         concentration per mass instead of concentration per volume.
        !         It is divided by air density (rho) in order to achieve
        !         that value.  That also explains why rho is removed from
        !         the denominator of the factor [(4*PI*rho_lw)/(3*rho)]^(2/3).
        !         The factor becomes [(4*PI*rho_lw)/(3)]^(2/3).
        !
        !         For our code using the local version, we use KK(2000)
        !         equation 22 exactly, as listed above.  We do not mess around
        !         with changing anything.
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION mean_volume_radius( rrm, Nrm, rhot )

        USE constants
        USE model_flags

        ! Input variables.
        REAL, INTENT(IN):: rrm    ! Grid-box average rrm
!        REAL, INTENT(IN):: rrp2   ! Grid-box rr variance
        REAL, INTENT(IN):: Nrm    ! Grid-box average Nrm
!        REAL, INTENT(IN):: Nrp2   ! Grid-box Nr variance
        REAL, INTENT(IN):: rhot   ! Grid-box average density (t-level)

        ! Output variables.
        REAL:: mean_volume_radius

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of rr
        REAL:: beta_exp   ! Exponent of Nr

        ! Original terms.
        REAL:: mu_rr       ! Grid-box average of rr
        REAL:: sigma_rr    ! Grid-box standard deviation of rr
        REAL:: mu_Nr       ! Grid-box average of Nr
        REAL:: sigma_Nr    ! Grid-box standard deviation of Nr
        REAL:: corr_rrNr   ! Correlation of rr and Nr

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              mean_volume_radius = 
     .           (  ((4.0*pi*rho_lw)/(3.0*rhot))**(-1.0/3.0)  )
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
              sigma_rr = rrm * SQRT(rrp2_rrm2)  ! rr'2/rrm^2 = 1.2

              ! Nr is distributed Lognormally.
              mu_Nr = (Nrm/rhot)
!              sigma_Nr = SQRT(Nrp2)
              sigma_Nr = (Nrm/rhot) * SQRT(Nrp2_Nrm2)  ! Nr'2/Nrm^2 = 1.0

              ! Correlations.
              corr_rrNr = corr_rrNr_LL  ! corr_rcNr_LL =  0.85


              mean_volume_radius =
     .           (  ((4.0*pi*rho_lw)/3.0)**(-1.0/3.0)  )
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
        !      * rr^(1/3) * Nr^(2/3) * S
        !
        ! where S = (e/es) - 1
        !
        ! since S = [ ( 1 + Beta_T * rs ) / rs ] * s
        !
        ! (where s is the extended liquid water specific humidity),
        !
        ! the equation becomes:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * rr^(1/3) * Nr^(2/3) * [ ( 1 + Beta_T * rs ) / rs ] * s
        !
        ! This turns into:
        !
        ! (drr/dt)cond = 
        !    3 * Cevap * G(T,p) * [(4*PI*rho_lw)/(3*rho)]^(2/3) 
        !      * [ ( 1 + Beta_T * rs ) / rs ]
        !      * s^alpha * rr^beta * Nr^gamma
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
        ! Y3 is used for Nr.  It is distributed as a log-normal with a
        ! domain from zero to infinity.
        !
        ! Notes:  For our PDF-version (non-local), we turn Nrm into
        !         concentration per mass instead of concentration per volume.
        !         It is divided by air density (rho) in order to achieve
        !         that value.  That also explains why rho is removed from
        !         the denominator of the factor [(4*PI*rho_lw)/(3*rho)]^(2/3).
        !         The factor becomes [(4*PI*rho_lw)/(3)]^(2/3).
        !
        !         For our code using the local version, we use KK(2000)
        !         equation 22 exactly, as listed above.  We do not mess around
        !         with changing anything.
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION cond_evap_rrm( rrm, Nrm, s1, ss1, s2, ss2,
     .                          thl1, thl2, rc1, rc2, a,
     .                          press, rhot, Temp, Supsat  )

        USE constants
        USE model_flags

        ! External functions

        real rsat
        external rsat

        ! Input variables.
        REAL, INTENT(IN):: rrm    ! Grid-box average rrm
!        REAL, INTENT(IN):: rrp2   ! Grid-box rr variance
        REAL, INTENT(IN):: Nrm    ! Grid-box average Nrm
!        REAL, INTENT(IN):: Nrp2   ! Grid-box Nr variance
        REAL, INTENT(IN):: s1     ! Plume 1 average s
        REAL, INTENT(IN):: ss1    ! Plume 1 sigma s1 (not sigma^2 s1)
        REAL, INTENT(IN):: s2     ! Plume 2 average s
        REAL, INTENT(IN):: ss2    ! Plume 2 sigma s2 (not sigma^2 s2)
        REAL, INTENT(IN):: thl1   ! Plume 1 average theta-l
        REAL, INTENT(IN):: thl2   ! Plume 2 average theta-l
        REAL, INTENT(IN):: rc1    ! Plume 1 average rc
        REAL, INTENT(IN):: rc2    ! Plume 2 average rc
        REAL, INTENT(IN):: a      ! Relative weight of each individual
                                  ! Gaussian "plume."
        REAL, INTENT(IN):: press  ! Grid-box average pressure
        REAL, INTENT(IN):: rhot   ! Grid-box average density (t-level)
        REAL, INTENT(IN):: Temp   ! Grid-box average Temperature
        REAL, INTENT(IN):: Supsat ! Grid-box average Supersaturation

        ! Output variables.
        REAL:: cond_evap_rrm

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s
        REAL:: beta_exp   ! Exponent of rr
        REAL:: gamma_exp  ! Exponent of Nr

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s
        REAL:: sigma_s1    ! Plume 1 standard deviation of s
        REAL:: mu_s2       ! Plume 2 average of s
        REAL:: sigma_s2    ! Plume 2 standard deviation of s
        REAL:: mu_rr       ! Grid-box average of rr
        REAL:: sigma_rr    ! Grid-box standard deviation of rr
        REAL:: mu_Nr       ! Grid-box average of Nr
        REAL:: sigma_Nr    ! Grid-box standard deviation of Nr
        REAL:: corr_rcrr   ! Correlation of rc and rr
        REAL:: corr_rcNr   ! Correlation of rc and Nr
        REAL:: corr_rrNr   ! Correlation of rr and Nr

        REAL:: Tl_1, Tl_2, T_1, T_2, rsl_1, rsl_2, Beta_T1, Beta_T2
        REAL:: plume_1_constants, plume_2_constants

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rrm > rr_tol .AND. Nrm > Nr_tol ) THEN

              cond_evap_rrm =
     .           3.0 * C_evap * G_T_p( Temp, press )
     .               * ( ((4.0*pi*rho_lw)/(3.0*rhot))**(2.0/3.0) )
     .               * (rrm**(1.0/3.0)) * (Nrm**(2.0/3.0)) * Supsat

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
              Tl_1 = thl1 * ((press/p0)**kappa)

              T_1 = Tl_1 + (Lv/Cp)*rc1

              rsl_1 = rsat(press, Tl_1)

              Beta_T1 = (Rd/Rv) * ( Lv/(Rd*Tl_1) ) 
     .                          * ( Lv/(Cp*Tl_1) )

              plume_1_constants = 
     .           3.0 * C_evap * G_T_p( T_1, press )
     .               * ( ((4.0*pi*rho_lw)/3.0)**(2.0/3.0) )
     .               * ( (1.0 + Beta_T1*rsl_1) / rsl_1 )

              ! Gaussian "plume" 2
              Tl_2 = thl2 * ((press/p0)**kappa)

              T_2 = Tl_2 + (Lv/Cp)*rc2

              rsl_2 = rsat(press, Tl_2)

              Beta_T2 = (Rd/Rv) * ( Lv/(Rd*Tl_2) ) 
     .                          * ( Lv/(Cp*Tl_2) )

              plume_2_constants = 
     .           3.0 * C_evap * G_T_p( T_2, press )
     .               * ( ((4.0*pi*rho_lw)/3.0)**(2.0/3.0) )
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
              sigma_rr = rrm * SQRT(rrp2_rrm2)  ! rr'2/rrm^2 = 1.2

              ! Nr is distributed Lognormally.
              mu_Nr = (Nrm/rhot)
!              sigma_Nr = SQRT(Nrp2)
              sigma_Nr = (Nrm/rhot) * SQRT(Nrp2_Nrm2)  ! Nr'2/Nrm^2 = 1.0

              ! Correlations.
              corr_rcrr = corr_rcrr_NL  ! corr_rcrr_NL =  0.10
              corr_rcNr = corr_rcNr_NL  ! corr_rcNr_NL = -0.15
              corr_rrNr = corr_rrNr_LL  ! corr_rrNr_LL =  0.85


              cond_evap_rrm = 
     .         ( a ) * plume_1_constants
     .         * PDF_TRIVAR_2G_LN_LN ( mu_s1, mu_rr, mu_Nr,
     .                                 sigma_s1, sigma_rr, sigma_Nr,
     .                                 corr_rcrr, corr_rcNr, corr_rrNr,
     .                                 alpha_exp, beta_exp, gamma_exp )
     .        +
     .         (1-a) * plume_2_constants
     .         * PDF_TRIVAR_2G_LN_LN ( mu_s2, mu_rr, mu_Nr,
     .                                 sigma_s2, sigma_rr, sigma_Nr,
     .                                 corr_rcrr, corr_rcNr, corr_rrNr,
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

        REAL, INTENT(IN):: cond_rrm, Nrm, rrm
        REAL:: cond_evap_Nrm

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
        ! (drr/dt)auto = 1350 * rc^2.47 * Nc^-1.79
        !
        ! where Nc is in concentration per cm^3.  Since there are 10^-6
        ! cubic meters in one cubic centimenter, the equation becomes:
        !
        ! (drr/dt)auto = 1350 * rc^2.47 * ((10^-6)Nc)^-1.79
        !
        ! simplified:
        !
        ! (drr/dt)auto = 1350 * (10^-6)^-1.79 * rc^2.47 * Nc^-1.79
        !
        ! This turns into:
        !
        ! (drr/dt)auto = 1350 * 10^10.74 * rc^2.47 * Nc^-1.79
        !
        ! with Nc in concentration per m^3.
        !
        ! For our PDF version, we want to express Nc in concentration per mass
        ! instead of concentration per volume.  We end up with:
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
        ! Notes:  For our code using the local version, we use KK(2000)
        !         equation 33, only modifying it to express Nc in
        !         concentration per m^3.  We do not mess around with 
        !         changing anything else.
        !
        ! Brian Griffin.  November 14, 2006.
        !
        !-----------------------------------------------------------------------
        
        FUNCTION autoconv_rrm( rcm, Ncm, s1, ss1, s2, ss2, a, rhot )

        USE constants
        USE model_flags

        ! Input variables.
        REAL, INTENT(IN):: rcm    ! Grid-box average rcm
!        REAL, INTENT(IN):: rcp2   ! Grid-box rc variance
        REAL, INTENT(IN):: Ncm    ! Grid-box average Ncm
!        REAL, INTENT(IN):: Ncp2   ! Grid-box Nc variance
        REAL, INTENT(IN):: s1     ! Plume 1 average s
        REAL, INTENT(IN):: ss1    ! Plume 1 sigma s1 (not sigma^2 s1)
        REAL, INTENT(IN):: s2     ! Plume 2 average s
        REAL, INTENT(IN):: ss2    ! Plume 2 sigma s2 (not sigma^2 s2)
        REAL, INTENT(IN):: a      ! Relative weight of each individual
                                  ! Gaussian "plume."
        REAL, INTENT(IN):: rhot   ! Grid-box average density (t-level)

        ! Output variables.
        REAL:: autoconv_rrm

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s1
        REAL:: beta_exp   ! Exponent of Nc

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s1
        REAL:: sigma_s1    ! Plume 1 standard deviation of s1
        REAL:: mu_s2       ! Plume 2 average of s2
        REAL:: sigma_s2    ! Plume 2 standard deviation of s2
        REAL:: mu_Nc       ! Grid-box average of Nc
        REAL:: sigma_Nc    ! Grid-box standard deviation of Nc
        REAL:: corr_rcNc   ! Correlation of rc and Nc

!-------------------------------------------------------------------------------

        IF ( local_kk ) THEN

           ! Tolerance values are used instead of 0 in order to prevent
           ! numerical error.
           IF ( rcm > rc_tol .AND. Ncm > Nc_tol ) THEN

              autoconv_rrm = 7.4188E13 * rcm**2.47 * Ncm**(-1.79)

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
              mu_Nc = (Ncm/rhot)
!              sigma_Nc = SQRT(Ncp2)
              sigma_Nc = (Ncm/rhot) * SQRT(Ncp2_Ncm2)  ! Nc'2/Ncm^2 = 0.07

              ! Correlations.
              corr_rcNc = corr_rcNc_NL  ! corr_rcNc_NL =  0.45


              autoconv_rrm = 7.4188E13 * rhot**(-1.79) * (
     .             ( a ) 
     .           * PDF_BIVAR_2G_LN ( mu_s1, mu_Nc, sigma_s1, sigma_Nc,
     .                               corr_rcNc, alpha_exp, beta_exp )
     .         +   (1-a)
     .           * PDF_BIVAR_2G_LN ( mu_s2, mu_Nc, sigma_s2, sigma_Nc,
     .                               corr_rcNc, alpha_exp, beta_exp )
     .                                                   )

           ELSE

              ! If either rcm or Ncm are 0.
              autoconv_rrm = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION autoconv_rrm

!===============================================================================

        FUNCTION autoconv_Nrm( auto_rrm, rhot )

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: auto_rrm, rhot
        REAL:: autoconv_Nrm

        autoconv_Nrm = auto_rrm /
     .                 ( ((4.0*pi*rho_lw)/(3.0*rhot)) * (r_0**3) )

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
        
        FUNCTION accretion_rrm( rcm, rrm, s1, ss1, s2, ss2, a )

        USE constants
        USE model_flags

        ! Input variables.
        REAL, INTENT(IN):: rcm    ! Grid-box average rcm
!        REAL, INTENT(IN):: rcp2   ! Grid-box rc variance
        REAL, INTENT(IN):: rrm    ! Grid-box average rrm
!        REAL, INTENT(IN):: rrp2   ! Grid-box rr variance
        REAL, INTENT(IN):: s1     ! Plume 1 average s
        REAL, INTENT(IN):: ss1    ! Plume 1 sigma s1 (not sigma^2 s1)
        REAL, INTENT(IN):: s2     ! Plume 2 average s
        REAL, INTENT(IN):: ss2    ! Plume 2 sigma s2 (not sigma^2 s2)
        REAL, INTENT(IN):: a      ! Relative weight of each individual
                                  ! Gaussian "plume."

        ! Output variables.
        REAL:: accretion_rrm

        ! Exponential terms.
        REAL:: alpha_exp  ! Exponent of s
        REAL:: beta_exp   ! Exponent of rr

        ! Original terms.
        REAL:: mu_s1       ! Plume 1 average of s
        REAL:: sigma_s1    ! Plume 1 standard deviation of s
        REAL:: mu_s2       ! Plume 2 average of s
        REAL:: sigma_s2    ! Plume 2 standard deviation of s
        REAL:: mu_rr       ! Grid-box average of rr
        REAL:: sigma_rr    ! Grid-box standard deviation of rr
        REAL:: corr_rcrr   ! Correlation of rc and rr

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
              sigma_rr = rrm * SQRT(rrp2_rrm2)  ! rr'2/rrm^2 = 1.2

              ! Correlations.
              corr_rcrr = corr_rcrr_NL  ! corr_rcrr_NL =  0.10


              accretion_rrm = 67.0 * (
     .             ( a ) 
     .           * PDF_BIVAR_2G_LN ( mu_s1, mu_rr, sigma_s1, sigma_rr,
     .                               corr_rcrr, alpha_exp, beta_exp )
     .         +   (1-a)
     .           * PDF_BIVAR_2G_LN ( mu_s2, mu_rr, sigma_s2, sigma_rr,
     .                               corr_rcrr, alpha_exp, beta_exp )
     .                               )

           ELSE

              ! If either rcm or rrm are 0.
              accretion_rrm = 0.0

           ENDIF

        ENDIF

        RETURN

        END FUNCTION accretion_rrm

!===============================================================================

        FUNCTION G_T_p( Temp, press )

        USE constants

        IMPLICIT NONE


        ! Here we compute G(T,p) as in KK (17) 
        !      and Eq. (7.17) of Rogers and Yau 1989

        ! External functions

        real, external :: rsat

        ! Input      
        REAL, INTENT(IN):: Temp, press

        ! Ouput
        REAL:: G_T_p

        ! Internal
        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: es, rs, Celsius

        Celsius = Temp - 273.16

        Ka = (5.69 + 0.017*Celsius)*0.00001  ! Ka in cal./(cm.*sec.*C)
        Ka = 4.1868*100.0*Ka  ! Ka in J./(m.*sec.*K)

        Dv = 0.221*((Temp/273.16)**1.94)*(101325.0/press)  
                                ! Dv in (cm.^2)/sec.  ! .221 is correct.
        Dv = Dv/10000.0  ! Dv in (m.^2)/sec.

        rs = rsat(press, Temp)
        es = (press*rs)/(ep + rs)

        Fk = (Lv/(Rv*Temp) - 1.0)*(Lv*rho_lw)/(Ka*Temp)
        Fd = (rho_lw*Rv*Temp)/(Dv*es)

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
        ! Gaussian.
        !
        ! P_1(X1,X2,X3) and P_2(X1,X2,X3) are simply the equation for a
        ! Trivariate Gaussian Distribution.  The only difference between the
        ! two is the dependence on each individual plume for the plume average
        ! and plume variance for the X1 variable.
        !
        ! The Probability Density Function for a Trivariate Normal Distribution:
        !
        ! P_i(X1,X2,X3) = 
        !
        ! 1 /
        !     {
        !       (2*PI)^(3/2) * sigmaX1i * sigmaX2 * sigmaX3 
        !       * SQRT[ 1 - ( corrX1X2^2 + corrX1X3^2 + corrX2X3^2 ) 
        !                 + 2*corrX1X2*corrX1X3*corrX2X3 ]
        !     }
        ! * EXP{ -(1/2) * phi }
        !
        ! where,
        !
        ! phi =
        !
        ! 1 /
        !     [
        !       1 - ( corrX1X2^2 + corrX1X3^2 + corrX2X3^2 ) 
        !         + 2*corrX1X2*corrX1X3*corrX2X3
        !     ]
        ! * {
        !       [ (1-corrX2X3^2) / sigmaX1i^2 ] * ( X1i - muX1i )^2
        !     + [ (1-corrX1X3^2) / sigmaX2^2 ] * ( X2 - muX2 )^2
        !     + [ (1-corrX1X2^2) / sigmaX3^2 ] * ( X3 - muX3 )^2
        !     + [ 2*(corrX1X3*corrX2X3-corrX1X2) / (sigmaX1i*sigmaX2) ]
        !       * ( X1i - muX1i ) * ( X2 - muX2 )
        !     + [ 2*(corrX1X2*corrX2X3-corrX1X3) / (sigmaX1i*sigmaX3) ]
        !       * ( X1i - muX1i ) * ( X3 - muX3 )
        !     + [ 2*(corrX1X2*corrX1X3-corrX2X3) / (sigmaX2*sigmaX3) ]
        !       * ( X2 - muX2 ) * ( X3 - muX3 )
        !   }
        !
        ! Definitions:
        !
        ! muX1i:  The plume average of the X1 (or Y1) variable.
        ! muX2:  The average of the X2 (or ln Y2) variable.
        ! muX3:  The average of the X3 (or ln Y3) variable.
        !
        ! sigmaX1i:  The plume standard deviation (square-root of variance)
        !            of the X1 (or Y1) variable.
        ! sigmaX2:  The standard deviation (square-root of variance) of the
        !           X2 (or ln Y2) variable.
        ! sigmaX3:  The standard deviation (square-root of variance) of the 
        !           X3 (or ln Y3) variable.
        !
        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
        !            Y1 and ln Y2).
        ! corrX1X3:  Intra-Gaussian correlation between X1 and X3 (or between
        !            Y1 and ln Y3).
        ! corrX2X3:  Intra-Gaussian correlation between X2 and X3 (or between
        !            ln Y2 and ln Y3).
        !
        ! After careful (and long) triple-integration of the above equation,
        ! one gets the result for each plume of the trivariate normal
        ! -- log-normal -- log-normal equation.
        !
        ! The result for each plume is (Equation #1):
        !
        !   [ 1 / SQRT(2*PI) ] * ( -sigmaX1i )^alpha 
        ! * EXP{ muX2*beta + muX3*gamma }
        ! * EXP{ (1/2) * [  (1-corrX1X2^2)*sigmaX2^2*beta^2 
        !                 + (1-corrX1X3^2)*sigmaX3^2*gamma^2 
        !                 + 2*(corrX2X3-corrX1X2*corrX1X3)
        !                    *sigmaX2*beta*sigmaX3*gamma
        !                ]
        !      }
        ! * EXP{  (1/4) * [  ( muX1i / sigmaX1i ) 
        !                  + corrX1X2*sigmaX2*beta
        !                  + corrX1X3*sigmaX3*gamma
        !                 ]^2
        !       - ( muX1i / sigmaX1i ) *
        !                 [  ( muX1i / sigmaX1i ) 
        !                  + corrX1X2*sigmaX2*beta
        !                  + corrX1X3*sigmaX3*gamma
        !                 ]
        !       + (1/2) * ( muX1i^2 / sigmaX1i^2 )
        !      }
        ! * GAMMA_FNC { alpha + 1 }
        ! * PARAB_CYL_FNC_[-(alpha+1)] {  ( muX1i / sigmaX1i ) 
        !                               + corrX1X2*sigmaX2*beta
        !                               + corrX1X3*sigmaX3*gamma
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
        ! cylinder function is made.
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
        ! * EXP{ (1/2) * [  (1-corrX1X2^2)*sigmaX2^2*beta^2 
        !                 + (1-corrX1X3^2)*sigmaX3^2*gamma^2 
        !                 + 2*(corrX2X3-corrX1X2*corrX1X3)
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
     .                                 corr_rcrr, corr_rcNr, corr_rrNr,
     .                                 alpha_exp, beta_exp, gamma_exp )

        USE constants
        USE parabolic, ONLY: gamma

        IMPLICIT NONE

        ! Input variables.
        REAL, INTENT(IN):: mu_si     ! Plume (i) average of s
        REAL, INTENT(IN):: mu_rr     ! Average of rr
        REAL, INTENT(IN):: mu_Nr     ! Average of Nr
        REAL, INTENT(IN):: sigma_si  ! Plume (i) standard deviation of s
        REAL, INTENT(IN):: sigma_rr  ! Standard deviation of rr
        REAL, INTENT(IN):: sigma_Nr  ! Standard deviation of Nr
        REAL, INTENT(IN):: corr_rcrr ! Intra-Gaussian correlation between rc, rr
        REAL, INTENT(IN):: corr_rcNr ! Intra-Gaussian correlation between rc, Nr
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
        REAL:: corr_rcrrG  ! Intra-Gaussian correlation between rc(G), rr(G)
        REAL:: corr_rcNrG  ! Intra-Gaussian correlation between rc(G), Nr(G)
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

        ! corr_rcrr is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_rcrrG, which is a correlation between
        ! two Gaussians.
        corr_rcrrG = corr_rcrr
     .                   * SQRT( EXP(sigma_rrG**2.0) - 1.0 ) / sigma_rrG
        ! corr_rcNr is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_rcNrG, which is a correlation between
        ! two Gaussians.
        corr_rcNrG = corr_rcNr
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
     .           + corr_rcrrG*sigma_rrG*beta_exp 
     .           + corr_rcNrG*sigma_NrG*gamma_exp

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

           PDF_TRIVAR_2G_LN_LN =
     .       (1.0/SQRT(2.0*pi)) * (-sigma_siG)**alpha_exp 
     .     * EXP( mu_rrG*beta_exp + mu_NrG*gamma_exp )
     .     * EXP( (1.0/2.0) * (
     .                            ( 1.0 - corr_rcrrG**2.0 )
     .                           *(sigma_rrG**2.0)*(beta_exp**2.0)
     .                         +  ( 1.0 - corr_rcNrG**2.0 )
     .                           *(sigma_NrG**2.0)*(gamma_exp**2.0)
     .                         + 2.0
     .                           *( corr_rrNrG - corr_rcrrG*corr_rcNrG )
     .                           *sigma_rrG*beta_exp*sigma_NrG*gamma_exp
     .                        )
     .          )
     .     * EXP(   (1.0/4.0) * ( scci )**2.0
     .            - (mu_siG/sigma_siG) * ( scci )
     .            + (1.0/2.0) * ( (mu_siG**2.0) / (sigma_siG**2.0) )
     .          )
     .     * GAMMA( gamma_fnc_input )
     .     * Dv_fnc( parab_cyl_fnc_ord, scci )

        ELSEIF ( .not. use_eq_1 .AND. mu_siG < 0.0 ) THEN

           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is less than 0.

           PDF_TRIVAR_2G_LN_LN =
     .       (mu_siG)**alpha_exp 
     .     * EXP( mu_rrG*beta_exp + mu_NrG*gamma_exp )
     .     * EXP( (1.0/2.0) * (
     .                            ( 1.0 - corr_rcrrG**2.0 )
     .                           *(sigma_rrG**2.0)*(beta_exp**2.0)
     .                         +  ( 1.0 - corr_rcNrG**2.0 )
     .                           *(sigma_NrG**2.0)*(gamma_exp**2.0)
     .                         + 2.0
     .                           *( corr_rrNrG - corr_rcrrG*corr_rcNrG )
     .                           *sigma_rrG*beta_exp*sigma_NrG*gamma_exp
     .                        )
     .          )

        ELSEIF ( .not. use_eq_1 .AND. mu_siG >=  0.0 ) THEN
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
        ! Gaussian.
        !
        ! P_1(X1,X2) and P_2(X1,X2) are simply the equation for a
        ! Bivariate Gaussian Distribution.  The only difference between the
        ! two is the dependence on each individual plume for the plume average
        ! and plume variance for the X1 variable.
        !
        ! The Probability Density Function for a Bivariate Normal Distribution:
        !
        ! P_i(X1,X2) = 
        !
        ! 1 / { 2 * PI * sigmaX1i * sigmaX2 * SQRT[ 1 - corrX1X2^2 ] }
        ! * EXP{ -(1/2) * phi }
        !
        ! where,
        !
        ! phi =
        !
        ! 1 / [ 1 - corrX1X2^2 ]
        ! * {
        !       [ 1 / sigmaX1i^2 ] * ( X1i - muX1i )^2
        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
        !     + [ 2*corrX1X2 / (sigmaX1i*sigmaX2) ]
        !       * ( X1i - muX1i ) * ( X2 - muX2 )
        !   }
        !
        ! Definitions:
        !
        ! muX1i:  The plume average of the X1 (or Y1) variable.
        ! muX2:  The average of the X2 (or ln Y2) variable.
        !
        ! sigmaX1i:  The plume standard deviation (square-root of variance)
        !            of the X1 (or Y1) variable.
        ! sigmaX2:  The standard deviation (square-root of variance) of the
        !           X2 (or ln Y2) variable.
        !
        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
        !            Y1 and ln Y2).
        !
        ! After careful (and long) double-integration of the above equation,
        ! one gets the result for each plume of the bivariate normal
        ! -- log-normal equation.
        !
        ! The result for each plume is (Equation #1):
        !
        !   [ 1 / SQRT(2*PI) ] * ( sigmaX1i )^alpha 
        ! * EXP{    muX2*beta + (1/2)*sigmaX2^2*beta^2
        !        - (1/4) * [ ( muX1i / sigmaX1i ) + corrX1X2*sigmaX2*beta ]^2
        !      }
        ! * GAMMA_FNC { alpha + 1 }
        ! * PARAB_CYL_FNC_[-(alpha+1)] { - [  ( muX1i / sigmaX1i )
        !                                   + corrX1X2*sigmaX2*beta ]
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
     .                             corr_rcxx, alpha_exp, beta_exp )

        USE constants
        USE parabolic, ONLY: gamma

        IMPLICIT NONE

        ! Input variables.
        REAL, INTENT(IN):: mu_si     ! Plume (i) average of s
        REAL, INTENT(IN):: mu_xx     ! Average of xx
        REAL, INTENT(IN):: sigma_si  ! Plume (i) standard deviation of s
        REAL, INTENT(IN):: sigma_xx  ! Standard deviation of Y2
        REAL, INTENT(IN):: corr_rcxx ! Intra-Gaussian correlation between rc, xx
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
        REAL:: corr_rcxxG  ! Intra-Gaussian correlation between rc(G), xx(G)

        DOUBLE PRECISION:: gamma_fnc_input
        DOUBLE PRECISION:: parab_cyl_fnc_ord
        DOUBLE PRECISION:: sci
        DOUBLE PRECISION:: test
        LOGICAL:: use_eq_1
        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308
        REAL, PARAMETER:: sigs_tol = 10.0**(-18)

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

        ! corr_rcxx is a correlation between a Gaussian and a Log-normal.
        ! It must be converted to corr_rcxxG, which is a correlation between
        ! two Gaussians.
        corr_rcxxG = corr_rcxx
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
     .          + corr_rcxxG*sigma_xxG*beta_exp 

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

           PDF_BIVAR_2G_LN =
     .       (1.0/SQRT(2.0*pi)) * (sigma_siG)**alpha_exp 
     .     * EXP(   mu_xxG*beta_exp
     .            + (1.0/2.0)*(sigma_xxG**2.0)*(beta_exp**2.0)
     .            - (1.0/4.0)*( sci )**2.0
     .          )
     .     * GAMMA( gamma_fnc_input )
     .     * Dv_fnc( parab_cyl_fnc_ord, -sci )

        ELSEIF ( .not. use_eq_1 .AND. mu_siG > 0.0 ) THEN

           ! Case where sigma_siG is equal to 0 (or sigma_siG is sufficiently
           ! small enough to cause the parabolic cylinder function to produce a
           ! result greater than what the computer can represent numerically)
           ! and mu_siG is greater than 0.

           PDF_BIVAR_2G_LN =
     .       (mu_siG)**alpha_exp 
     .     * EXP(   mu_xxG*beta_exp
     .            + (1.0/2.0)*(sigma_xxG**2.0)*(beta_exp**2.0) )

        ELSEIF ( .not. use_eq_1 .AND. mu_siG <=  0.0 ) THEN
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

        USE constants

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
!        ! Gaussian.
!        !
!        ! P_1(X1,X2,X3) and P_2(X1,X2,X3) are simply the equation for a
!        ! Trivariate Gaussian Distribution.  The only difference between the
!        ! two is the dependence on each individual plume for the plume average
!        ! and plume variance for the X1 variable.
!        !
!        ! The Probability Density Function for a Trivariate Normal Distribution:
!        !
!        ! P_i(X1,X2,X3) = 
!        !
!        ! 1 /
!        !     {
!        !       (2*PI)^(3/2) * sigmaX1i * sigmaX2 * sigmaX3 
!        !       * SQRT[ 1 - ( corrX1X2^2 + corrX1X3^2 + corrX2X3^2 ) 
!        !                 + 2*corrX1X2*corrX1X3*corrX2X3 ]
!        !     }
!        ! * EXP{ -(1/2) * phi }
!        !
!        ! where,
!        !
!        ! phi =
!        !
!        ! 1 /
!        !     [
!        !       1 - ( corrX1X2^2 + corrX1X3^2 + corrX2X3^2 ) 
!        !         + 2*corrX1X2*corrX1X3*corrX2X3
!        !     ]
!        ! * {
!        !       [ (1-corrX2X3^2) / sigmaX1i^2 ] * ( X1i - muX1i )^2
!        !     + [ (1-corrX1X3^2) / sigmaX2^2 ] * ( X2 - muX2 )^2
!        !     + [ (1-corrX1X2^2) / sigmaX3^2 ] * ( X3 - muX3 )^2
!        !     + [ 2*(corrX1X3*corrX2X3-corrX1X2) / (sigmaX1i*sigmaX2) ]
!        !       * ( X1i - muX1i ) * ( X2 - muX2 )
!        !     + [ 2*(corrX1X2*corrX2X3-corrX1X3) / (sigmaX1i*sigmaX3) ]
!        !       * ( X1i - muX1i ) * ( X3 - muX3 )
!        !     + [ 2*(corrX1X2*corrX1X3-corrX2X3) / (sigmaX2*sigmaX3) ]
!        !       * ( X2 - muX2 ) * ( X3 - muX3 )
!        !   }
!        !
!        ! Definitions:
!        !
!        ! muX1i:  The plume average of the X1 (or Y1) variable.
!        ! muX2:  The average of the X2 (or ln Y2) variable.
!        ! muX3:  The average of the X3 (or ln Y3) variable.
!        !
!        ! sigmaX1i:  The plume standard deviation (square-root of variance)
!        !            of the X1 (or Y1) variable.
!        ! sigmaX2:  The standard deviation (square-root of variance) of the
!        !           X2 (or ln Y2) variable.
!        ! sigmaX3:  The standard deviation (square-root of variance) of the 
!        !           X3 (or ln Y3) variable.
!        !
!        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
!        !            Y1 and ln Y2).
!        ! corrX1X3:  Intra-Gaussian correlation between X1 and X3 (or between
!        !            Y1 and ln Y3).
!        ! corrX2X3:  Intra-Gaussian correlation between X2 and X3 (or between
!        !            ln Y2 and ln Y3).
!        !
!        ! After careful (and long) triple-integration of the above equation,
!        ! one gets the result for each plume of the trivariate normal
!        ! -- log-normal -- log-normal equation.
!        !
!        ! The result for each plume is (Equation #1):
!        !
!        !   [ 1 / SQRT(2*PI) ] * ( -sigmaX1i )^alpha 
!        ! * EXP{ muX2*beta + muX3*gamma }
!        ! * EXP{ (1/2) * [  (1-corrX1X2^2)*sigmaX2^2*beta^2 
!        !                 + (1-corrX1X3^2)*sigmaX3^2*gamma^2 
!        !                 + 2*(corrX2X3-corrX1X2*corrX1X3)
!        !                    *sigmaX2*beta*sigmaX3*gamma
!        !                ]
!        !      }
!        ! * EXP{  (1/4) * [  ( muX1i / sigmaX1i ) 
!        !                  + corrX1X2*sigmaX2*beta
!        !                  + corrX1X3*sigmaX3*gamma
!        !                 ]^2
!        !       - ( muX1i / sigmaX1i ) *
!        !                 [  ( muX1i / sigmaX1i ) 
!        !                  + corrX1X2*sigmaX2*beta
!        !                  + corrX1X3*sigmaX3*gamma
!        !                 ]
!        !       + (1/2) * ( muX1i^2 / sigmaX1i^2 )
!        !      }
!        ! * GAMMA_FNC { alpha + 1 }
!        ! * PARAB_CYL_FNC_[-(alpha+1)] {  ( muX1i / sigmaX1i ) 
!        !                               + corrX1X2*sigmaX2*beta
!        !                               + corrX1X3*sigmaX3*gamma
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
!        ! * EXP{ (1/2) * [  (1-corrX1X2^2)*sigmaX2^2*beta^2 
!        !                 + (1-corrX1X3^2)*sigmaX3^2*gamma^2 
!        !                 + 2*(corrX2X3-corrX1X2*corrX1X3)
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
!     .                                 corrY1Y2, corrY1Y3, corrY2Y3,
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
!        REAL, INTENT(IN):: corrY1Y2  ! Intra-Gaussian correlation between Y1,Y2
!        REAL, INTENT(IN):: corrY1Y3  ! Intra-Gaussian correlation between Y1,Y3
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
!        REAL:: corrX1X2  ! Intra-Gaussian correlation between X1,X2
!        REAL:: corrX1X3  ! Intra-Gaussian correlation between X1,X3
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
!        ! corrY1Y2 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1X2, which is a correlation between
!        ! two Gaussians.
!        corrX1X2 = corrY1Y2 
!     .                   * SQRT( EXP(sigmaX2**2.0) - 1.0 ) / sigmaX2
!        ! corrY1Y3 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1X3, which is a correlation between
!        ! two Gaussians.
!        corrX1X3 = corrY1Y3
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
!     .                          + corrX1X2*sigmaX2*beta_exp 
!     .                          + corrX1X3*sigmaX3*gamma_exp
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
!     .                            ( 1.0 - corrX1X2**2.0 )
!     .                           *(sigmaX2**2.0)*(beta_exp**2.0)
!     .                         +  ( 1.0 - corrX1X3**2.0 )
!     .                           *(sigmaX3**2.0)*(gamma_exp**2.0)
!     .                         + 2.0*( corrX2X3 - corrX1X2*corrX1X3 )
!     .                              *sigmaX2*beta_exp*sigmaX3*gamma_exp
!     .                        )
!     .          )
!     .     * EXP(   (1.0/4.0) * (
!     .                             (muX1i/sigmaX1i)
!     .                           + corrX1X2*sigmaX2*beta_exp 
!     .                           + corrX1X3*sigmaX3*gamma_exp
!     .                          )**2.0
!     .            - (muX1i/sigmaX1i) * (
!     .                                    (muX1i/sigmaX1i)
!     .                                  + corrX1X2*sigmaX2*beta_exp 
!     .                                  + corrX1X3*sigmaX3*gamma_exp
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
!     .                            ( 1.0 - corrX1X2**2.0 )
!     .                           *(sigmaX2**2.0)*(beta_exp**2.0)
!     .                         +  ( 1.0 - corrX1X3**2.0 )
!     .                           *(sigmaX3**2.0)*(gamma_exp**2.0)
!     .                         + 2.0*( corrX2X3 - corrX1X2*corrX1X3 )
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
!        ! Gaussian.
!        !
!        ! P_1(X1,X2) and P_2(X1,X2) are simply the equation for a
!        ! Bivariate Gaussian Distribution.  The only difference between the
!        ! two is the dependence on each individual plume for the plume average
!        ! and plume variance for the X1 variable.
!        !
!        ! The Probability Density Function for a Bivariate Normal Distribution:
!        !
!        ! P_i(X1,X2) = 
!        !
!        ! 1 / { 2 * PI * sigmaX1i * sigmaX2 * SQRT[ 1 - corrX1X2^2 ] }
!        ! * EXP{ -(1/2) * phi }
!        !
!        ! where,
!        !
!        ! phi =
!        !
!        ! 1 / [ 1 - corrX1X2^2 ]
!        ! * {
!        !       [ 1 / sigmaX1i^2 ] * ( X1i - muX1i )^2
!        !     + [ 1 / sigmaX2^2 ] * ( X2 - muX2 )^2
!        !     + [ 2*corrX1X2 / (sigmaX1i*sigmaX2) ]
!        !       * ( X1i - muX1i ) * ( X2 - muX2 )
!        !   }
!        !
!        ! Definitions:
!        !
!        ! muX1i:  The plume average of the X1 (or Y1) variable.
!        ! muX2:  The average of the X2 (or ln Y2) variable.
!        !
!        ! sigmaX1i:  The plume standard deviation (square-root of variance)
!        !            of the X1 (or Y1) variable.
!        ! sigmaX2:  The standard deviation (square-root of variance) of the
!        !           X2 (or ln Y2) variable.
!        !
!        ! corrX1X2:  Intra-Gaussian correlation between X1 and X2 (or between
!        !            Y1 and ln Y2).
!        !
!        ! After careful (and long) double-integration of the above equation,
!        ! one gets the result for each plume of the bivariate normal
!        ! -- log-normal equation.
!        !
!        ! The result for each plume is (Equation #1):
!        !
!        !   [ 1 / SQRT(2*PI) ] * ( sigmaX1i )^alpha 
!        ! * EXP{    muX2*beta + (1/2)*sigmaX2^2*beta^2
!        !        - (1/4) * [ ( muX1i / sigmaX1i ) + corrX1X2*sigmaX2*beta ]^2
!        !      }
!        ! * GAMMA_FNC { alpha + 1 }
!        ! * PARAB_CYL_FNC_[-(alpha+1)] { - [  ( muX1i / sigmaX1i )
!        !                                   + corrX1X2*sigmaX2*beta ]
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
!     .                             corrY1Y2, alpha_exp, beta_exp )
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
!        REAL, INTENT(IN):: corrY1Y2  ! Intra-Gaussian correlation between Y1,Y2
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
!        REAL:: corrX1X2  ! Intra-Gaussian correlation between X1,X2
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
!        ! corrY1Y2 is a correlation between a Gaussian and a Log-normal.
!        ! It must be converted to corrX1X2, which is a correlation between
!        ! two Gaussians.
!        corrX1X2 = corrY1Y2 
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
!     .                              + corrX1X2*sigmaX2*beta_exp )
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
!     .                         + corrX1X2*sigmaX2*beta_exp )**2.0
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
        use parabolic, only: gamma, parab
        use constants, only: pi_dp

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
     .  argument ! Order 'x' of Dv(a,x)         [-]

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
