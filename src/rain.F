        
! This subroutine advances drizzle (rrm and Nrm) one timestep.
! It uses implicit discretization and formulas by 
!     Khairoutdinov and Kogan (2000)

        SUBROUTINE rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .                   rcm, rtm, rrm, Nrm, cond_rrm, auto_rrm,
     .                   accr_rrm, cond_Nrm, auto_Nrm, sed_rrm, 
     .                   diff_rrm, subs_rrm, Vr, mean_vol_rad,
! LH Kessler autoconversion.  Vince Larson 22 May 2005:
     .                   AKm_est, AKm, 
     .                   runtype )


        USE grid_class
        USE constants

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
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm    ! cloud water mixing ratio
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rtm    ! total water mixing ratio
        INTEGER, INTENT(IN):: runtype

        ! Input/Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: rrm ! rain water mixing ratio
        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: Nrm ! rain droplet number conc.

        ! Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: cond_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: auto_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: accr_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: cond_Nrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: auto_Nrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: sed_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: diff_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: subs_rrm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Vr
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: mean_vol_rad

        ! Parameter
!        LOGICAL, PARAMETER:: local_kk = .false.  ! in constants.F

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: rvm          ! water vapor mixing ratio; kg/kg
        REAL, DIMENSION(1:gr%nnzp):: conc_cloud   ! concentration of cloud drops; 
                                                  !   number/m^3
        REAL, DIMENSION(1:gr%nnzp):: Temp         ! Temperature; K
        REAL, DIMENSION(1:gr%nnzp):: es           ! Saturation vapor pressure; Pa
        REAL, DIMENSION(1:gr%nnzp):: Supsat       ! Supersaturation
        REAL, DIMENSION(1:gr%nnzp):: VNr          ! Number conc. sedimentation 
                                                  !   velocity; m/s 
        REAL, DIMENSION(1:gr%nnzp):: aa, bb, cc, dd, oo, pp, qq, rr
! Vince Larson made rs into an array.  20 April 2005
        REAL, DIMENSION(1:gr%nnzp):: rs        ! Saturation mixing ratio; kg/kg



        REAL:: e             ! Vapor pressure; Pa
        REAL:: thm           ! Potential temperature; K
!        REAL:: rs            ! Saturation mixing ratio; kg/kg
! End Vince Larson's change

! Vince Larson added a new eddy diffusivity coeff for rain.  20 Apr 2005
        REAL:: nu_r
! End Vince Larson's addition

        REAL:: a, thl1, thl2, rc1, rc2, s1, s2, ss1, ss2

        INTEGER:: k
        INTEGER:: km1, kp1

        REAL :: Beta_T


        ! Function Calls
        REAL:: cond_evap_rrm, autoconv_rrm, accretion_rrm
        REAL:: cond_evap_Nrm, autoconv_Nrm

! Logical flag: LH_on = .true. if Latin hypercube sampling is on.
! Vince Larson 22 May 2005
! Now occurs in constants.F
!        LOGICAL, INTENT(IN) :: LH_on
! Latin hypercube estimate of Kessler autoconversion   Vince Larson 22 May 2005
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN) :: AKm_est
! Analytic calculation of Kessler autoconversion       Vince Larson 22 May 2005
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN) :: AKm


! DESCRIPTION
!
! This subroutine computes the new rain water mixing ratio.  It computes it
! according to the equation:
!
! drr/dt = -w*(drr/dz) + Vr*(drr/dz) + nu2*(d2rr/dz2)
!                      + (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! where Vr is the mean sedimentation velocity of the rain drops.  It
! is determined by Vr = 0.012*rvr - 0.2, with Vr in m/s and rvr, the
! rain drop mean volume radius, in um.  In turn the rain drop mean
! volume radius is given by the following equation:
!
! rvr = [(4*PI*rho_lw)/(3*rho)]^(-1/3) * rr^(1/3) * Nr^(-1/3);
!
! where Nr is the rain drop concentration (num/m^3).  rvr is given in
! meters.
!
! It is important to note that a one-sided upwind advection scheme is used
! for sedimentation velocity, with Vr being from the level above the level
! which is being dealt with.
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
! drr/dt + w*(drr/dz) - Vr*(drr/dz) - nu2*(d2rr/dz2)
!                      = (drr/dt)cond + (drr/dt)auto + (drr/dt)accr
!
! The equation is then broken down as following:
! (n is the current timestep; n+1 is the next timestep; 
! k is the vertical level).
!
! rr(n+1,k)/dt
! + w(n,k)*(rr(n+1,k+1)-rr(n+1,k-1))/(del_zm(n+1,k)-del_zm(n+1,k-1))
! - Vr(n,k+1)*(rr(n+1,k+1)-rr(n+1,k))/(del_zm(n+z,k))
! + nu2*(1/del_zt(n+1,k)) 
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

        DO k = 1, gr%nnzp-1, 1

           ! Find the important conditions

           rvm(k) = rtm(k) - rcm(k)   ! rvm is water vapor mixing ratio.

           ! Find current vapor pressure.
           e = (p(k)*rvm(k))/(ep + rvm(k))

           ! Concentration of cloud drops at this altitude.
           IF (rcm(k) > 0.0) THEN
              IF ( runtype == 0 ) THEN
                 conc_cloud(k) = 50000000.0  ! Average maritime distribution.
!                 conc_cloud(k) = 150000000.0 ! Average continental distribution.
              ELSEIF ( runtype == 8 ) THEN
                 conc_cloud(k) = 25000000.0  ! ASTEX A209 case (peak value).
              ELSEIF ( runtype == 9 ) THEN
                 conc_cloud(k) = 65000000.0  ! DYCOMS RF02 case.
              ELSE
                 conc_cloud(k) = 50000000.0  ! Average maritime distribution.
!                 conc_cloud(k) = 150000000.0 ! Average continental distribution.
              ENDIF
           ELSE
              conc_cloud(k) = 0.0
           ENDIF

           ! Get theta from theta-l
           thm = thlm(k) + (Lv/(Cp*exner(k)))*rcm(k)

! Vince Larson made rs into an array.  20 April 2005
!           ! Saturation mixing ratio
!           rs = rsat(p(k),thm*exner(k))
!           ! Find temperature and saturation vapor pressure.
!           Temp(k) = thm/((p0/p(k))**kappa)
!           es(k) = (p(k)*rs)/(ep + rs)

           ! Saturation mixing ratio
           rs(k) = rsat(p(k),thm*exner(k))
           ! Find temperature and saturation vapor pressure.
           Temp(k) = thm/((p0/p(k))**kappa)
           es(k) = (p(k)*rs(k))/(ep + rs(k))
! End Vince Larson's change

           ! Supersaturation
!           Supsat(k) = (e/es(k)) - 1.0

           CALL sed_velocity( rhot(k), Nrm(k), rrm(k),
     .                        Vr(k), VNr(k), mean_vol_rad(k) )

        ENDDO

        ! This is not looped over the top level because the pdf
        ! parameters are not looped over the top level.

        ! Eddy diffusivity for turbulent transport of 
        !     rain mix ratio and # conc.
        ! Use nu_r=10 m^2/s for 10-m vertical grid spacing.
        ! Use nu_r=70 m^2/s for 100-m vertical grid spacing.
        nu_r = 10.0  ! 10-m grid spacing
!        nu_r = 70.0  ! 100-m grid spacing

        DO k = 2, gr%nnzp-1, 1

           km1 = MAX(k-1,1)
           kp1 = MIN(k+1,gr%nnzp)

           ! LHS terms for rain water mixing ratio, rrm.

           ! Find subdiagonal aa

           aa(k) = 
     .            - wmt(k) / (1.0/gr%dzm(k) + 1.0/gr%dzm(km1))   ! advection
     .            - nu_r*gr%dzt(k)*gr%dzm(km1)                  ! diffusion
     
          ! Find main diagonal bb

           bb(k) = 
     .              (1.0/dt)                                     ! time tendency
     .            + Vr(kp1)*gr%dzm(k)                            ! sedimentation
     .            + nu_r*gr%dzt(k)*(gr%dzm(k) + gr%dzm(km1))    ! diffusion

           ! Find superdiagonal cc

           cc(k) = 
     .              wmt(k) / (1.0/gr%dzm(k) + 1.0/gr%dzm(km1))   ! advection
     .            - Vr(kp1)*gr%dzm(k)                            ! sedimentation
     .            - nu_r*gr%dzt(k)*gr%dzm(k)                    ! diffusion

           ! LHS terms for rain droplet number concentration, Nrm.

           ! Find subdiagonal oo

           oo(k) = 
     .            - wmt(k) / (1.0/gr%dzm(k) + 1.0/gr%dzm(km1))   ! advection
     .            - nu_r*gr%dzt(k)*gr%dzm(km1)                    ! diffusion

           ! Find main diagonal pp

           pp(k) = 
     .              (1.0/dt)                                     ! time tendency
     .            + VNr(kp1)*gr%dzm(k)                           ! sedimentation
     .            + nu_r*gr%dzt(k)*(gr%dzm(k) + gr%dzm(km1))      ! diffusion

           ! Find superdiagonal qq

           qq(k) = 
     .              wmt(k) / (1.0/gr%dzm(k) + 1.0/gr%dzm(km1))   ! advection
     .            - VNr(kp1)*gr%dzm(k)                           ! sedimentation
     .            - nu_r*gr%dzt(k)*gr%dzm(k)                      ! diffusion

           ! PDF parameters from pdf_closure_new for use in

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
           Beta_T = (R/Rv) * ( Lv/(R*Temp(k)) ) * ( Lv/(Cp*Temp(k)) )

           Supsat(k) = (a*s1+(1-a)*s2)*((1.0 + Beta_T*rs(k))/rs(k))

           ! Now find the elements that make up the right-hand side of the
           ! equation, dd, for rain water mixing ratio, rrm.

           cond_rrm(k) = cond_evap_rrm( Temp(k), p(k), rhot(k),
     .                     a, thl1, thl2, rc1, rc2, s1, s2,
     .                     ss1, ss2, es(k), Supsat(k), Nrm(k), rrm(k) )

! Vince Larson added option to call LH sampled Kessler autoconversion.
! 22 May 2005
!           auto_rrm(k) = autoconv_rrm( rcm(k), conc_cloud(k), rhot(k),
!     .                        a, s1, s2, ss1, ss2, rc1, rc2 )
           if ( LH_on ) then

!              auto_rrm(k) = AKm_est(k)
              auto_rrm(k) = AKm(k)

           else

              auto_rrm(k) = autoconv_rrm( rcm(k), conc_cloud(k), 
     .               rhot(k), a, s1, s2, ss1, ss2, rc1, rc2 )

           endif
! End Vince Larson's addition

           accr_rrm(k) = accretion_rrm( rcm(k), rrm(k),
     .                             a, s1, s2, ss1, ss2, rc1, rc2 )

           ! Enter the entire dd array, which is the RHS of the equation.
           dd(k) = (1.0/dt)*rrm(k)
     .                     + cond_rrm(k) + auto_rrm(k) + accr_rrm(k)

           ! Now find the elements that make up the right-hand side of the
           ! equation, rr, for rain droplet number concentration, Nrm.

           cond_Nrm(k) = cond_evap_Nrm( cond_rrm(k), Nrm(k), rrm(k) )

           auto_Nrm(k) = autoconv_Nrm( auto_rrm(k), rhot(k) )

           ! Enter the entire rr array, which is the RHS of the equation.
           rr(k) = (1.0/dt)*Nrm(k) + cond_Nrm(k) + auto_Nrm(k)

        ENDDO

!       set up boundary conditions

        ! rrm
        aa(1) = aa(2)
        aa(gr%nnzp) = aa(gr%nnzp-1)
        bb(1) = bb(2)
        bb(gr%nnzp) = bb(gr%nnzp-1)
        cc(1) = cc(2)
        cc(gr%nnzp) = cc(gr%nnzp-1)
        dd(1) = dd(2)
        dd(gr%nnzp) = dd(gr%nnzp-1)

        ! Nrm
        oo(1) = oo(2)
        oo(gr%nnzp) = oo(gr%nnzp-1)
        pp(1) = pp(2)
        pp(gr%nnzp) = pp(gr%nnzp-1)
        qq(1) = qq(2)
        qq(gr%nnzp) = qq(gr%nnzp-1)
        rr(1) = rr(2)
        rr(gr%nnzp) = rr(gr%nnzp-1)


!       solve tridiagonal systems

        call tridag(aa, bb, cc, dd, rrm, gr%nnzp)

        call tridag(oo, pp, qq, rr, Nrm, gr%nnzp)


!       boundary conditions on results

        ! rrm
        rrm(1) = rrm(2)
        rrm(gr%nnzp) = rrm(gr%nnzp-1)

        ! Nrm
        Nrm(1) = Nrm(2)
        Nrm(gr%nnzp) = Nrm(gr%nnzp-1)

        DO k = 1, gr%nnzp, 1

           ! Value of rain water mixing ratio cannot fall below 0.
           IF (rrm(k) < 0.0) THEN
              rrm(k) = 0.0
           ENDIF

           ! Value of rain droplet number concentration cannot fall below 0.
           IF (Nrm(k) < 0.0) THEN
              Nrm(k) = 0.0
           ENDIF

           ! Breakdown for output purposes.

           km1 = MAX(k-1,1)
           kp1 = MIN(k+1,gr%nnzp)

           sed_rrm(k) = Vr(kp1) * ( (rrm(kp1) - rrm(k)) * gr%dzm(k) )

           diff_rrm(k) = nu2 * gr%dzt(k)
     .                  * (  (rrm(kp1) - rrm(k))*gr%dzm(k)
     .                     - (rrm(k) - rrm(km1))*gr%dzm(km1)  )

           subs_rrm(k) = -wmt(k) * ( (rrm(kp1) - rrm(km1))
     .                            / (1.0/gr%dzm(k) + 1.0/gr%dzm(km1)) )


        ENDDO

        END SUBROUTINE rain

!------------------------------------------------------------------------

        SUBROUTINE sed_velocity( rhot, Nrm, rrm, 
     .                           Vr, VNr, mean_vol_rad )

        USE grid_class
        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: rhot, Nrm, rrm
!        LOGICAL, INTENT(IN):: local_kk
        REAL, INTENT(OUT):: Vr, VNr, mean_vol_rad

!       Added to constants.F
!        REAL, PARAMETER:: rho_lw = 1000.0 ! density of liquid water; kg/m^3

        ! Lognormal Distribution Variables.
        REAL:: Nrm_kg
        REAL:: rrp2_rrm2, Nrp2_Nrm2
        REAL:: u_rr_G, sig2_rr_G, u_Nr_G, sig2_Nr_G
        REAL:: corr_rrNr_G, bivar_term

        ! Vince Larson: 0.75 agrees with COAMPS DYCOMS2 RF02
        !     in cloud and below cloud.  7 Apr 2005
        REAL, PARAMETER:: corr_rrNr_L = 0.75 

        ! Find sedimentation velocity, Vr.

        ! First, find mean volume radius of rain drops (in m).
        IF (rrm > 0.0 .AND. Nrm > 0.0) THEN

           IF ( local_kk ) THEN
        
              mean_vol_rad = ( ( (4.0*pi*rho_lw)
     .                      /(3.0*rhot) )**(-1.0/3.0) )
     .                      *(rrm**(1.0/3.0))
     .                      *(Nrm**(-1.0/3.0))

           ELSEIF ( .not. local_kk ) THEN

              ! The following code is to find the grid-box average mean
              ! volume radius for a drizzle drop based on a lognormal 
              ! distribution.  This takes the original Khairoutdinov and
              ! Kogan code and turns it into a calculation of grid-box average
              ! based on a bivariate lognormal distribution.  This is done
              ! because the orignal code works for a local formula for 3-D 
              ! models with high resolution, but this code is for a 1-D model.
              ! The equation was produced by Vince Larson (2005).

              ! Note:  In order to make this equation able to be turned into
              ! a bivariate lognormal distribution, the first term must be 
              ! constant.  Therefore, rho is taken out of it.  In order to make
              ! up for that, Nr is divided by rho to become concentration per kg.
              !
              ! rvr = [(4*pi*rho_lw)/3]^(-1/3) * rr^(1/3) * Nr^(-1/3)
              !
              ! For our purposes, we need the mean rvr:
              !
              ! mean(rvr) = [(4*pi*rho_lw)/3]^(-1/3) * mean(rr^(1/3) * Nr^(-1/3))

              Nrm_kg = Nrm/rhot

              ! rrp2_rrm2=1 agrees with COAMPS DYCOMS2 RF02 in and below cloud.  
              !     7 Apr 2005
!              rrp2_rrm2 = 2 * rcp2 / MAX(rcm**2.0, 1.0E-10)
              rrp2_rrm2 = 1.0

              Nrp2_Nrm2 = 0.75

              u_rr_G = LOG( rrm 
     .                     * ((1.0 + rrp2_rrm2)**(-1.0/2.0)) )
              sig2_rr_G = LOG(1.0 + rrp2_rrm2)

              u_Nr_G = LOG( Nrm_kg 
     .                     * ((1.0 + Nrp2_Nrm2)**(-1.0/2.0)) )
              sig2_Nr_G = LOG(1.0 + Nrp2_Nrm2)

              ! Correlation coefficient for rr and Nr.
              corr_rrNr_G = (  LOG( 1.0 + corr_rrNr_L
     .                       * ((EXP(sig2_rr_G) - 1.0)**(1.0/2.0))
     .                       * ((EXP(sig2_Nr_G) - 1.0)**(1.0/2.0)) )  )
     .                     / ( SQRT(sig2_rr_G)*SQRT(sig2_Nr_G) )

              ! This term is:  mean(rr^(1/3) * Nr^(-1/3))
              bivar_term = EXP(
     .                       (1.0/3.0)*u_rr_G + (-1.0/3.0)*u_Nr_G
     .                     + (1.0/2.0) * ((1.0/3.0)**2.0) * sig2_rr_G
     .                     + (1.0/2.0) * ((-1.0/3.0)**2.0) * sig2_Nr_G
     .                     + corr_rrNr_G
     .                     * (1.0/3.0) * SQRT(sig2_rr_G)
     .                     * (-1.0/3.0) * SQRT(sig2_Nr_G)        )

              ! The mean volume radius
              mean_vol_rad = ( ((4.0*pi*rho_lw)/3.0)**(-1.0/3.0) )
     .                      * bivar_term

           ENDIF

        ELSE

           ! If there isn't any rain, then there isn't any drop rad.
           mean_vol_rad = 0.0

        ENDIF
        
        ! Khairoutdinov and Kogan Sedimentation Velocity Calculation.
        ! sedimentation velocity of rain drops (in m/s).
        ! mean volume radius converted to um.
        ! Note 1: positive sedimentation velocity means downwards.
        ! Note 2: Vr(m/s) = 0.012*rvr(um) - 0.2
        ! If rvr(um) is too small, the equation will result in an
        ! upwards (negative) sedimentation velocity.  The limiter on
        ! the first line is put in there for that reason.

        IF (1000000.0*mean_vol_rad > (0.2/0.012) ) THEN
           Vr = 0.012*(1000000.0*mean_vol_rad) - 0.2
        ELSE
           Vr = 0.0
        ENDIF

        IF (1000000.0*mean_vol_rad > (0.1/0.007) ) THEN
           VNr = 0.007*(1000000.0*mean_vol_rad) - 0.1
        ELSE
           VNr = 0.0
        ENDIF

        END SUBROUTINE sed_velocity

!------------------------------------------------------------------------

        FUNCTION cond_evap_rrm(Temp, press, rhot, a, 
     .                         thl1, thl2, rc1, rc2, s1, s2, ss1, ss2,
     .                         es, Supsat, Nrm, rrm)

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: Temp, press, rhot
        REAL, INTENT(IN):: a, thl1, thl2, rc1, rc2, s1, s2, ss1, ss2
        REAL, INTENT(IN):: es, Supsat, Nrm, rrm
!        LOGICAL, INTENT(IN):: local_kk
        REAL:: cond_evap_rrm

        REAL:: Nrm_kg
        REAL:: rrp2_rrm2, Nrp2_Nrm2
        REAL:: u_rr_G, u_Nr_G
        REAL:: sig2_rr_G, sig2_Nr_G
        REAL:: corr_rcrr_G, corr_rcNr_G, corr_rrNr_G
!       Now in constants.F
!        REAL, PARAMETER:: rho_lw = 1000.0 ! density of liquid water; kg/m^3
        REAL, PARAMETER:: C_evap = 0.86*0.2 ! COAMPS value of
                                            ! ratio of drizzle drop mean 
                                            ! geometric radius to drizzle drop
                                            ! mean volume radius. Khairoutdinov 
                                            ! and Kogan (2000); p. 233.
!        REAL, PARAMETER:: C_evap = 0.86    ! KK 2000, Eq. 22 value. 

        ! corr_rcrr_L=0.1 agrees w/ COAMPS DYCOMS2 RF02 below cloud.  
        REAL, PARAMETER:: corr_rcrr_L = 0.1
        ! corr_rcNr_L=-0.15 agrees w/ COAMPS DYCOMS2 RF02 below cloud.  
        REAL, PARAMETER:: corr_rcNr_L = -0.15
        ! corr_rrNr_L=0.85 agrees w/ COAMPS DYCOMS2 RF02 below cloud.  
        REAL, PARAMETER:: corr_rrNr_L = 0.85

        ! Function call
        REAL:: G_t_p
        REAL:: NLL_PDF


        IF (rrm > 0.0 .AND. Nrm > 0.0 ) THEN

           IF (local_kk ) THEN

              IF ( Supsat < 0.0 ) THEN

                 ! change in rain mixing ratio over time due to condensation.
                 cond_evap_rrm = 3.0 * C_evap * G_t_p(Temp, press)
     .                   * ( ((4.0*pi*rho_lw)/(3.0*rhot))**(2.0/3.0) )
     .                   * (rrm**(1.0/3.0)) * (Nrm**(2.0/3.0))
     .                   * Supsat

              ELSE  ! layer contains cloud

                 cond_evap_rrm = 0.0

              ENDIF   

           ELSEIF (.not. local_kk ) THEN

              ! Convert to number per unit mass
              Nrm_kg = Nrm/rhot

              ! rrp2_rrm2=1.2 below cloud acc. to COAMPS DYCOMS2 RF02.
!              rrp2_rrm2 = 2 * rcp2 / MAX(rcm**2.0, 1.0E-10)
              rrp2_rrm2 = 1.2

              ! Nrp2_Nrm2=1.0 below cloud acc. to COAMPS DYCOMS2 RF02.
              Nrp2_Nrm2 = 1.0

              u_rr_G = LOG( rrm 
     .                     * ((1.0 + rrp2_rrm2)**(-1.0/2.0)) )
              sig2_rr_G = LOG(1.0 + rrp2_rrm2)

              u_Nr_G = LOG( Nrm_kg 
     .                     * ((1.0 + Nrp2_Nrm2)**(-1.0/2.0)) )
              sig2_Nr_G = LOG(1.0 + Nrp2_Nrm2)

              ! Correlation coefficient for rc and rr.
              corr_rcrr_G = corr_rcrr_L
     .                     * ((EXP(sig2_rr_G) - 1.0)**(1.0/2.0)) 
     .                     / SQRT(sig2_rr_G)

              ! Correlation coefficient for rc and Nr.
              corr_rcNr_G = corr_rcNr_L
     .                     * ((EXP(sig2_Nr_G) - 1.0)**(1.0/2.0)) 
     .                     / SQRT(sig2_Nr_G)

              ! Correlation coefficient for rr and Nr.
              corr_rrNr_G = (  LOG( 1.0 + corr_rrNr_L
     .                       * ((EXP(sig2_rr_G) - 1.0)**(1.0/2.0))
     .                       * ((EXP(sig2_Nr_G) - 1.0)**(1.0/2.0)) )  )
     .                     / ( SQRT(sig2_rr_G)*SQRT(sig2_Nr_G) )

              cond_evap_rrm =    a    * NLL_PDF(
     .                          thl1, press, rc1, s1, ss1,
     .                          1.0, 1.0/3.0, 2.0/3.0, 
     .                          u_rr_G, u_Nr_G, sig2_rr_G, sig2_Nr_G,
     .                          corr_rcrr_G, corr_rcNr_G, corr_rrNr_G )
     .                      + (1.0-a) * NLL_PDF(
     .                          thl2, press, rc2, s2, ss2,
     .                          1.0, 1.0/3.0, 2.0/3.0, 
     .                          u_rr_G, u_Nr_G, sig2_rr_G, sig2_Nr_G,
     .                          corr_rcrr_G, corr_rcNr_G, corr_rrNr_G )


           ENDIF

        ELSE  ! Nrm or rrm is zero, i.e. there is no rain

           cond_evap_rrm = 0.0

        ENDIF

        END FUNCTION cond_evap_rrm

!------------------------------------------------------------------------

        FUNCTION cond_evap_Nrm(cond_rrm, Nrm, rrm)

        IMPLICIT NONE

        REAL, INTENT(IN):: cond_rrm, Nrm, rrm
        REAL:: cond_evap_Nrm

        IF (rrm > 0.0 .AND. Nrm > 0.0) THEN

           cond_evap_Nrm = ( Nrm / rrm ) * cond_rrm

        ELSE

           cond_evap_Nrm = 0.0

        ENDIF

        END FUNCTION cond_evap_Nrm

!------------------------------------------------------------------------

        FUNCTION NLL_PDF( thl_pl, press, rc_pl, s, ss,
     .                    alpha_exp, beta_exp, gamma_exp, 
     .                    u_rr_G, u_Nr_G, sig2_rr_G, sig2_Nr_G,
     .                    corr_rcrr_G, corr_rcNr_G, corr_rrNr_G )

        USE constants

        IMPLICIT NONE

        ! External functions

        real rsat
        external rsat

        REAL, INTENT(IN):: thl_pl, press, rc_pl, s, ss
        REAL, INTENT(IN):: alpha_exp, beta_exp, gamma_exp
        REAL, INTENT(IN):: u_rr_G, u_Nr_G, sig2_rr_G, sig2_Nr_G
        REAL, INTENT(IN):: corr_rcrr_G, corr_rcNr_G, corr_rrNr_G
        REAL:: NLL_PDF

        REAL:: tmpl_pl, tmp_pl, rsl_pl, Beta_T

        DOUBLE PRECISION:: scc, gamma_in, gamma_out, Dv_order, test

!       Now in constants.F
!        REAL, PARAMETER:: rho_lw = 1000.0 ! density of liquid water; kg/m^3
        REAL, PARAMETER:: C_evap = 0.86*0.2 ! COAMPS value of
                                            ! ratio of drizzle drop mean 
                                            ! geometric radius to drizzle drop
                                            ! mean volume radius. Khairoutdinov 
                                            ! and Kogan (2000); p. 233.
!        REAL, PARAMETER:: C_evap = 0.86    ! KK 2000, Eq. 22 value. 
!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308.0
!       Found that above is not valid on most compilers -dschanen
        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308

        ! Make function D_fnc results a Double Precision number.
        DOUBLE PRECISION:: D_fnc

        ! Function call
        REAL:: G_t_p

        tmpl_pl =  thl_pl * ((press/p0)**kappa)

        tmp_pl = tmpl_pl + (Lv/Cp)*rc_pl

        rsl_pl = rsat(press, tmpl_pl)

        Beta_T = (R/Rv) * ( Lv/(R*tmpl_pl) ) * ( Lv/(Cp*tmpl_pl) )

        ! The following IF statement is for a normal situation where
        ! there is variance in liquid water, meaning that ss > 0.
        ! This is Larson's eq. number 46.


        ! Vince Larson skipped IF to simplify subroutine.  3 May 2005.  
        ! Then only the final ELSEIF and below is called.
!        IF (ss > 0.0) THEN
        IF (-1.0 > 0.0) THEN

           scc = (s/ss) + corr_rcrr_G*beta_exp*SQRT(sig2_rr_G)
     .                  + corr_rcNr_G*gamma_exp*SQRT(sig2_Nr_G)

           Dv_order = -alpha_exp - 1.0

           test = D_fnc( Dv_order, scc )


           ! The following if statement makes sure that the result
           ! of the Parabolic Cylinder function falls within the range
           ! of a Double Precision number for the given alpha and scc.

           IF ( test >= 0.0d0 .AND. test < limit ) THEN

              gamma_in = alpha_exp + 1.0

              CALL GAMMA( gamma_in, gamma_out )

              NLL_PDF = 
                 ! Factor 1
     .           3.0 * C_evap * G_t_p(tmp_pl, press)
     .           * ( ((4.0*pi*rho_lw)/3.0)**(2.0/3.0) )
     .           * ( (1.0 + Beta_T*rsl_pl)/rsl_pl )
                 ! Factor 2
     .           * (1.0/SQRT(2.0*pi))*(-ss**alpha_exp)
     .           * EXP(beta_exp*u_rr_G + gamma_exp*u_Nr_G)
                 ! Factor 3
     .           * EXP( (1.0/2.0) * (
     .              (1.0 - corr_rcrr_G**2.0)*(beta_exp**2.0)*sig2_rr_G
     .            + (1.0 - corr_rcNr_G**2.0)*(gamma_exp**2.0)*sig2_Nr_G
     .            + 2.0*(corr_rrNr_G - corr_rcNr_G*corr_rcrr_G)
     .              *beta_exp*SQRT(sig2_rr_G)*gamma_exp*SQRT(sig2_Nr_G)
     .                              )     )
                 ! Factor 4 
     .            * EXP( (1.0/4.0)*scc**2.0 - (1.0/2.0)*(s/ss)*scc )
     .            * gamma_out * D_fnc( Dv_order, scc )


           ! The following ELSEIF statement is for a situation where there
           ! is such small variance (ss) that the result of the Parabolic
           ! Cylinder function (D_fnc) is higher than can be represented
           ! by a Double Precision number.  Therefore, it is treated as a 
           ! spike with ss = 0.  This is Larson's eq. number 45.
           ! Also, for this particular situation, s must be less than 0, 
           ! meaning that the air is subsaturated and evaporation occurs.
           ! If this is not the case, then the air would be saturated or 
           ! supersaturated and condensation would occur.  However, the HOC
           ! code requires that all supersaturation be removed from the air
           ! and made into cloud water.  None of it is left for rain water
           ! condensation and therefore we must make this factor 0 if 
           ! a saturation or a supersaturation is present. 

           ELSEIF (s < 0.0) THEN

              NLL_PDF = 
                 ! Factor 1
     .           3.0 * C_evap * G_t_p(tmp_pl, press)
     .           * ( ((4.0*pi*rho_lw)/3.0)**(2.0/3.0) )
     .           * ( (1.0 + Beta_T*rsl_pl)/rsl_pl )
                 ! Factor 2
     .           * (s**alpha_exp)
     .           * EXP(beta_exp*u_rr_G + gamma_exp*u_Nr_G)
                 ! Factor 3
     .           * EXP( (1.0/2.0) * (
     .              (1.0 - corr_rcrr_G**2.0)*(beta_exp**2.0)*sig2_rr_G
     .            + (1.0 - corr_rcNr_G**2.0)*(gamma_exp**2.0)*sig2_Nr_G
     .            + 2.0*(corr_rrNr_G - corr_rcNr_G*corr_rcrr_G)
     .              *beta_exp*SQRT(sig2_rr_G)*gamma_exp*SQRT(sig2_Nr_G)
     .                              )     )


           ! The following ELSE statement is for a situation where 
           ! s >= 0 and ss = 0, which means no clear air.  Therefore,
           ! no evaporation can occur.

           ELSE

              NLL_PDF = 0.0

           ENDIF


        ! The following ELSEIF statement is for a situation where there
        ! is no variance (ss) so that the result of the Parabolic
        ! Cylinder function (D_fnc) is infinite.  Therefore, it is a 
        ! spike with ss = 0.  This is Larson's eq. number 45.
        ! Also, for this particular situation, s must be less than 0, 
        ! meaning that the air is subsaturated and evaporation occurs.
        ! If this is not the case, then the air would be saturated or 
        ! supersaturated and condensation would occur.  However, the HOC
        ! code requires that all supersaturation be removed from the air
        ! and made into cloud water.  None of it is left for rain water
        ! condensation and therefore we must make this factor 0 if 
        ! a saturation or a supersaturation is present. 

        ELSEIF (s < 0.0) THEN

           NLL_PDF = 
              ! Factor 1
     .        3.0 * C_evap * G_t_p(tmp_pl, press)
     .        * ( ((4.0*pi*rho_lw)/3.0)**(2.0/3.0) )
     .        * ( (1.0 + Beta_T*rsl_pl)/rsl_pl )
              ! Factor 2
     .        * (s**alpha_exp)
     .        * EXP(beta_exp*u_rr_G + gamma_exp*u_Nr_G)
              ! Factor 3
     .        * EXP( (1.0/2.0) * (
     .           (1.0 - corr_rcrr_G**2.0)*(beta_exp**2.0)*sig2_rr_G
     .         + (1.0 - corr_rcNr_G**2.0)*(gamma_exp**2.0)*sig2_Nr_G
     .         + 2.0*(corr_rrNr_G - corr_rcNr_G*corr_rcrr_G)
     .           *beta_exp*SQRT(sig2_rr_G)*gamma_exp*SQRT(sig2_Nr_G)
     .                           )     )

        ! The following ELSE statement is for a situation where 
        ! s >= 0 and ss = 0, which means no clear air.  Therefore,
        ! no evaporation can occur.

        ELSE   ! No clear air

           NLL_PDF = 0.0

        ENDIF

        END FUNCTION NLL_PDF

!------------------------------------------------------------------------

        FUNCTION G_t_p(Temp, press)

        USE constants

        IMPLICIT NONE


        ! Here we compute G(T,p) as in KK (17) 
        !      and Eq. (7.17) of Rogers and Yau 1989

        ! External functions

        real rsat
        external rsat

        REAL, INTENT(IN):: Temp, press
        REAL:: G_t_p

        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: es, rs, Celsius
!       Now in constants.F
!        REAL, PARAMETER:: rho_lw = 1000.0 ! density of liquid water; kg/m^3

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

        G_t_p = 1.0/(Fk + Fd)

        END FUNCTION G_t_p

!------------------------------------------------------------------------

        FUNCTION autoconv_rrm( rcm, conc_cloud, rhot,
     .                           a, s1, s2, ss1, ss2, rc1, rc2 )
        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: rcm, conc_cloud, rhot
!        LOGICAL, INTENT(IN):: local_kk
        REAL, INTENT(IN):: a, s1, s2, ss1, ss2, rc1, rc2
        REAL:: autoconv_rrm

        ! Single Normal-Lognormal PDF Variables.
        REAL:: conc_cloud_kg
        REAL:: Ncp2_Ncm2
        REAL:: u_Nc_G, sig2_Nc_G
        REAL:: corr_rcNc_G

        ! corr_rcNc_L=0.45 agrees with DYCOMS II RF02 in cloud.
        REAL, PARAMETER:: corr_rcNc_L = 0.45

        ! Function call
        REAL:: NL_PDF

        ! Change in rain water mixing ratio over time due to autoconversion.

        IF (rcm > 0.0 .AND. conc_cloud > 0.0) THEN

           IF ( local_kk ) THEN

              ! Khairoutdinov and Kogan local formula
              ! cloud droplet concentration converted to num/cm^3.

              autoconv_rrm = 1350.0*(rcm**2.47)
     .                      *((conc_cloud/1000000.0)**(-1.79))

           ELSEIF ( .not. local_kk ) THEN

              ! Larson's formulation for subgrid variability.
              ! Note:  For this version the cloud droplet concentration
              ! stays as num/m^3 and Larson has changed the value of the
              ! coefficient.

              conc_cloud_kg = conc_cloud/rhot

              ! Ncp2_Ncm2=0.07 agrees with COAMPS DYCOMS2 RF02 (in cloud). 
              Ncp2_Ncm2 = 0.07

              u_Nc_G = LOG( conc_cloud_kg 
     .                     * ((1.0 + Ncp2_Ncm2)**(-1.0/2.0)) )
              sig2_Nc_G = LOG(1.0 + Ncp2_Ncm2)

              corr_rcNc_G = corr_rcNc_L
     .                     * ((EXP(sig2_Nc_G) - 1.0)**(1.0/2.0)) 
     .                     / SQRT(sig2_Nc_G)

              autoconv_rrm = 7.419E13 * (
     .                       a    * NL_PDF(s1, ss1, 2.47, u_Nc_G, 
     .                              sig2_Nc_G, -1.79, corr_rcNc_G, rc1)
     .                  + (1.0-a) * NL_PDF(s2, ss2, 2.47, u_Nc_G, 
     .                              sig2_Nc_G, -1.79, corr_rcNc_G, rc2)
     .                                  )

           ENDIF

        ELSE

           autoconv_rrm = 0.0

        ENDIF

        END FUNCTION autoconv_rrm

!------------------------------------------------------------------------

        FUNCTION autoconv_Nrm(auto_rrm, rhot)

        USE constants


        IMPLICIT NONE

        REAL, INTENT(IN):: auto_rrm, rhot
        REAL:: autoconv_Nrm

!        now in constants.F
!        REAL, PARAMETER:: rho_lw = 1000.0 ! density of liquid water; kg/m^3
! Vince Larson set r_0=28mum to agree with COAMPS-LES formula. 15 April 2005
!        REAL, PARAMETER:: r_0 = 25.0e-6   ! Assumed radius of all new drops; m.
        REAL, PARAMETER:: r_0 = 28.0e-6   ! Assumed radius of all new drops; m.
! End Vince Larson's change.

        autoconv_Nrm = auto_rrm /
     .                 ( ((4.0*pi*rho_lw)/(3.0*rhot)) * (r_0**3.0) )

        END FUNCTION autoconv_Nrm

!------------------------------------------------------------------------

        FUNCTION accretion_rrm( rcm, rrm,
     .                          a, s1, s2, ss1, ss2, rc1, rc2 )

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: rcm, rrm
!        LOGICAL, INTENT(IN):: local_kk
        REAL, INTENT(IN):: a, s1, s2, ss1, ss2, rc1, rc2

        REAL:: accretion_rrm

        ! Single Normal-Lognormal PDF Variables.
        REAL:: rrp2_rrm2
        REAL:: u_rr_G, sig2_rr_G
        REAL:: corr_rcrr_G

        ! corr_rcrr_L = 0.35 agrees with DYCOMS II RF02 in cloud.
        REAL, PARAMETER:: corr_rcrr_L = 0.35

        ! Function call
        REAL:: NL_PDF

        ! Change in rain water mixing ratio over time due to accretion.

        IF (rcm > 0.0 .AND. rrm > 0.0) THEN

           IF ( local_kk ) THEN

              ! Khairoutdinov and Kogan local formula
              accretion_rrm = 67.0*((rcm*rrm)**1.15)

           ELSEIF ( .not. local_kk ) THEN

              ! Larson's formulation for subgrid variability.

              ! rrp2_rrm2=0.4 agrees with COAMPS DYCOMS2 RF02 in cloud.  
!              rrp2_rrm2 = rcp2 / MAX(rcm**2.0, 1.0E-10)
              rrp2_rrm2 = 0.4

              u_rr_G = LOG( rrm 
     .                     * ((1.0 + rrp2_rrm2)**(-1.0/2.0)) )
              sig2_rr_G = LOG(1.0 + rrp2_rrm2)

              corr_rcrr_G = corr_rcrr_L
     .                     * ((EXP(sig2_rr_G) - 1.0)**(1.0/2.0)) 
     .                     / SQRT(sig2_rr_G)

              accretion_rrm = 67.0 * (
     .                       a    * NL_PDF(s1, ss1, 1.15, u_rr_G, 
     .                              sig2_rr_G, 1.15, corr_rcrr_G, rc1)
     .                  + (1.0-a) * NL_PDF(s2, ss2, 1.15, u_rr_G, 
     .                              sig2_rr_G, 1.15, corr_rcrr_G, rc2)
     .                               )

           ENDIF

        ELSE

           accretion_rrm = 0.0

        ENDIF

        END FUNCTION accretion_rrm

!------------------------------------------------------------------------

        FUNCTION NL_PDF( s, ss, alpha_exp, u, sig2, 
     .                   beta_exp, corr, rc_pl )

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: s, ss, alpha_exp, u, sig2
        REAL, INTENT(IN):: beta_exp, corr, rc_pl
        REAL:: NL_PDF

        DOUBLE PRECISION:: sc, gamma_in, gamma_out, Dv_order, test

!        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308.0d0
!       Found that above is not valid on most compilers -dschanen
        DOUBLE PRECISION, PARAMETER:: limit = 10.0d0**308

        ! Make function D_fnc results a Double Precision number.
        DOUBLE PRECISION:: D_fnc

        ! DESCRIPTION
        !
        ! This function solves the equation for AUTO_NL or ACCR_NL,
        ! which is:
        !
        ! ( (1/SQRT(2*pi))*(ss^alpha) )
        ! * EXP( beta*u + (1/2)*(beta^2)*sig2 - (1/4)*sc^2 )
        ! * GAMMA(alpha + 1.0) * D[order: -alpha - 1](-sc )
        !
        ! where: sc = (s/ss) + corr*beta*SQRT(sig2)
        !
        ! In a case where there is no variability in s throughout
        ! the entire grid box, ss (sigma s) = 0.  In that case, sc
        ! would be infinite.  This would give the parabolic cylinder
        ! function and the EXP term both infinite values, and the 1st
        ! term a 0 value, leading to an answer that could be anything.
        ! Therefore, we have a special equation for this condition.
        ! It is:
        !
        ! ( rc_pl**alpha ) * EXP( beta*u + (1.0/2.0)*(beta**2.0)*sig2 )
        !
        ! rc_pl is the liquid water mixing ratio for that specific plume.
        ! s and rc_pl are the exact same parameter, except s can be negative
        ! in order to represent a subsaturation.  In those cases, rc_pl is 0.
        !
        ! There is one more thing which further complicates this equation.
        ! When the ratio of s to ss is very large (low, but not quite 0 
        ! variability), the value of sc becomes large.  When the value of
        ! sc becomes large and positive, taking the value of the parabolic
        ! cylinder function of -sc results in a very large number.
        ! Mathematically, this number is offset the 1st and 2nd terms of the 
        ! equation.  Also mathematically, as the variability in s gets lower
        ! and lower, the value of the 1st equation (the one with the Parabolic
        ! Cylinder function) gets closer and closer to the value of the 2nd
        ! equation (the one without the Parabolic Cylinder function). However,
        ! the computer induces errors into these mathematical calculations due
        ! to the fact that there is a limit to the size of the values that can
        ! be represented.  If the value of the parabolic cylinder function
        ! gets above 10^308 for a Double Precision data type, it is 
        ! represented as infinite by the computer.  Therefore, anything
        ! that results in that value being two high must be treated as
        ! though it has constant variability.  This is very close to the
        ! actual value anyway, so no accuracy is lost.


        ! The following IF statement is for a case where there is variability
        ! in liquid water (ss > 0).  This is for Larson's eq. number 24 and 36.

        IF (ss > 0.0) THEN

           sc = (s/ss) + corr*beta_exp*SQRT(sig2)

           Dv_order = -alpha_exp - 1.0

           test = D_fnc( Dv_order, -sc )


           ! The following IF statement is to test whether the result of the
           ! Parabolic Cylinder function is within the limits of a 
           ! Double Precision number.  If it is not, the Parabolic Cylinder
           ! function will produce an infinite result.  If that is the case,
           ! a different equation (the equation for a spike) must be used.

           IF ( test >= 0.0d0 .AND. test < limit ) THEN

              gamma_in = alpha_exp + 1.0

              CALL GAMMA( gamma_in, gamma_out )

              NL_PDF = ( (1.0/SQRT(2.0*pi))*(ss**alpha_exp) )
     .                * EXP( beta_exp*u
     .                      + (1.0/2.0)*(beta_exp**2.0)*sig2 
     .                      - (1.0/4.0)*sc**2.0 )
     .                * gamma_out * D_fnc( Dv_order, -sc )


           ! The following ELSEIF statement is for a case where the variability
           ! of liquid water (ss) is so low that the Parabolic Cylinder function
           ! becomes infinite.  When this happens, the equation for a spike is 
           ! used.  If there is a spike and sc > 0, then the following equation
           ! is used.  If sc <= 0, then there is no cloud and NL_PDF = 0.  The
           ! following equation is Larson's eq. number 23 and 35.

           ELSEIF (sc > 0.0) THEN

              NL_PDF = ( rc_pl**alpha_exp )
     .                * EXP( beta_exp*u 
     .                      + (1.0/2.0)*(beta_exp**2.0)*sig2 )


           ! The following ELSE statement is for a case where ss = 0 and
           ! sc <= 0, which means only clear air.  Therefore, NL_PDF = 0.

           ELSE

              NL_PDF = 0.0

           ENDIF


        ! The following ELSEIF statement is for a case where the variability
        ! of liquid water (ss) is 0.  When this happens, the equation for a 
        ! spike is used.  If there is a spike and sc > 0, then the following
        ! equation is used.  Since sc = s/ss + ... and ss = 0, then sc is 
        ! either positive or negative infinity based on if s is positive or
        ! negative.  If s <= 0 (and therefore sc <= 0), then there is no 
        ! cloud and NL_PDF = 0.  The following equation is Larson's eq. number
        ! 23 and 35.

        ELSEIF (s > 0.0) THEN 

           NL_PDF = ( rc_pl**alpha_exp )
     .             * EXP( beta_exp*u
     .                   + (1.0/2.0)*(beta_exp**2.0)*sig2 )

        ELSE

           NL_PDF = 0.0

        ENDIF

        END FUNCTION NL_PDF

!------------------------------------------------------------------------

        ! DESCRIPTION
        !
        ! The following function is a function which computes the Parabolic
        ! Cylinder function.  It's inputs are V (the order of the function)
        ! and X (the value at which the function is evaluated at).  Both
        ! inputs HAVE to be declared as Double Precision numbers in any code
        ! that calls this particular function.  The output of this function
        ! is also a Double Precision number.  As part of the computing the 
        ! Parabolic Cylinder function, there is a Subroutine GAMMA which
        ! computes the Gamma function.  It can be called from a separate
        ! area of code, as well.  In the subroutine call, the first argument
        ! is the value at which the function is evaluated at.  It must be
        ! declared as a Double Precision number in the code that calls that
        ! subroutine.  The second argument is the output of the gamma
        ! function.  It must also be declared as a Double Precision number
        ! in the code that calls Subroutine GAMMA.

        DOUBLE PRECISION FUNCTION D_fnc(V,X)

!       Reference: http://jin.ece.uiuc.edu/routines/mpbdv.for
!       COMPUTATION OF SPECIAL FUNCTIONS
!          by Shanjie Zhang and Jianming Jin
!       Copyright 1996 by John Wiley & Sons, Inc.
!
!       ====================================================
!       Purpose: Compute parabolic cylinder function Dv(x).
!                This is the D function of Whittaker,
!                Eq. 19.3.7 Abramowitz and Stegun.
!       Input:   x --- Argument of Dv(x)
!                v --- Order of Dv(x)
!       Output:  D_fnc = PDF --- Dv(x)
!       Extras:  DV(na) --- Dn+v0(x)
!                DP(na) --- Dn+v0'(x)
!                ( na = |n|, v0 = v-n, |v0| < 1,
!                  n = 0,ñ,ñ,úú)
!                PDD --- Dv'(x)
!       Routines called:
!             (1) DVSA for computing Dv(x) for small |x|
!             (2) DVLA for computing Dv(x) for large |x|
!       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DV(0:100),DP(0:100)
        XA=DABS(X)
        VH=V
        V=V+DSIGN(1.0D0,V)
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        EP=DEXP(-.25D0*X*X)
        IF (NA.GE.1) JA=1
        IF (V.GE.0.0) THEN
           IF (V0.EQ.0.0) THEN
              PD0=EP
              PD1=X*EP
           ELSE
              DO 10 L=0,JA
                 V1=V0+L
                 IF (XA.LE.5.8) CALL DVSA(V1,X,PD1)
                 IF (XA.GT.5.8) CALL DVLA(V1,X,PD1)
                 IF (L.EQ.0) PD0=PD1
10            CONTINUE
           ENDIF
           DV(0)=PD0
           DV(1)=PD1
           DO 15 K=2,NA
              PDF=X*PD1-(K+V0-1.0D0)*PD0
              DV(K)=PDF
              PD0=PD1
15            PD1=PDF
        ELSE  ! (V < 0.0)
           IF (X.LE.0.0) THEN
              IF (XA.LE.5.8D0)  THEN
                 CALL DVSA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVSA(V1,X,PD1)
              ELSE
                 CALL DVLA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVLA(V1,X,PD1)
              ENDIF
              DV(0)=PD0
              DV(1)=PD1
              DO 20 K=2,NA
                 PD=(-X*PD1+PD0)/(K-1.0D0-V0)
                 DV(K)=PD
                 PD0=PD1
20               PD1=PD
           ELSE IF (X.LE.2.0) THEN
              V2=NV+V0
              IF (NV.EQ.0) V2=V2-1.0D0
              NK=INT(-V2)
              CALL DVSA(V2,X,F1)
              V1=V2+1.0D0
              CALL DVSA(V1,X,F0)
              DV(NK)=F1
              DV(NK-1)=F0
              DO 25 K=NK-2,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 DV(K)=F
                 F1=F0
25               F0=F
           ELSE ! (X > 2.0)
              IF (XA.LE.5.8) CALL DVSA(V0,X,PD0)
              IF (XA.GT.5.8) CALL DVLA(V0,X,PD0)
              DV(0)=PD0
              M=100+NA
              F1=0.0D0
              F0=1.0D-30
              DO 30 K=M,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 IF (K.LE.NA) DV(K)=F
                 F1=F0
30               F0=F
              S0=PD0/F
              DO 35 K=0,NA
35               DV(K)=S0*DV(K)
           ENDIF ! X categories
        ENDIF  ! V categories
!        DO 40 K=0,NA-1
!           V1=ABS(V0)+K
!           IF (V.GE.0.0D0) THEN
!              DP(K)=0.5D0*X*DV(K)-DV(K+1)
!           ELSE
!              DP(K)=-0.5D0*X*DV(K)-V1*DV(K+1)
!           ENDIF
!40      CONTINUE
!        PDF=DV(NA-1)
!        PDD=DP(NA-1)
        V=VH

        D_fnc=DV(NA-1)

        RETURN
        END

!------------------------------------------------------------------------

        SUBROUTINE DVSA(VA,X,PD)

!       ===================================================
!       Purpose: Compute parabolic cylinder function Dv(x)
!                for small argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PD --- Dv(x)
!       Routine called: GAMMA for computing âx)
!       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        SQ2=DSQRT(2.0D0)
        EP=DEXP(-.25D0*X*X)
        VA0=0.5D0*(1.0D0-VA)
        IF (VA.EQ.0.0) THEN
           PD=EP
        ELSE
           IF (X.EQ.0.0) THEN
              IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0)) THEN
                 PD=0.0D0
              ELSE
                 CALL GAMMA(VA0,GA0)
                 PD=DSQRT(PI)/(2.0D0**(-.5D0*VA)*GA0)
              ENDIF
           ELSE
              CALL GAMMA(-VA,G1)
              A0=2.0D0**(-0.5D0*VA-1.0D0)*EP/G1
              VT=-.5D0*VA
              CALL GAMMA(VT,G0)
              PD=G0
              R=1.0D0
              DO 10 M=1,250
                 VM=.5D0*(M-VA)
                 CALL GAMMA(VM,GM)
                 R=-R*SQ2*X/M
                 R1=GM*R
                 PD=PD+R1
                 IF (DABS(R1).LT.DABS(PD)*EPS) GO TO 15
10            CONTINUE
15            PD=A0*PD
           ENDIF
        ENDIF
        RETURN
        END

!------------------------------------------------------------------------

        SUBROUTINE DVLA(VA,X,PD)

!       ====================================================
!       Purpose: Compute parabolic cylinder functions Dv(x)
!                for large argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PD --- Dv(x)
!       Routines called:
!             (1) VVLA for computing Vv(x) for large |x|
!             (2) GAMMA for computing âx)
!       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        EP=DEXP(-.25*X*X)
        A0=DABS(X)**VA*EP
        R=1.0D0
        PD=1.0D0
        DO 10 K=1,16
           R=-0.5D0*R*(2.0*K-VA-1.0)*(2.0*K-VA-2.0)/(K*X*X)
           PD=PD+R
           IF (DABS(R/PD).LT.EPS) GO TO 15
10      CONTINUE
15      PD=A0*PD
        IF (X.LT.0.0D0) THEN
            X1=-X
            CALL VVLA(VA,X1,VL)
            CALL GAMMA(-VA,GL)
            PD=PI*VL/GL+DCOS(PI*VA)*PD
        ENDIF
        RETURN
        END

!------------------------------------------------------------------------

        SUBROUTINE VVLA(VA,X,PV)

!       ===================================================
!       Purpose: Compute parabolic cylinder function Vv(x)
!                for large argument
!       Input:   x  --- Argument
!                va --- Order
!       Output:  PV --- Vv(x)
!       Routines called:
!             (1) DVLA for computing Dv(x) for large |x|
!             (2) GAMMA for computing âx)
!       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        QE=DEXP(0.25*X*X)
        A0=DABS(X)**(-VA-1.0D0)*DSQRT(2.0D0/PI)*QE
        R=1.0D0
        PV=1.0D0
        DO 10 K=1,18
           R=0.5D0*R*(2.0*K+VA-1.0)*(2.0*K+VA)/(K*X*X)
           PV=PV+R
           IF (DABS(R/PV).LT.EPS) GO TO 15
10      CONTINUE
15      PV=A0*PV
        IF (X.LT.0.0D0) THEN
           X1=-X
           CALL DVLA(VA,X1,PDL)
           CALL GAMMA(-VA,GL)
           DSL=DSIN(PI*VA)*DSIN(PI*VA)
           PV=DSL*GL/PI*PDL-DCOS(PI*VA)*PV
        ENDIF
        RETURN
        END

!------------------------------------------------------------------------

        SUBROUTINE GAMMA(X,GA)

!       ==================================================
!       Purpose: Compute gamma function âx)
!       Input :  x  --- Argument of âx)
!                       ( x is not equal to 0,-1,-2,úú
!       Output:  GA --- âx)
!       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     .          -0.6558780715202538D0, -0.420026350340952D-1,
     .          0.1665386113822915D0,-.421977345555443D-1,
     .          -.96219715278770D-2, .72189432466630D-2,
     .          -.11651675918591D-2, -.2152416741149D-3,
     .          .1280502823882D-3, -.201348547807D-4,
     .          -.12504934821D-5, .11330272320D-5,
     .          -.2056338417D-6, .61160950D-8,
     .          .50020075D-8, -.11812746D-8,
     .          .1043427D-9, .77823D-11,
     .          -.36968D-11, .51D-12,
     .          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

!------------------------------------------------------------------------
