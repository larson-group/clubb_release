MODULE module_bl_newhoc

! This module contains four subroutines and functions:
!
! 1) SUBROUTINE HOC_BL_SCHEME
!    This subroutine is called every timestep in order to provide results
!    from the HOC parameterization for use in WRF code.
!
! 2) FUNCTION LIN_EXT
!    This is a simple linear extension computation to use for any variable.
!
! 3) FUNCTION LIN_INT
!    This is a simple linear interpolation computation to use for any variable.
!    It is not currently being referenced in the code.
!
! 4) SUBROUTINE HOC_BL_INIT
!    This is a subroutine that is called once in the initialization part of
!    WRF.  It is important because it provides some essential initialization
!    information for some HOC variables (namely tau and wp2) that need to be
!    initialized before they are first called in parameterization_timestep
!    from SUBROUTINE HOC_BL_SCHEME.
!
!==============================================================================
!
!!!!  WRF VS. HOC vertical grid comparison.
!
!
!          WRF VERTICAL GRID           !          HOC VERTICAL GRID
!          #################           !          #################
!                                      !
!                                      !
!     TOP ========== Momentum (top)    !     TOP ========= Momentum (top)
!                                      !
!                                      !
!         ------- Thermodynamic (top-1)!         --------- Thermodynamic (top)
!                                      !
!                                      !
!         ========== Momentum (top-1)  !         ========= Momentum (top-1)
!              .                       !              .
!              .                       !              .
!              .                       !              .
!              .                       !              .
!              .                       !              .
!         ========== Momentum (k+1)    !         ========== Momentum (k+1)
!                                      !
!                                      !
!                                      !
!         ---------- Thermodynamic (k) !         ---------- Thermodynamic (k+1)
!                                      !
!                                      !
!                                      !
!         ========== Momentum (k)      !         ========== Momentum (k)
!              .                       !              .
!              .                       !              .
!              .                       !              .
!              .                       !              .
!              .                       !              .
!         ========== Momentum (2)      !         ========== Momentum(2)
!          *    |                      !          *    |
!          *    | 1/2 delta zw         !          *    | 1/2 delta zw
! delta zw *    |                      ! delta zw *    |
!         -*-------- Thermodynamic (1) !         -*-------- Thermodynamic (2)
!          *    |                      !          *    |
!          *    | 1/2 delta zw         !          *    | 1/2 delta zw
!          *    |                      !          *    |
! SURFACE ========== Momentum (1)      ! SURFACE ========== Momentum (1)
!         / / / / /                    !         / / / / /
!        / / / / /                     !        / / / / /
!                                      !
!                                      !         ---------- Thermodynamic (1)
!
!==============================================================================
!
!!!! Note on vertical advection in WRF-HOC interface.
!
! The vertical advection equation:  - w * dQ/dz; where Q is some property.
!
! When it is by itself, HOC handles this term for all the mean and the
! turbulent quantities.
!
! In the WRF-HOC interface, HOC handles the advection of all the turbulent
! quantities, for they don't appear in WRF.  However, this is not the case
! for the mean quantities.
!
! w = w_lss + w_WRF
!
! vertical velocity is equal to the imposed large-scale subsidence (w_lss)
! specified for each individual test case, plus the WRF resolved vertical
! velocity (w_WRF or w_2(i,k,j) in the WRF code).
!
! In it's code, WRF already handles the vertical advection of the mean 
! quantities due to it's own resolved vertical velocity.  In addition, the
! vertical advection of mean quantities due to specified large-scale
! subsidence in handled by the cumulus physics scheme created for each
! individual test case.
!
! Therefore, for mean quantities, the HOC code needs to show:
! 
! w = 0
!
! and for the turbulent quantities, the HOC code needs to show:
!
! w = w_lss + w_WRF
!
!==============================================================================
!  ############################################################################
!                         SET HOC CONSTANTS HERE!!!!!

REAL, PARAMETER:: C1_in          = 1.450000 ! Low Skewness in C1 Skewness Function.
REAL, PARAMETER:: C1b_in         = 3.500000 ! High Skewness in C1 Skewness Function.
REAL, PARAMETER:: C1c_in         = 1.000000 ! Degree of Slope of C1 Skewness Function.
REAL, PARAMETER:: C4_in          = 5.200000 ! Currently Not Used in the Model.
REAL, PARAMETER:: C5_in          = 0.300000 ! Constant from Chris Golaz.
REAL, PARAMETER:: C6rt_in        = 2.250000 ! Low Skewness in C6rt Skewness Function.
REAL, PARAMETER:: C6rtb_in       = 4.000000 ! High Skewness in C6rt Skewness Function.
REAL, PARAMETER:: C6rtc_in       = 1.000000 ! Degree of Slope of C6rt Skewness Function.
REAL, PARAMETER:: C6thl_in       = 2.250000 ! Low Skewness in C6thl Skewness Function.
REAL, PARAMETER:: C6thlb_in      = 4.000000 ! High Skewness in C6thl Skewness Function.
REAL, PARAMETER:: C6thlc_in      = 1.000000 ! Degree of Slope of C6thl Skewness Function.
REAL, PARAMETER:: C7_in          = 0.500000 ! Low Skewness in C7 Skewness Function.
REAL, PARAMETER:: C7b_in         = 0.750000 ! High Skewness in C7 Skewness Function.
REAL, PARAMETER:: C7c_in         = 1.000000 ! Degree of Slope of C7 Skewness Function.
REAL, PARAMETER:: C8_in          = 1.500000 ! Coefficient #1 in C8 Skewness Equation.
REAL, PARAMETER:: C8b_in         = 0.005000 ! Coefficient #2 in C8 Skewness Equation.  From Golaz.
REAL, PARAMETER:: C10_in         = 3.300000 ! Currently Not Used in the Model.
REAL, PARAMETER:: C11_in         = 0.900000 ! Low Skewness in C11 Skewness Function.
REAL, PARAMETER:: C11b_in        = 0.350000 ! High Skewness in C11 Skewness Function.  From Golaz.
REAL, PARAMETER:: C11c_in        = 0.500000 ! Degree of Slope of C11 Skewness Function. From Golaz.
REAL, PARAMETER:: C12_in         = 1.000000 ! Constant from Chris Golaz.
REAL, PARAMETER:: nu1_in         = 20.00000 ! Previously determined eddy viscosity.
REAL, PARAMETER:: nu2_in         = 5.000000 ! Previously determined eddy viscosity.
REAL, PARAMETER:: nu6_in         = 5.000000 ! Previously determined eddy viscosity.
REAL, PARAMETER:: nu8_in         = 20.00000 ! Previously determined eddy viscosity.
REAL, PARAMETER:: nu_r_in        = 3.000000 ! Tuned to LES from Surface Value of DYCOMS2 RF02 DS.
REAL, PARAMETER:: beta_in        = 1.750000 ! Constant from Chris Golaz.
REAL, PARAMETER:: gamma_coef_in  = 0.320000 ! Low Skewness in gamma coefficient Skewness Function.
REAL, PARAMETER:: gamma_coefb_in = 0.320000 ! High Skewness in gamma coefficient Skewness Function.
REAL, PARAMETER:: gamma_coefc_in = 5.000000 ! Degree of Slope of gamma coefficient Skewness Function.
REAL, PARAMETER:: c_K_in         = 0.548000 ! Constant from Chris Golaz.
REAL, PARAMETER:: taumin_in      = 90.00000 ! Previously determined value.
REAL, PARAMETER:: taumax_in      = 3600.000 ! Previously determined value.
REAL, PARAMETER:: lmin_coef_in   = 0.500000 ! Previously determined value.
REAL, PARAMETER:: mu_in          = 6.000E-4 ! Constant from Chris Golaz.

!  HOC Flags
LOGICAL, PARAMETER:: lbugsrad      = .false.
LOGICAL, PARAMETER:: kk_rain       = .false.
LOGICAL, PARAMETER:: cloud_sed     = .false.
LOGICAL, PARAMETER:: bottom_at_sfc = .true.

!  ############################################################################
!==============================================================================

CONTAINS

SUBROUTINE HOC_BL_SCHEME ( ITIMESTEP, DT, z_at_w, dz8w, z, config_flags,       &
                           u_phy, v_phy, rho, th_phy, p_phy, pi_phy, t_phy,    &
                           w_2, f, moist, n_moist, P_QV, P_QC, P_QR,           &
                           wpthlp_sgs, wprtp_sgs, wp2_sgs, rtp2_sgs,           &
                           thlp2_sgs, rtpthlp_sgs, upwp_sgs, vpwp_sgs,         &
                           Sc_3d, tau_3d,                                      &
                           wp3_sgs, Ncm_HOC, rcm_HOC,                          &
                           rrm_HOC, Nrm_HOC, cf_HOC,                           &
                           RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN,               &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           its, ite, jts, jte, kts, kte                        )

!------------------------------------------------------------------------------

! HOC modules that need to be referenced by this subroutine.
USE hoc_parameterization_interface
USE grid_class
#ifdef STATS
USE statistics
#endif /*STATS*/
! WRF modules that might need to be referenced by this subroutine.
USE module_configure
USE module_model_constants
USE module_state_description
! WRF module for surface fluxes.
USE module_diffusion_em
! WRF modules for case-specific forcings.
USE module_cu_GCSS

IMPLICIT NONE

TYPE(grid_config_rec_type), INTENT(IN):: config_flags

INTEGER, INTENT (IN):: ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       n_moist, P_QV, P_QC, P_QR

! Current time step of model run.
INTEGER, INTENT (IN):: ITIMESTEP

! Length of Time step.
REAL, INTENT (IN):: DT

! Heights of WRF Grid Points.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN)::  z_at_w, &  ! height at full-level of z.
                                                          dz8w, &  ! delta-z between full-levels of z.
                                                             z     ! height at half-level of z.

! Standard WRF variables to be brought into this subroutine.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme, n_moist), INTENT(IN)::   moist  ! water species mixing ratio.

! Standard WRF variables to be brought into this subroutine.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN)::   u_phy, &  ! u-wind at half-level of x.
                                                         v_phy, &  ! v-wind at half-level of y.
                                                           rho, &  ! density (dry air + vapor)
                                                        th_phy, &  ! theta
                                                         p_phy, &  ! pressure
                                                        pi_phy, &  ! exner function (p/p0)^(Rd/cp)
                                                         t_phy, &  ! Temperature
                                                           w_2     ! w-wind at full-level of z.

! Standard WRF variables to be brought into this subroutine.
REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: f  ! Coriolis paramter

! Standard WRF tendency variables to be taken out of this subroutine.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT)::   RUBLTEN, &  ! du/dt due to BL scheme.
                                                          RVBLTEN, &  ! dv/dt due to BL scheme.
                                                         RTHBLTEN, &  ! dth/dt due to BL scheme.
                                                         RQVBLTEN     ! dqv/dt due to BL scheme.

! HOC variables stored for entire 3-D WRF Field.
! Momentum-level variables.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT)::  wpthlp_sgs, &  ! overbar[w'thl'] on subgrid-scale.
                                                            wprtp_sgs, &  ! overbar[w'rt'] on subgrid-scale.
                                                              wp2_sgs, &  ! overbar[w'^2] on subgrid-scale.
                                                             rtp2_sgs, &  ! overbar[rt'^2] on subgrid-scale.
                                                            thlp2_sgs, &  ! overbar[thl'^2] on subgrid-scale.
                                                          rtpthlp_sgs, &  ! overbar[rt'thl'] on subgrid-scale.
                                                             upwp_sgs, &  ! overbar[u'w'] on subgrid-scale.
                                                             vpwp_sgs, &  ! overbar[v'w'] on subgrid-scale.
                                                                Sc_3d, &  ! Sc (sigma^2).
                                                               tau_3d     ! Time-scale Tau.

! HOC variables stored for entire 3-D WRF Field.
! Thermodynamic-level variables.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT)::  wp3_sgs, &  ! overbar[w'^3] on subgrid-scale.
                                                           rcm_HOC, &  ! HOC's mean cloud water mixing ratio.
                                                           Ncm_HOC, &  ! HOC's mean cloud droplet concentration.
                                                           rrm_HOC, &  ! HOC's mean rain water mixing ratio.
                                                           Nrm_HOC, &  ! HOC's mean rain droplet concentration.
                                                            cf_HOC     ! HOC's grid-box cloud fraction

! Local variables.

! HOC variables translated into 1-D HOC Field.
! Momentum-level variables.
REAL, DIMENSION(kms:kme)::  wpthlp, &  ! overbar[w'thl'] on subgrid-scale.
                             wprtp, &  ! overbar[w'rt'] on subgrid-scale.
                               wp2, &  ! overbar[w'^2] on subgrid-scale.
                              rtp2, &  ! overbar[rt'^2] on subgrid-scale.
                             thlp2, &  ! overbar[thl'^2] on subgrid-scale.
                           rtpthlp, &  ! overbar[rt'thl'] on subgrid-scale.
                              upwp, &  ! overbar[u'w'] on subgrid-scale.
                              vpwp, &  ! overbar[v'w'] on subgrid-scale.
                               Scm, &  ! Sc (sigma^2).
                              taum, &  ! Time-scale Tau.
                               wmm, &
                  momentum_heights     

! HOC variables translated into 1-D HOC Field.
! Thermodynamic-level variables.
REAL, DIMENSION(kms:kme)::  wp3, &  ! overbar[w'^3] on subgrid-scale.
                            rcm, &  ! HOC's mean cloud water mixing ratio.
                            Ncm, &  ! HOC's mean cloud droplet concentration.
                            rrm, &  ! HOC's mean rain water mixing ratio.
                            Nrm, &  ! HOC's mean rain droplet concentration.
                             cf, &  ! HOC's grid-box cloud fraction
                            wmt, &
                 thermo_heights, &
                        w_at_tl

! HOC variables translated into 1-D HOC Field.
! Thermodynamic-level variables.
REAL, DIMENSION(kms:kme)::  thetam, &
                              thlm, &
                               rtm, &
                                 p, &
                              rhot, &
                              rhom, &
                             exner, &
                                um, &
                                vm

! Interface with WRF Surface Fluxes.
REAL, DIMENSION(ims:ime,jms:jme):: sfc_wind, u_star
REAL, DIMENSION(ims:ime,jms:jme):: heat_flux, latent_flux
! Relevant HOC variables.
REAL:: wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc

! Interface with WRF case-specific forcings.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme):: w_sub, th_tend, rt_tend, r_vapor
INTEGER:: ra_lw_physics
! Relevant HOC variables.
REAL, DIMENSION(kms:kme):: w_lss, thlm_forcing, rtm_forcing

! Other variables needed for HOC import.
REAL:: deltaz, fcor
INTEGER:: nzmax

#ifdef SCALARS
!REAL, DIMENSION(kms:kme,sclrm_dimension):: sclrm, sclrm_forcing, edsclrm
!REAL, DIMENSION(kms:kme,sclrm_dimension):: edsclrmt
REAL, DIMENSION(kms:kme,2):: sclrm, sclrm_forcing, edsclrm
REAL, DIMENSION(kms:kme,2):: edsclrmt
#endif /*SCALARS*/

REAL:: dtclosure, dt_hoc, time

REAL, PARAMETER:: time_initial = 0.0
REAL, PARAMETER:: eps = 1.0E-10

INTEGER:: niterlong, iter

LOGICAL:: isValid = .true.

INTEGER:: i_start, i_end, j_start, j_end, k_start, k_end, i, j, k

#ifdef STATS
LOGICAL:: stdout = .true.
INTEGER:: x_center, y_center
LOGICAL:: sample_column = .false.
#endif /*STATS*/

!------------------------------------------------------------------------------

i_start = its
i_end   = MIN(ide-1, ite)
j_start = jts
j_end   = MIN(jde-1, jte)
k_start = kts
k_end   = MIN(kde-1, kte)


!  HOC uses theta-l and rt.  WRF uses theta, rv, rc, ...
!  We are just treating WRF's theta as theta-l and WRF's rv as rt.
!  ***We are not using any of WRF's microphysics***

!###############################################################################
!
! We need to import information about the surface heat flux, surface latent
! heat flux, and the value of u*.
!
! We need this information in order to compute the following fluxes at the 
! surface:  w'thl', w'rt', u'w', and v'w'.
!
! HOC uses these fluxes as input at the lowest level (surface level) for each
! of the above listed profiles.
!
! Brian:  I have previously created a subroutine in WRF's module_diffusion_em.F
!         which outputs surface fluxes (w'thl' and w'rt'), as well as the value
!         of u* (used for u'w' and v'w'), for each of the GCSS cases we have
!         previously set up in WRF.

CALL set_sfc_fluxes ( th_phy, p_phy, u_phy, v_phy, z,     &
                          rho, moist, n_moist, dt, itimestep, &
                          u_star, heat_flux, latent_flux,     &
                          config_flags,                       &
                          ids, ide, jds, jde, kds, kde,       &
                          ims, ime, jms, jme, kms, kme,       &
                          its, ite, jts, jte, kts, kte        )

!###############################################################################
!
! NOTES ON INTERFACE WITH GCSS SUBROUTINES
!
! We need to import information about the imposed large-scale subsidence for
! each test case to this subroutine.  This information is needed in order to 
! properly apply the vertical velocity, w, to each of the turbulent parameters
! in the HOC parameterization (such as w'rt', w'^2, w'thetal', etc.)
!
! By the way, I have previously set up the WRF code for all of the test cases.
! That portion of the WRF code (which is in the cumulus parameterization 
! section of WRF, referenced from module_cumulus_driver.F) computes the 
! GCSS-specific forcings.  This includes the effects of large-scale subsidence,
! other large-scale forcings, and radiative cooling.  Due to the structure of 
! these subroutines, we import other large-scale forcings (th_tend, rt_tend),
! as well as large-scale subsidence (w_sub).  However, the HOC code is set up 
! so that the information on the other large-scale forcings is not needed.
!
! If available, we also want information on cloud droplet concentration (Ncm).
!
! It should also be noted that information on cloud water mixing ratio,
! gathered from the HOC parameterization, is used by the GCSS subroutines in
! order to compute radiative cooling.

IF ( config_flags%cu_physics > 0 ) THEN

   SELECT CASE ( config_flags%cu_physics )

   ! BOMEX case
   CASE ( BOMEXSCHEME )
      CALL BOMEX_forcings( z, w_sub, th_tend, rt_tend,   &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte  )

   ! ARM case
   CASE ( ARMSCHEME )
      CALL ARM_forcings( z, dt, itimestep, th_tend, rt_tend, &
                         ids, ide, jds, jde, kds, kde,       &
                         ims, ime, jms, jme, kms, kme,       &
                         its, ite, jts, jte, kts, kte        )

   ! FIRE case
   CASE ( FIRESCHEME )

      DO i = i_start, i_end, 1
         DO j = j_start, j_end, 1
            DO k = k_start, k_end, 1
               r_vapor(i,k,j) = moist(i,k,j,P_QV) - rcm_HOC(i,k,j)
            ENDDO
         ENDDO
      ENDDO

      ra_lw_physics = config_flags%ra_lw_physics

      CALL FIRE_forcings( z, dz8w, rho, th_phy, r_vapor,     &
                          rcm_HOC, pi_phy, p_phy, w_sub,     &
                          th_tend, ra_lw_physics,            &
                          ids, ide, jds, jde, kds, kde,      &
                          ims, ime, jms, jme, kms, kme,      &
                          its, ite, jts, jte, kts, kte       )

   ! ATEX case
   CASE ( ATEXSCHEME )

      DO i = i_start, i_end, 1
         DO j = j_start, j_end, 1
            DO k = k_start, k_end, 1
               r_vapor(i,k,j) = moist(i,k,j,P_QV) - rcm_HOC(i,k,j)
            ENDDO
         ENDDO
      ENDDO

      CALL ATEX_forcings( rho, pi_phy, z, dz8w, r_vapor,         &
                          rcm_HOC, p_phy, th_phy, dt, itimestep, &
                          w_sub, th_tend, rt_tend,               &
                          ids, ide, jds, jde, kds, kde,          &
                          ims, ime, jms, jme, kms, kme,          &
                          its, ite, jts, jte, kts, kte           )

   CASE DEFAULT

      DO i = i_start, i_end, 1
         DO j = j_start, j_end, 1
            DO k = k_start, k_end, 1
               w_sub(i,k,j) = 0.0
               th_tend(i,k,j) = 0.0
               rt_tend(i,k,j) = 0.0
            ENDDO
         ENDDO
      ENDDO


   END SELECT

ENDIF

!###############################################################################

! CALL THE HOC PARAMETERIZATION FOR EVERY COLUMN IN THE 3-D WRF FIELD.
DO i = i_start, i_end, 1
   DO j = j_start, j_end, 1

      !#########################################################################
      !#########################################################################
      !#########################################################################

      ! Setup momentum-level grid heights.  This is a one-dimensional array
      ! used for HOC based on the heights from the three-dimensional array
      ! used by WRF.
      DO k = k_start, kte, 1
         momentum_heights(k) = z_at_w(i,k,j)
      ENDDO

      ! Setup thermodynamic-level grid heights.  This is also the same as
      ! listed above, but for thermodynamic levels.  One should note that
      ! thermodynamic levels are always located exactly halfway between
      ! the momentum levels that encase them.  The indices of WRF's
      ! thermodynamics levels and HOC's thermodynamic levels differ by one
      ! due to the fact that HOC's first thermodynamic level is below the
      ! surface (see above diagram).
      DO k = k_start+1, kte, 1
         thermo_heights(k) = z(i,k-1,j)
      ENDDO
      thermo_heights(k_start) = -thermo_heights(k_start+1)

      ! Number of vertical levels (momentum).
      nzmax = kte
      ! deltaz is for HOC as a model by itself.  It is useless when HOC
      ! is inside of a host model.  Since it needs to be sent in anyway
      ! as part of the argument list, I will define it as the average
      ! height of the vertical layers.
      deltaz =   ( momentum_heights(kte) - momentum_heights(k_start) )    &
               / ( nzmax - 1 )


      CALL parameterization_setup                                         &
           ( nzmax, deltaz, C1_in, C1b_in, C1c_in,                        &
             C4_in, C5_in,                                                &
             C6rt_in, C6rtb_in, C6rtc_in, C6thl_in,                       &
             C6thlb_in, C6thlc_in, C7_in, C7b_in, C7c_in,                 &
             C8_in, C8b_in, C10_in, C11_in, C11b_in, C11c_in, C12_in,     &
             nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,                     &
             gamma_coef_in, gamma_coefb_in, gamma_coefc_in,               &
             mu_in, c_K_in, beta_in, lmin_coef_in,                        &
             taumin_in, taumax_in, fcor,                                  &
             lbugsrad, kk_rain, cloud_sed, bottom_at_sfc,                 &
             .true., momentum_heights, thermo_heights                     )
      ! HOC grid is now setup and HOC grid operators can now be used.
      ! All of the HOC variables in HOC's prognostic arrays have been
      ! allocated and initialized.  We must overwrite that initialization
      ! below.

      !------------------------------------------------------------------------

      ! Surface heat and moisture fluxes.
      wpthlp_sfc = heat_flux(i,j)
      wprtp_sfc  = latent_flux(i,j)

      ! Surface momentum fluxes.
      sfc_wind(i,j) = MAX( SQRT( u_phy(i,kts,j)**2.0 + v_phy(i,kts,j)**2.0 ), &
                           eps )
      upwp_sfc = - ( u_star(i,j)**2.0 ) * ( u_phy(i,kts,j) / sfc_wind(i,j) )
      vpwp_sfc = - ( u_star(i,j)**2.0 ) * ( v_phy(i,kts,j) / sfc_wind(i,j) )

      ! Large-scale subsidence and large-scale forcings.
      ! Remember that HOC thermodynamic-level indices are offset
      ! from WRF thermodynamic-level indices by one.  See above diagram.
      DO k = k_start+1, kte, 1
         w_lss(k)        = w_sub(i,k-1,j)
!         thlm_forcing(k) = th_tend(i,k-1,j)
!         rtm_forcing(k)  = rt_tend(i,k-1,j)
         ! Note from Brian:  theta-l and rt forcings in this model are 
         ! done in the cumulus physics schemes for the different test
         ! cases.  This is done in a different part of WRF.  In my case,
         ! they do not need to be brought into HOC.
         thlm_forcing(k) = 0.0
         rtm_forcing(k)  = 0.0
      ENDDO
      ! The lowest HOC thermodynamic-level is below ground.  Set forcings
      ! to 0.0 there.
      w_lss(k_start)        = 0.0
      thlm_forcing(k_start) = 0.0
      rtm_forcing(k_start)  = 0.0

      ! Cloud droplet concentration (if applicable).
      DO k = k_start+1, kte, 1
         Ncm(k) = 0.0
      ENDDO
      Ncm(k_start) = 0.0

      ! Coriolis force
      fcor = f(i,j)

      !------------------------------------------------------------------------

      ! Set up the 1-D HOC variables for every grid column.
      ! Momentum-level variables.
      DO k = k_start, kte, 1

         wpthlp(k)  = wpthlp_sgs(i,k,j)
         wprtp(k)   = wprtp_sgs(i,k,j)
         wp2(k)     = wp2_sgs(i,k,j)
         rtp2(k)    = rtp2_sgs(i,k,j)
         thlp2(k)   = thlp2_sgs(i,k,j)
         rtpthlp(k) = rtpthlp_sgs(i,k,j)
         upwp(k)    = upwp_sgs(i,k,j)
         vpwp(k)    = vpwp_sgs(i,k,j)
         Scm(k)     = Sc_3d(i,k,j)
         taum(k)    = tau_3d(i,k,j)

      ENDDO

      ! Thermodynamic-level variables.
      DO k = k_start+1, kte, 1

         wp3(k)     = wp3_sgs(i,k-1,j)
         rcm(k)     = rcm_HOC(i,k-1,j)
         Ncm(k)     = Ncm_HOC(i,k-1,j)
         rrm(k)     = rrm_HOC(i,k-1,j)
         Nrm(k)     = Nrm_HOC(i,k-1,j)
         cf(k)      = cf_HOC(i,k-1,j)

         ! The following thermodynamic-level variables are passed from WRF
         ! to HOC as INPUT ONLY!  They are not directly altered by HOC.
         thetam(k)  = th_phy(i,k-1,j)
         rtm(k)     = moist(i,k-1,j,P_QV) !+ moist(i,k-1,j,P_QC)
         um(k)      = u_phy(i,k-1,j)
         vm(k)      = v_phy(i,k-1,j)
         p(k)       = p_phy(i,k-1,j)
         exner(k)   = pi_phy(i,k-1,j)
         rhot(k)    = rho(i,k-1,j)

         thlm(k)    = thetam(k) !- (  XLV / ( cp * exner(k) )  ) * rcm(k)

         ! Note on rtm and thlm:  WRF uses theta and rvm.  However, I am
         ! turning off WRF's internal microphysics (mp_physics = 0), allowing
         ! WRF's theta and rvm to be treated as theta-l and rtm (which are the
         ! variables that HOC uses).  Therefore, there is no need to convert
         ! from theta to theta-l and from rvm to rtm.  Cloud water mixing ratio
         ! is determined from HOC's internal subgrid microphysics.  Rain water
         ! mixing ratio can also be determined from HOC's subgrid microphysics 
         ! (based on scheme by Khairoutdinov and Kogan (2000)).  In order to do
         ! so, the kk_rain logical flag at the top must be set to .true.  Brian 

      ENDDO

      ! Lowest-level Thermodynamic HOC setup.
      wp3(k_start)    = 0.0
      rcm(k_start)    = 0.0
      Ncm(k_start)    = 0.0
      rrm(k_start)    = rrm(k_start+1)
      Nrm(k_start)    = Nrm(k_start+1)
      cf(k_start)     = 0.0

      thetam(k_start) = thetam(k_start+1)
      rtm(k_start)    = rtm(k_start+1)
      um(k_start)     = um(k_start+1)
      vm(k_start)     = vm(k_start+1)
      ! Pressure always increases as you go down, so use linear extension.
!      p(k_start)      = LIN_EXT( p(k_start+2), p(k_start+1),                 & 
!                        gr%zt(k_start+2), gr%zt(k_start+1), gr%zt(k_start) )
      ! To stay consistent with HOC code, we record pressure at the lowest
      ! thermodynamic level (which is below the surface) to be the pressure
      ! at the surface itself.  The same is true for exner and density.
      p(k_start)      = LIN_EXT( p(k_start+2), p(k_start+1),                 & 
                        gr%zt(k_start+2), gr%zt(k_start+1), gr%zm(k_start) )
      ! The exner function at any level is always equal to (p/p0)^(Rd/Cp)
      exner(k_start)  = (p(k_start)/p1000mb)**rcp
      ! Computation of density:
      ! Temp = Theta * Exner
      ! rho = p / (Rd * Temp)
      ! therefore, rho = p / (Rd * Theta * Exner)
      rhot(k_start)   = p(k_start) / ( r_d * thetam(k_start) * exner(k_start) )

      thlm(k_start)   = thlm(k_start+1)

      ! Interpolate rho to momentum levels.
      ! rho is a thermodynamic-level variable and the
      ! HOC variable rhot stores these values.  However,
      ! we also want to know rho at momentum levels.  The
      ! HOC variable rhom stores these values.
      rhom(k_start) = rhot(k_start)
!      DO k = k_start, kte, 1
      DO k = k_start+1, kte, 1
         rhom(k) = zt2zm( rhot, k )
      ENDDO

      ! ##### Set up vertical velocity profile #####

      ! WRF's resolved vertical velocity (w_2) is located on
      ! momentum levels.  It will be interpolated to thermodynamic
      ! levels.  The thermodynamic levels are always located exactly
      ! halfway between the momentum levels that sandwich them.
      ! This is using WRF grid indices.
      DO k = k_start, k_end, 1
         w_at_tl(k) = (1.0/2.0) * ( w_2(i,k+1,j) + w_2(i,k,j) )
      ENDDO

      ! Large-scale subsidence has been brought in on thermodynamic levels.
      ! It will be added to the WRF resolved vertical velocity at
      ! thermodynamic levels.
      ! Since wmt and w_lss are HOC grid variables and w_at_tl is a WRF grid
      ! variables, the grid indices will be offset by one (since we are 
      ! dealing with thermodynamic levels).
      DO k = k_start+1, kte, 1
         wmt(k) = w_lss(k) + w_at_tl(k-1)
      ENDDO
      ! The lowest HOC grid level is always below ground.
      wmt(k_start) = 0.0

      ! Find the vertical velocity at momentum levels.
      DO k = k_start, kte, 1
         wmm(k) = zt2zm( wmt, k )
      ENDDO
      ! Vertical velocity will be zero at the ground surface.
      wmm(k_start) = 0.0


      ! ##### Setup HOC Parameterization Timestep Information #####
      dtclosure = DT

!!     Compute number of iterations for closure loop.
!      IF ( time > time_spinup ) THEN
!         niterlong = 1
!         dt_hoc    = DT
!      ELSE
         niterlong = INT( DT / dtclosure )
         dt_hoc    = dtclosure
!      ENDIF

      ! Pre-HOC-loop time
      time = time_initial + (ITIMESTEP-1) * DT

#ifdef STATS

         ! Only record stats for one column,
         ! so choose central or a border-central column.
         x_center = INT( ( i_end + i_start )/2.0 )
         y_center = INT( ( j_end + j_start )/2.0 )

         IF ( i == x_center .AND. j == y_center ) THEN

            sample_column = .true.

!           When this time step is over, the time will be time + DT
            CALL stats_begin_timestep( time + DT, DT )
            if ( lstats_last .and. stdout ) then
               write(*,'(a,i7,a,f8.1)') 'iteration = ',                &
                                         ITIMESTEP, '; time = ', time
            endif

         ELSE

            sample_column = .false.

         ENDIF

#endif /*STATS*/

      ! CALL HOC PARAMETERIZATION TIMESTEP.
      ! This can be done once or multiple times.
      DO iter = 1, niterlong, 1

         CALL parameterization_timestep                          &
               ( iter, niterlong, time_initial, time, dt_hoc,    &
                 thlm_forcing, rtm_forcing, wmm, wmt,            &
                 wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,      &
                 p, rhom, rhot, exner,                           &
                 um, vm, upwp, vpwp,                             &
                 thlm, rtm, wprtp, wpthlp, wp2, wp3,             &
                 rtp2, thlp2, rtpthlp,                           &
                 Scm, taum, rcm, Ncm, rrm, Nrm, cf,              &
                 isValid, .true.                                 &
#ifdef STATS
                 ,sample_column                                  &
#endif /*STATS*/
#ifdef SCALARS
                 ,sclrm, sclrm_forcing,                          &
                 edsclrm, edsclrmt                               &
#endif /*SCALARS*/
               )


#ifdef STATS

         IF ( i == x_center .AND. j == y_center ) THEN

            CALL stats_end_timestep( time + DT, DT )

         ENDIF

#endif /*STATS*/

         ! SET TIME (during HOC-loop)
         ! Advance time here, not in hoc_closure_timestep,
         ! in order to facilitate use of stats.
         ! A host model, e.g. WRF, would advance time outside
         ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
         IF ( iter < niterlong ) THEN
            time = time_initial + (ITIMESTEP-1) * DT + iter * dt_hoc
         ELSEIF ( iter == niterlong ) THEN
            time = time_initial + ITIMESTEP * DT
         ENDIF

         if ( .not. isValid ) exit

      ENDDO 

      !########## SAVE OUTPUT FROM HOC PARAMETERIZATION ##########

      ! Overwrite the WRF Storage variables with the 1-D HOC results for every
      ! grid column.
      ! NOTE:  All of the momentum-level variables listed below were specially
      !        created 3-D variables in the WRF Registry.  They are being used
      !        to save all of the HOC output (from every single column) for use
      !        by HOC during the next timestep.  The same is true for most of
      !        the thermodynamic-level variables listed below, namely wp3_sgs,
      !        rcm_HOC, Ncm_HOC, rrm_HOC, Nrm_HOC, and cf_HOC.  However, the 
      !        thermodynamic-level variables moist(i,k,j,P_QC) and
      !        moist(i,k,j,P_QR) are WRF's own variables for cloud water mixing
      !        ratio and rain water mixing ratio, respectively.  Usually, they
      !        are updated in the microphysics portion of WRF, while various
      !        moisture tendencies are saved along the way, only to be entered
      !        into the equation in the microphysics portion of WRF.  However,
      !        all WRF microphysics are turned off in order to run HOC.  HOC 
      !        uses its own intrinsic microphysics instead.  We make WRF track
      !        total moisture mixing ratio using its moist(i,k,j,P_QV), which 
      !        is its normal variable for water vapor mixing ratio.  However, 
      !        when GCSS cases are being run, values need to be saved in the WRF
      !        moist(P_QC and P_QR) arrays.  This is because those arrays allow
      !        WRF to take the data to its cumulus_driver portions.  There, 
      !        GCSS subroutines use the data to compute such things as radiative
      !        cooling.  There are no other subroutines that will alter the data
      !        between its output here and its input into the GCSS subroutines.
      !        Once the data has been used there, it does not matter what else
      !        happens to it.  It is not used again for HOC purposes.  It will
      !        simply be overwritten here at the next timestep.  However, the
      !        only way to turn off WRF's microphysics is to use the passiveqv
      !        microphysics option (mp_physics = 0 in the WRF namelist).  This
      !        allows for advective and diffusional effects on the water vapor
      !        (in our case, total moisture) field, but does not effect it with
      !        microphysics.  The unfortunate thing with WRF is that using this
      !        option means that only the moist(P_QV) variable will be
      !        available.  In order to have the same effects while using the
      !        moist(P_QC and P_QR) variables, one has to create a brand new
      !        package in the WRF Registry.  In addition to WRF's own passiveqv
      !        package (where only P_QV is available), I have also created
      !        passiveqvqc (where P_QV and P_QC are available) and
      !        passiveqvqcqr (where P_QV, P_QC, and P_QR are available).  Brian

      ! Momentum-level variables.
      DO k = k_start, kte, 1

         wpthlp_sgs(i,k,j)  = wpthlp(k)
         wprtp_sgs(i,k,j)   = wprtp(k)
         wp2_sgs(i,k,j)     = wp2(k)
         rtp2_sgs(i,k,j)    = rtp2(k)
         thlp2_sgs(i,k,j)   = thlp2(k)
         rtpthlp_sgs(i,k,j) = rtpthlp(k)
         upwp_sgs(i,k,j)    = upwp(k)
         vpwp_sgs(i,k,j)    = vpwp(k)
         Sc_3d(i,k,j)       = Scm(k)
         tau_3d(i,k,j)      = taum(k)

      ENDDO

      ! Thermodynamic-level variables.
      ! Once again, note that WRF's vertical grid indices and HOC's vertical
      ! grid indices differ by 1.  See above diagram.
      DO k = k_start, k_end, 1

         wp3_sgs(i,k,j)     = wp3(k+1)
         rcm_HOC(i,k,j)     = rcm(k+1)
         moist(i,k,j,P_QC)  = rcm(k+1)
         Ncm_HOC(i,k,j)     = Ncm(k+1)
         rrm_HOC(i,k,j)     = rrm(k+1)
         moist(i,k,j,P_QR)  = rrm(k+1)
         Nrm_HOC(i,k,j)     = Nrm(k+1)
         cf_HOC(i,k,j)      = cf(k+1)

      ENDDO

      ! Note:  HOC's results for the first thermodynamic-level are thrown out,
      !        for they are located below the surface level and WRF does not
      !        have a level there.

      ! Free memory
      CALL parameterization_cleanup( .true. )
      ! All of HOC's prognostic variable arrays have now been deallocated.
      ! They have been eliminated from memory.

      !#########################################################################
      !#########################################################################
      !#########################################################################

   ENDDO  ! OUTER j-loop
ENDDO     ! OUTER i-loop



! Tendencies to send back into WRF:
!
! (dTH/dt)_BL_HOC = - d/dz [w'thl']
! (dQV/dt)_BL_HOC = - d/dz [w'rt']
! (dU/dt)_BL_HOC  = - d/dz [u'w']
! (dV/dt)_BL_HOC  = - d/dz [v'w']

DO i = i_start, i_end, 1
   DO j = j_start, j_end, 1
      DO k = k_start, k_end, 1
         RTHBLTEN(i,k,j) = - ( wpthlp_sgs(i,k+1,j) - wpthlp_sgs(i,k,j) ) / dz8w(i,k,j)
         RQVBLTEN(i,k,j) = - ( wprtp_sgs(i,k+1,j) - wprtp_sgs(i,k,j) ) / dz8w(i,k,j)
         RUBLTEN(i,k,j)  = - ( upwp_sgs(i,k+1,j) - upwp_sgs(i,k,j) ) / dz8w(i,k,j)
         RVBLTEN(i,k,j)  = - ( vpwp_sgs(i,k+1,j) - vpwp_sgs(i,k,j) ) / dz8w(i,k,j)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE HOC_BL_SCHEME

!==============================================================================

FUNCTION LIN_EXT( var_high, var_low, height_high, height_low, height_ext )

! This function computes a linear extension of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height outside of those two height levels 
! (rather than a height between those two height levels) is computed.
!
! Here is a diagram:
!
!  -------------------------------- Height to be extended to; linear extension
!
!  ################################ Height high, know variable value
!
!
!
!  ################################ Height low, know variable value
!
!
!
!  -------------------------------- Height to be extended to; linear extension
!
!
! FORMULA:
!
! variable(@ Height extension) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height extension - Height high)  +  variable(@ Height high)

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_ext

! Output Variable
REAL:: lin_ext

lin_ext = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_ext - height_high ) + var_high


END FUNCTION LIN_EXT

!==============================================================================

FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: lin_int

lin_int = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT

!==============================================================================

SUBROUTINE HOC_BL_INIT ( RUBLTEN, RVBLTEN, RTHBLTEN, RQVBLTEN,                 &
                         restart, config_flags, DT,                            &
                         wpthlp_sgs, wprtp_sgs, wp2_sgs, rtp2_sgs,             &
                         thlp2_sgs, rtpthlp_sgs, upwp_sgs, vpwp_sgs,           &
                         Sc_3d, tau_3d,                                        &
                         wp3_sgs, rcm_HOC, Ncm_HOC, rrm_HOC,                   &
                         Nrm_HOC, cf_HOC,                                      &
                         u_2, v_2,                                             &
                         t_2, pp, pb, php, phb, tke_2, f,                      &
                         moist, n_moist, P_QV, P_QC,                           &
                         ids, ide, jds, jde, kds, kde,                         &
                         ims, ime, jms, jme, kms, kme,                         &
                         its, ite, jts, jte, kts, kte                          )

! HOC modules that need to be referenced by this subroutine.
USE hoc_parameterization_interface
USE grid_class
#ifdef STATS
USE statistics
#endif /*STATS*/
! WRF modules that might need to be referenced by this subroutine.
USE module_configure
USE module_model_constants

TYPE(grid_config_rec_type), INTENT(IN):: config_flags

INTEGER, INTENT (IN):: ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       n_moist, P_QV, P_QC

LOGICAL, INTENT(IN):: restart

REAL, DIMENSION(ims:ime,kms:kme,jms:jme,n_moist), INTENT(IN)::  moist

REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN)::  t_2, pp, pb, php, phb, tke_2

REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN):: f

REAL, INTENT(IN):: DT

REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN)::  u_2, v_2

! Standard WRF tendency variables to be initialized in this subroutine.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT)::   RUBLTEN, &  ! du/dt due to BL scheme.
                                                          RVBLTEN, &  ! dv/dt due to BL scheme.
                                                         RTHBLTEN, &  ! dth/dt due to BL scheme.
                                                         RQVBLTEN     ! dqv/dt due to BL scheme.

! HOC variables initialized for entire 3-D WRF Field.
! Momentum-level variables.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT)::  wpthlp_sgs, &  ! overbar[w'thl'] on subgrid-scale.
                                                          wprtp_sgs, &  ! overbar[w'rt'] on subgrid-scale.
                                                            wp2_sgs, &  ! overbar[w'^2] on subgrid-scale.
                                                           rtp2_sgs, &  ! overbar[rt'^2] on subgrid-scale.
                                                          thlp2_sgs, &  ! overbar[thl'^2] on subgrid-scale.
                                                        rtpthlp_sgs, &  ! overbar[rt'thl'] on subgrid-scale.
                                                           upwp_sgs, &  ! overbar[u'w'] on subgrid-scale.
                                                           vpwp_sgs, &  ! overbar[v'w'] on subgrid-scale.
                                                              Sc_3d, &  ! Sc (sigma^2).
                                                             tau_3d     ! Time-scale Tau.

! HOC variables initialized for entire 3-D WRF Field.
! Thermodynamic-level variables.
REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT)::  wp3_sgs, &  ! overbar[w'^3] on subgrid-scale.
                                                         rcm_HOC, &  ! HOC's mean cloud water mixing ratio.
                                                         Ncm_HOC, &  ! HOC's mean cloud droplet concentration.
                                                         rrm_HOC, &  ! HOC's mean rain water mixing ratio.
                                                         Nrm_HOC, &  ! HOC's mean rain droplet concentration.
                                                          cf_HOC     ! HOC's grid-box cloud fraction

! Local variables.
REAL, PARAMETER:: emin = 1.0E-6
REAL, PARAMETER:: taumin = 90.0
REAL, PARAMETER:: taumax = 3600.0

!
REAL, DIMENSION(kms:kme):: momentum_heights, thermo_heights
REAL, DIMENSION(kms:kme):: tke_HOC_t, em
REAL, DIMENSION(kms:kme):: thetam, p, rtm, rcm, exner, thlm, thvm
REAL, DIMENSION(kms:kme):: Lscale, lup, ldown
REAL, DIMENSION(kms:kme):: um, vm

INTEGER:: nzmax
REAL:: deltaz, fcor

INTEGER:: i_start, i_end, j_start, j_end, k_start, k_end, i, j, k

#ifdef STATS
CHARACTER(50):: runfile
INTEGER:: day, month, year
REAL:: rlat, rlon, time, time_initial
INTEGER:: x_center, y_center
#endif /*STATS*/



!------------------------------------------------------------------------------

i_start = its
i_end   = MIN(ide-1, ite)
j_start = jts
j_end   = MIN(jde-1, jte)
k_start = kts
k_end   = MIN(kde-1, kte)

IF (.not. restart) THEN

   DO i = i_start, i_end, 1
      DO j = j_start, j_end, 1

         ! Initialize WRF 3-D Storage for HOC variables.
         DO k = kts, kte, 1

            RUBLTEN(i,k,j)     = 0.0
            RVBLTEN(i,k,j)     = 0.0
            RTHBLTEN(i,k,j)    = 0.0
            RQVBLTEN(i,k,j)    = 0.0

            wpthlp_sgs(i,k,j)  = 0.0
            wprtp_sgs(i,k,j)   = 0.0
            wp2_sgs(i,k,j)     = (2.0/3.0)*emin
            rtp2_sgs(i,k,j)    = 0.0
            thlp2_sgs(i,k,j)   = 0.0
            rtpthlp_sgs(i,k,j) = 0.0
            upwp_sgs(i,k,j)    = 0.0
            vpwp_sgs(i,k,j)    = 0.0
            Sc_3d(i,k,j)       = 0.0
            tau_3d(i,k,j)      = 0.0

            wp3_sgs(i,k,j)     = 0.0
            rcm_HOC(i,k,j)     = 0.0
            Ncm_HOC(i,k,j)     = 0.0
            rrm_HOC(i,k,j)     = 0.0
            Nrm_HOC(i,k,j)     = 0.0
            cf_HOC(i,k,j)      = 0.0

         ENDDO

         ! Set value of the Coriolis parameter.
         fcor = f(i,j)

         ! Setup momentum-level grid heights.  This is a one-dimensional array
         ! used for HOC based on the heights from the three-dimensional array
         ! used by WRF.
         DO k = k_start, kte, 1
            momentum_heights(k) = ( phb(i,k,j) + php(i,k,j) ) / g
         ENDDO

         ! Setup thermodynamic-level grid heights.  This is also the same as
         ! listed above, but for thermodynamic levels.  One should note that
         ! thermodynamic levels are always located exactly halfway between
         ! the momentum levels that encase them.  The indices of WRF's
         ! thermodynamics levels and HOC's thermodynamic levels differ by one
         ! due to the fact that HOC's first thermodynamic level is below the
         ! surface (see above diagram).
         DO k = k_start+1, kte, 1
            thermo_heights(k) =  (1.0/2.0)                                     &
                               * ( momentum_heights(k) + momentum_heights(k-1) )
         ENDDO
         thermo_heights(k_start) = - thermo_heights(k_start+1)

         ! Number of vertical levels (momentum).
         nzmax = kte
         ! deltaz is for HOC as a model by itself.  It is useless when HOC
         ! is inside of a host model.  Since it needs to be sent in anyway
         ! as part of the argument list, I will define it as the average
         ! height of the vertical layers.
         deltaz =   ( momentum_heights(kte) - momentum_heights(k_start) )    &
                  / ( nzmax - 1 )


         ! Set up HOC interface.
         CALL parameterization_setup                                         &
              ( nzmax, deltaz, C1_in, C1b_in, C1c_in,                        &
                C4_in, C5_in,                                                &
                C6rt_in, C6rtb_in, C6rtc_in, C6thl_in,                       &
                C6thlb_in, C6thlc_in, C7_in, C7b_in, C7c_in,                 &
                C8_in, C8b_in, C10_in, C11_in, C11b_in, C11c_in, C12_in,     &
                nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,                     &
                gamma_coef_in, gamma_coefb_in, gamma_coefc_in,               &
                mu_in, c_K_in, beta_in, lmin_coef_in,                        &
                taumin_in, taumax_in, fcor,                                  &
                lbugsrad, kk_rain, cloud_sed, bottom_at_sfc,                 &
                .true., momentum_heights, thermo_heights                     )
         ! HOC grid is now setup and HOC grid operators can now be used.
         ! All of the HOC variables in HOC's prognostic arrays have been
         ! allocated and initialized.  We must overwrite that initialization
         ! below.


         ! ########## TKE Initialization ##########

         ! Initial subgrid TKE profiles for all cases that are setup in WRF
         ! are located in WRF's module_initialize_quarter_ss.F.  There, the
         ! initial value of subgrid TKE is set as tke_2.
         ! WRF's TKE is located on thermodynamic levels, while HOC's TKE
         ! is located on momentum levels.
         ! Start by using the value of TKE from WRF's thermodynamic levels
         ! for HOC's equivalent thermodynamic levels.
         DO k = kts+1, kte, 1
            tke_HOC_t(k) = tke_2(i,k-1,j)
         ENDDO
         ! Find the value of TKE for HOC at level kts (the below-ground level)
         ! from a linear extension of the values of TKE at the first two
         ! above-ground levels.
         tke_HOC_t(kts) = LIN_EXT( tke_HOC_t(kts+2), tke_HOC_t(kts+1),       &
                          gr%zt(kts+2), gr%zt(kts+1), gr%zt(kts) )

         ! HOC's value of TKE on it's own thermodynamic levels now needs to be
         ! interpolated to it's momentum levels.
         DO k = kts, kte, 1
            em(k) = zt2zm( tke_HOC_t, k )
            em(k) = MAX( em(k), emin )
         ENDDO

         ! Initial value of subgrid-scale wp2 is equal to 2/3 of 
         ! subgrid-scale TKE
         DO k = kts, kte, 1
            wp2_sgs(i,k,j) = (2.0/3.0) * em(k)
         ENDDO

         ! ########## Initialization of time-scale Tau ##########

         ! First, the length-scale (Lscale) must be found.
         ! Set up thermodynamic-level variables.
         DO k = k_start+1, kte, 1 

            um(k)     = (1.0/2.0) * ( u_2(i,k-1,j) + u_2(i+1,k-1,j) )
            vm(k)     = (1.0/2.0) * ( v_2(i,k-1,j) + v_2(i,k-1,j+1) )

            thetam(k) = t_2(i,k-1,j) + t0
            p(k)      = pp(i,k-1,j) + pb(i,k-1,j)
            rtm(k)    = moist(i,k-1,j,P_QV) !+ moist(i,k-1,j,P_QC)
            rcm(k)    = rcm_HOC(i,k-1,j)
!            rcm(k)    = moist(i,k-1,j,P_QC)
            exner(k)  = (p(k)/p1000mb)**rcp

            thlm(k)   = thetam(k) !- ( XLV / (cp*exner(k)) ) * rcm(k)
            thvm(k)   = thlm(k)                                    &
                       + EP_1 * t0 * rtm(k)                         &
                       + ( XLV / (cp*exner(k)) - rvovrd * t0 ) * rcm(k)

         ENDDO

         ! Set up lowest-level (below ground in HOC) thermodynamic variables.
         um(k_start)     = um(k_start+1)
         vm(k_start)     = vm(k_start+1)
         thetam(k_start) = thetam(k_start+1)
         p(k_start)      = p(k_start+1)
         rtm(k_start)    = rtm(k_start+1)
         rcm(k_start)    = 0.0
!         rcm(k_start)    = 0.0
         exner(k_start)  = exner(k_start+1)
         thlm(k_start)   = thlm(k_start+1)
         thvm(k_start)   = thvm(k_start+1)

         ! Call HOC subroutine to compute length-scale.
         CALL length( thvm, thlm, rtm, rcm, em, p, exner, &
                      Lscale, lup, ldown )

         ! Finally, compute time-scale tau
         DO k = k_start, kte, 1       

            IF ( wp2_sgs(i,k,j) < 0.005 ) THEN
               tau_3d(i,k,j) = taumin
            ELSE
               tau_3d(i,k,j) = zt2zm(Lscale,k) / SQRT( em(k) )
!               tau_3d(i,k,j) = MAX( MIN( tau_3d(i,k,j), taumax ), taumin )
               tau_3d(i,k,j) = MIN( tau_3d(i,k,j), taumax )
            ENDIF

         ENDDO

#ifdef STATS

         ! Only initialize stats one time; so it needs to be done for only
         ! one column, chose central or a border-central column.
         x_center = INT( ( i_end + i_start )/2.0 )
         y_center = INT( ( j_end + j_start )/2.0 )

         IF ( i == x_center .AND. j == y_center ) THEN

            IF ( config_flags%cu_physics > 0 ) THEN

               SELECT CASE ( config_flags%cu_physics )

               ! BOMEX case
               CASE ( BOMEXSCHEME )
                 runfile      = "bomex_hoc.in"
                 day          = 22
                 month        = 6
                 year         = 1969
                 rlat         = 15.0
                 rlon         = 0.0
                 time_initial = 0.0

               ! ARM case
               CASE ( ARMSCHEME )
                 runfile      = "arm_hoc.in"
                 day          = 21
                 month        = 6
                 year         = 1997
                 rlat         = 36.62
                 rlon         = -97.5
                 time_initial = 41400.0

               ! FIRE case
               CASE ( FIRESCHEME )
                 runfile      = "fire_hoc.in"
                 day          = 7
                 month        = 7
                 year         = 1987
                 rlat         = 43.23
                 rlon         = -125.0
                 time_initial = 0.0

               ! ATEX case
               CASE ( ATEXSCHEME )
                 runfile      = "atex_hoc.in"
                 day          = 7
                 month        = 2
                 year         = 1969
                 rlat         = 14.9406
                 rlon         = -35.0
                 time_initial = 0.0

               CASE DEFAULT
                 runfile      = "hoc.in"
                 day          = 1
                 month        = 1
                 year         = 2000
                 rlat         = 45.0
                 rlon         = 0.0
                 time_initial = 0.0

               END SELECT

            ENDIF

            time = time_initial

!           Initialize statistics output
            CALL stats_init( 50, runfile, gr%nz, gr%zt, gr%zm,          &
                             day, month, year, rlat, rlon, time, DT )

         ENDIF

#endif /*STATS*/

         ! Free memory
         CALL parameterization_cleanup( .true. )
         ! All of HOC's prognostic variable arrays have now been deallocated.
         ! They have been eliminated from memory.

      ENDDO
   ENDDO

ENDIF

END SUBROUTINE HOC_BL_INIT

!==============================================================================

END MODULE module_bl_newhoc
