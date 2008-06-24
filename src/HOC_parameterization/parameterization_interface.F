!-----------------------------------------------------------------------
! $Id: parameterization_interface.F,v 1.98 2008-06-24 18:13:58 faschinj Exp $

!  Directives for scalars;  
!  Theta and rt are arbitrary and used to test the code.

#ifdef SCALARS
#define SCLR_THETA 1
#define SCLR_RT 2
#endif

!-----------------------------------------------------------------------
        module hoc_parameterization_interface

!       Description:
!       The module containing the `core' of the HOC model.

!       References:
!       None
!-----------------------------------------------------------------------

        implicit none

        public :: 
     .  parameterization_setup,
     .  parameterization_timestep,
     .  parameterization_cleanup

        private :: latin_hypercube_sampling

        private ! Default Scope

        contains

!-----------------------------------------------------------------------
       subroutine parameterization_timestep
     .            ( iter, dt, fcor,
     .              thlm_forcing, rtm_forcing, wmm, wmt,
     .              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .              p, rhom, rhot, exner,
     .              um, vm, upwp, vpwp, up2, vp2,
     .              thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .              rtp2, thlp2, rtpthlp,
     .              Scm, taum, rcm, cf,
     .              err_code, implemented
#ifdef SCALARS
     .              ,wpsclrp_sfc, wpedsclrp_sfc,   ! Optional
     .              sclrm, sclrm_forcing, edsclrm, ! Optional
     .              wpsclrp
#endif
     .            )

!      Description:
!      Subroutine to advance the model one timestep

!      References:
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!-----------------------------------------------------------------------

       ! Modules to be included
       use constants, only:
     .     wtol, ! Variable(s)
     .     emin,
     .     thltol,
     .     rttol,
     .     ep2,
     .     Cp,
     .     Lv,
     .     ep1,
     .     eps,
     .     sclr_dim
       use parameters, only:
     .     gamma_coefc, ! Variable(s)
     .     gamma_coefb,
     .     gamma_coef,
     .     T0,
     .     taumax,
     .     c_K,
     .     ts_nudge
       use model_flags, only:
     .     LH_on, ! Variable(s)
     .     lKhm_aniso,
     .     luv_nudge
       use grid_class, only:
     .     gr, ! Variable(s)
     .     zm2zt, ! Procedure(s)
     .     zt2zm,
     .     ddzm
       use numerical_check, only:
     .     parameterization_check ! Procedure(s)
       use diagnostic_variables, only:
     .     Skwt, ! Varible(s)
     .     Skwm,
     .     Sct,
     .     wp4,
     .     wpthvp,
     .     thlpthvp,
     .     rtpthvp,
     .     wprcp,
     .     rtprcp,
     .     thlprcp,
     .     pdf_parms,
     .     rcp2,
     .     rsm,
     .     shear,
     .     ustar,
     .     Kht,
     .     wprtp2,
     .     wp2rtp,
     .     wpthlp2,
     .     wp2thlp,
     .     wprtpthlp,
     .     wp2thvp,
     .     wp2rcp,
     .     thvm,
     .     em,
     .     Lscale,
     .     taut,
     .     Khm,
     .     umt,
     .     vmt,
     .     vg,
     .     ug,
     .     um_ref,
     .     vm_ref,
     .     wp2zt,
     .     thlp2zt,
     .     wpthlpzt,
     .     wprtpzt,
     .     rtp2zt,
     .     rtpthlpzt
#ifdef SCALARS
     .    ,edsclrmt,
     .     wpedsclrp,
     .     sclrpthvp,  ! sclr'th_v'
     .     sclrprtp,   ! sclr'rt'
     .     sclrp2,     ! sclr'^2
     .     sclrpthlp,  ! sclr'th_l'
     .     sclrprcp,   ! sclr'rc'
     .     wp2sclrp,   ! w'^2 sclr'
     .     wpsclrp2,   ! w'sclr'^2
     .     wpsclrprtp, ! w'sclr'rt'
     .     wpsclrpthlp ! w'sclr'thl'
#endif

       use mixing, only:
           ! Variable(s) 
     .     timestep_mixing          ! Compute mean/flux terms
       use diagnose_variances, only:
           ! Variable(s) 
     .     diag_var     ! Computes variance terms
       use surface_var, only: 
     .     sfc_var ! Procedure
       use pdf_closure, only: 
           ! Procedure 
     .     pdf_closure_new     ! Prob. density function
       use mixing_length, only:
     .     compute_length ! Procedure
       use compute_um_edsclrm_mod, only: 
     .     compute_um_edsclrm, ! Procedure(s)
     .     compute_uv_tndcy
       use saturation, only: 
           ! Procedure
     .     rsat ! Saturation mixing ratio
       use wp23, only: 
     .     timestep_wp23 ! Procedure
       use stats_prec, only: 
     .     time_prec ! Variable(s)
       use error_code, only : 
     .     clubb_var_equals_NaN,! Variable(s)
     .     lapack_error, ! Procedure(s)
     .     clubb_at_debug_level
       use Skw, only: 
     .     Skw_func ! Procedure
       use permute_height_time_mod, only: 
     .     permute_height_time ! Procedure
#ifdef STATS
       use stats_variables, only:
     .     zm, ! Variable(s)
     .     sfc,
     .     lstats_samp,
     .     ish,
     .     ilh,
     .     iustar,
     .     iup2_plus_vp2
       use stats_subs, only:
     .     stats_accumulate ! Procedure
       use stats_type, only:
     .     stat_update_var, stat_update_var_pt ! Procedure(s)
#endif /*STATS*/

       implicit none

       intrinsic :: sqrt, min, max, exp, mod

       ! Input
       integer, intent(in) :: 
     . iter      ! Closure iteration number

       logical, intent(in) :: 
     . implemented ! Is this part of a larger host model (T/F) ?

       ! Note on dt, dmain, and dtclosure: since being moved out of
       ! hoc.F, all subroutines within parameterization_timestep now use
       ! dt for time dependent calculations.  The old dt is noted in
       ! each section of the code -dschanen 20 April 2006 
       real(kind=time_prec), intent(in) :: 
     . dt            ! Current timestep size    [s]

       real, intent(in) :: 
     . fcor          ! Coriolis forcing         [s^-1]

       real, intent(in), dimension(gr%nnzp) :: 
     . thlm_forcing,  ! theta_l forcing.        [K/s]
     . rtm_forcing,   ! r_t forcing.            [(kg/kg)/s] 
     . wmm,           ! wm on moment. grid.     [m/s]
     . wmt,           ! wm on thermo. grid.     [m/s]
     . p,             ! Pressure.               [Pa] 
     . rhom,          ! Density on moment. grid [kg/m^3]
     . rhot,          ! Density on thermo. grid [kg/m^3] 
     . exner          ! Exner function.         [-]

       real, intent(in) :: 
     . wpthlp_sfc,  ! w' theta_l' at surface.   [(m K)/s]
     . wprtp_sfc,   ! w' r_t' at surface.       [(kg m)/( kg s)]
     . upwp_sfc,    ! u'w' at surface.          [m^2/s^2]
     . vpwp_sfc     ! v'w' at surface.          [m^2/s^2]

       ! Input/Output
       ! These are prognostic or are planned to be in the future
       real, intent(inout), dimension(gr%nnzp) :: 
     . um,      ! u wind.                       [m/s]
     . upwp,    ! u'w'.                         [m^2/s^2]
     . vm,      ! v wind.                       [m/s]
     . vpwp,    ! u'w'.                         [m^2/s^2]
     . up2,     ! u'^2                          [m^2/s^2]
     . vp2,     ! v'^2                          [m^2/s^2]
     . rtm,     ! r_t Total water mixing ratio. [kg/kg]
     . wprtp,   ! w' r_t'.                      [(m kg)/(s kg)]
     . thlm,    ! th_l Liquid potential temp.   [K]
     . wpthlp,  ! w' th_l'.                     [(m K)/s]
     . wp2,     ! w'^2.                         [m^2/s^2]
     . wp3,     ! w'^3.                         [m^3/s^3]
     . Scm,     ! Sc on moment. grid.           [-]
     . rtp2,    ! r_t'^2.                       [(kg/kg)^2]
     . thlp2,   ! th_l'^2.                      [K^2]
     . rtpthlp, ! r_t' th_l'.                   [(kg K)/kg]
     . taum,    ! Tau on moment. grid.          [s]
     . rcm      ! Liquid water mixing ratio.    [kg/kg]

       ! Needed for output for host models
       real, intent(inout), dimension(gr%nnzp) :: 
     . cf ! Cloud fraction.     [%]

       ! Diagnostic, for if some calculation goes amiss.
       integer, intent(inout) :: err_code
     
#ifdef SCALARS

       ! Optional Input Variables
       real, intent(in),  dimension(sclr_dim) :: 
     . wpsclrp_sfc,  ! Scalar flux at surface           [units m/s]
     . wpedsclrp_sfc ! Eddy-Scalar flux at surface      [units m/s]

       ! Optional Input/Output Variables
       real, intent(inout), dimension(gr%nnzp,sclr_dim) ::
     . sclrm,          ! Passive scalar mean.           [units vary]
     . sclrm_forcing,  ! Passive scalar forcing.        [units vary/s]
     . edsclrm,        ! Eddy passive scalar mean.      [units vary]
     . wpsclrp         ! w'sclr'                        [units vary m/s]

#endif /*SCALARS*/

       ! Local Variables
       integer :: i, j, k, km1, kp1

       logical :: ltmp

       real, dimension(gr%nnzp) ::
     . tmp1, gamma_Skw_fnc
 
#ifdef SCALARS
       real, dimension(gr%nnzp,sclr_dim) ::
     . sclr_tmp1, sclr_tmp2, sclr_tmp3, sclr_tmp4 ! for PDF closure
#endif /* SCALARS */

!------- Local variables for Latin Hypercube sampling ------------------

       integer i_rmd 

! Number of variables to sample
       integer, parameter :: d_variables = 5

! n = number of calls to microphysics per timestep (normally=2)
       integer, parameter :: n_micro_call = 12

! sequence_length = nt/n = number of timesteps before sequence repeats.
       integer, parameter :: sequence_length = 1

! nt = number of random samples before sequence of repeats (normally=10)
       integer, parameter :: nt_repeat = n_micro_call * sequence_length

! A true/false flag that determines whether
!     the PDF allows us to construct a sample
       logical sample_flag

       integer, dimension(gr%nnzp, nt_repeat, d_variables+1) 
     . :: p_height_time ! matrix of rand ints

! coeffs of s from pdf_closure_new
       real :: crt1, crt2, cthl1, cthl2   

!-------- End Latin hypercube section ----------------------------------


!-------- Test input variables ----------------------------------------
       if ( clubb_at_debug_level( 2 ) ) then
       call parameterization_check
     .      ( thlm_forcing, rtm_forcing, wmm, wmt, p, rhom,
     .        rhot, exner, wpthlp_sfc, wprtp_sfc,
     .        upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,
     .        up2, vp2, rtm, wprtp, thlm,
     .        wpthlp, wp2, wp3, Scm, rtp2, thlp2,
     .        rtpthlp, taum, rcm, cf, "beginning of " 
#ifdef SCALARS
     .        ,wpsclrp_sfc, wpedsclrp_sfc,
     .        sclrm, sclrm_forcing, edsclrm
#endif /*SCALARS*/
     .       )
       end if
!-----------------------------------------------------------------------

       !----------------------------------------------------------------
       ! Interpolate wp2 & wp3, and then compute Skw for m & t grid
       !----------------------------------------------------------------

       do k = 1, gr%nnzp, 1

         Skwt(k) = Skw_func( zm2zt(wp2,k), wp3(k), wtol )
         Skwm(k) = Skw_func( wp2(k), zt2zm(wp3,k), wtol )

       enddo

       ! SET SURFACE VALUES OF FLUXES (BROUGHT IN)
       wpthlp(1) = wpthlp_sfc
       wprtp(1)  = wprtp_sfc
       upwp(1)   = upwp_sfc
       vpwp(1)   = vpwp_sfc

#ifdef SCALARS
       ! Set fluxes for passive scalars (if enabled)
       wpsclrp(1,1:sclr_dim)   = wpsclrp_sfc(1:sclr_dim)
       wpedsclrp(1,1:sclr_dim) = wpedsclrp_sfc(1:sclr_dim)
#endif
       !----------------------------------------------------------------
       ! Set Surface variances
       !----------------------------------------------------------------

!      Surface variances should be set here, before the call to diag_var.
!      The reasons that surface variances can be set here are because the
!      only variables that are the input into surface variances are the
!      surface values of wpthlp, wprtp, upwp, and vpwp.  The surface values
!      of all those variables are set in the surface forcings section of the
!      GCSS cases listed in the main timestep above.  Even if they weren't
!      set there, the updates to wpthlp, wprtp, upwp, and vpwp are at the
!      end of the closure loop, right before the code loops back around to
!      this point at the top of the closure loop.
!      Surface variances need to be set here for two reasons.  One reason is
!      that the values of rtp2, thlp2, and rtpthlp at the surface will be
!      used to find the diffusional term and the mean advection term in each
!      predictive equation for those respective terms.  The other reason is
!      that if the correct surface variances are not set here and diag_var
!      outputs it's own value for them, it will results in a faulty value for
!      Scm at the surface.  Brian Griffin.  December 18, 2005.

!      Surface effects should not be included with any case where the lowest
!      level is not the ground level.  Brian Griffin.  December 22, 2005.
       IF ( gr%zm(1) == 0.0 ) THEN
          call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1), 
     .                  ustar, wp2(1), up2(1), vp2(1), 
     .                  thlp2(1), rtp2(1), rtpthlp(1), err_code
#ifdef SCALARS
     .                 ,wpsclrp(1,1:sclr_dim), sclrp2(1,1:sclr_dim),
     .                  sclrprtp(1,1:sclr_dim), sclrpthlp(1,1:sclr_dim)
#endif /* SCALARS */
     .                )
          ! Subroutine may produce NaN values, and if so, exit
          ! gracefully.
          ! Joshua Fasching March 2008
          if( err_code == clubb_var_equals_NaN ) return
    
       ELSE
          ! Variances for cases where the lowest level is not at the surface.
          ! Eliminate surface effects on lowest level variances.
          wp2(1)     = (2.0/3.0) * emin
          up2(1)     = (2.0/3.0) * emin
          vp2(1)     = (2.0/3.0) * emin
          thlp2(1)   = 0.0
          rtp2(1)    = 0.0
          rtpthlp(1) = 0.0
#ifdef SCALARS
          DO i = 1, sclr_dim, 1
             sclrp2(1,i)    = 0.0
             sclrprtp(1,i)  = 0.0
             sclrpthlp(1,i) = 0.0
          END DO
#endif /* SCALARS */
       END IF

       !----------------------------------------------------------------
       ! Diagnose variances
       !----------------------------------------------------------------

       ! We also found that certain cases require a time tendency to run
       ! at shorter timesteps.
       ! This requires us to store in memory Scm and taum between timesteps.

       ! We found that if we call diag_var first, we can use a longer timestep.
       call diag_var( taum, wmm, rtm, wprtp, 
     .                thlm, wpthlp, wpthvp, um, vm,
     .                wp2, wp3, upwp, vpwp, Scm, Skwm, Kht,
     .                rtp2, thlp2, rtpthlp,
     .                up2, vp2,
! Vince Larson used prognostic timestepping of variances 
!    in order to increase numerical stability.  17 Jul 2007
!     .                .false., dt, isValid
     .                .true., dt, err_code
! End Vince Larson's change
#ifdef SCALARS
     .                ,sclrm, wpsclrp, sclrp2,    ! optional vars
     .                sclrprtp, sclrpthlp         !   "       "
#endif /* SCALARS */
     .               )

       ! Check stability
       ! Changed from a logical flag to an integer indicating nature of
       ! error.
       ! Joshua Fasching March 2008
       if ( lapack_error(err_code)) return

#ifndef DISABLE_GAMMA_SKW
       !----------------------------------------------------------------
       ! Compute gamma as a function of Skw  - 14 April 06 dschanen
       !----------------------------------------------------------------
       ! Use -DDISABLE_GAMMA_SKW with the fortran preprocessor
       ! to bring back the use of gamma_coef the value of gamma 
       ! without Skw dependence

        gamma_Skw_fnc 
     .  = gamma_coefb + (gamma_coef-gamma_coefb)
     .    *exp( -(1.0/2.0) * (Skwm/gamma_coefc)**2 )

#else
        gamma_Skw_fnc = gamma_coef

#endif /* DISABLE */

       !----------------------------------------------------------------
       ! Compute Sc with new formula from Vince
       !----------------------------------------------------------------

       Scm = gamma_Skw_fnc * ( 1.0 - 
     .       min(
     .       max( ( wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) 
     .              + 0.01 * wtol * thltol ) )**2,
     .            ( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) 
     .              + 0.01 * wtol * rttol ) )**2 ),
     .       1.0      ) )

       Sct = max( zm2zt( Scm ), 0.0 )   ! Positive definite quantity

!    Latin hypercube sample generation
!    Generate p_height_time, an nnzp x nt_repeat x d_variables array of random integers
       if ( LH_on ) then
         i_rmd = mod( iter-1, sequence_length )
         if ( i_rmd == 0) then
           call permute_height_time( gr%nnzp, nt_repeat, 
     .                               d_variables+1, p_height_time )
         end if
       end if
!    End Latin hypercube sample generation

!       print*, 'hoc.F: i_rmd=', i_rmd

       !----------------------------------------------------------------
       ! Call closure scheme
       !----------------------------------------------------------------

#ifdef SCALARS
       ! Put passive scalar input on the t grid for the PDF
       do i = 1, sclr_dim
          sclr_tmp1(:,i) = zm2zt( wpsclrp(:,i) ) 
          sclr_tmp2(:,i) = zm2zt( sclrprtp(:,i) ) 
          sclr_tmp3(:,i) = max( zm2zt( sclrp2(:,i) ), 0.0 )   ! Pos. def. quantity
          sclr_tmp4(:,i) = zm2zt( sclrpthlp(:,i) ) 
       end do ! i = 1, sclr_dim
#endif /* SCALARS */

       do k = 2, gr%nnzp, 1
         call pdf_closure_new
     .   ( p(k), exner(k), wmt(k), zm2zt(wp2, k), wp3(k), Sct(k),
     .     rtm(k), zm2zt(rtp2, k), zm2zt( wprtp, k ),
     .     thlm(k), zm2zt( thlp2, k ), zm2zt( wpthlp, k ),
     .     zm2zt(rtpthlp, k), wp4(k), wprtp2(k), wp2rtp(k),
     .     wpthlp2(k), wp2thlp(k), wprtpthlp(k),
     .     cf(k), rcm(k), wpthvp(k), wp2thvp(k), rtpthvp(k), 
     .     thlpthvp(k), wprcp(k), wp2rcp(k), rtprcp(k), thlprcp(k),
     .     rcp2(k), pdf_parms(k, :), crt1, crt2, cthl1, cthl2, err_code
#ifdef SCALARS
     .     ,sclrm(k,:), sclr_tmp1(k,:),
     .     wpsclrprtp(k,:), wpsclrp2(k,:), sclr_tmp2(k,:), 
     .     sclr_tmp3(k,:), sclr_tmp4(k,:), sclrpthvp(k,:), 
     .     wpsclrpthlp(k,:), sclrprcp(k,:), wp2sclrp(k,:)
#endif /*SCALARS*/
     .   )

        ! Subroutine may produce NaN values, and if so, exit
        ! gracefully.
        ! Joshua Fasching March 2008
         
        if( err_code == clubb_var_equals_NaN ) then
            print *, "At grid level = ",k
            return
         endif
         
         !--------------------------------------------------------------
         ! Latin hypercube sampling
         !--------------------------------------------------------------
         if ( LH_on ) then 
!          call latin_hypercube_sampling
!    .     ( k, n_micro_call, d_variables, nt_repeat, i_rmd, 
!    .       crt1, crt2, cthl1, cthl2, hydromet(:,1), 
!    .       cf, gr%nnzp, sample_flag, p_height_time )
         end if

       end do ! k = 2, nz-1

!            print*, 'hoc.F: AKm=', AKm
!            print*, 'hoc.F: AKm_est=', AKm_est

!      Interpolate momentum variables back to momentum grid.
!      Since top momentum level is higher than top thermo level,
!      set variables at top momentum level to 0.

       wp4               = max( zt2zm( wp4 ), 0.0 )   ! Pos. def. quantity
       wp4(gr%nnzp)      = 0.0
       wpthvp            = zt2zm( wpthvp )
       wpthvp(gr%nnzp)   = 0.0
       thlpthvp          = zt2zm( thlpthvp )
       thlpthvp(gr%nnzp) = 0.0
       rtpthvp           = zt2zm( rtpthvp )
       rtpthvp(gr%nnzp)  = 0.0
       wprcp             = zt2zm( wprcp )
       wprcp(gr%nnzp)    = 0.0
       rtprcp            = zt2zm( rtprcp )
       rtprcp(gr%nnzp)   = 0.0
       thlprcp           = zt2zm( thlprcp )
       thlprcp(gr%nnzp)  = 0.0
       rcp2              = max( zt2zm( rcp2 ), 0.0 )   ! Pos. def. quantity
       rcp2(gr%nnzp)     = 0.0
#ifdef SCALARS
!      Interpolate passive scalars back onto the m grid
       do i = 1, sclr_dim
          sclrpthvp(:,i)       = zt2zm( sclrpthvp(:,i) )
          sclrpthvp(gr%nnzp,i) = 0.0
          sclrprcp(:,i)        = zt2zm( sclrprcp(:,i) )
          sclrprcp(gr%nnzp,i)  = 0.0
       end do ! i=1, sclr_dim
#endif /* SCALARS */

       !----------------------------------------------------------------
       ! Compute thvm
       !----------------------------------------------------------------

       thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

       !----------------------------------------------------------------
       ! Compute tke
       !----------------------------------------------------------------

       if ( .not. lKhm_aniso ) then
         ! tke is assumed to be 3/2 of wp2
         em = 1.5 * wp2
       else
         em = 0.5 * ( wp2 + vp2 + up2 )
       end if

       !----------------------------------------------------------------
       ! Compute mixing length
       !----------------------------------------------------------------

       call compute_length
     .      ( thvm, thlm, rtm, rcm, em, p, exner, ! Input
     .        Lscale, err_code )  ! Output
       
       ! Subroutine may produce NaN values, and if so, exit
       ! gracefully.
       ! Joshua Fasching March 2008       
       if( err_code == clubb_var_equals_NaN ) return


       !----------------------------------------------------------------
       ! Dissipation time
       !----------------------------------------------------------------
! Vince Larson replaced the cutoff of emin by wtol**2.  7 Jul 2007
!     This is to prevent tau from being too large (producing little damping)
!     in stably stratified layers with little turbulence.
!       tmp1 = SQRT( MAX( emin, zm2zt( em ) ) )
!       taut = MIN( Lscale / tmp1, taumax )
!       taum 
!     . = MIN( ( zt2zm( Lscale ) / SQRT( MAX( emin, em ) ) ), taumax )
       tmp1 = SQRT( MAX( wtol**2, zm2zt( em ) ) )
       taut = MIN( Lscale / tmp1, taumax )
       taum 
     . = MIN( ( MAX( zt2zm( Lscale ), 0.0 ) 
     .         / SQRT( MAX( wtol**2, em ) ) ), taumax )
! End Vince Larson's replacement.

       ! Modification to damp noise in stable region
! Vince Larson commented out because it may prevent turbulence from 
!    initiating in unstable regions.  7 Jul 2007
!       do k = 1, gr%nnzp
!         if ( wp2(k) <= 0.005 ) then
!           taut(k) = taumin
!           taum(k) = taumin
!         end if
!       end do
! End Vince Larson's commenting.

       !----------------------------------------------------------------
       ! Eddy diffusivity coefficient
       !----------------------------------------------------------------
       ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)

       Kht = c_K * Lscale * tmp1
       Khm = c_K * max( zt2zm( Lscale ), 0.0 ) 
     .           * sqrt( max( em, emin ) )
 
!#######################################################################
!############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
!#######################################################################


       ! Store the saturation mixing ratio for output purposes.  Brian
       do k = 1, gr%nnzp, 1
         rsm(k) =
     .   rsat( p(k), (thlm(k)+(Lv/(Cp*exner(k)))*rcm(k))*exner(k) )
       end do

       !----------------------------------------------------------------
       ! Advance rtm/wprtp and thlm/wpthlp one time step
       !----------------------------------------------------------------

        call timestep_mixing( dt, Scm, wmm, wmt, wp2, wp3,
     .                        Kht, taum, Skwm, rtpthvp,
     .                        rtm_forcing, thlpthvp,
     .                        thlm_forcing, rtp2, thlp2,
     .                        rtm, wprtp, thlm, wpthlp,
     .                        implemented, err_code
#ifdef SCALARS
     .                        ,sclrpthvp, sclrm_forcing    ! optional
     .                        ,sclrp2, sclrm, wpsclrp      ! optional 
#endif /* SCALARS */
     .                      )

       ! Wrapped LAPACK procedures may report errors, and if so, exit
       ! gracefully.
       ! Joshua Fasching March 2008

       if ( lapack_error( err_code ) ) return

        ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
        ! This code won't work unless rtm >= 0 !!!
        do k = 1, gr%nnzp
           if ( rtm(k) < rcm(k) ) then
              if ( clubb_at_debug_level(1) ) then
              print*, 'rtm < rcm in timestep_mixing at k=', k, '.',
     .                '  Clipping rcm.'
              end if
              rcm(k) = max( 0.0, rtm(k) - eps )
           end if
        end do

       !----------------------------------------------------------------
       ! Advance wp2/wp3 one timestep
       !----------------------------------------------------------------

       call timestep_wp23
     .      ( dt, Scm, wmm, wmt, wpthvp, wp2thvp,
     .        um, vm, upwp, vpwp, up2, vp2, Khm, Kht,
     .        taum, taut, Skwm, Skwt, pdf_parms(:, 13),
     .        wp2, wp3, err_code )
 
       ! Wrapped LAPACK procedures may report errors, and if so, exit
       ! gracefully.
       ! Joshua Fasching March 2008
       if ( lapack_error( err_code ) ) return

#ifdef SCALARS
       !----------------------------------------------------------------
       ! Compute Eddy-diff. Passive Scalars
       !----------------------------------------------------------------

       do i=1, sclr_dim

         edsclrmt(1:gr%nnzp,i) = - wmt * ddzm( zt2zm( edsclrm(:,i) ) )

         call compute_um_edsclrm( "edsclr", wpedsclrp(1,i), 
     .                            edsclrmt(:,i), Khm, dt, edsclrm(:,i), 
     .                            wpedsclrp(:,i), err_code )


       end do

       ! Set boundary condition as in rt
       edsclrm(1,1:sclr_dim) = edsclrm(2,1:sclr_dim)
       if ( lapack_error(err_code) ) return
#endif /*SCALARS*/


       !----------------------------------------------------------------
       ! Update winds
       !----------------------------------------------------------------

       call compute_uv_tndcy
     .    ( "um", um, wmt, fcor, vm, vg, implemented, umt )

       ! dtmain to dt -dschanen
       call compute_um_edsclrm( "um", upwp(1), umt,
     .                          Khm, dt, um, upwp, err_code )

       um(1)       = ( ( um(3)-um(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + um(2)
       um(gr%nnzp) = ( ( um(gr%nnzp-1)-um(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + um(gr%nnzp-1)
       
       ! Wrapped LAPACK procedures may report errors, and if so, exit
       ! gracefully.
       ! Joshua Fasching March 2008
       if ( lapack_error(err_code) ) return

       call compute_uv_tndcy
     .    ("vm", vm, wmt, fcor, um, ug, implemented, vmt )

       call compute_um_edsclrm( "vm", vpwp(1), vmt,
     .                          Khm, dt, vm, vpwp, err_code )

       vm(1)       = ( ( vm(3)-vm(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + vm(2)
       vm(gr%nnzp) = ( ( vm(gr%nnzp-1)-vm(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + vm(gr%nnzp-1)
       
       ! Wrapped LAPACK procedures may report errors, and if so, exit
       ! gracefully.
       ! Joshua Fasching March 2008    
       if ( lapack_error( err_code ) ) return

       if ( luv_nudge ) then
         um(1:gr%nnzp) = um(1:gr%nnzp) 
     .      - ((um(1:gr%nnzp) - um_ref(1:gr%nnzp)) * (dt/ts_nudge))
         vm(1:gr%nnzp) = vm(1:gr%nnzp) 
     .      - ((vm(1:gr%nnzp) - vm_ref(1:gr%nnzp)) * (dt/ts_nudge))
       end if

       ! Compute Shear Production  -Brian
       do k = 1, gr%nnzp-1, 1
         shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k)
     .              -vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
       end do
       shear(gr%nnzp) = 0.0

!#######################################################################
!#############            ACCUMULATE STATISTICS            #############
!#######################################################################

#ifdef STATS
!      Store values of surface fluxes for statistics
       if (lstats_samp) then
         call stat_update_var_pt( ish, 1, wpthlp(1)*rhom(1)*Cp, sfc )

         call stat_update_var_pt( ilh, 1, wprtp(1)*rhom(1)*Lv, sfc )

         call stat_update_var_pt( iustar, 1, ustar, sfc )

         ! Michael Falk added this for up2_plus_vp2, 23 Jan 2007.
         call stat_update_var( iup2_plus_vp2, up2+vp2, zm )
         ! eMFc
       end if 

       ! Added to allow tuning without using the zm stats file
       wp2zt     = max( zm2zt( wp2 ), 0.0 )   ! Positive definite quantity
       thlp2zt   = max( zm2zt( thlp2 ), 0.0 )   ! Positive definite quantity
       wpthlpzt  = zm2zt( wpthlp )
       wprtpzt   = zm2zt( wprtp )
       rtp2zt    = max( zm2zt( rtp2 ), 0.0 )   ! Positive definite quantity
       rtpthlpzt = zm2zt( rtpthlp )

       call stats_accumulate
     .           ( um, vm, upwp, vpwp, up2, vp2, thlm,
     .             rtm, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp,
     .             p, exner, rhot, rhom, thlm_forcing, rtm_forcing, wmm,
     .             wmt, Scm, taum, rcm, cf
#ifdef SCALARS
     .             ,sclrm, edsclrm, sclrm_forcing, wpsclrp
#endif /*SCALARS*/
     .            )

#endif /*STATS*/

       if( clubb_at_debug_level( 2 ) ) then
       call parameterization_check
     .      ( thlm_forcing, rtm_forcing, wmm, wmt, p, rhom,
     .        rhot, exner, wpthlp_sfc, wprtp_sfc,
     .        upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,
     .        up2, vp2, rtm, wprtp, thlm,
     .        wpthlp, wp2, wp3, Scm, rtp2, thlp2,
     .        rtpthlp, taum, rcm, cf, "end of "
#ifdef SCALARS
     .        ,wpsclrp_sfc, wpedsclrp_sfc,
     .        sclrm, sclrm_forcing, edsclrm
#endif /*SCALARS*/
     . )
       endif
!-----------------------------------------------------------------------

       return
       end subroutine parameterization_timestep


!-----------------------------------------------------------------------
        subroutine latin_hypercube_sampling
     .             ( k, n, dvar, nt, i_rmd,
     .               crt1, crt2, cthl1, cthl2,
     .               rrm, cf, grid, sflag, p_height_time )
!       Description:
!       Estimate using Latin Hypercubes.  This is usually disabled by default.
!       The actual generation of a random matrix is done in a call from the
!       subroutine hoc_initialize to permute_height_time()
!       References:
!-----------------------------------------------------------------------

        use diagnostic_variables, only: 
     .      pdf_parms, ! Variable(s) 
     .      AKm_est, 
     .      AKm,
     .      AKstd,
     .      AKstd_cld,
     .      AKm_rcm,
     .      AKm_rcc,
     .      rcm_est

        use lh_sampler_mod, only:
     .      lh_sampler ! Procedure
        use micro_calcs_mod, only:
     .      micro_calcs ! Procedure
        
        implicit none  


        ! Input Variables 
        integer, intent(in) :: k  ! index
        integer, intent(in) :: n, dvar, i_rmd, nt, grid
        logical, intent(out) :: sflag

        ! coeffs of s from pdf_closure_new
        real, intent(in) :: crt1, crt2, cthl1, cthl2

        real, dimension(grid), intent(in) :: 
     .  rrm, ! Rain water mixing ratio  [kg/kg]
     .  cf   ! Cloud fraction           [%]

        integer, dimension(1:grid, 1:nt, 1:(dvar+1) ), intent(in) ::
     .  p_height_time ! matrix of rand ints

        ! Local Variables

        integer :: p_matrix(n, dvar+1)
        ! Sample drawn from uniform distribution
        double precision, dimension(1:n,1:(dvar+1)) :: X_u

        ! Sample that is transformed ultimately to normal-lognormal
        double precision, dimension(1:n,1:dvar) :: X_nl

        ! Choose which rows of LH sample to feed into closure.
        p_matrix(1:n,1:(dvar+1)) =
     .  p_height_time( k,n*i_rmd+1:n*i_rmd+n, 1:(dvar+1) )

!       print*, 'hoc.F: got past p_matrix'

        ! Generate LH sample, represented by X_u and X_nl, for level k
        call lh_sampler( n, nt, dvar, p_matrix,
     .                   cf(k), pdf_parms(k, :),
     .                   crt1, crt2, cthl1, cthl2,
     .                   rrm(k), X_u, X_nl, sflag )

!       print *, 'hoc.F: got past lh_sampler'

        ! Perform LH and analytic microphysical calculations
        call micro_calcs( n, dvar, X_u, X_nl, sflag,
     .                    pdf_parms(k,:),
     .                    AKm_est(k), AKm(k), AKstd(k), AKstd_cld(k),
     .                    AKm_rcm(k), AKm_rcc(k), rcm_est(k) )

!       print*, 'k, AKm_est=', k, AKm_est(k)
!       print*, 'k, AKm=', k, AKm(k)

        return
        end subroutine latin_hypercube_sampling
!-----------------------------------------------------------------------
        subroutine parameterization_setup
     .             ( nzmax, T0_in, ts_nudge_in, params,
     .               lbugsrad, kk_rain, licedfs, lcoamps_micro,
     .               cloud_sed, luv_nudge, lKhm_aniso, 
     .               implemented, grid_type, deltaz, zm_init,
     .               momentum_heights, thermodynamic_heights, 
     .               host_dx, host_dy, err_code )

        use grid_class, only:
     .      gridsetup ! Procedure
        use param_index, only: 
     .      nparams ! Variable(s)
        use parameters, only:
     .      setup_parameters ! Procedure
        use diagnostic_variables, only:
     .      setup_diagnostic_variables ! Procedure
        use prognostic_variables, only:
     .      setup_prognostic_variables ! Procedure
        use constants, only: 
     .      fstderr, ! Variable(s)
     .      Lscale_max
        use error_code, only: 
     .      clubb_var_out_of_bounds ! Variable(s)
        implicit none

        ! External
        external :: setup_model_flags

        ! Input

        ! Grid definition
        integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

        ! Flag to see if CLUBB is running on it's own,
        ! or if it's implemented as part of a host model.
        logical, intent(in) :: implemented   ! (T/F)

        ! If CLUBB is running on it's own, this option determines
        ! if it is using:
        ! 1) an evenly-spaced grid,
        ! 2) a stretched (unevenly-spaced) grid entered on the
        !    thermodynamic grid levels (with momentum levels set
        !    halfway between thermodynamic levels), or
        ! 3) a stretched (unevenly-spaced) grid entered on the
        !    momentum grid levels (with thermodynamic levels set
        !    halfway between momentum levels).
        integer, intent(in) :: grid_type

        ! If the CLUBB model is running by itself, and is using an
        ! evenly-spaced grid (grid_type = 1), it needs the vertical
        ! grid spacing and momentum-level starting altitude as input.
        real, intent(in) ::
     .  deltaz,  ! Change in altitude per level           [m]
     .  zm_init  ! Initial grid altitude (momentum level) [m]

        ! If the CLUBB parameterization is implemented in a host model,
        ! it needs to use the host model's momentum level altitudes
        ! and thermodynamic level altitudes.
        ! If the CLUBB model is running by itself, but is using a
        ! stretched grid entered on thermodynamic levels (grid_type = 2),
        ! it needs to use the thermodynamic level altitudes as input.
        ! If the CLUBB model is running by itself, but is using a
        ! stretched grid entered on momentum levels (grid_type = 3),
        ! it needs to use the momentum level altitudes as input.
        real, intent(in), dimension(nzmax) ::
     .  momentum_heights,     ! Momentum level altitudes (input)      [m]
     .  thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

        ! Host model horizontal grid spacing, if part of host model.
        real, intent(in) ::
     .  host_dx, ! East-West horizontal grid spacing     [m]
     .  host_dy  ! North-South horizontal grid spacing   [m]

        ! Model parameters
        real, intent(in) :: 
     .  T0_in, ts_nudge_in

        real, intent(in), dimension(nparams) ::
     .  params  ! Including C1, nu1, nu2, etc.

        ! Flags
        logical, intent(in) :: 
     .  lbugsrad,     ! BUGSrad interactive radiation scheme
     .  kk_rain,      ! K & K rain microphysics
     .  licedfs,      ! Simplified ice scheme
     .  lcoamps_micro,! COAMPS microphysics scheme
     .  cloud_sed,    ! Cloud Sedimentation
     .  luv_nudge,    ! Wind nudging for mpace_b case
     .  lKhm_aniso    ! Whether to use anisotropic Khm. - Michael Falk 2 Feb 2007

        ! Output variables
        integer, intent(out) ::
     .  err_code   ! Diagnostic for a problem with the setup

        ! Setup flags

        call setup_model_flags
     .       ( lbugsrad, kk_rain, cloud_sed,
     .         licedfs, lcoamps_micro, luv_nudge,
     .         lKhm_aniso )

        ! Define model constant parameters

        call setup_parameters
     .       ( deltaz, T0_in, ts_nudge_in, params, err_code )

!       Error Report
!       Joshua Fasching February 2008
        if ( err_code == clubb_var_out_of_bounds ) then
                
           write(fstderr,*) "Error in parameterization_setup"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "deltaz = ", deltaz
           write(fstderr,*) "zm_init = ", zm_init
           write(fstderr,*) "momentum_heights = ", momentum_heights
           write(fstderr,*) "thermodynamic_heights = ", 
     .                       thermodynamic_heights
           write(fstderr,*) "T0_in = ", T0_in
           write(fstderr,*) "ts_nudge_in = ", ts_nudge_in
           write(fstderr,*) "params = ", params 

           return

        endif

!        if ( .not. implemented ) then
!          call setup_diagnostic_variables( nzmax )
!        end if

        ! Both prognostic variables and diagnostic variables need to be
        ! declared, allocated, initialized, and deallocated whether HOC
        ! is part of a larger model or not.
        call setup_prognostic_variables( nzmax )
        call setup_diagnostic_variables( nzmax )

        ! Setup grid
        call gridsetup( nzmax, implemented, grid_type,
     .                  deltaz, zm_init, momentum_heights,
     .                  thermodynamic_heights )

        ! Determine the maximum allowable value for Lscale.
        if ( implemented ) then
           Lscale_max = 0.25 * min(host_dx,host_dy)
        else
           Lscale_max = 1.0e5
        endif

        return
        end subroutine parameterization_setup

!-----------------------------------------------------------------------
        subroutine parameterization_cleanup( implemented )

        use diagnostic_variables, only:
     .      cleanup_diagnostic_variables ! Procedure
        use prognostic_variables, only:
     .      cleanup_prognostic_variables ! Procedure

        implicit none

        ! whether this is part of larger model (T/F)
        logical, intent(in) :: implemented 
  
!        if ( .not. implemented ) then
!          call cleanup_diagnostic_variables( )
!        end if

        ! Both prognostic variables and diagnostic variables need to be
        ! declared, allocated, initialized, and deallocated whether HOC
        ! is part of a larger model or not.
        call cleanup_prognostic_variables( )
        call cleanup_diagnostic_variables( )

        return
        end subroutine parameterization_cleanup

      end module hoc_parameterization_interface
