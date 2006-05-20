!-----------------------------------------------------------------------
! $Id: parameterization_interface.F,v 1.14 2006-05-20 12:51:14 griffinb Exp $


!  Module hoc_parameterization_interface

!  The `core' of the HOC model 

!  Directives for scalars;  
!  Theta and rt are arbitrary and used to test the code.

#ifdef SCALARS
#define SCLR_THETA 1
#define SCLR_RT 2
#endif

        module hoc_parameterization_interface
        implicit none

        public parameterization_setup
        public parameterization_timestep
        public parameterization_cleanup

        private latin_hypercube_sampling

        contains

!-----------------------------------------------------------------------
!    SUBROUTINE parameterization_timestep

!    Calls:
!    subroutine diag_var
!    subroutine length
!    subroutine pdf_closure_new
!    subroutine rain (when kk_rain is true)
!    subroutine timestep_mixing
!    subroutine timestep_wp23
!    subroutine update_ed_var

!    function rsat
!    function zm2zt (from grid_class)
!    function zt2zm (from grid_class)

!    If compiled with STATS:
!    subroutine stats_accumulate

!    If LH_on is true:
!    permute_height_time (external)
!    latin_hypercube_sampling (Internal to module)
!-----------------------------------------------------------------------
       subroutine parameterization_timestep
     .            ( iter, niterlong, time_initial, time_current, dt,
     .              thlm_forcing, rtm_forcing, wmm, wmt,
     .              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .              p, rhom, rhot, exner,
     .              um, vm, upwp, vpwp,
     .              thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .              rtp2, thlp2, rtpthlp,
     .              Scm, taum, Ncm, rrm, Nrm,
     .              isValid, implemented,
     .              sclrm, sclrm_forcing, ! optional
     .              edsclrm, edsclrmt )   ! optional

!      Modules to be included
       use constants
       use model_flags
       use grid_class
       use diagnostic_variables
!      use prognostic_variables
       use mixing
       use diagnose_variances
       use surface
       use pdf_closure
       use rain_equations
#ifdef STATS
       use statistics
#endif /*STATS*/

       implicit none

!      External
       real, external :: rsat
       intrinsic sqrt, min, max, size, exp

!      Input
       integer, intent(in) :: iter      ! closure iteration number
       integer, intent(in) :: niterlong ! total # of closure iterations

       real, intent(in) :: time_initial, time_current ! initial time, current time
       ! Note on dt, dmain, and dtclosure: since being moved out of
       ! hoc.F, all subroutines within parameterization_timestep now use
       ! dt for time dependent calculations.  The old dt is noted in
       ! each section of the code -dschanen 20 April 2006 
       real, intent(in) :: dt ! current timestep size

       real, intent(in), dimension(gr%nnzp) :: thlm_forcing, rtm_forcing
       real, intent(in), dimension(gr%nnzp) :: wmm, wmt
       real, intent(in), dimension(gr%nnzp) :: p, rhom, rhot, exner

!      Input/Output
!      These are prognostic or are planned to be in the future
       real, intent(inout), dimension(gr%nnzp) :: um, vm, upwp, vpwp
       real, intent(inout), dimension(gr%nnzp) :: thlm, rtm
       real, intent(inout), dimension(gr%nnzp) :: wprtp, wpthlp
       real, intent(inout), dimension(gr%nnzp) :: wp2, wp3
       real, intent(inout), dimension(gr%nnzp) :: rtp2, thlp2, rtpthlp
       real, intent(inout), dimension(gr%nnzp) :: Scm, taum
       real, intent(inout), dimension(gr%nnzp) :: Ncm

!      You can safely leave these uninitialized for a non-drizzle case
       real, intent(inout), dimension(gr%nnzp) :: rrm, Nrm

       real, intent(inout), dimension(gr%nnzp,sclrm_dimension),
     . optional :: sclrm, sclrm_forcing, edsclrm

!      Output

       real, intent(out), dimension(gr%nnzp,sclrm_dimension),
     . optional :: edsclrmt

       real:: wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc

       logical, intent(in) :: implemented

       logical, intent(inout) :: isValid

!      Internal

       integer i, j, k, km1, kp1
       real, dimension(1:gr%nnzp) :: tmp1
       real, dimension(1:gr%nnzp) :: Skw_arr_m
       real, dimension(1:gr%nnzp) :: gamma_Skw_fnc

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

       integer, dimension(1:gr%nnzp, 1:nt_repeat, 1:(d_variables+1) ) 
     . :: p_height_time ! matrix of rand ints

! coeffs of s from pdf_closure_new
       real :: crt1, crt2, cthl1, cthl2   

!-------- End Latin hypercube section ----------------------------------

#ifdef SCALARS
       real, dimension(1:gr%nnzp,1:sclrm_dimension) ::
     . sclr_tmp1, sclr_tmp2, sclr_tmp3, sclr_tmp4 ! for pdf closure
#endif /* SCALARS */

!      SET SURFACE VALUES OF FLUXES (BROUGHT IN)
       wpthlp(1) = wpthlp_sfc
       wprtp(1)  = wprtp_sfc
       upwp(1)   = upwp_sfc
       vpwp(1)   = vpwp_sfc


!      Set Surface variances
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
       IF ( bottom_at_sfc ) THEN
          call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .                  wp2(1), thlp2(1), rtp2(1), rtpthlp(1)
#ifdef SCALARS
     .                  ,wpsclrp(1,:), sclrp2(1,:),
     .                  sclrprtp(1,:), sclrpthlp(1,:)
#endif /* SCALARS */
     .                )
       ELSE
          ! Variances for cases where the lowest level is not at the surface.
          ! Eliminate surface effects on lowest level variances.
          wp2(1)     = (2.0/3.0) * emin
          thlp2(1)   = 0.0
          rtp2(1)    = 0.0
          rtpthlp(1) = 0.0
#ifdef SCALARS
          DO i = 1, sclrm_dimension, 1
             sclrp2(1,i)    = 0.0
             sclrprtp(1,i)  = 0.0
             sclrpthlp(1,i) = 0.0
          ENDDO
#endif /* SCALARS */
       ENDIF

!      Diagnose variances
!      We found that if we call diag_var first, we can use a longer timestep.
!      This requires us to store in memory Scm and taum between timesteps.
       call diag_var( taum, wmm, rtm, wprtp, 
     .                thlm, wpthlp, wp2, wp3, Scm, 
     .                rtp2, thlp2, rtpthlp,
     .                .false., dt, isValid ! dtclosure to dt -dschanen
#ifdef SCALARS
     .                ,sclrm, wpsclrp, sclrp2,    ! optional vars
     .                sclrprtp, sclrpthlp         !   "       "
#endif /* SCALARS */
     .               )
       if ( .not. isValid ) return

!      Set Boundary Conditions on thlp2, rtp2, and rtpthlp output.
       thlp2(gr%nnzp)   = 0.0
       rtp2(gr%nnzp)    = 0.0
       rtpthlp(gr%nnzp) = 0.0
#ifdef SCALARS
!      Set Boundary Conditions for optional passive scalars
       do i = 1, sclrm_dimension
          sclrprtp(gr%nnzp,i)  = 0.0
          sclrp2(gr%nnzp,i)    = 0.0
          sclrpthlp(gr%nnzp,i) = 0.0
       enddo ! i = 1, sclrm_dimension
#endif /* SCALARS */

#ifndef DISABLE_GAMMA_SKW
!      Compute gamma as a function of Skw  - 14 April 06 dschanen
       do k = 1, gr%nnzp, 1
        km1 = MAX(k-1,1)
        kp1 = MIN(k+1,gr%nnzp)

        Skw_arr_m(k) = (1.0/2.0) * ( wp3(kp1)+wp3(k) ) /
     .                 MAX( wp2(k), wtol**2.0 )**(3.0/2.0)
       enddo

        gamma_Skw_fnc 
     .  = gamma_coefb + (gamma_coef-gamma_coefb)
     .    *exp( -(1.0/2.0) * (Skw_arr_m/gamma_coefc)**2.0 )

#else
        gamma_Skw_fnc = gamma_coef

#endif /* DISABLE */

       ! Brian added this line of code which can be uncommented
       ! in order to bring back the use of gamma_coef by overwriting
       ! the value of C2_Skw_fnc before the computing the equation
       ! below.
!       gamma_Skw_fnc = 0.3200   ! gamma_coef

!      Compute Sc with new formula from Vince

       Scm = gamma_Skw_fnc * ( 1.0 - 
     .       min(
     .       max( ( wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) 
     .              + 0.01 * wtol * thltol ) )**2,
     .            ( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) 
     .              + 0.01 * wtol * rttol ) )**2 ),
     .       1.0      ) )

       Sct = zm2zt( Scm )


!    Latin hypercube sample generation
!    Generate p_height_time, an nnzp x nt_repeat x d_variables array of random integers
       if ( LH_on ) then
         i_rmd = mod( iter-1, sequence_length )
         if ( i_rmd == 0) then
           call permute_height_time( gr%nnzp, nt_repeat, 
     .                               d_variables+1, p_height_time )
         endif
       endif
!    End Latin hypercube sample generation

!       print*, 'hoc.F: i_rmd=', i_rmd


!      Call closure scheme

#ifdef SCALARS
!      Put closure scheme input on the m grid
       do i = 1, sclrm_dimension
          sclr_tmp1(:,i) = zm2zt( wpsclrp(:,i) ) 
          sclr_tmp2(:,i) = zm2zt( sclrprtp(:,i) ) 
          sclr_tmp3(:,i) = zm2zt( sclrp2(:,i) ) 
          sclr_tmp4(:,i) = zm2zt( sclrpthlp(:,i) ) 
       enddo ! i = 1, sclrm_dimension
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
     .     pdf_parms(k, :), crt1, crt2, cthl1, cthl2
#ifdef SCALARS
     .     ,sclrm(k,:), sclr_tmp1(k,:),
     .     wpsclrprtp(k,:), wpsclrp2(k,:), sclr_tmp2(k,:), 
     .     sclr_tmp3(k,:), sclr_tmp4(k,:), sclrpthvp(k,:), 
     .     wpsclrpthlp(k,:), sclrprcp(k,:), wp2sclrp(k,:)
#endif /*SCALARS*/
     .   )
!        Latin hypercube sampling.
         if ( LH_on ) then 
           call latin_hypercube_sampling
     .     ( k, n_micro_call, d_variables, 
     .       nt_repeat, i_rmd, crt1, crt2, cthl1, cthl2, 
     .       rrm, gr%nnzp, sample_flag, p_height_time )
         endif

       enddo ! k = 2, nz-1

!            print*, 'hoc.F: AKm=', AKm
!            print*, 'hoc.F: AKm_est=', AKm_est

!      Interpolate momentum variables back to momentum grid.
!      Since top momentum level is higher than top thermo level,
!      set variables at top momentum level to 0.

       wp4               = zt2zm( wp4 )
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
       thlprcp           = 0.0
#ifdef SCALARS
!      Interpolate passive scalars back onto the m grid
       do i = 1, sclrm_dimension
          sclrpthvp(:,i)       = zt2zm( sclrpthvp(:,i) )
          sclrpthvp(gr%nnzp,i) = 0.0
          sclrprcp(:,i)        = zt2zm( sclrprcp(:,i) )
          sclrprcp(gr%nnzp,i)  = 0.0
       enddo ! i=1, sclrm_dimension
#endif /* SCALARS */

!      Alternate approach for lower wp4 boundary condition:
!      call pdf scheme at the 1st momentum level, zm(1)

#ifdef SCALARS
       do i = 1, sclrm_dimension
          sclr_tmp1(:,i) = zt2zm( sclrm(:,i) )
       enddo
#endif /* SCALARS */

        call pdf_closure_new
     .      ( p(1), exner(1),
     .        wmm(1), wp2(1), zt2zm( wp3, 1 ), Scm(1),
     .        zt2zm( rtm, 1 ), rtp2(1), wprtp(1),
     .        zt2zm( thlm, 1 ), thlp2(1), wpthlp(1),
     .        rtpthlp(1), wp4(1), wprtp2(1), 
     .        wp2rtp(1), wpthlp2(1), wp2thlp(1),
     .        wprtpthlp(1), cf(1), rcm(1),
     .        wpthvp(1), wp2thvp(1), rtpthvp(1), thlpthvp(1),
     .        wprcp(1), wp2rcp(1), rtprcp(1), thlprcp(1),
     .        pdf_parms(1,:), crt1, crt2, cthl1, cthl2
#ifdef SCALARS
     .        ,sclr_tmp1(1,:), wpsclrp(1,:),
     .        wpsclrprtp(1,:), wpsclrp2(1,:),
     .        sclrprtp(1,:), sclrp2(1,:),
     .        sclrpthlp(1,:), sclrpthvp(1,:),
     .        wpsclrpthlp(1,:), sclrprcp(1,:),
     .        wp2sclrp(1,:)
#endif /* SCALARS */
     .      )

!      Set other specific boundary conditions.
       wpthvp(1)         = wpthlp( 1 )
       thlpthvp(1)       = thlp2( 1 )
       wp2thvp(gr%nnzp)  = 0.0
#ifdef SCALARS
       sclrpthvp(1,SCLR_THETA) = sclrpthlp(1,SCLR_THETA)
#endif /* SCALARS */

!      Compute thvm

       thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

!      tke is assumed to be 3/2 of wp2

       em = 1.5 * wp2

!      Compute mixing length

       call length( thvm, thlm, rtm, rcm, 
     .              em, p, exner, Lscale, lup, ldown )

!      Dissipation time

       tmp1 = SQRT( MAX( emin, zm2zt( em ) ) )
       taut = MIN( Lscale / tmp1, taumax )
       taum = MIN( ( zt2zm( Lscale ) / SQRT( max(em, emin) ) ), taumax )

!      Modification to damp noise in stable region
       do k = 1, gr%nnzp
         if ( wp2(k) <= 0.005 ) then
           taut(k) = taumin
           taum(k) = taumin
         endif
       enddo

!      Eddy diffusivity coefficient
!      c_K is 0.548 usually (Duynkerke and Driedonks 1987)

       Kht = c_K * Lscale * tmp1
       Khm = c_K * zt2zm( Lscale ) * SQRT( max( em, emin ) )


!      call Khairoutdinov and Kogan (2000) scheme for rain microphysics.

       IF ( kk_rain ) THEN

          CALL rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .               rtm, rcm, Ncm, rrm, Nrm, cond_rrm,
     .               auto_rrm, accr_rrm, Vrr, mean_vol_rad,
     .               AKm_est, AKm, isValid )

          ! Rainfall rate (mm/day) should be defined on thermodynamic
          ! levels.  Brian.
          rain_rate = ( rrm * zm2zt(Vrr) ) * ( rhot / rho_lw )   ! rainfall rate (m/s)
     .                                     * ( 86400.0 * 1000.0 )

          ! Precipitation Flux (W/m^2) should be defined on
          ! momentum levels.  Brian.
          Fprec = ( zt2zm(rrm) * Vrr ) * ( rhom / rho_lw )
     .                                 * rho_lw * Lv

       ENDIF ! kk_rain

       if ( .not. isValid ) return

!      Store the saturation mixing ratio for output purposes.  Brian
       do k = 1, gr%nnzp, 1
         rsm(k) =
     .   rsat(p(k), (thlm(k)+(Lv/(Cp*exner(k)))*rcm(k))*exner(k))
       enddo

!#######################################################################
!############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
!#######################################################################


!      Advance rtm/wprtp and thlm/wpthlp one time step

       call timestep_mixing
     .      ( dt, taum, wmm, wmt, wp2, wp3, Scm, 
     .        rtm, wprtp, rtpthvp, rtm_forcing,
     .        thlm, wpthlp, thlpthvp, thlm_forcing, 
     .        rhot, rhom, exner, rcm, Ncm, 
     .        cond_rrm, auto_rrm, accr_rrm,             
     .        Fcsed, sed_rcm,
     .        cloud_sed, kk_rain, isValid, implemented  
#ifdef SCALARS
     .        ,sclrm, wpsclrp, sclrpthvp, sclrm_forcing   ! optional
#endif /* SCALARS */
     .      )
       if ( .not. isValid ) return

!      Set Boundary Conditions on wpthlp and wprtp output.
       wpthlp(gr%nnzp)  = 0.0
       wprtp(gr%nnzp)   = 0.0
#ifdef SCALARS
!      Set Boundary Conditions for optional passive scalars
       do i = 1, sclrm_dimension
          wpsclrp(gr%nnzp,i)   = 0.0
       enddo ! i = 1, sclrm_dimension
#endif /* SCALARS */


!      Advance wp2/wp3 one timestep

       call timestep_wp23
     .      ( dt, Lscale, taum, taut, Scm, ! dtclosure to dt -dschanen
     .        pdf_parms(:, 13), Kht, wmt, wmm,
     .        um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .        wp2, wp3, isValid )

       if ( .not. isValid ) return

!      Set Boundary Conditions on wp2 and wp3 output.
       wp2(gr%nnzp)     = 0.0
       wp3(gr%nnzp)     = 0.0


#ifdef SCALARS
!      Eddy Diffusivity Code
       do i=1, sclrm_dimension
         edsclrmt(1:gr%nnzp,i) = - wmt * ddzm( zt2zm( edsclrm(:,i) ) )
         if ( i == SCLR_THETA) then
          ! Note changed from dtmain=>dt -dschanen
           call update_ed_var( edsclrm(:,i), edsclrmt(:,i), wpthlp(1), 
     .                         Khm, dt, wpedsclrp(:,i), isValid )
         else ! i == SCLR_RT
           call update_ed_var( edsclrm(:,i), edsclrmt(:,i), wprtp(1),
     .                         Khm, dt, wpedsclrp(:,i), isValid )
         end if
       end do
       ! Set boundary condition as in rt
       edsclrm(1,1:sclrm_dimension) = edsclrm(2,1:sclrm_dimension)
       if ( .not. isValid ) return
#endif /*SCALARS*/


!      Update winds

       umt = + fcor * ( vm - vg )
       IF ( .NOT. implemented ) THEN
          umt = umt - ( wmt * ddzm( zt2zm( um ) ) )
       ENDIF
       ! dtmain to dt -dschanen
       call update_ed_var( um, umt, upwp(1), Khm, dt, upwp, isValid)
       um(1)       = ( ( um(3)-um(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + um(2)
       um(gr%nnzp) = ( ( um(gr%nnzp-1)-um(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + um(gr%nnzp-1)

       if ( .not. isValid ) return

       vmt = - fcor * ( um - ug )
       IF ( .NOT. implemented ) THEN
          vmt = vmt - ( wmt * ddzm( zt2zm( vm ) ) )
       ENDIF
       ! dtmain to dt -dschanen
       call update_ed_var( vm, vmt, vpwp(1), Khm, dt, vpwp, isValid)
       vm(1)       = ( ( vm(3)-vm(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + vm(2)
       vm(gr%nnzp) = ( ( vm(gr%nnzp-1)-vm(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + vm(gr%nnzp-1)
           
       if ( .not. isValid ) return


!      Shear production.  Brian
       do k = 1, gr%nnzp-1, 1
         shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k)
     .              -vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
       enddo
       shear(gr%nnzp) = 0.0

!#######################################################################
!#############            ACCUMULATE STATISTICS            #############
!#######################################################################

#ifdef STATS
!      Store values of surface fluxes for statistics
       if (lstats_samp) then
         if ( ish > 0 ) then
           sfc%x(1,ish) = sfc%x(1,ish) + wpthlp(1) * rhom(1) * cp
           sfc%n(1,ish) = sfc%n(1,ish) + 1
         endif
         if ( ilh > 0 ) then
           sfc%x(1,ilh) = sfc%x(1,ilh) + wprtp(1) * rhom(1) * Lv
           sfc%n(1,ilh) = sfc%n(1,ilh) + 1
         endif
         if ( iustar > 0 ) then   ! Added by Brian
           sfc%x(1,iustar) = sfc%x(1,iustar) + ustar
           sfc%n(1,iustar) = sfc%n(1,iustar) + 1
         endif
         if ( irain > 0 .and. kk_rain ) then   ! Added by Brian
           sfc%x(1,irain) = sfc%x(1,irain) + rain_rate(2)
           sfc%n(1,irain) = sfc%n(1,irain) + 1
         endif
         if ( ipflux > 0 .and. kk_rain ) then   ! Added by Brian
           sfc%x(1,ipflux) = sfc%x(1,ipflux) + Fprec(1)
           sfc%n(1,ipflux) = sfc%n(1,ipflux) + 1
         endif
       endif 

!      Added to allow tuning without using the zm stats file

       wp2zt     = zm2zt( wp2 )
       thlp2zt   = zm2zt( thlp2 )
       wpthlpzt  = zm2zt( wpthlp )
       wprtpzt   = zm2zt( wprtp )
       rtp2zt    = zm2zt( rtp2 )
       rtpthlpzt = zm2zt( rtpthlp )

!      call stats_accumulate( time_current + dtmain, dtmain )
       call stats_accumulate( )

#endif /*STATS*/

       return
       end subroutine parameterization_timestep
!-----------------------------------------------------------------------
!       Subroutine latin_hypercube_sampling()

!       Estimate using Latin Hypercubes.  This is usually disabled by default.
!       The actual generation of a random matrix is done in a call from the
!       subroutine hoc_initialize to permute_height_time()

!       Calls:
!       subroutine lh_sampler
!       subroutine micro_calcs
!-----------------------------------------------------------------------
        subroutine latin_hypercube_sampling( k, n, dvar, i_rmd, nt,
     .                                       crt1, crt2, cthl1, cthl2,
     .                                       rrm, grid, 
     .                                       sflag, p_height_time )
        !use model_arrays
        use diagnostic_variables, only: cf, pdf_parms, AKm_est, 
     .  AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, rcm_est

        implicit none  

!       Input 
        integer, intent(in) :: k  ! index
        integer, intent(in) :: n, dvar, i_rmd, nt, grid
        logical, intent(in) :: sflag

        ! coeffs of s from pdf_closure_new
        real, intent(in) :: crt1, crt2, cthl1, cthl2

        real, dimension(grid), intent(in) :: rrm

        integer, dimension(1:grid, 1:nt, 1:(dvar+1) ) ::
     .  p_height_time ! matrix of rand ints

!       Internal

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
     .             ( nzmax, deltaz, C1, C1b, C1c, C2, C2b, C2c, C4, C5, 
     .               C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .               C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12,
     .               nu1, nu2, nu6, nu8, nu_r, 
     .               gamma_coef, gamma_coefb, gamma_coefc,
     .               mu, c_K, beta, lmin_coef, taumin, taumax, fcor,
     .               lbugsrad, kk_rain, cloud_sed, bottom_at_sfc,
     .               implemented, host_momentum, host_thermodynamic )

        use grid_class
        use diagnostic_variables
        use prognostic_variables

        implicit none

!       External
        external setup_constants, setup_model_flags

!       Input

        ! Grid definition
        integer, intent(in) :: nzmax     ! extent in levels
        real, intent(in)    :: deltaz    ! change in altitude in m per level

        ! Grid setup if part of larger (host) model.
        real, dimension(1:nzmax):: host_momentum
        real, dimension(1:nzmax):: host_thermodynamic

        ! Constant parameters
        real, intent(in) :: C1, C1b, C1c
        real, intent(in) :: C2, C2b, C2c
        real, intent(in) :: C4, C5
        real, intent(in) :: C6rt, C6rtb, C6rtc
        real, intent(in) :: C6thl, C6thlb, C6thlc
        real, intent(in) :: C7, C7b, C7c
        real, intent(in) :: C8, C8b, C10
        real, intent(in) :: C11, C11b, C11c
        real, intent(in) :: C12
        real, intent(in) :: nu1, nu2, nu6, nu8, nu_r
        real, intent(in) :: beta
        real, intent(in) :: lmin_coef, taumin, taumax
        real, intent(in) :: gamma_coef, gamma_coefb, gamma_coefc
        real, intent(in) :: mu, c_K
        real, intent(in) :: fcor

        ! Flags
        logical, intent(in) :: lbugsrad  ! BUGSrad interactive radiation scheme
        logical, intent(in) :: kk_rain   ! K & K rain microphysics
        logical, intent(in) :: cloud_sed ! Cloud Sedimentation
        ! Whether the lowest model level is at the surface or not (T/F).
        logical, intent(in) :: bottom_at_sfc


        ! Whether this is part of larger model (T/F)
        logical, intent(in) :: implemented 

!       Setup flags

        call setup_model_flags( lbugsrad, kk_rain, cloud_sed,
     .                          bottom_at_sfc )

        if ( .not. implemented ) then
          call setup_diagnostic_variables( nzmax )
        end if

        call setup_prognostic_variables( nzmax )

!       Setup grid

        call gridsetup( nzmax, deltaz, implemented,
     .                  host_momentum, host_thermodynamic )

!       Define constants

        call setup_constants
     .       ( deltaz, C1, C1b, C1c, C2, C2b, C2c, C4, C5, 
     .         C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .         C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12,
     .         nu1, nu2, nu6, nu8, nu_r, 
     .         gamma_coef, gamma_coefb, gamma_coefc,
     .         mu, c_K, beta, lmin_coef, taumin, taumax, fcor )

        return
        end subroutine parameterization_setup

!-----------------------------------------------------------------------
        subroutine parameterization_cleanup( implemented )
        use diagnostic_variables
        use prognostic_variables
#ifdef STATS
        use statistics
#endif

        implicit none

        ! whether this is part of larger model (T/F)
        logical, intent(in) :: implemented 
  
        if ( .not. implemented ) then
          call cleanup_diagnostic_variables( )
        end if

        call cleanup_prognostic_variables( )

        return
        end subroutine parameterization_cleanup

        end module hoc_parameterization_interface
