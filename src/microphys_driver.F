!-----------------------------------------------------------------------
! $Id: microphys_driver.F,v 1.7 2008-06-06 17:57:50 griffinb Exp $
        module microphys_driver

!       Description:
!       Call a microphysical scheme to compute hydrometeor species,
!       and advect, sediment, & diffuse using a tridiagonal system.

!       References:
!       None
!-----------------------------------------------------------------------
        implicit none

        public :: advance_microphys

        ! Subroutines
        private :: microphys_lhs, microphys_solve
        private :: adj_microphys_tndcy

        ! Functions
        private :: sedimentation

        private ! Default

        contains

!-----------------------------------------------------------------------
        subroutine advance_microphys
     .             ( runtype, dt, time_current, time_initial, 
     .               thlm, p, exner, rhot, rhom, rtm, rcm, Ncm, 
     .               pdf_parms, wmt, wmm, Khm, AKm_est, AKm, 
     .               Ncnm, Nim,
     .               hydromet,
     .               rtm_forcing, thlm_forcing, err_code )

!       Description:
!       Compute pristine ice, snow, graupel, & rain hydrometeor fields.
!       Uses implicit discretization.

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class, only:
     .      gr, ! Variable(s)
     .      zm2zt, ! Procedure(s)
     .      zt2zm

        use rain_equations, only:
     .      kk_microphys ! Procedure(s)

        use ice_dfsn_mod, only:
     .      ice_dfsn ! Procedure(s)

        use model_flags, only:
     .      lcoamps_micro, ! Procedure(s)
     .      kk_rain,
     .      licedfs

        use parameters, only:
     .      c_Krrm, ! Procedure(s) 
     .      nu_r

        use constants, only: 
     .      hydromet_dim, ! Variable(s)
     .      Lv, 
     .      Cp, 
     .      rho_lw, 
     .      fstderr

        use stats_prec, only: 
     .      time_prec ! Variable(s)

        use error_code, only: 
     .      lapack_error ! Procedure

        use coamps_micro_driver_mod, only: 
     .      coamps_micro_driver ! Procedure

        use temp_in_K_mod, only: thlm2temp_in_K


#ifdef STATS
        use stats_hoc, only:
     .      iVrr, ! Variable(s)
     .      iVnr,
     .      iVsnow,
     .      iVice,
     .      iVgraupel,
     .      ithlm_mc,
     .      irtm_mc,
     .      irain_rate,
     .      iFprec,
     .      zm,
     .      irrm_bt,
     .      irrm_mc,
     .      irrm_cond_adj,
     .      irrm_cl,
     .      iNrm_bt,
     .      iNrm_mc,
     .      iNrm_cond_adj,
     .      iNrm_cl,
     .      iNcm,
     .      iNim,
     .      iNcnm,
     .      irrm,
     .      iNrm,
     .      irsnowm,
     .      iricem,
     .      irgraupelm,
     .      iricem_bt,
     .      iricem_mc,
     .      iricem_cl,
     .      irgraupelm_bt,
     .      irgraupelm_mc,
     .      irgraupelm_cl,
     .      irsnowm_bt,
     .      irsnowm_mc,
     .      irsnowm_cl,
     .      zt,
     .      irain,
     .      ipflux,
     .      irrm_sfc,
     .      sfc,
     .      lstats_samp

#endif
        implicit none

        ! Constant Parameters


        character(len=*), intent(in) ::
     .  runtype ! Name of the run, for case specific effects.

        real(kind=time_prec), intent(in) :: 
     .  dt           ! Timestep         [s]

        real(kind=time_prec), intent(in) :: 
     .  time_current,! Current time     [s]
     .  time_initial ! Initial time     [s]

        real, dimension(gr%nnzp), intent(in) ::
     .  thlm,   ! Liquid potential temp.                 [K]
     .  p,      ! Pressure                               [Pa]
     .  exner,  ! Exner function                         [-]
     .  rhot,   ! Density on thermo. grid                [kg/m^3]
     .  rhom,   ! Density on moment. grid                [kg/m^3]
     .  rtm,    ! Total water mixing ratio               [kg/kg]
     .  rcm,    ! Liquid water mixing ratio              [kg/kg]
     .  wmt,    ! w wind on moment. grid                 [m/s]
     .  wmm,    ! w wind on thermo. grid                 [m/s]
     .  Khm,    ! Kh Eddy diffusivity on momentum grid   [m^2/s]
     .  Akm_est,! Analytic Kessler ac                    [kg/kg]
     .  Akm     ! Analytic Kessler estimate              [kg/kg]

        ! Note:
        ! K & K only uses Ncm, while for COAMPS Ncnm is initialized
        ! and Nim & Ncm are computed within subroutine adjtg.
        real, dimension(gr%nnzp), intent(inout) ::
     .  Ncm,    ! Cloud drop number concentration       [count/kg]
     .  Ncnm,   ! Cloud nuclei number concentration     [count/m^3]
     .  Nim     ! Ice crystal number concentration      [count/m^3]

        real, target, dimension(gr%nnzp,26), intent(in) ::
     .  pdf_parms     ! PDF parameters

        real, target, dimension(gr%nnzp,hydromet_dim), intent(inout) ::
     .  hydromet      ! Array of rain, prist. ice, graupel, etc. [units vary]

        real, dimension(gr%nnzp), intent(inout) ::
     .  rtm_forcing, ! Imposed contributions to total water        [kg/kg/s]
     .  thlm_forcing ! Imposed contributions to liquid potential temp. [K/s]

        integer, intent(out) :: err_code

        ! Local Variables
        real, dimension(3,gr%nnzp) ::
     .  lhs ! Left hand side of tridiagonal matrix

        real, dimension(gr%nnzp) ::
     .  Vrr,     ! Rain mixing ratio sedimentation velocity    [m/s]
     .  VNr,     ! Rain number conc. sedimentation velocity    [m/s]
     .  Vice,    ! Ice mixing ratio sedimentation velocity     [m/s]
     .  Vsnow,   ! Snow mixing ratio sedimentation velocity    [m/s]
     .  Vgraupel ! Graupel mixing ratio sedimentation velocity [m/s]

        real, dimension(gr%nnzp) ::
     .  rtm_mc,  ! Change in total water due to microphysics    [(kg/kg)/s]
     .  thlm_mc  ! Change in liquid potential temperature 
                 ! due to microphysics                          [K/s]

        real, dimension(gr%nnzp) ::
     .  rrm_tndcy,      ! Rain mixing ratio tendency      [(kg/kg)/s]
     .  Nrm_tndcy,      ! Rain number conc. tendency      [(count/kg)/s]
     .  rsnowm_tndcy,   ! Snow mixing ratio tendency      [(kg/kg)/s]
     .  ricem_tndcy,    ! Ice mixing ratio tendency       [(kg/kg)/s]
     .  rgraupelm_tndcy ! Graupel mixing ratio tendency   [(kg/kg)/s]  

        real, dimension(gr%nnzp) ::
     .  T_in_K  ! Temperature   [K]

        real, dimension(1,1,gr%nnzp) ::
     .  cond ! COAMPS stat for condesation/evap of rcm

        real, pointer, dimension(:) ::
     .  rrm,      ! Pointer for rain water mixing ratio   [kg/kg]
     .  Nrm,      ! Pointer for rain droplet number conc. [count/kg]
     .  rsnowm,   ! Pointer for snow mixing ratio         [kg/kg]
     .  ricem,    ! Pointer for ice mixing ratio          [kg/kg]
     .  rgraupelm ! Pointer for graupel mixing ratio      [kg/kg]

        ! Various PDF parameters needed for Brian's K&K microphysics
        real, pointer, dimension(:) ::
     .  a,
     .  thl1, thl2,
     .  s1, s2,
     .  ss1, ss2,
     .  rc1, rc2

        ! Eddy diffusivity for rain and rain drop concentration.
        ! It is also used for the other hydrometeor variables.
        ! Kr = Constant * Khm; Constant is named c_Krrm.
        real, dimension(gr%nnzp) :: Kr   ! [m^2/s]

        ! Variable needed to handle correction to rtm and thlm microphysics
        ! tendency arrays, as well rrm_cond and Nrm_cond statistical
        ! tendency arrays, due to a negative result being produced by 
        ! over-evaporation of rain water over the course of a timestep.
        ! Brian Griffin.  April 14, 2007.
        real :: overevap_rate ! Absolute value of negative evap. rate.

        integer :: k, km1, kp1 ! Array indices

!-----------------------------------------------------------------------

        ! Assign pointers to each field
        rrm       => hydromet(1:gr%nnzp,1)
        Nrm       => hydromet(1:gr%nnzp,2)
        rsnowm    => hydromet(1:gr%nnzp,3)
        ricem     => hydromet(1:gr%nnzp,4)
        rgraupelm => hydromet(1:gr%nnzp,5)

        ! Assign pointers to pdf_parms
        thl1 => pdf_parms(1:gr%nnzp,9)
        thl2 => pdf_parms(1:gr%nnzp,10)
        a    => pdf_parms(1:gr%nnzp,13)
        rc1  => pdf_parms(1:gr%nnzp,14)
        rc2  => pdf_parms(1:gr%nnzp,15)
        s1   => pdf_parms(1:gr%nnzp,20)
        s2   => pdf_parms(1:gr%nnzp,21)
        ss1  => pdf_parms(1:gr%nnzp,22)
        ss2  => pdf_parms(1:gr%nnzp,23)

        ! Solve for the value of Kr, the hydrometeor eddy diffusivity.
        do k = 1, gr%nnzp, 1
           Kr(k) = c_Krrm * Khm(k)
        end do

        ! Determine temperature in K for the microphysics
        T_in_K = thlm2temp_in_K( thlm, exner, rcm )

        ! Begin by calling either Brian Griffin's implementation of the
        ! Khairoutdinov and Kogan microphysical scheme or 
        ! alternatively the Rutlege and Hobbes scheme from COAMPS(R).
        ! Note: COAMPS appears to have some K&K elements to it as well.

        if ( lcoamps_micro ) then

           call coamps_micro_driver
     .          ( runtype, gr%nnzp-1, gr%nnzp, time_current, dt,
     .            rtm, wmm, p, exner, rhot, T_in_K,
     .            thlm, ricem, rrm, rgraupelm, rsnowm,
     .            rcm, Ncm, Nrm, Ncnm, Nim,
     .            cond, Vsnow, Vice, Vrr, VNr, Vgraupel,
     .            ricem_tndcy, rrm_tndcy, rgraupelm_tndcy, 
     .            rsnowm_tndcy, Nrm_tndcy,
     .            rtm_mc, thlm_mc )

#ifdef STATS
          if ( lstats_samp ) then
              
           ! Sedimentation velocity for rrm
           if ( iVrr > 0 ) then
             zm%x(:,iVrr) = zm%x(:,iVrr) + Vrr(:)
             zm%n(:,iVrr) = zm%n(:,iVrr) + 1
           end if

           ! Sedimentation velocity for Nrm
           if ( iVNr > 0 ) then
             zm%x(:,iVNr) = zm%x(:,iVNr) + VNr(:)
             zm%n(:,iVNr) = zm%n(:,iVNr) + 1
           end if

           ! Sedimentation velocity for snow
           if ( iVsnow > 0 ) then
             zm%x(:,iVsnow) = zm%x(:,iVsnow) + Vsnow(:)
             zm%n(:,iVsnow) = zm%n(:,iVsnow) + 1
           end if

           ! Sedimentation velocity for pristine ice
           if ( iVice > 0 ) then
             zm%x(:,iVice) = zm%x(:,iVice) + Vice(:)
             zm%n(:,iVice) = zm%n(:,iVice) + 1
           end if

           ! Sedimentation velocity for graupel
           if ( iVgraupel > 0 ) then
             zm%x(:,iVgraupel) = zm%x(:,iVgraupel) + Vgraupel(:)
             zm%n(:,iVgraupel) = zm%n(:,iVgraupel) + 1
           end if

           ! Sum total of rrm microphysics
           if ( irrm_mc > 0 ) then
             zt%x(:,irrm_mc) = zt%x(:,irrm_mc) + rrm_tndcy(:)
             zt%n(:,irrm_mc) = zt%n(:,irrm_mc) + 1
           end if

           ! Sum total of Nrm microphysics
           if ( iNrm_mc > 0 ) then
             zt%x(:,iNrm_mc) = zt%x(:,iNrm_mc) + Nrm_tndcy(:)
             zt%n(:,iNrm_mc) = zt%n(:,iNrm_mc) + 1
           end if

           ! Sum total of pristine ice microphysics
           if ( iricem_mc > 0 ) then
             zt%x(:,iricem_mc) = zt%x(:,iricem_mc) + ricem_tndcy(:)
             zt%n(:,iricem_mc) = zt%n(:,iricem_mc) + 1
           end if

           ! Sum total of graupel microphysics
           if ( irgraupelm_mc > 0 ) then
             zt%x(:,irgraupelm_mc) = zt%x(:,irgraupelm_mc) 
     .                             + rgraupelm_tndcy(:)
             zt%n(:,irgraupelm_mc) = zt%n(:,irgraupelm_mc) + 1
           end if

           ! Sum total of snow microphysical processeses
           if ( irsnowm_mc > 0 ) then
             zt%x(:,irsnowm_mc) = zt%x(:,irsnowm_mc) + rsnowm_tndcy(:)
             zt%n(:,irsnowm_mc) = zt%n(:,irsnowm_mc) + 1
           end if

         end if ! lstats_samp
#endif /*STATS*/
        else if ( kk_rain ) then

          ! Note that Ncm for DYCOMS II RF02 is specified in the GCSS,
          ! so new cases will need this computed beforehand as well.

          call kk_microphys
     .         ( thlm, T_in_K, p, exner, rhot, 
     .           thl1, thl2, a, rc1, rc2,
     .           s1, s2, ss1, ss2, rtm, rcm, Ncm, rrm, Nrm,
     .           .true., AKm_est, AKm, 
     .           rrm_tndcy, Nrm_tndcy, 
     .           rtm_mc, thlm_mc,
     .           Vrr, VNr )

#ifdef STATS
          if ( lstats_samp ) then
              
           ! Sedimentation velocity for rrm
           if ( iVrr > 0 ) then
             zm%x(:,iVrr) = zm%x(:,iVrr) + Vrr(:)
             zm%n(:,iVrr) = zm%n(:,iVrr) + 1
           end if

           ! Sedimentation velocity for Nrm
           if ( iVNr > 0 ) then
             zm%x(:,iVNr) = zm%x(:,iVNr) + VNr(:)
             zm%n(:,iVNr) = zm%n(:,iVNr) + 1
           end if

           ! Sum total of rrm microphysics (auto + accr + cond)
           if ( irrm_mc > 0 ) then
             zt%x(:,irrm_mc) = zt%x(:,irrm_mc) + rrm_tndcy(:)
             zt%n(:,irrm_mc) = zt%n(:,irrm_mc) + 1
           end if

           ! Sum total of Nrm microphysics (auto + cond)
           if ( iNrm_mc > 0 ) then
             zt%x(:,iNrm_mc) = zt%x(:,iNrm_mc) + Nrm_tndcy(:)
             zt%n(:,iNrm_mc) = zt%n(:,iNrm_mc) + 1
           end if

         end if ! lstats_samp
#endif /*STATS*/

        end if ! coamps micro or KK rain. 

        ! Compute rrm and Nrm for either case
        if ( lcoamps_micro .or. kk_rain ) then
          ! --------------------------------------------------------------------
          ! Rain water mixing ratio (rrm)
          ! --------------------------------------------------------------------
          call microphys_lhs
     .         ( "rrm", .true., dt, Kr, nu_r, wmt, Vrr, rrm, lhs )

       
          call microphys_solve
     .         ( "rrm", dt, lhs, rrm_tndcy, rrm, err_code )

          ! Handle over-evaporation of rrm and adjust rt and theta-l hydrometeor
          ! tendency arrays accordingly.
          do k = 1, gr%nnzp, 1
            if ( rrm(k) < 0.0 ) then


              call adj_microphys_tndcy( rrm_tndcy, wmt, Vrr, Kr, nu_r, 
     .                                  dt, k, .true., "rrm",
     .                                  rrm, overevap_rate )

              ! overevap_rate is defined as positive.
              ! It is a correction factor.
              rtm_mc(k)  = rtm_mc(k) - overevap_rate

              thlm_mc(k) = thlm_mc(k)
     .                  + ( Lv / ( Cp*exner(k) ) ) * overevap_rate

            ! Moved from adj_microphys_tndcy  
#ifdef STATS
              if( lstats_samp .and. irrm_cond_adj > 0 ) then
                zt%x(k,irrm_cond_adj) = zt%x(k,irrm_cond_adj) 
     .            + overevap_rate
              end if
#endif /*STATS*/

            end if ! rrm(k) < 0

#ifdef STATS
            if ( lstats_samp .and. irrm_cond_adj > 0 ) then
              zt%n(k,irrm_cond_adj) = zt%n(k,irrm_cond_adj) + 1
            end if
             ! Joshua Faschinj December 2007
#endif /*STATS*/

          end do ! k=1..gr%nnzp

#ifdef STATS
          ! Store computed value before clipping
          if (lstats_samp .and. irrm_cl > 0 ) then
            zt%x(:,irrm_cl) = zt%x(:,irrm_cl) - rrm / dt
          end if
#endif

          ! Value of rain water mixing ratio cannot fall below 0, 
          ! so we clip accordingly
          where ( rrm < 0.0 ) rrm = 0.0

#ifdef STATS
          if (lstats_samp .and. irrm_cl > 0 ) then
            zt%x(:,irrm_cl) = zt%x(:,irrm_cl) + rrm / dt
            zt%n(:,irrm_cl) = zt%n(:,irrm_cl) + 1
          end if
#endif

#ifdef STATS
          ! Finalize time tendency
          if ( lstats_samp .and. irrm_bt > 0 ) then
            zt%x(:,irrm_bt) = zt%x(:,irrm_bt) + ( rrm / dt )
            zt%n(:,irrm_bt) = zt%n(:,irrm_bt) + 1
          end if
#endif

          !---------------------------------------------------------------------
          ! Rain drop number concentration (Nrm)
          !---------------------------------------------------------------------

          call microphys_lhs
     .         ( "Nrm", .true., dt, Kr, nu_r, wmt, VNr, Nrm, lhs )

          call microphys_solve
     .         ( "Nrm", dt, lhs, Nrm_tndcy, Nrm, err_code )


          ! Handle over-evaporation similar to rrm.  However, in the case of Nrm
          ! there is no effect on rtm or on thlm.  Brian Griffin.  April 14, 2007.
          do k = 1, gr%nnzp, 1
            if ( Nrm(k) < 0.0 ) then

              call adj_microphys_tndcy( Nrm_tndcy, wmt, VNr, Kr, nu_r,
     .                                  dt, k, .true., "Nrm",
     .                                  Nrm, overevap_rate )

              ! Moved from adj_microphys_tndcy
#ifdef STATS
              if ( lstats_samp .and. iNrm_cond_adj > 0 ) then
                   zt%x(k,iNrm_cond_adj) = zt%x(k,iNrm_cond_adj) 
     .                 + overevap_rate
              end if
#endif /*STATS*/

            end if

#ifdef STATS
            if ( lstats_samp .and. iNrm_cond_adj > 0 ) then
              zt%n(k,iNrm_cond_adj) = zt%n(k,iNrm_cond_adj) + 1
            end if
#endif /*STATS*/
              ! Joshua Fasching December 2007
          end do

#ifdef STATS
          ! Store computed value before clipping
          if (lstats_samp .and. iNrm_cl > 0 ) then
            zt%x(:,iNrm_cl) = zt%x(:,iNrm_cl) - Nrm / dt
          end if
#endif

          ! Value of rain water mixing ratio cannot fall below 0, 
          ! so we clip accordingly
          where ( Nrm < 0.0 ) Nrm = 0.0

#ifdef STATS
          if (lstats_samp .and. iNrm_cl > 0 ) then
            zt%x(:,iNrm_cl) = zt%x(:,iNrm_cl) + Nrm / dt
            zt%n(:,iNrm_cl) = zt%n(:,iNrm_cl) + 1
          end if
#endif

#ifdef STATS
          ! Finalize time tendency
          if ( lstats_samp .and. iNrm_bt > 0 ) then
            zt%x(:,iNrm_bt) = zt%x(:,iNrm_bt) + ( Nrm / dt )
            zt%n(:,iNrm_bt) = zt%n(:,iNrm_bt) + 1
          end if
#endif
        end if ! kk_rain or lcoamps_micro

        ! Compute rsnowm, rgraupelm, and ricem for COAMPS microphysics only.

        if ( lcoamps_micro ) then

          ! --------------------------------------------------------------------
          ! Snow mixing ratio (rsnowm)
          ! --------------------------------------------------------------------

          call microphys_lhs
     .         ( "rsnowm", .true., dt, Kr, nu_r, wmt, Vsnow, 
     .           rsnowm, lhs )

          call microphys_solve
     .         ( "rsnowm", dt, lhs, rsnowm_tndcy, rsnowm, err_code )

          ! Both sublimation and/or melting can contribute to the depletion of
          ! any frozen water particle.  In addition, any type of frozen water
          ! particle can be depleted by turning into another type of frozen
          ! water particle.  For example, pristine ice can turn into graupel
          ! through heavy riming.  It can also turn into snow after combining
          ! with more ice particles so it is heavy enough to fall.
          !
          ! In any place where a negative result is obtained from excess 
          ! depletion over the course of a timestep, it isn't automatically 
          ! clear which process is the most responsible.  Therefore, any
          ! negative result will be clipped and the negative tendency will be
          ! stored in the clipping array for statistical purposes.

#ifdef STATS
          ! Store computed value before clipping
          if (lstats_samp .and. irsnowm_cl > 0 ) then
            zt%x(:,irsnowm_cl) = zt%x(:,irsnowm_cl) - rsnowm / dt
          end if
#endif

          ! Value of snow mixing ratio cannot fall below 0, 
          ! so we clip accordingly
          where ( rsnowm < 0.0 ) rsnowm = 0.0

#ifdef STATS
          if (lstats_samp .and. irsnowm_cl > 0 ) then
            zt%x(:,irsnowm_cl) = zt%x(:,irsnowm_cl) + rsnowm / dt
            zt%n(:,irsnowm_cl) = zt%n(:,irsnowm_cl) + 1
          end if
#endif

#ifdef STATS
          ! Finalize time tendency
          if ( lstats_samp .and. irsnowm_bt > 0 ) then
            zt%x(:,irsnowm_bt) = zt%x(:,irsnowm_bt) + ( rsnowm / dt )
            zt%n(:,irsnowm_bt) = zt%n(:,irsnowm_bt) + 1
          end if
#endif

          ! --------------------------------------------------------------------
          ! Graupel mixing ratio (rgraupelm)
          ! --------------------------------------------------------------------

          call microphys_lhs
     .         ( "rgraupelm", .true., dt, Kr, nu_r, wmt, Vgraupel, 
     .            rgraupelm, lhs )

          call microphys_solve
     .         ( "rgraupelm", dt, lhs, rgraupelm_tndcy,
     .           rgraupelm, err_code )

          ! Handled the same as rsnowm, as listed above.

#ifdef STATS
          ! Store computed value before clipping
          if (lstats_samp .and. irgraupelm_cl > 0 ) then
            zt%x(:,irgraupelm_cl) = zt%x(:,irgraupelm_cl) 
     .                              - rgraupelm / dt
          end if
#endif

          ! Value of graupel mixing ratio cannot fall below 0, 
          ! so we clip accordingly
          where ( rgraupelm < 0.0 ) rgraupelm = 0.0

#ifdef STATS
          if (lstats_samp .and. irgraupelm_cl > 0 ) then
            zt%x(:,irgraupelm_cl) = zt%x(:,irgraupelm_cl) 
     .                              + rgraupelm / dt
            zt%n(:,irgraupelm_cl) = zt%n(:,irgraupelm_cl) + 1
          end if
#endif

#ifdef STATS
          ! Finalize time tendency
          if ( lstats_samp .and. irgraupelm_bt > 0 ) then
            zt%x(:,irgraupelm_bt) = zt%x(:,irgraupelm_bt) 
     .                              + ( rgraupelm / dt )
            zt%n(:,irgraupelm_bt) = zt%n(:,irgraupelm_bt) + 1
          end if
#endif

          ! --------------------------------------------------------------------
          ! Pristine ice mixing ratio (ricem)
          ! --------------------------------------------------------------------

          ! Note: Ice sediments in new COAMPS

          call microphys_lhs
     .         ( "ricem", .true., dt, Kr, nu_r, wmt, Vice, ricem, lhs )

          call microphys_solve
     .         ( "ricem", dt, lhs, ricem_tndcy, ricem, err_code )

          ! Handled the same as rsnowm and rgraupelm, as listed above.

#ifdef STATS
          ! Store computed value before clipping
          if (lstats_samp .and. iricem_cl > 0 ) then
            zt%x(:,iricem_cl) = zt%x(:,iricem_cl) - ricem / dt
          end if
#endif

          ! Value of ice mixing ratio cannot fall below 0, 
          ! so we clip accordingly
          where ( ricem < 0.0 ) ricem = 0.0

#ifdef STATS
          if (lstats_samp .and. iricem_cl > 0 ) then
            zt%x(:,iricem_cl) = zt%x(:,iricem_cl) + ricem / dt
            zt%n(:,iricem_cl) = zt%n(:,iricem_cl) + 1
          end if
#endif

#ifdef STATS
          ! Finalize time tendency
          if ( lstats_samp .and. iricem_bt > 0 ) then
            zt%x(:,iricem_bt) = zt%x(:,iricem_bt) + ( ricem / dt )
            zt%n(:,iricem_bt) = zt%n(:,iricem_bt) + 1
          end if
#endif

        end if ! lcoamps_micro

        ! Call the ice diffusion scheme
        if ( licedfs ) then
          call ice_dfsn( dt, T_in_K, rcm, p, exner, rhot, rtm_mc )
          thlm_mc = - ( Lv/(Cp*exner) ) * rtm_mc
        end if

#ifdef STATS
        if ( lstats_samp ) then

          if ( iNcm > 0 ) then           ! Brian
            zt%x(:,iNcm) = zt%x(:,iNcm) + Ncm
            zt%n(:,iNcm) = zt%n(:,iNcm) + 1
          end if
          if ( iNcnm > 0 ) then
            zt%x(:,iNcnm) = zt%x(:,iNcnm) + Ncnm
            zt%n(:,iNcnm) = zt%n(:,iNcnm) + 1
          end if
          if ( iNim > 0 ) then
            zt%x(:,iNim) = zt%x(:,iNim) + Nim
            zt%n(:,iNim) = zt%n(:,iNim) + 1
          end if

          if ( irrm > 0 ) then     ! Brian
            zt%x(:,irrm) = zt%x(:,irrm) + hydromet(:,1)
            zt%n(:,irrm) = zt%n(:,irrm) + 1
          end if

          if ( iNrm > 0 ) then     ! Brian
            zt%x(:,iNrm) = zt%x(:,iNrm) + hydromet(:,2)
            zt%n(:,iNrm) = zt%n(:,iNrm) + 1
          end if

          if ( irsnowm > 0 ) then
            zt%x(:,irsnowm) = zt%x(:,irsnowm) + hydromet(:,3)
            zt%n(:,irsnowm) = zt%n(:,irsnowm) + 1
          end if

          if ( iricem > 0 ) then
            zt%x(:,iricem) = zt%x(:,iricem) + hydromet(:,4)
            zt%n(:,iricem) = zt%n(:,iricem) + 1
          end if

          if ( irgraupelm > 0 ) then
            zt%x(:,irgraupelm) = zt%x(:,irgraupelm) + hydromet(:,5)
            zt%n(:,irgraupelm) = zt%n(:,irgraupelm) + 1
          end if

          if ( ithlm_mc > 0 ) then
            zt%x(:,ithlm_mc) 
     .      = zt%x(:,ithlm_mc) + thlm_mc(:)
            zt%n(:,ithlm_mc) = zt%n(:,ithlm_mc) + 1
          end if

          if ( irtm_mc > 0 ) then
            zt%x(:,irtm_mc) 
     .      = zt%x(:,irtm_mc) + rtm_mc(:)
            zt%n(:,irtm_mc) = zt%n(:,irtm_mc) + 1
          end if
        end if
#endif /*STATS*/

        rtm_forcing  = rtm_forcing + rtm_mc
        thlm_forcing = thlm_forcing + thlm_mc


#ifdef STATS
        if ( lstats_samp .and. ( lcoamps_micro .or. kk_rain ) ) then
          ! Rainfall rate (mm/day) should be defined on thermodynamic
          ! levels.  -Brian
          ! The absolute value of Vrr is taken because rainfall rate
          ! is a scalar quantity, and is therefore positive.
          if ( irain_rate > 0 ) then
            ! Rainfall rate (mm/day)
            zt%x(:,irain_rate) = zt%x(:,irain_rate) + 
     .        ( rrm * zm2zt( abs( Vrr ) ) ) * ( rhot / rho_lw )
     .        * ( 86400.0 * 1000.0 )
            zt%n(:,irain_rate) = zt%n(:,irain_rate) + 1
          end if

          ! Precipitation Flux (W/m^2) should be defined on
          ! momentum levels.  -Brian
          ! Normally, a flux is a vector quantity.  Since rain obviously
          ! falls downward, the sign of the flux would normally be negative.
          ! However, it is generally a convention in meteorology to show
          ! Precipitation Flux as a positive downward quantity.  Thus, the
          ! absolute value of vrr is taken.
          if ( iFprec > 0 ) then    ! Brian
            zm%x(:,iFprec) = zm%x(:,iFprec) +
     .      ( zt2zm( rrm ) * abs( Vrr ) ) 
     .      * ( rhom / rho_lw ) * rho_lw * Lv
            zm%n(:,iFprec) = zm%n(:,iFprec) + 1
          end if

          ! Store values of surface fluxes for statistics
          ! See notes above.
          if ( irain > 0 ) then   ! Added by Brian
            sfc%x(1,irain) = sfc%x(1,irain) +
     .        ( rrm(2) * abs( zm2zt( Vrr, 2 ) ) ) * ( rhot(2) / rho_lw )
     .        * ( 86400.0 * 1000.0 ) ! Rainfall rate (mm/day)
            sfc%n(1,irain) = sfc%n(1,irain) + 1
          end if
          if ( ipflux > 0 ) then   ! Added by Brian
            sfc%x(1,ipflux) 
     .      = sfc%x(1,ipflux) + ( zt2zm( rrm, 1 )
     .           * abs( Vrr(1) ) ) * ( rhom(1) / rho_lw ) * rho_lw * Lv
            sfc%n(1,ipflux) = sfc%n(1,ipflux) + 1
          end if
          ! Also store the value of surface rain water mixing ratio.
          if ( irrm_sfc > 0 ) then   ! Added by Brian
            sfc%x(1,irrm_sfc)
     .         = sfc%x(1,irrm_sfc) + ( zt2zm( rrm, 1 ) )
            sfc%n(1,irrm_sfc) = sfc%n(1,irrm_sfc) + 1
          end if

        end if ! lstats_samp
#endif /*STATS*/

!       Error Report
!       Joshua Fasching Feb 2008
        
        if ( lapack_error(err_code) ) then
                
           write(fstderr,*) "Error in advance_microphys"
           
           write(fstderr,*) "Intent(in)"
           
           write(fstderr,*) "thlm = ", thlm
           write(fstderr,*) "p = ", p
           write(fstderr,*) "exner = ", exner
           write(fstderr,*) "rhot = ", rhot
           write(fstderr,*) "rhom = ", rhom
           write(fstderr,*) "rtm = ", rtm
           write(fstderr,*) "rcm = ", rcm
           write(fstderr,*) "wmt = ", wmt
           write(fstderr,*) "wmm = ", wmm
           write(fstderr,*) "Khm = ", Khm
           write(fstderr,*) "Akm_est = ", Akm_est
           write(fstderr,*) "Akm = ", Akm
           write(fstderr,*) "pdf_parms = ", pdf_parms
           
           write(fstderr,*) "Intent(inout)"
           
           write(fstderr,*) "Ncm = ", Ncm
           write(fstderr,*) "Ncnm = ", Ncnm
           write(fstderr,*) "Nim = ", Nim
           write(fstderr,*) "hydromet = ", hydromet
           write(fstderr,*) "Intent(out)"
           write(fstderr,*) "rtm_mc = ", rtm_mc
           write(fstderr,*) "thlm_mc = ", thlm_mc
           
        endif
        
        return

        end subroutine advance_microphys
!-----------------------------------------------------------------------
        subroutine microphys_solve( solve_type, dt, lhs,
     .                              xrm_tndcy, xrm, err_code )

!       Description:

!       References:
!-----------------------------------------------------------------------
        use grid_class, only:
     .      gr ! Variable(s)
        use stats_prec, only: 
     .      time_prec ! Variable(s)
        use lapack_wrap, only: 
     .      tridag_solve ! Procedure(s)
!    .     ,band_solve
#ifdef STATS
        use stats_hoc, only:
     .      zt, ! Variable(s)
     .      irrm_ma,
     .      irrm_sd,
     .      irrm_dff,
     .      iNrm_ma,
     .      iNrm_sd,
     .      iNrm_dff,
     .      iricem_ma,
     .      iricem_sd,
     .      iricem_dff,
     .      irsnowm_ma,
     .      irsnowm_sd,
     .      irsnowm_dff,
     .      irgraupelm_ma,
     .      irgraupelm_sd,
     .      irgraupelm_dff,
     .      lstats_samp,
     .      ztscr01,
     .      ztscr02,
     .      ztscr03,
     .      ztscr04,
     .      ztscr05,
     .      ztscr06,
     .      ztscr07,
     .      ztscr08,
     .      ztscr09
#endif

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real(kind=time_prec), intent(in) :: dt ! Timestep     [s]

        ! Tendency computed tendency from COAMPS routine adjtq 
        ! or Brian Griffin's K & K microphysics implementation
        real, intent(in), dimension(gr%nnzp) ::
     .  xrm_tndcy !                                     [units/s]

        ! Input/Output Variables
        real, intent(inout), dimension(3,gr%nnzp) ::
     .  lhs ! Left hand side

        real, intent(inout), dimension(gr%nnzp) ::
     .  xrm ! Hydrometeor being solved for              [units vary]

        ! Output Variables
        integer, intent(out) :: err_code

        ! Local Variables
        real, dimension(gr%nnzp) ::
     .  rhs ! Right hand side

        integer :: k, kp1, km1 ! Array indices

#ifdef STATS
        integer ::
     .  ixrm_ma, ! Mean advection budget stats toggle
     .  ixrm_sd, ! Sedimentation budget stats toggle
     .  ixrm_dff ! Diffusion budget stats toggle

        select case( solve_type )
        case( "rrm" )
          ixrm_ma  = irrm_ma
          ixrm_sd  = irrm_sd
          ixrm_dff = irrm_dff
        case( "Nrm" )
          ixrm_ma  = iNrm_ma
          ixrm_sd  = iNrm_sd
          ixrm_dff = iNrm_dff
        case( "ricem" )
          ixrm_ma  = iricem_ma
          ixrm_sd  = iricem_sd
          ixrm_dff = iricem_dff
        case( "rsnowm" )
          ixrm_ma  = irsnowm_ma
          ixrm_sd  = irsnowm_sd
          ixrm_dff = irsnowm_dff
        case( "rgraupelm" )
          ixrm_ma  = irgraupelm_ma
          ixrm_sd  = irgraupelm_sd
          ixrm_dff = irgraupelm_dff
        case default
          ixrm_ma  = 0
          ixrm_sd  = 0
          ixrm_dff = 0
        end select

#endif /*STATS*/


        ! RHS of equation, following Brian's method from 
        ! the rain subroutine
        rhs(2:gr%nnzp-1) 
     .  = (xrm(2:gr%nnzp-1) / dt ) ! Time tendency
     .  + xrm_tndcy(2:gr%nnzp-1)

        ! Boundary condition on the RHS
        rhs(1) = xrm(1) / dt
        rhs(gr%nnzp) = ( xrm(gr%nnzp) / dt ) + xrm_tndcy(gr%nnzp-1)

        ! Solve system using tridag_solve. This uses LAPACK sgtsv,
        ! which relies on Gaussian elimination to decompose the matrix.
        call tridag_solve
     .       ( solve_type, gr%nnzp, 1, lhs(1,:), lhs(2,:), lhs(3,:),
     .         rhs, xrm, err_code )

        ! Alternative: Use LU decomposition instead.
        ! This doesn't seem to change the answer too much.
!       call band_solve
!    .       ( solve_type, 1, 1, gr%nnzp, 1, 
!    .         lhs, rhs, xrm, isValid )

#ifdef STATS
         if ( lstats_samp ) then
           do k = 1, gr%nnzp, 1

             km1 = max( k-1, 1 )
             kp1 = min( k+1, gr%nnzp )

             ! Finalize implicit contributions
             if ( ixrm_ma > 0 ) then
               zt%x(k,ixrm_ma) = zt%x(k,ixrm_ma)
     .           + ztscr01(k) * xrm(km1)
     .           + ztscr02(k) * xrm(k)
     .           + ztscr03(k) * xrm(kp1)
               zt%n(k,ixrm_ma) = zt%n(k,ixrm_ma) + 1
             end if

             if ( ixrm_sd > 0 ) then
               zt%x(k,ixrm_sd) = zt%x(k,ixrm_sd)
     .           + ztscr04(k) * xrm(km1)
     .           + ztscr05(k) * xrm(k)
     .           + ztscr06(k) * xrm(kp1)
               zt%n(k,ixrm_sd) = zt%n(k,ixrm_sd) + 1
             end if

             if ( ixrm_dff > 0 ) then
               zt%x(k,ixrm_dff) = zt%x(k,ixrm_dff)
     .           + ztscr07(k) * xrm(km1)
     .           + ztscr08(k) * xrm(k)
     .           + ztscr09(k) * xrm(kp1)
               zt%n(k,ixrm_dff) = zt%n(k,ixrm_dff) + 1
             end if

           end do ! 1..gr%nnzp
        end if ! lstats_samp

#endif /*STATS*/

        ! Boundary conditions on results
        !xrm(1) = xrm(2)
! Michael Falk, 7 Sep 2007, made this change to eliminate problems
! with anomalous rain formation at the top boundary.
!        xrm(gr%nnzp) = 0
        !xrm(gr%nnzp) = xrm(gr%nnzp-1)
! eMFc

        return
        end subroutine microphys_solve

!-----------------------------------------------------------------------
        subroutine microphys_lhs
     .             ( solve_type, lsed, dt, Kr, nu, wmt, v_t, xrm, lhs )

!       Description:
!       Setup the matrix of implicit contributions to a term
!       Includes the effects of sedimentation, diffusion, and advection.
!
!       Notes:
!       Setup for tridiagonal system and boundary conditions should be
!       the same as the original rain subroutine code.
!-----------------------------------------------------------------------
        use grid_class, only: 
     .      gr, ! Variable(s)
     .      zm2zt ! Procedure(s)
        use stats_prec, only: 
     .      time_prec ! Variable(s)
        use diffusion, only: 
     .      diffusion_zt_lhs ! Procedure(s)
        use mean_adv, only: 
     .      term_ma_zt_lhs ! Procedure(s)
#ifdef STATS
        use stats_hoc, only:
     .      irrm_bt,  ! Variable(s)
     .      irrm_ma,
     .      irrm_sd,
     .      irrm_dff,
     .      iNrm_bt, 
     .      iNrm_ma,
     .      iNrm_sd,
     .      iNrm_dff,
     .      iricem_bt, 
     .      iricem_ma,
     .      iricem_sd,
     .      iricem_dff,
     .      irsnowm_bt, 
     .      irsnowm_ma,
     .      irsnowm_sd,
     .      irsnowm_dff,
     .      irgraupelm_bt, 
     .      irgraupelm_ma,
     .      irgraupelm_sd,
     .      irgraupelm_dff,
     .      zt,
     .      ztscr01,
     .      ztscr02,
     .      ztscr03,
     .      ztscr04,
     .      ztscr05,
     .      ztscr06,
     .      ztscr07,
     .      ztscr08,
     .      ztscr09,
     .      lstats_samp
#endif
        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        logical, intent(in) :: 
     .  lsed ! Whether to add a sedimentation term

        real(kind=time_prec), intent(in) :: 
     .  dt  ! Timestep                  [s]

        real, intent(in) :: 
     .  nu  ! Diffusion coefficient     [m^2/s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  xrm, ! Value of hydrometeor as prev. time       [units vary]
     .  wmt, ! w wind on t                              [m/s]
     .  v_t, ! Terminal velocity                        [m/s]
     .  Kr   ! Eddy diffusivity for hydromet. (m-lev)   [m^2/s]

        real, intent(out), dimension(3,gr%nnzp) ::
     .  lhs ! Left hand side of tridiagonal matrix

        ! Local Variables
        real, dimension(3) :: tmp
        integer :: k, kp1, km1

#ifdef STATS
        integer ::
     .  ixrm_bt, ! Change in xrm over time stats toggle
     .  ixrm_ma, ! Mean advection budget stats toggle
     .  ixrm_sd, ! Sedimentation budget stats toggle
     .  ixrm_dff ! Diffusion budget stats toggle

        select case( solve_type )
        case( "rrm" )
          ixrm_bt  = irrm_bt 
          ixrm_ma  = irrm_ma
          ixrm_sd  = irrm_sd
          ixrm_dff = irrm_dff
        case( "Nrm" )
          ixrm_bt  = iNrm_bt 
          ixrm_ma  = iNrm_ma
          ixrm_sd  = iNrm_sd
          ixrm_dff = iNrm_dff
        case( "ricem" )
          ixrm_bt  = iricem_bt 
          ixrm_ma  = iricem_ma
          ixrm_sd  = iricem_sd
          ixrm_dff = iricem_dff
        case( "rsnowm" )
          ixrm_bt  = irsnowm_bt 
          ixrm_ma  = irsnowm_ma
          ixrm_sd  = irsnowm_sd
          ixrm_dff = irsnowm_dff
        case( "rgraupelm" )
          ixrm_bt  = irgraupelm_bt 
          ixrm_ma  = irgraupelm_ma
          ixrm_sd  = irgraupelm_sd
          ixrm_dff = irgraupelm_dff
        case default
          ixrm_bt  = 0
          ixrm_ma  = 0
          ixrm_sd  = 0
          ixrm_dff = 0
        end select

        if ( lstats_samp .and. ixrm_bt > 0 ) then
          zt%x(:,ixrm_bt) = zt%x(:,ixrm_bt) - ( xrm / dt )
        end if
#endif /*STATS*/

        ! Reset LHS Matrix for current timestep.
        lhs = 0.0

        ! Setup LHS Matrix
        do k = 1, gr%nnzp, 1
          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          ! Main diagonal
          ! Time Tendency
          lhs(kdiag,k) = lhs(kdiag,k) + ( 1.0 / dt )

          ! All diagonals

          ! Diffusion
          lhs(kp1diag:km1diag,k)
     .    = lhs(kp1diag:km1diag,k)
     .    + diffusion_zt_lhs( Kr(k), Kr(km1), nu, 
     .                        gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

          ! Mean Advection
          lhs(kp1diag:km1diag,k)
     .    = lhs(kp1diag:km1diag,k)
     .    + term_ma_zt_lhs( wmt(k), gr%dzt(k), k )

          ! Sedimentation
          ! Note: originally pristine ice did not sediment, so it was
          ! setup to be disabled as needed, but now pristine ice does
          ! sediment in the COAMPS case. -dschanen 12 Feb 2007
          if ( lsed ) then
             lhs(kp1diag:km1diag,k)
     .       = lhs(kp1diag:km1diag,k)
     .       + sedimentation( v_t(k), v_t(km1), zm2zt(v_t,k), 
     .                        gr%dzt(k), k )
          endif

         ! Implicit contributions to xrm
#ifdef STATS
          if ( lstats_samp ) then

            if ( ixrm_ma > 0 ) then
              tmp(1:3) = term_ma_zt_lhs( wmt(k), gr%dzt(k), k )
              ztscr01(k) = -tmp(3)
              ztscr02(k) = -tmp(2)
              ztscr03(k) = -tmp(1)
            end if

            if ( ixrm_sd > 0 ) then
              tmp(1:3) = sedimentation( v_t(k), v_t(km1), zm2zt(v_t,k), 
     .                                  gr%dzt(k), k )
              ztscr04(k) = -tmp(3)
              ztscr05(k) = -tmp(2)
              ztscr06(k) = -tmp(1)
            end if

            if ( ixrm_dff > 0 ) then
              tmp(1:3)
     .        = diffusion_zt_lhs( Kr(k), Kr(km1), nu, 
     .                            gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
              ztscr07(k) = -tmp(3)
              ztscr08(k) = -tmp(2)
              ztscr09(k) = -tmp(1)
            end if

          end if ! lstats_samp
#endif /*STATS*/

        end do ! 1..gr%nnzp

        return
        end subroutine microphys_lhs

!===============================================================================
        pure function sedimentation( V_hm, V_hmm1, V_hmzt, dzt, level )
     .  result( lhs ) 

!       Description:
!       Sedimentation of a hydrometeor:  implicit portion of the code.
!
!       The variable "r_hm" stands for one of the five hydrometeor 
!       variables currently in the code:  mean rain mixing ratio (rrm),
!       mean rain drop concentration (Nrm), mean ice mixing ratio 
!       (ricem), mean snow mixing ratio (rsnowm), or mean graupel mixing
!       ratio (rgraupelm).  The variable "V_hm" stands for the 
!       sedimentation velocity of the appropriate hydrometeor.
!
!       The d(r_hm)/dt equation contains a sedimentation term:
!
!       - d(V_hm*r_hm)/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - d( V_hm(t) * r_hm(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next
!       timestep, which is being advanced to in solving the d(r_hm)/dt
!       equation.
!
!       This term is discretized as follows:
!
!       The values of r_hm are found on the thermodynamic levels, while 
!       the values of V_hm are found on the momentum levels.  The 
!       variable r_hm is interpolated to the intermediate momentum 
!       levels.  At the intermediate momentum levels, the interpolated 
!       values of r_hm are multiplied by the values of V_hm.  Then, the
!       derivative of (r_hm*V_hm) is taken over the central thermodynamic
!       level.
!
!       -----r_hmp1---------------------------------------------- t(k+1)
!
!       ==============r_hm(interp)====V_hm======================= m(k)
!
!       -----r_hm--------------------------------d(V_hm*r_hm)/dz- t(k)
!
!       ==============r_hm(interp)====V_hmm1===================== m(k-1)
!
!       -----r_hmm1---------------------------------------------- t(k-1)
!
!       The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1)
!       correspond with altitudes zt(k+1), zm(k), zt(k), zm(k-1),
!       and zt(k-1), respectively.  The letter "t" is used for
!       thermodynamic levels and the letter "m" is used for momentum
!       levels.
!
!       dzt(k) = 1 / ( zm(k) - zm(k-1) )
!
!
!       Conservation Properties:
!
!       When a hydrometeor is sedimented to the ground (or to the lower
!       boundary of the model), it is removed from the atmosphere (or
!       from the model domain).  Thus, the quantity of the hydrometeor 
!       over the entire vertical domain should not be conserved due to 
!       the process of sedimentation.  Thus, not all of the column 
!       totals in the left-hand side matrix should be equal to 0.
!
!       Left-hand side matrix contributions from eddy diffusion term;
!       first four vertical levels:
!
!            ------------------------------------------------------------------------------------------>
!       k=1 |            0                            0                                0
!           |
!       k=2 | -dzt(k)*V_hm(k-1)*D(k)       +dzt(k)*[ V_hm(k)*B(k)
!           |                                       -V_hm(k-1)*C(k) ]       +dzt(k)*V_hm(k)*A(k)
!           |
!       k=3 |            0                 -dzt(k)*V_hm(k-1)*D(k)           +dzt(k)*[ V_hm(k)*B(k)
!           |                                                                        -V_hm(k-1)*C(k) ]
!           |
!       k=4 |            0                            0                     -dzt(k)*V_hm(k-1)*D(k)
!           |
!          \ /
!
!       The variables A(k), B(k), C(k), and D(k) are weights of interpolation 
!       around the central thermodynamic level (k), such that:
!       A(k) = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ),
!       B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ],
!       C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
!       D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ].
!
!       Note:  The superdiagonal term from level 3 and both the main
!              diagonal and superdiagonal terms from level 4 are not
!              shown on this diagram.

!       References:
!       None

!       Notes:
!       Both COAMPS Microphysics and Brian Griffin's implementation use
!       Khairoutdinov and Kogan (2000) for the calculation of rain mixing
!       ratio and rain droplet number concentration sedimentation
!       velocities.
!-----------------------------------------------------------------------

        use grid_class, only: 
     .      gr ! Variable(s)

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1_tdiag = 1,   ! Thermodynamic superdiagonal index.
     .  k_tdiag   = 2,   ! Thermodynamic main diagonal index.
     .  km1_tdiag = 3    ! Thermodynamic subdiagonal index.

        integer, parameter ::
     .  t_above = 1,   ! Index for upper thermodynamic level grid weight.
     .  t_below = 2    ! Index for lower thermodynamic level grid weight.

        ! Input Variables
        real, intent(in) ::
     .  V_hm,   ! Sedimentation velocity of hydrometeor (k)                [m/s]
     .  V_hmm1, ! Sedimentation velocity of hydrometeor (k-1)              [m/s]
     .  V_hmzt, ! Sed. vel. of hydromet. interpolated to thermo. level (k) [m/s]
     .  dzt     ! Inverse of grid spacing (k)                              [m] 

        integer, intent(in) :: 
     .  level ! Central thermodynamic level (on which calculation occurs).

        ! Return Variable
        real, dimension(3) :: lhs

        ! Local Variables
        integer ::
     .  mk,   ! Momentum level directly above central thermodynamic level.
     .  mkm1  ! Momentum level directly below central thermodynamic level.

        ! Momentum level (k) is between thermodynamic level (k+1)
        ! and thermodynamic level (k).
        mk   = level
        ! Momentum level (k-1) is between thermodynamic level (k)
        ! and thermodynamic level (k-1).
        mkm1 = level - 1

        ! Note:  The code has now been changed so that Vxx has been pulled 
        !        inside of the derivative.  The sedimentation term is now of 
        !        the form -d(Vxx*xx)/dt, rather than of the form -Vxx d(xx)/dt.
        !        The term has been re-discretized in a conservative manner and
        !        the results are listed below, with the old code commented out.

        if ( level == 1 ) then

           ! k = 1 (bottom level); lower boundary level; no effects.

           ! Thermodynamic superdiagonal: [ x r_hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) = 0.0

           ! Thermodynamic main diagonal: [ x r_hm(k,<t+1>) ]
           lhs(k_tdiag)   = 0.0

           ! Thermodynamic subdiagonal: [ x r_hm(k-1,<t+1>) ]
           lhs(km1_tdiag) = 0.0


        elseif ( level > 1 .and. level < gr%nnzp ) then

           ! Most of the interior model; normal conditions.

        ! Vince Larson pulled v_t inside derivative to make conservative. 
        ! 13 Dec 2007
!
!           ! Thermodynamic superdiagonal: [ x r_hm(k+1,<t+1>) ]
!           lhs(kp1_tdiag) 
!     .     = + V_hmzt * dzt * gr%weights_zt2zm(t_above,mk)
!
!           ! Thermodynamic main diagonal: [ x r_hm(k,<t+1>) ]
!           lhs(k_tdiag) 
!     .     = + V_hmzt * dzt * (   gr%weights_zt2zm(t_below,mk)
!     .                          - gr%weights_zt2zm(t_above,mkm1)  )
!
!           ! Thermodynamic subdiagonal: [ x r_hm(k-1,<t+1>) ]
!           lhs(km1_tdiag) 
!     .     = - V_hmzt * dzt * gr%weights_zt2zm(t_below,mkm1)

           ! Thermodynamic superdiagonal: [ x r_hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) 
     .     = + dzt * V_hm * gr%weights_zt2zm(t_above,mk)

           ! Thermodynamic main diagonal: [ x r_hm(k,<t+1>) ]
           lhs(k_tdiag) 
     .     = + dzt * (   V_hm * gr%weights_zt2zm(t_below,mk)
     .                 - V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

           ! Thermodynamic subdiagonal: [ x r_hm(k-1,<t+1>) ]
           lhs(km1_tdiag) 
     .     = - dzt * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)

        !  End Vince Larson change


        elseif ( level == gr%nnzp ) then

           ! k = gr%nnzp (top level); upper boundary level; no flux.

        ! Vince Larson pulled v_t inside derivative to make conservative. 
        ! 13 Dec 2007
!
!           ! Thermodynamic superdiagonal: [ x r_hm(k+1,<t+1>) ]
!           lhs(kp1_tdiag) 
!     .     = + V_hmzt * dzt * gr%weights_zt2zm(t_above,mk)
!
!           ! Thermodynamic main diagonal: [ x r_hm(k,<t+1>) ]
!           lhs(k_tdiag) 
!     .     = + V_hmzt * dzt * (   gr%weights_zt2zm(t_below,mk)
!     .                          - gr%weights_zt2zm(t_above,mkm1)  )
!
!           ! Thermodynamic subdiagonal: [ x r_hm(k-1,<t+1>) ]
!           lhs(km1_tdiag) 
!     .     = - V_hmzt * dzt * gr%weights_zt2zm(t_below,mkm1)

           ! Thermodynamic superdiagonal: [ x r_hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) 
     .     = + dzt * V_hm * gr%weights_zt2zm(t_above,mk)

           ! Thermodynamic main diagonal: [ x r_hm(k,<t+1>) ]
           lhs(k_tdiag) 
     .     = + dzt * (   V_hm * gr%weights_zt2zm(t_below,mk)
     .                 - V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

           ! Thermodynamic subdiagonal: [ x r_hm(k-1,<t+1>) ]
           lhs(km1_tdiag) 
     .     = - dzt * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)

        !  End Vince Larson change


        endif

        return
        end function sedimentation

!===============================================================================
        subroutine adj_microphys_tndcy( xrm_tndcy, wmt, v_t, Kr, nu,
     .                                  dt, level, lsed, solve_type,
     .                                  xrm, overevap_rate )

        ! DESCRIPTION:  Correction for the over-evaporation of a hydrometeor.
        !
        ! If a small amount of a hydrometeor (such as rain water) gets 
        ! diffused into an area that is very dry (such as right above the 
        ! cloud top), the hydrometeor (rain water) will have a very high rate
        ! of evaporation and will evaporate entirely in a short amount of time.
        ! However, the evaporation rate is computed instantaneously at a given
        ! moment in time.  This rate is then projected over the entire length 
        ! of the given timestep.  Therefore, a high-enough rate of evaporation
        ! combined with a small-enough amount of the hydrometeor (rain water)
        ! and a long-enough timestep will cause the hydrometeor value (rain 
        ! water mixing ratio) to be negative by the end of the timestep.  
        ! Therefore, a correction factor needs to be imposed on the 
        ! evaporation rate so that the amount of the hydrometeor (rain water 
        ! mixing ratio) does not fall below 0.
        !
        ! Besides over-evaporation of a hydrometeor, other factors may come
        ! into play that cause the value of a hydrometeor to fall below 0.
        ! These factors are due to the nature of implicit discretization
        ! and numerical errors.  In a nutshell, the eddy diffusion parameter
        ! used currently in this model smooths out the entire hydrometeor
        ! profile as a whole at every timestep.  This smoothing may cause
        ! negative values at certain levels.  Also, mean advection and 
        ! hydrometeor sedimentation can cause negative values to occur in
        ! the hydrometeor.  This can happen in places where the profile
        ! abruptly goes from a large positive value to 0 (such as at cloud 
        ! top).  The nature of the discretization of taking a derivative at
        ! these levels may cause negative values of a hydrometeor.
        !
        ! This subroutine is called only if a hydrometeor at a certain level 
        ! contains a negative value.  First, this subroutine uses the same
        ! methods that the model statistical code uses in computing budget terms
        ! in order to determine what factors effected the value of the given
        ! hydrometeor during the timestep that was just solved for.  The
        ! mean advection, sedimentation, and diffusion budget terms are all 
        ! computed.  These three terms are then added together to make up the
        ! total transport and sedimentation tendency.  This tendency is then
        ! added to the total microphysical tendency to find the overall
        ! hydrometeor tendency.  The overall hydrometeor tendency is then
        ! multiplied by the timestep length to find the net change in the
        ! hydrometeor over the last timestep.  This net change is then added
        ! to the current value of the hydrometeor in order to find the value
        ! of the hydrometeor at the previous timestep.  This method has been
        ! well tested and produces accurate results.
        !
        ! Once the value of the hydrometeor at the previous timestep has been
        ! found, the net change in the hydrometeor due to ONLY mean advection,
        ! diffusion, and sedimentation is calculated.  This net change is
        ! added to the value of the hydrometeor at the previous timestep.  If
        ! the new value is below zero, then the negative value of the
        ! hydrometeor was caused by the mean advection, diffusion, and 
        ! sedimentation terms.  The microphysical terms (evaporation) did not
        ! cause the negative value.  There was no over-evaporation and the
        ! evaporation rate can be set to 0.  However, if the new value of
        ! the hydrometeor is greater than or equal to 0, then the microphysical
        ! tendencies (evaporation) did cause the hydrometeor array to have
        ! negative values.  The amount of hydrometeor evaporated is set equal
        ! to the amount that was left-over after the transport and 
        ! sedimentation effects were added in.  The evaporation rate is that
        ! amount divided by the timestep.  This can be viewed as the
        ! timestep-average evaporation rate, whereas the rate previously
        ! calculated can be viewed as the instantaneous evaporation rate.
        ! The amount of the hydrometeor that was over-evaporated is the amount
        ! of the hydrometeor that is negative.  The over-evaporation rate is
        ! that amount divided by the timestep.
        !
        ! It should be noted that this is important because the rain water 
        ! mixing ratio time tendency (drr/dt) due to microphysics at every 
        ! level is incorporated into the total water mixing ratio (rtm) and 
        ! liquid water potential temperature (thlm) equations.  Any artificial
        ! excess in evaporation will artificially increase water vapor, and 
        ! thus rtm, and artificially decrease thlm (due to evaporative 
        ! cooling).  This may result in an artificial increase in cloud water.
        !
        ! rrm_mc_tndcy = rrm_cond + rrm_auto + rrm_accr
        ! rtm_mc  = - rrm_mc_tndcy
        ! thlm_mc = ( Lv / (Cp*exner) ) * rrm_mc_tndcy
        !
        ! anyplace where rrm drops below zero due to microphysics, there is 
        ! too much evaporation rate for the timestep, so rrm_cond is too 
        ! negative.  We must add in the over-evaporated amount of rrm/dt to 
        ! make the rate accurate.  The over-evaporated amount is being defined
        ! as a positive scalar, so that:  overevap_rrm = -rrm (where rrm < 0)
        ! -- this makes overevap_rrm positive.
        !
        ! New cond/evap rate = rrm_cond + overevap_rrm/dt
        ! (overevap_rate = overevap_rrm/dt)
        ! -- since rrm_cond can only be negative (we don't allow rain droplets
        !    to grow by condensation) and overevap_rrm/dt can only be positive
        !    (we define it that way), the new cond/evap rate will be less
        !    negative, which is what we want.
        !
        ! To update the effects of microphysics on rtm and thl:
        !
        ! rtm_mc = rtm_mc - overevap_rate
        ! thlm_mc = thlm_mc + ( Lv / (Cp*exner) ) * overevap_rate
        !
        ! This is done in the subroutine which calls this one.
        ! 
        ! If the hydrometeor is negative due to reasons besides 
        ! over-evaporation, the value is clipped.  This is statistically
        ! stored in the clipping array.  This is also done in the subroutine
        ! which calls this one.
        !
        ! Brian Griffin.

        use grid_class, only: 
     .      gr, ! Variable(s) 
     .      zm2zt ! Procedure(s)
        use stats_prec, only:
     .      time_prec ! Variable(s)
        use diffusion, only: 
     .      diffusion_zt_lhs ! Procedure(s)
        use mean_adv, only: 
     .      term_ma_zt_lhs ! Procedure(s)

!#ifdef STATS
!        use stats_hoc
!#endif /*STATS*/
! Joshua Fasching 2007

        ! Input variables.

        ! Hydrometeor microphysical tendency.
        real, dimension(gr%nnzp), intent(in) :: xrm_tndcy  ! [hm_units/s]
        ! Vertical velocity (thermo. levels).
        real, dimension(gr%nnzp), intent(in) :: wmt        ! [m/s]
        ! Sedimentation velocity (interpolated to thermo. levels).
        real, dimension(gr%nnzp), intent(in) :: v_t        ! [m/s]
        ! Eddy diffusivity for hydrometeors (m-lev).
        real, dimension(gr%nnzp), intent(in) :: Kr         ! [m^2/s]

        real, intent(in) :: nu  ! Diffusion coefficient      [m^2/s]

        real(kind=time_prec), intent(in) :: dt  ! Timestep   [s]

        integer, intent(in) :: level  ! Vertical grid index

        logical, intent(in) :: lsed   ! Whether to add a sedimentation term

        character(len=*), intent(in) :: solve_type

        ! Input/output variable.

        ! Hydrometeor.
        real, dimension(gr%nnzp), intent(inout) :: xrm     ! [hm_units]

        ! Output variable.

        ! Excess evaporation rate.
        real, intent(out) :: overevap_rate                 ! [hm_units/s]

        ! Local variables.
        real :: ma_subdiag   ! Term to be multiplied by xrm(k-1) in m.a. eq.
        real :: ma_maindiag  ! Term to be multiplied by xrm(k) in m.a. eq.
        real :: ma_supdiag   ! Term to be multiplied by xrm(k+1) in m.a. eq.
        real :: sd_subdiag   ! Term to be multiplied by xrm(k-1) in sed. eq.
        real :: sd_maindiag  ! Term to be multiplied by xrm(k) in sed. eq.
        real :: sd_supdiag   ! Term to be multiplied by xrm(k+1) in sed. eq.
        real :: df_subdiag   ! Term to be multiplied by xrm(k-1) in diff. eq.
        real :: df_maindiag  ! Term to be multiplied by xrm(k) in diff. eq.
        real :: df_supdiag   ! Term to be multiplied by xrm(k+1) in diff. eq.

        real :: ma_tndcy     ! Mean advection tendency  [hm_units/s]
        real :: sd_tndcy     ! Sedimentation tendency   [hm_units/s]
        real :: df_tndcy     ! Diffusion tendency       [hm_units/s]

        real :: trnsprt_sed_tndcy ! Total transport and sedimentation tendency.
        real :: tot_tndcy         ! Overall hydrometeor total tendency.
        real :: xrm_chge          ! Total change in hydrometeor over last t.s.
        real :: xrm_old           ! Value of hydrometeor at previous timestep.
        real :: xrm_chge_trsed    ! Net change in hm. due to only transport/sed.
        real :: xrm_trsed_only    ! New hm. val. due only to transport/sed.

        real :: evap_amt          ! The actual evaporation amount over the t.s.
        real :: evap_rate         ! The time-averaged rate.
        real :: overevap_amt      ! The amount of h.m. that was over-evap.

        real, dimension(1:3) :: tmp

        integer :: k, km1, kp1

!#ifdef STATS
!       integer ::
!    .  ixrm_cond_adj  ! Adjustment to xrm evaporation rate due to over-evap.

!       select case( solve_type )
!       case( "rrm" )
!         ixrm_cond_adj  = irrm_cond_adj 
!       case( "Nrm" )
!         ixrm_cond_adj  = iNrm_cond_adj 
!       end select
!#endif /*STATS*/
! Joshua Fasching 2007

        k = level
        km1 = max(k-1,1)
        kp1 = min(k+1,gr%nnzp)


        ! Mean advection tendency component

        ! The implicit (LHS) value of the mean advection component of the
        ! equation used during the timestep that was just solved for.
        tmp(1:3) = term_ma_zt_lhs( wmt(k), gr%dzt(k), k )

        ma_subdiag  = -tmp(3) ! subdiagonal
        ma_maindiag = -tmp(2) ! main diagonal
        ma_supdiag  = -tmp(1) ! superdiagonal

        ma_tndcy =
     .     + ma_subdiag  * xrm(km1)
     .     + ma_maindiag * xrm(k)
     .     + ma_supdiag  * xrm(kp1)


        ! Sedimentation tendency component

        if ( lsed ) then

           ! The implicit (LHS) value of the sedimentation component of the
           ! equation used during the timestep that was just solved for.
           tmp(1:3) = sedimentation( v_t(k), v_t(km1), zm2zt(v_t,k), 
     .                               gr%dzt(k), k )

           sd_subdiag  = -tmp(3) ! subdiagonal
           sd_maindiag = -tmp(2) ! main diagonal
           sd_supdiag  = -tmp(1) ! superdiagonal

           sd_tndcy = 
     .        + sd_subdiag  * xrm(km1)
     .        + sd_maindiag * xrm(k)
     .        + sd_supdiag  * xrm(kp1)

        else

           sd_tndcy = 0.0

        endif


        ! Diffusion tendency component

        ! The implicit (LHS) value of the diffusion component of the
        ! equation used during the timestep that was just solved for.
        tmp(1:3)
     .     = diffusion_zt_lhs( Kr(k), Kr(km1), nu, 
     .                         gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

        df_subdiag  = -tmp(3) ! subdiagonal
        df_maindiag = -tmp(2) ! main diagonal
        df_supdiag  = -tmp(1) ! superdiagonal

        df_tndcy = 
     .     + df_subdiag  * xrm(km1)
     .     + df_maindiag * xrm(k)
     .     + df_supdiag  * xrm(kp1)


        ! Total transport and sedimentation tendency
        trnsprt_sed_tndcy = ma_tndcy + df_tndcy + sd_tndcy

        ! Overall hydrometeor tendency
        tot_tndcy = trnsprt_sed_tndcy + xrm_tndcy(k)

        ! The net amount of change in the hydrometeor over the last timestep.
        xrm_chge = tot_tndcy * dt

        ! The value of xrm at the previous timestep.
        xrm_old = xrm(k) - xrm_chge

        ! The net amount of change in the hydrometeor due to only the
        ! transport (mean advection and diffusion) and sedimentation terms.
        xrm_chge_trsed = trnsprt_sed_tndcy * dt

        ! The new value of the hydrometeor at this timestep due to only
        ! the transport and sedimentation terms.
        xrm_trsed_only = xrm_old + xrm_chge_trsed

        if ( xrm_trsed_only >= 0.0 ) then
           ! The negative value of hydrometeor (xrm) is due ONLY to 
           ! microphysical tendencies, namely the over-evaporation of xrm.
           ! Find the actual amount of the hydrometeor that evaporated during
           ! the timestep to make the value of xrm go to 0.
           evap_amt = -xrm_trsed_only
           ! Divide by the timestep to find the actual evaporation rate.
           evap_rate = evap_amt / dt
           ! The amount of the hydrometeor that was artificially excessively
           ! evaporated.  Define as positive.
           overevap_amt = -xrm(k)
           ! Divide by the timestep to find the over-evaporation rate.
           ! Define as positive.  This rate should also be the difference
           ! between the computed evaporation rate (xrm_tndcy) and the
           ! actual evaporation rate (evap_rate).
           overevap_rate = overevap_amt / dt
           ! Reset the value of the hydrometeor (xrm) to 0.
           xrm(k) = 0.0
        else
           ! The negative value of hydrometeor (xrm) is due to transport 
           ! (mean advection and diffusion) and sedimentation.
           ! Find the actual amount of the hydrometeor that evaporated during
           ! the timestep to make the value of xrm go to 0.  Even though the
           ! microphysical tendency portion of the code may have computed an 
           ! evaporation rate, we figure that the transport and sedimentation 
           ! terms made the value of the hydrometeor negative, so we say that 
           ! the evaporation amount and rate is 0.
           evap_amt = 0.0
           evap_rate = 0.0
           ! The amount of the hydrometeor that was artificially excessively
           ! evaporated.  Define as positive.  In this case, any evaporation
           ! that was computed is considered to be over-evaporation.
           ! Define as positive.
           overevap_amt = -xrm_tndcy(k) * dt
           overevap_rate = -xrm_tndcy(k)
           ! Currently reset xrm to xrm_trsed_only.
           ! This is done to make the statistical budget for xrm balance
           ! correctly.  The value of xrm(k) will still be negative at this
           ! this point.  However, it will be less negative because it has been
           ! adjusted for over-evaporation.  The remaining negative value of
           ! hydrometeor xrm, which is due to transport and sedimentation, will
           ! be zeroed out in clipping in the subroutine that calls this one.
           xrm(k) = xrm_trsed_only
        endif

        end subroutine adj_microphys_tndcy

!===============================================================================

        end module microphys_driver
