!-----------------------------------------------------------------------
! $Id: microphysics.F,v 1.9 2007-06-19 21:06:36 dschanen Exp $
        module microphysics

!       Description:
!       Call a microphysical scheme to compute hydrometeor species,
!       and advect, sediment, & diffuse using a tridiagonal system.

!       References:
!       None
!-----------------------------------------------------------------------
        implicit none

        public :: timestep_microphys

        ! Subroutines
        private :: microphys_lhs, microphys_solve

        ! Functions
        private :: sedimentation, advection, diffusion

        private ! Default

        contains

!-----------------------------------------------------------------------
        subroutine timestep_microphys
     .             ( dt, time_current, time_initial, 
     .               thlm, p, exner, rhot, rhom, rtm, rcm, Ncm, 
     .               pdf_parms, wmt, wmm, AKm_est, AKm, 
     .               Ncnm, Nim,
     .               hydromet,
     .               hm_rt_tndcy, hm_thl_tndcy, isValid )

!       Description:
!       Compute pristine ice, snow, graupel, & rain hydrometeor fields.
!       Uses implicit discretization.

!       References:
!       None
!-----------------------------------------------------------------------
        use grid_class
        use rain_equations, only: kk_microphys
        use model_flags
        use parameters, only: nu_r
        use constants, only: hydromet_dim, Lv, Cp, rho_lw
#ifdef STATS
        use hoc_stats
#endif
        implicit none

        ! External
        external :: coamps_micro_driver

        ! Constant Parameters

        ! COAMPS parameters
        integer, parameter :: 
     .  nrdamp = 0, ! Number of points in the damping layer (COAMPS)
     .  icond  = 3  ! Autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none

        logical, parameter ::
     .  lice = .true. ! Whether to produce ice in COAMPS.
                      ! According to Jerry Schmidt of NRL,
                      !   if lice = .true., then we should
                      !   set ldrizzle = .false.
                      !   because collection of drizzle by ice
                      !   is not implemented yet.

        real, intent(in) :: 
     .  dt,          ! Timestep         [s]
     .  time_current,! Current time     [s]
     .  time_initial ! Initial time     [s]

        real, dimension(gr%nnzp), intent(in) ::
     .  thlm,   ! Liquid potential temp.        [K]
     .  p,      ! Pressure                      [Pa]
     .  exner,  ! Exner function                [-]
     .  rhot,   ! Density on thermo. grid       [kg/m^3]
     .  rhom,   ! Density on moment. grid       [kg/m^3]
     .  rtm,    ! Total water mixing ratio      [kg/kg]
     .  rcm,    ! Liquid water mixing ratio     [kg/kg]
     .  wmt,    ! w wind on moment. grid        [m/s]
     .  wmm,    ! w wind on thermo. grid        [m/s]
     .  Akm_est,! Analytic Kessler ac           [kg/kg]
     .  Akm     ! Analytic Kessler estimate     [kg/kg]

        ! Note:
        ! K & K only uses Ncm, while for COAMPS Ncnm is initialized
        ! and Nim & Ncm are computed within subroutine adjtg.
        real, dimension(gr%nnzp), intent(inout) ::
     .  Ncm,    ! Rain droplet number concentration     [count/m^3]
     .  Ncnm,   ! Cloud nuclei number concentration     [count/m^3]
     .  Nim     ! Ice crystal number concentration      [count/m^3]

        real, target, dimension(gr%nnzp,26), intent(in) ::
     .  pdf_parms     ! PDF parameters

        real, target, dimension(gr%nnzp,hydromet_dim), intent(inout) ::
     .  hydromet      ! Array of rain, prist. ice, graupel, etc. [units vary]

        real, dimension(gr%nnzp), intent(out) ::
     .  hm_rt_tndcy,  ! Change in total water due to hydrometeors    [(kg/kg)/s]
     .  hm_thl_tndcy  ! Change in liquid potential temperature 
                      ! due to hydrometeors                          [K/s]

        logical, intent(out) :: isValid

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
     .  rrm_tndcy,      ! Rain mixing ratio tendency      [(kg/kg)/s]
     .  Nrm_tndcy,      ! Rain number conc. tendency      [(count/m^3)/s]
     .  rsnowm_tndcy,   ! Snow mixing ratio tendency      [(kg/kg)/s]
     .  ricem_tndcy,    ! Ice mixing ratio tendency       [(kg/kg)/s]
     .  rgraupelm_tndcy ! Graupel mixing ratio tendency   [(kg/kg)/s]  

        real, dimension(1,1,gr%nnzp) ::
     .  cond ! COAMPS stat for condesation/evap of rcm

        real, pointer, dimension(:) ::
     .  rrm,      ! Pointer for rain water mixing ratio   [kg/kg]
     .  Nrm,      ! Pointer for rain droplet number conc. [count/m^3]
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

        ! Variables needed to handle correction to rtm and thlm microphysics
        ! tendency arrays, as well rrm_cond and Nrm_cond statistical
        ! tendency arrays, due to a negative result being produced by 
        ! over-evaporation of rain water over the course of a timestep.
        ! Brian Griffin.  April 14, 2007.
        real :: overevap_rrm, ! Absolute value of negative rrm [kg/kg]
     .          overevap_Nrm  ! Absolute value of negative Nrm [count/m^3]

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

        ! Begin by calling either Brian Griffin's implementation of the
        ! Khairoutdinov and Kogan microphysical scheme or 
        ! alternatively the Rutlege and Hobbes scheme from COAMPS(R).
        ! Note: COAMPS appears to have some K&K elements to it as well.

        if ( lcoamps_micro ) then

           call coamps_micro_driver
     .          ( gr%nnzp-1, gr%nnzp, time_current, dt,
     .            rtm, wmm, p, exner, rhot,
     .            thlm, ricem, rrm, rgraupelm, rsnowm,
     .            1, nrdamp, 1, 1, icond, 
     .            lice,
     .            rcm, Ncm, Nrm, Ncnm, Nim,
     .            cond, Vsnow, Vice, Vrr, VNr, Vgraupel,
     .            ricem_tndcy, rrm_tndcy, rgraupelm_tndcy, 
     .            rsnowm_tndcy, Nrm_tndcy,
     .            hm_rt_tndcy, hm_thl_tndcy )

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
     .         ( thlm, p, exner, rhot, thl1, thl2, a, rc1, rc2,
     .           s1, s2, ss1, ss2, rtm, rcm, Ncm, rrm, Nrm,
     .           .true., AKm_est, AKm,
     .           rrm_tndcy, Nrm_tndcy, 
     .           hm_rt_tndcy, hm_thl_tndcy,
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

        end if ! kk_rain

        ! Compute rrm and Nrm for either case

        ! ----------------------------------------------------------------------
        ! Rain water mixing ratio (rrm)
        ! ----------------------------------------------------------------------

        call microphys_lhs
     .       ( "rrm", .true., dt, nu_r, wmt, Vrr, rrm, lhs )

        call microphys_solve
     .       ( "rrm", dt, lhs, rrm_tndcy, rrm, isValid )

        ! Correction for the over-evaporation of rain.
        !
        ! If a small amount of rain water gets diffused into an area that is
        ! very dry (such as right above the cloud top), the rain water will
        ! have a very high rate of evaporation and will evaporate entirely in
        ! a short amount of time.  However, the evaporation rate is computed
        ! instantaneously at a given moment in time.  This rate is then
        ! projected over the entire length of the given timestep.  Therefore,
        ! a high-enough rate of evaporation combined with a small-enough amount
        ! of rain water and a long-enough timestep will cause the rain water
        ! mixing ratio to be negative by the end of the timestep.  Therefore,
        ! a correction factor needs to be imposed on the evaporation rate so 
        ! that the amount of rain water mixing ratio does not fall below 0.
        !
        ! This is accomplished by adjusting the evaporation rate at any level
        ! where the value of rain water mixing ratio has fallen below 0.  
        ! The negative amount of rain water is divided by the timestep length
        ! and then added back into the evaporation rate.  This is important
        ! because the rain water mixing ratio time tendency (drr/dt) at every
        ! level is incorporated into the total water mixing ratio (rtm) and 
        ! liquid water potential temperature (thlm) equations.  Any artificial
        ! excess in evaporation will artificially increase water vapor, and 
        ! thus rtm, and artificially decrease thlm (due to evaporative 
        ! cooling).  This may result in an artificial increase in cloud water.
        !
        ! rrm_mc_tndcy = rrm_cond + rrm_auto + rrm_accr
        ! hm_rt_tndcy  = - rrm_mc_tndcy
        ! hm_thl_tndcy = ( Lv / (Cp*exner) ) * rrm_mc_tndcy
        !
        ! anyplace where rrm drops below zero, there is too much evaporation
        ! rate for the timestep, so rrm_cond is too negative.  We must add in 
        ! the overevaporated amount of rrm/dt to make the rate accurate.  
        ! The overevaporated amount is being defined as a positive scalar, so
        ! that:  overevap_rrm = -rrm (where rrm < 0) 
        ! -- this makes overevap_rrm positive.
        !
        ! New cond/evap rate = rrm_cond + overevap_rrm/dt
        ! -- since rrm_cond can only be negative (we don't allow rain droplets
        !    to grow by condensation) and overevap_rrm/dt can only be positive
        !    (we define it that way), the new cond/evap rate will be less
        !    negative, which is what we want.
        !
        ! To update the effects of microphysics on rtm and thl:
        !
        ! hm_rt_tndcy = hm_rt_tndcy - overevap_rrm/dt
        ! hm_thl_tndcy = hm_thl_tndcy + ( Lv / (Cp*exner) ) * overevap_rrm/dt
        !
        ! Brian Griffin.  April 14, 2007.

        do k = 1, gr%nnzp, 1
           if ( rrm(k) < 0.0 ) then
!              if ( hm_rt_tndcy(k) > 0.0 ) then
!                 ! A positive hm_rt_tndcy means a negative rrm tendency due
!                 ! to microphysics.  The only way this can happen is if
!                 ! evaporation takes place.
                 overevap_rrm = -rrm(k)
!                 rrm(k) = 0.0
                 hm_rt_tndcy(k) = hm_rt_tndcy(k) - (overevap_rrm/dt)
                 hm_thl_tndcy(k) = hm_thl_tndcy(k)
     .                + ( Lv / ( Cp*exner(k) ) ) * (overevap_rrm/dt)
!#ifdef STATS
!                 ! Adjustment to condensation/evaporation
!                 if ( lstats_samp .and. irrm_cond > 0 ) then
!                    zt%x(k,irrm_cond) = zt%x(k,irrm_cond) 
!     .                                  + (overevap_rrm/dt)
!                    ! This is part of the same sample that was done in rain.F.
!                    ! Since it is only an adjustment, no new sample addition
!                    ! is needed.
!                 endif
!#endif
!              endif
           endif
        enddo

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
        endif
#endif

        ! ----------------------------------------------------------------------
        ! Rain droplet number concentration (Nrm)
        ! ----------------------------------------------------------------------

        call microphys_lhs
     .       ( "Nrm", .true., dt, nu_r, wmt, VNr, Nrm, lhs )

        call microphys_solve
     .       ( "Nrm", dt, lhs, Nrm_tndcy, Nrm, isValid )

        ! Handle over-evaporation similar to rrm.  However, in the case of Nrm
        ! there is no effect on rtm or on thlm.  Brian Griffin.  April 14, 2007.

        do k = 1, gr%nnzp, 1
           if ( Nrm(k) < 0.0 ) then
              overevap_Nrm = -Nrm(k)
!              Nrm(k) = 0.0
!#ifdef STATS
!              ! Adjustment to condensation/evaporation
!              if ( lstats_samp .and. iNrm_cond > 0 ) then
!                 zt%x(k,iNrm_cond) = zt%x(k,iNrm_cond) 
!     .                               + (overevap_Nrm/dt)
!                 ! This is part of the same sample that was done in rain.F.
!                 ! Since it is only an adjustment, no new sample addition
!                 ! is needed.
!              endif
!#endif
           endif
        enddo

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
          if (lstats_samp .and. irrm_cl > 0 ) then
            zt%x(:,iNrm_cl) = zt%x(:,iNrm_cl) + Nrm / dt
            zt%n(:,iNrm_cl) = zt%n(:,iNrm_cl) + 1
          end if
#endif

#ifdef STATS
        ! Finalize time tendency
        if ( lstats_samp .and. iNrm_bt > 0 ) then
          zt%x(:,iNrm_bt) = zt%x(:,iNrm_bt) + ( Nrm / dt )
          zt%n(:,iNrm_bt) = zt%n(:,iNrm_bt) + 1
        endif
#endif

        ! Compute rsnowm, rgraupelm, and ricem for COAMPS microphysics only.

        if ( lcoamps_micro ) then

          ! --------------------------------------------------------------------
          ! Snow mixing ratio (rsnowm)
          ! --------------------------------------------------------------------

          call microphys_lhs
     .         ( "rsnowm", .true., dt, nu_r, wmt, Vsnow, rsnowm, lhs )

          call microphys_solve
     .         ( "rsnowm", dt, lhs, rsnowm_tndcy, rsnowm, isValid )

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
     .         ( "rgraupelm", .true., dt, nu_r, wmt, Vgraupel, 
     .            rgraupelm, lhs )

          call microphys_solve
     .         ( "rgraupelm", dt, lhs, rgraupelm_tndcy,
     .           rgraupelm, isValid )

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
     .         ( "ricem", .true., dt, nu_r, wmt, Vice, ricem, lhs )

          call microphys_solve
     .         ( "ricem", dt, lhs, ricem_tndcy, ricem, isValid )

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

#ifdef STATS
        if ( lstats_samp ) then
          ! Rainfall rate (mm/day) should be defined on thermodynamic
          ! levels.  -Brian
          ! The absolute value of Vrr is taken because rainfall rate
          ! is a scalar quantity, and is therefore positive.
          if ( irain_rate > 0 ) then
            ! Rainfall rate (mm/day)
            zt%x(:,irain_rate) = zt%x(:,irain_rate) + 
     .        ( rrm * zm2zt( abs(Vrr) ) ) * ( rhot / rho_lw )
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

        end if ! lstats_samp
#endif /*STATS*/

        return
        end subroutine timestep_microphys
!-----------------------------------------------------------------------
        subroutine microphys_solve( solve_type, dt, lhs,
     .                              xrm_tndcy, xrm, isValid )

!       Description:

!       References:
!-----------------------------------------------------------------------
        use grid_class
        use lapack_wrap, only: tridag_solve !,band_solve
#ifdef STATS
        use hoc_stats
#endif

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real, intent(in) :: dt ! Timestep               [s]

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
        logical, intent(out) :: isValid

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
     .         rhs, xrm, isValid )

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
        xrm(1) = xrm(2)
        xrm(gr%nnzp) = xrm(gr%nnzp-1)

        return
        end subroutine microphys_solve

!-----------------------------------------------------------------------
        subroutine microphys_lhs
     .             ( solve_type, lsed, dt, nu, wmt, v_t, xrm, lhs )

!       Description:
!       Setup the matrix of implicit contributions to a term
!       Includes the effects of sedimentation, diffusion, and advection.
!
!       Notes:
!       Setup for tridiagonal system and boundary conditions should be
!       the same as the original rain subroutine code.
!-----------------------------------------------------------------------
        use grid_class
#ifdef STATS
        use hoc_stats
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

        real, intent(in) :: 
     .  dt, ! Timestep                  [s]
     .  nu  ! Diffusion coefficient     [-]

        real, intent(in), dimension(gr%nnzp) :: 
     .  xrm, ! Value of hydrometeor as prev. time       [units vary]
     .  wmt, ! w wind on t                              [m/s]
     .  v_t  ! Terminal velocity                        [m/s]

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
     .    + diffusion( nu, gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          ! Mean Advection
          lhs(kp1diag:km1diag,k)
     .    = lhs(kp1diag:km1diag,k)
     .    + advection( wmt(k), gr%dzt(k), k )
          ! Sedimentation
          ! Note: originally pristine ice did not sediment, so it was
          ! setup to be disabled as needed, but now pristine ice does
          ! sediment in the COAMPS case. -dschanen 12 Feb 2007
          if ( lsed ) then
             lhs(kp1diag:km1diag,k)
     .       = lhs(kp1diag:km1diag,k)
     .       + sedimentation( v_t(k), v_t(km1), gr%dzt(k), k )
          endif

         ! Implicit contributions to xrm
#ifdef STATS
          if ( lstats_samp ) then

            if ( ixrm_ma > 0 ) then
              tmp(1:3) = advection( wmt(k), gr%dzt(k), k )
              ztscr01(k) = -tmp(3)
              ztscr02(k) = -tmp(2)
              ztscr03(k) = -tmp(1)
            end if

            if ( ixrm_sd > 0 ) then
              tmp(1:3) = sedimentation( v_t(k), v_t(km1), gr%dzt(k), k )
              ztscr04(k) = -tmp(3)
              ztscr05(k) = -tmp(2)
              ztscr06(k) = -tmp(1)
            end if

            if ( ixrm_dff > 0 ) then
              tmp(1:3)
     .        = diffusion( nu, gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
              ztscr07(k) = -tmp(3)
              ztscr08(k) = -tmp(2)
              ztscr09(k) = -tmp(1)
            end if

          end if ! lstats_samp
#endif /*STATS*/

        end do ! 1..gr%nnzp

        return
        end subroutine microphys_lhs

!-----------------------------------------------------------------------
        pure function diffusion( nu, dzmm1, dzm, dzt, level ) 
     .  result( lhs ) 

!       Description:
!       Vertical diffusion of a hydrometeor using an implicit scheme.

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        ! Input Variables
        real, intent(in) :: nu ! Diffusion coefficient

        real, intent(in) :: 
     .  dzmm1, dzm, dzt ! Inverse of the grid spacing [m]

        integer, intent(in) :: level ! thermo. level that the calc. is for.

        ! Output
        real, dimension(3) :: lhs

        ! DESCRIPTION
        !
        ! This function handles the portion of the hydrometeor equation that
        ! is concerned with eddy diffusion:
        !
        ! d(xx)/dt = ... + nu * [ d^2(xx)/(dz)^2 ] + ...
        !
        ! The normal discretization is to use the values of the hydrometeor, xx,
        ! which is located on thermodynamic levels, from three levels -- the 
        ! central level and one-level-above and one-level-below the central 
        ! level.  These three values lead to two first derivatives centered at 
        ! the intermediate momentum levels.  In turn, these two first 
        ! derivatives lead to one second derivative centered at the central 
        ! thermodynamic level.  This is illustrated in the diagram below:
        !
        ! ---xx(k+1)----------------------------------  t(k+1)
        !
        ! ===============d(xx)/dz=====================  m(k)
        !
        ! ---xx( k )------------------d^2(xx)/(dz)^2--  t(k)
        !
        ! ===============d(xx)/dz=====================  m(k-1)
        !
        ! ---xx(k-1)----------------------------------  t(k-1)
        !
        ! In equation form, this looks like:
        !
        ! nu * dzt(k)*[dzm(k)*(xx(k+1)-xx(k)) - dzm(k-1)*(xx(k)-xx(k-1))];
        ! 
        ! where dzt(k) and dzm(k) are 1/(delta z) at level k.
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! However, there are special conditions at the boundaries, as any 
        ! quantity cannot be diffused outside of the boundaries, but rather
        ! builds up at the boundaries and then re-diffuses back into the rest of
        ! the domain.  We define a "phantom-point" just outside the boundary to 
        ! have the same value of the quantity as the point just inside the 
        ! boundary.  This results in a zero derivative at the point of the 
        ! boundary.  From a purely mathematical standpoint, any point that has a 
        ! first derivative of zero will also have a second derivative of zero.  
        ! Therefore, there will be no change to the quantity at the point of the 
        ! boundary due to diffusional effects.  This is illustrated in the 
        ! diagram of the lower boundary below:
        !
        ! --------------------xx(3)--------------------  t(3)
        !
        ! =============================================  m(2)
        !
        ! --------------------xx(2)--------------------  t(2)
        !
        ! ===SURFACE OR LOWER BOUNDARY (IMPERMEABLE)===  m(1)
        !
        ! ----"phantom point" xx(1) = xx(2)------------  t(1)
        !
        ! This is still discretized the same way as before, yielding:
        !
        ! nu * dzt(2)*[dzm(2)*(xx(3)-xx(2)) - dzm(1)*(xx(2)-xx(1))],
        !
        ! but since xx(1) = xx(2):
        !
        ! nu * dzt(2)*[dzm(2)*(xx(3)-xx(2))]
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! The upper boundary is treated in the exact same way as the lower
        ! boundary.  In HOC, the first thermodynamic level is below the surface,
        ! so there are no diffusional effects at level 1.  The boundary levels
        ! are level 2 and the top level.  Everything in between operates as 
        ! normal.
        !
        ! Finally, it should be pointed out that this technique of discretizing
        ! the equations leads to conservative differencing (for the diffusional
        ! term), as long as the grid levels are equally spaced.  When 
        ! conservative differencing is in place, the column totals in the LHS
        ! matrix (for the diffusional term) should be equal to 0.  This insures
        ! that the total amount of the quantity over the entire domain is being
        ! conserved, meaning that nothing is lost due to diffusional effects.
        ! This is illustrated in the matrix setup below:
        !
        ! Left-hand side matrix, diffusional term only; first five levels:
        !
        !   ----->
        ! 0 | 0            0                     0                              0
        !   |
        !   | 0   +nu*dzt(k)*dzm(k)     -nu*dzt(k)*dzm(k)                       0
        !   |
        !   | 0   -nu*dzt(k)*dzm(k-1)   +nu*dzt(k)*(dzm(k)+dzm(k-1))   -nu*dzt(k)*dzm(k)
        !   |
        !   | 0            0            -nu*dzt(k)*dzm(k-1)            +nu*dzt(k)*(dzm(k)+dzm(k-1))
        !   |
        !   | 0            0                     0                     -nu*dzt(k)*dzm(k-1)
        !   |
        !  \ /
        !
        ! Brian Griffin.  February 27, 2007.

        if ( level == 1 ) then
           ! k = 1; below surface; no effects.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = 0.0
           lhs(km1diag) = 0.0

        elseif ( level == 2 ) then
           ! k = 2; lower boundary level.
           lhs(kp1diag) = - nu * dzt * dzm
           lhs(kdiag)   = + nu * dzt * dzm
           lhs(km1diag) = 0.0

        elseif ( level > 2 .and. level < gr%nnzp ) then
           ! Most of the interior model; normal conditions.
           lhs(kp1diag) = - nu * dzt * dzm
           lhs(kdiag)   = + nu * dzt * ( dzm + dzmm1 )
           lhs(km1diag) = - nu * dzt * dzmm1

        elseif ( level == gr%nnzp ) then
           ! k = gr%nnzp (top level); upper boundary level.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = + nu * dzt * dzmm1
           lhs(km1diag) = - nu * dzt * dzmm1

        endif

        return
        end function diffusion

!-----------------------------------------------------------------------
        pure function advection( wmt, dzt, level )
     .  result( lhs )

!       Description:
!       Mean advection term

!       References:
!       None
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        real, intent(in) ::
     .  wmt, ! w wind on t                      [m/s]
     .  dzt  ! Inverse of the grid spacing      [m] 

        integer, intent(in) :: level ! thermo. level that the calc. is for.

        ! Output
        real, dimension(3) :: lhs

        ! DESCRIPTION
        !
        ! This function handles the portion of the hydrometeor equation that
        ! is concerned with mean advection in the vertical:
        !
        ! d(xx)/dt = ... -w * d(xx)/dz  + ...
        !
        ! The normal discretization is to use the values of the hydrometeor, xx,
        ! which is located on thermodynamic levels, from three levels -- the 
        ! central level and one-level-above and one-level-below the central 
        ! level.  These three values lead to two averaged values centered at 
        ! the intermediate momentum levels.  In turn, these two averaged 
        ! values lead to one first derivative centered at the central 
        ! thermodynamic level.  This is illustrated in the diagram below:
        !
        ! ---xx(k+1)----------------------------------  t(k+1)
        !
        ! ===============ave(xx)m=====================  m(k)
        !
        ! ---xx( k )------------------d(xx)/dz--------  t(k)
        !
        ! ===============ave(xx)m=====================  m(k-1)
        !
        ! ---xx(k-1)----------------------------------  t(k-1)
        !
        ! The value for vertical wind velocity, w, is taken from the central
        ! thermodynamic level (w is normally a momentum-level variable, but
        ! HOC also interpolates its value to thermodynamic levels).
        !
        ! In equation form, this looks like:
        !
        ! -w_t(k) * [dzt(k)*( (xx(k+1)+xx(k))/2 - (xx(k)+xx(k-1))/2 )];
        !
        ! which is reduced to:
        ! 
        ! -w_t(k) * [(1/2)*dzt(k)*(xx(k+1)-xx(k-1))];
        ! 
        ! where dzt(k) and dzm(k) are 1/(delta z) at level k.
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! However, there are special conditions at the boundaries.  We need 
        ! three points to take a centered derivative, but a boundary yields only
        ! two points.  We define a "phantom-point" just outside the boundary to 
        ! have the same value of the quantity as the point just inside the 
        ! boundary.  This results in a zero derivative at the point of the 
        ! boundary, although the level just above the boundary, which is where 
        ! the effects are taking place, will have a derivative.  This is 
        ! illustrated in the diagram of the lower boundary below:
        !
        ! --------------------xx(3)--------------------  t(3)
        !
        ! =============================================  m(2)
        !
        ! --------------------xx(2)--------------------  t(2)
        !
        ! ===SURFACE OR LOWER BOUNDARY (IMPERMEABLE)===  m(1)
        !
        ! ----"phantom point" xx(1) = xx(2)------------  t(1)
        !
        ! This is still discretized the same way as before, yielding:
        !
        ! -w_t(2) * [dzt(2)*( (xx(3)+xx(2))/2 - (xx(2)+xx(1))/2 )];
        !
        ! but since xx(1) = xx(2):
        !
        ! -w_t(2) * [(1/2)*dzt(2)*(xx(3)-xx(2))];
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! The upper boundary is treated in the exact same way as the lower
        ! boundary.  In HOC, the first thermodynamic level is below the surface,
        ! so there are no advective effects at level 1.  The boundary levels
        ! are level 2 and the top level.  Everything in between operates as 
        ! normal.
        !
        ! Finally, it should be pointed out that this technique of discretizing
        ! is not conservative.  Column totals at any level do not necessarily
        ! add up to 0, as w may be different at any level.
        !
        ! Brian Griffin.  March 7, 2007.

        if ( level == 1 ) then
           ! k = 1; below surface; no effects.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = 0.0
           lhs(km1diag) = 0.0

        elseif ( level == 2 ) then
           ! k = 2; lower boundary level.
           lhs(kp1diag) = + wmt * (1.0/2.0) * dzt
           lhs(kdiag)   = - wmt * (1.0/2.0) * dzt
           lhs(km1diag) = 0.0

        elseif ( level > 2 .and. level < gr%nnzp ) then
           ! Most of the interior model; normal conditions.
           lhs(kp1diag) = + wmt * (1.0/2.0) * dzt
           lhs(kdiag)   = 0.0
           lhs(km1diag) = - wmt * (1.0/2.0) * dzt

        elseif ( level == gr%nnzp ) then
           ! k = gr%nnzp (top level); upper boundary level.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = + wmt * (1.0/2.0) * dzt
           lhs(km1diag) = - wmt * (1.0/2.0) * dzt

        endif

        return
        end function advection

!-----------------------------------------------------------------------
        pure function sedimentation( v_t, v_tm1, dzt, level )
     .  result( lhs ) 

!       Description:
!       Sediment a hydrometer using v_t, the calculated terminal 
!       velocity.

!       References:
!       None

!       Notes:
!       Both COAMPS Microphysics and Brian Griffin's implementation use
!       Khairoutdinov and Kogan (2000) for the calculation of rain mixing
!       ratio and rain droplet number concentration sedimentation
!       velocities.
!-----------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  kdiag   = 2,
     .  km1diag = 3

        real, intent(in) ::
     .  v_t,   ! Terminal velocity <k>          [m/s]
     .  v_tm1, ! Terminal velocity <k-1>        [m/s]
     .  dzt    ! Inverse of the grid spacing <k>[m] 

        integer, intent(in) :: level ! thermo. level that the calc. is for.

        ! Output
        real, dimension(3) :: lhs

        ! DESCRIPTION
        !
        ! This function handles the portion of the hydrometeor equation that
        ! is concerned with drop or particle sedimentation:
        !
        ! d(xx)/dt = ... -Vxx * d(xx)/dz  + ...
        !
        ! The normal discretization is to use the values of the hydrometeor, xx,
        ! which is located on thermodynamic levels, from three levels -- the 
        ! central level and one-level-above and one-level-below the central 
        ! level.  These three values lead to two averaged values centered at 
        ! the intermediate momentum levels.  In turn, these two averaged 
        ! values lead to one first derivative centered at the central 
        ! thermodynamic level.  Meanwhile, the values for sedimentation 
        ! velocity, Vxx, are located at the intermediate momentum levels.  
        ! The two values lead to one averaged value centered at the central 
        ! thermodynamic level.  This is illustrated in the diagram below:
        !
        ! ---xx(k+1)------------------------------------  t(k+1)
        !
        ! ===Vxx( k )=====ave(xx)m======================  m(k)
        !
        ! ---xx( k )------ave(Vxx)t-----d(xx)/dz--------  t(k)
        !
        ! ===Vxx(k-1)=====ave(xx)m======================  m(k-1)
        !
        ! ---xx(k-1)------------------------------------  t(k-1)
        !
        ! In equation form, this looks like:
        !
        ! -( (Vxx(k)+Vxx(k-1))/2 ) * 
        !      [dzt(k)*( (xx(k+1)+xx(k))/2 - (xx(k)+xx(k-1))/2 )];
        !
        ! which is reduced to:
        ! 
        ! -( (Vxx(k)+Vxx(k-1))/2 ) * [(1/2)*dzt(k)*(xx(k+1)-xx(k-1))];
        !
        ! and then rearranged to:
        !
        ! -(1/4) * (Vxx(k)+Vxx(k-1)) * dzt(k) * (xx(k+1)-xx(k-1));
        ! 
        ! where dzt(k) and dzm(k) are 1/(delta z) at level k.
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! However, there are special conditions at the boundaries.  We need 
        ! three points to take a centered derivative, but a boundary yields only
        ! two points.  We define a "phantom-point" just outside the boundary to 
        ! have the same value of the quantity as the point just inside the 
        ! boundary.  This is illustrated in the diagram of the lower boundary 
        ! below:
        !
        ! --------------------xx(3)--------------------  t(3)
        !
        ! =============================================  m(2)
        !
        ! --------------------xx(2)--------------------  t(2)
        !
        ! ===SURFACE OR LOWER BOUNDARY (IMPERMEABLE)===  m(1)
        !
        ! ----"phantom point" xx(1) = xx(2)------------  t(1)
        !
        ! This is still discretized the same way as before, yielding:
        !
        ! -( (Vxx(2)+Vxx(1))/2 ) * 
        !      [dzt(2)*( (xx(3)+xx(2))/2 - (xx(2)+xx(1))/2 )];
        !
        ! but since xx(1) = xx(2):
        !
        ! -(1/4) * (Vxx(2)+Vxx(1)) * dzt(2) * (xx(3)-xx(2));
        !
        ! This is then moved to the left-hand side of the equation and 
        ! rearranged to be placed onto a tridiagonal matrix.
        !
        ! The upper boundary is treated in the exact same way as the lower
        ! boundary.  In HOC, the first thermodynamic level is below the surface,
        ! so there are no sedimentary effects at level 1.  The boundary levels
        ! are level 2 and the top level.  Everything in between operates as 
        ! normal.
        !
        ! Finally, it should be pointed out that this technique of discretizing
        ! is not conservative.  Column totals at any level do not necessarily
        ! add up to 0, as Vxx may be different at any level.  We want this, as
        ! we should be losing some of the hydrometeor through the lower boundary
        ! as it sediments.
        !
        ! Brian Griffin. March 7, 2007.

        if ( level == 1 ) then
           ! k = 1; below surface; no effects.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = 0.0
           lhs(km1diag) = 0.0

        elseif ( level == 2 ) then
           ! k = 2; lower boundary level.
           lhs(kp1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           lhs(kdiag)   = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           lhs(km1diag) = 0.0
           ! For a positive downward sedimentation velocity
!           lhs(kp1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
!           lhs(kdiag)   = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt
!           lhs(km1diag) = 0.0

        elseif ( level > 2 .and. level < gr%nnzp ) then
           ! Most of the interior model; normal conditions.
           lhs(kp1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           lhs(kdiag)   = 0.0
           lhs(km1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           ! For a positive downward sedimentation velocity
!           lhs(kp1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
!           lhs(kdiag)   = 0.0
!           lhs(km1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt

        elseif ( level == gr%nnzp ) then
           ! k = gr%nnzp (top level); upper boundary level.
           lhs(kp1diag) = 0.0
           lhs(kdiag)   = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           lhs(km1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
           ! For a positive downward sedimentation velocity
!           lhs(kp1diag) = 0.0
!           lhs(kdiag)   = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
!           lhs(km1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt

        endif

        return
        end function sedimentation

        end module microphysics
