!-----------------------------------------------------------------------
! $Id: microphysics.F,v 1.4 2007-02-21 19:59:21 mjfalk Exp $
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
        use constants, only: nu_r, hydromet_dim, Lv, Cp, rho_lw
#ifdef STATS
        use statistics
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
     .  lice = .true. ! Whether to produce ice in COAMPS

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
          ! Note that Ncm for DYCOMS II RF 02 is specified in the GCSS,
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

        ! Rain water mixing ratio
        call microphys_lhs
     .       ( "rrm", .true., dt, nu_r, wmt, Vrr, rrm, lhs )

        call microphys_solve
     .       ( "rrm", dt, lhs, rrm_tndcy, rrm, isValid )

        ! Rain droplet number concentration
        call microphys_lhs
     .       ( "Nrm", .true., dt, nu_r, wmt, VNr, Nrm, lhs )

        call microphys_solve
     .       ( "Nrm", dt, lhs, Nrm_tndcy, Nrm, isValid )


        if ( lcoamps_micro ) then

          ! Snow mixing ratio
          call microphys_lhs
     .         ( "rsnowm", .true., dt, nu_r, wmt, Vsnow, rsnowm, lhs )

          call microphys_solve
     .         ( "rsnowm", dt, lhs, rsnowm_tndcy, rsnowm, isValid )

          ! Graupel mixing ratio
          call microphys_lhs
     .         ( "rgraupelm", .true., dt, nu_r, wmt, Vgraupel, 
     .            rgraupelm, lhs )

          call microphys_solve
     .         ( "rgraupelm", dt, lhs, rgraupelm_tndcy,
     .           rgraupelm, isValid )

          ! Pristine ice mixing ratio
          ! Note: Ice sediments in new COAMPS
          call microphys_lhs
     .         ( "ricem", .true., dt, nu_r, wmt, Vice, ricem, lhs )

          call microphys_solve
     .         ( "ricem", dt, lhs, ricem_tndcy, ricem, isValid )

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
        use statistics
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
     .  ixrm_bt, ! Change in xrm over time stats toggle
     .  ixrm_ma, ! Mean advection budget stats toggle
     .  ixrm_sd, ! Sedimentation budget stats toggle
     .  ixrm_dff,! Diffusion budget stats toggle
     .  ixrm_cl  ! Clipping term stats toggle

        select case( solve_type )
        case( "rrm" )
          ixrm_bt  = irrm_bt 
          ixrm_ma  = irrm_ma
          ixrm_sd  = irrm_sd
          ixrm_dff = irrm_dff
          ixrm_cl  = irrm_cl
        case( "Nrm" )
          ixrm_bt  = iNrm_bt 
          ixrm_ma  = iNrm_ma
          ixrm_sd  = iNrm_sd
          ixrm_dff = iNrm_dff
          ixrm_cl  = iNrm_cl
        case( "ricem" )
          ixrm_bt  = iricem_bt 
          ixrm_ma  = iricem_ma
          ixrm_sd  = iricem_sd
          ixrm_dff = iricem_dff
          ixrm_cl  = iricem_cl
        case( "rsnowm" )
          ixrm_bt  = irsnowm_bt 
          ixrm_ma  = irsnowm_ma
          ixrm_sd  = irsnowm_sd
          ixrm_dff = irsnowm_dff
          ixrm_cl  = irsnowm_cl
        case( "rgraupelm" )
          ixrm_bt  = irgraupelm_bt 
          ixrm_ma  = irgraupelm_ma
          ixrm_sd  = irgraupelm_sd
          ixrm_dff = irgraupelm_dff
          ixrm_cl  = irgraupelm_cl
        case default
          ixrm_bt  = 0
          ixrm_ma  = 0
          ixrm_sd  = 0
          ixrm_dff = 0
          ixrm_cl  = 0
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
     .         + ztscr01(k) * xrm(km1)
     .         + ztscr02(k) * xrm(kp1)
               zt%n(k,ixrm_ma) = zt%n(k,ixrm_ma) + 1
             end if

             if ( ixrm_sd > 0 ) then
               zt%x(k,ixrm_sd) = zt%x(k,ixrm_sd)
     .           + ztscr03(k) * xrm(km1)
     .           + ztscr04(k) * xrm(kp1)
               zt%n(k,ixrm_sd) = zt%n(k,ixrm_sd) + 1
             end if

             if ( ixrm_dff > 0 ) then
               zt%x(k,ixrm_dff) = zt%x(k,ixrm_dff)
     .           + ztscr05(k) * xrm(km1)
     .           + ztscr06(k) * xrm(k)
     .           + ztscr07(k) * xrm(kp1)
               zt%n(k,ixrm_dff) = zt%n(k,ixrm_dff) + 1
             end if

           end do ! 1..gr%nnzp
        end if ! lstats_samp

#endif /*STATS*/

#ifdef STATS
        ! Store computed value before clipping
        if (lstats_samp .and. ixrm_cl > 0 ) then
          zt%x(:,ixrm_cl) = zt%x(:,ixrm_cl) - xrm / dt
        end if
#endif
        ! Value of mixing ratios / number concentrations cannot fall
        ! below 0, so we clip accordingly
        where ( xrm < 0.0 ) xrm = 0.0

#ifdef STATS
        if (lstats_samp .and. ixrm_cl > 0 ) then
          zt%x(:,ixrm_cl) = zt%x(:,ixrm_cl) + xrm / dt
          zt%n(:,ixrm_cl) = zt%n(:,ixrm_cl) + 1
        end if
#endif

        ! Boundary conditions on results
        xrm(1) = xrm(2)
        xrm(gr%nnzp) = xrm(gr%nnzp-1)

#ifdef STATS
        ! Finalize time tendency
        if ( lstats_samp .and. ixrm_bt > 0 ) then
          zt%x(:,ixrm_bt) = zt%x(:,ixrm_bt) + ( xrm / dt )
          zt%n(:,ixrm_bt) = zt%n(:,ixrm_bt) + 1
        end if
#endif

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
        use statistics
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

        ! Setup LHS Matrix
        do k = 2, gr%nnzp-1, 1
          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          ! All diagonals
          lhs(kp1diag:km1diag,k)
     .    = diffusion( nu, gr%dzm(km1), gr%dzm(k), gr%dzt(k) )

          ! Main diagonal
          lhs(kdiag,k) = lhs(kdiag,k) + ( 1.0 / dt ) ! Time tendency

          ! Super and sub diagonal
          lhs((/kp1diag,km1diag/),k) 
     .    = lhs((/kp1diag,km1diag/),k) 
     .    + advection( wmt(k), gr%dzt(k) )

          ! Super and sub diagonal
          ! Note: originally pristine ice did not sediment, so it was
          ! setup to be disabled as needed, but now pristine ice does
          ! sediment in the COAMPS case. -dschanen 12 Feb 2007
          if ( lsed ) then
            lhs((/kp1diag,km1diag/),k) 
     .      = lhs((/kp1diag,km1diag/),k) 
     .      + sedimentation( v_t(k), v_t(km1), gr%dzt(k) )
          end if

         ! Implicit contributions to xrm
#ifdef STATS
          if ( lstats_samp ) then

            if ( ixrm_ma > 0 ) then
              tmp(1:2) = advection( wmt(k), gr%dzt(k) )
              ztscr01(k) = -tmp(2)
              ztscr02(k) = -tmp(1)
            end if

            if ( ixrm_sd > 0 ) then
              tmp(1:2) = sedimentation( v_t(k), v_t(km1), gr%dzt(k) )
              ztscr03(k) = -tmp(2)
              ztscr04(k) = -tmp(1)
            end if

            if ( ixrm_dff > 0 ) then
              tmp(1:3)
     .        = diffusion( nu, gr%dzm(km1), gr%dzm(k), gr%dzt(k) )
              ztscr05(k) = -tmp(3)
              ztscr06(k) = -tmp(2)
              ztscr07(k) = -tmp(1)
            end if

          end if ! lstats_samp
#endif /*STATS*/

        end do ! 2..gr%nnzp-1

        ! Boundary conditions

        ! These have been chosen so that the column totals in the LHS matrix
        ! of the diffusional term are equal to zero.  This leads to 
        ! "conservative differencing" of the diffusion term for equally spaced grids.     
        ! The ground is impermeable to diffusion of rrm, Nrm.   Brian.

        ! Lower boundary conditions

        ! Not referenced by solver
        !lhs(km1diag,1) = 0.0   ! sub(1) term

        lhs(kdiag,1)
     .  =  + (1.0/dt)                      ! main(1) time tendency
     .  - ( - nu * gr%dzt(2) * gr%dzm(1) ) ! minus sub(2) diffusion

        lhs(kp1diag,1)
     .  = + wmt(1) * (1.0/2.0) * gr%dzt(1)                 ! super(1) advection
     .  + (1.0/4.0) * ( v_t(1) + v_t(1) ) * gr%dzt(1)      ! super(1) sediment.
     .  - ( + nu * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) ) ) ! minus main(2) diffusion
     .  - ( - nu * gr%dzt(3) * gr%dzm(2) )                 ! minus sub(3) diffusion

        ! Upper boundary conditions
        lhs(km1diag,gr%nnzp)
     .  = - ( nu*gr%dzt(gr%nnzp-1)
     .        *( gr%dzm(gr%nnzp-1) + gr%dzm(gr%nnzp-2) ) ) ! main(gr%nnzp-1) diffusion
     .    - ( - nu*gr%dzt(gr%nnzp-2)*gr%dzm(gr%nnzp-2) )   ! super(gr%nnzp-2) diffusion
     .    - wmt(gr%nnzp) * (1.0/2.0) * gr%dzt(gr%nnzp)     ! sub(gr%nnzp) advection
                                                           ! sub(gr%nnzp) sediment = 0.

        lhs(kdiag,gr%nnzp)
     .  = - ( - nu*gr%dzt(gr%nnzp-1)
     .      *gr%dzm(gr%nnzp-1) ) ! super(gr%nnzp-1) diffusion
     .    + (1.0/dt) ! main(gr%nnzp) time tendency

        ! Not referenced
        !lhs(kp1diag,gr%nnzp) = 0.0 ! super(gr%nnzp) term

        ! Implicit contributions (as in subroutine rain)
#ifdef STATS

        if ( lstats_samp ) then

         ! Lower boundary
         if ( ixrm_ma > 0 ) then
           ztscr01(1) = 0.0
           ztscr02(1) = - wmt(1) * (1.0/2.0) * gr%dzt(1)
         endif

         if ( ixrm_sd > 0 ) then
           ztscr03(1) = 0.0
           ztscr04(1) = - (1.0/4.0) * ( v_t(1) + v_t(1) ) * gr%dzt(1)
         endif

         if ( ixrm_dff > 0 ) then
           ztscr05(1) = 0.0
           ztscr06(1) = - nu * gr%dzt(2) * gr%dzm(1)
           ztscr07(1) = + nu * gr%dzt(2) * ( gr%dzm(2) + gr%dzm(1) )
     .                  - nu * gr%dzt(3) * gr%dzm(2)
         endif

         ! Upper boundary
         if ( ixrm_ma > 0 ) then
           ztscr01(gr%nnzp) = + wmt(gr%nnzp) 
     .                          * (1.0/2.0) * gr%dzt(gr%nnzp)
           ztscr02(gr%nnzp) = 0.0
         end if

         if ( ixrm_sd > 0 ) then
           ztscr03(gr%nnzp) = 0.0
           ztscr04(gr%nnzp) = 0.0
         end if

         if ( ixrm_dff > 0 ) then
           ztscr05(gr%nnzp) = 
     .         + (   nu*gr%dzt(gr%nnzp-1)
     .              *(gr%dzm(gr%nnzp-1) + gr%dzm(gr%nnzp-2)) )
     .         + ( - nu*gr%dzt(gr%nnzp-2)*gr%dzm(gr%nnzp-2) )

           ztscr06(gr%nnzp) 
     .     = + ( - nu*gr%dzt(gr%nnzp-1)*gr%dzm(gr%nnzp-1) )

           ztscr07(gr%nnzp) = 0.0 
         end if

        end if ! lstats_samp

#endif /*STATS*/

        return
        end subroutine microphys_lhs

!-----------------------------------------------------------------------
        pure function diffusion( nu, dzmm1, dzm, dzt ) 
     .  result( lhs ) 

!       Description:
!       Vertical diffusion of a hydrometeor using an implicit scheme.

!       References:
!       None
!-----------------------------------------------------------------------

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

        ! Output
        real, dimension(3) :: lhs

        lhs(kp1diag) = - nu * dzt * dzm
        lhs(kdiag)   = + nu * dzt * ( dzm + dzmm1 )
        lhs(km1diag) = - nu * dzt * dzmm1

        return
        end function diffusion

!-----------------------------------------------------------------------
        pure function advection( wmt, dzt )
     .  result( lhs )

!       Description:
!       Mean advection term

!       References:
!       None
!-----------------------------------------------------------------------

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  km1diag = 2

        real, intent(in) ::
     .  wmt, ! w wind on t                      [m/s]
     .  dzt  ! Inverse of the grid spacing      [m] 

        ! Output
        real, dimension(2) :: lhs

        lhs(kp1diag) = + wmt * (1.0/2.0) * dzt
        lhs(km1diag) = - wmt * (1.0/2.0) * dzt

        return
        end function advection

!-----------------------------------------------------------------------
        pure function sedimentation( v_t, v_tm1, dzt )
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

        implicit none

        ! Constant parameters
        integer, parameter ::
     .  kp1diag = 1,
     .  km1diag = 2

        real, intent(in) ::
     .  v_t,   ! Terminal velocity <k>          [m/s]
     .  v_tm1, ! Terminal velocity <k-1>        [m/s]
     .  dzt    ! Inverse of the grid spacing <k>[m] 

        ! Output
        real, dimension(2) :: lhs

        ! For a positive downward sedimentation velocity
!       lhs(kp1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt
!       lhs(km1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt

        ! For a negative downward sedimentation velocity
        lhs(kp1diag) = + (1.0/4.0) * ( v_t + v_tm1 ) * dzt
        lhs(km1diag) = - (1.0/4.0) * ( v_t + v_tm1 ) * dzt

        return
        end function sedimentation

        end module microphysics
