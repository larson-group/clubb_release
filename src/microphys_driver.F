!-----------------------------------------------------------------------
! $Id: microphys_driver.F,v 1.23 2008-07-21 20:16:18 faschinj Exp $
        module microphys_driver

!       Description:
!       Call a microphysical scheme to compute hydrometeor species,
!       and advect, sediment, & diffuse using a tridiagonal system.

!       References:
!       None
!-----------------------------------------------------------------------
        implicit none

        public :: advance_microphys, init_microphys

        ! Subroutines
        private :: microphys_lhs, microphys_solve
        private :: adj_microphys_tndcy

        ! Functions
        private :: sedimentation

        private ! Default Scope

        character(len=10), dimension(:), allocatable ::
     .  hydromet_list

        contains

!-----------------------------------------------------------------------
        subroutine init_microphys
     .             ( kk_rain, lcoamps_micro, licedfs, hydromet_dim )

!       Description:
!       Set indices to the various hydrometeor species and define
!       hydromet_dim for the purposes of allocating memory.

!       References:
!       None
!-----------------------------------------------------------------------
        use array_index, only: 
     .      iirrm, iiNrm, iirsnowm, iiricem, iirgraupelm

        implicit none

        logical, intent(in) ::
     .  kk_rain,        ! K&K microphysics
     .  lcoamps_micro,  ! COAMPS microphysics
     .  licedfs         ! Simplified ice

        integer, intent(out) ::
     .  hydromet_dim ! Number of hydrometeor fields.

        ! The location of the fields in the hydromet array are arbitrary,
        ! and don't need to be set consistently among schemes so long as
        ! the 'i' indices point to the correct parts of the array.
        if ( lcoamps_micro ) then
          iirrm       = 1
          iiNrm       = 2
          iirsnowm    = 3
          iiricem     = 4
          iirgraupelm = 5

          hydromet_dim = 5

          allocate( hydromet_list(hydromet_dim) )

          hydromet_list(iirrm)       = "rrm"
          hydromet_list(iiNrm)       = "Nrm"
          hydromet_list(iirsnowm)    = "rsnowm"
          hydromet_list(iiricem)     = "ricem"
          hydromet_list(iirgraupelm) = "rgraupelm"

        else if ( kk_rain ) then
          iirrm       = 1
          iiNrm       = 2
          iirsnowm    = -1
          iiricem     = -1
          iirgraupelm = -1

          hydromet_dim = 2

          allocate( hydromet_list(hydromet_dim) )

          hydromet_list(iirrm) = "rrm"
          hydromet_list(iiNrm) = "Nrm"

        else if ( licedfs ) then
          iirrm       = -1
          iiNrm       = -1
          iirsnowm    = -1
          iiricem     = -1
          iirgraupelm = -1

          hydromet_dim = 0

        else
          iirrm       = -1
          iiNrm       = -1
          iirsnowm    = -1
          iiricem     = -1
          iirgraupelm = -1

          hydromet_dim = 0

        end if

        return
        end subroutine init_microphys

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
     .      c_Krrm, ! Variable(s) 
     .      hydromet_dim, 
     .      nu_r

        use constants, only: 
     .      Lv,  ! Constant(s)
     .      Cp, 
     .      rho_lw, 
     .      fstderr,
     .      sec_per_day

        use stats_precision, only: 
     .      time_precision ! Variable(s)

        use error_code, only: 
     .      lapack_error ! Procedure

        use coamps_micro_driver_mod, only: 
     .      coamps_micro_driver ! Procedure

        use T_in_K_mod, only: thlm2T_in_K ! Procedure(s)

        use array_index, only: 
     .      iirrm, iiNrm, iirsnowm, iiricem, iirgraupelm


#ifdef STATS
        use stats_variables, only:
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

        use stats_type, only: 
     .      stat_update_var, stat_update_var_pt, ! Procedure(s)
     .      stat_begin_update, stat_end_update 

#endif
        implicit none

        ! Constant Parameters
        logical, parameter :: lsed = .true.

        character(len=*), intent(in) ::
     .  runtype ! Name of the run, for case specific effects.

        real(kind=time_precision), intent(in) :: 
     .  dt           ! Timestep         [s]

        real(kind=time_precision), intent(in) :: 
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

        real, target,dimension(gr%nnzp,26),intent(in) ::
     .  pdf_parms     ! PDF parameters

        real, dimension(gr%nnzp,hydromet_dim), intent(inout) ::
     .  hydromet      ! Array of rain, prist. ice, graupel, etc. [units vary]

        real, dimension(gr%nnzp), intent(inout) ::
     .  rtm_forcing, ! Imposed contributions to total water        [kg/kg/s]
     .  thlm_forcing ! Imposed contributions to liquid potential temp. [K/s]

        integer, intent(out) :: err_code ! Exit code returned from subroutine

        ! Local Variables
        real, dimension(3,gr%nnzp) ::
     .  lhs ! Left hand side of tridiagonal matrix

        real, dimension(gr%nnzp,hydromet_dim) ::
     .  hydromet_vel ! Contains vel. of the hydrometeors        [m/s]
        !  Can contain:
        !  Vrr      ! Rain mixing ratio sedimentation velocity     [m/s]
        !  VNr      ! Rain number conc. sedimentation velocity     [m/s]
        !  Vice     ! Ice mixing ratio sedimentation velocity      [m/s]
        !  Vsnow    ! Snow mixing ratio sedimentation velocity     [m/s]
        !  Vgraupel ! Graupel mixing ratio sedimentation velocity  [m/s]

        real, dimension(gr%nnzp) ::
     .  rtm_mc,  ! Change in total water due to microphysics    [(kg/kg)/s]
     .  thlm_mc  ! Change in liquid potential temperature 
                 ! due to microphysics                          [K/s]

        real, dimension(gr%nnzp,hydromet_dim) ::
     .  hydromet_mc
     .  ! Rain mixing ratio tendency      [(kg/kg)/s]
     .  ! Rain number conc. tendency      [(count/kg)/s]
     .  ! Snow mixing ratio tendency      [(kg/kg)/s]
     .  ! Ice mixing ratio tendency       [(kg/kg)/s]
     .  ! Graupel mixing ratio tendency   [(kg/kg)/s]  

        real, dimension(gr%nnzp) ::
     .  T_in_K  ! Temperature   [K]

        real, dimension(1,1,gr%nnzp) ::
     .  cond ! COAMPS stat for condesation/evap of rcm

!       real, pointer, dimension(:) ::
!    .  rrm,      ! Pointer for rain water mixing ratio   [kg/kg]
!    .  Nrm,      ! Pointer for rain droplet number conc. [count/kg]
!    .  rsnowm,   ! Pointer for snow mixing ratio         [kg/kg]
!    .  ricem,    ! Pointer for ice mixing ratio          [kg/kg]
!    .  rgraupelm ! Pointer for graupel mixing ratio      [kg/kg]

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

        integer :: i, k ! Array index
#ifdef STATS
        integer :: ixrm_cl, ixrm_bt
#endif

!-----------------------------------------------------------------------

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
        T_in_K = thlm2T_in_K( thlm, exner, rcm )

        ! Begin by calling either Brian Griffin's implementation of the
        ! Khairoutdinov and Kogan microphysical scheme or 
        ! alternatively the Rutlege and Hobbes scheme from COAMPS(R).
        ! Note: COAMPS appears to have some K&K elements to it as well.

        if ( lcoamps_micro ) then

           call coamps_micro_driver
     .          ( runtype, time_current, dt,
     .            rtm, wmm, p, exner, rhot, T_in_K,
     .            thlm, hydromet(:,iiricem), hydromet(:,iirrm), 
     .            hydromet(:,iirgraupelm), hydromet(:,iirsnowm),
     .            rcm, Ncm, hydromet(:,iiNrm), Ncnm, Nim, cond,
     .            hydromet_vel(:,iirsnowm), hydromet_vel(:,iiricem),
     .            hydromet_vel(:,iirrm), hydromet_vel(:,iiNrm), 
     .            hydromet_vel(:,iirgraupelm),
     .            hydromet_mc(:,iiricem), hydromet_mc(:,iirrm),
     .            hydromet_mc(:,iirgraupelm), hydromet_mc(:,iirsnowm),
     .            hydromet_mc(:,iiNrm),
     .            rtm_mc, thlm_mc )

#ifdef STATS
          if ( lstats_samp ) then
              
           ! Sedimentation velocity for rrm
           call stat_update_var(iVrr, hydromet_vel(:,iirrm), zm)

           ! Sedimentation velocity for Nrm
           call stat_update_var(iVNr, hydromet_vel(:,iiNrm), zm )

           ! Sedimentation velocity for snow
           call stat_update_var(iVsnow, hydromet_vel(:,iirsnowm), zm )

           ! Sedimentation velocity for pristine ice
           call stat_update_var( iVice, hydromet_vel(:,iiricem), zm )

           ! Sedimentation velocity for graupel
           call stat_update_var( iVgraupel, 
     .                           hydromet_vel(:,iirgraupelm), zm )

           ! Sum total of rrm microphysics
           call stat_update_var( irrm_mc, hydromet_mc(:,iirrm), zt )

           ! Sum total of Nrm microphysics
           call stat_update_var( iNrm_mc, hydromet_mc(:,iiNrm), zt )
           
           ! Sum total of pristine ice microphysics
           call stat_update_var( iricem_mc, hydromet_mc(:,iiricem), zt )

           ! Sum total of graupel microphysics
           call stat_update_var( irgraupelm_mc,
     .                           hydromet_mc(:,iirgraupelm), zt )

           ! Sum total of snow microphysical processeses
           call stat_update_var( irsnowm_mc, 
     .                           hydromet_mc(:,iirsnowm), zt )

         end if ! lstats_samp
#endif /*STATS*/
        else if ( kk_rain ) then

          ! Note that Ncm for DYCOMS II RF02 is specified in the
          ! dycoms_rf02_tndcy subroutine, so new cases will need 
          ! Ncm computed beforehand as well.

          call kk_microphys
     .         ( thlm, T_in_K, p, exner, rhot, 
     .           thl1, thl2, a, rc1, rc2, s1, 
     .           s2, ss1, ss2, rtm, rcm, Ncm, 
     .           hydromet(:,iirrm), hydromet(:,iiNrm),
     .           .true., AKm_est, AKm, 
     .           hydromet_mc(:,iirrm), hydromet_mc(:,iiNrm), 
     .           rtm_mc, thlm_mc,
     .           hydromet_vel(:,iirrm), hydromet_vel(:,iiNrm) )

#ifdef STATS
          if ( lstats_samp ) then
              
           ! Sedimentation velocity for rrm
           call stat_update_var( iVrr, hydromet_vel(:,iirrm), zm )

           ! Sedimentation velocity for Nrm
           call stat_update_var( iVNr, hydromet_vel(:,iiNrm), zm )

           ! Sum total of rrm microphysics (auto + accr + cond)
           call stat_update_var( irrm_mc, hydromet_mc(:,iirrm), zt )

           ! Sum total of Nrm microphysics (auto + cond)
           call stat_update_var( iNrm_mc, hydromet_mc(:,iiNrm), zt )

         end if ! lstats_samp
#endif /*STATS*/

        end if ! coamps micro or KK rain.

!-----------------------------------------------------------------------
!       Loop over all hydrometeor species and apply sedimentation,
!       advection and diffusion.
!-----------------------------------------------------------------------
        if ( hydromet_dim > 0 ) then

          do i = 1, hydromet_dim
#ifdef STATS
            select case( trim( hydromet_list(i) ) )
            case( "rrm" )
              ixrm_bt = irrm_bt
              ixrm_cl = irrm_cl
            case( "Nrm" )
              ixrm_bt = iNrm_bt
              ixrm_cl = iNrm_cl
            case( "ricem" )
              ixrm_bt = iricem_bt
              ixrm_cl = iricem_cl
            case( "rsnowm" )
              ixrm_bt = irsnowm_bt
              ixrm_cl = irsnowm_cl
            case( "rgraupelm" )
              ixrm_bt = irgraupelm_bt
              ixrm_cl = irgraupelm_cl
            case default
              ixrm_bt = 0
              ixrm_cl = 0
            end select

           if ( lstats_samp ) then
            call stat_begin_update
     .           ( ixrm_bt, real(hydromet(:,i) / dt), zt )
           end if
        
#endif /*STATS*/

            call microphys_lhs
     .           ( trim( hydromet_list(i) ), lsed, dt, Kr, nu_r, wmt,
     .             hydromet_vel(:,i), hydromet(:,i), lhs )

            call microphys_solve
     .           ( trim( hydromet_list(i) ), dt, lhs, hydromet_mc(:,i), 
     .             hydromet(:,i), err_code )

            if ( i == iirrm ) then
            ! Handle over-evaporation of rrm and adjust rt and theta-l 
            ! hydrometeor tendency arrays accordingly.
              do k = 1, gr%nnzp, 1
                if ( hydromet(k,i) < 0.0 ) then

                call adj_microphys_tndcy
     .               ( hydromet_mc(:,i), wmt, hydromet_vel(:,i), 
     .                 Kr, nu_r, dt, k, .true., "rrm",
     .                 hydromet(:,i), overevap_rate )

              ! overevap_rate is defined as positive.
              ! It is a correction factor.
              rtm_mc(k)  = rtm_mc(k) - overevap_rate

              thlm_mc(k) = thlm_mc(k)
     .                  + ( Lv / ( Cp*exner(k) ) ) * overevap_rate

              ! Moved from adj_microphys_tndcy  
#ifdef STATS
              if ( lstats_samp ) then
                call stat_update_var_pt( irrm_cond_adj, k, 
     .                                   overevap_rate, zt )
              end if

              else
#endif /*STATS*/

#ifdef STATS
              if ( lstats_samp ) then
                call stat_update_var_pt( irrm_cond_adj, k, 
     .                                  0.0, zt )
              end if
             ! Joshua Faschinj December 2007
#endif /*STATS*/

              end if 
            end do ! k=1..gr%nnzp

            else if ( i == iiNrm ) then
            ! Handle over-evaporation similar to rrm.  However, in the case 
            ! of Nrm there is no effect on rtm or on thlm.  
            ! Brian Griffin.  April 14, 2007.
              do k = 1, gr%nnzp, 1
                if ( hydromet(k,i) < 0.0 ) then

                  call adj_microphys_tndcy
     .               ( hydromet_mc(:,i), wmt, hydromet_vel(:,i), 
     .                 Kr, nu_r, dt, k, .true., "Nrm",
     .                 hydromet(:,i), overevap_rate )

              ! Moved from adj_microphys_tndcy
#ifdef STATS
              call stat_update_var_pt( iNrm_cond_adj, k, 
     .                                 overevap_rate, zt )

#endif /*STATS*/

#ifdef STATS
               else
                 if( lstats_samp ) then      
                   call stat_update_var_pt( iNrm_cond_adj,k, 0.0, zt )
                 end if
#endif /*STATS*/
              end if ! ! Nrm(k) < 0
              ! Joshua Fasching December 2007
            end do

            end if
#ifdef STATS
          if ( lstats_samp ) then
            call stat_begin_update
     .           ( ixrm_cl, real( hydromet(:,i) / dt ), zt )
          end if
#endif /*STATS*/
          ! Clip to zero
          where ( hydromet(:,i) < 0.0 ) hydromet(:,i) = 0.0

#ifdef STATS
          if ( lstats_samp ) then

            ! Effects of clipping
            call stat_end_update
     .           ( ixrm_cl, real( hydromet(:,i) / dt ), zt )

           ! Total time tendency
            call stat_end_update
     .           ( ixrm_bt, real(hydromet(:,i) / dt), zt )

          end if ! lstats_samp
#endif /*STATS*/

          end do ! i=1..hydromet_dim
        end if ! hydromet_dim > 0


        ! Call the ice diffusion scheme
        if ( licedfs ) then
          call ice_dfsn( dt, T_in_K, rcm, p, exner, rhot, rtm_mc )
          thlm_mc = - ( Lv/(Cp*exner) ) * rtm_mc
        end if

#ifdef STATS
        if ( lstats_samp ) then
          call stat_update_var( iNcm, Ncm, zt )

          call stat_update_var( iNcnm, Ncnm, zt )

          call stat_update_var( iNim, Nim, zt )

          if ( iirrm > 0 ) then
            call stat_update_var( irrm, hydromet(:,iirrm), zt )
          end if

          if ( iiNrm > 0 ) then
            call stat_update_var( iNrm, hydromet(:,iiNrm), zt )
          end if
          
          if ( iirsnowm > 0 ) then
            call stat_update_var( irsnowm, hydromet(:,iirsnowm), zt )
          end if

          if ( iiricem > 0 ) then
            call stat_update_var( iricem, hydromet(:,iiricem), zt )
          end if

          if ( iirgraupelm > 0 ) then
            call stat_update_var( irgraupelm, 
     .                            hydromet(:,iirgraupelm), zt )
          end if

          call stat_update_var( ithlm_mc, thlm_mc(:), zt )

          call stat_update_var( irtm_mc, rtm_mc(:), zt )

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
          call stat_update_var( irain_rate, 
     .       ( hydromet(:,iirrm) 
     .         * zm2zt( abs( hydromet_vel(:,iirrm) ) ) )
     .          * ( rhot / rho_lw )
     .        * real( sec_per_day * 1000.0 ), zt )

          ! Precipitation Flux (W/m^2) should be defined on
          ! momentum levels.  -Brian
          ! Normally, a flux is a vector quantity.  Since rain obviously
          ! falls downward, the sign of the flux would normally be negative.
          ! However, it is generally a convention in meteorology to show
          ! Precipitation Flux as a positive downward quantity.  Thus, the
          ! absolute value of vrr is taken.
          call stat_update_var( iFprec, 
     .      ( zt2zm( hydromet(:,iirrm) ) 
     .        * abs( hydromet_vel(:,iirrm) ) )
     .      * ( rhom / rho_lw ) * rho_lw * Lv, zm )

          ! Store values of surface fluxes for statistics
          ! See notes above.
          call stat_update_var_pt( irain, 1, 
     .        ( hydromet(2,iirrm)
     .          * abs( zm2zt( hydromet_vel(:,iirrm), 2 ) ) )
     .            * ( rhot(2) / rho_lw )
     .        * real( sec_per_day * 1000.0 ), sfc )

          call stat_update_var_pt( ipflux, 1,
     .       ( zt2zm( hydromet(:,iirrm), 1 ) 
     .         * abs( hydromet_vel(1,iirrm) ) ) * ( rhom(1) / rho_lw ) 
     .         * rho_lw * Lv, sfc )

          ! Also store the value of surface rain water mixing ratio.
          call stat_update_var_pt( irrm_sfc, 1, 
     .           ( zt2zm( hydromet(:,iirrm), 1 ) ), sfc )

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

        use stats_precision, only: 
     .      time_precision ! Variable(s)

        use lapack_wrap, only: 
     .      tridag_solve ! Procedure(s)
!    .     ,band_solve

#ifdef STATS
        use stats_variables, only:
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

        use stats_type, only: stat_update_var_pt ! Procedure(s)
#endif

        implicit none

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        real(kind=time_precision), intent(in) :: dt ! Timestep     [s]

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
     .  = real((xrm(2:gr%nnzp-1) / dt ) ! Time tendency
     .  + xrm_tndcy(2:gr%nnzp-1))

        ! Boundary condition on the RHS
        rhs(1) = real( xrm(1) / dt )
        rhs(gr%nnzp) = 
     .     real( ( xrm(gr%nnzp) / dt ) + xrm_tndcy(gr%nnzp-1) )

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
             call stat_update_var_pt( ixrm_ma, k,
     .             ztscr01(k) * xrm(km1)
     .           + ztscr02(k) * xrm(k)
     .           + ztscr03(k) * xrm(kp1), zt)
             
             call stat_update_var_pt( ixrm_sd, k,
     .             ztscr04(k) * xrm(km1)
     .           + ztscr05(k) * xrm(k)
     .           + ztscr06(k) * xrm(kp1), zt )
             
             call stat_update_var_pt( ixrm_dff, k,
     .             ztscr07(k) * xrm(km1)
     .           + ztscr08(k) * xrm(k)
     .           + ztscr09(k) * xrm(kp1), zt )

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

!===============================================================================
        subroutine microphys_lhs
     .             ( solve_type, lsed, dt, Kr, nu, wmt, V_hm, xrm, lhs )

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
        
        use stats_precision, only: 
     .      time_precision ! Variable(s)
        
        use diffusion, only: 
     .      diffusion_zt_lhs ! Procedure(s)
        
        use mean_adv, only: 
     .      term_ma_zt_lhs ! Procedure(s)
#ifdef STATS 
        use constants, only: sec_per_day ! Variable(s)
      
        use stats_variables, only:
     .      irrm_ma,  ! Variable(s)
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
     .  kp1_tdiag = 1,   ! Thermodynamic superdiagonal index.
     .  k_tdiag   = 2,   ! Thermodynamic main diagonal index.
     .  km1_tdiag = 3    ! Thermodynamic subdiagonal index.

        ! Input Variables
        character(len=*), intent(in) :: solve_type

        logical, intent(in) :: 
     .  lsed ! Whether to add a sedimentation term

        real(kind=time_precision), intent(in) :: 
     .  dt    ! Timestep                                                 [s]

        real, intent(in) :: 
     .  nu    ! Background diffusion coefficient                         [m^2/s]

        real, intent(in), dimension(gr%nnzp) :: 
     .  xrm,  ! Hydrometeor (at current timestep) (thermodynamic levels) [units vary]
     .  wmt,  ! w wind component on thermodynamic levels                 [m/s]
     .  V_hm, ! Sedimentation velocity of hydrometeor (momentum levels)  [m/s]
     .  Kr    ! Eddy diffusivity for hydrometeor on momentum levels      [m^2/s]

        real, intent(out), dimension(3,gr%nnzp) ::
     .  lhs ! Left hand side of tridiagonal matrix

        ! Local Variables
        real, dimension(3) :: tmp

        ! Array indices
        integer :: k, km1

!       integer kp1

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

        ! Reset LHS Matrix for current timestep.
        lhs = 0.0

        ! Setup LHS Matrix
        do k = 1, gr%nnzp, 1

          km1 = max( k-1, 1 )
!          kp1 = min( k+1, gr%nnzp )

          ! Main diagonal
          ! Time Tendency
          lhs(k_tdiag,k) = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

          ! All diagonals

          ! Diffusion
          lhs(kp1_tdiag:km1_tdiag,k)
     .    = lhs(kp1_tdiag:km1_tdiag,k)
     .    + diffusion_zt_lhs( Kr(k), Kr(km1), nu, 
     .                        gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

          ! Mean Advection
          lhs(kp1_tdiag:km1_tdiag,k)
     .    = lhs(kp1_tdiag:km1_tdiag,k)
     .    + term_ma_zt_lhs( wmt(k), gr%dzt(k), k )

          ! Sedimentation
          ! Note: originally pristine ice did not sediment, so it was
          ! setup to be disabled as needed, but now pristine ice does
          ! sediment in the COAMPS case. -dschanen 12 Feb 2007
          if ( lsed ) then
             lhs(kp1_tdiag:km1_tdiag,k)
     .       = lhs(kp1_tdiag:km1_tdiag,k)
     .       + sedimentation( V_hm(k), V_hm(km1), zm2zt(V_hm,k), 
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
              tmp(1:3) = sedimentation( V_hm(k), V_hm(km1), 
     .                                  zm2zt(V_hm,k), gr%dzt(k), k )
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
!       The variable "hm" stands for one of the five hydrometeor 
!       variables currently in the code:  mean rain mixing ratio (rrm),
!       mean rain drop concentration (Nrm), mean ice mixing ratio 
!       (ricem), mean snow mixing ratio (rsnowm), or mean graupel mixing
!       ratio (rgraupelm).  The variable "V_hm" stands for the 
!       sedimentation velocity of the appropriate hydrometeor.
!
!       The d(hm)/dt equation contains a sedimentation term:
!
!       - d(V_hm*hm)/dz.
!
!       This term is solved for completely implicitly, such that:
!
!       - d( V_hm(t) * hm(t+1) )/dz.
!
!       Note:  When the term is brought over to the left-hand side, the
!              sign is reversed and the leading "-" in front of the
!              term is changed to a "+".
!
!       Timestep index (t) stands for the index of the current timestep,
!       while timestep index (t+1) stands for the index of the next
!       timestep, which is being advanced to in solving the d(hm)/dt
!       equation.
!
!       This term is discretized as follows:
!
!       The values of hm are found on the thermodynamic levels, while 
!       the values of V_hm are found on the momentum levels.  The 
!       variable hm is interpolated to the intermediate momentum 
!       levels.  At the intermediate momentum levels, the interpolated 
!       values of hm are multiplied by the values of V_hm.  Then, the
!       derivative of (hm*V_hm) is taken over the central thermodynamic
!       level.
!
!       -----hmp1------------------------------------------------ t(k+1)
!
!       =============hm(interp)=====V_hm========================= m(k)
!
!       -----hm--------------------------------d(V_hm*hm)/dz----- t(k)
!
!       =============hm(interp)=====V_hmm1======================= m(k-1)
!
!       -----hmm1------------------------------------------------ t(k-1)
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
!       When a hydrometeor is sedimented to the ground (or out the lower
!       boundary of the model), it is removed from the atmosphere (or
!       from the model domain).  Thus, the quantity of the hydrometeor 
!       over the entire vertical domain should not be conserved due to 
!       the process of sedimentation.  Thus, not all of the column 
!       totals in the left-hand side matrix should be equal to 0. 
!       Instead, the sum of all the column totals should equal the flux 
!       of hm out the bottom (zm(1) level) of the domain,
!       -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) ), where the factor in 
!       parentheses is the interpolated value of hm at the zm(1) level.
!       Furthermore, most of the individual column totals should sum 
!       to 0, but the 1st and 2nd (from the left) columns should combine
!       to sum to the flux out the bottom of the domain.
!
!       To see that this modified conservation law is satisfied, 
!       compute the sedimentation of hm and integrate vertically.  In 
!       discretized matrix notation (where "i" stands for the matrix 
!       column and "j" stands for the matrix row):
!
!        -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) )
!           =
!        Sum_j Sum_i ( 1/dzt )_i ( d (V_hm * weights_hm) / dz )_ij hm_j.
!
!       The left-hand side matrix, ( d (V_hm * weights_hm) / dz )_ij, is
!       partially written below.  The sum over i in the above equation 
!       removes dzt everywhere from the matrix below.  The sum over j
!       leaves the column totals and the flux at zm(1) that are desired.
!
!       Left-hand side matrix contributions from the sedimentation term
!       (only); first four vertical levels:
!
!            ------------------------------------------------------------------------------------------>
!       k=1 |            0                            0                                0
!           |
!       k=2 | -dzt(k)*V_hm(k-1)*D(k)       +dzt(k)*[ V_hm(k)*B(k)           +dzt(k)*V_hm(k)*A(k)
!           |                                       -V_hm(k-1)*C(k) ]       
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
!       B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
!            = 1 - A(k);
!       C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
!       D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
!            = 1 - C(k).
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

        ! Note:  The code is now written so that V_hm has been pulled 
        !        inside of the derivative.  The sedimentation term is 
        !        now of the form -d(V_hm*hm)/dz, rather than of the form
        !        -V_hm d(hm)/dz.  The term has been re-discretized in a 
        !        conservative manner and the results are listed below, 
        !        with the old code commented out.

        if ( level == 1 ) then

           ! k = 1 (bottom level); lower boundary level; no effects.

           ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) = 0.0

           ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
           lhs(k_tdiag)   = 0.0

           ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
           lhs(km1_tdiag) = 0.0


        elseif ( level > 1 .and. level < gr%nnzp ) then

           ! Most of the interior model; normal conditions.

        ! Vince Larson pulled V_hm inside derivative to make conservative. 
        ! 13 Dec 2007
!
!           ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
!           lhs(kp1_tdiag) 
!     .     = + V_hmzt * dzt * gr%weights_zt2zm(t_above,mk)
!
!           ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
!           lhs(k_tdiag) 
!     .     = + V_hmzt * dzt * (   gr%weights_zt2zm(t_below,mk)
!     .                          - gr%weights_zt2zm(t_above,mkm1)  )
!
!           ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
!           lhs(km1_tdiag) 
!     .     = - V_hmzt * dzt * gr%weights_zt2zm(t_below,mkm1)

           ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) 
     .     = + dzt * V_hm * gr%weights_zt2zm(t_above,mk)

           ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
           lhs(k_tdiag) 
     .     = + dzt * (   V_hm * gr%weights_zt2zm(t_below,mk)
     .                 - V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

           ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
           lhs(km1_tdiag) 
     .     = - dzt * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)

        !  End Vince Larson change


        elseif ( level == gr%nnzp ) then

           ! k = gr%nnzp (top level); upper boundary level; no flux.

           ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
           lhs(kp1_tdiag) = 0.0

           ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
           lhs(k_tdiag)   = 0.0

           ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
           lhs(km1_tdiag) = 0.0


        endif

        return
        end function sedimentation

!===============================================================================
        subroutine adj_microphys_tndcy( xrm_tndcy, wmt, V_hm, Kr, nu,
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
        
        use stats_precision, only:
     .      time_precision ! Variable(s)
        
        use diffusion, only: 
     .      diffusion_zt_lhs ! Procedure(s)
        
        use mean_adv, only: 
     .      term_ma_zt_lhs ! Procedure(s)

        implicit none
!#ifdef STATS
!        use stats_variables
!#endif /*STATS*/
! Joshua Fasching 2007

        ! Input variables.

        ! Hydrometeor microphysical tendency.
        real, dimension(gr%nnzp), intent(in) :: xrm_tndcy  ! [hm_units/s]
        ! Vertical velocity (thermo. levels).
        real, dimension(gr%nnzp), intent(in) :: wmt        ! [m/s]
        ! Sedimentation velocity (interpolated to thermo. levels).
        real, dimension(gr%nnzp), intent(in) :: V_hm       ! [m/s]
        ! Eddy diffusivity for hydrometeors (m-lev).
        real, dimension(gr%nnzp), intent(in) :: Kr         ! [m^2/s]

        real, intent(in) :: nu  ! Diffusion coefficient      [m^2/s]

        real(kind=time_precision), intent(in) :: dt  ! Timestep   [s]

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

!        real :: evap_amt          ! The actual evaporation amount over the t.s.
!        real :: evap_rate         ! The time-averaged rate.
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
           tmp(1:3) = sedimentation( V_hm(k), V_hm(km1), zm2zt(V_hm,k),
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
        xrm_chge = real( tot_tndcy * dt )

        ! The value of xrm at the previous timestep.
        xrm_old = xrm(k) - xrm_chge

        ! The net amount of change in the hydrometeor due to only the
        ! transport (mean advection and diffusion) and sedimentation terms.
        xrm_chge_trsed = real( trnsprt_sed_tndcy * dt )

        ! The new value of the hydrometeor at this timestep due to only
        ! the transport and sedimentation terms.
        xrm_trsed_only = xrm_old + xrm_chge_trsed

        if ( xrm_trsed_only >= 0.0 ) then
           ! The negative value of hydrometeor (xrm) is due ONLY to 
           ! microphysical tendencies, namely the over-evaporation of xrm.
           ! Find the actual amount of the hydrometeor that evaporated during
           ! the timestep to make the value of xrm go to 0.
!           evap_amt = -xrm_trsed_only
           ! Divide by the timestep to find the actual evaporation rate.
!           evap_rate = evap_amt / dt
           ! The amount of the hydrometeor that was artificially excessively
           ! evaporated.  Define as positive.
           overevap_amt = -xrm(k)
           ! Divide by the timestep to find the over-evaporation rate.
           ! Define as positive.  This rate should also be the difference
           ! between the computed evaporation rate (xrm_tndcy) and the
           ! actual evaporation rate (evap_rate).
           overevap_rate = real( overevap_amt / dt )
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
!           evap_amt = 0.0
!           evap_rate = 0.0
           ! The amount of the hydrometeor that was artificially excessively
           ! evaporated.  Define as positive.  In this case, any evaporation
           ! that was computed is considered to be over-evaporation.
           ! Define as positive.
           overevap_amt = real( -xrm_tndcy(k) * dt )
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
