module prescribe_forcings_module
!
!  Description:
!    This module is responsible for managing the reading in and
!    storage of time dependent information for a case.
!
!  References:
!    None
!--------------------------------------------------------------------------------------------------

  implicit none

  public :: prescribe_forcings

  contains 
  
  !----------------------------------------------------------------------
  subroutine prescribe_forcings( gr, nzm, nzt, ngrdcol, &
                                 sclr_dim, edsclr_dim, sclr_idx, &
                                 runtype, sfctype, &
                                 time_current, time_initial, dt, &
                                 um, vm, thlm, &
                                 p_in_Pa, exner, rho, rho_zm, thvm, &
                                 veg_T_in_K, &
                                 l_modify_bc_for_cnvg_test, &
                                 saturation_formula, &
                                 stats,         &
                                 l_add_dycore_grid, &
                                 grid_remap_method, &
                                 total_idx_rho_lin_spline, rho_lin_spline_vals, &
                                 rho_lin_spline_levels, &
                                 gr_dycore, &
                                 rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref, &
                                 thlm_forcing, rtm_forcing, um_forcing, &
                                 vm_forcing, wprtp_forcing, wpthlp_forcing, &
                                 rtp2_forcing, thlp2_forcing, rtpthlp_forcing, &
                                 wpsclrp, sclrm_forcing, edsclrm_forcing, &
                                 wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
                                 T_sfc, p_sfc, sens_ht, latent_ht, &
                                 wpsclrp_sfc, wpedsclrp_sfc, err_info )

    ! Description:
    !   Calculate tendency and surface variables
    ! References:
    !   None
    !----------------------------------------------------------------------

    ! Modules to be included

    use grid_class, only: grid ! Type

    use grid_class, only: zt2zm_api, zm2zt_api !---------------------- Procedure(s)

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    use constants_clubb, only: &
      Cp, Lv, kappa, p0, & !---------------------------------- Variable(s)
      zero, fstderr

    use clubb_precision, only: &
      core_rknd, &   !-------------------- Variable(s)
      time_precision

    use time_dependent_input, only: &
      apply_time_dependent_forcings, &
      apply_time_dependent_forcings_from_dycore, &
      l_t_dependent, &
      l_ignore_forcings

    use array_index, only: &
      sclr_idx_type

    ! Case specific modules
    use arm, only: arm_sfclyr !------------------------------ Procedure(s)

    use arm_0003, only: arm_0003_sfclyr !-------------------- Procedure(s)

    use arm_3year, only: arm_3year_sfclyr !------------------ Procedure(s)

    use astex_a209, only: astex_a209_sfclyr !---------------- astex_a209_tndcy ! Procedure(s)

    use atex, only: atex_tndcy, atex_sfclyr !---------------- Procedure(s)

    use atex_long, only: atex_long_tndcy, atex_long_sfclyr !- Procedure(s)

    use arm_97, only: arm_97_sfclyr !------------------------ Procedure(s)

    use bomex, only: bomex_tndcy, bomex_sfclyr !------------- Procedure(s)

    use clex9_nov02, only: clex9_nov02_read_t_dependent !---- Procedure(s)

    use clex9_oct14, only: clex9_oct14_read_t_dependent !---- Procedure(s)

    use cloud_feedback, only: cloud_feedback_sfclyr !-------- Procedure(s)

    use cobra, only: cobra_sfclyr !-------------------------- Procedure(s)

    use dycoms2_rf01, only: &           !---------------- Procedure(s)
        dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

    use dycoms2_rf02, only: &
        dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr !------------ Procedure(s)

    use fire, only: &
      fire_sfclyr !------------------------------------------ Procedure(s)

    use gabls2, only: gabls2_tndcy, gabls2_sfclyr !---------- Procedure(s)

    use gabls3, only: gabls3_sfclyr !------------------------ Procedures(s)

    use gabls3_night, only: gabls3_night_sfclyr

    use jun25, only: jun25_altocu_read_t_dependent !--------- Procedure(s)
    
    use lba, only: lba_tndcy, lba_sfclyr !------------------- Procedure(s)

    use mpace_a, only: mpace_a_tndcy, mpace_a_sfclyr !------- Procedure(s)

    use mpace_b, only: mpace_b_tndcy, mpace_b_sfclyr !------- Procedure(s)

    use nov11, only: nov11_altocu_rtm_adjust, nov11_altocu_read_t_dependent ! Procedure(s)

    use rico, only: rico_tndcy, rico_sfclyr !---------------- Procedure(s)

    use neutral_case, only: neutral_case_sfclyr   ! Procedure(s)

    use ekman, only: ekman_sfclyr                 ! Procedure(s)

    use twp_ice, only: twp_ice_sfclyr !---------------------- Procedure(s)

    use wangara, only: wangara_tndcy, wangara_sfclyr !------- Procedure(s)

    use sfc_flux, only: &   !--------------------------- Procedure(s)
      compute_momentum_flux, &
      compute_ubar,          &
      set_sclr_sfc_rtm_thlm

    use err_info_type_module, only: &
      err_info_type        ! Type

    use clubb_api_module, only: &
        clubb_fatal_error,              & ! Constant
        clubb_at_least_debug_level_api    ! Procedure

    implicit none

    ! Input Variables
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol, &
      sclr_dim, &
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (grid), intent(in) :: &
      gr
      
    character(len=50), intent(in) ::  & 
      runtype ! String identifying the model case; e.g. bomex

    integer, intent(in) :: &
      sfctype

    real(kind = time_precision ), intent(in) :: & 
      time_initial, & ! Time of start of simulation     [s]
      time_current   !  Current time of simulation      [s]

    real(kind=core_rknd), intent(in) :: &
      dt         ! Model timestep         [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      um,           & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,           & ! northward grid-mean wind component (thermo. levs.) [m/s]
      thlm,         & ! liq. water pot. temp., th_l (thermo. levels)       [K]
      p_in_Pa,      & ! Air pressure (thermodynamic levels)                [Pa]
      exner,        & ! Exner function (thermodynamic levels)              [-]
      rho,          & ! Air density on thermodynamic levels                [kg/m^3]
      thvm            ! Virtual potential temperature                      [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      rho_zm          ! Air density on momentum levels                     [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      veg_T_in_K    

    ! Flag to activate modifications on boundary condition for convergence test
    ! (surface fluxes computed at fixed 25 m height).
    logical, intent(in) :: &
      l_modify_bc_for_cnvg_test

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type(stats_type), intent(inout) :: &
      stats

    logical, intent(in) :: &
      l_add_dycore_grid

    integer, intent(in) :: &
      grid_remap_method

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of indices for the linear spline definition arrays

    real( kind = core_rknd ), dimension(ngrdcol,total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, & ! rho values at the given altitudes
      rho_lin_spline_levels  ! altitudes for the given rho values
    ! Note: both these arrays need to be sorted from low to high altitude

    type( grid ), intent(in) :: &
      gr_dycore ! only allocated if l_add_dycore_grid = .true.

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      rtm,             & ! total water mixing ratio, r_t (thermo. levs.)        [kg/kg]
      wm_zt,           & ! vertical mean wind comp. on thermo. levs             [m/s]
      ug,              & ! u geostrophic wind                                   [m/s]
      vg,              & ! v geostrophic wind                                   [m/s]
      um_ref,          & ! Initial u wind                                       [m/s]
      vm_ref,          & ! Initial v wind                                       [m/s]
      thlm_forcing,    & ! liquid potential temp. forcing (thermodynamic levels)[K/s]
      rtm_forcing,     & ! total water forcing (thermodynamic levels)           [(kg/kg)/s]
      um_forcing,      & ! eastward wind forcing (thermodynamic levels)         [m/s/s]
      vm_forcing         ! northward wind forcing (thermodynamic levels)        [m/s/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wm_zm,           & ! vertical mean wind comp. on momentum levs     [m/s]
      wprtp_forcing,   & ! total water turbulent flux forcing (momentum levels) [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (momentum levels)     [m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (momentum levels)       [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (momentum levels)      [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (momentum levels)     [K(kg/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(inout) :: &
      wpsclrp          ! w'sclr' (momentum levels)       [{units vary} m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(inout) :: &
      sclrm_forcing    ! Passive scalar forcing          [{units vary}/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,edsclr_dim), intent(inout) :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      wpthlp_sfc, & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,  & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,   & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc,   & ! v'w' at surface          [m^2/s^2]
      T_sfc,      & ! surface temperature      [K]
      p_sfc         ! surface pressure         [Pa]

    real( kind = core_rknd ), intent(inout) :: &
      sens_ht,    & ! sensible heat flux       [K m/s]
      latent_ht     ! latent heat flux         [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,sclr_dim), intent(inout) :: &
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), dimension(ngrdcol,edsclr_dim), intent(inout) :: &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface [{un vary}m/s]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      ustar,    & ! Average value of friction velocity [m/s]
      ubar        ! mean sfc wind speed   [m/s]

    ! Flags to help avoid code duplication
    logical :: &
      l_compute_momentum_flux, &
      l_set_sclr_sfc_rtm_thlm, &
      l_fixed_flux            

    ! Variables to store the values at fixed model height
    ! used to calculate surface fluxes in order to implement the
    ! modified boundary conditions for convergence test
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      um_bot, &
      vm_bot, &
      rtm_bot, &
      thlm_bot, &
      rho_bot, &
      exner_bot, &
      z_bot

    integer :: i, k, sclr, edsclr

!-----------------------------------------------------------------------

    !$acc enter data create( um_bot, vm_bot, rtm_bot, thlm_bot, rho_bot, exner_bot, z_bot, ustar, ubar )

!-----------------------------------------------------------------------
!                    FIND ALL DIAGNOSTIC VARIABLES
!-----------------------------------------------------------------------
    l_compute_momentum_flux = .false.
    l_set_sclr_sfc_rtm_thlm = .false.
    l_fixed_flux            = .false.

    !----------------------------------------------------------------
    ! Set vertical velocity, w, and compute large-scale forcings
    !----------------------------------------------------------------

    ! These lines were added to reset the forcing arrays to 0 each iteration.
    ! This was previously done in the <case>_tndcy subroutine.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        rtm_forcing(i,k)  = zero
        thlm_forcing(i,k) = zero
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        wprtp_forcing(i,k)   = zero
        wpthlp_forcing(i,k)  = zero
        rtp2_forcing(i,k)    = zero
        thlp2_forcing(i,k)   = zero
        rtpthlp_forcing(i,k) = zero
      end do
    end do

    if ( l_t_dependent .and. .not. l_ignore_forcings ) then

      ! This should include the following:
      ! "cloud_feedback_s6", "cloud_feedback_s6_p2k",
      !   "cloud_feedback_s11", "cloud_feedback_s11_p2k",
      !   "cloud_feedback_s12", "cloud_feedback_s12_p2k",
      !   "gabls3_night", "arm_97", "gabls3", "twp_ice",
      !   "arm", "arm_0003", "arm_3year", "astex_a209", & "cobra".
      if ( l_add_dycore_grid ) then
        ! if we want to simulate forcings coming in from the host model on the dycore grid,
        ! we have the forcings stored on the dycore grid, so we need to remap them everytime
        ! to the current physics grid
        call apply_time_dependent_forcings_from_dycore( &
                ngrdcol, gr%nzm, gr%nzt, &
                sclr_dim, edsclr_dim, sclr_idx, &
                gr, gr_dycore, time_current, rtm, rho, exner,  &
                grid_remap_method, &
                total_idx_rho_lin_spline, rho_lin_spline_vals, &
                rho_lin_spline_levels, &
                p_sfc(1), &
                thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing, &
                wm_zt, wm_zm,  ug, vg, &
                sclrm_forcing, edsclrm_forcing )
      else
        call apply_time_dependent_forcings( &
              ngrdcol, gr%nzm, gr%nzt, &
              sclr_dim, edsclr_dim, sclr_idx, & ! In
              gr, time_current, rtm, rho, exner, & ! In
              thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing, & ! In/Out
              wm_zt, wm_zm, ug, vg, & ! In/Out
              sclrm_forcing, edsclrm_forcing ) ! In/Out
      end if

      ! Vince Larson set forcing to zero at the top point so that we don't need
      ! so much sponge damping, which is associated with sawtooth noise
      ! in the cloud_feedback cases.  I don't know how it will affect
      ! the other cases.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        rtm_forcing(i,nzt) = zero
        thlm_forcing(i,nzt) = zero
      end do
      ! End Vince Larson's addition

    else ! Legacy method of setting the forcings

      select case ( runtype )

      !  case ( "astex_a209" ) ! ASTEX Sc case for K & K
      !    call astex_a209_tndcy( sclr_dim, edsclr_dim, sclr_idx, &     ! Intent(in)
      !                           wm_zt, wm_zm, &                       ! Intent(out)
      !                           thlm_forcing, rtm_forcing , &         ! Intent(out)
      !                           sclrm_forcing, edsclrm_forcing )      ! Intent(out)

      case ( "atex" ) ! ATEX case

        call atex_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &      ! Intent(in)
                         gr, time_current, time_initial, &      ! Intent(in)
                         rtm, &                                 ! Intent(in)
                         l_add_dycore_grid, &                   ! Intent(in)
                         grid_remap_method, &                   ! Intent(in)
                         gr_dycore, &                           ! Intent(in)
                         total_idx_rho_lin_spline, &            ! Intent(in)
                         rho_lin_spline_vals, &                 ! Intent(in)
                         rho_lin_spline_levels, &               ! Intent(in)
                         p_sfc(1), &                            ! Intent(in)
                         err_info, &                            ! Intent(inout)
                         wm_zt, wm_zm, &                        ! Intent(out)
                         thlm_forcing, rtm_forcing, &           ! Intent(out)
                         sclrm_forcing, edsclrm_forcing )       ! Intent(out)

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
            write(fstderr, *) err_info%err_header_global
            write(fstderr, *) "Calling atex_tndcy in prescribe_forcings"
            return
          end if
        end if

      case ( "atex_long" ) ! Long ATEX case

        call atex_long_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                              gr, time_current, &                        ! Intent(in)
                              wm_zt, wm_zm, &                            ! Intent(out)
                              thlm_forcing, rtm_forcing, &               ! Intent(out)
                              sclrm_forcing, edsclrm_forcing )           ! Intent(out)

      case ( "bomex" ) ! BOMEX Cu case

        call bomex_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &     ! Intent(in)
                          gr, rtm, &                            ! Intent(in)
                          thlm_forcing, rtm_forcing, &          ! Intent(out)
                          sclrm_forcing, edsclrm_forcing )      ! Intent(out)

      case ( "dycoms2_rf01" ) ! DYCOMS2 RF01 case

        call dycoms2_rf01_tndcy( ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                                  thlm_forcing, rtm_forcing, &                  ! Intent(out)
                                  sclrm_forcing, edsclrm_forcing )              ! Intent(out)

      case ( "dycoms2_rf02" ) ! DYCOMS2 RF02 case

        call dycoms2_rf02_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, & ! Intent(in)
                                 wm_zt, wm_zm, &                                ! Intent(inout)
                                 thlm_forcing, rtm_forcing, &                   ! Intent(out)
                                 sclrm_forcing, edsclrm_forcing )               ! Intent(out)

      case ( "fire", "generic", "coriolis_test", "ekman" ) ! FIRE Sc case

        ! Analytic radiation is computed elsewhere
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            thlm_forcing(i,k) = 0._core_rknd
            rtm_forcing(i,k)  = 0._core_rknd
          end do
        end do

      case ( "gabls2" ) ! GABLS 2 case

        call gabls2_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &      ! Intent(in)
                           gr, time_current, time_initial, &      ! Intent(in)
                           wm_zt, wm_zm, thlm_forcing, &          ! Intent(out)
                           rtm_forcing, &                         ! Intent(out)
                           sclrm_forcing, edsclrm_forcing )       ! Intent(out)

      case ( "lba" )

        call lba_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx,  & ! Intent(in)
                        gr, thlm_forcing, rtm_forcing, &            ! Intent(out)
                        sclrm_forcing, edsclrm_forcing )            ! Intent(out)

      case ( "mpace_a" ) ! mpace_a arctic stratus case

        !$acc update host( p_in_Pa )
        call mpace_a_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                            gr, time_current, p_in_Pa, &               ! Intent(in)
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                            um_ref, vm_ref, &                          ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )           ! Intent(out)
        !$acc update device( wm_zt, wm_zm, thlm_forcing, rtm_forcing, um_ref, vm_ref, sclrm_forcing, edsclrm_forcing )

      case ( "mpace_b" ) ! mpace_b arctic stratus case

        call mpace_b_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                            gr, p_in_Pa, thvm, &                       ! Intent(in)
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )           ! Intent(out)

      case ( "rico" ) ! RICO case

        call rico_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                         gr, rtm, exner, &                          ! Intent(in)
                         thlm_forcing, rtm_forcing, &               ! Intent(out)
                         sclrm_forcing, edsclrm_forcing )           ! Intent(out)

      case ( "neutral" )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            thlm_forcing(i,k) = 0.0_core_rknd
            rtm_forcing(i,k) = 0.0_core_rknd
          end do
        end do

        if ( sclr_dim > 0 ) then
          !$acc parallel loop gang vector collapse(2) default(present)
          do sclr = 1, sclr_dim
            do k = 1, gr%nzt
              do i = 1, ngrdcol
                sclrm_forcing(i,k,sclr) = 0.0_core_rknd
              end do
            end do
          end do
        end if

        if ( edsclr_dim > 0 ) then
          !$acc parallel loop gang vector collapse(2) default(present)
          do edsclr = 1, edsclr_dim
            do k = 1, gr%nzt
              do i = 1, ngrdcol
                edsclrm_forcing(i,k,edsclr) = 0.0_core_rknd
              end do
            end do
          end do
        end if

      case ( "wangara" ) ! Wangara dry CBL

          ! compute_momentum
        call wangara_tndcy( ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx, &  ! Intent(in)
                            wm_zt, wm_zm, &                                 ! Intent(out) 
                            thlm_forcing, rtm_forcing, &                    ! Intent(out)
                            sclrm_forcing, edsclrm_forcing )                ! Intent(out)

      case default

        write(unit=fstderr,fmt=*) &
           "prescribe_forcings: Don't know how to handle " &
           //"LS forcing for runtype: "//trim( runtype )
        error stop

      end select

    end if ! l_t_dependent



    !----------------------------------------------------------------
    ! Compute Surface Fluxes
    !----------------------------------------------------------------

    ! A new subrutine is added here to derive the physical quantities at
    ! the bottom model level that are used for computing the surface fluxes
    ! (i.e. boundary conditions)
    call read_surface_var_for_bc( gr, ngrdcol,                      & ! Intent
                                  um, vm, rtm, thlm, rho_zm, exner, & ! Intent(in)
                                  p_sfc, l_modify_bc_for_cnvg_test, & ! Intent(in)
                                  z_bot, um_bot, vm_bot, rtm_bot,   & ! Intent (out)
                                  thlm_bot, rho_bot, exner_bot )      ! Intent (out)

    ! Boundary conditions for the second order moments
    call compute_ubar( ngrdcol, um_bot, vm_bot, &
                       ubar )

    select case ( trim( runtype ) )

      case ( "rico" )

        l_set_sclr_sfc_rtm_thlm = .true.
        call rico_sfclyr( ngrdcol, time_current, um_bot, vm_bot, thlm_bot, rtm_bot, & ! Intent(in)
                            ! 299.8_core_rknd K is the RICO T_sfc;
                            ! 101540 Pa is the sfc pressure.
                            !gr%zt(1,2), 299.8_core_rknd, 101540._core_rknd, &        ! Intent(in)
                          z_bot, p_sfc, exner_bot, &                      ! Intent(in)
                          saturation_formula, &                           ! Intent(in)
                          upwp_sfc, vpwp_sfc, wpthlp_sfc, &               ! Intent(out)
                          wprtp_sfc, ustar, T_sfc )                       ! Intent(out)

      case ( "gabls3" )

        l_compute_momentum_flux = .true.
        call gabls3_sfclyr( ngrdcol, ubar, veg_T_in_K,           & ! Intent(in)
                            thlm_bot, rtm_bot, z_bot, exner_bot, & ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar )         ! Intent(out)

      case ( "gabls3_night" )

        call gabls3_night_sfclyr( ngrdcol, time_current, um_bot, vm_bot,  & ! Intent(in)
                                  thlm_bot, rtm_bot, z_bot,               & ! Intent(in)
                                  upwp_sfc, vpwp_sfc,                     & ! Intent(out)
                                  wpthlp_sfc, wprtp_sfc, ustar )            ! Intent(out)

      case ( "jun25_altocu" )
        ! There are no surface momentum or heat fluxes
        ! for the Jun. 25 Altocumulus case.

        ! Ensure ustar(i) is set
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          ustar(i) = 0._core_rknd
        end do

        ! Read in time dependent inputs
        call jun25_altocu_read_t_dependent( time_current, &       ! Intent(in)
                                            sens_ht, latent_ht )  ! Intent(out)

      case ( "cobra" )

        l_compute_momentum_flux = .true.

        call cobra_sfclyr( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, & ! Intent(in)
                          time_current, z_bot, rho_bot,            & ! Intent(in)
                          thlm_bot, ubar,                          & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar,            & ! Intent(out)
                          wpsclrp_sfc, wpedsclrp_sfc, T_sfc )        ! Intent(out)

      case ( "clex9_nov02" )
        ! There are no surface momentum or heat fluxes
        ! for the CLEX-9: Nov. 02 Altocumulus case.

        ! Ensure ustar is set
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          ustar(i) = 0._core_rknd
        end do

        ! Read in time dependent inputs
        call clex9_nov02_read_t_dependent( time_current, &      ! Intent(in)
                                          sens_ht, latent_ht ) ! Intent(out)

      case ( "clex9_oct14" )
        ! There are no surface momentum or heat fluxes
        ! for the CLEX-9: Oct. 14 Altocumulus case.

        ! Ensure ustar is set.
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          ustar(i) = 0._core_rknd
        end do

        ! Read in time dependent inputs
        call clex9_oct14_read_t_dependent( time_current, &      ! Intent(in)
                                          sens_ht, latent_ht ) ! Intent(out)

      case ( "astex_a209" )

        l_compute_momentum_flux = .true.
        call astex_a209_sfclyr( ngrdcol, time_current, ubar, rtm_bot, & ! Intent(in)
                                thlm_bot, z_bot, exner_bot, p_sfc,    & ! Intent(in)
                                saturation_formula,                   & ! Intent(in)
                                wpthlp_sfc, wprtp_sfc, ustar, T_sfc )   ! Intent(out)

      case ( "nov11_altocu" )
        ! There are no surface momentum or heat fluxes
        ! for the Nov. 11 Altocumulus case.

        ! Ensure ustar is set
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          ustar(i) = 0._core_rknd
        end do

        ! However, the Nov. 11 Altocumulus case has a one-time adjustment
        ! of rtm at t=3600s after the start of the simulation.
        ! As the nov11_altocu_tndcy subroutine is now obsolete, this was
        ! moved to a separate subroutine, nov11_altocu_rtm_adjust.
        ! This subroutine is called here, as the surface momentum/heat fluxes
        ! are called every timestep.
        ! ~EIHoppe/20110104
        call nov11_altocu_rtm_adjust( ngrdcol, gr,                    & ! (in)
                                      time_current, time_initial, dt, & ! (in)
                                      rtm )                             ! (inout)

        ! Read in time dependent inputs
        call nov11_altocu_read_t_dependent( time_current, &      ! Intent(in)
                                            sens_ht, latent_ht ) ! Intent(out)

      case ( "fire", "generic" )  ! Generic setup, and GCSS FIRE

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        l_fixed_flux            = .true.
        call fire_sfclyr( ngrdcol, time_current, ubar, p_sfc,  & ! Intent(in)
                          thlm_bot, rtm_bot, exner_bot,        & ! Intent(in)
                          saturation_formula,                  & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc )  ! Intent(out)

      case ( "cloud_feedback_s6", "cloud_feedback_s6_p2k",  &
            "cloud_feedback_s11", "cloud_feedback_s11_p2k", &
            "cloud_feedback_s12", "cloud_feedback_s12_p2k", &
            "cgils_s6", "cgils_s6_p2k", "cgils_s11",        &
            "cgils_s11_p2k", "cgils_s12", "cgils_s12_p2k"  ) ! Cloud Feedback cases

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        l_fixed_flux            = .true.
        call cloud_feedback_sfclyr( ngrdcol, time_current, sfctype, & ! Intent(in)
                                    thlm_bot, rtm_bot, z_bot,       & ! Intent(in)
                                    ubar, p_sfc, T_sfc,             & ! Intent(in)
                                    saturation_formula,             & ! Intent(in)
                                    wpthlp_sfc, wprtp_sfc, ustar)     ! Intent(out)

      case ( "arm" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call arm_sfclyr( ngrdcol, time_current, z_bot,    & ! Intent(in)
                          thlm_bot, ubar,                 & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar )    ! Intent(out)

      case ( "arm_0003" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call arm_0003_sfclyr( ngrdcol, time_current, z_bot,   & ! Intent(in)
                              rho_bot, thlm_bot, ubar,        & ! Intent(in)
                              wpthlp_sfc, wprtp_sfc, ustar )    ! Intent(out)

      case ( "arm_3year" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call arm_3year_sfclyr( ngrdcol, time_current, z_bot, rho_bot, & ! Intent(in)
                              thlm_bot, ubar,                         & ! Intent(in)
                              wpthlp_sfc, wprtp_sfc, ustar )           ! Intent(out)

      case ( "arm_97", "mc3e" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call arm_97_sfclyr( ngrdcol, time_current, z_bot, rho_bot, & ! Intent(in)
                            thlm_bot, ubar,                        & ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar )           ! Intent(out)

      case ( "atex" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call atex_sfclyr( ngrdcol, time_current, ubar,  &       ! Intent(in)
                          thlm_bot, rtm_bot, exner_bot, &       ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc ) ! Intent(out)

      case ( "atex_long" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call atex_long_sfclyr( ngrdcol, time_current, ubar,  &          ! Intent(in)
                               thlm_bot, rtm_bot, exner_bot, rho_bot, & ! Intent(in)
                               wpthlp_sfc, wprtp_sfc, ustar, T_sfc )    ! Intent(out)

      case ( "bomex" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call bomex_sfclyr( ngrdcol, time_current, rtm_bot,  & ! Intent(in)
                          wpthlp_sfc, wprtp_sfc, ustar )      ! Intent(out)

      case ( "dycoms2_rf01" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call dycoms2_rf01_sfclyr( ngrdcol, time_current, sfctype, p_sfc,  & ! Intent(in)
                                  exner_bot, ubar,                        & ! Intent(in)
                                  thlm_bot, rtm_bot, rho_bot,             & ! Intent(in)
                                  saturation_formula,                     & ! Intent(in)
                                  wpthlp_sfc, wprtp_sfc, ustar, T_sfc )     ! Intent(out)
      case ( "dycoms2_rf02" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call dycoms2_rf02_sfclyr( ngrdcol, time_current, & ! Intent(in)
                                  wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)

      case ( "gabls2" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call gabls2_sfclyr( ngrdcol, time_current, time_initial,    & ! Intent(in)
                            z_bot, p_sfc,                           & ! Intent(in)
                            ubar, thlm_bot, rtm_bot, exner_bot,     & ! Intent(in)
                            saturation_formula,                     & ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )     ! Intent(out)

      case ( "lba" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call lba_sfclyr( ngrdcol, time_current, time_initial, & ! Intent(in)
                        z_bot, rho_bot, thlm_bot, ubar, &       ! Intent(in)
                        wpthlp_sfc, wprtp_sfc, ustar )          ! Intent(out)

      case ( "mpace_a" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        !$acc update host( rho_bot )
        call mpace_a_sfclyr( ngrdcol, time_current, rho_bot,  & ! Intent(in)
                             wpthlp_sfc, wprtp_sfc, ustar )     ! Intent(out)
        !$acc update device( wpthlp_sfc, wprtp_sfc, ustar )

      case ( "mpace_b" )
        
        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call mpace_b_sfclyr( ngrdcol, time_current, rho_bot, &        ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)

      case ( "neutral" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call neutral_case_sfclyr( ngrdcol, time_current,         & ! Intent(in)
                                  !  z_bot(i), thlm_bot(i),      & ! Intent(in)
                                  um_bot, vm_bot, ubar,          & ! Intent(in)
                                  upwp_sfc, vpwp_sfc,            & ! Intent(out)
                                  wpthlp_sfc, wprtp_sfc, ustar )   ! Intent(out)

      case ( "ekman" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call ekman_sfclyr( ngrdcol, z_bot,                & ! Intent(in)
                           um_bot, vm_bot, ubar,          & ! Intent(in)
                           upwp_sfc, vpwp_sfc,            & ! Intent(out)
                           wpthlp_sfc, wprtp_sfc, ustar )   ! Intent(out)

      case ( "twp_ice" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call twp_ice_sfclyr( ngrdcol, time_current, z_bot, exner_bot, & ! Intent(in)
                            thlm_bot, ubar, rtm_bot, p_sfc,          & ! Intent(in)
                            saturation_formula,                      & ! Intent(in)
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )      ! Intent(out)

      case ( "wangara" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .true.
        call wangara_sfclyr( ngrdcol, time_current, &                 ! Intent(in)
                              wpthlp_sfc, wprtp_sfc, ustar ) ! Intent(out)

      case ( "coriolis_test" )

        l_compute_momentum_flux = .true.
        l_set_sclr_sfc_rtm_thlm = .false.
        l_fixed_flux            = .true.

        ! Ensure ustar is set
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          ustar(i) = 0._core_rknd
        end do

      case default

        write(unit=fstderr,fmt=*)  &
          "Invalid value of runtype = ", runtype
        error stop

    end select ! runtype

    ! These have been placed here to help avoid repetition in the cases
    if( l_compute_momentum_flux ) then
      call compute_momentum_flux( ngrdcol, um_bot, vm_bot, ubar, ustar, & ! Intent(in)
                                  upwp_sfc, vpwp_sfc )           ! Intent(out)
    end if

    if( l_set_sclr_sfc_rtm_thlm ) then
      call set_sclr_sfc_rtm_thlm( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                                  wpthlp_sfc, wprtp_sfc, &      ! Intent(in)
                                  wpsclrp_sfc, wpedsclrp_sfc )  ! Intent(out)
    end if

    ! If the surface type is 0, use fixed fluxes
    if ( sfctype == 0 .and. l_fixed_flux ) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wpthlp_sfc(i) = sens_ht
        wprtp_sfc(i)  = latent_ht
      end do

      if ( sclr_idx%iisclr_thl > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            wpsclrp(i,k,sclr_idx%iisclr_thl) = sens_ht
          end do
        end do
      end if

      if ( sclr_idx%iisclr_rt > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            wpsclrp(i,k,sclr_idx%iisclr_rt)   = latent_ht
          end do
        end do
      end if

    end if

    ! Store values of surface fluxes for statistics
    if ( stats%l_sample ) then
      ! Keep host copies in sync before surface-flux stats updates.
      ! $acc update host( wpthlp_sfc, rho_zm, wprtp_sfc, upwp_sfc, vpwp_sfc, ustar, T_sfc )
      !$acc update host( wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, ustar, T_sfc, rho_zm )
      call stats_update( "sh", wpthlp_sfc(:)*rho_zm(:,1)*Cp, stats )
      call stats_update( "lh", wprtp_sfc(:)*rho_zm(:,1)*Lv, stats )
      call stats_update( "wpthlp_sfc", wpthlp_sfc(:), stats )
      call stats_update( "wprtp_sfc", wprtp_sfc(:), stats )
      call stats_update( "upwp_sfc", upwp_sfc(:), stats )
      call stats_update( "vpwp_sfc", vpwp_sfc(:), stats )
      call stats_update( "ustar", ustar(:), stats )
      call stats_update( "T_sfc", T_sfc(:), stats )
    end if

    !$acc exit data delete( um_bot, vm_bot, rtm_bot, thlm_bot, rho_bot, exner_bot, z_bot, ustar, ubar )

    return

  end subroutine prescribe_forcings

  !-------------------------------------------------------------------------------
  subroutine read_surface_var_for_bc( gr, ngrdcol,                      & ! Intent(in)
                                      um, vm, rtm, thlm, rho_zm, exner, & ! Intent(in)
                                      p_sfc, l_modify_bc_for_cnvg_test, & ! Intent(in)
                                      z_bot, um_bot, vm_bot, rtm_bot,   & ! Intent (out)
                                      thlm_bot, rho_bot, exner_bot )      ! Intent (out)

    ! Description:
    ! Derives the physical quantities at the bottom model level for calculating
    ! surface fluxes (boundary conditions). The default option is to use the
    ! quantities at first/second model level. When l_modify_bc_for_cnvg_test =
    ! .true., the quantities at a fixed model height (25m) is obtained via
    ! vertical interpolation and used for calculating the surface fluxes. The
    ! purpose is to eleminate the space-dependence of quantities in default option
    ! when model is refined vertically, which results in a space-dependence of
    ! surface fluxes. The modified option is found to be correct treatment for
    ! evaluating space-time convergence in CLUBB-SCM.
    !
    ! Author: Shixuan Zhang (Shixuan.Zhang@pnnl.gov).

    use clubb_precision, only: &
        core_rknd,    & !------------------- Constants
        time_precision

    use interpolation, only: &
        mono_cubic_interp  ! Procedure(s)

    use grid_class, only: &
        grid  ! Type

    use grid_class, only: &
        zt2zm_api,  & ! Procedure(s)
        zm2zt_api

    use constants_clubb, only: &
        fstderr, & ! Constant
        !fstdout, &
        p0,      & ! Reference pressure of 100000 Pa             [Pa]
        kappa      ! Rd/Cp                                       [-]

    use clubb_api_module, only: &
        clubb_at_least_debug_level_api ! Error indicator

    implicit none

    integer, intent(in) :: &
      ngrdcol

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      um,           & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,           & ! northward grid-mean wind component (thermo. levs.) [m/s]
      rtm,          & ! total water mixing ratio, r_t (thermo. levs.) [kg/kg]
      thlm,         & ! liq. water pot. temp., th_l (thermo. levels)       [K]
      exner           ! Exner function (thermodynamic levels)              [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      rho_zm          ! Air density on momentum levels [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc    ! Surface pressure         [Pa]

    logical, intent(in) :: &
      l_modify_bc_for_cnvg_test ! Flag to activate modifications on boundary condition for
                                ! convergence test (surface fluxes computed at fixed 25m height)

    ! the variable at a fixed model hight for the derivation of surface fluxes
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
      z_bot,        & ! height at bottom model level [m]
      um_bot,       & ! um at bottom model level (thermo. levs.) [m/s]
      vm_bot,       & ! vm at bottom model level (thermo. levs.) [m/s]
      rtm_bot,      & ! rtm at bottom model level (thermo. levs.) [kg/kg]
      thlm_bot,     & ! thlm at bottom model level (thermo. levels) [K]
      rho_bot,      & ! rho at bottom model level (momentum levels) [kg/m^3]
      exner_bot       ! exner at bottom model level (thermodynamic levels) [-]

    ! Options for finding the fixed model height
    integer, parameter :: &
      constant_height_option   = 2 ! option 1: find the nearest level
                                   ! option 2: interpolate to the constant
                                   ! height level

    ! Local variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      um_zm,   &
      vm_zm,   &
      exner_zm,&
      rtm_zm,  &
      thlm_zm

    real( kind = core_rknd ) :: &
      min_val

    integer, dimension(ngrdcol) :: &
      k_min

    integer :: km1,kp1,kp2,k00, i, k

    !$acc data create( um_zm, vm_zm, exner_zm, rtm_zm, thlm_zm, k_min )

    if ( .not. l_modify_bc_for_cnvg_test ) then

      ! Default model setup in CLUBB-SCM
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        z_bot(i)     = gr%zt(i,1)
        um_bot(i)    = um(i,1)
        vm_bot(i)    = vm(i,1)
        rtm_bot(i)   = rtm(i,1)
        thlm_bot(i)  = thlm(i,1)
        rho_bot(i)   = rho_zm(i,1)
        exner_bot(i) = ( p_sfc(i) / p0 )**kappa
      end do

      !write(unit=fstdout, fmt='(a,f10.4)')'Surface fluxes calculated at height
      !of = ', z_bot
      !write(unit=fstdout, fmt='(a,f10.4)')'The nearest zt-levels to z_bot = ',
      !gr%zt(1,2)
      !fstdout constant is commented out
    else

      ! Modified option which find the values of physical quantities
      ! at a fixed model height (25m)
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        z_bot(i)  = 25.0_core_rknd !user-specified

        min_val  = abs( gr%zt(i,1) - z_bot(i) )
        k_min(i) = 1

        ! find the neareast level to the constant model height
        do k = 2, gr%nzt
          if ( abs(gr%zt(i,k) - z_bot(i) ) < min_val ) then
            min_val  = abs( gr%zt(i,k) - z_bot(i) )
            k_min(i) = k
          end if
        end do
      end do

      if ( clubb_at_least_debug_level_api( 1 ) ) then

        !$acc update host( k_min )
      
        do i = 1, ngrdcol

          if ( (k_min(i) < 1) .or. (k_min(i) > gr%nzt) ) then
            write(fstderr,*) "Sanity check failed! constant model height is not properly set"
            write(fstderr,*) "in get_fixed_height_values at i = ", i
          end if
          
        end do

      end if

      if (constant_height_option == 1) then ! option 1 (non-interpolation)

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          um_bot(i)    = um(i,k_min(i))
          vm_bot(i)    = vm(i,k_min(i))
          rtm_bot(i)   = rtm(i,k_min(i))
          thlm_bot(i)  = thlm(i,k_min(i))
          rho_bot(i)   = rho_zm(i,k_min(i))
          exner_bot(i) = exner(i,k_min(i))
        end do

      else ! option 2 (interpolation)

        um_zm       = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, um )
        vm_zm       = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, vm )
        thlm_zm     = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, thlm )
        rtm_zm      = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, rtm )
        exner_zm    = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, exner)

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol

          exner_zm(i,1) = ( p_sfc(i) / p0 )**kappa

          ! use the mono cubic interpolation to get the values
          if ( k_min(i) == 1 ) then
            km1 = 1
            k00 = 1
            kp1 = 2
            kp2 = 3
          else
            km1 = k_min(i)-1
            k00 = k_min(i)
            kp1 = k_min(i)+1
            kp2 = k_min(i)+2
          end if
          
          rho_bot(i)   = rho_zm(i,k_min(i))

          um_bot(i)    = mono_cubic_interp( z_bot(i), km1, k00, kp1, kp2, &
                                            gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                                            um_zm(i,km1), um_zm(i,k00), um_zm(i,kp1), um_zm(i,kp2) )

          vm_bot(i)    = mono_cubic_interp( z_bot(i), km1, k00, kp1, kp2, &
                                            gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                                            vm_zm(i,km1), vm_zm(i,k00), vm_zm(i,kp1), vm_zm(i,kp2) )

          exner_bot(i) = mono_cubic_interp( z_bot(i), km1, k00, kp1, kp2, &
                                            gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                                            exner_zm(i,km1), exner_zm(i,k00), exner_zm(i,kp1), exner_zm(i,kp2) )

          thlm_bot(i)  = mono_cubic_interp( z_bot(i), km1, k00, kp1, kp2, &
                                            gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                                            thlm_zm(i,km1), thlm_zm(i,k00), thlm_zm(i,kp1), thlm_zm(i,kp2) )

          rtm_bot(i)   = mono_cubic_interp( z_bot(i), km1, k00, kp1, kp2, &
                                            gr%zm(i,km1), gr%zm(i,k00), gr%zm(i,kp1), gr%zm(i,kp2), &
                                            rtm_zm(i,km1), rtm_zm(i,k00), rtm_zm(i,kp1), rtm_zm(i,kp2) )
        end do

      end if

      !write(unit=fstdout, fmt='(a,f10.4)')'Surface fluxes calculated at height
      !of = ', z_bot
      !write(unit=fstdout, fmt='(a,f10.4)')'The nearest zt-levels to z_bot = ',
      !gr%zt(1,kk)
      !fstdout constant is commented out

    end if

    !$acc end data

    return

  end subroutine read_surface_var_for_bc

end module prescribe_forcings_module
