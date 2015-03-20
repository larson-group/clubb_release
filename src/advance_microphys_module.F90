! $Id$
!===============================================================================
module advance_microphys_module

  ! Description:
  ! Advance mean precipitating hydrometeors and mean cloud droplet concentration
  ! one time step.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  ! Subroutines
  public :: advance_microphys

  private :: advance_hydrometeor, &
             advance_Ncm, &
             microphys_solve, &
             microphys_lhs, &
             microphys_rhs

  ! Functions
  private :: sed_centered_diff_lhs, &
             sed_upwind_diff_lhs, &
             term_turb_sed_lhs, &
             term_turb_sed_rhs

  private ! Default Scope

  contains

  !=============================================================================
  subroutine advance_microphys( dt, time_current, wm_zt, wp2, &            ! In
                                exner, rho, rho_zm, rcm, &                 ! In
                                cloud_frac, Kh_zm, Skw_zm, &               ! In
                                rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &   ! In
                                hydromet_mc, Ncm_mc, Lscale, &             ! In
                                hydromet_vel_covar_zt_impc, &              ! In
                                hydromet_vel_covar_zt_expc, &              ! In
                                hydromet, hydromet_vel_zt, hydrometp2, &   ! Inout
                                K_hm, Ncm, Nc_in_cloud, rvm_mc, thlm_mc, & ! Inout
                                wphydrometp, wpNcp, err_code )             ! Out

    ! Description:
    ! Advance mean precipitating hydrometeors and mean cloud droplet
    ! concentration one model time step, and calculate some statistics.

    ! References:
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

    use parameters_tunable, only: & 
        c_K_hm,  & ! Variable(s) 
        c_K_hmb, &  
        K_hm_min_coef    

    use parameters_model, only: & 
        hydromet_dim   ! Integer

    use constants_clubb, only: & 
        one,         & ! Constant(s)
        zero,        &
        Lv,          &
        rho_lw,      & 
        fstderr,     &
        sec_per_day, &
        mm_per_m,    &
        eps

    use parameters_microphys, only: &
        l_predict_Nc,          & ! Predict cloud droplet number conc (Morrison)
        microphys_scheme,         & ! The microphysical scheme in use
        microphys_start_time    ! When to start the microphysics [s]

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use error_code, only:  & 
        fatal_error,                & ! Procedure(s)
        clubb_at_least_debug_level, &
        clubb_no_error                ! Constant(s)

    use array_index, only:  & 
        hydromet_list, & ! Names of the hydrometeor species
        hydromet_tol,  & ! Tolerance values for hydrometeor species
        iirrm            ! Variable(s)

    use stats_variables, only: & 
        stats_zt,  & ! Variable(s)
        stats_zm,  &
        stats_sfc, &
        l_stats_samp

    !use stats_variables, only: & 
        !iVrrprrp, & ! Variable(s)
        !iVNrpNrp

    use stats_variables, only: & 
        iprecip_rate_zt, & ! Variable(s) 
        iFprec, &
        iprecip_rate_sfc, & 
        irain_flux_sfc, & 
        irrm_sfc

    use stats_variables, only: & 
        iNcm,         & ! Variable(s)
        iK_hm,        & 
        iNc_in_cloud

    use stats_type_utilities, only: & 
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_clubb_utilities, only: &
        stats_accumulate_hydromet  ! Procedure(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt           ! Model timestep duration         [s]

    real( kind = time_precision ), intent(in) ::  & 
      time_current   ! Current time     [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zt,      & ! w wind component on thermodynamic levels        [m/s]
      wp2,        & ! Variance of vertical velocity (momentum levels) [m^2/s^2]
      exner,      & ! Exner function                                  [-]
      rho,        & ! Density on thermodynamic levels                 [kg/m^3]
      rho_zm,     & ! Density on momentum levels                      [kg/m^3]
      rcm,        & ! Mean cloud water mixing ratio                   [kg/kg]
      cloud_frac, & ! Cloud fraction                                  [-]
      Kh_zm,      & ! Kh Eddy diffusivity on momentum grid            [m^2/s]
      Skw_zm        ! Skewness of w on momentum levels                [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: & 
      hydromet_mc    ! Microphysics tendency for mean hydrometeors  [units/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Ncm_mc, & ! Microphysics tendency for Ncm                     [num/kg/s]
      Lscale    ! Length-scale                                      [m]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(inout) :: &
      hydromet,        & ! Hydrometeor mean, <h_m> (thermo. levels)    [units]
      hydromet_vel_zt, & ! Mean hydrometeor sed. vel. on thermo. levs. [m/s]
      hydrometp2,      & ! Variance of hydrometeor (overall) (m-levs.) [units^2]
      K_hm               ! hm eddy diffusivity on momentum grid        [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water          [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp.   [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)   [(m/s)units]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      wpNcp          ! Covariance < w'N_c' > (momentum levels)   [(m/s)(num/kg)]

    integer, intent(out) :: &
      err_code  ! Exit code (used to check for errors)

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      hydromet_vel    ! Mean hydrometeor sedimentation velocity, <V_xx> [m/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      hydromet_vel_covar,    & ! Covariance of V_xx & x_x (m-levs)  [units(m/s)]
      hydromet_vel_covar_zt    ! Covariance of V_xx & x_x (t-levs)  [units(m/s)]

    ! Turbulent advection for hydrometeors -- down-gradient approximation for
    ! covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
    ! <w'hm'> = - K_hm * d<hm>/dz; and
    ! <w'Nc'> = - K_Nc & d<Nc>/dz;
    ! where the coefficients of diffusion, K_hm and K_Nc are variable and depend
    ! on multiple factors.
    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim) :: &
      K_gamma ! Non-local factor of diffusion (turb. adv.) for hydrometeors [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      K_Nc    ! Coefficient of diffusion (turb. adv.) for Nc                [m^2/s]

    integer :: k, kp1, i    ! Loop indices

    integer, dimension(hydromet_dim) :: &
      err_code_hydromet    ! Exit code (used to check for errors) for hydromet

    integer :: &
      err_code_Ncm    ! Exit code (used to check for errors) for Ncm

    logical, parameter :: &
      l_use_non_local_diff_fac = .false. ! Use a non-local factor for eddy-diffusivity 
                                         ! applied to hydrometeors 

    ! Initialize intent(out) variables -- covariances <w'hm'> (for any
    ! hydrometeor, hm) and <w'Nc'>.
    if ( hydromet_dim > 0 ) then
       wphydrometp = zero
    endif
    wpNcp = zero

    ! Initialize intent(out) variables -- the error code.
    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Return if there is delay between the model start time and start of the
    ! microphysics
    if ( time_current < microphys_start_time ) return

    ! Turbulent advection for hydrometeors -- down-gradient approximation for
    ! covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
    ! <w'hm'> = - K_hm * d<hm>/dz; and
    ! <w'Nc'> = - K_Nc & d<Nc>/dz;
    ! where the coefficients of diffusion, K_hm and K_Nc are variable and depend
    ! on multiple factors.
    if ( hydromet_dim > 0 ) then

       ! Solve for the value of K_hm, the coefficient of diffusion for
       ! hydrometeors.
       do i = 1, hydromet_dim, 1
          do k = 1, gr%nz, 1

             kp1 = min( k+1, gr%nz )

             K_hm(k,i) &
             = c_K_hm * Kh_zm(k) &
               * ( sqrt( hydrometp2(k,i) ) &
                   / max( zt2zm( hydromet(:,i), k ), hydromet_tol(i) ) ) &
               * ( one + abs( Skw_zm(k) ) ) 

             if ( l_use_non_local_diff_fac ) then
               K_gamma(k,i) &
               = one -  c_K_hmb * ( ( zt2zm( Lscale(:) , k ) &
                 / max( zt2zm( hydromet(:,i), k ),hydromet_tol(i) ) )  &
                 * ( gr%invrs_dzm(k) * ( hydromet(kp1,i)  - hydromet(k,i) ) ) ) 

               K_hm(k,i) = K_hm(k,i) * max(K_gamma(k,i), K_hm_min_coef)
             endif

             if ( abs( gr%invrs_dzm(k) &
                       * ( hydromet(kp1,i) - hydromet(k,i) ) ) > eps ) then
                ! Ensure the abs( correlation ) between w and hydromet does not have
                ! a value greater than one.
                K_hm(k,i) &
                = min( K_hm(k,i), &
                       ( sqrt( wp2(k) ) * sqrt( hydrometp2(k,i) ) ) &
                       / abs( gr%invrs_dzm(k) &
                              * ( hydromet(kp1,i) - hydromet(k,i) ) ) )

             endif ! | d<hm>/dz | > 0

          enddo ! k = 1, gr%nz, 1
       enddo ! i = hydromet_dim, 1

    endif ! hydromet_dim > 0

    if ( l_stats_samp ) then
      do i = 1, hydromet_dim, 1
        do k = 1, gr%nz, 1
          call stat_update_var_pt( iK_hm(i), k, K_hm(k,i), stats_zm )
        end do
      end do ! i = hydromet_dim, 1
    end if 


    if ( l_predict_Nc ) then

       ! Solve for the value of K_Nc, the coefficient of diffusion for cloud
       ! droplet concentration.
       do k = 1, gr%nz, 1
          K_Nc(k) = c_K_hm * Kh_zm(k)
       enddo ! k = 1, gr%nz, 1

    endif ! l_predict_Nc

    !------------------------------------------------------------------------
    ! Advance predictive mean precipitating hydrometeors one model timestep.
    ! Loop over all hydrometeor species, apply mean advection, turbulent
    ! advection (diffusion), mean sedimentation, turbulent sedimentation, and
    ! microphysics tendency terms.
    !------------------------------------------------------------------------

    if ( hydromet_dim > 0 ) then

       call advance_hydrometeor( dt, wm_zt, exner, cloud_frac, K_hm, &
                                 rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                                 hydromet_mc, hydromet_vel_covar_zt_impc, &
                                 hydromet_vel_covar_zt_expc, &
                                 hydromet, hydromet_vel_zt, &
                                 hydrometp2, rvm_mc, thlm_mc, &
                                 wphydrometp, hydromet_vel, &
                                 hydromet_vel_covar, hydromet_vel_covar_zt, &
                                 err_code_hydromet )

    endif ! hydromet_dim > 0

    !-----------------------------------------------------------------------
    ! When mean cloud droplet concentration, Ncm, is predicted, apply
    ! sedimentation, advection, and diffusion, and advance Ncm one model
    ! timestep.
    !-----------------------------------------------------------------------

    if ( l_predict_Nc ) then

       ! Nc is predicted.

       call advance_Ncm( dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm, &
                         rho_ds_zt, invrs_rho_ds_zt, Ncm_mc, &
                         Ncm, Nc_in_cloud, &
                         wpNcp, err_code_Ncm )

    else

       ! Nc is prescribed.
       ! The in-cloud mean of cloud droplet concentration, Nc_in_cloud, is
       ! constant in this scenario.  The overall mean of cloud droplet
       ! concentration, Ncm, depends of the in-cloud mean and cloud fraction.

       Ncm = Nc_in_cloud * cloud_frac

    endif ! l_predict_Nc

    if ( l_stats_samp ) then
      call stat_update_var( iNcm, Ncm, stats_zt )
      call stat_update_var( iNc_in_cloud, Nc_in_cloud, stats_zt )
    endif ! l_stats_samp

    if ( l_stats_samp .and. iirrm > 0 ) then

      ! Rainfall rate (mm/day) is defined on thermodynamic levels.  The rainfall
      ! rate is given by the equation:
      ! Rainfall rate (mm/day) = - V_rr * r_r * ( rho_a / rho_lw )
      !                            * ( 86400 s/day ) * ( 1000 mm/m ).
      ! The level average rainfall rate is given by:
      ! < Rainfall rate (mm/day) > = - < V_rr * r_r > * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m );
      ! which can also be written as:
      ! < Rainfall rate (mm/day) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                                * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m ).
      ! Rainfall rate is defined as positive.  Since V_rr is always negative,
      ! the minus (-) sign is necessary.
      call stat_update_var( iprecip_rate_zt,  & 
                            max( - ( hydromet(:,iirrm) &
                                     * hydromet_vel_zt(:,iirrm) &
                                     + hydromet_vel_covar_zt(:,iirrm) ), &
                                 zero ) &
                            * ( rho / rho_lw ) & 
                            * sec_per_day &
                            * mm_per_m, stats_zt )

      ! Precipitation Flux (W/m^2) is defined on momentum levels.  The
      ! precipitation flux is given by the equation:
      ! Precip. flux (W/m^2) = - V_rr * r_r * rho_a * L_v.
      ! The level average precipitation flux is given by:
      ! < Precip. flux (W/m^2) > = - < V_rr * r_r > * rho_a * L_v;
      ! which can also be written as:
      ! < Precip. flux (W/m^2) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                              * rho_a * L_v.
      ! It is generally a convention in meteorology to show Precipitation Flux
      ! as a positive downward quantity, so the minus (-) sign is necessary.
      call stat_update_var( iFprec,  & 
                            max( - ( zt2zm( hydromet(:,iirrm) )  & 
                                     * hydromet_vel(:,iirrm) &
                                     + hydromet_vel_covar(:,iirrm) ), &
                                 zero ) &
                            * rho_zm * Lv, stats_zm )

      ! Store values of surface fluxes for statistics
      ! See notes above.

      if ( trim( microphys_scheme ) /= "morrison" ) then
        call stat_update_var_pt( iprecip_rate_sfc, 1,  & 
                                 max( - ( hydromet(2,iirrm) &
                                          * hydromet_vel_zt(2,iirrm) &
                                          + hydromet_vel_covar_zt(2,iirrm) ), &
                                      zero ) &
                                  * ( rho(2) / rho_lw ) & 
                                  * sec_per_day &
                                  * mm_per_m, stats_sfc )
      endif ! microphys_scheme /= "morrison"

      call stat_update_var_pt( irain_flux_sfc, 1, & 
                               max( - ( zt2zm( hydromet(:,iirrm), 1 )  & 
                                        * hydromet_vel(1,iirrm) &
                                        + hydromet_vel_covar(1,iirrm) ), &
                                    zero ) &
                               * rho_zm(1) * Lv, stats_sfc )

      ! Also store the value of surface rain water mixing ratio.
      call stat_update_var_pt( irrm_sfc, 1,  & 
                               ( zt2zm( hydromet(:,iirrm), 1 ) ), stats_sfc )

    endif ! l_stats_samp and iirrm > 0

    call stats_accumulate_hydromet( hydromet, rho_ds_zt )

    ! Perform error checking.
    do i = 1, hydromet_dim, 1

       if ( fatal_error( err_code_hydromet(i) ) ) then

          if ( clubb_at_least_debug_level(1) ) then

             write(fstderr,*) "Error in hydrometeor field " &
                              // trim( hydromet_list(i) )

             write(fstderr,*) trim( hydromet_list(i) ) // " = ", hydromet(:,i)

          endif ! clubb_at_least_debug_level(1)

          err_code = err_code_hydromet(i)

       endif ! fatal_error( err_code_hydromet(i) )

    enddo ! i = 1, hydromet_dim, 1


    if ( l_predict_Nc ) then

       if ( fatal_error( err_code_Ncm ) ) then

          if ( clubb_at_least_debug_level(1) ) then

             write(fstderr,*) "Error in Ncm"

             write(fstderr,*) "Ncm = ", Ncm

          endif ! clubb_at_least_debug_level(1)

          err_code = err_code_Ncm

       endif ! fatal_error( err_code_Ncm )

    endif ! l_predict_Nc

!       Error Report
!       Joshua Fasching Feb 2008

    if ( fatal_error( err_code ) .and.  &
         clubb_at_least_debug_level( 1 ) ) then

       write(fstderr,*) "Error in advance_microphys"

       write(fstderr,*) "Intent(in)"

       write(fstderr,*) "wm_zt = ", wm_zt
       write(fstderr,*) "exner = ", exner
       write(fstderr,*) "rho = ", rho
       write(fstderr,*) "rho_zm = ", rho_zm
       write(fstderr,*) "cloud_frac = ", cloud_frac
       write(fstderr,*) "Kh_zm = ", Kh_zm
       write(fstderr,*) "rho_ds_zm = ", rho_ds_zm
       write(fstderr,*) "rho_ds_zt = ", rho_ds_zt
       write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt

       write(fstderr,*) "hydromet_mc = ", hydromet_mc

       write(fstderr,*) "Ncm_mc = ", Ncm_mc

       write(fstderr,*) "hydromet_vel_covar_zt_impc = ", &
                        hydromet_vel_covar_zt_impc
       write(fstderr,*) "hydromet_vel_covar_zt_expc = ", &
                        hydromet_vel_covar_zt_expc

       write(fstderr,*) "Intent(inout)"

       write(fstderr,*) "hydromet = ", hydromet
       write(fstderr,*) "hydromet_vel_zt = ", hydromet_vel_zt

       write(fstderr,*) "Ncm = ", Ncm
       write(fstderr,*) "Nc_in_cloud = ", Nc_in_cloud

       write(fstderr,*) "rvm_mc = ", rvm_mc
       write(fstderr,*) "thlm_mc = ", thlm_mc

       write(fstderr,*) "Intent(out)" 

       write(fstderr,*) "wphydrometp = ", wphydrometp

       write(fstderr,*) "wpNcp = ", wpNcp

    endif


    return

  end subroutine advance_microphys

  !=============================================================================
  subroutine advance_hydrometeor( dt, wm_zt, exner, cloud_frac, K_hm, &
                                  rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                                  hydromet_mc, hydromet_vel_covar_zt_impc, &
                                  hydromet_vel_covar_zt_expc, &
                                  hydromet, hydromet_vel_zt, &
                                  hydrometp2, rvm_mc, thlm_mc, &
                                  wphydrometp, hydromet_vel, &
                                  hydromet_vel_covar, hydromet_vel_covar_zt, &
                                  err_code_hydromet )

    ! Description:
    !   Advance each hydrometeor (precipitating hydrometeor) one model time step.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

    use constants_clubb, only: & 
        one_half,       & ! Constant(s)
        zero,           &
        zero_threshold, &
        fstderr

    use parameters_model, only: & 
        hydromet_dim   ! Variable(s)

    use advance_windm_edsclrm_module, only : &
        xpwp_fnc  ! Procedure(s)

    use fill_holes, only: &
        fill_holes_driver,   & ! Procedure(s)
        setup_stats_indices

    use parameters_tunable, only: & 
        nu_hm_vert_res_dep  ! Variable(s)

    use parameters_microphys, only: &
        l_hydromet_sed,    & ! Variable(s)
        l_upwind_diff_sed

    use error_code, only:  & 
        clubb_at_least_debug_level, & ! Procedure(s)
        clubb_no_error

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_update_var,      & ! Procedure(s)
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_variables, only: & 
        stats_zt,           &  ! Variable(s)
        stats_zm,           &
        l_stats_samp

    use array_index, only:  & 
        hydromet_list, & ! Variable(s)
        hydromet_tol,  &
        iirrm,         &
        iiNrm

    use stats_variables, only: & 
        iVrrprrp,     & ! Variable(s)
        iVNrpNrp,     &
        ihydrometp2,  &
        iwphydrometp

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Duration of one model time step         [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zt,      & ! mean w wind component on thermodynamic levels  [m/s]
      exner,      & ! Exner function, (p/p1000mb)^(Rd/Cp)            [-]
      cloud_frac    ! Cloud fraction                                 [-]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      K_hm    ! Coefficient of diffusion (turb. adv.) for hydrometeors [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: & 
      hydromet_mc     ! Change in hydrometeors due to microphysics  [units/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_hm'h_m'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_hm'h_m'> t-levs [units(m/s)]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(inout) :: &
      hydromet,        & ! Hydrometeor mean, <h_m> (thermodynamic levs.) [units]
      hydromet_vel_zt, & ! Mean hydrometeor sed. velocity on thermo. levs. [m/s]
      hydrometp2         ! Variance of hydrometeor (overall) (m-levs.) [units^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water          [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp.   [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: &
      wphydrometp,  & ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]
      hydromet_vel    ! Mean hydrometeor sedimentation velocity, <V_hm>    [m/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(out) :: &
      hydromet_vel_covar,    & ! Covariance of V_hm & h_m (m-levs)  [units(m/s)]
      hydromet_vel_covar_zt    ! Covariance of V_hm & h_m (t-levs)  [units(m/s)]

    integer, dimension(hydromet_dim), intent(out) :: &
      err_code_hydromet    ! Exit code (used to check for errors) for hydromet

    ! Local Variables
    real( kind = core_rknd ), dimension(3,gr%nz) :: & 
      lhs    ! Left hand side array

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      rhs    ! Right hand side vector

    character(len=10) :: hydromet_name

    real( kind = core_rknd ) :: &
      max_velocity    ! Maximum sedimentation velocity         [m/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      hydromet_zm     ! Mean of hydrometeor interp. to m-levs. [units vary]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      ratio_hmp2_on_hmm2    ! Value of <hm'^2> / <hm>^2    [-]

    logical :: l_fill_holes_hm = .true.

    integer :: i, k    ! Loop indices

    ! Stat indices
    integer :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
               ixrm_bt, ixrm_mc


    ! Initialize error code
    err_code_hydromet = clubb_no_error

    ! Loop over each type of precipitating hydrometeor and advance one model
    ! time step.
    do i = 1, hydromet_dim

       ! Set up the stats indices for hydrometeor index i
       call setup_stats_indices( i,                           & ! Intent(in)
                                 ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                 ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                 max_velocity )                 ! Intent(inout)

       if ( l_stats_samp ) then

          ! Update explicit contributions to the hydrometeor species
          call stat_update_var( ixrm_mc, hydromet_mc(:,i), stats_zt )

          ! Save prior value of the hydrometeors for determining total time
          ! tendency.
          if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
             call stat_begin_update( ixrm_bt, &
                                     hydromet(:,i) &
                                     / dt, stats_zt )
          else
             ! The mean hydrometeor field at thermodynamic level k = 1 is simply
             ! set to the value of the mean hydrometeor field at k = 2.
             ! Don't include budget stats for level k = 1.
             do k = 2, gr%nz, 1
                call stat_begin_update_pt( ixrm_bt, k, &
                                           hydromet(k,i) &
                                           / dt, stats_zt )
             enddo
          endif

       endif ! l_stats_samp

       ! Set realistic limits on sedimentation velocities, following the
       ! numbers in the Morrison microphysics.
       do k = 1, gr%nz
          if ( clubb_at_least_debug_level( 1 ) ) then
            ! Print a warning if the velocity has a large magnitude or the
            ! velocity is in the wrong direction.
             if ( hydromet_vel_zt(k,i) < max_velocity .or. &
                  hydromet_vel_zt(k,i) > zero_threshold ) then

                write(fstderr,*) trim( hydromet_list(i) )// &
                                 " velocity at k = ", k, " = ", &
                                 hydromet_vel_zt(k,i), "m/s"
             endif
          endif
          hydromet_vel_zt(k,i) &
          = min( max( hydromet_vel_zt(k,i), max_velocity ), zero_threshold )
       enddo ! k = 1, gr%nz, 1

       ! Interpolate velocity to the momentum grid for a centered difference
       ! approximation of the sedimenation term.
       hydromet_vel(:,i) = zt2zm( hydromet_vel_zt(:,i) )
       hydromet_vel(gr%nz,i) = zero ! Upper boundary condition

       ! Calculate the value of <hm'^2> / <hm>^2.  This will be used to update
       ! <hm'^2> after <hm> has been advanced one model timestep.  This method
       ! is being used because CLUBB does not currently have a predictive
       ! equation for <hm'^2> (hydrometp2).
       hydromet_zm = zt2zm( hydromet(:,i) )

       do k = 1, gr%nz, 1

          if ( hydromet_zm(k) > hydromet_tol(i) ) then

             ! Calculate the ratio of the overall variance of the hydrometeor
             ! to the overall mean of the hydrometeor squared.
             ratio_hmp2_on_hmm2(k,i) = hydrometp2(k,i) / hydromet_zm(k)**2

          else  ! hydromet_zm(k) <= hydromet_tol(i)

             ! The overall mean of the hydrometeor is 0 or is treated as 0.  The
             ! overall variance of the hydrometeor is also 0 in this scenario.
             ! The ratio is undefined, but will be assigned a value of 0.
             ratio_hmp2_on_hmm2(k,i) = zero

          endif ! hydromet_zm(k) > hydromet_tol(i)

       enddo ! k = 1, gr%nz, 1

       ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
       ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
       ! A Crank-Nicholson time-stepping scheme is used for this variable.
       ! This is the portion of the calculation using < h_m > from timestep t. 
       wphydrometp(1:gr%nz-1,i) &
       = - one_half &
           * xpwp_fnc( K_hm(1:gr%nz-1,i)+nu_hm_vert_res_dep(1:gr%nz-1), &
                       hydromet(1:gr%nz-1,i), hydromet(2:gr%nz,i), &
                       gr%invrs_dzm(1:gr%nz-1) )

       ! A zero-flux boundary condition at the top of the model is used for
       ! hydrometeors.
       wphydrometp(gr%nz,i) = zero

       ! Add implicit terms to the LHS matrix
       call microphys_lhs( trim( hydromet_list(i) ), l_hydromet_sed(i), & ! In
                           dt, K_hm(:,i), nu_hm_vert_res_dep, wm_zt,    & ! In
                           hydromet_vel(:,i), hydromet_vel_zt(:,i),     & ! In
                           hydromet_vel_covar_zt_impc(:,i),             & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,       & ! In
                           lhs )                                          ! Out

       ! Set up explicit term in the RHS vector
       call microphys_rhs( trim( hydromet_list(i) ), dt, l_hydromet_sed(i), &
                           hydromet(:,i), hydromet_mc(:,i), &
                           K_hm(:,i), nu_hm_vert_res_dep, cloud_frac, &
                           hydromet_vel_covar_zt_expc(:,i), &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                           rhs )

       !!!!! Advance hydrometeor one time step.
       call microphys_solve( trim( hydromet_list(i) ), l_hydromet_sed(i), &
                             cloud_frac, &
                             lhs, rhs, hydromet(:,i), err_code_hydromet(i) )

    enddo ! i = 1, hydromet_dim, 1

    ! Now that all precipitating hydrometeors have been advanced, fill holes in
    ! hydromet profiles.
    call fill_holes_driver( gr%nz, dt, hydromet_dim,     & ! Intent(in)
                            l_fill_holes_hm,             & ! Intent(in)
                            rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
                            thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

    ! Loop over each type of precipitating hydrometeor and calculate hydrometeor
    ! covariances (<w'hm'> and <V_hm'hm'>) and other quantities requiring the
    ! value of hydromet (<hm>) from the (t+1) timestep.
    do i = 1, hydromet_dim

       ! Set up the stats indices for hydrometeor at index i
       call setup_stats_indices( i,                           & ! Intent(in)
                                 ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                 ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                 max_velocity )                 ! Intent(inout)

       ! Print warning message if any hydrometeor species has a value < 0.
       if ( any( hydromet(:,i) < zero_threshold ) ) then

          hydromet_name = hydromet_list(i)

          if ( clubb_at_least_debug_level( 1 ) ) then
             do k = 1, gr%nz
                if ( hydromet(k,i) < zero_threshold ) then
                   write(fstderr,*) trim( hydromet_name ) //" < ", &
                                    zero_threshold, &
                                    " in advance_microphys at k= ", k
                endif
             enddo
          endif

       endif ! hydromet(:,i) < 0

       ! Lower boundary condition
       ! Hydrometeors that are below the model lower boundary level have
       ! sedimented out of the model domain, and is not conserved.
       if ( hydromet(1,i) < hydromet_tol(i) ) then
          hydromet(1,i) = zero_threshold
       endif

       ! Calculate the value of <hm'^2> (hydrometp2).  This is based on the
       ! ratio of <hm'^2> / <hm>^2 (ratio_hmp2_on_hmm2) that was saved before
       ! hydrometeors were updated.  This method is being used because CLUBB
       ! does not currently have a predictive equation for <hm'^2>.
       hydromet_zm = zt2zm( hydromet(:,i) )

       do k = 1, gr%nz, 1
          hydrometp2(k,i) = ratio_hmp2_on_hmm2(k,i) * hydromet_zm(k)**2
       enddo ! k = 1, gr%nz, 1

       ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
       ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
       ! A Crank-Nicholson time-stepping scheme is used for this variable.
       ! This is the portion of the calculation using < h_m > from timestep t+1.
       wphydrometp(1:gr%nz-1,i) &
       = wphydrometp(1:gr%nz-1,i) &
         - one_half &
           * xpwp_fnc( K_hm(1:gr%nz-1,i)+nu_hm_vert_res_dep(1:gr%nz-1), &
                       hydromet(1:gr%nz-1,i), hydromet(2:gr%nz,i), &
                       gr%invrs_dzm(1:gr%nz-1) )

       ! A zero-flux boundary condition at the top of the model is used for
       ! hydrometeors.
       wphydrometp(gr%nz,i) = zero

       !!! Calculate the covariance of hydrometeor sedimentation velocity and
       !!! the hydrometeor, which is solved semi-implicitly on thermodynamic
       !!! levels.
       hydromet_vel_covar_zt(:,i) &
       = hydromet_vel_covar_zt_impc(:,i) * hydromet(:,i) &
         + hydromet_vel_covar_zt_expc(:,i)

       ! Boundary conditions for < V_hm'hm' >|_zt.
       hydromet_vel_covar_zt(1,i) = hydromet_vel_covar_zt(2,i)

       !!! Calculate the covariance of hydrometeor sedimentation velocity and
       !!! the hydrometeor, < V_hm'h_m' >, by interpolating the thermodynamic
       !!! level results to momentum levels.
       hydromet_vel_covar(:,i) = zt2zm( hydromet_vel_covar_zt(:,i) )

       ! Boundary conditions for < V_hm'hm' >.
       hydromet_vel_covar(1,i)     = hydromet_vel_covar_zt(2,i)
       hydromet_vel_covar(gr%nz,i) = zero

       ! Statistics for all covariances involving hydrometeors:  < w'h_m' >,
       ! <V_rr'r_r'>, and <V_Nr'N_r'>.
       if ( l_stats_samp ) then

          if ( ihydrometp2(i) > 0 ) then

             ! Covariance of vertical velocity and the hydrometeor.
             call stat_update_var( ihydrometp2(i), hydrometp2(:,i), stats_zm )

          endif

          if ( iwphydrometp(i) > 0 ) then

             ! Covariance of vertical velocity and the hydrometeor.
             call stat_update_var( iwphydrometp(i), wphydrometp(:,i), stats_zm )

          endif

          if ( trim( hydromet_list(i) ) == "rrm" .and. iVrrprrp > 0 ) then

             ! Covariance of sedimentation velocity of r_r and r_r.
             call stat_update_var( iVrrprrp, hydromet_vel_covar(:,iirrm), &
                                   stats_zm )

          elseif ( trim( hydromet_list(i) ) == "Nrm" .and. iVNrpNrp > 0 ) then

             ! Covariance of sedimentation velocity of N_r and N_r.
             call stat_update_var( iVNrpNrp, hydromet_vel_covar(:,iiNrm), stats_zm )

          endif

       endif ! l_stats_samp

       if ( l_stats_samp ) then

          ! Total time tendency
          if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
             call stat_end_update( ixrm_bt, &
                                   hydromet(:,i) &
                                   / dt, stats_zt )
          else
             ! The mean hydrometeor field at thermodynamic level k = 1 is simply
             ! set to the value of the mean hydrometeor field at k = 2.  Don't
             ! include budget stats for level k = 1.
             do k = 2, gr%nz, 1
                call stat_end_update_pt( ixrm_bt, k, &
                                         hydromet(k,i) &
                                         / dt, stats_zt )
             enddo
          endif

       endif ! l_stats_samp

    enddo ! i = 1, hydromet_dim, 1


    return

  end subroutine advance_hydrometeor

  !=============================================================================
  subroutine advance_Ncm( dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm, &
                          rho_ds_zt, invrs_rho_ds_zt, Ncm_mc, &
                          Ncm, Nc_in_cloud, &
                          wpNcp, err_code_Ncm )

    ! Description:
    ! Advance cloud droplet concentration (Ncm) one model time step.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        pi,              & ! Constant(s)
        four_thirds,     &
        one,             &
        one_half,        &
        zero,            &
        rho_lw,          &
        mvr_cloud_max,   &
        cloud_frac_min,  &
        Nc_in_cloud_min, &
        fstderr

    use advance_windm_edsclrm_module, only : &
        xpwp_fnc  ! Procedure(s)

    use parameters_tunable, only: & 
        nu_hm_vert_res_dep  ! Variable(s)

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use error_code, only:  & 
        clubb_at_least_debug_level, & ! Procedure(s)
        clubb_no_error

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_update_var,      & ! Procedure(s)
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_variables, only: & 
        stats_zt,           & ! Variable(s)
        stats_zm,           &
        l_stats_samp, &
        iNcm_bt,      &
        iNcm_mc,      &
        iNcm_cl,      &
        iwpNcp

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Duration of one model time step         [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zt,      & ! mean w wind component on thermodynamic levels  [m/s]
      cloud_frac, & ! Cloud fraction                                 [-]
      K_Nc,       & ! Coefficient of diffusion (turb. adv.) for Nc   [m^2/s]
      rcm           ! Mean cloud water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      Ncm_mc     ! Change in Ncm due to microphysics  [num/kg/s]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Ncm,         & ! Mean cloud droplet conc., <N_c> (thermo. levs.)  [num/kg]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration      [num/kg]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      wpNcp    ! Covariance < w'N_c' > (momentum levels)    [(m/s)(num/kg)]

    integer, intent(out) :: &
      err_code_Ncm    ! Exit code (used to check for errors) for Ncm

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Ncm_vel_covar_zt_impc, & ! Imp. comp. <V_Nc'N_c'> t-levs [m/s]
      Ncm_vel_covar_zt_expc    ! Exp. comp. <V_Nc'N_c'> t-levs [(num/kg)(m/s)]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Ncm_vel,    & ! Mean N_c sedimentation velocity, <V_Nc> (m-levs.) [m/s]
      Ncm_vel_zt    ! Mean N_c sedimenation velocity on thermo. levs.   [m/s]

    real( kind = core_rknd ), dimension(3,gr%nz) :: & 
      lhs    ! Left hand side array

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Ncm_mvr_min, & ! Min. allowable Ncm due to max. cloud droplet mvr [num/kg]
      Ncm_min,     & ! Minimum allowable mean cloud droplet conc.       [num/kg]
      Ncic_min       ! Min. allowable in-cloud mean cloud droplet conc. [num/kg]

    ! Sedimentation velocity of cloud droplet concentration is not considered in
    ! the solution of N_c.
    logical, parameter :: &
      l_Ncm_sed = .false. ! Flag to sediment Ncm with sedimentation vel. Ncm_vel

    integer :: k    ! Loop index


    ! Initialize error code
    err_code_Ncm = clubb_no_error

    ! The mean sedimentation velocity of cloud droplet concentration, < V_Nc >,
    ! and the covariance of N_c sedimentation velocity with N_c, < V_Nc'Nc' >,
    ! are not used in the advancement of < N_c >.  Sedimentation velocity of N_c
    ! is not considered in the solution.  However, these variables need to be
    ! initialized to 0 and passed into microphys_lhs and microphys_rhs.
    Ncm_vel_covar_zt_impc = zero
    Ncm_vel_covar_zt_expc = zero
    Ncm_vel    = zero
    Ncm_vel_zt = zero

    if ( l_stats_samp ) then

       ! Update explicit contributions to cloud droplet concentration.
       call stat_update_var( iNcm_mc, Ncm_mc, stats_zt )

       ! Save prior value of Ncm for determining total time tendency.
       ! The value of Ncm at thermodynamic level k = 1 is simply set to the
       ! value of Ncm at k = 2.  Don't include budget stats for level k = 1.
       do k = 2, gr%nz, 1
          if ( l_in_cloud_Nc_diff ) then
             call stat_begin_update_pt( iNcm_bt, k, &
                                        ( Nc_in_cloud(k) &
                                         * max(cloud_frac(k),cloud_frac_min) ) &
                                        / dt, stats_zt )
          else
             call stat_begin_update_pt( iNcm_bt, k, &
                                        Ncm(k) / dt, &
                                        stats_zt )
          endif
       enddo

    endif ! l_stats_samp

    ! Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
    ! a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
    ! A Crank-Nicholson time-stepping scheme is used for this variable.
    ! This is the portion of the calculation using < N_c > from timestep t. 
    wpNcp(1:gr%nz-1) &
    = - one_half * xpwp_fnc( K_Nc(1:gr%nz-1)+nu_hm_vert_res_dep(1:gr%nz-1), &
                             Ncm(1:gr%nz-1), Ncm(2:gr%nz), &
                             gr%invrs_dzm(1:gr%nz-1) )

    ! A zero-flux boundary condition at the top of the model is used for N_c.
    wpNcp(gr%nz) = zero

    ! Add implicit terms to the LHS array
    call microphys_lhs( "Ncm", l_Ncm_sed, & ! In
                        dt, K_Nc, nu_hm_vert_res_dep, wm_zt, &  ! In
                        Ncm_vel, Ncm_vel_zt, & ! In
                        Ncm_vel_covar_zt_impc, & ! In
                        rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, & ! In
                        lhs ) ! Out

    ! Set up explicit term in the RHS vector
    if ( l_in_cloud_Nc_diff ) then

       call microphys_rhs( "Ncm", dt, l_Ncm_sed, &
                           Nc_in_cloud, &
                           Ncm_mc / max( cloud_frac, cloud_frac_min ), &
                           K_Nc, nu_hm_vert_res_dep, cloud_frac, &
                           Ncm_vel_covar_zt_expc, &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                           rhs )

    else

       call microphys_rhs( "Ncm", dt, l_Ncm_sed, &
                           Ncm, Ncm_mc, &
                           K_Nc, nu_hm_vert_res_dep, cloud_frac, &
                           Ncm_vel_covar_zt_expc, &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                           rhs )

    endif


    !!!!! Advance Ncm one time step.
    if ( l_in_cloud_Nc_diff ) then

       call microphys_solve( "Ncm", l_Ncm_sed, &
                             cloud_frac, &
                             lhs, rhs, Nc_in_cloud, err_code_Ncm )

       Ncm = Nc_in_cloud * max( cloud_frac, cloud_frac_min )

    else

       call microphys_solve( "Ncm", l_Ncm_sed, &
                             cloud_frac, &
                             lhs, rhs, Ncm, err_code_Ncm )

       Nc_in_cloud = Ncm / max( cloud_frac, cloud_frac_min )

    endif


    ! Clipping for mean cloud droplet concentration, <Nc>.

    ! When mean cloud water mixing ratio, <rc>, is found at a grid level, mean
    ! cloud droplet concentration must be at least a minimum value so that
    ! average cloud droplet mean volume radius stays within an upper bound.

    ! The minimum mean rain drop concentration is given by:
    !
    ! <Nc> = <rc> / ( (4/3) * pi * rho_lw * mvr_cloud_max^3 ).

    Ncm_mvr_min &
    = ( one / ( four_thirds * pi * rho_lw * mvr_cloud_max**3 ) ) * rcm

    ! Calculate the minimum threshold value for mean cloud droplet
    ! concentration.  This is the greater of the minimum allowable Ncm based on
    ! the maximum cloud droplet mean volume radius and Nc_in_cloud_min
    ! multiplied by max( cloud_frac, cloud_frac_min ).
    Ncm_min &
    = max( Nc_in_cloud_min * max( cloud_frac, cloud_frac_min ), Ncm_mvr_min )

    ! Print warning message if Ncm (at any level) has a value < Ncm_min.
    if ( any( Ncm < Ncm_min ) ) then

       if ( clubb_at_least_debug_level( 1 ) ) then
          do k = 1, gr%nz
             if ( Ncm(k) < Ncm_min(k) ) then
                write(fstderr,*) "Ncm < ", Ncm_min(k), &
                                 " in advance_microphys at k = ", k
             endif ! Ncm(k) < Ncm_min(k)
          enddo ! k = 1, gr%nz, 1
       endif ! clubb_at_least_debug_level( 1 )

    endif ! Ncm < Ncm_min

    ! Store the previous value of Ncm for the effect of clipping.
    if ( l_stats_samp ) then
       call stat_begin_update( iNcm_cl, &
                               Ncm / dt, stats_zt )
    endif

    if ( any( Ncm < Ncm_min ) ) then

       ! Clip any values of Ncm below the minimum threshold value of Ncm_min to
       ! the minimum threshold value.
       where ( Ncm < Ncm_min )
          Ncm = Ncm_min
       end where

    endif ! Ncm < Ncm_min

    ! Enter the new value of Ncm for the effect of clipping.
    if ( l_stats_samp ) then
       call stat_end_update( iNcm_cl, &
                             Ncm / dt, stats_zt )
    endif

    ! Calculate the minimum threshold value for in-cloud mean cloud droplet
    ! concentration.  This is the greater of the minimum allowable Nc_in_cloud
    ! based on the maximum cloud droplet mean volume radius (Ncm_mvr_min)
    ! divided by max( cloud_frac, cloud_frac_min ) and Nc_in_cloud_min.
    Ncic_min &
    = max( Nc_in_cloud_min, Ncm_mvr_min /  max( cloud_frac, cloud_frac_min ) )

    ! Clip any negative values of Nc_in_cloud to the minimum threshold value of
    ! Ncic_min.
    where ( Nc_in_cloud < Ncic_min )
       Nc_in_cloud = Ncic_min
    end where

    ! Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
    ! a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
    ! A Crank-Nicholson time-stepping scheme is used for this variable.
    ! This is the portion of the calculation using < N_c > from timestep t+1. 
    wpNcp(1:gr%nz-1) &
    = wpNcp(1:gr%nz-1) &
      - one_half * xpwp_fnc( K_Nc(1:gr%nz-1)+nu_hm_vert_res_dep(1:gr%nz-1), &
                             Ncm(1:gr%nz-1), Ncm(2:gr%nz), &
                             gr%invrs_dzm(1:gr%nz-1) )

    ! A zero-flux boundary condition at the top of the model is used for N_c.
    wpNcp(gr%nz) = zero

    ! Statistics for all covariances involving N_c:  < w'N_c' >.
    if ( l_stats_samp ) then

       if ( iwpNcp > 0 ) then

          ! Covariance of vertical velocity and N_c.
          call stat_update_var( iwpNcp, wpNcp, stats_zm )

       endif

    endif ! l_stats_samp

    if ( l_stats_samp ) then

       ! Total time tendency
       ! The value of Ncm at thermodynamic level k = 1 is simply set to the
       ! value of Ncm at k = 2.  Don't include budget stats for level k = 1.
       do k = 2, gr%nz, 1
          call stat_end_update_pt( iNcm_bt, k, &
                                   Ncm(k) / dt, stats_zt )
       enddo ! k = 2, gr%nz, 1

    endif ! l_stats_samp


    return

  end subroutine advance_Ncm

  !=============================================================================
  subroutine microphys_solve( solve_type, l_sed, &
                              cloud_frac, &
                              lhs, rhs, hmm, err_code )

    ! Description:
    ! Solve the tridiagonal system for hydrometeor variable.

    ! References:
    !  None
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only: &
        cloud_frac_min  ! Constant(s)

    use clubb_precision, only:  & 
        core_rknd    ! Variable(s)

    use lapack_wrap, only:  & 
        tridag_solve    ! Procedure(s)

    use error_code, only: &
        clubb_no_error ! Constant

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use stats_variables, only: & 
        stats_zt,  & ! Variable(s)
        irrm_ma, & 
        irrm_sd, & 
        irrm_ts, & 
        irrm_ta, & 
        irim_ma, & 
        irim_sd, & 
        irim_ta, & 
        irsm_ma, & 
        irsm_sd, & 
        irsm_ta, & 
        irgm_ma, & 
        irgm_sd, & 
        irgm_ta, & 
        l_stats_samp, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, &
        ztscr10, &
        ztscr11, &
        ztscr12

    use stats_variables, only: & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_ts, & 
        iNrm_ta, & 
        iNim_ma, & 
        iNim_sd, & 
        iNim_ta, & 
        iNsm_ma, & 
        iNsm_sd, & 
        iNsm_ta, & 
        iNgm_ma, & 
        iNgm_sd, & 
        iNgm_ta

    use stats_variables, only: & 
        iNcm_ma, & 
        iNcm_ta

    use stats_type_utilities, only: &
        stat_update_var_pt, & ! Procedure(s)
        stat_end_update_pt

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      cloud_frac    ! Cloud fraction (thermodynamic levels)        [-]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(3,gr%nz) :: & 
      lhs    ! Left hand side

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
      hmm    ! Mean value of hydrometeor (thermodynamic levels)    [units vary]

    ! Output Variables
    integer, intent(out) :: err_code

    ! Local Variables
    integer :: k, km1, kp1  ! Array indices

    integer :: & 
      ihmm_ma,  & ! Mean advection budget stats toggle
      ihmm_ta,  & ! Turbulent advection budget stats toggle
      ihmm_sd,  & ! Mean sedimentation budget stats toggle
      ihmm_ts     ! Turbulent sedimentation budget stats toggle


    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Initializing ihmm_ma, ihmm_ta, ihmm_sd, ihmm_ts, and in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case( "rrm" )
      ihmm_ma = irrm_ma
      ihmm_ta = irrm_ta
      ihmm_sd = irrm_sd
      ihmm_ts = irrm_ts
    case( "rim" )
      ihmm_ma = irim_ma
      ihmm_ta = irim_ta
      ihmm_sd = irim_sd
      ihmm_ts = 0
    case( "rsm" )
      ihmm_ma = irsm_ma
      ihmm_ta = irsm_ta
      ihmm_sd = irsm_sd
      ihmm_ts = 0
    case( "rgm" )
      ihmm_ma = irgm_ma
      ihmm_ta = irgm_ta
      ihmm_sd = irgm_sd
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ma = iNcm_ma
      ihmm_ta = iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nrm" )
      ihmm_ma = iNrm_ma
      ihmm_ta = iNrm_ta
      ihmm_sd = iNrm_sd
      ihmm_ts = iNrm_ts
    case( "Nim" )
      ihmm_ma = iNim_ma
      ihmm_ta = iNim_ta
      ihmm_sd = iNim_sd
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ma = iNsm_ma
      ihmm_ta = iNsm_ta
      ihmm_sd = iNsm_sd
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ma = iNgm_ma
      ihmm_ta = iNgm_ta
      ihmm_sd = iNgm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select

    ! Solve system using tridag_solve. This uses LAPACK sgtsv,
    ! which relies on Gaussian elimination to decompose the matrix.
    call tridag_solve( solve_type, gr%nz, 1, lhs(1,:), lhs(2,:), lhs(3,:), & 
                       rhs, hmm, err_code )


    ! Statistics
    if ( l_stats_samp ) then

       do k = 1, gr%nz, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nz )

          ! Finalize implicit contributions

          ! hmm term ma is completely implicit; call stat_update_var_pt.
          if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

             ! For Ncm, we divide by cloud_frac when entering the subroutine,
             ! but do not multiply until we return from the subroutine, so we
             ! must account for this here for the budget to balance.
             call stat_update_var_pt( ihmm_ma, k, & 
               ztscr01(k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr02(k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr03(k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                      stats_zt )

          else

             call stat_update_var_pt( ihmm_ma, k, & 
                                      ztscr01(k) * hmm(km1) & 
                                      + ztscr02(k) * hmm(k) & 
                                      + ztscr03(k) * hmm(kp1), stats_zt)

          endif

          ! hmm term sd is completely implicit; call stat_update_var_pt.
          if ( l_sed ) then
             call stat_update_var_pt( ihmm_sd, k, & 
                                      ztscr04(k) * hmm(km1) & 
                                      + ztscr05(k) * hmm(k) & 
                                      + ztscr06(k) * hmm(kp1), stats_zt )
          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( l_sed .and. k > 1 ) then
             call stat_end_update_pt( ihmm_ts, k, & 
                                      ztscr07(k) * hmm(km1) & 
                                      + ztscr08(k) * hmm(k) & 
                                      + ztscr09(k) * hmm(kp1), stats_zt )
          endif

          ! hmm term ta has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( k > 1 ) then

             if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

                ! For Ncm, we divide by cloud_frac when entering the subroutine,
                ! but do not multiply until we return from the subroutine, so we
                ! must account for this here for the budget to balance.
                call stat_end_update_pt( ihmm_ta, k, & 
               ztscr10(k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr11(k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr12(k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                         stats_zt )

             else

                call stat_end_update_pt( ihmm_ta, k, & 
                                         ztscr10(k) * hmm(km1) & 
                                         + ztscr11(k) * hmm(k) & 
                                         + ztscr12(k) * hmm(kp1), stats_zt )

             endif

          endif ! k > 1

       enddo ! 1..gr%nz

    endif ! l_stats_samp


    ! Boundary conditions on results
    !hmm(1) = hmm(2)
    ! Michael Falk, 7 Sep 2007, made this change to eliminate problems
    ! with anomalous rain formation at the top boundary.
    !        hmm(gr%nz) = 0
    !hmm(gr%nz) = hmm(gr%nz-1)
    ! eMFc

    return

  end subroutine microphys_solve

  !=============================================================================
  subroutine microphys_lhs & 
             ( solve_type, l_sed, dt, K_hm, nu, wm_zt, &
               V_hm, V_hmt, &
               Vhmphmp_zt_impc, &
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
               lhs )

    ! Description:
    ! Setup the matrix of implicit contributions to a term.
    ! Can include the effects of sedimentation, diffusion, and advection.
    ! The Morrison microphysics has an explicit sedimentation code, which is
    ! handled elsewhere.
    !
    ! Notes:
    ! Setup for tridiagonal system and boundary conditions should be the same as
    ! the original rain subroutine code.
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm    ! Procedure(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zt_lhs ! Procedure(s)

    use constants_clubb, only: &
        one,         & ! Constant(s)
        one_half,    &
        zero

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use stats_variables, only: & 
        irrm_ma, & ! Variable(s)
        irrm_sd, & 
        irrm_ts, & 
        irrm_ta, & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_ts, & 
        iNrm_ta, & 
        iNcm_ma, &
        iNcm_ta, &
        irim_ma, & 
        irim_sd, & 
        irim_ta, & 
        irsm_ma, & 
        irsm_sd, & 
        irsm_ta, & 
        irgm_ma, & 
        irgm_sd, & 
        irgm_ta, & 
        iNim_ma, & 
        iNim_sd, & 
        iNim_ta, & 
        iNsm_ma, & 
        iNsm_sd, & 
        iNsm_ta, & 
        iNgm_ma, & 
        iNgm_sd, & 
        iNgm_ta

    use stats_variables, only: & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        ztscr10, & 
        ztscr11, & 
        ztscr12, & 
        l_stats_samp

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real( kind = core_rknd ), intent(in) ::  & 
      dt       ! Model timestep                                          [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      nu       ! Background diffusion coefficient                        [m^2/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      wm_zt, & ! w wind component on thermodynamic levels                [m/s]
      V_hm,  & ! Sedimentation velocity of hydrometeor (momentum levels) [m/s]
      V_hmt, & ! Sedimentation velocity of hydrometeor (thermo. levels)  [m/s]
      K_hm     ! Coefficient of diffusion (turb. adv.) for hydrometeor   [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      Vhmphmp_zt_impc, & ! Implicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), intent(out), dimension(3,gr%nz) :: & 
      lhs      ! Left hand side of tridiagonal matrix.

    ! Local Variables
    real( kind = core_rknd ), dimension(3) :: tmp

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Vhmphmp_impc    ! Implicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    integer :: k, km1, kp1  ! Array indices
    integer :: diff_k_in

    integer :: & 
      ihmm_ma, & ! Mean advection budget stats toggle
      ihmm_ta, & ! Turbulent advection budget stats toggle
      ihmm_sd, & ! Mean sedimentation budget stats toggle
      ihmm_ts    ! Turbulent sedimentation budget stats toggle


    ! Initializing ihmm_ma, ihmm_sd, ihmm_ts, and ihmm_ta in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrm" )
      ihmm_ma = irrm_ma
      ihmm_ta = irrm_ta
      ihmm_sd = irrm_sd
      ihmm_ts = irrm_ts
    case ( "Nrm" )
      ihmm_ma = iNrm_ma
      ihmm_ta = iNrm_ta
      ihmm_sd = iNrm_sd
      ihmm_ts = iNrm_ts
    case ( "rim" )
      ihmm_ma = irim_ma
      ihmm_ta = irim_ta
      ihmm_sd = irim_sd
      ihmm_ts = 0
    case ( "rsm" )
      ihmm_ma = irsm_ma
      ihmm_ta = irsm_ta
      ihmm_sd = irsm_sd
      ihmm_ts = 0
    case ( "rgm" )
      ihmm_ma = irgm_ma
      ihmm_ta = irgm_ta
      ihmm_sd = irgm_sd
      ihmm_ts = 0
    case ( "Ncm" )
      ihmm_ma = iNcm_ma
      ihmm_ta = iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ma = iNim_ma
      ihmm_ta = iNim_ta
      ihmm_sd = iNim_sd
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ma = iNsm_ma
      ihmm_ta = iNsm_ta
      ihmm_sd = iNsm_sd
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ma = iNgm_ma
      ihmm_ta = iNgm_ta
      ihmm_sd = iNgm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select


    ! Interpolate the implicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_impc = zt2zm( Vhmphmp_zt_impc )

    ! Reset LHS Matrix for current timestep.
    lhs = zero

    ! Setup LHS Matrix
    do k = 2, gr%nz, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nz )

       ! LHS time tendency.
       lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( one / dt )

       ! LHS turbulent advection term.
       ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
       ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
       !        the turbulent advection term is solved as an eddy-diffusion
       !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
       ! A Crank-Nicholson time-stepping scheme is used for this term.
       if ( k == 2 ) then
          ! The lower boundary condition needs to be applied here at level 2.
          ! The lower boundary condition is a zero-flux boundary condition.
          ! A hydrometeor is not allowed to be fluxed through the model lower
          ! boundary by the processes of mean or turbulent advection.  Only
          ! mean or turbulent sedimentation can flux a hydrometeor through the
          ! lower boundary.  Subroutine diffusion_zt_lhs is set-up to apply a
          ! zero-flux boundary condition at thermodynamic level 1.  In order to
          ! apply the same boundary condition code here at level 2, an adjuster
          ! needs to be used to tell diffusion_zt_lhs to use the code at level 2
          ! that it normally uses at level 1.
          diff_k_in = 1
       else
          diff_k_in = k
       endif
       lhs(kp1_tdiag:km1_tdiag,k) & 
       = lhs(kp1_tdiag:km1_tdiag,k) &
         + one_half &
           * invrs_rho_ds_zt(k) &
           * diffusion_zt_lhs( rho_ds_zm(k) * K_hm(k), &
                               rho_ds_zm(km1) * K_hm(km1), nu, & 
                               gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                               gr%invrs_dzt(k), diff_k_in )

       ! LHS mean advection term.
       lhs(kp1_tdiag:km1_tdiag,k) & 
       = lhs(kp1_tdiag:km1_tdiag,k) & 
         + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
                           gr%invrs_dzm(km1) )

       if ( l_sed ) then

          ! LHS mean sedimentation term.
          ! Note: the Morrison microphysics has its own sedimentation code,
          ! which is applied through the _mc terms to each hydrometeor species.
          ! Therefore, l_sed will always be false when the Morrison microphysics
          ! is enabled.  -dschanen 24 Jan 2011
          if ( .not. l_upwind_diff_sed ) then

             ! Sedimentation (both mean and turbulent) uses centered
             ! differencing.  This is the default method.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_centered_diff_lhs( V_hm(k), V_hm(km1), rho_ds_zm(k), &
                                        rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                        gr%invrs_dzt(k), k )

          else

             ! Sedimentation (both mean and turbulent) uses "upwind"
             ! differencing.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                      rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                      gr%invrs_dzm(k), k )

          endif

          ! LHS turbulent sedimentation term.
          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
             + term_turb_sed_lhs( Vhmphmp_impc(k), Vhmphmp_impc(km1), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(k), rho_ds_zm(km1), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                  invrs_rho_ds_zt(k), k )

       endif ! l_sed

       if ( l_stats_samp ) then

          ! Statistics:  implicit contributions to hydrometeor hmm.

          if ( ihmm_ma > 0 ) then
             tmp(1:3) &
             = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
                               gr%invrs_dzm(km1) )

             ztscr01(k) = -tmp(3)
             ztscr02(k) = -tmp(2)
             ztscr03(k) = -tmp(1)
          endif

          if ( ihmm_sd > 0 .and. l_sed ) then
             if ( .not. l_upwind_diff_sed ) then
                tmp(1:3) &
                = sed_centered_diff_lhs( V_hm(k), V_hm(km1), rho_ds_zm(k), &
                                         rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                         gr%invrs_dzt(k), k )
             else
                tmp(1:3) &
                = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                       rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                       gr%invrs_dzm(k), k )
             endif

             ztscr04(k) = -tmp(3)
             ztscr05(k) = -tmp(2)
             ztscr06(k) = -tmp(1)

          endif

          if ( ihmm_ts > 0 .and. l_sed ) then
             tmp(1:3) &
             = term_turb_sed_lhs( Vhmphmp_impc(k), Vhmphmp_impc(km1), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(k), rho_ds_zm(km1), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                  invrs_rho_ds_zt(k), k )
             ztscr07(k) = -tmp(3)
             ztscr08(k) = -tmp(2)
             ztscr09(k) = -tmp(1)
          endif

          if ( ihmm_ta > 0 ) then
             tmp(1:3) &
             = one_half &
               * invrs_rho_ds_zt(k) & 
               * diffusion_zt_lhs( rho_ds_zm(k) * K_hm(k), &
                                   rho_ds_zm(km1) * K_hm(km1), nu,  & 
                                   gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                   gr%invrs_dzt(k), diff_k_in )
             ztscr10(k) = -tmp(3)
             ztscr11(k) = -tmp(2)
             ztscr12(k) = -tmp(1)
          endif

       endif ! l_stats_samp

    enddo ! 2..gr%nz-1


    ! Boundary Conditions

    ! Lower Boundary
    k = 1

    ! This is set so that < h_m > at thermodynamic level k = 1, which is below
    ! the model lower boundary, is equal to < h_m > at k = 2.
    lhs(k_tdiag,k)   = one
    lhs(kp1_tdiag,k) = -one


    return

  end subroutine microphys_lhs

  !=============================================================================
  subroutine microphys_rhs( solve_type, dt, l_sed, &
                            hmm, hmm_tndcy, &
                            K_hm, nu, cloud_frac, &
                            Vhmphmp_zt_expc, &
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                            rhs )

    ! Description:
    ! Compute RHS vector for a given hydrometeor.
    ! This subroutine computes the explicit portion of the predictive equation
    ! for a given hydrometeor.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one_half,       & ! Constant(s)
        zero,           &
        cloud_frac_min

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff  ! Use in cloud values of Nc for diffusion

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_variables, only: & 
        irrm_ta, & ! Variable(s)
        irrm_ts, &
        iNrm_ta, &
        iNrm_ts, &
        stats_zt, &
        l_stats_samp

    use stats_variables, only: &
        iNim_ta, &
        iNsm_ta, &
        iNgm_ta, &
        iNcm_ta, &
        irim_ta, &
        irsm_ta, &
        irgm_ta

    use stats_type_utilities, only: &
        stat_begin_update_pt ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Duration of model timestep     [s]

    logical, intent(in) :: &
      l_sed    ! Flag for hydrometeor sedimentation

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      hmm,             & ! Mean value of hydrometeor (t-levs.)      [units]
      hmm_tndcy,       & ! Microphysics tendency (thermo. levels)   [units/s]
      K_hm,            & ! Coef. of diffusion for hydrometeor       [m^2/s]
      nu,              & ! Background diffusion coefficient         [m^2/s]
      cloud_frac,      & ! Cloud fraction                           [-]
      Vhmphmp_zt_expc, & ! Explicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      rhs   ! Right hand side

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Vhmphmp_expc    ! Explicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    real( kind = core_rknd ), dimension(3) :: &
      rhs_diff    ! For use in Crank-Nicholson eddy diffusion

    integer :: k, kp1, km1  ! Array indices
    integer :: diff_k_in

    integer :: ihmm_ta, & ! Turbulent advection budget toggle.
               ihmm_ts    ! Turbulent sedimentation budget toggle.
 

    ! Initializing ihmm_ta and ihmm_ts in order to avoid compiler warnings.
    ihmm_ta = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrm" )
      ihmm_ta = irrm_ta
      ihmm_ts = irrm_ts
    case ( "Nrm" )
      ihmm_ta = iNrm_ta
      ihmm_ts = iNrm_ts
    case( "rim" )
      ihmm_ta = irim_ta
      ihmm_ts = 0
    case( "rsm" )
      ihmm_ta = irsm_ta
      ihmm_ts = 0
    case( "rgm" )
      ihmm_ta = irgm_ta
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ta = iNcm_ta
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ta = iNim_ta
      ihmm_ts = 0
    case( "Nsm" )
      ihmm_ta = iNsm_ta
      ihmm_ts = 0
    case( "Ngm" )
      ihmm_ta = iNgm_ta
      ihmm_ts = 0
    case default
      ihmm_ta = 0
      ihmm_ts = 0
    end select


    ! Interpolate the explicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_expc = zt2zm( Vhmphmp_zt_expc )

    ! Initialize right-hand side vector to 0.
    rhs = zero

    ! Hydrometeor right-hand side (explicit portion of the code).
    do k = 2, gr%nz, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nz )

       ! RHS time tendency.
       rhs(k) = hmm(k) / dt

       ! RHS microphysics tendency term (autoconversion, accretion, evaporation,
       ! etc.).
       rhs(k) = rhs(k) + hmm_tndcy(k)

       ! RHS turbulent advection term.
       ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
       ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
       !        the turbulent advection term is solved as an eddy-diffusion
       !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
       ! A Crank-Nicholson time-stepping scheme is used for this term.
       if ( k == 2 ) then
          ! The lower boundary condition needs to be applied here at level 2.
          ! The lower boundary condition is a zero-flux boundary condition.
          ! A hydrometeor is not allowed to be fluxed through the model lower
          ! boundary by the processes of mean or turbulent advection.  Only
          ! mean or turbulent sedimentation can flux a hydrometeor through the
          ! lower boundary.  Subroutine diffusion_zt_lhs is set-up to apply a
          ! zero-flux boundary condition at thermodynamic level 1.  In order to
          ! apply the same boundary condition code here at level 2, an adjuster
          ! needs to be used to tell diffusion_zt_lhs to use the code at level 2
          ! that it normally uses at level 1.
          diff_k_in = 1
       else
          diff_k_in = k
       endif
       rhs_diff(1:3) &
       = one_half &
         * invrs_rho_ds_zt(k) &
         * diffusion_zt_lhs( rho_ds_zm(k) * K_hm(k), &
                             rho_ds_zm(km1) * K_hm(km1), nu, & 
                             gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                             gr%invrs_dzt(k), diff_k_in )

       rhs(k) &
       = rhs(k) &
         - rhs_diff(3) * hmm(km1) &
         - rhs_diff(2) * hmm(k) &
         - rhs_diff(1) * hmm(kp1)

       ! RHS turbulent sedimentation term.
       if ( l_sed ) then
          rhs(k) &
          = rhs(k) &
            + term_turb_sed_rhs( Vhmphmp_expc(k), Vhmphmp_expc(km1), &
                                 Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                 rho_ds_zm(k), rho_ds_zm(km1), &
                                 rho_ds_zt(kp1), rho_ds_zt(k), &
                                 gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                 invrs_rho_ds_zt(k), k )
       endif

       
       if ( l_stats_samp ) then

          ! Statistics: explicit contributions for the hydrometeor.

          ! hmm term ta has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on the right-hand side
          ! turbulent advection component.
          if ( ihmm_ta > 0 ) then

             if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

                ! For Ncm, we divide by cloud_frac when entering the subroutine,
                ! but do not multiply until we return from the subroutine, so we
                ! must account for this here for the budget to balance.
                call stat_begin_update_pt( ihmm_ta, k, & 
              rhs_diff(3) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
              + rhs_diff(2) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
              + rhs_diff(1) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                           stats_zt )

             else

                call stat_begin_update_pt( ihmm_ta, k, & 
                                           rhs_diff(3) * hmm(km1) &
                                           + rhs_diff(2) * hmm(k)   &
                                           + rhs_diff(1) * hmm(kp1), stats_zt )

             endif

          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_update_var_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_turb_sed_rhs.
          if ( ihmm_ts > 0 .and. l_sed ) then
             call stat_begin_update_pt( ihmm_ts, k, &
                 -term_turb_sed_rhs( Vhmphmp_expc(k), Vhmphmp_expc(km1), &
                                     Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                     rho_ds_zm(k), rho_ds_zm(km1), &
                                     rho_ds_zt(kp1), rho_ds_zt(k), &
                                     gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                     invrs_rho_ds_zt(k), k ), &
                                        stats_zt )
          endif ! ihmm_ts > 0 and l_sed

       endif ! l_stats_samp

    enddo ! k = 2, gr%nz, 1
    

    ! Lower boundary conditions on the RHS

    ! This is set so that < h_m > at thermodynamic level k = 1, which is below
    ! the model lower boundary, is equal to < h_m > at k = 2.
    rhs(1) = zero


    return

  end subroutine microphys_rhs

  !=============================================================================
  pure function sed_centered_diff_lhs( V_hm, V_hmm1, rho_ds_zm, &
                                       rho_ds_zmm1, invrs_rho_ds_zt, &
                                       invrs_dzt, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the centered difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the centered-difference
    ! approximation:
    !
    ! The values of <hm> are found on the thermodynamic levels, while the values
    ! of <V_hm> are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The variable <hm>
    ! is interpolated to the intermediate momentum levels.  At the intermediate
    ! momentum levels, the interpolated values of <hm> are multiplied by the
    ! values of <V_hm> and the values of rho_ds_zm.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken over the central thermodynamic level, where
    ! it is multiplied by invrs_rho_ds_zt.
    !
    ! -----hmp1------------------------------------------------ t(k+1)
    !
    ! =============hm(interp)=====V_hm=====rho_ds_zm=========== m(k)
    !
    ! -----hm--------invrs_rho_ds_zt----d(rho_ds*V_hm*hm)/dz--- t(k)
    !
    ! =============hm(interp)=====V_hmm1===rho_ds_zmm1========= m(k-1)
    !
    ! -----hmm1------------------------------------------------ t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1),
    ! respectively.  The letter "t" is used for thermodynamic levels and the
    ! letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level) of the domain,
    ! -rho_ds_zm(1) * V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) ), where the factor in
    ! parentheses is the interpolated value of hm at the zm(1) level.
    ! Furthermore, most of the individual column totals should sum to 0, but the
    ! 1st and 2nd (from the left) columns should combine to sum to the flux out
    ! the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zm(1) * V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) )
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzt )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzt everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zm(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first four vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 |           0                     0                       0
    !    |
    !k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
    !    |    *rho_ds_zm(k-1)      *[ rho_ds_zm(k)        *rho_ds_zm(k)
    !    |    *V_hm(k-1)*D(k)         *V_hm(k)*B(k)       *V_hm(k)*A(k)
    !    |                           -rho_ds_zm(k-1)
    !    |                            *V_hm(k-1)*C(k) ]
    !    |
    !k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzt(k)          *invrs_dzt(k)
    !    |                         *rho_ds_zm(k-1)        *[ rho_ds_zm(k)
    !    |                         *V_hm(k-1)*D(k)           *V_hm(k)*B(k)
    !    |                                                  -rho_ds_zm(k-1)
    !    |                                                   *V_hm(k-1)*C(k) ]
    !    |
    !k=4 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzt(k)
    !    |                                                *rho_ds_zm(k-1)
    !    |                                                *V_hm(k-1)*D(k)
    !    |
    !   \ /
    !
    ! The variables A(k), B(k), C(k), and D(k) are weights of interpolation
    ! around the central thermodynamic level (k), such that:
    !
    ! A(k) = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ),
    ! B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !      = 1 - A(k);
    ! C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
    ! D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !      = 1 - C(k).
    !
    ! Furthermore, for all intermediate thermodynamic grid levels (as long as
    ! k /= gr%nz and k /= 1), the four weighting factors have the following
    ! relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    !
    ! Note:  The superdiagonal term from level 3 and both the main diagonal
    !        and superdiagonal terms from level 4 are not shown on this
    !        diagram.

    ! References:
    ! None

    ! Notes:
    !   Both COAMPS Microphysics and Brian Griffin's implementation use
    !   Khairoutdinov and Kogan (2000) for the calculation of rain
    !   mixing ratio and rain droplet number concentration sedimentation
    !   velocities, but COAMPS has only the local parameterization.
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hm,            & ! Sedimentation velocity of hydrometeor (k)    [m/s]
      V_hmm1,          & ! Sedimentation velocity of hydrometeor (k-1)  [m/s]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)    [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k-1)  [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                  [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,   & ! Momentum level directly above central thermodynamic level.
      mkm1    ! Momentum level directly below central thermodynamic level.

    ! ---- Begin Code ----

    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk   = level
    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    if ( level == 1 ) then

       ! k = 1 (bottom level); lower boundary level; no effects.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    elseif ( level > 1 .and. level < gr%nz ) then

       ! Most of the interior model; normal conditions.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag)  & 
       = + invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zm * V_hm * gr%weights_zt2zm(t_above,mk)

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)  & 
       = + invrs_rho_ds_zt &
           * invrs_dzt &
           * ( rho_ds_zm * V_hm * gr%weights_zt2zm(t_below,mk) & 
               - rho_ds_zmm1 * V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag)  & 
       = - invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zmm1 * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)


    elseif ( level == gr%nz ) then

       ! k = gr%nz (top level); upper boundary level; no flux.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    endif


    return

  end function sed_centered_diff_lhs

  !=============================================================================
  pure function sed_upwind_diff_lhs( V_hmt, V_hmtp1, rho_ds_zt, &
                                     rho_ds_ztp1, invrs_rho_ds_zt, &
                                     invrs_dzm, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the "upwind" difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the upwind-difference
    ! approximation:
    !
    ! The values of <hm> and the values of V_hmt are found on the thermodynamic
    ! levels.  Additionally, the values of rho_ds_zt and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of <hm> are multiplied by the values of
    ! V_hmt and the values of rho_ds_zt.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken between the thermodynamic level above the
    ! central thermodynamic level and the central thermodynamic level.  The
    ! derivative is multiplied by invrs_rho_ds_zt.
    !
    ! --hmp1--V_hmtp1--rho_ds_ztp1--------------------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! --hm----V_hmt----rho_ds_zt--invrs_rho_ds_zt--d(rho_ds*V_hm*hm)/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level, approximated by the zt(1) level for the
    ! "upwind" sedimentation option) of the domain,
    ! -rho_ds_zt(1) * V_hmt(1) * hm(1).  Furthermore, most of the individual
    ! column totals should sum to 0, but the 2nd (from the left) column should
    ! be equal to the flux out the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zt(1) * V_hmt(1) * hm(1)
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzm )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzm everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zt(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first three vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
    !    |  *invrs_dzm(k)          *invrs_dzm(k)
    !    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |  *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzm(k)          *invrs_dzm(k)
    !    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |                         *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=3 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzm(k)
    !    |                                                *rho_ds_zt(k)
    !    |                                                *V_hmt(k)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 3 is not shown on this diagram.

    ! References:
    ! None

    ! Notes:
    ! Both COAMPS Microphysics and Brian Griffin's implementation use
    ! Khairoutdinov and Kogan (2000) for the calculation of rain
    ! mixing ratio and rain droplet number concentration sedimentation
    ! velocities, but COAMPS has only the local parameterization.
    !
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hmtp1,         & ! Sed. velocity of hydrometeor at t-lev (k+1)  [m/s]
      V_hmt,           & ! Sed. velocity of hydrometeor at t-lev (k)    [m/s]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzm          ! Inverse of grid spacing over m-lev. (k)      [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! ---- Begin Code ----

    ! Sedimention is always a downward process, so we omit the upward case
    ! (i.e. the V_hmt variable will always be negative).
    if ( level == gr%nz ) then

       ! k = gr%nz (top level); upper boundary level; no flux.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    else

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = + invrs_rho_ds_zt * invrs_dzm * rho_ds_ztp1 * V_hmtp1

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = - invrs_rho_ds_zt * invrs_dzm * rho_ds_zt * V_hmt

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    endif


    return

  end function sed_upwind_diff_lhs

  !=============================================================================
  pure function term_turb_sed_lhs( Vhmphmp_impcm, Vhmphmp_impcm1, &
                                   Vhmphmp_zt_impcp1, Vhmphmp_zt_impc, &
                                   rho_ds_zm, rho_ds_zmm1, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzm, &
                                   invrs_rho_ds_zt, level ) &
    result( lhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  implicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The implicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < hm >(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the d[ ] / dz term is
    !        changed to a "+".
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    ! 
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variables < hm > and Vhmphmp_zt_impc are both
    ! interpolated to the intermediate momentum levels.  At the momentum levels,
    ! the values of interpolated < hm > and interpolated Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zm.  The mathematical expression F is the product of these three
    ! variables at momentum levels.  Then, the derivative dF/dz is taken over
    ! the central thermodynamic level, where it is multiplied by
    ! invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zm * Vhmphmp_impc(interp) * hmm(interp).
    !
    !
    ! ----hmmp1--------Vhmphmp_zt_impcp1--------------------------------- t(k+1)
    !
    ! =====hmm(interp)=====Vhmphmp_impc(interp)=======rho_ds_zm========== m(k)
    !
    ! ----hmm----------Vhmphmp_zt_impc-----invrs_rho_ds_zt-----dF/dz----- t(k)
    !
    ! =====hmm(interp)=====Vhmphmp_impcm1(interp)=====rho_ds_zmm1======== m(k-1)
    !
    ! ----hmmm1--------Vhmphmp_zt_impcm1--------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of < hm > and Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zt.  The mathematical expression F is the product of these three
    ! variables at thermodynamic levels.  Then, the derivative dF/dz is taken
    ! between the thermodynamic level above the central thermodynamic level and
    ! the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_impc * hmm.
    !
    !
    ! --hmmp1---Vhmphmp_zt_impcp1---rho_ds_ztp1-------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! --hmm-----Vhmphmp_zt_impc-----rho_ds_zt---invrs_rho_ds_zt---dF/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) ).

    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_impcm,     & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_impcm1,    & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k-1) [vary]
      Vhmphmp_zt_impcp1, & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_impc,   & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_zmm1,       & ! Dry, static density at moment. lev (k-1)   [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzm,         & ! Inverse of grid spacing over m-levs. (k)   [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,   & ! Momentum level directly above central thermodynamic level.
      mkm1    ! Momentum level directly below central thermodynamic level.


    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk   = level
    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1


    ! LHS (implicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Implicit component:
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < h_m > ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zm * Vhmphmp_impcm * gr%weights_zt2zm(t_above,mk)

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * ( rho_ds_zm * Vhmphmp_impcm &
                            * gr%weights_zt2zm(t_below,mk) & 
                  - rho_ds_zmm1 * Vhmphmp_impcm1 &
                                * gr%weights_zt2zm(t_above,mkm1)  )

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag)  & 
          = - invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zmm1 * Vhmphmp_impcm1 * gr%weights_zt2zm(t_below,mkm1)


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       endif  ! level


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) &
          = + invrs_rho_ds_zt &
              * invrs_dzm * rho_ds_ztp1 * Vhmphmp_zt_impcp1

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzm * rho_ds_zt * Vhmphmp_zt_impc

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       endif  ! level


    endif


    return

  end function term_turb_sed_lhs

  !=============================================================================
  pure function term_turb_sed_rhs( Vhmphmp_expcm, Vhmphmp_expcm1, &
                                   Vhmphmp_zt_expcp1, Vhmphmp_zt_expc, &
                                   rho_ds_zm, rho_ds_zmm1, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzm, &
                                   invrs_rho_ds_zt, level ) &
    result( rhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  explicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The explicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variable Vhmphmp_zt_expc is interpolated to the
    ! intermediate momentum levels.  At the momentum levels, the values of
    ! interpolated Vhmphmp_zt_expc are multiplied by the values of rho_ds_zm.
    ! Then, the derivative d(rho_ds*Vhmphmp_zt_expc)/dz is taken over the
    ! central thermodynamic level, where it is multiplied by invrs_rho_ds_zt.
    !
    ! ---Vhmphmp_zt_expcp1----------------------------------------------- t(k+1)
    !
    ! ======Vhmphmp_expc(interp)=========rho_ds_zm======================= m(k)
    !
    ! ---Vhmphmp_zt_expc--invrs_rho_ds_zt--d(rho_ds*Vhmphmp_zt_expc)/dz-- t(k)
    !
    ! ======Vhmphmp_expcm1(interp)=======rho_ds_zmm1===================== m(k-1)
    !
    ! ---Vhmphmp_zt_expcm1----------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of Vhmphmp_zt_expc are multiplied by the
    ! values of rho_ds_zt.  The mathematical expression F is the product of
    ! these variables at thermodynamic levels.  Then, the derivative dF/dz is
    ! taken between the thermodynamic level above the central thermodynamic
    ! level and the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_expc.
    !
    !
    ! -----Vhmphmp_zt_expcp1---rho_ds_ztp1------------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! -----Vhmphmp_zt_expc-----rho_ds_zt----invrs_rho_ds_zt----dF/dz----- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) ).
    
    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use parameters_microphys, only: &
        l_upwind_diff_sed  ! Use "upwind" differencing approx. for sedimentation

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_expcm,     & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_expcm1,    & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k-1) [vary]
      Vhmphmp_zt_expcp1, & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_expc,   & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_zmm1,       & ! Dry, static density at moment. lev (k-1)   [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzm,         & ! Inverse of grid spacing over m-levs. (k)   [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    ! Initialize Right-hand side term.
    rhs = zero

    ! RHS (explicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Explicit component:  - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.
          rhs = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzt * ( rho_ds_zm * Vhmphmp_expcm &
                              - rho_ds_zmm1 * Vhmphmp_expcm1 )


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.
          rhs = zero


       endif


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.
          rhs = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzm * ( rho_ds_ztp1 * Vhmphmp_zt_expcp1 &
                              - rho_ds_zt * Vhmphmp_zt_expc )


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.
          rhs = zero


       endif


    endif


    return

  end function term_turb_sed_rhs

!===============================================================================

end module advance_microphys_module
