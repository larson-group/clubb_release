module rev_direction_grid_test

  implicit none

  contains

  !=============================================================================
  function rev_direction_grid_unit_test( )

    use grid_class, only: &
        flip,         & ! Procedure(s)
        cleanup_grid_api, &
        zt2zm2zt,     &
        zm2zt2zm,     &
        ddzm,         &
        ddzt

    use constants_clubb, only: &
        one,     & ! Constant(s)
        zero,    &
        fstdout

    use clubb_api_module, only: &
        grid,               & ! Type(s)
        setup_grid_api,     & ! Procedure(s)
        zt2zm_api,          &
        zm2zt_api,          &
        calc_derrived_params_api,   &
        core_rknd             ! Variable(s)

    use mean_adv, only: &
        term_ma_zt_lhs, & ! Procedure(s)
        term_ma_zm_lhs

    use diffusion, only: &
        diffusion_zt_lhs, & ! Procedure(s)
        diffusion_zm_lhs

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs,         & ! Procedure(s)
        xpyp_term_ta_pdf_lhs_godunov, &
        xpyp_term_ta_pdf_rhs,         &
        xpyp_term_ta_pdf_rhs_godunov

    use err_info_type_module, only: &
      err_info_type,                & ! Type
      init_default_err_info_api,    & ! Procedures
      cleanup_err_info_api

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use parameters_tunable, only: &
        init_clubb_params_api, & ! Procedure(s)
        nu_vertical_res_dep      ! Type(s)

    implicit none

    integer, parameter :: &
      t_above = 1, & ! Upper thermodynamic level index (gr%weights_zt2zm).
      t_below = 2, & ! Lower thermodynamic level index (gr%weights_zt2zm).
      m_above = 1, & ! Upper momentum level index (gr%weights_zm2zt).
      m_below = 2    ! Lower momentum level index (gr%weights_zm2zt).
  
    integer, parameter :: &
      nzm = 10, &
      nzt = nzm-1

    logical, parameter :: &
      l_implemented = .false.

    logical :: l_ascending_grid

    integer :: grid_type

    real( kind = core_rknd ) ::  &
      deltaz,        & ! Vertical grid spacing                  [m]
      zm_init,       & ! Initial grid altitude (momentum level) [m]
      zm_top,        & ! Top of the grid (momentum level)       [m]
      sfc_elevation    ! Surface elevation                      [m]

    real( kind = core_rknd ), dimension(nzm) ::  & 
      momentum_heights    ! Momentum level altitudes                [m]

    real( kind = core_rknd ), dimension(nzt) ::  & 
      thermodynamic_heights    ! Thermodynamic level altitudes      [m]

    real( kind = core_rknd ), dimension(1,nzm) ::  &
      var_zm,                & ! Variable defined on momentum levels [-]
      var_zm_flip,           & ! Variable defined on momentum levels [-]
      var_zm_interp_ascend,  & ! Variable defined on momentum levels [-]
      var_zm_interp_descend, & ! Variable defined on momentum levels [-]
      var_zm_smooth_ascend,  & ! Variable defined on momentum levels [-]
      var_zm_smooth_descend, & ! Variable defined on momentum levels [-]
      var_zm_deriv_ascend,   & ! Variable defined on momentum levels [-]
      var_zm_deriv_descend,  & ! Variable defined on momentum levels [-]
      wm_zm,                 & ! Vertical velocity on m levels       [m/s]
      wm_zm_flip,            & ! Vertical velocity on m levels       [m/s]
      K_zm,                  & ! Eddy diffusivity coef on m levels   [m^s/2]
      K_zm_flip,             & ! Eddy diffusivity coef on m levels   [m^2/s]
      rho_ds_zm,             & ! Base-state density of m levels      [kg/m^3]
      rho_ds_zm_flip,        & ! Base-state density of m levels      [kg/m^3]
      invrs_rho_ds_zm,       & ! Inverse base density of m levs      [m^3/kg]
      invrs_rho_ds_zm_flip     ! Inverse base density of m levs      [m^3/kg]

    real( kind = core_rknd ), dimension(1,nzm) ::  &
      coef_wpxpyp_implicit_zm,      & ! Implicit coef of wpxpyp term [m/s]
      term_wpxpyp_explicit_zm,      & ! Implicit coef of wpxpyp term [m/s] 
      coef_wpxpyp_implicit_zm_flip, & ! Explicit coef of wpxpyp term [units vary]
      term_wpxpyp_explicit_zm_flip, & ! Explicit coef of wpxpyp term [units vary]
      xpyp,                         & ! xpyp                         [units vary]
      sgn_turbulent_vel,            & ! Sign of turbulent velocity   [-]
      sgn_turbulent_vel_flip          ! Sign of turbulent velocity   [-]

    real( kind = core_rknd ), dimension(1,nzt) ::  & 
      var_zt,                & ! Variable defined on thermo. levels [-]
      var_zt_flip,           & ! Variable defined on thermo. levels [-]
      var_zt_interp_ascend,  & ! Variable defined on thermo. levels [-]
      var_zt_interp_descend, & ! Variable defined on thermo. levels [-]
      var_zt_smooth_ascend,  & ! Variable defined on thermo. levels [-]
      var_zt_smooth_descend, & ! Variable defined on thermo. levels [-]
      var_zt_deriv_ascend,   & ! Variable defined on thermo. levels [-]
      var_zt_deriv_descend,  & ! Variable defined on thermo. levels [-]
      wm_zt,                 & ! Vertical velocity on t levels      [m/s]
      wm_zt_flip,            & ! Vertical velocity on t levels      [m/s]
      K_zt,                  & ! Eddy diffusivity coef on t levels  [m^s/2]
      K_zt_flip,             & ! Eddy diffusivity coef on t levels  [m^2/s]
      rho_ds_zt,             & ! Base-state density of t levels     [kg/m^3]
      rho_ds_zt_flip,        & ! Base-state density of t levels     [kg/m^3]
      invrs_rho_ds_zt,       & ! Inverse base density of t levs     [m^3/kg]
      invrs_rho_ds_zt_flip     ! Inverse base density of t levs     [m^3/kg]

    real( kind = core_rknd ), dimension(1,nzt) ::  & 
      coef_wpxpyp_implicit,      & ! Implicit coef of wpxpyp term   [m/s]
      term_wpxpyp_explicit,      & ! Implicit coef of wpxpyp term   [m/s]
      coef_wpxpyp_implicit_flip, & ! Explicit coef of wpxpyp term   [units vary]
      term_wpxpyp_explicit_flip, & ! Explicit coef of wpxpyp term   [units vary]
      sgn_turbulent_vel_zt,      & ! Sign of turbulent velocity     [-]
      sgn_turbulent_vel_zt_flip    ! Sign of turbulent velocity     [-]

    real( kind = core_rknd ), dimension(1) ::  &
      nu    ! Constant background coef of eddy diffusivity    [m^2/s]

    type (grid), target :: &
      gr_ascending,  & ! Grid type variable for ascending grid
      gr_descending    ! Grid type variable for descending grid

    integer :: &
      seed_size

    integer, dimension(:), allocatable :: &
      seed_output

    real( kind = core_rknd ), dimension(nzt+2) :: &
      thermo_hghts_rand_ext

    real( kind = core_rknd ), dimension(nzm) :: &
      momentum_hghts_rand_ext

    real( kind = core_rknd ), dimension(3,1,nzt) ::  & 
      lhs_ma_zt_center_ascend,  & ! lhs zt mean advection array (ascending)
      lhs_ma_zt_center_descend, & ! lhs zt mean advection array (descending)
      lhs_ma_zt_upwind_ascend,  & ! lhs zt mean advection array (ascending)
      lhs_ma_zt_upwind_descend, & ! lhs zt mean advection array (descending)
      lhs_diff_zt_ascend,       & ! lhs zt diffusion array (ascending)
      lhs_diff_zt_descend         ! lhs zt diffusion array (descending)

    real( kind = core_rknd ), dimension(3,1,nzm) ::  & 
      lhs_ma_zm_ascend,       & ! lhs zm mean advection array (ascending)
      lhs_ma_zm_descend,      & ! lhs zm mean advection array (descending)
      lhs_diff_zm_ascend,     & ! lhs zm diffusion array (ascending)
      lhs_diff_zm_descend,    & ! lhs zm diffusion array (descending)
      lhs_ta_center_ascend,   & ! lhs zm turbulent advection array (ascending)
      lhs_ta_center_descend,  & ! lhs zm turbulent advection array (ascending)
      lhs_ta_upwind_ascend,   & ! lhs zm turbulent advection array (ascending)
      lhs_ta_upwind_descend,  & ! lhs zm turbulent advection array (ascending)
      lhs_ta_godunov_ascend,  & ! lhs zm turbulent advection array (ascending)
      lhs_ta_godunov_descend    ! lhs zm turbulent advection array (ascending)

    real( kind = core_rknd ), dimension(1,nzm) ::  & 
      rhs_ta_center_ascend,   & ! rhs zm turbulent advection array (ascending)
      rhs_ta_center_descend,  & ! rhs zm turbulent advection array (ascending)
      rhs_ta_upwind_ascend,   & ! rhs zm turbulent advection array (ascending)
      rhs_ta_upwind_descend,  & ! rhs zm turbulent advection array (ascending)
      rhs_ta_godunov_ascend,  & ! rhs zm turbulent advection array (ascending)
      rhs_ta_godunov_descend    ! rhs zm turbulent advection array (ascending)

    logical :: &
      l_upwind_xm_ma,              & ! Flag for upwind (zt) mean advection
      l_upwind_xpyp_turbulent_adv    ! Flag for upwind (zm pdf) turbulent advection

    integer, parameter :: iunit = 10

    character(len=13), parameter :: &
      namelist_filename = ""

    real( kind = core_rknd ), dimension(1,nparams) :: & 
      clubb_params  ! Array of the model constants

    ! Flag for using prescribed avg_deltaz in calc_derrived_params_api
    logical, parameter :: &
      l_prescribed_avg_deltaz = .false.

    type(nu_vertical_res_dep) :: &
      nu_vert_res_dep_ascend,  & ! Vertical resolution dependent nu (ascending)
      nu_vert_res_dep_descend    ! Vertical resolution dependent nu (descending)
  
    real( kind = core_rknd ) :: & 
      mixt_frac_max_mag, &
      lmin

    integer :: &
      rev_direction_grid_unit_test

    type(err_info_type) :: &
      err_info_dummy        ! err_info struct containing err_code and err_header

    integer :: k

  !--------------------------------- Begin Code ---------------------------------

    ! Init err_info with dummy values.
    ! Latitude and Longitude are set to zero
    call init_default_err_info_api(1, err_info_dummy)

    write(fstdout,*) ""
    write(fstdout,*) "Performing reverse direction grid unit test"
    write(fstdout,*) "==========================================="
    write(fstdout,*) ""

    call random_seed(size=seed_size)
    allocate(seed_output(seed_size))
    !seed used =  -1102123834  1480475276  -248259584   525496705  1614213535 -1267052688     6290033   -78466652
    !seed_output = (/ -1102123834,  1480475276,  -248259584,   525496705,  1614213535, -1267052688,     6290033,   -78466652 /)
    !call random_seed(put=seed_output)
    call random_seed(get=seed_output)

    ! Initialize the output variable to 0.
    ! This value will remain at 0 if all tests are successful.
    rev_direction_grid_unit_test = 0

    ! Read in model parameter values for the call to calc_derrived_params_api
    call init_clubb_params_api( 1, iunit, namelist_filename, &
                                clubb_params )
      
    ! Loop over each grid type
    do grid_type = 1, 3, 1

       if ( grid_type == 1 ) then

          write(fstdout,*) ""
          write(fstdout,*) "Grid type 1 -- evenly-spaced grid"
          write(fstdout,*) "---------------------------------"
          write(fstdout,*) ""

          ! Evenly-spaced grid
          zm_init = 1000.0_core_rknd
          zm_top = 1900.0_core_rknd
          deltaz = 100.0_core_rknd
          momentum_heights = 0.0_core_rknd
          thermodynamic_heights = 0.0_core_rknd

       elseif ( grid_type == 2 ) then

          write(fstdout,*) ""
          write(fstdout,*) "Grid type 2 -- declared on thermodynamic levels"
          write(fstdout,*) "-----------------------------------------------"
          write(fstdout,*) ""

          ! Thermodynamic-level grid
          call random_number( thermo_hghts_rand_ext )
          deltaz = 0.0_core_rknd
          momentum_heights = 0.0_core_rknd
          ! Declare grid on thermodynamic grid levels
          zm_init = 1000.0_core_rknd + thermo_hghts_rand_ext(1)
          zm_top = 3500.0_core_rknd + thermo_hghts_rand_ext(11)
          thermodynamic_heights(1) = 1020.0_core_rknd + thermo_hghts_rand_ext(2)
          thermodynamic_heights(2) = 1075.0_core_rknd + thermo_hghts_rand_ext(3)
          thermodynamic_heights(3) = 1150.0_core_rknd + thermo_hghts_rand_ext(4)
          thermodynamic_heights(4) = 1250.0_core_rknd + thermo_hghts_rand_ext(5)
          thermodynamic_heights(5) = 1500.0_core_rknd + thermo_hghts_rand_ext(6)
          thermodynamic_heights(6) = 1750.0_core_rknd + thermo_hghts_rand_ext(7)
          thermodynamic_heights(7) = 2000.0_core_rknd + thermo_hghts_rand_ext(9)
          thermodynamic_heights(8) = 2500.0_core_rknd + thermo_hghts_rand_ext(9)
          thermodynamic_heights(9) = 3000.0_core_rknd + thermo_hghts_rand_ext(10)

       else ! grid_type = 3

          write(fstdout,*) ""
          write(fstdout,*) "Grid type 3 -- declared on momentum levels"
          write(fstdout,*) "------------------------------------------"
          write(fstdout,*) ""

          ! Momentum-level grid
          call random_number( momentum_hghts_rand_ext )
          deltaz = 0.0_core_rknd
          thermodynamic_heights = 0.0_core_rknd
          ! Declare grid on momentum grid levels
          zm_init = 1000.0_core_rknd + momentum_hghts_rand_ext(1)
          zm_top = 4500.0_core_rknd + momentum_hghts_rand_ext(10)
          momentum_heights(1)  = zm_init
          momentum_heights(2)  = 1050.0_core_rknd + momentum_hghts_rand_ext(2)
          momentum_heights(3)  = 1150.0_core_rknd + momentum_hghts_rand_ext(3)
          momentum_heights(4)  = 1300.0_core_rknd + momentum_hghts_rand_ext(4)
          momentum_heights(5)  = 1500.0_core_rknd + momentum_hghts_rand_ext(5)
          momentum_heights(6)  = 1750.0_core_rknd + momentum_hghts_rand_ext(6)
          momentum_heights(7)  = 2200.0_core_rknd + momentum_hghts_rand_ext(7)
          momentum_heights(8)  = 2750.0_core_rknd + momentum_hghts_rand_ext(8)
          momentum_heights(9)  = 3500.0_core_rknd + momentum_hghts_rand_ext(9)
          momentum_heights(10) = zm_top

       endif

       sfc_elevation = zm_init

       ! Call setup_grid_api for an ascending grid.
       l_ascending_grid = .true.

       call setup_grid_api( nzm, sfc_elevation, &
                            l_implemented, l_ascending_grid, &
                            grid_type, deltaz, zm_init, zm_top, &
                            momentum_heights, thermodynamic_heights, &
                            gr_ascending, err_info_dummy )

       ! Call setup_grid_api for a descending grid.
       l_ascending_grid = .false.
       thermodynamic_heights = flip( thermodynamic_heights, nzt )
       momentum_heights = flip( momentum_heights, nzm )

       call setup_grid_api( nzm, sfc_elevation, &
                            l_implemented, l_ascending_grid, &
                            grid_type, deltaz, zm_init, zm_top, &
                            momentum_heights, thermodynamic_heights, &
                            gr_descending, err_info_dummy )

       !!!!! Grid Test 1 -- Test symmetry for grid variables.

       !!! Momentum level grid variables.
       ! nzm
       if ( abs( gr_ascending%nzm - gr_descending%nzm ) > 0 ) then
          write(fstdout,*) "The variable nzm is not the same between " &
                           // "the ascending and descending grids."
          write(fstdout,*) "Ascending grid nzm = ", gr_ascending%nzm
          write(fstdout,*) "Descending grid nzm = ", gr_descending%nzm
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       do k = 1, nzm
          ! zm
          if ( abs( gr_ascending%zm(1,k) &
                    - gr_descending%zm(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The variable zm is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid zm = ", gr_ascending%zm(1,k)
             write(fstdout,*) "Descending grid zm = ", gr_descending%zm(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! dzm (note: dzm on the ascending grid is defined to be positive
          ! while dzm on the descending grid is defined to be negative.)
          if ( abs( gr_ascending%dzm(1,k) &
                    + gr_descending%dzm(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The variable dzm is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid dzm = ", gr_ascending%dzm(1,k)
             write(fstdout,*) "Descending grid dzm = ", gr_descending%dzm(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! invrs_dzm (note: invrs_dzm on the ascending grid is defined to be positive
          ! while invrs_dzm on the descending grid is defined to be negative.)
          if ( abs( gr_ascending%invrs_dzm(1,k) &
                    + gr_descending%invrs_dzm(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The variable invrs_dzm is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid invrs_dzm = ", gr_ascending%invrs_dzm(1,k)
             write(fstdout,*) "Descending grid invrs_dzm = ", gr_descending%invrs_dzm(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! weights_zt2zm (note: the t_above weight on the ascending grid should equal
          ! the t_below weight on the descending grid.)
          if ( abs( gr_ascending%weights_zt2zm(1,k,t_above) &
                    - gr_descending%weights_zt2zm(1,nzm-k+1,t_below) ) > zero ) then
             write(fstdout,*) "The variable weights_zt2zm is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid weights_zt2zm (t_above) = ", &
                              gr_ascending%weights_zt2zm(1,k,t_above)
             write(fstdout,*) "Descending grid weights_zt2zm (t_below) = ", &
                              gr_descending%weights_zt2zm(1,nzm-k+1,t_below)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! weights_zt2zm (note: the t_below weight on the ascending grid should equal
          ! the t_above weight on the descending grid.)
          if ( abs( gr_ascending%weights_zt2zm(1,k,t_below) &
                    - gr_descending%weights_zt2zm(1,nzm-k+1,t_above) ) > zero ) then
             write(fstdout,*) "The variable weights_zt2zm is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid weights_zt2zm (t_below) = ", &
                              gr_ascending%weights_zt2zm(1,k,t_below)
             write(fstdout,*) "Descending grid weights_zt2zm (t_above) = ", &
                              gr_descending%weights_zt2zm(1,nzm-k+1,t_above)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       !!! Thermodynamic level grid variables.
       ! nzt
       if ( abs( gr_ascending%nzt - gr_descending%nzt ) > 0 ) then
          write(fstdout,*) "The variable nzt is not the same between " &
                          // "the ascending and descending grids."
          write(fstdout,*) "Ascending grid nzt = ", gr_ascending%nzt
          write(fstdout,*) "Descending grid nzt = ", gr_descending%nzt
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       do k = 1, nzt
          ! zt
          if ( abs( gr_ascending%zt(1,k) &
                    - gr_descending%zt(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The variable zt is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid zt = ", gr_ascending%zt(1,k)
             write(fstdout,*) "Descending grid zt = ", gr_descending%zt(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! dzt (note: dzt on the ascending grid is defined to be positive
          ! while dzt on the descending grid is defined to be negative.)
          if ( abs( gr_ascending%dzt(1,k) &
                    + gr_descending%dzt(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The variable dzt is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid dzt = ", gr_ascending%dzt(1,k)
             write(fstdout,*) "Descending grid dzt = ", gr_descending%dzt(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! invrs_dzt (note: invrs_dzt on the ascending grid is defined to be positive
          ! while invrs_dzt on the descending grid is defined to be negative.)
          if ( abs( gr_ascending%invrs_dzt(1,k) &
                    + gr_descending%invrs_dzt(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The variable invrs_dzt is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid invrs_dzt = ", gr_ascending%invrs_dzt(1,k)
             write(fstdout,*) "Descending grid invrs_dzt = ", gr_descending%invrs_dzt(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! weights_zm2zt (note: the m_above weight on the ascending grid should equal
          ! the m_below weight on the descending grid.)
          if ( abs( gr_ascending%weights_zm2zt(1,k,m_above) &
                    - gr_descending%weights_zm2zt(1,nzt-k+1,m_below) ) > zero ) then
             write(fstdout,*) "The variable weights_zm2zt is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid weights_zm2zt (m_above) = ", &
                              gr_ascending%weights_zm2zt(1,k,m_above)
             write(fstdout,*) "Descending grid weights_zm2zt (m_below) = ", &
                              gr_descending%weights_zm2zt(1,nzt-k+1,m_below)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! weights_zm2zt (note: the m_below weight on the ascending grid should equal
          ! the m_above weight on the descending grid.)
          if ( abs( gr_ascending%weights_zm2zt(1,k,m_below) &
                    - gr_descending%weights_zm2zt(1,nzt-k+1,m_above) ) > zero ) then
             write(fstdout,*) "The variable weights_zm2zt is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid weights_zm2zt (m_below) = ", &
                              gr_ascending%weights_zm2zt(1,k,m_below)
             write(fstdout,*) "Descending grid weights_zm2zt (m_above) = ", &
                              gr_descending%weights_zm2zt(1,nzt-k+1,m_above)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       !!!!! Grid Test 2 -- Test symmetry for grid functions.

       ! Set up a variable, var_zm, that is defined on momentum levels
       ! on an ascending grid.
       call random_number( var_zm )
       ! Set up a variable, var_zt, that is defined on thermodynamic levels
       ! on an ascending grid.
       call random_number( var_zt )
       ! Store the flipped array of var_zm to be on momentum levels
       ! on a descending grid.
       var_zm_flip(1,:) = flip( var_zm(1,:), nzm )
       ! Store the flipped array of var_zt to be on thermodynamic levels
       ! on a descending grid.
       var_zt_flip(1,:) = flip( var_zt(1,:), nzt )

       ! Interpolate to thermodynamic levels using zm2zt
       var_zt_interp_ascend = zm2zt_api( nzm, nzt, 1, gr_ascending, var_zm )
       var_zt_interp_descend = zm2zt_api( nzm, nzt, 1, gr_descending, var_zm_flip )

       ! Smooth the thermodynamic-level variables using zt2zm2zt
       var_zt_smooth_ascend = zt2zm2zt( nzm, nzt, 1, gr_ascending, var_zt )
       var_zt_smooth_descend = zt2zm2zt( nzm, nzt, 1, gr_descending, var_zt_flip )

       ! Take the derivative of the momentum-level variable over the
       ! intermediate thermodynamic grid level.
       var_zt_deriv_ascend = ddzm( nzm, nzt, 1, gr_ascending, var_zm )
       var_zt_deriv_descend = ddzm( nzm, nzt, 1, gr_descending, var_zm_flip )

       ! Compare the symmetry of the results
       do k = 1, nzt
          ! zm2zt standard interpolation
          if ( abs( var_zt_interp_ascend(1,k) &
                    - var_zt_interp_descend(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The standard zm2zt interp is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid zm2zt var_zt = ", var_zt_interp_ascend(1,k)
             write(fstdout,*) "Descending grid zm2zt var_zt = ", var_zt_interp_descend(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! zt2zm2zt smoothing
          if ( abs( var_zt_smooth_ascend(1,k) &
                    - var_zt_smooth_descend(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The zt2zm2zt smoothing is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid zt2zm2zt var_zt = ", var_zt_smooth_ascend(1,k)
             write(fstdout,*) "Descending grid zt2zm2zt var_zt = ", var_zt_smooth_descend(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! ddzm derivative
          if ( abs( var_zt_deriv_ascend(1,k) &
                    - var_zt_deriv_descend(1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The ddzm derivative is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid ddzm var_zt = ", var_zt_deriv_ascend(1,k)
             write(fstdout,*) "Descending grid ddzm var_zt = ", var_zt_deriv_descend(1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       ! Interpolate to momentum levels using zt2zm
       var_zm_interp_ascend = zt2zm_api( nzm, nzt, 1, gr_ascending, var_zt )
       var_zm_interp_descend = zt2zm_api( nzm, nzt, 1, gr_descending, var_zt_flip )

       ! Smooth the momentum-level variables using zm2zt2zm
       var_zm_smooth_ascend = zm2zt2zm( nzm, nzt, 1, gr_ascending, var_zm )
       var_zm_smooth_descend = zm2zt2zm( nzm, nzt, 1, gr_descending, var_zm_flip )

       ! Take the derivative of the thermodynamic-level variable over the
       ! intermediate momentum grid level.
       var_zm_deriv_ascend = ddzt( nzm, nzt, 1, gr_ascending, var_zt )
       var_zm_deriv_descend = ddzt( nzm, nzt, 1, gr_descending, var_zt_flip )

       ! Compare the symmetry of the results
       do k = 1, nzm
          ! zt2zm standard interpolation
          if ( abs( var_zm_interp_ascend(1,k) &
                    - var_zm_interp_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The standard zt2zm interp is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid zt2zm var_zm = ", var_zm_interp_ascend(1,k)
             write(fstdout,*) "Descending grid zt2zm var_zm = ", var_zm_interp_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! zm2zt2zm smoothing
          if ( abs( var_zm_smooth_ascend(1,k) &
                    - var_zm_smooth_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The zm2zt2zm smoothing is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid zm2zt2zm var_zm = ", var_zm_interp_ascend(1,k)
             write(fstdout,*) "Descending grid zm2zt2zm var_zm = ", var_zm_interp_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! ddzt derivative
          if ( abs( var_zm_deriv_ascend(1,k) &
                    - var_zm_deriv_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The ddzt derivative is not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid ddzt var_zm = ", var_zm_deriv_ascend(1,k)
             write(fstdout,*) "Descending grid ddzt var_zm = ", var_zm_deriv_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       !!!!! Grid Test 3 -- Test symmetry for some functions or subroutines
       !!!!!                used in CLUBB.

       ! LHS mean advection code
       call random_number( wm_zm(1,:) ) ! wm_zm values between 0 and 1 m/s
       ! wm_zm values beween -1 m/s and 1 m/s
       wm_zm(1,:) = 2.0_core_rknd * ( wm_zm(1,:) - 0.5_core_rknd )
       wm_zt = zm2zt_api( nzm, nzt, 1, gr_ascending, wm_zm )
       wm_zm_flip(1,:) = flip( wm_zm(1,:), nzm )
       wm_zt_flip(1,:) = flip( wm_zt(1,:), nzt )
       ! zt -- upwinding off
       l_upwind_xm_ma = .false.
       call term_ma_zt_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, wm_zt, &
                            gr_ascending%weights_zt2zm, &
                            gr_ascending%invrs_dzt, gr_ascending%invrs_dzm, &
                            l_upwind_xm_ma, gr_ascending%grid_dir, &
                            lhs_ma_zt_center_ascend )
       call term_ma_zt_lhs( gr_descending%nzm, gr_descending%nzt, 1, wm_zt_flip, &
                            gr_descending%weights_zt2zm, &
                            gr_descending%invrs_dzt, gr_descending%invrs_dzm, &
                            l_upwind_xm_ma, gr_descending%grid_dir, &
                            lhs_ma_zt_center_descend )
       ! zt -- upwinding on
       l_upwind_xm_ma = .true.
       call term_ma_zt_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, wm_zt, &
                            gr_ascending%weights_zt2zm, &
                            gr_ascending%invrs_dzt, gr_ascending%invrs_dzm, &
                            l_upwind_xm_ma, gr_ascending%grid_dir, &
                            lhs_ma_zt_upwind_ascend )
       call term_ma_zt_lhs( gr_descending%nzm, gr_descending%nzt, 1, wm_zt_flip, &
                            gr_descending%weights_zt2zm, &
                            gr_descending%invrs_dzt, gr_descending%invrs_dzm, &
                            l_upwind_xm_ma, gr_descending%grid_dir, &
                            lhs_ma_zt_upwind_descend )
       ! zm
       call term_ma_zm_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, wm_zm, &
                            gr_ascending%invrs_dzm, gr_ascending%weights_zm2zt, &
                            lhs_ma_zm_ascend )
       call term_ma_zm_lhs( gr_descending%nzm, gr_descending%nzt, 1, wm_zm_flip, &
                            gr_descending%invrs_dzm, gr_descending%weights_zm2zt, &
                            lhs_ma_zm_descend )

       ! LHS diffusion code
       call random_number( K_zm ) ! K_zm is between 0 and 1
       K_zm = 5.0_core_rknd * K_zm ! K_zm is between 0 and 5
       K_zt = zm2zt_api( nzm, nzt, 1, gr_ascending, K_zm )
       K_zm_flip(1,:) = flip( K_zm(1,:), nzm )
       K_zt_flip(1,:) = flip( K_zt(1,:), nzt )
       call random_number( nu ) ! nu is between 0 and 1
       nu = 20.0_core_rknd * nu ! nu is between 0 and 20
       call random_number( rho_ds_zm ) ! rho_ds_zm is between 0 and 1
       ! Set rho_ds_zm to be between 0.75 and 1.
       rho_ds_zm = min( 0.25_core_rknd * rho_ds_zm + 0.75_core_rknd, 1.0_core_rknd )
       ! Set rho_ds_zm at the surface to 1.2 kg/m^3.
       rho_ds_zm(:,1) = 1.2_core_rknd
       ! Set the value of rho_ds_zm so that its value at each ascending grid
       ! level is less than its value at the grid level below. This will
       ! result in a profile that decreases with altitude but remains positive.
       do k = 2, nzm, 1
          rho_ds_zm(:,k) = rho_ds_zm(:,k) * rho_ds_zm(:,k-1)
       enddo
       rho_ds_zt = zm2zt_api( nzm, nzt, 1, gr_ascending, rho_ds_zm )
       invrs_rho_ds_zm = one / rho_ds_zm
       invrs_rho_ds_zt = one / rho_ds_zt
       rho_ds_zm_flip(1,:) = flip( rho_ds_zm(1,:), nzm )
       rho_ds_zt_flip(1,:) = flip( rho_ds_zt(1,:), nzt )
       invrs_rho_ds_zm_flip(1,:) = flip( invrs_rho_ds_zm(1,:), nzm )
       invrs_rho_ds_zt_flip(1,:) = flip( invrs_rho_ds_zt(1,:), nzt )
       ! zt
       call diffusion_zt_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, &
                              K_zm, K_zt, nu, &
                              invrs_rho_ds_zt, rho_ds_zm, &
                              lhs_diff_zt_ascend ) 
       call diffusion_zt_lhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, &
                              K_zm_flip, K_zt_flip, nu, &
                              invrs_rho_ds_zt_flip, rho_ds_zm_flip, &
                              lhs_diff_zt_descend ) 
       ! zm
       call diffusion_zm_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, &
                              K_zt, K_zm, nu, &
                              invrs_rho_ds_zm, rho_ds_zt, &
                              lhs_diff_zm_ascend )
       call diffusion_zm_lhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, &
                              K_zt_flip, K_zm_flip, nu, &
                              invrs_rho_ds_zm_flip, rho_ds_zt_flip, &
                              lhs_diff_zm_descend )

       ! xpyp turbulent advection (PDF)
       call random_number( coef_wpxpyp_implicit ) ! coef_wpxpyp_implicit between 0 and 1
       ! Set coef_wpxpyp_implicit between -1 and 1 for this test since it can
       ! have a positive or negative value.
       coef_wpxpyp_implicit = 2.0_core_rknd * ( coef_wpxpyp_implicit - 0.5_core_rknd )
       call random_number( term_wpxpyp_explicit ) ! term_wpxpyp_explicit between 0 and 1
       ! Set term_wpxpyp_explicit between -1 and 1 for this test since it can
       ! have a positive or negative value.
       term_wpxpyp_explicit = 2.0_core_rknd * ( term_wpxpyp_explicit - 0.5_core_rknd )
       coef_wpxpyp_implicit_zm = zt2zm_api( nzm, nzt, 1, gr_ascending, coef_wpxpyp_implicit )
       term_wpxpyp_explicit_zm = zt2zm_api( nzm, nzt, 1, gr_ascending, term_wpxpyp_explicit )
       call random_number( xpyp ) ! xpyp between 0 and 1
       ! Let's have some fun with sgn_turbulent_vel (momentum level variable) by making
       ! xpyp range between -1 and 1 too, since it can be positive or negative.
       xpyp = 2.0_core_rknd * ( xpyp - 0.5_core_rknd )
       sgn_turbulent_vel = sign( one, coef_wpxpyp_implicit_zm * xpyp + term_wpxpyp_explicit_zm )
       sgn_turbulent_vel_zt = zm2zt_api( nzm, nzt, 1, gr_ascending, sgn_turbulent_vel )
       coef_wpxpyp_implicit_flip(1,:) = flip( coef_wpxpyp_implicit(1,:), nzt )
       term_wpxpyp_explicit_flip(1,:) = flip( term_wpxpyp_explicit(1,:), nzt )
       coef_wpxpyp_implicit_zm_flip(1,:) = flip( coef_wpxpyp_implicit_zm(1,:), nzm )
       term_wpxpyp_explicit_zm_flip(1,:) = flip( term_wpxpyp_explicit_zm(1,:), nzm )
       sgn_turbulent_vel_flip(1,:) = flip( sgn_turbulent_vel(1,:), nzm )
       sgn_turbulent_vel_zt_flip(1,:) = flip( sgn_turbulent_vel_zt(1,:), nzt )
       ! lhs -- centered differencing
       l_upwind_xpyp_turbulent_adv = .false.
       call xpyp_term_ta_pdf_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                  coef_wpxpyp_implicit, & ! In
                                  rho_ds_zt, rho_ds_zm, & ! In
                                  invrs_rho_ds_zm, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel, & ! In
                                  coef_wpxpyp_implicit_zm, & ! In
                                  lhs_ta_center_ascend ) ! Out
       call xpyp_term_ta_pdf_lhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                  coef_wpxpyp_implicit_flip, & ! In
                                  rho_ds_zt_flip, rho_ds_zm_flip, & ! In
                                  invrs_rho_ds_zm_flip, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel_flip, & ! In
                                  coef_wpxpyp_implicit_zm_flip, & ! In
                                  lhs_ta_center_descend ) ! Out
       ! lhs -- upwind differencing
       l_upwind_xpyp_turbulent_adv = .true.
       call xpyp_term_ta_pdf_lhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                  coef_wpxpyp_implicit, & ! In
                                  rho_ds_zt, rho_ds_zm, & ! In
                                  invrs_rho_ds_zm, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel, & ! In
                                  coef_wpxpyp_implicit_zm, & ! In
                                  lhs_ta_upwind_ascend ) ! Out
       call xpyp_term_ta_pdf_lhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                  coef_wpxpyp_implicit_flip, & ! In
                                  rho_ds_zt_flip, rho_ds_zm_flip, & ! In
                                  invrs_rho_ds_zm_flip, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel_flip, & ! In
                                  coef_wpxpyp_implicit_zm_flip, & ! In
                                  lhs_ta_upwind_descend ) ! Out
       ! lhs -- godunov type
       call xpyp_term_ta_pdf_lhs_godunov( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                          coef_wpxpyp_implicit, & ! In
                                          invrs_rho_ds_zm, rho_ds_zm,  & ! In
                                          lhs_ta_godunov_ascend ) ! Out
       call xpyp_term_ta_pdf_lhs_godunov( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                          coef_wpxpyp_implicit_flip, & ! In
                                          invrs_rho_ds_zm_flip, rho_ds_zm_flip,  & ! In
                                          lhs_ta_godunov_descend ) ! Out
       ! rhs -- centered differencing
       l_upwind_xpyp_turbulent_adv = .false.
       call xpyp_term_ta_pdf_rhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                  term_wpxpyp_explicit, & ! In
                                  rho_ds_zt, rho_ds_zm, & ! In
                                  invrs_rho_ds_zm, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel, & ! In
                                  term_wpxpyp_explicit_zm, & ! In
                                  rhs_ta_center_ascend ) ! Out
       call xpyp_term_ta_pdf_rhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                  term_wpxpyp_explicit_flip, & ! In
                                  rho_ds_zt_flip, rho_ds_zm_flip, & ! In
                                  invrs_rho_ds_zm_flip, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel_flip, & ! In
                                  term_wpxpyp_explicit_zm_flip, & ! In
                                  rhs_ta_center_descend ) ! Out
       ! rhs -- upwind differencing
       l_upwind_xpyp_turbulent_adv = .true.
       call xpyp_term_ta_pdf_rhs( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                  term_wpxpyp_explicit, & ! In
                                  rho_ds_zt, rho_ds_zm, & ! In
                                  invrs_rho_ds_zm, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel, & ! In
                                  term_wpxpyp_explicit_zm, & ! In
                                  rhs_ta_upwind_ascend ) ! Out
       call xpyp_term_ta_pdf_rhs( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                  term_wpxpyp_explicit_flip, & ! In
                                  rho_ds_zt_flip, rho_ds_zm_flip, & ! In
                                  invrs_rho_ds_zm_flip, & ! In
                                  l_upwind_xpyp_turbulent_adv, & ! In
                                  sgn_turbulent_vel_flip, & ! In
                                  term_wpxpyp_explicit_zm_flip, & ! In
                                  rhs_ta_upwind_descend ) ! Out
       ! rhs -- godunov type
       call xpyp_term_ta_pdf_rhs_godunov( gr_ascending%nzm, gr_ascending%nzt, 1, gr_ascending, & ! In
                                          term_wpxpyp_explicit_zm, & ! In
                                          invrs_rho_ds_zm, & ! In
                                          sgn_turbulent_vel_zt, & ! In
                                          rho_ds_zm, & ! In
                                          rhs_ta_godunov_ascend ) ! Out
       call xpyp_term_ta_pdf_rhs_godunov( gr_descending%nzm, gr_descending%nzt, 1, gr_descending, & ! In
                                          term_wpxpyp_explicit_zm_flip, & ! In
                                          invrs_rho_ds_zm_flip, & ! In
                                          sgn_turbulent_vel_zt_flip, & ! In
                                          rho_ds_zm_flip, & ! In
                                          rhs_ta_godunov_descend ) ! Out

       ! calc_derrived_params_api (called from within setup_parameters in CLUBB)
       call calc_derrived_params_api( gr_ascending, grid_type, deltaz, & ! Intent(in)
                              clubb_params(1,:),               & ! Intent(in)
                              l_prescribed_avg_deltaz,         & ! Intent(in)
                              nu_vert_res_dep_ascend, lmin,    & ! intent(inout)
                              mixt_frac_max_mag )                ! intent(inout)
       call calc_derrived_params_api( gr_descending, grid_type, deltaz, & ! Intent(in)
                              clubb_params(1,:),                & ! Intent(in)
                              l_prescribed_avg_deltaz,          & ! Intent(in)
                              nu_vert_res_dep_descend, lmin,    & ! intent(inout)
                              mixt_frac_max_mag )                 ! intent(inout)

       ! Compare the symmetry of the results -- zt level variables
       do k = 1, nzt
          ! term_ma_zt_lhs (centered differencing)
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ma_zt_center_ascend(1,1,k) &
                    - lhs_ma_zt_center_descend(3,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (centered) superdiag = ", &
                              lhs_ma_zt_center_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (centered) subdiag = ", &
                              lhs_ma_zt_center_descend(3,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ma_zt_center_ascend(2,1,k) &
                    - lhs_ma_zt_center_descend(2,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (centered) main diag = ", &
                              lhs_ma_zt_center_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (centered) main diag = ", &
                              lhs_ma_zt_center_descend(2,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ma_zt_center_ascend(3,1,k) &
                    - lhs_ma_zt_center_descend(1,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (centered) subdiag = ", &
                              lhs_ma_zt_center_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (centered) superdiag = ", &
                              lhs_ma_zt_center_descend(1,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! term_ma_zt_lhs (upwind differencing)
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ma_zt_upwind_ascend(1,1,k) &
                    - lhs_ma_zt_upwind_descend(3,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (upwind) superdiag = ", &
                              lhs_ma_zt_upwind_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (upwind) subdiag = ", &
                              lhs_ma_zt_upwind_descend(3,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ma_zt_upwind_ascend(2,1,k) &
                    - lhs_ma_zt_upwind_descend(2,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (upwind) main diag = ", &
                              lhs_ma_zt_upwind_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (upwind) main diag = ", &
                              lhs_ma_zt_upwind_descend(2,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ma_zt_upwind_ascend(3,1,k) &
                    - lhs_ma_zt_upwind_descend(1,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zt (upwind) subdiag = ", &
                              lhs_ma_zt_upwind_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zt (upwind) superdiag = ", &
                              lhs_ma_zt_upwind_descend(1,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! diffusion_zt_lhs
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_diff_zt_ascend(1,1,k) &
                    - lhs_diff_zt_descend(3,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zt superdiag = ", &
                              lhs_diff_zt_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zt subdiag = ", &
                              lhs_diff_zt_descend(3,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_diff_zt_ascend(2,1,k) &
                    - lhs_diff_zt_descend(2,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zt main diag = ", &
                              lhs_diff_zt_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zt main diag = ", &
                              lhs_diff_zt_descend(2,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_diff_zt_ascend(3,1,k) &
                    - lhs_diff_zt_descend(1,1,nzt-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zt_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzt-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zt subdiag = ", &
                              lhs_diff_zt_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zt superdiag = ", &
                              lhs_diff_zt_descend(1,1,nzt-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       ! Compare the symmetry of the results -- zm level variables
       do k = 1, nzm
          ! term_ma_zm_lhs
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ma_zm_ascend(1,1,k) &
                    - lhs_ma_zm_descend(3,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zm superdiag = ", &
                              lhs_ma_zm_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zm subdiag = ", &
                              lhs_ma_zm_descend(3,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ma_zm_ascend(2,1,k) &
                    - lhs_ma_zm_descend(2,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zm main diag = ", &
                              lhs_ma_zm_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zm main diag = ", &
                              lhs_ma_zm_descend(2,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ma_zm_ascend(3,1,k) &
                    - lhs_ma_zm_descend(1,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of term_ma_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ma_zm subdiag = ", &
                              lhs_ma_zm_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ma_zm superdiag = ", &
                              lhs_ma_zm_descend(1,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! diffusion_zm_lhs
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_diff_zm_ascend(1,1,k) &
                    - lhs_diff_zm_descend(3,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zm superdiag = ", &
                              lhs_diff_zm_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zm subdiag = ", &
                              lhs_diff_zm_descend(3,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_diff_zm_ascend(2,1,k) &
                    - lhs_diff_zm_descend(2,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zm main diag = ", &
                              lhs_diff_zm_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zm main diag = ", &
                              lhs_diff_zm_descend(2,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_diff_zm_ascend(3,1,k) &
                    - lhs_diff_zm_descend(1,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of diffusion_zm_lhs are not symmetric between " &
                              // "the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_diff_zm subdiag = ", &
                              lhs_diff_zm_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_diff_zm superdiag = ", &
                              lhs_diff_zm_descend(1,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! xpyp_term_ta_pdf_lhs
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ta_center_ascend(1,1,k) &
                    - lhs_ta_center_descend(3,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_center superdiag = ", &
                              lhs_ta_center_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ta_center subdiag = ", &
                              lhs_ta_center_descend(3,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ta_center_ascend(2,1,k) &
                    - lhs_ta_center_descend(2,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_center main diag = ", &
                              lhs_ta_center_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ta_center main diag = ", &
                              lhs_ta_center_descend(2,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ta_center_ascend(3,1,k) &
                    - lhs_ta_center_descend(1,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_center subdiag = ", &
                              lhs_ta_center_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ta_center superdiag = ", &
                              lhs_ta_center_descend(1,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ta_upwind_ascend(1,1,k) &
                    - lhs_ta_upwind_descend(3,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_upwind superdiag = ", &
                              lhs_ta_upwind_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ta_upwind subdiag = ", &
                              lhs_ta_upwind_descend(3,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ta_upwind_ascend(2,1,k) &
                    - lhs_ta_upwind_descend(2,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_upwind main diag = ", &
                              lhs_ta_upwind_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ta_upwind main diag = ", &
                              lhs_ta_upwind_descend(2,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ta_upwind_ascend(3,1,k) &
                    - lhs_ta_upwind_descend(1,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_upwind subdiag = ", &
                              lhs_ta_upwind_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ta_upwind superdiag = ", &
                              lhs_ta_upwind_descend(1,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! xpyp_term_ta_pdf_lhs_godunov
          ! Note: the superdiagonal term on the ascending grid should exactly
          ! match the subdiagonal term on the descending grid.
          if ( abs( lhs_ta_godunov_ascend(1,1,k) &
                    - lhs_ta_godunov_descend(3,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs_godunov are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_godunov superdiag = ", &
                              lhs_ta_godunov_ascend(1,1,k)
             write(fstdout,*) "Descending grid lhs_ta_godunov subdiag = ", &
                              lhs_ta_godunov_descend(3,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the main diagonal term on the ascending grid should exactly
          ! match the main diagonal term on the descending grid.
          if ( abs( lhs_ta_godunov_ascend(2,1,k) &
                    - lhs_ta_godunov_descend(2,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs_godunov are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_godunov main diag = ", &
                              lhs_ta_godunov_ascend(2,1,k)
             write(fstdout,*) "Descending grid lhs_ta_godunov main diag = ", &
                              lhs_ta_godunov_descend(2,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! Note: the subdiagonal term on the ascending grid should exactly
          ! match the superdiagonal term on the descending grid.
          if ( abs( lhs_ta_godunov_ascend(3,1,k) &
                    - lhs_ta_godunov_descend(1,1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_lhs_godunov are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid lhs_ta_godunov subdiag = ", &
                              lhs_ta_godunov_ascend(3,1,k)
             write(fstdout,*) "Descending grid lhs_ta_godunov superdiag = ", &
                              lhs_ta_godunov_descend(1,1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! xpyp_term_ta_pdf_rhs
          if ( abs( rhs_ta_center_ascend(1,k) &
                    - rhs_ta_center_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_rhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid rhs_ta_center superdiag = ", &
                              rhs_ta_center_ascend(1,k)
             write(fstdout,*) "Descending grid rhs_ta_center subdiag = ", &
                              rhs_ta_center_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          if ( abs( rhs_ta_upwind_ascend(1,k) &
                    - rhs_ta_upwind_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_rhs are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid rhs_ta_upwind superdiag = ", &
                              rhs_ta_upwind_ascend(1,k)
             write(fstdout,*) "Descending grid rhs_ta_upwind subdiag = ", &
                              rhs_ta_upwind_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
          ! xpyp_term_ta_pdf_rhs_godunov
          if ( abs( rhs_ta_godunov_ascend(1,k) &
                    - rhs_ta_godunov_descend(1,nzm-k+1) ) > zero ) then
             write(fstdout,*) "The results of xpyp_term_ta_pdf_rhs_godunov are not " &
                              // "symmetric between the ascending and descending grids."
             write(fstdout,*) "Ascending grid level = ", k
             write(fstdout,*) "Descending grid level = ", nzm-k+1
             write(fstdout,*) "Ascending grid rhs_ta_godunov superdiag = ", &
                              rhs_ta_godunov_ascend(1,k)
             write(fstdout,*) "Descending grid rhs_ta_godunov subdiag = ", &
                              rhs_ta_godunov_descend(1,nzm-k+1)
             rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
          endif
       enddo

       ! Compare results for nu_vert_res_dep between ascending and descending grids.
       ! % nu1
       if ( abs( nu_vert_res_dep_ascend%nu1(1) &
                 - nu_vert_res_dep_descend%nu1(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu1 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu1 = ", &
                           nu_vert_res_dep_ascend%nu1(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu1 = ", &
                           nu_vert_res_dep_descend%nu1(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu2
       if ( abs( nu_vert_res_dep_ascend%nu2(1) &
                 - nu_vert_res_dep_descend%nu2(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu2 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu2 = ", &
                           nu_vert_res_dep_ascend%nu2(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu2 = ", &
                           nu_vert_res_dep_descend%nu2(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu6
       if ( abs( nu_vert_res_dep_ascend%nu6(1) &
                 - nu_vert_res_dep_descend%nu6(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu6 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu6 = ", &
                           nu_vert_res_dep_ascend%nu6(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu6 = ", &
                           nu_vert_res_dep_descend%nu6(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu8
       if ( abs( nu_vert_res_dep_ascend%nu8(1) &
                 - nu_vert_res_dep_descend%nu8(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu8 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu8 = ", &
                           nu_vert_res_dep_ascend%nu8(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu8 = ", &
                           nu_vert_res_dep_descend%nu8(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu9
       if ( abs( nu_vert_res_dep_ascend%nu9(1) &
                 - nu_vert_res_dep_descend%nu9(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu9 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu9 = ", &
                           nu_vert_res_dep_ascend%nu9(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu9 = ", &
                           nu_vert_res_dep_descend%nu9(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu10
       if ( abs( nu_vert_res_dep_ascend%nu10(1) &
                 - nu_vert_res_dep_descend%nu10(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu10 are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu10 = ", &
                           nu_vert_res_dep_ascend%nu10(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu10 = ", &
                           nu_vert_res_dep_descend%nu10(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif
       ! % nu_hm
       if ( abs( nu_vert_res_dep_ascend%nu_hm(1) &
                 - nu_vert_res_dep_descend%nu_hm(1) ) > zero ) then
          write(fstdout,*) "The results of nu_vert_res_dep for nu_hm are not symmetric " &
                           // "between the ascending and descending grids."
          write(fstdout,*) "Ascending grid nu_vert_res_dep nu_hm = ", &
                           nu_vert_res_dep_ascend%nu_hm(1)
          write(fstdout,*) "Descending grid nu_vert_res_dep nu_hm = ", &
                           nu_vert_res_dep_descend%nu_hm(1)
          rev_direction_grid_unit_test = rev_direction_grid_unit_test + 1
       endif

       ! Deallocate gr variables
       ! (They are allocated during the call to setup_grid_api).
       call cleanup_grid_api( gr_ascending )
       call cleanup_grid_api( gr_descending )

    enddo ! iter = 1, 3, 1


    ! Print results and return exit code.
    if ( rev_direction_grid_unit_test == 0 ) then
       ! Exit Code = 0, Success!
       write(fstdout,'(1x,A)') "Success!"
       write(fstdout,*) ""
    else ! rev_direction_grid_unit_test > 0
       ! Exit Code > 0, Fail
       write(fstdout,'(1x,A,I4,A)') "There were ", rev_direction_grid_unit_test, &
                                    " failed test sets."
       write(fstdout,*) "seed used = ", seed_output
       write(fstdout,*) ""
    endif ! rev_direction_grid_unit_test = 0

    ! Deallocate err_info
    call cleanup_err_info_api( err_info_dummy )

    return

  end function rev_direction_grid_unit_test

  !=============================================================================

end module rev_direction_grid_test 
