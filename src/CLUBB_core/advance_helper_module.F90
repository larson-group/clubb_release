
module advance_helper_module

! Description:
!   This module contains helper methods for the advance_* modules.
!------------------------------------------------------------------------

  implicit none

  public :: &
    set_boundary_conditions_lhs, &
    set_boundary_conditions_rhs, &
    calc_stability_correction,   &
    calc_brunt_vaisala_freq_sqd, &
    compute_Cx_fnc_Richardson, &
    wp2_term_splat_lhs, &
    wp3_term_splat_lhs, &
    smooth_min, smooth_max, &
    smooth_heaviside_peskin, &
    calc_xpwp, &
    vertical_avg, &
    vertical_integral, &
    Lscale_width_vert_avg

  interface calc_xpwp
    module procedure calc_xpwp_1D
    module procedure calc_xpwp_2D
  end interface

  private ! Set Default Scope

!===============================================================================
  interface smooth_min

    ! These functions smooth the output of the min function 
    ! by introducing a varyingly steep path between the two input variables.
    ! The degree to which smoothing is applied depends on the value of 'smth_coef'.
    ! If 'smth_coef' goes toward 0, the output of the min function will be 
    !        0.5 * ((a+b) - abs(a-b))
    ! If a > b, then this comes out to be b. Likewise, if a < b, abs(a-b)=b-a so we get a.
    ! Increasing the smoothing coefficient will lead to a greater degree of smoothing
    ! in the smooth min and max functions. Generally, the coefficient should roughly scale
    ! with the magnitude of data in the data structure that is to be smoothed, in order to
    ! obtain a sensible degree of smoothing (not too much, not too little).

    module procedure smooth_min_scalar_array
    module procedure smooth_min_array_scalar
    module procedure smooth_min_arrays
    module procedure smooth_min_scalars

  end interface

!===============================================================================
  interface smooth_max

    ! These functions smooth the output of the max functions 
    ! by introducing a varyingly steep path between the two input variables.
    ! The degree to which smoothing is applied depends on the value of 'smth_coef'.
    ! If 'smth_coef' goes toward 0, the output of the max function will be 
    !        0.5 * ((a+b) + abs(a-b))
    ! If a > b, then this comes out to be a. Likewise, if a < b, abs(a-b)=b-a so we get b.
    ! Increasing the smoothing coefficient will lead to a greater degree of smoothing
    ! in the smooth min and max functions. Generally, the coefficient should roughly scale
    ! with the magnitude of data in the data structure that is to be smoothed, in order to
    ! obtain a sensible degree of smoothing (not too much, not too little).

    module procedure smooth_max_scalar_array
    module procedure smooth_max_array_scalar
    module procedure smooth_max_arrays
    module procedure smooth_max_scalars

  end interface

!===============================================================================
  contains

  !---------------------------------------------------------------------------
  subroutine set_boundary_conditions_lhs( diag_index, low_bound, high_bound, &
                                          lhs, &
                                          diag_index2, low_bound2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a left-hand side LAPACK matrix.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)
        
    use constants_clubb, only: &
        one, zero

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables
    integer, intent(in) :: &
      diag_index, low_bound, high_bound ! boundary indexes for the first variable

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(:,:), intent(inout) :: &
      lhs ! left hand side of the LAPACK matrix equation

    ! Optional Input Variables
    integer, intent(in), optional :: &
      diag_index2, low_bound2, high_bound2 ! boundary indexes for the second variable

    ! --------------------- BEGIN CODE ----------------------

    if ( ( present( low_bound2 ) .or. present( high_bound2 ) ) .and. &
         ( .not. present( diag_index2 ) ) ) then

      error stop "Boundary index provided without diag_index."

    end if

    ! Set the lower boundaries for the first variable
    lhs(:,low_bound) = zero
    lhs(diag_index,low_bound) = one

    ! Set the upper boundaries for the first variable
    lhs(:,high_bound) = zero
    lhs(diag_index,high_bound) = one

    ! Set the lower boundaries for the second variable, if it is provided
    if ( present( low_bound2 ) ) then

      lhs(:,low_bound2) = zero
      lhs(diag_index2,low_bound2) = one

    end if

    ! Set the upper boundaries for the second variable, if it is provided
    if ( present( high_bound2 ) ) then

      lhs(:,high_bound2) = zero
      lhs(diag_index2,high_bound2) = one

    end if

    return
  end subroutine set_boundary_conditions_lhs

  !--------------------------------------------------------------------------
  subroutine set_boundary_conditions_rhs( &
               low_value, low_bound, high_value, high_bound, &
               rhs, &
               low_value2, low_bound2, high_value2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a right-hand side LAPACK vector.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables

    ! The values for the first variable
    real( kind = core_rknd ), intent(in) :: low_value, high_value

    ! The bounds for the first variable
    integer, intent(in) :: low_bound, high_bound

    ! Input / Output Variables

    ! The right-hand side vector
    real( kind = core_rknd ), dimension(:), intent(inout) :: rhs

    ! Optional Input Variables

    ! The values for the second variable
    real( kind = core_rknd ), intent(in), optional :: low_value2, high_value2

    ! The bounds for the second variable
    integer, intent(in), optional :: low_bound2, high_bound2


    ! -------------------- BEGIN CODE ------------------------

    ! Stop execution if a boundary was provided without a value
    if ( (present( low_bound2 ) .and. (.not. present( low_value2 ))) .or. &
         (present( high_bound2 ) .and. (.not. present( high_value2 ))) ) then

      error stop "Boundary condition provided without value."

    end if

    ! Set the lower and upper bounds for the first variable
    rhs(low_bound) = low_value
    rhs(high_bound) = high_value

    ! If a lower bound was given for the second variable, set it
    if ( present( low_bound2 ) ) then
      rhs(low_bound2) = low_value2
    end if

    ! If an upper bound was given for the second variable, set it
    if ( present( high_bound2 ) ) then
      rhs(high_bound2) = high_value2
    end if

    return
  end subroutine set_boundary_conditions_rhs

  !===============================================================================
  subroutine calc_stability_correction( nz, ngrdcol, gr, &
                                        thlm, Lscale, em, &
                                        exner, rtm, rcm, &
                                        p_in_Pa, thvm, ice_supersat_frac, &
                                        lambda0_stability_coef, &
                                        bv_efold, &
                                        l_brunt_vaisala_freq_moist, &
                                        l_use_thvm_in_bv_freq, &
                                        stability_correction )
  !
  ! Description:
  !   Stability Factor
  !
  ! References:
  !
  !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero, one, three    ! Constant(s)

    use grid_class, only:  &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    implicit none

    ! ---------------- Input Variables ----------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      Lscale,          & ! Turbulent mixing length                   [m]
      em,              & ! Turbulent Kinetic Energy (TKE)            [m^2/s^2]
      thlm,            & ! th_l (thermo. levels)                     [K]
      exner,           & ! Exner function                            [-]
      rtm,             & ! total water mixing ratio, r_t             [kg/kg]
      rcm,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,         & ! Air pressure                              [Pa]
      thvm,            & ! Virtual potential temperature             [K]
      ice_supersat_frac

    real( kind = core_rknd ), intent(in) :: &
      lambda0_stability_coef, &     ! CLUBB tunable parameter lambda0_stability_coef
      bv_efold                      ! Control parameter for inverse e-folding of
                                    ! cloud fraction in the mixed Brunt Vaisala frequency

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq         ! Use thvm in the calculation of Brunt-Vaisala frequency

    ! ---------------- Output Variables ----------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz) :: &
      stability_correction
      
    ! ---------------- Local Variables ----------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_vaisala_freq_sqd, & !  []
      brunt_vaisala_freq_sqd_mixed, &
      brunt_vaisala_freq_sqd_dry, & !  []
      brunt_vaisala_freq_sqd_moist, &
      lambda0_stability, &
      Lscale_zm

    integer :: i, k

    !------------ Begin Code --------------

    !$acc enter data create( brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                    brunt_vaisala_freq_sqd_moist, brunt_vaisala_freq_sqd_dry, &
    !$acc                    lambda0_stability, Lscale_zm )

    call calc_brunt_vaisala_freq_sqd( nz, ngrdcol, gr, thlm, &          ! intent(in)
                                      exner, rtm, rcm, p_in_Pa, thvm, & ! intent(in)
                                      ice_supersat_frac, &              ! intent(in)
                                      l_brunt_vaisala_freq_moist, &     ! intent(in)
                                      l_use_thvm_in_bv_freq, &          ! intent(in)
                                      bv_efold, &                       ! intent(in)
                                      brunt_vaisala_freq_sqd, &         ! intent(out)
                                      brunt_vaisala_freq_sqd_mixed,&    ! intent(out)
                                      brunt_vaisala_freq_sqd_dry, &     ! intent(out)
                                      brunt_vaisala_freq_sqd_moist )    ! intent(out)

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        if ( brunt_vaisala_freq_sqd(i,k) > zero  ) then
          lambda0_stability(i,k) = lambda0_stability_coef
        else
          lambda0_stability(i,k) = zero
        end if
      end do
    end do
    !$acc end parallel loop

    Lscale_zm = zt2zm( nz, ngrdcol, gr, Lscale(:,:) )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        stability_correction(i,k) = one + min( lambda0_stability(i,k) * brunt_vaisala_freq_sqd(i,k) &
                                               * Lscale_zm(i,k)**2 / em(i,k), three )
      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                   brunt_vaisala_freq_sqd_moist, brunt_vaisala_freq_sqd_dry, &
    !$acc                   lambda0_stability, Lscale_zm )

    return

  end subroutine calc_stability_correction

  !===============================================================================
  subroutine calc_brunt_vaisala_freq_sqd(  nz, ngrdcol, gr, thlm, &
                                           exner, rtm, rcm, p_in_Pa, thvm, &
                                           ice_supersat_frac, &
                                           l_brunt_vaisala_freq_moist, &
                                           l_use_thvm_in_bv_freq, &
                                           bv_efold, &
                                           brunt_vaisala_freq_sqd, &
                                           brunt_vaisala_freq_sqd_mixed,&
                                           brunt_vaisala_freq_sqd_dry, &
                                           brunt_vaisala_freq_sqd_moist )

  ! Description:
  !   Calculate the Brunt-Vaisala frequency squared, N^2.

  ! References:
  !   ?
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Konstant

    use constants_clubb, only: &
        grav, & ! Constant(s)
        Lv, &
        Cp, &
        Rd, &
        ep, &
        one, &
        one_half, &
        zero_threshold

    use parameters_model, only: & 
        T0 ! Variable! 

    use grid_class, only: &
        grid, & ! Type
        ddzt,   &  ! Procedure(s)
        zt2zm

    use T_in_K_module, only: &
        thlm2T_in_K ! Procedure

    use saturation, only: &
        sat_mixrat_liq ! Procedure

    implicit none

    !---------------------------- Input Variables ----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,    &  ! th_l (thermo. levels)              [K]
      exner,   &  ! Exner function                     [-]
      rtm,     &  ! total water mixing ratio, r_t      [kg/kg]
      rcm,     &  ! cloud water mixing ratio, r_c      [kg/kg]
      p_in_Pa, &  ! Air pressure                       [Pa]
      thvm,    &  ! Virtual potential temperature      [K]
      ice_supersat_frac

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq         ! Use thvm in the calculation of Brunt-Vaisala frequency

    real( kind = core_rknd ), intent(in) :: &
      bv_efold                      ! Control parameter for inverse e-folding of
                                    ! cloud fraction in the mixed Brunt Vaisala frequency

    !---------------------------- Output Variables ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      brunt_vaisala_freq_sqd, & ! Brunt-Vaisala frequency squared, N^2 [1/s^2]
      brunt_vaisala_freq_sqd_mixed, &
      brunt_vaisala_freq_sqd_dry,&
      brunt_vaisala_freq_sqd_moist

    !---------------------------- Local Variables ----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      T_in_K, T_in_K_zm, rsat, rsat_zm, thm, thm_zm, ddzt_thlm, &
      ddzt_thm, ddzt_rsat, ddzt_rtm, thvm_zm, ddzt_thvm

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      stat_dry, stat_liq, ddzt_stat_liq, ddzt_stat_liq_zm, &
      stat_dry_virtual, stat_dry_virtual_zm,  ddzt_rtm_zm

    integer :: i, k

    !---------------------------- Begin Code ----------------------------

    !$acc data copyin( gr, gr%zt, &
    !$acc              thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac ) &
    !$acc      copyout( brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc               brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist ) &
    !$acc       create( T_in_K, T_in_K_zm, rsat, rsat_zm, thm, thm_zm, ddzt_thlm, &
    !$acc               ddzt_thm, ddzt_rsat, ddzt_rtm, thvm_zm, ddzt_thvm, stat_dry, &
    !$acc               stat_liq, ddzt_stat_liq, ddzt_stat_liq_zm, stat_dry_virtual, &
    !$acc               stat_dry_virtual_zm, ddzt_rtm_zm )

    ddzt_thlm = ddzt( nz, ngrdcol, gr, thlm )
    thvm_zm = zt2zm( nz, ngrdcol, gr, thvm )
    ddzt_thvm = ddzt( nz, ngrdcol, gr, thvm )

    ! Dry Brunt-Vaisala frequency
    if ( l_use_thvm_in_bv_freq ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_vaisala_freq_sqd(i,k) = ( grav / thvm_zm(i,k) ) * ddzt_thvm(i,k)
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_vaisala_freq_sqd(i,k) = ( grav / T0 ) * ddzt_thlm(i,k)
        end do
      end do
      !$acc end parallel loop

    end if

    T_in_K = thlm2T_in_K( nz, ngrdcol, thlm, exner, rcm )
    T_in_K_zm = zt2zm( nz, ngrdcol, gr, T_in_K )
    rsat = sat_mixrat_liq( nz, ngrdcol, p_in_Pa, T_in_K )
    rsat_zm = zt2zm( nz, ngrdcol, gr, rsat )
    ddzt_rsat = ddzt( nz, ngrdcol, gr, rsat )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        thm(i,k) = thlm(i,k) + Lv/(Cp*exner(i,k)) * rcm(i,k)
      end do
    end do
    !$acc end parallel loop

    thm_zm = zt2zm( nz, ngrdcol, gr, thm )
    ddzt_thm = ddzt( nz, ngrdcol, gr, thm )
    ddzt_rtm = ddzt( nz, ngrdcol, gr, rtm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        stat_dry(i,k)  =  Cp * T_in_K(i,k) + grav * gr%zt(i,k)
        stat_liq(i,k)  =  stat_dry(i,k) - Lv * rcm(i,k)
      end do
    end do
    !$acc end parallel loop

    ddzt_stat_liq    = ddzt( nz, ngrdcol, gr, stat_liq )
    ddzt_stat_liq_zm = zt2zm( nz, ngrdcol, gr, ddzt_stat_liq)

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        stat_dry_virtual(i,k) = stat_dry(i,k) + Cp * T_in_K(i,k) &
                                                *( 0.608 * ( rtm(i,k) - rcm(i,k) )- rcm(i,k) )
      end do
    end do
    !$acc end parallel loop

    stat_dry_virtual_zm = zt2zm( nz, ngrdcol, gr, stat_dry_virtual)
    ddzt_rtm_zm         = zt2zm( nz, ngrdcol, gr, ddzt_rtm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        brunt_vaisala_freq_sqd_dry(i,k) = ( grav / thm_zm(i,k) )* ddzt_thm(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        ! In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and
        ! Klemp (1982)
        brunt_vaisala_freq_sqd_moist(i,k) = &
          grav * ( ((one + Lv*rsat_zm(i,k) / (Rd*T_in_K_zm(i,k))) / &
          (one + ep*(Lv**2)*rsat_zm(i,k)/(Cp*Rd*T_in_K_zm(i,k)**2))) * &
          ( (one/thm_zm(i,k) * ddzt_thm(i,k)) + (Lv/(Cp*T_in_K_zm(i,k)))*ddzt_rsat(i,k)) - &
          ddzt_rtm(i,k) )
      end do
    end do ! k=1, gr%nz
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
         brunt_vaisala_freq_sqd_mixed(i,k) = &
             brunt_vaisala_freq_sqd_moist(i,k) + &
                 exp( - bv_efold * ice_supersat_frac(i,k) ) * &
                 ( brunt_vaisala_freq_sqd_dry(i,k) - brunt_vaisala_freq_sqd_moist(i,k) )
      end do
    end do
    !$acc end parallel loop

    if ( l_brunt_vaisala_freq_moist ) then

      brunt_vaisala_freq_sqd = brunt_vaisala_freq_sqd_moist

    end if

    !$acc end data

    return

  end subroutine calc_brunt_vaisala_freq_sqd

!===============================================================================
  subroutine compute_Cx_fnc_Richardson( nz, ngrdcol, gr, &
                                        thlm, um, vm, em, Lscale, exner, rtm, &
                                        rcm, p_in_Pa, thvm, rho_ds_zm, &
                                        ice_supersat_frac, &
                                        clubb_params, &
                                        l_brunt_vaisala_freq_moist, &
                                        l_use_thvm_in_bv_freq, &
                                        l_use_shear_Richardson, &
                                        l_modify_limiters_for_cnvg_test, & 
                                        stats_metadata, &
                                        stats_zm, & 
                                        Cx_fnc_Richardson )

  ! Description:
  !   Compute Cx as a function of the Richardson number

  ! References:
  !   cam:ticket:59
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Konstant

    use grid_class, only: &
        grid, & ! Type
        ddzt, & ! Procedure(s)
        zt2zm, & 
        zm2zt2zm

    use constants_clubb, only: &
        one, zero

    use interpolation, only: &
        linear_interp_factor ! Procedure

    use parameter_indices, only: &
        nparams,             & ! Variable(s)
        iCx_min,             &
        iCx_max,             &
        iRichardson_num_min, &
        iRichardson_num_max, &
        ibv_efold

    use stats_variables, only: &
        stats_metadata_type

    use stats_type_utilities, only: &
        stat_update_var      ! Procedure

    use stats_type, only: stats ! Type

    implicit none

    !------------------------------ Constant Parameters ------------------------------
    real( kind = core_rknd ), parameter :: &
      Richardson_num_divisor_threshold = 1.0e-6_core_rknd, &
      Cx_fnc_Richardson_below_ground_value = one

    logical, parameter :: &
      l_Cx_fnc_Richardson_vert_avg = .false.    ! Vertically average Cx_fnc_Richardson over a
                                                !  distance of Lscale

    real( kind = core_rknd ), parameter :: &
      min_max_smth_mag = 1.0e-9_core_rknd ! "base" smoothing magnitude before scaling 
                                          ! for the respective data structure. See
                                          ! https://github.com/larson-group/clubb/issues/965#issuecomment-1119816722
                                          ! for a plot on how output behaves with varying min_max_smth_mag

    !------------------------------ Input Variables ------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,      & ! th_l (liquid water potential temperature)      [K]
      um,        & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,        & ! v mean wind component (thermodynamic levels)   [m/s]
      em,        & ! Turbulent Kinetic Energy (TKE)                 [m^2/s^2]
      Lscale,    & ! Turbulent mixing length                        [m]
      exner,     & ! Exner function                                 [-]
      rtm,       & ! total water mixing ratio, r_t                  [kg/kg]
      rcm,       & ! cloud water mixing ratio, r_c                  [kg/kg]
      p_in_Pa,   & ! Air pressure                                   [Pa]
      thvm,      & ! Virtual potential temperature                  [K]
      rho_ds_zm, &  ! Dry static density on momentum levels          [kg/m^3]
      ice_supersat_frac  ! ice cloud fraction

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,      & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_use_shear_Richardson        ! Use shear in the calculation of Richardson number

    ! Flag to activate modifications on limiters for convergence test 
    ! (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
    ! (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
    ! (reduce threshold on limiters for sqrt_Ri_zm in mixing_length.F90)
    logical, intent(in) :: &
      l_modify_limiters_for_cnvg_test

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------------ InOut Variable ------------------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zm

    !------------------------------ Output Variable ------------------------------
    real( kind = core_rknd), dimension(ngrdcol,nz), intent(out) :: &
      Cx_fnc_Richardson

    !------------------------------ Local Variables ------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_vaisala_freq_sqd, &
      brunt_vaisala_freq_sqd_mixed,&
      brunt_vaisala_freq_sqd_dry, &
      brunt_vaisala_freq_sqd_moist, &
      fnc_Richardson, &
      fnc_Richardson_clipped, &
      fnc_Richardson_smooth, &
      Ri_zm, &
      ddzt_um, &
      ddzt_vm, &
      shear_sqd, &
      Lscale_zm, &
      Cx_fnc_interp, &
      Cx_fnc_Richardson_avg

    real ( kind = core_rknd ) :: &
      invrs_min_max_diff, &
      invrs_num_div_thresh

    real( kind = core_rknd ) :: &
      Richardson_num_max, & ! CLUBB tunable parameter Richardson_num_max
      Richardson_num_min, & ! CLUBB tunable parameter Richardson_num_min
      Cx_max,             & ! CLUBB tunable parameter max of Cx_fnc_Richardson
      Cx_min                ! CLUBB tunable parameter min of Cx_fnc_Richardson

    integer :: smth_type = 1

    integer :: i, k

    !------------------------------ Begin Code ------------------------------

    !$acc enter data create( brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                    brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    !$acc                    Cx_fnc_interp, &
    !$acc                    Ri_zm, ddzt_um, ddzt_vm, shear_sqd, Lscale_zm, &
    !$acc                    Cx_fnc_Richardson_avg, fnc_Richardson, &
    !$acc                    fnc_Richardson_clipped, fnc_Richardson_smooth )

    call calc_brunt_vaisala_freq_sqd( nz, ngrdcol, gr, thlm, &          ! intent(in)
                                      exner, rtm, rcm, p_in_Pa, thvm, & ! intent(in)
                                      ice_supersat_frac, &              ! intent(in)
                                      l_brunt_vaisala_freq_moist, &     ! intent(in)
                                      l_use_thvm_in_bv_freq, &          ! intent(in)
                                      clubb_params(ibv_efold), &        ! intent(in)
                                      brunt_vaisala_freq_sqd, &         ! intent(out)
                                      brunt_vaisala_freq_sqd_mixed,&    ! intent(out)
                                      brunt_vaisala_freq_sqd_dry, &     ! intent(out)
                                      brunt_vaisala_freq_sqd_moist )    ! intent(out)

    Richardson_num_max = clubb_params(iRichardson_num_max)
    Richardson_num_min = clubb_params(iRichardson_num_min)
    Cx_max = clubb_params(iCx_max)
    Cx_min = clubb_params(iCx_min)

    invrs_min_max_diff = one / ( Richardson_num_max - Richardson_num_min )
    invrs_num_div_thresh = one / Richardson_num_divisor_threshold

    Lscale_zm = zt2zm( nz, ngrdcol, gr, Lscale )

    ! Calculate shear_sqd
    ddzt_um = ddzt( nz, ngrdcol, gr, um )
    ddzt_vm = ddzt( nz, ngrdcol, gr, vm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        shear_sqd(i,k) = ddzt_um(i,k)**2 + ddzt_vm(i,k)**2
      end do
    end do
    !$acc end parallel loop

    if ( stats_metadata%l_stats_samp ) then
      !$acc update host(shear_sqd)
      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%ishear_sqd, shear_sqd(i,:), & ! intent(in)
                              stats_zm(i) )               ! intent(inout)
      end do
    end if

    if ( l_use_shear_Richardson ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Ri_zm(i,k) = max( 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_mixed(i,k) ) &
                       / max( shear_sqd(i,k), 1.0e-7_core_rknd )
        end do
      end do
      !$acc end parallel loop

    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Ri_zm(i,k) = brunt_vaisala_freq_sqd(i,k) * invrs_num_div_thresh
        end do
      end do
      !$acc end parallel loop
    end if

    ! Cx_fnc_Richardson is interpolated based on the value of Richardson_num
    ! The min function ensures that Cx does not exceed Cx_max, regardless of the
    !     value of Richardson_num_max.
    if ( l_modify_limiters_for_cnvg_test ) then 

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          fnc_Richardson(i,k) = ( Ri_zm(i,k) - Richardson_num_min ) * invrs_min_max_diff
        end do
      end do

      fnc_Richardson_clipped = smooth_min( nz, ngrdcol, one, &
                                           fnc_Richardson, &
                                           min_max_smth_mag )

      fnc_Richardson_smooth = smooth_max( nz, ngrdcol, zero, &
                                          fnc_Richardson_clipped, &
                                          min_max_smth_mag )

      ! use smoothed max amd min to achive smoothed profile and avoid discontinuities 
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Cx_fnc_interp(i,k) = fnc_Richardson_smooth(i,k) * ( Cx_max - Cx_min ) + Cx_min
        end do
      end do

      Cx_fnc_Richardson = zm2zt2zm( nz, ngrdcol, gr, Cx_fnc_interp )

    else ! default method 

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Cx_fnc_Richardson(i,k) = ( max(min(Richardson_num_max, Ri_zm(i,k)), Richardson_num_min) &
                                     - Richardson_num_min )  &
                                   * invrs_min_max_diff * ( Cx_max - Cx_min ) + Cx_min
        end do
      end do
      !$acc end parallel loop

    end if 

    if ( l_Cx_fnc_Richardson_vert_avg ) then
      Cx_fnc_Richardson = Lscale_width_vert_avg( nz, ngrdcol, gr, smth_type, &
                                                 Cx_fnc_Richardson, Lscale_zm, rho_ds_zm, &
                                                 Cx_fnc_Richardson_below_ground_value )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Cx_fnc_Richardson(i,k) = Cx_fnc_Richardson_avg(i,k)
        end do
      end do
      !$acc end parallel loop
    end if

    ! On some compilers, roundoff error can result in Cx_fnc_Richardson being
    ! slightly outside the range [0,1]. Thus, it is clipped here.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        Cx_fnc_Richardson(i,k) = max( zero, min( one, Cx_fnc_Richardson(i,k) ) )
      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                   brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    !$acc                   Cx_fnc_interp, Ri_zm, &
    !$acc                   ddzt_um, ddzt_vm, shear_sqd, Lscale_zm, &
    !$acc                   Cx_fnc_Richardson_avg, fnc_Richardson, &
    !$acc                   fnc_Richardson_clipped, fnc_Richardson_smooth )

    return

  end subroutine compute_Cx_fnc_Richardson
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  function Lscale_width_vert_avg( nz, ngrdcol, gr, smth_type, &
                                  var_profile, Lscale_zm, rho_ds_zm, &
                                  var_below_ground_value )&
  result (Lscale_width_vert_avg_output)

  ! Description:
  !   Averages a profile with a running mean of width Lscale_zm

  ! References:
  !   cam:ticket:59

    use clubb_precision, only: &
        core_rknd ! Precision

    use grid_class, only: &
        grid ! Type
        
    use constants_clubb, only: &
        zero

    implicit none

    !-------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol, &
      smth_type
      
    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      var_profile, &      ! Profile on momentum levels
      Lscale_zm, &        ! Lscale on momentum levels
      rho_ds_zm           ! Dry static energy on momentum levels!

    real( kind = core_rknd ), intent(in) :: &
      var_below_ground_value ! Value to use below ground

    !-------------------------- Result Variable --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Lscale_width_vert_avg_output ! Vertically averaged profile (on momentum levels)

    !-------------------------- Local Variables --------------------------
    integer :: &
        k, i,        & ! Loop variable
        k_avg_lower, &
        k_avg_upper, &
        k_avg

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      one_half_avg_width, &
      numer_terms, &
      denom_terms

    integer :: &
      n_below_ground_levels

    real( kind = core_rknd ) :: & 
      numer_integral, & ! Integral in the numerator (see description)
      denom_integral    ! Integral in the denominator (see description)

    !-------------------------- Begin Code --------------------------

    !$acc enter data create( one_half_avg_width, numer_terms, denom_terms )

    if ( smth_type == 1 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          one_half_avg_width(i,k) = max( Lscale_zm(i,k), 500.0_core_rknd )
        end do
      end do
    else if (smth_type == 2 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          one_half_avg_width(i,k) = 60.0_core_rknd
        end do
      end do
    endif

    ! Pre calculate numerator and denominator terms
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        numer_terms(i,k) = rho_ds_zm(i,k) * gr%dzm(i,k) * var_profile(i,k)
        denom_terms(i,k) = rho_ds_zm(i,k) * gr%dzm(i,k)
      end do
    end do

    ! For every grid level
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        !-----------------------------------------------------------------------
        ! Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
        ! 
        ! Note: Outdated explanation of version that improves CPU performance
        !       below. Reworked due to it requiring a k dependency. Now we
        !       begin looking for k_avg_upper and k_avg_lower starting at 
        !       the kth level.
        ! 
        ! Outdated but potentially useful note:
        !     k_avg_upper and k_avg_lower can be saved each loop iteration, this 
        !     reduces iterations beacuse the kth values are likely to be within
        !     one or two grid levels of the k-1th values. Less searching is required
        !     by starting the search at the previous values and incrementing or 
        !     decrement as needed.
        !-----------------------------------------------------------------------

        ! Determine if k_avg_upper needs to increment
        do k_avg_upper = k, nz-1
          if ( gr%zm(i,k_avg_upper+1) - gr%zm(i,k) > one_half_avg_width(i,k) ) then
            exit
          end if
        end do

        ! Determine if k_avg_lower needs to decrement
        do k_avg_lower = k, 2, -1
          if ( gr%zm(i,k) - gr%zm(i,k_avg_lower-1) > one_half_avg_width(i,k) ) then
            exit
          end if
        end do

        ! Compute the number of levels below ground to include.
        if ( k_avg_lower > 1 ) then

          ! k=1, the lowest "real" level, is not included in the average, so no
          ! below-ground levels should be included.
          n_below_ground_levels = 0

          numer_integral = zero
          denom_integral = zero

        else

          ! The number of below-ground levels included is equal to the distance
          ! below the lowest level spanned by one_half_avg_width(k)
          ! divided by the distance between vertical levels below ground; the
          ! latter is assumed to be the same as the distance between the first and
          ! second vertical levels.
          n_below_ground_levels = int( ( one_half_avg_width(i,k)-(gr%zm(i,k)-gr%zm(i,1)) ) / &
                                      ( gr%zm(i,2)-gr%zm(i,1) ) )

          numer_integral = n_below_ground_levels * denom_terms(i,1) * var_below_ground_value
          denom_integral = n_below_ground_levels * denom_terms(i,1)

        end if
            
        ! Add numerator and denominator terms for all above-ground levels
        do k_avg = k_avg_lower, k_avg_upper

          numer_integral = numer_integral + numer_terms(i,k_avg)
          denom_integral = denom_integral + denom_terms(i,k_avg)

        end do

        Lscale_width_vert_avg_output(i,k) = numer_integral / denom_integral

      end do
    end do

    !$acc exit data delete( one_half_avg_width, numer_terms, denom_terms )

    return

  end function Lscale_width_vert_avg

 !=============================================================================
  subroutine wp2_term_splat_lhs( nz, ngrdcol, gr, C_wp2_splat, &
                                 brunt_vaisala_freq_sqd_splat, &
                                 lhs_splat_wp2 )

    ! Description
    ! DESCRIBE TERM

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        grid, & ! Type
        zm2zt2zm

    use constants_clubb, only: &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! --------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      brunt_vaisala_freq_sqd_splat  ! Inverse time-scale tau at momentum levels  [1/s^2]

    real( kind = core_rknd ), intent(in) :: &
      C_wp2_splat    ! Model parameter C_wp2_splat             [ -]

    ! --------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      lhs_splat_wp2    ! LHS coefficient of wp2 splatting term  [1/s]

    ! --------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_vaisala_freq_splat_clipped, &
      brunt_vaisala_freq_splat_smooth

    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc enter data create( brunt_vaisala_freq_splat_clipped, brunt_vaisala_freq_splat_smooth )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        brunt_vaisala_freq_splat_clipped(i,k) &
                = sqrt( max( zero, brunt_vaisala_freq_sqd_splat(i,k) ) )
      end do
    end do
    !$acc end parallel loop
    
    brunt_vaisala_freq_splat_smooth = zm2zt2zm( nz, ngrdcol, gr, &
                                                brunt_vaisala_freq_splat_clipped )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        lhs_splat_wp2(i,k) = + C_wp2_splat * brunt_vaisala_freq_splat_smooth(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_vaisala_freq_splat_clipped, brunt_vaisala_freq_splat_smooth )

    return

  end subroutine wp2_term_splat_lhs

 !=============================================================================
  subroutine wp3_term_splat_lhs( nz, ngrdcol, gr, C_wp2_splat, &
                                 brunt_vaisala_freq_sqd_splat, &
                                 lhs_splat_wp3 )

    ! Description
    ! DESCRIBE TERM

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        grid, & ! Type
        zm2zt2zm

    use constants_clubb, only: &
        zero, &
        one_half, &
        three

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! --------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      brunt_vaisala_freq_sqd_splat  ! Inverse time-scale tau at momentum levels  [1/s^2]

    real( kind = core_rknd ), intent(in) :: &
      C_wp2_splat    ! Model parameter C_wp2_splat              [-]

    ! --------------------- Output Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      lhs_splat_wp3    ! LHS coefficient of wp3 splatting term [1/s]

    ! --------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_vaisala_freq_splat_clipped, &
      brunt_vaisala_freq_splat_smooth

    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc enter data create( brunt_vaisala_freq_splat_clipped, brunt_vaisala_freq_splat_smooth )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        brunt_vaisala_freq_splat_clipped(i,k) &
                = sqrt( max( zero, brunt_vaisala_freq_sqd_splat(i,k) ) )
      end do
    end do
    !$acc end parallel loop
    
    brunt_vaisala_freq_splat_smooth = zm2zt2zm( nz, ngrdcol, gr, &
                                                brunt_vaisala_freq_splat_clipped )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        lhs_splat_wp3(i,k) = + one_half * three * C_wp2_splat &
                               * brunt_vaisala_freq_splat_smooth(i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_vaisala_freq_splat_clipped, brunt_vaisala_freq_splat_smooth )

    return

  end subroutine wp3_term_splat_lhs

!===============================================================================
  function smooth_min_scalar_array( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function, using one scalar and
  !   one 1d array as inputs. For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    !----------------------------- Input Variables -----------------------------
    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var2          ! Units vary

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var2 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( (input_var1+input_var2(i,k)) - &
                                  sqrt((input_var1-input_var2(i,k))**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_min_scalar_array

!===============================================================================
  function smooth_min_array_scalar( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function, using one scalar and 
  !   one 1d array as inputs. For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)

    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1          ! Units vary

    real ( kind = core_rknd ), intent(in) :: &
      input_var2, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var1 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( (input_var1(i,k)+input_var2) - &
                                  sqrt((input_var1(i,k)-input_var2)**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_min_array_scalar

!===============================================================================
  function smooth_min_arrays( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function, using two 1d arrays as inputs.
  !   For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables-----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2          ! Units vary
      
    real ( kind = core_rknd ), intent(in) :: &
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var1, input_var2 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( (input_var1(i,k)+input_var2(i,k)) - &
                                  sqrt((input_var1(i,k)-input_var2(i,k))**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_min_arrays

!===============================================================================
  function smooth_min_scalars( input_var1, input_var2, smth_coef ) &
  result( output_var )
  !$acc routine

  ! Description:
  !   Computes a smoothed version of the min function, using two scalars as inputs.
  !   For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

  ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

  ! Output Variables
    real( kind = core_rknd ) :: &
      output_var          ! Same unit as input_var1 and input_var2

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) - &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_min_scalars

!===============================================================================
  function smooth_max_scalar_array( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function, using one scalar and 
  !   one 1d array as inputs. For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var2          ! Units vary

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var2 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( (input_var1+input_var2(i,k)) + &
                                  sqrt((input_var1-input_var2(i,k))**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_max_scalar_array

!===============================================================================
  function smooth_max_array_scalar( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function, using one scalar and 
  !   one 1d array as inputs. For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1          ! Units vary

    real ( kind = core_rknd ), intent(in) :: &
      input_var2, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var1 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( ( input_var1(i,k) + input_var2 ) + &
                                  sqrt(( input_var1(i,k) - input_var2 )**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_max_array_scalar

!===============================================================================
  function smooth_max_arrays( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function, using two 1d arrays as inputs.
  !   For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2          ! Units vary
      
    real( kind = core_rknd ), intent(in) :: &
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    !$acc data copyin( input_var1, input_var2 ) &
    !$acc     copyout( output_var )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        output_var(i,k) = one_half * ( (input_var1(i,k)+input_var2(i,k)) + &
                                  sqrt((input_var1(i,k)-input_var2(i,k))**2 + smth_coef**2) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end function smooth_max_arrays

!===============================================================================
  function smooth_max_scalars( input_var1, input_var2, smth_coef ) &
  result( output_var )
  !$acc routine

  ! Description:
  !   Computes a smoothed version of the max function, using two scalars as inputs.
  !   For more details, see the interface in this file.

  ! References:
  !   See clubb:ticket: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none

    !----------------------------- Input Variables -----------------------------
    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2, &       ! Units vary
      smth_coef           ! "intensity" of the smoothing. Should be of a similar magnitude to
                          ! that of the data structures input_var1 and input_var2

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ) :: &
      output_var          ! Same unit as input_var1 and input_var2

    !----------------------------- Local Variables -----------------------------
    integer :: i, k

    !----------------------------- Begin Code -----------------------------

    output_var = one_half * ( (input_var1+input_var2) + &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )
    return

  end function smooth_max_scalars

  function smooth_heaviside_peskin( nz, ngrdcol, input, smth_range ) &
    result( smth_output )
    
  ! Description:
  !   Computes a smoothed heaviside function as in 
  !       [Lin, Lee et al., 2005, A level set characteristic Galerkin 
  !       finite element method for free surface flows], equation (2)
  
  ! References:
  !   See clubb:ticket:965
  !----------------------------------------------------------------------
  
    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        pi, invrs_pi, one, one_half, zero

    implicit none

    !------------------------- Input Variables -------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real ( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      input    ! Units vary

    real ( kind = core_rknd ), intent(in) :: &
      smth_range  ! Smooth Heaviside function on [-smth_range, smth_range]

    !------------------------- Output Variables -------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      smth_output    ! Same units as input
    
    !------------------------- Local Variables -------------------------
    real ( kind = core_rknd ) :: &
      input_over_smth_range  ! input divided by smth_range

    integer :: i, k
      
    !------------------------- Begin Code -------------------------

    !$acc data copyin( input ) &
    !$acc     copyout( smth_output )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        if ( input(i,k) < -smth_range ) then 
          smth_output(i,k) = zero
        elseif ( input(i,k) > smth_range ) then
           smth_output(i,k) = one
        else 
          ! Note that this case will only ever be reached if smth_range != 0,
          ! so this division is fine and should not cause any issues
          input_over_smth_range = input(i,k) / smth_range
          smth_output(i,k) = one_half &
                             * (one + input_over_smth_range &
                               + invrs_pi * sin(pi * input_over_smth_range))
        end if

      end do
    end do
    !$acc end parallel loop

    !$acc end data
    
    return

  end function smooth_heaviside_peskin

  !===============================================================================
  subroutine calc_xpwp_1D( gr, Km_zm, xm, &
                           xpwp )

    ! Description:
    ! Compute x'w' from x<k>, x<k+1>, Kh and invrs_dzm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)
        
    use grid_class, only: &
      grid

    implicit none

    ! ----------------------- Input variables -----------------------
    type (grid), target, intent(in) :: gr
      
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Km_zm,     & ! Eddy diff. (k momentum level)                 [m^2/s]
      xm           ! x (k thermo level)                            [units vary]
      
    ! ----------------------- Output variable -----------------------
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      xpwp ! x'w'   [(units vary)(m/s)]
      
    integer :: k

    ! ----------------------- Begin Code -----------------------

    ! Solve for x'w' at all intermediate model levels.
    do k = 1, gr%nz-1
      xpwp(k) = Km_zm(k) * gr%invrs_dzm(1,k) * ( xm(k+1) - xm(k) )
    end do

    return
  end subroutine calc_xpwp_1D
  
  !===============================================================================
  subroutine calc_xpwp_2D( nz, ngrdcol, gr, &
                        Km_zm, xm, &
                        xpwp )

    ! Description:
    ! Compute x'w' from x<k>, x<k+1>, Kh and invrs_dzm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)
        
    use grid_class, only: &
      grid

    implicit none

    ! ----------------------- Input variables -----------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
      
    type (grid), target, intent(in) :: gr
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Km_zm,     & ! Eddy diff. (k momentum level)                 [m^2/s]
      xm           ! x (k thermo level)                            [units vary]
      
    ! ----------------------- Output variable -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      xpwp ! x'w'   [(units vary)(m/s)]
      
    integer :: i, k

    ! ----------------------- Begin Code -----------------------

    !$acc data copyin( gr, gr%invrs_dzm, Km_zm, xm ) &
    !$acc     copyout( xpwp )

    ! Solve for x'w' at all intermediate model levels.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz-1
      do i = 1, ngrdcol
        xpwp(i,k) = Km_zm(i,k) * gr%invrs_dzm(i,k) * ( xm(i,k+1) - xm(i,k) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine calc_xpwp_2D

  !=============================================================================
  function vertical_avg( total_idx, rho_ds, field, dz )

    ! Description:
    ! Computes the density-weighted vertical average of a field.
    !
    ! The average value of a function, f, over a set domain, [a,b], is
    ! calculated by the equation:
    !
    ! f_avg = ( INT(a:b) f*g ) / ( INT(a:b) g );
    !
    ! as long as f is continous and g is nonnegative and integrable.  Therefore,
    ! the density-weighted (by dry, static, base-static density) vertical
    ! average value of any model field, x, is calculated by the equation:
    !
    ! x_avg|_z = ( INT(z_bot:z_top) x rho_ds dz )
    !            / ( INT(z_bot:z_top) rho_ds dz );
    !
    ! where z_bot is the bottom of the vertical domain, and z_top is the top of
    ! the vertical domain.
    !
    ! This calculation is done slightly differently depending on whether x is a
    ! thermodynamic-level or a momentum-level variable.
    !
    ! Thermodynamic-level computation:
    
    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given thermodynamic level, x and rho_ds are
    ! both thermodynamic-level variables, and delta_z(k) = zm(k) - zm(k-1).  The
    ! indices k_bot and k_top are the indices of the respective lower and upper
    ! thermodynamic levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) thermodynamic level is below ground (or below the
    ! official lower boundary at the first momentum level), so it should not
    ! count in a vertical average, whether that vertical average is used for
    ! the hole-filling scheme or for statistical purposes. Begin no lower
    ! than level k=2, which is the first thermodynamic level above ground (or
    ! above the model lower boundary).
    !
    ! For cases where hole-filling over the entire (global) vertical domain
    ! is desired, or where statistics over the entire (global) vertical
    ! domain are desired, the lower (thermodynamic-level) index of k = 2 and
    ! the upper (thermodynamic-level) index of k = gr%nz, means that the
    ! overall vertical domain will be gr%zm(1,gr%nz) - gr%zm(1,1).
    !
    !
    ! Momentum-level computation:
    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given momentum level, x and rho_ds are both
    ! momentum-level variables, and delta_z(k) = zt(k+1) - zt(k).  The indices
    ! k_bot and k_top are the indices of the respective lower and upper momentum
    ! levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) momentum level is right at ground level (or right at
    ! the official lower boundary).  The momentum level variables that call
    ! the hole-filling scheme have set values at the surface (or lower
    ! boundary), and those set values should not be changed.  Therefore, the
    ! vertical average (for purposes of hole-filling) should not include the
    ! surface level (or lower boundary level).  For hole-filling purposes,
    ! begin no lower than level k=2, which is the second momentum level above
    ! ground (or above the model lower boundary).  Likewise, the value at the
    ! model upper boundary (k=gr%nz) is also set for momentum level
    ! variables.  That value should also not be changed.
    !
    ! However, this function is also used to keep track (for statistical
    ! purposes) of the vertical average of certain variables.  In that case,
    ! the vertical average needs to be taken over the entire vertical domain
    ! (level 1 to level gr%nz).
    !
    !
    ! In both the thermodynamic-level computation and the momentum-level
    ! computation, the numerator integral is divided by the denominator integral
    ! in order to find the average value (over the vertical domain) of x.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      total_idx ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds, & ! Dry, static density on either thermodynamic or momentum levels    [kg/m^3]
      field,  & ! The field (e.g. wp2) to be vertically averaged                    [Units vary]
      dz  ! Reciprocal of thermodynamic or momentum level thickness           [1/m]
                ! depending on whether we're on zt or zm grid.
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = 1.

    ! Output variable
    real( kind = core_rknd ) :: & 
      vertical_avg  ! Vertical average of field    [Units of field]

    ! Local variables
    real( kind = core_rknd ) :: & 
      numer_integral, & ! Integral in the numerator (see description)
      denom_integral    ! Integral in the denominator (see description)
      

    integer :: k

    !-----------------------------------------------------------------------
    
    ! Initialize variable
    numer_integral = 0.0_core_rknd
    denom_integral = 0.0_core_rknd

    ! Compute the numerator and denominator integral.
    ! Multiply rho_ds at level k by the level thickness
    ! at level k.  Then, sum over all vertical levels.
    do k=1, total_idx

        numer_integral = numer_integral + rho_ds(k) * dz(k) * field(k)
        denom_integral = denom_integral + rho_ds(k) * dz(k)

    end do

    ! Find the vertical average of 'field'.
    vertical_avg = numer_integral / denom_integral
    !vertical_avg = sum( rho_ds(:) * dz(:) * field(:) ) / sum( rho_ds(:) * dz(:) )

    return
  end function vertical_avg

  !=============================================================================
  function vertical_integral( total_idx, rho_ds, &
                                       field, dz )

    ! Description:
    ! Computes the vertical integral. rho_ds, field, and dz must all be
    ! of size total_idx and should all start at the same index.
    ! 
    
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      dz         ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = k_start.

    ! Local variables
    real( kind = core_rknd ) :: &
      vertical_integral ! Integral in the numerator (see description)

    !-----------------------------------------------------------------------

    !  Assertion checks: that k_start <= gr%nz - 1
    !                    that k_end   >= 2
    !                    that k_start <= k_end


    ! Initializing vertical_integral to avoid a compiler warning.
    vertical_integral = 0.0_core_rknd

    ! Compute the integral.
    ! Multiply the field at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    ! Note:  The values of the field and rho_ds are passed into this function
    !        so that field(1) and rho_ds(1) are actually the field and rho_ds
    !        at level k_start.
    vertical_integral = sum( field * rho_ds * dz )

    !print *, vertical_integral

    return
  end function vertical_integral


end module advance_helper_module
