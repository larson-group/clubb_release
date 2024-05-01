!-----------------------------------------------------------------------
! $Id: mixing_length.F90 8664 2018-05-10 20:21:35Z huebler@uwm.edu $
!===============================================================================
module mixing_length

  implicit none

  private ! Default Scope

  public :: calc_Lscale_directly,  &
            diagnose_Lscale_from_tau

  contains

  !=============================================================================
  subroutine compute_mixing_length( nz, ngrdcol, gr, thvm, thlm, &
                                    rtm, em, Lscale_max, p_in_Pa, &
                                    exner, thv_ds, mu, lmin, &
                                    saturation_formula, &
                                    l_implemented, &
                                    stats_metadata, &
                                    Lscale, Lscale_up, Lscale_down )

    ! Description:
    !   Larson's 5th moist, nonlocal length scale
    !
    ! References:
    !   Section 3b ( /Eddy length formulation/ ) of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !
    ! Notes:
    !
    !   The equation for the rate of change of theta_l and r_t of the parcel with
    !   respect to height, due to entrainment, is:
    !
    !           d(thl_par)/dz = - mu * ( thl_parcel - thl_environment );
    !
    !           d(rt_par)/dz = - mu * ( rt_parcel - rt_environment );
    !
    !   where mu is the entrainment rate,
    !   such that:
    !
    !           mu = (1/m)*(dm/dz);
    !
    !   where m is the mass of the parcel.  The value of mu is set to be a
    !   constant.
    !
    !   The differential equations are solved for given the boundary condition
    !   and given the fact that the value of thl_environment and rt_environment
    !   are treated as changing linearly for a parcel of air from one grid level
    !   to the next.
    !
    !   For the special case where entrainment rate, mu, is set to 0,
    !   thl_parcel and rt_parcel remain constant
    !
    !
    !   The equation for Lscale_up is:
    !
    !       INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
    !
    !   and for Lscale_down
    !
    !       INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
    !
    !   where thv_par is theta_v of the parcel, thvm is the mean
    !   environmental value of theta_v, z_i is the altitude that the parcel
    !   started from, and em is the mean value of TKE at
    !   altitude z_i (which gives the parcel its initial boost)
    !
    !   The increment of CAPE (convective air potential energy) for any two
    !   successive vertical levels is:
    !
    !       Upwards:
    !           CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
    !
    !       Downwards:
    !           CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz
    !
    !   Thus, the derivative of CAPE with respect to height is:
    !
    !           dCAPE/dz = g * ( thv_par - thvm ) / thvm.
    !
    !   A purely trapezoidal rule is used between levels, and is considered
    !   to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be
    !   of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
    !   where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )
    !
    !   The integral is evaluated to find the CAPE increment between two
    !   successive vertical levels.  The result either adds to or depletes
    !   from the total amount of energy that keeps the parcel ascending/descending.
    !
    !
    ! IMPORTANT NOTE:
    !   This subroutine has been optimized by adding precalculations, rearranging
    !   equations to avoid divides, and modifying the algorithm entirely.
    !       -Gunther Huebler
    !
    !   The algorithm previously used looped over every grid level, following a
    !   a parcel up from its initial grid level to its max. The very nature of this
    !   algorithm is an N^2
    !--------------------------------------------------------------------------------

    ! mu = (1/M) dM/dz > 0.  mu=0 for no entrainment.
    ! Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
    ! When mu was fixed, we used the value mu = 6.e-4

    use constants_clubb, only:  &  ! Variable(s)
        Cp,             & ! Dry air specific heat at constant pressure [J/kg/K]
        Rd,             & ! Dry air gas constant                       [J/kg/K]
        ep,             & ! Rd / Rv                                    [-]
        ep1,            & ! (1-ep)/ep                                  [-]
        ep2,            & ! 1/ep                                       [-]
        Lv,             & ! Latent heat of vaporiztion                 [J/kg/K]
        grav,           & ! Gravitational acceleration                 [m/s^2]
        fstderr,        &
        zero_threshold, &
        eps,            &
        one_half,       &
        one,            &
        two,            &
        zero

    use grid_class, only:  &
        grid, & ! Type
        zm2zt ! Procedure(s)

    use numerical_check, only:  &
        length_check ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use saturation, only:  &
        sat_mixrat_liq ! Procedure(s)
        
    use stats_variables, only: & 
        stats_metadata_type

    implicit none

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  &
      zlmin = 0.1_core_rknd, & ! Minimum value for Lscale [m]
      Lscale_sfclyr_depth = 500._core_rknd ! [m]

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol  

    type (grid), target, intent(in) :: &
      gr
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      thvm,    & ! Virtual potential temp. on themodynamic level  [K]
      thlm,    & ! Liquid potential temp. on themodynamic level   [K]
      rtm,     & ! Total water mixing ratio on themodynamic level [kg/kg]
      em,      & ! em = 3/2 * w'^2; on momentum level             [m^2/s^2]
      exner,   & ! Exner function on thermodynamic level          [-]
      p_in_Pa, & ! Pressure on thermodynamic level                [Pa]
      thv_ds     ! Dry, base-state theta_v on thermodynamic level [K]
    ! Note:  thv_ds used as a reference theta_l here

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      Lscale_max ! Maximum allowable value for Lscale             [m]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      mu      ! mu Fractional extrainment rate per unit altitude  [1/m]
      
    real( kind = core_rknd ), intent(in) :: &
      lmin    ! CLUBB tunable parameter lmin

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, intent(in) :: &
      l_implemented ! Flag for CLUBB being implemented in a larger model

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------------- Local Variables ---------------------------------

    integer :: i, j, k, start_index

    real( kind = core_rknd ) :: tke, CAPE_incr

    real( kind = core_rknd ) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    ! Temporary 2D arrays to store calculations to speed runtime
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
        exp_mu_dzm, &
        invrs_dzm_on_mu, &
        grav_on_thvm, &
        Lv_coef, &
        entrain_coef, &
        thl_par_j_precalc, &
        rt_par_j_precalc, &
        tl_par_1, &
        rt_par_1, &
        rsatl_par_1, &
        thl_par_1, &
        dCAPE_dz_1, &
        s_par_1, &
        rc_par_1, &
        CAPE_incr_1, &
        thv_par_1

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
        tke_i

    ! Minimum value for Lscale that will taper off with height
    real( kind = core_rknd ) :: lminh

    ! Parcel quantities at grid level j
    real( kind = core_rknd ) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

    ! Used in latent heating calculation
    real( kind = core_rknd ) :: tl_par_j, rsatl_par_j, s_par_j

    ! Variables to make L nonlocal
    real( kind = core_rknd ) :: Lscale_up_max_alt, Lscale_down_min_alt

    ! Variables used to precalculate values
    real( kind = core_rknd ) :: &
        Lv2_coef, &
        tl_par_j_sqd, &
        invrs_dCAPE_diff, &
        invrs_Lscale_sfclyr_depth

    ! --------------------------------- Begin Code ---------------------------------

    !$acc enter data create( exp_mu_dzm, invrs_dzm_on_mu, grav_on_thvm, Lv_coef, &
    !$acc                    entrain_coef, thl_par_j_precalc, rt_par_j_precalc, &
    !$acc                    tl_par_1, rt_par_1, rsatl_par_1, thl_par_1, dCAPE_dz_1, &
    !$acc                    s_par_1, rc_par_1, CAPE_incr_1, thv_par_1, tke_i )
 
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      if ( abs(mu(i)) < eps ) then
        err_code = clubb_fatal_error
        print *, "mu = ", mu(i)
      end if
    end do
    !$acc end parallel loop

    if ( err_code == clubb_fatal_error ) then
      write(fstderr,*) "Entrainment rate mu cannot be 0"
      error stop "Fatal error in subroutine compute_mixing_length"
    end if

    ! Calculate initial turbulent kinetic energy for each grid level
    tke_i = zm2zt( nz, ngrdcol, gr, em )
 
    ! Initialize arrays and precalculate values for computational efficiency
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 1, nz

        ! Initialize up and down arrays
        Lscale_up(i,k) = zlmin
        Lscale_down(i,k) = zlmin

        ! Precalculate values to avoid unnecessary calculations later
        exp_mu_dzm(i,k) = exp( -mu(i) * gr%dzm(i,k) )
        invrs_dzm_on_mu(i,k) = ( gr%invrs_dzm(i,k) ) / mu(i)
        grav_on_thvm(i,k) = grav / thvm(i,k)
        Lv_coef(i,k) = Lv / ( exner(i,k) * cp ) - ep2 * thv_ds(i,k)
        entrain_coef(i,k) = ( one - exp_mu_dzm(i,k) ) * invrs_dzm_on_mu(i,k)

      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! Avoid uninitialized memory (these values are not used in Lscale)
      Lscale_up(i,1)   = zero
      Lscale_down(i,1) = zero
    end do
    !$acc end parallel loop

    ! Precalculations of single values to avoid unnecessary calculations later
    Lv2_coef = ep * Lv**2 / ( Rd * cp )
    invrs_Lscale_sfclyr_depth = one / Lscale_sfclyr_depth


    ! ---------------- Upwards Length Scale Calculation ----------------

    ! Precalculate values for upward Lscale, these are useful only if a parcel can rise
    ! more than one level. They are used in the equations that calculate thl and rt
    ! recursively for a parcel as it ascends

    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol  
      do j = 2, nz-1

        thl_par_j_precalc(i,j) = thlm(i,j) - thlm(i,j-1) * exp_mu_dzm(i,j-1)  &
                               - ( thlm(i,j) - thlm(i,j-1) ) * entrain_coef(i,j-1)

        rt_par_j_precalc(i,j) = rtm(i,j) - rtm(i,j-1) * exp_mu_dzm(i,j-1)  &
                              - ( rtm(i,j) - rtm(i,j-1) ) * entrain_coef(i,j-1)
      end do
    end do
    !$acc end parallel loop

    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations.

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
     do j = 3, nz

        thl_par_1(i,j) = thlm(i,j) - ( thlm(i,j) - thlm(i,j-1) ) * entrain_coef(i,j-1)

        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)

        rt_par_1(i,j) = rtm(i,j) - ( rtm(i,j) - rtm(i,j-1) ) * entrain_coef(i,j-1)

      end do
    end do
    !$acc end parallel loop


    ! Caclculate initial rsatl for parcels at each grid level

    ! The entire pressure and temperature arrays are passed as 
    ! argument and the sub-arrays are choosen using 
    ! start_index. This workaround is used to solve 
    ! subarray issues with OpenACC.
    ! rsatl_par_1(i,3:) = sat_mixrat_liq_acc( nz-2, ngrdcol, p_in_Pa(i,3:), tl_par_1(i,3:) )
    ! since subarray 3:, the start_index is 3 and it is an optional argument
    start_index = 3
    rsatl_par_1 = sat_mixrat_liq( nz, ngrdcol, p_in_Pa, tl_par_1, saturation_formula, start_index )
    
    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do j = 3, nz

        tl_par_j_sqd = tl_par_1(i,j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
        !
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )

        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )

        ! theta_v of entraining parcel at grid level j
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) + Lv_coef(i,j) * rc_par_1(i,j)


        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) * gr%dzm(i,j-1)

      end do


      ! Calculate Lscale_up for each grid level. If the TKE from a parcel has not been completely
      ! exhausted by the initial change then continue the exhaustion calculations here for a single
      ! grid level at a time until the TKE is exhausted.

      Lscale_up_max_alt = zero    ! Set initial max value for Lscale_up to 0
      do k = 2, nz-2

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i,k) + CAPE_incr_1(i,k+1) > zero ) then

          ! Calculate new TKE for parcel
          tke = tke_i(i,k) + CAPE_incr_1(i,k+1)

          ! Set j to 2 levels above current Lscale_up level, this is because we've already
          ! determined that the parcel can rise at least 1 full level
          j = k + 2

          ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
          thl_par_j = thl_par_1(i,k+1)
          rt_par_j  = rt_par_1(i,k+1)
          dCAPE_dz_j_minus_1 = dCAPE_dz_1(i,k+1)


          ! Continue change in TKE calculations until it is exhausted or the max grid
          ! level has been reached. j is the next grid level above the level that can
          ! be reached for a parcel starting at level k. If TKE is exhausted in this loop
          ! that means the parcel starting at k cannot reach level j, but has reached j-1
          do while ( j < nz )

            ! thl, rt of parcel are conserved except for entrainment
            !
            ! The values of thl_env and rt_env are treated as changing linearly for a parcel
            ! of air ascending from level j-1 to level j

            ! theta_l of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
            thl_par_j = thl_par_j_precalc(i,j) + thl_par_j * exp_mu_dzm(i,j-1)


            ! r_t of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
            rt_par_j = rt_par_j_precalc(i,j) + rt_par_j * exp_mu_dzm(i,j-1)


            ! Include effects of latent heating on Lscale_up 6/12/00
            ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
            ! Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(i,j)

            rsatl_par_j = sat_mixrat_liq( p_in_Pa(i,j), tl_par_j, saturation_formula )

            tl_par_j_sqd = tl_par_j**2

            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )
            ! and SD's beta (eqn. 8),
            !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
            !
            ! Simplified by multiplying top and bottom by tl^2 to avoid a
            ! divide and precalculating ep * Lv**2 / ( Rd * cp )
            s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

            rc_par_j = max( s_par_j, zero_threshold )

            ! theta_v of entraining parcel at grid level j
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j  &
                        + Lv_coef(i,j) * rc_par_j

            ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )

            ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
            ! Trapezoidal estimate between grid levels j and j-1
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * gr%dzm(i,j-1)

            ! Exit loop early if tke has been exhaused between level j and j+1
            if ( tke + CAPE_incr <= zero ) then
                exit
            end if

            ! Save previous dCAPE value for next cycle
            dCAPE_dz_j_minus_1 = dCAPE_dz_j

            ! Caclulate new TKE and increment j
            tke = tke + CAPE_incr
            j = j + 1

          enddo


          ! Add full grid level thickness for each grid level that was passed without the TKE
          ! being exhausted, difference between starting level (k) and last level passed (j-1)
          Lscale_up(i,k) = Lscale_up(i,k) + gr%zt(i,j-1) - gr%zt(i,k)


          if ( j < nz ) then

            ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
            ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * 2 <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps ) then

              ! Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0
              ! Find the remaining distance z - z_0 that it takes to
              ! exhaust the remaining TKE

              Lscale_up(i,k) = Lscale_up(i,k) + ( - tke / dCAPE_dz_j )

            else

              ! Case used for most scenarios where dCAPE/dz|_(z_1) /= dCAPE/dz|_(z_0)
              ! Find the remaining distance z - z_0 that it takes to exhaust the
              ! remaining TKE (tke_i), using the quadratic formula (only the
              ! negative (-) root works in this scenario).
              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )

              Lscale_up(i,k) = Lscale_up(i,k) &
                             - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff * gr%dzm(i,j-1)  &
                             - sqrt( dCAPE_dz_j_minus_1**2 &
                                      - two * tke * gr%invrs_dzm(i,j-1) &
                                        * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                               * invrs_dCAPE_diff  * gr%dzm(i,j-1)
            endif

          end if

        else    ! TKE for parcel at level (k) was exhaused before one full grid level

          ! Find the remaining distance z - z_0 that it takes to exhaust the
          ! remaining TKE (tke_i), using the quadratic formula. Simplified
          ! since dCAPE_dz_j_minus_1 = 0.0
          Lscale_up(i,k) = Lscale_up(i,k) - sqrt( - two * tke_i(i,k) &
                                                * gr%dzm(i,k) * dCAPE_dz_1(i,k+1) ) &
                                        / dCAPE_dz_1(i,k+1)
        endif


        ! If a parcel at a previous grid level can rise past the parcel at this grid level
        ! then this one should also be able to rise up to that height. This feature insures
        ! that the profile of Lscale_up will be smooth, thus reducing numerical instability.
        if ( gr%zt(i,k) + Lscale_up(i,k) < Lscale_up_max_alt ) then

            ! A lower starting parcel can ascend higher than this one, set height to the max
            ! that any lower starting parcel can ascend to
            Lscale_up(i,k) = Lscale_up_max_alt - gr%zt(i,k)
        else

            ! This parcel can ascend higher than any below it, save final height
            Lscale_up_max_alt = Lscale_up(i,k) + gr%zt(i,k)
        end if


      end do
    end do
    !$acc end parallel loop

    ! ---------------- Downwards Length Scale Calculation ----------------

    ! Precalculate values for downward Lscale, these are useful only if a parcel can descend
    ! more than one level. They are used in the equations that calculate thl and rt
    ! recursively for a parcel as it descends
    !$acc parallel loop gang vector collapse(2) default(present)    
    do i = 1, ngrdcol
      do j = 2, nz-1

        thl_par_j_precalc(i,j) = thlm(i,j) - thlm(i,j+1) * exp_mu_dzm(i,j)  &
                               - ( thlm(i,j) - thlm(i,j+1) ) * entrain_coef(i,j)

        rt_par_j_precalc(i,j) = rtm(i,j) - rtm(i,j+1) * exp_mu_dzm(i,j)  &
                              - ( rtm(i,j) - rtm(i,j+1) ) * entrain_coef(i,j)
      end do
    end do
    !$acc end parallel loop

    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations.

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    !$acc parallel loop gang vector collapse(2) default(present)    
    do i = 1, ngrdcol
      do j = 2, nz-1

        thl_par_1(i,j) = thlm(i,j) - ( thlm(i,j) - thlm(i,j+1) )  * entrain_coef(i,j)

        tl_par_1(i,j) = thl_par_1(i,j) * exner(i,j)

        rt_par_1(i,j) = rtm(i,j) - ( rtm(i,j) - rtm(i,j+1) ) * entrain_coef(i,j)

      end do
    end do
    !$acc end parallel loop

    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental

    ! The entire pressure and temperature arrays are passed as 
    ! argument and the sub-arrays are choosen using 
    ! start_index. This workaround is used to solve 
    ! subarray issues with OpenACC.
    ! rsatl_par_1(i,2:) = sat_mixrat_liq_acc( nz-1, p_in_Pa(i,2:), tl_par_1(i,2:) )
    ! since subarray 2:, the start_index is 2 and it is an optional argument
    start_index = 2
    rsatl_par_1 = sat_mixrat_liq( nz, ngrdcol, p_in_Pa, tl_par_1, saturation_formula, start_index )

    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do j = 2, nz-1

        tl_par_j_sqd = tl_par_1(i,j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
        !
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )
        s_par_1(i,j) = ( rt_par_1(i,j) - rsatl_par_1(i,j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(i,j) )

        rc_par_1(i,j) = max( s_par_1(i,j), zero_threshold )

        ! theta_v of entraining parcel at grid level j
        thv_par_1(i,j) = thl_par_1(i,j) + ep1 * thv_ds(i,j) * rt_par_1(i,j) + Lv_coef(i,j) * rc_par_1(i,j)

        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(i,j) = grav_on_thvm(i,j) * ( thv_par_1(i,j) - thvm(i,j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(i,j) = one_half * dCAPE_dz_1(i,j) * gr%dzm(i,j)

      end do


      ! Calculate Lscale_down for each grid level. If the TKE from a parcel has not been completely
      ! exhausted by the initial change then continue the exhaustion calculations here for a single
      ! grid level at a time until the TKE is exhausted.

      Lscale_down_min_alt = gr%zt(i,nz)  ! Set initial min value for Lscale_down to max zt
      do k = nz, 3, -1

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i,k) - CAPE_incr_1(i,k-1) > zero ) then

          ! Calculate new TKE for parcel
          tke = tke_i(i,k) - CAPE_incr_1(i,k-1)

          ! Set j to 2 levels below current Lscale_down level, this is because we've already
          ! determined that the parcel can descend at least 1 full level
          j = k - 2

          ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
          thl_par_j = thl_par_1(i,k-1)
          rt_par_j = rt_par_1(i,k-1)
          dCAPE_dz_j_plus_1 = dCAPE_dz_1(i,k-1)


          ! Continue change in TKE calculations until it is exhausted or the min grid
          ! level has been reached. j is the next grid level below the level that can
          ! be reached for a parcel starting at level k. If TKE is exhausted in this loop
          ! that means the parcel starting at k cannot sink to level j, but can sink to j+1
          do while ( j >= 2 )

            ! thl, rt of parcel are conserved except for entrainment
            !
            ! The values of thl_env and rt_env are treated as changing linearly for a parcel
            ! of air descending from level j to level j-1

            ! theta_l of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
            thl_par_j = thl_par_j_precalc(i,j) + thl_par_j * exp_mu_dzm(i,j)


            ! r_t of the parcel starting at grid level k, and currenly
            ! at grid level j
            !
            ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
            rt_par_j = rt_par_j_precalc(i,j) + rt_par_j * exp_mu_dzm(i,j)


            ! Include effects of latent heating on Lscale_up 6/12/00
            ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
            ! Probably should use properties of bump 1 in Gaussian, not mean!!!

            tl_par_j = thl_par_j*exner(i,j)

            rsatl_par_j = sat_mixrat_liq( p_in_Pa(i,j), tl_par_j, saturation_formula )

            tl_par_j_sqd = tl_par_j**2

            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
            !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )
            ! and SD's beta (eqn. 8),
            !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
            !
            ! Simplified by multiplying top and bottom by tl^2 to avoid a
            ! divide and precalculating ep * Lv**2 / ( Rd * cp )
            s_par_j = (rt_par_j - rsatl_par_j) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

            rc_par_j = max( s_par_j, zero_threshold )

            ! theta_v of entraining parcel at grid level j
            thv_par_j = thl_par_j + ep1 * thv_ds(i,j) * rt_par_j + Lv_coef(i,j) * rc_par_j

            ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
            dCAPE_dz_j = grav_on_thvm(i,j) * ( thv_par_j - thvm(i,j) )

            ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
            ! Trapezoidal estimate between grid levels j+1 and j
            CAPE_incr = one_half * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * gr%dzm(i,j)

            ! Exit loop early if tke has been exhaused between level j+1 and j
            if ( tke - CAPE_incr <= zero ) then
              exit
            endif

            ! Save previous dCAPE value for next cycle
            dCAPE_dz_j_plus_1 = dCAPE_dz_j

            ! Caclulate new TKE and increment j
            tke = tke - CAPE_incr
            j = j - 1

          enddo

          ! Add full grid level thickness for each grid level that was passed without the TKE
          ! being exhausted, difference between starting level (k) and last level passed (j+1)
          Lscale_down(i,k) = Lscale_down(i,k) + gr%zt(i,k) - gr%zt(i,j+1)


          if ( j >= 2 ) then

            ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
            ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

            if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * 2 <= &
                 abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps ) then

              ! Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0
              ! Find the remaining distance z_0 - z that it takes to
              ! exhaust the remaining TKE

              Lscale_down(i,k) = Lscale_down(i,k) + ( tke / dCAPE_dz_j )

            else

              ! Case used for most scenarios where dCAPE/dz|_(z_(-1)) /= dCAPE/dz|_(z_0)
              ! Find the remaining distance z_0 - z that it takes to exhaust the
              ! remaining TKE (tke_i), using the quadratic formula (only the
              ! negative (-) root works in this scenario) -- however, the
              ! negative (-) root is divided by another negative (-) factor,
              ! which results in an overall plus (+) sign in front of the
              ! square root term in the equation below).
              invrs_dCAPE_diff = one / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )

              Lscale_down(i,k) = Lscale_down(i,k) &
                               - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff * gr%dzm(i,j)  &
                               + sqrt( dCAPE_dz_j_plus_1**2 &
                                       + two * tke * gr%invrs_dzm(i,j)  &
                                         * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                                 * invrs_dCAPE_diff * gr%dzm(i,j)
            endif

          end if

        else    ! TKE for parcel at level (k) was exhaused before one full grid level

          ! Find the remaining distance z_0 - z that it takes to exhaust the
          ! remaining TKE (tke_i), using the quadratic formula. Simplified
          ! since dCAPE_dz_j_plus_1 = 0.0
          Lscale_down(i,k) = Lscale_down(i,k) + sqrt( two * tke_i(i,k) &
                                                  * gr%dzm(i,k-1) * dCAPE_dz_1(i,k-1) ) &
                                            / dCAPE_dz_1(i,k-1)
        endif

        ! If a parcel at a previous grid level can descend past the parcel at this grid level
        ! then this one should also be able to descend down to that height. This feature insures
        ! that the profile of Lscale_down will be smooth, thus reducing numerical instability.
        if ( gr%zt(i,k) - Lscale_down(i,k) > Lscale_down_min_alt ) then
          Lscale_down(i,k) = gr%zt(i,k) - Lscale_down_min_alt
        else
          Lscale_down_min_alt = gr%zt(i,k) - Lscale_down(i,k)
        end if

      end do
    end do
    !$acc end parallel loop

      ! ---------------- Final Lscale Calculation ----------------

    !$acc parallel loop gang vector default(present) 
    do i = 1, ngrdcol
      do k = 2, nz, 1

        ! Make lminh a linear function starting at value lmin at the bottom
        ! and going to zero at 500 meters in altitude.
        if( l_implemented ) then

          ! Within a host model, increase mixing length in 500 m layer above *ground*
          lminh = max( zero_threshold, Lscale_sfclyr_depth - ( gr%zt(i,k) - gr%zm(i,1) ) ) &
                  * lmin * invrs_Lscale_sfclyr_depth
        else

          ! In standalone mode, increase mixing length in 500 m layer above *mean sea level*
          lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i,k) ) &
                  * lmin * invrs_Lscale_sfclyr_depth
        end if

        Lscale_up(i,k)    = max( lminh, Lscale_up(i,k) )
        Lscale_down(i,k)  = max( lminh, Lscale_down(i,k) )

        ! When L is large, turbulence is weakly damped
        ! When L is small, turbulence is strongly damped
        ! Use a geometric mean to determine final Lscale so that L tends to become small
        ! if either Lscale_up or Lscale_down becomes small.
        Lscale(i,k) = sqrt( Lscale_up(i,k) * Lscale_down(i,k) )

      enddo

      ! Set the value of Lscale at the upper and lower boundaries.
      Lscale(i,1) = Lscale(i,2)
      Lscale(i,nz) = Lscale(i,nz-1)

      ! Vince Larson limited Lscale to allow host
      ! model to take over deep convection.  13 Feb 2008.
      Lscale(i,:) = min( Lscale(i,:), Lscale_max(i) )
      
    end do
    !$acc end parallel loop

    ! Ensure that no Lscale values are NaN
    if ( clubb_at_least_debug_level( 1 ) ) then

      !$acc update host( Lscale, Lscale_up, Lscale_down, &
      !$acc              thvm, thlm, rtm, em, exner, p_in_Pa, thv_ds )

      do i = 1, ngrdcol
        call length_check( nz, Lscale(i,:), Lscale_up(i,:), Lscale_down(i,:) ) ! intent(in)
      end do

      if ( err_code == clubb_fatal_error ) then

        write(fstderr,*) "Errors in compute_mixing_length subroutine"

        write(fstderr,*) "Intent(in)"

        write(fstderr,*) "thvm = ", thvm
        write(fstderr,*) "thlm = ", thlm
        write(fstderr,*) "rtm = ", rtm
        write(fstderr,*) "em = ", em
        write(fstderr,*) "exner = ", exner
        write(fstderr,*) "p_in_Pa = ", p_in_Pa
        write(fstderr,*) "thv_ds = ", thv_ds

        write(fstderr,*) "Intent(out)"

        write(fstderr,*) "Lscale = ", Lscale
        write(fstderr,*) "Lscale_up = ", Lscale_up
        write(fstderr,*) "Lscale_down = ", Lscale_down

      endif ! Fatal error

    end if

    !$acc exit data delete( exp_mu_dzm, invrs_dzm_on_mu, grav_on_thvm, Lv_coef, &
    !$acc                   entrain_coef, thl_par_j_precalc, rt_par_j_precalc, &
    !$acc                   tl_par_1, rt_par_1, rsatl_par_1, thl_par_1, dCAPE_dz_1, &
    !$acc                   s_par_1, rc_par_1, CAPE_incr_1, thv_par_1, tke_i )

    return

  end subroutine compute_mixing_length

!===============================================================================
  subroutine calc_Lscale_directly ( ngrdcol, nz, gr, &
                                    l_implemented, p_in_Pa, exner, rtm,    &
                                    thlm, thvm, newmu, rtp2, thlp2, rtpthlp, &
                                    pdf_params, em, thv_ds_zt, Lscale_max, lmin, &
                                    clubb_params, &
                                    saturation_formula, &
                                    l_Lscale_plume_centered, &
                                    stats_metadata, &
                                    stats_zt, &
                                    Lscale, Lscale_up, Lscale_down)

    use constants_clubb, only: &
        thl_tol,      &
        rt_tol,       &
        one_half,     &
        one_third,    &
        one,          &
        three,        &
        unused_var

    use parameter_indices, only: &
        nparams, &
        iLscale_mu_coef, &
        iLscale_pert_coef

    use grid_class, only: &
        grid ! Type

    use clubb_precision, only: &
        core_rknd

    use stats_variables, only: &
        stats_metadata_type

    use pdf_parameter_module, only: &
        pdf_parameter

    use stats_type_utilities, only:   &
        stat_update_var

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use constants_clubb, only:  &
        fstderr  ! Variable(s)

    use stats_type, only: &
        stats ! Type

    implicit none

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: &
      gr

    logical, intent(in) ::  &
      l_implemented ! True if CLUBB is being run within a large-scale hostmodel,
                    !   rather than a standalone single-column model.

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      rtp2,      &
      thlp2,     &
      rtpthlp,   &
      thlm,      &
      thvm,      &
      rtm,       &
      em,        &
      p_in_Pa,   & ! Air pressure (thermodynamic levels)       [Pa]
      exner,     &
      thv_ds_zt

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      newmu, &
      Lscale_max
      
    real( kind = core_rknd ), intent(in) ::  &
      lmin

    type (pdf_parameter), intent(in) :: &
      pdf_params    ! PDF Parameters  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, intent(in) :: &
      l_Lscale_plume_centered    ! Alternate that uses the PDF to compute the perturbed values

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------------- InOut Variables ---------------------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt

    !--------------------------------- Output Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    !--------------------------------- Local Variables ---------------------------------
    integer :: k, i

    logical, parameter :: &
      l_avg_Lscale = .false.  ! Lscale is calculated in subroutine compute_mixing_length
                              ! if l_avg_Lscale is true, compute_mixing_length is called two additional
                              ! times with
                              ! perturbed values of rtm and thlm.  An average value of Lscale
                              ! from the three calls to compute_mixing_length is then calculated.
                              ! This reduces temporal noise in RICO, BOMEX, LBA, and other cases.

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sign_rtpthlp,         & ! Sign of the covariance rtpthlp       [-]
      Lscale_pert_1, Lscale_pert_2, & ! For avg. calculation of Lscale  [m]
      thlm_pert_1, thlm_pert_2, &     ! For avg. calculation of Lscale  [K]
      rtm_pert_1, rtm_pert_2,   &     ! For avg. calculation of Lscale  [kg/kg]
      thlm_pert_pos_rt, thlm_pert_neg_rt, &     ! For avg. calculation of Lscale [K]
      rtm_pert_pos_rt, rtm_pert_neg_rt     ! For avg. calculation of Lscale [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      mu_pert_1, mu_pert_2, &
      mu_pert_pos_rt, mu_pert_neg_rt  ! For l_Lscale_plume_centered

    real( kind = core_rknd ) :: &
      Lscale_pert_coef

    !Lscale_weight Uncomment this if you need to use this vairable at some
    !point.

    !--------------------------------- Begin Code ---------------------------------

    !$acc enter data create( sign_rtpthlp, Lscale_pert_1, Lscale_pert_2, &
    !$acc                    thlm_pert_1, thlm_pert_2, rtm_pert_1, rtm_pert_2, &
    !$acc                    thlm_pert_pos_rt, thlm_pert_neg_rt, rtm_pert_pos_rt, &
    !$acc                    rtm_pert_neg_rt, &
    !$acc                    mu_pert_1, mu_pert_2, mu_pert_pos_rt, mu_pert_neg_rt )

    if ( clubb_at_least_debug_level( 0 ) ) then

      if ( l_Lscale_plume_centered .and. .not. l_avg_Lscale ) then
        write(fstderr,*) "l_Lscale_plume_centered requires l_avg_Lscale"
        write(fstderr,*) "Fatal error in advance_clubb_core"
        err_code = clubb_fatal_error
        return
      end if

    end if

    if ( l_avg_Lscale .and. .not. l_Lscale_plume_centered ) then

      ! Call compute length two additional times with perturbed values
      ! of rtm and thlm so that an average value of Lscale may be calculated.

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz, 1
        do  i = 1, ngrdcol
          sign_rtpthlp(i,k) = sign( one, rtpthlp(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz, 1
        do  i = 1, ngrdcol
          rtm_pert_1(i,k)  = rtm(i,k) + clubb_params(i,iLscale_pert_coef) &
                                        * sqrt( max( rtp2(i,k), rt_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz, 1
        do  i = 1, ngrdcol
          thlm_pert_1(i,k) = thlm(i,k) + sign_rtpthlp(i,k) * clubb_params(i,iLscale_pert_coef) &
                                         * sqrt( max( thlp2(i,k), thl_tol**2 ) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do  i = 1, ngrdcol
        mu_pert_1(i)   = newmu(i) / clubb_params(i,iLscale_mu_coef)
      end do 
      !$acc end parallel loop

      call compute_mixing_length( nz, ngrdcol, gr, thvm, thlm_pert_1,  & ! In
                    rtm_pert_1, em, Lscale_max, p_in_Pa,               & ! In
                    exner, thv_ds_zt, mu_pert_1, lmin,                 & ! In 
                    saturation_formula,                                & ! In
                    l_implemented,                                     & ! In
                    stats_metadata,                                    & ! In
                    Lscale_pert_1, Lscale_up, Lscale_down )              ! Out


      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz, 1
        do  i = 1, ngrdcol
          rtm_pert_2(i,k)  = rtm(i,k) - clubb_params(i,iLscale_pert_coef) &
                                        * sqrt( max( rtp2(i,k), rt_tol**2 ) )
        end do
      end do
      !$acc end parallel loop
      
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz, 1
        do  i = 1, ngrdcol
          thlm_pert_2(i,k) = thlm(i,k) - sign_rtpthlp(i,k) * clubb_params(i,iLscale_pert_coef) &
                               * sqrt( max( thlp2(i,k), thl_tol**2 ) )
        end do
      end do
      !$acc end parallel loop
           
      !$acc parallel loop gang vector default(present) 
      do  i = 1, ngrdcol
        mu_pert_2(i)   = newmu(i) * clubb_params(i,iLscale_mu_coef)
      end do 
      !$acc end parallel loop         

      call compute_mixing_length( nz, ngrdcol, gr, thvm, thlm_pert_2, & ! In
                    rtm_pert_2, em, Lscale_max, p_in_Pa,              & ! In
                    exner, thv_ds_zt, mu_pert_2, lmin,                & ! In 
                    saturation_formula,                               & ! In
                    l_implemented,                                    & ! In
                    stats_metadata,                                   & ! In
                    Lscale_pert_2, Lscale_up, Lscale_down )             ! Out

    else if ( l_avg_Lscale .and. l_Lscale_plume_centered ) then

      ! Take the values of thl and rt based one 1st or 2nd plume

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          sign_rtpthlp(i,k) = sign( one, rtpthlp(i,k) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol

          Lscale_pert_coef = clubb_params(i,iLscale_pert_coef)

          if ( pdf_params%rt_1(i,k) > pdf_params%rt_2(i,k) ) then

            rtm_pert_pos_rt(i,k) = pdf_params%rt_1(i,k) &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1(i,k), rt_tol**2 ) )

            thlm_pert_pos_rt(i,k) = pdf_params%thl_1(i,k) + ( sign_rtpthlp(i,k) * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_1(i,k), thl_tol**2 ) ) )

            thlm_pert_neg_rt(i,k) = pdf_params%thl_2(i,k) - ( sign_rtpthlp(i,k) * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_2(i,k), thl_tol**2 ) ) )

            rtm_pert_neg_rt(i,k) = pdf_params%rt_2(i,k) &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2(i,k), rt_tol**2 ) )

            !Lscale_weight = pdf_params%mixt_frac(i,k)

          else

            rtm_pert_pos_rt(i,k) = pdf_params%rt_2(i,k) &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2(i,k), rt_tol**2 ) )

            thlm_pert_pos_rt(i,k) = pdf_params%thl_2(i,k) + ( sign_rtpthlp(i,k) * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_2(i,k), thl_tol**2 ) ) )

            thlm_pert_neg_rt(i,k) = pdf_params%thl_1(i,k) - ( sign_rtpthlp(i,k) * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_1(i,k), thl_tol**2 ) ) )

            rtm_pert_neg_rt(i,k) = pdf_params%rt_1(i,k) &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1(i,k), rt_tol**2 ) )

            !Lscale_weight = 1.0_core_rknd - pdf_params%mixt_frac(i,k)

          end if

        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present) 
      do i = 1, ngrdcol
        mu_pert_pos_rt(i) = newmu(i) / clubb_params(i,iLscale_mu_coef)
        mu_pert_neg_rt(i) = newmu(i) * clubb_params(i,iLscale_mu_coef)
      end do
      !$acc end parallel loop

      ! Call length with perturbed values of thl and rt
      call compute_mixing_length( nz, ngrdcol, gr, thvm, thlm_pert_pos_rt,  & ! In
                rtm_pert_pos_rt, em, Lscale_max, p_in_Pa,                   & ! In
                exner, thv_ds_zt, mu_pert_pos_rt, lmin,                     & ! In 
                saturation_formula,                                         & ! In
                l_implemented,                                              & ! In
                stats_metadata,                                             & ! In
                Lscale_pert_1, Lscale_up, Lscale_down )                       ! Out

      call compute_mixing_length( nz, ngrdcol, gr, thvm, thlm_pert_neg_rt,  & ! In
                rtm_pert_neg_rt, em, Lscale_max, p_in_Pa,                   & ! In
                exner, thv_ds_zt, mu_pert_neg_rt, lmin,                     & ! In 
                saturation_formula,                                         & ! In
                l_implemented,                                              & ! In
                stats_metadata,                                             & ! In
                Lscale_pert_2, Lscale_up, Lscale_down )                       ! Out
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Lscale_pert_1(i,k) = unused_var ! Undefined
          Lscale_pert_2(i,k) = unused_var ! Undefined
        end do
      end do
      !$acc end parallel loop
    end if ! l_avg_Lscale

    if ( stats_metadata%l_stats_samp ) then
      !$acc update host( Lscale_pert_1, Lscale_pert_2 )
      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%iLscale_pert_1, Lscale_pert_1(i,:), & ! intent(in)
                              stats_zt(i) )                       ! intent(inout)
        call stat_update_var( stats_metadata%iLscale_pert_2, Lscale_pert_2(i,:), & ! intent(in)
                              stats_zt(i) )                       ! intent(inout)
      end do
    end if ! stats_metadata%l_stats_samp


    ! ********** NOTE: **********
    ! This call to compute_mixing_length must be last.  Otherwise, the values
    ! of
    ! Lscale_up and Lscale_down in stats will be based on perturbation length
    ! scales
    ! rather than the mean length scale.

    ! Diagnose CLUBB's turbulent mixing length scale.
    call compute_mixing_length( nz, ngrdcol, gr, thvm, thlm,            & ! In
                          rtm, em, Lscale_max, p_in_Pa,                 & ! In
                          exner, thv_ds_zt, newmu, lmin,                & ! In 
                          saturation_formula,                           & ! In
                          l_implemented,                                & ! In
                          stats_metadata,                               & ! In
                          Lscale, Lscale_up, Lscale_down )                ! Out

    if ( l_avg_Lscale ) then
      if ( l_Lscale_plume_centered ) then
        ! Weighted average of mean, pert_1, & pert_2
        !       Lscale = 0.5_core_rknd * ( Lscale + Lscale_weight*Lscale_pert_1 &
        !                                  + (1.0_core_rknd-Lscale_weight)*Lscale_pert_2
        !                                  )
        ! Weighted average of just the perturbed values
        !       Lscale = Lscale_weight*Lscale_pert_1 +
        !       (1.0_core_rknd-Lscale_weight)*Lscale_pert_2

        ! Un-weighted average of just the perturbed values
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            Lscale(i,k) = one_half *( Lscale_pert_1(i,k) + Lscale_pert_2(i,k) )
          end do
        end do
        !$acc end parallel loop
      else
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            Lscale(i,k) = one_third * ( Lscale(i,k) + Lscale_pert_1(i,k) + Lscale_pert_2(i,k) )
          end do
        end do
        !$acc end parallel loop
      end if
    end if

    !$acc exit data delete( sign_rtpthlp, Lscale_pert_1, Lscale_pert_2, &
    !$acc                   thlm_pert_1, thlm_pert_2, rtm_pert_1, rtm_pert_2, &
    !$acc                   thlm_pert_pos_rt, thlm_pert_neg_rt, rtm_pert_pos_rt, &
    !$acc                   rtm_pert_neg_rt, &
    !$acc                   mu_pert_1, mu_pert_2, mu_pert_pos_rt, mu_pert_neg_rt )

   return
   
 end subroutine  calc_Lscale_directly



!===============================================================================

 subroutine diagnose_Lscale_from_tau( nz, ngrdcol, gr, &
                        upwp_sfc, vpwp_sfc, um, vm, & !intent in
                        exner, p_in_Pa, & !intent in
                        rtm, thlm, thvm, & !intent in
                        rcm, ice_supersat_frac, &! intent in
                        em, sqrt_em_zt, & ! intent in
                        ufmin, tau_const, & ! intent in
                        sfc_elevation, Lscale_max, & ! intent in
                        clubb_params, & ! intent in
                        stats_metadata, & ! intent in
                        saturation_formula, & ! intent in
                        l_e3sm_config, & ! intent in
                        l_brunt_vaisala_freq_moist, & !intent in
                        l_use_thvm_in_bv_freq, &! intent in
                        l_smooth_Heaviside_tau_wpxp, & ! intent in
                        l_modify_limiters_for_cnvg_test, & ! intent in
                        stats_zm, & ! intent inout
                        brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, & ! intent out
                        brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, & ! intent out
                        Ri_zm, & ! intent out
                        invrs_tau_zt, invrs_tau_zm, & ! intent out
                        invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd, & ! intent out
                        invrs_tau_shear, invrs_tau_N2_iso, & ! intent out
                        invrs_tau_wp2_zm, invrs_tau_xp2_zm, & ! intent out
                        invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm, & ! intent out
                        tau_max_zm, tau_max_zt, tau_zm, tau_zt, & !intent out
                        Lscale, Lscale_up, Lscale_down)! intent out
! Description:
!     Diagnose inverse damping time scales (invrs_tau_...) and turbulent mixing length (Lscale)
! References:
!     Guo et al.(2021, JAMES)
!--------------------------------------------------------------------------------------------------

    use advance_helper_module, only: &
        calc_brunt_vaisala_freq_sqd, &
        smooth_heaviside_peskin, &
        smooth_min, smooth_max

    use stats_type_utilities, only:   &
        stat_update_var ! Procedure

    use constants_clubb, only: &
        one_fourth,     &
        one_half,       &
        vonk,           &
        zero,           &
        one,            & 
        two,            &
        em_min,         &
        zero_threshold, &
        eps

    use grid_class, only: &
        grid, & ! Type
        zt2zm, &
        zm2zt, &
        zm2zt2zm, &
        zt2zm2zt, &
        ddzt

    use clubb_precision, only: &
        core_rknd

    use parameter_indices, only: &
        nparams,                     & ! Variable(s)
        iC_invrs_tau_bkgnd,          &
        iC_invrs_tau_shear,          &
        iC_invrs_tau_sfc,            &
        iC_invrs_tau_N2,             &
        iC_invrs_tau_N2_wp2 ,        &
        iC_invrs_tau_N2_wpxp,        &
        iC_invrs_tau_N2_xp2,         &
        iC_invrs_tau_wpxp_N2_thresh, &
        iC_invrs_tau_N2_clear_wp3,   &
        iC_invrs_tau_wpxp_Ri,        &
        ialtitude_threshold,         &
        ibv_efold,                   &
        iwpxp_Ri_exp,                &
        iz_displace

    use error_code, only: &
      err_code, &
      clubb_fatal_error, &
      clubb_at_least_debug_level

    use stats_variables, only: &
      stats_metadata_type

    use stats_type, only: stats ! Type

    implicit none

    !--------------------------------- Input Variables ---------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      upwp_sfc,      &
      vpwp_sfc
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      um,                &
      vm,                &
      exner,             &
      p_in_Pa,           &
      rtm,               &
      thlm,              &
      thvm,              &
      rcm,               &
      ice_supersat_frac, &
      em,                &
      sqrt_em_zt

    real(kind = core_rknd), intent(in) :: &
      ufmin,         &
      tau_const
      
    real(kind = core_rknd), dimension(ngrdcol), intent(in) :: &
      sfc_elevation, &
      Lscale_max

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]
 
    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, intent(in) :: &
      l_e3sm_config,              &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq, &      ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_smooth_Heaviside_tau_wpxp   ! Use the smoothed Heaviside 'Peskin' function
                                    ! to compute invrs_tau_wpxp_zm

    ! Flag to activate modifications on limiters for convergence test 
    ! (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
    ! (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
    ! (reduce threshold on limiters for Ri_zm in mixing_length.F90)
    logical, intent(in) :: &
      l_modify_limiters_for_cnvg_test

    !--------------------------- Input/Output Variables ---------------------------
    type (stats), target, intent(inout), dimension(ngrdcol) :: &
      stats_zm

    !--------------------------------- Output Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      brunt_vaisala_freq_sqd,       &
      brunt_vaisala_freq_sqd_mixed, &
      brunt_vaisala_freq_sqd_dry,   &
      brunt_vaisala_freq_sqd_moist, &
      Ri_zm,                        &
      invrs_tau_zt,                 &
      invrs_tau_zm,                 &
      invrs_tau_sfc,                &
      invrs_tau_no_N2_zm,           &
      invrs_tau_bkgnd,              &
      invrs_tau_shear,              &
      invrs_tau_N2_iso,             &
      invrs_tau_wp2_zm,             &
      invrs_tau_xp2_zm,             &
      invrs_tau_wp3_zm,             &
      invrs_tau_wp3_zt,             &
      invrs_tau_wpxp_zm,            &
      tau_max_zm,                   &
      tau_max_zt,                   &
      tau_zm,                       &
      tau_zt,                       &
      Lscale,                       &
      Lscale_up,                    &
      Lscale_down

    !--------------------------------- Local Variables ---------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_freq_pos,               &
      brunt_vaisala_freq_sqd_smth,  & ! smoothed Buoyancy frequency squared, N^2     [s^-2]
      brunt_freq_out_cloud,         &
      smooth_thlm,                  & 
      bvf_thresh,                   & ! temporatory array  
      H_invrs_tau_wpxp_N2             ! Heaviside function for clippings of invrs_tau_wpxp_N2

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      ustar
      
    real( kind = core_rknd ) :: &
      C_invrs_tau_N2,             &
      C_invrs_tau_N2_wp2

    real( kind = core_rknd ), parameter :: &
      min_max_smth_mag = 1.0e-9_core_rknd, &  ! "base" smoothing magnitude before scaling 
                                              ! for the respective data structure. See
                                              ! https://github.com/larson-group/clubb/issues/965#issuecomment-1119816722
                                              ! for a plot on how output behaves with varying min_max_smth_mag
      heaviside_smth_range = 1.0e-0_core_rknd ! range where Heaviside function is smoothed
   
    logical, parameter :: l_smooth_min_max = .false.  ! whether to apply smooth min and max functions

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      ddzt_um, &
      ddzt_vm, &
      ddzt_umvm_sqd, &
      ddzt_umvm_sqd_clipped, &
      norm_ddzt_umvm, &
      smooth_norm_ddzt_umvm, &
      brunt_vaisala_freq_clipped, &
      ice_supersat_frac_zm, &
      invrs_tau_shear_smooth, &
      Ri_zm_clipped, &
      Ri_zm_smooth, &
      em_clipped, &
      tau_zm_unclipped, & 
      tau_zt_unclipped, &
      tmp_calc, &
      tmp_calc_max, &
      tmp_calc_min_max

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      tmp_calc_ngrdcol

    integer :: i, k

    !--------------------------------- Begin Code ---------------------------------

    !$acc enter data create( brunt_freq_pos, brunt_vaisala_freq_sqd_smth, brunt_freq_out_cloud, &
    !$acc                    smooth_thlm, bvf_thresh, H_invrs_tau_wpxp_N2, ustar, &
    !$acc                    ddzt_um, ddzt_vm, norm_ddzt_umvm, smooth_norm_ddzt_umvm, &
    !$acc                    brunt_vaisala_freq_clipped, &
    !$acc                    ice_supersat_frac_zm, invrs_tau_shear_smooth, &
    !$acc                    ddzt_umvm_sqd, tmp_calc_ngrdcol )

    !$acc enter data if( l_smooth_min_max .or. l_modify_limiters_for_cnvg_test ) &
    !$acc            create( Ri_zm_clipped, ddzt_umvm_sqd_clipped, &
    !$acc                    tau_zm_unclipped, tau_zt_unclipped, Ri_zm_smooth, em_clipped, &
    !$acc                    tmp_calc, tmp_calc_max, tmp_calc_min_max )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      if ( gr%zm(i,1) - sfc_elevation(i) + clubb_params(i,iz_displace) < eps ) then
        err_code = clubb_fatal_error
      end if
    end do
    !$acc end parallel loop

    if ( clubb_at_least_debug_level(0) ) then
      if ( err_code == clubb_fatal_error ) then
        error stop  "Lowest zm grid level is below ground in CLUBB."
      end if
    end if

    ! Smooth thlm by interpolating to zm then back to zt
    smooth_thlm = zt2zm2zt( nz, ngrdcol, gr, thlm )

    call calc_brunt_vaisala_freq_sqd( nz, ngrdcol, gr, smooth_thlm, & ! intent(in)
                                      exner, rtm, rcm, p_in_Pa, thvm, & ! intent(in)
                                      ice_supersat_frac, & ! intent(in)
                                      saturation_formula, & ! intent(in)
                                      l_brunt_vaisala_freq_moist, & ! intent(in)
                                      l_use_thvm_in_bv_freq, & ! intent(in)
                                      clubb_params(:,ibv_efold), & ! intent(in)
                                      brunt_vaisala_freq_sqd, & ! intent(out)
                                      brunt_vaisala_freq_sqd_mixed,& ! intent(out)
                                      brunt_vaisala_freq_sqd_dry, & ! intent(out)
                                      brunt_vaisala_freq_sqd_moist ) ! intent(out)

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      tmp_calc_ngrdcol(i) = ( upwp_sfc(i)**2 + vpwp_sfc(i)**2 )**one_fourth
    end do
    !$acc end parallel loop

    if ( l_smooth_min_max ) then

      ! tmp_calc_ngrdcol used here because openacc has an issue with
      !  a variable being both input and output of a function
      ustar = smooth_max( ngrdcol, tmp_calc_ngrdcol, ufmin, min_max_smth_mag )

    else 

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ustar(i) = max( tmp_calc_ngrdcol(i), ufmin )
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        invrs_tau_bkgnd(i,k) = clubb_params(i,iC_invrs_tau_bkgnd) / tau_const
      end do
    end do
    !$acc end parallel loop

    ddzt_um = ddzt( nz, ngrdcol, gr, um )
    ddzt_vm = ddzt( nz, ngrdcol, gr, vm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        ddzt_umvm_sqd(i,k) = ddzt_um(i,k)**2 + ddzt_vm(i,k)**2
        norm_ddzt_umvm(i,k) = sqrt( ddzt_umvm_sqd(i,k) )
      end do
    end do
    !$acc end parallel loop

    smooth_norm_ddzt_umvm = zm2zt2zm( nz, ngrdcol, gr, norm_ddzt_umvm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        invrs_tau_shear_smooth(i,k) = clubb_params(i,iC_invrs_tau_shear) *  smooth_norm_ddzt_umvm(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Enforce that invrs_tau_shear is positive
    invrs_tau_shear = smooth_max( nz, ngrdcol, invrs_tau_shear_smooth, &
                                  zero_threshold, min_max_smth_mag )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        invrs_tau_sfc(i,k) = clubb_params(i,iC_invrs_tau_sfc) &
                             * ( ustar(i) / vonk ) &
                               / ( gr%zm(i,k) - sfc_elevation(i) + clubb_params(i,iz_displace) )
         !C_invrs_tau_sfc * ( wp2 / vonk /ustar ) / ( gr%zm(1,:) -sfc_elevation + z_displace )
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        invrs_tau_no_N2_zm(i,k) = invrs_tau_bkgnd(i,k) + invrs_tau_sfc(i,k) + invrs_tau_shear(i,k)
      end do
    end do
    !$acc end parallel loop

    !The min function below smooths the slope discontinuity in brunt freq
    !  and thereby allows tau to remain large in Sc layers in which thlm may
    !  be slightly stably stratified.
    if ( l_modify_limiters_for_cnvg_test ) then 

      !Remove the limiters to improve the solution convergence
      brunt_vaisala_freq_sqd_smth = zm2zt2zm( nz,ngrdcol,gr, brunt_vaisala_freq_sqd_mixed )

    else  ! default method  

      if ( l_smooth_min_max ) then

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            tmp_calc(i,k) = 1.e8_core_rknd * abs(brunt_vaisala_freq_sqd_mixed(i,k))**3
          end do
        end do
        !$acc end parallel loop

        brunt_vaisala_freq_clipped = smooth_min( nz, ngrdcol, &
                                                 brunt_vaisala_freq_sqd_mixed, &
                                                 tmp_calc, &
                                                 1.0e-4_core_rknd * min_max_smth_mag)

        brunt_vaisala_freq_sqd_smth = zm2zt2zm( nz, ngrdcol, gr, brunt_vaisala_freq_clipped )

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            brunt_vaisala_freq_clipped(i,k) = min( brunt_vaisala_freq_sqd_mixed(i,k), &
                                                   1.e8_core_rknd * abs(brunt_vaisala_freq_sqd_mixed(i,k))**3)
          end do
        end do
        !$acc end parallel loop

        brunt_vaisala_freq_sqd_smth = zm2zt2zm( nz, ngrdcol, gr, brunt_vaisala_freq_clipped )

      end if

    end if

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( brunt_vaisala_freq_sqd_smth )

      do i = 1, ngrdcol
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_smth, brunt_vaisala_freq_sqd_smth(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
      end do
    end if

    if ( l_modify_limiters_for_cnvg_test ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Ri_zm_clipped(i,k) = max( 0.0_core_rknd, brunt_vaisala_freq_sqd_smth(i,k) ) &
                                  / max( ddzt_umvm_sqd(i,k), 1.0e-12_core_rknd )
        end do
      end do
      !$acc end parallel loop

      Ri_zm = zm2zt2zm( nz, ngrdcol, gr, Ri_zm_clipped )

    else ! default method 

      if ( l_smooth_min_max ) then

        brunt_vaisala_freq_clipped = smooth_max( nz, ngrdcol, 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_smth, &
                                                 1.0e-4_core_rknd * min_max_smth_mag )

        ddzt_umvm_sqd_clipped = smooth_max( nz, ngrdcol, ddzt_umvm_sqd, 1.0e-7_core_rknd, &
                                        1.0e-6_core_rknd * min_max_smth_mag )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            Ri_zm(i,k) = brunt_vaisala_freq_clipped(i,k) / ddzt_umvm_sqd_clipped(i,k)
          end do
        end do
        !$acc end parallel loop

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            Ri_zm(i,k) = max( 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_smth(i,k) ) &
                            / max( ddzt_umvm_sqd(i,k), 1.0e-7_core_rknd )
          end do
        end do
        !$acc end parallel loop

      end if

    end if

    if ( l_smooth_min_max ) then

      brunt_vaisala_freq_clipped = smooth_max( nz, ngrdcol, zero_threshold, &
                                               brunt_vaisala_freq_sqd_smth, &
                                               1.0e-4_core_rknd * min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_freq_pos(i,k) = sqrt( brunt_vaisala_freq_clipped(i,k) )
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_freq_pos(i,k) = sqrt( max( zero_threshold, brunt_vaisala_freq_sqd_smth(i,k) ) )
        end do
      end do
      !$acc end parallel loop

    end if

    ice_supersat_frac_zm = zt2zm( nz, ngrdcol, gr, ice_supersat_frac, zero_threshold )

    if ( l_smooth_min_max ) then

      ! roll this back as well once checks have passed
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tmp_calc(i,k) = one - ice_supersat_frac_zm(i,k) / 0.001_core_rknd
        end do
      end do
      !$acc end parallel loop

      tmp_calc_max = smooth_max( nz, ngrdcol, zero_threshold, tmp_calc, &
                                 min_max_smth_mag)

      tmp_calc_min_max = smooth_min( nz, ngrdcol, one, tmp_calc_max, &
                                     min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_freq_out_cloud(i,k) =  brunt_freq_pos(i,k) * tmp_calc_min_max(i,k)
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          brunt_freq_out_cloud(i,k) &
            = brunt_freq_pos(i,k) &
                * min(one, max(zero_threshold, &
                               one - ( ( ice_supersat_frac_zm(i,k) / 0.001_core_rknd) )))
        end do
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        if ( gr%zt(i,k) < clubb_params(i,ialtitude_threshold) ) then
          brunt_freq_out_cloud(i,k) = zero
        end if
      end do
    end do
    !$acc end parallel loop

    ! Write both bv extra terms to invrs_taus to disk
    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( brunt_freq_pos, brunt_freq_out_cloud )

      do i = 1, ngrdcol
        call stat_update_var(stats_metadata%ibrunt_freq_pos, brunt_freq_pos(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
        call stat_update_var(stats_metadata%ibrunt_freq_out_cloud, brunt_freq_out_cloud(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
      end do
    end if

    ! This time scale is used optionally for the return-to-isotropy term. It
    ! omits invrs_tau_sfc based on the rationale that the isotropization
    ! rate shouldn't be enhanced near the ground.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        C_invrs_tau_N2     = clubb_params(i,iC_invrs_tau_N2)
        C_invrs_tau_N2_wp2 = clubb_params(i,iC_invrs_tau_N2_wp2)

        invrs_tau_N2_iso(i,k) = invrs_tau_bkgnd(i,k) + invrs_tau_shear(i,k) &
                                + C_invrs_tau_N2_wp2 * brunt_freq_pos(i,k)

        invrs_tau_wp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                C_invrs_tau_N2 * brunt_freq_pos(i,k) + &
                                C_invrs_tau_N2_wp2 * brunt_freq_out_cloud(i,k)

        invrs_tau_zm(i,k) = invrs_tau_no_N2_zm(i,k) + C_invrs_tau_N2 * brunt_freq_pos(i,k)
      end do
    end do
    !$acc end parallel loop


    if ( l_e3sm_config ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          invrs_tau_zm(i,k) = one_half * invrs_tau_zm(i,k)
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          invrs_tau_xp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) &
                                  + clubb_params(i,iC_invrs_tau_N2_xp2) * brunt_freq_pos(i,k) & ! 0
                                  + clubb_params(i,iC_invrs_tau_sfc) * two &
                                  * sqrt(em(i,k)) &
                                  / ( gr%zm(i,k) - sfc_elevation(i) + clubb_params(i,iz_displace) )  ! small
        end do
      end do
      !$acc end parallel loop

      if ( l_smooth_min_max ) then

        brunt_vaisala_freq_clipped = smooth_max( nz, ngrdcol, 1.0e-7_core_rknd, &
                                                 brunt_vaisala_freq_sqd_smth, &
                                                 1.0e-4_core_rknd * min_max_smth_mag )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            tmp_calc(i,k) = sqrt( ddzt_umvm_sqd(i,k) / brunt_vaisala_freq_clipped(i,k) )
          end do
        end do
        !$acc end parallel loop

        tmp_calc_max = smooth_max( nz, ngrdcol, tmp_calc, &
                                   0.3_core_rknd, 0.3_core_rknd * min_max_smth_mag )

        tmp_calc_min_max = smooth_min( nz, ngrdcol, tmp_calc_max, &
                                       one, min_max_smth_mag )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            invrs_tau_xp2_zm(i,k) =  tmp_calc_min_max(i,k) * invrs_tau_xp2_zm(i,k)
          end do
        end do
        !$acc end parallel loop

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nz
          do i = 1, ngrdcol
            invrs_tau_xp2_zm(i,k) &
              = min( max( sqrt( ddzt_umvm_sqd(i,k) &
                                / max( 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_smth(i,k) ) ), &
                          0.3_core_rknd ), one ) &
                     * invrs_tau_xp2_zm(i,k)
          end do
        end do
        !$acc end parallel loop

      end if

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          invrs_tau_wpxp_zm(i,k) = two * invrs_tau_zm(i,k) &
                                   + clubb_params(i,iC_invrs_tau_N2_wpxp) * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

    else ! l_e3sm_config = false

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          invrs_tau_xp2_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                  clubb_params(i,iC_invrs_tau_N2) * brunt_freq_pos(i,k) + &
                                  clubb_params(i,iC_invrs_tau_N2_xp2) * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

      ice_supersat_frac_zm = zt2zm( nz, ngrdcol, gr, ice_supersat_frac, zero_threshold )

!      !$acc parallel loop gang vector collapse(2) default(present)
!      do k = 1, nz
!        do i = 1, ngrdcol
!          if ( ice_supersat_frac_zm(i,k) <= 0.01_core_rknd &
!               .and. invrs_tau_xp2_zm(i,k)  >= 0.003_core_rknd ) then
!            invrs_tau_xp2_zm(i,k) = 0.003_core_rknd
!          end if
!        end do
!      end do
!      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          invrs_tau_wpxp_zm(i,k) = invrs_tau_no_N2_zm(i,k) + &
                                   clubb_params(i,iC_invrs_tau_N2) * brunt_freq_pos(i,k) + &
                                   clubb_params(i,iC_invrs_tau_N2_wpxp) * brunt_freq_out_cloud(i,k)
        end do
      end do
      !$acc end parallel loop

    end if ! l_e3sm_config

    if ( l_smooth_Heaviside_tau_wpxp ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          bvf_thresh(i,k) = brunt_vaisala_freq_sqd_smth(i,k) &
                            / clubb_params(i,iC_invrs_tau_wpxp_N2_thresh) - one
          end do
      end do
      !$acc end parallel loop

      H_invrs_tau_wpxp_N2 = smooth_heaviside_peskin( nz, ngrdcol, bvf_thresh, heaviside_smth_range )

    else ! l_smooth_Heaviside_tau_wpxp = .false.

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          if ( brunt_vaisala_freq_sqd_smth(i,k) > clubb_params(i,iC_invrs_tau_wpxp_N2_thresh) ) then
            H_invrs_tau_wpxp_N2(i,k) = one
          else
            H_invrs_tau_wpxp_N2(i,k) = zero
          end if
        end do
      end do
      !$acc end parallel loop

    end if ! l_smooth_Heaviside_tau_wpxp

    if ( l_smooth_min_max ) then

      Ri_zm_smooth = smooth_max( nz, ngrdcol, Ri_zm, zero, &
                                  12.0_core_rknd * min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tmp_calc(i,k) = clubb_params(i,iC_invrs_tau_wpxp_Ri) * Ri_zm_smooth(i,k)**clubb_params(i,iwpxp_Ri_exp) 
        end do 
      end do
      !$acc end parallel loop

      Ri_zm_smooth = smooth_min( nz, ngrdcol, tmp_calc, &
                                 12.0_core_rknd, 12.0_core_rknd * min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol

          if ( gr%zt(i,k) > clubb_params(i,ialtitude_threshold) ) then
             invrs_tau_wpxp_zm(i,k) = invrs_tau_wpxp_zm(i,k) &
                                      * ( one + H_invrs_tau_wpxp_N2(i,k) &
                                          * Ri_zm_smooth(i,k) )

          end if
        end do 
      end do
      !$acc end parallel loop

    else ! l_smooth_min_max

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          if ( gr%zt(i,k) > clubb_params(i,ialtitude_threshold) ) then
             invrs_tau_wpxp_zm(i,k) = invrs_tau_wpxp_zm(i,k) &
                                      * ( one  + H_invrs_tau_wpxp_N2(i,k) &
                                      * min( clubb_params(i,iC_invrs_tau_wpxp_Ri) &
                                      * max( Ri_zm(i,k), zero)**clubb_params(i,iwpxp_Ri_exp), 12.0_core_rknd ) )
          end if
        end do
      end do
      !$acc end parallel loop

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        invrs_tau_wp3_zm(i,k) = invrs_tau_wp2_zm(i,k) &
                                + clubb_params(i,iC_invrs_tau_N2_clear_wp3) * brunt_freq_out_cloud(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Calculate the maximum allowable value of time-scale tau,
    ! which depends of the value of Lscale_max.
    if ( l_smooth_min_max ) then

      em_clipped = smooth_max( nz, ngrdcol, em, em_min, min_max_smth_mag )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tau_max_zt(i,k) = Lscale_max(i) / sqrt_em_zt(i,k)
          tau_max_zm(i,k) = Lscale_max(i) / sqrt( em_clipped(i,k) )
        end do
      end do
      !$acc end parallel loop

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tau_max_zt(i,k) = Lscale_max(i) / sqrt_em_zt(i,k)
          tau_max_zm(i,k) = Lscale_max(i) / sqrt( max( em(i,k), em_min ) )
        end do
      end do
      !$acc end parallel loop

    end if

    if ( l_smooth_min_max ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tau_zm_unclipped(i,k) = one / invrs_tau_zm(i,k)
        end do
      end do
      !$acc end parallel loop

      tau_zm = smooth_min( nz, ngrdcol, tau_zm_unclipped, &
                           tau_max_zm, 1.0e3_core_rknd * min_max_smth_mag )

      tau_zt_unclipped = zm2zt( nz, ngrdcol, gr, tau_zm )

      tau_zt = smooth_min( nz, ngrdcol, tau_zt_unclipped, tau_max_zt, 1.0e3_core_rknd * min_max_smth_mag )

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tau_zm(i,k) = min( one / invrs_tau_zm(i,k), tau_max_zm(i,k) )
        end do
      end do
      !$acc end parallel loop

      tau_zt = zm2zt( nz, ngrdcol, gr, tau_zm )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          tau_zt(i,k) = min( tau_zt(i,k), tau_max_zt(i,k) )
        end do
      end do
      !$acc end parallel loop

    end if

    invrs_tau_zt     = zm2zt( nz, ngrdcol, gr, invrs_tau_zm )
    invrs_tau_wp3_zt = zm2zt( nz, ngrdcol, gr, invrs_tau_wp3_zm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        Lscale(i,k) = tau_zt(i,k) * sqrt_em_zt(i,k)

        ! Lscale_up and Lscale_down aren't calculated with this option.
        ! They are set to 0 for stats output.
        Lscale_up(i,k) = zero
        Lscale_down(i,k) = zero

      end do
    end do
    !$acc end parallel loop

    !$acc exit data delete( brunt_freq_pos, brunt_vaisala_freq_sqd_smth, brunt_freq_out_cloud, &
    !$acc                   smooth_thlm, bvf_thresh, H_invrs_tau_wpxp_N2, ustar, &
    !$acc                   ddzt_um, ddzt_vm, norm_ddzt_umvm, smooth_norm_ddzt_umvm, &
    !$acc                   brunt_vaisala_freq_clipped, &
    !$acc                   ice_supersat_frac_zm, invrs_tau_shear_smooth, &
    !$acc                   ddzt_umvm_sqd, tmp_calc_ngrdcol )

    !$acc exit data if( l_smooth_min_max .or. l_modify_limiters_for_cnvg_test ) &
    !$acc           delete( Ri_zm_clipped, ddzt_umvm_sqd_clipped, &
    !$acc                   tau_zm_unclipped, tau_zt_unclipped, Ri_zm_smooth, em_clipped, &
    !$acc                   tmp_calc, tmp_calc_max, tmp_calc_min_max )

    return
    
  end subroutine diagnose_Lscale_from_tau

end module mixing_length
