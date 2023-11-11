!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module sfc_varnce_module

  implicit none

  private ! Default to private

  public :: calc_sfc_varnce

  contains

  !=============================================================================
  subroutine calc_sfc_varnce( nz, ngrdcol, gr, dt, sfc_elevation, & 
                              upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc, & 
                              um, vm, Lscale_up, wpsclrp_sfc, & 
                              lhs_splat_wp2, tau_zm, &
                              !wp2_splat_sfc, tau_zm_sfc, &
                              l_vary_convect_depth, &
                              clubb_params, &
                              stats_metadata, &
                              stats_zm, &
                              wp2, up2, vp2, & 
                              thlp2, rtp2, rtpthlp, & 
                              sclrp2, sclrprtp,  sclrpthlp )

    ! Description:
    ! This subroutine computes estimate of the surface thermodynamic and wind
    ! component second order moments.

    ! References:
    ! Andre, J. C., G. De Moor, P. Lacarrere, G. Therry, and R. Du Vachat, 1978.
    !   Modeling the 24-Hour Evolution of the Mean and Turbulent Structures of
    !   the Planetary Boundary Layer.  J. Atmos. Sci., 35, 1861 -- 1883.
    !-----------------------------------------------------------------------

    use grid_class, only: &
      grid

    use parameters_model, only:  & 
        T0 ! Variable(s)

    use constants_clubb, only: &
        four,       & ! Variable(s)
        two,        &
        one,        &
        two_thirds, &
        one_third,  &
        one_fourth, &
        zero,       &
        grav,       &
        eps,        &
        w_tol_sqd,  &
        thl_tol,    &
        rt_tol,     &
        max_mag_correlation_flux, &
        fstderr,    &
        wp2_max

    use parameters_model, only: & 
        sclr_dim  ! Variable(s)

    use parameter_indices, only: &
        nparams, &
        ia_const, &
        iup2_sfc_coef

    use numerical_check, only: & 
        sfc_varnce_check ! Procedure

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use array_index, only: &
        iisclr_rt, & ! Index for a scalar emulating rt
        iisclr_thl   ! Index for a scalar emulating thetal

    use stats_type, only: &
      stats ! Type

    use stats_type_utilities, only: & 
        stat_end_update_pt,   & ! Procedure(s)
        stat_begin_update_pt, &
        stat_update_var_pt

    use stats_variables, only: &
        stats_metadata_type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters

    ! Logical for Andre et al., 1978 parameterization.
    logical, parameter :: l_andre_1978 = .false.

    real( kind = core_rknd ), parameter ::  & 
      z_const = one, & ! Defined height of 1 meter                [m]
        ! Vince Larson increased ufmin to stabilize arm_97.  24 Jul 2007
        ! ufmin = 0.0001_core_rknd, &
      ufmin = 0.01_core_rknd, & ! Minimum allowable value of u*   [m/s]
        ! End Vince Larson's change.
        ! Vince Larson changed in order to make correlations between [-1,1].  31 Jan 2008.
        ! sclr_var_coef = 0.25_core_rknd, & ! This value is made up! - Vince Larson 12 Jul 2005
      sclr_var_coef = 0.4_core_rknd,  & ! This value is made up! - Vince Larson 12 Jul 2005
        ! End Vince Larson's change
        ! Vince Larson reduced surface spike in scalar variances associated
        ! w/ Andre et al. 1978 scheme
      reduce_coef   = 0.2_core_rknd

    !-------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ) ::  & 
      dt

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      sfc_elevation, &
      upwp_sfc,         & ! Surface u momentum flux, <u'w'>|_sfc   [m^2/s^2]
      vpwp_sfc,         & ! Surface v momentum flux, <v'w'>|_sfc   [m^2/s^2]
      wprtp_sfc           ! Surface moisture flux, <w'rt'>|_sfc    [kg/kg m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      wpthlp,        & ! Surface thetal flux, <w'thl'>|_sfc     [K m/s]
      um,            & ! Surface u wind component, <u>          [m/s]
      vm,            & ! Surface v wind component, <v>          [m/s]
      Lscale_up,     & ! Upward component of Lscale at surface  [m]
      lhs_splat_wp2, & 
      !wp2_splat,    & ! Tendency of <w'^2> due to splatting of eddies at zm(1) [m^2/s^3]
      tau_zm           ! Turbulent dissipation time at level zm(1)  [s]
      

    real( kind = core_rknd ), dimension(ngrdcol,sclr_dim), intent(in) ::  & 
      wpsclrp_sfc    ! Passive scalar flux, <w'sclr'>|_sfc   [units m/s]

    logical, intent(in) :: &
      l_vary_convect_depth

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !-------------------------- InOut Variables --------------------------
    type (stats), target, intent(inout), dimension(ngrdcol) :: &
      stats_zm

    !-------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) ::  & 
      wp2,     & ! Surface variance of w, <w'^2>|_sfc            [m^2/s^2]
      up2,     & ! Surface variance of u, <u'^2>|_sfc            [m^2/s^2]
      vp2,     & ! Surface variance of v, <v'^2>|_sfc            [m^2/s^2]
      thlp2,   & ! Surface variance of theta-l, <thl'^2>|_sfc    [K^2]
      rtp2,    & ! Surface variance of rt, <rt'^2>|_sfc          [(kg/kg)^2]
      rtpthlp    ! Surface covariance of rt and theta-l          [kg K/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(inout) ::  & 
      sclrp2,    & ! Surface variance of passive scalar            [units^2]
      sclrprtp,  & ! Surface covariance of pssv scalar and rt  [units kg/kg]
      sclrpthlp    ! Surface covariance of pssv scalar and theta-l [units K]

    !-------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol) ::  & 
      uf,               & ! Surface friction vel., u*, in older version   [m/s]
      depth_pos_wpthlp, & ! Thickness of the layer near the surface with wpthlp > 0 [m]
      min_wp2_sfc_val     ! Minimum value of wp2_sfc that guarantees  [m^2/s^2]
                          ! correlation of (w,rt) and (w,thl) is within (-1,1)
      
    ! Variables for Andre et al., 1978 parameterization.
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      um_sfc_sqd, & ! Surface value of <u>^2                           [m^2/s^2]
      vm_sfc_sqd, & ! Surface value of <v>^2                           [m^2/s^2]
      usp2_sfc,   & ! u_s (vector oriented w/ mean sfc. wind) variance [m^2/s^2]
      vsp2_sfc,   & ! v_s (vector perpen. to mean sfc. wind) variance  [m^2/s^2]
      ustar,      & ! Surface friction velocity, u*                             [m/s]
      zeta,       & ! Dimensionless height z_const/Lngth, where z_const = 1 m.  [-]
      wp2_splat_sfc_correction  ! Reduction in wp2_sfc due to splatting [m^2/s^2]

    real( kind = core_rknd ) ::  & 
      ustar2, & ! Square of surface friction velocity, u*       [m^2/s^2]
      wstar     ! Convective velocity, w*                       [m/s]

    real( kind = core_rknd ) :: &
      Lngth    ! Monin-Obukhov length [m]

    real( kind = core_rknd ) :: &
      up2_sfc_coef,   & ! CLUBB tunable parameter up2_sfc_coef   [-]
      a_const           ! Coefficient in front of wp2, up2, and vp2

    integer :: i, k, sclr ! Loop index

    !-------------------------- Begin Code --------------------------

    !$acc enter data create( uf, depth_pos_wpthlp, min_wp2_sfc_val, &
    !$acc                    um_sfc_sqd, vm_sfc_sqd, usp2_sfc, vsp2_sfc, &
    !$acc                    ustar, zeta, wp2_splat_sfc_correction )

    up2_sfc_coef = clubb_params(iup2_sfc_coef)

    ! Reflect surface varnce changes in budget
    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( wp2, up2, vp2, thlp2, rtp2, rtpthlp )

      do i = 1, ngrdcol
        call stat_begin_update_pt( stats_metadata%ithlp2_sf, 1,      & ! intent(in)
                                   thlp2(i,1) / dt,   & ! intent(in)
                                   stats_zm(i) )        ! intent(inout)

        call stat_begin_update_pt( stats_metadata%irtp2_sf, 1,       & ! intent(in)
                                   rtp2(i,1) / dt,    & ! intent(in)
                                   stats_zm(i) )        ! intent(inout)

        call stat_begin_update_pt( stats_metadata%irtpthlp_sf, 1,    & ! intent(in)
                                   rtpthlp(i,1) / dt, & ! intent(in)
                                   stats_zm(i) )        ! intent(inout)

        call stat_begin_update_pt( stats_metadata%iup2_sf, 1,        & ! intent(in)
                                   up2(i,1) / dt,     & ! intent(in)
                                   stats_zm(i) )        ! intent(inout)

        call stat_begin_update_pt( stats_metadata%ivp2_sf, 1,        & ! intent(in)
                                   vp2(i,1) / dt,     & ! intent(in)
                                   stats_zm(i) ) ! intent(inout)

        call stat_begin_update_pt( stats_metadata%iwp2_sf, 1,        & ! intent(in)
                                   wp2(i,1) / dt,     & ! intent(in)
                                   stats_zm(i) )        ! intent(inout)
      end do
    end if

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      ! Find thickness of layer near surface with positive heat flux.
      ! This is used when l_vary_convect_depth=.true. in order to determine wp2.
      if ( wpthlp(i,1) <= zero ) then
         depth_pos_wpthlp(i) = one ! When sfc heat flux is negative, set depth to 1 m.
      else ! When sfc heat flux is positive, march up sounding until wpthlp 1st becomes negative.
         k = 1
         do while ( wpthlp(i,k) > zero .and. (gr%zm(i,k)-sfc_elevation(i)) < 1000._core_rknd )
            k = k + 1
        end do
        depth_pos_wpthlp(i) = max( one, gr%zm(i,k)-sfc_elevation(i) )
      end if
    end do
    !$acc end parallel loop

    ! a_const used to be set here; now it is a tunable parameter
    !if ( .not. l_vary_convect_depth ) then
    !   a_const = 1.8_core_rknd
    !else
    !   a_const = 0.6_core_rknd
    !end if

    a_const = clubb_params(ia_const) 

    if ( l_andre_1978 ) then

      ! Calculate <u>^2 and <v>^2.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        um_sfc_sqd(i) = um(i,2)**2
        vm_sfc_sqd(i) = vm(i,2)**2
      end do
      !$acc end parallel loop

      ! Calculate surface friction velocity, u*.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ustar(i) = max( ( upwp_sfc(i)**2 + vpwp_sfc(i)**2 )**(one_fourth), ufmin )
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        if ( abs(wpthlp(i,1)) > eps) then

          ! Find Monin-Obukhov Length (Andre et al., 1978, p. 1866).
          Lngth = - ( ustar(i)**3 ) &
                    / ( 0.35_core_rknd * (one/T0) * grav * wpthlp(i,1) )

          ! Find the value of dimensionless height zeta
          ! (Andre et al., 1978, p. 1866).
          zeta(i) = z_const / Lngth

       else ! wpthlp = 0

          ! The value of Monin-Obukhov length is +inf when ustar < 0 and -inf
          ! when ustar > 0.  Either way, zeta = 0.
          zeta(i) = zero

        endif ! wpthlp /= 0
      end do
      !$acc end parallel loop

      ! Andre et al, 1978, Eq. 29.
      ! Notes:  1) "reduce_coef" is a reduction coefficient intended to make
      !            the values of rtp2, thlp2, and rtpthlp smaller at the
      !            surface.
      !         2) With the reduction coefficient having a value of 0.2, the
      !            surface correlations of both w & rt and w & thl have a value
      !            of about 0.845.  These correlations are greater if zeta < 0.
      !            The correlations have a value greater than 1 if
      !            zeta <= -0.212.
      !         3) The surface correlation of rt & thl is 1.
      ! Brian Griffin; February 2, 2008.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        if ( zeta(i) < zero ) then

          thlp2(i,1)   = reduce_coef  & 
                        * ( wpthlp(i,1)**2 / ustar(i)**2 ) & 
                        * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

          rtp2(i,1)    = reduce_coef  & 
                        * ( wprtp_sfc(i)**2 / ustar(i)**2 ) & 
                        * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

          rtpthlp(i,1) = reduce_coef  & 
                        * ( wprtp_sfc(i) * wpthlp(i,1) / ustar(i)**2 ) & 
                        * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

          wp2(i,1)    = ( ustar(i)**2 ) & 
                          * ( 1.75_core_rknd + two * (-zeta(i))**(two_thirds) )

       else

          thlp2(i,1)   = reduce_coef  & 
                        * four * ( wpthlp(i,1)**2 / ustar(i)**2 )

          rtp2(i,1)    = reduce_coef  & 
                        * four * ( wprtp_sfc(i)**2 / ustar(i)**2 )

          rtpthlp(i,1) = reduce_coef  & 
                        * four * ( wprtp_sfc(i) * wpthlp(i,1) / ustar(i)**2 )

          wp2(i,1)     = 1.75_core_rknd * ustar(i)**2

        endif
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        thlp2(i,1) = max( thl_tol**2, thlp2(i,1) )
        rtp2(i,1) = max( rt_tol**2, rtp2(i,1) )
      end do
      !$acc end parallel loop

      ! Calculate wstar following Andre et al., 1978, p. 1866.
      ! w* = ( ( 1 / T0 ) * g * <w'thl'>|_sfc * z_i )^(1/3);
      ! where z_i is the height of the mixed layer.  The value of CLUBB's
      ! upward component of mixing length, Lscale_up, at the surface will be
      ! used as z_i.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        wstar = ( (one/T0) * grav * wpthlp(i,1) * Lscale_up(i,2) )**(one_third)

        ! Andre et al., 1978, Eq. 29.
        ! Andre et al. (1978) defines horizontal wind surface variances in terms
        ! of orientation with the mean surface wind.  The vector u_s is the wind
        ! vector oriented with the mean surface wind.  The vector v_s is the wind
        ! vector oriented perpendicular to the mean surface wind.  Thus, <u_s> is
        ! equal to the mean surface wind (both in speed and direction), and <v_s>
        ! is 0.  Equation 29 gives the formula for the variance of u_s, which is
        ! <u_s'^2> (usp2_sfc in the code), and the formula for the variance of
        ! v_s, which is <v_s'^2> (vsp2_sfc in the code).
        if ( wpthlp(i,1) > zero ) then

          usp2_sfc(i) = four * ustar(i)**2 + 0.3_core_rknd * wstar**2

          vsp2_sfc(i) = 1.75_core_rknd * ustar(i)**2 + 0.3_core_rknd * wstar**2

        else

          usp2_sfc(i) = four * ustar(i)**2

          vsp2_sfc(i) = 1.75_core_rknd * ustar(i)**2

        endif
      end do
      !$acc end parallel loop

      ! Add effect of vertical compression of eddies on horizontal gustiness.
      ! First, ensure that wp2 does not make the correlation 
      !   of (w,thl) or (w,rt) outside (-1,1).
      ! Perhaps in the future we should also ensure that the correlations 
      !   of (w,u) and (w,v) are not outside (-1,1).
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        min_wp2_sfc_val(i) &
               = max( w_tol_sqd, &
                      wprtp_sfc(i)**2 / ( rtp2(i,1) * max_mag_correlation_flux**2 ), &
                      wpthlp(i,1)**2 / ( thlp2(i,1) * max_mag_correlation_flux**2 ) )
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        if ( wp2(i,1) - tau_zm(i,1) * lhs_splat_wp2(i,1) * wp2(i,1) &
             < min_wp2_sfc_val(i) ) then
        !if ( wp2(i,1) + tau_zm(i,1) * wp2_splat(i,1) < min_wp2_sfc_val(i) ) then 
                              ! splatting correction drives wp2 to overly small value
          wp2_splat_sfc_correction(i) = -wp2(i,1) + min_wp2_sfc_val(i)
          wp2(i,1) = min_wp2_sfc_val(i)
        else
          wp2_splat_sfc_correction(i) = tau_zm(i,1) * lhs_splat_wp2(i,1) * wp2(i,1)
          !wp2_splat_sfc_correction(i) = tau_zm(i,1) * wp2_splat(i,1)
          wp2(i,1) = wp2(i,1) + wp2_splat_sfc_correction(i)
        end if
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        usp2_sfc(i) = usp2_sfc(i) - 0.5_core_rknd * wp2_splat_sfc_correction(i)
        vsp2_sfc(i) = vsp2_sfc(i) - 0.5_core_rknd * wp2_splat_sfc_correction(i)
      end do
      !$acc end parallel loop

      ! Variance of u, <u'^2>, at the surface can be found from <u_s'^2>,
      ! <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
      !    <u'^2>|_sfc = <u_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
      !                  + <v_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
      ! where <u>^2 + <v>^2 /= 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        up2(i,1) = usp2_sfc(i) * ( um_sfc_sqd(i) / max( um_sfc_sqd(i) + vm_sfc_sqd(i) , eps ) )  &
                   + vsp2_sfc(i) * ( vm_sfc_sqd(i) / max( um_sfc_sqd(i) + vm_sfc_sqd(i) , eps ) )
      end do
      !$acc end parallel loop

      ! Variance of v, <v'^2>, at the surface can be found from <u_s'^2>,
      ! <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
      !    <v'^2>|_sfc = <v_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
      !                  + <u_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
      ! where <u>^2 + <v>^2 /= 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol 
        vp2(i,1) = vsp2_sfc(i) * ( um_sfc_sqd(i) / max( um_sfc_sqd(i) + vm_sfc_sqd(i) , eps ) )  &
                   + usp2_sfc(i) * ( vm_sfc_sqd(i) / max( um_sfc_sqd(i) + vm_sfc_sqd(i) , eps ) )
      end do
      !$acc end parallel loop

      ! Passive scalars
      if ( sclr_dim > 0 ) then

        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol

            ! Notes:  1) "reduce_coef" is a reduction coefficient intended to
            !            make the values of sclrprtp, sclrpthlp, and sclrp2
            !            smaller at the surface.
            !         2) With the reduction coefficient having a value of 0.2,
            !            the surface correlation of w & sclr has a value of
            !            about 0.845.  The correlation is greater if zeta < 0.
            !            The correlation has a value greater than 1 if
            !            zeta <= -0.212.
            !         3) The surface correlations of both rt & sclr and
            !            thl & sclr are 1.
            ! Brian Griffin; February 2, 2008.
            if ( zeta(i) < zero ) then

              sclrprtp(i,1,sclr)  & 
              = reduce_coef  & 
                * ( wpsclrp_sfc(i,sclr) * wprtp_sfc(i) / ustar(i)**2 ) & 
                * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

              sclrpthlp(i,1,sclr)  & 
              = reduce_coef  & 
                * ( wpsclrp_sfc(i,sclr) * wpthlp(i,1) / ustar(i)**2 ) & 
                * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

              sclrp2(i,1,sclr)  & 
              = reduce_coef   & 
                * ( wpsclrp_sfc(i,sclr)**2 / ustar(i)**2 ) & 
                * four * ( one - 8.3_core_rknd * zeta(i) )**(-two_thirds)

            else

              sclrprtp(i,1,sclr)  & 
              = reduce_coef  & 
                * four * ( wpsclrp_sfc(i,sclr) * wprtp_sfc(i) / ustar(i)**2 )

              sclrpthlp(i,1,sclr)  & 
              = reduce_coef  & 
                * four * ( wpsclrp_sfc(i,sclr) * wpthlp(i,1) / ustar(i)**2 )

              sclrp2(i,1,sclr)  & 
              = reduce_coef & 
                * four * ( wpsclrp_sfc(i,sclr)**2 / ustar(i)**2 )

            endif

          end do
        enddo ! i = 1, sclr_dim
        !$acc end parallel loop

      endif


    else ! Previous code.

      ! Compute ustar^2
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ustar2 = sqrt( upwp_sfc(i) * upwp_sfc(i) + vpwp_sfc(i) * vpwp_sfc(i) )

        ! Compute wstar following Andre et al., 1976
        if ( wpthlp(i,1) > zero ) then
          if ( .not. l_vary_convect_depth ) then
             wstar = ( one/T0 * grav * wpthlp(i,1) * z_const )**(one_third)
          else
             wstar = ( one/T0 * grav * wpthlp(i,1) * 0.2_core_rknd * depth_pos_wpthlp(i) )**(one_third)
          end if
        else
          wstar = zero
        endif

        ! Surface friction velocity following Andre et al. 1978
        if ( .not. l_vary_convect_depth ) then
          uf(i) = sqrt( ustar2 + 0.3_core_rknd * wstar * wstar )
        else
          uf(i) = sqrt( ustar2 + wstar * wstar )
        end if

        uf(i) = max( ufmin, uf(i) )
      end do
      !$acc end parallel loop

      ! Compute estimate for surface second order moments
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wp2(i,1) = a_const * uf(i)**2
        up2(i,1) = up2_sfc_coef * a_const * uf(i)**2  ! From Andre, et al. 1978
        vp2(i,1) = up2_sfc_coef * a_const * uf(i)**2  ! "  "
      end do
      !$acc end parallel loop

      ! Notes:  1) With "a" having a value of 1.8, the surface correlations of
      !            both w & rt and w & thl have a value of about 0.878.
      !         2) The surface correlation of rt & thl is 0.5.
      ! Brian Griffin; February 2, 2008.

      if ( .not. l_vary_convect_depth )  then
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          thlp2(i,1)   = 0.4_core_rknd * a_const * ( wpthlp(i,1) / uf(i) )**2
          rtp2(i,1)    = 0.4_core_rknd * a_const * ( wprtp_sfc(i) / uf(i) )**2
          rtpthlp(i,1) = 0.2_core_rknd * a_const * ( wpthlp(i,1) / uf(i) ) &
                                                 * ( wprtp_sfc(i) / uf(i) )
        end do
        !$acc end parallel loop
      else
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          thlp2(i,1)   = ( wpthlp(i,1) / uf(i) )**2 / ( max_mag_correlation_flux**2 * a_const )
          rtp2(i,1)    = ( wprtp_sfc(i) / uf(i) )**2 / ( max_mag_correlation_flux**2 * a_const )
          rtpthlp(i,1) = max_mag_correlation_flux * sqrt( thlp2(i,1) * rtp2(i,1) )
        end do
        !$acc end parallel loop
      end if

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        thlp2(i,1) = max( thl_tol**2, thlp2(i,1) )
        rtp2(i,1)  = max( rt_tol**2, rtp2(i,1) )
      end do
      !$acc end parallel loop

      ! Add effect of vertical compression of eddies on horizontal gustiness.
      ! First, ensure that wp2 does not make the correlation 
      !   of (w,thl) or (w,rt) outside (-1,1).
      ! Perhaps in the future we should also ensure that the correlations 
      !   of (w,u) and (w,v) are not outside (-1,1).
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        min_wp2_sfc_val(i) &
            = max( w_tol_sqd, &
                   wprtp_sfc(i)**2 / ( rtp2(i,1) * max_mag_correlation_flux**2 ), &
                   wpthlp(i,1)**2 / ( thlp2(i,1) * max_mag_correlation_flux**2 ) )
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        if ( wp2(i,1) - tau_zm(i,1) * lhs_splat_wp2(i,1) * wp2(i,1) &
             < min_wp2_sfc_val(i) ) then
        !if ( wp2(i,1) + tau_zm(i,1) * wp2_splat(i,1) < min_wp2_sfc_val(i) ) then 
                           ! splatting correction drives wp2 to overly small values
          wp2_splat_sfc_correction(i) = -wp2(i,1) + min_wp2_sfc_val(i)
          wp2(i,1) = min_wp2_sfc_val(i)
        else
         wp2_splat_sfc_correction(i) = tau_zm(i,1) * lhs_splat_wp2(i,1) * wp2(i,1)
         !wp2_splat_sfc_correction(i) = tau_zm(i,1) * wp2_splat(i,1)
         wp2(i,1) = wp2(i,1) + wp2_splat_sfc_correction(i)
        end if

        up2(i,1) = up2(i,1) - 0.5_core_rknd * wp2_splat_sfc_correction(i)
        vp2(i,1) = vp2(i,1) - 0.5_core_rknd * wp2_splat_sfc_correction(i)
      end do
      !$acc end parallel loop

      ! Passive scalars
      if ( sclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol
            ! Vince Larson changed coeffs to make correlations between [-1,1].
            ! 31 Jan 2008
            ! sclrprtp(i) &
            ! = a * (wprtp_sfc / uf) * (wpsclrp(i) / uf)
            ! sclrpthlp(i) &
            ! = a * (wpthlp / uf) * (wpsclrp(i) / uf)
            ! sclrp2(i) &
            ! = sclr_var_coef * a * ( wpsclrp(i) / uf )**2
            ! Notes:  1) With "a" having a value of 1.8 and "sclr_var_coef"
            !            having a value of 0.4, the surface correlation of
            !            w & sclr has a value of about 0.878.
            !         2) With "sclr_var_coef" having a value of 0.4, the
            !            surface correlations of both rt & sclr and
            !            thl & sclr are 0.5.
            ! Brian Griffin; February 2, 2008.

            ! We use the following if..then's to make sclr_rt and sclr_thl
            ! close to the actual thlp2/rtp2 at the surface.
            ! -dschanen 25 Sep 08
            if ( sclr == iisclr_rt ) then
              ! If we are trying to emulate rt with the scalar, then we
              ! use the variance coefficient from above
              sclrprtp(i,1,sclr) = 0.4_core_rknd * a_const &
                                * ( wprtp_sfc(i) / uf(i) ) * ( wpsclrp_sfc(i,sclr) / uf(i) )
            else
              sclrprtp(i,1,sclr) = 0.2_core_rknd * a_const &
                                * ( wprtp_sfc(i) / uf(i) ) * ( wpsclrp_sfc(i,sclr) / uf(i) )
            endif

            if ( sclr == iisclr_thl ) then
              ! As above, but for thetal
              sclrpthlp(i,1,sclr) = 0.4_core_rknd * a_const &
                                 * ( wpthlp(i,1) / uf(i) ) &
                                 * ( wpsclrp_sfc(i,sclr) / uf(i) )
            else
              sclrpthlp(i,1,sclr) = 0.2_core_rknd * a_const &
                                 * ( wpthlp(i,1) / uf(i) ) &
                                 * ( wpsclrp_sfc(i,sclr) / uf(i) )
            endif

            sclrp2(i,1,sclr) = sclr_var_coef * a_const &
                               * ( wpsclrp_sfc(i,sclr) / uf(i) )**2

            ! End Vince Larson's change.

          end do
        enddo ! 1,...sclr_dim
        !$acc end parallel loop

      endif ! sclr_dim > 0

    endif ! l_andre_1978

    ! Clip wp2 at wp2_max, same as in advance_wp2_wp3
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wp2(i,1) = min( wp2(i,1), wp2_max )
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      if ( abs(gr%zm(i,1)-sfc_elevation(i)) > abs(gr%zm(i,1)+sfc_elevation(i))*eps/2 ) then

        ! Variances for cases where the lowest level is not at the surface.
        ! Eliminate surface effects on lowest level variances.
        wp2(i,1)     = w_tol_sqd
        up2(i,1)     = w_tol_sqd
        vp2(i,1)     = w_tol_sqd
        thlp2(i,1)   = thl_tol**2
        rtp2(i,1)    = rt_tol**2
        rtpthlp(i,1) = 0.0_core_rknd

      end if ! gr%zm(1,1) == sfc_elevation
    end do
    !$acc end parallel loop

    if ( sclr_dim > 0 ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do sclr = 1, sclr_dim
        do i = 1, ngrdcol
          if ( abs(gr%zm(i,1)-sfc_elevation(i)) > abs(gr%zm(i,1)+sfc_elevation(i))*eps/2 ) then
            
            ! Variances for cases where the lowest level is not at the surface.
            ! Eliminate surface effects on lowest level variances.
            sclrp2(i,1,sclr)    = 0.0_core_rknd
            sclrprtp(i,1,sclr)  = 0.0_core_rknd
            sclrpthlp(i,1,sclr) = 0.0_core_rknd
          end if
        end do
      end do
    end if

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( wp2, up2, vp2, thlp2, rtp2, rtpthlp )

      do i = 1, ngrdcol
        call stat_end_update_pt( stats_metadata%ithlp2_sf, 1,    & ! intent(in)
                                 thlp2(i,1) / dt, & ! intent(in)
                                 stats_zm(i) )      ! intent(inout)

        call stat_end_update_pt( stats_metadata%irtp2_sf, 1,     & ! intent(in)
                                 rtp2(i,1) / dt,  & ! intent(in)
                                 stats_zm(i) )      ! intent(inout)

        call stat_end_update_pt( stats_metadata%irtpthlp_sf, 1,    & ! intent(in)
                                 rtpthlp(i,1) / dt, & ! intent(in)
                                 stats_zm(i) )        ! intent(inout)

        call stat_end_update_pt( stats_metadata%iup2_sf, 1,    & ! intent(in)
                                 up2(i,1) / dt, & ! intent(in)
                                 stats_zm(i) )    ! intent(inout)

        call stat_end_update_pt( stats_metadata%ivp2_sf, 1,    & ! intent(in)
                                 vp2(i,1) / dt, & ! intent(in)
                                 stats_zm(i) )    ! intent(inout)

        call stat_end_update_pt( stats_metadata%iwp2_sf, 1,    & ! intent(in)
                                 wp2(i,1) / dt, & ! intent(in)
                                 stats_zm(i) )    ! intent(inout)
      end do
    end if

    if ( clubb_at_least_debug_level( 2 ) ) then 

      !$acc update host( wp2, up2, vp2, thlp2, rtp2, rtpthlp, &
      !$acc              sclrp2, sclrprtp, sclrpthlp, &
      !$acc              upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc )

      do i = 1, ngrdcol
        call sfc_varnce_check( wp2(i,1), up2(i,1), vp2(i,1),                    & ! intent(in)
                               thlp2(i,1), rtp2(i,1), rtpthlp(i,1),             & ! intent(in)
                               sclrp2(i,1,:), sclrprtp(i,1,:), sclrpthlp(i,1,:) ) ! intent(in)
      end do

      if ( err_code == clubb_fatal_error ) then

        write(fstderr,*) "Error in calc_sfc_varnce"
        write(fstderr,*) "Intent(in)"

        write(fstderr,*) "upwp_sfc = ", upwp_sfc
        write(fstderr,*) "vpwp_sfc = ", vpwp_sfc
        write(fstderr,*) "wpthlp = ", wpthlp
        write(fstderr,*) "wprtp_sfc = ", wprtp_sfc

        if ( sclr_dim > 0 ) then
           write(fstderr,*) "wpsclrp = ", wpsclrp_sfc
        endif

        write(fstderr,*) "Intent(out)"

        write(fstderr,*) "wp2 = ", wp2
        write(fstderr,*) "up2 = ", up2
        write(fstderr,*) "vp2 = ", vp2
        write(fstderr,*) "thlp2 = ", thlp2
        write(fstderr,*) "rtp2 = ", rtp2
        write(fstderr,*) "rtpthlp = ", rtpthlp

        if ( sclr_dim > 0 ) then
           write(fstderr,*) "sclrp2 = ", sclrp2
           write(fstderr,*) "sclrprtp = ", sclrprtp
           write(fstderr,*) "sclrpthlp = ", sclrpthlp
        endif

      endif ! err_code == clubb_fatal_error

    endif ! clubb_at_least_debug_level ( 2 )! Update surface stats

    !$acc exit data delete( uf, depth_pos_wpthlp, min_wp2_sfc_val, &
    !$acc                   um_sfc_sqd, vm_sfc_sqd, usp2_sfc, vsp2_sfc, &
    !$acc                   ustar, zeta, wp2_splat_sfc_correction )

    return

  end subroutine calc_sfc_varnce

!===============================================================================

end module sfc_varnce_module
