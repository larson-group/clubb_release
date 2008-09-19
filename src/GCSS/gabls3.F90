#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: gabls3.F90 2839 2008-09-04 14:50:14Z faschinj $
module gabls3

!       Description:
!       Contains subroutines for the GABLS3.
!----------------------------------------------------------------------

  implicit none

  public :: gabls3_tndcy, gabls3_sfclyr

  private :: time_select ! Defualt Scope

  private

  contains

!----------------------------------------------------------------------
  subroutine gabls3_tndcy( time, rtm, exner, p_in_Pa, thvm, &
                           wm_zt, wm_zm, thlm_forcing, rtm_forcing )
!       Description:
!       Subroutine to set thetal and total water tendencies for GABLS3 case

!       References:
!       None
!----------------------------------------------------------------------

    use grid_class, only: gr,zt2zm ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use error_code, only: clubb_debug ! Procedure(s)

    use constants, only: Cp, Lv, grav, Rd ! Procedure(s)

    use interpolation, only:factor_interp, lin_int ! Procedure(s)

    implicit none

    ! Prescribed Geostrophic wind
    real, parameter, dimension(6) :: &
        ug_sfc = (/-7.8,-7.8, -6.5, -5.0, -5.0, -6.5/), &
        vg_sfc = (/0.0, 0.0, 4.5, 4.5, 4.5, 2.5/)
    !t_geo = [43200.0 64800.0 82800.0 97200.0 108000.0];


    ! Prescribed Vertical Dynamic tendency at specified time %
    real, parameter, dimension(4) :: omega_t = (/0.12, 0.12, 0.00, 0.00/),&
    ! Times specified for omega
    omega_time = (/43200.0, 61200.0, 68400.0, 129600.0/);

    ! Prescribed Temperature dynamic tendency
    real, parameter, dimension(6) :: &
      T_adv_t = (/-2.5E-5, -2.5E-5, 7.5E-5, 7.5E-5, 0.0E+0, 0.0E+0/), &
      t_adv_time = (/43200.0, 90000.0, 90000.0, 108000.0, 108000.0, 129600.0/);

    ! Prescribed Specific Humidity dynamic tendency
    real, parameter, dimension(10) :: &
      q_adv_t = (/0.0E+0, 0.0E+0, 8.0E-8, 8.0E-8, 0.0E+0, 0.0E+0, &
                 -8.0E-8, -8.0E-8, 0.0E+0, 0.0E+0/), &
      q_adv_time = (/43200, 75600, 75600, 86400, 86400, & 
                     93600, 93600, 97200, 97200, 129600/);

    ! Prescribed winds
    real, parameter, dimension(8) :: & 
      u_adv_t = (/0.0E+0, 0.0E+0, -1.5E-4, -1.5E-4, 5.0E-4, &
                  5.0E-4, 0.0E+0, 0.0E+0/),&
      v_adv_t = (/0.0E+0, 0.0E+0, 1.0E-4, 1.0E-4, 0.0E+0, & 
                  0.0E+0, 0.0E+0, 0.0E+0/);
    !t_adv_wind = [43200 64800 64800 82800 82800 97200 97200 129600];

    ! Input Variables

    real(kind=time_precision), intent(in) :: time ! Model time [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      rtm,  &     ! Total water mixing ratio                        [kg/kg]
      exner,&     ! Exner Function function = (p/p0 ** kappa)       [-]
      p_in_Pa, &  ! Pressure (Pa) on thermodynamic points           [Pa]
      thvm


    ! Output Variables

    real, intent(out), dimension(gr%nnzp) ::  & 
      wm_zt,         &
      wm_zm,         &
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing    ! Total water mixing ratio tendency            [kg/kg/s]

    real, dimension(gr%nnzp) :: velocity_omega
    real :: time_frac, T_t_interp, q_t_interp, omega_t_interp
    integer :: i1, i2,i

    call time_select( time, T_adv_time, 6, i1, i2 )

    time_frac = real((time - T_adv_time(i1)) /  &          ! at the first time a=0;
            (T_adv_time(i2) - T_adv_time(i1)))             ! at the second time a=1.

    T_t_interp = factor_interp(time_frac, T_adv_t(i2), T_adv_t(i1))

    call time_select( time, q_adv_time, 10, i1, i2)
    time_frac = real((time - q_adv_time(i1)) /  &          ! at the first time a=0;
            (q_adv_time(i2) - q_adv_time(i1)))             ! at the second time a=1.

    q_t_interp = factor_interp( time_frac, q_adv_t(i2), q_adv_t(i1))

    call time_select( time, omega_time, 4, i1, i2)
    time_frac = real((time - omega_time(i1)) /  &          ! at the first time a=0;
            (omega_time(i2) - omega_time(i1)))             ! at the second time a=1.

    omega_t_interp = factor_interp( time_frac, omega_t(i2), omega_t(i1))



    do i = 1, gr%nnzp,1
      select case (int (gr%zt(i)))
      case (1:199)
        thlm_forcing(i) = lin_int(gr%zt(i), 200., 0., T_t_interp, 0. )
        rtm_forcing(i)  =  lin_int(gr%zt(i), 200., 0., q_t_interp, 0.)
        velocity_omega(i)  =  lin_int(gr%zt(i), 200., 0., omega_t_interp, 0.)
      case (200:999)
        thlm_forcing(i) = T_t_interp ! Convert!
        rtm_forcing(i) = q_t_interp
        velocity_omega(i) = omega_t_interp
      case (1000:1499)
        thlm_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0., T_t_interp )
        rtm_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0. ,q_t_interp )
        velocity_omega(i) = lin_int( gr%zt(i), 1500., 1000., 0. ,omega_t_interp )
      case default
        thlm_forcing(i) = 0
        rtm_forcing(i) = 0
        velocity_omega(i) = 0
      end select
      !print *, "zt (",i,") =", gr%zt(i)
      !print *, "T_t_interp (",i,") =", T_t_interp
      !print *, "q_t_interp (",i,") =", q_t_interp
      !print *, "Thlm_forcing (",i,") =", thlm_forcing(i)
      !print *, "rtm_forcing (",i,") =", rtm_forcing(i)
    end do
!    rtm_forcing = rtm_forcing/(1. - rtm_forcing)
    rtm_forcing = rtm_forcing * ( 1. + rtm )**2
!    thlm_forcing = ( thlm_forcing - Lv * rcm/Cp ) / exner
    thlm_forcing = thlm_forcing / exner
   ! Compute vertical motion
    do i=2,gr%nnzp
      wm_zt(i) = -velocity_omega(i) * Rd * thvm(i) / p_in_Pa(i) / grav

    end do
 
   ! Boundary condition
   wm_zt(1) = 0.0        ! Below surface
 
   ! Interpolation
   wm_zm = zt2zm( wm_zt )
 
   ! Boundary condition
   wm_zm(1) = 0.0        ! At surface
   wm_zm(gr%nnzp) = 0.0  ! Model top



!-----------------------------------------------------------------------


  end subroutine gabls3_tndcy
!----------------------------------------------------------------------
  subroutine gabls3_sfclyr( time, z, rho0, & 
                            thlm_sfc, um_sfc, vm_sfc,  & 
                            upwp_sfc, vpwp_sfc, & 
                            wpthlp_sfc, wprtp_sfc, ustar, & 
                            wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

    use constants, only: Cp, Lv, grav ! Variable(s)

    use parameters_tunable, only: sclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_mod, only: diag_ustar ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl

    use interpolation, only: factor_interp

    implicit none

    intrinsic :: max, sqrt, present

    real, parameter ::  & 
      ubmin = 0.25, & ! Minimum value for ubar 
      z0    = 0.035   ! ARM Cu mom. roughness height


! Input Variables
    real(time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real, intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho0,      & ! Density at zm=1     [kg/m^3] 
      um_sfc,    & ! um at (2)           [m/s]
      vm_sfc,    & ! vm at (2)           [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

! Output variables
    real, intent(out) ::  & 
      upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
      vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

! Output variables (optional)
    real, intent(out), optional, dimension(sclr_dim) ::  & 
      wpsclrp_sfc,   & ! Passive scalar surface flux      [units m/s] 
      wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]


    return
  end subroutine gabls3_sfclyr
!----------------------------------------------------------------------
  subroutine time_select( time, time_array, nvar, left_time, right_time )
    use stats_precision, only: time_precision ! Variable(s)

    use interpolation, only: binary_search
    implicit none
    real(kind=time_precision), intent(in) :: time
    real, dimension(nvar), intent(in) :: time_array
    integer, intent(in) :: nvar
    integer, intent(out) :: right_time, left_time
    integer :: k

    ! Interpolate in time for Temperature tendency
    if( time <= time_array(1)) then
      left_time = 1
      right_time = 2
    else if ( time >= time_array(nvar) ) then
      left_time = nvar
      right_time = nvar - 1
    else
      do k=1,nvar-1
        if ((time > time_array(k)) .AND. &
         (time <= time_array(k+1))) then
          left_time = k
          right_time = k+1
        end if
      end do
    endif
    !print *, "Left time = ", left_time
    !print *, "Right time = ", right_time
  end subroutine time_select
end module gabls3

