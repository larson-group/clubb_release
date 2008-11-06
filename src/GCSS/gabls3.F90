!----------------------------------------------------------------------
! $Id: gabls3.F90 2839 2008-09-04 14:50:14Z faschinj $
module gabls3

  !       Description:
  !       Contains subroutines for the GABLS3.
  !----------------------------------------------------------------------

  implicit none

  public :: gabls3_tndcy, gabls3_sfclyr

  private :: time_select ! Defualt Scope

  real, private, dimension(2) :: veg_T_in_K, sfc_soil_T_in_K, deep_soil_T_in_K

  private

  contains

  !----------------------------------------------------------------------
  subroutine gabls3_tndcy( time, rtm, exner, p_in_Pa, thvm, &
                           wm_zt, wm_zm, thlm_forcing, rtm_forcing,&
                           um_forcing, vm_forcing, ug, vg )
    !       Description:
    !       Subroutine to set thetal and total water tendencies for GABLS3 case

    !       References:
    !       None
    !----------------------------------------------------------------------

    use grid_class, only: gr, zt2zm ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use constants, only: Cp, Lv, grav, Rd ! Procedure(s)

    use interpolation, only:factor_interp, lin_int ! Procedure(s)

    implicit none

    ! Prescribed Geostrophic wind
    real, parameter, dimension(6) :: &
        ug_sfc = (/-7.8,-7.8, -6.5, -5.0, -5.0, -6.5/), &
        vg_sfc = (/0.0, 0.0, 4.5, 4.5, 4.5, 2.5/), &
        geo_time = (/43200.0, 64800.0, 82800.0, 97200.0, 108000.0, 129600.0/)


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
                     93600, 93600, 104400, 104400, 129600/);

    ! Prescribed winds
    real, parameter, dimension(8) :: & 
      u_adv_t = (/0.0E+0, 0.0E+0, -1.5E-4, -1.5E-4, 5.0E-4, &
                  5.0E-4, 0.0E+0, 0.0E+0/),&
      v_adv_t = (/0.0E+0, 0.0E+0, 1.0E-4, 1.0E-4, 0.0E+0, & 
                  0.0E+0, 0.0E+0, 0.0E+0/),&
      uv_adv_time = (/43200, 64800, 64800, 82800, 82800, 97200, 97200, 129600/);

    ! Input Variables

    real(kind=time_precision), intent(in) :: time ! Model time [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      rtm,  &     ! Total water mixing ratio                        [kg/kg]
      exner,&     ! Exner Function function = (p/p0 ** kappa)       [-]
      p_in_Pa, &  ! Pressure (Pa) on thermodynamic points           [Pa]
      thvm        ! Virtual Potential Temperature


    ! Output Variables

    real, intent(out), dimension(gr%nnzp) ::  & 
      wm_zt,        &
      wm_zm,        &
      thlm_forcing, &  ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing,  &  ! Total water mixing ratio tendency            [kg/kg/s]
      um_forcing,   &
      vm_forcing,   &
      ug,           &
      vg

    real, dimension(gr%nnzp) :: velocity_omega, T_in_K_forcing, sp_humidity_forcing

    real :: time_frac, T_t_interp, q_t_interp, omega_t_interp,&
            u_t_interp, v_t_interp, ug_interp,vg_interp

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

    call time_select( time, uv_adv_time, 8, i1, i2)
    time_frac = real((time - uv_adv_time(i1)) /  &          ! at the first time a=0;
            (uv_adv_time(i2) - uv_adv_time(i1)))             ! at the second time a=1.

    u_t_interp = factor_interp( time_frac, u_adv_t(i2), u_adv_t(i1) )
    v_t_interp = factor_interp( time_frac, v_adv_t(i2), v_adv_t(i1) )

    call time_select( time, geo_time, 6, i1, i2)

    time_frac = real((time - geo_time(i1)) /  &          ! at the first time a=0;
            (geo_time(i2) - geo_time(i1)))             ! at the second time a=1.

    ug_interp = factor_interp( time_frac, ug_sfc(i2), ug_sfc(i1) )
    vg_interp = factor_interp( time_frac, vg_sfc(i2), vg_sfc(i1) )

    do i = 1, gr%nnzp,1
      select case (int (gr%zt(i)))
      case (1:199)
        T_in_K_forcing(i) = lin_int(gr%zt(i), 200., 0., T_t_interp, 0. )
        sp_humidity_forcing(i)  =  lin_int(gr%zt(i), 200., 0., q_t_interp, 0.)
        velocity_omega(i)  =  lin_int(gr%zt(i), 200., 0., omega_t_interp, 0.)
        um_forcing(i)  =  lin_int(gr%zt(i), 200., 0., u_t_interp, 0.)
        vm_forcing(i)  =  lin_int(gr%zt(i), 200., 0., v_t_interp, 0.)
      case (200:999)
        T_in_K_forcing(i) = T_t_interp ! Convert!
        sp_humidity_forcing(i) = q_t_interp
        velocity_omega(i) = omega_t_interp
        um_forcing(i) = u_t_interp
        vm_forcing(i) = v_t_interp
      case (1000:1499)
        T_in_K_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0., T_t_interp )
        sp_humidity_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0. ,q_t_interp )
        velocity_omega(i) = lin_int( gr%zt(i), 1500., 1000., 0. ,omega_t_interp )
        um_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0. , u_t_interp )
        vm_forcing(i) = lin_int( gr%zt(i), 1500., 1000., 0. , v_t_interp )
      case default
        T_in_K_forcing(i) = 0
        sp_humidity_forcing(i) = 0
        velocity_omega(i) = 0
        um_forcing(i) = 0
        vm_forcing(i) = 0
      end select

      if( gr%zt(i) < 2000 ) then
        ! Winterpolate
        ug(i) = lin_int( gr%zt(i), 2000., 0., -2., ug_interp )
        vg(i) = lin_int( gr%zt(i), 2000., 0., 2., vg_interp )
      else
        ug(i) = -2.0
        vg(i) = 2.0
      endif
      !print *, "zt (",i,") =", gr%zt(i)
      !print *, "T_t_interp (",i,") =", T_t_interp
      !print *, "q_t_interp (",i,") =", q_t_interp
      !print *, "Thlm_forcing (",i,") =", thlm_forcing(i)
      !print *, "rtm_forcing (",i,") =", rtm_forcing(i)
    end do
!    rtm_forcing = rtm_forcing/(1. - rtm_forcing)
    rtm_forcing = sp_humidity_forcing * ( 1. + rtm )**2
!    thlm_forcing = ( thlm_forcing - Lv * rcm/Cp ) / exner
    thlm_forcing = T_in_K_forcing / exner
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
    i = 1

  end subroutine gabls3_tndcy

  !-----------------------------------------------------------------------
  subroutine gabls3_sfclyr( time_start, time_current, dt, rho_sfc, um_sfc, vm_sfc, &
                            thlm_sfc, rtm_sfc, lowest_level, psfc,& 
                            Frad_SW_up_sfc, Frad_SW_down_sfc, &
                            Frad_LW_up_sfc, Frad_LW_down_sfc, &
                            upwp_sfc, vpwp_sfc,  & 
                            wpthlp_sfc, wprtp_sfc, ustar, & 
                            wpsclrp_sfc, wpedsclrp_sfc )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ATEX specifications

    !       References:

    !----------------------------------------------------------------------

    use constants, only: kappa, grav, Rd, Cp, p0,lv ! Variable(s)

    use parameters_model, only: sclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl ! Variable(s)

    use diag_ustar_mod, only: diag_ustar ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use stats_variables, only: l_stats_samp, sfc, iveg_t_sfc, it_sfc, ideep_T_sfc ! Variables

    use stats_type, only: stat_update_var_pt ! Procedure(s)

    use interpolation, only: lin_int ! Procedure(s)

    use surface, only: prognose_soil_T_in_K ! Procedure(s)

    implicit none

    ! Constants

    real, parameter ::  & 
      ubmin = 0.25, & 
     ! ustar = 0.3,
!     C_10  = 0.0013, & !ATEX value
!     C_10  = 0.013, & ! Fudged value
!      C_10  = 0.0039, & ! Fudged value
      C_10 = 0.00195, &
      z0 = 0.15

    real, parameter, dimension(25) :: sst_given = (/300., 300.8, 300.9, 301.,300.9, &
                                        300.5, 300., 298.5, 297., 296., 295.,&
                                        294., 293.5, 292.5, 291.5, 291.,&
                                        290.5, 292.5, 294.5, 296.5, 298.,&
                                        298.5, 300.5, 301.5, 301./)
    real, parameter, dimension(25) :: sst_time = (/43200., 46800., 50400., 54000., 57600., &
                                        61200., 64800., 68400., 72000., 75600.,&
                                        79200., 82800., 86400., 90000., 93600., &
                                        97200., 100800., 104400., 108000., 111600.,&
                                        115200., 118800., 122400., 126000., 129600./)


    ! Input variables

    real(kind=time_precision), intent(in) :: time_start, time_current,dt

    real, intent(in) ::  & 
      um_sfc,     & ! um at zt(2)           [m/s]
      vm_sfc,     & ! vm at zt(2)           [m/s]
      thlm_sfc,   & ! Theta_l at zt(2)      [K]
      rtm_sfc,    & ! rt at zt(2)           [kg/kg]
      rho_sfc,    &
      lowest_level, &
      psfc,&
      Frad_SW_up_sfc, &
      Frad_SW_down_sfc,&
      Frad_LW_up_sfc, &
      Frad_LW_down_sfc

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,    & ! u'w' at surface           [m^2/s^2]
      vpwp_sfc,    & ! v'w' at surface           [m^2/s^2]
      ustar          ! surface friction velocity [m/s]

    real, intent(inout):: &
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc      ! w'rt' surface flux        [(m kg)/(kg s)]

    ! Output variables (optional)

    real, dimension(sclr_dim), intent(out) ::  & 
      wpsclrp_sfc,    & ! Passive scalar surface flux      [units m/s]
      wpedsclrp_sfc     ! Passive eddy-scalar surface flux [units m/s]

    ! Local Variables
    real :: ubar, veg_theta_in_K, bflx, wpthep

    integer :: i1, i2
    ! Compute heat and moisture fluxes
    ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

    ! Set SST by time in lieu of a surface scheme
    !call time_select( time_current, sst_time, 25, i1, i2)
    !sfc_soil_T_in_K(1) = lin_int( real(time_current), sst_time(i2), &
    !                   sst_time(i1), sst_given(i2), sst_given(i1) )
    !----- Experimental Code -------------
    ! Turbulent Flux of equivalent potential temperature
    wpthep = wpthlp_sfc + (Lv/Cp) * ((p0/psfc)**kappa) * wprtp_sfc

    call prognose_soil_T_in_K( time_start, time_current, real( dt), 1, 2, rho_sfc, &
                               Frad_SW_down_sfc-Frad_SW_up_sfc, Frad_SW_down_sfc,&
                               Frad_LW_down_sfc, Frad_LW_up_sfc, wpthep, &
                               veg_T_in_K, sfc_soil_T_in_K, deep_soil_T_in_K )

    if(l_stats_samp) then
      call stat_update_var_pt( iveg_t_sfc, 1, veg_T_in_K(1), sfc )
      call stat_update_var_pt( it_sfc, 1, sfc_soil_T_in_K(1), sfc )
      call stat_update_var_pt( ideep_t_sfc, 1, deep_soil_T_in_K(1), sfc )
    end if

    sfc_soil_T_in_K(1) = sfc_soil_T_in_K(2)
    veg_T_in_K(1) = veg_T_in_K(2)
    deep_soil_T_in_K(1) = deep_soil_T_in_K(2)
    !-------------------------------------

    !sstheta = sfc_soil_T_in_K(1) * (( p0 / psfc )**(Rd/Cp))
    veg_theta_in_K = veg_T_in_K(1) * (( p0 / psfc )**(Rd/Cp))

    wpthlp_sfc = -C_10 * ubar * ( thlm_sfc - veg_theta_in_K )
    wprtp_sfc  = -C_10 * ubar * ( rtm_sfc - 7.5e-3 )

    ! Let passive scalars be equal to rt and theta_l for now
    if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

    if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

    ! Compute momentum fluxes
    bflx = wpthlp_sfc * grav / veg_theta_in_K
    ustar = diag_ustar( lowest_level, bflx, ubar, z0)

    upwp_sfc = -um_sfc * ustar**2 / ubar
    vpwp_sfc = -vm_sfc * ustar**2 / ubar

    return
  end subroutine gabls3_sfclyr

  !----------------------------------------------------------------------
  subroutine time_select( time, time_array, nvar, left_time, right_time )
    !
    !   Description: This subroutine determines which indexes of the given
    !                time_array should be used when interpolating a value 
    !                at the specified time.
    !
    !
    !----------------------------------------------------------------------

    use stats_precision, only: time_precision ! Variable(s)

    use interpolation, only: binary_search    ! Procedure(s)
  
    implicit none
  
    real(kind=time_precision), intent(in) :: time       ! Target time              [s]
  
    real, dimension(nvar), intent(in) :: time_array     ! Array of times           [s]
  
    integer, intent(in) :: nvar                         ! Number of array elements []
  
    integer, intent(out) :: &
      right_time, &                                     ! Index of a time later
      !                                                   than the target time []
      left_time                                         ! Index of time before 
      !                                                   the target time       []
  
    integer :: k

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

