!----------------------------------------------------------------------
! $Id$
module gabls3

  !       Description:
  !       Contains subroutines for the GABLS3.
  !----------------------------------------------------------------------

  implicit none

  public :: gabls3_tndcy, gabls3_sfclyr

  private

  contains

  !----------------------------------------------------------------------
  subroutine gabls3_tndcy( time, rtm, exner, rho, &
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

    use interpolation, only:factor_interp, lin_int,binary_search ! Procedure(s)

    use time_dependant_input, only: l_t_dependant, &
                                    time_select, &
                                    thlm_f_given, &
                                    rtm_f_given,&
                                    wm_given,&
                                    um_f_given,&
                                    vm_f_given,&
                                    ug_given,&
                                    vg_given,&
                                    time_f_given
    implicit none

    ! Input Variables

    real(kind=time_precision), intent(in) :: time ! Model time [s]

    real, intent(in), dimension(gr%nnzp) ::  & 
      rtm,  &     ! Total water mixing ratio                        [kg/kg]
      exner,&     ! Exner Function function = (p/p0 ** kappa)       [-]
      rho         ! Air density at surface                          [kg/m^3]


    ! Output Variables

    real, intent(out), dimension(gr%nnzp) ::  & 
      wm_zt,        &  ! w on the thermodynamic levels                [m/s]
      wm_zm,        &  ! w on the momentum levels                     [m/s]
      thlm_forcing, &  ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing,  &  ! Total water mixing ratio tendency            [kg/kg/s]
      um_forcing,   &  ! u wind forcing                               [m/s/s]
      vm_forcing,   &  ! v wind forcing                               [m/s/s]
      ug,           &  ! u geostrophic wind                           [m/s]
      vg               ! v geostrophic wind                           [m/s]

    real, dimension(gr%nnzp) :: velocity_omega, T_in_K_forcing, sp_humidity_forcing

    real :: time_frac

    integer :: i1, i2
 
    if( l_t_dependant ) then
      call time_select( time, size(time_f_given), time_f_given, i1, i2 )

      time_frac = real((time - time_f_given(i1)) /  &          ! at the first time a=0;
              (time_f_given(i2) - time_f_given(i1)))             ! at the second time a=1.

      T_in_K_forcing = factor_interp(time_frac, thlm_f_given(:,i2), thlm_f_given(:,i1))


      sp_humidity_forcing = factor_interp( time_frac, rtm_f_given(:,i2), rtm_f_given(:,i1))


      velocity_omega  = factor_interp( time_frac, wm_given(:,i2), wm_given(:,i1))

      um_forcing  = factor_interp( time_frac, um_f_given(:,i2), um_f_given(:,i1) )
      vm_forcing = factor_interp( time_frac, vm_f_given(:,i2), vm_f_given(:,i1) )

      ug = factor_interp( time_frac, ug_given(:,i2), ug_given(:,i1) )
      vg = factor_interp( time_frac, vg_given(:,i2), vg_given(:,i1) )

      rtm_forcing = sp_humidity_forcing * ( 1. + rtm )**2

      thlm_forcing = T_in_K_forcing / exner

      wm_zt = -velocity_omega /( rho * grav );

      ! Boundary condition
      wm_zt(1) = 0.0        ! Below surface

      ! Interpolation
      wm_zm = zt2zm( wm_zt )

      ! Boundary condition
      wm_zm(1) = 0.0        ! At surface
      wm_zm(gr%nnzp) = 0.0  ! Model top

    end if
  end subroutine gabls3_tndcy

  !-----------------------------------------------------------------------
  subroutine gabls3_sfclyr( um_sfc, vm_sfc, veg_t_in_K, &
                            thlm_sfc, rtm_sfc, lowest_level, exner_sfc, & 
                            upwp_sfc, vpwp_sfc, &
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ATEX specifications

    !       References:

    !----------------------------------------------------------------------

    use constants, only: kappa, grav, Rd, Cp, p0, Lv ! Variable(s)

    use diag_ustar_mod, only: diag_ustar ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use surface_flux, only: compute_momentum_flux, &
                            compute_wpthlp_sfc, &
                            compute_wprtp_sfc

    !use surface, only: prognose_soil_T_in_K ! Procedure(s)

    implicit none

    ! Constants

    real, parameter ::  & 
      ubmin = 0.25, & 
     ! ustar = 0.3,
     ! C_10  = 0.0013, & !ATEX value
     ! C_10  = 0.013, & ! Fudged value
     ! C_10  = 0.0049, & ! Fudged value
     ! C_10  = 0.0039, & ! Fudged value
      C_10 = 0.00195, &
     ! C_10 = 0.001, &
     ! C_10 = 0.003, &
      z0 = 0.15

!    real, parameter, dimension(25) :: sst_given = (/300., 300.8, 300.9, 301.,300.9, &
!                                        300.5, 300., 298.5, 297., 296., 295.,&
!                                        294., 293.5, 292.5, 291.5, 291.,&
!                                        290.5, 292.5, 294.5, 296.5, 298.,&
!                                        298.5, 300.5, 301.5, 301./)
!    real, parameter, dimension(25) :: sst_time = (/43200., 46800., 50400., 54000., 57600., &
!                                        61200., 64800., 68400., 72000., 75600.,&
!                                        79200., 82800., 86400., 90000., 93600., &
!                                        97200., 100800., 104400., 108000., 111600.,&
!                                        115200., 118800., 122400., 126000., 129600./)


    ! Input variables


    real, intent(in) ::  & 
      um_sfc,       & ! um at zt(2)            [m/s]
      vm_sfc,       & ! vm at zt(2)            [m/s]
      thlm_sfc,     & ! Theta_l at zt(2)       [K]
      rtm_sfc,      & ! rt at zt(2)            [kg/kg]
      veg_T_in_K,   & ! Vegetation temperature [K]
      lowest_level, & ! gr%zt(2)               [m]
      exner_sfc       ! Exner function         [-]

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,    & ! u'w' at surface           [m^2/s^2]
      vpwp_sfc,    & ! v'w' at surface           [m^2/s^2]
      ustar          ! surface friction velocity [m/s]

    real, intent(inout):: &
      wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
      wprtp_sfc      ! w'rt' surface flux        [(m kg)/(kg s)]

    ! Local Variables
    real :: ubar, veg_theta_in_K, bflx

    ! Compute heat and moisture fluxes
    ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

    veg_theta_in_K = veg_T_in_K / exner_sfc

    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, veg_T_in_K, exner_sfc)

    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, 9.9e-3 ) * 10
  
    ! Compute momentum fluxes
    bflx = wpthlp_sfc * grav / veg_theta_in_K

    ustar = diag_ustar( lowest_level, bflx, ubar, z0)

    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    return
  end subroutine gabls3_sfclyr

end module gabls3
