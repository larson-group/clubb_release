!----------------------------------------------------------------------
! $Id$
module arm_97

  !       Description:
  !       Contains subroutines for the July 26-30 1997 ARM IOP A case.
  !       
  !       References:
  !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
  !----------------------------------------------------------------------

  implicit none

  public :: arm_97_sfclyr

  private ! Defualt Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_97_sfclyr( ngrdcol, time, z, rho_sfc, & 
                            thlm_sfc, ubar,  & 
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !
    !       References:
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
    !----------------------------------------------------------------------

    use constants_clubb, only: grav ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use sfc_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

    use time_dependent_input, only: &
      time_select, &
      time_sfc_given, &
      latent_ht_given, &
      sens_ht_given, &
      l_t_dependent

    implicit none

    intrinsic :: max, sqrt, present

    real( kind = core_rknd ), parameter ::  & 
      z0    = 0.035_core_rknd   ! ARM Cu mom. roughness height

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real(time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho_sfc,      & ! Density at zm=1     [kg/m^3] 
      ubar,      & ! mean sfc wind speed [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables

    real( kind = core_rknd ) :: bflx, heat_flx, moisture_flx, time_frac

    integer :: before_time, after_time, i

    !----------------------------------------------------------------------

    if ( l_t_dependent ) then

      time_frac = -1.0_core_rknd ! Default initialization

      call time_select( time, size(time_sfc_given), time_sfc_given, &
                                   before_time, after_time, time_frac )

      heat_flx = linear_interp_factor( time_frac, sens_ht_given(after_time), &
                                       sens_ht_given(before_time) )
      moisture_flx = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                           latent_ht_given(before_time) )

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol

        ! Convert W/m^2 into w'thl' w'rt' units
        wpthlp_sfc(i) = convert_sens_ht_to_km_s( heat_flx, rho_sfc(i) )     ! (K m/s)
        wprtp_sfc(i)  = convert_latent_ht_to_m_s( moisture_flx, rho_sfc(i) ) ! (kg m/ kg s)

        ! Compute momentum fluxes using ARM Cu formulae

        bflx = grav / thlm_sfc(i) * wpthlp_sfc(i)

        ! Compute ustar
        ustar(i) = diag_ustar( z(i), bflx, ubar(i), z0 )

      end do

    end if

    return

  end subroutine arm_97_sfclyr
  !----------------------------------------------------------------------
end module arm_97
