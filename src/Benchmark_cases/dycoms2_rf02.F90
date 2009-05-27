!----------------------------------------------------------------------
! $Id$
module dycoms2_rf02

  !       Description:
  !       Contains subroutines for the DYCOMS II RF02 case.
  !----------------------------------------------------------------------

  implicit none

  public :: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine dycoms2_rf02_tndcy( rho, & 
                                 rho_zm, rtm, rcm, exner,  & 
                                 err_code, &
                                 thlm_forcing, rtm_forcing,  & 
                                 Frad, radht, sclrm_forcing, &
                                 edsclrm_forcing )
    ! Description:
    !   Compute thlm_ls, rtm_ls, radiative heating rate, and cloud
    !   droplet number concentration as needed.

    ! References:
    !  ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
    !  Stevens, Bjorn, et al., (2003)
    !  Bull. Amer. Meteorol. Soc., 84, 579-593.
    !----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants, only: fstderr, Cp, rc_tol ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use parameters_radiation, only: rad_scheme  ! Variable(s)

    use error_code, only: clubb_rtm_level_not_found ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use stats_type, only: stat_update_var, stat_update_var_pt ! Procedure(s)

    use stats_variables, only:  & 
        iradht_LW, izi, sfc, zt, l_stats_samp ! Variable(s)

    use interpolation, only: lin_int

    use parameters_radiation, only: &
      F0,  & ! Variable(s)
      F1,  &
      kappa

    implicit none

    ! Constant parameters
    real, parameter ::  & 
      ls_div = 3.75e-6
!     kap    = 85.0       ! [m^2/kg]
!     F0     = 70.0,    & ! [W/m^2]
!     F1     = 22.0       ! [W/m^2]

    ! Input Variables

    real, intent(in), dimension(gr%nnzp) :: & 
      rho,    & ! Density on thermo. grid        [kg/m^3] 
      rho_zm, & ! Density on moment. grid        [kg/m^3]
      rtm,    & ! Total water mixing ratio       [kg/kg]
      rcm,    & ! Cloud water mixing ratio       [kg/kg]
      exner     ! Exner function.                [-]

    ! Input/Output Variables
    integer, intent(inout) :: err_code

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      thlm_forcing, & ! theta_l forcing                [K/s]
      rtm_forcing,  & ! r_t forcing                    [(kg/kg)/s] 
      Frad,         & ! Radiative flux                 [W/m^2]
      radht           ! Radiative heating rate         [K/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing    ! Passive scalar tendency        [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing  ! Eddy-passive scalar tendency   [units/s]

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      LWP,      & ! Liquid water path
      Heaviside

    real :: z_i

    integer :: k  ! Loop index

    if ( trim( rad_scheme ) == "simplified" ) then

      ! Radiation

      ! Compute liquid water path from top of the model
      ! We define liquid water path on momentum levels

      LWP(gr%nnzp) = 0.0
      do k = gr%nnzp-1, 1, -1
        LWP(k) = LWP(k+1) + rho(k+1)*rcm(k+1)/gr%dzt(k+1)
      end do  ! k = gr%nnzp..1

      ! Find the height of the isotherm rtm = 8.0 g/kg.

      k = 2
      do while ( k <= gr%nnzp .and. rtm(k) > 8.0e-3 )
        k = k + 1
      end do
      if ( k == gr%nnzp+1 .or. k == 2 ) then
        write(fstderr,*) "Identification of 8.0 g/kg level failed"
        write(fstderr,*) "Subroutine: dycoms2_rf02_tndcy. " & 
          // "File: dycoms2_rf02.F"
        write(fstderr,*) "k = ", k
        write(fstderr,*) "rtm(k) = ", rtm(k)
        err_code = clubb_rtm_level_not_found
        return
      end if
!  z_i = ( (gr%zt(k)-gr%zt(k-1))/(rtm(k)-rtm(k-1)) ) &
!      * (8.0e-3-rtm(k-1)) + gr%zt(k-1)

      z_i = lin_int( 8.0e-3, rtm(k), rtm(k-1), gr%zt(k), gr%zt(k-1) )
!         Compute the Heaviside step function for z - z_i.

      do k = 1, gr%nnzp, 1
        if ( gr%zm(k) - z_i  <  0.0 ) then
          Heaviside(k) = 0.0
        else if ( gr%zm(k) - z_i  ==  0.0 ) then
          Heaviside(k) = 0.5
        else if ( gr%zm(k) - z_i  >  0.0 ) then
          Heaviside(k) = 1.0
        end if
      end do

      ! Compute radiative flux profile (Frad).
      ! Radiative flux is defined on momentum levels.

      do k = 1, gr%nnzp, 1

        Frad(k) = F0 * EXP( -kappa * LWP(k) ) & 
                + F1 * EXP( -kappa * (LWP(1) - LWP(k)) )

        if ( Heaviside(k) > 0.0 ) then
          Frad(k) = Frad(k) & 
                  + rho_zm(k) * Cp * ls_div * Heaviside(k) & 
                    * ( 0.25 * ((gr%zm(k)-z_i)**(4.0/3.0)) & 
                  + z_i * ((gr%zm(k)-z_i)**(1.0/3.0)) )
        end if

      end do

      ! Compute the radiative heating rate.
      ! The radiative heating rate is defined on thermodynamic levels.

      do k = 2, gr%nnzp, 1
        radht(k) = ( 1.0 / exner(k) ) * ( -1.0/(Cp*rho(k)) ) & 
                 * ( Frad(k) - Frad(k-1) ) * gr%dzt(k)
      end do
      radht(1) = radht(2)

      if ( l_stats_samp ) then
        call stat_update_var( iradht_LW, radht, zt )
      end if

      thlm_forcing(1:gr%nnzp) = radht(1:gr%nnzp)

    else
      ! Compute heating rate elsewhere
      radht(1:gr%nnzp)        = 0.0
      thlm_forcing(1:gr%nnzp) = 0.0

    end if ! simplified

    ! Enter the final rtm tendency

    rtm_forcing(1:gr%nnzp) = 0.0

    ! Update surface statistics
    if ( l_stats_samp ) then

      call stat_update_var_pt( izi, 1, z_i, sfc )

    end if

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine dycoms2_rf02_tndcy


!----------------------------------------------------------------------

  subroutine dycoms2_rf02_sfclyr( um_sfc, vm_sfc,  & 
                                  upwp_sfc, vpwp_sfc,  & 
                                  wpthlp_sfc, wprtp_sfc, ustar, & 
                                  wpsclrp_sfc, wpedsclrp_sfc )

    use constants, only: Cp, Lv ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use surface_flux, only: &
        compute_ubar, compute_momentum_flux

    implicit none

    ! External
    intrinsic :: sqrt

    ! Constant parameters
    real, parameter ::  & 
      SH = 16.0, & 
      LH = 93.0

    ! Input Variables
    real, intent(in) ::  & 
      um_sfc,  & ! um(2) [m/s]
      vm_sfc  ! vm(2) [m/s]

    ! Output
    real, intent(out) ::  & 
      upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
      vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    real, intent(out), dimension(sclr_dim) ::  & 
      wpsclrp_sfc       ! w' scalar at surface [units m/s]

    real, intent(out), dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc     ! w' eddy-scalar at surface [units m/s]

    ! Local Variables
    real :: ubar  ! ? [m^2/s^2]?

    ! Declare the value of ustar.
    ustar = 0.25

    ubar = compute_ubar( um_sfc, vm_sfc )

    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    wpthlp_sfc = SH / (1.21 * Cp)
    wprtp_sfc  = LH / (1.21 * Lv)

    ! Let passive scalars be equal to rt and theta_l for now
    if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

    if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
    if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

    return
  end subroutine dycoms2_rf02_sfclyr

end module dycoms2_rf02
