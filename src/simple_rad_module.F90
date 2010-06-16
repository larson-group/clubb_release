! $Id$
module simple_rad_module

  implicit none

  public :: simple_rad, simple_rad_bomex

  private :: liq_water_path

  private

  contains

!-------------------------------------------------------------------------------
  subroutine simple_rad( rho, rho_zm, rtm, rcm, exner,  & 
                         err_code, Frad, radht )
! Description:
!   A simplified radiation driver
! References:
!   None
!-------------------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants, only: fstderr, Cp, rc_tol ! Variable(s)

    use error_code, only: clubb_rtm_level_not_found ! Variable(s)

    use stats_type, only: stat_update_var, stat_update_var_pt ! Procedure(s)

    use stats_variables, only:  & 
        iradht_LW, izi, sfc, zt, l_stats_samp ! Variable(s)

    use interpolation, only: lin_int ! Procedure(s)

    use parameters_radiation, only: &
      F0,  & ! Variable(s)
      F1,  &
      l_rad_above_cloud, &
      kappa

    implicit none

    ! External
    intrinsic :: exp

    ! Constant parameters

    real, parameter ::  & 
      ls_div = 3.75e-6

    ! Input Variables

    real, intent(in), dimension(gr%nnzp) :: & 
      rho,    & ! Density on thermodynamic grid  [kg/m^3] 
      rho_zm, & ! Density on momentum grid       [kg/m^3]
      rtm,    & ! Total water mixing ratio       [kg/kg]
      rcm,    & ! Cloud water mixing ratio       [kg/kg]
      exner     ! Exner function.                [-]

    integer, intent(inout) :: err_code

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      Frad,         & ! Radiative flux                 [W/m^2]
      radht           ! Radiative heating rate         [K/s]

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      LWP,      & ! Liquid water path
      Heaviside

    real :: z_i

    integer :: k

    ! ---- Begin Code ----

    lwp(1:gr%nnzp) = liq_water_path( gr%nnzp, rho, rcm, gr%invrs_dzt )

    do k = 1, gr%nnzp, 1

      if ( F1 /= 0 ) then
        Frad(k) = F0 * exp( -kappa * LWP(k) ) & 
                + F1 * exp( -kappa * (LWP(1) - LWP(k)) )

      else ! Mathematically equivalent to the above, but computationally cheaper
        Frad(k) = F0 * exp( -kappa * LWP(k) )

      end if ! F1 /= 0

    end do 

    if ( l_rad_above_cloud ) then
      ! Find the height of the isotherm rtm = 8.0 g/kg.

      k = 2
      do while ( k <= gr%nnzp .and. rtm(k) > 8.0e-3 )
        k = k + 1
      end do
      if ( k == gr%nnzp+1 .or. k == 2 ) then
        write(fstderr,*) "Identification of 8.0 g/kg level failed"
        write(fstderr,*) "Subroutine: simple_rad. " & 
          // "File: simple_rad_module.F90"
        write(fstderr,*) "k = ", k
        write(fstderr,*) "rtm(k) = ", rtm(k)
        err_code = clubb_rtm_level_not_found
        return
      end if

      z_i = lin_int( 8.0e-3, rtm(k), rtm(k-1), gr%zt(k), gr%zt(k-1) )

      ! Compute the Heaviside step function for z - z_i.
      do k = 1, gr%nnzp, 1
        if ( gr%zm(k) - z_i  <  0.0 ) then
          Heaviside(k) = 0.0
        else if ( gr%zm(k) - z_i  ==  0.0 ) then
          Heaviside(k) = 0.5
        else if ( gr%zm(k) - z_i  >  0.0 ) then
          Heaviside(k) = 1.0
        end if
      end do

      do k = 1, gr%nnzp, 1
        if ( Heaviside(k) > 0.0 ) then
          Frad(k) = Frad(k) & 
                  + rho_zm(k) * Cp * ls_div * Heaviside(k) & 
                    * ( 0.25 * ((gr%zm(k)-z_i)**(4.0/3.0)) & 
                  + z_i * ((gr%zm(k)-z_i)**(1.0/3.0)) )
        end if
      end do ! k=1..gr%nnzp

      ! Update surface statistics
      if ( l_stats_samp ) then

        call stat_update_var_pt( izi, 1, z_i, sfc )

      end if

    end if ! l_rad_above_cloud

    ! Compute the radiative heating rate.
    ! The radiative heating rate is defined on thermodynamic levels.

    do k = 2, gr%nnzp, 1
      radht(k) = ( 1.0 / exner(k) ) * ( -1.0/(Cp*rho(k)) ) & 
               * ( Frad(k) - Frad(k-1) ) * gr%invrs_dzt(k)
    end do
    radht(1) = radht(2)

    if ( l_stats_samp ) then
      call stat_update_var( iradht_LW, radht, zt )
    end if

    return
  end subroutine simple_rad

!-------------------------------------------------------------------------------
  subroutine simple_rad_bomex( radht )
! Description:
!   Compute radiation as per the GCSS BOMEX specification.
! References:
!   <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!-------------------------------------------------------------------------------
    use grid_class, only: gr ! Type

    implicit none

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) :: & 
      radht  ! Radiative heating rate [K/s]

    ! Local Variables
    integer :: k

    ! ---- Begin Code ----

    ! Radiative theta-l tendency
    do k = 2, gr%nnzp

      if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
        radht(k) = -2.315e-5
      else if ( gr%zt(k) >= 1500. .and. gr%zt(k) < 2500. ) then
        radht(k) & 
          = - 2.315e-5  & 
            + 2.315e-5  & 
              * ( gr%zt(k) - 1500. ) / ( 2500. - 1500. )
      else
        radht(k) = 0.
      end if

    end do ! k=2..gr%nnzp

    ! Boundary condition
    radht(1) = 0.0

    return
  end subroutine simple_rad_bomex
!-------------------------------------------------------------------------------
  pure function liq_water_path( nnzp, rho, rcm, dzt )

! Description:
!   Compute liquid water path
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: nnzp

    real, intent(in), dimension(nnzp) :: &
      rho, & ! Air Density                      [kg/m^3]
      rcm, & ! Cloud water mixing ratio         [kg/kg]
      dzt    ! Inverse of distance per level    [1/m]

    ! Output Variables
    real, dimension(nnzp) :: &
      liq_water_path ! Liquid water path

    integer :: k

    ! ---- Begin Code ----

    liq_water_path(nnzp) = 0.0

    ! Liquid water path is defined on the intermediate model levels between the
    ! rcm and rho levels (i.e. the momentum levels in CLUBB).
    do k = nnzp-1, 1, -1
       liq_water_path(k) = liq_water_path(k+1) + rcm(k+1)*rho(k+1) / dzt(k+1)
    end do ! k = nnzp..1

    return
  end function liq_water_path

end module simple_rad_module
