!----------------------------------------------------------------------
! $Id$
module bomex

!       Description:
!       Contains subroutines for the GCSS BOMEX case.
!
!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

  implicit none

  public :: bomex_tndcy, bomex_sfclyr

  private ! Default Scope

  contains

!----------------------------------------------------------------------
  subroutine bomex_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                          gr, rtm, & 
                          thlm_forcing, rtm_forcing, & 
                          sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------


    use grid_class, only: &
        zt2zm, & ! Procedure(s)
        grid     ! Type

    use spec_hum_to_mixing_ratio, only: &
        force_spec_hum_to_mixing_ratio ! Procedure(s)

    use array_index, only: &
        sclr_idx_type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzt) :: &
      rtm    ! Total water mixing ratio (thermodynamic levels)        [kg/kg]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing        [units vary/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing [units vary/s]

    !--------------------- Local Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      qtm_forcing  ! Specified total water spec. humidity tendency    [kg/kg/s]

    integer :: i, k

    !--------------------- Begin Code ---------------------

    ! Large scale advective moisture tendency
    ! The BOMEX specifications give large-scale advective moisture tendency in
    ! terms of total water specific humidity.
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        if ( gr%zt(i,k) >= 0._core_rknd .and. gr%zt(i,k) < 300._core_rknd ) then
          qtm_forcing(i,k) = -1.2e-8_core_rknd
        else if ( gr%zt(i,k) >= 300._core_rknd .and. gr%zt(i,k) < 500._core_rknd ) then
          qtm_forcing(i,k)  & 
            = - 1.2e-8_core_rknd  & 
                * ( 1._core_rknd - ( gr%zt(i,k) - 300._core_rknd )/ &
                ( 500._core_rknd - 300._core_rknd ) ) !Known magic number
        else
          qtm_forcing(i,k) = 0._core_rknd
        end if

      end do
    end do

    ! Convert forcings from terms of total water specific humidity to terms of
    ! total water mixing ratio.
    call force_spec_hum_to_mixing_ratio( ngrdcol, gr%nzt, rtm, qtm_forcing, rtm_forcing )

    ! Test scalars with thetal and rt if desired
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        ! Radiative theta-l tendency
        thlm_forcing(i,k) = 0.0_core_rknd

        if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_forcing(i,k)

        if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_forcing(i,k)
      end do
    end do

    return
  end subroutine bomex_tndcy

!----------------------------------------------------------------------
  subroutine bomex_sfclyr( ngrdcol, time, rtm_sfc, & 
                           wpthlp_sfc, wprtp_sfc, ustar )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS BOMEX specifications

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

    use spec_hum_to_mixing_ratio, only: &
        flux_spec_hum_to_mixing_ratio ! Procedure(s)

    use time_dependent_input, only: &
        time_select, &  ! Procedure(s)
        wpthlp_sfc_given, wpqtp_sfc_given, &
        time_sfc_given ! Variable(s)

    use interpolation, only: &
        linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real(time_precision), intent(in) ::  & 
      time    ! the current time [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      rtm_sfc    ! rtm(2) [kg/kg]

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      wpqtp_sfc     ! w'q_t' at (1)         [(m kg)/(s kg)]   

    real( kind = core_rknd ) :: &
      time_frac     ! The time fraction used for interpolation.

    integer :: before_time, after_time, i ! The time bounds used for interpolation


    ! Compute heat and moisture fluxes

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                      before_time, after_time, time_frac)

    do i = 1, ngrdcol

      ! Declare the value of ustar.
      ustar(i) = 0.28_core_rknd

      wpthlp_sfc(i) = linear_interp_factor( time_frac, wpthlp_sfc_given(after_time), &
                                            wpthlp_sfc_given(before_time) )

    ! The BOMEX specifications give surface moisture flux in terms of total water
    ! specific humidity.
      wpqtp_sfc(i)  = linear_interp_factor( time_frac, wpqtp_sfc_given(after_time), &
                                            wpqtp_sfc_given(before_time) )
    end do


    ! Convert flux from terms of total water specific humidity to terms of total
    ! water mixing ratio.
    call flux_spec_hum_to_mixing_ratio( ngrdcol, rtm_sfc, wpqtp_sfc, wprtp_sfc )

    return
  end subroutine bomex_sfclyr

end module bomex
