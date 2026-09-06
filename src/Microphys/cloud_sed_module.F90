!---------------------------------------------------------------------------
!$Id$
!===============================================================================
module cloud_sed_module

  implicit none

  public :: cloud_drop_sed

  private ! Default Scope

  !=====================================================================
  ! NOTE:  ADDITION OF RAIN EFFECTS AND CLOUD WATER SEDIMENTATION ON
  !        RTM AND THLM.
  !=====================================================================
  !
  ! Equations:  rtm = rvm + rcm;
  !             thlm = thm - ( Lv / (Cp*exner) ) * rcm
  !
  ! When water condenses, latent heat is given off and theta (thm)
  ! increases by a factor of ( Lv / (Cp*exner) ) * rcm(condensed).
  ! The opposite effect occurs with evaporation.
  !
  !======================================================================
  !||     Effect      |  rvm  |  rcm  |  rtm  |  rrm  |  thm  |  thlm  ||
  !|====================================================================|
  !|| Sedimentation   |       |       |       |       |       |        ||
  !|| Effects of      | stays | incr. | incr. | stays | stays | decr.  ||
  !|| Cloud Water.    | same  |       |       | same  | same  |        ||
  !|| sed_rcm > 0     |       |       |       |       |       |        ||
  !|====================================================================|
  !|| Evaporation     |       |       |       |       |       |        ||
  !|| of rain to      | incr. | stays | incr. | decr. | decr. | decr.  ||
  !|| water vapor.    |       | same  |       |       |       |        ||
  !|| cond_rrainm < 0 |       |       |       |       |       |        ||
  !|--------------------------------------------------------------------|
  !|| Autoconversion  |       |       |       |       |       |        ||
  !|| of cloud water  | stays | decr. | decr. | incr. | stays | incr.  ||
  !|| to rain water.  | same  |       |       |       | same  |        ||
  !|| auto_rrainm > 0 |       |       |       |       |       |        ||
  !|--------------------------------------------------------------------|
  !|| Accretion of    |       |       |       |       |       |        ||
  !|| cloud water by  | stays | decr. | decr. | incr. | stays | incr.  ||
  !|| rain water.     | same  |       |       |       | same  |        ||
  !|| accr_rrainm > 0 |       |       |       |       |       |        ||
  !======================================================================
  !
  ! Note: In HOC, cond_rrainm will always be either negative or zero.
  !
  ! Overall effects of rain and cloud water sedimentation:
  !
  ! (drtm/dt)|_t  = (drtm/dt)|_0 + sed_rcm
  !                 - cond_rrainm - auto_rrainm - accr_rrainm
  !
  ! (dthlm/dt)|_t = (dthlm/dt)|_0 - ( Lv / (Cp*exner) ) * sed_rcm
  !                 + ( Lv / (Cp*exner) )
  !                   * ( cond_rrainm + auto_rrainm + accr_rrainm )

contains

  !=============================================================================
  subroutine cloud_drop_sed( gr, ngrdcol, rcm, Ncm, &
                             rho_zm, rho, exner, sigma_g, &
                             stats, rcm_mc, thlm_mc )

    ! Description:
    ! Account for cloud droplet sedimentation in all model columns.
    !
    ! Sedimentation flux of cloud droplets should be treated by assuming a
    ! log-normal size distribution of droplets falling in a Stoke's regime, in
    ! which the sedimentation flux (Fcsed) is given by:
    !
    ! Sedimentation Flux = constant
    !                      * [ ( 3 / ( 4 * pi * rho_lw * NcV ) )^(2/3) ]
    !                      * [ ( rho * rc )^(5/3) ]
    !                      * EXP[ 5 * ( ( LOG( sigma_g ) )^2 ) ];
    !
    ! where constant = 1.19 x 10^8 (m^-1 s^-1) and sigma_g is the geometric
    ! standard deviation of cloud droplets falling in a Stoke's regime.
    !
    ! When written for a mass-dependent cloud droplet concentration, Nc:
    !
    ! Sedimentation Flux = constant
    !                      * [ ( 3 / ( 4 * pi * rho_lw * Nc * rho ) )^(2/3) ]
    !                      * [ ( rho * rc )^(5/3) ]
    !                      * EXP[ 5 * ( ( LOG( sigma_g ) )^2 ) ].
    !
    ! According to the above equation, sedimentation flux, Fcsed, is defined
    ! positive downwards.  Therefore, 
    !
    ! (drc/dt)|_Fcsed = (1.0/rho) * d(Fcsed)/dz.

    ! References:
    ! http://journals.ametsoc.org/doi/abs/10.1175/2008MWR2582.1
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid,       & ! Type
        zt2zm_api,  & ! Procedure(s)
        ddzm

    use constants_clubb, only: &
        five,       & ! Constant(s)
        four,       &
        three,      &
        two_thirds, &
        one,        &
        zero,       &
        rho_lw,     &
        pi,         &
        Cp,         &
        Lv
 
    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    implicit none

    ! External
    intrinsic :: exp, log

    !--------------------------- Input Variables ---------------------------
    type (grid), intent(in) :: &
      gr                       ! Model grid [-]

    integer, intent(in) :: &
      ngrdcol                  ! Number of model columns [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzt) :: &
      rcm,    & ! Mean cloud water mixing ratio        [kg/kg]
      Ncm,    & ! Mean cloud droplet concentration     [num/kg]
      rho,    & ! Density on thermodynamic levels      [kg/m^3]
      exner     ! Exner function                       [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzm) :: &
      rho_zm    ! Density on momentum levels           [kg/m^3]

    real( kind = core_rknd ), intent(in) :: &
      sigma_g   ! Geometric standard deviation of cloud droplets   [-]

    !------------------------ Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: &
      stats     ! Statistics runtime context [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,gr%nzt) :: &
      rcm_mc,  & ! r_c tendency due to microphysics     [kg/kg)/s] 
      thlm_mc    ! thlm tendency due to microphysics    [K/s] 

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      Fcsed      ! Cloud water sedimentation flux         [kg/(m^2 s)]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      sed_rcm    ! d(rcm)/dt due to cloud sedimentation   [kg/(m^2 s)]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      rcm_zm, & ! Cloud water mixing ratio on momentum levels [kg/kg]
      Ncm_zm    ! Cloud droplet concentration on momentum levels [num/kg]

    integer :: &
      i, & ! Column loop index [-]
      k    ! Vertical-level loop index [-]

    !--------------------------- Begin Code ---------------------------

    rcm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, rcm )
    Ncm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, Ncm )

    ! Define cloud water sedimentation flux on momentum levels.
    do k = 2, gr%nzm-1
      do i = 1, ngrdcol
        if ( rcm_zm(i,k) > zero .and. Ncm_zm(i,k) > zero ) then
          Fcsed(i,k) = 1.19E8_core_rknd &
            * ( ( three / ( four * pi * rho_lw &
                * Ncm_zm(i,k) * rho_zm(i,k) ) )**two_thirds ) &
            * ( ( rho_zm(i,k) * rcm_zm(i,k) )**(five/three) ) &
            * exp( five * ( log( sigma_g )**2 ) ) ! See Ackerman - eq. no. 7
        else
          Fcsed(i,k) = zero
        endif
      enddo
    enddo

    ! Boundary conditions.
    do i = 1, ngrdcol
      Fcsed(i,1) = zero
      Fcsed(i,gr%nzm) = zero
    enddo

    ! Find drc/dt due to cloud water sedimentation flux.
    ! This value is defined on thermodynamic levels.
    ! Fcsed units:  kg (liquid) / [ m^2 * s ]
    ! Multiply by Lv for units of W / m^2.
    ! sed_rcm units:  [ kg (liquid) / kg (air) ] / s
    !$acc data copyin( gr, gr%invrs_dzt, Fcsed ) copyout( sed_rcm )
    sed_rcm = ddzm( gr%nzm, gr%nzt, ngrdcol, gr, Fcsed )
    !$acc end data

    sed_rcm = ( one / rho ) * sed_rcm

    if ( stats%l_sample ) then
      call stats_update( "sed_rcm", sed_rcm, stats )
      call stats_update( "Fcsed", Fcsed, stats )
    end if

    ! + thlm/rtm_microphysics -- cloud water sedimentation.
    ! Code addition by Brian for cloud water sedimentation.

    do k = 1, gr%nzt
      do i = 1, ngrdcol
        rcm_mc(i,k) = rcm_mc(i,k) + sed_rcm(i,k)
        thlm_mc(i,k) = thlm_mc(i,k) - ( Lv / ( Cp * exner(i,k) ) ) * sed_rcm(i,k)
      enddo
    enddo


    return

  end subroutine cloud_drop_sed

!===============================================================================

end module cloud_sed_module
