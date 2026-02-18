!$Id$
module soil_vegetation

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: advance_soil_veg, initialize_soil_veg

  logical, public :: l_soil_veg

!$omp threadprivate(l_soil_veg)

  private

  contains

  !----------------------------------------------------------------------
  subroutine advance_soil_veg( ngrdcol, dt, rho_sfc, &
                               Frad_SW_up_sfc, Frad_SW_down_sfc, &
                               Frad_LW_down_sfc, &
                               wpthlp_sfc, wprtp_sfc, p_sfc, &
                               stats,         &
                               deep_soil_T_in_K, sfc_soil_T_in_K, &
                               veg_T_in_K )
    !
    !     Description:
    !
    !     This subroutine updates the surface and soil temp, soil heat flux, 
    !     while assuming that net radiation and turbulent heat fluxes are 
    !     already available from another subroutine.
    !
    !     The surface temperature (sfc_soil_T_in_K) is calculated  from
    !     the surface energy budget.
    !
    ! ************************************** 2
    !               *                   *
    !               *                   *
    !             --*         +         *--  1
    !               *                   *
    ! *** surface ************+************* 0
    !                        sfc_soil_T_in_K              l (evel)
    !
    !
    !     Heat conduction (in a homogeneous medium) can be described by
    !     the equation:
    !
    !              d(sfc_soil_T_in_K)/dt=ks*d(d(sfc_soil_T_in_K)/dz)/dz
    !
    !     In which ks is the soil thermal diffusity.
    !     We consider a semi half infinite medium, initially at the
    !     constant temperature deep_soil_T_in_K. if we vary the surface temperature
    !     sfc_soil_T_in_K sinussodally in time we can deduce a relation between the
    !     surface temperature, the soil heat flux (soil_heat_flux)  and deep_soil_T_in_K:
    !
    !              d(sfc_soil_T_in_K)/dt=c1*soil_heat_flux - c2*(sfc_soil_T_in_K-deep_soil_T_in_K)
    !
    !     However in reality the temperature deep_soil_T_in_K also varies in time, it
    !     may be calculated from:
    !
    !              d(deep_soil_T_in_K)/dt= c2*soil_heat_flux
    !
    !     The equations given above are analogous to those used by
    !     Deardorff  (1978).
    ! 
    !     Reference:
    !     Duynkerke, Peter G.  "Radiation Fog: A Comparison of Model 
    !     Simulation with Detailed Observations"
    !     (February 1991) _Monthly Weather Review_ Vol 119, p. 324-341.
    !     
    !
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Constant

    use constants_clubb, only: &
      pi, &
      Cp, &
      Lv, &
      kappa, &
      p0, &
      stefan_boltzmann ! Variable(s)

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    implicit none

    ! This subroutine does not produce any output variables. Instead the module
    ! variables listed below are updated.
    !
    ! veg_T_in_K_in, &            ! Temperature of vegetation layer [K]
    ! sfc_soil_T_in_K_in, &       ! Temperature of surface soil layer [K]
    ! deep_soil_T_in_K_in         ! Temperature of deep soil layer [K]

    ! External

    intrinsic :: sqrt, exp

    ! Input variables
    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), intent(in) :: &
      dt ! Current model timestep (Must be < 60s) [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      rho_sfc, &             ! Air density at the surface   [kg/m^3]
      Frad_SW_up_sfc, &      ! SW Net                       [W/m^3]
      Frad_SW_down_sfc, &    ! SW radiative upwelling flux  [W/m^2]
      Frad_LW_down_sfc, &    ! LW downwelling flux          [W/m^2]
      wpthlp_sfc, &
      wprtp_sfc, &
      p_sfc

    type(stats_type), intent(inout) :: &
      stats
    
    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      deep_soil_T_in_K, &
      sfc_soil_T_in_K, &
      veg_T_in_K

    ! Local variables

    real( kind = core_rknd ) :: &
      cs, &  ! soil heat capacity              [Jg/K]
      ks, &  ! soil heat diffusivity           [m^2/s]
      rs, &  ! soil density                    [g/m^3]
      c1, &  ! coefficient in force restore 1
      c2, &  ! coefficient in force restore 2
      c3, &  ! coefficient in force restore 3
      d1, &
      veg_heat_flux, &                         
      Frad_LW_up_sfc, & ! LW upwelling flux [W/m2]
      wpthep            ! Turbulent Flux of equivalent potential temperature   [K]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      soil_heat_flux ! Soil Heat flux [W/m^2]

    integer :: i

    !--------------------------------- Begin Code ---------------------------------

    !----------------------------
    !  Soil parameters
    !---------------------------

    cs = 2.00e3_core_rknd  ! cs
    rs = 1.00e3_core_rknd  ! ps 
    ks = 2.00e-7_core_rknd ! as
    d1 = sqrt( ks*3600.e0_core_rknd*24.e0_core_rknd ) ! Known magic number
    c1 = 2.e0_core_rknd*sqrt( pi )/(rs*cs*d1) 
    c2 = 2.e0_core_rknd*pi/(3600.e0_core_rknd*24.e0_core_rknd) ! Omega - known magic number
    c3 = sqrt(pi*2.e0_core_rknd)/(exp( pi/4.e0_core_rknd )*rs*cs* &
                     sqrt( ks*3600.e0_core_rknd*24.e0_core_rknd* &
                     365.e0_core_rknd )) ! Known magic number

    if ( stats%l_sample ) then
      call stats_update( "veg_T_in_K", veg_T_in_K, stats )
      call stats_update( "sfc_soil_T_in_K", sfc_soil_T_in_K, stats )
      call stats_update( "deep_soil_T_in_K", deep_soil_T_in_K, stats )
    end if

    do i = 1, ngrdcol 
    
      Frad_LW_up_sfc = stefan_boltzmann * (veg_T_in_K(i)**4)

      wpthep = wpthlp_sfc(i) + (Lv/Cp) * ((p0/p_sfc(i))**kappa) * wprtp_sfc(i)

      ! Calculate net radiation minus turbulent heat flux
      veg_heat_flux = Frad_LW_down_sfc(i) - Frad_LW_up_sfc &
                      - wpthep * rho_sfc(i) * Cp + ( Frad_SW_down_sfc(i) - Frad_SW_up_sfc(i) )

      ! Calculate soil heat flux
      ! Duynkerke (1991) used a coefficient of 3.0 W/m^2*K, not 10.0 W/m^2*K
      !
      ! Equation 19 p.328
      
      soil_heat_flux(i) = 10.0_core_rknd * ( veg_T_in_K(i) - sfc_soil_T_in_K(i) ) &
                          + 0.05_core_rknd * Frad_SW_down_sfc(i) ! Known magic number

      ! Update surf veg temp
      veg_T_in_K(i) = veg_T_in_K(i) + dt * 5.e-5_core_rknd &
                      * ( veg_heat_flux - soil_heat_flux(i) ) ! Known magic number

      ! Update soil temp
      sfc_soil_T_in_K(i) = sfc_soil_T_in_K(i) + dt &
                           * ( c1 * soil_heat_flux(i) - c2 &
                                    * ( sfc_soil_T_in_K(i) - deep_soil_T_in_K(i) ) )

      ! Update deep soil temp
      deep_soil_T_in_K(i) = deep_soil_T_in_K(i) + dt * c3 * soil_heat_flux(i)
    end do

    if ( stats%l_sample ) then
      do i = 1, ngrdcol
          call stats_update( "soil_heat_flux", soil_heat_flux(i), stats, i )
      end do
    end if

    return

  end subroutine advance_soil_veg

  !-----------------------------------------------------------------------------
  subroutine initialize_soil_veg( ngrdcol, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K )
  ! Description:
  !   Sets some default values for the soil scheme
  ! References:
  !   None
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Constant

    implicit none

    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
      deep_soil_T_in_K, &
      sfc_soil_T_in_K, &
      veg_T_in_K       

    ! ---- Begin Code ----

    ! These default values are the values for gabls3
    deep_soil_T_in_K = 288.58_core_rknd
    sfc_soil_T_in_K  = 300._core_rknd
    veg_T_in_K       = 300._core_rknd

    ! Disable this for most cases
    l_soil_veg       = .false.

    return
  end subroutine initialize_soil_veg

end module soil_vegetation
