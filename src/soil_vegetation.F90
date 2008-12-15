!$Id$
module soil_vegetation

  implicit none

  public :: advance_soil_veg, initialize_soil_veg, get_veg_T_in_K

  integer, private, parameter :: leaf = 2

  real, private, dimension(leaf) :: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

  private

  contains

  !----------------------------------------------------------------------
  subroutine advance_soil_veg( dt, rho_sfc, &
                                   Frad_SW_net, Frad_SW_down_sfc, &
                                   Frad_LW_down_sfc, wpthep )
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

    use stats_precision, only: time_precision ! Variable(s)

    use stats_variables, only: l_stats_samp, sfc, &
                                iveg_T_in_K, isfc_soil_T_in_K, ideep_soil_T_in_K ! Variables

    use stats_type, only: stat_update_var_pt ! Procedure(s)

    use constants, only: pi, Cp, stefan_boltzmann ! Variable(s)

    implicit none

    ! This subroutine does not produce any output variables. Instead the module
    ! variables listed below are updated.
    !
    ! veg_T_in_K_in, &            ! Temperature of vegetation layer [K]
    ! sfc_soil_T_in_K_in, &       ! Temperature of surface soil layer [K]
    ! deep_soil_T_in_K_in         ! Temperature of deep soil layer [K]

    ! Input variables

    real, intent(in) :: dt ! Current model timestep (Must be < 60s) [s]

    real, intent(in) :: &
    rho_sfc, &             ! Air density at the surface [kg/m^3]
    Frad_SW_net, &         ! SW Net                     [W/m^3]
    Frad_SW_down_sfc, &    ! SW downwelling flux        [W/m^3]
    Frad_LW_down_sfc, &    ! LW downwelling flux        [W/m^2]
    wpthep                 ! Turbulent Flux of equivalent potential temperature   [K]

    ! Local variables

    real cs, &  ! soil heat capacity
         ks, &  ! soil heat diffusivity
         rs, &  ! soil density
         c1, &  ! coefficient in force restore 1
         c2, &  ! coefficient in force restore 2
         c3, &  ! coefficient in force restore 3
         d1, &
         soil_heat_flux, & ! Soil Heat Flux
         veg_heat_flux,&
         Frad_LW_up_sfc ! LW upwelling flux [W/m2]

    integer :: &
    itf, & ! Current Timestep [-]
    itl    ! Next Timestep    [-]

    !----------------------------
    !  Soil parameters
    !---------------------------

    cs=2.00e3
    rs=1.00e3
    ks=2.00e-7
    d1=sqrt(ks*3600.e0*24.e0)
    c1=2.e0*sqrt(pi)/(rs*cs*d1)
    c2=2.e0*pi/(3600.e0*24.e0)
    c3=sqrt(pi*2.e0)/(exp(pi/4.e0)*rs*cs*sqrt(ks*3600.e0*24.e0* &
                     365.e0))

    itl = 2
    itf = 1


    Frad_LW_up_sfc = stefan_boltzmann * (veg_T_in_K(itf)**4)

    ! Calculate net radiation minus turbulent heat flux
    veg_heat_flux = Frad_LW_down_sfc - Frad_LW_up_sfc - wpthep * rho_sfc * Cp + Frad_SW_net

    ! Calculate soil heat flux
    ! Duynkerke (1990) used a coefficient of 3.0, not 10.0
    
    soil_heat_flux = 10.0 * ( veg_T_in_K(itf) - sfc_soil_T_in_K(itf) ) + 0.05 * Frad_SW_down_sfc

    ! Update surf veg temp
    veg_T_in_K(itl) = veg_T_in_K(itf) + dt * 5.e-5 * ( veg_heat_flux - soil_heat_flux )

    ! Update soil temp
    sfc_soil_T_in_K(itl) = sfc_soil_T_in_K(itf) & 
      + dt * ( c1 * soil_heat_flux - c2 * ( sfc_soil_T_in_K(itf)-deep_soil_T_in_K(itf) ) )

    ! Update deep soil temp
    deep_soil_T_in_K(itl) = deep_soil_T_in_K(itf) + dt * c3 * soil_heat_flux

    if( l_stats_samp ) then
      call stat_update_var_pt( iveg_T_in_K, 1, veg_T_in_K(1), sfc )
      call stat_update_var_pt( isfc_soil_T_in_K, 1, sfc_soil_T_in_K(1), sfc )
      call stat_update_var_pt( ideep_soil_T_in_K, 1, deep_soil_T_in_K(1), sfc )
    end if

    veg_T_in_K = veg_T_in_K(itl)
    sfc_soil_T_in_K = sfc_soil_T_in_K(itl)
    deep_soil_T_in_K = deep_soil_T_in_K(itl)

    return
  end subroutine advance_soil_veg

  !------------------------------------------------------------------------------------------
  subroutine initialize_soil_veg( initial_veg_T_in_K_in, initial_sfc_soil_T_in_K_in, &
                                 initial_deep_soil_T_in_K_in )
    !
    !       Description: This subroutine sets the initial state of the temperature
    !       of the ground at vegetation, surface soil, and deep soil levels.
    !
    !-----------------------------------------------------------------------------------------
    implicit none

    ! Input variables
    real, intent(in) :: &
    initial_veg_T_in_K_in, &            ! Initial temperature of vegetation layer [K]
    initial_sfc_soil_T_in_K_in, &       ! Initial temperature of surface soil layer [K]
    initial_deep_soil_T_in_K_in         ! Initial temperature of deep soil layer [K]

    !-----------------------------------------------------------------------------------------

    veg_T_in_K = initial_veg_T_in_K_in
    sfc_soil_T_in_K = initial_sfc_soil_T_in_K_in
    deep_soil_T_in_K = initial_deep_soil_T_in_K_in

  end subroutine initialize_soil_veg

  !--------------------------------------------------------------------------------------------
  real function get_veg_T_in_K()
    !
    !     Description: This function serves as accessor to sfc_soil_T_in_K. It
    !     abstracts the information that sfc_soil_T_in_K is stored as an array.
    !
    !--------------------------------------------------------------------------------------------
    implicit none

    get_veg_T_in_K = veg_T_in_K(1)

    return

  end function
end module soil_vegetation
