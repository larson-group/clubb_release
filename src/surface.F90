!$Id$
module surface

  implicit none

  public :: prognose_soil_t_in_k, initialize_surface, get_veg_T_in_K

  integer, private, parameter :: leaf = 2

  real, private, dimension(leaf) :: deep_soil_t_in_k, sfc_soil_t_in_k, veg_t_in_k

  private

  contains

  !----------------------------------------------------------------------
  subroutine prognose_soil_t_in_k( dt, rho_sfc, &
                                   Frad_SW_down, fsin, &
                                   Frad_LW_down_sfc, wpthep )
    !
    !     Description:
    !
    !     This subroutine updates the surface and soil temp, soil heat flux, while assuming that net radiation
    !     and turbulent heat fluxes are already available from another subroutine.
    !
    !     The surface temperature (sfc_soil_t_in_k) is calculated  from
    !     the surface energy budget.
    !
    ! ************************************** 2
    !               *                   *
    !               *                   *
    !             --*         +         *--  1
    !               *                   *
    ! *** surface ************+************* 0
    !                        sfc_soil_t_in_k              l (evel)
    !
    !
    !     Heat conduction (in a homogeneous medium) can be described by
    !     the equation:
    !
    !              d(sfc_soil_t_in_k)/dt=ks*d(d(sfc_soil_t_in_k)/dz)/dz
    !
    !     In which ks is the soil thermal diffusity.
    !     We consider a semi half infinite medium, initially at the
    !     constant temperature deep_soil_t_in_k. if we vary the surface temperature
    !     sfc_soil_t_in_k sinussodally in time we can deduce a relation between the
    !     surface temperature, the soil heat flux (shf)  and deep_soil_t_in_k:
    !
    !              d(sfc_soil_t_in_k)/dt=c1*shf - c2*(sfc_soil_t_in_k-deep_soil_t_in_k)
    !
    !     However in reality the temperature deep_soil_t_in_k also varies in time, it
    !     may be calculated from:
    !
    !              d(deep_soil_t_in_k)/dt= c2*shf
    !
    !     The equations given above are analogous to those used by
    !     Deardorff  (1978).
    !
    !-----------------------------------------------------------------------

    use stats_precision, only: time_precision ! Variable(s)

    use stats_variables, only: l_stats_samp, sfc, iveg_t_sfc, it_sfc, ideep_T_sfc ! Variables

    use stats_type, only: stat_update_var_pt ! Procedure(s)

    use constants, only: pi, Cp, stefan_boltzmann ! Variable(s)

    implicit none

    ! Input variables

    real, intent(in) :: dt ! Current model timestep (Must be < 60s) [s]

    real, intent(in) :: &
    rho_sfc, &             ! Air density at the surface [kg/m^3]
    Frad_SW_down, &        ! SW downwelling flux        [W/m^3]
    fsin, &                ! Net Solar radiation        [W/m^3]
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
         shf, & ! Soil Heat Flux
         shfs,&
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


    Frad_LW_up_sfc = stefan_boltzmann * (veg_t_in_k(itf)**4)

    ! Calculate net radiation minus turbulent heat flux
    shfs = Frad_LW_down_sfc - Frad_LW_up_sfc - wpthep * rho_sfc * Cp + Frad_SW_down

    ! Calculate soil heat flux
    shf = 10.0 * ( veg_t_in_k(itf) - sfc_soil_t_in_k(itf) ) + 0.05 * fsin

    ! Update surf veg temp
    veg_t_in_k(itl) = veg_t_in_k(itf) + dt * 5.e-5 * ( shfs - shf )

    ! Update soil temp
    sfc_soil_t_in_k(itl) = sfc_soil_t_in_k(itf) & 
      + dt * ( c1 * shf - c2 * ( sfc_soil_t_in_k(itf)-deep_soil_t_in_k(itf) ) )

    ! Update deep soil temp
    deep_soil_t_in_k(itl) = deep_soil_t_in_k(itf) + dt * c3 * shf

    if( l_stats_samp ) then
      call stat_update_var_pt( iveg_t_sfc, 1, veg_T_in_K(1), sfc )
      call stat_update_var_pt( it_sfc, 1, sfc_soil_T_in_K(1), sfc )
      call stat_update_var_pt( ideep_t_sfc, 1, deep_soil_T_in_K(1), sfc )
    end if

    veg_t_in_k = veg_t_in_k(itl)
    sfc_soil_t_in_k = sfc_soil_t_in_k(itl)
    deep_soil_t_in_k = deep_soil_t_in_k(itl)

    return
  end subroutine prognose_soil_t_in_k

  !------------------------------------------------------------------------------------------
  subroutine initialize_surface( initial_veg_t_in_k_in, initial_sfc_soil_t_in_k_in, &
                                 initial_deep_soil_t_in_k_in )
    !
    !       Description: This subroutine sets the initial state of the temperature
    !       of the ground at vegetation, surface soil, and deep soil levels.
    !
    !-----------------------------------------------------------------------------------------
    implicit none

    ! Input variables
    real, intent(in) :: &
    initial_veg_t_in_k_in, &            ! Initial temperature of vegetation layer [K]
    initial_sfc_soil_t_in_k_in, &       ! Initial temperature of surface soil layer [K]
    initial_deep_soil_t_in_k_in         ! Initial temperature of deep soil layer [K]

    !-----------------------------------------------------------------------------------------

    veg_t_in_k = initial_veg_t_in_k_in
    sfc_soil_t_in_k = initial_sfc_soil_t_in_k_in
    deep_soil_t_in_k = initial_deep_soil_t_in_k_in

  end subroutine initialize_surface

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
end module surface
