!$Id$
module surface

  implicit none

  public :: prognose_soil_t_in_k, initialize_surface

  real, private :: initial_veg_t_in_k, initial_sfc_soil_t_in_k, initial_deep_soil_t_in_k

  private


  contains

  !----------------------------------------------------------------------
  subroutine prognose_soil_t_in_k( time_start, time_current, dt, itf, itl, rho_sfc, &
                                   fs, fsin, &
                                   Frad_LW_down_sfc, wpthep, &
                                   veg_t_in_k, sfc_soil_t_in_k, deep_soil_t_in_k )
    !
    !
    !     In this subroutine tq and qw at the surface are calculated as a
    !     function of the surface temperature (sfc_soil_t_in_k). sfc_soil_t_in_k is calculated  from
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

    use constants, only: pi, Cp, stefan_boltzmann ! Variable(s)

    implicit none

    integer, parameter :: leaf = 2

    ! Input variables

    real(kind=time_precision), intent(in) :: &
    time_start, &  ! Time of start of simulation [s]
    time_current   ! Current time of simulation  [s]

    real, intent(in) :: dt ! Current model timestep (Must be < 60s) [s]

    integer, intent(in) :: &
    itf, & ! Current Timestep [-]
    itl    ! Next Timestep    [-]

    real, intent(in) :: &
    rho_sfc, & 
    fs, &
    fsin, &
    Frad_LW_down_sfc, &
    wpthep

    ! Input/Output variables
    real, dimension(leaf), intent(inout):: deep_soil_t_in_k, sfc_soil_t_in_k ,veg_t_in_k

    ! Local variables

    real cs,ks,rs,c1,c2,c3,d1,shf,shfs,Frad_LW_up_sfc

    integer :: it

    ! Functionality of the code:
    ! Update surface and soil temp, soil heat flux, while assuming that net radiation
    ! and turbulent heat fluxes are already available from another subroutine.
    !
    !
    !     meaning of variables
    !        Frad_LW_down_sfc= downwelling longwaverad
    !        Frad_LW_up_sfc= upwelling shortwaverad
    !        shf= soil heatflux
    !        rho_sfc(0)= air density at surface
    !        kappa= von karman constant
    !        Cpd= air heatcapacity at constant pressure
    !        fsin= net solar radiation
    !        fs= net solar radiation array
    !        wpthep= wet bulb equavelnt heatflux (is sum of sensible and latent surfaceheat flux)
    !        lv = latent heat of vaporization
    !        ix= index for horizontal dimension (always 0 here)
    !        iz= index for vertical height
    !        time_current = time
    !        time_start= starttime
    !        dt= model timestep (should be smaller than 60 sec)
    !        c1= coefficient in force restore 1
    !        c2= coefficient in force restore 2
    !        c3= coefficient in force restore 3
    !        cs= soil heat capacity
    !        rs= soil density
    !        ks= soil heat diffusivity
    !        pi = 3.1415927.....
    !        itf= current time step
    !        itl= next time step
    !        hvel= vector wind speed at first model level
    !        pes0te,t0tet,atet,btet= coefficiensfc_soil_t_in_k in tetens formula
    !        g = gravity acceleration
    !        rd = gas constant dry air
    !        rv = gas constant water vapor
    !        cpl = heat capacity liquid water
    !        qwa= intermediate variable for water vapor
    !        qw = water vapor array
    !        qvs= saturated water vapor pressure at surface
    !        deqw = difference between first model level and surface
    !        fact= resistance factor between first level and surface humidity (see duynkerke 1991)
    !        rcan= canopy resistance (rcan about 100 sm-1)
    !        kgz= surface layer exchange coefficient.


    !
    !       soil parameters
    !
    cs=2.00e3
    rs=1.00e3
    ks=2.00e-7
    d1=sqrt(ks*3600.e0*24.e0)
    c1=2.e0*sqrt(pi)/(rs*cs*d1)
    c2=2.e0*pi/(3600.e0*24.e0)
    c3=sqrt(pi*2.e0)/(exp(pi/4.e0)*rs*cs*sqrt(ks*3600.e0*24.e0* &
                     365.e0))

    ! Initialize surface vegetation temp (veg_t_in_k), soil surf.
    ! temp (sfc_soil_t_in_k), deep soil temp (deep_soil_t_in_k)

    if (time_current .eq. time_start)  then
      do it=1,leaf
        veg_t_in_k(it) = initial_veg_t_in_k
        sfc_soil_t_in_k(it) = initial_sfc_soil_t_in_k
        deep_soil_t_in_k(it) = initial_deep_soil_t_in_k
        shf=0.
      end do
    else
      Frad_LW_up_sfc = stefan_boltzmann * (veg_t_in_k(itf)**4)

      ! Calculate net radiation minus turbulent heat flux
      shfs = Frad_LW_down_sfc - Frad_LW_up_sfc - wpthep * rho_sfc * Cp + fs

      ! Calculate soil heat flux
      shf = 10.0 * ( veg_t_in_k(itf) - sfc_soil_t_in_k(itf) ) + 0.05 * fsin

      ! Update surf veg temp
      veg_t_in_k(itl) = veg_t_in_k(itf) + dt * 5.e-5 * ( shfs - shf )

      ! Update soil temp
      sfc_soil_t_in_k(itl) = sfc_soil_t_in_k(itf) & 
        + dt * ( c1 * shf - c2 * ( sfc_soil_t_in_k(itf)-deep_soil_t_in_k(itf) ) )

      ! Update deep soil temp
      deep_soil_t_in_k(itl) = deep_soil_t_in_k(itf) + dt * c3 * shf
    endif

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

    initial_veg_t_in_k = initial_veg_t_in_k_in
    initial_sfc_soil_t_in_k = initial_sfc_soil_t_in_k_in
    initial_deep_soil_t_in_k = initial_deep_soil_t_in_k_in

  end subroutine initialize_surface

end module surface
