!-------------------------------------------------------------------------
! $Id$ 
!===============================================================================
module T_in_K_module

  implicit none

  private ! Default scope

  public :: thlm2T_in_K, T_in_K2thlm

  interface thlm2T_in_K
    module procedure thlm2T_in_K_k    ! Works over a single vertical level
    module procedure thlm2T_in_K_1D   ! Works over all vertical levels 
    module procedure thlm2T_in_K_2D   ! Works over all vertical levels and columns
  end interface thlm2T_in_K

  contains

  !-------------------------------------------------------------------------------
  ! Wrapped in interface thlm2T_in_K
  function thlm2T_in_K_k( thlm, exner, rcm )  & 
                  result( T_in_K )

  ! Description:
  !   Calculates absolute temperature from liquid water potential
  !   temperature.  (Does not include ice.)

  ! References: 
  !   Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51). 
  !-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
        Cp,  & ! Dry air specific heat at constant p [J/kg/K]
        Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! ---------------------- Input ----------------------
    real( kind = core_rknd ), intent(in) :: & 
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: & 
      T_in_K ! Result temperature [K]

    ! ---------------------- Begin Code ----------------------
    T_in_K = thlm * exner + Lv * rcm / Cp

    return
  end function thlm2T_in_K_k

  !-------------------------------------------------------------------------------
  ! Wrapped in interface thlm2T_in_K
  function thlm2T_in_K_1D( nz, thlm, exner, rcm )  & 
                   result( T_in_K )

  ! Description:
  !   Calculates absolute temperature from liquid water potential
  !   temperature.  (Does not include ice.)

  ! References: 
  !   Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51). 
  !-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
        Cp,  & ! Dry air specific heat at constant p [J/kg/K]
        Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! ---------------------- Input ----------------------
    integer, intent(in) :: &
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) :: & 
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(nz) :: & 
      T_in_K ! Result temperature [K]

    integer :: k

    ! ---------------------- Begin Code ----------------------
    do k = 1, nz
      T_in_K(k) = thlm(k) * exner(k) + Lv * rcm(k) / Cp
    end do

    return
  end function thlm2T_in_K_1D

  !-------------------------------------------------------------------------------
  ! Wrapped in interface thlm2T_in_K
  function thlm2T_in_K_2D( nz, ngrdcol, thlm, exner, rcm )  & 
                   result( T_in_K )

  ! Description:
  !   Calculates absolute temperature from liquid water potential
  !   temperature.  (Does not include ice.)

  ! References: 
  !   Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51). 
  !-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
        Cp,  & ! Dry air specific heat at constant p [J/kg/K]
        Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input 
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      T_in_K ! Result temperature [K]

    integer :: i, k

    ! ---- Begin Code ----

    !$acc data copyin( thlm, exner, rcm ) &
    !$acc     copyout( T_in_K )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        T_in_K(i,k) = thlm(i,k) * exner(i,k) + Lv * rcm(i,k) / Cp
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return
  end function thlm2T_in_K_2D

!-------------------------------------------------------------------------------
  elemental function T_in_K2thlm( T_in_K, exner, rcm )  & 
    result( thlm )

! Description:
!   Calculates liquid water potential temperature from absolute temperature 

! References: 
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
        Cp,  & ! Dry air specific heat at constant p [J/kg/K]
        Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input 
    real( kind = core_rknd ), intent(in) :: & 
      T_in_K, &! Result temperature [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: & 
      thlm    ! Liquid potential temperature  [K]

    ! ---- Begin Code ----

    thlm = ( T_in_K - Lv/Cp * rcm ) / exner 

    return
  end function T_in_K2thlm
!-------------------------------------------------------------------------------

end module T_in_K_module
