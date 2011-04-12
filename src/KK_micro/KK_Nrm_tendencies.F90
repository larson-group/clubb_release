! $Id$
!===============================================================================
module KK_Nrm_tendencies

  private

  public :: KK_Nrm_evap, &
            KK_Nrm_auto

  contains

  !=============================================================================
  function KK_Nrm_evap( KK_rrm_evap_tndcy, Nrm, rrm )

    ! Description:
    ! This function calculates the evaporation rate of rain drop concentration.
    ! The KK equation for evaporation rate of rain drop concentration relates
    ! the rate of change of rain concentration to the rate of change of rain
    ! water content according to:
    !
    ! ( Delta N_rV / N_rV )_evap = ( Delta r_r / r_r )_evap^nu;
    !
    ! where N_rV has units of num/m^3, and nu is a tuned parameter.  However, KK
    ! suggest nu = 1 is reasonable choice.  This equation can be rearranged and
    ! written in terms of N_r (in num/kg):
    !
    ! ( Delta N_r / Delta r_r )_evap = N_r / r_r;
    !
    ! and then can be written as:
    !
    ! ( dN_r / dr_r )_evap = N_r / r_r.
    !
    ! The rate of change of N_r with respect to time can be written as:
    !
    ! ( dN_r / dt )_evap = ( dN_r / dr_r )_evap * ( dr_r / dt )_evap;
    !
    ! which becomes:
    !
    ! ( dN_r / dt )_evap = ( N_r / r_r ) * ( dr_r / dt )_evap;
    !
    ! where ( dr_r / dt )_evap is the KK rain water evaporation rate.

    ! References:
    !  Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    !    Parameterization in a Large-Eddy Simulation Model of Marine
    !    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    !  -- Eq. 23.
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: &
      KK_rrm_evap_tndcy, & ! Mean < dr_r/dt > due to evaporation    [(kg/kg)/s]
      Nrm,               & ! Mean rain drop concentration, < N_r >  [num/kg]
      rrm                  ! Mean rain water mixing ratio, < r_r >  [kg/kg]

    ! Return Variables
    real :: &
      KK_Nrm_evap  ! KK rain drop concentration evaporation rate    [(kg/kg)/s]
 

    ! Evaporation of N_r.
    KK_Nrm_evap = ( Nrm / rrm ) * KK_rrm_evap_tndcy


    return

  end function KK_Nrm_evap

  !=============================================================================
  function KK_Nrm_auto( KK_rrm_auto_tndcy )

    ! Description:
    ! This function calculates the production rate of rain drop concentration
    ! from the process of autoconversion of cloud droplets into rain drops.
    ! The KK equation for the source of rain drop concentration due to
    ! autoconversion is:
    !
    ! ( dN_rV / dt )_auto
    !    = ( dr_r / dt )_auto / ( ( 4*pi*rho_lw / (3*rho_a) ) * r_0^3 ).
    !
    ! Since Nr (in num/kg) = N_rV (in num/m^3) / rho_a, the equation becomes:
    !
    ! ( dN_r / dt )_auto
    !    = ( dr_r / dt )_auto / ( ( 4*pi*rho_lw / 3 ) * r_0^3 ).
    !
    ! CLUBB uses N_r specified in terms of num/kg.

    ! References:
    !  Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    !    Parameterization in a Large-Eddy Simulation Model of Marine
    !    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    !  -- Eq. 32.
    !-----------------------------------------------------------------------

    use constants_clubb, only: & 
        rho_lw,  & ! Constant(s)
        pi

    use parameters_microphys, only: &
        r_0 ! Constant(s)

    implicit none

    ! Input Variables
    real, intent(in) :: &
      KK_rrm_auto_tndcy  ! KK rain water autoconversion rate   [(kg/kg)/s]

    ! Return Variable
    real :: &
      KK_Nrm_auto  ! KK rain drop concentration autoconversion rate [(kg/kg)/s]


    ! Production of N_r through autoconversion.
    KK_Nrm_auto = KK_rrm_auto_tndcy / ( (4.0/3.0) * pi * rho_lw * r_0**3.0 )


    return

  end function KK_Nrm_auto

!===============================================================================

end module KK_Nrm_tendencies
