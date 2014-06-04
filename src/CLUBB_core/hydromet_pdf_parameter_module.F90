!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module hydromet_pdf_parameter_module

  ! Description:
  ! This module defines the derived type hydromet_pdf_parameter.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd  ! Variable(s)

  implicit none

  private ! Default scope

  public :: hydromet_pdf_parameter,   & ! Variable type
            init_hydromet_pdf_params    ! Procedure

  integer, parameter, private :: &
    max_hydromet_dim = 8

  type hydromet_pdf_parameter

    real( kind = core_rknd ), dimension(max_hydromet_dim) :: &
      hm1,        & ! Mean of hydrometeor, hm (1st PDF component)   [units vary]
      hm2,        & ! Mean of hydrometeor, hm (2nd PDF component)   [units vary]
      mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip) [units vary]
      mu_hm_2,    & ! Mean of hm (2nd PDF component) ip             [units vary]
      sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip [un. vary]
      sigma_hm_2    ! Standard deviation of hm (2nd PDF component) ip [un. vary]

    real( kind = core_rknd ) :: &
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2    ! Standard deviation of Ncn (2nd PDF component)    [num/kg]

    real( kind = core_rknd ) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

  end type hydromet_pdf_parameter

contains

  !=============================================================================
  subroutine init_hydromet_pdf_params( hydromet_pdf_params )

    ! Description:
    ! Initialize the elements of hydromet_pdf_params.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr  ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    implicit none

    ! Output Variable
    type(hydromet_pdf_parameter), dimension(gr%nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]

    ! Local Variable
    integer :: k  ! Loop index


    ! Initialize hydromet_pdf_params.
    do k = 1, gr%nz, 1

       hydromet_pdf_params(k)%hm1 = zero
       hydromet_pdf_params(k)%hm2 = zero
       hydromet_pdf_params(k)%mu_hm_1 = zero
       hydromet_pdf_params(k)%mu_hm_2 = zero
       hydromet_pdf_params(k)%sigma_hm_1 = zero
       hydromet_pdf_params(k)%sigma_hm_2 = zero

       hydromet_pdf_params(k)%mu_Ncn_1 = zero
       hydromet_pdf_params(k)%mu_Ncn_2 = zero
       hydromet_pdf_params(k)%sigma_Ncn_1 = zero
       hydromet_pdf_params(k)%sigma_Ncn_2 = zero

       hydromet_pdf_params(k)%precip_frac = zero
       hydromet_pdf_params(k)%precip_frac_1 = zero
       hydromet_pdf_params(k)%precip_frac_2 = zero

    enddo ! k = 1, gr%nz, 1


    return

  end subroutine init_hydromet_pdf_params

!===============================================================================

end module hydromet_pdf_parameter_module
