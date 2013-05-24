! $Id: pdf_parameter_module.F90 5668 2012-01-29 03:40:28Z bmg2@uwm.edu $
module hydromet_pdf_parameter_module
! Description:
!   This module defines the derived type pdf_parameter.
! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default scope

  public :: hydromet_pdf_parameter

  type hydromet_pdf_parameter
    real( kind = core_rknd ) :: &
      mu_w_1,      & ! Mean of w (1st PDF component)                       [m/s]
      mu_w_2,      & ! Mean of w (2nd PDF component)                       [m/s]
      mu_s_1,      & ! Mean of s (1st PDF component)                     [kg/kg]
      mu_s_2,      & ! Mean of s (2nd PDF component)                     [kg/kg]
      mu_t_1,      & ! Mean of t (1st PDF component)                     [kg/kg]
      mu_t_2,      & ! Mean of t (2nd PDF component)                     [kg/kg]
      mu_rr_1,     & ! Mean of rr (1st PDF component) in-precip (ip)     [kg/kg]
      mu_rr_2,     & ! Mean of rr (2nd PDF component) ip                 [kg/kg]
      mu_Nr_1,     & ! Mean of Nr (1st PDF component) ip                [num/kg]
      mu_Nr_2,     & ! Mean of Nr (2nd PDF component) ip                [num/kg]
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_w_1,   & ! Standard deviation of w (1st PDF component)         [m/s]
      sigma_w_2,   & ! Standard deviation of w (2nd PDF component)         [m/s]
      sigma_s_1,   & ! Standard deviation of s (1st PDF component)       [kg/kg]
      sigma_s_2,   & ! Standard deviation of s (2nd PDF component)       [kg/kg]
      sigma_t_1,   & ! Standard deviation of t (1st PDF component)       [kg/kg]
      sigma_t_2,   & ! Standard deviation of t (2nd PDF component)       [kg/kg]
      sigma_rr_1,  & ! Standard deviation of rr (1st PDF component) ip   [kg/kg]
      sigma_rr_2,  & ! Standard deviation of rr (2nd PDF component) ip   [kg/kg]
      sigma_Nr_1,  & ! Standard deviation of Nr (1st PDF component) ip  [num/kg]
      sigma_Nr_2,  & ! Standard deviation of Nr (2nd PDF component) ip  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2    ! Standard deviation of Ncn (2nd PDF component)    [num/kg]
    end type hydromet_pdf_parameter

end module hydromet_pdf_parameter_module
