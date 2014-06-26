!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_parameter_module
! Description:
!   This module defines the derived type pdf_parameter.
! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default scope

  public :: pdf_parameter

  type pdf_parameter
    real( kind = core_rknd ) :: &
      w1,          & ! Mean of w (1st PDF component)                       [m/s]
      w2,          & ! Mean of w (2nd PDF component)                       [m/s]
      varnce_w1,   & ! Variance of w (1st PDF component)               [m^2/s^2]
      varnce_w2,   & ! Variance of w (2nd PDF component)               [m^2/s^2]
      rt1,         & ! Mean of r_t (1st PDF component)                   [kg/kg]
      rt2,         & ! Mean of r_t (2nd PDF component)                   [kg/kg]
      varnce_rt1,  & ! Variance of r_t (1st PDF component)           [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t (2nd PDF component)           [kg^2/kg^2]
      thl1,        & ! Mean of th_l (1st PDF component)                      [K]
      thl2,        & ! Mean of th_l (2nd PDF component)                      [K]
      varnce_thl1, & ! Variance of th_l (1st PDF component)                [K^2]
      varnce_thl2, & ! Variance of th_l (2nd PDF component)                [K^2]
      rrtthl,      & ! Correlation between r_t and th_l (both components)    [-]
      alpha_thl,   & ! Factor relating to normalized variance for th_l       [-]
      alpha_rt,    & ! Factor relating to normalized variance for r_t        [-]
      crt1,        & ! Coef. on r_t in s/t eqns. (1st PDF comp.)             [-]
      crt2,        & ! Coef. on r_t in s/t eqns. (2nd PDF comp.)             [-]
      cthl1,       & ! Coef. on th_l in s/t eqns. (1st PDF comp.)    [(kg/kg)/K]
      cthl2,       & ! Coef. on th_l in s/t eqns. (2nd PDF comp.)    [(kg/kg)/K]
      chi_1,          & ! Mean of s (1st PDF component)                     [kg/kg]
      chi_2,          & ! Mean of s (2nd PDF component)                     [kg/kg]
      stdev_chi_1,    & ! Standard deviation of s (1st PDF component)       [kg/kg]
      stdev_chi_2,    & ! Standard deviation of s (2nd PDF component)       [kg/kg]
      stdev_eta_1,    & ! Standard deviation of t (1st PDF component)       [kg/kg]
      stdev_eta_2,    & ! Standard deviation of t (2nd PDF component)       [kg/kg]
      covar_chi_eta_1,  & ! Covariance of s and t (1st PDF component)     [kg^2/kg^2]
      covar_chi_eta_2,  & ! Covariance of s and t (2nd PDF component)     [kg^2/kg^2]
      corr_chi_eta_1,   & ! Correlation between s and t (1st PDF component)       [-]
      corr_chi_eta_2,   & ! Correlation between s and t (2nd PDF component)       [-]
      rsl1,        & ! Mean of r_sl (1st PDF component)                  [kg/kg]
      rsl2,        & ! Mean of r_sl (2nd PDF component)                  [kg/kg]
      rc1,         & ! Mean of r_c (1st PDF component)                   [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)                   [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)                    [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)                    [-]
      mixt_frac      ! Weight of 1st PDF component (Sk_w dependent)          [-]
  end type pdf_parameter

#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */

  public :: pack_pdf_params, unpack_pdf_params

  integer, public, parameter :: num_pdf_params = 36

  !-------
  contains
  !-------

  subroutine pack_pdf_params(pdf_params, nz, r_param_array)
    implicit none
    ! Input a pdf_parameter array with nz instances of pdf_parameter
    integer, intent(in) :: nz ! Num Vert Model Levs
    type (pdf_parameter), dimension(nz), intent(in) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(out) :: &
       r_param_array  

    ! Local Loop vars
    integer :: k, p    

    do k = 1,nz
       do p = 1,num_pdf_params
 
	r_param_array(k,p) = get_param_at_ind(pdf_params(k), p)
       
       end do ! p
    end do ! k

  end subroutine pack_pdf_params

  subroutine unpack_pdf_params(r_param_array, nz, pdf_params)
    implicit none
    ! Input a two dimensional real array with pdf values
    integer, intent(in) :: nz ! Num Vert Model Levs
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
       r_param_array 

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), dimension(nz), intent(out) :: pdf_params

    ! Local Loop vars
    integer :: k, p
    ! temp var
    real (kind = core_rknd) :: value
    
    do k = 1,nz
       do p = 1,num_pdf_params

       	  value = r_param_array(k,p)
          call set_param_at_ind(pdf_params(k), p, value)       

       end do ! p
    end do ! k

  end subroutine unpack_pdf_params

  real function get_param_at_ind(pp_struct, ind)
    implicit none
    type (pdf_parameter), intent(in) :: pp_struct
    integer, intent(in) :: ind

    SELECT CASE (ind)
      CASE (1)
      	   get_param_at_ind = pp_struct%w1
      CASE (2)
      	   get_param_at_ind = pp_struct%w2
      CASE (3)
      	   get_param_at_ind = pp_struct%varnce_w1
      CASE (4)
      	   get_param_at_ind = pp_struct%varnce_w2
      CASE (5)
      	   get_param_at_ind = pp_struct%rt1
      CASE (6)
      	   get_param_at_ind = pp_struct%rt2
      CASE (7)
      	   get_param_at_ind = pp_struct%varnce_rt1
      CASE (8)
      	   get_param_at_ind = pp_struct%varnce_rt2
      CASE (9)
      	   get_param_at_ind = pp_struct%thl1
      CASE (10)
      	   get_param_at_ind = pp_struct%thl2
      CASE (11)
      	   get_param_at_ind = pp_struct%varnce_thl1
      CASE (12)
      	   get_param_at_ind = pp_struct%varnce_thl2
      CASE (13)
      	   get_param_at_ind = pp_struct%rrtthl
      CASE (14)
      	   get_param_at_ind = pp_struct%alpha_thl
      CASE (15)
      	   get_param_at_ind = pp_struct%alpha_rt
      CASE (16)
      	   get_param_at_ind = pp_struct%crt1
      CASE (17)
      	   get_param_at_ind = pp_struct%crt2
      CASE (18)
      	   get_param_at_ind = pp_struct%cthl1
      CASE (19)
      	   get_param_at_ind = pp_struct%cthl2
      CASE (20)
      	   get_param_at_ind = pp_struct%chi_1
      CASE (21)
      	   get_param_at_ind = pp_struct%chi_2
      CASE (22)
      	   get_param_at_ind = pp_struct%stdev_chi_1
      CASE (23)
      	   get_param_at_ind = pp_struct%stdev_chi_2
      CASE (24)
      	   get_param_at_ind = pp_struct%stdev_eta_1
      CASE (25)
      	   get_param_at_ind = pp_struct%stdev_eta_2
      CASE (26)
      	   get_param_at_ind = pp_struct%covar_chi_eta_1
      CASE (27)
      	   get_param_at_ind = pp_struct%covar_chi_eta_2
      CASE (28)
      	   get_param_at_ind = pp_struct%corr_chi_eta_1
      CASE (29)
      	   get_param_at_ind = pp_struct%corr_chi_eta_2
      CASE (30)
      	   get_param_at_ind = pp_struct%rsl1
      CASE (31)
      	   get_param_at_ind = pp_struct%rsl2
      CASE (32)
      	   get_param_at_ind = pp_struct%rc1
      CASE (33)
      	   get_param_at_ind = pp_struct%rc2
      CASE (34)
      	   get_param_at_ind = pp_struct%cloud_frac1
      CASE (35)
      	   get_param_at_ind = pp_struct%cloud_frac2
      CASE (36)
      	   get_param_at_ind = pp_struct%mixt_frac
      CASE DEFAULT
      	   get_param_at_ind = 0.0/0.0 !NaN - watch out!
    END SELECT

    RETURN
  end function get_param_at_ind

  subroutine set_param_at_ind(pp_struct, ind, val)
    implicit none
    type (pdf_parameter), intent(inout) :: pp_struct
    integer, intent(in) :: ind
    real (kind = core_rknd), intent(in) :: val

    SELECT CASE (ind)
      CASE (1)
      	   pp_struct%w1 = val
      CASE (2)
      	   pp_struct%w2 = val
      CASE (3)
      	   pp_struct%varnce_w1 = val
      CASE (4)
      	   pp_struct%varnce_w2 = val
      CASE (5)
      	   pp_struct%rt1 = val
      CASE (6)
      	   pp_struct%rt2 = val
      CASE (7)
      	   pp_struct%varnce_rt1 = val
      CASE (8)
      	   pp_struct%varnce_rt2 = val
      CASE (9)
      	   pp_struct%thl1 = val
      CASE (10)
      	   pp_struct%thl2 = val
      CASE (11)
      	   pp_struct%varnce_thl1 = val
      CASE (12)
      	   pp_struct%varnce_thl2 = val
      CASE (13)
      	   pp_struct%rrtthl = val
      CASE (14)
      	   pp_struct%alpha_thl = val
      CASE (15)
      	   pp_struct%alpha_rt = val
      CASE (16)
      	   pp_struct%crt1 = val
      CASE (17)
      	   pp_struct%crt2 = val
      CASE (18)
      	   pp_struct%cthl1 = val
      CASE (19)
      	   pp_struct%cthl2 = val
      CASE (20)
      	   pp_struct%chi_1 = val
      CASE (21)
      	   pp_struct%chi_2 = val
      CASE (22)
      	   pp_struct%stdev_chi_1 = val
      CASE (23)
      	   pp_struct%stdev_chi_2 = val
      CASE (24)
      	   pp_struct%stdev_eta_1 = val
      CASE (25)
      	   pp_struct%stdev_eta_2 = val
      CASE (26)
      	   pp_struct%covar_chi_eta_1 = val
      CASE (27)
      	   pp_struct%covar_chi_eta_2 = val
      CASE (28)
      	   pp_struct%corr_chi_eta_1 = val
      CASE (29)
      	   pp_struct%corr_chi_eta_2 = val
      CASE (30)
      	   pp_struct%rsl1 = val
      CASE (31)
      	   pp_struct%rsl2 = val
      CASE (32)
      	   pp_struct%rc1 = val
      CASE (33)
      	   pp_struct%rc2 = val
      CASE (34)
      	   pp_struct%cloud_frac1 = val
      CASE (35)
      	   pp_struct%cloud_frac2 = val
      CASE (36)
      	   pp_struct%mixt_frac = val
      CASE DEFAULT
      	   ! do nothing !
    END SELECT

  end subroutine set_param_at_ind

#endif

end module pdf_parameter_module
