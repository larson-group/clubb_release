! $Id$
module corr_array_fixed_indices

  ! Prescribed correlation array indices corresponding to the *_corr_array_cloud/below.in files
  integer, parameter, public :: &
    iicorr_s_mellor = 1, &
    iicorr_t_mellor = 2, &
    iicorr_w        = 3

  integer, parameter, public :: &
    iicorr_rrain    = 5, &
    iicorr_rice     = 7, &
    iicorr_rsnow    = 9, &
    iicorr_rgraupel = 11

  integer, parameter, public :: &
    iicorr_Ncn      = 4, &
    iicorr_Nr       = 6, &
    iicorr_Ni       = 8, &
    iicorr_Nsnow    = 10, &
    iicorr_Ngraupel = 12

  ! Number of variables in the *_corr_array_cloud/below.in files
  integer, parameter, public :: &
    n_variables = 12

  implicit none

  private

end module corr_array_fixed_indices
