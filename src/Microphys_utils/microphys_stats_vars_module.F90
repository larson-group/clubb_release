! $Id$
module microphys_stats_vars_module

! Description:
!   This module contains the derived type microphys_stats_vars, which is used to
!   feed output variables out of microphysics schemes

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Set Default Scope

  public :: microphys_stats_vars

  type microphys_stats_vars

    ! Number of output variables from microphysics
    integer :: num_vars

    ! An array of statistics indices corresponding to the output variables
    integer, dimension(:), pointer :: &
      stats_indices

    ! Values of the output variables
    real( kind = core_rknd ), pointer :: &
      output_values

  end type microphys_stats_vars

end module microphys_stats_vars_module
