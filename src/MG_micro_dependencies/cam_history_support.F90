! $Id$
module cam_history_support
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  type, public :: field_info
    real(r8) :: fillvalue              ! fillvalue for this variable, set to default if not explicit in addfld

  end type field_info

  real(r8), parameter, public :: fillvalue = 1.e36_r8     ! fill value for netcdf fields

end module cam_history_support

