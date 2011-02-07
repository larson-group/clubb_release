! $Id$
module rad_constituents
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: rad_cnst_get_clim_info, rad_cnst_get_clim_aer_props
  
  ! These variables are not used anywhere in MG, they are just imported. Because of this we
  ! are setting them to dummy values
  integer :: rad_cnst_get_clim_info = 0, rad_cnst_get_clim_aer_props = 0

end module rad_constituents
