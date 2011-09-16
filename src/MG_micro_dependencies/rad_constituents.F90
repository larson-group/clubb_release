! $Id$
module rad_constituents
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: rad_cnst_get_info, rad_cnst_get_aer_props
  
  contains
  
  subroutine rad_cnst_get_info( naero, aernames, aersources, aerindices, &
                                   ngas,  gasnames, gassources, gasindices, &
                                   use_data_o3, diagnosticindex )
     !
     !  Description: The original subroutine returns info about aerosol climate list.
     !               In our case, we don't directly use this, so this
     !               subroutine does nothing.
     !
     !---------------------------------------------------------------------------------
     
     implicit none

     ! Arguments
     integer,           optional, intent(out) :: naero
     integer,           optional, intent(out) :: ngas
     character(len=64), optional, intent(out) :: aernames(:)
     character(len=64), optional, intent(out) :: gasnames(:)
     character(len=1),  optional, intent(out) :: aersources(:)
     character(len=1),  optional, intent(out) :: gassources(:)
     integer,           optional, intent(out) :: aerindices(:)
     integer,           optional, intent(out) :: gasindices(:)
     logical,           optional, intent(out) :: use_data_o3
     integer,           optional, intent(in)  :: diagnosticindex
     
     ! Dummy values
     if (present(naero)) then
       naero = 1
     end if
     
     if (present(ngas)) then
       ngas = 1
     end if
     
     if (present(aernames)) then
       aernames(:) = ' '
     end if
     
     if (present(gasnames)) then
       gasnames(:) = ' '
     end if
     
     if (present(aersources)) then
       aersources(:) = ' '
     end if
     
     if (present(gassources)) then
       gassources(:) = ' '
     end if
     
     if (present(aerindices)) then
       aerindices(:) = 1
     end if
     
     if (present(gasindices)) then
       gasindices(:) = 1
     end if
     
     if (present(use_data_o3)) then
       use_data_o3 = .false.
     end if
     
     return
     
  end subroutine rad_cnst_get_info
  
  subroutine rad_cnst_get_aer_props( &
     list_idx, diagnosticindex, &
     sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
     sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
     sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
     refindex_real_aer_sw, refindex_im_aer_sw, refindex_real_aer_lw, refindex_im_aer_lw, &
     refindex_real_water_sw, refindex_im_water_sw, refindex_real_water_lw, refindex_im_water_lw, &
     r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
     aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, num_to_mass_aer)
    !
    !  Description: The original subroutine returns requested properties for specified
    !               climate aerosol. In our case, we don't directly use this, so this
    !               subroutine does nothing.
    !
    !---------------------------------------------------------------------------------
   
    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none
     
     integer,                     intent(in)  :: list_idx
     integer,           optional, intent(in)  :: diagnosticindex ! index to rad diagnostic call
     real(r8),          optional, pointer     :: sw_hygro_ext(:,:) 
     real(r8),          optional, pointer     :: sw_hygro_ssa(:,:) 
     real(r8),          optional, pointer     :: sw_hygro_asm(:,:) 
     real(r8),          optional, pointer     :: lw_hygro_ext(:,:)         
     real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
     real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
     real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
     real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
     real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
     real(r8),          optional, pointer     :: lw_ext(:)         
     real(r8),          optional, pointer     :: refindex_real_aer_sw(:)
     real(r8),          optional, pointer     :: refindex_im_aer_sw(:)
     real(r8),          optional, pointer     :: refindex_real_aer_lw(:)
     real(r8),          optional, pointer     :: refindex_im_aer_lw(:)
     real(r8),          optional, pointer     :: refindex_real_water_sw(:)
     real(r8),          optional, pointer     :: refindex_im_water_sw(:)
     real(r8),          optional, pointer     :: refindex_real_water_lw(:)
     real(r8),          optional, pointer     :: refindex_im_water_lw(:)
     character(len=20), optional, intent(out) :: aername           
     real(r8),          optional, intent(out) :: density_aer
     real(r8),          optional, intent(out) :: hygro_aer
     real(r8),          optional, intent(out) :: dryrad_aer        
     real(r8),          optional, intent(out) :: dispersion_aer    
     real(r8),          optional, intent(out) :: num_to_mass_aer   

     real(r8),          optional, pointer     :: r_sw_ext(:,:)         
     real(r8),          optional, pointer     :: r_sw_scat(:,:)         
     real(r8),          optional, pointer     :: r_sw_ascat(:,:)         
     real(r8),          optional, pointer     :: r_lw_abs(:,:)         
     real(r8),          optional, pointer     :: mu(:)
     
     ! Dummy values
     if (present(aername)) then
       aername = 'a'
     end if
     
     if (present(density_aer)) then
       density_aer = 1._r8
     end if
     
     if (present(hygro_aer)) then
       hygro_aer = 1._r8
     end if
     
     if (present(dryrad_aer)) then
       dryrad_aer = 1._r8
     end if
     
     if (present(dispersion_aer)) then
       dispersion_aer = 10._r8
     end if
     
     if (present(num_to_mass_aer)) then
       num_to_mass_aer = 1._r8
     end if
     
     return
     
  end subroutine rad_cnst_get_aer_props

end module rad_constituents
