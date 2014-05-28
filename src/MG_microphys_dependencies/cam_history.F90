! $Id$
module cam_history
!
! Dummy module for importing variables into morrison_gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: addfld, add_default, phys_decomp, outfld
  
  integer :: phys_decomp = 0 ! Not used in CLUBB
  
  contains

!================================================================================================
  subroutine addfld (fname, units, numlev, avgflag, long_name, &
                     decomp_type, flag_xyfill, flag_isccplev, &
                     flag_cospprstaulev, flag_cospprstaumodislev, flag_cosphtdbzelev, &
                     flag_cosphtsrlev, flag_cosphtmlscollev, &
                     flag_cosphtmisrtaulev,flag_cospht,flag_cospscol,&
                     flag_cospsza,sampling_seq)
    !
    !  Description: The original subroutine initiates variables to output history
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------

      character(len=*), intent(in) :: fname      ! field name--should be "max_fieldname_len"
                                                 ! characters long or less
      character(len=*), intent(in) :: units      ! units of fname--should be 8 chars
      character(len=1), intent(in) :: avgflag    ! averaging flag
      character(len=*), intent(in) :: long_name  ! long name of field
      
      integer, intent(in) :: numlev              ! number of vertical levels (dimension and loop)
      integer, intent(in) :: decomp_type         ! decomposition type

      logical, intent(in), optional :: flag_xyfill
      logical, intent(in), optional :: flag_isccplev            ! levels are ISCCP levels not vertical
      logical, intent(in), optional :: flag_cospprstaulev       ! COSP prstau output dimension
      logical, intent(in), optional :: flag_cospprstaumodislev
      logical, intent(in), optional :: flag_cosphtdbzelev       ! COSP htdbze levels output dimension
      logical, intent(in), optional :: flag_cosphtsrlev         ! COSP htsr levels output dimension
      logical, intent(in), optional :: flag_cosphtmlscollev
      logical, intent(in), optional :: flag_cosphtmisrtaulev
      logical, intent(in), optional :: flag_cospht              ! COSP ht output dimension
      logical, intent(in), optional :: flag_cospscol            ! COSP scol output dimension
      logical, intent(in), optional :: flag_cospsza             ! COSP sza output dimension

      character(len=*), intent(in), optional :: sampling_seq

  end subroutine addfld
  
!================================================================================================
  subroutine add_default (name, tindex, flag)
    !
    !  Description: The original subroutine initiates variables to output history for default
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------
  
      character(len=*), intent(in) :: name  ! field name
      character(len=1), intent(in) :: flag  ! averaging flag

      integer, intent(in) :: tindex         ! history tape index
      
  end subroutine add_default
  
!================================================================================================
   subroutine outfld (fname, field, idim, c)
    !
    !  Description: The original subroutine is used for saving variables to history
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8 

      character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

      integer, intent(in) :: idim           ! Longitude dimension of field array
      integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

      real(kind=r8), intent(in) :: field(idim,*)     ! Array containing field values
      
   end subroutine outfld
   
!================================================================================================

end module cam_history
