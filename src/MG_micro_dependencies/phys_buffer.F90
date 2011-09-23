module phys_buffer

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Implement a physics buffer to hold arrays that must persist
!   across timesteps or between calls to different physics packages.
!
!   Current implementation only supports one buffer.
!
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols
   use abortutils,   only: endrun

   implicit none
   private
   save

! Public methods 

   public ::&
      pbuf_init,          &! initialize physics buffer
      pbuf_add,           &! add field to physics buffer
      pbuf_get_fld_idx,   &! get index of specified field in the physics buffer
      pbuf_old_tim_idx,   &! return the index for the oldest time
      pbuf_next_tim_idx,  &! return the index for the next time
      pbuf_allocate,      &! allocate memory for physics buffer fields
      pbuf_deallocate,    &! deallocate memory for physics buffer fields
      pbuf_setval          ! set value for a field in the buffer

! Public types and data

   type, public :: pbuf_fld
      character(len=16)                       :: name
      character(len=16)                       :: scope
      integer                                 :: fdim, mdim, ldim
      real(r8), pointer, dimension(:,:,:,:,:) :: fld_ptr
   end type pbuf_fld


   integer, public, parameter :: pbuf_size_max=1000
   type(pbuf_fld), public, dimension(pbuf_size_max) :: pbuf


! Private module data

   integer :: pbuf_size = 0
   integer :: old_time_idx = 1
!=========================================================================================
contains
!=========================================================================================
subroutine pbuf_init()

! Initialize physics buffer.

   implicit none
!-----------------------------------------------------------------------------------------
   integer :: i

   do i=1,pbuf_size_max
      nullify(pbuf(i)%fld_ptr)
   end do


end subroutine pbuf_init

!=========================================================================================
subroutine pbuf_add(name, fdim, mdim, ldim)

! Add a field to the physics buffer

   implicit none

   character(len=*), intent(in)  :: name   ! field name (case insensitive)
   integer,          intent(in)  :: fdim   ! first generic field dimension
   integer,          intent(in)  :: mdim   ! middle generic field dimension
   integer,          intent(in)  :: ldim   ! last generic field dimension

! Local variables
   character(len=*), parameter :: sub = 'pbuf_add'
   integer :: i
   character(len=len(name)) :: uname
   integer :: index
!-----------------------------------------------------------------------------------------

   if ( pbuf_size >= pbuf_size_max ) then
      call endrun (sub//': max number physics buffer fields exceeded. Increase pbuf_size_max in phys_buffer.F90')
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%name == name ) then
         call endrun (sub//': ERROR: field name '//name//' is already in use.')
      end if
   end do

   pbuf_size = pbuf_size + 1
   index = pbuf_size
   pbuf(index)%name = name
   pbuf(index)%fdim = fdim
   pbuf(index)%mdim = mdim
   pbuf(index)%ldim = ldim

end subroutine pbuf_add

!=========================================================================================
function pbuf_get_fld_idx(name, failcode)

! Get index of specified field in the physics buffer.  String matching is case insensitive.
! Call endrun if name not found

   implicit none

   character(len=*), intent(in)  :: name   ! field name 
   integer, intent(in), optional :: failcode

! Return value
   integer :: pbuf_get_fld_idx

! Local variables
   integer :: i
   character(len=len(name)) :: Uname
!-----------------------------------------------------------------------------------------

!
!  Search for specified field in physics buffer, assuming that case of
!  argument "name" matches definition in pbuf structure.
!

   do i = 1, pbuf_size
      if ( pbuf(i)%name == name ) then
         pbuf_get_fld_idx = i
         return
      end if
   end do

   if ( present(failcode)) then
      pbuf_get_fld_idx = failcode
      return
   else
      call endrun ('PBUF_GET_FLD_IDX: index not found for '//name)
   endif

end function pbuf_get_fld_idx

!=========================================================================================
function pbuf_old_tim_idx()

! Return index of oldest time sample in the physics buffer.

   implicit none

! Return value
   integer :: pbuf_old_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_old_tim_idx = 1

end function pbuf_old_tim_idx

!=========================================================================================
function pbuf_next_tim_idx(idx)

! Return index of next time sample in the physics buffer.

   implicit none

   integer, intent(in) :: idx

! Return value
   integer :: pbuf_next_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_next_tim_idx = 1

end function pbuf_next_tim_idx

!=========================================================================================
subroutine pbuf_allocate()

! N.B. This routine must be called after phys_grid_init because that's
!      where begchunk and endchunk are set

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'pbuf_allocate'
   integer :: i, fdim, mdim, ldim, istat
   logical :: allocate_all
!-----------------------------------------------------------------------------------------

   do i = 1, pbuf_size
     fdim = pbuf(i)%fdim
     mdim = pbuf(i)%mdim
     ldim = pbuf(i)%ldim

     allocate(pbuf(i)%fld_ptr(fdim,pcols,mdim,1,ldim), stat=istat)
     if ( istat /= 0 ) then
        call endrun (sub//': ERROR: allocate failed for '//pbuf(i)%name)
     end if
     pbuf(i)%fld_ptr = -9999.99
         
   end do

end subroutine pbuf_allocate

!=========================================================================================
subroutine pbuf_deallocate()

! Deallocate storage for fields in the physics buffer with the specified scope.
! If global_allocate_all=.true. then storage for both global and physpkg scope 
! is deallocated just once, when scope='global'.

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'pbuf_deallocate'
   integer :: i, fdim, mdim, ldim
   logical :: deallocate_all
!-----------------------------------------------------------------------------------------

   do i = 1, pbuf_size
     if (associated(pbuf(i)%fld_ptr)) then
        deallocate(pbuf(i)%fld_ptr)
     else
        call endrun (sub//': ERROR: '//pbuf(i)%name//' is not allocated')
     end if
   end do

end subroutine pbuf_deallocate

!=========================================================================================
subroutine pbuf_setval(name, value)

! Set a value for a field in the physics buffer.

   implicit none

   character(len=*), intent(in)  :: name   ! field name 
   real(r8), dimension(pbuf(1)%mdim), intent(in)  :: value

! Local variables
   character(len=*), parameter :: sub = 'pbuf_setval'
   integer :: idx, lchnk, ncols
   integer :: failcode = -1
!-----------------------------------------------------------------------------------------

   ! get field index
   idx = pbuf_get_fld_idx(name, failcode)
   ! check return value
   if (idx == failcode) then
      call endrun (sub//': ERROR: name not found:'//name)
   end if

   ! check that field pointer is associated -- if so then set the value
   if (associated(pbuf(idx)%fld_ptr)) then
      pbuf(idx)%fld_ptr(1,1,:,1,1) = value
   else
      call endrun (sub//': ERROR: field '//name//' is not allocated')
   end if

end subroutine pbuf_setval

end module phys_buffer
