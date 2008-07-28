!-----------------------------------------------------------------------
! $Id: stats_type.F90,v 1.3 2008-07-28 19:16:20 faschinj Exp $
module stats_type
#ifdef STATS 
!     Description:
!     Contains derived data type 'stats'.
!     Used for storing output statistics to disk.
!-----------------------------------------------------------------------

  use outputfile_class, only: & 
      outputfile ! Type

  use stats_precision, only: & 
      stat_rknd,  & ! Variable(s)
      stat_nknd
  
  implicit none

  private ! Set Default Scope

  public :: stats, stat_assign, stat_update_var, & 
            stat_update_var_pt, stat_begin_update, & 
            stat_end_update, stat_begin_update_pt, & 
            stat_end_update_pt, stat_modify, stat_modify_pt


  ! Derived data types to store GrADS/netCDF statistics
  type stats

    ! Number of fields to sample
    integer :: nn

    ! Vertical extent of variable
    integer :: kk

    ! Vertical levels
    real, pointer, dimension(:) :: z

!   Array to store sampled fields

    real(kind=stat_rknd), pointer, dimension(:,:) :: x

    integer(kind=stat_nknd), pointer, dimension(:,:) :: n
    
    ! Tracks if a field is in the process of an update
    logical, pointer, dimension(:,:) :: in_update

!   Data for GrADS output

    type (outputfile) f

  end type stats


  contains


!-----------------------------------------------------------------------
   subroutine stat_assign( var_index, var_name,  & 
                     var_description, var_units, grid_kind )
!  Description: Assigns pointers for statistics variables in grid
!    
!-----------------------------------------------------------------------
   implicit none
     
   ! Input Variables
       
   integer,intent(in) :: var_index                   ! Variable index       [#]
   character(len = *), intent(in) :: var_name        ! Variable name        []
   character(len = *), intent(in) :: var_description ! Variable description []
   character(len = *), intent(in) :: var_units       ! Variable units       []

   ! Output Variable

   ! Which grid the variable is located on (zt, zm, or sfc )
   type(stats), intent(out) :: grid_kind

   grid_kind%f%var(var_index)%ptr => grid_kind%x(:,var_index)
   grid_kind%f%var(var_index)%name = var_name
   grid_kind%f%var(var_index)%description = var_description
   grid_kind%f%var(var_index)%units = var_units

       !Example of the old format 
       !changed by Joshua Fasching 23 August 2007
       
       !zt%f%var(ithlm)%ptr => zt%x(:,k)
       !zt%f%var(ithlm)%name = "thlm"
       !zt%f%var(ithlm)%description = "thetal (K)"
       !zt%f%var(ithlm)%units = "K"

   return

   end subroutine stat_assign

!---------------------------------------------------------------------
   subroutine stat_update_var( var_index, value, grid_kind )
!
!  Description: This updates the value of a statistics variable 
!         located at var_index associated with grid type 'grid_kind'
!         at a specific level.
!
!---------------------------------------------------------------------            
   use grid_class, only: gr ! Variable(s)
   
   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index ! The index at which the variable is stored  []
            
   real, dimension(gr%nnzp), intent(in) :: & 
   value ! Value of field being added to the statistic    [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc )

   if ( var_index > 0 ) then
      grid_kind%x(:,var_index) =  & 
           grid_kind%x(:,var_index) + value
      grid_kind%n(:,var_index) =  & 
           grid_kind%n(:,var_index) + 1
   end if

   end subroutine stat_update_var
   
!---------------------------------------------------------------------
   subroutine stat_update_var_pt( var_index, grid_level, value, grid_kind )
!
!  Description: This updates the value of a statistics variable 
!     located at var_index associated with grid type 'grid_kind' 
!      at a specific grid_level.
!
!---------------------------------------------------------------------            
  
   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index,       & ! The index at which the variable is stored           []
   grid_level      ! The level at which the variable is to be modified   []
   
   real, intent(in) :: & 
   value ! Value of field being added to the statistic         [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).
  
   if ( var_index > 0 ) then
      grid_kind%x( grid_level, var_index ) =  & 
        grid_kind%x( grid_level, var_index ) + value

      grid_kind%n( grid_level, var_index ) =  & 
         grid_kind%n( grid_level, var_index ) + 1
   end if

   end subroutine stat_update_var_pt
       
!---------------------------------------------------------------------
   subroutine stat_begin_update( var_index, value, grid_kind )
!
!
!        Description: This begins an update of the value of a 
!        statistics variable located at var_index on the (zt, zm, or sfc) grid.
!        Commonly this is used for beginning a budget term calculation. 
!           In this type of stats calculation, we first subtract the field
!        (e.g. rtm) from the statistic, then update rtm by a term (e.g. clip rtm), 
!        and then re-add rtm to the statistic. 
!            
!        Example:
!          call stat_begin_update( irtm_bt, real(rtm / dt), zt ) 
!
!          !!! Perform clipping of rtm !!!
!          
!          call stat_end_update( irtm_bt, real(rtm / dt), zt )
!---------------------------------------------------------------------                    

   use grid_class, only: gr  ! Variable(s)        

   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index      ! The index at which the variable is stored           []
   
   real, dimension(gr%nnzp), intent(in) :: & 
   value          ! Value of field being added to the statistic         [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).

   integer :: i

   do i = 1, gr%nnzp 

     call stat_begin_update_pt & 
             ( var_index, i, value(i), grid_kind )

   enddo

   end subroutine stat_begin_update



!---------------------------------------------------------------------
   subroutine stat_begin_update_pt & 
             ( var_index, grid_level, value, grid_kind )
!
!        Description: This begins an update of the value of a 
!         statistics variable located at var_index associated with the grid type (grid_kind) 
!         at a specific grid_level.
!
!         Commonly this is used for beginning a budget. See stat_begin_update 
!---------------------------------------------------------------------            
 
   use error_code, only: clubb_debug ! Procedure(s)

   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index,    & ! The index at which the variable is stored           []
   grid_level      ! The level at which the variable is to be modified   []
   
   real, intent(in) :: & 
   value ! Value of field being added to the statistic                [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).

   ! Local Variable
   character(len=3) :: char_index
  
   if ( var_index > 0 ) then  ! Are we storing this variable?

       if ( .not. grid_kind%in_update( grid_level, var_index ) ) & 
         then ! Can we begin an update?  

         grid_kind%x(grid_level, var_index) =  & 
                 grid_kind%x(grid_level, var_index) - value

         grid_kind%in_update(grid_level, var_index) = .true.  ! Start Record

       else
        write(char_index(1:3),'(i3.3)') var_index
        call clubb_debug( 1, & 
        "Beginning an update before finishing previous.  "// & 
        "Var index = "// char_index   )

       end if

   end if

   end subroutine stat_begin_update_pt

!---------------------------------------------------------------------
   subroutine stat_end_update( var_index, value, grid_kind )
!
!        Description: This ends the update of a statistics variable 
!         located at var_index associated with the grid type (grid_kind). 
!
!         Commoly used to end the storage of a budget term. See stat_begin_update
!
!---------------------------------------------------------------------            

   use grid_class, only: gr ! Variable(s)

   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index ! The index at which the variable is stored           []
           
   real, dimension(gr%nnzp), intent(in) :: & 
   value ! Value of field being added to the statistic             [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).

   integer :: i 
           
   do i = 1,gr%nnzp             
        call stat_end_update_pt & 
                 ( var_index, i, value(i), grid_kind )
   end do


   end subroutine stat_end_update

!---------------------------------------------------------------------
   subroutine stat_end_update_pt & 
                ( var_index, grid_level, value, grid_kind )
!
!        Description: This ends an update of the value of a 
!         statistics variable located at var_index associated with the grid type (grid_kind) 
!         at a specific grid_level.
!
!         Commonly this is used for beginning a budget. See stat_begin_update 
!---------------------------------------------------------------------            

   use error_code, only: clubb_debug ! Procedure(s)

   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index,   & ! The index at which the variable is stored           []
   grid_level   ! The level at which the variable is to be modified   []     
   
   real, intent(in) :: & 
   value       ! Value of field being added to the statistic         [Units Vary]

   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).

   if ( var_index > 0 ) then ! Are we storing this variable?
           
       if ( grid_kind%in_update(grid_level, var_index) ) then ! Can we end an update? 
               
          call stat_update_var_pt & 
                   ( var_index, grid_level, value, grid_kind ) 
        
          grid_kind%in_update( grid_level, var_index ) = .false. ! End Record

        else

          call clubb_debug( 1, "Ending before beginning update. For variable "// &
          grid_kind%f%var(var_index)%name )

        endif

   endif
    
   end subroutine stat_end_update_pt
   
!---------------------------------------------------------------------
   subroutine stat_modify( var_index, value, grid_kind )
!
!        Description: This modifies the value of a statistics variable 
!        located at var_index in the grid. It does not increment the
!        sampling count.
!
!---------------------------------------------------------------------            
            
   use grid_class, only: gr ! Variable(s)

   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index ! The index at which the variable is stored           []
   
   real, dimension(gr%nnzp), intent(in) :: & 
   value     ! Value of field being added to the statistic         [Units Vary]


   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).

   integer :: i

   do i = 1, gr%nnzp 

    call stat_modify_pt( var_index, i, value(i), grid_kind )

   enddo

   end subroutine stat_modify

!---------------------------------------------------------------------
   subroutine stat_modify_pt & 
             ( var_index, grid_level, value, grid_kind )
!
!        Description: This modifies the value of a statistics variable 
!        located at var_index in the grid at a specific point. It does
!        not increment the sampling count.
!
!---------------------------------------------------------------------            
 
   implicit none

   ! Input Variables(s)

   integer, intent(in) ::  & 
   var_index ! The index at which the variable is stored            []
   
   
   real, intent(in) :: & 
   value      ! Value of field being added to the statistic         [Units Vary]

   integer, intent(in) ::  & 
   grid_level ! The level at which the variable is to be modified   []
   
   ! Input/Output Variable(s)
   type(stats), intent(inout) ::  & 
   grid_kind ! Which grid the variable is located on (zt, zm, or sfc).
  
   if ( var_index > 0 ) then 

     grid_kind%x( grid_level, var_index )  & 
         = grid_kind%x( grid_level, var_index ) + value

   end if
  
   end subroutine stat_modify_pt
   
#endif /*STATS*/
end module stats_type
