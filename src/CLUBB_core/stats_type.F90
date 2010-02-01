!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_type
  
  ! Description:
  ! Contains derived data type 'stats'.
  ! Used for storing output statistics to disk.
  !-----------------------------------------------------------------------

  use stat_file_module, only: & 
      stat_file ! Type

  use stats_precision, only: & 
      stat_rknd,  & ! Variable(s)
      stat_nknd
  
  implicit none

  private ! Set Default Scope

  public :: stats, &
            stat_assign, &
            stat_update_var, &
            stat_update_var_pt, &
            stat_begin_update, &
            stat_begin_update_pt, &
            stat_end_update, &
            stat_end_update_pt, &
            stat_modify, &
            stat_modify_pt

  ! Derived data types to store GrADS/netCDF statistics
  type stats

    ! Number of fields to sample
    integer :: nn

    ! Vertical extent of variable
    integer :: kk

    ! Vertical levels
    real, pointer, dimension(:) :: z

    ! Array to store sampled fields

    real(kind=stat_rknd), pointer, dimension(:,:,:,:) :: x

    integer(kind=stat_nknd), pointer, dimension(:,:,:,:) :: n
    
    ! Tracks if a field is in the process of an update
    logical, pointer, dimension(:,:,:,:) :: l_in_update

    ! Data for GrADS / netCDF output

    type (stat_file) f

  end type stats

  contains

  !=============================================================================
  subroutine stat_assign( var_index, var_name,  & 
                     var_description, var_units, grid_kind )

    ! Description: Assigns pointers for statistics variables in grid.
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
    type(stats), intent(inout) :: grid_kind

    grid_kind%f%var(var_index)%ptr => grid_kind%x(:,:,:,var_index)
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

  !=============================================================================
  subroutine stat_update_var( var_index, value, grid_kind )

    ! Description:
    ! This updates the value of a statistics variable located at var_index 
    ! associated with grid type 'grid_kind' (zt, zm, or sfc).
    !
    ! This subroutine is used when a statistical variable needs to be updated
    ! only once during a model timestep.
    !
    ! In regards to budget terms, this subroutine is used for variables that
    ! are either completely implicit (e.g. wprtp_ma) or completely explicit 
    ! (e.g. wp2_pr3).  For completely implicit terms, once the variable has been
    ! solved for, the implicit contribution can be finalized.  The finalized
    ! implicit contribution is sent into stat_update_var_pt.  For completely
    ! explicit terms, the explicit contribution is sent into stat_update_var_pt
    ! once it has been calculated.
    !---------------------------------------------------------------------            
   
    implicit none

    ! Input Variables(s)

    integer, intent(in) ::  & 
      var_index ! The index at which the variable is stored  []

    ! Input/Output Variable(s)
    type(stats), intent(inout) ::  & 
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc )

    ! Input Variable(s) NOTE: Due to the implicit none above, these must
    ! be declared below to allow the use of grid_kind

    real, dimension(grid_kind%kk), intent(in) :: & 
      value ! Value of field being added to the statistic    [Units Vary]  

    integer :: k

    if ( var_index > 0 ) then
       do k = 1, grid_kind%kk
         grid_kind%x(1,1,k,var_index) =  & 
              grid_kind%x(1,1,k,var_index) + value(k)
         grid_kind%n(1,1,k,var_index) =  & 
              grid_kind%n(1,1,k,var_index) + 1
       end do
    endif

  end subroutine stat_update_var
   
  !=============================================================================
  subroutine stat_update_var_pt( var_index, grid_level, value, grid_kind )

    ! Description:
    ! This updates the value of a statistics variable located at var_index 
    ! associated with grid type 'grid_kind' at a specific grid_level.
    !
    ! See the description of stat_update_var for more details.
    !---------------------------------------------------------------------            
  
    implicit none

    ! Input Variables(s)

    integer, intent(in) ::  & 
      var_index,    & ! The index at which the variable is stored           []
      grid_level      ! The level at which the variable is to be modified   []
   
    real, intent(in) :: & 
      value ! Value of field being added to the statistic         [Units Vary]

    ! Input/Output Variable(s)
    type(stats), intent(inout) ::  & 
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).
  
    if ( var_index > 0 ) then

       grid_kind%x(1,1,grid_level,var_index) = grid_kind%x(1,1,grid_level,var_index) + value

       grid_kind%n(1,1,grid_level,var_index) = grid_kind%n(1,1,grid_level,var_index) + 1

    endif

  end subroutine stat_update_var_pt
       
  !=============================================================================
  subroutine stat_begin_update( var_index, value, &
                                grid_kind )

    ! Description:
    ! This begins an update of the value of a statistics variable located at 
    ! var_index on the (zt, zm, or sfc) grid.  It is used in conjunction with
    ! subroutine stat_end_update.
    !
    ! This subroutine is used when a statistical variable needs to be updated
    ! more than one time during a model timestep.  Commonly, this is used for 
    ! beginning a budget term calculation.
    !
    ! In this type of stats calculation, we first subtract the field 
    ! (e.g. rtm / dt ) from the statistic, then update rtm by a term 
    ! (e.g. clip rtm), and then re-add the field (e.g. rtm / dt) to the 
    ! statistic. 
    !            
    ! Example:
    !
    !  call stat_begin_update( irtm_bt, real(rtm / dt), zt ) 
    !
    !  !!! Perform clipping of rtm !!!
    !          
    !  call stat_end_update( irtm_bt, real(rtm / dt), zt )
    !
    ! This subroutine is often used with stats budget terms for variables that
    ! have both implicit and explicit components (e.g. wp3_ta).  The explicit
    ! component is sent into stat_begin_update_pt (with the sign reversed 
    ! because stat_begin_update_pt automatically subtracts the value sent into 
    ! it).  Then, once the variable has been solved for, the implicit 
    ! statistical contribution can be finalized.  The finalized implicit 
    ! component is sent into stat_end_update_pt.
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
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).

    integer :: i

    do i = 1, gr%nnzp 

       call stat_begin_update_pt & 
             ( var_index, i, value(i), grid_kind )

    enddo

  end subroutine stat_begin_update

  !=============================================================================
  subroutine stat_begin_update_pt & 
             ( var_index, grid_level, value, grid_kind )

    ! Description:
    ! This begins an update of the value of a statistics variable located at 
    ! var_index associated with the grid type (grid_kind) at a specific 
    ! grid_level.  It is used in conjunction with subroutine stat_end_update_pt.
    !
    ! Commonly this is used for beginning a budget.  See the description of
    ! stat_begin_update for more details.
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
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).

    ! Local Variable
    character(len=3) :: char_index
  
    if ( var_index > 0 ) then  ! Are we storing this variable?

       if ( .not. grid_kind%l_in_update(1,1,grid_level,var_index) ) then ! Can we begin an update?

          grid_kind%x(1,1,grid_level, var_index) =  & 
                  grid_kind%x(1,1,grid_level, var_index) - value

          grid_kind%l_in_update(1,1,grid_level, var_index) = .true.  ! Start Record

       else

          write(char_index(1:3),'(i3.3)') var_index
          call clubb_debug( 1, & 
          "Beginning an update before finishing previous.  "// & 
          "Var index = "// char_index   )

       endif

    endif

  end subroutine stat_begin_update_pt

  !=============================================================================
  subroutine stat_end_update( var_index, value, grid_kind )

    ! Description:
    ! This ends an update of the value of a statistics variable located at 
    ! var_index on the (zt, zm, or sfc) grid.  It is used in conjunction with
    ! subroutine stat_begin_update.
    !
    ! This subroutine is used when a statistical variable needs to be updated
    ! more than one time during a model timestep.  Commonly, this is used for 
    ! finishing a budget term calculation.
    !
    ! In this type of stats calculation, we first subtract the field 
    ! (e.g. rtm / dt ) from the statistic, then update rtm by a term 
    ! (e.g. clip rtm), and then re-add the field (e.g. rtm / dt) to the 
    ! statistic. 
    !            
    ! Example:
    !
    !  call stat_begin_update( irtm_bt, real(rtm / dt), zt ) 
    !
    !  !!! Perform clipping of rtm !!!
    !          
    !  call stat_end_update( irtm_bt, real(rtm / dt), zt )
    !
    ! This subroutine is often used with stats budget terms for variables that
    ! have both implicit and explicit components (e.g. wp3_ta).  The explicit
    ! component is sent into stat_begin_update_pt (with the sign reversed 
    ! because stat_begin_update_pt automatically subtracts the value sent into 
    ! it).  Then, once the variable has been solved for, the implicit 
    ! statistical contribution can be finalized.  The finalized implicit 
    ! component is sent into stat_end_update_pt.
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
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).

    integer :: i 
           
    do i = 1,gr%nnzp             
        call stat_end_update_pt & 
                 ( var_index, i, value(i), grid_kind )
    enddo

   end subroutine stat_end_update

  !=============================================================================
  subroutine stat_end_update_pt & 
                ( var_index, grid_level, value, grid_kind )

    ! Description:
    ! This ends an update of the value of a statistics variable located at 
    ! var_index associated with the grid type (grid_kind) at a specific 
    ! grid_level.  It is used in conjunction with subroutine 
    ! stat_begin_update_pt.
    !
    ! Commonly this is used for finishing a budget.  See the description of
    ! stat_end_update for more details.
    !---------------------------------------------------------------------            

    use error_code, only: clubb_debug ! Procedure(s)

    implicit none

    ! Input Variables(s)

    integer, intent(in) ::  & 
      var_index,   & ! The index at which the variable is stored           []
      grid_level     ! The level at which the variable is to be modified   []     
   
    real, intent(in) :: & 
      value       ! Value of field being added to the statistic         [Units Vary]

    ! Input/Output Variable(s)
    type(stats), intent(inout) ::  & 
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).

    if ( var_index > 0 ) then ! Are we storing this variable?
           
       if ( grid_kind%l_in_update(1,1,grid_level,var_index) ) then ! Can we end an update? 
               
          call stat_update_var_pt & 
                   ( var_index, grid_level, value, grid_kind ) 
        
          grid_kind%l_in_update(1,1,grid_level,var_index) = .false. ! End Record

       else

          call clubb_debug( 1, "Ending before beginning update. For variable "// &
          grid_kind%f%var(var_index)%name )

       endif

    endif
    
   end subroutine stat_end_update_pt
   
  !=============================================================================
  subroutine stat_modify( var_index, value, & 
                          grid_kind )

    ! Description:
    ! This modifies the value of a statistics variable located at var_index on 
    ! the (zt, zm, or sfc) grid.  It does not increment the sampling count.
    !
    ! This subroutine is normally used when a statistical variable needs to be 
    ! updated more than twice during a model timestep.  Commonly, this is used 
    ! if a budget term calculation needs an intermediate modification between
    ! stat_begin_update and stat_end_update.
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
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).

    integer :: i

    do i = 1, gr%nnzp 

       call stat_modify_pt( var_index, i, value(i), grid_kind )

    enddo

  end subroutine stat_modify

  !=============================================================================
  subroutine stat_modify_pt( var_index, grid_level, value, & 
                             grid_kind )

    ! Description:
    ! This modifies the value of a statistics variable located at var_index on 
    ! the grid at a specific point. It does not increment the sampling count.
    !
    ! Commonly this is used for intermediate updates to a budget.  See the 
    ! description of stat_modify for more details.
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
      grid_kind ! Which grid the variable is located on (zt, zm, rad, or sfc).
  
    if ( var_index > 0 ) then 

       grid_kind%x(1,1,grid_level,var_index )  & 
          = grid_kind%x(1,1,grid_level,var_index ) + value

    endif
  
  end subroutine stat_modify_pt
   
!===============================================================================
 
end module stats_type
