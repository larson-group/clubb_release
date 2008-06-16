!-----------------------------------------------------------------------
! $Id: stats_type.F,v 1.9 2008-06-16 22:19:11 faschinj Exp $
      module stats_type
#ifdef STATS
 
!     Description:
!     Contains derived data type 'stats'
!-----------------------------------------------------------------------

        use outputfile_class, only:
     .      outputfile ! Type

        use stats_prec, only:
     .      stat_rknd, ! Variable(s)
     .      stat_nknd
        
        implicit none

        private ! Set Default Scope

        public :: stats, stat_assign, stat_update_var,
     .            stat_update_var_pt, stat_begin_update,
     .            stat_end_update, stat_begin_update_pt,
     .            stat_end_update_pt, stat_modify_pt


        ! Derived data types to store GrADS/netCDF statistics
        type stats

          ! Number of fields to sample
          integer :: nn

          ! Vertical extent of variable
          integer :: kk

          ! Vertical levels
          real, pointer, dimension(:) :: z

!         Array to store sampled fields

          real(kind=stat_rknd), pointer, dimension(:,:) :: x

          integer(kind=stat_nknd), pointer, dimension(:,:) :: n
          
          ! Tracks if a field is in the process of an update
          logical, pointer, dimension(:,:) :: in_update

!         Data for GrADS output

          type (outputfile) f

        end type stats


        contains
!-----------------------------------------------------------------------
         subroutine stat_assign( var_index, name, 
     .                     description, units, grid )

!        Description: Assigns pointers for statistics variables in grid
!    
!-----------------------------------------------------------------------
         implicit none
           
         ! Input Variables
             
         integer,intent(in) :: var_index               ! Variable index
         character(len = *), intent(in) :: name        ! Variable name
         character(len = *), intent(in) :: description ! Variable description
         character(len = *), intent(in) :: units       ! Variable units

         ! Output Variable
         type(stats), intent(out) :: grid

         grid%f%var(var_index)%ptr => grid%x(:,var_index)
         grid%f%var(var_index)%name = name
         grid%f%var(var_index)%description = description
         grid%f%var(var_index)%units = units

             !Example of the old format 
             !changed by Joshua Fasching 23 August 2007
             
             !zt%f%var(ithlm)%ptr => zt%x(:,k)
             !zt%f%var(ithlm)%name = "thlm"
             !zt%f%var(ithlm)%description = "thetal (K)"
             !zt%f%var(ithlm)%units = "K"
      
         return

         end subroutine stat_assign

!---------------------------------------------------------------------
         subroutine stat_update_var( var_index, value, grid )
!
!        Description: This updates the value of a statistics variable 
!        located at var_index in the grid
!
!---------------------------------------------------------------------            
         use grid_class, only: gr ! Variable(s)
   
         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
                  
         real, dimension(gr%nnzp), intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid

         if ( var_index > 0 ) then
            grid%x(:,var_index) = grid%x(:,var_index) + value
            grid%n(:,var_index) = grid%n(:,var_index) + 1
         end if

         end subroutine stat_update_var
         
!---------------------------------------------------------------------
         subroutine stat_update_var_pt
     .             ( var_index, point, value, grid )
!
!        Description: This updates the value of a statistics variable 
!        located at var_index in the grid at a specific point.
!
!---------------------------------------------------------------------            
        
         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index, ! The index at which the variable is stored
     .   point      ! The level the variable is to be modified at
         
         real, intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid
        
         if ( var_index > 0 ) then
            grid%x( point, var_index ) = 
     .        grid%x( point ,var_index) + value
            grid%n( point, var_index ) = 
     .         grid%n( point, var_index ) + 1
         end if

         end subroutine stat_update_var_pt
             
!---------------------------------------------------------------------
         subroutine stat_begin_update( var_index, value, grid )
!
!        Description: This begins an update of the value of a 
!        statistics variable located at var_index in the grid.
!        Commonly this is used for beginning a budget. 
!       
!         Ex.
!          call stat_begin_update( irtm_bt, real(rtm / dt), zt ) 
!
!          !!! Perform clipping of statistic !!!
!          
!          call stat_end_update( irtm_bt, real(rtm / dt), zt )
!---------------------------------------------------------------------            
         use grid_class, only: gr          

         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         real, dimension(gr%nnzp), intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid

         integer :: i

         do i = 1, gr%nnzp 

           call stat_begin_update_pt( var_index, i, value(i), grid )

         enddo

         end subroutine stat_begin_update

!---------------------------------------------------------------------
         subroutine stat_begin_update_pt
     .             ( var_index, i, value, grid )
!
!        Description: This begins an update of the value of a 
!        statistics variable located at var_index in the grid 
!        at a specific point i.
!        Commonly this is used for beginning a budget. See stat_begin_update 
!---------------------------------------------------------------------            
 
         use error_code, only: clubb_debug

         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         
         real, intent(in) ::
     .   value ! The value being stored

         integer, intent(in) :: i
         
         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid
        
         if ( var_index > 0 ) then
             if ( .not. grid%in_update(i,var_index) ) then  
               grid%x(i,var_index) = 
     .                 grid%x(i,var_index) - value
               grid%in_update(i,var_index) = .true.
             else
               print *, var_index
              call clubb_debug( 2,
     .        "Warning! Beginning an update before finishing previous")
             end if
         end if
         end subroutine stat_begin_update_pt

!---------------------------------------------------------------------
         subroutine stat_end_update( var_index, value, grid )
!
!        Description: This ends the update of a statistics variable 
!        located at var_index in the grid. Commoly used to end the
!        storage of a budget term. See stat_begin_update
!
!---------------------------------------------------------------------            

         use grid_class, only: gr ! Variable(s)

         use error_code, only: clubb_debug ! Procedure(s)

         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         
         real, dimension(gr%nnzp), intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid

         integer :: i 
                 
         do i = 1, gr%nnzp            
              call stat_end_update_pt( var_index, i, value(i), grid )
         end do


         end subroutine stat_end_update

!---------------------------------------------------------------------
         subroutine stat_end_update_pt( var_index, i, value, grid )
!
!        Description: This ends an update of the value of a 
!        statistics variable located at var_index in the grid 
!        at a specific point i.
!        Commonly this is used for ending a budget. See stat_begin_update 
!
!---------------------------------------------------------------------            

         use error_code, only: clubb_debug ! Procedure(s)

         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         
         real, intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid

         integer, intent(in) :: i 

         if ( var_index > 0 ) then
                 
             if ( grid%in_update(i, var_index) ) then
                     
                call stat_update_var_pt( var_index, i, value, grid )
              
                grid%in_update( i, var_index ) = .false.

              else

                call clubb_debug(2,"Ending before beginning update")

              endif

         endif
    
         end subroutine stat_end_update_pt
         
!---------------------------------------------------------------------
         subroutine stat_modify( var_index, value, grid )
!
!        Description: This modifies the value of a statistics variable 
!        located at var_index in the grid. It does not increment the
!        sampling count.
!
!---------------------------------------------------------------------            
                  
         use grid_class, only: gr ! Variable(s)

         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         real, dimension(gr%nnzp), intent(in) ::
     .   value ! The value being stored

         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid

         integer :: i

         do i = 1, gr%nnzp 

          call stat_modify_pt( var_index, i, value(i), grid )

         enddo

         end subroutine stat_modify

!---------------------------------------------------------------------
         subroutine stat_modify_pt
     .             ( var_index, i, value, grid )
!
!        Description: This modifies the value of a statistics variable 
!        located at var_index in the grid at a specific point. It does
!        not increment the sampling count.
!
!---------------------------------------------------------------------            
 
         implicit none

         ! Input Variables(s)

         integer, intent(in) :: 
     .   var_index ! The index at which the variable is stored
         
         
         real, intent(in) ::
     .   value ! The value being stored

         integer, intent(in) :: i
         
         ! Input/Output Variable(s)
         type(stats), intent(inout) :: grid
        
         if ( var_index > 0 ) then 
           grid%x(i,var_index) = grid%x(i,var_index) + value
         end if
        
         end subroutine stat_modify_pt
         
#endif /*STATS*/
      end module stats_type
