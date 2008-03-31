!-----------------------------------------------------------------------
! $Id: stats_type.F,v 1.4 2008-03-31 18:30:12 faschinj Exp $
#ifdef STATS
      module stats_type
 
!     Description:
!     Contains derived data type 'stats'
!-----------------------------------------------------------------------

        use outputfile_class
        
        implicit none

        private
        public :: stats, stat_assign

        ! Derived data types to store GrADS statistics
        type stats 

          ! Number of fields to sample
          integer :: nn

          ! Vertical extent of variable
          integer :: kk

          ! Vertical levels
          real, pointer, dimension(:) :: z

!         Array to store sampled fields

          real, pointer, dimension(:,:) :: x

          integer, pointer, dimension(:,:) :: n

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
           !input
             integer,intent(in) :: var_index               ! Variable index
             character(len = *), intent(in) :: name        ! Variable name
             character(len = *), intent(in) :: description ! Variable description
             character(len = *), intent(in) :: units       ! Variable units

           !output
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

      end module stats_type
#endif /*STATS*/
