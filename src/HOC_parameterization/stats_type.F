!-----------------------------------------------------------------------
! $Id: stats_type.F,v 1.1 2007-09-08 20:57:57 faschinj Exp $
#ifdef STATS
      module stats_type

!     Description:
!     Contains derived data type 'stats'
!-----------------------------------------------------------------------

        use outputfile_class

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
      end module stats_type
#endif /*STATS*/
