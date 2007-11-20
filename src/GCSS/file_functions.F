      module file_functions
      implicit none

      public :: file_read_1d, file_read_2d
      contains

      subroutine file_read_1d (device,file_path,file_dimension1,
     .                       file_per_line,variable)
      
      integer, intent(in) ::
     . device,
     . file_dimension1,
     . file_per_line

      character(*), intent(in) ::
     . file_path

      real, dimension(file_dimension1), intent(out) ::
     . variable

      integer i, k

      open(device,file=file_path,action='read')

! Michael Falk wrote this routine to read data files in a particular format for mpace_a.
! Each line has a specific number of values, until the last line in the file, which
! has the last few values and then ends.  This reads the correct number of values on
! each line.  24 September 2007

      do k=0,((file_dimension1/file_per_line)-1)			 ! reads file_per_line values from each line of the input file
        read(device,*) (variable((k*file_per_line)+i),i=1,file_per_line) ! does that for every full line
      end do
      read (device,*) (variable((k*file_per_line)+i),
     .     i=1,(mod(file_dimension1,file_per_line)))                     ! reads the partial line

      close (device)

      end subroutine file_read_1d





      subroutine file_read_2d (device,file_path,file_dimension1,
     .                       file_dimension2,file_per_line,variable)
      
      integer, intent(in) ::
     . device,
     . file_dimension1,
     . file_dimension2,
     . file_per_line

      character(*), intent(in) ::
     . file_path

      real, dimension(file_dimension1,file_dimension2), intent(out) ::
     . variable

      integer i, j, k

      open(device,file=file_path,action='read')

! Michael Falk wrote this routine to read data files in a particular format for mpace_a.
! The 2d mpace_a files list the (file_dimension2) values on a given vertical level, then
! moves to the next level to list its values.
! Each line has a specific number of values, until the last line on a level, which
! is short-- it has the last few values and then a line break.  The next line, beginning
! the next level, is full-sized again.  24 September 2007

      do k=1,(file_dimension1)                                     ! For each level in the data file,
        do j=0,((file_dimension2/file_per_line)-1)
          read(device,*) (variable(k,(j*file_per_line)+i),         ! read file_per_line values in,
     .        i=1,file_per_line)
        end do
        read (device,*) (variable(k,(j*file_per_line)+i),          ! then read the partial line
     .        i=1,(mod(file_dimension2,file_per_line)))
      end do                                                       ! and then start over at the next level.

      close (device)

      end subroutine file_read_2d

      end module file_functions
