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

      do k=0,((file_dimension1/file_per_line)-1)
        read(device,*) (variable((k*file_per_line)+i),i=1,file_per_line)
      end do
      read (device,*) (variable((k*file_per_line)+i),
     .     i=1,(mod(file_dimension1,file_per_line)))

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

!      do k=0,(file_dimension2-1)
      do k=1,(file_dimension1)
        do j=0,((file_dimension2/file_per_line)-1)
          read(device,*) (variable(k,(j*file_per_line)+i),
     .        i=1,file_per_line)
        end do
        read (device,*) (variable(k,(j*file_per_line)+i),
     .        i=1,(mod(file_dimension2,file_per_line)))
      end do

      close (device)

      end subroutine file_read_2d

      end module file_functions
