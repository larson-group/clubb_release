!-----------------------------------------------------------------------
! $Id$

module output_writer

! Description: Outputs information to the screen and optionally to a file
!-----------------------------------------------------------------------

  implicit none

  public :: write_output

  private ! Default to private

  interface write_output
    module procedure write_output_text, write_output_real, write_output_real_array, &
      write_output_integer, write_output_logical
  end interface

  contains

  !----------------------------------------------------------------------
  subroutine write_output_text( text, l_write_to_file, iunit )
    ! Description:
    !   Outputs a string
    ! References:
    !   None
    !--------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text
    logical, intent(in) :: l_write_to_file
    integer, intent(in) :: iunit

    print *, trim( text )
    
    if( l_write_to_file ) then
      write(iunit, *) trim( text )
    end if

  end subroutine write_output_text

  !----------------------------------------------------------------------
  subroutine write_output_real( text, value, l_write_to_file, iunit )
    ! Description:
    !   Outputs a string and a real
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text
    real, intent(in) :: value
    logical, intent(in) :: l_write_to_file
    integer, intent(in) :: iunit

    print *, text, value
    
    if( l_write_to_file ) then
      write(iunit, *) text, value
    end if

  end subroutine write_output_real

  !----------------------------------------------------------------------
  subroutine write_output_real_array( text, value, l_write_to_file, iunit )
    ! Description:
    !   Outputs a string and a real array
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text
    real, dimension(:), intent(in) :: value
    logical, intent(in) :: l_write_to_file
    integer, intent(in) :: iunit

    print *, text, value
    
    if( l_write_to_file ) then
      write(iunit, *) text, value
    end if

  end subroutine write_output_real_array

  !----------------------------------------------------------------------
  subroutine write_output_integer( text, value, l_write_to_file, iunit )
    ! Description:
    !   Outputs a string and an integer
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text
    integer, intent(in) :: value
    logical, intent(in) :: l_write_to_file
    integer, intent(in) :: iunit

    print *, text, value
    
    if( l_write_to_file ) then
      write(iunit, *) text, value
    end if

  end subroutine write_output_integer

  !----------------------------------------------------------------------
  subroutine write_output_logical( text, value, l_write_to_file, iunit )
    ! Description:
    !   Outputs a string an a logical
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text
    logical, intent(in) :: value
    logical, intent(in) :: l_write_to_file
    integer, intent(in) :: iunit

    print *, text, value
    
    if( l_write_to_file ) then
      write(iunit, *) text, value
    end if

  end subroutine write_output_logical

end module output_writer
