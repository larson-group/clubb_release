!-----------------------------------------------------------------------
! $Id$

module output_writer

! Description: Outputs information to the screen and optionally to a file
!-----------------------------------------------------------------------

  implicit none

  public :: write_output, write_date

  private ! Default to private

  interface write_output
    module procedure write_output_text, write_output_real, write_output_real_array, &
      write_output_integer, write_output_logical
  end interface

  contains

  !----------------------------------------------------------------------
  subroutine write_output_text( text, l_write_to_file, iunit, disp_format )
    ! Description:
    !   Outputs a string
    ! References:
    !   None
    !--------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text ! The text to write/display
    logical, intent(in) :: l_write_to_file ! Whether or not to write to a file
    integer, intent(in) :: iunit           ! The file to write to
    character(len = *), optional, intent(in) :: disp_format ! The way to format the text

    if( present( disp_format ) ) then
      write(unit=*, fmt=disp_format) trim( text )
    else
      print *, trim( text )
    endif

    if( l_write_to_file ) then
      if( present( disp_format ) ) then
        write(unit=iunit, fmt=disp_format) trim( text )
      else
        write(iunit, *) trim( text )
      endif
    end if

  end subroutine write_output_text

  !----------------------------------------------------------------------
  subroutine write_output_real( text, value, l_write_to_file, iunit, disp_format )
    ! Description:
    !   Outputs a string and a real
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text ! The text to write/display
    real, intent(in) :: value              ! The value to write/display
    logical, intent(in) :: l_write_to_file ! Whether or not to write to a file
    integer, intent(in) :: iunit           ! The file to write to
    character(len = *), optional, intent(in) :: disp_format ! The way to format the text

    if( present( disp_format ) ) then
      write(unit=*, fmt=disp_format) text, value
    else
      print *, text, value
    endif

    if( l_write_to_file ) then
      if( present( disp_format ) ) then
        write(unit=iunit, fmt=disp_format) text, value
      else
        write(iunit, *) text, value
      endif
    end if

  end subroutine write_output_real

  !----------------------------------------------------------------------
  subroutine write_output_real_array( text, value, l_write_to_file, iunit, disp_format )
    ! Description:
    !   Outputs a string and a real array
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text  ! The text to write/display
    real, dimension(:), intent(in) :: value ! The value to write/display
    logical, intent(in) :: l_write_to_file  ! Whether or not to write to a file
    integer, intent(in) :: iunit            ! The file to write to
    character(len = *), optional, intent(in) :: disp_format ! The way to format the text

    if( present( disp_format ) ) then
      write(unit=*, fmt=disp_format) text, value
    else
      print *, text, value
    endif

    if( l_write_to_file ) then
      if( present( disp_format ) ) then
        write(unit=iunit, fmt=disp_format) text, value
      else
        write(iunit, *) text, value
      endif
    end if
  end subroutine write_output_real_array

  !----------------------------------------------------------------------
  subroutine write_output_integer( text, value, l_write_to_file, iunit, disp_format )
    ! Description:
    !   Outputs a string and an integer
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text ! The text to write/display
    integer, intent(in) :: value           ! The value to write/display
    logical, intent(in) :: l_write_to_file ! Whether or not to write to a file
    integer, intent(in) :: iunit           ! The file to write to
    character(len = *), optional, intent(in) :: disp_format ! The way to format the text

    if( present( disp_format ) ) then
      write(unit=*, fmt=disp_format) text, value
    else
      print *, text, value
    endif

    if( l_write_to_file ) then
      if( present( disp_format ) ) then
        write(unit=iunit, fmt=disp_format) text, value
      else
        write(iunit, *) text, value
      endif
    end if
  end subroutine write_output_integer

  !----------------------------------------------------------------------
  subroutine write_output_logical( text, value, l_write_to_file, iunit, disp_format )
    ! Description:
    !   Outputs a string an a logical
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    character(len = *), intent(in) :: text ! The text to write/display
    logical, intent(in) :: value           ! The value to write/display
    logical, intent(in) :: l_write_to_file ! Whether or not to write to a file
    integer, intent(in) :: iunit           ! The file to write to
    character(len = *), optional, intent(in) :: disp_format ! The way to format the text

    if( present( disp_format ) ) then
      write(unit=*, fmt=disp_format) text, value
    else
      print *, text, value
    endif

    if( l_write_to_file ) then
      if( present( disp_format ) ) then
        write(unit=iunit, fmt=disp_format) text, value
      else
        write(iunit, *) text, value
      endif
    end if
  end subroutine write_output_logical

  !----------------------------------------------------------------------
  subroutine write_date( l_write_to_file, iunit )
    ! Description:
    !   Outputs the current date and time in this format:
    !     YYYY/MM/DD HH:MM:SS
    ! References:
    !   None
    !----------------------------------------------------------------------

    implicit none

    logical, intent(in) :: l_write_to_file   ! Whether or not to write to a file
    integer, intent(in) :: iunit             ! The file to write to
    character(len=8)    :: current_date      ! Current date string (ccyymmdd)
    character(len=10)   :: current_time      ! Current time string (hhmmss.sss)

    call date_and_time( current_date, current_time )

    write(unit=*, fmt="(A4,A1,4(A2,A1),A2)") &
      current_date(1:4), "/", current_date(5:6), "/", current_date(7:8), " ", &
      current_time(1:2), ":", current_time(3:4), ":", current_time(5:6)
 
    if( l_write_to_file ) then
      write(unit=iunit, fmt="(A4,A1,4(A2,A1),A2)") &
        current_date(1:4), "/", current_date(5:6), "/", current_date(7:8), " ", &
        current_time(1:2), ":", current_time(3:4), ":", current_time(5:6)
    endif

  end subroutine write_date

end module output_writer
