!-------------------------------------------------------------------------------
! $Id$
!-------------------------------------------------------------------------------

module error_code

! Description:
!   Since f90/95 lacks enumeration, we're stuck numbering each
!   error code by hand like this.

!   We are "enumerating" error codes to be used with CLUBB. Adding
!   additional codes is as simple adding an additional integer
!   parameter. The error codes are ranked by severity, the higher
!   number being more servere. When two errors occur, assign the
!   most servere to the output.

!   This code also handles subroutines related to debug_level. See
!   the 'set_clubb_debug_level' description for more detail.

! References:
!   None
!-------------------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: & 
    reportError,  & 
    fatal_error, & 
    lapack_error,     & 
    clubb_at_least_debug_level,  & 
    set_clubb_debug_level, & 
    clubb_debug, &
    getCrcCheckSum1D

  private :: clubb_debug_level

  ! Model-Wide Debug Level
  integer, save :: clubb_debug_level = 0

!$omp threadprivate(clubb_debug_level)

  ! Error Code Values
  integer, parameter, public :: & 
    clubb_no_error                 =  0, & 
    clubb_var_less_than_zero       =  1, & 
    clubb_var_equals_NaN           =  2, & 
    clubb_singular_matrix          =  3, & 
    clubb_bad_lapack_arg           =  4, & 
    clubb_rtm_level_not_found      =  5, & 
    clubb_var_out_of_bounds        =  6, &
    clubb_var_out_of_range         =  7

  contains

!-------------------------------------------------------------------------------
  subroutine reportError( err_code )
!
! Description: 
!   Reports meaning of error code to console.
!
!-------------------------------------------------------------------------------

    use constants_clubb, only: & 
        fstderr ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! ---- Begin Code ----

    select case ( err_code )

    case ( clubb_no_error )
      write(fstderr,*) "No errors reported."

    case ( clubb_var_less_than_zero )
      write(fstderr,*) "Variable in CLUBB is less than zero."

    case ( clubb_singular_matrix )
      write(fstderr,*) "Singular Matrix in CLUBB."

    case ( clubb_var_equals_NaN )
      write(fstderr,*) "Variable in CLUBB is NaN."

    case ( clubb_bad_lapack_arg )
      write(fstderr,*) "Argument passed to a LAPACK procedure is invalid."

    case ( clubb_rtm_level_not_found )
      write(fstderr,*) "rtm level not found"

    case ( clubb_var_out_of_bounds )
      write(fstderr,*) "Input variable is out of bounds."

    case ( clubb_var_out_of_range )
      write(fstderr,*) "A CLUBB variable had a value outside the valid range."

    case default
      write(fstderr,*) "Unknown error: ", err_code

    end select

    return
  end subroutine reportError
!-------------------------------------------------------------------------------
  elemental function lapack_error( err_code )
!
! Description: 
!   Checks to see if the err_code is equal to one
!   caused by an error encountered using LAPACK.
! Reference:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! Input variable
    integer,intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: lapack_error

    ! ---- Begin Code ----

    lapack_error = (err_code == clubb_singular_matrix .or. & 
        err_code == clubb_bad_lapack_arg )

    return
  end function lapack_error

!-------------------------------------------------------------------------------
  elemental function fatal_error( err_code )
!
! Description: Checks to see if the err_code is one that usually
!   causes an exit in other parts of CLUBB.
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: fatal_error

    ! ---- Begin Code ----

    fatal_error = err_code /= clubb_no_error .and. & 
                  err_code /= clubb_var_less_than_zero
    return
  end function fatal_error

!------------------------------------------------------------------	
  logical function clubb_at_least_debug_level( level )
!
! Description:
!   Checks to see if clubb has been set to a specified debug level
!------------------------------------------------------------------
    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    ! ---- Begin Code ----

    clubb_at_least_debug_level = ( level <= clubb_debug_level )

    return
  end function clubb_at_least_debug_level

!-------------------------------------------------------------------------------
  subroutine set_clubb_debug_level( level )
!
!  Description:
!    Accessor for clubb_debug_level
!
!   0 => Print no debug messages to the screen
!   1 => Print lightweight debug messages, e.g. print statements
!   2 => Print debug messages that require extra testing,
!        e.g. checks for NaNs and spurious negative values.
!  References:
!    None
!-------------------------------------------------------------------------------
    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    ! ---- Begin Code ----

    clubb_debug_level = level

    return
  end subroutine set_clubb_debug_level

!-------------------------------------------------------------------------------
  subroutine clubb_debug( level, str )
!
! Description:
!   Prints a message to file unit fstderr if the level is greater
!   than or equal to the current debug level.
!-------------------------------------------------------------------------------
    use constants_clubb, only: & 
        fstderr ! Variable(s)

    implicit none

    ! Input Variable(s)

    character(len=*), intent(in) :: str ! The message being reported

    ! The debug level being checked against the current setting
    integer, intent(in) :: level

    ! ---- Begin Code ----

    if ( level <= clubb_debug_level ) then
      write(fstderr,*) str
    end if

    return
  end subroutine clubb_debug


!-------------------------------------------------------------------------------
function getCrcCheckSum1D(realArray,len) result(crc)
!
! Description:
!   This function computes an integer cyclic redundancy check sum (CRC)
!   for an real valued 1D array using the routine icrc()
!   from Numerical Recipes.
!
! Input:
!   realArray ... array of real values for which a CRC should be generated
!   len       ... length of realArray  
!
! Output:
!   crc       ... CRC checksum of realArray
!
!-------------------------------------------------------------------------------
  
  implicit none

  ! Input Variable(s)
  integer, intent(in)                  :: len
  real, dimension(1:len), intent(in)   :: realArray 

  ! Output Variable(s)                                        
  integer                              :: checksum

  ! Internal Variable(s)
  character(20), dimension(len)        :: strArray
  integer                              :: i
  character(1), dimension(2)           :: crc
  


  ! ---- Begin Code ----
  
  ! Convert real array to string array
  do i=1,len
    write(strArray(i),*) realArray(i) 
  end do
 
  ! calculate checksum
  checksum = icrc(crc,strArray,len,1,1)
  return

end function getCrcCheckSum1D


!-------------------------------------------------------------------------------
function icrc(crc,bufptr,len,jinit,jrev)

! Description:
!   This routine is taken from Numerical Recipes in Fortran 77. Please see
!   http://www.haoli.org/nr/bookfpdf.html or  our trac system 
!   wrf:ticket:19#comment:26.
!   It is written in Fortran77 and I didn't rewrite it. 
!
!   Original description: 
!   Computes a 16-bit Cyclic Redundancy Check for an array bufptr of length len 
!   bytes, using any of several conventions as determined by the settings of jinit
!   and jrev (see accompanying table). The result is returned both as an integer
!   icrc and as a 2-byte array crc. If jinit is negative, then crc is used on
!   input to initialize the remainder register, in effect concatenating bufptr to
!   the previous call.

!-------------------------------------------------------------------------------
  

  INTEGER icrc,jinit,jrev,len
  CHARACTER*1 bufptr(*),crc(2)


  INTEGER ich,init,ireg,j,icrctb(0:255),it(0:15),icrc1,ib1,ib2,ib3
  CHARACTER*1 creg(4),rchr(0:255)
  SAVE icrctb,rchr,init,it,ib1,ib2,ib3

  EQUIVALENCE (creg,ireg) 

  !Table of 4-bit bit-reverses, and flag for initialization.
  DATA it/0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15/, init /0/
  
  if (init.eq.0) then
    init=1
    ireg=256*(256*ichar('3')+ichar('2'))+ichar('1')
    do j=1,4
      if (creg(j).eq.'1') ib1=j
      if (creg(j).eq.'2') ib2=j
      if (creg(j).eq.'3') ib3=j
    enddo
    
    do j=0,255
      ireg=j*256
      icrctb(j)=icrc1(creg,char(0),ib1,ib2,ib3)
      ich=it(mod(j,16))*16+it(j/16)
      rchr(j)=char(ich)
    enddo
  endif

  if (jinit.ge.0) then
    crc(1)=char(jinit)
    crc(2)=char(jinit)
  else if (jrev.lt.0) then
    ich=ichar(crc(1))
    crc(1)=rchr(ichar(crc(2)))
    crc(2)=rchr(ich)
  endif

  do j=1,len
    ich=ichar(bufptr(j))
    if(jrev.lt.0)ich=ichar(rchr(ich))
    ireg=icrctb(ieor(ich,ichar(crc(2))))
    crc(2)=char(ieor(ichar(creg(ib2)),ichar(crc(1))))
    crc(1)=creg(ib1)
  enddo

  if (jrev.ge.0) then
    creg(ib1)=crc(1)
    creg(ib2)=crc(2)
  else
    creg(ib2)=rchr(ichar(crc(1)))
    creg(ib1)=rchr(ichar(crc(2)))
    crc(1)=creg(ib1)
    crc(2)=creg(ib2)
  endif

  icrc=ireg
  return
end function icrc







!-------------------------------------------------------------------------------
function icrc1(crc,onech,ib1,ib2,ib3)

! Description:
!   This routine is taken from Numerical Recipes in Fortran 77. Please see
!   http://www.haoli.org/nr/bookfpdf.html or  our trac system 
!   wrf:ticket:19#comment:26.
!   It is written in Fortran77 and I didn't rewrite it. 
!
!   Original description: 
!   Given a remainder up to now, return the new CRC after one character is added. 
!   This routine is functionally equivalent to icrc(,,1,-1,1) , but slower. It is
!   used by icrc to initialize its table.

!-------------------------------------------------------------------------------

  INTEGER icrc1,ib1,ib2,ib3

  INTEGER i,ichr,ireg
  CHARACTER*1 onech,crc(4),creg(4)
  EQUIVALENCE (creg,ireg)
  ireg=0
  creg(ib1)=crc(ib1)
  !Here is where the character is folded into the register.
  creg(ib2)=char(ieor(ichar(crc(ib2)),ichar(onech)))
  do i=1,8
    !Here is where 8 one-bit shifts, and some XORs with the gen-
    !erator polynomial, are done.
    ichr=ichar(creg(ib2))
    ireg=ireg+ireg
    creg(ib3)=char(0)

    if(ichr.gt.127)ireg=ieor(ireg,4129)
  enddo
  
  icrc1=ireg
  return
end function icrc1


end module error_code
!-------------------------------------------------------------------------------
