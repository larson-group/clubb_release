! $Id: gradsaverage.F90,v 1.3 2008-07-29 16:44:00 nielsenb Exp $
module grads_common

implicit none

public :: grads_average, grads_average_interval, grads_zlvl

private ! Default Scope

contains

!----------------------------------------------------------------------
function grads_average( filename, nz, t1, t2, variable_name,  & 
                        npower, error )
!       Description:
!       Average a GrADS file variable over the interval t1 to t2

!       References:
!       None
!----------------------------------------------------------------------

use constants, only: fstderr ! Variable(s)

use inputfile_class, only: inputgrads ! Type(s)

use inputfile_class, only: open_grads_read, get_var,  & ! Procedures
                     close_grads_read

implicit none

! Input Variables
character(len=*), intent(in) ::  & 
  filename ! Name of the file

integer, intent(in) ::  & 
  nz,  & ! Number of vertial levels in the GrADS file.
  t1,  & ! Beginning timestep to look at
  t2  ! Ending timestep to look at

character(len=*), intent(in) ::  & 
  variable_name ! Name of the variable to read in

integer, intent(in) ::  & 
  npower ! Exponent operator, must be 1 or 2

! Output Variable
logical, intent(out) ::  & 
  error ! error status for this function

! Return Variable for function
real, dimension(nz) :: grads_average

! Local Variables
type (inputgrads) :: faverage ! Data file derived type

real, dimension(nz) :: grads_temp ! Temporary variable

integer ::  & 
  i,  & ! les_array loop index 
  t  ! timestep loop index

integer ::  & 
  num_timesteps ! steps between t1 and t2 

!-----------------------------------------------------------------------

! Initialize variables
num_timesteps = ( t2 - t1 ) + 1
grads_average = 0.0

! Open grads file
call open_grads_read( 10, filename, faverage )

! Read in variables from GrADS file
do t = t1, t2
  call get_var( faverage, variable_name, t,  & 
                grads_temp(1:nz), error ) 
  if ( error ) then
     write(fstderr,*) "grads_average: get_var failed for "  & 
       //trim( variable_name )//" in "//trim( filename ) & 
       //" at time=", t
     return 
  end if

  if ( npower == 1 ) then
    grads_average(1:nz)  & 
    = grads_average(1:nz) + grads_temp(1:nz)

  else if ( npower == 2 ) then
    grads_average(1:nz)  & 
    = grads_average(1:nz) + grads_temp(1:nz)*grads_temp(1:nz)

  else
    write(fstderr,*) "gradsaverage: invalid npower = ", npower
    error = .true.
    return

  end if ! npower

end do ! t = t1, t2

! Close GrADS file
call close_grads_read( faverage )

! Take average over num_timesteps
grads_average(1:nz)  & 
= grads_average(1:nz) / real( num_timesteps )

return
end function grads_average

!-------------------------------------------------------------------------
function grads_average_interval & 
         ( filename, nz, t, variable_name, & 
           npower, error )

!       Description:
!       Reads in GrADS data from a file and then takes several averages 
!       over an interval.

!       References:
!       None

!       Notes:
!       The variable t is assumed size, which needs to be used with
!       caution.
!-------------------------------------------------------------------------
use constants, only: fstderr ! Variable(s)

implicit none

! Constant Parameters
integer, parameter ::  & 
  tmax = huge( 1 ) ! Sanity check for huge t array

! Input Variables
character(len=*), intent(in) ::  & 
  filename ! Name of the file

integer, intent(in) ::  & 
  nz ! Number of vertical grid levels

integer, dimension(:), intent(in) ::  & 
  t ! Timesteps to use for taking an average over an interval

character(len=*), intent(in) ::  & 
  variable_name  ! Name of the variable to read in

integer, intent(in) ::  & 
  npower ! exponent applied to data retrieved from the file( 1 or 2) 

! Output Variables
logical, intent(out) ::  & 
  error ! status of this function

! Return Variables
real, dimension(nz) ::  & 
  grads_average_interval

! Local Variables 
real, dimension(nz) :: grads_temp

integer ::  & 
  i,       & ! Loop variable 
  tdim,    & ! Dimension to read over for t variable
  divisor

!-------------------------------------------------------------------------

! Sanity check
if ( size( t ) > tmax .or. size( t ) < 2 ) then
  write(unit=fstderr,fmt=*)  & 
    "grads_average_interval: Invalid time interval"
  error = .true.
  return
end if

do i = 1, size( t )
  if ( t( i ) == 0 ) exit
  tdim = i
end do

grads_average_interval & 
= grads_average & 
  ( filename, nz, t(1), t(2), variable_name, npower, error )  & 
  * ( t(2) - t(1) )

divisor = t(2) - t(1)

if ( error ) return

do i=3, tdim, 2 
  grads_temp = grads_average & 
               ( filename, nz, t(i), t(i+1),  & 
                 variable_name, npower, error )
  grads_average_interval  & 
  = grads_average_interval + grads_temp * ( t(i+1) - t(i) )
  divisor = divisor + ( t(i+1) - t(i) )
end do

grads_average_interval(1:nz)  & 
= grads_average_interval(1:nz) / real( divisor )

return
end function grads_average_interval

!-------------------------------------------------------------------------
integer function grads_zlvl( filename )

!       Description:
!       Returns a integer for the number of vertical levels in a file

!       References:
!       None

!       Notes:
!       Somewhat inefficient, since it reads in the entire .ctl file to
!       determine the number of levels
!-------------------------------------------------------------------------

use inputfile_class, only: inputgrads ! Type(s)

use inputfile_class, only: open_grads_read, close_grads_read ! Procedure(s)

implicit none

! Input Variables
character(len=*), intent(in) ::  & 
  filename ! File name

! Local Variables
type (inputgrads) :: fz            ! Data file

! Read in the control file
call open_grads_read( 10, filename, fz )

! Set return variable
grads_zlvl = fz%iz

! Close file
call close_grads_read( fz )

return
end function grads_zlvl

end module grads_common
