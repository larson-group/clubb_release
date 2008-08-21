!----------------------------------------------------------------------
!  $Id: compare_runs.F90,v 1.3 2008-07-29 16:44:00 nielsenb Exp $

program compare_runs 

!       Description:
!       Tests selected GrADS file variables between runs

!       References:
!       None
!----------------------------------------------------------------------


use grads_common, only: grads_zlvl, grads_average ! Procedure(s)   

implicit none
 
integer, parameter ::  & 
  nvar = 14, & 
  max_times = 10


real(kind=4), allocatable, dimension(:) ::  & 
  fvar1, & 
  fvar2, & 
  fvar3

character(len=10), dimension(nvar) ::  & 
  hoc_variables, & 
  les_variables

real(kind=8) ::  & 
  percent_mean_sqr, & 
  reference

real(kind=4) :: minmax

character(len=80) ::  & 
  file1, & ! Always a HOC file
  file2, & ! Either a HOC or LES file
  file3 ! Always a HOC file

integer, dimension(max_times) ::  & 
  t ! times to analyze

integer ::  & 
  nz,        & ! smallest of the 3 grids
  nz_file1,  & ! grid levels
  nz_file2,  & ! grid levels
  nz_file3  ! grid levels

integer :: & 
  n, i ! loop

integer :: tmax  ! max time to analyze 

logical :: les_comparison

logical :: error ! set, never used

namelist /compare/ les_comparison, file1, file2, file3, t

!-----------------------------------------------------------------------

open(unit=10,file="compare_runs.in",status='old')
read(unit=10,nml=compare)

hoc_variables(1)  = "rcm"
hoc_variables(2)  = "rtm" 
hoc_variables(3)  = "cf" 
hoc_variables(4)  = "thlm"
hoc_variables(5)  = "um"
hoc_variables(6)  = "vm"
hoc_variables(7)  = "wp3"
hoc_variables(8)  = "wp2rtp"
hoc_variables(9)  = "wp2rcp"
hoc_variables(10) = "wp2zt"
hoc_variables(11) = "thlp2zt"
hoc_variables(12) = "wpthlpzt"
hoc_variables(13) = "wprtpzt"
hoc_variables(14) = "rtp2zt"

if ( les_comparison ) then
  les_variables(1)  = "qcm"
  les_variables(2)  = "qtm" 
  les_variables(3)  = "cf" 
  les_variables(4)  = "thlm"
  les_variables(5)  = "um"
  les_variables(6)  = "vm"
  les_variables(7)  = "wp3"
  les_variables(8)  = "wp2qtp"
  les_variables(9)  = "wp2qlp"
  les_variables(10) = "wp2"
  les_variables(11) = "thlp2"
  les_variables(12) = "wpthlp"
  les_variables(13) = "wpqtp"
  les_variables(14) = "qtp2"
else
  les_variables(1:nvar) = hoc_variables(1:nvar)
end if ! les_comparison

do i = 1, max_times
  if ( t(i) == 0 ) exit
end do

tmax = i - 1
 
write(unit=*,fmt=*) "file1: ", file1
write(unit=*,fmt=*) "file2: ", file2
write(unit=*,fmt=*) "file3: ", file3

do i = 1, tmax, 2
  write(unit=*,fmt=*) "time: ", t(i), " to ", t(i+1)
  write(unit=*,fmt='(1a11)') "Comparison:"
  write(unit=*,fmt='(1a85)')  & 
 "( (file1-file2)/(max-min) )^2   ( (file2-file3)/(max-min) )^2"
  do n = 1, nvar, 1

    ! HOC reference file
    nz_file3 = grads_zlvl( file3 )
    allocate( fvar3(nz_file3) )

    fvar3 = grads_average( file3, nz_file3, t(i), t(i+1),  & 
                           hoc_variables(n), 1, error )

    ! Second HOC file or the LES file
    nz_file2 = grads_zlvl( file2 )
    allocate( fvar2(nz_file2) )

    fvar2 = grads_average( file2, nz_file2, t(i), t(i+1),  & 
                           les_variables(n), 1, error )

    ! First hoc file
    nz_file1 = grads_zlvl( file1 )
    allocate( fvar1(nz_file1) )

    fvar1 = grads_average( file1, nz_file1, t(i), t(i+1),  & 
                           hoc_variables(n), 1, error )

    ! Choose smallest number of vertical gridpoints.
    nz = min( nz_file3, nz_file2, nz_file1 )

    minmax = maxval( fvar2 ) - minval( fvar2 )
    if ( minmax == 0.0 ) then
      write(unit=*,fmt=*) "maxval - minval = 0. for ",  & 
         les_variables(n)
      minmax = 1.0 ! no normalization
    end if

    if ( .not. les_comparison ) then
      percent_mean_sqr =  & 
        sum( (( fvar1(1:nz)-fvar2(1:nz) )/minmax )**2, 1 ) / nz
    else
      call les_interpolate( )
    end if

    reference = & 
      sum( (( fvar2(1:nz)-fvar3(1:nz) )/minmax )**2, 1 ) / nz

    write(unit=*,fmt='(a30,2f18.6)') "Pct mean-square-diff "// & 
      hoc_variables(n)//" = ", percent_mean_sqr, reference

    deallocate( fvar1, fvar2, fvar3 )

  end do   ! n = 1, nvar, 1

end do     ! i = 1, tmax, 2

contains
!-----------------------------------------------------------------------

subroutine les_interpolate( ) 

!       Description:
!       This is more a "shift" than an actual interpolation of grid
!       points, done on a case by case basis.  At some point perhaps
!       this could be changed to do an actual generalized interpolation.

!       References:
!       None
!-----------------------------------------------------------------------

implicit none

select case ( nz )
case( 48, 75, 110, 147, 150, 160 )
  ! Due to HOC's lower starting point, we can only use
  ! (total number of z-levels) - 1 (a maximum of 74 
  ! for BOMEX).  The code below assumes the LES data 
  ! are on an evenly spaced grid.
  ! fvar1 = HOC profile; fvar2 = LES profile.
  ! (Need to interpolate hoc to LES' levels.  Right now we just
  ! compare adjacent z levels.  Vince Larson 12 Jan 2005)
  percent_mean_sqr =  & 
    sum( (( fvar1(2:nz)-fvar2(1:nz-1) )/minmax)**2, 1 ) / (nz-1)
case( 132 ) !  the DYCOMS II RF01 case
  percent_mean_sqr =  & 
    sum( (( fvar1(3:nz)-fvar2(1:nz-2) )/minmax)**2, 1 )/ (nz-2)
case( 50 ) !  the Wangara case
  percent_mean_sqr = & 
    sum( (( fvar1(1:nz)-fvar2(1:nz) )/minmax)**2, 1 ) / nz

case default !
  stop "Not able to handle specified number of HOC z-levels"
end select ! nz

return
end subroutine les_interpolate

end program compare_runs 
!----------------------------------------------------------------------
