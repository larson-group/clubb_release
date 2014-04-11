program sounding_convert
!
!       Simply this program takes sounding information from old
!       <casename>_model.in and converts it to the "SAM-like" format.
!       Casenames of files to be converted go into a file called 'toconvert'.
!       This program must be run in  the directory where the model files are
!       stored.
!
!
!
!--------------------------------------------------------------------------------

! Arbitrary max number of sounding points
  integer, parameter :: nmaxsnd = 600

! Name of the case being converted.
  character(len=200) :: casename

  integer :: i, j
  integer :: nlevels


! Sounding Information
double precision, dimension(nmaxsnd) :: &
      z,      & ! Altitude                               [m]
      theta,  & ! Liquid potential temperature sounding  [K]
      rt,     & ! Total water mixing ratio sounding      [kg/kg]
      u,      & ! u wind sounding                        [m/s] 
      v,      & ! v wind sounding                        [m/s]
      ug,     & ! u geostrophic wind sounding            [m/s]
      vg     ! v geostrophic wind sounding            [m/s]


  integer, parameter :: iunit = 43
  integer, parameter :: ounit = 44

  namelist /sounding/ nlevels, z, theta, rt, u, v, ug, vg

  open(unit = ounit, file="toconvert", status="old")

  do while( .true. ) ! While there are names in toconvert list
    read(ounit, end=55, fmt=*) casename ! Get the casename to convert
    open(unit = iunit, file = trim(casename)//'_model.in', status = 'old') ! Open that file
    read(unit = iunit, nml = sounding) ! Read in its namelist
    close(iunit)

    ! Write out namelist in new format
    open(unit = iunit, file = trim(casename)//'_sounding.in', status = 'new')
    write(iunit, '(1X, 7A12)' ) "z[m]", "theta[K]", "rt[kg\kg]", "u[m\s]" ,"v[m\s]", "ug[m\s]", "vg[m\s]"
    do i=1, nlevels
      write(iunit, *) z(i), theta(i), rt(i), u(i), v(i), ug(i), vg(i)
    end do

    close(iunit)

  end do
  55 close(ounit)

end program sounding_convert
