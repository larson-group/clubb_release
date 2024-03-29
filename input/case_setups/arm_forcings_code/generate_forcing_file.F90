!-------------------------------------------------------------------------------
program generate_forcing_file

! Description:
!   Generate a forcing file for the ARM case
!-------------------------------------------------------------------------------
  implicit none

  integer, parameter :: &
    nzmax   = 5, & ! Number of forcing levels
    ntimes = 6    ! Number of times

  real, parameter :: &
    time_initial = 41400., & ! Initial time        [s]
    sec_per_hr   = 3600.,  & ! Seconds per hour
    g_per_kg     = 1000.     ! Grams per kg

  ! From Table A.3 in A.R. Brown, et al.
  real, parameter, dimension(ntimes) ::  & 
    times = (/ 0.0, 3., 6., 9., 12., 14.5/), &  ! hrs
    atheta = (/ 0.000, 0.000,  0.000, -0.080, -0.160, -0.160/), & ! K/h
    rtheta = (/-0.125, 0.000,  0.000,  0.000,  0.000, -0.100/), & ! K/h
    art    = (/ 0.080, 0.020, -0.040, -0.100, -0.160, -0.300/)    ! g/kg/h

  real, dimension(nzmax) :: &
    heights, &     ! [m]
    rtm_forcing, & ! [kg/kg/s]
    thlm_forcing   ! [K/s]

  real :: time, theta_tmp, rad_tmp, rt_tmp, b

  integer :: i, k

  ! --- Begin Code ---

  heights(1)  = 0.   ! Starting point   [m]
  heights(2)  = 999.
  heights(3)  = 1000.
  heights(4)  = 2999.
  heights(5)  = 3000.

  ! Write the header for the forcing file
  write(6,*) "! Forcing file generated by generate_forcing_file for ARM"
  write(6,*) "z[m]        'thlm_f[K/s]'     'rtm_f[kg/kg/s]'   'um_ref[m/s]'   "// &
    "'vm_ref[m/s]'   'um_f[m/s^2]'   'vm_f[m/s^2]'   'w[m/s]'   'ug[m/s]'   'vg[m/s]'"

  do i = 1, ntimes
    ! Convert to MKS units

    time = times(i) * sec_per_hr + time_initial

    theta_tmp = atheta(i) / sec_per_hr
    rad_tmp   = rtheta(i) / sec_per_hr
    rt_tmp    = art(i) / ( sec_per_hr * g_per_kg )

    ! Interpolate with respect to height above 1000m
    ! This is more elaborate than in needs to be
    do k = 1, nzmax
      select case( int( heights(k) ) )
      case ( 0:999 )
        rtm_forcing(k)  = rt_tmp
        thlm_forcing(k) = theta_tmp + rad_tmp

      case ( 1000:2999 )
        b               = 1. - ( heights(k) - 1000. ) / 2000.
        rtm_forcing(k)  = b * rt_tmp
        thlm_forcing(k) = b * ( theta_tmp + rad_tmp )

      case default ! > 3000
        rtm_forcing(k)  = 0.0
        thlm_forcing(k) = 0.0

      end select
    end do ! k=1..nzmax

    ! Write time data
    write(6,'(f10.1,i3)') time, nzmax
    do k = 1, nzmax
      write(6,'(f10.1,9(g12.4))') heights(k), thlm_forcing(k), rtm_forcing(k), &
        -999.9, -999.9, -999.9, -999.9, 0.0, -999.0, -999.0
    end do
  end do ! i = 1, ntimes

end program generate_forcing_file
