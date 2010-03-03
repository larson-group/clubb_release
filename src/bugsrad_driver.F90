!-------------------------------------------------------------------------------
! $Id$
module bugsrad_driver

  implicit none

  public :: compute_bugsrad_radiation

  ! Constant parameters
  integer, private, parameter :: &
    nlen = 1, &   ! Length of the total domain
    slen = 1      ! Length of the sub domain

  private ! Default Scope

  contains

  subroutine compute_bugsrad_radiation &
             ( alt, nz, lin_int_buffer,             &
               extend_atmos_range_size,             &
               extend_atmos_bottom_level,           &
               extend_atmos_top_level,              &
               amu0, &
               thlm, rcm, rtm, rsnowm, rim,& 
               cloud_frac, p_in_Pa, p_in_Pam,       &
               exner, rho_zm,                       &
               radht, Frad,                         &
               Frad_SW_up, Frad_LW_up,              &
               Frad_SW_down, Frad_LW_down )          

! Description:
!   Does the necessary operations to interface the CLUBB model with
!   the BUGS radition scheme.

! Grid Layout:
!
! /////////////// Layers from U.S. Standard Atmosphere or sounding   ///////////////
! ///////////////       Dimension: extend_atmos_range_size           ///////////////
!                                  .
!                                  .
! ---------------         Interpolated Layers                        ---------------
! ---------------         Dimension: lin_int_buffer                  ---------------
!                                  .
!                                  .
! ///////////////          Top of CLUBB Grid                         ///////////////
! ///////////////          Dimension: nz                             ///////////////

! References:
! Stevens, et al., (2001) _Journal of Atmospheric Science_, Vol 58, p.3391-3409
! McClatchey, et al., (1972) _Environmental Research Papers_, No. 411, p.94

! Contact for information on BUGSrad (other than this routine)
!   Norm Wood <norm@atmos.colostate.edu>

! All code external to this based on the BUGSrad source from 2004/7/10
!-------------------------------------------------------------------------------

    use constants, only: fstderr, grav, Cp ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use T_in_K_module, only: thlm2T_in_K ! Procedure(s)

    use error_code, only: clubb_at_least_debug_level ! Procedure(s)

    use stats_type, only: stat_update_var ! Procedure(s)
  
    use extend_atmosphere_module, only: &
      extend_atmos_dim, extend_alt, extend_pinmb, & ! Variable(s)
      extend_T_in_K, extend_sp_hmdty, extend_o3l

    use stats_variables, only: zt, zm, rad_zt, rad_zm, l_stats_samp, & ! Variable(s)
      iT_in_K_rad, ircil_rad, io3l_rad, irsnowm_rad, ircm_in_cloud_rad, &
      icloud_frac_rad, iFrad_SW, iFrad_SW_rad, iFrad_LW, iFrad_LW_rad, &
      iradht_rad, iradht_SW, iradht_SW_rad, iradht_LW, iradht_LW_rad, &
      iFrad_SW_up, iFrad_SW_up_rad, iFrad_LW_up, iFrad_LW_up_rad, &
      iFrad_SW_down, iFrad_SW_down_rad, iFrad_LW_down, iFrad_LW_down_rad, &
      l_output_rad_files

    use parameters_radiation, only: &
      sol_const, &
      alvdr, &
      alvdf, &
      alndr, &
      alndf, &
      slr

    implicit none

    intrinsic :: dble, real

    ! Input Variables
    double precision, intent(in) :: &
      amu0  ! Cosine of the solar zenith angle  [-]

    integer, intent(in) :: &
      nz ! Vertical extent;  i.e. nnzp in the grid class

    ! Number of levels to interpolate from the bottom of extend_atmos to the top
    ! of the CLUBB profile, hopefully enough to eliminate cooling spikes, etc.
    integer, intent(in) :: lin_int_buffer, extend_atmos_range_size

    real, intent(in), dimension(nz) :: &
      alt,          & ! Altitudes of the model              [m]
      thlm,         & ! Liquid potential temp.              [K]
      rcm,          & ! Liquid water mixing ratio           [kg/kg]
      rsnowm,       & ! Snow water mixing ratio             [kg/kg]
      rim,          & ! Ice water mixing ratio              [kg/kg]
      rtm,          & ! Total water mixing ratio            [kg/kg]
      rho_zm,       & ! Density                             [kg/m^3]
      cloud_frac,   & ! Cloud fraction                      [-]
      p_in_Pa,      & ! Pressure on the t grid              [Pa]
      p_in_Pam,     & ! Pressure on the m grid              [Pa]
      exner           ! Exner function                      [-]

    
    integer,intent(in) ::&
      extend_atmos_bottom_level, &
      extend_atmos_top_level

    ! Output Variables
    real, intent(out), dimension(nz) :: &
      Frad,         & ! Total radiative flux          [W/m^2]
      Frad_SW_up,   & ! SW radiative upwelling flux   [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux   [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux [W/m^2]
      Frad_LW_down, & ! LW radiative downwelling flux [W/m^2]
      radht           ! Total heating rate            [K/s]

    ! Local Variables
    real, dimension(nz) :: &
      Frad_SW, & ! SW radiative flux [W/m^2]
      Frad_LW, & ! LW radiative flux [W/m^2]
      radht_SW,& ! SW heating rate   [K/s]
      radht_LW   ! LW heating rate   [K/s]

    real, dimension(nz) :: &
      rcm_in_cloud  ! Liquid water mixing ratio in cloud  [kg/kg]

    ! Altered 3 Oct 2005 to be buffer levels higher
    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size) :: &
      T_in_K,& ! Temperature        [K]
      rcil, &  ! Ice mixing ratio    [kg/kg]
      o3l      ! Ozone mixing ratio  [kg/kg]

    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size+1) :: &
      Frad_uLW, & ! LW upwelling flux         [W/m^2]
      Frad_dLW, & ! LW downwelling flux       [W/m^2]
      Frad_uSW, & ! SW upwelling flux         [W/m^2]
      Frad_dSW    ! SW downwelling flux       [W/m^2]

    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size) :: &
      sp_humidity, & ! Specific humidity      [kg/kg]
      pinmb          ! Pressure in millibars  [hPa]

    ! Pressure in millibars for layers (calculated as an average of pinmb)
    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size+1) :: &
      playerinmb ! [hPa]

    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size) :: &
      dpl, &                                  ! Difference in pressure levels       [hPa]
      rsnowm_2d, rcm_in_cloud_2d, cloud_frac_2d ! Two-dimensional copies of the input parameters

    double precision, dimension(nlen,(nz-1)+lin_int_buffer+extend_atmos_range_size) :: &
      radht_SW_2d, & ! SW Radiative heating rate        [W/m^2]
      radht_LW_2d    ! LW Radiative heating rate        [W/m^2]

    double precision, dimension(nlen) :: &
      ts  ! Surface temperature [K]
    
    double precision :: z1_fact, z2_fact, tmp ! Temp storage

    integer :: i, z, z1, z2  ! Loop indices

    integer :: buffer ! The sum of the two buffers

    !character(len=40) :: time_char

    !-------------------------------------------------------------------------------

    buffer = lin_int_buffer + extend_atmos_range_size

    ! Convert to millibars
    pinmb(1,1:(nz-1))  = dble( p_in_Pa(2:nz) / 100.0 ) ! t grid in CLUBB

    playerinmb(1,1:nz) = dble( p_in_Pam / 100.0 ) ! m grid in CLUBB

    ! Determine rcm in cloud
    rcm_in_cloud = rcm / max( cloud_frac, 0.01 )

    ! Convert theta_l to temperature

    T_in_K(1,1:(nz-1)) = thlm2T_in_K( thlm(2:nz), exner(2:nz), rcm(2:nz) )

    ! Derive Specific humidity from rc & rt.
    do z = 2, nz
      if ( rtm(z) < rcm(z) ) then
        sp_humidity(1,z-1) = 0.0d0
        if ( clubb_at_least_debug_level(1) ) then
          write(fstderr,*) "rvm < 0 at ", z, " before BUGSrad, specific humidity set to 0."
        endif
      else
        sp_humidity(1,z-1) &
          = dble( rtm(z) - rcm(z) ) / dble( 1.0+rtm(z) )
      end if
    end do

    ! Setup miscellaneous variables
    
    ! Ozone at < 1 km = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976.
    !   Convert from g to kg.
    o3l(1,1:(nz-1)) = dble( ( 5.4e-5 / rho_zm(1:(nz-1)) ) * 0.001 )

    ! Convert and transpose as needed
    rcil(1,buffer+1:(nz-1)+buffer)            = flip( dble( rim(2:nz) ), nz-1 )
    rsnowm_2d(1,buffer+1:(nz-1)+buffer)       = flip( dble( rsnowm(2:nz) ), nz-1 )
    rcm_in_cloud_2d(1,buffer+1:(nz-1)+buffer) = flip( dble( rcm_in_cloud(2:nz) ), nz-1 )
    cloud_frac_2d(1,buffer+1:(nz-1)+buffer)   = flip( dble( cloud_frac(2:nz) ), nz-1 )

    T_in_K(1,buffer+1:(nz-1)+buffer) = flip( T_in_K(1,1:(nz-1)), nz-1 )

    sp_humidity(1,buffer+1:(nz-1)+buffer) = flip( sp_humidity(1,1:(nz-1)), nz-1 )

    pinmb(1,(buffer+1):(nz-1+buffer))        = flip( pinmb(1,1:(nz-1)), nz-1 )
    playerinmb(1,(buffer+1):(nz-1+buffer+1)) = flip( playerinmb(1,1:nz), nz )

    o3l(1,buffer+1:(nz-1)+buffer) = flip( o3l(1,1:(nz-1)), nz-1 )

    ! Assume these are all zero above the CLUBB profile
    rsnowm_2d(1,1:buffer)       = 0.0d0
    rcil(1,1:buffer)            = 0.0d0
    rcm_in_cloud_2d(1,1:buffer) = 0.0d0
    cloud_frac_2d(1,1:buffer)   = 0.0d0

    if ( alt(nz) > extend_alt(extend_atmos_dim) ) then

      write(fstderr,*) "The CLUBB model grid (for zm levels) contains an ",  &
                       "altitude above the top of the extended atmosphere ",  &
                       "profile."
      write(fstderr,*) "Top of CLUBB model zm grid =", alt(nz), "m."
      write(fstderr,*) "Top of extended atmosphere profile =",  &
                       extend_alt(extend_atmos_dim), "m."
      write(fstderr,*) "Reduce the vertical extent of the CLUBB model grid."
      ! CLUBB zm grid exceeds a 50 km altitude
      stop "compute_bugsrad_radiation: cannot handle this altitude"

    else

    end if

    ! Add the extended atmospheric profile above the linear interpolation
     T_in_K(1,1:extend_atmos_range_size) = &
                flip( extend_T_in_K( extend_atmos_bottom_level:extend_atmos_top_level), &
                      extend_atmos_range_size )

     sp_humidity(1,1:extend_atmos_range_size) = &
                flip( extend_sp_hmdty( extend_atmos_bottom_level:extend_atmos_top_level ), &
                                       extend_atmos_range_size )

     o3l(1,1:extend_atmos_range_size) = &
                flip( extend_o3l( extend_atmos_bottom_level:extend_atmos_top_level ), &
                                  extend_atmos_range_size )

     pinmb(1,1:extend_atmos_range_size) = &
                flip( extend_pinmb( extend_atmos_bottom_level:extend_atmos_top_level ), &
                    extend_atmos_range_size )
    
    ! Do a linear interpolation to produce the levels between the extended
    ! atmospheric levels and the CLUBB levels;
    ! These levels should number the lin_int_buffer parameter
    z1 = buffer + 1
    z2 = extend_atmos_range_size
    do z = buffer, extend_atmos_range_size+1, -1
      z1_fact = dble( z2 - z ) / dble( z2 - z1 )
      z2_fact = dble( z - z1 ) / dble( z2 - z1 )

      T_in_K(1,z) = z1_fact * T_in_K(1,z1) + z2_fact * T_in_K(1,z2)

      sp_humidity(1,z) = z1_fact * sp_humidity(1,z1) + z2_fact * sp_humidity(1,z2)
      o3l(1,z) = z1_fact * o3l(1,z1) + z2_fact * o3l(1,z2)

      pinmb(1,z) = z1_fact * pinmb(1,z1) + z2_fact * pinmb(1,z2)
    end do

    ! Do a linear interpolation to find playerinmb.  Since this interpolation
    ! occurs at levels above the top of the CLUBB model, the CLUBB zt2zm function
    ! or CLUBB weighted averages do not apply.  The variable playerinmb is being
    ! defined on momentum levels above the top of the CLUBB model, which are
    ! being defined here at points half-way inbetween the thermodynamic levels
    ! above the top of the CLUBB model.  Brian Griffin; May 13, 2008.
    playerinmb(1,2:buffer+1) = ( pinmb(1,1:buffer) + pinmb(1,2:buffer+1) ) / 2.

    ! Do a linear extension to find playerinmb at the uppermost standard
    ! atmosphere momentum level.  The grid is evenly-spaced at these points.
    ! Brian Griffin; May 13, 2008.
    tmp = 2. * playerinmb(1,2) - playerinmb(1,3)
    if ( tmp > 0. ) then
      playerinmb(1,1) = tmp
    else ! Assuming a linear extension didn't work
      playerinmb(1,1) = .5 * playerinmb(1,2)
    end if

    ! Calculate the difference in pressure layers (including buffer levels)
    do i = 1, (nz-1)+buffer
      dpl(1,i) = playerinmb(1,i+1) - playerinmb(1,i)
    end do

    ts(1) = T_in_K(1,(nz-1)+buffer)
!  Write a profile for Kurt's driver program for debugging purposes
!  write(time_char ,*) time
!  time_char =adjustl(time_char)
 !open(10, file="profile"//trim(time_char)//"dat")
 !write(10,'(2i4,a10)') nlen, (nz-1)+buffer, "TROPICAL"
 !do i=1, (nz-1)+buffer
  ! write(10,'(i4,9f12.6)') i, pinmb(1,i), playerinmb(1,i),T_in_K(1,i), &
   !sp_humidity(1,i), 100000.0*o3l(1,i), rcm_in_cloud_2d(1,i), &
   !rcil(1,i), cloud_frac_2d(1,i), dpl(1,i)
 !end do
 !write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+buffer), ts(1), amu0
 !close(10)
 

!  print *, "playerinmb = ", playerinmb
!  print *, "sp_humidity = ", sp_humidity

    call bugs_rad( nlen, slen, (nz-1)+buffer, playerinmb,              &
                   pinmb, dpl, T_in_K, sp_humidity,                    &
                   rcm_in_cloud_2d, rcil, rsnowm_2d, o3l,              &
                   ts, amu0, slr, alvdf,                               &
                   alndf, alvdr, alndr, sol_const,                     &
                   dble( grav ), dble( Cp ), radht_SW_2d, radht_LW_2d, &
                   Frad_dSW, Frad_uSW, Frad_dLW, Frad_uLW,             &
                   cloud_frac_2d )

    ! Michael pointed out that this was a temperature tendency, not a theta_l
    ! tendency.  The 2nd line should fix both.  -dschanen 28 July 2006
    radht_SW(2:nz) = real( flip( radht_SW_2d(1,buffer+1:(nz-1)+buffer), nz-1 ) ) &
                     * ( 1.0 / exner(2:nz) )

    radht_LW(2:nz) = real( flip( radht_LW_2d(1,buffer+1:(nz-1)+buffer), nz-1 ) ) &
                     * ( 1.0 / exner(2:nz) )

    ! No radiative heating below ground
    radht_SW(1) = 0.0
    radht_LW(1) = 0.0

    radht = radht_SW + radht_LW

    ! These are on the m grid, and require no adjusting
    Frad_SW_up = real( flip( Frad_uSW(1,buffer+1:nz+buffer), nz ) )

    Frad_LW_up = real( flip( Frad_uLW(1,buffer+1:nz+buffer), nz ) )

    Frad_SW_down = real( flip( Frad_dSW(1,buffer+1:nz+buffer), nz ) )

    Frad_LW_down = real( flip( Frad_dLW(1,buffer+1:nz+buffer), nz ) )

    Frad_SW(1:nz) = Frad_SW_up - Frad_SW_down

    Frad_LW(1:nz) = Frad_LW_up - Frad_LW_down

    Frad(1:nz) = Frad_SW(1:nz) + Frad_LW(1:nz)

    if ( l_stats_samp ) then

      call stat_update_var( iradht_LW, radht_LW, zt )

      call stat_update_var( iradht_SW, radht_SW, zt )

      call stat_update_var( iFrad_SW, Frad_SW, zm )

      call stat_update_var( iFrad_LW, Frad_LW, zm )

      call stat_update_var( iFrad_SW_up, Frad_SW_up, zm )

      call stat_update_var( iFrad_LW_up, Frad_LW_up, zm )

      call stat_update_var( iFrad_SW_down, Frad_SW_down, zm )

      call stat_update_var( iFrad_LW_down, Frad_LW_down, zm )

      if ( l_output_rad_files ) then
        call stat_update_var( iT_in_K_rad, real(T_in_K(1,:)), rad_zt )

        call stat_update_var( ircil_rad, real(rcil(1,:)), rad_zt )

        call stat_update_var( io3l_rad, real(o3l(1,:)), rad_zt )

        call stat_update_var( irsnowm_rad, real(rsnowm_2d(1,:)), rad_zt )

        call stat_update_var( ircm_in_cloud_rad, real(rcm_in_cloud_2d(1,:)), rad_zt )

        call stat_update_var( icloud_frac_rad, real(cloud_frac_2d(1,:)), rad_zt )

        call stat_update_var( iradht_rad, real(radht_SW_2d(1,:) + radht_LW_2d(1,:)), rad_zt )

        call stat_update_var( iradht_LW_rad, real(radht_LW_2d(1,:)), rad_zt )

        call stat_update_var( iradht_SW_rad, real(radht_SW_2d(1,:)), rad_zt )

        call stat_update_var( iFrad_SW_rad, real(Frad_uSW(1,:) - Frad_dSW(1,:)), rad_zm )

        call stat_update_var( iFrad_LW_rad, real(Frad_uLW(1,:) - Frad_dLW(1,:)), rad_zm )

        call stat_update_var( iFrad_SW_up_rad, real(Frad_uSW(1,:)), rad_zm )

        call stat_update_var( iFrad_LW_up_rad, real(Frad_uLW(1,:)), rad_zm )

        call stat_update_var( iFrad_SW_down_rad, real(Frad_dSW(1,:)), rad_zm )

        call stat_update_var( iFrad_LW_down_rad, real(Frad_dLW(1,:)), rad_zm )
      end if ! l_output_rad_files

    end if ! lstats_samp

    return
  end subroutine compute_bugsrad_radiation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  function flip( x, xdim )
! Description:
! Flips a single dimension array (i.e. a vector), so the first element
! becomes the last and vice versa for the whole column.  This is a
! necessary part of the code because BUGSrad and CLUBB store altitudes in
! reverse order
!-------------------------------------------------------------------------------
    implicit none

    ! Input
    integer, intent(in) :: xdim

    double precision, dimension(xdim), intent(in) :: x

    ! Output
    double precision, dimension(xdim) :: flip

    ! Internal
    double precision, dimension(xdim) :: tmp
    integer :: indx

    do indx = 1, xdim, 1
      tmp(indx) = x((xdim+1) - (indx))
    end do

    flip = tmp

    return
  end function flip
!-------------------------------------------------------------------------------

end module bugsrad_driver
