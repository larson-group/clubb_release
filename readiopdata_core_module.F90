! Added by Steffen Domke (2019):
! Added module statement in order to avoid conversion errors from real(4) to real
! regarding  daysnd(i) in calcdate
module readiopdata_core_module

!bloss #include <params.h>
!bloss #include <max.h>
!------------------------------------------------------------------------
! File: readiopdata_core_module.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id$
!
!------------------------------------------------------------------------

  use netcdf, only: &
        NF90_OPEN, &
        NF90_NOWRITE, &
        NF90_NOERR, &
        NF90_CLOSE, &
        NF90_INQ_VARID, &
        NF90_INQ_DIMID, &
        NF90_GET_VAR, &
        NF90_MAX_NAME, &
        NF90_MAX_VAR_DIMS, &
        NF90_REAL, &
        NF90_FLOAT, &
        NF90_DOUBLE, &
        NF90_INT, &
        NF90_INQUIRE_VARIABLE, &
        NF90_INQUIRE_DIMENSION

  implicit none

  real, parameter, private :: missing_value = -99999.
  ! Variables used in previously `contained` subroutines
  logical, private :: masterproc_
contains

! Allocates and reads dimensions from a netcdf dataset
!
! masterproc, iopfilepath,      : SAM inputs
! scamiop_from_global_cam       |
! ntime, nlev, nlat, nlon       : Length of dimension
! tsec, dplevs, lat_in, lon_in  : Dimension data
! bdate                         : Start date as 6-digit int
subroutine readiopdata_core_dims( masterproc, iopfilepath, scamiop_from_global_cam, &
                                 ntime, nlev, nlat, nlon, tsec, bdate, dplevs, lat_in, lon_in )
    implicit none

    logical, intent(in) :: masterproc, scamiop_from_global_cam
    character(*), intent(in) :: iopfilepath

    integer, intent(out) :: ntime, nlev, nlat, nlon

    integer, allocatable, intent(out) :: tsec(:)
    integer, intent(out) :: bdate
    real, allocatable, intent(out) :: dplevs(:), lat_in(:), lon_in(:)

    integer :: ncid, status, varID

    logical :: use_NF90_real

   ! Assign module-level variable
   masterproc_ = masterproc

!     
!     Open IOP dataset
!     
   if(masterproc) write(*,*) 'Opening ', iopfilepath
   STATUS = NF90_OPEN( iopfilepath, NF90_NOWRITE, NCID )
   if ( STATUS /= NF90_NOERR ) then
      if(masterproc) write( 6,* ) &
           'ERROR(readiopdata_core_module.f90):Cant open iop dataset: ' ,iopfilepath
      call task_abort() 
   end if

!=====================================================================
!     Read time variables     
!
! Modified by Steffen Domke (2019):
! The SCAM forcing files include both a 'time' dimension and a 'tsec' variable
! So we need to read dimlength from 'time'.
! However, the correct time values are stored in tsec, as 'time' is specified in day units
#ifdef UWM_MISC
    ! Get dimlength from time dimension, usually called 'time' 
    call get_netcdf_dimlength(ncid, 'time', ntime, status, .false.)
    if ( STATUS /= NF90_NOERR )  then
       ! If it is not called 'time', try 'tsec'
       call get_netcdf_dimlength(ncid, 'tsec', ntime, status, .true.)
    end if
    ! If both fail, abort
    if ( status /= NF90_NOERR ) then
       if (masterproc) write(6,*)'Error(readiopdata_core_module.f90):Could not read length of time dimension.'
       status = NF90_CLOSE ( NCID )
       call task_abort()
    end if
    
    ! Get variable ID for time variable. Try 'tsec' first, as we want the time values to be in units of seconds
    STATUS = NF90_INQ_VARID( NCID, 'tsec', varID )
    if ( status /= NF90_NOERR ) then
       ! If failed, try 'time'
       STATUS = NF90_INQ_VARID( NCID, 'time', varID )
    end if
#else
    call get_netcdf_dimlength(ncid, 'time', ntime, status, .false.)
    if ( STATUS /= NF90_NOERR )  then
       call get_netcdf_dimlength(ncid, 'tsec', ntime, status, .true.)
       STATUS = NF90_INQ_VARID( NCID, 'tsec', varID )
    else
       STATUS = NF90_INQ_VARID( NCID, 'time', varID )
       if ( STATUS /= NF90_NOERR ) then
          ! you might end up here if the dimension is called time, but the variable is tsec.
          STATUS = NF90_INQ_VARID( NCID, 'tsec', varID )
       end if
    end if
#endif /*UWM_MISC*/
    
    if ( STATUS /= NF90_NOERR ) then
      write( 6,* )'ERROR(readiopdata_core_module.f90):Could not get variable ID for tsec'
      STATUS = NF90_CLOSE ( NCID )
      call task_abort()
    end if

    ALLOCATE(tsec(ntime),STAT=status)
    if(status/=0) then
       write(6,*) 'Could not allocate tsec in readiopdata'
       call task_abort()
    end if

    STATUS = NF90_GET_VAR( NCID, varID, tsec )
    if ( STATUS /= NF90_NOERR ) then
      write( 6,* )'ERROR(readiopdata_core_module.f90):Could not read variable tsec'
      STATUS = NF90_CLOSE ( NCID )
      call task_abort()
    end if

    ! Load date
    STATUS = NF90_INQ_VARID( NCID, 'nbdate', varID )
    if ( STATUS /= NF90_NOERR ) then
       STATUS = NF90_INQ_VARID( NCID, 'bdate', varID )
       if ( STATUS /= NF90_NOERR ) then
          write( 6,* )'ERROR(readiopdata_core_module.f90):Could not find variable ID for bdate'
          STATUS = NF90_CLOSE ( NCID )
          call task_abort()
       end if
    end if

    STATUS = NF90_GET_VAR( NCID, varID, bdate )
    if ( STATUS /= NF90_NOERR )then
      write( 6,* )'ERROR(readiopdata_core_module.f90):Could not find variable bdate'
      STATUS = NF90_CLOSE ( NCID )
      call task_abort()
    end if

!     
!======================================================
!     read level data
!     
    call get_netcdf_dimlength(ncid, 'lev', nlev, status, .true.)

    ALLOCATE(dplevs(nlev+1),STAT=status)
    if(status/=0) then
       write(6,*) 'Could not allocate dplevs in readiopdata'
       call task_abort()
    end if
    dplevs(:) = 0

    ! get pressure levels (in Pascal)
    call get_netcdf_var1d_real( NCID, 'lev', dplevs, use_NF90_REAL, status,.true. )

#ifdef UWM_MISC
    if ( scamiop_from_global_cam) then
    ! CAM pressure levels are output in mb, convert to Pa
      dplevs = dplevs*100.
    endif
#endif /* UWM_MISC */  

!     
!======================================================
!     read lat/lon data
!     
    call get_netcdf_dimlength(ncid, 'lat', nlat, status, .true.)
    call get_netcdf_dimlength(ncid, 'lon', nlon, status, .true.)

    ALLOCATE(lat_in(nlat),lon_in(nlon),STAT=status)
    if(status/=0) then
       write(6,*) 'Could not allocate lat/lon in readiopdata'
       call task_abort()
    end if

    ! get latitude
    call get_netcdf_var1d_real( NCID, 'lat',lat_in,use_NF90_REAL,status,.false. )

    ! get longitude
    call get_netcdf_var1d_real( NCID, 'lon',lon_in,use_NF90_REAL,status,.false. )

    STATUS = NF90_CLOSE( NCID )

end subroutine readiopdata_core_dims

! Allocates and reads sounding data from netcdf dataset
!
! Inputs:
! masterproc, iopfilepath,      : SAM inputs
! scamiop_from_global_cam       |
! ntime, nlev, dplevs           : Dimension info
! rgas, cp, fac_cond, fac_sub   : Required constants
!
! Outputs:
! zsnd                          : Empty sounding for z
! usnd, vsnd, ugls, vgls        : Soundings for u, v, ug, vg
! tsnd, qsnd, psnd, wgls        : Soundings for T, q, p, omega
! get_add_surface_data,         : Loading flags
! have_omega, wgls_holds_omega  |
subroutine readiopdata_core_snd( masterproc, iopfilepath, scamiop_from_global_cam, &
                               ntime, nlev, dplevs, rgas, cp, fac_cond, fac_sub, &
                               zsnd, usnd, vsnd, tsnd, qsnd, psnd, ugls, vgls, wgls, &
                               get_add_surface_data, have_omega, wgls_holds_omega &
#ifdef UWM_MISC
                               , have_w, wgls_holds_w, subs_type &
#endif
                               )
!-----------------------------------------------------------------------
!     
!     Open and read netCDF file containing initial IOP  conditions
!     
!---------------------------Code history--------------------------------
!     
!     Written by J.  Truesdale    August, 1996, revised January, 1998
!     
!     Modified by Peter Blossey (pblossey@u.washington.edu) July 2006
!       Reorganization and some changes that enable SAM to read 
!       initial soundings and forcings from SCAM IOP input datasets.
!
!     Modified by Peter Blossey (pblossey@u.washington.edu) March 2008
!       With SAM6.5, sounding and forcing data no longer need to be
!       interpolated to SAM vertical grid.  These changes should not
!       make much effect on results, but seem to clean up the code a bit.
!
!     Modified by Peter Blossey (pblossey@u.washington.edu) April 2008
!       For release in SAM6.7, have made further changes:
!         - added logical get_add_surface_data to eliminate
!             the appending of surface data to soundings/forcings.
!             If there is surface data in the netcdf file, it is read
!             in.  Otherwise, it is extrapolated.
!
!         - added input of CLDLIQ and CLDICE variables, in case cloud is
!              explicitly included in SCAM IOP initial condition.  Note
!              that the cloud liquid and ice is immediately added to water
!              vapor and the temperature is modified to account for the
!              latent heat release.  The cloud will arise with the
!              initial calls to the microphysics.
!
!     ===================================================================
!     NOTE: If the initial sounding has cloud, I (Peter) recommend using:
!             - Tli = Tabs - (L/Cp)*ql - (Ls/Cp)*qi in place of Tabs, and
!             - qtot = qv + ql + qi in place of qv
!       in the netcdf SCAM input file.  SAM will perform a saturation
!       adjustment of the initial profile that will convert the excess
!       vapor to cloud in the initial condition and give the correct
!       initial absolute temperature profile.
!     ===================================================================
!
!-----------------------------------------------------------------------

   implicit none

   ! Outputs
   real, allocatable, intent(out) :: usnd(:,:), vsnd(:,:), tsnd(:,:), qsnd(:,:), &
        psnd(:,:), zsnd(:,:), ugls(:,:), vgls(:,:), wgls(:,:)

#ifdef UWM_MISC
   character(len=*), intent(out) :: subs_type

   ! Optional outputs
   logical, optional, intent(out) :: have_w, wgls_holds_w
#endif

   logical, optional, intent(out) :: get_add_surface_data, have_omega, wgls_holds_omega

   ! Inputs
   integer, intent(in) :: ntime, nlev
   real, intent(in) :: dplevs(:)
   logical, intent(in) :: masterproc, scamiop_from_global_cam
   character(*), intent(in) :: iopfilepath
   real, intent(in) :: fac_cond, fac_sub, rgas, cp

   ! Local variables
   integer :: nlev_in, nsnd, nzsnd, nlsf, nzlsf
   integer :: STATUS, i, n, nlat, nlon
   integer :: varID, NCID

   real :: coef
   real, allocatable :: tmp_srf(:), Ts_in(:), T_in(:,:), q_in(:,:), u_in(:,:), v_in(:,:), &
        ug_in(:,:), vg_in(:,:), omega_in(:,:), cldliq_in(:,:), cldice_in(:,:), Tg_in(:), Ps_in(:)
#ifdef UWM_MISC
   real, allocatable :: w_in(:,:)
#endif
   real :: levs(nlev+1)

   logical :: use_NF90_real, have_tsair, have_tg, have_geostrophic_wind

#ifdef UWM_MISC
   ! This is manually defined since sam_clubb doesn't have input_names file while clubb does
   character(len=6)  :: wm_name = 'w[m/s]'
   character(len=11) :: omega_name = 'omega[Pa/s]'
#endif

   ! Assign module-level variable
   masterproc_ = masterproc

!     
!     Open IOP dataset
!     
   if(masterproc) write(*,*) 'Opening ', iopfilepath
   STATUS = NF90_OPEN( iopfilepath, NF90_NOWRITE, NCID )
   if ( STATUS /= NF90_NOERR ) then
      if(masterproc) write( 6,* ) &
           'ERROR(readiopdata_core_module.f90):Cant open iop dataset: ' ,iopfilepath
      call task_abort() 
   end if

!
!======================================================
!     allocate surface variables
!     
   ALLOCATE(Tg_in(ntime), Ts_in(ntime), Ps_in(ntime), tmp_srf(ntime), STAT=status)

   Tg_in(:) = missing_value
   Ts_in(:) = missing_value
   Ps_in(:) = missing_value

   if(status/=0) then
      write(6,*) 'Could not allocate surface variables in readiopdata'
      call task_abort()
   end if

!
!======================================================
!     read surface variables
!     
! note that the last argument is whether the run should die if the variable
!   is not present in the netcdf file.
!
   ! surface air temperature
   call get_netcdf_var1d_real( ncid, 'Tsair', Ts_in,use_NF90_REAL,status,.false.)
   have_tsair = .true.
   if ( STATUS /= NF90_NOERR ) have_tsair = .false.

#ifdef UWM_MISC
!In sam_clubb:ticket:87, we noticed that Tsair is present for our supplied
!CAM-IOP file, however it was all zero.
!
!Check to see if any zeros (unphysical) are present. If so, ignore Tsair
   if( scamiop_from_global_cam ) then
     if( ANY(Ts_in == 0.0) ) then
       have_tsair = .false.
     endif
   endif
#endif
   ! ground/sea surface temperature
   call get_netcdf_var1d_real( ncid, 'Tg', Tg_in, use_NF90_REAL, status,.false.)
   have_tg = .true.
   if ( STATUS /= NF90_NOERR ) have_tg = .false.

   ! surface pressure
   call get_netcdf_var1d_real( ncid, 'Ps', Ps_in, use_NF90_REAL, status, .true.)

!         
!====================================================================
!     check whether surface pressure exceeds largest pressure
!       in pressure sounding (dplevs)
!     
   if(MINVAL(Ps_in)<=MAXVAL(dplevs(1:nlev))) then
      ! Surface pressure is included in dplevs.
      ! Do not bother with adding surface data to soundings/forcings.
      get_add_surface_data = .false.
      ! do not leave room for surface data in input datasets
      nlev_in = nlev
   else
      ! Surface pressure exceeds max pressure in dplevs
      ! Get/add surface data to soundings/forcings.
      get_add_surface_data = .true.
      ! leave room for surface data in input datasets
      nlev_in = nlev+1
   end if

   !convert pressures to millibar
   levs(1:nlev) = dplevs(1:nlev)/100.
   Ps_in(:) = Ps_in(:)/100. ! convert to millibar

!         
!====================================================================
!     allocate variables with pressure and time dependence (q,T,etc.)
!     
   Allocate(T_in(ntime,nlev_in), q_in(ntime,nlev_in), &
        u_in(ntime,nlev_in), v_in(ntime,nlev_in), &
        ug_in(ntime,nlev_in), vg_in(ntime,nlev_in), &
        omega_in(ntime,nlev_in), cldliq_in(ntime,nlev_in), &
        cldice_in(ntime,nlev_in), &
#ifdef UWM_MISC
        w_in(ntime,nlev_in), &
#endif
        STAT=status)

   if(status/=0) then
      write(6,*) 'Could not allocate surface variables in readiopdata'
      call task_abort()
   end if

!
!====================================================================
!     read variables with pressure and time dependence (q,T,etc.)
!     
   ! Temperature
   T_in(:,:) = missing_value
#ifdef UWM_MISC
   if ( scamiop_from_global_cam ) then
   ! CAM outputs temperature as 't' instead of 'T'
     call get_netcdf_var2d_real( ncid,'t',ntime,nlev,T_in, &
                                 use_NF90_REAL,status,.true.)
   else
#endif
   call get_netcdf_var2d_real( ncid,'T',ntime,nlev,T_in, &
        use_NF90_REAL,status,.true.)
#ifdef UWM_MISC
   endif
#endif

   !==================
   ! Moisture
   q_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'q',ntime,nlev,q_in,use_NF90_REAL,status,.true.)

   cldliq_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'CLDLIQ',ntime,nlev,cldliq_in,use_NF90_REAL,status,.false.)
   if( STATUS == NF90_NOERR ) then
      ! if CLDLIQ is present, add cloud liquid water to water vapor and
      !   modify initial temperature to reflect release of latent heat.
      !   SAM does not support a specified initial cloud layer at this
      !   point.  The initial cloud will arise with the first saturation adjustment.
      q_in(:,1:nlev) = q_in(:,1:nlev) + cldliq_in(:,1:nlev)
      t_in(:,1:nlev) = t_in(:,1:nlev) - fac_cond*cldliq_in(:,1:nlev)
   end if

   cldice_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'CLDICE',ntime,nlev,cldice_in,use_NF90_REAL,status,.false.)
   if( STATUS == NF90_NOERR ) then
      ! if CLDICE is present, add cloud ice water to water vapor 
      !   modify initial temperature to reflect release of latent heat
      !   SAM does not support a specified initial cloud layer at this
      !   point.  The initial cloud will arise with the first saturation adjustment.
      q_in(:,1:nlev) = q_in(:,1:nlev) + cldice_in(:,1:nlev)
      t_in(:,1:nlev) = t_in(:,1:nlev) - fac_sub*cldice_in(:,1:nlev)
   end if

   !==================================
   ! omega or w: vertical pressure velocity
   
   omega_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'omega',ntime,nlev,omega_in, &
         use_NF90_REAL,status,.false.)
   have_omega = .true.
   if( STATUS /= NF90_NOERR ) have_omega = .false.

#ifdef UWM_MISC
   if( .not. have_omega ) then
      ! If omega fails, try to get w (vertical velocity) instead
      w_in(:,:) = missing_value
      call get_netcdf_var2d_real( ncid,'w',ntime,nlev,w_in, &
         use_NF90_REAL,status,.false.)
      have_w = .true.
      if( STATUS /= NF90_NOERR ) have_w = .false.
   end if
   
   if( have_omega ) then
      subs_type = omega_name
   else if( have_w ) then
      subs_type = wm_name
   end if


   if( scamiop_from_global_cam .and. ( have_omega .or. have_w ) ) then
   ! In sam_clubb:ticket:87, we noticed that SAM is imposing subsidence since
   ! the CAM IOP files have omega. Since there is no horizontal (only)
   ! divergence term for temperature and moisture in the IOP file and SAM is
   ! computing ther vertical component- the large scale horizontal advective
   ! tendencies are being ignored.

   ! This is a fix until we can generate IOP files on our own. Therefore, we
   ! will set omega to zero and use the IOP supplied horiz + vert tendencies of
   ! temperature and moisture. weberjk(UWM)
     omega_in(:,:) = 0.
     w_in(:,:) = 0.
     have_omega = .false.
     have_w = .false.
   endif
#endif /* UWM_MISC */

   !==================================
   ! horizontal wind
   u_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'u',ntime,nlev,u_in, &
        use_NF90_REAL,status,.true.)

   v_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'v',ntime,nlev,v_in, &
        use_NF90_REAL,status,.true.)

   !==================================
   ! geostrophic wind (not native to SCAM, but useful in SAM)
   ug_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'ug',ntime,nlev,ug_in, &
        use_NF90_REAL,status,.false.)
   have_geostrophic_wind = .true.
   if( STATUS /= NF90_NOERR ) have_geostrophic_wind = .false.

   vg_in(:,:) = missing_value
   call get_netcdf_var2d_real( ncid,'vg',ntime,nlev,vg_in, &
        use_NF90_REAL,status,.false.)
   if( STATUS /= NF90_NOERR ) have_geostrophic_wind = .false.

   !==================================
   !==================================
   ! READ IN SURFACE DATA AND PUT INTO FORCINGS/SOUNDINGS IF SURFACE
   !   PRESSURE IS BIGGER THAN MAX PRESSURE IN SOUNDING.
   if(get_add_surface_data) then
      ! temperature
      if(have_tsair) T_in(:,nlev+1) = Ts_in(:)

      ! surface pressure tendency --> surface omega
      call get_netcdf_var1d_real(ncid,'Ptend',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) then
#ifdef UWM_MISC
         if( have_w ) then
            w_in(:,nlev+1) = tmp_srf(:)
         end if

         if( have_omega ) then
            omega_in(:,nlev+1) = tmp_srf(:) 
         end if
#else
         omega_in(:,nlev+1) = tmp_srf(:) 
#endif
      end if

      ! winds
      call get_netcdf_var1d_real(ncid,'usrf',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) then
         u_in(:,nlev+1) = tmp_srf(:)
      end if

      call get_netcdf_var1d_real(ncid,'vsrf',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) then
         v_in(:,nlev+1) = tmp_srf(:)
      end if

      ! moisture
      call get_netcdf_var1d_real( ncid,'qsrf',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) q_in(:,nlev+1) = tmp_srf(:)

      call get_netcdf_var1d_real(ncid,'ugsrf',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) ug_in(:,nlev+1) = tmp_srf(:)

      call get_netcdf_var1d_real(ncid,'vgsrf',tmp_srf,use_NF90_REAL,status,.false.)
      if( STATUS == NF90_NOERR ) vg_in(:,nlev+1) = tmp_srf(:)
   endif

   !=========================================
   ! fix ground temperature if not in nc file
   if(.not.have_tg) then
      if ( have_tsair ) then
         write(6,*) 'In readiopdata(): Using Tsair for Tground'
         Tg_in = Ts_in(:)
      else
         write(6,*) 'In readiopdata(): Using lowest level T for Tground'
         Tg_in = T_in(:,nlev)
      end if
   end if


   STATUS = NF90_CLOSE( NCID )
   !===========================================================================
   ! END LOADING VALUES
   !===========================================================================

   !bloss (10Mar2008): Modify to handle SAM's new forcing setup (as of ~6.5).
   !      This means that sounding and forcing data can sit on its own
   !      grid, rather than needing to be interpolated to the model grid.
   nsnd = ntime
   nzsnd = nlev_in
   nlsf = ntime
   nzlsf = nlev_in
   if(masterproc) print*,'sounding data: nsnd=',nsnd,'  nzsnd=',nzsnd
   if(masterproc)print*,'forcing data: nlsf=',nlsf,'  nzlsf=',nzlsf
   allocate(zsnd(nzsnd,nsnd),usnd(nzsnd,nsnd),vsnd(nzsnd,nsnd), &
        tsnd(nzsnd,nsnd),qsnd(nzsnd,nsnd), &
        psnd(nzsnd,nsnd), &
        ugls(nzlsf,nlsf),vgls(nzlsf,nlsf), wgls(nzlsf,nlsf), &
        STAT=STATUS)
   if(STATUS/=0) then
      if(masterproc) then
         write(*,*) 'Error in allocating snd/lsf/rad/sfc vars in readiopdata'
      end if
      call task_abort()
   end if

   ! change absolute temperature sounding to potential temperature
   do i = 1,ntime
      do n = 1,nlev
         T_in(i,n) = T_in(i,n)*(1000./levs(n))**(rgas/cp) ! T --> theta
      end do
   end do

   if(get_add_surface_data) then
      ! either use surface data from netcdf file or interpolate it to surface.

      ! fix surface theta
      do i = 1,ntime
         if(T_in(i,nlev+1)/=missing_value) then
            ! set iop dataset surface pressure for this time
            n = nlev+1
            levs(n) = Ps_in(i) 
            T_in(i,n) = T_in(i,n)*(1000./levs(n))**(rgas/cp) ! T --> theta
         end if
      end do

      ! interpolate to fill surface data/forcings if necessary
      do i = 1,ntime

         ! extrapolate to surface if iop data at surface is not available.
         levs(nlev+1) = Ps_in(i) 
         coef = (levs(nlev+1)-levs(nlev-1))/(levs(nlev)-levs(nlev-1))

         ! fill in surface sounding data if missing
         if(T_in(i,nlev+1)==missing_value) &
              T_in(i,nlev+1) = (1-coef)*T_in(i,nlev-1) + coef*T_in(i,nlev)
         if(q_in(i,nlev+1)==missing_value) &
              q_in(i,nlev+1) = (1-coef)*q_in(i,nlev-1) + coef*q_in(i,nlev)
         if(u_in(i,nlev+1)==missing_value) &
              u_in(i,nlev+1) = (1-coef)*u_in(i,nlev-1) + coef*u_in(i,nlev)
         if(v_in(i,nlev+1)==missing_value) &
              v_in(i,nlev+1) = (1-coef)*v_in(i,nlev-1) + coef*v_in(i,nlev)

#ifdef UWM_MISC
         if(have_omega) then
            if(omega_in(i,nlev+1)==missing_value) omega_in(i,nlev+1) = &
                 (1-coef)*omega_in(i,nlev-1) + coef*omega_in(i,nlev)
         end if
         if(have_w) then
            if(w_in(i,nlev+1)==missing_value) w_in(i,nlev+1) = &
                 (1-coef)*w_in(i,nlev-1) + coef*w_in(i,nlev)
         end if
#else
         if(omega_in(i,nlev+1)==missing_value) omega_in(i,nlev+1) = &
              (1-coef)*omega_in(i,nlev-1) + coef*omega_in(i,nlev)
#endif
         if(ug_in(i,nlev+1)==missing_value) &
              ug_in(i,nlev+1) = (1-coef)*ug_in(i,nlev-1) + coef*ug_in(i,nlev)
         if(vg_in(i,nlev+1)==missing_value) &
              vg_in(i,nlev+1) = (1-coef)*vg_in(i,nlev-1) + coef*vg_in(i,nlev)
      end do

   end if !if(get_add_surface_data)

   ! now fill in profiles for snd
   zsnd(:,:) = -999.
   do i = 1,nsnd
      psnd(1,i) = Ps_in(i)
      psnd(2:nzsnd,i) = levs(nzsnd-1:1:-1)
      tsnd(1:nzsnd,i) = T_in(i,nzsnd:1:-1)
      qsnd(1:nzsnd,i) = q_in(i,nzsnd:1:-1)
      usnd(1:nzsnd,i) = u_in(i,nzsnd:1:-1) 
      vsnd(1:nzsnd,i) = v_in(i,nzsnd:1:-1) 
   end do

   ! convert qsnd to g/kg to be consistent with Marat's implementation
   qsnd(:,:) = 1.e3*qsnd(:,:)

   wgls(:,:) = 0.

   if(have_omega) then
      ! use omega for large-scale vertical advection if it exists
      do i = 1,nlsf
         ! NOTE: Here, we are putting omega into wgls (large-scale w)
         ! THIS WILL BE CONVERTED INTO w IN forcing()
         wgls(1:nzlsf,i) = omega_in(i,nzlsf:1:-1)
      end do
      wgls_holds_omega = .true.
#ifdef UWM_MISC
   else if(have_w) then
      ! use w for large-scale vertical advection if it exists
      do i = 1,nlsf
         wgls(1:nzlsf,i) = w_in(i,nzlsf:1:-1)
      end do
      wgls_holds_w = .true.
   else
      ! no large-scale vertical advection is present in the dataset
      wgls_holds_omega = .false.
      wgls_holds_w = .false.
      write(6,*) 'No large-scale vertical advection is present in the dataset'
      call task_abort()
#endif
   end if

   if(have_geostrophic_wind) then
      do i = 1,nlsf
         ugls(1:nzlsf,i) = ug_in(i,nzlsf:1:-1)
         vgls(1:nzlsf,i) = vg_in(i,nzlsf:1:-1)
      end do
   else
      do i = 1,nlsf ! use wind sounding as geostrophic wind
         ugls(1:nzlsf,i) = u_in(i,nzlsf:1:-1)
         vgls(1:nzlsf,i) = v_in(i,nzlsf:1:-1)
      end do
   end if

   deallocate(Tg_in, Ps_in, Ts_in, STAT=STATUS)
   if(status/=0) then
      write(6,*) 'Could not de-allocate surface data arrays in readiopdata'
      call task_abort()
   end if

   deallocate(T_in, q_in, u_in, v_in, ug_in, vg_in, omega_in, cldliq_in, cldice_in, &
#ifdef UWM_MISC
               w_in, &
#endif
               STAT=STATUS)

   if(status/=0) then
      write(6,*) 'Could not de-allocate sounding/forcing arrays in readiopdata'
      call task_abort()
   end if

   return
end subroutine readiopdata_core_snd

! Allocate and read surface data
!
! Inputs:
! masterproc, iopfilepath,      : SAM inputs
! scamiop_from_global_cam,      |
! SFC_FLX_FXD,                  |
! get_add_surface_data          |
! ntime                         : Number of timesteps
!
! Outputs:
! sstsfc,                       : Surface temperature
! shsfc, lhsfc, tausfc          : SAM outputs
subroutine readiopdata_core_sfc( masterproc, iopfilepath, scamiop_from_global_cam, SFC_FLX_FXD, &
                               get_add_surface_data, ntime, &
                               sstsfc, shsfc, lhsfc, tausfc )
    implicit none

    logical, intent(in) :: masterproc, scamiop_from_global_cam, SFC_FLX_FXD, get_add_surface_data
    character(*), intent(in) :: iopfilepath
    integer, intent(in) :: ntime

    real, allocatable, intent(out) :: sstsfc(:), shsfc(:), lhsfc(:), tausfc(:)

    real, allocatable :: shf_in(:), lhf_in(:), tausrf_in(:), tmp_srf(:), Tg_in(:)
    logical :: have_shflx, have_lhflx, have_tausrf, use_nf90_real
    integer :: ncid, status, nsfc, nlsf, i

    ! Assign module-level variable
    masterproc_ = masterproc

!     
!     Open IOP dataset
!     
   if(masterproc) write(*,*) 'Opening ', iopfilepath
   STATUS = NF90_OPEN( iopfilepath, NF90_NOWRITE, NCID )
   if ( STATUS /= NF90_NOERR ) then
      if(masterproc) write( 6,* ) &
           'ERROR(readiopdata_core_module.f90):Cant open iop dataset: ' ,iopfilepath
      call task_abort() 
   end if

!
!======================================================
!     allocate surface variables
!     
   ALLOCATE(Tg_in(ntime), shf_in(ntime), lhf_in(ntime), tmp_srf(ntime), tausrf_in(ntime), &
        STAT=status)

   shf_in(:) = missing_value
   lhf_in(:) = missing_value
   tausrf_in(:) = missing_value

   if(status/=0) then
      write(6,*) 'Could not allocate surface variables in readiopdata'
      call task_abort()
   end if

!
!======================================================
!     read surface variables
!     
! note that the last argument is whether the run should die if the variable
!   is not present in the netcdf file.
!
   ! sensible heat flux
   call get_netcdf_var1d_real( ncid, 'shflx', shf_in, use_NF90_REAL, &
        status,.false.)
   have_shflx = .true.
   if ( STATUS /= NF90_NOERR ) then
      ! old name - backwards compatibility
      call get_netcdf_var1d_real( ncid, 'sh', shf_in, use_NF90_REAL, &
           status,.false.)
      if ( STATUS /= NF90_NOERR ) have_shflx = .false.
   end if

   ! latent heat flux
   call get_netcdf_var1d_real( ncid, 'lhflx', lhf_in, use_NF90_REAL, &
        status,.false.)
   have_lhflx = .true.
   if ( STATUS /= NF90_NOERR ) then
      ! old name - backwards compatibility
      call get_netcdf_var1d_real( ncid, 'lh', lhf_in, use_NF90_REAL, &
           status,.false.)
      if ( STATUS /= NF90_NOERR ) have_lhflx = .false.
   end if

   ! abort if surface fluxes are required, but are not present in netcdf file
   if(SFC_FLX_FXD.and.(.NOT.have_lhflx.OR..NOT.have_shflx)) then
      if(masterproc) then
         write(*,*) 'ERROR(readiopdata_core_module.f90): If SFC_FLX_FXD is true, '
         write(*,*) '          shflx and lhflx needed in SCAM iop netcdf file.'
      end if
      call task_abort()
   end if

   ! ground/sea surface temperature
   call get_netcdf_var1d_real( ncid, 'Tg', Tg_in, use_NF90_REAL, status,.false.)


   ! winds
   call get_netcdf_var1d_real(ncid,'usrf',tmp_srf,use_NF90_REAL,status,.false.)
   have_tausrf = .false.
   if( STATUS == NF90_NOERR ) then
      tausrf_in(:) = tmp_srf(:)**2
      have_tausrf = .true.
   end if

   call get_netcdf_var1d_real(ncid,'vsrf',tmp_srf,use_NF90_REAL,status,.false.)
   if( STATUS == NF90_NOERR ) then
      if(have_tausrf) then
         tausrf_in(:) = tausrf_in(:) + tmp_srf(:)**2
      end if
   else
      have_tausrf = .false.
   end if

   status = nf90_close( ncid )

   nsfc = ntime
   nlsf = ntime
   if(masterproc)print*,'surface forcing data: nsfc=',nsfc
   allocate(sstsfc(nsfc),shsfc(nsfc),lhsfc(nsfc), tausfc(nsfc), STAT=status)
   if(status/=0) then
      if(masterproc) then
         write(*,*) 'Error in allocating snd/lsf/rad/sfc vars in readiopdata'
      end if
      call task_abort()
   end if

   !set up sfc stuff (surface forcings)
   do i = 1,ntime
      sstsfc(i) = Tg_in(i)
      shsfc(i)   = shf_in(i)
      lhsfc(i)  = lhf_in(i)
      if(have_tausrf) then
         tausfc(i) = sqrt(tausrf_in(i)) !!!!! FIX THIS !!!!!
      else
         tausfc(i) = 0.
      end if
   end do

   deallocate(Tg_in, shf_in,lhf_in,tmp_srf,STAT=status)
   if(status/=0) then
      write(6,*) 'Could not de-allocate surface data arrays in readiopdata'
      call task_abort()
   end if


end subroutine readiopdata_core_sfc

subroutine readiopdata_core_frc( masterproc, iopfilepath, get_add_surface_data, have_omega, &
                               dplevs, nlev, ntime, &
                               dtls, dqls )
    implicit none

    ! Inputs
    logical, intent(in) :: masterproc, get_add_surface_data, have_omega
    character(*), intent(in) :: iopfilepath
    real, intent(in) :: dplevs(:)
    integer, intent(in) :: nlev, ntime

    ! Outputs
    real, allocatable, intent(out) :: dtls(:,:), dqls(:,:)

    ! Local variables
    integer :: ncid, status, nzlsf, nlsf, i, nlev_in
    real :: coef
    real :: levs(nlev+1)
    logical :: use_nf90_real
    ! soundings, omega, advective tendencies (only function of time, lev here)
    real, allocatable :: divT_in(:,:), vertdivT_in(:,:), divT3d_in(:,:), &
        divq_in(:,:), vertdivq_in(:,:), divq3d_in(:,:), tmp_srf(:)
    logical :: have_divq, have_vertdivq, have_divq3d, have_divT, have_vertdivT, have_divT3d

    ! Assign module-level variable
    masterproc_ = masterproc

    if(get_add_surface_data) then
        nlev_in = nlev+1
    else
        nlev_in = nlev
    endif

    levs(1:nlev) = dplevs(1:nlev)/100.

!     
!     Open IOP dataset
!     
    if(masterproc) write(*,*) 'Opening ', iopfilepath
    STATUS = NF90_OPEN( iopfilepath, NF90_NOWRITE, NCID )
    if ( STATUS /= NF90_NOERR ) then
       if(masterproc) write( 6,* ) &
            'ERROR(readiopdata_core_module.f90):Cant open iop dataset: ' ,iopfilepath
       call task_abort() 
    end if

    allocate(tmp_srf(ntime))

!         
!====================================================================
!     allocate variables with pressure and time dependence (q,T,etc.)
!     
    Allocate(divT_in(ntime,nlev_in), divq_in(ntime,nlev_in), &
         divT3d_in(ntime,nlev_in), divq3d_in(ntime,nlev_in), &
         vertdivT_in(ntime,nlev_in), vertdivq_in(ntime,nlev_in), &
         STAT=status)
 
    if(status/=0) then
       write(6,*) 'Could not allocate surface variables in readiopdata'
       call task_abort()
    end if
!
!====================================================================
!     read variables with pressure and time dependence (q,T,etc.)
!     
    ! Horizontal Advective Temperature Forcing
    divT_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'divT',ntime,nlev,divT_in, &
         use_NF90_REAL,status,.false.)
    have_divT = .true.
    if( STATUS /= NF90_NOERR ) have_divT = .false.

    ! Vertical Advective Temperature Forcing
    vertdivT_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'vertdivT',ntime,nlev,vertdivT_in, &
         use_NF90_REAL,status,.false.)
    have_vertdivT = .true.
    if( STATUS /= NF90_NOERR ) have_vertdivT = .false.

    ! Three-dimensional Advective Temperature Forcing
    divT3d_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'divT3d',ntime,nlev,divT3d_in, &
         use_NF90_REAL,status,.false.)
    have_divT3d = .true.
    if( STATUS /= NF90_NOERR ) have_divT3d = .false.

    ! Horizontal Advective Moisture Forcing
    divq_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'divq',ntime,nlev,divq_in, &
         use_NF90_REAL,status,.false.)
    have_divq = .true.
    if( STATUS /= NF90_NOERR ) have_divq = .false.

    ! Vertical Advective Moisture Forcing
    vertdivq_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'vertdivq',ntime,nlev,vertdivq_in, &
         use_NF90_REAL,status,.false.)
    have_vertdivq = .true.
    if( STATUS /= NF90_NOERR ) have_vertdivq = .false.

    ! Three-dimensional Advective Moisture Forcing
    divq3d_in(:,:) = missing_value
    call get_netcdf_var2d_real( ncid,'divq3d',ntime,nlev,divq3d_in, &
         use_NF90_REAL,status,.false.)
    have_divq3d = .true.
    if( STATUS /= NF90_NOERR ) have_divq3d = .false.

    !==================================
    !==================================
    ! READ IN SURFACE DATA AND PUT INTO FORCINGS/SOUNDINGS IF SURFACE
    !   PRESSURE IS BIGGER THAN MAX PRESSURE IN SOUNDING.
    if(get_add_surface_data) then
       ! temperature forcing
       call get_netcdf_var1d_real( ncid,'divTsrf',tmp_srf, &
            use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) divT_in(:,nlev+1) = tmp_srf(:)

       call get_netcdf_var1d_real(ncid,'vertdivTsrf', &
            tmp_srf,use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) vertdivT_in(:,nlev+1) = tmp_srf(:)

       call get_netcdf_var1d_real(ncid,'divT3dsrf', &
            tmp_srf,use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) divT3d_in(:,nlev+1) = tmp_srf(:)

       ! moisture forcing
       call get_netcdf_var1d_real( ncid,'divqsrf',tmp_srf, &
            use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) divq_in(:,nlev+1) = tmp_srf(:)

       call get_netcdf_var1d_real(ncid,'vertdivqsrf', &
            tmp_srf,use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) vertdivq_in(:,nlev+1) = tmp_srf(:)

       ! get surface data if available
       call get_netcdf_var1d_real(ncid,'divq3dsrf', &
            tmp_srf,use_NF90_REAL,status,.false.)
       if( STATUS == NF90_NOERR ) divq3d_in(:,nlev+1) = tmp_srf(:)
    end if ! if(get_add_surface_data)

    status = nf90_close( ncid )

    if(get_add_surface_data) then
       ! either use surface data from netcdf file or interpolate it to surface.
 
       ! interpolate to fill surface data/forcings if necessary
        do i = 1,ntime
            coef = (levs(nlev+1)-levs(nlev-1))/(levs(nlev)-levs(nlev-1))

            ! fill in surface large-scale forcing data if missing
            if(divT_in(i,nlev+1)==missing_value) divT_in(i,nlev+1) = &
                 (1-coef)*divT_in(i,nlev-1) + coef*divT_in(i,nlev)
            if(vertdivT_in(i,nlev+1)==missing_value) vertdivT_in(i,nlev+1) = &
                 (1-coef)*vertdivT_in(i,nlev-1) + coef*vertdivT_in(i,nlev)
            if(divT3d_in(i,nlev+1)==missing_value) divT3d_in(i,nlev+1) = &
                 (1-coef)*divT3d_in(i,nlev-1) + coef*divT3d_in(i,nlev)
 
            if(divq_in(i,nlev+1)==missing_value) divq_in(i,nlev+1) = &
                 (1-coef)*divq_in(i,nlev-1) + coef*divq_in(i,nlev)
            if(vertdivq_in(i,nlev+1)==missing_value) vertdivq_in(i,nlev+1) = &
                 (1-coef)*vertdivq_in(i,nlev-1) + coef*vertdivq_in(i,nlev)
            if(divq3d_in(i,nlev+1)==missing_value) divq3d_in(i,nlev+1) = &
                 (1-coef)*divq3d_in(i,nlev-1) + coef*divq3d_in(i,nlev)
        end do
 
    end if !if(get_add_surface_data)


    nzlsf = nlev
    nlsf = ntime
    allocate(dtls(nzlsf,nlsf),dqls(nzlsf,nlsf), stat=STATUS)

    dtls(:,:) = 0.
    dqls(:,:) = 0.
 
    if(have_omega) then
       ! use large-scale horizontal advection w/omega if it exists.
       if(have_divT) then
          do i = 1,nlsf
             dtls(1:nzlsf,i) = divT_in(i,nzlsf:1:-1)
          end do
       end if
       if(have_divq) then
          do i = 1,nlsf
             dqls(1:nzlsf,i) = divq_in(i,nzlsf:1:-1)
          end do
       end if
 
    else
 
       ! if no omega in dataset, use 3d or vert+horiz forcing.
       if(have_divT3d) then
          do i = 1,nlsf
             dtls(1:nzlsf,i) = divT3d_in(i,nzlsf:1:-1)
          end do
       elseif(have_vertdivT.and.have_divT) then
          do i = 1,nlsf
             dtls(1:nzlsf,i) = divT_in(i,nzlsf:1:-1) + vertdivT_in(i,nzlsf:1:-1)
          end do
       elseif(have_divT) then
          do i = 1,nlsf
             dtls(1:nzlsf,i) = divT_in(i,nzlsf:1:-1)
          end do
       end if
 
       if(have_divq3d) then
          do i = 1,nlsf
             dqls(1:nzlsf,i) = divq3d_in(i,nzlsf:1:-1)
          end do
       elseif(have_vertdivq.and.have_divq) then
          do i = 1,nlsf
             dqls(1:nzlsf,i) = divq_in(i,nzlsf:1:-1) + vertdivq_in(i,nzlsf:1:-1)
          end do
       elseif(have_divq) then
          do i = 1,nlsf
             dqls(1:nzlsf,i) = divq_in(i,nzlsf:1:-1)
          end do
       end if
 
    end if

    deallocate(divT_in,divT3d_in,vertdivT_in, &
         vertdivq_in,divq3d_in,divq_in,STAT=status)
    if(status/=0) then
       write(6,*) 'Could not de-allocate sounding/forcing arrays in readiopdata'
       call task_abort()
    end if

end subroutine readiopdata_core_frc



!=====================================================================
subroutine get_netcdf_dimlength( NCID, dimName, dimlength, status, required)
  implicit none

  ! input/output variables
  integer, intent(in)   :: NCID
  character, intent(in) :: dimName*(*)
  logical, intent(in) :: required

  integer, intent(out) :: status, dimlength

  ! local variables
  integer :: dimID

  ! get variable ID
  STATUS = NF90_INQ_DIMID( NCID, dimName, dimID )
  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find dimension ID for ', &
             dimName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note(readiopdata_core_module.f90): No dimension ID for ', dimName
        return
     endif
  endif

  STATUS = NF90_INQUIRE_DIMENSION( NCID, dimID, len=dimlength )
  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find length of ',dimName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note - readiopdata_core_module.f90 : Could not find length of ',&
             dimName
     endif
  endif

end subroutine get_netcdf_dimlength

!=====================================================================
subroutine get_netcdf_var1d_real( NCID, varName, var, use_NF90_REAL, &
     status, required)
  implicit none

  ! input/output variables
  integer, intent(in)   :: NCID
  character, intent(in) :: varName*(*)
  logical, intent(in) :: required, use_NF90_REAL

  integer, intent(out) :: status
  real(4), intent(inout) :: var(:)

  ! local variables
  integer :: varID, dimlen, ndims, n
  integer, allocatable :: dimids(:), start(:), dimlength(:)

  ! get variable ID
  STATUS = NF90_INQ_VARID( NCID, varName, varID )
  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find variable ID for ', &
             varName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note(readiopdata_core_module.f90): Optional variable ', varName,&
             ' not found'! in ', TRIM(iopfile) !Remove reference to iopfile
        return
     endif
  endif

  ! Determine the dimensions of the variable, and use the number and length of
  ! dimensions to set the values of start (start) and count (dimlength) for
  ! the call to nf90_get_var, respectively
  STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, ndims=ndims )
  allocate(dimids(ndims))
  allocate(start(ndims))
  allocate(dimlength(ndims))
  STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, dimids=dimids(:) )
  ! Loop through the dimensions, setting start to 1 and dimlength to the length
  ! for every dimension
  do n=1,ndims
    STATUS = NF90_INQUIRE_DIMENSION( NCID, dimids(n), len=dimlen )
    start(n) = 1
    dimlength(n) = dimlen
  end do

  ! Load all values connected with varID into var, using start and count
  STATUS = NF90_GET_VAR( NCID, varID, var(:), start=start(:), count=dimlength(:) )

  ! Deallocate local variables
  deallocate(dimids, start, dimlength)

  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find variable ', varName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note (readiopdata_core_module.f90) : Could not find ', varName
     endif
  endif

end subroutine get_netcdf_var1d_real

!=====================================================================
subroutine get_netcdf_var2d_real( NCID, varName, ntime, nlev, &
     var, use_NF90_REAL, status, required)
  !based on John Truesdale's getncdata_real_1d
  implicit none

  ! input/output variables
  integer, intent(in)   :: NCID, ntime, nlev
  character, intent(in) :: varName*(*)
  logical, intent(in) :: required, USE_NF90_REAL

  integer, intent(out) :: status
  real(4), intent(inout) :: var(:,:)

  ! local variables
  integer :: varID
  character     dim_name*( NF90_MAX_NAME )
  integer     var_dimIDs( NF90_MAX_VAR_DIMS )
  integer     var_ndims, dim_size, dims_set, i, n, var_type
  integer :: dimlen, ndims
  integer :: row_index, col_index, reshape_index
  integer, allocatable :: dimids(:), start(:), dimlength(:)
  real(4), allocatable :: var_tmp(:,:)

  ! get variable ID
  STATUS = NF90_INQ_VARID( NCID, varName, varID )
  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find variable ID for ',&
             varName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note(readiopdata_core_module.f90): Optional variable ', varName,&
             ' not found'! in ', TRIM(iopfile) ! Remove unnecessary reference to iopfile
        return
     endif
  endif


  ! Check the var variable's information with what we are expecting
  ! it to be.

  STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, xtype=var_type, ndims=var_ndims, dimids=var_dimIDs )
  if ( var_ndims > 4 ) then
     if(masterproc_) write(6,*) &
          'ERROR - getncdata.f90: The input var',varName, &
          'has', var_ndims, 'dimensions'
     STATUS = -1
     return
  endif

  if ( var_type /= NF90_FLOAT .and. var_type /= NF90_DOUBLE .and. &
       var_type /= NF90_INT ) then
     if(masterproc_) write(6,*) &
          'ERROR - getncdata.f90: The input var',varName, &
          'has unknown type', var_type
     STATUS = -1
     return
  endif

  if ( STATUS /= NF90_NOERR ) then
     if(masterproc_) write(6,*) &
          'ERROR - getncdata.f90:Cant get dimension IDs for', varName
     return
  endif

  ! I'm unsure if the above code is necessary after the change, but wil leave
  ! it there

  ! Allocate an array with the same size as var
  allocate(var_tmp(size(var,1), size(var,2)))

  ! Determine the dimensions of the variable, and use the number and length of
  ! dimensions to set the values of start (start) and count (dimlength) for
  ! the call to nf90_get_var, respectively
  STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, ndims=ndims )
  allocate(dimids(ndims))
  allocate(start(ndims))
  allocate(dimlength(ndims))
  STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, dimids=dimids(:) )
  ! Loop through the dimensions, setting start to 1 and dimlength to the length
  ! for every dimension
  do n=1,ndims
    STATUS = NF90_INQUIRE_DIMENSION( NCID, dimids(n), len=dimlen )
    start(n) = 1
    dimlength(n) = dimlen
  end do

  ! Load all values connected with varID into var, using start and count
  STATUS = NF90_GET_VAR( NCID, varID, var_tmp(:,:), start=start(:), count=dimlength(:) )

  ! Shift the values around in the array to match the organization in the file
  ! For some reason, the above call seems to result in incorrect order of
  ! elements
  !
  ! This section converts to the correct order:
  ! ---------     ---------
  ! |var_tmp|     |  var  |
  ! |-------|     |-------|
  ! | a | b |     | a | e |
  ! | c | d | --> | b | f |
  ! | e | f |     | c | g |
  ! | g | h |     | d | h |
  ! ---------     ---------
  !
  ! Conversion is done as if it is zero indexed
  ! The array is flattened to a single index using  `reshape = row * size(col) + col`
  ! This is converted back to a 2d index using `var(row,col) = var_tmp(reshape % size(row), reshape / size(row))`
  do row_index=0,size(var_tmp,1)-1
    do col_index=0,size(var_tmp,2)-1
      reshape_index=row_index*size(var_tmp,2)+col_index
      ! Set the values of var with a transpose from the values of var_tmp
      var(row_index+1, col_index+1) = var_tmp(modulo(reshape_index,size(var_tmp,1))+1,reshape_index/size(var_tmp,1)+1)
    end do
  end do

  ! Deallocate local variables
  deallocate(dimids,start,dimlength,var_tmp)

  if (STATUS /= NF90_NOERR ) then
     if(required) then
        if(masterproc_) write(6,*) &
             'ERROR(readiopdata_core_module.f90):Could not find variable ', &
             varName
        STATUS = NF90_CLOSE( NCID )
        call task_abort()
     else
        if(masterproc_) write(6,*) &
             'Note (readiopdata_core_module.f90) : Could not find ', varName
     endif
  endif

end subroutine get_netcdf_var2d_real

!------------------------------------------------------------------------
! File: calcdate.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id$
!
! Modified by Peter Blossey (pblossey@u.washington.edu)
!    to handle leap years and output calendar day.
!------------------------------------------------------------------------
!bloss #include <params.h>
subroutine calcdate(inDate, inSecs,  outDate, outSecs, outCalday)
!-----------------------------------------------------------------------
!  calcdate           Calculate Date from base date plus seconds
!
! INPUTS:
!
!	inDate	       Base date as YYMMDD.
!       inSecs         number of seconds the model has run
!
! OUTPUTS:
!       outDate        Current date as YYMMDD
!       outSecs        number of seconds into current date
!       outCalday      calendar day of current year (Jan-01 00:00:00 == 1.)
!
!
!-----------------------------------------------------------------------
! Computational notes: 
!
! 86400 is the number of seconds in 1 day.
!
! Dividing an integer by 10**n has the effect of right-shifting the           
! decimal digits n positions (ex: 861231/100 = 008612).
!
! mod(integer,10**n) has the effect of extracting the rightmost     
! n decimal digits of the integer (ex: mod(861231,10000) = 1231).
!
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: inDate
   integer, intent(in) :: inSecs       

!         
! Output arguments       
!                        
   integer, intent(out) :: outSecs
   integer, intent(out) :: outDate
#ifdef UWM_MISC
   real, intent(out) :: outCalday
#else
   real(4), intent(out) :: outCalday
#endif /*UWM_MISC*/
!
!---------------------------Local workspace-----------------------------
!
   integer     YY
   integer     MM
   integer     DD
   integer     i
   integer     iyear !bloss: added to deal with leap years
   integer     byear !bloss: added to deal with leap years
   integer     bmnth
   integer     bday
   integer     jday
   integer     jsec
   integer     jdcon(12)
   integer     ndm(12)
   integer     secs_per_year
   integer     secs_this_year  !bloss: added to deal with leap years
   integer     days_this_month  !bloss: added to deal with leap years

   data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
   data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
!
!-----------------------------------------------------------------------
!
! Check validity of input data
!
   byear = inDate/10000 
   if(byear<100) byear = 1900+byear !bloss: if 2 digit year, add 1900
   bmnth = mod(inDate,10000)/100
   bday =  mod(inDate,100)
   if (bmnth<1 .or. bmnth>12) then
      write(6,*)' CALCDATE: Invalid base month input:',bmnth
      call task_abort()
   end if
   if (bday<1 .or. bday>ndm(bmnth)) then
      write(6,*)' CALCDATE: Invalid base day of base date input:',bday
      call task_abort()
   end if
!
!
!
   jday = jdcon(bmnth) + bday

   !bloss: add a day if past February and this is a leap year.
   if( (bmnth>2).and. isLeapYear(byear) ) jday = jday + 1

   jsec  = (jday-1) * 86400 + insecs

   secs_per_year = 86400 * 365

   !bloss: count through years until jsec is less than a year.
   do iyear = byear,byear+insecs/secs_per_year
      secs_this_year = secs_per_year
      if(isLeapYear(iyear)) secs_this_year = secs_this_year + 86400 ! leap day
      if(jsec<secs_this_year) EXIT ! break from loop is jsec < one year
      jsec = jsec - secs_this_year
   end do

   YY  = mod(iyear,100) !bloss reduce year to two digit YY -- breaks 19YY

   outCalday = 1. + float(jsec)/86400. !bloss: compute calendar day

   !bloss: count through months until jsec is less than the next month
   do i=1, 12
      MM = i
      days_this_month = ndm(i)
      if(i==2.and.isLeapYear(iyear)) days_this_month = days_this_month + 1
      if(jsec<86400*days_this_month) EXIT
      jsec = jsec - 86400*days_this_month
   end do

   DD = jsec/86400 +1

   outSecs = mod(jsec,86400)

   outDate = YY*10000 + MM*100 + DD

!      write( *,* )'date =' , outDate
!
   return

contains

  logical function isLeapYear(iyear)
    implicit none
    integer iyear

    isLeapYear = mod(iyear,4)==0.and.mod(iyear,100)/=0 &
                          .or.mod(iyear,400)==0
  end function isLeapYear

end subroutine calcdate

end module readiopdata_core_module
