!-----------------------------------------------------------------------
! $Id$

! Module inputfields

!  This exists because I wanted to keep the grads_reader code 
!  generalized and bypass having to pass the datafile as a parameter,
!  since there may be situations where the fields are be calculated 
!  analytically or by calling a different model without reading a datafile.
!  Using a module also saves the trouble of writing an interface definition
!  within the hoc_inputfields code.
!===============================================================================
module inputfields
  implicit none

!----- Run information--------------------------------------------------
  character(len=80), public :: datafile

  character(len=100), public ::  & 
  datafilet, datafilem

  character(3), public      :: input_type

  logical, public :: input_um, input_vm, input_rtm, input_thlm, & 
                     input_wp2, input_wprtp, input_wpthlp,  & 
                     input_wp3, input_rtp2, input_thlp2,  & 
                     input_rtpthlp, input_upwp, input_vpwp, & 
                     input_ug, input_vg, input_rcm,  & 
                     input_wm_zt, input_exner, input_em, & 
                     input_p, input_rho, input_rho_zm, & 
                     input_Lscale, input_Lscale_up, input_Lscale_down, & 
                     input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
                     input_thvm, input_rrainm,input_Nrm,  & 
                     input_rsnowm, input_ricem, input_rgraupelm,  & 
                     input_thlm_forcing, input_rtm_forcing, & 
                     input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
                     input_Ncnm, input_Nim, input_cf, input_sigma_sqd_w_zt, &
                     input_veg_T_in_K, input_deep_soil_T_in_K, &
                     input_sfc_soil_T_in_K 


  public  :: grads_fields_reader, compute_timestep, set_filenames

  private :: lin_ext_zm_bottom, lin_ext_zt_bottom

  private ! Default Scope

  contains

!===============================================================================


!-----------------------------------------------------------------------
  subroutine set_filenames( )
!       Description: Set the names of the GrADS files to be used.
!       Used by hoc_inputfields.
!-----------------------------------------------------------------------

  implicit none

  select case ( input_type )
  case ( "les", "rf1" )
    datafilet = trim( datafile )//"_coamps_sm.ctl"
    datafilem = trim( datafile )//"_coamps_sw.ctl"
  case ( "hoc" )
    datafilet = trim( datafile )//"_zt.ctl"
    datafilem = trim( datafile )//"_zm.ctl"
  case default
    write(0,*) "Don't know how to handle input_type = "// & 
      input_type
    stop
  end select

  return
  end subroutine set_filenames
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  subroutine grads_fields_reader( timestep )
!       Description:
!       Reads in variables for the model from GrADS data

!       Calls:
!       subroutine open_grads_read
!       subroutine get_var
!       subroutine close_grads_read
!-----------------------------------------------------------------------

  use variables_prognostic_module, only: & 
      um,  & ! Variable(s)
      vm, & 
      rtm, & 
      thlm, & 
      wp2, & 
      wp3, & 
      wprtp, & 
      wpthlp, & 
      rtp2, & 
      thlp2, & 
      rtpthlp, & 
      upwp, & 
      vpwp, & 
      p_in_Pa, & 
      exner, & 
      rcm, & 
      wm_zt, & 
      rho, & 
      rho_zm, & 
      thlm_forcing, & 
      rtm_forcing, & 
      cf, & 
      tau_zm, & 
      up2, & 
      vp2, & 
      sigma_sqd_w

  use variables_diagnostic_module, only: & 
      hydromet,  & ! Variable(s)
      tau_zt, & 
      ug, & 
      vg, & 
      Lscale, & 
      Lscale_up, & 
      Lscale_down, & 
      Kh_zt, & 
      Kh_zm, & 
      thvm, & 
      Ncm, & 
      Ncnm, & 
      Nim, & 
      sigma_sqd_w_zt, & 
      em

  use grid_class, only: & 
      gr,  & ! Variable(s)
      zt2zm ! Procedure(s)

  use constants, only:  &
      rttol,   & ! Variable(s)
      wtol_sqd

  use array_index, only:  & 
      iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm

  use inputfile_class, only: & 
      inputgrads,  & ! Type
      get_var,  & ! Procedure(s)
      open_grads_read, & 
      close_grads_read

  use soil_vegetation, only: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

  implicit none

  ! Arguments
  integer, intent(in) :: timestep

  ! Local Variables
  logical :: l_error

  type (inputgrads) :: fread_var

  real, dimension(gr%nnzp+1) :: tmp1

  integer :: k


  select case( input_type )

  case( "hoc" )

    ! NOTE:  The code is not set up to compensate for grid 
    !        discrepancies between the CLUBB GrADS file that is 
    !        having its variable values passed in and the CLUBB 
    !        inputfields run that is using those values as 
    !        variable inputs.  Therefore, CLUBB should be set up
    !        to match the number of grid levels and altitude of 
    !        each grid level found in the CLUBB GrADS zt and zm
    !        files.

    !  Thermo grid - zt file 
    call open_grads_read( 15, trim( datafile )//"_zt.ctl",  & 
                          fread_var )

    if ( input_um ) then
      call get_var( fread_var, "um", timestep, & 
                    um(1:gr%nnzp),  l_error )
    endif

    if ( input_vm ) then 
      call get_var( fread_var, "vm", timestep, & 
                    vm(1:gr%nnzp),  l_error )
    endif

    if ( input_rtm ) then 
      call get_var( fread_var, "rtm", timestep, & 
                    rtm(1:gr%nnzp),  l_error )
    endif

    if ( input_thlm ) then 
      call get_var( fread_var, "thlm",  & 
                    timestep, & 
                    thlm(1:gr%nnzp),  l_error )
    endif

    if ( input_wp3 ) then 
      call get_var( fread_var, "wp3", timestep, & 
                    wp3(1:gr%nnzp),  l_error )
    endif
    if ( input_tau_zt ) then 
      call get_var( fread_var, "tau_zt", timestep, & 
                    tau_zt(1:gr%nnzp),  l_error )
    endif
    if ( input_rrainm ) then 
      call get_var( fread_var, "rrainm", timestep, & 
                    hydromet(1:gr%nnzp,iirrainm),  l_error )
    endif
    if ( input_rsnowm ) then 
      call get_var( fread_var, "rsnowm", timestep, & 
                    hydromet(1:gr%nnzp,iirsnowm),  l_error )
    endif
    if ( input_ricem ) then 
      call get_var( fread_var, "ricem", timestep, & 
                    hydromet(1:gr%nnzp,iiricem),  l_error )
    endif
    if ( input_rgraupelm ) then 
      call get_var( fread_var, "rgraupelm", timestep, & 
                    hydromet(1:gr%nnzp,iirgraupelm),  l_error )
    endif

!--------------------------------------------------------
! Added variables for hoc_restart
    if ( input_p ) then
      call get_var( fread_var, "p_in_Pa", timestep, & 
                    p_in_Pa(1:gr%nnzp),  l_error )
    endif
    if ( input_exner) then
      call get_var( fread_var , "exner", timestep, & 
                    exner(1:gr%nnzp),  l_error)
    endif
    if ( input_ug) then
      call get_var( fread_var , "ug", timestep, & 
                    ug(1:gr%nnzp),  l_error)
    endif
    if ( input_vg) then
      call get_var( fread_var , "vg", timestep, & 
                    vg(1:gr%nnzp),  l_error)
    endif
    if ( input_rcm) then
      call get_var( fread_var , "rcm", timestep, & 
                    rcm(1:gr%nnzp),  l_error)
    endif
    if ( input_wm_zt) then
      call get_var( fread_var , "wm", timestep, & 
                    wm_zt(1:gr%nnzp),  l_error)
    endif
    if ( input_rho) then
      call get_var( fread_var , "rho", timestep, & 
                    rho(1:gr%nnzp), l_error)
    endif
    if ( input_Lscale) then
      call get_var( fread_var , "lscale", timestep, & 
                    Lscale(1:gr%nnzp), l_error)
    endif
    if ( input_Lscale_up) then
      call get_var( fread_var , "Lscale_up", timestep, & 
                    Lscale_up(1:gr%nnzp), l_error)
    endif
    if ( input_Lscale_down) then
      call get_var( fread_var , "Lscale_down", timestep, & 
                    Lscale_down(1:gr%nnzp), l_error)
    endif
    if ( input_Kh_zt) then
      call get_var( fread_var , "Kh_zt", timestep, & 
                    Kh_zt(1:gr%nnzp), l_error)
    endif
    if ( input_thvm) then
      call get_var( fread_var , "thvm", timestep, & 
                    thvm(1:gr%nnzp), l_error)
    endif
    if ( input_thlm_forcing ) then
      call get_var( fread_var , "thlm_f", timestep, & 
                    thlm_forcing(1:gr%nnzp), l_error)
    endif
    if ( input_rtm_forcing ) then
      call get_var( fread_var , "rtm_f", timestep, & 
                    rtm_forcing(1:gr%nnzp), l_error)
    endif
    if ( input_Ncm) then
      call get_var( fread_var , "Ncm", timestep, & 
                    Ncm(1:gr%nnzp), l_error)
    endif
    if ( input_Ncnm) then
      call get_var( fread_var , "Ncnm", timestep, & 
                    Ncnm(1:gr%nnzp), l_error)
    endif
    if ( input_Nim) then
      call get_var( fread_var , "Nim", timestep, & 
                    Nim(1:gr%nnzp), l_error)
    endif
    if ( input_cf) then
      call get_var( fread_var , "cf", timestep, & 
                    cf(1:gr%nnzp), l_error)
    endif
    if ( input_Nrm ) then
      call get_var( fread_var , "Nrm", timestep, & 
                    hydromet(1:gr%nnzp,iiNrm), l_error)
    endif
    if ( input_sigma_sqd_w_zt ) then
      call get_var( fread_var , "sigma_sqd_w_zt", timestep, & 
                    sigma_sqd_w_zt(1:gr%nnzp), l_error)
    endif

!--------------------------------------------------------
    call close_grads_read( fread_var )

!         zm file
    call open_grads_read( 15, trim(datafile)//"_zm.ctl", & 
                          fread_var )

    if ( input_wp2) then 
      call get_var( fread_var, "wp2", timestep, & 
                    wp2(1:gr%nnzp),  l_error )
    endif

    if ( input_wprtp) then 
      call get_var( fread_var, "wprtp",  & 
                    timestep, wprtp(1:gr%nnzp), & 
                    l_error )
    endif

    if ( input_wpthlp) then 
      call get_var( fread_var, "wpthlp",  & 
                    timestep,  & 
                    wpthlp(1:gr%nnzp),  & 
                    l_error )
    endif

    if ( input_rtp2) then 
       call get_var( fread_var, "rtp2",  & 
                     timestep, & 
                     rtp2(1:gr%nnzp), l_error )
    endif

    if ( input_thlp2) then 
       call get_var( fread_var, "thlp2",  & 
                     timestep, & 
                     thlp2(1:gr%nnzp), l_error )
    endif

    if ( input_rtpthlp) then 
       call get_var( fread_var, "rtpthlp",  & 
                     timestep,  & 
                     rtpthlp(1:gr%nnzp), & 
                     l_error )
    endif

    if ( input_upwp) then 
       call get_var( fread_var, "upwp",  & 
                     timestep, & 
                     upwp(1:gr%nnzp), l_error )
    endif

    if ( input_vpwp) then 
       call get_var( fread_var, "vpwp",  & 
                     timestep, & 
                     vpwp(1:gr%nnzp), l_error )
    endif
!-----------------------------------------------------------
   if ( input_em) then
      call get_var( fread_var, "em", & 
                    timestep, & 
                    em(1:gr%nnzp), l_error )
   endif
   if ( input_rho_zm) then
      call get_var( fread_var, "rho_zm", & 
                    timestep, & 
                    rho_zm(1:gr%nnzp), l_error )
   endif
   if ( input_Kh_zm) then
      call get_var( fread_var, "Kh_zm", & 
                    timestep, & 
                    Kh_zm(1:gr%nnzp), l_error )
   endif
   if ( input_tau_zm) then
      call get_var( fread_var, "tau_zm", & 
                    timestep, & 
                    tau_zm(1:gr%nnzp), l_error )
   endif
   if ( input_up2) then
      call get_var( fread_var, "up2", & 
                    timestep, & 
                    up2(1:gr%nnzp), l_error )
   endif
   if ( input_vp2) then
      call get_var( fread_var, "vp2", & 
                    timestep, & 
                    vp2(1:gr%nnzp), l_error )
   endif
   if ( input_sigma_sqd_w ) then
      call get_var( fread_var, "sigma_sqd_w", & 
                    timestep, & 
                    sigma_sqd_w(1:gr%nnzp), l_error )
   endif

!-----------------------------------------------------------
    call close_grads_read( fread_var )

    call open_grads_read( 15, trim( datafile )//"_sfc.ctl",  & 
                          fread_var )

   if ( input_veg_T_in_K ) then
      call get_var( fread_var, "veg_T_in_K", & 
                    timestep, & 
                    tmp1, l_error )
      veg_T_in_K = tmp1(1)
      print *, "Veg T = ", veg_T_in_K
   endif
   if ( input_deep_soil_T_in_K ) then
      call get_var( fread_var, "deep_soil_T_in_", & 
                    timestep, & 
                    tmp1, l_error )
        deep_soil_T_in_K = tmp1(1)
      print *,"Deep soil = ",deep_soil_T_in_K
   endif
   if ( input_sfc_soil_T_in_K ) then
      call get_var( fread_var, "sfc_soil_T_in_K", & 
                    timestep, & 
                    tmp1, l_error )
        sfc_soil_T_in_K = tmp1(1)
        print *,"surface_soil = ", sfc_soil_T_in_K
   endif


    if ( l_error ) stop "oops, get_var failed in field_reader"

    call close_grads_read( fread_var )


  case( "rf1" )   ! special case for COAMPS DYCOMS-II RF01

    ! NOTE:  The code is not set up to compensate for discrepancies
    !        in thermodynamic level altitudes between CLUBB and 
    !        COAMPS LES (other than for thermodynamic level 
    !        indices).  Therefore, CLUBB should be set up so that 
    !        CLUBB thermodynamic level altitudes from thermodynamic 
    !        level indices 3 to gr%nnzp match COAMPS thermodynamic 
    !        level altitudes from thermodynamic level indices 1 to 
    !        gr%nnzp-2.

!         stats_sm
    call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl",  & 
                          fread_var )

    if ( input_um) then
      call get_var( fread_var, "um", timestep, & 
                    tmp1(1:gr%nnzp-2), l_error )
      ! tmp1 is the value of um from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.
      um(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
      ! Use the values of um at thermodynamic levels 4 and 3 
      ! to find the value at thermodynamic level 2 through the use
      ! of a linear extension.  Then, use the values of um at
      ! thermodynamic levels 3 and 2 to find the value at 
      ! thermodynamic level 1 through the use of a linear extension.
      um(2)  & 
      = lin_ext_zt_bottom( um(4), um(3), & 
                           gr%zt(4), gr%zt(3), gr%zt(2) )
      um(1)  & 
      = lin_ext_zt_bottom( um(3), um(2), & 
                           gr%zt(3), gr%zt(2), gr%zt(1) )
    endif

    if ( input_vm ) then 
      call get_var( fread_var, "vm", timestep, & 
                    tmp1(1:gr%nnzp-2), l_error )
      ! tmp1 is the value of vm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.
      vm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
      ! Use the values of vm at thermodynamic levels 4 and 3 
      ! to find the value at thermodynamic level 2 through the use
      ! of a linear extension.  Then, use the values of vm at
      ! thermodynamic levels 3 and 2 to find the value at 
      ! thermodynamic level 1 through the use of a linear extension.
      vm(2)  & 
      = lin_ext_zt_bottom( vm(4), vm(3), & 
                           gr%zt(4), gr%zt(3), gr%zt(2) )
      vm(1)  & 
      = lin_ext_zt_bottom( vm(3), vm(2), & 
                           gr%zt(3), gr%zt(2), gr%zt(1) )
    endif

    if ( input_rtm) then 
      call get_var( fread_var, "qtm", timestep, & 
                    tmp1(1:gr%nnzp-2), l_error )
      ! tmp1 is the value of rtm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.
      rtm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
      ! Set values of rtm at thermodynamic levels 2 and 1 to the
      ! value at thermodynamic level 3.
      rtm(1:2) = rtm(3)
    endif

    if ( input_thlm) then
      call get_var( fread_var, "thlm",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp-2), l_error )
      ! tmp1 is the value of thlm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.
      thlm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
      ! Set values of thlm at thermodynamic levels 2 and 1 to the
      ! value at thermodynamic level 3.
      thlm(1:2) = thlm(3)
    endif

    if ( input_wp3) then 
      call get_var( fread_var, "wp3", timestep, & 
                    tmp1(1:gr%nnzp-2), l_error )
      ! tmp1 is the value of wp3 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.
      wp3(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
      ! Set values of wp3 at thermodynamic levels 2 and 1 to 0.
      wp3(1:2) = 0.
    endif

    if ( input_wprtp) then
      call get_var( fread_var, "wpqtp",  & 
                    timestep, tmp1(1:gr%nnzp-1), & 
                    l_error )
      ! tmp1 is the value of wprtp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.  Interpolate the read-in values of 
      ! wprtp to their appropriate places on the momentum levels.
      wprtp(3:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-2) )
      ! Use the values of wprtp at momentum levels 4 and 3 to 
      ! find the value at momentum level 2 through the use of a
      ! linear extension.  Then, use the values of wprtp at
      ! momentum levels 3 and 2 to find the value at momentum
      ! level 1 through the use of a linear extension.  It should 
      ! be pointed out that the boundary flux is usually solved in
      ! LES or hoc via a subroutine like sfc_var.
      wprtp(2)  & 
      = lin_ext_zm_bottom( wprtp(4), wprtp(3), & 
                           gr%zm(4), gr%zm(3), gr%zm(2) )
      wprtp(1)  & 
      = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( input_wpthlp) then 
      call get_var( fread_var, "wpthlp",  & 
                    timestep, tmp1(1:gr%nnzp-1),  & 
                    l_error )
      ! tmp1 is the value of wpthlp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.  Interpolate the read-in values of 
      ! wpthlp to their appropriate places on the momentum levels.
      wpthlp(3:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-2) )
      ! Use the values of wpthlp at momentum levels 4 and 3 to 
      ! find the value at momentum level 2 through the use of a
      ! linear extension.  Then, use the values of wpthlp at
      ! momentum levels 3 and 2 to find the value at momentum
      ! level 1 through the use of a linear extension.  It should 
      ! be pointed out that the boundary flux is usually solved in
      ! LES or hoc via a subroutine like sfc_var.
      wpthlp(2)  & 
      = lin_ext_zm_bottom( wpthlp(4), wpthlp(3), & 
                           gr%zm(4), gr%zm(3), gr%zm(2) )
      wpthlp(1)  & 
      = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( input_rtp2) then 
      call get_var( fread_var, "qtp2",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp-1), l_error )
      ! tmp1 is the value of rtp2 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.  Interpolate the read-in values of 
      ! rtp2 to their appropriate places on the momentum levels.
      rtp2(3:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-2) )
      ! Using a linear extension here resulted in negatives.
      rtp2(1:2) =  rtp2(3)
      if ( any ( rtp2(1:gr%nnzp) < rttol**2 ) ) then
! %% debug
!              print *, "Some values of rtp2 are negative, compensating."
! %% debug
        do k=1, gr%nnzp
          rtp2(k) = max(rtp2(k), rttol**2)
        enddo
      endif
    endif

    if ( input_thlp2 ) then 
      call get_var( fread_var, "thlp2",  & 
                    timestep, tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of thlp2 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.  Interpolate the read-in values of 
      ! thlp2 to their appropriate places on the momentum levels.
      thlp2(3:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-2) )
      ! Using a linear extension here resulted in negatives.
      thlp2(1:2) = thlp2(3)
    endif

    if ( input_rtpthlp) then 
      call get_var( fread_var, "qtpthlp",  & 
                    timestep, tmp1(1:gr%nnzp-1), & 
                    l_error )
      ! tmp1 is the value of rtpthlp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! thermodynamic level 3 on the CLUBB grid for this 
      ! particular case.  Interpolate the read-in values of 
      ! rtpthlp to their appropriate places on the momentum levels.
      rtpthlp(3:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-2) )
      ! Use the values of rtpthlp at momentum levels 4 and 3 to 
      ! find the value at momentum level 2 through the use of a
      ! linear extension.  Then, use the values of rtpthlp at
      ! momentum levels 3 and 2 to find the value at momentum
      ! level 1 through the use of a linear extension.  It should 
      ! be pointed out that the boundary flux is usually solved in
      ! LES or hoc via a subroutine like sfc_var.
      rtpthlp(2)  & 
      = lin_ext_zm_bottom( rtpthlp(4), rtpthlp(3), & 
                           gr%zm(4), gr%zm(3), gr%zm(2) )
      rtpthlp(1)  & 
      = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( l_error ) stop "oops, get_var failed in field_reader"

    call close_grads_read( fread_var )


  case( "les" )   ! COAMPS LES -- all other cases.

    ! NOTE:  The code is not set up to compensate for discrepancies
    !        in thermodynamic level altitudes between CLUBB and 
    !        COAMPS LES (other than for thermodynamic level 
    !        indices).  Therefore, CLUBB should be set up so that 
    !        CLUBB thermodynamic level altitudes from thermodynamic 
    !        level indices 2 to gr%nnzp match COAMPS thermodynamic 
    !        level altitudes from thermodynamic level indices 1 to 
    !        gr%nnzp-1.

!         stats_sm
    call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl",  & 
                          fread_var )

    if ( input_um) then
      call get_var( fread_var, "um", timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of um from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid). 
      um(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
      ! Use the values of um at thermodynamic levels 3 and 2 
      ! to find the value at thermodynamic level 1 through the use
      ! of a linear extension.  
      um(1)  & 
      = lin_ext_zt_bottom( um(3), um(2), & 
                           gr%zt(3), gr%zt(2), gr%zt(1) )
    endif

    if ( input_vm) then 
      call get_var( fread_var, "vm", timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of vm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid). 
      vm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
      ! Use the values of um at thermodynamic levels 3 and 2 
      ! to find the value at thermodynamic level 1 through the use
      ! of a linear extension.  
      vm(1)  & 
      = lin_ext_zt_bottom( vm(3), vm(2), & 
                           gr%zt(3), gr%zt(2), gr%zt(1) )
    endif

    if ( input_rtm) then 
      call get_var( fread_var, "qtm", timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of rtm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid). 
      rtm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
      ! Set values of rtm at thermodynamic level 1 to the value 
      ! at thermodynamic level 2, as it is done in mixing.F.
      rtm(1) = rtm(2)
    endif

    if ( input_thlm) then
      call get_var( fread_var, "thlm",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of thlm from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid). 
      thlm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
      ! Set values of thlm at thermodynamic level 1 to the value 
      ! at thermodynamic level 2, as it is done in mixing.F.
      thlm(1) = thlm(2)
    endif

    if ( input_wp3) then 
      call get_var( fread_var, "wp3", timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of wp3 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid). 
      wp3(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
      ! Set values of wp3 at thermodynamic level 1 to 0, as it is
      ! done in wp23.F.
      wp3(1) = 0.  ! Computed as in hoc.F
    endif

!          if ( ( wp2 )) then 
!            call get_var( fread_var, "wp2", timestep,
!     .                    tmp1(1:gr%nnzp), l_error )
!            ! tmp1 is the value of wp2 from the LES GrADS file.  
!            ! It has been output onto thermodynamic levels starting at 
!            ! the first level above ground (thermodynamic level 2 on 
!            ! the CLUBB grid).  Interpolate the read-in values of 
!            ! wp2 to their appropriate places on the momentum levels.
!            wp2(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
!            ! Use the values of wp2 at momentum levels 3 and 2 to 
!            ! find the value at momentum level 1 through the use of a
!            ! linear extension.  It should be pointed out that the 
!            ! boundary flux is usually solved in LES or hoc via a 
!            ! subroutine like sfc_var.
!            wp2(1) 
!     .      = lin_ext_zm_bottom( wp2(3), wp2(2),
!     .                           gr%zm(3), gr%zm(2), gr%zm(1) )
!          endif

    if ( input_wprtp) then
      call get_var( fread_var, "wpqtp",  & 
                    timestep, tmp1(1:gr%nnzp), & 
                    l_error )
      ! tmp1 is the value of wprtp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid).  Interpolate the read-in values of 
      ! wprtp to their appropriate places on the momentum levels.
      wprtp(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
      ! Use the values of wprtp at momentum levels 3 and 2 to 
      ! find the value at momentum level 1 through the use of a
      ! linear extension.  It should be pointed out that the 
      ! boundary flux is usually solved in LES or hoc via a 
      ! subroutine like sfc_var.
      wprtp(1)  & 
      = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( input_wpthlp) then 
      call get_var( fread_var, "wpthlp",  & 
                    timestep, tmp1(1:gr%nnzp),  & 
                    l_error )
      ! tmp1 is the value of wpthlp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid).  Interpolate the read-in values of 
      ! wpthlp to their appropriate places on the momentum levels.
      wpthlp(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
      ! Use the values of wpthlp at momentum levels 3 and 2 to 
      ! find the value at momentum level 1 through the use of a
      ! linear extension.  It should be pointed out that the 
      ! boundary flux is usually solved in LES or hoc via a 
      ! subroutine like sfc_var.
      wpthlp(1)  & 
      = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( input_rtp2) then 
      call get_var( fread_var, "qtp2",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of rtp2 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid).  Interpolate the read-in values of 
      ! rtp2 to their appropriate places on the momentum levels.
      rtp2(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
      ! Using a linear extension here resulted in negatives.
      rtp2(1) =  rtp2(2)
      if ( any ( rtp2(1:gr%nnzp) < rttol**2 ) ) then
! %% debug
!              print *, "Some values of rtp2 are negative, compensating."
! %% debug
        do k=1, gr%nnzp
          rtp2(k) = max(rtp2(k), rttol**2)
        enddo
      endif
    endif

    if ( input_thlp2) then 
      call get_var( fread_var, "thlp2",  & 
                    timestep, tmp1(1:gr%nnzp), l_error )
      ! tmp1 is the value of thlp2 from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid).  Interpolate the read-in values of 
      ! thlp2 to their appropriate places on the momentum levels.
      thlp2(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
      ! Using a linear extension here resulted in negatives.
      thlp2(1) = thlp2(2)
    endif

    if ( input_rtpthlp) then 
      call get_var( fread_var, "qtpthlp",  & 
                    timestep,  & 
                    tmp1(1:gr%nnzp), & 
                    l_error )
      ! tmp1 is the value of rtpthlp from the LES GrADS file.  
      ! It has been output onto thermodynamic levels starting at 
      ! the first level above ground (thermodynamic level 2 on 
      ! the CLUBB grid).  Interpolate the read-in values of 
      ! rtpthlp to their appropriate places on the momentum levels.
      rtpthlp(2:gr%nnzp) = zt2zm( tmp1(1:gr%nnzp-1) )
      ! Use the values of rtpthlp at momentum levels 3 and 2 to 
      ! find the value at momentum level 1 through the use of a
      ! linear extension.  It should be pointed out that the 
      ! boundary flux is usually solved in LES or hoc via a 
      ! subroutine like sfc_var.
      rtpthlp(1)  & 
      = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                           gr%zm(3), gr%zm(2), gr%zm(1) )
    endif

    if ( l_error ) stop "oops, get_var failed in field_reader"

    call close_grads_read( fread_var )


  end select


  select case( input_type )

  case( "les", "rf1" )

    ! NOTE:  The code is not set up to compensate for discrepancies
    !        in momentum level altitudes between CLUBB and COAMPS 
    !        LES.  Therefore, CLUBB should be set up so that CLUBB 
    !        momentum level altitudes from momentum level indices 
    !        1 to gr%nnzp match COAMPS momentum level altitudes 
    !        from momentum level indices 1 to gr%nnzp.

    ! stats_sw
    call open_grads_read( 15, trim(datafile)//"_coamps_sw.ctl",  & 
                          fread_var )
   ! no interpolation is required, however, the stats_sw files have
   ! an extra top z-level, and wpup_sgs must be added to make the
   ! u'w' and v'w' terms as they are in CLUBB.

    if ( input_upwp) then 
      call get_var( fread_var, "wpup",  & 
                    timestep, tmp1(1:gr%nnzp+1), l_error )
      upwp(1:gr%nnzp) = tmp1(1:gr%nnzp) 

      call get_var( fread_var, "wpup_sgs",  & 
                    timestep, tmp1(1:gr%nnzp+1), l_error )
      upwp(1:gr%nnzp) = tmp1(1:gr%nnzp) + upwp(1:gr%nnzp)
    endif

    if ( l_error ) stop "get_var failed for upwp in field_reader"

    if ( input_vpwp) then
      call get_var( fread_var, "wpvp",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp+1), l_error )
      vpwp(1:gr%nnzp) = tmp1(1:gr%nnzp) 
      call get_var( fread_var, "wpvp_sgs",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp+1), l_error )
      vpwp(1:gr%nnzp) = tmp1(1:gr%nnzp) + vpwp(1:gr%nnzp)
    endif
    if ( l_error ) stop "get_var failed for vpwp in field_reader"

    if ( input_wp2 ) then 
      call get_var( fread_var, "wp2",  & 
                    timestep, & 
                    tmp1(1:gr%nnzp+1), l_error )
      wp2(1:gr%nnzp) = tmp1(1:gr%nnzp)
      if ( any ( wp2(1:gr%nnzp) < wtol_sqd ) ) then
! %% debug
!              print *, "Some values of wp2 are negative, compensating."
! %% debug
        do k=1, gr%nnzp
          wp2(k) = max(wp2(k), wtol_sqd)
        end do
      end if
    end if
    if ( l_error ) stop "get_var failed for wp2 in field_reader"

    call close_grads_read( fread_var )


  end select


  return
  end subroutine grads_fields_reader

!===============================================================================
  pure function lin_ext_zm_bottom( var_zmp2, var_zmp1, & 
                                   zmp2, zmp1, zm ) & 
  result( var_zm )

!       Description:
!       This function computes the value of a momentum-level variable
!       at a bottom grid level by using a linear extension of the values
!       of the variable at the two levels immediately above the level
!       where the result value is needed.

!-------------------------------------------------------------------------

  implicit none


  ! Input Variables
  real, intent(in) :: & 
  var_zmp2,    & ! Momentum level variable at level (k+2)   [units vary]
  var_zmp1,    & ! Momentum level variable at level (k+1)   [units vary]
  zmp2,        & ! Altitude at momentum level (k+2)         [m]
  zmp1,        & ! Altitude at momentum level (k+1)         [m]
  zm          ! Altitude at momentum level (k)           [m]

  ! Return Variable
  real :: var_zm  ! Momentum level variable at level (k) [units vary]

  var_zm = ( ( var_zmp2 - var_zmp1 ) / ( zmp2 - zmp1 ) ) & 
           * ( zm - zmp1 ) + var_zmp1

  return
  end function lin_ext_zm_bottom

!===============================================================================
  pure function lin_ext_zt_bottom( var_ztp2, var_ztp1, & 
                                   ztp2, ztp1, zt ) & 
  result( var_zt )

!       Description:
!       This function computes the value of a thermodynamic-level 
!       variable at a bottom grid level by using a linear extension of 
!       the values of the variable at the two levels immediately above 
!       the level where the result value is needed.

!-------------------------------------------------------------------------

  implicit none

  ! Input Variables
  real, intent(in) :: & 
  var_ztp2,    & ! Thermodynamic level variable at level (k+2)   [units vary]
  var_ztp1,    & ! Thermodynamic level variable at level (k+1)   [units vary]
  ztp2,        & ! Altitude at thermodynamic level (k+2)         [m]
  ztp1,        & ! Altitude at thermodynamic level (k+1)         [m]
  zt          ! Altitude at thermodynamic level (k)           [m]

  ! Return Variable
  real :: var_zt  ! Thermodynamic level variable at level (k) [units vary]

  var_zt = ( ( var_ztp2 - var_ztp1 ) / ( ztp2 - ztp1 ) ) & 
           * ( zt - ztp1 ) + var_ztp1

  return
  end function lin_ext_zt_bottom

!===============================================================================
  subroutine compute_timestep( iunit, filename, l_restart, & 
                               time, nearest_timestep )
!
!       Description: Given a time 'time', determines the closest 
!       output time in a GrADS file
!
!-------------------------------------------------------------------------

    use inputfile_class, only: & 
        inputgrads,  & ! Type
        open_grads_read,  & ! Procedure(s)
        close_grads_read
    use constants, only:  & 
        sec_per_min ! Variable(s)

    use stats_precision, only:  & 
        time_precision

    implicit none
    
    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) ::filename ! Path to the file and its name
    
    logical, intent(in) :: l_restart ! Whether this is a restart run

    real(kind=time_precision), intent(in) ::  & 
      time ! Time near which we want to find GrADS output,
           ! e.g. time_restart     [s]


    ! Output Variable(s)
    integer, intent(out) ::  & 
      nearest_timestep ! Nearest GrADS output time to time [min]
    
    ! Local Variables
    type (inputgrads) :: fread_var         

    real(kind=time_precision) :: delta_time   ! In seconds
    
    call open_grads_read( iunit, trim( filename ), fread_var )

    ! (restart time) - (initial time) 
    delta_time = time - (fread_var%time - fread_var%dtwrite)
    
    !    Joshua Fasching March 2008
!     .        time - fread_var%time
    
    ! Reporting
    if ( l_restart ) then
      print *, "Initial time of GrADS reference file ", & 
               "[seconds since midnight]: ",  & 
               fread_var%time
      print *, "Model restart time [s]: ", time
      print *, "Elapsed time between ", & 
               "initial time of ref file and restart time [s]: ",  & 
               delta_time
      print *, "GrADS file output time interval [s]: ",  & 
               fread_var%dtwrite

      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8 ) .or.  & 
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8 ) ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of the reference GrADS output time interval."
        print*, "Elapsed time [s] = ", delta_time
        print*, "GrADS output time interval = ", fread_var%dtwrite
        stop
      end if 

      if ( mod( delta_time , sec_per_min ) > 1e-8 & 
            .or. mod( delta_time, sec_per_min ) < -1e-8 ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of one minute."
        print*, "Elapsed time [s] = ", delta_time
        stop
      end if

    end if ! l_restart

    ! Determines the closest recorded timestep to the restart
    ! time.
    nearest_timestep = nint( delta_time / sec_per_min ) 

    if ( l_restart ) then 
      print *, "Elapsed time between ", & 
               "initial time of ref file and restart time ", & 
               "rounded to nearest minute: ",  & 
               nearest_timestep

      ! Print the actual record being recalled.
      ! Joshua Fasching March 2008
      print *, "Nearest GrADS output time iteration [ ]: ", & 
               nint( nearest_timestep /  & 
                     (fread_var%dtwrite/sec_per_min) ) - 1
    end if ! l_restart
    
    call close_grads_read( fread_var )
    
  end subroutine compute_timestep
  
!===============================================================================

end module inputfields
