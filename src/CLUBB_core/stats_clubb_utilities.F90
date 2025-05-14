!-----------------------------------------------------------------------
!  $Id$
!===============================================================================
module stats_clubb_utilities

  implicit none

  private ! Set Default Scope

  public :: stats_init, stats_init_w_diff_output_gr, &
    stats_begin_timestep, stats_end_timestep, &
    stats_end_timestep_w_diff_output_gr, &
    stats_accumulate, stats_finalize, stats_accumulate_hydromet, &
    stats_accumulate_lh_tend

  private :: stats_zero, stats_avg, stats_check_num_samples, stats_init_helper, &
             stats_end_timestep_helper

  contains

  !-----------------------------------------------------------------------
  subroutine stats_init_helper( iunit, fname_prefix, fdir, l_stats_in, &
                                stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                                hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, &
                                hydromet_list, l_mix_rat_hm, &
                                nzmax, ngrdcol, nlon, nlat, gzt, gzm, nnrad_zt, &
                                grad_zt, nnrad_zm, grad_zm, day, month, year, &
                                lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                                clubb_params, &
                                l_uv_nudge, &
                                l_tke_aniso, &
                                l_standard_term_ta, &
                                l_different_output_grid, &
                                output_gr_nzm, &
                                output_gr_zt, &
                                output_gr_zm, &
                                stats_metadata, &
                                stats_zt, stats_zm, stats_sfc, &
                                stats_lh_zt, stats_lh_sfc, &
                                stats_rad_zt, stats_rad_zm, &
                                err_code )
    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: &
        stats_metadata_type

    use parameter_indices, only: &
        nparams

    use clubb_precision, only: &
        time_precision, & ! Constant(s)
        core_rknd

    use output_grads, only: &
        open_grads ! Procedure

#ifdef NETCDF
    use output_netcdf, only: &
        open_netcdf_for_writing, & ! Procedure
        first_write
#endif

    use stats_zm_module, only: &
        nvarmax_zm, & ! Constant(s)
        stats_init_zm ! Procedure(s)

    use stats_zt_module, only: &
        nvarmax_zt, & ! Constant(s)
        stats_init_zt ! Procedure(s)

    use stats_lh_zt_module, only: &
        nvarmax_lh_zt, & ! Constant(s)
        stats_init_lh_zt ! Procedure(s)

    use stats_lh_sfc_module, only: &
        nvarmax_lh_sfc, & ! Constant(s)
        stats_init_lh_sfc ! Procedure(s)

    use stats_rad_zt_module, only: &
        nvarmax_rad_zt, & ! Constant(s)
        stats_init_rad_zt ! Procedure(s)

    use stats_rad_zm_module, only: &
        nvarmax_rad_zm, & ! Constant(s)
        stats_init_rad_zm ! Procedure(s)

    use stats_sfc_module, only: &
        nvarmax_sfc, & ! Constant(s)
        stats_init_sfc ! Procedure(s)

    use constants_clubb, only: &
        fstdout, fstderr, var_length ! Constants

    use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        clubb_fatal_error               ! Constant

    use stats_type, only: stats ! Type

    implicit none

    ! Local Constants
    integer, parameter :: &
      silhs_num_importance_categories = 8

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real( kind = core_rknd ), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      ngrdcol,  &      ! Number of columns
      nlon,     &      ! Number of points in the X direction [-]
      nlat,     &      ! Number of points in the Y direction [-]
      nzmax,    &      ! Grid points in the vertical         [-]
      output_gr_nzm    ! Grid points in the vertical for the grid that gets written to file;
                       ! can vary from nzmax since output is remapped to dycore grid if
                       ! we adapt grid and want to simulate forcings from dycore grid,
                       ! since we need a common grid for all results if the grid is adapted
                       ! over time, and the dycore grid can have a different number of
                       ! levels as the physics grid          [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax-1) ::  & 
      gzt       ! Thermodynamic levels           [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) ::  & 
      gzm       ! Momentum levels                [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,output_gr_nzm-1) ::  & 
      output_gr_zt  ! Thermodynamic levels of the grid that gets written to file    [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,output_gr_nzm) ::  & 
      output_gr_zm  ! Momentum levels of the grid that gets written to file         [m]
    
    ! Note: output_gr_zt and output_gr_zm are only different from gzt and gzm if we adapt the grid
    !       over time, since then we need a common grid we can remap all the results to,
    !       that were calculated on different grids

    integer, intent(in) :: &
      nnrad_zt,     & ! Grid points in the radiation grid [count]
      hydromet_dim, &
      sclr_dim,     &
      edsclr_dim

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    character(len=10), dimension(hydromet_dim), intent(in) :: & 
      hydromet_list

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt, stats_lh_sfc)  [boolean]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,          &   ! For wind speed nudging
      l_tke_aniso,         &   ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                               ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta,  &   ! Use the standard discretization for the turbulent advection terms.
                               ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                               ! derivative in advance_wp2_wp3_module.F90 and in
                               ! advance_xp2_xpyp_module.F90.
      l_different_output_grid  ! use different grid to output values to file

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code    ! Error code catching and relaying any errors occurring in this subroutine

    ! Local Variables
    logical :: l_error

    character(len=200) :: fname

    integer :: i, ivar, ntot, read_status

    ! Namelist Variables

    character(len=10) :: stats_fmt  ! File storage convention

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_lh_zt) ::  & 
      vars_lh_zt  ! Latin Hypercube variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_lh_sfc) ::  & 
      vars_lh_sfc  ! Latin Hypercube variables at the surface

    character(len=var_length), dimension(nvarmax_zm) ::  & 
      vars_zm  ! Variables on the momentum levels

    character(len=var_length), dimension(nvarmax_rad_zt) ::  & 
      vars_rad_zt  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_rad_zm) ::  & 
      vars_rad_zm  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_sfc) ::  &
      vars_sfc ! Variables at the model surface

    namelist /statsnl/ & 
      vars_zt, & 
      vars_zm, &
      vars_lh_zt, &
      vars_lh_sfc, &
      vars_rad_zt, &
      vars_rad_zm, & 
      vars_sfc

    ! ------------------- Begin Code -------------------

    ! Initialize
    l_error = .false.

    ! Set stats_variables variables with inputs from calling subroutine
    stats_metadata%l_stats = l_stats_in

    stats_metadata%stats_tsamp = stats_tsamp_in
    stats_metadata%stats_tout  = stats_tout_in
    stats_fmt   = trim( stats_fmt_in )
    stats_metadata%l_silhs_out = l_silhs_out_in

    if ( .not. stats_metadata%l_stats ) then
      stats_metadata%l_stats_samp  = .false.
      stats_metadata%l_stats_last  = .false.
      return
    end if

    ! Initialize namelist variables

    vars_zt  = ''
    vars_zm  = ''
    vars_lh_zt = ''
    vars_lh_sfc = ''
    vars_rad_zt = ''
    vars_rad_zm = ''
    vars_sfc = ''

    ! Reads list of variables that should be output to GrADS/NetCDF (namelist &statsnl)

    open(unit=iunit, file=fnamelist)
    read(unit=iunit, nml=statsnl, iostat=read_status, end=100)
    if ( read_status /= 0 ) then
      if ( read_status > 0 ) then
        write(fstderr,*) "Error reading stats namelist in file ",  &
                         trim( fnamelist )
      else ! Read status < 0
        write(fstderr,*) "End of file marker reached while reading stats namelist in file ", &
          trim( fnamelist )
      end if
      write(fstderr,*) "One cause is having more statistical variables ",  &
                       "listed in the namelist for var_zt, var_zm, or ",  &
                       "var_sfc than allowed by nvarmax_zt, nvarmax_zm, ",  &
                       "or nvarmax_sfc, respectively."
      write(fstderr,*) "Maximum variables allowed for var_zt = ", nvarmax_zt
      write(fstderr,*) "Maximum variables allowed for var_zm = ", nvarmax_zm
      write(fstderr,*) "Maximum variables allowed for var_rad_zt = ", nvarmax_rad_zt
      write(fstderr,*) "Maximum variables allowed for var_rad_zm = ", nvarmax_rad_zm
      write(fstderr,*) "Maximum variables allowed for var_sfc = ", nvarmax_sfc
      write(fstderr,*) "stats_init: Error reading stats namelist."
      err_code = clubb_fatal_error
      close(unit=iunit)
      return
    end if ! read_status /= 0

    close(unit=iunit)

    if ( clubb_at_least_debug_level( 1 ) ) then
      write(fstdout,*) "--------------------------------------------------"

      write(fstdout,*) "Statistics"

      write(fstdout,*) "--------------------------------------------------"
      write(fstdout,*) "vars_zt = "
      ivar = 1
      do while ( vars_zt(ivar) /= '' )
        write(fstdout,*) vars_zt(ivar)
        ivar = ivar + 1
      end do

      write(fstdout,*) "vars_zm = "
      ivar = 1
      do while ( vars_zm(ivar) /= '' )
        write(fstdout,*) vars_zm(ivar)
        ivar = ivar + 1
      end do

      if ( stats_metadata%l_silhs_out ) then
        write(fstdout,*) "vars_lh_zt = "
        ivar = 1
        do while ( vars_lh_zt(ivar) /= '' )
          write(fstdout,*) vars_lh_zt(ivar)
          ivar = ivar + 1
        end do

        write(fstdout,*) "vars_lh_sfc = "
        ivar = 1
        do while ( vars_lh_sfc(ivar) /= '' )
          write(fstdout,*) vars_lh_sfc(ivar)
          ivar = ivar + 1
        end do
      end if ! l_silhs_out

      if ( stats_metadata%l_output_rad_files ) then
        write(fstdout,*) "vars_rad_zt = "
        ivar = 1
        do while ( vars_rad_zt(ivar) /= '' )
          write(fstdout,*) vars_rad_zt(ivar)
          ivar = ivar + 1
        end do

        write(fstdout,*) "vars_rad_zm = "
        ivar = 1
        do while ( vars_rad_zm(ivar) /= '' )
          write(fstdout,*) vars_rad_zm(ivar)
          ivar = ivar + 1
        end do
      end if ! l_output_rad_files

      write(fstdout,*) "vars_sfc = "
      ivar = 1
      do while ( vars_sfc(ivar) /= '' )
        write(fstdout,*) vars_sfc(ivar)
        ivar = ivar + 1
      end do

      write(fstdout,*) "--------------------------------------------------"
    end if ! clubb_at_least_debug_level 1

    ! Determine file names for GrADS or NetCDF files
    stats_metadata%fname_zt  = trim( fname_prefix )//"_zt"
    stats_metadata%fname_zm  = trim( fname_prefix )//"_zm"
    stats_metadata%fname_lh_zt  = trim( fname_prefix )//"_lh_zt"
    stats_metadata%fname_lh_sfc  = trim( fname_prefix )//"_lh_sfc"
    stats_metadata%fname_rad_zt  = trim( fname_prefix )//"_rad_zt"
    stats_metadata%fname_rad_zm  = trim( fname_prefix )//"_rad_zm"
    stats_metadata%fname_sfc = trim( fname_prefix )//"_sfc"

    ! Parse the file type for stats output.  Currently only GrADS and
    ! netCDF > version 3.5 are supported by this code.
    select case ( trim( stats_fmt ) )
    case ( "GrADS", "grads", "gr" )
      stats_metadata%l_netcdf = .false.
      stats_metadata%l_grads  = .true.

    case ( "NetCDF", "netcdf", "nc" )
      stats_metadata%l_netcdf = .true.
      stats_metadata%l_grads  = .false.

    case default
      write(fstderr,*) "In module stats_clubb_utilities subroutine stats_init: "
      write(fstderr,*) "Invalid stats output format "//trim( stats_fmt )
      err_code = clubb_fatal_error
      return

    end select

    ! Check sampling and output frequencies

    ! The model time step length, delt (which is dt_main), should multiply
    ! evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_metadata%stats_tsamp/delt &
              - real( floor( stats_metadata%stats_tsamp/delt ), kind=core_rknd ) )  & 
           > 1.e-8_core_rknd) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                       'delt (which is dt_main).  Check the appropriate ',  &
                       'model.in file.'
      write(fstderr,*) 'stats_tsamp = ', stats_metadata%stats_tsamp
      write(fstderr,*) 'delt = ', delt
    end if

    ! The statistical sampling time step length, stats_tsamp, should multiply
    ! evenly into the statistical output time step length, stats_tout.
    if ( abs( stats_metadata%stats_tout/stats_metadata%stats_tsamp &
           - real( floor( stats_metadata%stats_tout    &
                          / stats_metadata%stats_tsamp ), kind=core_rknd) ) & 
         > 1.e-8_core_rknd) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tout should be an even multiple of ',  &
                       'stats_tsamp.  Check the appropriate model.in file.'
      write(fstderr,*) 'stats_tout = ', stats_metadata%stats_tout
      write(fstderr,*) 'stats_tsamp = ', stats_metadata%stats_tsamp
    end if

    ! Initialize zt (mass points)

    ivar = 1
    do while ( ichar(vars_zt(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_zt(ivar)) /= 0 & 
               .and. ivar <= nvarmax_zt )
      ivar = ivar + 1
    end do
    ntot = ivar - 1

    if ( any( vars_zt == "hm_i" ) ) then
       ! Correct for number of variables found under "hm_i".
       ! Subtract "hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_hm_i" ) ) then
       ! Correct for number of variables found under "mu_hm_i".
       ! Subtract "mu_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_Ncn_i" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i".
       ! Subtract "mu_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "mu_hm_i_n" ) ) then
       ! Correct for number of variables found under "mu_hm_i_n".
       ! Subtract "mu_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i_n".
       ! Subtract "mu_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "sigma_hm_i" ) ) then
       ! Correct for number of variables found under "sigma_hm_i".
       ! Subtract "sigma_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "sigma_Ncn_i" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i".
       ! Subtract "sigma_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "sigma_hm_i_n" ) ) then
       ! Correct for number of variables found under "sigma_hm_i_n".
       ! Subtract "sigma_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "sigma_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i_n".
       ! Subtract "sigma_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif

    if ( any( vars_zt == "corr_w_hm_i" ) ) then
       ! Correct for number of variables found under "corr_w_hm_i".
       ! Subtract "corr_w_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_w_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_w_Ncn_i".
       ! Subtract "corr_w_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_chi_hm_i" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i".
       ! Subtract "corr_chi_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i".
       ! Subtract "corr_chi_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_eta_hm_i" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i".
       ! Subtract "corr_eta_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i".
       ! Subtract "corr_eta_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i" ) ) then
       ! Correct for number of variables found under "corr_Ncn_hm_i".
       ! Subtract "corr_Ncn_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i" ) ) then
       ! Correct for number of variables found under "corr_hmx_hmy_i".
       ! Subtract "corr_hmx_hmy_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of correlations of two hydrometeors, which is found by:
       ! (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zt statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 )
    endif

    if ( any( vars_zt == "corr_w_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_w_hm_i_n".
       ! Subtract "corr_w_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_w_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_w_Ncn_i_n".
       ! Subtract "corr_w_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_chi_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i_n".
       ! Subtract "corr_chi_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i_n".
       ! Subtract "corr_chi_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_eta_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i_n".
       ! Subtract "corr_eta_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i_n".
       ! Subtract "corr_eta_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_Ncn_hm_i_n".
       ! Subtract "corr_Ncn_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i_n" ) ) then
       ! Correct for number of variables found under "corr_hmx_hmy_i_n".
       ! Subtract "corr_hmx_hmy_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of normal space correlations of two hydrometeors, which is
       ! found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zt statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 )
    endif

    if ( any( vars_zt == "hmp2_zt" ) ) then
       ! Correct for number of variables found under "hmp2_zt".
       ! Subtract "hmp2_zt" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zt statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zt == "wp2hmp" ) ) then
       ! Correct for number of variables found under "wp2hmp".
       ! Subtract "wp2hmp" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zt statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zt == "sclrm" ) ) then
       ! Correct for number of variables found under "sclrm".
       ! Subtract "sclrm" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + sclr_dim
    endif   

    if ( any( vars_zt == "sclrm_f" ) ) then
       ! Correct for number of variables found under "sclrm_f".
       ! Subtract "sclrm_f" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + sclr_dim
    endif

    if ( any( vars_zt == "edsclrm" ) ) then
       ! Correct for number of variables found under "edsclrm".
       ! Subtract "edsclrm" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + edsclr_dim
    endif

    if ( any( vars_zt == "edsclrm_f" ) ) then
       ! Correct for number of variables found under "edsclrm_f".
       ! Subtract "edsclrm_f" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + edsclr_dim
    endif

    if ( ntot >= nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      write(fstderr,*) "number of variables in vars_zt = ", ntot
      write(fstderr,*) "stats_init:  number of zt statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return
    end if

    ! Allocate and initialize for all columns, but we only write for 1
    do i = 1, ngrdcol
      stats_zt(i)%num_output_fields = ntot
      stats_zt(i)%kk = nzmax - 1
      stats_zt(i)%ii = nlon
      stats_zt(i)%jj = nlat

      allocate( stats_zt(i)%z( stats_zt(i)%kk ) )
      stats_zt(i)%z = gzt(i,:)

      allocate( stats_zt(i)%accum_field_values( stats_zt(i)%ii, stats_zt(i)%jj, &
        stats_zt(i)%kk, stats_zt(i)%num_output_fields ) )
      allocate( stats_zt(i)%accum_num_samples( stats_zt(i)%ii, stats_zt(i)%jj, &
        stats_zt(i)%kk, stats_zt(i)%num_output_fields ) )
      allocate( stats_zt(i)%l_in_update( stats_zt(i)%ii, stats_zt(i)%jj, stats_zt(i)%kk, &
        stats_zt(i)%num_output_fields ) )
      call stats_zero( stats_zt(i)%ii, stats_zt(i)%jj, stats_zt(i)%kk, & !In
                       stats_zt(i)%num_output_fields, & ! In
                       stats_zt(i)%accum_field_values, stats_zt(i)%accum_num_samples, & !Out
                       stats_zt(i)%l_in_update ) ! Out

      allocate( stats_zt(i)%file%grid_avg_var( stats_zt(i)%num_output_fields ) )
      if ( l_different_output_grid ) then
        allocate( stats_zt(i)%file%z( output_gr_nzm-1 ) )
      else
        allocate( stats_zt(i)%file%z( stats_zt(i)%kk ) )
      end if

      ! Default initialization for array indices for zt

      call stats_init_zt( hydromet_dim, sclr_dim, edsclr_dim, & ! intent(in)
                          hydromet_list, l_mix_rat_hm,        & ! intent(in)
                          vars_zt,                            & ! intent(in)
                          l_error,                            & ! intent(inout)
                          stats_metadata, stats_zt(i) )         ! intent(inout)
    end do

    fname = trim( stats_metadata%fname_zt )

    if ( stats_metadata%l_grads ) then

      ! Open GrADS file
      if ( l_different_output_grid ) then
        call open_grads( iunit, fdir, fname,  &  ! In
                         1, output_gr_nzm-1, nlat, nlon, output_gr_zt(1,:), & ! In 
                         day, month, year, lat_vals, lon_vals, &  ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_zt(1)%num_output_fields, & ! In
                         stats_zt(1)%file ) ! intent(inout)
      else
        call open_grads( iunit, fdir, fname,  &  ! In
                         1, stats_zt(1)%kk, nlat, nlon, stats_zt(1)%z, & ! In 
                         day, month, year, lat_vals, lon_vals, &  ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_zt(1)%num_output_fields, & ! In
                         stats_zt(1)%file ) ! intent(inout)
      end if

    else ! Open NetCDF file
#ifdef NETCDF
      if ( l_different_output_grid ) then
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, output_gr_nzm-1, &
                                      output_gr_zt(1,:), day, month, year, lat_vals, lon_vals, & ! In
                                      time_current, stats_metadata%stats_tout, & ! In
                                      stats_zt(1)%num_output_fields, & ! In
                                      stats_zt(1)%file, err_code ) ! InOut
      else
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_zt(1)%kk, &
                                      stats_zt(1)%z, day, month, year, lat_vals, lon_vals, & ! In
                                      time_current, stats_metadata%stats_tout, & ! In
                                      stats_zt(1)%num_output_fields, & ! In
                                      stats_zt(1)%file, err_code ) ! InOut
      end if

      ! Finalize the variable definitions
      call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                        l_uv_nudge, & ! intent(in)
                        l_tke_aniso, & ! intent(in)
                        l_standard_term_ta, & ! intent(in)
                        stats_zt(1)%file, err_code ) ! intent(inout)

      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Fatal error setting up stats_zt"
        return
      end if
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif

    end if


    ! Setup output file for stats_lh_zt (Latin Hypercube stats)

    if ( stats_metadata%l_silhs_out ) then

      ivar = 1
      do while ( ichar(vars_lh_zt(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_lh_zt(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_lh_zt )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( any( vars_lh_zt == "silhs_variance_category" ) ) then
        ! Correct for number of variables found under "silhs_variance_category".
        ! Subtract "silhs_variance_category" from the number of lh_zt statistical
        ! variables.
        ntot = ntot - 1
        ! Add 1 for each SILHS category to the number of lh_zt statistical variables
        ntot = ntot + silhs_num_importance_categories
      end if

      if ( any( vars_lh_zt == "lh_samp_frac_category" ) ) then
        ! Correct for number of variables found under "lh_samp_frac_category".
        ! Subtract "lh_samp_frac_category" from the number of lh_zt statistical
        ! variables.
        ntot = ntot - 1
        ! Add 1 for each SILHS category to the number of lh_zt statistical variables
        ntot = ntot + silhs_num_importance_categories
      end if

      if ( ntot == nvarmax_lh_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_lh_zt than allowed for by nvarmax_lh_zt."
        write(fstderr,*) "Check the number of variables listed for vars_lh_zt ",  &
                         "in the stats namelist, or change nvarmax_lh_zt."
        write(fstderr,*) "nvarmax_lh_zt = ", nvarmax_lh_zt
        write(fstderr,*) "number of variables in vars_lh_zt = ", ntot
        write(fstderr,*) "stats_init:  number of lh_zt statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      ! Allocate and initialize for all columns, but we only write for 1
      do i = 1, ngrdcol
        stats_lh_zt(i)%num_output_fields = ntot
        stats_lh_zt(i)%kk = nzmax - 1
        stats_lh_zt(i)%ii = nlon
        stats_lh_zt(i)%jj = nlat

        allocate( stats_lh_zt(i)%z( stats_lh_zt(i)%kk ) )
        stats_lh_zt(i)%z = gzt(i,:)

        allocate( stats_lh_zt(i)%accum_field_values( stats_lh_zt(i)%ii, stats_lh_zt(i)%jj, &
          stats_lh_zt(i)%kk, stats_lh_zt(i)%num_output_fields ) )
        allocate( stats_lh_zt(i)%accum_num_samples( stats_lh_zt(i)%ii, stats_lh_zt(i)%jj, &
          stats_lh_zt(i)%kk, stats_lh_zt(i)%num_output_fields ) )
        allocate( stats_lh_zt(i)%l_in_update( stats_lh_zt(i)%ii, stats_lh_zt(i)%jj, &
          stats_lh_zt(i)%kk, stats_lh_zt(i)%num_output_fields ) )
        call stats_zero( stats_lh_zt(i)%ii, stats_lh_zt(i)%jj, stats_lh_zt(i)%kk, & ! intent(in)
          stats_lh_zt(i)%num_output_fields, & ! intent(in)
          stats_lh_zt(i)%accum_field_values, stats_lh_zt(i)%accum_num_samples, & ! intent(out)
          stats_lh_zt(i)%l_in_update ) ! intent(out)

        allocate( stats_lh_zt(i)%file%grid_avg_var( stats_lh_zt(i)%num_output_fields ) )
        allocate( stats_lh_zt(i)%file%z( stats_lh_zt(i)%kk ) )

        call stats_init_lh_zt( vars_lh_zt,                    & ! intent(in)
                              l_error,                        & ! intent(inout)
                              stats_metadata, stats_lh_zt(i) )  ! intent(inout)
      end do


      fname = trim( stats_metadata%fname_lh_zt )

      if ( stats_metadata%l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & ! In
                         1, stats_lh_zt(1)%kk, nlat, nlon, stats_lh_zt(1)%z, & ! In
                         day, month, year, lat_vals, lon_vals, &  ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_lh_zt(1)%num_output_fields, & ! In
                         stats_lh_zt(1)%file ) ! In/Out

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_lh_zt(1)%kk, & ! In
                          stats_lh_zt(1)%z, day, month, year, lat_vals, lon_vals, &  ! In
                          time_current, stats_metadata%stats_tout, & ! In
                          stats_lh_zt(1)%num_output_fields, & ! In
                          stats_lh_zt(1)%file, err_code ) ! InOut

        ! Finalize the variable definitions
        call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                          l_uv_nudge, & ! intent(in)
                          l_tke_aniso, & ! intent(in)
                          l_standard_term_ta, & ! intent(in)
                          stats_lh_zt(1)%file, err_code ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error setting up stats_lh_zt"
          return
        end if
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

      ivar = 1
      do while ( ichar(vars_lh_sfc(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_lh_sfc(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_lh_sfc )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_lh_sfc ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_lh_sfc than allowed for by nvarmax_lh_sfc."
        write(fstderr,*) "Check the number of variables listed for vars_lh_sfc ",  &
                         "in the stats namelist, or change nvarmax_lh_sfc."
        write(fstderr,*) "nvarmax_lh_sfc = ", nvarmax_lh_sfc
        write(fstderr,*) "number of variables in vars_lh_sfc = ", ntot
        write(fstderr,*) "stats_init:  number of lh_sfc statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      ! Allocate and initialize for all columns, but we only write for 1
      do i = 1, ngrdcol
        stats_lh_sfc(i)%num_output_fields = ntot
        stats_lh_sfc(i)%kk = 1
        stats_lh_sfc(i)%ii = nlon
        stats_lh_sfc(i)%jj = nlat

        allocate( stats_lh_sfc(i)%z( stats_lh_sfc(i)%kk ) )
        stats_lh_sfc(i)%z = gzm(i,1)

        allocate( stats_lh_sfc(i)%accum_field_values( stats_lh_sfc(i)%ii, stats_lh_sfc(i)%jj, &
          stats_lh_sfc(i)%kk, stats_lh_sfc(i)%num_output_fields ) )
        allocate( stats_lh_sfc(i)%accum_num_samples( stats_lh_sfc(i)%ii, stats_lh_sfc(i)%jj, &
          stats_lh_sfc(i)%kk, stats_lh_sfc(i)%num_output_fields ) )
        allocate( stats_lh_sfc(i)%l_in_update( stats_lh_sfc(i)%ii, stats_lh_sfc(i)%jj, &
          stats_lh_sfc(i)%kk, stats_lh_sfc(i)%num_output_fields ) )

        call stats_zero( stats_lh_sfc(i)%ii, stats_lh_sfc(i)%jj, stats_lh_sfc(i)%kk, & ! intent(in)
            stats_lh_sfc(i)%num_output_fields, & ! intent(in)
            stats_lh_sfc(i)%accum_field_values, & ! intent(out)
            stats_lh_sfc(i)%accum_num_samples, stats_lh_sfc(i)%l_in_update ) ! intent(out)

        allocate( stats_lh_sfc(i)%file%grid_avg_var( stats_lh_sfc(i)%num_output_fields ) )
        allocate( stats_lh_sfc(i)%file%z( stats_lh_sfc(i)%kk ) )

        call stats_init_lh_sfc( vars_lh_sfc,                    & ! intent(in)
                                l_error,                        & ! intent(inout)
                                stats_metadata, stats_lh_sfc(i) ) ! intent(inout)
      end do

      fname = trim( stats_metadata%fname_lh_sfc )

      if ( stats_metadata%l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & ! In
                         1, stats_lh_sfc(1)%kk, nlat, nlon, stats_lh_sfc(1)%z, & ! In 
                         day, month, year, lat_vals, lon_vals, &  ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_lh_sfc(1)%num_output_fields, & ! In
                         stats_lh_sfc(1)%file ) ! intent(inout)

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_lh_sfc(1)%kk, &  ! In
                          stats_lh_sfc(1)%z, day, month, year, lat_vals, lon_vals, &  ! In
                          time_current, stats_metadata%stats_tout, & ! In
                          stats_lh_sfc(1)%num_output_fields, & ! In
                          stats_lh_sfc(1)%file, err_code ) ! InOut

        ! Finalize the variable definitions
        call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                          l_uv_nudge, & ! intent(in)
                          l_tke_aniso, & ! intent(in)
                          l_standard_term_ta, & ! intent(in)
                          stats_lh_sfc(1)%file, err_code ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error setting up stats_lh_sfc"
          return
        end if

#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

    end if ! l_silhs_out

    ! Initialize stats_zm (momentum points)

    ivar = 1
    do while ( ichar(vars_zm(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_zm(ivar)) /= 0 & 
               .and. ivar <= nvarmax_zm )
      ivar = ivar + 1
    end do
    ntot = ivar - 1

    if ( any( vars_zm == "hydrometp2" ) ) then
       ! Correct for number of variables found under "hydrometp2".
       ! Subtract "hydrometp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "wphydrometp" ) ) then
       ! Correct for number of variables found under "wphydrometp".
       ! Subtract "wphydrometp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "rtphmp" ) ) then
       ! Correct for number of variables found under "rtphmp".
       ! Subtract "rtphmp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "thlphmp" ) ) then
       ! Correct for number of variables found under "thlphmp".
       ! Subtract "thlphmp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "hmxphmyp" ) ) then
       ! Correct for number of variables found under "hmxphmyp".
       ! Subtract "hmxphmyp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add the number of overall covariances of two hydrometeors, which is
       ! found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zm statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 ) / 2
    endif

    if ( any( vars_zm == "K_hm" ) ) then
       ! Correct for number of variables found under "K_hm".
       ! Subtract "K_hm" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "sclrprtp" ) ) then
       ! Correct for number of variables found under "sclrprtp".
       ! Subtract "sclrprtp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif

    if ( any( vars_zm == "sclrp2" ) ) then
       ! Correct for number of variables found under "sclrp2".
       ! Subtract "sclrp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrpthvp" ) ) then
       ! Correct for number of variables found under "sclrpthvp".
       ! Subtract "sclrpthvp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrpthlp" ) ) then
       ! Correct for number of variables found under "sclrpthlp".
       ! Subtract "sclrpthlp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrprcp" ) ) then
       ! Correct for number of variables found under "sclrprcp".
       ! Subtract "sclrprcp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrp" ) ) then
       ! Correct for number of variables found under "wpsclrp".
       ! Subtract "wpsclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrp2" ) ) then
       ! Correct for number of variables found under "wpsclrp2".
       ! Subtract "wpsclrp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wp2sclrp" ) ) then
       ! Correct for number of variables found under "wp2sclrp".
       ! Subtract "wp2sclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrprtp" ) ) then
       ! Correct for number of variables found under "wpsclrprtp".
       ! Subtract "wpsclrprtp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrpthlp" ) ) then
       ! Correct for number of variables found under "wpsclrpthlp".
       ! Subtract "wpsclrpthlp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpedsclrp" ) ) then
       ! Correct for number of variables found under "wpedsclrp".
       ! Subtract "wpedsclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + edsclr_dim
    endif



    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      write(fstderr,*) "number of variables in vars_zm = ", ntot
      write(fstderr,*) "stats_init:  number of zm statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return
    end if

    ! Allocate and initialize for all columns, but we only write for 1
    do i = 1, ngrdcol
      stats_zm(i)%num_output_fields = ntot
      stats_zm(i)%kk = nzmax
      stats_zm(i)%ii = nlon
      stats_zm(i)%jj = nlat

      allocate( stats_zm(i)%z( stats_zm(i)%kk ) )
      stats_zm(i)%z = gzm(i,:)

      allocate( stats_zm(i)%accum_field_values( stats_zm(i)%ii, stats_zm(i)%jj, &
        stats_zm(i)%kk, stats_zm(i)%num_output_fields ) )
      allocate( stats_zm(i)%accum_num_samples( stats_zm(i)%ii, stats_zm(i)%jj, &
        stats_zm(i)%kk, stats_zm(i)%num_output_fields ) )
      allocate( stats_zm(i)%l_in_update( stats_zm(i)%ii, stats_zm(i)%jj, stats_zm(i)%kk, &
        stats_zm(i)%num_output_fields ) )

      call stats_zero( stats_zm(i)%ii, stats_zm(i)%jj, stats_zm(i)%kk, & ! In
                       stats_zm(i)%num_output_fields, & ! In
                       stats_zm(i)%accum_field_values, stats_zm(i)%accum_num_samples, & ! Out
                       stats_zm(i)%l_in_update ) ! Out

      allocate( stats_zm(i)%file%grid_avg_var( stats_zm(i)%num_output_fields ) )
      if ( l_different_output_grid ) then
        allocate( stats_zm(i)%file%z( output_gr_nzm ) )
      else
        allocate( stats_zm(i)%file%z( stats_zm(i)%kk ) )
      end if

      call stats_init_zm( hydromet_dim, sclr_dim, edsclr_dim, & ! intent(in)
                          hydromet_list, l_mix_rat_hm,        & ! intent(in)
                          vars_zm,                            & ! intent(in)
                          l_error,                            & ! intent(inout)
                          stats_metadata, stats_zm(i) )         ! intent(inout)
    end do

    fname = trim( stats_metadata%fname_zm )
    if ( stats_metadata%l_grads ) then

      ! Open GrADS files
      if ( l_different_output_grid ) then
        call open_grads( iunit, fdir, fname,  & ! In
                         1, output_gr_nzm, nlat, nlon, output_gr_zm(1,:), & ! In
                         day, month, year, lat_vals, lon_vals, & ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_zm(1)%num_output_fields, & ! In
                         stats_zm(1)%file ) ! intent(inout)
      else
        call open_grads( iunit, fdir, fname,  & ! In
                         1, stats_zm(1)%kk, nlat, nlon, stats_zm(1)%z, & ! In
                         day, month, year, lat_vals, lon_vals, & ! In
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_zm(1)%num_output_fields, & ! In
                         stats_zm(1)%file ) ! intent(inout)
      end if

    else ! Open NetCDF file
#ifdef NETCDF
      if ( l_different_output_grid ) then
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, output_gr_nzm, & !In
                                      output_gr_zm(1,:), day, month, year, lat_vals, lon_vals, & !In
                                      time_current, stats_metadata%stats_tout, & ! In
                                      stats_zm(1)%num_output_fields, & ! In
                                      stats_zm(1)%file, err_code ) ! InOut
      else
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_zm(1)%kk, & !In
                                      stats_zm(1)%z, day, month, year, lat_vals, lon_vals, & ! In
                                      time_current, stats_metadata%stats_tout, & ! In
                                      stats_zm(1)%num_output_fields, & ! In
                                      stats_zm(1)%file, err_code ) ! InOut
      end if
      
      ! Finalize the variable definitions
      call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                        l_uv_nudge, & ! intent(in)
                        l_tke_aniso, & ! intent(in)
                        l_standard_term_ta, & ! intent(in)
                        stats_zm(1)%file, err_code ) ! intent(inout)

      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Fatal error setting up stats_zm"
        return
      end if
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    ! Initialize stats_rad_zt (radiation points)

    if ( stats_metadata%l_output_rad_files ) then

      ivar = 1
      do while ( ichar(vars_rad_zt(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zt(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_rad_zt )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        write(fstderr,*) "number of variables in vars_rad_zt = ", ntot
        write(fstderr,*) "stats_init:  number of rad_zt statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      ! Allocate and initialize for all columns, but we only write for 1
      do i = 1, ngrdcol
        stats_rad_zt(i)%num_output_fields = ntot
        stats_rad_zt(i)%kk = nnrad_zt
        stats_rad_zt(i)%ii = nlon
        stats_rad_zt(i)%jj = nlat
        allocate( stats_rad_zt(i)%z( stats_rad_zt(i)%kk ) )
        stats_rad_zt(i)%z = grad_zt

        allocate( stats_rad_zt(i)%accum_field_values( stats_rad_zt(i)%ii, stats_rad_zt(i)%jj, &
          stats_rad_zt(i)%kk, stats_rad_zt(i)%num_output_fields ) )
        allocate( stats_rad_zt(i)%accum_num_samples( stats_rad_zt(i)%ii, stats_rad_zt(i)%jj, &
          stats_rad_zt(i)%kk, stats_rad_zt(i)%num_output_fields ) )
        allocate( stats_rad_zt(i)%l_in_update( stats_rad_zt(i)%ii, stats_rad_zt(i)%jj, &
          stats_rad_zt(i)%kk, stats_rad_zt(i)%num_output_fields ) )

        call stats_zero( stats_rad_zt(i)%ii, stats_rad_zt(i)%jj, stats_rad_zt(i)%kk, & ! intent(in)
                      stats_rad_zt(i)%num_output_fields, & ! intent(in)
                      stats_rad_zt(i)%accum_field_values, & ! intent(out)
                      stats_rad_zt(i)%accum_num_samples, stats_rad_zt(i)%l_in_update )! intent(out)

        allocate( stats_rad_zt(i)%file%grid_avg_var( stats_rad_zt(i)%num_output_fields ) )
        allocate( stats_rad_zt(i)%file%z( stats_rad_zt(i)%kk ) )

        call stats_init_rad_zt( vars_rad_zt,                    & ! intent(in)
                                l_error,                        & ! intent(inout)
                                stats_metadata, stats_rad_zt(i) ) ! intent(inout)
      end do

      fname = trim( stats_metadata%fname_rad_zt )
      if ( stats_metadata%l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & ! In
                         1, stats_rad_zt(1)%kk, nlat, nlon, stats_rad_zt(1)%z, & ! In
                         day, month, year, lat_vals, lon_vals, & 
                         time_current+real(stats_metadata%stats_tout, kind=time_precision), & ! In
                         stats_metadata, stats_rad_zt(1)%num_output_fields, & ! In
                         stats_rad_zt(1)%file ) ! In/Out

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname,  & ! intent(in)
                          1, stats_rad_zt(1)%kk, stats_rad_zt(1)%z, & ! intent(in)
                          day, month, year, lat_vals, lon_vals, & ! intent(in)
                          time_current, stats_metadata%stats_tout, & ! intent(in)
                          stats_rad_zt(1)%num_output_fields, & ! intent(in)
                          stats_rad_zt(1)%file, err_code ) ! intent(inout)

        ! Finalize the variable definitions
        call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                          l_uv_nudge, & ! intent(in)
                          l_tke_aniso, & ! intent(in)
                          l_standard_term_ta, & ! intent(in)
                          stats_rad_zt(1)%file, err_code ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error setting up stats_rad_zt"
          return
        end if
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

      ! Initialize stats_rad_zm (radiation points)

      ivar = 1
      do while ( ichar(vars_rad_zm(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zm(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_rad_zm )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        write(fstderr,*) "number of variables in vars_rad_zm = ", ntot
        write(fstderr,*) "stats_init:  number of rad_zm statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      ! Allocate and initialize for all columns, but we only write for 1
      do i = 1, ngrdcol
        stats_rad_zm(i)%num_output_fields = ntot
        stats_rad_zm(i)%kk = nnrad_zm
        stats_rad_zm(i)%ii = nlon
        stats_rad_zm(i)%jj = nlat

        allocate( stats_rad_zm(i)%z( stats_rad_zm(i)%kk ) )
        stats_rad_zm(i)%z = grad_zm

        allocate( stats_rad_zm(i)%accum_field_values( stats_rad_zm(i)%ii, stats_rad_zm(i)%jj, &
          stats_rad_zm(i)%kk, stats_rad_zm(i)%num_output_fields ) )
        allocate( stats_rad_zm(i)%accum_num_samples( stats_rad_zm(i)%ii, stats_rad_zm(i)%jj, &
          stats_rad_zm(i)%kk, stats_rad_zm(i)%num_output_fields ) )
        allocate( stats_rad_zm(i)%l_in_update( stats_rad_zm(i)%ii, stats_rad_zm(i)%jj, &
          stats_rad_zm(i)%kk, stats_rad_zm(i)%num_output_fields ) )

        call stats_zero( stats_rad_zm(i)%ii, stats_rad_zm(i)%jj, stats_rad_zm(i)%kk, & ! intent(in)
                    stats_rad_zm(i)%num_output_fields, & ! intent(in)
                    stats_rad_zm(i)%accum_field_values, & ! intent(out)
                    stats_rad_zm(i)%accum_num_samples, stats_rad_zm(i)%l_in_update ) ! intent(out)

        allocate( stats_rad_zm(i)%file%grid_avg_var( stats_rad_zm(i)%num_output_fields ) )
        allocate( stats_rad_zm(i)%file%z( stats_rad_zm(i)%kk ) )

        call stats_init_rad_zm( vars_rad_zm,                    & ! intent(in)
                                l_error,                        & ! intent(inout)
                                stats_metadata, stats_rad_zm(i) ) ! intent(inout)
      end do

      fname = trim( stats_metadata%fname_rad_zm )
      if ( stats_metadata%l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & ! In
                         1, stats_rad_zm(1)%kk, nlat, nlon, stats_rad_zm(1)%z, & ! In
                         day, month, year, lat_vals, lon_vals, & 
                         time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                         stats_metadata, stats_rad_zm(1)%num_output_fields, & ! In
                         stats_rad_zm(1)%file ) ! In/Out

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname,  & ! intent(in)
                          1, stats_rad_zm(1)%kk, stats_rad_zm(1)%z, & ! intent(in)
                          day, month, year, lat_vals, lon_vals, & ! intent(in)
                          time_current, stats_metadata%stats_tout, & ! intent(in)
                          stats_rad_zm(1)%num_output_fields, & ! intent(in)
                          stats_rad_zm(1)%file, err_code ) ! intent(inout)

        ! Finalize the variable definitions
        call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                          l_uv_nudge, & ! intent(in)
                          l_tke_aniso, & ! intent(in)
                          l_standard_term_ta, & ! intent(in)
                          stats_rad_zm(1)%file, err_code ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error setting up stats_rad_zm"
          return
        end if

#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

    end if ! l_output_rad_files


    ! Initialize stats_sfc (surface point)

    ivar = 1
    do while ( ichar(vars_sfc(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_sfc(ivar)) /= 0 & 
               .and. ivar <= nvarmax_sfc )
      ivar = ivar + 1
    end do
    ntot = ivar - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      write(fstderr,*) "number of variables in vars_sfc = ", ntot
      write(fstderr,*) "stats_init:  number of sfc statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return

    end if

    ! Allocate and initialize for all columns, but we only write for 1
    do i = 1, ngrdcol
      stats_sfc(i)%num_output_fields = ntot
      stats_sfc(i)%kk = 1
      stats_sfc(i)%ii = nlon
      stats_sfc(i)%jj = nlat

      allocate( stats_sfc(i)%z( stats_sfc(i)%kk ) )
      stats_sfc(i)%z = gzm(i,1)

      allocate( stats_sfc(i)%accum_field_values( stats_sfc(i)%ii, stats_sfc(i)%jj, &
        stats_sfc(i)%kk, stats_sfc(i)%num_output_fields ) )
      allocate( stats_sfc(i)%accum_num_samples( stats_sfc(i)%ii, stats_sfc(i)%jj, &
        stats_sfc(i)%kk, stats_sfc(i)%num_output_fields ) )
      allocate( stats_sfc(i)%l_in_update( stats_sfc(i)%ii, stats_sfc(i)%jj, &
        stats_sfc(i)%kk, stats_sfc(i)%num_output_fields ) )

      call stats_zero( stats_sfc(i)%ii, stats_sfc(i)%jj, stats_sfc(i)%kk, & ! In
                       stats_sfc(i)%num_output_fields, & ! In
                       stats_sfc(i)%accum_field_values, stats_sfc(i)%accum_num_samples, & ! Out
                       stats_sfc(i)%l_in_update ) ! Out

      allocate( stats_sfc(i)%file%grid_avg_var( stats_sfc(i)%num_output_fields ) )
      allocate( stats_sfc(i)%file%z( stats_sfc(i)%kk ) )

      call stats_init_sfc( vars_sfc,                    & ! intent(in)
                          l_error,                      & ! intent(inout)
                          stats_metadata, stats_sfc(i) )  ! intent(inout)
    end do

    fname = trim( stats_metadata%fname_sfc )

    if ( stats_metadata%l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & ! In
                       1, stats_sfc(1)%kk, nlat, nlon, stats_sfc(1)%z, & ! In
                       day, month, year, lat_vals, lon_vals, & ! In
                       time_current+real(stats_metadata%stats_tout,kind=time_precision), & ! In
                       stats_metadata, stats_sfc(1)%num_output_fields, & ! In
                       stats_sfc(1)%file ) ! In/Out

    else ! Open NetCDF files
#ifdef NETCDF
      call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_sfc(1)%kk, & ! In
                        stats_sfc(1)%z, day, month, year, lat_vals, lon_vals, & ! In
                        time_current, stats_metadata%stats_tout, & ! In
                        stats_sfc(1)%num_output_fields, & ! In
                        stats_sfc(1)%file, err_code ) ! InOut

      ! Finalize the variable definitions
      call first_write( clubb_params(1,:), sclr_dim, sclr_tol, & ! intent(in)
                        l_uv_nudge, & ! intent(in)
                        l_tke_aniso, & ! intent(in)
                        l_standard_term_ta, & ! intent(in)
                        stats_sfc(1)%file, err_code ) ! intent(inout)

      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Fatal error setting up stats_sfc"
        return
      end if
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      err_code = clubb_fatal_error
      return
    endif

    return

    ! If namelist was not found in input file, turn off statistics

    100 continue
    write(fstderr,*) 'Error with statsnl, statistics is turned off'
    stats_metadata%l_stats       = .false.
    stats_metadata%l_stats_samp  = .false.
    stats_metadata%l_stats_last  = .false.

    if ( err_code == clubb_fatal_error ) error stop

    return

  end subroutine stats_init_helper

  !-----------------------------------------------------------------------
  subroutine stats_init( iunit, fname_prefix, fdir, l_stats_in, &
                         stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                         hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, &
                         hydromet_list, l_mix_rat_hm, &
                         nzmax, ngrdcol, nlon, nlat, gzt, gzm, nnrad_zt, &
                         grad_zt, nnrad_zm, grad_zm, day, month, year, &
                         lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                         clubb_params, &
                         l_uv_nudge, &
                         l_tke_aniso, &
                         l_standard_term_ta, &
                         stats_metadata, &
                         stats_zt, stats_zm, stats_sfc, &
                         stats_lh_zt, stats_lh_sfc, &
                         stats_rad_zt, stats_rad_zm, &
                         err_code )
    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: &
        stats_metadata_type

    use parameter_indices, only: &
        nparams

    use clubb_precision, only: &
        time_precision, & ! Constant(s)
        core_rknd

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real( kind = core_rknd ), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      ngrdcol,  & ! Number of columns
      nlon,     & ! Number of points in the X direction [-]
      nlat,     & ! Number of points in the Y direction [-]
      nzmax       ! Grid points in the vertical         [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax-1) ::  & 
      gzt       ! Thermodynamic levels           [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) ::  & 
      gzm       ! Momentum levels                [m]

    integer, intent(in) :: &
      nnrad_zt,     & ! Grid points in the radiation grid [count]
      hydromet_dim, &
      sclr_dim,     &
      edsclr_dim

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    character(len=10), dimension(hydromet_dim), intent(in) :: & 
      hydromet_list

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt, stats_lh_sfc)  [boolean]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code    ! Error code catching and relaying any errors occurring in this subroutine

    ! Local Variables
    logical :: l_different_output_grid

    integer :: output_gr_nzm_placeholder

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      output_gr_zt_placeholder, &
      output_gr_zm_placeholder

    ! ------------------- Begin Code -------------------

    l_different_output_grid = .false.
    output_gr_nzm_placeholder = 2
    allocate( output_gr_zt_placeholder(ngrdcol,output_gr_nzm_placeholder-1) )
    allocate( output_gr_zm_placeholder(ngrdcol,output_gr_nzm_placeholder) )

    call stats_init_helper( iunit, fname_prefix, fdir, l_stats_in, &
                            stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                            hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, &
                            hydromet_list, l_mix_rat_hm, &
                            nzmax, ngrdcol, nlon, nlat, gzt, gzm, nnrad_zt, &
                            grad_zt, nnrad_zm, grad_zm, day, month, year, &
                            lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                            clubb_params, &
                            l_uv_nudge, &
                            l_tke_aniso, &
                            l_standard_term_ta, &
                            l_different_output_grid, &
                            output_gr_nzm_placeholder, &
                            output_gr_zt_placeholder, &
                            output_gr_zm_placeholder, &
                            stats_metadata, &
                            stats_zt, stats_zm, stats_sfc, &
                            stats_lh_zt, stats_lh_sfc, &
                            stats_rad_zt, stats_rad_zm, &
                            err_code )

    deallocate( output_gr_zt_placeholder )
    deallocate( output_gr_zm_placeholder )


  end subroutine stats_init

  !-----------------------------------------------------------------------
  subroutine stats_init_w_diff_output_gr( iunit, fname_prefix, fdir, l_stats_in, &
                                          stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                                          hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, &
                                          hydromet_list, l_mix_rat_hm, &
                                          nzmax, ngrdcol, nlon, nlat, gzt, gzm, nnrad_zt, &
                                          grad_zt, nnrad_zm, grad_zm, day, month, year, &
                                          lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                                          clubb_params, &
                                          l_uv_nudge, &
                                          l_tke_aniso, &
                                          l_standard_term_ta, &
                                          output_gr_nzm, &
                                          output_gr_zt, &
                                          output_gr_zm, &
                                          stats_metadata, &
                                          stats_zt, stats_zm, stats_sfc, &
                                          stats_lh_zt, stats_lh_sfc, &
                                          stats_rad_zt, stats_rad_zm, &
                                          err_code )
    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: &
        stats_metadata_type

    use parameter_indices, only: &
        nparams

    use clubb_precision, only: &
        time_precision, & ! Constant(s)
        core_rknd

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real( kind = core_rknd ), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      ngrdcol,  &    ! Number of columns
      nlon,     &    ! Number of points in the X direction [-]
      nlat,     &    ! Number of points in the Y direction [-]
      nzmax,    &    ! Grid points in the vertical         [-]
      output_gr_nzm  ! Grid points in the vertical for the grid that gets written to file;
                     ! can vary from nzmax since output is remapped to dycore grid if
                     ! we adapt grid and want to simulate forcings from dycore grid,
                     ! since we need a common grid for all results if the grid is adapted
                     ! over time, and the dycore grid can have a different number of
                     ! levels as the physics grid          [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax-1) ::  & 
      gzt       ! Thermodynamic levels           [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) ::  & 
      gzm       ! Momentum levels                [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,output_gr_nzm-1) ::  & 
      output_gr_zt  ! Thermodynamic levels of the grid that gets written to file    [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,output_gr_nzm) ::  & 
      output_gr_zm  ! Momentum levels of the grid that gets written to file         [m]
    
    ! Note: output_gr_zt and output_gr_zm are only different from gzt and gzm if we adapt the grid
    !       over time, since then we need a common grid we can remap all the results to,
    !       that were calculated on different grids

    integer, intent(in) :: &
      nnrad_zt,     & ! Grid points in the radiation grid [count]
      hydromet_dim, &
      sclr_dim,     &
      edsclr_dim

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    character(len=10), dimension(hydromet_dim), intent(in) :: & 
      hydromet_list

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt, stats_lh_sfc)  [boolean]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code    ! Error code catching and relaying any errors occurring in this subroutine

    ! Local variables
    logical :: &
      l_different_output_grid ! use different grid to output values to file

    ! ------------------- Begin Code -------------------

    l_different_output_grid = .true.

    call stats_init_helper( iunit, fname_prefix, fdir, l_stats_in, &
                            stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                            hydromet_dim, sclr_dim, edsclr_dim, sclr_tol, &
                            hydromet_list, l_mix_rat_hm, &
                            nzmax, ngrdcol, nlon, nlat, gzt, gzm, nnrad_zt, &
                            grad_zt, nnrad_zm, grad_zm, day, month, year, &
                            lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                            clubb_params, &
                            l_uv_nudge, &
                            l_tke_aniso, &
                            l_standard_term_ta, &
                            l_different_output_grid, &
                            output_gr_nzm, &
                            output_gr_zt, &
                            output_gr_zm, &
                            stats_metadata, &
                            stats_zt, stats_zm, stats_sfc, &
                            stats_lh_zt, stats_lh_sfc, &
                            stats_rad_zt, stats_rad_zm, &
                            err_code )

  end subroutine stats_init_w_diff_output_gr

  !-----------------------------------------------------------------------
  subroutine stats_zero( ii, jj, kk, nn, &
                         x, n, l_in_update )

    ! Description:
    !   Initialize stats to zero
    ! References:
    !   None
    !-----------------------------------------------------------------------
    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: ii, jj, kk, nn

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(ii,jj,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(ii,jj,kk,nn), intent(out) :: n
    logical, dimension(ii,jj,kk,nn), intent(out) :: l_in_update

    ! Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0_stat_rknd
      n(:,:,:,:) = 0_stat_nknd
      l_in_update(:,:,:,:) = .false.
    end if

    return
  end subroutine stats_zero

  !-----------------------------------------------------------------------
  subroutine stats_avg( ii, jj, kk, nn, n, &
                        x )

    ! Description:
    !   Compute the average of stats fields
    ! References:
    !   None
    !-----------------------------------------------------------------------
    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    use stat_file_module, only: &
        clubb_i, clubb_j ! Variable(s)

    implicit none

    ! External
    intrinsic :: real

    ! Input Variable(s)
    integer, intent(in) :: &
      ii, & ! Number of points in X (i.e. latitude) dimension
      jj, & ! Number of points in Y (i.e. longitude) dimension
      kk, & ! Number of levels in vertical (i.e. Z) dimension
      nn    ! Number of variables being output to disk (e.g. cloud_frac, rain rate, etc.)

    integer(kind=stat_nknd), dimension(ii,jj,kk,nn), intent(in) :: &
      n ! n is the number of samples for each of the nn fields 
        ! and each of the kk vertical levels

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(ii,jj,kk,nn), intent(inout) :: &
      x ! The variable x contains the cumulative sums of n sample values of each of
        ! the nn output fields (e.g. the sum of the sampled rain rate values)

    ! ---- Begin Code ----

    ! Compute averages
    where ( n(1,1,1:kk,1:nn) > 0 )
      x(clubb_i,clubb_j,1:kk,1:nn) = x(clubb_i,clubb_j,1:kk,1:nn) &
         / real( n(clubb_i,clubb_j,1:kk,1:nn), kind=stat_rknd )
    end where

    return
  end subroutine stats_avg

  !-----------------------------------------------------------------------
  subroutine stats_begin_timestep( itime, stats_nsamp, stats_nout, &
                                   stats_metadata )

    !     Description:
    !       Given the elapsed time, set flags determining specifics such as
    !       if this time set should be sampled or if this is the first or
    !       last time step.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
      stats_metadata_type

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) ::  & 
      itime, &       ! Elapsed model time       [timestep]
      stats_nsamp, & ! Stats sampling interval  [timestep]
      stats_nout     ! Stats output interval    [timestep]

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    if ( .not. stats_metadata%l_stats ) return

    ! Only sample time steps that are multiples of "stats_tsamp"
    ! in a case's "model.in" file to shorten length of run
    if ( mod( itime, stats_nsamp ) == 0 ) then
      stats_metadata%l_stats_samp = .true.
    else
      stats_metadata%l_stats_samp = .false.
    end if

    ! Indicates the end of the sampling time period. Signals to start writing to the file
    if ( mod( itime, stats_nout ) == 0 ) then
      stats_metadata%l_stats_last = .true.
    else
      stats_metadata%l_stats_last = .false.
    end if

    return

  end subroutine stats_begin_timestep

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep_helper( stats_metadata,                & ! intent(in)
                                        l_different_output_grid,       & ! intent(in)
                                        gr_source, gr_target,          & ! intent(in)
                                        total_idx_rho_lin_spline,      & ! intent(in)
                                        rho_lin_spline_vals,           & ! intent(in)
                                        rho_lin_spline_levels,         & ! intent(in)
                                        p_sfc,                         & ! intent(in)
                                        grid_remap_method,             & ! intent(in)
                                        stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                                        stats_lh_zt, stats_lh_sfc,     & ! intent(inout)
                                        stats_rad_zt, stats_rad_zm,    & ! intent(inout)
                                        err_code                       & ! intent(inout)
                                      )

    ! Description:
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Constant(s)

    use stats_variables, only: & 
        stats_metadata_type

    use output_grads, only: &
        write_grads ! Procedure(s)

    use stat_file_module, only: &
        clubb_i, & ! Variable(s)
        clubb_j

#ifdef NETCDF
    use output_netcdf, only: &
        write_netcdf_w_diff_output_gr, &
        write_netcdf ! Procedure(s)
#endif

    use error_code, only : &
        clubb_fatal_error   ! Constant

    use stats_type, only: stats ! Type

    use grid_class, only: grid ! Type

    implicit none

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    logical, intent(in) :: &
      l_different_output_grid

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of points in the rho spline

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, &  ! the rho values for constructing the spline for remapping;
                              ! only used if l_different_output_gr is .true.
      rho_lin_spline_levels   ! the levels at which the rho values are given;
                              ! only used if l_different_output_gr is .true.

    type( grid ), intent(in) :: &
      gr_source, & ! the grid where the values are currently given on;
                   ! only used if l_different_output_gr is .true.
      gr_target    ! the grid where the values should be remapped to;
                   ! only used if l_different_output_gr is .true.

    real( kind = core_rknd ), dimension(1), intent(in) :: &
      p_sfc

    integer, intent(in) :: &
      grid_remap_method

    type (stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code  ! Error code catching and relaying any errors occurring in this subroutine

    ! External
    intrinsic :: floor

    ! Local Variables

    logical :: l_error, &
               l_zt_variable

    ! ------------------- Begin Code -------------------

    ! Check if it is time to write to file

    if ( .not. stats_metadata%l_stats_last ) return

    ! Initialize
    l_error = .false.

    call stats_check_num_samples( stats_zt, stats_metadata,  & ! intent(in)
                                  l_error )                    ! intent(inout)
    call stats_check_num_samples( stats_zm, stats_metadata,  & ! intent(in)
                                  l_error )                    ! intent(inout)
    call stats_check_num_samples( stats_sfc, stats_metadata, & ! intent(in)
                                  l_error )                    ! intent(inout)
    if ( stats_metadata%l_silhs_out ) then
      call stats_check_num_samples( stats_lh_zt, stats_metadata,  & ! intent(in)
                                    l_error )                       ! intent(inout)
      call stats_check_num_samples( stats_lh_sfc, stats_metadata, & ! intent(in)
                                    l_error )                       ! intent(inout)
    end if
    if ( stats_metadata%l_output_rad_files ) then
      call stats_check_num_samples( stats_rad_zt, stats_metadata, & ! intent(in)
                                    l_error )                       ! intent(inout)
      call stats_check_num_samples( stats_rad_zm, stats_metadata, & ! intent(in)
                                    l_error )                       ! intent(inout)
    end if

    ! Return if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      write(fstderr,*) 'stats_end_timestep:  error(s) found'
      err_code = clubb_fatal_error
      return
    end if ! l_error

    ! Compute averages
    call stats_avg( stats_zt%ii, stats_zt%jj, stats_zt%kk, stats_zt%num_output_fields, & ! In
                    stats_zt%accum_num_samples, & ! intent(in)
                    stats_zt%accum_field_values ) ! intent(inout)
    call stats_avg( stats_zm%ii, stats_zm%jj, stats_zm%kk, stats_zm%num_output_fields, & ! In
                    stats_zm%accum_num_samples, & ! intent(in)
                    stats_zm%accum_field_values ) ! intent(inout)
    if ( stats_metadata%l_silhs_out ) then
      call stats_avg( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, & ! intent(in)
                      stats_lh_zt%num_output_fields, stats_lh_zt%accum_num_samples, & ! intent(in)
                      stats_lh_zt%accum_field_values ) ! intent(inout)
      call stats_avg( stats_lh_sfc%ii, stats_lh_sfc%jj, stats_lh_sfc%kk, & ! intent(in)
                      stats_lh_sfc%num_output_fields, stats_lh_sfc%accum_num_samples, & ! intent(in)
                      stats_lh_sfc%accum_field_values ) ! intent(inout)
    end if
    if ( stats_metadata%l_output_rad_files ) then
      call stats_avg( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zt%kk, & ! intent(in)
                      stats_rad_zt%num_output_fields, & ! intent(in)
                      stats_rad_zt%accum_num_samples, & ! intent(in)
                      stats_rad_zt%accum_field_values ) ! intent(inout)
      call stats_avg( stats_rad_zm%ii, stats_rad_zm%jj, stats_rad_zm%kk, & ! intent(in)
                      stats_rad_zm%num_output_fields, & ! intent(in)
                      stats_rad_zm%accum_num_samples, & ! intent(in)
                      stats_rad_zm%accum_field_values ) ! intent(inout)
    end if
    call stats_avg( stats_sfc%ii, stats_sfc%jj, stats_sfc%kk, stats_sfc%num_output_fields, & ! In
                    stats_sfc%accum_num_samples, & ! intent(in)
                    stats_sfc%accum_field_values ) ! intent(inout)

    ! Only write to the file and zero out the stats fields if we've reach the horizontal
    ! limits of the domain (this is always true in the single-column case because it's 1x1).
    if ( clubb_i == stats_zt%ii .and. clubb_j == stats_zt%jj ) then
      ! Write to file
      if ( stats_metadata%l_grads ) then
        call write_grads( stats_zt%file ) ! intent(inout)
        call write_grads( stats_zm%file ) ! intent(inout)
        if ( stats_metadata%l_silhs_out ) then
          call write_grads( stats_lh_zt%file ) ! intent(inout)
          call write_grads( stats_lh_sfc%file ) ! intent(inout)
        end if
        if ( stats_metadata%l_output_rad_files ) then
          call write_grads( stats_rad_zt%file ) ! intent(inout)
          call write_grads( stats_rad_zm%file ) ! intent(inout)
        end if
        call write_grads( stats_sfc%file ) ! intent(inout)
      else ! l_netcdf

#ifdef NETCDF
        if ( l_different_output_grid ) then
          l_zt_variable = .true.
          call write_netcdf_w_diff_output_gr( gr_source, gr_target, &      ! Intent(in)
                                              l_zt_variable, &             ! Intent(in)
                                              total_idx_rho_lin_spline, &  ! Intent(in)
                                              rho_lin_spline_vals, &       ! Intent(in)
                                              rho_lin_spline_levels, &     ! Intent(in)
                                              p_sfc, &
                                              grid_remap_method, &
                                              stats_zt%file, err_code )    ! Intent(inout)
          
          l_zt_variable = .false.
          call write_netcdf_w_diff_output_gr( gr_source, gr_target, &      ! Intent(in)
                                              l_zt_variable, &             ! Intent(in)
                                              total_idx_rho_lin_spline, &  ! Intent(in)
                                              rho_lin_spline_vals, &       ! Intent(in)
                                              rho_lin_spline_levels, &     ! Intent(in)
                                              p_sfc, &
                                              grid_remap_method, &
                                              stats_zm%file, err_code )    ! Intent(inout)
        else
          call write_netcdf( stats_zt%file, err_code ) ! intent(inout)

          call write_netcdf( stats_zm%file, err_code ) ! intent(inout)
        endif

        if ( stats_metadata%l_silhs_out ) then

          call write_netcdf( stats_lh_zt%file, err_code ) ! intent(inout)

          call write_netcdf( stats_lh_sfc%file, err_code ) ! intent(inout)

        end if
        if ( stats_metadata%l_output_rad_files ) then

          call write_netcdf( stats_rad_zt%file, err_code ) ! intent(inout)

          call write_netcdf( stats_rad_zm%file, err_code ) ! intent(inout)

        end if

        call write_netcdf( stats_sfc%file, err_code ) ! intent(inout)
        
        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This program was not compiled with netCDF support"
#endif /* NETCDF */
      end if ! l_grads

      ! Reset sample fields
      call stats_zero( stats_zt%ii, stats_zt%jj, stats_zt%kk, stats_zt%num_output_fields, & ! In
      stats_zt%accum_field_values, stats_zt%accum_num_samples, stats_zt%l_in_update ) ! out
      call stats_zero( stats_zm%ii, stats_zm%jj, stats_zm%kk, stats_zm%num_output_fields, & ! In
        stats_zm%accum_field_values, stats_zm%accum_num_samples, stats_zm%l_in_update ) ! Out
      if ( stats_metadata%l_silhs_out ) then
        call stats_zero( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, & ! intent(in)
          stats_lh_zt%num_output_fields, & ! intent(in)
          stats_lh_zt%accum_field_values, & ! intent(out)
          stats_lh_zt%accum_num_samples, stats_lh_zt%l_in_update ) ! intent(out)
        call stats_zero( stats_lh_sfc%ii, stats_lh_sfc%jj, stats_lh_sfc%kk, & ! intent(in)
          stats_lh_sfc%num_output_fields, & ! intent(in)
          stats_lh_sfc%accum_field_values, & ! intent(out)
          stats_lh_sfc%accum_num_samples, stats_lh_sfc%l_in_update ) ! intent(out)
      end if
      if ( stats_metadata%l_output_rad_files ) then
        call stats_zero( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zt%kk, & ! intent(in)
          stats_rad_zt%num_output_fields, & ! intent(in)
          stats_rad_zt%accum_field_values, & ! intent(out)
          stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update ) ! intent(out)
        call stats_zero( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zm%kk, & ! intent(in)
          stats_rad_zm%num_output_fields, & ! intent(in)
          stats_rad_zm%accum_field_values, & ! intent(out)
          stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update ) ! intent(out)
      end if
      call stats_zero( stats_sfc%ii, stats_sfc%jj, stats_sfc%kk, stats_sfc%num_output_fields, & !IN
        stats_sfc%accum_field_values, & ! intent(out)
        stats_sfc%accum_num_samples, stats_sfc%l_in_update ) ! intent(out)

    end if ! clubb_i = stats_zt%ii .and. clubb_j == stats_zt%jj

    return
  end subroutine stats_end_timestep_helper

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep( stats_metadata,                & ! intent(in)
                                 stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                                 stats_lh_zt, stats_lh_sfc,     & ! intent(inout)
                                 stats_rad_zt, stats_rad_zm,    & ! intent(inout)
                                 err_code                       & ! intent(inout)
                               )

    ! Description:
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        stats_metadata_type

    use stats_type, only: stats ! Type

    use grid_class, only: grid ! Type

    use clubb_precision, only: &
        core_rknd ! Constant(s)

    implicit none

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type (stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code  ! Error code catching and relaying any errors occurring in this subroutine

    ! External
    intrinsic :: floor

    ! Local Variables

    logical :: l_different_output_grid

    integer :: &
      total_idx_rho_lin_spline_placeholder, & ! number of points in the rho spline
      grid_remap_method_placeholder

    real( kind = core_rknd ), dimension(:), allocatable :: &
      rho_lin_spline_vals_placeholder, &  ! the rho values for constructing the spline for
                                          ! remapping; only used if l_different_output_gr is .true.
      rho_lin_spline_levels_placeholder   ! the levels at which the rho values are given;
                                          ! only used if l_different_output_gr is .true.

    type( grid ) :: &
      gr_source_placeholder, & ! the grid where the values are currently given on;
                               ! only used if l_different_output_gr is .true.
      gr_target_placeholder    ! the grid where the values should be remapped to;
                               ! only used if l_different_output_gr is .true.

    real( kind = core_rknd ), dimension(1) :: &
      p_sfc_placeholder

    ! ------------------- Begin Code -------------------
    l_different_output_grid = .false.
    total_idx_rho_lin_spline_placeholder = 1
    p_sfc_placeholder = -9999.0
    grid_remap_method_placeholder = 0

    allocate( rho_lin_spline_vals_placeholder(total_idx_rho_lin_spline_placeholder) )
    allocate( rho_lin_spline_levels_placeholder(total_idx_rho_lin_spline_placeholder) )

    call stats_end_timestep_helper( stats_metadata,                               & ! intent(in)
                                    l_different_output_grid,                      & ! intent(in)
                                    gr_source_placeholder, gr_target_placeholder, & ! intent(in)
                                    total_idx_rho_lin_spline_placeholder,         & ! intent(in)
                                    rho_lin_spline_vals_placeholder,              & ! intent(in)
                                    rho_lin_spline_levels_placeholder,            & ! intent(in)
                                    p_sfc_placeholder,                            & ! intent(in)
                                    grid_remap_method_placeholder,                & ! intent(in)
                                    stats_zt, stats_zm, stats_sfc,                & ! intent(inout)
                                    stats_lh_zt, stats_lh_sfc,                    & ! intent(inout)
                                    stats_rad_zt, stats_rad_zm,                   & ! intent(inout)
                                    err_code                                      & ! intent(inout)
                                  )

    deallocate( rho_lin_spline_vals_placeholder )
    deallocate( rho_lin_spline_levels_placeholder )
    
  end subroutine stats_end_timestep

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep_w_diff_output_gr( stats_metadata,                & ! intent(in)
                                                  gr_source, gr_target,          & ! intent(in)
                                                  total_idx_rho_lin_spline,      & ! intent(in)
                                                  rho_lin_spline_vals,           & ! intent(in)
                                                  rho_lin_spline_levels,         & ! intent(in)
                                                  p_sfc,                         & ! intent(in)
                                                  grid_remap_method,             & ! intent(in)
                                                  stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                                                  stats_lh_zt, stats_lh_sfc,     & ! intent(inout)
                                                  stats_rad_zt, stats_rad_zm,    & ! intent(inout)
                                                  err_code                       & ! intent(inout)
                                                )

    ! Description:
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        stats_metadata_type

    use stats_type, only: stats ! Type

    use grid_class, only: grid ! Type

    use clubb_precision, only: &
        core_rknd ! Constant(s)

    implicit none

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    integer, intent(in) :: &
      total_idx_rho_lin_spline ! number of points in the rho spline

    real( kind = core_rknd ), dimension(total_idx_rho_lin_spline), intent(in) :: &
      rho_lin_spline_vals, &  ! the rho values for constructing the spline for remapping;
                              ! only used if l_different_output_gr is .true.
      rho_lin_spline_levels   ! the levels at which the rho values are given;
                              ! only used if l_different_output_gr is .true.

    type( grid ), intent(in) :: &
      gr_source, & ! the grid where the values are currently given on;
                   ! only used if l_different_output_gr is .true.
      gr_target    ! the grid where the values should be remapped to;
                   ! only used if l_different_output_gr is .true.

    real( kind = core_rknd ), dimension(1), intent(in) :: &
      p_sfc

    integer, intent(in) :: &
      grid_remap_method

    type (stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    integer, intent(inout) :: &
      err_code  ! Error code catching and relaying any errors occurring in this subroutine

    ! External
    intrinsic :: floor

    ! Local Variables

    logical :: l_different_output_grid

    ! ------------------- Begin Code -------------------
    l_different_output_grid = .true.

    call stats_end_timestep_helper( stats_metadata,                 & ! intent(in)
                                    l_different_output_grid,        & ! intent(in)
                                    gr_source, gr_target,           & ! intent(in)
                                    total_idx_rho_lin_spline,       & ! intent(in)
                                    rho_lin_spline_vals,            & ! intent(in)
                                    rho_lin_spline_levels,          & ! intent(in)
                                    p_sfc,                          & ! intent(in)
                                    grid_remap_method,              & ! intent(in)
                                    stats_zt, stats_zm, stats_sfc,  & ! intent(inout)
                                    stats_lh_zt, stats_lh_sfc,      & ! intent(inout)
                                    stats_rad_zt, stats_rad_zm,     & ! intent(inout)
                                    err_code                        & ! intent(inout)
                                  )
    
  end subroutine stats_end_timestep_w_diff_output_gr

  !----------------------------------------------------------------------
  subroutine stats_accumulate( &
                     nzm, nzt, sclr_dim, edsclr_dim, &
                     invrs_dzm, zt, dzm, dzt, dt, &
                     um, vm, upwp, vpwp, up2, vp2, &
                     thlm, rtm, wprtp, wpthlp, &
                     wp2, wp3, rtp2, rtp3, thlp2, thlp3, &
                     rtpthlp, &
                     wpthvp, wp2thvp, rtpthvp, thlpthvp, &
                     p_in_Pa, exner, rho, rho_zm, &
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                     wm_zt, wm_zm, rcm, wprcp, rc_coef, &
                     rc_coef_zm, &
                     rcm_zm, rtm_zm, thlm_zm, cloud_frac, &
                     ice_supersat_frac, &
                     cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
                     cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
                     thvm, ug, vg, Lscale, wpthlp2, wp2thlp, &
                     wprtp2, wp2rtp, &
                     Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
                     wprtpthlp, sigma_sqd_w_zt, rsat, wp2_zt, &
                     thlp2_zt, &
                     wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, &
                     up2_zt, &
                     vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, & 
                     wp2up2, wp2vp2, wp4, &
                     tau_zm, Kh_zm, thlprcp, &
                     rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
                     wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
                     w_up_in_cloud, w_down_in_cloud, &
                     cloudy_updraft_frac, cloudy_downdraft_frac, &
                     pdf_params, pdf_params_zm, &
                     sclrm, sclrp2, &
                     sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp, &
                     wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, &
                     wpsclrprtp, &
                     wpsclrpthlp, wpedsclrp, edsclrm, &
                     edsclrm_forcing, &
                     saturation_formula, &
                     stats_metadata, &
                     stats_zt, stats_zm, stats_sfc )

    ! Description:
    !   Accumulate those stats variables that are preserved in CLUBB from timestep to
    !   timestep, but not those stats that are not, (e.g. budget terms, longwave and
    !   shortwave components, etc.)
    !
    ! References:
    !   None
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        cloud_frac_min, &  ! Constant
        eps

    use pdf_utilities, only: &
        compute_variance_binormal    ! Procedure

    use stats_variables, only: & 
        stats_metadata_type

    use pdf_parameter_module, only: & 
        pdf_parameter ! Type

    use T_in_K_module, only: & 
        thlm2T_in_K ! Procedure

    use constants_clubb, only: & 
        rc_tol, fstderr    ! Constant(s)

    use stats_type_utilities, only: & 
        stat_update_var,  & ! Procedure(s)
        stat_update_var_pt

    use advance_helper_module, only: &
        vertical_avg, &     ! Procedure(s)
        vertical_integral

    use interpolation, only: & 
        lin_interpolate_two_points             ! Procedure

    use saturation, only: &
        sat_mixrat_ice ! Procedure

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      nzm, &
      nzt, &
      sclr_dim, &
      edsclr_dim
    
    real( kind = core_rknd ), intent(in), dimension(nzm) :: & 
      invrs_dzm, & ! The inverse spacing between thermodynamic grid
                   ! levels; centered over momentum grid levels.
      dzm          ! Spacing between thermodynamic grid levels; centered over
                   ! momentum grid levels

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      zt,        & ! Thermo grid
      dzt          ! Spcaing between momentum grid levels; centered over
                   ! thermodynamic grid levels

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Model timestep                        [s]

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      um,       & ! u wind (thermodynamic levels)          [m/s]
      vm,       & ! v wind (thermodynamic levels)          [m/s]
      thlm,     & ! liquid potential temperature (t levs.) [K]
      rtm,      & ! total water mixing ratio (t levs.)     [kg/kg]
      wp3,      & ! < w'^3 > (thermodynamic levels)        [m^3/s^3]
      rtp3,     & ! < r_t'^3 > (thermodynamic levels)      [(kg/kg)^3]
      thlp3,    & ! < th_l'^3 > (thermodynamic levels)     [K^3]
      wp2thvp     ! < w'^2 th_v' > (thermodynamic levels)  [m^2/s^2 K]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: & 
      upwp,     & ! vertical u momentum flux (m levs.)     [m^2/s^2]
      vpwp,     & ! vertical v momentum flux (m levs.)     [m^2/s^2]
      up2,      & ! < u'^2 > (momentum levels)             [m^2/s^2]
      vp2,      & ! < v'^2 > (momentum levels)             [m^2/s^2]
      wprtp,    & ! < w' r_t' > (momentum levels)          [m/s kg/kg]
      wpthlp,   & ! < w' th_l' > (momentum levels)         [m/s K]
      wp2,      & ! < w'^2 > (momentum levels)             [m^2/s^2]
      rtp2,     & ! < r_t'^2 > (momentum levels)           [(kg/kg)^2]
      thlp2,    & ! < th_l'^2 > (momentum levels)          [K^2]
      rtpthlp,  & ! < r_t' th_l' > (momentum levels)       [kg/kg K]
      wpthvp,   & ! < w' th_v' > (momentum levels)         [kg/kg K]
      rtpthvp,  & ! < r_t' th_v' > (momentum levels)       [kg/kg K]
      thlpthvp    ! < th_l' th_v' > (momentum levels)      [K^2]

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
      exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
      rho,          & ! Density (thermodynamic levels)           [kg/m^3]
      rho_ds_zt,    & ! Dry, static density (thermo. levs.)      [kg/m^3]
      thv_ds_zt,    & ! Dry, base-state theta_v (thermo. levs.)  [K]
      wm_zt           ! w on thermodynamic levels                [m/s]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: & 
      rho_zm,       & ! Density on momentum levels               [kg/m^3]
      rho_ds_zm,    & ! Dry, static density (momentum levels)    [kg/m^3]
      thv_ds_zm,    & ! Dry, base-state theta_v (momentum levs.) [K]
      wm_zm           ! w on momentum levels                     [m/s]

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      rcm,                  & ! Cloud water mixing ratio (t levs.)       [kg/kg]
      rc_coef,              & ! Coefficient of X'r_c' (t-levs.)      [K/(kg/kg)]
      cloud_frac,           & ! Cloud fraction (thermodynamic levels)    [-]
      ice_supersat_frac,    & ! Ice cloud fracion (thermodynamic levels) [-]
      rcm_in_layer,         & ! Cloud water mixing ratio in cloud layer  [kg/kg]
      cloud_cover,          & ! Cloud cover                              [-]
      rcm_supersat_adj        ! rcm adjustment due to supersaturation    [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: & 
      rcm_zm,               & ! Cloud water mixing ratio on m levs.      [kg/kg]
      rtm_zm,               & ! Total water mixing ratio on m levs.      [kg/kg]
      thlm_zm,              & ! Liquid potential temperature on m levs.  [K]
      wprcp,                & ! < w' r_c' > (momentum levels)            [m/s kg/kg]
      rc_coef_zm,           & ! Coefficient of X'r_c' on m-levs.     [K/(kg/kg)]
      cloud_frac_zm,        & ! Cloud fraction on zm levels              [-]
      ice_supersat_frac_zm    ! Ice cloud fraction on zm levels          [-]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: &
      sigma_sqd_w    ! PDF width parameter (momentum levels)    [-]

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      thvm,           & ! Virtual potential temperature        [K]
      ug,             & ! u geostrophic wind                   [m/s]
      vg,             & ! v geostrophic wind                   [m/s]
      Lscale,         & ! Length scale                         [m]
      wpthlp2,        & ! w'thl'^2                             [m K^2/s]
      wp2thlp,        & ! w'^2 thl'                            [m^2 K/s^2]
      wprtp2,         & ! w'rt'^2                              [m/s kg^2/kg^2]
      wp2rtp,         & ! w'^2rt'                              [m^2/s^2 kg/kg]
      Lscale_up,      & ! Length scale (upwards component)     [m]
      Lscale_down,    & ! Length scale (downwards component)   [m]
      tau_zt,         & ! Eddy diss. time scale; thermo. levs. [s]
      Kh_zt,          & ! Eddy diff. coef. on thermo. levels   [m^2/s]
      wp2rcp,         & ! w'^2 rc'                             [m^2/s^2 kg/kg]
      wprtpthlp,      & ! w'rt'thl'                            [m/s kg/kg K]
      sigma_sqd_w_zt, & ! PDF width parameter (thermo. levels) [-]
      rsat              ! Saturation mixing ratio              [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(nzt) :: & 
      wp2_zt,                & ! w'^2 on thermo. grid                  [m^2/s^2]
      thlp2_zt,              & ! thl'^2 on thermo. grid                [K^2]
      wpthlp_zt,             & ! w'thl' on thermo. grid                [m K/s]
      wprtp_zt,              & ! w'rt' on thermo. grid                 [m kg/(kg s)]
      rtp2_zt,               & ! rt'^2 on therm. grid                  [(kg/kg)^2]
      rtpthlp_zt,            & ! rt'thl' on thermo. grid               [kg K/kg]
      up2_zt,                & ! u'^2 on thermo. grid                  [m^2/s^2]
      vp2_zt,                & ! v'^2 on thermo. grid                  [m^2/s^2]
      upwp_zt,               & ! u'w' on thermo. grid                  [m^2/s^2]
      vpwp_zt,               & ! v'w' on thermo. grid                  [m^2/s^2]
      wpup2,                 & ! w'u'^2 (thermodynamic levels)         [m^3/s^3]
      wpvp2,                 & ! w'v'^2 (thermodynamic levels)         [m^3/s^3]
      a3_coef_zt,            & ! The a3 coef. interp. to the zt grid   [-]
      wp3_on_wp2_zt,         & ! w'^3 / w'^2 on the zt grid            [m/s]
      w_up_in_cloud,         & ! Mean cloudy updraft speed             [m/s]
      w_down_in_cloud,       & ! Mean cloudy downdraft speed           [m/s]
      cloudy_updraft_frac,   & ! Cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! Cloudy downdraft fraction             [-]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: & 
      wp2up2,                & ! < w'^2u'^2 > (momentum levels)        [m^4/s^4]
      wp2vp2,                & ! < w'^2v'^2 > (momentum levels)        [m^4/s^4]
      wp4,                   & ! < w'^4 > (momentum levels)            [m^4/s^4]
      tau_zm,                & ! Eddy diss. time scale; momentum levs. [s]
      Kh_zm,                 & ! Eddy diff. coef. on momentum levels   [m^2/s]
      thlprcp,               & ! thl'rc'                               [K kg/kg]
      rtprcp,                & ! rt'rc'                                [kg^2/kg^2]
      rcp2,                  & ! rc'^2                                 [kg^2/kg^2]
      em,                    & ! Turbulent Kinetic Energy (TKE)        [m^2/s^2]
      a3_coef,               & ! The a3 coefficient from CLUBB eqns    [-]
      wp3_zm,                & ! w'^3 interpolated to momentum levels  [m^3/s^3]
      wp3_on_wp2,            & ! w'^3 / w'^2 on the zm grid            [m/s]
      Skw_velocity             ! Skewness velocity                     [m/s]

    type(pdf_parameter), intent(in) :: & 
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    real( kind = core_rknd ), intent(in), dimension(nzt,sclr_dim) :: & 
      sclrm,         & ! High-order passive scalar            [units vary]
      sclrm_forcing, & ! Large-scale forcing of scalar        [units/s]
      wp2sclrp,      & ! w'^2 sclr'                           [units vary]
      wpsclrp2,      & ! w'sclr'^2                            [units vary]
      wpsclrprtp,    & ! w'sclr'rt'                           [units vary]
      wpsclrpthlp      ! w'sclr'thl'      [units vary]

    real( kind = core_rknd ), intent(in), dimension(nzm,sclr_dim) :: & 
      sclrprcp,  & ! sclr'rc'                                [units vary]
      sclrp2,    & ! High-order passive scalar variance      [units^2]
      sclrprtp,  & ! High-order passive scalar covariance    [units kg/kg]
      sclrpthlp, & ! High-order passive scalar covariance    [units K]
      sclrpthvp, & ! High-order passive scalar covariance    [units K]
      wpsclrp      ! w'sclr'                                 [units m/s]

    real( kind = core_rknd ), intent(in), dimension(nzt,edsclr_dim) :: & 
      edsclrm,         & ! Eddy-diff passive scalar         [units vary] 
      edsclrm_forcing    ! Large-scale forcing of edscalar  [units vary]

    real( kind = core_rknd ), intent(in), dimension(nzm,edsclr_dim) :: & 
      wpedsclrp          ! w'edsclr'                        [units vary]

    integer, intent(in) :: &
      saturation_formula

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type (stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! Local Variables

    integer :: sclr, edsclr, k
    integer :: grid_level = 1  ! grid level for stats where there is only one sensible level 
                               ! (eg timeseries)

    real( kind = core_rknd ), dimension(nzt) :: &
      T_in_K,        &  ! Absolute temperature         [K]
      rsati,         &  ! Saturation w.r.t ice         [kg/kg]
      chi,           &  ! Mellor's 's'                 [kg/kg]
      chip2,         &  ! Variance of Mellor's 's'     [kg/kg]
      rcm_in_cloud      ! rcm in cloud                 [kg/kg]

    real( kind = core_rknd ), dimension(nzm) :: &
      shear       ! Wind shear production term   [m^2/s^3]

    real( kind = core_rknd ) :: xtmp

    ! ---- Begin Code ----

    ! Sample fields

    if ( stats_metadata%l_stats_samp ) then

      ! stats_zt variables


      if ( stats_metadata%iT_in_K > 0 .or. stats_metadata%irsati > 0 ) then
        T_in_K = thlm2T_in_K( nzt, thlm, exner, rcm )
      else
        T_in_K = -999._core_rknd
      end if

      call stat_update_var( stats_metadata%iT_in_K, T_in_K, & ! intent(in)
                            stats_zt ) ! intent(inout)
 
      call stat_update_var( stats_metadata%ithlm, thlm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithvm, thvm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irtm, rtm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ircm, rcm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ium, um, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivm, vm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwm_zt, wm_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwm_zm, wm_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( stats_metadata%iug, ug, & ! intent(in) 
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivg, vg, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icloud_frac, cloud_frac, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iice_supersat_frac, ice_supersat_frac, & ! intent(in)
                            stats_zt) ! intent(inout)
      call stat_update_var( stats_metadata%ircm_in_layer, rcm_in_layer, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icloud_cover, cloud_cover, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ircm_supersat_adj, rcm_supersat_adj, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ip_in_Pa, p_in_Pa, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iexner, exner, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irho_ds_zt, rho_ds_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithv_ds_zt, thv_ds_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iLscale, Lscale, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwpup2, wpup2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwpvp2, wpvp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp3, wp3, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwpthlp2, wpthlp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp2thlp, wp2thlp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwprtp2, wprtp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp2rtp, wp2rtp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iLscale_up, Lscale_up, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iLscale_down, Lscale_down, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%itau_zt, tau_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iKh_zt, Kh_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp2thvp, wp2thvp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp2rcp, wp2rcp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iw_up_in_cloud, w_up_in_cloud, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iw_down_in_cloud, w_down_in_cloud, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icld_updr_frac, cloudy_updraft_frac, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icld_downdr_frac, cloudy_downdraft_frac, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwprtpthlp, wprtpthlp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irc_coef, rc_coef, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%isigma_sqd_w_zt, sigma_sqd_w_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irho, rho, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irsat, rsat, & ! intent(in)
                            stats_zt ) ! intent(inout)
      if ( stats_metadata%irsati > 0 ) then
        do k = 1, nzt
          rsati(k) = sat_mixrat_ice( p_in_Pa(k), T_in_K(k), saturation_formula )
        end do
        call stat_update_var( stats_metadata%irsati, rsati, & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      call stat_update_var( stats_metadata%imixt_frac, pdf_params%mixt_frac(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iw_1, pdf_params%w_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iw_2, pdf_params%w_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_w_1, pdf_params%varnce_w_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_w_2, pdf_params%varnce_w_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithl_1, pdf_params%thl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithl_2, pdf_params%thl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_thl_1, pdf_params%varnce_thl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_thl_2, pdf_params%varnce_thl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irt_1, pdf_params%rt_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irt_2, pdf_params%rt_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_rt_1, pdf_params%varnce_rt_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivarnce_rt_2, pdf_params%varnce_rt_2(1,:), & ! In
                            stats_zt ) ! intent(inout )
      call stat_update_var( stats_metadata%irc_1, pdf_params%rc_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irc_2, pdf_params%rc_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irsatl_1, pdf_params%rsatl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irsatl_2, pdf_params%rsatl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icloud_frac_1, pdf_params%cloud_frac_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icloud_frac_2, pdf_params%cloud_frac_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ichi_1, pdf_params%chi_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ichi_2, pdf_params%chi_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%istdev_chi_1, pdf_params%stdev_chi_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%istdev_chi_2, pdf_params%stdev_chi_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%istdev_eta_1, pdf_params%stdev_eta_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%istdev_eta_2, pdf_params%stdev_eta_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icovar_chi_eta_1, pdf_params%covar_chi_eta_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icovar_chi_eta_2, pdf_params%covar_chi_eta_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_chi_1, pdf_params%corr_w_chi_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_chi_2, pdf_params%corr_w_chi_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_eta_1, pdf_params%corr_w_eta_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_eta_2, pdf_params%corr_w_eta_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_chi_eta_1, pdf_params%corr_chi_eta_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_chi_eta_2, pdf_params%corr_chi_eta_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_rt_1, pdf_params%corr_w_rt_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_rt_2, pdf_params%corr_w_rt_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_thl_1, pdf_params%corr_w_thl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_w_thl_2, pdf_params%corr_w_thl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_rt_thl_1, pdf_params%corr_rt_thl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icorr_rt_thl_2, pdf_params%corr_rt_thl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icrt_1, pdf_params%crt_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icrt_2, pdf_params%crt_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icthl_1, pdf_params%cthl_1(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%icthl_2, pdf_params%cthl_2(1,:), & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp2_zt, wp2_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithlp2_zt, thlp2_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ithlp3, thlp3, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwpthlp_zt, wpthlp_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwprtp_zt, wprtp_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irtp2_zt, rtp2_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irtp3, rtp3, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%irtpthlp_zt, rtpthlp_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iup2_zt, up2_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivp2_zt, vp2_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iupwp_zt, upwp_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ivpwp_zt, vpwp_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%ia3_coef_zt, a3_coef_zt, & ! In
                            stats_zt ) ! intent(inout)
      call stat_update_var( stats_metadata%iwp3_on_wp2_zt, wp3_on_wp2_zt, & ! In
                            stats_zt ) ! In/Out

      if ( stats_metadata%ichi > 0 ) then
        ! Determine 's' from Mellor (1977) (extended liquid water)
        chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
                    + (1.0_core_rknd-pdf_params%mixt_frac(1,:)) * pdf_params%chi_2(1,:)
        call stat_update_var( stats_metadata%ichi, chi, & ! In
                             stats_zt ) ! In/Out
      end if 

      ! Calculate variance of chi
      if ( stats_metadata%ichip2 > 0 ) then
        chip2 = compute_variance_binormal( chi, pdf_params%chi_1(1,:), pdf_params%chi_2(1,:), &
                                         pdf_params%stdev_chi_1(1,:), pdf_params%stdev_chi_2(1,:), &
                                         pdf_params%mixt_frac(1,:) )
        call stat_update_var( stats_metadata%ichip2, chip2, & ! In
                              stats_zt ) ! In/Out
      end if

      if ( sclr_dim > 0 ) then
        do sclr=1, sclr_dim
          call stat_update_var( stats_metadata%isclrm(sclr), sclrm(:,sclr), & ! In
                                stats_zt ) ! In/Out
          call stat_update_var( stats_metadata%isclrm_f(sclr), sclrm_forcing(:,sclr),  & ! In
                                stats_zt ) ! In/Out
        end do
      end if

      if ( edsclr_dim > 0 ) then
        do edsclr = 1, edsclr_dim
          call stat_update_var( stats_metadata%iedsclrm(edsclr), edsclrm(:,edsclr), & ! In
                                stats_zt ) ! In/Out
          call stat_update_var( stats_metadata%iedsclrm_f(edsclr), edsclrm_forcing(:,edsclr), & ! In
                                stats_zt ) ! In/Out
        end do
      end if

      ! Calculate rcm in cloud
      if ( stats_metadata%ircm_in_cloud > 0 ) then
        where ( cloud_frac(:) > cloud_frac_min )
            rcm_in_cloud(:) = rcm / cloud_frac
        elsewhere
            rcm_in_cloud(:) = rcm
        endwhere

        call stat_update_var( stats_metadata%ircm_in_cloud, rcm_in_cloud, & ! intent(in)
                              stats_zt ) ! In/Out
      end if

      ! stats_zm variables

      call stat_update_var( stats_metadata%iwp2, wp2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp3_zm, wp3_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irtp2, rtp2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ithlp2, thlp2, & ! In
                            stats_zm ) ! In/Out 
      call stat_update_var( stats_metadata%irtpthlp, rtpthlp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwprtp, wprtp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwpthlp, wpthlp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp2up2, wp2up2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp2vp2, wp2vp2, &  ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp4, wp4, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwpthvp, wpthvp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irtpthvp, rtpthvp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ithlpthvp, thlpthvp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%itau_zm, tau_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iKh_zm, Kh_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwprcp, wprcp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irc_coef_zm, rc_coef_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ithlprcp, thlprcp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irtprcp, rtprcp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ircp2, rcp2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iupwp, upwp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ivpwp, vpwp, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ivp2, vp2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iup2, up2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irho_zm, rho_zm, & ! In
                            stats_zm ) ! In/Out 
      call stat_update_var( stats_metadata%isigma_sqd_w, sigma_sqd_w, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irho_ds_zm, rho_ds_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ithv_ds_zm, thv_ds_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iem, em, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iSkw_velocity, Skw_velocity, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ia3_coef, a3_coef, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp3_on_wp2, wp3_on_wp2, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iwp3_on_wp2_cfl_num, wp3_on_wp2 * dt / dzm, & ! In
                            stats_zm ) ! In/Out

      call stat_update_var( stats_metadata%icloud_frac_zm, cloud_frac_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iice_supersat_frac_zm, ice_supersat_frac_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ircm_zm, rcm_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%irtm_zm, rtm_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ithlm_zm, thlm_zm, & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iw_1_zm, pdf_params_zm%w_1(1,:), & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%iw_2_zm, pdf_params_zm%w_2(1,:), & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ivarnce_w_1_zm, pdf_params_zm%varnce_w_1(1,:), & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%ivarnce_w_2_zm, pdf_params_zm%varnce_w_2(1,:), & ! In
                            stats_zm ) ! In/Out
      call stat_update_var( stats_metadata%imixt_frac_zm, pdf_params_zm%mixt_frac(1,:), & ! In
                            stats_zm ) ! In/Out

      if ( sclr_dim > 0 ) then
        do sclr=1, sclr_dim
          call stat_update_var( stats_metadata%isclrp2(sclr), sclrp2(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%isclrprtp(sclr), sclrprtp(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%isclrpthvp(sclr), sclrpthvp(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%isclrpthlp(sclr), sclrpthlp(:,sclr), & ! In
                                 stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%isclrprcp(sclr), sclrprcp(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%iwpsclrp(sclr), wpsclrp(:,sclr), & ! In
                               stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%iwp2sclrp(sclr), wp2sclrp(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%iwpsclrp2(sclr), wpsclrp2(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%iwpsclrprtp(sclr), wpsclrprtp(:,sclr), & ! In
                                stats_zm ) ! In/Out
          call stat_update_var( stats_metadata%iwpsclrpthlp(sclr), wpsclrpthlp(:,sclr), & ! In
                                stats_zm ) ! In/Out
        end do
      end if
      if ( edsclr_dim > 0 ) then
        do edsclr = 1, edsclr_dim
          call stat_update_var( stats_metadata%iwpedsclrp(edsclr), wpedsclrp(:,edsclr), & ! In
                                stats_zm ) ! In/Out
        end do
      end if

      ! Calculate shear production
      if ( stats_metadata%ishear > 0 ) then
        do k = 2, nzm-1, 1
          shear(k) = - upwp(k) * ( um(k) - um(k-1) ) * invrs_dzm(k)  &
                     - vpwp(k) * ( vm(k) - vm(k-1) ) * invrs_dzm(k)
        enddo
        shear(1)   = 0.0_core_rknd
        shear(nzm) = 0.0_core_rknd
      end if
      call stat_update_var( stats_metadata%ishear, shear, & ! intent(in)
                            stats_zm ) ! intent(inout)

      ! stats_sfc variables

      ! Cloud cover
      call stat_update_var_pt( stats_metadata%icc, grid_level, & ! intent(in)
                               maxval( cloud_frac(1:nzt) ),    & ! intent(in)
                               stats_sfc )                       ! intent(inout)

      ! Cloud base
      if ( stats_metadata%iz_cloud_base > 0 ) then

        k = 1
        do while ( rcm(k) < rc_tol .and. k < nzt )
          k = k + 1
        enddo

        if ( k == 1 ) then

          ! Set the cloud base to the height of the lowest thermodynamic 
          ! grid level above the surface.
          call stat_update_var_pt( stats_metadata%iz_cloud_base, grid_level, & ! intent(in)
                                   zt(1),                                    & ! intent(in)
                                   stats_sfc )                                 ! intent(inout)

        elseif ( k > 1 .and. k < nzt ) then

          ! Use linear interpolation to find the exact height of the
          ! rc_tol kg/kg level.
          call stat_update_var_pt( stats_metadata%iz_cloud_base, grid_level, & ! intent(in)
                                   lin_interpolate_two_points( rc_tol, rcm(k), & ! intent(in)
                                   rcm(k-1), zt(k), zt(k-1) ), & ! intent(in)
                                   stats_sfc ) ! intent(inout)

        else

          ! Set the cloud base output to -10m, if it's clear. 
          ! Known magic number
          call stat_update_var_pt( stats_metadata%iz_cloud_base, grid_level, & ! intent(in)
                                   -10.0_core_rknd , & ! intent(in)
                                   stats_sfc ) ! intent(inout)
 
        end if

      end if ! stats_metadata%iz_cloud_base > 0

      ! Liquid Water Path
      if ( stats_metadata%ilwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( nzt, rho_ds_zt(1:nzt), &
                 rcm(1:nzt), dzt(1:nzt) )

        call stat_update_var_pt( stats_metadata%ilwp, grid_level, & ! intent(in)
                                 xtmp,                            & ! intent(in)
                                 stats_sfc )                        ! intent(inout)

      end if

      ! Vapor Water Path (Precipitable Water)
      if ( stats_metadata%ivwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( nzt, rho_ds_zt(1:nzt), &
                 ( rtm(1:nzt) - rcm(1:nzt) ), dzt(1:nzt) )

        call stat_update_var_pt( stats_metadata%ivwp, grid_level, & ! intent(in)
                                 xtmp,                            & ! intent(in)
                                 stats_sfc )                        ! intent(inout)

      end if


      ! Vertical average of thermodynamic level variables.

      ! Find the vertical average of thermodynamic level variables, averaged
      ! from level 1 through level nzt (the top of the model).  Use the vertical
      ! averaging function found in advance_helper_module.F90.

      ! Vertical average of thlm.
      call stat_update_var_pt( stats_metadata%ithlm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzt, rho_ds_zt(1:nzt), & ! intent(in)
                         thlm(1:nzt), dzt(1:nzt) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of rtm.
      call stat_update_var_pt( stats_metadata%irtm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzt, rho_ds_zt(1:nzt), & ! intent(in)
                         rtm(1:nzt), dzt(1:nzt) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of um.
      call stat_update_var_pt( stats_metadata%ium_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzt, rho_ds_zt(1:nzt), & ! intent(in)
                         um(1:nzt), dzt(1:nzt) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of vm.
      call stat_update_var_pt( stats_metadata%ivm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzt, rho_ds_zt(1:nzt), & ! intent(in)
                         vm(1:nzt), dzt(1:nzt) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of momentum level variables.

      ! Find the vertical average of momentum level variables, averaged over the
      ! entire vertical profile (level 1 through level nz).  Use the vertical
      ! averaging function found in advance_helper_module.F90.

      ! Vertical average of wp2.
      call stat_update_var_pt( stats_metadata%iwp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzm, rho_ds_zm(1:nzm), & ! intent(in)
                         wp2(1:nzm), dzm(1:nzm) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of up2.
      call stat_update_var_pt( stats_metadata%iup2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzm, rho_ds_zm(1:nzm), & ! intent(in)
                         up2(1:nzm), dzm(1:nzm) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of vp2.
      call stat_update_var_pt( stats_metadata%ivp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzm, rho_ds_zm(1:nzm), & ! intent(in)
                         vp2(1:nzm), dzm(1:nzm) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of rtp2.
      call stat_update_var_pt( stats_metadata%irtp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzm, rho_ds_zm(1:nzm), & ! intent(in)
                         rtp2(1:nzm), dzm(1:nzm) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of thlp2.
      call stat_update_var_pt( stats_metadata%ithlp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( nzm, rho_ds_zm(1:nzm), & ! intent(in)
                         thlp2(1:nzm), dzm(1:nzm) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)
      
      
      if (stats_metadata%itot_vartn_normlzd_rtm > 0) then
        if (abs(rtm(nzt) - rtm(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_rtm tried to divide by zero denominator ", &
                            "(surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(rtm(2 : nzt) - rtm(1 : nzt-1)) / abs(rtm(nzt) - rtm(1)))
        end if
        
        call stat_update_var_pt( stats_metadata%itot_vartn_normlzd_rtm, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
     
      if (stats_metadata%itot_vartn_normlzd_thlm > 0) then
        if (abs(thlm(nzt) - thlm(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_thlm tried to divide by zero ", &
                            "denominator (surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(thlm(2 : nzt) - thlm(1 : nzt-1)) / abs(thlm(nzt) - thlm(1)))
        end if
        
        call stat_update_var_pt( stats_metadata%itot_vartn_normlzd_thlm, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
     
      if (stats_metadata%itot_vartn_normlzd_wprtp > 0) then
        if (abs(wprtp(nzm) - wprtp(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_wprtp tried to divide by zero ", &
                            "denominator (surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(wprtp(2 : nzm) - wprtp(1 : nzm-1)) / abs(wprtp(nzm) - wprtp(1)))
        end if
        
        call stat_update_var_pt( stats_metadata%itot_vartn_normlzd_wprtp, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
    end if ! stats_metadata%l_stats_samp


    return
  end subroutine stats_accumulate
!------------------------------------------------------------------------------
  subroutine stats_accumulate_hydromet( gr, hydromet_dim, hm_metadata, & ! intent(in)
                                        hydromet, rho_ds_zt,           & ! intent(in)
                                        stats_metadata,                & ! intent(in)
                                        stats_zt, stats_sfc )            ! intent(inout)
! Description:
!   Compute stats related the hydrometeors

! References:
!   None
!------------------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type
        
    use corr_varnce_module, only: &
        hm_metadata_type

    use advance_helper_module, only: &
        vertical_integral ! Procedure(s)

    use stats_type_utilities, only: & 
        stat_update_var, & ! Procedure(s)
        stat_update_var_pt

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: & 
        stats_metadata_type

    implicit none

    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type (stats), intent(inout) :: &
      stats_zt, &
      stats_sfc

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    ! Local Variables
    real(kind=core_rknd) :: xtmp
    
    integer :: grid_level = 1

    ! ---- Begin Code ----

    if ( stats_metadata%l_stats_samp ) then

      if ( hm_metadata%iirr > 0 ) then
        call stat_update_var( stats_metadata%irrm, hydromet(:,hm_metadata%iirr), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iirs > 0 ) then
        call stat_update_var( stats_metadata%irsm, hydromet(:,hm_metadata%iirs), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if 

      if ( hm_metadata%iiri > 0 ) then 
        call stat_update_var( stats_metadata%irim, hydromet(:,hm_metadata%iiri), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iirg > 0 ) then
        call stat_update_var( stats_metadata%irgm,  &  ! intent(in)
                              hydromet(:,hm_metadata%iirg), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iiNi > 0 ) then
        call stat_update_var( stats_metadata%iNim, hydromet(:,hm_metadata%iiNi), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iiNr > 0 ) then
        call stat_update_var( stats_metadata%iNrm, hydromet(:,hm_metadata%iiNr), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iiNs > 0 ) then
        call stat_update_var( stats_metadata%iNsm, hydromet(:,hm_metadata%iiNs), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( hm_metadata%iiNg > 0 ) then
        call stat_update_var( stats_metadata%iNgm, hydromet(:,hm_metadata%iiNg), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      ! Snow Water Path
      if ( stats_metadata%iswp > 0 .and. hm_metadata%iirs > 0 ) then

        ! Calculate snow water path
        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iirs), gr%dzt(1,:) )

        call stat_update_var_pt( stats_metadata%iswp, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if ! stats_metadata%iswp > 0 .and. iirs > 0

      ! Ice Water Path
      if ( stats_metadata%iiwp > 0 .and. hm_metadata%iiri > 0 ) then

        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iiri), gr%dzt(1,:) )

        call stat_update_var_pt( stats_metadata%iiwp, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if

      ! Rain Water Path
      if ( stats_metadata%irwp > 0 .and. hm_metadata%iirr > 0 ) then

        xtmp &
        = vertical_integral &
               ( gr%nzt, rho_ds_zt, &
                 hydromet(:,hm_metadata%iirr), gr%dzt(1,:) )

        call stat_update_var_pt( stats_metadata%irwp, grid_level, & ! intent(in)
                                 xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
 
      end if ! stats_metadata%irwp > 0 .and. stats_metadata%irrm > 0
    end if ! stats_metadata%l_stats_samp

    return
  end subroutine stats_accumulate_hydromet
!------------------------------------------------------------------------------
  subroutine stats_accumulate_lh_tend( gr, hydromet_dim, hm_metadata, &
                                       lh_hydromet_mc, lh_Ncm_mc, &
                                       lh_thlm_mc, lh_rvm_mc, lh_rcm_mc, &
                                       lh_AKm, AKm, AKstd, AKstd_cld, &
                                       lh_rcm_avg, AKm_rcm, AKm_rcc, &
                                       stats_metadata, &
                                       stats_lh_zt )

! Description:
!   Compute stats for the tendency of latin hypercube sample points.

! References:
!   None
!------------------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use stats_type_utilities, only: & 
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        stats_metadata_type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_varnce_module, only: &
        hm_metadata_type

    use stats_type, only: &
        stats ! Type

    implicit none

    !----------------------- Input Variables -----------------------
    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      lh_hydromet_mc ! Tendency of hydrometeors except for rvm/rcm  [units vary]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      lh_Ncm_mc,  & ! Tendency of cloud droplet concentration  [num/kg/s]
      lh_thlm_mc, & ! Tendency of liquid potential temperature [kg/kg/s]
      lh_rcm_mc,  & ! Tendency of cloud water                  [kg/kg/s]
      lh_rvm_mc     ! Tendency of vapor                        [kg/kg/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      lh_AKm,     & ! Kessler ac estimate                 [kg/kg/s]
      AKm,        & ! Exact Kessler ac                    [kg/kg/s]
      AKstd,      & ! St dev of exact Kessler ac          [kg/kg/s]
      AKstd_cld,  & ! Stdev of exact w/in cloud ac        [kg/kg/s]
      lh_rcm_avg, & ! Monte Carlo rcm estimate            [kg/kg]
      AKm_rcm,    & ! Kessler ac based on rcm             [kg/kg/s]
      AKm_rcc       ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]
      
    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !----------------------- InOut Variables -----------------------
    type (stats), intent(inout) :: &
      stats_lh_zt

    !----------------------- Begin Code -----------------------

    if ( stats_metadata%l_stats_samp ) then

      call stat_update_var( stats_metadata%ilh_thlm_mc, lh_thlm_mc, & ! In
                            stats_lh_zt ) ! In/Out
      call stat_update_var( stats_metadata%ilh_rcm_mc, lh_rcm_mc, & ! In
                            stats_lh_zt ) ! In/Out
      call stat_update_var( stats_metadata%ilh_rvm_mc, lh_rvm_mc, & ! In
                            stats_lh_zt ) ! In/Out

      call stat_update_var( stats_metadata%ilh_Ncm_mc, lh_Ncm_mc, & ! In
                            stats_lh_zt ) ! In/Out

      if ( hm_metadata%iirr > 0 ) then
        call stat_update_var( stats_metadata%ilh_rrm_mc, lh_hydromet_mc(:,hm_metadata%iirr), & ! In
                              stats_lh_zt ) ! In/Out
      end if

      if ( hm_metadata%iirs > 0 ) then
        call stat_update_var( stats_metadata%ilh_rsm_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iirs),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if 

      if ( hm_metadata%iiri > 0 ) then
        call stat_update_var( stats_metadata%ilh_rim_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iiri),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if

      if ( hm_metadata%iirg > 0 ) then
        call stat_update_var( stats_metadata%ilh_rgm_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iirg),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if

      if ( hm_metadata%iiNi > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nim_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iiNi),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if

      if ( hm_metadata%iiNr > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nrm_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iiNr),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if

      if ( hm_metadata%iiNs > 0 ) then
        call stat_update_var( stats_metadata%ilh_Nsm_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iiNs),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if

      if ( hm_metadata%iiNg > 0 ) then
        call stat_update_var( stats_metadata%ilh_Ngm_mc,              & ! In
                              lh_hydromet_mc(:,hm_metadata%iiNg),  & ! In
                              stats_lh_zt )                             ! In/Out
      end if 

      call stat_update_var( stats_metadata%iAKm, AKm, & ! In
                            stats_lh_zt )               ! In/Out

      call stat_update_var( stats_metadata%ilh_AKm, lh_AKm, & ! In
                            stats_lh_zt)                      ! In/Out

      call stat_update_var( stats_metadata%ilh_rcm_avg, lh_rcm_avg, & ! In
                            stats_lh_zt )                             ! In/Out

      call stat_update_var( stats_metadata%iAKstd, AKstd, & ! In
                            stats_lh_zt )                   ! In/Out

      call stat_update_var( stats_metadata%iAKstd_cld, AKstd_cld, & ! In
                            stats_lh_zt )                       ! In/Out

      call stat_update_var( stats_metadata%iAKm_rcm, AKm_rcm, & ! In
                            stats_lh_zt)                        ! In/Out

      call stat_update_var( stats_metadata%iAKm_rcc, AKm_rcc, & ! In
                            stats_lh_zt )                       ! In/Out

    end if ! stats_metadata%l_stats_samp

    return
  end subroutine stats_accumulate_lh_tend
    
  !-----------------------------------------------------------------------
  subroutine stats_finalize( ngrdcol, stats_metadata, &
                             stats_zt, stats_zm, stats_sfc, &
                             stats_lh_zt, stats_lh_sfc, &
                             stats_rad_zt, stats_rad_zm )

    !     Description:
    !     Close NetCDF files and deallocate scratch space and
    !     stats file structures.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        stats_metadata_type

#ifdef NETCDF
    use output_netcdf, only:  & 
        close_netcdf ! Procedure
#endif

    use stats_type, only: stats ! Type

    implicit none

    integer, intent(in) :: &
      ngrdcol

    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    integer :: i

    if ( stats_metadata%l_stats .and. stats_metadata%l_netcdf ) then
#ifdef NETCDF

      ! We only created the files for the first column (see stats_init)
      ! so we only close for the first column as well

      call close_netcdf( stats_zt(1)%file ) ! intent(inout)
      call close_netcdf( stats_sfc(1)%file ) ! intent(inout)
      call close_netcdf( stats_zm(1)%file ) ! intent(inout)

      if ( stats_metadata%l_silhs_out ) then
        call close_netcdf( stats_lh_zt(1)%file ) ! intent(inout)
        call close_netcdf( stats_lh_sfc(1)%file ) ! intent(inout)
      end if

      if ( stats_metadata%l_output_rad_files ) then
        call close_netcdf( stats_rad_zt(1)%file ) ! intent(inout)
        call close_netcdf( stats_rad_zm(1)%file ) ! intent(inout)
      end if
#else
      error stop "This program was not compiled with netCDF support"
#endif
    end if

      
    if ( stats_metadata%l_stats ) then

      do i = 1, ngrdcol
        
        ! De-allocate all stats_zt variables
        if (allocated(stats_zt(i)%z)) then
          deallocate( stats_zt(i)%z )

          deallocate( stats_zt(i)%accum_field_values )
          deallocate( stats_zt(i)%accum_num_samples )
          deallocate( stats_zt(i)%l_in_update )

          deallocate( stats_zt(i)%file%grid_avg_var )
          deallocate( stats_zt(i)%file%z )
                
          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated( stats_zt(i)%file%lat_vals ) ) then
            deallocate( stats_zt(i)%file%lat_vals )
          end if
          if ( allocated( stats_zt(i)%file%lon_vals ) ) then
            deallocate( stats_zt(i)%file%lon_vals )
          end if

        end if

        if ( stats_metadata%l_silhs_out .and. allocated(stats_lh_zt(i)%z) ) then
          ! De-allocate all stats_lh_zt variables
          deallocate( stats_lh_zt(i)%z )

          deallocate( stats_lh_zt(i)%accum_field_values )
          deallocate( stats_lh_zt(i)%accum_num_samples )
          deallocate( stats_lh_zt(i)%l_in_update )

          deallocate( stats_lh_zt(i)%file%grid_avg_var )
          deallocate( stats_lh_zt(i)%file%z )
          
          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated(stats_lh_zt(i)%file%lat_vals ) ) then
            deallocate( stats_lh_zt(i)%file%lat_vals )
          end if
          if ( allocated(stats_lh_zt(i)%file%lon_vals ) ) then
            deallocate( stats_lh_zt(i)%file%lon_vals )
          end if

          ! De-allocate all stats_lh_sfc variables
          deallocate( stats_lh_sfc(i)%z )

          deallocate( stats_lh_sfc(i)%accum_field_values )
          deallocate( stats_lh_sfc(i)%accum_num_samples )
          deallocate( stats_lh_sfc(i)%l_in_update )

          deallocate( stats_lh_sfc(i)%file%grid_avg_var )
          deallocate( stats_lh_sfc(i)%file%z )
              
          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated( stats_lh_sfc(i)%file%lat_vals ) ) then
            deallocate( stats_lh_sfc(i)%file%lat_vals )
          end if
          if ( allocated( stats_lh_sfc(i)%file%lon_vals ) ) then
            deallocate( stats_lh_sfc(i)%file%lon_vals )
          end if
        end if ! l_silhs_out

        ! De-allocate all stats_zm variables
        if (allocated(stats_zm(i)%z)) then
          deallocate( stats_zm(i)%z )

          deallocate( stats_zm(i)%accum_field_values )
          deallocate( stats_zm(i)%accum_num_samples )
          deallocate( stats_zm(i)%l_in_update )

          deallocate( stats_zm(i)%file%grid_avg_var )
          deallocate( stats_zm(i)%file%z )
              
          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated( stats_zm(i)%file%lat_vals ) ) then
            deallocate( stats_zm(i)%file%lat_vals )
          end if
          if ( allocated( stats_zm(i)%file%lon_vals ) ) then
            deallocate( stats_zm(i)%file%lon_vals )
          end if
          
        end if

        if ( stats_metadata%l_output_rad_files ) then
          ! De-allocate all stats_rad_zt variables
          if (allocated(stats_rad_zt(i)%z)) then
            deallocate( stats_rad_zt(i)%z )

            deallocate( stats_rad_zt(i)%accum_field_values )
            deallocate( stats_rad_zt(i)%accum_num_samples )
            deallocate( stats_rad_zt(i)%l_in_update )

            deallocate( stats_rad_zt(i)%file%grid_avg_var )
            deallocate( stats_rad_zt(i)%file%z )
                
            ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
            if ( allocated( stats_rad_zt(i)%file%lat_vals ) ) then
              deallocate( stats_rad_zt(i)%file%lat_vals )
            end if
            if ( allocated( stats_rad_zt(i)%file%lon_vals ) ) then
              deallocate( stats_rad_zt(i)%file%lon_vals )
            end if

            ! De-allocate all stats_rad_zm variables
            deallocate( stats_rad_zm(i)%z )

            deallocate( stats_rad_zm(i)%accum_field_values )
            deallocate( stats_rad_zm(i)%accum_num_samples )
            deallocate( stats_rad_zm(i)%l_in_update )

            deallocate( stats_rad_zm(i)%file%grid_avg_var )
            deallocate( stats_rad_zm(i)%file%z )

            ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
            if ( allocated( stats_rad_zm(i)%file%lat_vals ) ) then
              deallocate( stats_rad_zm(i)%file%lat_vals )
            end if
            if ( allocated( stats_rad_zm(i)%file%lon_vals ) ) then
              deallocate( stats_rad_zm(i)%file%lon_vals )
            end if

          end if

        end if ! l_output_rad_files

        ! De-allocate all stats_sfc variables
        if (allocated(stats_sfc(i)%z)) then
          deallocate( stats_sfc(i)%z )

          deallocate( stats_sfc(i)%accum_field_values )
          deallocate( stats_sfc(i)%accum_num_samples )
          deallocate( stats_sfc(i)%l_in_update )

          deallocate( stats_sfc(i)%file%grid_avg_var )
          deallocate( stats_sfc(i)%file%z )
        end if
              
        ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
        if ( allocated( stats_sfc(i)%file%lat_vals ) ) then
          deallocate( stats_sfc(i)%file%lat_vals )
        end if
        if ( allocated( stats_sfc(i)%file%lon_vals ) ) then
          deallocate( stats_sfc(i)%file%lon_vals )
        end if

      end do

      ! De-allocate scalar indices
      if (allocated(stats_metadata%isclrm)) then
        deallocate( stats_metadata%isclrm )
        deallocate( stats_metadata%isclrm_f )
        deallocate( stats_metadata%iedsclrm )
        deallocate( stats_metadata%iedsclrm_f )
        deallocate( stats_metadata%isclrprtp )
        deallocate( stats_metadata%isclrp2 )
        deallocate( stats_metadata%isclrpthvp )
        deallocate( stats_metadata%isclrpthlp )
        deallocate( stats_metadata%isclrprcp )
        deallocate( stats_metadata%iwpsclrp )
        deallocate( stats_metadata%iwp2sclrp )
        deallocate( stats_metadata%iwpsclrp2 )
        deallocate( stats_metadata%iwpsclrprtp )
        deallocate( stats_metadata%iwpsclrpthlp )
        deallocate( stats_metadata%iwpedsclrp )
      end if

      ! De-allocate hyderometeor statistical variables
      if (allocated(stats_metadata%ihm_1)) then
        deallocate( stats_metadata%ihm_1 )
        deallocate( stats_metadata%ihm_2 )
        deallocate( stats_metadata%imu_hm_1 )
        deallocate( stats_metadata%imu_hm_2 )
        deallocate( stats_metadata%imu_hm_1_n )
        deallocate( stats_metadata%imu_hm_2_n )
        deallocate( stats_metadata%isigma_hm_1 )
        deallocate( stats_metadata%isigma_hm_2 )
        deallocate( stats_metadata%isigma_hm_1_n )
        deallocate( stats_metadata%isigma_hm_2_n )
        deallocate( stats_metadata%icorr_w_hm_1 )
        deallocate( stats_metadata%icorr_w_hm_2 )
        deallocate( stats_metadata%icorr_chi_hm_1 )
        deallocate( stats_metadata%icorr_chi_hm_2 )
        deallocate( stats_metadata%icorr_eta_hm_1 )
        deallocate( stats_metadata%icorr_eta_hm_2 )
        deallocate( stats_metadata%icorr_Ncn_hm_1 )
        deallocate( stats_metadata%icorr_Ncn_hm_2 )
        deallocate( stats_metadata%icorr_hmx_hmy_1 )
        deallocate( stats_metadata%icorr_hmx_hmy_2 )
        deallocate( stats_metadata%icorr_w_hm_1_n )
        deallocate( stats_metadata%icorr_w_hm_2_n )
        deallocate( stats_metadata%icorr_chi_hm_1_n )
        deallocate( stats_metadata%icorr_chi_hm_2_n )
        deallocate( stats_metadata%icorr_eta_hm_1_n )
        deallocate( stats_metadata%icorr_eta_hm_2_n )
        deallocate( stats_metadata%icorr_Ncn_hm_1_n )
        deallocate( stats_metadata%icorr_Ncn_hm_2_n )
        deallocate( stats_metadata%icorr_hmx_hmy_1_n )
        deallocate( stats_metadata%icorr_hmx_hmy_2_n )
        deallocate( stats_metadata%ihmp2_zt )
        deallocate( stats_metadata%iwp2hmp )
        deallocate( stats_metadata%ihydrometp2 )
        deallocate( stats_metadata%iwphydrometp )
        deallocate( stats_metadata%irtphmp )
        deallocate( stats_metadata%ithlphmp )
        deallocate( stats_metadata%ihmxphmyp )
        deallocate( stats_metadata%iK_hm )
      end if

      if ( allocated( stats_metadata%isilhs_variance_category ) ) then
        deallocate( stats_metadata%isilhs_variance_category )
        deallocate( stats_metadata%ilh_samp_frac_category )
      end if

    end if ! l_stats

    return
  end subroutine stats_finalize

!===============================================================================

!-----------------------------------------------------------------------
subroutine stats_check_num_samples( stats_grid, stats_metadata, &
                                    l_error )

! Description:
!   Ensures that each variable in a stats grid is sampled the correct
!   number of times.
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    use error_code, only: &
        clubb_at_least_debug_level   ! Procedure

    implicit none

  ! Input Variables
    type (stats), intent(in) :: &
      stats_grid               ! Grid type              [grid]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

  ! Input/Output Variables
    logical, intent(inout) :: &
      l_error                  ! Indicates an error     [boolean]

  ! Local Variables
    integer :: ivar, kvar      ! Loop variable          [index]

    logical :: l_proper_sample

!-----------------------------------------------------------------------

  !----- Begin Code -----

  ! Look for errors by checking the number of sampling points
  ! for each variable in the statistics grid at each vertical level.
  do ivar = 1, stats_grid%num_output_fields
    do kvar = 1, stats_grid%kk

      l_proper_sample = ( stats_grid%accum_num_samples(1,1,kvar,ivar) == 0 .or. &
                          stats_grid%accum_num_samples(1,1,kvar,ivar) &
                          == floor(  stats_metadata%stats_tout  &
                                   / stats_metadata%stats_tsamp ) )

      if ( .not. l_proper_sample ) then

        l_error = .true.  ! This will stop the run

        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) 'Possible sampling error for variable ',  &
                           trim(stats_grid%file%grid_avg_var(ivar)%name), ' in stats_grid ',  &
                           'at k = ', kvar,  &
                           '; stats_grid%accum_num_samples(',kvar,',',ivar,') = ', &
                            stats_grid%accum_num_samples(1,1,kvar,ivar)
        end if ! clubb_at_lest_debug_level 1


      end if ! .not. l_proper_sample

    end do ! kvar = 1 .. stats_grid%kk
  end do ! ivar = 1 .. stats_grid%num_output_fields

  return
end subroutine stats_check_num_samples
!-----------------------------------------------------------------------

end module stats_clubb_utilities
