!----------------------------------------------------------------------
! $Id$

program clubb_tuner 

!     Description:
!     ``Tunes'' constants in clubb so that the output matches LES output.
!     Uses amoeba or amebsa to calculate the min of (f_les - f_clubb)^2
! 
!     References:
!     _Numerical Recipes in Fortran 90_ (Chapter 10) 
!     (Amoeba & Amebsa subroutine)
!
!----------------------------------------------------------------------
  use error, only:  & 
    tuner_init, min_les_clubb_diff,                & ! Subroutines 
    write_results,                                 & ! Subroutine
    output_nml_standalone, output_nml_tuner,       & ! Subroutines
    param_vals_matrix, anneal_temp,                & ! Variables
    l_results_stdout, l_save_tuning_run,           & ! Variables
    l_results_file, tune_type, f_tol, ndim,        & ! Variables
    tuning_filename, file_unit                       ! Variable

  use error, only:  & 
    iamoeba, & ! Constants
    iamebsa, &
    iesa, &
    iflags, &
    iploops

  use constants_clubb, only: & 
    fstdout ! Variable

  use text_writer, only: &
    write_text ! Subroutine

  implicit none

  ! External
  external :: enhanced_simann_driver, amoeba_driver, amebsa_driver, &
    logical_flags_driver

  ! Variables
  character(len=10) :: current_time  ! Current time string (no seconds)
  character(len=8)  :: current_date  ! Current date string
  character(len=75) :: results_f     ! Results file

  character(len=3)  :: user_response ! Simple Y/N query

  !----------------- Begin Code -------------------

  ! Determine date and time to be used with file names
  call date_and_time( current_date, current_time )

  ! Create file to save tuning results if specified
  if ( l_save_tuning_run ) then
    ! File where tuning run results are written
    tuning_filename = "../input/tuning_run_results_"// &
      current_date//'_'//current_time(1:4)//".log"
    open(unit=file_unit, file=tuning_filename, action="write")
    write(file_unit,*) "Tuning..."
    close(unit=file_unit)
  end if ! l_save_tuning_run

  ! Read in namelists and define parameters
  call tuner_init( l_read_files=.true. )

  ! Attempt to find the optimal parameter set
  do
    select case ( tune_type )

        case ( iamoeba )
          call amoeba_driver( )

        case ( iamebsa )
          call amebsa_driver( )

        case ( iesa )
          call enhanced_simann_driver( )

        case ( iflags )
          call logical_flags_driver( current_date, current_time )
          stop "Program exited normally"

        case( iploops )
          call param_loops_driver( )
          write(fstdout,*) "All parameter sets have been run"

        case default
           stop "Unknown tuning type"
    end select

    ! Print to stdout if specified
    if ( l_results_stdout ) call write_results( fstdout )


    ! Save results in file if specified
    if ( l_save_tuning_run ) then
      open(unit=file_unit, file=tuning_filename, action="write", position="append")
      call write_results( file_unit )
    end if ! l_save_tuning_run


    ! Print completion message
    call write_text( "Run Complete.", l_save_tuning_run, file_unit )


    ! Query to see if we should exit the loop
    if ( tune_type /= iploops ) then
       ! Prompt the user to re-run with new parameters.
       write(unit=fstdout,fmt='(A)', advance='no')  & 
         "Re-run with new parameters?(y/n) "
       read(*,*) user_response
    else
       ! For the parameter loops tuner, automatically enter a response of "no".
       ! This allows the tuner to be run in the background.
       user_response = "n"
    endif


    ! Save tuning results in file if specified
    if( l_save_tuning_run ) then
      write(file_unit,*) "Re-run with new parameters?(y/n) ", user_response
      close(unit=file_unit)
    end if ! l_save_tuning_run

    
    ! what does a continue without a line number do? apparently what it says, continue :o
    select case ( trim( user_response ) )
        case ( "Yes", "yes", "y", "Y" )
          continue
        case default
          exit
    end select


    ! Save tuning results in file if specified
    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then

      if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
        action="write", position="append")

      call write_text( "Current f_tol= ", f_tol, l_save_tuning_run, file_unit )

      write(fstdout,fmt='(A)', advance='no') "Enter new f_tol=   "

      read(*,*) f_tol

      ! Save tuning results in file if specified
      if ( l_save_tuning_run ) write(file_unit,*) "Enter new f_tol=   ", f_tol

      if ( l_save_tuning_run ) close(unit=file_unit)

    end if


    if ( tune_type == iesa .or. tune_type == iamebsa ) then
      write(fstdout,fmt='(A)', advance='no') "New annealing temp =   "
      read(*,*) anneal_temp
    end if


    call tuner_init( l_read_files=.false. )


  end do ! user_response /= 'y', 'Y' or 'yes'

  ! Final namelist file output 

  if ( l_results_file ) then 

    ! Tuner namelist
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "Generating new error.in file...", l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    results_f = "../input/error_"//current_date//'_' & 
      //current_time(1:4)//".in" 

    ! Note:
    ! The first column of param_vals_matrix is the optimized result, which
    ! is swapped in by amoeba.
    call output_nml_tuner( results_f, param_vals_matrix(1,1:ndim) )
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "New filename is: "//results_f, l_save_tuning_run, file_unit )

    ! Parameters namelist
    call write_text( "Generating new tunable_parameters.in file...", l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    results_f = "../input/tunable_parameters/tunable_parameters_"//current_date//'_' & 
      //current_time(1:4)//".in" 

    call output_nml_standalone( results_f,  & 
                                param_vals_matrix(1,1:ndim) )
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "New filename is: "//results_f, l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

  end if

  ! Print message if tuning results were saved in a file
  if ( l_save_tuning_run ) then
    print*, "***The tuning results have been saved in the file: "
    print*, "  "//tuning_filename
  end if ! l_save_tuning_run

  ! Exit Program

  stop "Program exited normally"

  end program clubb_tuner

  !------------------------------------------------------------------------
  subroutine amoeba_driver

  !     Description:
  !     Simple interface for the amoeba minimization algorithm

  !     References:
  !     _Numerical Recipes in Fortran 90_.  See full citation above.
  !------------------------------------------------------------------------
#ifdef TUNER
  use nr, only:  & 
      amoeba ! Procedure(s)

  use error, only:  & 
    ! Variable(s)
    ndim,                                & ! Array dimensions
    param_vals_matrix, cost_fnc_vector,  & ! The 'p' matrix and 'y' vector resp.
    f_tol,                                & ! Tolerance of tuning run
    iter,                                & ! Iteration number
    min_les_clubb_diff,                  & ! Cost function
    min_err                                ! Minimum value of the cost function

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  real, dimension(1:ndim+1, 1:ndim) :: &
    param_vals_matrix_r4

  real, dimension(1:ndim+1) :: &
    cost_fnc_vector_r4

  ! ---- Begin Code ----

  param_vals_matrix_r4 = real(param_vals_matrix(1:ndim+1,1:ndim))
  cost_fnc_vector_r4 = real(cost_fnc_vector(1:ndim+1))

  call amoeba( param_vals_matrix_r4,  & 
               cost_fnc_vector_r4,  & 
               real(f_tol), min_les_clubb_diff, iter)

  param_vals_matrix(1:ndim+1,1:ndim) = real(param_vals_matrix_r4, kind = core_rknd)
  cost_fnc_vector(1:ndim+1) = real(cost_fnc_vector_r4, kind = core_rknd)


  ! Note:
  ! Amoeba will make the optimal cost result the first element of
  ! cost_fnc_vector and the optimal parameter set the first column of 
  ! param_vals_matrix, where param_vals_matrix is 'p' in the NR subroutine
  ! and cost_fnc_vector is 'y'.

  min_err = cost_fnc_vector(1)

  return
#else
  stop "Numerical recipes subroutines were disabled at compile time"
#endif
  end subroutine amoeba_driver

  !-----------------------------------------------------------------------

  subroutine amebsa_driver

  ! Description:
  !   Interface for the amoeba simulated annealing minimization algorithm
  !   At the end of the subroutine, the param_vals_matrix's first row gets the
  !   optimal values assigned to it.
  !
  ! References:
  !   None
  !-----------------------------------------------------------------------
#ifdef TUNER
  use nr, only:  & 
      amebsa ! Procedure(s)

  use nrtype, only:  & 
      SP ! Variable(s)

  use error, only: & 
      param_vals_matrix,  & ! Variable(s)
      anneal_temp, & 
      anneal_iter, & 
      iter, & 
      ndim, & 
      cost_fnc_vector, & 
      f_tol, & 
      min_err, & 
      min_les_clubb_diff ! Procedure(s)

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  ! Note:
  ! Most of these are taken from xamebsa, so I only have a vague idea
  ! what the their purpose is and can't pick less obscure names.
  ! dschanen 1 Apr 05
   
  ! Local Variables
  integer ::  & 
    iiter, jiter, & ! Loop variables
    nit ! ???

  real( kind = core_rknd ), dimension(ndim) ::  & 
    pb ! ???

  real( kind = core_rknd ) ::  & 
    ybb,   & ! ???
    yb,    & ! ???
    tmptr ! ???

  ! Temporary variables for passing into subroutines that require reals
  real, dimension(1:ndim+1, 1:ndim) :: &
    param_vals_matrix_r4

  real, dimension(1:ndim+1) :: &
    cost_fnc_vector_r4

  real, dimension(1:ndim) :: &
    pb_r4

  real :: &
    yb_r4

  ! ---- Begin Code ----

  ybb   = 1.0e30_core_rknd
  yb    = 1.0e30_core_rknd
  nit   = 0
  iiter = 0
  tmptr = anneal_temp ! anneal_temp taken from /stat/ namelist

  param_vals_matrix_r4 = real(param_vals_matrix(1:ndim+1,1:ndim))
  cost_fnc_vector_r4 = real(cost_fnc_vector(1:ndim+1))

  do jiter = 1, anneal_iter ! anneal_iter taken from /stat/ namelist


    iter  = iiter
    tmptr = tmptr * 0.8_core_rknd

    pb_r4 = real(pb(1:ndim))
    yb_r4 = real(yb)

    call amebsa & 
         ( param_vals_matrix_r4,  & 
           cost_fnc_vector_r4, & 
           pb_r4, yb_r4, real(f_tol), min_les_clubb_diff, iter, real(tmptr) )

    pb(1:ndim) = real(pb_r4, kind = core_rknd)
    yb = real(yb_r4, kind = core_rknd)

    nit = nit + iiter - iter

    if ( yb < ybb ) then
      ybb = yb
    end if

    if ( iter > 0 ) exit
  end do

  param_vals_matrix(1:ndim+1,1:ndim) = real(param_vals_matrix_r4, kind = core_rknd)
  cost_fnc_vector(1:ndim+1) = real(cost_fnc_vector_r4, kind = core_rknd)

  param_vals_matrix(1,1:ndim) = pb(1:ndim)
  min_err = ybb


  return

#else
  stop "Numerical recipes subroutines were disabled at compile time"
#endif
end subroutine amebsa_driver
!----------------------------------------------------------------------
subroutine enhanced_simann_driver

  ! Description:
  !   Wrapper subroutine for the ESA driver

  ! References: 
  !   ``Enhanced Simulated Annealing for Many Globally Minimized Functions
  !   of Many-Continuous Variables'', Siarry, et al. ACMS TOMS Vol. 23,
  !   No. 2, June 1997, pp. 209--228.
  !------------------------------------------------------------------------

  use enhanced_simann, only: &
        esa_driver, & ! Procedure(s)
        esa_driver_siarry

  use error, only: & ! Variable(s)
    ndim,                 & ! Array dimensions
    param_vals_matrix,    & ! The parameters to tune matrix
    param_vals_max,       & ! The maximum values for the parameters
    anneal_temp,          & ! Start annealing temperature
    max_final_temp,       & ! Maximum final annealing temperature
    min_err,              & ! Minimum value of the cost function
    stp_adjst_center_in,  & 
    stp_adjst_spread_in,  &
    iter,                 &
    tuning_filename,      &
    file_unit,            &
    f_tol

  use error, only:  & ! Procedure(s)
    min_les_clubb_diff  ! Cost function

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  logical, parameter :: l_esa_siarry = .false.    ! Toggle between old and new code

  ! Local Variables

  real( kind = core_rknd ), dimension(ndim) :: &
    xinit,  & ! Initial values for the tunable parameters
    xmin,  & ! Minimum value for the tunable parameters
    xmax,  & ! Maximum value for the tunable parameters
    xopt      ! Final values for the tunable parameters

  real( kind = core_rknd ) :: enopt ! Optimal cost

  ! ---- Begin Code ----

  xinit = param_vals_matrix(1,1:ndim)

  ! set constraints, this should be read in from file
  xmin = 0._core_rknd 
  xmax = param_vals_max

  if ( l_esa_siarry ) then 
    call esa_driver_siarry( xinit, xmin, xmax, anneal_temp, min_les_clubb_diff, xopt, enopt )
  else 
    call esa_driver( xinit, xmin, xmax,                         & ! intent(in)
                     anneal_temp, max_final_temp,               & ! intent(inout)
                     xopt, enopt,                               & ! intent(out)
                     min_les_clubb_diff,                        & ! procedure    
                     stp_adjst_center_in, stp_adjst_spread_in,  & ! optional(in)
                     f_tol, tuning_filename, file_unit,         & ! ^
                     iter                                       ) ! optional(inout) 
  end if

                   
  param_vals_matrix(1,1:ndim) = xopt(1:ndim)
  min_err = enopt

  return
end subroutine enhanced_simann_driver
!-------------------------------------------------------------------------------
subroutine logical_flags_driver( current_date, current_time )

! Description:
!   While not a true search algorithm in same sense as the simulated annealing
!   or the downhill simplex method is, this driver will try all permutations.
!
! References:
!   None
!-------------------------------------------------------------------------------
  use error, only: &
    param_vals_matrix, & ! Variable(s)
    model_flags_array, &
    iter, &
    l_results_file, &
    l_results_stdout

  use error, only: &
    min_les_clubb_diff ! Procedure(s)

  use constants_clubb, only: &
    fstdout ! Constant(s)

  use quicksort, only: &
    Qsort_flags ! Procedure(s)

  use model_flags, only: &
    set_default_clubb_config_flags, & ! Procedure(s)
    l_quintic_poly_interp ! Variable(s)

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! External
  intrinsic :: btest, selected_int_kind, int

  ! Constant parameters
  integer, parameter :: &
    i8 = selected_int_kind( 15 )

  integer, parameter :: &
    ndim = 10, & ! Temporarily hardwired for a fixed number of flags
    two_ndim = 2**ndim, &
    iunit = 10

  ! Input Variables
  character(len=10), intent(in) :: current_time  ! Current time string (no seconds)
  character(len=8), intent(in)  :: current_date  ! Current date string

  ! Local Variables

  real( kind = core_rknd ), dimension(two_ndim) :: &
    cost_function  ! Values from the cost function

  real( kind = core_rknd ), dimension(ndim) :: &
    cost_func_sum_true,  & ! Averaged cost function when the flag is true
    cost_func_sum_false, & ! Averaged cost function when the flag is false
    cost_func_avg          ! Averaged cost function true - false.

  logical, dimension(ndim) :: model_flags_default

  character(len=128) :: &
    filename_nml, &  ! Namelist file name
    filename_csv     ! Comma seperated values filename

  integer(kind=i8) :: bit_string, bit_iter
  real( kind = core_rknd) :: cost_func_default

  integer :: i, j

  integer :: &
    iiPDF_type,          & ! Selected option for the two-component normal
                           ! (double Gaussian) PDF type to use for the w, rt,
                           ! and theta-l (or w, chi, and eta) portion of
                           ! CLUBB's multivariate, two-component PDF.
    ipdf_call_placement    ! Selected option for the placement of the call to
                           ! CLUBB's PDF.

  logical :: &
    l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                    ! precipitation fraction is automatically set to 1 when this
                                    ! flag is turned off.
    l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                    ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                    ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                    ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                    ! approximated by eddy diffusivity when <u> and <v> are
                                    ! advanced in subroutine advance_windm_edsclrm.
    l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                    ! the overall correlation of w and x (w and rt, as well as w
                                    ! and theta-l) within the limits of -max_mag_correlation_flux
                                    ! to max_mag_correlation_flux.
    l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                    ! thlp2) on keeping the overall correlation of w and x within
                                    ! the limits of -max_mag_correlation_flux to
                                    ! max_mag_correlation_flux.
    l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                    ! turbulent dissipation coefficient, C2.
    l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
    l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
    l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
    l_upwind_wpxp_ta,             & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered
                                    ! differencing for turbulent or mean advection terms. It
                                    ! affects wprtp, wpthlp, & wpsclrp.
    l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered
                                    ! differencing for turbulent or mean advection terms. It
                                    ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                    ! sclrpthlp.
    l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered
                                    ! differencing for turbulent or mean advection terms. It
                                    ! affects rtm, thlm, sclrm, um and vm.
    l_uv_nudge,                   & ! For wind speed nudging.
    l_rtm_nudge,                  & ! For rtm nudging
    l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                    ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
    l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                    ! compute the varibles that are output from high order
                                    ! closure
    l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                    ! thermodynamic-level variables output from pdf_closure.
    l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                    ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                    ! output from pdf_closure.
    l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                    ! subroutine pdf_closure twice.  If true, pdf_closure is
                                    ! called first on thermodynamic levels and then on momentum
                                    ! levels so that each variable is computed on its native
                                    ! level.  If false, pdf_closure is only called on
                                    ! thermodynamic levels, and variables which belong on
                                    ! momentum levels are interpolated.
    l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                    ! terms.  Setting to .false. means that a_1 and a_3 are
                                    ! pulled outside of the derivative in
                                    ! advance_wp2_wp3_module.F90 and in
                                    ! advance_xp2_xpyp_module.F90.
    l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                    ! and rcm to help increase cloudiness at coarser grid
                                    ! resolutions.
    l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
    l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
    l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
    l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
    l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                    ! correction
    l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                    ! -(2/3)*em/tau_zm, as in Bougeault (1981)
    l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
    l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
    l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                    ! mixing length scale as Lscale = tau * tke
    l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
    l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
    l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
    l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
    l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
    l_single_C2_Skw,              & ! Use a single Skewness dependent C2 for rtp2, thlp2, and
                                    ! rtpthlp
    l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
    l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
    l_update_pressure               ! Flag for having CLUBB update pressure and exner

  namelist /configurable_clubb_flags_nl/ &
    iiPDF_type, ipdf_call_placement, &
    l_upwind_wpxp_ta, l_upwind_xpyp_ta, l_upwind_xm_ma, l_quintic_poly_interp, &
    l_tke_aniso, l_vert_avg_closure, l_single_C2_Skw, l_standard_term_ta, &
    l_use_cloud_cover, l_rcm_supersat_adj, l_damp_wp3_Skw_squared, &
    l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, l_C2_cloud_frac, &
    l_predict_upwp_vpwp, l_diag_Lscale_from_tau, l_stability_correct_tau_zm, &
    l_damp_wp2_using_em, l_use_C7_Richardson, l_use_precip_frac, l_do_expldiff_rtm_thlm, &
    l_use_C11_Richardson, l_prescribed_avg_deltaz, l_diffuse_rtm_and_thlm, &
    l_stability_correct_Kh_N2_zm, l_trapezoidal_rule_zt, l_trapezoidal_rule_zm, &
    l_call_pdf_closure_twice, l_Lscale_plume_centered, &
    l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, l_update_pressure

  ! ---- Begin Code ----

  call set_default_clubb_config_flags( iiPDF_type, & ! Intent(out)
                                       ipdf_call_placement, & ! Intent(out)
                                       l_use_precip_frac, & ! Intent(out)
                                       l_predict_upwp_vpwp, & ! Intent(out)
                                       l_min_wp2_from_corr_wx, & ! Intent(out)
                                       l_min_xp2_from_corr_wx, & ! Intent(out)
                                       l_C2_cloud_frac, & ! Intent(out)
                                       l_diffuse_rtm_and_thlm, & ! Intent(out)
                                       l_stability_correct_Kh_N2_zm, & ! Intent(out)
                                       l_calc_thlp2_rad, & ! Intent(out)
                                       l_upwind_wpxp_ta, & ! Intent(out)
                                       l_upwind_xpyp_ta, & ! Intent(out)
                                       l_upwind_xm_ma, & ! Intent(out)
                                       l_uv_nudge, & ! Intent(out)
                                       l_rtm_nudge, & ! Intent(out)
                                       l_tke_aniso, & ! Intent(out)
                                       l_vert_avg_closure, & ! Intent(out)
                                       l_trapezoidal_rule_zt, & ! Intent(out)
                                       l_trapezoidal_rule_zm, & ! Intent(out)
                                       l_call_pdf_closure_twice, & ! Intent(out)
                                       l_standard_term_ta, & ! Intent(out)
                                       l_use_cloud_cover, & ! Intent(out)
                                       l_diagnose_correlations, & ! Intent(out)
                                       l_calc_w_corr, & ! Intent(out)
                                       l_const_Nc_in_cloud, & ! Intent(out)
                                       l_fix_w_chi_eta_correlations, & ! Intent(out)
                                       l_stability_correct_tau_zm, & ! Intent(out)
                                       l_damp_wp2_using_em, & ! Intent(out)
                                       l_do_expldiff_rtm_thlm, & ! Intent(out)
                                       l_Lscale_plume_centered, & ! Intent(out)
                                       l_diag_Lscale_from_tau, & ! Intent(out)
                                       l_use_C7_Richardson, & ! Intent(out)
                                       l_use_C11_Richardson, & ! Intent(out)
                                       l_brunt_vaisala_freq_moist, & ! Intent(out)
                                       l_use_thvm_in_bv_freq, & ! Intent(out)
                                       l_rcm_supersat_adj, & ! Intent(out)
                                       l_single_C2_Skw, & ! Intent(out)
                                       l_damp_wp3_Skw_squared, & ! Intent(out)
                                       l_prescribed_avg_deltaz, & ! Intent(out)
                                       l_update_pressure ) ! Intent(out)

  ! Determine the current flags
  model_flags_default(1) = l_upwind_wpxp_ta
  model_flags_default(2) = l_upwind_xpyp_ta
  model_flags_default(3) = l_upwind_xm_ma
  model_flags_default(4) = l_quintic_poly_interp
  model_flags_default(5) = l_vert_avg_closure
  model_flags_default(6) = l_single_C2_Skw
  model_flags_default(7) = l_standard_term_ta
  model_flags_default(8) = l_tke_aniso
  model_flags_default(9) = l_use_cloud_cover
  model_flags_default(10) = l_rcm_supersat_adj

  ! This should always be 1.0; it's here as a sanity check
  cost_func_default = real( min_les_clubb_diff( real(param_vals_matrix(1,:)) ), kind = core_rknd )

  allocate( model_flags_array(two_ndim,ndim) )
  bit_string = 0_i8 ! Initialize bits to 00 ... 00
  do i = 1, two_ndim
    do j = 1, ndim
     ! This loop sets 1:n logicals using individual bits, i.e. 0 means
     ! false and 1 means true for the purposes of trying all possibilities
     bit_iter = int( j, i8 )
     model_flags_array(i,j) = btest( bit_string, bit_iter-1_i8 ) 
    end do
    bit_string = bit_string + 1_i8 ! Increment the binary adder
  end do

  ! Compute the cost function with new set of flags.  The model_flags array is
  ! indexed using the iter variable in min_les_clubb_diff to avoid having to
  ! modify the Numerical Recipes code.
  do iter = 1, two_ndim
    cost_function(iter) = real( min_les_clubb_diff( real(param_vals_matrix(1,:)) ), &
         kind = core_rknd )
  end do

  ! Compute a metric of false cost function - true cost function
  cost_func_sum_true = 0.0_core_rknd
  cost_func_sum_false = 0.0_core_rknd
  do i = 1, two_ndim
    do j = 1, ndim
      if ( model_flags_array(i,j) ) then ! Flag is true
        cost_func_sum_true(j) = cost_func_sum_true(j) + cost_function(i)
      else ! Flag is false
        cost_func_sum_false(j) = cost_func_sum_false(j) + cost_function(i)
      end if
    end do
  end do
  cost_func_avg(:) = ( cost_func_sum_false(:) / real( two_ndim / 2, kind = core_rknd ) ) &
                   - ( cost_func_sum_true(:) / real( two_ndim / 2, kind = core_rknd ) )

  ! Sort flags and the cost function in ascending order
  call Qsort_flags( model_flags_array, cost_function )

  if ( l_results_stdout ) then
    ! Output results to the terminal
    write(fstdout,'(A30)') "Default flags:                "
    do j = 1, ndim
      write(fstdout,'(L6,4X)',advance='no') model_flags_default(j)
    end do
    write(fstdout,'(G10.3)') cost_func_default
    write(fstdout,'(A)') "Results from trying all permutations of the flags: "
    do i = 1, ndim
      write(fstdout,'(I6,4X)',advance='no') i
    end do
    write(fstdout,'(A10)') " Cost func" 
    do i = 1, ndim
      write(fstdout,'(G10.3)',advance='no') cost_func_avg(i)
    end do
    write(fstdout,*) ""
    do i = 1, two_ndim
      do j = 1, ndim
        write(fstdout,'(L6,4X)',advance='no') model_flags_array(i,j)
      end do
      write(fstdout,'(G10.3)') cost_function(i)
    end do

    write(fstdout,'(A30)') &
      "Column 1 = upwind_wpxp_ta     ", &
      "Column 2 = upwind_xpyp_ta     ", &
      "Column 3 = upwind_xm_ma       ", &
      "Column 4 = quintic_poly_interp", &
      "Column 5 = vert_avg_closure   ", &
      "Column 6 = single_C2_Skw      ", &
      "Column 7 = standard_term_ta   ", &
      "Column 8 = tke_aniso          ", &
      "Column 9 = use_cloud_cover    "
  end if ! l_results_stdout

  ! Generate CSV file of the results
  filename_csv = "../output/clubb_model_flags_"//current_date//"_"//current_time//".csv"
  open(unit=iunit,file=trim( filename_csv ))
  write(iunit,'(10A20)') "upwind_wpxp_ta, ", "upwind_xpyp_ta, ", "upwind_xm_ma, ", &
    "quintic_poly_interp, ", "vert_avg_closure, ", &
    "single_C2_Skw, ", "standard_term_ta, ", "tke_aniso, ", "use_cloud_cover, ", "Cost func."
  write(iunit,'(A30)') "Default flags:               ,"
  do j = 1, ndim
    write(iunit,'(L20,A2)',advance='no') model_flags_default(j), ", "
  end do
  write(iunit,'(G20.6,A2)') cost_func_default, ", "
  do i = 1, ndim
    write(iunit,'(G20.6,A2)',advance='no') cost_func_avg(i), ", "
  end do
  write(iunit,'(A20)') "Avg false - Avg true ,"
  do i = 1, two_ndim
    do j = 1, ndim
      write(iunit,'(L20,A2)',advance='no') model_flags_array(i,j), ", "
    end do
    write(iunit,'(G20.6,A2)') cost_function(i), ", "
  end do
  close(unit=iunit)
  write(fstdout,*) "Results of tuning model flags written to: ", trim( filename_csv )

  ! Generate namelist file of the optimal result
  if ( l_results_file ) then
    l_upwind_wpxp_ta = model_flags_array(1,1)
    l_upwind_xpyp_ta = model_flags_array(1,2)
    l_upwind_xm_ma = model_flags_array(1,3)
    l_quintic_poly_interp = model_flags_array(1,4)
    l_vert_avg_closure = model_flags_array(1,5)
    l_single_C2_Skw = model_flags_array(1,6)
    l_standard_term_ta = model_flags_array(1,7)
    l_tke_aniso = model_flags_array(1,8)
    l_use_cloud_cover = model_flags_array(1,9)
    l_rcm_supersat_adj = model_flags_array(1,10)

    if ( l_vert_avg_closure ) then
      l_trapezoidal_rule_zt    = .true.
      l_trapezoidal_rule_zm    = .true.
      l_call_pdf_closure_twice = .true.
    else
      l_trapezoidal_rule_zt    = .false.
      l_trapezoidal_rule_zm    = .false.
      l_call_pdf_closure_twice = .false.
    end if

    filename_nml = "../input/tunable_parameters/configurable_model_flags_"//current_date//'_' & 
      //current_time(1:4)//".in"

    ! Write namelist to file
    open(unit=iunit, file=trim( filename_nml ), status='unknown', action='write')
    write(unit=iunit, nml=configurable_clubb_flags_nl)
    close(unit=iunit)

    write(fstdout,*) "New namelist of tuning model flags written to: ", trim( filename_nml )
  end if

  deallocate( model_flags_array )

  return
end subroutine logical_flags_driver
  !------------------------------------------------------------------------
  subroutine param_loops_driver( )

    ! Description:
    ! The parameters loops tuner allows the user to select (through hard-wired
    ! coding) a few tunable parameters (up to 4) and define values to loop over
    ! for each parameters.  This allows the tuner to loop through all of
    ! parameter space to look for the lowest value of the cost function.
    !
    ! This tuner does not use random numbers or "downhill" tuning, etc.  It
    ! loops over all of the defined values in parameter space, calculates the
    ! cost function for each parameter set, and exits.  It does not do further
    ! tuning.
    !
    ! The number of values and the values for each selected parameter are set
    ! in subroutine tuner_init in error.F90.

    !-----------------------------------------------------------------------

    use error, only: &
        ndim,              &
        param_vals_matrix, &
        cost_fnc_vector,   &
        min_err

    implicit none

    ! Local Variables
    integer :: &
      min_idx    ! Index for parameter set with minimum value of cost function


    min_idx = minloc( cost_fnc_vector, 1 )

    ! Place the parameter set with the minimum value of the cost function in
    ! the 1st row of the parameter array.
    param_vals_matrix(1,1:ndim) = param_vals_matrix(min_idx,1:ndim)

    min_err = cost_fnc_vector(min_idx)


    return

  end subroutine param_loops_driver
