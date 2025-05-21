# Reproduce results from 'Adaptive grid redistribution for a 1D model of turbulence and clouds'

To reproduce the results and plots from the paper, this branch needs to be checked out with the tag: `TODO`.

Then the code needs to be compiled using `./compile/compile.bash`.

## Reproducing the grid adaptation results
1. Run `./run_scripts/run_scm.bash arm`, `./run_scripts/run_scm.bash astex_a209` and `./run_scripts/run_scm.bash gabls2`

## Reproducing the dycore results
1. Change `grid_adapt_in_time_method` to `0` in `./input/tunable_parameters/configurable_model_flags.in`
2. Run `./run_scripts/run_scm.bash arm`, `./run_scripts/run_scm.bash astex_a209` and `./run_scripts/run_scm.bash gabls2`

## Reproducing the hi-res/benchmark results
1. Change `grid_adapt_in_time_method` to `0` in `./input/tunable_parameters/configurable_model_flags.in`
2. Set the following files accordingly  
`./input/case_setups/arm_model.in`
```
! $Id$
! Parameter file for GCSS ARM Cu over land case

&model_setting

! Model Settings

runtype = "arm" ! Case name.

! For grid_type = 1 nzmax can be any nunber
nzmax = 426       ! Number of vertical levels.
!nzmax = 73       ! Number of vertical levels.

grid_type = 1     ! Select grid type:
                  ! 1 ==> evenly-spaced grid levels.
                  ! 2 ==> stretched (unevenly-spaced) grid
                  !       entered on thermodynamic grid levels;
                  !       mom. levels halfway between thermo. levels
                  !       (style of SAM stretched grid).
                  ! 3 ==> stretched (unevenly-spaced) grid
                  !       entered on momentum grid levels;
                  !       thermo. levels halfway between mom. levels
                  !       (style of WRF stretched grid).

! Evenly-spaced grid (grid_type = 1).
! Note:  For this option, both zt_grid_fname and zm_grid_fname must be ''.
deltaz_nl  = 20.0    ! Distance between grid levels on evenly-spaced grid.      [m]
zm_init_nl = 0.0 ! Minimum Altitude of lowest momentum level on any grid. [m]
!zm_top_nl = 5350.0 ! Maximum Altitude of highest momentum level on any grid. [m]
zm_top_nl = 8500.0 ! Maximum Altitude of highest momentum level on any grid. [m]

! Stretched grid on thermodynamic levels (grid_type = 2).
! Path (from ../run_scripts/) and filename of data file that contains
! thermodynamic level altitudes (in meters).
! Note:  For this option, zm_grid_fname must be ''.
zt_grid_fname = ''

! Stretched grid on momentum levels (grid_type = 3).
! Path (from ../run_scripts/) and filename to data file that contains
! momentum level altitudes (in meters).
! Note:  For this option, zt_grid_fname must be ''.
zm_grid_fname = ''
!zm_grid_fname = '../input/grid/dycore.grd'

day   = 21    ! Day of model start (1 to 31).
month = 6     ! Month of model start (1 to 12).
year  = 1997  ! Year of model start.

lat_vals = 36.62  ! Latitude (0 to 90 for north lat.; 0 to -90 for south lat.) [degrees]
lon_vals = -97.5  ! Longitude (0 to 180 for east long.; 0 to -180 for west long.) [degrees]

sfc_elevation_nl = 0.0     ! Elevation of ground level [meters above mean sea level]

time_initial = 41400.0  ! Model start time [seconds since midnight UTC on start date]
time_final = 93600.0    ! Model end time [seconds since midnight UTC on start date]

dt_main    = 60.0   ! Model timestep [s]
dt_rad = 60.0   ! Radiation timestep [s]

l_t_dependent = .true. ! Read values from _surface.in file. Also read
                       ! _forcings.in file if l_ignore_forcings = .false.


sfctype = 1      ! Select surface scheme (cases without specified surface fluxes).
                 ! 0 ==> fixed surface sensible and latent heat fluxes,
                 !       as defined by sens_ht and latent_ht, respectively.
                 ! 1 ==> determine surface sensible and latent heat fluxes
                 !       through a bulk formula that uses the given surface
                 !       temperature (T_sfc) and assumes over ocean.

p_sfc_nl = 970.e2    ! Pressure at model base              [Pa]

fcor_nl = 8.5e-5    ! Coriolis parameter                  [s^-1]
T0   = 300.      ! Reference temperature (usually 300) [K]



! Sponge Damping
!
! tau_sponge_damp_min : Minimum damping time-scale ( at the top0 [s]
! tau_sponge_damp_max : Maximum damping time-scale (base of damping layer) [s]
! sponge_damp_depth   : damping depth as a fraction of domain height [-]
!
thlm_sponge_damp_settings%l_sponge_damping = .false.
thlm_sponge_damp_settings%tau_sponge_damp_min = 60. 
thlm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
thlm_sponge_damp_settings%sponge_damp_depth = 0.25

rtm_sponge_damp_settings%l_sponge_damping = .false.
rtm_sponge_damp_settings%tau_sponge_damp_min = 60. 
rtm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
rtm_sponge_damp_settings%sponge_damp_depth = 0.25

uv_sponge_damp_settings%l_sponge_damping = .false.
uv_sponge_damp_settings%tau_sponge_damp_min = 60. 
uv_sponge_damp_settings%tau_sponge_damp_max = 1800. 
uv_sponge_damp_settings%sponge_damp_depth = 0.25

l_restart      = .false.  ! Flag for whether this is a restart run.
! Path (from ../run_scripts/) and filename of GrADS control files
! (without '_zt.ctl' or '_zm.ctl').
restart_path_case  = "restart/arm"
time_restart  = 0.0      ! Time of model restart.
                         ! [seconds since start time listed in GrADS control file]

! Select debug level:
!        0 => Print no debug messages to the screen
!        1 => Print lightweight debug messages, e.g. print statements
!        2 => Print debug messages that require extra testing,
!                e.g. checks for NaNs and spurious negative values.
debug_level = 2

! Add scalars that are computed using higher-order closure
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_sclr_sounding.in 
!	that will be read in and tracked
sclr_dim   = 0  ! Total number of scalars

iisclr_rt  = -1  ! Location in the array to place a scalar like total water
iisclr_thl = -2 ! Location in the array to place a scalar like thetal

sclr_tol_nl  = 1.e-2, 1.e-8 ! Tolerances below which we consider scalar to be 0.

! Add scalars that are computed using simple eddy diffusivity
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_edsclr_sounding.in 
!	that will be read in and tracked
edsclr_dim  = 0

iiedsclr_rt     = -1 ! Location in the edsclrm array to place total water
iiedsclr_thl    = -2 ! Location in the edsclrm array to place a scalar like thetal

l_input_fields = .false.
/

&gfdl_activation_setting

/

&microphysics_setting
! The schemes are: "khairoutdinov_kogan", "morrison", 
! "coamps", "simplified_ice", "none"
microphys_scheme = "none" 
/

&radiation_setting
! The schemes are: "bugsrad", simplified", "simplified_bomex", or "none"
rad_scheme = "none",
l_calc_thlp2_rad = .false.,
/

&stats_setting

! Statistics Settings

l_stats      = .true.   ! Flag for statistical output.
fname_prefix = "arm"    ! Prefix of output filename.
stats_tsamp  = 60.      ! Frequency of statistical sampling [s]
                        ! For most complete sampling, let 
                        ! stats_tsamp = dt_main
                        ! (stats_tsamp must multiply evenly into stats_tout)
stats_tout   = 60.      ! Frequency of statistical output.  [s]
stats_fmt    = "netcdf"  ! Type of statistical output file ("grads" or "netcdf").

/

&setfields
datafile = "../sam_benchmark_runs/ARM_96x96x110/GCSSARM_96x96x110_67m_40m_1s.nc"
input_type = "SAM"  ! input_type = "clubb", "coamps_les", or "SAM"
l_input_um = .true. 
l_input_vm = .true.
l_input_rtm = .true. 
l_input_thlm = .true.
l_input_wp2 = .true.
l_input_wprtp = .true.
l_input_wpthlp = .true.
l_input_wp3 = .true.
l_input_rtp2 = .true.
l_input_thlp2 = .true.
l_input_rtpthlp = .true.
l_input_upwp = .true.
l_input_vpwp = .true.
l_input_rcm = .true.
l_input_em = .true.
l_input_p = .true.
l_input_up2 = .true.
l_input_vp2 = .true.
/
```

`./input/case_setups/astex_a209_model.in`
```
! $Id$
! Parameter file for ASTEX intermediate regime

&model_setting

! Model Settings

runtype = "astex_a209"  ! Case name.

! For grid_type = 1 nzmax can be any nunber
!nzmax = 73       ! Number of vertical levels.
nzmax = 451       ! Number of vertical levels.

grid_type = 1     ! Select grid type:
                  ! 1 ==> evenly-spaced grid levels.
                  ! 2 ==> stretched (unevenly-spaced) grid
                  !       entered on thermodynamic grid levels;
                  !       mom. levels halfway between thermo. levels
                  !       (style of SAM stretched grid).
                  ! 3 ==> stretched (unevenly-spaced) grid
                  !       entered on momentum grid levels;
                  !       thermo. levels halfway between mom. levels
                  !       (style of WRF stretched grid).

! Evenly-spaced grid (grid_type = 1).
! Note:  For this option, both zt_grid_fname and zm_grid_fname must be ''.
deltaz_nl  = 20.0    ! Distance between grid levels on evenly-spaced grid.      [m]
zm_init_nl = 0.0 ! Minimum Altitude of lowest momentum level on any grid. [m]
zm_top_nl = 9000.0 ! Maximum Altitude of highest momentum level on any grid. [m]

! Stretched grid on thermodynamic levels (grid_type = 2).
! Path (from ../run_scripts/) and filename of data file that contains
! thermodynamic level altitudes (in meters).
! Note:  For this option, zm_grid_fname must be ''.
!zt_grid_fname = '../input/grid/deep_convection_128lev_27km_zt_grid.grd'
zt_grid_fname = ''

! Stretched grid on momentum levels (grid_type = 3).
! Path (from ../run_scripts/) and filename to data file that contains
! momentum level altitudes (in meters).
! Note:  For this option, zt_grid_fname must be ''.
!zm_grid_fname = '../input/grid/dycore.grd'
zm_grid_fname = ''

day   = 13     ! Day of model start (1 to 31).
month = 6     ! Month of model start (1 to 12).
year  = 1992  ! Year of model start.

lat_vals = 34. ! Latitude (0 to 90 for north lat.; 0 to -90 for south lat.) [degrees]
lon_vals = -25.0   ! Longitude (0 to 180 for east long.; 0 to -180 for west long.) [degrees]

sfc_elevation_nl = 0.0     ! Elevation of ground level [meters above mean sea level]

time_initial = 0.0    ! Model start time [seconds since midnight on start date]
time_final = 144000.0  ! Model end time [seconds since midnight on start date]

dt_main = 60.0   ! Model timestep [s]
dt_rad  = 60.0   ! Radiation timestep [s]

sfctype = 1

p_sfc_nl = 1029.e2     ! Pressure at model base              [Pa]

fcor_nl = 8.1553388e-5    ! Coriolis parameter                  [s^-1]
T0   = 300.        ! Reference temperature (usually 300) [K]

l_uv_nudge     = .true.   ! Flag for horizontal wind speed nudging.

l_t_dependent = .true. ! Read values from _surface.in file. Also read
                       ! _forcings.in file if l_ignore_forcings = .false.

l_modify_bc_for_cnvg_test = .true.  ! Modify surface boundary conditions for convergence test

! Sponge Damping
!
! tau_sponge_damp_min : Minimum damping time-scale ( at the top0 [s]
! tau_sponge_damp_max : Maximum damping time-scale (base of damping layer) [s]
! sponge_damp_depth   : damping depth as a fraction of domain height [-]
!
thlm_sponge_damp_settings%l_sponge_damping = .false.
thlm_sponge_damp_settings%tau_sponge_damp_min = 1800. 
thlm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
thlm_sponge_damp_settings%sponge_damp_depth = 0.25

rtm_sponge_damp_settings%l_sponge_damping = .false.
rtm_sponge_damp_settings%tau_sponge_damp_min = 1800. 
rtm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
rtm_sponge_damp_settings%sponge_damp_depth = 0.25

uv_sponge_damp_settings%l_sponge_damping = .false.
uv_sponge_damp_settings%tau_sponge_damp_min = 1800. 
uv_sponge_damp_settings%tau_sponge_damp_max = 1800. 
uv_sponge_damp_settings%sponge_damp_depth = 0.25


l_restart      = .false.  ! Flag for whether this is a restart run.
! Path (from ../run_scripts/) and filename of GrADS control files
! (without '_zt.ctl' or '_zm.ctl').
restart_path_case  = "restart/astex_a209"
time_restart  = 0.0      ! Time of model restart.
                         ! [seconds since start time listed in GrADS control file]

! Select debug level:
!        0 => Print no debug messages to the screen
!        1 => Print lightweight debug messages, e.g. print statements
!        2 => Print debug messages that require extra testing,
!                e.g. checks for NaNs and spurious negative values.
debug_level = 2

! Add scalars that are computed using higher-order closure
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_sclr_sounding.in 
!	that will be read in and tracked
sclr_dim   = 0  ! Total number of scalars

iisclr_rt  = -1  ! Location in the array to place a scalar like total water
iisclr_thl = -2 ! Location in the array to place a scalar like thetal

sclr_tol_nl  = 1.e-8 ! Tolerances below which we consider scalar to be 0.

! Add scalars that are computed using simple eddy diffusivity
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_edsclr_sounding.in 
!	that will be read in and tracked
edsclr_dim  = 0

iiedsclr_rt     = -1 ! Location in the edsclrm array to place total water
iiedsclr_thl    = -2 ! Location in the edsclrm array to place a scalar like thetal

l_diagnose_correlations = .false.
l_calc_w_corr = .false.
/

&microphysics_setting
! The schemes are: "khairoutdinov_kogan", "morrison", 
! "coamps", "simplified_ice", "none"
!microphys_scheme = "khairoutdinov_kogan",
microphys_scheme = "none",
Nc0_in_cloud  = 100.e6,
hmp2_ip_on_hmm2_ip_intrcpt%rr = 1.0
hmp2_ip_on_hmm2_ip_intrcpt%Nr = 1.0
hmp2_ip_on_hmm2_ip_slope%rr = 0.0
hmp2_ip_on_hmm2_ip_slope%Nr = 0.0
Ncnp2_on_Ncnm2                               = 2.2

C_evap = 0.5,
r_0 = 100.0e-6,

l_cloud_sed = .true.,
! This parameter is only used with l_cloud_sed = .true.
sigma_g = 1.2,
l_var_covar_src = .true.
/

&gfdl_activation_setting

/

&radiation_setting
! The schemes are: "bugsrad", "simplified", "simplified_bomex", or "none"
rad_scheme = "bugsrad"
!rad_scheme = "none",
!l_fix_cos_solar_zen = .true.
l_sw_radiation = .true.
!cos_solar_zen_values= 0.629
!cos_solar_zen_times = 144001. ! Use the same value all simulation
!F0 = 74.0    ! Coefficient for cloud top heating (see Stevens) [W/m^2]
!F1 = 0.      ! Coefficient for cloud base heating (see Stevens)[W/m^2]
!kappa = 130. ! [m^2/kg]
slr = 1.0
sol_const = 1376.0
radiation_top       = 30000 ! Height of profile used for radiation [m]
alvdf = 0.07 ! Visible diffuse surface albedo [-]
alndr = 0.07 ! Near-IR direct surface albedo  [-]
alndf = 0.07 ! Near-IR diffuse surface albedo [-]
alvdr = 0.07 ! Visible direct surface albedo  [-]

! Flag indicating if U.S. Std Atmosphere should 
! be used to extend the grid in radiation.
! If false, sounding data will be used.
l_use_default_std_atmosphere = .true. 

/

&stats_setting

! Statistics Settings

l_stats       = .true.  ! Flag for statistical output.
fname_prefix = "astex_a209"   ! Prefix of output filename.
stats_tsamp  = 60.      ! Frequency of statistical sampling [s]
                        ! For most complete sampling, let 
                        ! stats_tsamp = dt_main
                        ! (stats_tsamp must multiply evenly into stats_tout)
stats_tout   = 300.      ! Frequency of statistical output.  [s]
stats_fmt    = "netcdf"  ! Type of statistical output file ("grads" or "netcdf").

/
```

`./input/case_setups/gabls2_model.in`
```
! $Id$
! Parameter file for GABLS2 case
! Edited by Michael Falk, 28 December 2006

&model_setting

! Model Settings

runtype = "gabls2" ! Case name.

! For grid_type = 1 nzmax can be any nunber
!nzmax = 73       ! Number of vertical levels.
nzmax = 426       ! Number of vertical levels.

grid_type = 1     ! Select grid type:
                  ! 1 ==> evenly-spaced grid levels.
                  ! 2 ==> stretched (unevenly-spaced) grid
                  !       entered on thermodynamic grid levels;
                  !       mom. levels halfway between thermo. levels
                  !       (style of SAM stretched grid).
                  ! 3 ==> stretched (unevenly-spaced) grid
                  !       entered on momentum grid levels;
                  !       thermo. levels halfway between mom. levels
                  !       (style of WRF stretched grid).

! Evenly-spaced grid (grid_type = 1).
! Note:  For this option, both zt_grid_fname and zm_grid_fname must be ''.
!deltaz_nl  = 153.0    ! Distance between grid levels on evenly-spaced grid.      [m]
deltaz_nl  = 20.0    ! Distance between grid levels on evenly-spaced grid.      [m]
!deltaz_nl  = 10.0    ! Distance between grid levels on evenly-spaced grid.      [m]
!deltaz_nl  = 2.0    ! Distance between grid levels on evenly-spaced grid.      [m]
zm_init_nl = 0.0 ! Minimum Altitude of lowest momentum level on any grid. [m]
!zm_top_nl = 4000.0 ! Maximum Altitude of highest momentum level on any grid. [m]
zm_top_nl = 8500.0 ! Maximum Altitude of highest momentum level on any grid. [m]

! Stretched grid on thermodynamic levels (grid_type = 2).
! Path (from ../run_scripts/) and filename of data file that contains
! thermodynamic level altitudes (in meters).
! Note:  For this option, zm_grid_fname must be ''.
zt_grid_fname = ''

! Stretched grid on momentum levels (grid_type = 3).
! Path (from ../run_scripts/) and filename to data file that contains
! momentum level altitudes (in meters).
! Note:  For this option, zt_grid_fname must be ''.
!zm_grid_fname = '../input/grid/dycore.grd'
zm_grid_fname = ''

day   = 22    ! Day of model start (1 to 31).
month = 10    ! Month of model start (1 to 12).
year  = 1999  ! Year of model start.

lat_vals = 37.6   ! Latitude (0 to 90 for north lat.; 0 to -90 for south lat.) [degrees]
lon_vals = -96.7  ! Longitude (0 to 180 for east long.; 0 to -180 for west long.) [degrees]

sfc_elevation_nl = 0.0     ! Elevation of ground level [meters above mean sea level]

time_initial = 0.0       ! Model start time [seconds since midnight UTC on start date]
time_final   = 212400.0  ! Model end time [seconds since midnight UTC on start date]


dt_main    = 60.0   ! Model timestep [s]
dt_rad = 60.0   ! Radiation timestep [s] 

sfctype = 0      ! Select surface scheme (cases without specified surface fluxes).
                 ! 0 ==> fixed surface sensible and latent heat fluxes,
                 !       as defined by sens_ht and latent_ht, respectively.
                 ! 1 ==> determine surface sensible and latent heat fluxes
                 !       through a bulk formula that uses the given surface
                 !       temperature (T_sfc) and assumes over ocean.

p_sfc_nl = 972.e2        ! Pressure at model base              [Pa]

fcor_nl = 8.8983571e-5  ! Coriolis parameter                  [s^-1]
T0   = 283.15        ! Reference temperature (usually 300) [K]

l_modify_bc_for_cnvg_test = .true.  ! Modify surface boundary conditions for convergence test


! Sponge Damping
!
! tau_sponge_damp_min : Minimum damping time-scale ( at the top0 [s]
! tau_sponge_damp_max : Maximum damping time-scale (base of damping layer) [s]
! sponge_damp_depth   : damping depth as a fraction of domain height [-]
!
thlm_sponge_damp_settings%l_sponge_damping = .false.
thlm_sponge_damp_settings%tau_sponge_damp_min = 60. 
thlm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
thlm_sponge_damp_settings%sponge_damp_depth = 0.25

rtm_sponge_damp_settings%l_sponge_damping = .false.
rtm_sponge_damp_settings%tau_sponge_damp_min = 60. 
rtm_sponge_damp_settings%tau_sponge_damp_max = 1800. 
rtm_sponge_damp_settings%sponge_damp_depth = 0.25

uv_sponge_damp_settings%l_sponge_damping = .false.
uv_sponge_damp_settings%tau_sponge_damp_min = 60. 
uv_sponge_damp_settings%tau_sponge_damp_max = 1800. 
uv_sponge_damp_settings%sponge_damp_depth = 0.25

l_restart      = .false.  ! Flag for whether this is a restart run.
! Path (from ../run_scripts/) and filename of GrADS control files
! (without '_zt.ctl' or '_zm.ctl').
restart_path_case  = "restart/gabls2"
time_restart  = 0.0      ! Time of model restart.
                         ! [seconds since start time listed in GrADS control file]

l_input_fields = .false.,

! Select debug level:
!        0 => Print no debug messages to the screen
!        1 => Print lightweight debug messages, e.g. print statements
!        2 => Print debug messages that require extra testing,
!                e.g. checks for NaNs and spurious negative values.
debug_level   = 2

! Add scalars that are computed using higher-order closure
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_sclr_sounding.in
!	that will be read in and tracked 
sclr_dim   = 2  ! Total number of scalars

iisclr_thl = 1 ! Location in the array to place a scalar like thetal
iisclr_rt  = 2  ! Location in the array to place a scalar like total water
iisclr_CO2     = -1  ! Location in the array to place CO2

sclr_tol_nl  = 1.e-2, 1.e-8 ! Tolerances below which we consider scalar to be 0.

! Add scalars that are computed using simple eddy diffusivity
! 0 to shut off all scalar computations
! > 0 the specifies the number of columns in <runtype>_edsclr_sounding.in 
!	that will be read in and tracked
edsclr_dim   = 2

iiedsclr_thl = 1 ! Location in the edsclrm array to place a scalar like thetal

iiedsclr_rt  = 2 ! Location in the edsclrm array to place a scalar like total water 

iiedsclr_CO2 = -1


/

&gfdl_activation_setting

/

&microphysics_setting
! The schemes are: "khairoutdinov_kogan", "morrison", 
! "coamps", "simplified_ice", "none"
microphys_scheme = "none"
/

&radiation_setting
! The schemes are: "bugsrad", simplified", "simplified_bomex", or "none"
rad_scheme = "none",
l_calc_thlp2_rad = .false.,
/

&stats_setting

! Statistics Settings

l_stats       = .true.    ! Flag for statistical output.
fname_prefix = "gabls2"   ! Prefix of output filename.
stats_tsamp  = 60.        ! Frequency of statistical sampling [s]
                          ! For most complete sampling, let 
                          ! stats_tsamp = dt_main
                          ! (stats_tsamp must multiply evenly into stats_tout)
stats_tout   = 600.       ! Frequency of statistical output.  [s]
stats_fmt    = "netcdf"    ! Type of statistical output file ("grads" or "netcdf").
/

! Inputfields Settings
! You will not need to modify these if l_input_fields is false.
&setfields
datafile      = "../les_data/gabls2"
input_type    = "coamps_les"  ! input_type = "clubb" or "coamps_les"
l_input_um = .true. 
l_input_vm = .true.
l_input_rtm = .true. 
l_input_thlm = .true.
l_input_wp2 = .true.
l_input_wprtp = .true.
l_input_wpthlp = .true.
l_input_wp3 = .true.
l_input_rtp2 = .true.
l_input_thlp2 = .true.
l_input_rtpthlp = .true.
l_input_upwp = .true.
l_input_vpwp = .true.
l_input_ug = .false.
l_input_vg = .false.
l_input_rcm = .true.
l_input_wm_zt = .true.
l_input_exner = .true.
l_input_em = .true.
l_input_p = .true.
l_input_rho = .true.
l_input_rho_zm = .true.
l_input_Lscale = .false.
l_input_Lscale_up = .false.
l_input_Lscale_down = .false.
l_input_Kh_zt = .true.
l_input_Kh_zm = .false.
l_input_tau_zm = .false.
l_input_tau_zt = .false.
l_input_wpthvp = .true.
l_input_thl_1 = .false.
l_input_thl_2 = .false.
l_input_mixt_frac = .false.
l_input_s1 = .false.
l_input_s2 = .false.
l_input_stdev_s1 = .false.
l_input_stdev_s2 = .false.
l_input_rc_1 = .false.
l_input_rc_2 = .false.
l_input_thvm = .true.
l_input_rrm = .false.
l_input_Nrm = .false.
l_input_Ncm = .false.
l_input_rsm = .false.
l_input_rim = .false.
l_input_rgm = .false.
l_input_Ncnm = .false.
l_input_Nim = .false.
l_input_thlm_forcing = .false.
l_input_rtm_forcing = .false.
l_input_up2 = .true.
l_input_vp2 = .true.
l_input_sigma_sqd_w = .false.
l_input_cloud_frac = .false.
l_input_sigma_sqd_w_zt = .false.
l_input_veg_T_in_K = .false.
l_input_deep_soil_T_in_K = .false.
l_input_sfc_soil_T_in_K = .false.
/
```
3. Run `./run_scripts/run_scm.bash arm`, `./run_scripts/run_scm.bash astex_a209` and `./run_scripts/run_scm.bash gabls2`


## Reproduce the plots
The scripts mentioned before to reproduce the results produce `netcdf` files as output, these can be used to create the plots shwon in the paper.
The script `generate_paper_plots.py` can be used for that.

1. Generate the the netcdf files with the scripts before and each time copy the output into a different directory. Give the directory the name you want to have in the legend of the plots for each of the three different results (grid adaptation, dycore, hi-res).
2. Go into the `generate_paper_plots.py` script and set the variables `READ_DIR_HI_RES`, `READ_DIR_NO_ADAPT` and `READ_DIR_ADAPT` to the paths, where the directories are with the resulting netcdf files. Then set the output directory, where the plots should be stored, with `OUTPUT_DIR`.
3. Run the script with `python generate_paper_plots.py`
