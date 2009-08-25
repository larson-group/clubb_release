$Id$

INPUT/CASE_SETUPS DIRECTORY OVERVIEW
====================================

The <CASE NAME>_model.in files contain the following namelists:

1. &model_setting:

  Name        | Data type  
------------------------------------------

runtype | character 
  This a simulation type, e.g. "bomex".

nzmax | integer
 Number of vertical levels.

grid_type | integer
  1 ==> evenly spaced grid levels, 2 ==> unevenly-spaced thermo. grid levels,
  3 ==> unevenly-spaced momentum grid levels.

deltaz | real, default precision 
  Distance between grid levels for grid_type = 1 [m].

zm_init | real, default precision 
  Min Altitude of lowest mom. level on any grid  [m].

zm_top | real, default precision 
  Max Altitude of highest mom. level on any grid [m].

zt_grid_fname | character
  Stretched grid on thermodynamic levels (grid_type = 2).
  Path (from ../run_scripts/) and filename of data file that contains
  thermodynamic level altitudes (in meters).

zm_grid_fname | character
  Stretched grid on momentum levels (grid_type = 3).
  Path (from ../run_scripts/) and filename to data file that contains
  momentum level altitudes (in meters).
  Note:  For this option, zt_grid_fname must be ''.

day | integer
  Day of model start (1 to 31).

month | integer
  Month of model start (1 to 12).

year | integer
  Year of model start.

rlat | real, default precision
  Latitude (0 to 90 for north lat.; 0 to -90 for south lat.) [degrees].

rlon | real, default precision
  Longitude (0 to 180 for east long.; 0 to -180 for west long.) [degrees].

time_initial | real, minimum of 12 digits of precision
  Model start time [seconds since midnight on start date].

time_final | real, minimum of 12 digits of precision
  Model end time [seconds since midnight on start date].

time_spinup | real, minimum of 12 digits of precision
  End of spin-up period [seconds since midnight on start date].

dtmain | real, minimum of 12 digits of precision 
  Model timestep [s].

dtclosure | real, minimum of 12 digits of precision
  Closure timestep [s] (must multiply evenly into dtmain).

Tsfc | real, default precision
  Temperature at model base [K].

psfc | real, default precision
  Pressure at model base [Pa].

SE | real, default precision
  Fixed surface sensible heat flux  [K m/s].

LE ! real, default precision
  Fixed surface latent heat flux [kg/kg m/s].

fcor | real, default precision
  Coriolis parameter [s^-1].

T0 | real, default precision
  Reference temperature (usually 300) [K].

ts_nudge | real, default precision 
  Timescale of u/v nudging [s].

l_t_dependent | logical
  Flag used to indicate that this case uses a forcings.in and surface.in file for
  time dependent input.

l_tke_aniso | logical
  Flag for anisotropic tke.

l_uv_nudge | logical
  Flag for horizontal wind speed nudging.

l_use_default_std_atmosphere | logical 
  Flag indicating if U.S. Std Atmosphere should be used to extend the grid in radiation.
  If false, sounding data will be used.

rtm/thlm/um/vm_sponge_damp_settings%l_sponge_damping | logical
  Enable sponge damping.

rtm/thlm/um/vm_sponge_damp_settings%tau_sponge_damp_min | real, default precision
  Minimum damping time-scale ( at the top ) [s].

rtm/thlm/um/vm_sponge_damp_settings%tau_sponge_damp_max | real, default precision
  Maximum damping time-scale (base of damping layer) [s].

rtm/thlm/um/vm_sponge_damp_settings%sponge_damp_depth   | real, default precision
   damping depth as a fraction of domain height [-].

l_restart | logical 
  Flag for whether this is a restart run.

restart_path_case | character 
  Path (from ../run_scripts/) and filename of GrADS/netCDF files 
  sans '_zt.ctl', '_zm.ctl', or '_zt.nc', '_zm.nc'.

time_restart | real, minimum of 12 digits of precision
  Time of model restart.  [seconds since start time in stat file]

debug_level | integer
  0 => Print no debug messages to the screen
  1 => Print lightweight debug messages, e.g. print statements
  2 => Print debug messages that require extra testing,
       e.g. checks for NaNs and spurious negative values.

sclr_dim  | integer
  Adds passive scalars that are computed using higher-order closure.
  0 to shut off all high-order scalar computations.
  > 0 the specifies the number of columns in <runtype>_sclr_sounding.in 
      that will be read in and tracked.

iisclr_rt  | integer
  Location in the sclr arrays to place a scalar emulating total water.

iisclr_thl | integer 
  Location in the sclr arrays to place a scalar emulating thetal.

iisclr_CO2 | integer
  Location in the sclr arrays to place a scalar for CO2.

sclr_tol | real, default precision, array 
  Tolerances below which we consider the scalar to be 0.

edsclr_dim  | integer
  Adds scalars that are computed using simple eddy-diffusivity.
  0 to shut off all eddy-scalar computations.
  > 0 the specifies the number of columns in <runtype>_edsclr_sounding.in 
      that will be read in and tracked.

iiedsclr_rt | integer
  Location in the edsclr arrays to place a scalar emulating total water.

iiedsclr_thl | integer
  Location in the edsclr arrays to place a scalar emulating thetal.

iiedsclr_CO2 | integer
  Location in the edsclr arrays to place a scalar for CO2.

2. &microphys_setting:

  Name        | Data type  
------------------------------------------
TODO

3. &radiation_setting:

  Name        | Data type  
------------------------------------------

rad_scheme | character
  Currently only "bugsrad" or "simplified".

sol_const | double precision
  The solar constant [W/m^2].

radiation_top | real, default precision
  The top of the atmosphere fed into a radiation scheme [m].
  The computational grid should be extended to reach this altitude.

alvdr | double precision
  Visible direct surface albedo   [-].

alndr | double precision
  Near-IR direct surface albedo   [-].

alvdf | double precision
  Visible diffuse surface albedo  [-].

alndf| double precision
  Near-IR diffuse surface albedo  [-].

F0 | real, default precision
  Coefficient for cloud top heating (see Stevens) [W/m^2].

F1 | real, default precision
  Coefficient for cloud base heating (see Stevens) [W/m^2].

kappa | real, default precision
  A constant (Duynkerke eqn. 5) [m^2/kg].

gc | real, default precision
  Asymmetry parameter, "g" in Duynkerke [-].

omega |  real, default precision
  Single-scattering albedo [-].

slr | double precision
  Fraction of daylight (usually 1.0)  [-].

amu0 | double precision
  Cosine of the solar zenith angle [-].

eff_drop_radius | real, default precision
  Effective droplet radius  [m].

4. &stats_setting

  Name        | Data type  
------------------------------------------

l_stats | logical
  Enables statistics output

fname_prefix | character
  Prefix of the output filenames (a _zt, _zm, and _sfc will be generated).

stats_tsamp | real, minimum 12 digits of precision
  Frequency to sample statistics [s].

stats_tout | real, minimum 12 digits of precision
  Frequency to output statistics [s]

stats_fmt | character
  Either "grads" or "netcdf".

Other files in the case_setups directory:

<Case Name>_sounding.in
  Sounding information for initial CLUBB fields.

<Case Name>_sclr_sounding.in
  Sounding information for passive scalars.

<Case Name>_forcings.in
  Forcing information.

<Case Name>_surface.in
  Surface flux information.
