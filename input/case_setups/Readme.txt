$Id$

INPUT/CASE_SETUPS DIRECTORY OVERVIEW
====================================

The <CASE NAME>_model.in files contain the following namelists:

--------------------------------------------------------------------------------
1. &model_setting:
--------------------------------------------------------------------------------

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

LE | real, default precision
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

l_input_xpwp_sfc | logical
  Flag used to determine whether or not to read in the surface momentum
  fluxes, upwp_sfc and vpwp_sfc

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

--------------------------------------------------------------------------------
2. &microphys_setting:
--------------------------------------------------------------------------------

  Name        | Data type  
------------------------------------------
micro_scheme | character
  The microphysics scheme to use.  Either khairoutdinov_kogan, coamps,
  morrison, or none

l_cloud_sed | logical
  Cloud water sedimentation (K&K or no microphysics)

l_ice_micro | logical
  Compute ice (COAMPS / Morrison)

l_graupel | logical
  Compute graupel (COAMPS / Morrison)

l_hail | logical
  See module_mp_graupel for a description (Morrison)

l_seifert_beheng | logical
  Use Seifert and Beheng (2001) warm drizzle (Morrison)

l_predictnc | logical
  Predict cloud droplet number conc (Morrison)

l_specify_aerosol | logical
  Specify aerosol (Morrison)

l_subgrid_w | logical
  Use subgrid-scale w for cloud droplet activation (Morrison)

l_arctic_nucl | logical
  Use MPACE observations (Morrison)

l_cloud_edge_activation | logical
  Activate on cloud edges (Morrison)

l_fix_pgam | logical
  Fix pgam (Morrison)

l_latin_hypercube_sampling | logical
  Use latin sypercube sampling (K&K).

LH_microphys_calls | integer
  Number of latin hypercube samples to call the microphysics with (K&K).

LH_sequence_length | integer
  Number of timesteps before the latin hypercube seq. repeats (K&K).

l_local_kk | logical
  Use the local formulas for K&K, rather than the analytic formulas.

microphys_start_time | real, minimum 12 digits of precision
  Model time to start calling the microphysics scheme [s]

Ncm_initial 
  Initial value for Ncm in cc/m^3 (K&K, l_cloud_sed, Morrison)

rrp2_on_rrainm2_cloud | real, default precision
rrp2_on_rrainm2_below | real, default precision
  Variance of rain water mixing ratio divided by the mean squared (K&K).

Nrp2_on_Nrm2_cloud | real, default precision
Nrp2_on_Nrm2_below | real, default precision
  Variance of rain water num. conc. divided by the mean squared (K&K).

Ncp2_on_Ncm2_cloud | real, default precision
Ncp2_on_Ncm2_below | real, default precision
  Variance of cloud water num. conc. divided by the mean squared (K&K).

corr_rrNr_LL_cloud | real, default precision
corr_rrNr_LL_below | real, default precision
  Correlation between rain water mixing ratio and number conc. (K&K)

corr_srr_NL_cloud | real, default precision
corr_srr_NL_below | real, default precision
  Correlation between rain water mixing ratio and 's' in Mellor (K&K).

corr_sNr_NL_cloud | real, default precision
corr_sNr_NL_below | real, default precision
  Correlation between rain water num. conc. ratio and 's' in Mellor (K&K).

corr_sNc_NL_cloud | real, default precision
corr_sNc_NL_below | real, default precision
  Correlation between cloud droplet number conc. and 's' in Mellor (K&K).

C_evap | real, default precision
  Khairoutdinov and Kogan (2000) ratio of drizzle drop mean geometric radius to
  mean volume radius (K&K)

r_0 | real, default precision
  Assumed radius of all new drops (K&K) [m].

ccnconst | real, default precision
  Parameter for powerlaw CCN (Morrison) [#/cc].

ccnexpnt | real, default precision
  Exponent for powerlaw CCN (Morrison).

aer_rm1 | real, default precision 
aer_rm2 | real, default precision
  Mean geometric radius (Morrison) [μ].

aer_n1 | real, default precision
aer_n2 | real, default precision
  Aerosol concentration (Morrison) [#/cc].

aer_sig1 | real, default precision
aer_sig2 | real, default precision
  Standard deviation of size distribution (Morrison) [-]

pgam_fixed | real, default precision
  Value for fixed pgam (Morrison).

--------------------------------------------------------------------------------
3. &radiation_setting:
--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------
4. &stats_setting
--------------------------------------------------------------------------------

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
====================================

<Case Name>_sounding.in
  Sounding information for initial CLUBB fields.

<Case Name>_sclr_sounding.in
  Sounding information for passive scalars.

<Case Name>_forcings.in
  Forcing information.

<Case Name>_surface.in
  Surface flux information.
