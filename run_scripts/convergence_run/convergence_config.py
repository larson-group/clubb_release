#!/usr/bin/python
#######################################################################
# $Id$
#
# Description:
# Script to configure and run convergence test simulations.   
# The script only well tested for BOMEX, RICO, DYCOMS2_RF02 and Wangara cases 
# original script developer: Chris Vogl (vogl2@llnl.gov)
# revised by Shixuan Zhang (shixuan.zhang@pnnl.gov) 
#
#######################################################################
import argparse
import numpy as np
import os
import sys
import shutil
from convergence_function import modify_ic_profile

# check that Python 3 is being used
if (sys.version_info.major < 3):
  sys.exit('must use Python 3 instead of {}'.format(sys.version))

# get CLUBB root directory (assumed to be one above the CWD) and set directories
clubb_dir = os.path.join(os.getcwd(),'..')
output_dir = os.path.join(clubb_dir, 'output')
case_dir = os.path.join(clubb_dir, 'input', 'case_setups')
grid_dir = os.path.join(clubb_dir, 'input', 'grid')
tunable_dir = os.path.join(clubb_dir, 'input', 'tunable_parameters')
stats_dir = os.path.join(clubb_dir, 'input', 'stats')
bin_dir = os.path.join(clubb_dir, 'bin')
res_dir = os.path.join(clubb_dir, 'restart_ic') 

# remind user to load modules if this is Quartz
if ('quartz' in os.uname().nodename):
  print('########################################################')
  print('  make sure to load the mkl and netcdf-fortran modules  ')
  print('########################################################')

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('case',
  help='name of case to run', choices=('rico', 'bomex', 'wangara', 'dycoms2_rf02_nd'))
parser.add_argument ('-dt', metavar='seconds',
  help='timestep size')
parser.add_argument('-dto', metavar='seconds', dest='dt_output',
  help='time interval between output (default is max[600,dt] if dt specified)')
parser.add_argument('-dz', metavar='meters',
  help='use a uniform grid with specified dz')
parser.add_argument('-Tsfc', dest='Tsfc', action='store_true',
  help='note that the Tscf in the sounding file has been modified (user must have modified file)')
parser.add_argument('-ti', metavar='seconds', dest='tinitial',
  help='time (in seconds) to simulate')
parser.add_argument('-tf', metavar='seconds', dest='tfinal',
  help='time (in seconds) to simulate')
parser.add_argument('-ref', metavar='num', dest='refine',
  help='use a stretched grid with specified number of refinement')
parser.add_argument('-output-name',dest='output_name',
  help='use user specified naming convention instead of one generated from other flags')
parser.add_argument('-micro-off', dest='turn_off_microphysics', action='store_true',
  help='turn off microphysics instead of using default case setups')
parser.add_argument('-rad-off', dest='turn_off_radiation', action='store_true',
  help='turn off radiation scheme instead of using default case setup')
parser.add_argument('-binary-out', dest='default_format', action='store_true',
  help='use output format in default case setup (binary) instead of netcdf')
parser.add_argument('-new-bc', dest='modified_bc', action='store_true',
  help='use boundary condition at fixed height instead of method in default case setup (space dependent)')
parser.add_argument('-new-ic', dest='modified_ic', action='store_true',
  help='use smoothed initial condition instead of method in default case setup (linear interpolation)')
parser.add_argument('-fix-fc', dest='fixed_forcing', action='store_true',
  help='use non-time-dependent forcing instead of method in default case setup (time-dependent)')
parser.add_argument('-splat-off', dest='turn_off_splat', action='store_true',
  help='use nonzero spatting terms in default case setup instead of setting them to zero') 
parser.add_argument('-standard-aterms', dest='standard_aterms', action='store_true',
  help='use aterms in default case setup instead of keeping them in the derivative')
parser.add_argument('-smooth-tau', dest='smoothed_tau', action='store_true',
  help='use formula of tau (for wpxp) in default case setup without smoothed Heavidise function')
parser.add_argument('-lin-diff', dest='linear_diffusion', action='store_true',
  help='use nonlinear diffusion in default case setup instead of linear diffusion')
parser.add_argument('-new-lim', dest='modified_BVF_Ri_limiter', action='store_true',
  help='use modified limiters on BVF and Ri instead of default case setup')
parser.add_argument('-new-wp3cl', dest='modified_wp3_clip', action='store_true',
  help='use modified skewness clippings on wp3 instead of default case setup') 
parser.add_argument('-godwp3', metavar='name', dest='godunov_wp3',
  help='use specified Godunov-like scheme instead of central difference (for wp3 equation)')
parser.add_argument('-godxpyp', metavar='name', dest='godunov_xpyp',
  help='use specified Godunov-like scheme instead of central difference (for xpyp equation)')
parser.add_argument('-skip-check', action='store_true',
  help='skip the pause in the script that allows user to inspect configuration changes')
parser.add_argument('-warm-init', dest='restart_run', action='store_true',
  help='do restart run instead of cold initialization using an existing simulation')

# create dictionary of parameters to change (remove None and False vals)
parameters = vars(parser.parse_args(sys.argv[1:]))
unset_keys = []
for key, val in parameters.items():
  if (val is None or val is False):
    unset_keys.append(key)
for key in unset_keys:
  del parameters[key]

# concatenate all other input namelist files to one file
line_collections = []
# "parameter" file
input_file = open(os.path.join(tunable_dir, 'tunable_parameters.in'), 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "SILHS_PARAMS" file
input_file = open(os.path.join(tunable_dir, 'silhs_parameters.in'), 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "FLAGS" file
input_file = open(os.path.join(tunable_dir, 'configurable_model_flags.in'), 'r')
input_lines = input_file.readlines()
input_file.close()
for line in input_lines:
  if (not line.strip().startswith('!')):
    ind = line.find('!')
    if (ind != -1):
      line = line[0:ind] + '\n'
    line_collections.append(line)
    if (line.startswith('tridiag_solve_method')):
      if ('smoothed_tau' in parameters):
        if (not any('l_smooth_Heaviside_tau_wpxp' in tmpstr for tmpstr in input_lines)):
          lnew = 'l_smooth_Heaviside_tau_wpxp = .false.,' + '\n'
          line_collections.append(lnew)
          print('append line '+lnew)
      if ('linear_diffusion' in parameters):
        if (not any('l_linear_diffusion' in tmpstr for tmpstr in input_lines)):
          lnew = 'l_linear_diffusion = .false.,' + '\n'
          line_collections.append(lnew)
          print('append line '+lnew)
      if ('modified_BVF_Ri_limiter' in parameters):
        if (not any('l_use_modify_limiters' in tmpstr for tmpstr in input_lines)):
          lnew = 'l_use_modify_limiters = .false.,' + '\n'
          line_collections.append(lnew)
          print('append line '+lnew)
      if ('modified_wp3_clip' in parameters):
        if (not any('l_smooth_Heaviside_wp3_lim' in tmpstr for tmpstr in input_lines)):
          lnew = 'l_smooth_Heaviside_wp3_lim = .false.,' + '\n'
          line_collections.append(lnew)
          print('append line '+lnew)
# "MOD_MODEL" file
model_file_name = parameters['case'] + '_model.in'
input_file = open(os.path.join(case_dir, model_file_name), 'r')
input_lines = input_file.readlines()
input_file.close()
for line in input_lines:
  if (not line.strip().startswith('!')):
    ind = line.find('!')
    if (ind != -1):
      line = line[0:ind] + '\n'
    line_collections.append(line)
    if (line.startswith('runtype')):
      if (not any('nzmax' in tmpstr for tmpstr in line_collections)):
        lnew = 'nzmax = 999' + '\n'
        line_collections.append(lnew)
        print('append line '+lnew)
      if ('modified_bc' in parameters):
        if (not any('l_modify_bc_for_cnvg_test' in tmpstr for tmpstr in input_lines)):
          lnew = 'l_modify_bc_for_cnvg_test = .false.,' + '\n'
          line_collections.append(lnew)
          print('append line '+lnew)
    elif (line.startswith('rad_scheme')):
      if (not any('l_calc_thlp2_rad' in tmpstr for tmpstr in input_lines)):
        lnew = 'l_calc_thlp2_rad = .true.,' + '\n'
        if (line.split()[2] == '"none",'):
          lnew = 'l_calc_thlp2_rad = .false.,' + '\n'
        line_collections.append(lnew)
        print('append line '+lnew)
      if (not any('l_rad_above_cloud' in tmpstr for tmpstr in input_lines)):
        lnew = 'l_rad_above_cloud = .true.,' + '\n'
        if (line.split()[2] == '"none",'):
          lnew = 'l_rad_above_cloud = .false.,' + '\n'
        line_collections.append(lnew)
        print('append line '+lnew)
    elif (line.startswith('microphys_scheme')):
      if (not any('l_var_covar_src' in tmpstr for tmpstr in input_lines)):
        lnew = 'l_var_covar_src = .true.' + '\n'
        if (line.split()[2] == '"none"'):
          lnew = 'l_var_covar_src = .false.' + '\n'
        line_collections.append(lnew)
        print('append line '+lnew)
      if (not any('l_cloud_sed' in tmpstr for tmpstr in input_lines)):
        lnew = 'l_cloud_sed = .true.' + '\n'
        if (line.split()[2] == '"none"'):
          lnew = 'l_cloud_sed = .false.' + '\n'
        line_collections.append(lnew)
        print('append line '+lnew)
# "stats" file
input_file = open(os.path.join(stats_dir,'standard_stats.in'), 'r')
line_collections.append(input_file.readlines())
input_file.close()
# write out to namelist file
out_file_name = parameters['case'] + '_default.in' 
namelist_file = open(out_file_name, 'w')
for line_collection in line_collections:
  for line in line_collection:
    ind = line.find('!')
    if (ind != -1):
      line = line[0:ind] + '\n'
    namelist_file.write(line)
namelist_file.close()

# read the default namelist setup and modify it for convergence test 
model_file      = open(out_file_name, 'r')
model_config    = model_file.readlines()
model_file.close()
# create a file to save the modified namelist 
# use concatenated model file name from specified parameters
# unless user specified their own name
if ('output_name' not in parameters):
  model_file_name = model_file_name.replace('_model.in', '')
  for parameter in parameters:
    if (parameter == 'case' or parameter == 'skip_check'):
      continue
    if (isinstance(parameters[parameter], bool) and parameters[parameter]):
      model_file_name += '_' + parameter.replace('_','-')
    elif (parameters[parameter]):
      model_file_name += '_' + parameter + '-' + parameters[parameter]
else:
  model_file_name = model_file_name.replace('_model.in', '')
  for parameter in parameters:
    if (parameter != 'dt' and parameter != 'dz' and parameter != 'refine' and parameter != 'godunov_wp3' and 
        parameter != 'godunov_xpyp' and parameter != 'output_name'):
      continue
    if (isinstance(parameters[parameter], bool) and parameters[parameter]):
      model_file_name += '_' + parameter.replace('_','-')
    elif (parameters[parameter]):
      model_file_name += '_' + parameter + '-' + parameters[parameter]

# create modified configuration model file
model_file_name += '.in'
model_file = open(model_file_name, 'w')

# create dictionary of configuration strings from parameters
config_strings = {}
config_strings['fname_prefix'] = '"' + model_file_name.replace('.in','') + '"'

# set output to netcdf unless user has specified something else
if ('default_format' not in parameters):
  config_strings['stats_fmt'] = '"netcdf"'

# set microphysics to "none" unless user has specified something else
if ('turn_off_microphysics' in parameters):
  config_strings['microphys_scheme'] = '"none"'
  config_strings['l_var_covar_src'] = '.false.'
  config_strings['l_cloud_sed'] = '.false.'

# turn off radiation unless user has specified something else
if ('turn_off_radiation' in parameters):
  config_strings['rad_scheme']        = '"none",'
  config_strings['l_calc_thlp2_rad']  = '.false.,'
  config_strings['l_rad_above_cloud'] = '.false.,'

# set l_standard_term_ta to true unless user has specified something else
if ('standard_aterms' in parameters):
  config_strings['l_standard_term_ta'] = '.true.'

# fix flux computation height unless user has specified otherwise
if ('modified_ic' in parameters):
  #config_strings['l_modify_ic_with_cubic_int'] = '.true.,'
  if ('dycoms2_rf02' in parameters['case']):
    case_name = 'dycoms2_rf02'
  else: 
    case_name = parameters['case'] 
  if ('dz' in parameters):
    case_dz = 1.0
    case_ref  = -9999
  elif ('refine' in parameters):
    case_dz  = -9999.0
    case_ref = 7 
  else:
    sys.exit('must specify grid spacing or refinement level info.....') 
  modify_ic_profile(clubb_dir, case_dir, grid_dir, case_name, case_dz, case_ref)

# used revised boundary condition if user has specified 
if ('modified_bc' in parameters):
  config_strings['l_modify_bc_for_cnvg_test'] = '.true.,'

#set time-dependent forcing to false
if ('fixed_forcing' in parameters):
  config_strings['l_ignore_forcings'] = '.true.'

# set no splatting unless user has specified something else
if ('turn_off_splat' in parameters):
  config_strings['C_wp2_splat'] = '0.0'

# set l_smooth_Heaviside_tau_wpxp to true unless unless user has specified otherwise
if ('smoothed_tau' in parameters):
  config_strings['l_smooth_Heaviside_tau_wpxp'] = '.true.,'

# use linear diffusion instead of nonlinear diffusion unless user has specified otherwise
if ('linear_diffusion' in parameters):
  config_strings['l_linear_diffusion'] = '.true.,'
  #specific setup for coefficients in linear diffusion 
  #(used in clubb convergence paper) 
  config_strings['c_K1'] = '0.000000'
  config_strings['c_K2'] = '0.000000'
  config_strings['c_K6'] = '0.000000'
  config_strings['c_K8'] = '0.000000'
  config_strings['c_K9'] = '0.000000'
  config_strings['nu1']  = '21.00000'
  config_strings['nu2']  = '1.000000'
  config_strings['nu6']  = '5.000000'
  config_strings['nu8']  = '36.00000'
  config_strings['nu9']  = '10.00000'
  #slightly increase the dissipation term for wpxp 
  config_strings['C_invrs_tau_N2_wpxp'] = '0.02'

# use modified setup for limiters on BVF and Ri unless user has specified otherwise
if ('modified_BVF_Ri_limiter' in parameters):
  config_strings['l_use_modify_limiters'] = '.true.,'

if ('modified_wp3_clip' in parameters): 
  config_strings['l_smooth_Heaviside_wp3_lim'] = '.true.' 

# set restart run information
if ('restart_run' in parameters):
  config_strings['l_restart'] = '.true.'
  config_strings['restart_path'] = '"'+ os.path.join(res_dir,parameters['case']) +'"'
  if ('tinitial' in parameters):
    config_strings['time_restart'] = str(float(parameters['tinitial']) + 3600.0)
  else: 
    config_strings['time_restart'] = '3600.0'

# match dt_output to max of dt or 60s unless it is specifically provided
if ('dt' in parameters):
  config_strings['dt_main'] = parameters['dt']
  config_strings['dt_rad'] = parameters['dt']
  if ('dt_output' not in parameters):
    parameters['dt_output'] = str(max(float(parameters['dt']), 600.0))
if ('dt_output' in parameters):
  config_strings['stats_tsamp'] = parameters['dt_output']
  config_strings['stats_tout'] = parameters['dt_output']

# if dz specified, set other associated options
if ('dz' in parameters):
  config_strings['deltaz'] = parameters['dz']
  config_strings['grid_type'] = '1'
  config_strings['zt_grid_fname'] = "''"

# if refinement specified, create grid file and set appropriate file name
if ('refine' in parameters):
  # read namelist file to obtain model height in case setup 
  config_strings['grid_type'] = '2'
  tmp_file_name = parameters['case'] + '_model.in' 
  input_file = open(os.path.join(case_dir, tmp_file_name), 'r')
  input_lines = input_file.readlines()
  input_file.close()
  height = -999.0 
  for line in input_lines:
    if ('zm_top' in line):
      height = float(line.split()[2])
  if (height == -999.0): 
    sys.exit('Model height are not prescribed in namelist, space refinement failed...') 
  #generate refined grid 
  refine = int(parameters['refine'])
  grid128 = np.loadtxt(os.path.join(grid_dir,'deep_convection_128lev_27km_zt_grid.grd'))
  ind = np.where(grid128 > height)[0][0]
  grid128 = grid128[:ind]
  grid128[0] = 0.0
  size = (len(grid128)-1)*2**refine + 1
  grid = np.empty((size))
  coarse_ind = np.arange(len(grid128))
  refine_ind = np.arange(0,size,2**refine)
  grid[refine_ind] = grid128[coarse_ind]
  for level in range(refine):
    coarse_ind = np.arange(0, size, 2**(refine-level))
    refine_ind = coarse_ind[:-1] + 2**(refine-level-1)
    grid[refine_ind] = 0.5*(grid[coarse_ind[1:]] + grid[coarse_ind[:-1]])
  grid[0] = -grid[1]
  #specify the levels of refined grid in namelist file 
  config_strings['nzmax'] = str(len(grid))
  #save grid file and specify refined grid in namelist file 
  grid_filename = '{}_convection_{}lev_27km_zt_grid.grd'.format(parameters['case'],config_strings['nzmax'])
  np.savetxt(grid_filename, grid)
  config_strings['zt_grid_fname'] = f"'{grid_filename}'"

# set initial time if specified
if ('tinitial'in parameters):
  config_strings['time_initial'] = parameters['tinitial']
# set final time if specified
if ('tfinal' in parameters):
  config_strings['time_final'] = parameters['tfinal']

# set Godunov upwinding, if specified
if ('godunov_wp3' in parameters):
  config_strings['l_godunov_upwind_wp3_ta'] = '.true.'

if ('godunov_xpyp' in parameters):
  config_strings['l_godunov_upwind_xpyp_ta'] = '.true.'
  config_strings['l_upwind_xpyp_ta'] = '.false.'

# warn about Tsfc parameter
if ('Tsfc' in parameters):
  print("WARNING: Specifying Tsfc doesn't change model parammeters...")
  print("         this needs to be done in " + parameters['case'] + "_sounding.in")
  input("Press Enter to acknowledge")

modified_lines = []
for line in model_config:
  modified = False
  # Apply changes in model configurations 
  if (not line.strip().startswith('!')):
    for key, val in config_strings.items():
      # adjust values
      if (line.startswith(key)):
        if (line.split()[0] == key): 
          default_value = line.split()[2]
          line = line.replace(default_value, val)
          modified = True
          print(f'Setting {key} to {val}')
          print(line)

  # write line to model configuration (and record if modified)
  if (not line.strip().startswith('!')):
    model_file.write(line)
  if (modified):
    modified_lines.append(line)

model_file.close()

print('Wrote ' + model_file_name + ' file with the following modified lines:\n')
for line in modified_lines:
  print(line),

# execute CLUBB
shutil.copy(model_file_name, 'clubb.in')
if ('skip_check' not in parameters):
  input('Press Enter to run CLUBB...')
os.system(os.path.join(bin_dir, 'clubb_standalone'))
