#!/usr/bin/python

import argparse
import numpy as np
import os
import sys

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
res_dir = '../clubb_baseline/ic_file'

# remind user to load modules if this is Quartz
if ('quartz' in os.uname().nodename):
  print('########################################################')
  print('  make sure to load the mkl and netcdf-fortran modules  ')
  print('########################################################')

# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('case',
  help='name of case to run', choices=('rico', 'bomex', 'wangara', 'dycoms2_rf01', 'dycoms2_rf02'))
parser.add_argument ('-dt', metavar='seconds',
  help='timestep size')
parser.add_argument('-dto', metavar='seconds', dest='dt_output',
  help='time interval between output (default is max[60,dt] if dt specified)')
parser.add_argument('-dz', metavar='meters',
  help='use a uniform grid with specified dz')
parser.add_argument('-Tsfc', dest='Tsfc', action='store_true',
  help='note that the Tscf in the sounding file has been modified (user must have modified file)')
parser.add_argument('-god', metavar='name', dest='godunov',
  help='use specified Godunov-like scheme instead of central difference', choices=('scalarwp3'))
parser.add_argument('-ti', metavar='seconds', dest='tinitial',
  help='time (in seconds) to simulate')
parser.add_argument('-tf', metavar='seconds', dest='tfinal',
  help='time (in seconds) to simulate')
parser.add_argument('-ref', metavar='num', dest='refine',
  help='use a stretched grid with specified number of refinement')
parser.add_argument('-spl', metavar='val', dest='splat',
  help='specify strength of splatting terms (default is 0.0)')
parser.add_argument('-output-name',
  help='use user specified naming convention instead of one generated from other flags')
parser.add_argument('-default-micro', action='store_true',
  help='use microphysics in default case setups instead of turning them off')
parser.add_argument('-default-format', action='store_true',
  help='use output format in default case setup instead of netcdf')
parser.add_argument('-default-bc-height', action='store_true',
  help='use bc height method in default case setup instead of fixed height')
parser.add_argument('-default-aterms', action='store_true',
  help='use aterms in default case setup instead of keeping them in the derivative')
parser.add_argument('-default-rad', action='store_true',
  help='use radiation scheme in default case setup instead of turning it off')
parser.add_argument('-skip-check', action='store_true',
  help='skip the pause in the script that allows user to inspect configuration changes')
parser.add_argument('-warm-run', action='store_true',
  help='do restart run instead of cold initialization via souding profiles')
parser.add_argument ('-twarm', metavar='seconds',
  help='warm start time')

# create dictionary of parameters to change (remove None and False vals)
parameters = vars(parser.parse_args(sys.argv[1:]))
unset_keys = []
for key, val in parameters.items():
  if (val is None or val is False):
    unset_keys.append(key)
for key in unset_keys:
  del parameters[key]

# set case name
model_file_name = parameters['case'] + '_model.in'

# read in default configuration from model file
model_file = open(os.path.join(case_dir, model_file_name), 'r')
model_config = model_file.readlines()
model_file.close()

# create concatenated model file name from specified parameters
# (unless user specified their own name)
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
    if (parameter != 'dt' and parameter != 'dz' and parameter != 'refine' and parameter != 'godunov' and parameter != 'output_name'):
      continue
    if (isinstance(parameters[parameter], bool) and parameters[parameter]):
      model_file_name += '_' + parameter.replace('_','-')
    elif (parameters[parameter]):
      model_file_name += '_' + parameter + '-' + parameters[parameter]

model_file_name += '.in'

# create modified configuration model file
model_file = open(model_file_name, 'w')

# create dictionary of configuration strings from parameters
config_strings = {}
config_strings['fname_prefix'] = '"' + model_file_name.replace('.in','') + '"'

# set microphysics to "none" unless user has specified something else
if ('default_micro' not in parameters):
  config_strings['microphys_scheme'] = '"none"'

# set restart run information
if ('warm_run' in parameters):
  config_strings['l_restart'] = '.true.'
  if ('twarm' in parameters):
    config_strings['time_restart'] = parameters['twarm']
    if ('tinitial' in parameters):
      config_strings['time_restart'] = str(float(parameters['tinitial']) + float(parameters['twarm']))
  else:
    config_strings['time_restart'] = '600.0'
  config_strings['restart_path_case'] = '"'+ os.path.join(res_dir,parameters['case']) +'"'

# set output to netcdf unless user has specified something else
if ('default_format' not in parameters):
  config_strings['stats_fmt'] = '"netcdf"'

# set l_standard_term_ta to true unless user has specified something else
if ('default_aterms' not in parameters):
  config_strings['l_standard_term_ta'] = '.true.'

# fix flux computation height unless user has specified otherwise
if ('default_bc_height' not in parameters):
  config_strings['l_bc_at_constant_height'] = '.true.'

# turn off radiation unless user has specified something else
if ('default_rad' not in parameters):
  config_strings['rad_scheme'] = '"none"'
  config_strings['l_calc_thlp2_rad'] = '.false.'
elif ('radiation' not in parameters):
  parameters['radiation'] = '"none"'
  if parameters['radiation'].strip('\"') == 'none':
    parameters['l_calc_thlp2_rad'] = '.false.'
  else:
    parameters['l_calc_thlp2_rad'] = '.true.'

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
  if (parameters['case'] == 'bomex'):
    height = 3400.0
    config_strings['grid_type'] = '2'
  elif (parameters['case'] == 'rico'):
    height = 9900.0
  elif (parameters['case'] == 'wangara'):
    config_strings['grid_type'] = '2'
    height = 2200.0
  elif (parameters['case'] == 'dycoms2_rf02'):
    config_strings['grid_type'] = '2'
    height = 1500.0
  elif (parameters['case'] == 'dycoms2_rf01'):
    config_strings['grid_type'] = '2'
    height = 5000.0
  else:
    sys.exit('refinement only implemented for bomex, rico, dycoms_rf01, dycoms_rf02 and wangara cases')
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
  config_strings['nzmax'] = str(len(grid))
  grid_filename = '{}_convection_{}lev_27km_zt_grid.grd'.format(parameters['case'],config_strings['nzmax'])
  np.savetxt(grid_filename, grid)
  config_strings['zt_grid_fname'] = f"'{grid_filename}'"

# set initial time if specified
if ('tinitial'in parameters):
  config_strings['time_initial'] = parameters['tinitial']
# set final time if specified
if ('tfinal' in parameters):
  config_strings['time_final'] = parameters['tfinal']

# set no splatting unless user has specified something else
if ('splat' in parameters):
  config_strings['C_wp2_splat'] = parameters['splat']
else:
  config_strings['C_wp2_splat'] = '0.0'

# set Godunov upwinding, if specified
if ('godunov' in parameters):
  config_strings['l_godunov_upwind_wp3_ta'] = '.true.'
  config_strings['l_godunov_upwind_wpxp_ta'] = '.true.'
  config_strings['l_godunov_upwind_xpyp_ta'] = '.true.'
  config_strings['l_upwind_xpyp_ta'] = '.false.'

#set perturbation growth test, if specified 
if ('pergro' in parameters):
  config_strings['l_perturb_IC_at_rounding_level'] = '.true.'
  config_strings['perturb_factor'] = '10.0'
else:
  config_strings['l_perturb_IC_at_rounding_level'] = '.false.'
  config_strings['perturb_factor'] = '0.0'

# warn about Tsfc parameter
if ('Tsfc' in parameters):
  print("WARNING: Specifying Tsfc doesn't change model parammeters...")
  print("         this needs to be done in " + parameters['case'] + "_sounding.in")
  input("Press Enter to acknowledge")

modified_lines = []
for line in model_config:
  modified = False
  if (not line.strip().startswith('!')):
    for key, val in config_strings.items():
      # adjust values
      if (line.startswith(key)):
        default_value = line.split()[2]
        line = line.replace(default_value, val)
        modified = True

  # write line to model configuration (and record if modified)
  model_file.write(line)
  if (modified):
    modified_lines.append(line)

model_file.close()

print('Wrote ' + model_file_name + ' file with the following modified lines:\n')
for line in modified_lines:
  print(line),

# concatenate all other input namelist files to one file
line_collections = []
# "parameter" file
input_file = open(os.path.join(tunable_dir, 'tunable_parameters.in'), 'r')
input_lines = input_file.readlines()
input_file.close()
line_collection = []
for line in input_lines:
  if (line.startswith('C_wp2_splat')):
    default_value = line.split()[2]
    line = line.replace(default_value, config_strings['C_wp2_splat'])
    print('Setting splat coefficient to {}'.format(config_strings['C_wp2_splat']))
    print(line)
  line_collection.append(line)
line_collections.append(line_collection)
# "SILHS_PARAMS" file
input_file = open(os.path.join(tunable_dir, 'silhs_parameters.in'), 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "FLAGS" file
input_file = open(os.path.join(tunable_dir, 'configurable_model_flags.in'), 'r')
input_lines = input_file.readlines()
input_file.close()
line_collection = []
for line in input_lines:
  for key, val in config_strings.items():
    if (line.startswith(key)):
      default_value = line.split()[2]
      line = line.replace(default_value, val)
      print(f'Setting {key} flag to {val}')
      print(line)
  line_collection.append(line)
line_collections.append(line_collection)
# "MOD_MODEL" file
input_file = open(model_file_name, 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "stats" file
input_file = open(os.path.join(stats_dir,'standard_stats.in'), 'r')
line_collections.append(input_file.readlines())
input_file.close()
# write out to namelist file
namelist_file = open('clubb.in', 'w')
for line_collection in line_collections:
  for line in line_collection:
    ind = line.find('!')
    if (ind != -1):
      line = line[0:ind] + '\n'
    namelist_file.write(line)
namelist_file.close()

# execute CLUBB
if ('skip_check' not in parameters):
  input('Press Enter to run CLUBB...')
os.system(os.path.join(bin_dir, 'clubb_standalone'))
