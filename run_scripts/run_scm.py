#!/usr/bin/python

import numpy as np
import os
import sys

modifiable_parameters = ['dt', 'dt_output', 'microphysics', 'format', 'prefix', 'dz',
                         'Tsfc', 'godunov', 'aterms', 'levels', 'tfinal', 'refine',
                         'splat']

# check that Python 3 is being used
if (sys.version_info.major < 3):
  print('must use Python 3 instead of {}'.format(sys.version))
  sys.exit(1)

# remind user to load modules if this is Quartz
if ('quartz' in os.uname().nodename):
  print('########################################################')
  print('  make sure to load the mkl and netcdf-fortran modules  ')
  print('########################################################')

# TODO: check that this is being run from the run_scripts directory
os.chdir('../output')

if (len(sys.argv) < 2):
  print('usage: ./run_scm.py model_name ... option_name option value ...')
  print('valid options:')
  print(modifiable_parameters)
  sys.exit(1)

model_file_name = sys.argv[1] + '_model.in'

# read in default configuration from model file
model_file = open('../input/case_setups/' + model_file_name, 'r')
model_config = model_file.readlines()
model_file.close()

# read in parameters to modify from command line
parameters = {}
for n in range(2,len(sys.argv)-1,2):
  parameters[sys.argv[n]] = sys.argv[n+1]

# create concatenated model file name (and check parameters are valid)
model_file_name = model_file_name.replace('_model.in', '')
for parameter in parameters:
  model_file_name += '_' + parameter + '-' + parameters[parameter]
  if (parameter not in modifiable_parameters):
    print(parameter + ' is not modifiable')
    print('valid options:')
    print(modifiable_parameters)
    sys.exit(1)
model_file_name += '.in'

# create modified configuration model file
model_file = open(model_file_name, 'w')

# set microphysics to "none" unless user has specified something else
if ('microphysics' not in parameters):
  parameters['microphysics'] = '"none"'

# set output to netcdf unless user has specified something else
if ('format' not in parameters):
  parameters['format'] = '"netcdf"'

# set file prefix to modified model name unless user has specified something else
if ('prefix' not in parameters):
  parameters['prefix'] = '"' + model_file_name.replace('.in','') + '"'

# match dt_output to max of dt or 60s unless it is specifically provided
if ('dt' in parameters and 'dt_output' not in parameters):
  parameters['dt_output'] = str(max(float(parameters['dt']), 60.0))

# if dz specified, set other associated options
if ('dz' in parameters):
  parameters['grid_type'] = '1'
  parameters['zt_filename'] = "''"

# if levels specified, expand number into file name
if ('levels' in parameters):
  parameters['zt_filename'] = "'../input/grid/deep_convection_{}lev_27km_zt_grid.grd'".format(parameters['levels'])

# if refinement specified, create grid file and set appropriate file name
if ('refine' in parameters):
  height = 10000.0
  refine = int(parameters['refine'])
  grid128 = np.loadtxt('../input/grid/deep_convection_128lev_27km_zt_grid.grd')
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
  parameters['levels'] = str(len(grid))
  filename = 'deep_convection_{}lev_27km_zt_grid.grd'.format(parameters['levels'])
  np.savetxt(filename, grid)
  parameters['zt_filename'] = "'{}'".format(filename)

# set l_standard_term_ta to true unless user has specified something else
if ('aterms' not in parameters):
  parameters['aterms'] = 'unmodified'

# set time_final to 3 hours unless user has specified something else
if ('tfinal' not in parameters):
  parameters['tfinal'] = str(3*3600.0)

# set centered difference unless user has specified something else
if ('godunov' not in parameters):
  parameters['godunov'] = 'none'

# set no splatting unless user has specified something else
if ('splat' not in parameters):
  parameters['splat'] = '0.0'

# set new boundary conditions
parameters['newBC'] = '.true.'
parameters['scaleBC'] = '5.0'

# warn about Tsfc parameter
if ('Tsfc' in parameters):
  print("WARNING: Specifying Tsfc doesn't change model parammeters...")
  print("         this needs to be done in " + sys.argv[1] + "_sounding.in")
  input("Press Enter to acknowledge")

modified_lines = []
for line in model_config:
  modified = False
  if (not line.strip().startswith('!')):
    for parameter in parameters:
      # adjust values
      if (parameter == 'dt' and (line.startswith('dt_main') or line.startswith('dt_rad'))
        or parameter == 'dt_output' and (line.startswith('stats_tsamp') or line.startswith('stats_tout'))
        or parameter == 'microphysics' and line.startswith('microphys_scheme')
        or parameter == 'prefix' and line.startswith('fname_prefix')
        or parameter == 'format' and line.startswith('stats_fmt')
        or parameter == 'grid_type' and line.startswith('grid_type')
        or parameter == 'dz' and line.startswith('deltaz')
        or parameter == 'zmax' and line.startswith('zm_top')
        or parameter == 'zt_filename' and line.startswith('zt_grid_fname')
        or parameter == 'levels' and line.startswith('nzmax')
        or parameter == 'tfinal' and line.startswith('time_final')
        or parameter == 'newBC' and line.startswith('l_fixed_level_for_surflx')
        or parameter == 'scaleBC' and line.startswith('f_scale_surflx')):
        default_value = line.split()[2]
        line = line.replace(default_value, parameters[parameter])
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
input_file = open('../input/tunable_parameters/tunable_parameters.in', 'r')
input_lines = input_file.readlines()
input_file.close()
line_collection = []
for line in input_lines:
  if (line.startswith('C_wp2_splat')):
    default_value = line.split()[2]
    line = line.replace(default_value, parameters['splat'])
    print('Setting splat coefficient to {}'.format(parameters['splat']))
    print(line)
  line_collection.append(line)
line_collections.append(line_collection)
# "SILHS_PARAMS" file
input_file = open('../input/tunable_parameters/silhs_parameters.in', 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "FLAGS" file
input_file = open('../input/tunable_parameters/configurable_model_flags.in', 'r')
input_lines = input_file.readlines()
input_file.close()
line_collection = []
for line in input_lines:
  if (parameters['godunov'] == 'scalarwp3' and line.startswith('l_upwind_wp3_ta')):
    line = line.replace('false','true')
    print('Setting l_upwind_wp3_ta flag to true')
    print(line)
  elif (parameters['aterms'] == 'unmodified' and line.startswith('l_standard_term_ta')):
    line = line.replace('false','true')
    print('Setting l_standard_term_ta flag to true')
    print(line)
  line_collection.append(line)
line_collections.append(line_collection)
# "MOD_MODEL" file
input_file = open(model_file_name, 'r')
line_collections.append(input_file.readlines())
input_file.close()
# "stats" file
input_file = open('../input/stats/standard_stats.in', 'r')
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
input('Press Enter to run CLUBB...')
os.system('../bin/clubb_standalone')
