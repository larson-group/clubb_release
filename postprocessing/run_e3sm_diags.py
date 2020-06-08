# source /lcrc/soft/climate/e3sm-unified/base/etc/profile.d/conda.sh
# conda activate e3sm_diags_env
# python run_e3sm_diags.py
import os
from acme_diags.parameter.core_parameter import CoreParameter
from acme_diags.run import runner

param = CoreParameter()

param.reference_data_path = '/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
param.test_data_path = '/lcrc/group/acme/ac.griffin/E3SM_simulations/default.alpha22.F2010.clubb_c_K10h_0p35.ne30pg2_r05_oECv3.anvil/run/climatology_yrs_1_7/regridded_climo/'
param.test_name = 'default.alpha22.F2010.clubb_c_K10h_0p35.ne30pg2_r05_oECv3.anvil'
param.seasons = ["ANN"]   #all seasons ["ANN","DJF", "MAM", "JJA", "SON"] will run,if comment out"

prefix = '/lcrc/group/acme/public_html/diagnostic_output/griffinb/'
param.results_dir = os.path.join(prefix, 'clubb_c_K10h_0p35')
# Use the following if running in parallel:
#param.multiprocessing = True
#param.num_workers = 32

# Use below to run all core sets of diags:
#runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
# Use below to run lat_lon map only:
runner.sets_to_run = ['lat_lon']
runner.run_diags([param])
