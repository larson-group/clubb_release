"""
:author: Nicolas Strike
:date: Early 2019

This file is mostly a definition of Cases. Each case is defined in the following format
using python dictionaries (values surrounded with < > must have the < > removed to be valid).

.. code-block:: python
    :linenos:

    CASENAME = {'name': 'casename', 'start_time': <numeric value>, 'end_time': <numeric value>,
       'height_min_value': <numeric value>, 'height_max_value': <numeric value>,
       
       'blacklisted_vars': ['list', 'of', 'variable', 'names', 'to', 'exclude', 'from', 'plotting'],
       'sam_file': <path to sam file>",
       'coamps_file': <path to coamps file>,
       'r408_file': {'zm': <path to r408 file>,
                     'zt': <path to r408 file>,
                     'sfc': <path to r408 file>},
       'var_groups': [VariableGroupBase, <other variable groups to plot>]}

**Important note**:
When creating a new case, add it to the ALL_CASES list at the bottom of the file. Additionally, please add it in
alphabetical order.

**Case Definition values explained**:

    *name*: must be the same as the filename without the extention.
        E.g. to use lba_zt.nc and lba_zm.nc the case's name must be 'lba'. Extensions are determined
        by the last instance of _

    *start_time*: An integer value representing which timestep to begin the time-averaging interval.
        Valid options are from 1 -> list minute value. Give in terms of clubb minutes.

    *end_time*: An integer value representing which timestep to end the time-averaging interval.
        Valid options are from 1 -> list minute value. Give in terms of clubb minutes.
        Also used to determine where to stop timeseries plots

    *height_min_value*: The elevation to begin height plots at

    *height_max_value*: The elevation to end height plots at

    *blacklisted_vars*: List of variables to avoid plotting for this case. Names must use the clubb-name version

    *sam_file*: Path to the SAM .nc file for this case (please use the SAM_OUTPUT_ROOT variable as a base).

    *coamps_file*: dict containing paths to the coamps .nc files for this case (please use the LES_OUTPUT_ROOT variable as a base).

    *r408_file*: dict containing paths to the r408 .nc files for this case (please use the R408_OUTPUT_ROOT variable as a base).

    *var_groups*: list of python class names, where the classes use the naming scheme VariableGroup____.py and define a variable group.

"""

import os

from config.VariableGroupBase import VariableGroupBase
from config.VariableGroupCorrelations import VariableGroupCorrelations
from config.VariableGroupIceMP import VariableGroupIceMP
from config.VariableGroupKKMP import VariableGroupKKMP
from config.VariableGroupLiquidMP import VariableGroupLiquidMP
from config.VariableGroupWs import VariableGroupWs

# ---------------------------
# DO NOT MODIFY THESE PATHS
BENCHMARK_OUTPUT_ROOT = os.path.dirname(os.path.realpath(__file__)) + "/../les_and_clubb_benchmark_runs/"
SAM_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "sam_benchmark_runs"
LES_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "les_runs"
ARCHIVED_CLUBB_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "archived_clubb_runs"
R408_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + ""
# ---------------------------

"""
To plot only a subset of cases, reguardless of what output exists
in the input folder, uncomment the last line of this file and
fill that array with the cases you'd like to plot. This overwrites the
ALL_CASES variable such that pyplotgen will only know about cases in that
list and ignore all others. The name must match the python variable name
below (all caps).

For example, to plot only bomex and fire:

# ALL_CASES = [BOMEX, FIRE]
"""

ARM = {'name': 'arm', 'start_time': 481, 'end_time': 540, 'height_min_value': 0, 'height_max_value': 3500,
       
       'blacklisted_vars': ['radht'],
       'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/ARM_96x96x110/GCSSARM_96x96x110_67m_40m_1s.nc",
       'coamps_file': None,
       'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_zm.nc',
                     'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_zt.nc',
                     'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_sfc.nc'},
       'var_groups': [VariableGroupBase]}

ARM_97 = {'name': 'arm_97', 'start_time': 4321, 'end_time': 5580, 'height_min_value': 0, 'height_max_value': 18000,
          
          'blacklisted_vars': ['rtp3','Skrt_zt', 'Skthl_zt', 'thlp3', 'rtpthvp', 'thlpthvp'],
          'sam_file': SAM_OUTPUT_ROOT + "/ARM97_r1315_128x128x128_1km_Morrison/ARM9707.nc",
          'coamps_file': None,
          'r408_file': None,
          'var_groups': [VariableGroupBase, VariableGroupIceMP]}

ASTEX_A209 = {'name': 'astex_a209', 'start_time': 2340, 'end_time': 2400, 'height_min_value': 0,
              'height_max_value': 6000, 
              'blacklisted_vars': [],
              'sam_file': None,
              'coamps_file': None,
              'r408_file': None,
              'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                             VariableGroupKKMP]}

ATEX = {'name': 'atex', 'start_time': 421, 'end_time': 480, 'height_min_value': 0, 'height_max_value': 2500,
        
        'blacklisted_vars': [],
        'sam_file': None,
        'coamps_file': {'sm': LES_OUTPUT_ROOT + "/atex_coamps_sm.nc",
                        'sw': LES_OUTPUT_ROOT + "/atex_coamps_sw.nc"},
        'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_zm.nc',
                      'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_zt.nc',
                      'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_sfc.nc'},
        'var_groups': [VariableGroupBase, VariableGroupWs]}

BOMEX = {'name': 'bomex', 'start_time': 181, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 2500,
         
         'blacklisted_vars': [],
         'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc",
         'coamps_file': None,
         'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zm.nc',
                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zt.nc',
                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_sfc.nc'},
         'var_groups': [VariableGroupBase, VariableGroupWs]}

CGILS_S6 = {'name': 'cgils_s6', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0, 'height_max_value': 5950,
            
            'blacklisted_vars': ['Ngm', 'rgm', 'Skrt_zt', 'Skthl_zt', 'thlp3', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
            'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s6/ctl_s6_96x96x128_100m_DRZ_N100_tqndg.nc",
            'coamps_file': None,
            'r408_file': None,
            'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S11 = {'name': 'cgils_s11', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0,
             'height_max_value': 5950, 
             'blacklisted_vars': ['Ngm', 'rgm', 'Skthl_zt', 'Skrt_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s11/ctl_s11_96x96x320_50m_DRZ_N100_ref.nc",
             'coamps_file': None,
             'r408_file': None,
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S12 = {'name': 'cgils_s12', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0,
             'height_max_value': 5950, 
             'blacklisted_vars': ['Ngm', 'rgm', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s12/ctl_s12_96x96x192_25m_DRZ_N100_fixnudge.nc",
             'coamps_file': None,
             'r408_file': None,
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_NOV02 = {'name': 'clex9_nov02', 'start_time': 181, 'end_time': 240, 'height_min_value': 3072,
               'height_max_value': 6072, 
               'blacklisted_vars': [],
               'sam_file': None,
               'coamps_file': {'sm': LES_OUTPUT_ROOT + "/clex9_nov02_coamps_sm.nc",
                               'sw': LES_OUTPUT_ROOT + "/clex9_nov02_coamps_sw.nc"},
               'r408_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_OCT14 = {'name': 'clex9_oct14', 'start_time': 181, 'end_time': 240, 'height_min_value': 2188,
               'height_max_value': 6688, 
               'blacklisted_vars': [],
               'sam_file': None,
               'coamps_file': {'sm': LES_OUTPUT_ROOT + "/clex9_oct14_coamps_sm.nc",
                               'sw': LES_OUTPUT_ROOT + "/clex9_oct14_coamps_sw.nc"},
               'r408_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

DYCOMS2_RF01 = {'name': 'dycoms2_rf01', 'start_time': 181, 'end_time': 240, 'height_min_value': 0,
                'height_max_value': 1200, 
                'blacklisted_vars': [],
                'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF01_96x96x320/DYCOMS_RF01_96x96x320.nc",
                'coamps_file': None,
                'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_zm.nc',
                              'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_zt.nc',
                              'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_sfc.nc'},
                'var_groups': [VariableGroupBase, VariableGroupWs]}

DYCOMS2_RF01_FIXED_SST = {'name': 'dycoms2_rf01_fixed_sst', 'start_time': 2520, 'end_time': 2700, 'height_min_value': 0,
                          'height_max_value': 1200,
                          
                          'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp'],
                          'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/DYCOMS_RF01_fixed_sst/DYCOMS_RF01_96x96x320_LES_fixed_sst.nc",
                          'coamps_file': None,
                          'r408_file': None,
                          'var_groups': [VariableGroupBase]}

DYCOMS2_RF02_DO = {'name': 'dycoms2_rf02_do', 'start_time': 301, 'end_time': 360, 'height_min_value': 0,
                   'height_max_value': 1200, 
                   'blacklisted_vars': [],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_dr_nosed/DYCOMS_RF02_128x128x96_dr_nosed.nc",
                   'coamps_file': None,
                   'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_zm.nc',
                                 'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_zt.nc',
                                 'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_sfc.nc'},
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                                  VariableGroupKKMP]}

DYCOMS2_RF02_DS = {'name': 'dycoms2_rf02_ds', 'start_time': 301, 'end_time': 360, 'height_min_value': 0,
                   'height_max_value': 1200, 
                   'blacklisted_vars': [],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_dr_sed/DYCOMS_RF02_128x128x96_dr_sed.nc",
                   'coamps_file': None,
                   'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_zm.nc',
                                 'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_zt.nc',
                                 'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_sfc.nc'},
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                                  VariableGroupKKMP]}

DYCOMS2_RF02_ND = {'name': 'dycoms2_rf02_nd', 'start_time': 301, 'end_time': 360, 'height_min_value': 0,
                   'height_max_value': 1200, 
                   'blacklisted_vars': ['wprrp', 'wpNrp', 'corr_w_rr_1', 'corr_w_Nr_1'],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_nosed/DYCOMS_RF02_128x128x96_nodr_nosed.nc",
                   'coamps_file': None,
                   'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_zm.nc',
                                 'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_zt.nc',
                                 'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_sfc.nc'},
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

DYCOMS2_RF02_SO = {'name': 'dycoms2_rf02_so', 'start_time': 301, 'end_time': 360, 'height_min_value': 0,
                   'height_max_value': 1200, 
                   'blacklisted_vars': ['wprrp', 'wpNrp'],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_sed/DYCOMS_RF02_128x128x96_nodr_sed.nc",
                   'coamps_file': None,
                   'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_zm.nc',
                                 'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_zt.nc',
                                 'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_sfc.nc'},
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

FIRE = {'name': 'fire', 'start_time': 61, 'end_time': 120, 'height_min_value': 0, 'height_max_value': 1000,
        
        'blacklisted_vars': [],
        'sam_file': None,
        'coamps_file': { 'sm' : LES_OUTPUT_ROOT + "/fire_coamps_sm.nc",
                         'sw' : LES_OUTPUT_ROOT + "/fire_coamps_sw.nc"},
        'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_zm.nc',
                      'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_zt.nc',
                      'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_sfc.nc'},
        'var_groups': [VariableGroupBase, VariableGroupWs]}

# No budgets
GABLS2 = {'name': 'gabls2', 'start_time': 2101, 'end_time': 2160, 'height_min_value': 0, 'height_max_value': 2500,
          
          'blacklisted_vars': ['tau_zm', 'radht', 'Skw_zt', 'Skrt_zt', 'Skthl_zt', 'corr_w_chi_1', 'corr_chi_eta_1',
                               'rcp2', 'thlpthvp', 'rtpthvp'],
          'sam_file': None,
          'coamps_file': {'sm': LES_OUTPUT_ROOT + "/gabls2_coamps_sm.nc",
                          'sw': LES_OUTPUT_ROOT + "/gabls2_coamps_sw.nc"},
          'r408_file': None,
          'var_groups': [VariableGroupBase]}

GABLS3 = {'name': 'gabls3', 'start_time': 1081, 'end_time': 1200, 'height_min_value': 0, 'height_max_value': 4970,
          
          'blacklisted_vars': [],
          'sam_file': None,
          'coamps_file': None,
          'r408_file': None,
          'var_groups': [VariableGroupBase]}

GABLS3_NIGHT = {'name': 'gabls3_night', 'start_time': 421, 'end_time': 480, 'height_min_value': 0,
                'height_max_value': 800, 
                'blacklisted_vars': [],
                'sam_file': None,
                'coamps_file': None,
                'r408_file': None,
                'var_groups': [VariableGroupBase]}

JUN25_ALTOCU = {'name': 'jun25_altocu', 'start_time': 181, 'end_time': 240, 'height_min_value': 4808,
                'height_max_value': 7308, 
                'blacklisted_vars': ['Ngm', 'wprrp', 'wpNrp'],
                'sam_file': None,
                'coamps_file': {'sm': LES_OUTPUT_ROOT + "/jun25_altocu_qc3_coamps_sm.nc",
                                'sw': LES_OUTPUT_ROOT + "/jun25_altocu_qc3_coamps_sw.nc"},
                'r408_file': None,
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

LBA = {'name': 'lba', 'start_time': 300, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 12000,
       
       'blacklisted_vars': ['wprrp', 'wpNrp'],
       'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/LBA_128kmx128kmx128_1km_Morrison/LBA_128kmx128kmx128_1km_Morrison.nc",
       'coamps_file': None,
       'r408_file': None,
       'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP, VariableGroupWs]}

MC3E = {'name': 'mc3e', 'start_time': 1, 'end_time': 64800, 'height_min_value': 0, 'height_max_value': 18000,
        
        'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'Ngm', 'wprrp', 'wpNrp'],
        'sam_file': SAM_OUTPUT_ROOT + "/MC3E_r1359_128x128x128_1km_Morrison/MC3E.nc",
        'coamps_file': None,
        'r408_file': None,
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_A = {'name': 'mpace_a', 'start_time': 4141, 'end_time': 4320, 'height_min_value': 0, 'height_max_value': 10000,
           
           'blacklisted_vars': ['Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'Ngm', 'wpNrp'],
           'sam_file': None,
           'coamps_file': None,
           'r408_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B = {'name': 'mpace_b', 'start_time': 541, 'end_time': 720, 'height_min_value': 0, 'height_max_value': 2750,
           
           'blacklisted_vars': ['Ngm', 'wpNrp'],
           'sam_file': None,
           'coamps_file': {'sm': LES_OUTPUT_ROOT + "/mpace_b_coamps_sm.nc",
                           'sw': LES_OUTPUT_ROOT + "/mpace_b_coamps_sw.nc"},
           'r408_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B_SILHS = {'name': 'mpace_b_silhs', 'start_time': 541, 'end_time': 720, 'height_min_value': 0,
                 'height_max_value': 2750, 
                 'blacklisted_vars': ['Ngm', 'wpNrp'],
                 'sam_file': None,
                 'coamps_file': {'sm': LES_OUTPUT_ROOT + "/mpace_b_coamps_sm.nc",
                                 'sw': LES_OUTPUT_ROOT + "/mpace_b_coamps_sw.nc"},
                 'r408_file': None,
                 'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

NOV11_ALTOCU = {'name': 'nov11_altocu', 'start_time': 91, 'end_time': 150, 'height_min_value': 4150,
                'height_max_value': 6150, 
                'blacklisted_vars': ['Ngm'],
                'sam_file': None,
                'coamps_file':{'sm': LES_OUTPUT_ROOT + "/nov11_altocu_coamps_sm.nc",
                               'sw': LES_OUTPUT_ROOT + "/nov11_altocu_coamps_sw.nc"},
                'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_zm.nc',
                              'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_zt.nc',
                              'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_sfc.nc'},
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

RICO = {'name': 'rico', 'start_time': 4201, 'end_time': 4320, 'height_min_value': 0, 'height_max_value': 4000,
        
        'blacklisted_vars': ['wpNrp'],
        'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/RICO_256x256x100_drizzle/RICO_256x256x100_drizzle.nc",
        'coamps_file': None,
        'r408_file': None,
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupWs, VariableGroupCorrelations,
                       VariableGroupKKMP]}

TWP_ICE = {'name': 'twp_ice', 'start_time': 1, 'end_time': 9900, 'height_min_value': 0, 'height_max_value': 19000,
           
           'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'Ngm', 'wprrp', 'wpNrp'],
           'sam_file': SAM_OUTPUT_ROOT + "/TWP_ICE_r1315_128x128x128_1km_Morrison/TWP_ICE.nc",
           'coamps_file': None,
           'r408_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

WANGARA = {'name': 'wangara', 'start_time': 181, 'end_time': 240, 'height_min_value': 0, 'height_max_value': 1900,
           
           'blacklisted_vars': [],
           'sam_file': None,
           'coamps_file': None, # uses RAMS-LES /LES_files/wangara_rams.nc
           'r408_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_zm.nc',
                         'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_zt.nc',
                         'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_sfc.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX,
             BOMEX,
             CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14,
             DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST, DYCOMS2_RF02_DO, DYCOMS2_RF02_DS, DYCOMS2_RF02_ND, DYCOMS2_RF02_SO,
             FIRE,
             GABLS2, GABLS3, GABLS3_NIGHT,
             JUN25_ALTOCU,
             LBA,
             MC3E, MPACE_A, MPACE_B, MPACE_B_SILHS,
             NOV11_ALTOCU,
             RICO,
             TWP_ICE,
             WANGARA
             ]

# ALL_CASES = [JUN25_ALTOCU]