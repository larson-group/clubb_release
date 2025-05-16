"""
:author: Nicolas Strike
:date: Early 2019

This file is mostly a definition of Cases. Each case is defined in the following format
using python dictionaries (values surrounded with < > must have the < > removed to be valid).

.. code-block:: python
    :linenos:

    CASENAME = {'name': 'casename',
        'description': "",
        'start_time': <numeric value>, 'end_time': <numeric value>,
        'height_min_value': <numeric value>, 'height_max_value': <numeric value>,
        'blacklisted_vars': ['list', 'of', 'variable', 'names', 'to', 'exclude', 'from', 'plotting'],
        'sam_benchmark_file': <path to sam file>",
        'clubb_file': {'zm': <path to file>,
                      'zt': <path to file>,
                      'sfc': <path to file>},
        'coamps_benchmark_file': {'sm': <path to file>,
                          'sw': <path to file>},
        'clubb_r408_benchmark_file': {'zm': <path to file>,
                      'zt': <path to file>,
                      'sfc': <path to file>},
       'clubb_hoc_benchmark_file': {'zm': <path to file>',
                    'zt': <path to file>',
                    'sfc': <path to file>},
       'e3sm_file': <path to file>,
       'cam_file': <path to file>,
       'sam_file': <path to file>,
       'wrf_file': {'zm': <path to file>,
                    'zt': <path to file>,
                    'sfc': <path to file>},
        'var_groups': [VariableGroupBase, <other variable groups to plot>]}

**Important note**:
When creating a new case, add it to the CASES_TO_PLOT list at the bottom of the file. Additionally, please add it in
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

    *<model name>_file*: The path(s) to nc files for the given model.
        (please use the <model name>_OUTPUT_ROOT variables as the beginning of the path).

    *var_groups*: These are the groups of variables to be plotted for the given case. var_groups is defined as a
        list of python class names, where the classes use the naming scheme VariableGroup____.py and define a variable
        group. An example would be: 'var_groups': [VariableGroupBase, VariableGroupWs].
        The variables inside a VariableGroup can be found in the file with the same name,
        i.e. config/VariableGroupBase.py. An example would be thlm in VariableGroupBase.

"""

import os

from config.VariableGroupBase import VariableGroupBase
from config.VariableGroupCorrelations import VariableGroupCorrelations
from config.VariableGroupIceMP import VariableGroupIceMP
from config.VariableGroupKKMP import VariableGroupKKMP
from config.VariableGroupLiquidMP import VariableGroupLiquidMP
from config.VariableGroupSamProfiles import VariableGroupSamProfiles
from config.VariableGroupScalars import VariableGroupScalars
from config.VariableGroupWs import VariableGroupWs
from config.VariableGroupTaus import VariableGroupTaus
from config.VariableGroupNondimMoments import VariableGroupNondimMoments
from config.VariableGroupNormalizedVariations import VariableGroupNormalizedVariations

# ---------------------------
BENCHMARK_OUTPUT_ROOT = "/home/pub/les_and_clubb_benchmark_runs/"
if not os.path.isdir(BENCHMARK_OUTPUT_ROOT) and \
        not os.path.islink(BENCHMARK_OUTPUT_ROOT):
    print("Benchmark output was not found in " + BENCHMARK_OUTPUT_ROOT + ".\n\tChecking local location: " +
          os.path.dirname(os.path.realpath(__file__)) + "/../les_and_clubb_benchmark_runs/")
    BENCHMARK_OUTPUT_ROOT = os.path.dirname(os.path.realpath(__file__)) + "/../les_and_clubb_benchmark_runs/"
SAM_BENCHMARK_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "sam_benchmark_runs"
COAMPS_BENCHMARK_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "les_runs"
WRF_LASSO_BENCHMARK_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "wrf_lasso_runs"
ARCHIVED_CLUBB_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "archived_clubb_runs"
R408_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + ""
HOC_OUTPUT_ROOT = BENCHMARK_OUTPUT_ROOT + "HOC_20051217"

# This folder is passed in as a command line parameter
# It is not capitalized because it is not intended to
# be final, i.e. is changed depending on the cmd line arg
e3sm_output_root = ""
sam_output_root = ""
wrf_output_root = ""
cam_output_root = ""
clubb_output_root = ""
# ---------------------------

# These are all the names that represent the height variable within different models
HEIGHT_VAR_NAMES = ['z', 'Z3', 'altitude', 'lev', 'CSP_Zm', 'CSP_Z8Wm'] # CSP_* added for WRF-LASSO cases
TIME_VAR_NAMES = ['time', 'XTIME']
"""
To plot only a subset of cases, reguardless of what output exists
in the clubb folder, uncomment the last line of this file and
fill that array with the cases you'd like to plot. This overwrites the
CASES_TO_PLOT variable such that pyplotgen will only know about cases in that
list and ignore all others. The name must match the python variable name
below (all caps).

For example, to plot only bomex and fire:

CASES_TO_PLOT = [BOMEX, FIRE]
"""

ARM = {'name': 'arm',
       'description': "Output may differ from plotgen in some models (e.g. WRF) due to a difference in the time "
                      "averaging interval.",
       'start_time': 1, 'end_time': 870,
       #'height_min_value': 0, 'height_max_value': 3000,
       'height_min_value': 0, 'height_max_value': 8500,

       'blacklisted_vars': ['radht'],
       'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT + "/JULY_2017/ARM_96x96x110/GCSSARM_96x96x110_67m_40m_1s.nc"},
       'clubb_file': {'zm': clubb_output_root + '/arm_zm.nc',
                      'zt': clubb_output_root + '/arm_zt.nc',
                      'sfc': clubb_output_root + '/arm_sfc.nc'},
       'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/arm_coamps_sm.nc",
                                 'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/arm_coamps_sw.nc"},
       'wrf_benchmark_file': None,
       'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_zm.nc',
                           'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_zt.nc',
                           'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/arm_sfc.nc'},
       'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/arm_zm.nc',
                          'zt': HOC_OUTPUT_ROOT + '/arm_zt.nc',
                          'sfc': HOC_OUTPUT_ROOT + '/arm_sfc.nc'},
       'e3sm_file': { 'e3sm': e3sm_output_root + "/arm.nc"},
       'cam_file': None,
       'sam_file': {'sam': sam_output_root + "/GCSSARM_96x96x110_67m_40m_1s.nc"},
       'wrf_file': {'zm': wrf_output_root + "/arm_zm_wrf.nc",
                    'zt': wrf_output_root + "/arm_zt_wrf.nc",
                    'sfc': wrf_output_root + "/arm_sfc_wrf.nc"
                    },
       'var_groups': [VariableGroupBase, VariableGroupWs]}

ARM_97 = {'name': 'arm_97',
          'description': "",
          'start_time': 4321, 'end_time': 5580,
          'height_min_value': 0, 'height_max_value': 18000,

          'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'thlp3', 'rtpthvp', 'thlpthvp'],
          'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                  "/ARM97_r1315_128x128x128_1km_Morrison/ARM9707.nc"},
          'clubb_file': {'zm': clubb_output_root + '/arm_97_zm.nc',
                         'zt': clubb_output_root + '/arm_97_zt.nc',
                         'sfc': clubb_output_root + '/arm_97_sfc.nc',
                         'subcolumns': clubb_output_root + '/arm_97_nl_lh_sample_points_2D.nc'},
          'coamps_benchmark_file': None,
          'wrf_benchmark_file': None,
          'clubb_r408_benchmark_file': None,
          'clubb_hoc_benchmark_file': None,
          'e3sm_file': None,
          'cam_file': None,
          'sam_file': {'sam': sam_output_root + "/ARM9707_SAM_CLUBB.nc"},
          'wrf_file': None,
          'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupIceMP]}

ASTEX_A209 = {'name': 'astex_a209',
              'description': "",
              #'start_time': 2340, 'end_time': 2400,
              'start_time': 1, 'end_time': 2500,
              'height_min_value': 0, 'height_max_value': 8500,
              #'height_min_value': 0, 'height_max_value': 3000,
              #'height_min_value': 0, 'height_max_value': 2500,

              'blacklisted_vars': [],
              'sam_benchmark_file': None,
              'clubb_file': {'zm': clubb_output_root + '/astex_a209_zm.nc',
                             'zt': clubb_output_root + '/astex_a209_zt.nc',
                             'sfc': clubb_output_root + '/astex_a209_sfc.nc'},
              'coamps_benchmark_file': None,
              'wrf_benchmark_file': None,
              'clubb_r408_benchmark_file': None,
              'clubb_hoc_benchmark_file': None,
              'e3sm_file': None,
              'cam_file': None,
              'sam_file': None,
              'wrf_file': None,
              'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                             VariableGroupKKMP]}

ATEX = {'name': 'atex',
        'description': "",
        'start_time': 421, 'end_time': 480,
        'height_min_value': 0, 'height_max_value': 2500,

        'blacklisted_vars': [],
        'sam_benchmark_file': None,
        'clubb_file': {'zm': clubb_output_root + '/atex_zm.nc',
                       'zt': clubb_output_root + '/atex_zt.nc',
                       'sfc': clubb_output_root + '/atex_sfc.nc'},
        'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/atex_coamps_sm.nc",
                                  'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/atex_coamps_sw.nc"},
        'wrf_benchmark_file': None,
        'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_zm.nc',
                            'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_zt.nc',
                            'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/atex_sfc.nc'},
        'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/atex_zm.nc',
                           'zt': HOC_OUTPUT_ROOT + '/atex_zt.nc',
                           'sfc': HOC_OUTPUT_ROOT + '/atex_sfc.nc'},
        'e3sm_file': None,
        'cam_file': {'cam': cam_output_root + "/atex_cam.nc"},
        'sam_file': None,
        'wrf_file': {'zm': wrf_output_root + "/atex_zm_wrf.nc",
                     'zt': wrf_output_root + "/atex_zt_wrf.nc",
                     'sfc': wrf_output_root + "/atex_sfc_wrf.nc"
                     },
        'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupIceMP]}

BOMEX = {'name': 'bomex',
         'description': "",
         'start_time': 181, 'end_time': 360,
         'height_min_value': 0, 'height_max_value': 2500,

         'blacklisted_vars': [],
         'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                 "/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc"},
         'clubb_file': {'zm': clubb_output_root + '/bomex_zm.nc',
                        'zt': clubb_output_root + '/bomex_zt.nc',
                        'sfc': clubb_output_root + '/bomex_sfc.nc'},
         'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/bomex_coamps_sm.nc",
                                   'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/bomex_coamps_sw.nc"},
         'wrf_benchmark_file': None,
         'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zm.nc',
                             'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_zt.nc',
                             'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/bomex_sfc.nc'},
         'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/bomex_zm.nc',
                            'zt': HOC_OUTPUT_ROOT + '/bomex_zt.nc',
                            'sfc': HOC_OUTPUT_ROOT + '/bomex_sfc.nc'},
         'e3sm_file': { 'e3sm': e3sm_output_root + '/bomex.nc'},
         'cam_file': None,
         'sam_file': {'sam': sam_output_root + "/BOMEX_SAM_CLUBB.nc"},
         'wrf_file': {'zm': wrf_output_root + '/bomex_zm_wrf.nc',
                      'zt': wrf_output_root + '/bomex_zt_wrf.nc',
                      'sfc': wrf_output_root + '/bomex_sfc_wrf.nc'},
         'var_groups': [VariableGroupBase, VariableGroupWs]}

CGILS_S6 = {'name': 'cgils_s6',
            'description': "",
            'start_time': 12960, 'end_time': 14400,
            'height_min_value': 0, 'height_max_value': 5950,

            'blacklisted_vars': ['Ngm', 'rgm', 'Skrt_zt', 'Skthl_zt', 'thlp3',
                                 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
            'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                  "/SAM6.6/CLOUD_FEEDBACK_s6/ctl_s6_96x96x128_100m_DRZ_N100_tqndg.nc"},
            'clubb_file': {'zm': clubb_output_root + '/cgils_s6_zm.nc',
                           'zt': clubb_output_root + '/cgils_s6_zt.nc',
                           'sfc': clubb_output_root + '/cgils_s6_sfc.nc'},
            'coamps_benchmark_file': None,
            'wrf_benchmark_file': None,
            'clubb_r408_benchmark_file': None,
            'clubb_hoc_benchmark_file': None,
            'e3sm_file': None,
            'cam_file': None,
            'sam_file': None,
            'wrf_file': None,
            'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S11 = {'name': 'cgils_s11',
             'description': "",
             'start_time': 12960, 'end_time': 14400,
             'height_min_value': 0, 'height_max_value': 5950,

             'blacklisted_vars': ['Ngm', 'rgm', 'Skthl_zt', 'Skrt_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                   "/SAM6.6/CLOUD_FEEDBACK_s11/ctl_s11_96x96x320_50m_DRZ_N100_ref.nc"},
             'clubb_file': {'zm': clubb_output_root + '/cgils_s11_zm.nc',
                            'zt': clubb_output_root + '/cgils_s11_zt.nc',
                            'sfc': clubb_output_root + '/cgils_s11_sfc.nc'},
             'coamps_benchmark_file': None,
             'wrf_benchmark_file': None,
             'clubb_r408_benchmark_file': None,
             'clubb_hoc_benchmark_file': None,
             'e3sm_file': None,
             'cam_file': None,
             'sam_file': None,
             'wrf_file': None,
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S12 = {'name': 'cgils_s12',
             'description': "",
             'start_time': 12960, 'end_time': 14400,
             'height_min_value': 0, 'height_max_value': 5950,

             'blacklisted_vars': ['Ngm', 'rgm', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                   "/SAM6.6/CLOUD_FEEDBACK_s12/ctl_s12_96x96x192_25m_DRZ_N100_fixnudge.nc"},
             'clubb_file': {'zm': clubb_output_root + '/cgils_s12_zm.nc',
                            'zt': clubb_output_root + '/cgils_s12_zt.nc',
                            'sfc': clubb_output_root + '/cgils_s12_sfc.nc'},
             'coamps_benchmark_file': None,
             'wrf_benchmark_file': None,
             'clubb_r408_benchmark_file': None,
             'clubb_hoc_benchmark_file': None,
             'e3sm_file': None,
             'cam_file': None,
             'sam_file': None,
             'wrf_file': None,
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_NOV02 = {'name': 'clex9_nov02',
               'description': "",
               'start_time': 181, 'end_time': 240,
               'height_min_value': 4000, 'height_max_value': 6072,

               'blacklisted_vars': ['Ngm'],
               'sam_benchmark_file': None,
               'clubb_file': {'zm': clubb_output_root + '/clex9_nov02_zm.nc',
                              'zt': clubb_output_root + '/clex9_nov02_zt.nc',
                              'sfc': clubb_output_root + '/clex9_nov02_sfc.nc'},
               'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/clex9_nov02_coamps_sm.nc",
                                         'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/clex9_nov02_coamps_sw.nc"},
               'wrf_benchmark_file': None,
               'clubb_r408_benchmark_file': None,
               'clubb_hoc_benchmark_file': None,
               'e3sm_file': None,
               'cam_file': None,
               'sam_file': None,
               'wrf_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_OCT14 = {'name': 'clex9_oct14',
               'description': "",
               'start_time': 181, 'end_time': 240,
               'height_min_value': 2230, 'height_max_value': 6688,

               'blacklisted_vars': ['Ngm'],
               'sam_benchmark_file': None,
               'clubb_file': {'zm': clubb_output_root + '/clex9_oct14_zm.nc',
                              'zt': clubb_output_root + '/clex9_oct14_zt.nc',
                              'sfc': clubb_output_root + '/clex9_oct14_sfc.nc'},
               'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/clex9_oct14_coamps_sm.nc",
                                         'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/clex9_oct14_coamps_sw.nc"},
               'wrf_benchmark_file': None,
               'clubb_r408_benchmark_file': None,
               'clubb_hoc_benchmark_file': None,
               'e3sm_file': None,
               'cam_file': None,
               'sam_file': None,
               'wrf_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

DYCOMS2_RF01 = {'name': 'dycoms2_rf01',
                'description': "",
                'start_time': 181, 'end_time': 240,
                'height_min_value': 0, 'height_max_value': 1200,

                'blacklisted_vars': [],
                'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                        "/JULY_2017/DYCOMS_RF01_96x96x320/DYCOMS_RF01_96x96x320.nc"},
                'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf01_zm.nc',
                               'zt': clubb_output_root + '/dycoms2_rf01_zt.nc',
                               'sfc': clubb_output_root + '/dycoms2_rf01_sfc.nc'},
                'coamps_benchmark_file': None,
                'wrf_benchmark_file': None,
                'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_zm.nc',
                                    'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_zt.nc',
                                    'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf01_sfc.nc'},
                'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf01_zm.nc',
                                   'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf01_zt.nc',
                                   'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf01_sfc.nc'},
                'e3sm_file': { 'e3sm': e3sm_output_root + "/dycoms2_rf01.nc"},
                'cam_file': None,
                'sam_file': None,
                'wrf_file': None,
                'var_groups': [VariableGroupBase, VariableGroupWs]}

DYCOMS2_RF01_FIXED_SST = {'name': 'dycoms2_rf01_fixed_sst',
                          'description': "Copied from plotgen: Ran with a 5 min timestep and a 48-level grid",
                          'start_time': 2520, 'end_time': 2700,
                          'height_min_value': 0, 'height_max_value': 1200,

                          'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp'],
                          'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                "/SAM6.6/DYCOMS_RF01_fixed_sst/DYCOMS_RF01_96x96x320_LES_fixed_sst.nc"},
                          'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf01_fixed_sst_zm.nc',
                                         'zt': clubb_output_root + '/dycoms2_rf01_fixed_sst_zt.nc',
                                         'sfc': clubb_output_root + '/dycoms2_rf01_fixed_sst_sfc.nc'},
                          'coamps_benchmark_file': None,
                          'wrf_benchmark_file': None,
                          'clubb_r408_benchmark_file': None,
                          'clubb_hoc_benchmark_file': None,
                          'e3sm_file': None,
                          'cam_file': None,
                          'sam_file': None,
                          'wrf_file': None,
                          'var_groups': [VariableGroupBase]}

DYCOMS2_RF02_DO = {'name': 'dycoms2_rf02_do',
                   'description': "",
                   'start_time': 301, 'end_time': 360,
                   'height_min_value': 0, 'height_max_value': 1200,

                   'blacklisted_vars': [],
                   'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                         "/JULY_2017/DYCOMS_RF02_128x128x96_dr_nosed/DYCOMS_RF02_128x128x96_dr_nosed.nc"},
                   'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf02_do_zm.nc',
                                  'zt': clubb_output_root + '/dycoms2_rf02_do_zt.nc',
                                  'sfc': clubb_output_root + '/dycoms2_rf02_do_sfc.nc'},
                   'coamps_benchmark_file': None,
                   'wrf_benchmark_file': None,
                   'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_zm.nc',
                                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_zt.nc',
                                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_do_sfc.nc'},
                   'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf02_do_zm.nc',
                                      'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf02_do_zt.nc',
                                      'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf02_do_sfc.nc'},
                   'e3sm_file': None,
                   'cam_file': None,
                   'sam_file': {'sam': sam_output_root + "/DYCOMS_RF02_SAM_CLUBB.nc"},
                   'wrf_file': None,
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                                  VariableGroupKKMP]}

DYCOMS2_RF02_DS = {'name': 'dycoms2_rf02_ds',
                   'description': "",
                   'start_time': 301, 'end_time': 360,
                   'height_min_value': 0, 'height_max_value': 1200,

                   'blacklisted_vars': [],
                   'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                         "/JULY_2017/DYCOMS_RF02_128x128x96_dr_sed/DYCOMS_RF02_128x128x96_dr_sed.nc"},
                   'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf02_ds_zm.nc',
                                  'zt': clubb_output_root + '/dycoms2_rf02_ds_zt.nc',
                                  'sfc': clubb_output_root + '/dycoms2_rf02_ds_sfc.nc'},
                   'coamps_benchmark_file': None,
                   'wrf_benchmark_file': None,
                   'clubb_r408_benchmark_file': None,
                   'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_zm.nc',
                                      'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_zt.nc',
                                      'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_sfc.nc'},
                   'e3sm_file': {'e3sm': e3sm_output_root + "/dycoms2_rf02_ds.nc"},
                   'cam_file': None,
                   'sam_file': None,
                   'wrf_file': None,
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations,
                                  VariableGroupKKMP]}

DYCOMS2_RF02_ND = {'name': 'dycoms2_rf02_nd',
                   'description': "Copied from plotgen: ** Generated by doing a restart run after 7200 seconds. Note: "
                                  "t = 0 corresponds to start time of the restart run, not the original run. ** ",
                   'start_time': 301, 'end_time': 360,
                   'height_min_value': 0, 'height_max_value': 1200,

                   'blacklisted_vars': ['wprrp', 'wpNrp', 'corr_w_rr_1', 'corr_w_Nr_1'],
                   'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                         "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_nosed/DYCOMS_RF02_128x128x96_nodr_nosed.nc"},
                   'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf02_nd_zm.nc',
                                  'zt': clubb_output_root + '/dycoms2_rf02_nd_zt.nc',
                                  'sfc': clubb_output_root + '/dycoms2_rf02_nd_sfc.nc'},
                   'coamps_benchmark_file': None,
                   'wrf_benchmark_file': None,
                   'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_zm.nc',
                                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_zt.nc',
                                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_nd_sfc.nc'},
                   'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf02_nd_zm.nc',
                                      'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf02_nd_zt.nc',
                                      'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf02_nd_sfc.nc'},
                   'e3sm_file': None,
                   'cam_file': None,
                   'sam_file': None,
                   'wrf_file': None,
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

DYCOMS2_RF02_DS_RESTART = {'name': 'dycoms2_rf02_ds_restart',
                           'description': "Copied from plotgen: ** Uniform, coarse verticle grid spacing of 40 m. **",
                           'start_time': 181, 'end_time': 240,
                           'height_min_value': 0, 'height_max_value': 1200,

                           'blacklisted_vars': [],
                           'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                 "/JULY_2017/DYCOMS_RF02_128x128x96_dr_sed/DYCOMS_RF02_128x128x96_dr_sed.nc"},
                           'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf02_ds_restart_zm.nc',
                                          'zt': clubb_output_root + '/dycoms2_rf02_ds_restart_zt.nc',
                                          'sfc': clubb_output_root + '/dycoms2_rf02_ds_restart_sfc.nc'},
                           'coamps_benchmark_file': None,
                           'wrf_benchmark_file': None,
                           'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_zm.nc',
                                               'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_zt.nc',
                                               'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_ds_sfc.nc'},
                           'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_zm.nc',
                                              'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_zt.nc',
                                              'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf02_ds_sfc.nc'},
                           'e3sm_file': None,
                           'cam_file': None,
                           'sam_file': None,
                           'wrf_file': None,
                           'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP,
                                          VariableGroupCorrelations, VariableGroupKKMP]}

DYCOMS2_RF02_SO = {'name': 'dycoms2_rf02_so',
                   'description': "Copied from plotgen: " +
                                  "** WRF-type stretched (unevenly spaced) grid (grid_type = 3) ** ",
                   'start_time': 301, 'end_time': 360,
                   'height_min_value': 0, 'height_max_value': 1200,

                   'blacklisted_vars': ['wprrp', 'wpNrp'],
                   'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                         "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_sed/DYCOMS_RF02_128x128x96_nodr_sed.nc"},
                   'clubb_file': {'zm': clubb_output_root + '/dycoms2_rf02_so_zm.nc',
                                  'zt': clubb_output_root + '/dycoms2_rf02_so_zt.nc',
                                  'sfc': clubb_output_root + '/dycoms2_rf02_so_sfc.nc'},
                   'coamps_benchmark_file': None,
                   'wrf_benchmark_file': None,
                   'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_zm.nc',
                                       'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_zt.nc',
                                       'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/dycoms2_rf02_so_sfc.nc'},
                   'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/dycoms2_rf02_so_zm.nc',
                                      'zt': HOC_OUTPUT_ROOT + '/dycoms2_rf02_so_zt.nc',
                                      'sfc': HOC_OUTPUT_ROOT + '/dycoms2_rf02_so_sfc.nc'},
                   'e3sm_file': None,
                   'cam_file': None,
                   'sam_file': {'sam': sam_output_root + "/DYCOMS_RF02_SAM_CLUBB.nc"},
                   'wrf_file': None,
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

FIRE = {'name': 'fire',
        'description': "",
        'start_time': 61, 'end_time': 120,
        'height_min_value': 0, 'height_max_value': 1000,

        'blacklisted_vars': [],
        'sam_benchmark_file': None,
        'clubb_file': {'zm': clubb_output_root + '/fire_zm.nc',
                       'zt': clubb_output_root + '/fire_zt.nc',
                       'sfc': clubb_output_root + '/fire_sfc.nc'},
        'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/fire_coamps_sm.nc",
                                  'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/fire_coamps_sw.nc"},
        'wrf_benchmark_file': None,
        'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_zm.nc',
                            'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_zt.nc',
                            'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/fire_sfc.nc'},
        'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + "/fire_zm.nc",
                           'zt': HOC_OUTPUT_ROOT + '/fire_zt.nc',
                           'sfc': HOC_OUTPUT_ROOT + '/fire_sfc.nc'},
        'e3sm_file': None,
        'cam_file': None,
        'sam_file': None,
        'wrf_file': {'zm': wrf_output_root + "/fire_zm_wrf.nc",
                     'zt': wrf_output_root + "/fire_zt_wrf.nc",
                     'sfc': wrf_output_root + "/fire_sfc_wrf.nc"
                     },
        'var_groups': [VariableGroupBase, VariableGroupWs]}

# No budgets
GABLS2 = {'name': 'gabls2',
          'description': "",
          #'start_time': 2101, 'end_time': 2160,
          'start_time': 1, 'end_time': 3540,
          'height_min_value': 0, 'height_max_value': 8500,
          #'height_min_value': 0, 'height_max_value': 3500,
          #'height_min_value': 0, 'height_max_value': 1000,

          'blacklisted_vars': ['tau_zm', 'radht', 'Skw_zt', 'Skrt_zt', 'Skthl_zt', 'corr_w_chi_1', 'corr_chi_eta_1',
                               'rcp2', 'thlpthvp', 'rtpthvp'],
          'sam_benchmark_file': None,
          'clubb_file': {'zm': clubb_output_root + '/gabls2_zm.nc',
                         'zt': clubb_output_root + '/gabls2_zt.nc',
                         'sfc': clubb_output_root + '/gabls2_sfc.nc'},
          'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/gabls2_coamps_sm.nc",
                                    'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/gabls2_coamps_sw.nc",
                                    'sfc': COAMPS_BENCHMARK_OUTPUT_ROOT + "/gabls2_coamps_sfc.nc"},
          'wrf_benchmark_file': None,
          'clubb_r408_benchmark_file': None,
          'clubb_hoc_benchmark_file': None,
          'e3sm_file': None,
          'cam_file': None,
          'sam_file': None,
          'wrf_file': None,
          'var_groups': [VariableGroupBase]}

GABLS2_NIGHTLY = {'name': 'gabls2_nightly',
                  'description': "",
                  'start_time': 2101, 'end_time': 2160,
                  'height_min_value': 0, 'height_max_value': 2500,

                  'blacklisted_vars': [],
                  'sam_benchmark_file': None,
                  'clubb_file': {'zm': clubb_output_root + '/gabls2_zm.nc',
                                 'zt': clubb_output_root + '/gabls2_zt.nc',
                                 'sfc': clubb_output_root + '/gabls2_sfc.nc'},
                  'coamps_benchmark_file': None,
                  'wrf_benchmark_file': None,
                  'clubb_r408_benchmark_file': None,
                  'clubb_hoc_benchmark_file': None,
                  'e3sm_file': None,
                  'cam_file': None,
                  'sam_file': None,
                  'wrf_file': None,
                  'var_groups': [VariableGroupBase, VariableGroupScalars]}

GABLS3 = {'name': 'gabls3',
          'description': "",
          'start_time': 1081, 'end_time': 1200,
          'height_min_value': 0, 'height_max_value': 4970,

          'blacklisted_vars': [],
          'sam_benchmark_file': None,
          'clubb_file': {'zm': clubb_output_root + '/gabls3_zm.nc',
                         'zt': clubb_output_root + '/gabls3_zt.nc',
                         'sfc': clubb_output_root + '/gabls3_sfc.nc'},
          'coamps_benchmark_file': None,
          'wrf_benchmark_file': None,
          'clubb_r408_benchmark_file': None,
          'clubb_hoc_benchmark_file': None,
          'e3sm_file': None,
          'cam_file': None,
          'sam_file': None,
          'wrf_file': None,
          'var_groups': [VariableGroupBase]}

GABLS3_NIGHT = {'name': 'gabls3_night',
                'description': "Copied from plotgen: Uses a 5-min timestep with 48 levels",
                'start_time': 421, 'end_time': 480,
                'height_min_value': 0, 'height_max_value': 800,

                'blacklisted_vars': [],
                'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                        "/SAM6.6/GABLS3_NIGHT/gabls3_night.nc"},
                'clubb_file': {'zm': clubb_output_root + '/gabls3_night_zm.nc',
                               'zt': clubb_output_root + '/gabls3_night_zt.nc',
                               'sfc': clubb_output_root + '/gabls3_night_sfc.nc'},
                'coamps_benchmark_file': None,
                'wrf_benchmark_file': None,
                'clubb_r408_benchmark_file': None,
                'clubb_hoc_benchmark_file': None,
                'e3sm_file': None,
                'cam_file': None,
                'sam_file': None,
                'wrf_file': None,
                'var_groups': [VariableGroupBase]}

GATE_SHEAR_RLSF = {'name': 'gate_shear_rlsf',
                   'description': "",
                   'start_time': 540, 'end_time': 720,
                   'height_min_value': 0, 'height_max_value': 24000,

                   'blacklisted_vars': [],
                   'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                           "/SAM6.6/GATE_shear_rlsf/GATE_shear_rlsf_64x64x128_1km_5s.nc"},
                   'clubb_file': None,
                   'coamps_benchmark_file': None,
                   'wrf_benchmark_file': None,
                   'clubb_r408_benchmark_file': None,
                   'clubb_hoc_benchmark_file': None,
                   'e3sm_file': None,
                   'cam_file': None,
                   'sam_file': {'sam': sam_output_root + "/GATE_SAM_CLUBB.nc"},
                   'wrf_file': None,
                   'var_groups': [VariableGroupBase]}

# Use to plot IOP forced SAM runs
IOP = {'name': 'iop',
       'description': "",
       'start_time': 181, 'end_time': 1440,
       'height_min_value': 0, 'height_max_value': 27750,

       'blacklisted_vars': [],
       'clubb_datasets': None,
       'sam_benchmark_file': None,
       'clubb_file': None,
       'coamps_benchmark_file': None,
       'wrf_benchmark_file': None,
       'clubb_r408_benchmark_file': None,
       'clubb_hoc_benchmark_file': None,
       'e3sm_file': None,
       'cam_file': None,
       'var_groups': [VariableGroupBase, VariableGroupSamProfiles]}

JUN25_ALTOCU = {'name': 'jun25_altocu',
                'description': "",
                'start_time': 181, 'end_time': 240,
                'height_min_value': 4825, 'height_max_value': 7290,

                'blacklisted_vars': ['Ngm', 'wprrp', 'wpNrp'],
                'sam_benchmark_file': None,
                'clubb_file': {'zm': clubb_output_root + '/jun25_altocu_zm.nc',
                               'zt': clubb_output_root + '/jun25_altocu_zt.nc',
                               'sfc': clubb_output_root + '/jun25_altocu_sfc.nc'},
                'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/jun25_altocu_qc3_coamps_sm.nc",
                                          'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/jun25_altocu_qc3_coamps_sw.nc"},
                'wrf_benchmark_file': None,
                'clubb_r408_benchmark_file': None,
                'clubb_hoc_benchmark_file': None,
                'e3sm_file': None,
                'cam_file': None,
                'sam_file': None,
                'wrf_file': None,
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

LBA = {'name': 'lba',
       'description': "Note that sam-plotgen plots up to a height of 16000 not 12000.\n"
                      "Copied from plotgen: SAM-LES uses Morrison microphysics " +
                      "and CLUBB standalone uses COAMPS microphysics",
       'start_time': 300, 'end_time': 360,
       'height_min_value': 0, 'height_max_value': 14000,

       'blacklisted_vars': ['wprrp', 'wpNrp', 'Ngm'],
       'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                             "/JULY_2017/LBA_128kmx128kmx128_1km_Morrison/LBA_128kmx128kmx128_1km_Morrison.nc"},
       'clubb_file': {'zm': clubb_output_root + '/lba_zm.nc',
                      'zt': clubb_output_root + '/lba_zt.nc',
                      'sfc': clubb_output_root + '/lba_sfc.nc',
                      'subcolumns': clubb_output_root + '/lba_nl_lh_sample_points_2D.nc'},
       'coamps_benchmark_file': None,
       'wrf_benchmark_file': None,
       'clubb_r408_benchmark_file': None,
       'clubb_hoc_benchmark_file': None,
       'e3sm_file': None,
       'cam_file': None,
       'sam_file': {'sam': sam_output_root + "/LBA_SAM_CLUBB.nc"},
       'wrf_file': None,
       'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP, VariableGroupWs]}

MC3E = {'name': 'mc3e',
        'description': "",
        'start_time': 60, 'end_time': 64800,
        'height_min_value': 0, 'height_max_value': 18000,

        'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
        'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                "/MC3E_r1359_128x128x128_1km_Morrison/MC3E.nc"},
        'clubb_file': {'zm': clubb_output_root + '/mc3e_zm.nc',
                       'zt': clubb_output_root + '/mc3e_zt.nc',
                       'sfc': clubb_output_root + '/mc3e_sfc.nc',
                       'subcolumns': clubb_output_root + '/mc3e_nl_lh_sample_points_2D.nc'},
        'coamps_benchmark_file': None,
        'wrf_benchmark_file': None,
        'clubb_r408_benchmark_file': None,
        'clubb_hoc_benchmark_file': None,
        'e3sm_file': None,
        'cam_file': None,
        'sam_file': None,
        'wrf_file': None,
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_A = {'name': 'mpace_a',
           'description': "Copied from plotgen: SAM-LES and CLUBB standalone use Morrison microphysics",
           'start_time': 4141, 'end_time': 4320,
           'height_min_value': 0, 'height_max_value': 10000,

           'blacklisted_vars': ['Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'Ngm', 'wpNrp'],
           'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                   "/SAM6.6/MPACE_A/MPACE_A_128x128x69_morr_CEM.nc"},
           'clubb_file': {'zm': clubb_output_root + '/mpace_a_zm.nc',
                          'zt': clubb_output_root + '/mpace_a_zt.nc',
                          'sfc': clubb_output_root + '/mpace_a_sfc.nc'},
           'coamps_benchmark_file': None,
           'wrf_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B = {'name': 'mpace_b',
           'description': "Copied from plotgen: **The nightly simulation uses COAMPS microphysics**",
           'start_time': 541, 'end_time': 720,
           'height_min_value': 0, 'height_max_value': 2750,

           'blacklisted_vars': ['Ngm', 'wpNrp'],
           'sam_benchmark_file': None,
           'clubb_file': {'zm': clubb_output_root + '/mpace_b_zm.nc',
                          'zt': clubb_output_root + '/mpace_b_zt.nc',
                          'sfc': clubb_output_root + '/mpace_b_sfc.nc'},
           'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/mpace_b_coamps_sm.nc",
                                     'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/mpace_b_coamps_sw.nc",
                                     'sfc': COAMPS_BENCHMARK_OUTPUT_ROOT + "/mpace_b_coamps_sfc.nc"},
           'wrf_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B_SILHS = {'name': 'mpace_b_silhs',
                 'description': "",
                 'start_time': 541, 'end_time': 720,
                 'height_min_value': 0, 'height_max_value': 2750,

                 'blacklisted_vars': ['Ngm', 'wpNrp'],
                 'sam_benchmark_file': None,
                 'clubb_file': {'zm': clubb_output_root + '/mpace_b_silhs_zm.nc',
                                'zt': clubb_output_root + '/mpace_b_silhs_zt.nc',
                                'sfc': clubb_output_root + '/mpace_b_silhs_sfc.nc',
                                'subcolumns': clubb_output_root + '/mpace_b_silhs_nl_lh_sample_points_2D.nc'},
                 'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/mpace_b_coamps_sm.nc",
                                           'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/mpace_b_coamps_sw.nc"},
                 'wrf_benchmark_file': None,
                 'clubb_r408_benchmark_file': None,
                 'clubb_hoc_benchmark_file': None,
                 'e3sm_file': None,
                 'cam_file': None,
                 'sam_file': None,
                 'wrf_file': None,
                 'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

NOV11_ALTOCU = {'name': 'nov11_altocu',
                'description': "",
                'start_time': 91, 'end_time': 150,
                'height_min_value': 4160, 'height_max_value': 6150,

                'blacklisted_vars': ['Ngm'],
                'sam_benchmark_file': None,
                'clubb_file': {'zm': clubb_output_root + '/nov11_altocu_zm.nc',
                               'zt': clubb_output_root + '/nov11_altocu_zt.nc',
                               'sfc': clubb_output_root + '/nov11_altocu_sfc.nc'},
                'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/nov11_altocu_coamps_sm.nc",
                                          'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/nov11_altocu_coamps_sw.nc"},
                'wrf_benchmark_file': None,
                'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_zm.nc',
                                    'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_zt.nc',
                                    'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/nov11_altocu_sfc.nc'},
                'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/nov11_altocu_zm.nc',
                                   'zt': HOC_OUTPUT_ROOT + '/nov11_altocu_zt.nc',
                                   'sfc': HOC_OUTPUT_ROOT + '/nov11_altocu_sfc.nc'},
                'e3sm_file': None,
                'cam_file': None,
                'sam_file': None,
                'wrf_file': None,
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

RICO = {'name': 'rico',
        'description': "Cam output may differ from plotgen due to a difference in time averaging.",
        'start_time': 4201, 'end_time': 4320,
        'height_min_value': 0, 'height_max_value': 5000,

        'blacklisted_vars': [],
        'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                "/JULY_2017/RICO_256x256x100_drizzle/RICO_256x256x100_drizzle.nc"},
        'clubb_file': {'zm': clubb_output_root + '/rico_zm.nc',
                       'zt': clubb_output_root + '/rico_zt.nc',
                       'sfc': clubb_output_root + '/rico_sfc.nc'},
        'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/rico_coamps_sm.nc",
                                  'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/rico_coamps_sw.nc"},
        'wrf_benchmark_file': None,
        'clubb_r408_benchmark_file': None,
        'clubb_hoc_benchmark_file': None,
        'e3sm_file': {'e3sm': e3sm_output_root + "/rico.nc"},
        'cam_file': {'cam': cam_output_root + "/rico_cam.nc"},
        'sam_file': {'sam': sam_output_root + "/RICO_256x256x100_drizzle.nc"},
        'wrf_file': None,
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupWs, VariableGroupCorrelations,
                       VariableGroupKKMP]}

RICO_SILHS = {'name': 'rico_silhs',
              'description': "Copied from plotgen: CLUBB and SAM use Khairoutdinov-Kogan microphysics",
              'start_time': 4201, 'end_time': 4320,
              'height_min_value': 0, 'height_max_value': 4500,

              'blacklisted_vars': ['wpNrp'],
              'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                    "/JULY_2017/RICO_256x256x100_drizzle/RICO_256x256x100_drizzle.nc"},
              'clubb_file': {'zm': clubb_output_root + '/rico_silhs_zm.nc',
                             'zt': clubb_output_root + '/rico_silhs_zt.nc',
                             'sfc': clubb_output_root + '/rico_silhs_sfc.nc',
                             'subcolumns': clubb_output_root + '/rico_silhs_nl_lh_sample_points_2D.nc'},
              'coamps_benchmark_file': {'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/rico_coamps_sm.nc",
                                        'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/rico_coamps_sw.nc"},
              'wrf_benchmark_file': None,
              'clubb_r408_benchmark_file': None,
              'clubb_hoc_benchmark_file': None,
              'e3sm_file': None,
              'cam_file': None,
              'sam_file': None,
              'wrf_file': None,
              'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupWs, VariableGroupCorrelations,
                             VariableGroupKKMP]}

NEUTRAL = {'name': 'neutral',
          'description': "",
          'start_time': 181, 'end_time': 360,
          'height_min_value': 0, 'height_max_value': 1500,
  
          'blacklisted_vars': [],
          'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                  "/NEUTRAL/NEUTRAL_96x96x96_32m_10m_LES.nc"},
          'clubb_file': {'zm': clubb_output_root + '/neutral_zm.nc',
                         'zt': clubb_output_root + '/neutral_zt.nc',
                         'sfc': clubb_output_root + '/neutral_sfc.nc'},
          'coamps_benchmark_file': None,
          'wrf_benchmark_file': None,
          'clubb_r408_benchmark_file': None,
          'clubb_hoc_benchmark_file': None,
          'e3sm_file': None, 
          'cam_file': None,
          'sam_file': None,
          'wrf_file': None,
          'var_groups': [VariableGroupBase, VariableGroupWs]}

TWP_ICE = {'name': 'twp_ice',
           'description': "Copied from plotgen: Both vertical and horizontal fluxes applied to THLM and RTM for LES. "
                          "LES nudged U, V, RTM and THLM toward observed values. Forcings for LES derived from 10mb "
                          "forcing data.",
           'start_time': 60, 'end_time': 9900,
           'height_min_value': 0, 'height_max_value': 19000,

           'blacklisted_vars': ['rtp3', 'Skrt_zt', 'Skthl_zt', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
           'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                   "/TWP_ICE_r1315_128x128x128_1km_Morrison/TWP_ICE.nc"},
           'clubb_file': {'zm': clubb_output_root + '/twp_ice_zm.nc',
                          'zt': clubb_output_root + '/twp_ice_zt.nc',
                          'sfc': clubb_output_root + '/twp_ice_sfc.nc',
                          'subcolumns': clubb_output_root + '/twp_ice_nl_lh_sample_points_2D.nc'},
           'coamps_benchmark_file': None,
           'wrf_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_file': None,
           'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupIceMP]}

WANGARA = {'name': 'wangara',
           'description': "Note that COAMPS benchmark data is actually RAMS data by default.",
           'start_time': 181, 'end_time': 240,
           'height_min_value': 0, 'height_max_value': 1900,

           'blacklisted_vars': ['Ngm'],
           'sam_benchmark_file': {'sam_benchmark': SAM_BENCHMARK_OUTPUT_ROOT +
                                                   "/WANGARA/WANGARA_64x64x80_100m_40m_LES.nc"},
           'clubb_file': {'zm': clubb_output_root + '/wangara_zm.nc',
                          'zt': clubb_output_root + '/wangara_zt.nc',
                          'sfc': clubb_output_root + '/wangara_sfc.nc'},
           'coamps_benchmark_file': {'sw': COAMPS_BENCHMARK_OUTPUT_ROOT + "/wangara_rams.nc",
                                     'sm': COAMPS_BENCHMARK_OUTPUT_ROOT + "/wangara_rams.nc"},
           'wrf_benchmark_file': None,
           'clubb_r408_benchmark_file': {'zm': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_zm.nc',
                               'zt': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_zt.nc',
                               'sfc': R408_OUTPUT_ROOT + '/Chris_Golaz_best_ever/wangara_sfc.nc'},
           'clubb_hoc_benchmark_file': {'zm': HOC_OUTPUT_ROOT + '/wangara_zm.nc',
                              'zt': HOC_OUTPUT_ROOT + '/wangara_zt.nc',
                              'sfc': HOC_OUTPUT_ROOT + '/wangara_sfc.nc'},
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_file': {'zm': wrf_output_root + "/wangara_zm_wrf.nc",
                        'zt': wrf_output_root + "/wangara_zt_wrf.nc",
                        'sfc': wrf_output_root + "/wangara_sfc_wrf.nc"
                        },
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20170627 = {'name': 'lasso_20170627',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
              WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2017-06-27/wrf_lasso_stats_2017-06-27.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None, 
           'wrf_file': {'zm': clubb_output_root + '/lasso_2017-06-27_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2017-06-27_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2017-06-27_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2017-06-27_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20170717 = {'name': 'lasso_20170717',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2017-07-17/wrf_lasso_stats_2017-07-17.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2017-07-17_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2017-07-17_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2017-07-17_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2017-07-17_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20170728 = {'name': 'lasso_20170728',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2017-07-28/wrf_lasso_stats_2017-07-28.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2017-07-28_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2017-07-28_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2017-07-28_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2017-07-28_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20170923 = {'name': 'lasso_20170923',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2017-09-23/wrf_lasso_stats_2017-09-23.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2017-09-23_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2017-09-23_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2017-09-23_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2017-09-23_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20180911 = {'name': 'lasso_20180911',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2018-09-11/wrf_lasso_stats_2018-09-11.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2018-09-11_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2018-09-11_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2018-09-11_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2018-09-11_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20180917 = {'name': 'lasso_20180917',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2018-09-17/wrf_lasso_stats_2018-09-17.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None, 
           'wrf_file': {'zm': clubb_output_root + '/lasso_2018-09-17_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2018-09-17_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2018-09-17_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2018-09-17_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20180918 = {'name': 'lasso_20180918',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2018-09-18/wrf_lasso_stats_2018-09-18.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2018-09-18_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2018-09-18_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2018-09-18_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2018-09-18_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

LASSO_20181002 = {'name': 'lasso_20181002',
           'description': "Comparing WRF-CLUBB output to WRF-LASSO output.",
           'start_time': 301, 'end_time': 600,
           'height_min_value': 0, 'height_max_value': 4000,
           'blacklisted_vars': [],
           'e3sm_file': None,
           'cam_file': None,
           'sam_file': None,
           'wrf_benchmark_file': {'lasso_benchmark':
             WRF_LASSO_BENCHMARK_OUTPUT_ROOT + "/2018-10-02/wrf_lasso_stats_2018-10-02.nc"},
           'sam_benchmark_file': None,
           'coamps_benchmark_file': None,
           'clubb_r408_benchmark_file': None,
           'clubb_hoc_benchmark_file': None,
           'clubb_file': None,
           'wrf_file': {'zm': clubb_output_root + '/lasso_2018-10-02_zm_wrf.nc',
                          'zt': clubb_output_root + '/lasso_2018-10-02_zt_wrf.nc',
                          'sfc': clubb_output_root + '/lasso_2018-10-02_sfc_wrf.nc',
                          'subcolumns': clubb_output_root + '/lasso_2018-10-02_nl_lh_sample_points_2D.nc'},
           'var_groups': [VariableGroupBase, VariableGroupWs]}

# DO NOT EDIT THIS LIST UNLESS YOU ARE ADDING A NEW CASE. NEVER REMOVE CASES FROM THIS LIST.
# You may define a subset of cases at the end of this file.
ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX,
                 BOMEX,
                 CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14,
                 DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST, DYCOMS2_RF02_DO,
                 DYCOMS2_RF02_DS, DYCOMS2_RF02_DS_RESTART,
                 DYCOMS2_RF02_ND, DYCOMS2_RF02_SO,
                 FIRE,
                 GABLS2, GABLS2_NIGHTLY, GABLS3, GABLS3_NIGHT, GATE_SHEAR_RLSF,
                 # IOP,
                 JUN25_ALTOCU,
                 LBA,
                 MC3E, MPACE_A, MPACE_B, MPACE_B_SILHS,
                 NEUTRAL, NOV11_ALTOCU,
                 RICO, RICO_SILHS,
                 TWP_ICE,
                 WANGARA,
                 LASSO_20170627, LASSO_20170717, LASSO_20170728, LASSO_20170923,
                 LASSO_20180911, LASSO_20180917, LASSO_20180918, LASSO_20181002
                 ]

CASES_TO_PLOT = ALL_CASES
# If uncommented, this line will override the real CASES_TO_PLOT given above, forcing pyplotgen to only plot some cases.
# CASES_TO_PLOT = [ARM]
# CASES_TO_PLOT = CASES_TO_PLOT[:3]
