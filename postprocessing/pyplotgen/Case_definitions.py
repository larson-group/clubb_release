import os

from VariableGroupBase import VariableGroupBase
from VariableGroupCorrelations import VariableGroupCorrelations
from VariableGroupIceMP import VariableGroupIceMP
from VariableGroupKKMP import VariableGroupKKMP
from VariableGroupLiquidMP import VariableGroupLiquidMP
from VariableGroupWs import VariableGroupWs

SAM_OUTPUT_ROOT = os.path.dirname(os.path.realpath(__file__)) + "/sam_benchmark_runs"

# the 'name' parameter must be the same as the filename without the extention.
#   E.g. to use lba_zt.nc and lba_zm.nc the case's name must be 'lba'

ARM = {'name': 'arm', 'start_time': 481, 'end_time': 540, 'height_min_value': 0, 'height_max_value': 3500, 'enabled': True, 'disable_budgets': False,
       'blacklisted_vars': [], 'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/ARM_96x96x110/GCSSARM_96x96x110_67m_40m_1s.nc",
       'var_groups': [VariableGroupBase]}

ARM_97 = {'name': 'arm_97', 'start_time': 4321, 'end_time': 5580, 'height_min_value': 0, 'height_max_value': 18000, 'enabled': True, 'disable_budgets': False,
          'blacklisted_vars': ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp'],
          'sam_file': SAM_OUTPUT_ROOT + "/ARM97_r1315_128x128x128_1km_Morrison/ARM9707.nc",
          'var_groups': [VariableGroupBase]}

ASTEX_A209 = {'name': 'astex_a209', 'start_time': 2340, 'end_time': 2400, 'height_min_value': 0, 'height_max_value': 6000, 'enabled': True, 'disable_budgets': False,
              'blacklisted_vars': [],
              'sam_file': None,
              'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations, VariableGroupKKMP]}

ATEX = {'name': 'atex', 'start_time': 421, 'end_time': 480, 'height_min_value': 0, 'height_max_value': 2500, 'enabled': True, 'disable_budgets': False,
        'blacklisted_vars': [],
        'sam_file': None,
        'var_groups': [VariableGroupBase, VariableGroupWs]}

BOMEX = {'name': 'bomex', 'start_time': 181, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 2500, 'enabled': True, 'disable_budgets': False,
         'blacklisted_vars': [],
         'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/BOMEX_64x64x75/BOMEX_64x64x75_100m_40m_1s.nc",
         'var_groups': [VariableGroupBase, VariableGroupWs]}

CGILS_S6 = {'name': 'cgils_s6', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0, 'height_max_value': 5950, 'enabled': True, 'disable_budgets': False,
            'blacklisted_vars': ['Ngm', 'rgm', 'rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
            'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s6/ctl_s6_96x96x128_100m_DRZ_N100_tqndg.nc",
            'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S11 = {'name': 'cgils_s11', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0, 'height_max_value': 5950, 'enabled': True, 'disable_budgets': False,
             'blacklisted_vars': ['Ngm', 'rgm', 'rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s11/ctl_s11_96x96x320_50m_DRZ_N100_ref.nc",
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CGILS_S12 = {'name': 'cgils_s12', 'start_time': 12960, 'end_time': 14400, 'height_min_value': 0, 'height_max_value': 5950, 'enabled': True, 'disable_budgets': False,
             'blacklisted_vars': ['Ngm', 'rgm', 'rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'wprrp', 'wpNrp'],
             'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/CLOUD_FEEDBACK_s12/ctl_s12_96x96x192_25m_DRZ_N100_fixnudge.nc",
             'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_NOV02 = {'name': 'clex9_nov02', 'start_time': 181, 'end_time': 240, 'height_min_value': 3072, 'height_max_value': 6072, 'enabled': True, 'disable_budgets': False,
               'blacklisted_vars': [],
               'sam_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

CLEX9_OCT14 = {'name': 'clex9_oct14', 'start_time': 181, 'end_time': 240, 'height_min_value': 2188, 'height_max_value': 6688, 'enabled': True, 'disable_budgets': False,
               'blacklisted_vars': [],
               'sam_file': None,
               'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

DYCOMS2_RF01 = {'name': 'dycoms2_rf01', 'start_time': 181, 'end_time': 240, 'height_min_value': 0, 'height_max_value': 1200, 'enabled': True, 'disable_budgets': False,
                'blacklisted_vars': [],
                'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF01_96x96x320/DYCOMS_RF01_96x96x320.nc",
                'var_groups': [VariableGroupBase, VariableGroupWs]}

DYCOMS2_RF01_FIXED_SST = {'name': 'dycoms2_rf01_fixed_sst', 'start_time': 2520, 'end_time': 2700, 'height_min_value': 0, 'height_max_value': 1200,
                          'enabled': True, 'disable_budgets': False,
                          'blacklisted_vars': ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp'],
                          'sam_file': SAM_OUTPUT_ROOT + "/SAM6.6/DYCOMS_RF01_fixed_sst/DYCOMS_RF01_96x96x320_LES_fixed_sst.nc",
                          'var_groups': [VariableGroupBase]}

DYCOMS2_RF02_DO = {'name': 'dycoms2_rf02_do', 'start_time': 301, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 1200, 'enabled': True, 'disable_budgets': False,
                   'blacklisted_vars': [],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_dr_nosed/DYCOMS_RF02_128x128x96_dr_nosed.nc",
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations, VariableGroupKKMP]}

DYCOMS2_RF02_DS = {'name': 'dycoms2_rf02_ds', 'start_time': 301, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 1200, 'enabled': True, 'disable_budgets': False,
                   'blacklisted_vars': [],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_dr_sed/DYCOMS_RF02_128x128x96_dr_sed.nc",
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupCorrelations, VariableGroupKKMP]}

DYCOMS2_RF02_ND = {'name': 'dycoms2_rf02_nd', 'start_time': 301, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 1200, 'enabled': True, 'disable_budgets': False,
                   'blacklisted_vars': ['wprrp', 'wpNrp', 'corr_w_rr_1', 'corr_w_Nr_1'],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_nosed/DYCOMS_RF02_128x128x96_nodr_nosed.nc",
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

DYCOMS2_RF02_SO = {'name': 'dycoms2_rf02_so', 'start_time': 301, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 1200, 'enabled': True, 'disable_budgets': False,
                   'blacklisted_vars': ['wprrp', 'wpNrp'],
                   'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/DYCOMS_RF02_128x128x96_nodr_sed/DYCOMS_RF02_128x128x96_nodr_sed.nc",
                   'var_groups': [VariableGroupBase, VariableGroupWs, VariableGroupLiquidMP, VariableGroupKKMP]}

FIRE = {'name': 'fire', 'start_time': 61, 'end_time': 120, 'height_min_value': 0, 'height_max_value': 1000, 'enabled': True, 'disable_budgets': False,
        'blacklisted_vars': [],
        'sam_file': None,
        'var_groups': [VariableGroupBase, VariableGroupWs]}

# No budgets
GABLS2 = {'name': 'gabls2', 'start_time': 2101, 'end_time': 2160, 'height_min_value': 0, 'height_max_value': 2500, 'enabled': True, 'disable_budgets': True,
          'blacklisted_vars': ['tau_zm', 'radht', 'Skw_zt', 'Skrt_zt', 'Skthl_zt', 'corr_w_chi_1', 'corr_chi_eta_1', 'rcp2', 'thlpthvp', 'rtpthvp'],
          'sam_file': None,
          'var_groups': [VariableGroupBase]}

GABLS3 = {'name': 'gabls3', 'start_time': 1081, 'end_time': 1200, 'height_min_value': 0, 'height_max_value': 4970, 'enabled': True, 'disable_budgets': False,
          'blacklisted_vars': [],
          'sam_file': None,
          'var_groups': [VariableGroupBase]}

GABLS3_NIGHT = {'name': 'gabls3_night', 'start_time': 421, 'end_time': 480, 'height_min_value': 0, 'height_max_value': 800, 'enabled': True, 'disable_budgets': False,
                'blacklisted_vars': [],
                'sam_file': None,
                'var_groups': [VariableGroupBase]}

JUN25_ALTOCU = {'name': 'jun25_altocu', 'start_time': 181, 'end_time': 240, 'height_min_value': 4808, 'height_max_value': 7308, 'enabled': True, 'disable_budgets': False,
                'blacklisted_vars': ['Ngm', 'wprrp', 'wpNrp'],
                'sam_file': None,
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

LBA = {'name': 'lba', 'start_time': 300, 'end_time': 360, 'height_min_value': 0, 'height_max_value': 12000, 'enabled': True, 'disable_budgets': False,
       'blacklisted_vars': ['wprrp', 'wpNrp'],
       'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/LBA_128kmx128kmx128_1km_Morrison/LBA_128kmx128kmx128_1km_Morrison.nc",
       'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP, VariableGroupWs]}

MC3E = {'name': 'mc3e', 'start_time': 1, 'end_time': 64800, 'height_min_value': 0, 'height_max_value': 18000, 'enabled': True, 'disable_budgets': False,
        'blacklisted_vars': ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'Ngm', 'wprrp', 'wpNrp'],
        'sam_file': SAM_OUTPUT_ROOT + "/MC3E_r1359_128x128x128_1km_Morrison/MC3E.nc",
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_A = {'name': 'mpace_a', 'start_time': 4141, 'end_time': 4320, 'height_min_value': 0, 'height_max_value': 10000, 'enabled': True, 'disable_budgets': False,
           'blacklisted_vars': ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'Ngm', 'wpNrp'],
           'sam_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B = {'name': 'mpace_b', 'start_time': 541, 'end_time': 720, 'height_min_value': 0, 'height_max_value': 2750, 'enabled': True, 'disable_budgets': False,
           'blacklisted_vars': ['Ngm', 'wpNrp'],
           'sam_file': None,
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

MPACE_B_SILHS = {'name': 'mpace_b_silhs', 'start_time': 541, 'end_time': 720, 'height_min_value': 0, 'height_max_value': 2750, 'enabled': True, 'disable_budgets': False,
                 'blacklisted_vars': ['Ngm', 'wpNrp'],
                 'sam_file': None,
                 'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

NOV11_ALTOCU = {'name': 'nov11_altocu', 'start_time': 91, 'end_time': 150, 'height_min_value': 4150, 'height_max_value': 6150, 'enabled': True, 'disable_budgets': False,
                'blacklisted_vars': ['Ngm'],
                'sam_file': None,
                'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

RICO = {'name': 'rico', 'start_time': 4201, 'end_time': 4320, 'height_min_value': 0, 'height_max_value': 4000, 'enabled': True, 'disable_budgets': False,
        'blacklisted_vars': ['wpNrp'],
        'sam_file': SAM_OUTPUT_ROOT + "/JULY_2017/RICO_256x256x100_drizzle/RICO_256x256x100_drizzle.nc",
        'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupWs, VariableGroupCorrelations, VariableGroupKKMP]}

TWP_ICE = {'name': 'twp_ice', 'start_time': 1, 'end_time': 9900, 'height_min_value': 0, 'height_max_value': 19000, 'enabled': True, 'disable_budgets': False,
           'blacklisted_vars': ['rtp3', 'thlp3', 'rtpthvp', 'thlpthvp', 'Ngm', 'wprrp', 'wpNrp'],
           'sam_file': SAM_OUTPUT_ROOT + "/TWP_ICE_r1315_128x128x128_1km_Morrison/TWP_ICE.nc",
           'var_groups': [VariableGroupBase, VariableGroupLiquidMP, VariableGroupIceMP]}

WANGARA = {'name': 'wangara', 'start_time': 181, 'end_time': 240, 'height_min_value': 0, 'height_max_value': 1900, 'enabled': True, 'disable_budgets': False,
           'blacklisted_vars': [],
           'sam_file': None,
           'var_groups': [VariableGroupBase, VariableGroupWs]}

ALL_CASES = [ARM, ARM_97, ASTEX_A209, ATEX, BOMEX, CGILS_S6, CGILS_S11, CGILS_S12, CLEX9_NOV02, CLEX9_OCT14, DYCOMS2_RF01, DYCOMS2_RF01_FIXED_SST,
             DYCOMS2_RF02_DO, DYCOMS2_RF02_DS, DYCOMS2_RF02_ND, DYCOMS2_RF02_SO, FIRE, GABLS2, GABLS3, GABLS3_NIGHT, JUN25_ALTOCU, LBA, MC3E,
             MPACE_A, MPACE_B, MPACE_B_SILHS, NOV11_ALTOCU, RICO, TWP_ICE, WANGARA]

# ALL_CASES = [ARM_97]
