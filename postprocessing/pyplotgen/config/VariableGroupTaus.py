"""
:author: Ben Stephens
:date: December 2020
"""

from src.VariableGroup import VariableGroup


class VariableGroupTaus(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False,
                 background_rcm=False, background_rcm_folder=None):
        """
        """
        self.name = "tau variables"
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['bv_freq_sqd_smth'],
                'sam': ['BV_FREQ_SQD_SMTH'],
                'coamps': ['bv_freq_sqd_smth'],
                'r408': ['bv_freq_sqd_smth'],
                'hoc': ['bv_freq_sqd_smth'],
                'e3sm': ['bv_freq_sqd_smth'],
                'cam': ['bv_freq_sqd_smth'],
                'wrf': ['bv_freq_sqd_smth'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['bv_freq_sqd_mixed'],
                'sam': [],
                'coamps': [],
                'r408': ['bv_freq_sqd_mixed'],
                'hoc': ['bv_freq_sqd_mixed'],
                'e3sm': ['bv_freq_sqd_mixed'],
                'cam': ['bv_freq_sqd_mixed'],
                'wrf': ['bv_freq_sqd_mixed'],
                },
                'title': 'Mixed Brunt-Vaisala frequency squared',
                'axis_title': 'bv_freq_sqd_mixed [$\mathrm{1/s^2}$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['bv_freq_sqd_moist'],
                'sam': [],
                'coamps': [],
                'r408': ['bv_freq_sqd_moist'],
                'hoc': ['bv_freq_sqd_moist'],
                'e3sm': ['bv_freq_sqd_moist'],
                'cam': ['bv_freq_sqd_moist'],
                'wrf': ['bv_freq_sqd_moist'],
                },
                'title': 'Brunt-Vaisala frequency squared in moist air',
                'axis_title': 'bv_freq_sqd_moist [$\mathrm{1/s^2}$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['bv_freq_sqd_dry'],
                'sam': [],
                'coamps': [],
                'r408': ['bv_freq_sqd_dry'],
                'hoc': ['bv_freq_sqd_dry'],
                'e3sm': ['bv_freq_sqd_dry'],
                'cam': ['bv_freq_sqd_dry'],
                'wrf': ['bv_freq_sqd_dry'],
                },
                'title': 'Brunt-Vaisala frequency squared in dry air',
                'axis_title': 'bv_freq_sqd_dry [$\mathrm{1/s^2}$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['Richardson_num'],
                'sam': ['RICHARDSON_NUM'], # Dummy
                'coamps': ['Richardson_num'], # Dummy
                'r408': ['Richardson_num'], # Dummy
                'hoc': ['Richardson_num'], # Dummy
                'e3sm': ['Richardson_num'], # Dummy
                'cam': ['Richardson_num'], # Dummy
                'wrf': ['Richardson_num'], # Dummy
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['Kh_zm'],
                'sam': ['KH_ZM'],
                'coamps': ['Kh_zm'],
                'r408': ['Kh_zm'],
                'hoc': ['Kh_zm'],
                'e3sm': ['Kh_zm'],
                'cam': ['Kh_zm'],
                'wrf': ['Kh_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_zm'],
                'sam': ['INVRS_TAU_ZM'],
                'coamps': ['invrs_tau_zm'],
                'r408': ['invrs_tau_zm'],
                'hoc': ['invrs_tau_zm'],
                'e3sm': ['invrs_tau_zm'],
                'cam': ['invrs_tau_zm'],
                'wrf': ['invrs_tau_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_xp2_zm'],
                'sam': ['INVRS_TAU_XP2_ZM'],
                'coamps': ['invrs_tau_xp2_zm'],
                'r408': ['invrs_tau_xp2_zm'],
                'hoc': ['invrs_tau_xp2_zm'],
                'e3sm': ['invrs_tau_xp2_zm'],
                'cam': ['invrs_tau_xp2_zm'],
                'wrf': ['invrs_tau_xp2_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_wp2_zm'],
                'sam': ['INVRS_TAU_WP2_ZM'],
                'coamps': ['invrs_tau_wp2_zm'],
                'r408': ['invrs_tau_wp2_zm'],
                'hoc': ['invrs_tau_wp2_zm'],
                'e3sm': ['invrs_tau_wp2_zm'],
                'cam': ['invrs_tau_wp2_zm'],
                'wrf': ['invrs_tau_wp2_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_wpxp_zm'],
                'sam': ['INVRS_TAU_WPXP_ZM'],
                'coamps': ['invrs_tau_wpxp_zm'],
                'r408': ['invrs_tau_wpxp_zm'],
                'hoc': ['invrs_tau_wpxp_zm'],
                'e3sm': ['invrs_tau_wpxp_zm'],
                'cam': ['invrs_tau_wpxp_zm'],
                'wrf': ['invrs_tau_wpxp_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_wp3_zm'],
                'sam': ['INVRS_TAU_WP3_ZM'],
                'coamps': ['invrs_tau_wp3_zm'],
                'r408': ['invrs_tau_wp3_zm'],
                'hoc': ['invrs_tau_wp3_zm'],
                'e3sm': ['invrs_tau_wp3_zm'],
                'cam': ['invrs_tau_wp3_zm'],
                'wrf': ['invrs_tau_wp3_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['invrs_tau_no_N2_zm'],
                'sam': ['INVRS_TAU_NO_N2_ZM'],
                'coamps': ['invrs_tau_no_N2_zm'],
                'r408': ['invrs_tau_no_N2_zm'],
                'hoc': ['invrs_tau_no_N2_zm'],
                'e3sm': ['invrs_tau_no_N2_zm'],
                'cam': ['invrs_tau_no_N2_zm'],
                'wrf': ['invrs_tau_no_N2_zm'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['C6_term'],
                'sam': ['C6_TERM'], # Dummy
                'coamps': ['C6_term'], # Dummy
                'r408': ['C6_term'], # Dummy
                'hoc': ['C6_term'], # Dummy
                'e3sm': ['C6_term'], # Dummy
                'cam': ['C6_term'], # Dummy
                'wrf': ['C6_term'], # Dummy
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['C1_Skw_fnc'],
                'sam': ['C1_SKW_FNC'], # Dummy
                'coamps': ['C1_Skw_fnc'], # Dummy
                'r408': ['C1_Skw_fnc'], # Dummy
                'hoc': ['C1_Skw_fnc'], # Dummy
                'e3sm': ['C1_Skw_fnc'], # Dummy
                'cam': ['C1_Skw_fnc'], # Dummy
                'wrf': ['C1_Skw_fnc'], # Dummy
                },
                'title': 'C1_Skw_fnc',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['C6rt_Skw_fnc'],
                'sam': ['C6RT_SKW_FNC'], # Dummy
                'coamps': ['C6rt_Skw_fnc'], # Dummy
                'r408': ['C6rt_Skw_fnc'], # Dummy
                'hoc': ['C6rt_Skw_fnc'], # Dummy
                'e3sm': ['C6rt_Skw_fnc'], # Dummy
                'cam': ['C6rt_Skw_fnc'], # Dummy
                'wrf': ['C6rt_Skw_fnc'], # Dummy
                },
                'title': 'C6rt_Skw_fnc',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['C6thl_Skw_fnc'],
                'sam': ['C6THL_SKW_FNC'], # Dummy
                'coamps': ['C6thl_Skw_fnc'], # Dummy
                'r408': ['C6thl_Skw_fnc'], # Dummy
                'hoc': ['C6thl_Skw_fnc'], # Dummy
                'e3sm': ['C6thl_Skw_fnc'], # Dummy
                'cam': ['C6thl_Skw_fnc'], # Dummy
                'wrf': ['C6thl_Skw_fnc'], # Dummy
                },
                'title': 'C6thl_Skw_fnc',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['C7_Skw_fnc'],
                'sam': ['C7_SKW_FNC'], # Dummy
                'coamps': ['C7_Skw_fnc'], # Dummy
                'r408': ['C7_Skw_fnc'], # Dummy
                'hoc': ['C7_Skw_fnc'], # Dummy
                'e3sm': ['C7_Skw_fnc'], # Dummy
                'cam': ['C7_Skw_fnc'], # Dummy
                'wrf': ['C7_Skw_fnc'], # Dummy
                },
                'title': 'C7_Skw_fnc',
                'sci_scale': 0,
            },

        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset,wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)
