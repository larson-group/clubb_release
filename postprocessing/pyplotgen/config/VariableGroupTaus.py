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

        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset,wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)
