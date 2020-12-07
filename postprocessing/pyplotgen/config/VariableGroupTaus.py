"""
:author: Ben Stephens
:date: December 2020
"""

from src.VariableGroup import VariableGroup


class VariableGroupTaus(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False):
        """
        """
        self.name = "tau variables"
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['bv_freq_sqd'],
                'sam': ['BV_FREQ_SQD'],
                'coamps': ['bv_freq_sqd'],
                'r408': ['bv_freq_sqd'],
                'hoc': ['bv_freq_sqd'],
                'e3sm': ['bv_freq_sqd'],
                'cam': ['bv_freq_sqd'],
                'wrf': ['bv_freq_sqd'],
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
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars)
