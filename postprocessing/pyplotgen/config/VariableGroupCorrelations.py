"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupCorrelations(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None, r408_dataset=None,
                 wrf_benchmark_dataset=None, 
                 hoc_dataset=None, cam_datasets=None, e3sm_datasets=None, sam_datasets=None, wrf_datasets=None,
                 priority_vars=False, background_rcm=False):
        self.name = "corr variables"

        corr_w_rr_i_lines = [
            {'var_names': ['corr_w_rr_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_rr_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_w_Nr_i_lines = [
            {'var_names': ['corr_w_Nr_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_Nr_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_w_Ncn_i_lines = [
            {'var_names': ['corr_w_Ncn_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_Ncn_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_chi_rr_i_lines = [
            {'var_names': ['corr_chi_rr_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_chi_rr_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_chi_Nr_i_lines = [
            {'var_names': ['corr_chi_Nr_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_chi_Nr_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_chi_Ncn_i_lines = [
            {'var_names': ['corr_chi_Ncn_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_chi_Ncn_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_rr_Nr_i_lines = [
            {'var_names': ['corr_rr_Nr_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_rr_Nr_2'], 'legend_label': 'PDF comp. 2'},
        ]

        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['corr_w_rr_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_w_rr_i'],
                'hoc': ['corr_w_rr_i'],
                'e3sm': ['corr_w_rr_i'],
                'cam': ['corr_w_rr_i'],
                'wrf': [],
                },
                'lines': corr_w_rr_i_lines,
                'title': "Correlation (in-precip) of w and rr",
                'axis_title': "corr_w_rr_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_w_Nr_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_w_Nr_i'],
                'hoc': ['corr_w_Nr_i'],
                'e3sm': ['corr_w_Nr_i'],
                'cam': ['corr_w_Nr_i'],
                'wrf': [],
                },
                'lines': corr_w_Nr_i_lines,
                'title': "Correlation (in-precip) of w and Nr",
                'axis_title': "corr_w_Nr_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_w_Ncn_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_w_Ncn_i'],
                'hoc': ['corr_w_Ncn_i'],
                'e3sm': ['corr_w_Ncn_i'],
                'cam': ['corr_w_Ncn_i'],
                'wrf': [],
                },
                'lines': corr_w_Ncn_i_lines,
                'title': "Correlation (in-precip) of w and Ncn",
                'axis_title': "corr_w_Ncn_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_chi_rr_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_chi_rr_i'],
                'hoc': ['corr_chi_rr_i'],
                'e3sm': ['corr_chi_rr_i'],
                'cam': ['corr_chi_rr_i'],
                'wrf': [],
                },
                'lines': corr_chi_rr_i_lines,
                'title': "Correlation (in-precip) of chi and rr",
                'axis_title': "corr_chi_rr_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_chi_Nr_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_chi_Nr_i'],
                'hoc': ['corr_chi_Nr_i'],
                'e3sm': ['corr_chi_Nr_i'],
                'cam': ['corr_chi_Nr_i'],
                'wrf': [],
                },
                'lines': corr_chi_Nr_i_lines,
                'title': "Correlation (in-precip) of chi and Nr",
                'axis_title': "corr_chi_Nr_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_chi_Ncn_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_chi_Ncn_i'],
                'hoc': ['corr_chi_Ncn_i'],
                'e3sm': ['corr_chi_Ncn_i'],
                'cam': ['corr_chi_Ncn_i'],
                'wrf': [],
                },
                'lines': corr_chi_Ncn_i_lines,
                'title': "Correlation (in-precip) of chi and Ncn",
                'axis_title': "corr_chi_Ncn_i [-]"
            },
            {'var_names':
                {
                'clubb': ['corr_rr_Nr_i'],
                'sam': [],
                'coamps': [],
                'r408': ['corr_rr_Nr_i'],
                'hoc': ['corr_rr_Nr_i'],
                'e3sm': ['corr_rr_Nr_i'],
                'cam': ['corr_rr_Nr_i'],
                'wrf': [],
                },
                'lines': corr_rr_Nr_i_lines,
                'title': "Correlation (in-precip) of rr and Nr",
                'axis_title': "corr_rr_Nr_i [-]"
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset, coamps_benchmark_dataset=coamps_benchmark_dataset,
                         wrf_benchmark_dataset=wrf_benchmark_dataset, 
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm)
