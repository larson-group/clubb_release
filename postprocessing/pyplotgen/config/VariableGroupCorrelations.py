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
                 priority_vars=False, background_rcm=False, background_rcm_folder=None):
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
        corr_rt_thl_i_lines = [
            {'var_names': ['corr_rt_thl_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_rt_thl_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_w_rt_i_lines = [
            {'var_names': ['corr_w_rt_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_rt_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_w_thl_i_lines = [
            {'var_names': ['corr_w_thl_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_thl_2'], 'legend_label': 'PDF comp. 2'},
        ]
        stdev_chi_i_lines = [
            {'var_names': ['stdev_chi_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['stdev_chi_2'], 'legend_label': 'PDF comp. 2'},
        ]
        stdev_eta_i_lines = [
            {'var_names': ['stdev_eta_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['stdev_eta_2'], 'legend_label': 'PDF comp. 2'},
        ]
        corr_chi_eta_i_lines = [
            {'var_names': ['corr_chi_eta_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_chi_eta_2'], 'legend_label': 'PDF comp. 2'},
        ]

        self.variable_definitions = [
            {'var_names': # Correlation of chi and eta (multiple lines)
                {
                    'clubb': [''],
                    'sam': [], # no quotes here so no SAM lines plotted
                    'coamps': [], # no quotes here so no COAMPS lines plotted
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'lines': corr_chi_eta_i_lines,
                'title': "Correlation of chi and eta",
                'axis_title': "corr_chi_eta_i [$-$]",
                'sci_scale': 0,
            },
            {'var_names': # Correlation of rt and thl
                {
                    'clubb': ['corr_rt_thl_i'],
                    'sam': [], # no quotes here so no SAM lines plotted
                    'coamps': [], # no quotes here so no COAMPS lines plotted
                    'r408': ['corr_rt_thl_i'],
                    'hoc': ['corr_rt_thl_i'],
                    'e3sm': ['corr_rt_thl_i'],
                    'cam': ['corr_rt_thl_i'],
                    'wrf': ['corr_rt_thl_i'],
                },
                'lines': corr_rt_thl_i_lines,
                'title': "Correlation of rt and thl",
                'axis_title': "corr_rt_thl_i [$-$]",
                # 'sci_scale': 0,
            },
            {'var_names': # Correlation of w and rt
                {
                    'clubb': ['corr_w_rt_i'],
                    'sam': [], # no quotes here so no SAM lines plotted
                    'coamps': [], # no quotes here so no COAMPS lines plotted
                    'r408': ['corr_w_rt_i'],
                    'hoc': ['corr_w_rt_i'],
                    'e3sm': ['corr_w_rt_i'],
                    'cam': ['corr_w_rt_i'],
                    'wrf': ['corr_w_rt_i'],
                },
                'lines': corr_w_rt_i_lines,
                'title': "Correlation (in-precip) of w and rt",
                'axis_title': "corr_w_rt_i [$-$]",
                # 'sci_scale': 0,
            },
            {'var_names': # Correlation of w and thl
                {
                    'clubb': ['corr_w_thl_i'],
                    'sam': [], # no quotes here so no SAM lines plotted
                    'coamps': [], # no quotes here so no COAMPS lines plotted
                    'r408': ['corr_w_thl_i'],
                    'hoc': ['corr_w_thl_i'],
                    'e3sm': ['corr_w_thl_i'],
                    'cam': ['corr_w_thl_i'],
                    'wrf': ['corr_w_thl_i'],
                },
                'lines': corr_w_thl_i_lines,
                'title': "Correlation (in-precip) of w and thl",
                'axis_title': "corr_w_thl_i [$-$]",
                # 'sci_scale': 0,
            },
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
            {'var_names':
                {
                'clubb': ['stdev_chi_i'],
                'sam': [],
                'coamps': [],
                'r408': ['stdev_chi_i'],
                'hoc': ['stdev_chi_i'],
                'e3sm': ['stdev_chi_i'],
                'cam': ['stdev_chi_i'],
                'wrf': [],
                },
                'lines': stdev_chi_i_lines,
                'title': "Std dev chi",
                'axis_title': "stdev_chi_i [-]"
            },
            {'var_names':
                {
                'clubb': ['stdev_eta_i'],
                'sam': [],
                'coamps': [],
                'r408': ['stdev_eta_i'],
                'hoc': ['stdev_eta_i'],
                'e3sm': ['stdev_eta_i'],
                'cam': ['stdev_eta_i'],
                'wrf': [],
                },
                'lines': stdev_eta_i_lines,
                'title': "Std dev eta",
                'axis_title': "stdev_eta_i [-]"
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset, coamps_benchmark_dataset=coamps_benchmark_dataset,
                         wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)