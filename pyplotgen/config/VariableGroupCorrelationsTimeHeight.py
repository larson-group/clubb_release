"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup
from src.Panel import Panel


class VariableGroupCorrelationsTimeHeight(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None, r408_dataset=None,
                 wrf_benchmark_dataset=None,
                 hoc_dataset=None, cam_datasets=None, e3sm_datasets=None, sam_datasets=None, wrf_datasets=None,
                 priority_vars=False, background_rcm=False, background_rcm_folder=None):
        self.name = "corr time-height variables"

        self.variable_definitions = [
        {'var_names': # Correlation of chi and eta, 1st component
            {
                'clubb': ['corr_chi_eta_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of chi and eta, 1st component",
            'axis_title': "corr_chi_eta_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and eta, 2nd component
            {
                'clubb': ['corr_chi_eta_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of chi and eta, 2nd component",
            'axis_title': "corr_chi_eta_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and eta (ca), 1st component
            {
                'clubb': ['corr_chi_eta_1_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of chi and eta (ca), 1st component",
            'axis_title': "corr_chi_eta_1_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and eta (ca), 2nd component
            {
                'clubb': ['corr_chi_eta_2_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of chi and eta (ca), 2nd component",
            'axis_title': "corr_chi_eta_2_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and Ncn, 1st component
            {
                'clubb': ['corr_chi_Ncn_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Ncn, 1st component",
            'axis_title': "corr_chi_Ncn_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and Ncn, 2nd component
            {
                'clubb': ['corr_chi_Ncn_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Ncn, 2nd component",
            'axis_title': "corr_chi_Ncn_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of eta and Ncn, 1st component
            {
                'clubb': ['corr_eta_Ncn_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of eta and N_cn, 1st component",
            'axis_title': "corr_eta_Ncn_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of chi and Ncn, 2nd component
            {
                'clubb': ['corr_eta_Ncn_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of eta and N_cn, 2nd component",
            'axis_title': "corr_eta_Ncn_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of rt and thl, 1st component
            {
                'clubb': ['corr_rt_thl_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of rt and thl, 1st component",
            'axis_title': "corr_rt_thl_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of rt and thl, 2nd component
            {
                'clubb': ['corr_rt_thl_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of rt and thl, 2nd component",
            'axis_title': "corr_rt_thl_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and chi, 1st component
            {
                'clubb': ['corr_w_chi_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of w and chi, 1st component",
            'axis_title': "corr_w_chi_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and chi, 2nd component
            {
                'clubb': ['corr_w_chi_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            # 'lines': corr_w_chi_i_lines,
            'title': "Correlation of w and chi, 2nd component",
            'axis_title': "corr_w_chi_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and chi (ca), 1st component
            {
                'clubb': ['corr_w_chi_1_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of w and chi (ca), 1st component",
            'axis_title': "corr_w_chi_1_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and chi (ca), 2nd component
            {
                'clubb': ['corr_w_chi_2_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation of w and chi (ca), 2nd component",
            'axis_title': "corr_w_chi_2_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and Nr, 1st component
            {
                'clubb': ['corr_w_Nr_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': ['corr_w_Nr_1'],
                'hoc': ['corr_w_Nr_1'],
                'e3sm': ['corr_w_Nr_1'],
                'cam': ['corr_w_Nr_1'],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Nr, 1st component",
            'axis_title': "corr_w_Nr_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and Nr, 2nd component
            {
                'clubb': ['corr_w_Nr_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': ['corr_w_Nr_2'],
                'hoc': ['corr_w_Nr_2'],
                'e3sm': ['corr_w_Nr_2'],
                'cam': ['corr_w_Nr_2'],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Nr, 2nd component",
            'axis_title': "corr_w_Nr_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and Ncn, 1st component
            {
                'clubb': ['corr_w_Ncn_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': ['corr_w_Ncn_1'],
                'hoc': ['corr_w_Ncn_1'],
                'e3sm': ['corr_w_Ncn_1'],
                'cam': ['corr_w_Ncn_1'],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Ncn, 1st component",
            'axis_title': "corr_w_Ncn_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and Ncn, 2nd component
            {
                'clubb': ['corr_w_Ncn_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': ['corr_w_Ncn_2'],
                'hoc': ['corr_w_Ncn_2'],
                'e3sm': ['corr_w_Ncn_2'],
                'cam': ['corr_w_Ncn_2'],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and Ncn, 2nd component",
            'axis_title': "corr_w_Ncn_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and eta, 1st component
            {
                'clubb': ['corr_w_eta_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and eta, 1st component",
            'axis_title': "corr_w_eta_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and eta, 2nd component
            {
                'clubb': ['corr_w_eta_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and eta, 2nd component",
            'axis_title': "corr_w_eta_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and eta (ca), 1st component
            {
                'clubb': ['corr_w_eta_1_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and eta (ca), 1st component",
            'axis_title': "corr_w_eta_1_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and eta (ca), 2nd component
            {
                'clubb': ['corr_w_eta_2_ca'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and eta (ca), 2nd component",
            'axis_title': "corr_w_eta_2_ca [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and rt, 1st component
            {
                'clubb': ['corr_w_rt_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and rt, 1st component",
            'axis_title': "corr_w_rt_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and rt, 2nd component
            {
                'clubb': ['corr_w_rt_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and rt, 2nd component",
            'axis_title': "corr_w_rt_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and thl, 1st component
            {
                'clubb': ['corr_w_thl_1'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and thl, 1st component",
            'axis_title': "corr_w_thl_1 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        {'var_names': # Correlation of w and thl, 2nd component
            {
                'clubb': ['corr_w_thl_2'],
                'sam': [], # no quotes here so no SAM lines plotted
                'coamps': [], # no quotes here so no COAMPS lines plotted
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                'wrf': [''],
            },
            'title': "Correlation (in-precip) of w and thl, 2nd component",
            'axis_title': "corr_w_thl_2 [$-$]",
            'type': Panel.TYPE_TIMEHEIGHT,
            # 'sci_scale': 0,
        },
        # Correlations below were not found in netcdf output
        # Therefore, they are commented out for now
            # {'var_names':
            #     {
            #     'clubb': ['corr_w_rr_1'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_w_rr_1'],
            #     'hoc': ['corr_w_rr_1'],
            #     'e3sm': ['corr_w_rr_1'],
            #     'cam': ['corr_w_rr_1'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of w and rr, 1st component",
            #     'axis_title': "corr_w_rr_1 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_w_rr_2'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_w_rr_2'],
            #     'hoc': ['corr_w_rr_2'],
            #     'e3sm': ['corr_w_rr_2'],
            #     'cam': ['corr_w_rr_2'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of w and rr. 2nd component",
            #     'axis_title': "corr_w_rr_2 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_chi_rr_1'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_chi_rr_1'],
            #     'hoc': ['corr_chi_rr_1'],
            #     'e3sm': ['corr_chi_rr_1'],
            #     'cam': ['corr_chi_rr_1'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of chi and rr, 1st component",
            #     'axis_title': "corr_chi_rr_1 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_chi_rr_2'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_chi_rr_2'],
            #     'hoc': ['corr_chi_rr_2'],
            #     'e3sm': ['corr_chi_rr_2'],
            #     'cam': ['corr_chi_rr_2'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of chi and rr, 2nd component",
            #     'axis_title': "corr_chi_rr_2 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_chi_Nr_1'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_chi_Nr_1'],
            #     'hoc': ['corr_chi_Nr_1'],
            #     'e3sm': ['corr_chi_Nr_1'],
            #     'cam': ['corr_chi_Nr_1'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of chi and Nr, 1st component",
            #     'axis_title': "corr_chi_Nr_1 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_chi_Nr_2'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_chi_Nr_2'],
            #     'hoc': ['corr_chi_Nr_2'],
            #     'e3sm': ['corr_chi_Nr_2'],
            #     'cam': ['corr_chi_Nr_2'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of chi and Nr, 2nd component",
            #     'axis_title': "corr_chi_Nr_2 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_rr_Nr_1'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_rr_Nr_1'],
            #     'hoc': ['corr_rr_Nr_1'],
            #     'e3sm': ['corr_rr_Nr_1'],
            #     'cam': ['corr_rr_Nr_1'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of rr and Nr, 1st component",
            #     'axis_title': "corr_rr_Nr_1 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
            # {'var_names':
            #     {
            #     'clubb': ['corr_rr_Nr_2'],
            #     'sam': [],
            #     'coamps': [],
            #     'r408': ['corr_rr_Nr_2'],
            #     'hoc': ['corr_rr_Nr_2'],
            #     'e3sm': ['corr_rr_Nr_2'],
            #     'cam': ['corr_rr_Nr_2'],
            #     'wrf': [],
            #     },
            #     'title': "Correlation (in-precip) of rr and Nr, 2nd component",
            #     'axis_title': "corr_rr_Nr_2 [-]"
            #     'type': Panel.TYPE_TIMEHEIGHT,
            # },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset, coamps_benchmark_dataset=coamps_benchmark_dataset,
                         wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)