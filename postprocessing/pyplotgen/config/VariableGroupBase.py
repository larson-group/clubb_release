"""
:author: Nicolas Strike
:date: Mid 2019
"""
import numpy as np
from netCDF4 import Dataset

from src.Panel import Panel
from src.VariableGroup import VariableGroup

class VariableGroupBase(VariableGroup):
    """
    This is a panel group used for testing the functionality of pyplotgen.
    It contains a set of common panels being used for representing the majority
    of cases.
    """

    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False,
                 background_rcm=False, background_rcm_folder=None):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "base variables"

        corr_w_chi_i_lines = [
            {'var_names': ['corr_w_chi_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_w_chi_2'], 'legend_label': 'PDF comp. 2'},
        ]

        norm_grid_density = [
            {'var_names': ['norm_min_grid_dens'], 'legend_label': 'norm. min. grid dens.'},
            {'var_names': ['norm_grid_dens'], 'legend_label': 'norm. grid dens.'},
        ]

        norm_grid_density_plot = [
            {'var_names': ['norm_min_grid_dens'], 'legend_label': r"$g_\mathrm{min}(z)$"},
            {'var_names': ['norm_grid_dens'], 'legend_label': r"$g(z)$"},
        ]

        each_ref_crit_term = [
            {'var_names': ['alt_term'], 'legend_label': r"$3 \cdot g_z (z)$"},
            {'var_names': ['chi_term_time_avg'], 'legend_label': r"$100 \cdot \bar{g}_\chi (z)$"},
            {'var_names': ['brunt_term_time_avg'], 'legend_label': r"$7.5 \cdot \bar{g}_\mathrm{Ri} (z)$"},
        ]

        corr_chi_eta_i_lines = [
            {'var_names': ['corr_chi_eta_1'], 'legend_label': 'PDF comp. 1'},
            {'var_names': ['corr_chi_eta_2'], 'legend_label': 'PDF comp. 2'},
        ]

        tau_i_lines = [
            {'var_names': ['invrs_tau_zm'], 'legend_label': 'invrs_tau_zm'},
            {'var_names': ['invrs_tau_xp2_zm'], 'legend_label': 'invrs_tau_xp2_zm'},
            {'var_names': ['invrs_tau_wp2_zm'], 'legend_label': 'invrs_tau_wp2_zm'},
            {'var_names': ['invrs_tau_wpxp_zm'], 'legend_label': 'invrs_tau_wpxp_zm'},
            {'var_names': ['invrs_tau_wp3_zm'], 'legend_label': 'invrs_tau_wp3_zm'},
        ]

        tau_j_lines = [
            {'var_names': ['invrs_tau_no_N2_zm'], 'legend_label': 'invrs_tau_no_N2_zm'},
            {'var_names': ['invrs_tau_bkgnd'], 'legend_label': 'invrs_tau_bkgnd'},
            {'var_names': ['invrs_tau_sfc'], 'legend_label': 'invrs_tau_sfc'},
            {'var_names': ['invrs_tau_shear'], 'legend_label': 'invrs_tau_shear'},
        ]

        self.variable_definitions = [
            {'var_names':
                {
                    'clubb': ['thlm'],
                    'sam': [self.getThlmSamCalc, 'THETAL', 'THETA'],
                    'coamps': ['thlm'],
                    'r408': ['thlm'],
                    'hoc': ['thlm'],
                    'e3sm': ['thlm'],
                    'cam': ['thlm'],
                    'wrf': ['thlm','CSP_THL'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [self.getRtmSamCalc,'QT'],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [self.getRtmSamCalc,'QT'],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
                'title': 'GABLS2', # TODO change title here
                'axis_title': r"$\bar{r}_t ~~~[\mathrm{kg}/\mathrm{kg}]$",
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [self.getWpthlpSamCalc, 'WPTHLP'],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [self.getWpthlpSamCalc, 'WPTHLP'],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
                'title': 'ASTEX', # TODO change title here
                'axis_title': r"$\overline{w^{\prime}\theta^{\prime}_l} ~~~[\mathrm{m}\mathrm{s}^{-1}\mathrm{K}]$",
            },
            {'var_names':
                {
                    'clubb': ['wprtp'],
                    'sam': [self.getWprtpSamCalc, 'WPRTP'],
                    'coamps': ['wpqtp', 'wprtp'],
                    'r408': ['wprtp'],
                    'hoc': ['wprtp'],
                    'e3sm': ['wprtp'],
                    'cam': ['WPRTP_clubb', 'wprtp'],
                    'wrf': ['wprtp','CSP_WQT'],
                },
                'sci_scale': -4,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['cloud_frac'],
                    'sam': ['CLD_FRAC_CLUBB', 'CLD'],
                    'coamps': ['cf'],
                    'r408': ['cloud_frac', 'cf'],
                    'hoc': ['cloud_frac', 'cf'],
                    'e3sm': ['cloud_frac'],
                    'cam': ['CLOUD', 'cloud_frac'],
                    'wrf': ['cloud_frac'],
                },
                'title': 'Cloud fraction',
                'axis_title': 'cloud_frac [$-$]',
                #'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['cloud_frac'],
                    'sam': ['CLD_FRAC_CLUBB', 'CLD'],
                    'coamps': ['cf'],
                    'r408': ['cloud_frac', 'cf'],
                    'hoc': ['cloud_frac', 'cf'],
                    'e3sm': ['cloud_frac'],
                    'cam': ['CLOUD', 'cloud_frac'],
                    'wrf': ['cloud_frac'],
                },
                'title': 'ASTEX', # TODO change title here
                'axis_title': 'Cloud fraction    [$-$]',
                #'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': ['QCL', 'QC'],
                    'coamps': ['qcm', 'rcm'],
                    'r408': ['rcm'],
                    'hoc': ['rcm'],
                    'e3sm': ['rcm'],
                    'cam': ['CLDLIQ', 'rcm'],
                    'wrf': ['rcm','CSP_QC'],
                },
                'sam_conv_factor': 1 / 1000,
                #'sci_scale': -5,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': ['QCL', 'QC'],
                    'coamps': ['qcm', 'rcm'],
                    'r408': ['rcm'],
                    'hoc': ['rcm'],
                    'e3sm': ['rcm'],
                    'cam': ['CLDLIQ', 'rcm'],
                    'wrf': ['rcm','CSP_QC'],
                },
                'sam_conv_factor': 1 / 1000,
                #'sci_scale': -5,
                'priority': True,
                'title': 'ASTEX', # TODO change title here
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [self.get_wp2_sam_calc, 'W2', 'WP2'],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [self.get_wp2_sam_calc, 'W2', 'WP2'],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
                'title': 'ASTEX', # TODO change title here 
                'axis_title': r"$\overline{{w^{\prime}}^2} ~~~[\mathrm{m}^2/\mathrm{s}^2]$",
            },
            {'var_names':
                {
                    'clubb': ['wp3'],
                    'sam': [self.get_wp3_sam_calc, 'wp3', 'W3', 'WP3'],
                    'coamps': ['wp3', 'W3', 'WP3'],
                    'r408': ['wp3'],
                    'hoc': ['wp3'],
                    'e3sm': ['wp3'],
                    'cam': ['WP3_CLUBB', 'wp3'],
                    'wrf': ['wp3','CSP_W3'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['thlp2'],
                    'sam': [self.get_thlp2_sam_calc, 'THLP2', 'TL2'],
                    'coamps': ['thlp2'],
                    'r408': ['thlp2'],
                    'hoc': ['thlp2'],
                    'e3sm': ['thlp2'],
                    'cam': ['THLP2_CLUBB', 'thlp2'],
                    'wrf': ['thlp2','CSP_THL2'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtp2'],
                    'sam': [self.getRtp2SamCalc, 'RTP2', 'QT2'],
                    'coamps': ['qtp2'],
                    'r408': ['rtp2'],
                    'hoc': ['rtp2'],
                    'e3sm': ['rtp2'],
                    'cam': ['RTP2_CLUBB', 'rtp2'],
                    'wrf': ['rtp2'],
                },
                'sci_scale': -7,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtpthlp'],
                    'sam': ['RTPTHLP_SGS', 'RTPTHLP', 'TQ'],
                    'coamps': ['qtpthlp', 'rtpthlp'],
                    'r408': ['rtpthlp'],
                    'hoc': ['rtpthlp'],
                    'e3sm': ['rtpthlp'],
                    'cam': ['rtpthlp'],
                    'wrf': ['rtpthlp'],
                },
                'sci_scale': -4,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtp3'],
                    'sam': [ 'RTP3', self.getRtp3SamCalc],
                    'coamps': ['qtp3', 'rtp3'],
                    'r408': ['rtp3'],
                    'hoc': ['rtp3'],
                    'e3sm': ['rtp3'],
                    'cam': ['rtp3'],
                    'wrf': ['rtp3'],
                },
                'sci_scale': -9,
            },
            {'var_names':
                {
                    'clubb': ['thlp3'],
                    'sam': ['THLP3'],
                    'coamps': ['thlp3'],
                    'r408': ['thlp3'],
                    'hoc': ['thlp3'],
                    'e3sm': ['thlp3'],
                    'cam': ['thlp3'],
                    'wrf': ['thlp3'],
                },
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['Skw_zt'],
                    'sam': [self.getSkwZtLesCalc,'Skw_zt'],
                    'coamps': [self.getSkwZtLesCalc, 'Skw_zt'],
                    'r408': ['Skw_zt'],
                    'hoc': ['Skw_zt'],
                    'e3sm': ['Skw_zt'],
                    'cam': ['Skw_zt'],
                    'wrf': ['Skw_zt','CSP_WSKEW'],
                },
                'title': 'Skw, Skewness of vertical velocity w',
                'axis_title': 'Skw [$-$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['Skrt_zt'],
                    'sam': [self.getSkrtZtLesCalc,'Skrt_zt'],
                    'coamps': [self.getSkrtZtLesCalc, 'Skrt_zt'],
                    'r408': ['Skrt_zt'],
                    'hoc': ['Skrt_zt'],
                    'e3sm': ['Skrt_zt'],
                    'cam': ['Skrt_zt'],
                    'wrf': ['Skrt_zt'],
                },
                'title': 'Skrt, Skewness of total liquid water rt',
                'axis_title': 'Skrt [$-$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['Skthl_zt'],
                    'sam': [self.getSkthlZtLesCalc,'Skthl_zt'],
                    'coamps': [self.getSkthlZtLesCalc, 'Skthl_zt'],
                    'r408': ['Skthl_zt'],
                    'hoc': ['Skthl_zt'],
                    'e3sm': ['Skthl_zt'],
                    'cam': ['Skthl_zt'],
                    'wrf': ['Skthl_zt'],
                },
                'title': 'Skthl, Skewness of liq. water pot. temp. thl',
                'axis_title': 'Skthl [$-$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['wm', 'wlsm'],
                    'sam': ['WOBS', 'WM'],
                    'coamps': ['wlsm', 'wm'],
                    'r408': ['wm'],
                    'hoc': ['wm'],
                    'e3sm': ['wm'],
                    'cam': ['wm'], # -OMEGA /(9.81.*1)
                    'wrf': ['wm', 'wlsm','CSP_W'],
                },
                'sci_scale': -4,
            },
            {'var_names':
                {
                    'clubb': ['um'],
                    'sam': ['U'],
                    'coamps': ['um'],
                    'r408': ['um'],
                    'hoc': ['um'],
                    'e3sm': ['um'],
                    'cam': ['U','um'],
                    'wrf': ['um','CSP_U'],
                },
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['vm'],
                    'sam': ['V'],
                    'coamps': ['vm'],
                    'r408': ['vm'],
                    'hoc': ['vm'],
                    'e3sm': ['vm'],
                    'cam': ['V', 'vm'],
                    'wrf': ['vm','CSP_V'],
                },
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['upwp'],
                    'sam': [self.get_upwp_sam_calc, 'UW'],
                    'coamps': [self.getUwCoampsData, 'upwp'],
                    'r408': ['upwp'],
                    'hoc': ['upwp'],
                    'e3sm': ['upwp'],
                    'cam': ['UPWP_CLUBB', 'upwp'],
                    'wrf': ['upwp','CSP_UW'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['vpwp'],
                    'sam': [self.get_vpwp_sam_calc, 'VW'],
                    'coamps': [self.getVwCoampsData, 'vpwp'],
                    'r408': ['vpwp'],
                    'hoc': ['vpwp'],
                    'e3sm': ['vpwp'],
                    'cam': ['VPWP_CLUBB', 'vpwp'],
                    'wrf': ['vpwp','CSP_VW'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['up2'],
                    'sam': [self.get_up2_sam_calc, 'U2'],
                    'coamps': ['up2'],
                    'r408': ['up2'],
                    'hoc': ['up2'],
                    'e3sm': ['up2'],
                    'cam': ['UU', 'up2'],
                    'wrf': ['up2','CSP_U2'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['vp2'],
                    'sam': [self.get_vp2_sam_calc, 'V2'],
                    'coamps': ['vp2'],
                    'r408': ['vp2'],
                    'hoc': ['vp2'],
                    'e3sm': ['vp2'],
                    'cam': ['VV', 'vp2'],
                    'wrf': ['vp2','CSP_V2'],
                },
                'sci_scale': 0,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rcp2'],
                    'sam': ['QC2'],
                    'coamps': ['qcp2', 'rcp2', 'rlp2'],
                    'r408': ['rcp2'],
                    'hoc': ['rcp2'],
                    'e3sm': ['rcp2'],
                    'cam': ['rcp2'],
                    'wrf': ['rcp2'],
                },
                'sam_conv_factor': 1 / 10 ** 6,
                'sci_scale': -8,
            },
            {'var_names':
                {
                    'clubb': ['lwp'],
                    'sam': ['CWP'],
                    'coamps': ['lwp'],
                    'r408': ['lwp'],
                    'hoc': ['lwp'],
                    'e3sm': ['lwp'],
                    'cam': ['TGCLDLWP', 'lwp'],
                    'wrf': ['lwp','CST_LWP'],
                },
                'type': Panel.TYPE_TIMESERIES,
                'sam_conv_factor': 1 / 1000,
            },
            {'var_names':
                {
                    'clubb': ['thlm_vert_avg'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_TIMESERIES
            },
            {'var_names':
                {
                    'clubb': ['rtm_vert_avg'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_TIMESERIES
            },
            {'var_names':
                {
                    'clubb': ['wp2_vert_avg'],
                    'sam': ['W2_VERT_AVG'],
                    'coamps': ['wp2_vert_avg'],
                    'r408': ['wp2_vert_avg'],
                    'hoc': ['wp2_vert_avg'],
                    'e3sm': ['wp2_vert_avg'],
                    'cam': ['wp2_vert_avg'],
                    'wrf': ['wp2_vert_avg'],
                },
                'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names':
                {
                    'clubb': ['rtp2_vert_avg'],
                    'sam': [''],
                    'coamps': ['rtp2_vert_avg'],
                    'r408': ['rtp2_vert_avg'],
                    'hoc': ['rtp2_vert_avg'],
                    'e3sm': ['rtp2_vert_avg'],
                    'cam': ['rtp2_vert_avg'],
                    'wrf': ['rtp2_vert_avg'],
                },
                'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names':
                {
                    'clubb': ['thlp2_vert_avg'],
                    'sam': [''],
                    'coamps': ['thlp2_vert_avg'],
                    'r408': ['thlp2_vert_avg'],
                    'hoc': ['thlp2_vert_avg'],
                    'e3sm': ['thlp2_vert_avg'],
                    'cam': ['thlp2_vert_avg'],
                    'wrf': ['thlp2_vert_avg'],
                },
                'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names':
                {
                    'clubb': ['wpthlp_sfc'],
                    'sam': [''],
                    'coamps': ['wpthlp_sfc'],
                    'r408': ['wpthlp_sfc'],
                    'hoc': ['wpthlp_sfc'],
                    'e3sm': ['wpthlp_sfc'],
                    'cam': ['wpthlp_sfc'],
                    'wrf': ['wpthlp_sfc'],
                },
                'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names':
                {
                    'clubb': ['wprtp_sfc'],
                    'sam': [''],
                    'coamps': ['wprtp_sfc'],
                    'r408': ['wprtp_sfc'],
                    'hoc': ['wprtp_sfc'],
                    'e3sm': ['wprtp_sfc'],
                    'cam': ['wprtp_sfc'],
                    'wrf': ['wprtp_sfc'],
                },
                'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names': # First set of inverse tau variables (multiple lines)
                {
                    'clubb': [''],
                    'sam': [],  # no quotes here so no SAM lines plotted
                    'coamps': [], # no quotes here so no COAMPS lines plotted
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'lines': tau_i_lines,
                'title': 'CLUBB inverse time-scale, 1/tau (panel 1)',
                'axis_title': r'inverse tau [s$^{-1}$]',
                'sci_scale': 0,
            },
            {'var_names': # Second set of inverse tau variables (multiple lines)
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
                'lines': tau_j_lines,
                'title': 'CLUBB inverse time-scale, 1/tau (panel 2)',
                'axis_title': r'inverse tau [s$^{-1}$]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['Lscale'],
                    'sam': [],
                    'coamps': [],
                    'r408': ['Lscale'],
                    'hoc': ['Lscale'],
                    'e3sm': ['Lscale'],
                    'cam': ['Lscale'],
                    'wrf': ['Lscale'],
                },
                'title': 'CLUBB turbulent mixing length',
                'axis_title': 'Lscale [m]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['bv_freq_sqd'],
                'sam': [self.getBVSqdSamCalc, 'BV_FREQ_SQD'],
                'coamps': [],
                'r408': ['bv_freq_sqd'],
                'hoc': ['bv_freq_sqd'],
                'e3sm': ['bv_freq_sqd'],
                'cam': ['bv_freq_sqd'],
                'wrf': ['bv_freq_sqd'],
                },
                'title': 'Brunt-Vaisala frequency squared',
                'axis_title': 'bv_freq_sqd [$\mathrm{1/s^2}$]',
                #'sci_scale': 0,
                'sam_calc': self.getBVSqdSamCalc,
            },
            {'var_names':
                {
                'clubb': ['chi'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'Extended liquid water content',
                'axis_title': 'chi [$kg/kg$]',
            },
            {'var_names':
                {
                'clubb': ['ddzt_umvm_sqd'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'Square sum of horizontal shears, ddzt_umvm_sqd',
                'axis_title': '',
            },
            {'var_names':
                {
                'clubb': ['grid_density'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'Proposed non-normalized grid density',
                'axis_title': 'grid_density [1/m]',
            },
            {'var_names':
                {
                'clubb': [''],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'lines': norm_grid_density,
                'title': 'Proposed normalized grid density',
                'axis_title': 'norm_min_grid_dens [1/m]',
            },
            {'var_names':
                {
                'clubb': [''],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'lines': norm_grid_density_plot,
                'title': 'GABLS2', # TODO change title here
                'axis_title': 'grid density [1/m]',
            },
            {'var_names':
                {
                'clubb': [''],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'lines': each_ref_crit_term,
                'title': 'GABLS2', # TODO change title here
                'axis_title': 'grid density [1/m]',
            },
            {'var_names':
                {
                'clubb': ['alt_term'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': '$g_z(z)$',
                'axis_title': '$g_z$ [1/m]',
            },
            {'var_names':
                {
                'clubb': ['lscale_term'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'non-normalized inverse Lscale term of grid density',
                'axis_title': 'lscale_term [1/m]',
            },
            {'var_names':
                {
                'clubb': ['lscale_term_time_avg'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'non-normalized time avg. inverse Lscale term of grid density',
                'axis_title': 'lscale_term_time_avg [1/m]',
            },
            {'var_names':
                {
                'clubb': ['brunt_term_time_avg'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': '$\\bar{g}_\mathrm{Ri}(z)$',
                'axis_title': '$\\bar{g}_\mathrm{Ri}(z)$ [1/m]',
            },
            {'var_names':
                {
                'clubb': ['chi_term_time_avg'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': '$\\bar{g}_\chi(z)$',
                'axis_title': '$\\bar{g}_\chi(z)$ [1/m]',
            },
            {'var_names':
                {
                'clubb': ['chi_term'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'non-normalized chi term of grid density',
                'axis_title': 'chi_term [1/m]',
            },
            {'var_names':
                {
                'clubb': ['brunt_term'],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                'wrf': [],
                },
                'title': 'non-normalized brunt term of grid density',
                'axis_title': 'brunt_term [1/m]',
            },
            {'var_names':
                {
                'clubb': ['bv_freq_sqd_splat'],
                'sam': [],
                'coamps': [],
                'r408': ['bv_freq_sqd_splat'],
                'hoc': ['bv_freq_sqd_splat'],
                'e3sm': ['bv_freq_sqd_splat'],
                'cam': ['bv_freq_sqd_splat'],
                'wrf': ['bv_freq_sqd_splat'],
                },
                'title': 'Brunt-Vaisala frequency squared for splatting',
                'axis_title': 'bv_freq_sqd_splat [$\mathrm{1/s^2}$]',
                #'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['wpthvp'],
                    'sam': [self.getWpthvpSamCalc,'WPTHVP'],
                    'coamps': ['wpthvp'],
                    'r408': ['wpthvp'],
                    'hoc': ['wpthvp'],
                    'e3sm': ['wpthvp'],
                    'cam': ['WPTHVP_CLUBB', 'wpthvp'],
                    'wrf': ['wpthvp'],
                },
            },
            {'var_names':
                {
                    'clubb': ['radht'],
                    'sam': ['RADQR'],
                    'coamps': ['radht'],
                    'r408': ['radht'],
                    'hoc': ['radht'],
                    'e3sm': ['radht'],
                    'cam': ['radht'],
                    'wrf': ['radht'],
                },
                'sam_conv_factor': 1 / 86400,
            },
            {'var_names':
                {
                    'clubb': ['rtpthvp'],
                    'sam': ['RTPTHVP'],
                    'coamps': ['qtpthvp', 'rtpthvp'],
                    'r408': ['rtpthvp'],
                    'hoc': ['rtpthvp'],
                    'e3sm': ['rtpthvp'],
                    'cam': ['rtpthvp'],
                    'wrf': ['rtpthvp'],
                },
                'sci_scale': -5,
            },
            {'var_names': # Correlation of w and chi (multiple lines)
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
                'lines': corr_w_chi_i_lines,
                'title': "Correlation of w and chi",
                'axis_title': "corr_w_chi_i [$-$]",
                'sci_scale': 0,
            },
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
            {'var_names':
                {
                    'clubb': ['thlpthvp'],
                    'sam': ['thlpthvp', 'THLPTHVP'],
                    'coamps': ['thlpthvp'],
                    'r408': ['thlpthvp'],
                    'hoc': ['thlpthvp'],
                    'e3sm': ['thlpthvp'],
                    'cam': ['thlpthvp'],
                    'wrf': ['thlpthvp'],
                },
            },
            {'var_names':
                {
                    'clubb': [self.get_rc_coef_zm_X_wprcp_clubb_line, 'rc_coef_zm * wprcp'],
                    'sam': [self.get_rc_coef_zm_X_wprcp_sam_calc, 'rc_coef_zm * wprcp'],
                    'coamps': [self.get_rc_coef_zm_X_wprcp_coamps_calc, 'rc_coef_zm * wprcp'],
                    'r408': ['rc_coef_zm * wprcp'],
                    'hoc': ['rc_coef_zm * wprcp'],
                    'e3sm': ['rc_coef_zm * wprcp'],
                    'cam': ['rc_coef_zm * wprcp'],
                    'wrf': [self.get_rc_coef_zm_X_wprcp_wrf_line,'rc_coef_zm * wprcp'],
                },
                'title': 'Contribution of Cloud Water Flux to wpthvp',
                'axis_title': 'rc_coef_zm * wprcp [K m/s]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': [self.get_rc_coef_zm_X_thlprcp_clubb_calc, 'rc_coef_zm * thlprcp'],
                    'sam': [self.get_rc_coef_zm_X_thlprcp_sam_calc,'rc_coef_zm * thlprcp'],
                    'coamps': [self.get_rc_coef_zm_X_thlprcp_coamps_calc, 'rc_coef_zm * thlprcp'],
                    'r408': ['rc_coef_zm * thlprcp'],
                    'hoc': ['rc_coef_zm * thlprcp'],
                    'e3sm': ['rc_coef_zm * thlprcp'],
                    'cam': ['rc_coef_zm * thlprcp'],
                    'wrf': [self.get_rc_coef_zm_X_thlprcp_wrf_calc,'rc_coef_zm * thlprcp'],
                },
                'title': 'Contribution of Cloud Water Flux to thlpthvp',
                'axis_title': 'rc_coef_zm * thlprcp [K^2]',
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': [self.get_rc_coef_zm_X_rtprcp_clubb_calc, 'rc_coef_zm * rtprcp'],
                    'sam': [self.get_rc_coef_zm_X_rtprcp_sam_calc,'rc_coef_zm * rtprcp'],
                    'coamps': [self.get_rc_coef_zm_X_rtprcp_coamps_calc, 'rc_coef_zm * rtprcp'],
                    'r408': ['rc_coef_zm * rtprcp'],
                    'hoc': ['rc_coef_zm * rtprcp'],
                    'e3sm': ['rc_coef_zm * rtprcp'],
                    'cam': ['rc_coef_zm * rtprcp'],
                    'wrf': [self.get_rc_coef_zm_X_rtprcp_wrf_calc,'rc_coef_zm * rtprcp'],
                },
                'title': 'Contribution of Cloud Water Flux to rtpthvp',
                'axis_title': 'rc_coef_zm * rtprcp [kg/kg K]',
                'sci_scale': -4
            },
            {'var_names':
                {
                    'clubb': [self.get_rc_coef_X_wp2rcp_clubb_calc, 'rc_coef_zm * wp2rcp'],
                    'sam': [self.get_rc_coef_X_wp2rcp_sam_calc,'rc_coef_zm * wp2rcp'],
                    'coamps': [self.get_rc_coef_X_wp2rcp_coamps_calc, 'rc_coef_zm * wp2rcp'],
                    'r408': ['rc_coef_zm * wp2rcp'],
                    'hoc': ['rc_coef_zm * wp2rcp'],
                    'e3sm': ['rc_coef_zm * wp2rcp'],
                    'cam': ['rc_coef_zm * wp2rcp'],
                    'wrf': [self.get_rc_coef_X_wp2rcp_wrf_calc,'rc_coef_zm * wp2rcp'],
                },
                'title': 'Cloud water contribution to wp2thvp',
                'axis_title': 'rc_coef * wp2rcp [m^2/s^2 K]',
                'sci_scale': 0
            },
            {'var_names':
                {
                    'clubb': ['wp2thvp'],
                    'sam': ['WP2THVP'],
                    'coamps': ['wp2thvp'],
                    'r408': ['wp2thvp'],
                    'hoc': ['wp2thvp'],
                    'e3sm': ['wp2thvp'],
                    'cam': ['wp2thvp'],
                    'wrf': ['wp2thvp'],
                },
                'sci_scale': 0,
            },
            {'var_names':
                {
                    'clubb': ['em'],
                    'sam': [self.get_tke_sam_calc,'TKE'],
                    'coamps': ['em'],
                    'r408': ['em'],
                    'hoc': ['em'],
                    'e3sm': ['em'],
                    'cam': ['em'],
                    'wrf': ['em','CSP_TKE_RS'],
                },
            },
            {'var_names':
                {
                    'clubb':['sigma_sqd_w'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)

    def getThlmSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates thlm values from sam output using
        the following equation

        ``(THETAL + 2500.4.*(THETA/TABS).*(QI/1000))``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        # z,z, dataset = self.getVarForCalculations('z', self.sam_benchmark_dataset)
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        thetal, indep, dataset = self.getVarForCalculations('THETAL', dataset)
        theta, indep, dataset = self.getVarForCalculations('THETA', dataset)
        tabs, indep, dataset = self.getVarForCalculations('TABS', dataset)
        qi, indep, dataset = self.getVarForCalculations('QI', dataset)

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        return thlm, indep

    def getRtmSamCalc(self, dataset_override=None):
        """
         This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

       Calculates rtm values from sam output using
        the following equation

        ``(QT-QI) / 1000``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qt, indep, dataset = self.getVarForCalculations('QT', dataset)
        qi, indep, dataset = self.getVarForCalculations('QI', dataset)

        rtm = (qt - qi) / 1000
        if np.any(np.isnan(qi)):
            rtm = qt / 1000

        return rtm, indep

    def getSkwZtLesCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates Skw_zt values from sam output using
        the following equation

        ``WP3 / (WP2 + 1.6e-3)**1.5``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset
        if self.coamps_benchmark_dataset is not None:
            dataset = self.coamps_benchmark_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override

        wp3, indep, dataset = self.getVarForCalculations(['WP3', 'W3', 'wp3'], dataset)
        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        skw_zt = wp3 / (wp2 + 1.6e-3) ** 1.5

        return skw_zt, indep

    def getSkrtZtLesCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates Skrt_zt values from sam output using
        the following equation

         sam eqn
            ``RTP3 / (RTP2 + 4e-16)**1.5``
         coamps eqn
            ``qtp3 / (qtp2 + 4e-16)**1.5``
            ``rtp3 / (rtp2 + 4e-16)**1.5``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset

        if self.coamps_benchmark_dataset is not None:
            dataset = self.coamps_benchmark_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override
        rtp3, indep, dataset = self.getVarForCalculations(['RTP3', 'qtp3', 'rtp3'], dataset)
        rtp2, indep, dataset = self.getVarForCalculations(['RTP2', 'qtp2', 'rtp2', 'rlp2'], dataset)

        skrt_zt = rtp3 / (rtp2 + 4e-16) ** 1.5

        return skrt_zt, indep

    def getSkthlZtLesCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates Skthl_zt values from sam output using
        the following equation

        sam
        ``THLP3 / (THLP2 + 4e-4)**1.5``

        coamps
        ``thlp3 / (thlp2 + 4e-4)**1.5``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset

        if self.coamps_benchmark_dataset is not None:
            dataset = self.coamps_benchmark_dataset['sm']

        if dataset_override is not None:
            dataset = dataset_override


        thlp3, indep, dataset = self.getVarForCalculations(['THLP3', 'thlp3'], dataset)
        thlp2, indep, dataset = self.getVarForCalculations(['THLP2', 'thlp2'], dataset)

        skthl_zt = thlp3 / (thlp2 + 4e-4)**1.5
        return skthl_zt, indep

    def getWpthlpSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.

       ``WPTHLP = (TLFLUX) / (RHO * 1004) + WPTHLP_SGS``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        tlflux, indep, dataset = self.getVarForCalculations(['TLFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        wpthlp_sgs, indep, dataset = self.getVarForCalculations(['WPTHLP_SGS'], dataset)

        wpthlp = (tlflux / (rho * 1004))

        if not np.any(np.isnan(wpthlp_sgs)):
            wpthlp += wpthlp_sgs

        return wpthlp, indep

    def getWprtpSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.

        ``WPRTP = (QTFLUX) / (RHO * 2.5104e+6) + WPRTP_SGS``

        ``WPRTP = (QTFLUX) / (RHO * 2.5104e+6)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qtflux, indep, dataset = self.getVarForCalculations(['QTFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        wprtp_sgs, indep, dataset = self.getVarForCalculations(['WPRTP_SGS'], dataset)

        wprtp = qtflux / (rho * 2.5104e+6)

        if not np.any(np.isnan(wprtp_sgs)):
            wprtp += wprtp_sgs

        return wprtp, indep

    def getWpthvpSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.

        ``WPTHVP =  (TVFLUX) / ( RHO * 1004)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        tvflux, indep, dataset = self.getVarForCalculations(['TVFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        wpthvp = tvflux / (rho * 1004)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_benchmark_dataset)
        return wpthvp, indep

    def getRtp2SamCalc(self, dataset_override=None):
        """
         This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

       This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.

        ``(QT2 / 1e+6) + RTP2_SGS``

        ``(QT2 / 1e+6)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = None
        if self.sam_benchmark_dataset is not None:
            dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QT2, indep, dataset = self.getVarForCalculations(['QT2'], dataset)
        RTP2_SGS, indep, dataset = self.getVarForCalculations(['RTP2_SGS'], dataset)
        rtp2 = (QT2 / 1e+6) + RTP2_SGS
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_benchmark_dataset)
        return rtp2, indep

    def getRtp3SamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Caclulates Rtp3 output

        ``rc_coef_zm .* rtprcp``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        rtp3 = None
        z = None
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_benchmark_dataset
        if isinstance(dataset, Dataset):
            dataset = {'temp': dataset}
        for dataset in dataset.values():
            if 'rc_coef_zm' in dataset.variables.keys() and 'rtprcp' in dataset.variables.keys():
                rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
                rtprcp, indep, dataset = self.getVarForCalculations('rtprcp', dataset)
                rtp3 = rc_coef_zm * (rtprcp)

            elif 'QCFLUX' in dataset.variables.keys():
                QCFLUX, indep, dataset = self.getVarForCalculations('QCFLUX', dataset)
                RHO, indep, dataset = self.getVarForCalculations('RHO', dataset)
                PRES, indep, dataset = self.getVarForCalculations('PRES', dataset)
                THETAV, indep, dataset = self.getVarForCalculations('THETAV', dataset)
                rtp3 = ((QCFLUX) / (RHO * 2.5104e+6)) * (
                        2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return rtp3, indep

    def get_rc_coef_zm_X_wprcp_clubb_line(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the Contribution of Cloud Water Flux
        to wpthvp using the equation

        ``rc_coef_zm .* wprcp``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        wprcp, indep, dataset = self.getVarForCalculations('wprcp', dataset)

        output = rc_coef_zm * wprcp
        return output, indep

    def get_rc_coef_zm_X_wprcp_wrf_line(self, dataset_override=None):
        """
        Same as above function except used for WRF datasets.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.wrf_datasets['zm']
        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        wprcp, indep, dataset = self.getVarForCalculations('wprcp', dataset)

        output = rc_coef_zm * wprcp
        return output, indep

    def get_rc_coef_zm_X_wprcp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the Contribution of Cloud Water Flux
        to wpthvp for SAM using the equation

        sam eqn
        ``WPRCP                          * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)``

        ``((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_benchmark_dataset

        WPRCP, indep, dataset = self.getVarForCalculations('WPRCP', dataset)
        QCFLUX, indep, dataset = self.getVarForCalculations('QCFLUX', dataset)
        RHO, indep, dataset = self.getVarForCalculations('RHO', dataset)
        PRES, indep, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, indep, dataset = self.getVarForCalculations('THETAV', dataset)

        output1 = ((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        output2 = WPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)

        output = self.pickNonZeroOutput(output1, output2, favor_output=output2)

        return output, indep

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_clubb_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the Contribution of Cloud Water Flux
        to thlprcp using the equation

        ``rc_coef_zm * thlprcp``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        thlprcp, indep, dataset = self.getVarForCalculations('thlprcp', dataset)

        output = rc_coef_zm * thlprcp
        return output, indep

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_wrf_calc(self, dataset_override=None):
        """
        Same as above but for WRF datasets
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.wrf_datasets['zm']

        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        thlprcp, indep, dataset = self.getVarForCalculations('thlprcp', dataset)

        output = rc_coef_zm * thlprcp
        return output, indep

    def get_rc_coef_zm_X_rtprcp_clubb_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the Contribution of Cloud Water Flux
        to rtprcp using the equation

        `rc_coef_zm * rtprcp`

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        rtprcp, indep, dataset = self.getVarForCalculations('rtprcp', dataset)

        output = rc_coef_zm * rtprcp
        return output, indep

    def get_rc_coef_zm_X_rtprcp_wrf_calc(self, dataset_override=None):
        """
        Same as above except for WRF datasets
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.wrf_datasets['zm']

        rc_coef_zm, indep, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        rtprcp, indep, dataset = self.getVarForCalculations('rtprcp', dataset)

        output = rc_coef_zm * rtprcp
        return output, indep

    def getUwCoampsData(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        coamps eqn
        ``upwp = wpup + wpup_sgs``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
         :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_benchmark_dataset['sw']
        wpup, indep, dataset = self.getVarForCalculations('wpup', dataset)
        wpup_sgs, indep, dataset = self.getVarForCalculations('wpup_sgs', dataset)

        upwp = wpup + wpup_sgs
        return upwp, indep

    def getVwCoampsData(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        coamps eqn

        ``vpwp = wpvp + wpvp_sgs``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
         :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_benchmark_dataset['sw']
        wpvp, indep, dataset = self.getVarForCalculations('wpvp', dataset)
        wpvp_sgs, indep, dataset = self.getVarForCalculations('wpvp_sgs', dataset)

        vpwp = wpvp + wpvp_sgs
        return vpwp, indep

    def get_rc_coef_zm_X_wprcp_coamps_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        coamps eqn's

        ``thlpqcp * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)``

        ``wpqcp   * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)``

        ``wprlp   * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)``

        ``wprlp   * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return:  tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.coamps_benchmark_dataset['sw']
        if dataset_override is not None:
            dataset = dataset_override['sw']

        wprlp, indep, dataset = self.getVarForCalculations(['thlpqcp', 'wpqcp', 'wprlp'], dataset)
        ex0, indep, dataset = self.getVarForCalculations(['ex0'], dataset)
        p, indep, dataset = self.getVarForCalculations('p', dataset)
        thvm, indep, dataset = self.getVarForCalculations('thvm', dataset)
        output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
        output = self.pickNonZeroOutput(output1, output2)
        return output, indep

    def get_rc_coef_zm_X_thlprcp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        sam eqn
        ``THLPRCP .* (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override
        THLPRCP, indep, dataset = self.getVarForCalculations(['THLPRCP'], dataset)
        PRES, indep, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, indep, dataset = self.getVarForCalculations('THETAV', dataset)

        output = THLPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, indep

    def get_rc_coef_zm_X_thlprcp_coamps_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        coamps eqn
        .. code-block:: python
            :linenos:

            thlpqcp * (2.5e6 / (1004.67 * ex0)                             - 1.61*thvm)
            thlprlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61*thvm)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.coamps_benchmark_dataset['sw']
        if dataset_override is not None:
            dataset = dataset_override
        thlpqcp, indep, dataset = self.getVarForCalculations(['thlpqcp'], dataset)
        ex0, indep, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, indep, dataset = self.getVarForCalculations('thvm', dataset)

        thlprlp, indep, dataset = self.getVarForCalculations(['thlprlp'], dataset)
        p, indep, dataset = self.getVarForCalculations('p', dataset)

        output1 = thlpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = thlprlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61*thvm)
        output = self.pickNonZeroOutput(output1, output2)

        return output, indep

    def get_rc_coef_zm_X_rtprcp_coamps_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        coamps eqn
        ``qtpqcp * (2.5e6 / (1004.67*ex0)                           - 1.61*thvm)``

        ``qtpqcp * (2.5e6 / (1004.67*ex0)                           - 1.61*thvm)``

        ``rtprlp * (2.5e6 / (1004.67*((p/1.0e5)**(287.04/1004.67))) - 1.61*thvm)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.coamps_benchmark_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override['sm']

        qtpqcp, indep, dataset = self.getVarForCalculations(['qtpqcp', 'rtprcp'], dataset)
        ex0, indep, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, indep, dataset = self.getVarForCalculations('thvm', dataset)

        rtprlp, indep, dataset = self.getVarForCalculations('rtprlp', dataset)
        p, indep, dataset = self.getVarForCalculations('p', dataset)

        output1 = qtpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = rtprlp * (2.5e6 / (1004.67*((p/1.0e5)**(287.04/1004.67))) - 1.61*thvm)
        output = self.pickNonZeroOutput(output1, output2)

        return output, indep

    def get_rc_coef_zm_X_rtprcp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        sam eqn
        ``RTPRCP * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override

        RTPRCP, indep, dataset = self.getVarForCalculations('RTPRCP', dataset)
        PRES, indep, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, indep, dataset = self.getVarForCalculations('THETAV', dataset)

        output = RTPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, indep

    def get_rc_coef_X_wp2rcp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``WP2RCP * (2.5e6 / (1004.67*((PRES / 1000)^(287.04/1004.67))) - 1.61*THETAV)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override#['sam']
        WP2RCP, indep, dataset = self.getVarForCalculations('WP2RCP', dataset)
        PRES, indep, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, indep, dataset = self.getVarForCalculations('THETAV', dataset)

        output = WP2RCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, indep

    def get_rc_coef_X_wp2rcp_clubb_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``rc_coef * wrp2rcp``


        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        rc_coef, indep, dataset = self.getVarForCalculations('rc_coef', dataset)
        wp2rcp, indep, dataset = self.getVarForCalculations('wp2rcp', dataset)

        output = rc_coef * wp2rcp
        return output, indep

    def get_rc_coef_X_wp2rcp_wrf_calc(self, dataset_override=None):
        """
        Same as above except for WRF datasets
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.wrf_datasets['zm']
        rc_coef, indep, dataset = self.getVarForCalculations('rc_coef', dataset)
        wp2rcp, indep, dataset = self.getVarForCalculations('wp2rcp', dataset)

        output = rc_coef * wp2rcp
        return output, indep

    def get_rc_coef_X_wp2rcp_coamps_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``wp2qcp * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)``

        ``wp2rlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        wp2qcp, indep, dataset = self.getVarForCalculations('wp2qcp', dataset)
        ex0, indep, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, indep, dataset = self.getVarForCalculations('thvm', dataset)

        wp2rlp, indep, dataset = self.getVarForCalculations('wp2rlp', dataset)
        p, indep, dataset = self.getVarForCalculations('p', dataset)

        output1 = wp2qcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = wp2rlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)

        output = self.pickNonZeroOutput(output1, output2)

        return output, indep


    def get_wp2_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``WP2_SGS + W2``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        WP2_SGS, z, dataset = self.getVarForCalculations('WP2_SGS', dataset)
        W2, z, dataset = self.getVarForCalculations('W2', dataset)

        output = WP2_SGS + W2

        return output, z

    def get_wp3_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``WP3_SGS + W3``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        WP3_SGS, z, dataset = self.getVarForCalculations('WP3_SGS', dataset)
        W3, z, dataset = self.getVarForCalculations('W3', dataset)

        output = WP3_SGS + W3

        return output, z

    def get_thlp2_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``TL2 + THLP2_SGS``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        TL2, z, dataset = self.getVarForCalculations('TL2', dataset)
        THLP2_SGS, z, dataset = self.getVarForCalculations('THLP2_SGS', dataset)

        output = TL2 + THLP2_SGS

        return output, z

    def get_upwp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``UW + UPWP_SGS``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        UW, z, dataset = self.getVarForCalculations('UW', dataset)
        UPWP_SGS, z, dataset = self.getVarForCalculations('UPWP_SGS', dataset)

        output = UW + UPWP_SGS

        return output, z

    def get_vpwp_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``VW + VPWP_SGS``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        VW, z, dataset = self.getVarForCalculations('VW', dataset)
        VPWP_SGS, z, dataset = self.getVarForCalculations('VPWP_SGS', dataset)

        output = VW + VPWP_SGS

        return output, z

    def get_up2_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``UP2_SGS + U2``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        U2, z, dataset = self.getVarForCalculations('U2', dataset)
        UP2_SGS, z, dataset = self.getVarForCalculations('UP2_SGS', dataset)

        output = UP2_SGS + U2

        return output, z

    def get_vp2_sam_calc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        ``VP2_SGS + V2``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        V2, z, dataset = self.getVarForCalculations('V2', dataset)
        VP2_SGS, z, dataset = self.getVarForCalculations('VP2_SGS', dataset)

        output = VP2_SGS + V2

        return output, z

    def get_tke_sam_calc(self, dataset_override=None):
        """
        This function calculates TKE from SAM data by explicitly summing the squared
        components of the velocity.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        U2, z, dataset = self.getVarForCalculations('U2', dataset)
        UP2_SGS, z, dataset = self.getVarForCalculations('UP2_SGS', dataset)
        V2, z, dataset = self.getVarForCalculations('V2', dataset)
        VP2_SGS, z, dataset = self.getVarForCalculations('VP2_SGS', dataset)
        W2, z, dataset = self.getVarForCalculations('W2', dataset)
        WP2_SGS, z, dataset = self.getVarForCalculations('WP2_SGS', dataset)

        output = 0.5 * ( U2 + UP2_SGS + V2 + VP2_SGS + W2 + WP2_SGS )

        return output, z

    def getBVSqdSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates values of the squared Brunt-Väisälä frequency from sam output using
        the following equation

        ``N^2 = g/THETAV * dTHETAV/dz``
        where g is earth's gravitational acceleration constant (9.81 m/s^2)
        and dTHETAV/dz is the derivative of THETAV w.r.t. height (difference quotient here, for obv reasons)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        # z,z, dataset = self.getVarForCalculations('z', self.sam_benchmark_dataset)
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        thetav, z, dataset = self.getVarForCalculations('THETAV', dataset)

        # Calculate numeical derivative of theta by difference quotient.
        # Since the difference quotient "eliminates" one array entry,
        # the derivative array needs to be extended again.
        # Here, a 0 entry is added at the topmost level
        tmp = np.diff(thetav)/np.diff(z)
        dthetav_dz = np.zeros(thetav.shape)
        dthetav_dz[0:-1] = tmp
        dthetav_dz[-1] = tmp[-1]
        bv_sqd = 9.81/thetav * dthetav_dz
        return bv_sqd, z