"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup

class VariableGroupPaperPlots(VariableGroup):
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

        norm_grid_density_plot = [
            {'var_names': ['norm_min_grid_dens'], 'legend_label': r"$g_\mathrm{min}(z)$"},
            {'var_names': ['norm_grid_dens'], 'legend_label': r"$g(z)$"},
        ]

        each_ref_crit_term = [
            {'var_names': ['alt_term'], 'legend_label': r"$3 \cdot g_z (z)$"},
            {'var_names': ['chi_term_time_avg'], 'legend_label': r"$100 \cdot \bar{g}_\chi (z)$"},
            {'var_names': ['brunt_term_time_avg'], 'legend_label': r"$7.5 \cdot \bar{g}_\mathrm{Ri} (z)$"},
        ]

        RTM_VAR_MIN_ARM = 0
        RTM_VAR_MAX_ARM = 0.017
        RTM_VAR_MIN_ASTEX = 0
        RTM_VAR_MAX_ASTEX = 0.0119
        RTM_VAR_MIN_GABLS2 = 0
        RTM_VAR_MAX_GABLS2 = 0.00296

        RCM_VAR_MIN_ARM = 0
        RCM_VAR_MAX_ARM = 0.0002
        RCM_VAR_MIN_ASTEX = 0
        RCM_VAR_MAX_ASTEX = 0.0006
        RCM_VAR_MIN_DYCOMS2 = 0
        RCM_VAR_MAX_DYCOMS2 = 0.000067

        WPTHLP_VAR_MIN_ARM = -0.42337
        #WPTHLP_VAR_MAX_ARM = 0.12668
        WPTHLP_VAR_MAX_ARM = 0.15
        WPTHLP_VAR_MIN_ASTEX = -0.05
        #WPTHLP_VAR_MIN_ASTEX = -0.008
        #WPTHLP_VAR_MAX_ASTEX = 0.008
        WPTHLP_VAR_MAX_ASTEX = 0.05
        #WPTHLP_VAR_MIN_GABLS2 = -0.04
        WPTHLP_VAR_MIN_GABLS2 = -0.055
        WPTHLP_VAR_MAX_GABLS2 = 0.06
        WPTHLP_VAR_MIN_DYCOMS2_RF01 = -0.02
        WPTHLP_VAR_MAX_DYCOMS2_RF01 = 0.02

        CLOUD_FRAC_VAR_MIN_ARM = 0
        CLOUD_FRAC_VAR_MAX_ARM = 0.78
        CLOUD_FRAC_VAR_MIN_ASTEX = 0
        CLOUD_FRAC_VAR_MAX_ASTEX = 1.0
        CLOUD_FRAC_VAR_MIN_DYCOMS2_RF01 = 0
        CLOUD_FRAC_VAR_MAX_DYCOMS2_RF01 = 0.5

        WP2_VAR_MIN_ARM = 0
        WP2_VAR_MAX_ARM = 1.89
        WP2_VAR_MIN_ASTEX = 0
        WP2_VAR_MAX_ASTEX = 0.49
        WP2_VAR_MIN_GABLS2 = 0
        WP2_VAR_MAX_GABLS2 = 0.48
        WP2_VAR_MIN_DYCOMS2_RF01 = 0
        WP2_VAR_MAX_DYCOMS2_RF01 = 0.4

        self.variable_definitions = [
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
                'title': 'ARM',
                'file_identifier': 'norm_grid_dens',
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
                'lines': norm_grid_density_plot,
                'file_identifier': 'norm_grid_dens',
                'title': 'ASTEX',
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
                'lines': norm_grid_density_plot,
                'file_identifier': 'norm_grid_dens',
                'title': 'GABLS2',
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
                'lines': norm_grid_density_plot,
                'file_identifier': 'norm_grid_dens',
                'title': 'DYCOMS2',
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
                'file_identifier': 'each_ref_crit_term',
                'title': 'ARM',
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
                'file_identifier': 'each_ref_crit_term',
                'title': 'ASTEX',
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
                'file_identifier': 'each_ref_crit_term',
                'title': 'GABLS2',
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
                'file_identifier': 'each_ref_crit_term',
                'title': 'DYCOMS2',
                'axis_title': 'grid density [1/m]',
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
                'priority': True,
                'file_identifier': 'rcm',
                'title': 'ARM',
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
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
                'priority': True,
                'file_identifier': 'rcm',
                'title': 'ASTEX',
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
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
                'priority': True,
                'file_identifier': 'rcm',
                'title': 'DYCOMS2',
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
                'file_identifier': 'rtm',
                'title': 'ARM',
                'axis_title': r"$\bar{r}_t ~~~[\mathrm{kg}/\mathrm{kg}]$",
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
                'file_identifier': 'rtm',
                'title': 'ASTEX',
                'axis_title': r"$\bar{r}_t ~~~[\mathrm{kg}/\mathrm{kg}]$",
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
                'file_identifier': 'rtm',
                'title': 'GABLS2',
                'axis_title': r"$\bar{r}_t ~~~[\mathrm{kg}/\mathrm{kg}]$",
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm','CSP_QT'],
                },
                'sci_scale': -3,
                'priority': True,
                'file_identifier': 'rtm',
                'title': 'DYCOMS2',
                'axis_title': r"$\bar{r}_t ~~~[\mathrm{kg}/\mathrm{kg}]$",
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wpthlp',
                'title': 'ARM',
                'var_min': WPTHLP_VAR_MIN_ARM,
                'var_max': WPTHLP_VAR_MAX_ARM,
                'axis_title': r"$\overline{w^{\prime}\theta^{\prime}_l} ~~~[\mathrm{m}\mathrm{s}^{-1}\mathrm{K}]$",
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wpthlp',
                'title': 'ASTEX',
                'var_min': WPTHLP_VAR_MIN_ASTEX,
                'var_max': WPTHLP_VAR_MAX_ASTEX,
                'axis_title': r"$\overline{w^{\prime}\theta^{\prime}_l} ~~~[\mathrm{m}\mathrm{s}^{-1}\mathrm{K}]$",
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wpthlp',
                'title': 'GABLS2',
                'var_min': WPTHLP_VAR_MIN_GABLS2,
                'var_max': WPTHLP_VAR_MAX_GABLS2,
                'axis_title': r"$\overline{w^{\prime}\theta^{\prime}_l} ~~~[\mathrm{m}\mathrm{s}^{-1}\mathrm{K}]$",
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': [],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp','CSP_WTHL'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wpthlp',
                'title': 'DYCOMS2',
                'var_min': WPTHLP_VAR_MIN_DYCOMS2_RF01,
                'var_max': WPTHLP_VAR_MAX_DYCOMS2_RF01,
                'axis_title': r"$\overline{w^{\prime}\theta^{\prime}_l} ~~~[\mathrm{m}\mathrm{s}^{-1}\mathrm{K}]$",
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'file_identifier': 'rcm',
                'title': 'ARM',
                'var_min': RCM_VAR_MIN_ARM,
                'var_max': RCM_VAR_MAX_ARM,
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'file_identifier': 'rcm',
                'title': 'ASTEX',
                'var_min': RCM_VAR_MIN_ASTEX,
                'var_max': RCM_VAR_MAX_ASTEX,
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'file_identifier': 'rcm',
                'title': 'DYCOMS2',
                'var_min': RCM_VAR_MIN_DYCOMS2,
                'var_max': RCM_VAR_MAX_DYCOMS2,
                'axis_title': r"$\bar{r}_c \text{ [kg/kg]}$",
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
                'file_identifier': 'cloud_frac',
                'title': 'ARM',
                'var_min': CLOUD_FRAC_VAR_MIN_ARM,
                'var_max': CLOUD_FRAC_VAR_MAX_ARM,
                'axis_title': 'Cloud fraction    [$-$]',
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
                'file_identifier': 'cloud_frac',
                'title': 'ASTEX',
                'var_min': CLOUD_FRAC_VAR_MIN_ASTEX,
                'var_max': CLOUD_FRAC_VAR_MAX_ASTEX,
                'axis_title': 'Cloud fraction    [$-$]',
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
                'file_identifier': 'cloud_frac',
                'title': 'DYCOMS2',
                'var_min': CLOUD_FRAC_VAR_MIN_DYCOMS2_RF01,
                'var_max': CLOUD_FRAC_VAR_MAX_DYCOMS2_RF01,
                'axis_title': 'Cloud fraction    [$-$]',
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wp2',
                'title': 'ARM',
                'var_min': WP2_VAR_MIN_ARM,
                'var_max': WP2_VAR_MAX_ARM,
                'axis_title': r"$\overline{{w^{\prime}}^2} ~~~[\mathrm{m}^2/\mathrm{s}^2]$",
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wp2',
                'title': 'ASTEX',
                'var_min': WP2_VAR_MIN_ASTEX,
                'var_max': WP2_VAR_MAX_ASTEX,
                'axis_title': r"$\overline{{w^{\prime}}^2} ~~~[\mathrm{m}^2/\mathrm{s}^2]$",
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wp2',
                'title': 'GABLS2',
                'var_min': WP2_VAR_MIN_GABLS2,
                'var_max': WP2_VAR_MAX_GABLS2,
                'axis_title': r"$\overline{{w^{\prime}}^2} ~~~[\mathrm{m}^2/\mathrm{s}^2]$",
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': [],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2','CSP_W2'],
                },
                'sci_scale': 0,
                'priority': True,
                'file_identifier': 'wp2',
                'title': 'DYCOMS2',
                'var_min': WP2_VAR_MIN_DYCOMS2_RF01,
                'var_max': WP2_VAR_MAX_DYCOMS2_RF01,
                'axis_title': r"$\overline{{w^{\prime}}^2} ~~~[\mathrm{m}^2/\mathrm{s}^2]$",
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)