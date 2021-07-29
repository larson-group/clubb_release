"""
:author: Benjamin A. Stephens 
:date: April 2021
"""

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupBaseBudgetsSamStyle(VariableGroup):
    """

    """

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None, e3sm_datasets=None, wrf_datasets=None,
                 priority_vars=False):
        self.name = "base variables budgets"

#        rtm_budget_lines = [
#            {'var_names': ['rtm_bt'], 'legend_label': 'rtm_bt'},
#            {'var_names': ['rtm_ma'], 'legend_label': 'rtm_ma'},
#            {'var_names': ['rtm_ta'], 'legend_label': 'rtm_ta'},
#            {'var_names': ['rtm_mc'], 'legend_label': 'rtm_mc'},
#            {'var_names': ['rtm_clipping', self.getRtmClipping],
#                 'legend_label': 'rtm_clipping',
#             },
#            {'var_names': ['rtm_pd'],
#                'legend_label': 'rtm_pd'
#             },
#            {'var_names': ['ls_forcing', self.getRtmForcing],
#                 'legend_label': 'ls_forcing',
#             },
#            {'var_names': ['rtm_residual', self.getRtmResidual],
#                 'legend_label': 'rtm_residual',
#             },
#
#        ]

        wpthlp_budget_lines = [
            {'var_names': ['wpthlp_residual', self.getWpthlpResidual], 'legend_label': 'wpthlp_res'},
            {'var_names': ['wpthlp_tp'], 'legend_label': 'wpthlp_grad'},
            {'var_names': ['wpthlp_adv',self.calc_wpthlp_adv], 'legend_label': 'wpthlp_adv'},
            {'var_names': ['wpthlp_dp1'], 'legend_label': 'wpthlp_dfsn'},
            {'var_names': ['wpthlp_bp'], 'legend_label': 'wpthlp_buoy'},
            {'var_names': ['wpthlp_pres',self.calc_wpthlp_pres], 'legend_label': 'wpthlp_pres'},
            {'var_names': ['wpthlp_mc'], 'legend_label': 'wpthlp_prec'},
            {'var_names': ['wpthlp_rad',self.calc_wpthlp_rad], 'legend_label': 'wpthlp_rad'},
            {'var_names': ['wpthlp_forc',self.calc_wpthlp_forc], 'legend_label': 'wpthlp_forc'},
            {'var_names': ['wpthlp_bt'], 'legend_label': 'wpthlp_bt'},
            {'var_names': ['wpthlp_limiters',self.calc_wpthlp_limiters], 'legend_label': 'wpthlp_limiters'},
        ]

        wprtp_budget_lines = [
            {'var_names': ['wprtp_residual', self.getWprtpResidual], 'legend_label': 'wprtp_res'},
            {'var_names': ['wprtp_tp'], 'legend_label': 'wprtp_grad'},
            {'var_names': ['wprtp_adv',self.calc_wprtp_adv], 'legend_label': 'wprtp_adv'},
            {'var_names': ['wprtp_dp1'], 'legend_label': 'wprtp_dfsn'},
            {'var_names': ['wprtp_bp'], 'legend_label': 'wprtp_buoy'},
            {'var_names': ['wprtp_pres',self.calc_wprtp_pres], 'legend_label': 'wprtp_pres'},
            {'var_names': ['wprtp_mc'], 'legend_label': 'wprtp_prec'},
            {'var_names': ['wprtp_forc',self.calc_wprtp_forc], 'legend_label': 'wprtp_forc'},
            {'var_names': ['wprtp_bt'], 'legend_label': 'wprtp_bt'},
            {'var_names': ['wprtp_limiters',self.calc_wprtp_limiters], 'legend_label': 'wprtp_limiters'},
        ]

        wp2_budget_lines = [
            {'var_names': ['wp2_residual', self.getWp2Residual], 'legend_label': 'wp2_res'},
            {'var_names': ['wp2_adv',self.calc_wp2_adv], 'legend_label': 'wp2_adv'},
            {'var_names': ['wp2_pres',self.calc_wp2_pres], 'legend_label': 'wp2_pres'},
            {'var_names': ['wp2_pr1'], 'legend_label': 'wp2_redis'},
            {'var_names': ['wp2_bp'], 'legend_label': 'wp2_buoy'},
            {'var_names': ['wp2_dfsn',self.calc_wp2_dfsn], 'legend_label': 'wp2_dfsn'},
            {'var_names': ['wp2_sdmp'], 'legend_label': 'wp2_sdmp'},
            {'var_names': ['wp2_bt'], 'legend_label': 'wp2_bt'},
            {'var_names': ['wp2_limiters',self.calc_wp2_limiters], 'legend_label': 'wp2_limiters'},
            {'var_names': ['wp2_sf'], 'legend_label': 'wp2_sf'},
        ]

        wp3_budget_lines = [
            {'var_names': ['wp3_residual', self.getWp3Residual], 'legend_label': 'wp3_res'},
            {'var_names': ['wp3_pres',self.calc_wp3_pres], 'legend_label': 'wp3_pres'},
            {'var_names': ['wp3_adv',self.calc_wp3_adv], 'legend_label': 'wp3_adv'},
            {'var_names': ['wp3_bp1'], 'legend_label': 'wp3_buoy'},
            {'var_names': ['wp3_dp1'], 'legend_label': 'wp3_dfsn'},
            {'var_names': ['wp3_bt'], 'legend_label': 'wp3_bt'},
            {'var_names': ['wp3_cl'], 'legend_label': 'wp3_limiters'},
        ]

        # According to https://carson.math.uwm.edu/larson-group/internal/SAM_LES_BUDGET_PLOTS/,
        # the DISSIP and DIFFTR terms in SAM should be equivalent to CLUBB's dp1 and dp2 terms,
        # respectively.  The DISSIP term is calculated by (approximately), K * (dX/dz )^2, 
        # where K is eddy conductivity [m2/s].  (See e.g. around lines 486+ of SGS_TKE/sgs.F90.)
        # The DIFFTR term is calculated as (approximately), (X(curr)^2-X(prev)^2)/dt, 
        # (also from the sgs.F90 file) where a general diffusion subroutine is applied 
        # between the "curr" and "prev" values.  Because the DISSIP term more closely resembles
        # CLUBB's dp2 term, however, I've switched them here so that DISSIP = dp2 and DIFFTR = dp1.
        thlp2_budget_lines = [
            {'var_names': ['thlp2_residual', self.getThlp2Residual], 'legend_label': 'thlp2_res'},
            {'var_names': ['thlp2_adv',self.calc_thlp2_adv], 'legend_label': 'thlp2_adv'},
            {'var_names': ['thlp2_tp'], 'legend_label': 'thlp2_grad'},
            {'var_names': ['thlp2_dp2'], 'legend_label': 'thlp2_dissip'},
            {'var_names': ['thlp2_dp1'], 'legend_label': 'thlp2_difftr'},
            {'var_names': ['thlp2_mc'], 'legend_label': 'thlp2_prec'},
            {'var_names': ['thlp2_rad',self.calc_wpthlp_rad], 'legend_label': 'thlp2_rad'}, 
            {'var_names': ['thlp2_forcing',self.calc_thlp2_forc], 'legend_label': 'thlp2_forc'},
            {'var_names': ['thlp2_bt'], 'legend_label': 'thlp2_bt'},
            {'var_names': ['thlp2_limiters',self.calc_thlp2_limiters], 'legend_label': 'thlp2_limiters'},
            {'var_names': ['thlp2_sf'], 'legend_label': 'thlp2_sf'},
        ]

        # According to https://carson.math.uwm.edu/larson-group/internal/SAM_LES_BUDGET_PLOTS/,
        # the DISSIP and DIFFTR terms in SAM should be equivalent to CLUBB's dp1 and dp2 terms,
        # respectively.  The DISSIP term is calculated by (approximately), K * (dX/dz )^2, 
        # where K is eddy conductivity [m2/s].  (See e.g. around lines 486+ of SGS_TKE/sgs.F90.)
        # The DIFFTR term is calculated as (approximately), (X(curr)^2-X(prev)^2)/dt, 
        # (also from the sgs.F90 file) where a general diffusion subroutine is applied 
        # between the "curr" and "prev" values.  Because the DISSIP term more closely resembles
        # CLUBB's dp2 term, however, I've switched them here so that DISSIP = dp2 and DIFFTR = dp1.
        rtp2_budget_lines = [
            {'var_names': ['rtp2_residual', self.getRtp2Residual], 'legend_label': 'rtp2_res'},
            {'var_names': ['rtp2_adv',self.calc_rtp2_adv], 'legend_label': 'rtp2_adv'},
            {'var_names': ['rtp2_tp'], 'legend_label': 'rtp2_grad'},
            {'var_names': ['rtp2_dp2'], 'legend_label': 'rtp2_dissip'},
            {'var_names': ['rtp2_dp1'], 'legend_label': 'rtp2_difftr'},
            {'var_names': ['rtp2_mc'], 'legend_label': 'rtp2_prec'},
            {'var_names': ['rtp2_forcing',self.calc_rtp2_forc], 'legend_label': 'rtp2_forc'},
            {'var_names': ['rtp2_bt'], 'legend_label': 'rtp2_bt'},
            {'var_names': ['rtp2_limiters',self.calc_rtp2_limiters], 'legend_label': 'rtp2_limiters'},
            {'var_names': ['rtp2_sf'], 'legend_label': 'rtp2_sf'},
        ]

        # According to https://carson.math.uwm.edu/larson-group/internal/SAM_LES_BUDGET_PLOTS/,
        # the DISSIP and DIFFTR terms in SAM should be equivalent to CLUBB's dp1 and dp2 terms,
        # respectively.  The DISSIP term is calculated by (approximately), K * (dX/dz )^2, 
        # where K is eddy conductivity [m2/s].  (See e.g. around lines 486+ of SGS_TKE/sgs.F90.)
        # The DIFFTR term is calculated as (approximately), (X(curr)^2-X(prev)^2)/dt, 
        # (also from the sgs.F90 file) where a general diffusion subroutine is applied 
        # between the "curr" and "prev" values.  Because the DISSIP term more closely resembles
        # CLUBB's dp2 term, however, I've switched them here so that DISSIP = dp2 and DIFFTR = dp1.
        rtpthlp_budget_lines = [
            {'var_names': ['rtpthlp_residual', self.getRtpthlpResidual], 'legend_label': 'rtpthlp_res'},
            {'var_names': ['rtpthlp_adv',self.calc_rtpthlp_adv], 'legend_label': 'rtpthlp_adv'},
            {'var_names': ['rtpthlp_grad',self.calc_rtpthlp_grad], 'legend_label': 'rtpthlp_grad'},
            {'var_names': ['rtpthlp_dp2'], 'legend_label': 'rtpthlp_dissip'},
            {'var_names': ['rtpthlp_dp1'], 'legend_label': 'rtpthlp_difftr'},
            {'var_names': ['rtpthlp_mc'], 'legend_label': 'rtpthlp_prec'},
            {'var_names': ['rtpthlp_forcing',self.calc_rtpthlp_forc], 'legend_label': 'rtpthlp_forc'},
            {'var_names': ['rtpthlp_bt'], 'legend_label': 'rtpthlp_bt'},
            {'var_names': ['rtpthlp_pd'], 'legend_label': 'rtpthlp_limiters'},
            {'var_names': ['rtpthlp_sf'], 'legend_label': 'rtpthlp_sf'},
        ]

        upwp_budget_lines = [
            {'var_names': ['upwp_residual', self.getUpwpResidual], 'legend_label': 'upwp_res'},
            {'var_names': ['upwp_dp1'], 'legend_label': 'upwp_dfsn'},
            {'var_names': ['upwp_adv',self.calc_upwp_adv], 'legend_label': 'upwp_adv'},
            {'var_names': ['upwp_pres',self.calc_upwp_pres], 'legend_label': 'upwp_pres'},
            {'var_names': ['upwp_aniz',self.calc_upwp_aniz], 'legend_label': 'upwp_aniz'},
            {'var_names': ['upwp_bp'], 'legend_label': 'upwp_buoy'},
            {'var_names': ['upwp_tp'], 'legend_label': 'upwp_shear'},
            {'var_names': ['upwp_bt'], 'legend_label': 'upwp_bt'},
            {'var_names': ['upwp_limiters',self.calc_upwp_limiters], 'legend_label': 'upwp_limiters'},
        ]

        vpwp_budget_lines = [
            {'var_names': ['vpwp_residual', self.getVpwpResidual], 'legend_label': 'vpwp_res'},
            {'var_names': ['vpwp_dp1'], 'legend_label': 'vpwp_dfsn'},
            {'var_names': ['vpwp_adv',self.calc_vpwp_adv], 'legend_label': 'vpwp_adv'},
            {'var_names': ['vpwp_pres',self.calc_vpwp_pres], 'legend_label': 'vpwp_pres'},
            {'var_names': ['vpwp_aniz',self.calc_vpwp_aniz], 'legend_label': 'vpwp_aniz'},
            {'var_names': ['vpwp_bp'], 'legend_label': 'vpwp_buoy'},
            {'var_names': ['vpwp_tp'], 'legend_label': 'vpwp_shear'},         
            {'var_names': ['vpwp_bt'], 'legend_label': 'vpwp_bt'},
            {'var_names': ['vpwp_limiters',self.calc_vpwp_limiters], 'legend_label': 'vpwp_limiters'},
        ]

        up2_budget_lines = [
            {'var_names': ['up2_res',self.calc_up2_res], 'legend_label': 'up2_res'},
            {'var_names': ['up2_adv',self.calc_up2_adv], 'legend_label': 'up2_adv'},
            {'var_names': ['up2_tp'], 'legend_label': 'up2_shear'},
            {'var_names': ['up2_pr1'], 'legend_label': 'up2_redis'},
            {'var_names': ['up2_dfsn',self.calc_up2_dfsn], 'legend_label': 'up2_dfsn'},
            {'var_names': ['up2_bt'], 'legend_label': 'up2_bt'},
            {'var_names': ['up2_sdmp'], 'legend_label': 'up2_sdmp'},
            {'var_names': ['up2_limiters',self.calc_up2_limiters], 'legend_label': 'up2_limiters'},
            {'var_names': ['up2_sf'], 'legend_label': 'up2_sf'},
            {'var_names': ['up2_splat'], 'legend_label': 'up2_splat'}
        ]

        vp2_budget_lines = [
            {'var_names': ['vp2_res',self.calc_vp2_res], 'legend_label': 'vp2_res'},
            {'var_names': ['vp2_adv',self.calc_vp2_adv], 'legend_label': 'vp2_adv'},
            {'var_names': ['vp2_tp'], 'legend_label': 'vp2_shear'},
            {'var_names': ['vp2_pr1'], 'legend_label': 'vp2_redis'},
            {'var_names': ['vp2_dfsn',self.calc_vp2_dfsn], 'legend_label': 'vp2_dfsn'},
            {'var_names': ['vp2_bt'], 'legend_label': 'vp2_bt'},
            {'var_names': ['vp2_sdmp'], 'legend_label': 'vp2_sdmp'},
            {'var_names': ['vp2_limiters',self.calc_vp2_limiters], 'legend_label': 'vp2_limiters'},
            {'var_names': ['vp2_sf'], 'legend_label': 'vp2_sf'},
            {'var_names': ['vp2_splat'], 'legend_label': 'vp2_splat'}
        ]

        self.variable_definitions = [
#            {'var_names':
#                {
#                    'clubb': ['rtm'],
#                    'sam': ['rtm'],
#                    'coamps': ['rtm'],
#                    'r408': ['rtm'],
#                    'hoc': ['rtm'],
#                    'e3sm': ['rtm'],
#                    'cam': ['rtm'],
#                    'wrf': ['rtm'],
#                },
#                'lines': rtm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
#                'priority': True,
#            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': ['wpthlp'],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'],
                    'wrf': ['wpthlp'],
                },
                'lines': wpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wprtp'],
                    'sam': ['wprtp'],
                    'coamps': ['wprtp'],
                    'r408': ['wprtp'],
                    'hoc': ['wprtp'],
                    'e3sm': ['wprtp'],
                    'cam': ['wprtp'],
                    'wrf': ['wprtp'],
                },
                'lines': wprtp_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wp2'],
                    'sam': ['wp2'],
                    'coamps': ['wp2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['wp2'],
                    'wrf': ['wp2'],
                },
                'lines': wp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['wp3'],
                    'sam': ['wp3'],
                    'coamps': ['wp3'],
                    'r408': ['wp3'],
                    'hoc': ['wp3'],
                    'e3sm': ['wp3'],
                    'cam': ['wp3'],
                    'wrf': ['wp3'],
                },
                'lines': wp3_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['thlp2'],
                    'sam': ['thlp2'],
                    'coamps': ['thlp2'],
                    'r408': ['thlp2'],
                    'hoc': ['thlp2'],
                    'e3sm': ['thlp2'],
                    'cam': ['thlp2'],
                    'wrf': ['thlp2'],
                },
                'lines': thlp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtp2'],
                    'sam': ['rtp2'],
                    'coamps': ['rtp2'],
                    'r408': ['rtp2'],
                    'hoc': ['rtp2'],
                    'e3sm': ['rtp2'],
                    'cam': ['rtp2'],
                    'wrf': ['rtp2'],
                },
                'lines': rtp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['rtpthlp'],
                    'sam': ['rtpthlp'],
                    'coamps': ['rtpthlp'],
                    'r408': ['rtpthlp'],
                    'hoc': ['rtpthlp'],
                    'e3sm': ['rtpthlp'],
                    'cam': ['rtpthlp'],
                    'wrf': ['rtpthlp'],
                },
                'lines': rtpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['upwp'],
                    'sam': ['upwp'],
                    'coamps': ['upwp'],
                    'r408': ['upwp'],
                    'hoc': ['upwp'],
                    'e3sm': ['upwp'],
                    'cam': ['upwp'],
                    'wrf': ['upwp'],
                },
                'lines': upwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names':
                {
                    'clubb': ['vpwp'],
                    'sam': ['vpwp'],
                    'coamps': ['vpwp'],
                    'r408': ['vpwp'],
                    'hoc': ['vpwp'],
                    'e3sm': ['vpwp'],
                    'cam': ['vpwp'],
                    'wrf': ['vpwp'],
                },
                'lines': vpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names': {
                'clubb': ['up2'],
                'sam': ['up2'],
                'coamps': ['up2'],
                'r408': ['up2'],
                'hoc': ['up2'],
                'e3sm': ['up2'],
                'wrf': ['up2']
            },
                'lines': up2_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
            {'var_names': {
                'clubb': ['vp2'],
                'sam': ['vp2'],
                'coamps': ['vp2'],
                'r408': ['vp2'],
                'hoc': ['vp2'],
                'e3sm': ['vp2'],
                'wrf': ['vp2']
            },
                'lines': vp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
                'priority': True,
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars)

    def getRtmClipping(self, dataset_override=None):
        '''

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

        ``rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        rtm_mfl, indep, dataset = self.getVarForCalculations('rtm_mfl', dataset_override)
        rtm_cl, indep, dataset = self.getVarForCalculations('rtm_cl', dataset)
        rtm_tacl, indep, dataset = self.getVarForCalculations('rtm_tacl', dataset)
        rtm_sdmp, indep, dataset = self.getVarForCalculations('rtm_sdmp', dataset)

        output_data = rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp

        return output_data, indep

    def getRtmForcing(self, dataset_override=None):
        '''

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

        ``rtm_forcing - rtm_mc``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        rtm_mc, indep, dataset = self.getVarForCalculations('rtm_mc', dataset_override)
        rtm_forcing, indep, dataset = self.getVarForCalculations('rtm_forcing', dataset)

        output_data = rtm_forcing - rtm_mc

        return output_data, indep

    def getRtmResidual(self, dataset_override=None):
        '''

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

        ``rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtm_mfl, indep, dataset = self.getVarForCalculations('rtm_mfl', dataset_override)
        rtm_cl, indep, dataset = self.getVarForCalculations('rtm_cl', dataset)
        rtm_tacl, indep, dataset = self.getVarForCalculations('rtm_tacl', dataset)
        rtm_sdmp, indep, dataset = self.getVarForCalculations('rtm_sdmp', dataset)
        rtm_bt, indep, dataset = self.getVarForCalculations('rtm_bt', dataset)
        rtm_ta, indep, dataset = self.getVarForCalculations('rtm_ta', dataset)
        rtm_forcing, indep, dataset = self.getVarForCalculations('rtm_forcing', dataset)
        rtm_pd, indep, dataset = self.getVarForCalculations('rtm_pd', dataset)
        rtm_ma, indep, dataset = self.getVarForCalculations('rtm_ma', dataset)

        output_data = rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)

        return output_data, indep

    def calc_wpthlp_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_ma, indep, dataset = self.getVarForCalculations('wpthlp_ma', dataset_override)
        wpthlp_ta, indep, dataset = self.getVarForCalculations('wpthlp_ta', dataset)
        wpthlp_ac, indep, dataset = self.getVarForCalculations('wpthlp_ac', dataset)

        output_data = wpthlp_ma + wpthlp_ta + wpthlp_ac

        return output_data, indep

    def calc_wpthlp_pres(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_pr1, indep, dataset = self.getVarForCalculations('wpthlp_pr1', dataset_override)
        wpthlp_pr2, indep, dataset = self.getVarForCalculations('wpthlp_pr2', dataset)
        wpthlp_pr3, indep, dataset = self.getVarForCalculations('wpthlp_pr3', dataset)

        output_data = wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3

        return output_data, indep

    def calc_wpthlp_rad(self, dataset_override=None):
        '''
        This function currently outputs an array of zeros since CLUBB does not have a comparable term
        to SAM.  In the future one could replace " - wpthlp_forcing" below with all forcings except
        radiation.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_forcing, indep, dataset = self.getVarForCalculations('wpthlp_forcing', dataset_override)

        output_data = wpthlp_forcing - wpthlp_forcing 

        return output_data, indep

    def calc_wpthlp_forc(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_mc, indep, dataset = self.getVarForCalculations('wpthlp_mc', dataset_override)
        wpthlp_forcing, indep, dataset = self.getVarForCalculations('wpthlp_forcing', dataset)

        output_data = wpthlp_forcing - wpthlp_mc

        return output_data, indep

    def calc_wpthlp_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_mfl, indep, dataset = self.getVarForCalculations('wpthlp_mfl', dataset_override)
        wpthlp_cl, indep, dataset = self.getVarForCalculations('wpthlp_cl', dataset)
        wpthlp_sicl, indep, dataset = self.getVarForCalculations('wpthlp_sicl', dataset)

        output_data = wpthlp_mfl + wpthlp_cl + wpthlp_sicl

        return output_data, indep

    def getWpthlpResidual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 +
            wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wpthlp_mfl, indep, dataset = self.getVarForCalculations('wpthlp_mfl', dataset_override)
        wpthlp_cl, indep, dataset = self.getVarForCalculations('wpthlp_cl', dataset)
        wpthlp_tp, indep, dataset = self.getVarForCalculations('wpthlp_tp', dataset)
        wpthlp_ac, indep, dataset = self.getVarForCalculations('wpthlp_ac', dataset)
        wpthlp_pr1, indep, dataset = self.getVarForCalculations('wpthlp_pr1', dataset)
        wpthlp_pr3, indep, dataset = self.getVarForCalculations('wpthlp_pr3', dataset)
        wpthlp_pr2, indep, dataset = self.getVarForCalculations('wpthlp_pr2', dataset)
        wpthlp_dp1, indep, dataset = self.getVarForCalculations('wpthlp_dp1', dataset)
        wpthlp_sicl, indep, dataset = self.getVarForCalculations('wpthlp_sicl', dataset)
        wpthlp_bt, indep, dataset = self.getVarForCalculations('wpthlp_bt', dataset)
        wpthlp_ta, indep, dataset = self.getVarForCalculations('wpthlp_ta', dataset)
        wpthlp_forcing, indep, dataset = self.getVarForCalculations('wpthlp_forcing', dataset)
        wpthlp_bp, indep, dataset = self.getVarForCalculations('wpthlp_bp', dataset)
        wpthlp_ma, indep, dataset = self.getVarForCalculations('wpthlp_ma', dataset)

        output_data = wpthlp_bt - (
                wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 +
                wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)

        return output_data, indep

    def calc_wprtp_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wprtp_ma, indep, dataset = self.getVarForCalculations('wprtp_ma', dataset_override)
        wprtp_ta, indep, dataset = self.getVarForCalculations('wprtp_ta', dataset)
        wprtp_ac, indep, dataset = self.getVarForCalculations('wprtp_ac', dataset)

        output_data = wprtp_ma + wprtp_ta + wprtp_ac 

        return output_data, indep

    def calc_wprtp_pres(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wprtp_pr1, indep, dataset = self.getVarForCalculations('wprtp_pr1', dataset_override)
        wprtp_pr2, indep, dataset = self.getVarForCalculations('wprtp_pr2', dataset)
        wprtp_pr3, indep, dataset = self.getVarForCalculations('wprtp_pr3', dataset)

        output_data = wprtp_pr1 + wprtp_pr2 + wprtp_pr3

        return output_data, indep

    def calc_wprtp_forc(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wprtp_mc, indep, dataset = self.getVarForCalculations('wprtp_mc', dataset_override)
        wprtp_forcing, indep, dataset = self.getVarForCalculations('wprtp_forcing', dataset)

        output_data = wprtp_forcing - wprtp_mc

        return output_data, indep

    def calc_wprtp_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wprtp_mfl, indep, dataset = self.getVarForCalculations('wprtp_mfl', dataset_override)
        wprtp_cl, indep, dataset = self.getVarForCalculations('wprtp_cl', dataset)
        wprtp_sicl, indep, dataset = self.getVarForCalculations('wprtp_sicl', dataset)

        output_data = wprtp_mfl + wprtp_cl + wprtp_sicl

        return output_data, indep

    def getWprtpResidual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 +
            wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wprtp_mfl, indep, dataset = self.getVarForCalculations('wprtp_mfl', dataset_override)
        wprtp_cl, indep, dataset = self.getVarForCalculations('wprtp_cl', dataset)
        wprtp_tp, indep, dataset = self.getVarForCalculations('wprtp_tp', dataset)
        wprtp_ac, indep, dataset = self.getVarForCalculations('wprtp_ac', dataset)
        wprtp_pr1, indep, dataset = self.getVarForCalculations('wprtp_pr1', dataset)
        wprtp_pr3, indep, dataset = self.getVarForCalculations('wprtp_pr3', dataset)
        wprtp_pr2, indep, dataset = self.getVarForCalculations('wprtp_pr2', dataset)
        wprtp_dp1, indep, dataset = self.getVarForCalculations('wprtp_dp1', dataset)
        wprtp_sicl, indep, dataset = self.getVarForCalculations('wprtp_sicl', dataset)
        wprtp_bt, indep, dataset = self.getVarForCalculations('wprtp_bt', dataset)
        wprtp_ta, indep, dataset = self.getVarForCalculations('wprtp_ta', dataset)
        wprtp_forcing, indep, dataset = self.getVarForCalculations('wprtp_forcing', dataset)
        wprtp_bp, indep, dataset = self.getVarForCalculations('wprtp_bp', dataset)
        wprtp_ma, indep, dataset = self.getVarForCalculations('wprtp_ma', dataset)
        wprtp_pd, indep, dataset = self.getVarForCalculations('wprtp_pd', dataset)

        output_data = wprtp_bt - (
                wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 +
                wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_pd + wprtp_forcing)

        return output_data, indep

    def calc_wp2_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp2_ma, indep, dataset = self.getVarForCalculations('wp2_ma', dataset_override)
        wp2_ta, indep, dataset = self.getVarForCalculations('wp2_ta', dataset)
        wp2_ac, indep, dataset = self.getVarForCalculations('wp2_ac', dataset)

        output_data = wp2_ma + wp2_ta + wp2_ac

        return output_data, indep

    def calc_wp2_pres(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp2_pr2, indep, dataset = self.getVarForCalculations('wp2_pr2', dataset_override)
        wp2_pr3, indep, dataset = self.getVarForCalculations('wp2_pr3', dataset)
        wp2_splat, indep, dataset = self.getVarForCalculations('wp2_splat', dataset)
        wp2_pr_dfsn, indep, dataset = self.getVarForCalculations('wp2_pr_dfsn', dataset)

        output_data = wp2_pr2 + wp2_pr3 + wp2_splat + wp2_pr_dfsn

        return output_data, indep

    def calc_wp2_dfsn(self, dataset_override=None):
         '''
         '''
         # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
         wp2_dp1, indep, dataset = self.getVarForCalculations('wp2_dp1', dataset_override)
         wp2_dp2, indep, dataset = self.getVarForCalculations('wp2_dp2', dataset)

         output_data = wp2_dp1 + wp2_dp2

         return output_data, indep

    def calc_wp2_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp2_cl, indep, dataset = self.getVarForCalculations('wp2_cl', dataset_override)
        wp2_pd, indep, dataset = self.getVarForCalculations('wp2_pd', dataset)

        output_data = wp2_cl + wp2_pd

        return output_data, indep

    def getWp2Residual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            wp2_bt - (wp2_ma + wp2_ta + wp2_tp + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 +
            wp2_mfl + wp2_cl + wp2_sicl + wp2_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp2_sf, indep, dataset = self.getVarForCalculations('wp2_sf', dataset_override)
        wp2_cl, indep, dataset = self.getVarForCalculations('wp2_cl', dataset)
        wp2_ac, indep, dataset = self.getVarForCalculations('wp2_ac', dataset)
        wp2_pr1, indep, dataset = self.getVarForCalculations('wp2_pr1', dataset)
        wp2_pr3, indep, dataset = self.getVarForCalculations('wp2_pr3', dataset)
        wp2_pr2, indep, dataset = self.getVarForCalculations('wp2_pr2', dataset)
        wp2_pr_dfsn, indep, dataset = self.getVarForCalculations('wp2_pr_dfsn', dataset)
        wp2_dp1, indep, dataset = self.getVarForCalculations('wp2_dp1', dataset)
        wp2_dp2, indep, dataset = self.getVarForCalculations('wp2_dp2', dataset)
        wp2_bt, indep, dataset = self.getVarForCalculations('wp2_bt', dataset)
        wp2_ta, indep, dataset = self.getVarForCalculations('wp2_ta', dataset)
        wp2_splat, indep, dataset = self.getVarForCalculations('wp2_splat', dataset)
        wp2_bp, indep, dataset = self.getVarForCalculations('wp2_bp', dataset)
        wp2_ma, indep, dataset = self.getVarForCalculations('wp2_ma', dataset)
        wp2_pd, indep, dataset = self.getVarForCalculations('wp2_pd', dataset)

        output_data = wp2_bt - (
                wp2_ma + wp2_ta + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_dp2 +
                wp2_cl + wp2_pd + wp2_sf + wp2_splat + wp2_pr_dfsn )

        return output_data, indep

    def calc_wp3_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp3_ma, indep, dataset = self.getVarForCalculations('wp3_ma', dataset_override)
        wp3_ta, indep, dataset = self.getVarForCalculations('wp3_ta', dataset)
        wp3_tp, indep, dataset = self.getVarForCalculations('wp3_tp', dataset)
        wp3_ac, indep, dataset = self.getVarForCalculations('wp3_ac', dataset)

        output_data = wp3_ma + wp3_ta + wp3_tp + wp3_ac

        return output_data, indep
    
    def calc_wp3_pres(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp3_pr1, indep, dataset = self.getVarForCalculations('wp3_pr1', dataset_override)
        wp3_pr2, indep, dataset = self.getVarForCalculations('wp3_pr2', dataset)
        wp3_pr3, indep, dataset = self.getVarForCalculations('wp3_pr3', dataset)
        wp3_pr_tp, indep, dataset = self.getVarForCalculations('wp3_pr_tp', dataset)
        wp3_pr_turb, indep, dataset = self.getVarForCalculations('wp3_pr_turb', dataset)
        wp3_pr_dfsn, indep, dataset = self.getVarForCalculations('wp3_pr_dfsn', dataset)
        wp3_splat, indep, dataset = self.getVarForCalculations('wp3_splat', dataset)

        output_data = wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_pr_turb + wp3_pr_dfsn + wp3_splat + wp3_pr_tp

        return output_data, indep

    def getWp3Residual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_pr_turb + wp3_pr1 + wp3_pr2 + wp3_pr3 +
            wp3_dp1 + wp3_cl+wp3_splat)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp3_bp1, indep, dataset = self.getVarForCalculations('wp3_bp1', dataset_override)
        wp3_cl, indep, dataset = self.getVarForCalculations('wp3_cl', dataset)
        wp3_ac, indep, dataset = self.getVarForCalculations('wp3_ac', dataset)
        wp3_pr1, indep, dataset = self.getVarForCalculations('wp3_pr1', dataset)
        wp3_pr2, indep, dataset = self.getVarForCalculations('wp3_pr2', dataset)
        wp3_pr3, indep, dataset = self.getVarForCalculations('wp3_pr3', dataset)
        wp3_pr_turb, indep, dataset = self.getVarForCalculations('wp3_pr_turb', dataset)
        wp3_pr_dfsn, indep, dataset = self.getVarForCalculations('wp3_pr_dfsn', dataset)
        wp3_pr_tp, indep, dataset = self.getVarForCalculations('wp3_pr_tp', dataset)
        wp3_dp1, indep, dataset = self.getVarForCalculations('wp3_dp1', dataset)
        wp3_bt, indep, dataset = self.getVarForCalculations('wp3_bt', dataset)
        wp3_ta, indep, dataset = self.getVarForCalculations('wp3_ta', dataset)
        wp3_splat, indep, dataset = self.getVarForCalculations('wp3_splat', dataset)
        wp3_ma, indep, dataset = self.getVarForCalculations('wp3_ma', dataset)
        wp3_tp, indep, dataset = self.getVarForCalculations('wp3_tp', dataset)

        output_data = wp3_bt - (
                wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_pr_turb + wp3_pr_dfsn 
                + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_pr_tp
                + wp3_dp1 + wp3_cl + wp3_splat)

        return output_data, indep

    def calc_up2_res(self, dataset_override=None):

        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        up2_sf, indep, dataset = self.getVarForCalculations('up2_sf', dataset_override)
        up2_cl, indep, dataset = self.getVarForCalculations('up2_cl', dataset)
        up2_pd, indep, dataset = self.getVarForCalculations('up2_pd', dataset)
        up2_pr1, indep, dataset = self.getVarForCalculations('up2_pr1', dataset)
        up2_pr2, indep, dataset = self.getVarForCalculations('up2_pr2', dataset)
        up2_dp1, indep, dataset = self.getVarForCalculations('up2_dp1', dataset)
        up2_dp2, indep, dataset = self.getVarForCalculations('up2_dp2', dataset)
        up2_ta, indep, dataset = self.getVarForCalculations('up2_ta', dataset)
        up2_splat, indep, dataset = self.getVarForCalculations('up2_splat', dataset)
        up2_ma, indep, dataset = self.getVarForCalculations('up2_ma', dataset)
        up2_tp, indep, dataset = self.getVarForCalculations('up2_tp', dataset)
        up2_bt, indep, dataset = self.getVarForCalculations('up2_bt', dataset)

        output_data = up2_bt - (
                up2_ma + up2_ta + up2_tp + up2_dp1 + up2_dp2 + up2_pr1 + up2_pr2 + up2_pd +
                up2_sf + up2_cl + up2_splat)

        return output_data, indep

    def calc_up2_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        up2_ma, indep, dataset = self.getVarForCalculations('up2_ma', dataset_override)
        up2_ta, indep, dataset = self.getVarForCalculations('up2_ta', dataset)

        output_data = up2_ma + up2_ta

        return output_data, indep

    def calc_up2_dfsn(self, dataset_override=None):
         '''
         '''
         # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
         up2_dp1, indep, dataset = self.getVarForCalculations('up2_dp1', dataset_override)
         up2_dp2, indep, dataset = self.getVarForCalculations('up2_dp2', dataset)

         output_data = up2_dp1 + up2_dp2

         return output_data, indep

    def calc_up2_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        up2_cl, indep, dataset = self.getVarForCalculations('up2_cl', dataset_override)
        up2_pd, indep, dataset = self.getVarForCalculations('up2_pd', dataset)

        output_data = up2_cl + up2_pd

        return output_data, indep

    def calc_vp2_res(self, dataset_override=None):

        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vp2_sf, indep, dataset = self.getVarForCalculations('vp2_sf', dataset_override)
        vp2_cl, indep, dataset = self.getVarForCalculations('vp2_cl', dataset)
        vp2_pd, indep, dataset = self.getVarForCalculations('vp2_pd', dataset)
        vp2_pr1, indep, dataset = self.getVarForCalculations('vp2_pr1', dataset)
        vp2_pr2, indep, dataset = self.getVarForCalculations('vp2_pr2', dataset)
        vp2_dp1, indep, dataset = self.getVarForCalculations('vp2_dp1', dataset)
        vp2_dp2, indep, dataset = self.getVarForCalculations('vp2_dp2', dataset)
        vp2_ta, indep, dataset = self.getVarForCalculations('vp2_ta', dataset)
        vp2_splat, indep, dataset = self.getVarForCalculations('vp2_splat', dataset)
        vp2_ma, indep, dataset = self.getVarForCalculations('vp2_ma', dataset)
        vp2_tp, indep, dataset = self.getVarForCalculations('vp2_tp', dataset)
        vp2_bt, indep, dataset = self.getVarForCalculations('vp2_bt', dataset)

        output_data = vp2_bt - (
                vp2_ma + vp2_ta + vp2_tp + vp2_dp1 + vp2_dp2 + vp2_pr1 + vp2_pr2 + vp2_pd +
                vp2_sf + vp2_cl + vp2_splat)

        return output_data, indep

    def calc_vp2_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vp2_ma, indep, dataset = self.getVarForCalculations('vp2_ma', dataset_override)
        vp2_ta, indep, dataset = self.getVarForCalculations('vp2_ta', dataset)
        vp2_tp, indep, dataset = self.getVarForCalculations('vp2_tp', dataset)

        output_data = vp2_ma + vp2_ta + vp2_tp

        return output_data, indep

    def calc_vp2_dfsn(self, dataset_override=None):
         '''
         '''
         # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
         vp2_dp1, indep, dataset = self.getVarForCalculations('vp2_dp1', dataset_override)
         vp2_dp2, indep, dataset = self.getVarForCalculations('vp2_dp2', dataset)

         output_data = vp2_dp1 + vp2_dp2

         return output_data, indep

    def calc_vp2_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vp2_cl, indep, dataset = self.getVarForCalculations('vp2_cl', dataset_override)
        vp2_pd, indep, dataset = self.getVarForCalculations('vp2_pd', dataset)

        output_data = vp2_cl + vp2_pd

        return output_data, indep

    def calc_thlp2_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        thlp2_ma, indep, dataset = self.getVarForCalculations('thlp2_ma', dataset_override)
        thlp2_ta, indep, dataset = self.getVarForCalculations('thlp2_ta', dataset)

        output_data = thlp2_ma + thlp2_ta

        return output_data, indep

    def calc_thlp2_rad(self, dataset_override=None):
        '''
        This function currently outputs an array of zeros since CLUBB does not have a comparable term
        to SAM.  In the future one could replace " - thlp2_forcing" below with all forcings except
        radiation.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        thlp2_forcing, indep, dataset = self.getVarForCalculations('thlp2_forcing', dataset_override)

        output_data = thlp2_forcing - thlp2_forcing

        return output_data, indep

    def calc_thlp2_forc(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        thlp2_mc, indep, dataset = self.getVarForCalculations('thlp2_mc', dataset_override)
        thlp2_forcing, indep, dataset = self.getVarForCalculations('thlp2_forcing', dataset)

        output_data = thlp2_forcing - thlp2_mc

        return output_data, indep

    def calc_thlp2_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        thlp2_cl, indep, dataset = self.getVarForCalculations('thlp2_cl', dataset_override)
        thlp2_pd, indep, dataset = self.getVarForCalculations('thlp2_pd', dataset)

        output_data = thlp2_cl + thlp2_pd

        return output_data, indep

    def getThlp2Residual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd +
            thlp2_sf + thlp2_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        thlp2_cl, indep, dataset = self.getVarForCalculations('thlp2_cl', dataset_override)
        thlp2_dp2, indep, dataset = self.getVarForCalculations('thlp2_dp2', dataset)
        thlp2_forcing, indep, dataset = self.getVarForCalculations('thlp2_forcing', dataset)
        thlp2_sf, indep, dataset = self.getVarForCalculations('thlp2_sf', dataset)
        thlp2_dp1, indep, dataset = self.getVarForCalculations('thlp2_dp1', dataset)
        thlp2_bt, indep, dataset = self.getVarForCalculations('thlp2_bt', dataset)
        thlp2_ta, indep, dataset = self.getVarForCalculations('thlp2_ta', dataset)
        thlp2_pd, indep, dataset = self.getVarForCalculations('thlp2_pd', dataset)
        thlp2_ma, indep, dataset = self.getVarForCalculations('thlp2_ma', dataset)
        thlp2_tp, indep, dataset = self.getVarForCalculations('thlp2_tp', dataset)

        output_data = thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 +
                                  thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)

        return output_data, indep

    def calc_rtp2_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtp2_ma, indep, dataset = self.getVarForCalculations('rtp2_ma', dataset_override)
        rtp2_ta, indep, dataset = self.getVarForCalculations('rtp2_ta', dataset)

        output_data = rtp2_ma + rtp2_ta

        return output_data, indep

    def calc_rtp2_forc(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtp2_mc, indep, dataset = self.getVarForCalculations('rtp2_mc', dataset_override)
        rtp2_forcing, indep, dataset = self.getVarForCalculations('rtp2_forcing', dataset)

        output_data = rtp2_forcing - rtp2_mc

        return output_data, indep

    def calc_rtp2_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtp2_cl, indep, dataset = self.getVarForCalculations('rtp2_cl', dataset_override)
        rtp2_pd, indep, dataset = self.getVarForCalculations('rtp2_pd', dataset)

        output_data = rtp2_cl + rtp2_pd

        return output_data, indep

    def getRtp2Residual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtp2_cl, indep, dataset = self.getVarForCalculations('rtp2_cl', dataset_override)
        rtp2_dp2, indep, dataset = self.getVarForCalculations('rtp2_dp2', dataset)
        rtp2_forcing, indep, dataset = self.getVarForCalculations('rtp2_forcing', dataset)
        rtp2_sf, indep, dataset = self.getVarForCalculations('rtp2_sf', dataset)
        rtp2_dp1, indep, dataset = self.getVarForCalculations('rtp2_dp1', dataset)
        rtp2_bt, indep, dataset = self.getVarForCalculations('rtp2_bt', dataset)
        rtp2_ta, indep, dataset = self.getVarForCalculations('rtp2_ta', dataset)
        rtp2_pd, indep, dataset = self.getVarForCalculations('rtp2_pd', dataset)
        rtp2_ma, indep, dataset = self.getVarForCalculations('rtp2_ma', dataset)
        rtp2_tp, indep, dataset = self.getVarForCalculations('rtp2_tp', dataset)

        output_data = rtp2_bt - (
                rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)

        return output_data, indep

    def calc_rtpthlp_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtpthlp_ma, indep, dataset = self.getVarForCalculations('rtpthlp_ma', dataset_override)
        rtpthlp_ta, indep, dataset = self.getVarForCalculations('rtpthlp_ta', dataset)

        output_data = rtpthlp_ma + rtpthlp_ta

        return output_data, indep

    def calc_rtpthlp_grad(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtpthlp_tp1, indep, dataset = self.getVarForCalculations('rtpthlp_tp1', dataset_override)
        rtpthlp_tp2, indep, dataset = self.getVarForCalculations('rtpthlp_tp2', dataset)

        output_data = rtpthlp_tp1 + rtpthlp_tp2

        return output_data, indep    
    
    def calc_rtpthlp_forc(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtpthlp_mc, indep, dataset = self.getVarForCalculations('rtpthlp_mc', dataset_override)
        rtpthlp_forcing, indep, dataset = self.getVarForCalculations('rtpthlp_forcing', dataset)

        output_data = rtpthlp_forcing - rtpthlp_mc

        return output_data, indep

    def getRtpthlpResidual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl +
            rtpthlp_pd + rtpthlp_sf + rtpthlp_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rtpthlp_cl, indep, dataset = self.getVarForCalculations('rtpthlp_cl', dataset_override)
        rtpthlp_dp2, indep, dataset = self.getVarForCalculations('rtpthlp_dp2', dataset)
        rtpthlp_forcing, indep, dataset = self.getVarForCalculations('rtpthlp_forcing', dataset)
        rtpthlp_sf, indep, dataset = self.getVarForCalculations('rtpthlp_sf', dataset)
        rtpthlp_dp1, indep, dataset = self.getVarForCalculations('rtpthlp_dp1', dataset)
        rtpthlp_bt, indep, dataset = self.getVarForCalculations('rtpthlp_bt', dataset)
        rtpthlp_ta, indep, dataset = self.getVarForCalculations('rtpthlp_ta', dataset)
        rtpthlp_tp2, indep, dataset = self.getVarForCalculations('rtpthlp_tp2', dataset)
        rtpthlp_ma, indep, dataset = self.getVarForCalculations('rtpthlp_ma', dataset)
        rtpthlp_tp1, indep, dataset = self.getVarForCalculations('rtpthlp_tp1', dataset)

        output_data = rtpthlp_bt - (
                rtpthlp_ma + rtpthlp_ta + rtpthlp_tp1 + rtpthlp_tp2 + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl +
                rtpthlp_sf + rtpthlp_forcing)

        return output_data, indep

    def calc_upwp_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        upwp_ma, indep, dataset = self.getVarForCalculations('upwp_ma', dataset_override)
        upwp_ta, indep, dataset = self.getVarForCalculations('upwp_ta', dataset)
        upwp_ac, indep, dataset = self.getVarForCalculations('upwp_ac', dataset)

        output_data = upwp_ma + upwp_ta + upwp_ac

        return output_data, indep

    def calc_upwp_pres(self, dataset_override=None):
        '''
        This function currently outputs an array of zeros since CLUBB does not have a comparable term
        to SAM.  In the future one could replace this with an appropriate pressure term.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        upwp_pr2, indep, dataset = self.getVarForCalculations('upwp_pr2', dataset_override)

        output_data = upwp_pr2 - upwp_pr2

        return output_data, indep

    def calc_upwp_aniz(self, dataset_override=None):
        '''
        In SAM, upwp_aniz ends up being equal to (-w'*dp'/dx-u'*dp'/dz)-d(u'p')/dz.
        See around lines 990+ in statistics.f90 and the lines calculating presx.
        This is similar to the "scrambling" part of the pressure covariance,
        since it's the total pressure contribution minus a diffusion part.  (It's not 
        subtracting the total diffusion part, since there would also be d(w'p')/dz.)
        Hence I have set upwp_aniz equal to the sum of all current pressure contributions,
        since they all parameterize the pressure-scrambling part of the equation.
        As of now CLUBB doesn't explicitly model the "diffusion" part of the pressure
        covariance upwp_pres, so that part is set to zero.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        upwp_pr1, indep, dataset = self.getVarForCalculations('upwp_pr1', dataset_override)
        upwp_pr2, indep, dataset = self.getVarForCalculations('upwp_pr2', dataset)
        upwp_pr3, indep, dataset = self.getVarForCalculations('upwp_pr3', dataset)
        upwp_pr4, indep, dataset = self.getVarForCalculations('upwp_pr4', dataset)

        output_data = upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4

        return output_data, indep

    def calc_upwp_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        upwp_cl, indep, dataset = self.getVarForCalculations('upwp_cl', dataset_override)
        upwp_mfl, indep, dataset = self.getVarForCalculations('upwp_mfl', dataset)

        output_data = upwp_cl + upwp_mfl

        return output_data, indep

    def getUpwpResidual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_dp1 + upwp_dp2 + upwp_cl + upwp_pd + upwp_sf + upwp_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        upwp_cl, indep, dataset = self.getVarForCalculations('upwp_cl', dataset_override)
        upwp_tp, indep, dataset = self.getVarForCalculations('upwp_tp', dataset)
        upwp_ac, indep, dataset = self.getVarForCalculations('upwp_ac', dataset)
        upwp_bp, indep, dataset = self.getVarForCalculations('upwp_bp', dataset)
        upwp_dp1, indep, dataset = self.getVarForCalculations('upwp_dp1', dataset)
        upwp_bt, indep, dataset = self.getVarForCalculations('upwp_bt', dataset)
        upwp_ta, indep, dataset = self.getVarForCalculations('upwp_ta', dataset)
        upwp_pr1, indep, dataset = self.getVarForCalculations('upwp_pr1', dataset)
        upwp_pr2, indep, dataset = self.getVarForCalculations('upwp_pr2', dataset)
        upwp_pr3, indep, dataset = self.getVarForCalculations('upwp_pr3', dataset)
        upwp_pr4, indep, dataset = self.getVarForCalculations('upwp_pr4', dataset)
        upwp_mfl, indep, dataset = self.getVarForCalculations('upwp_mfl', dataset)
        upwp_ma, indep, dataset = self.getVarForCalculations('upwp_ma', dataset)

        output_data = upwp_bt - (
                upwp_ma + upwp_ta + upwp_tp + upwp_ac + upwp_bp + upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4 +
                upwp_dp1 + upwp_mfl + upwp_cl)

        return output_data, indep

    def calc_vpwp_adv(self, dataset_override=None):
        '''
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vpwp_ma, indep, dataset = self.getVarForCalculations('vpwp_ma', dataset_override)
        vpwp_ta, indep, dataset = self.getVarForCalculations('vpwp_ta', dataset)
        vpwp_ac, indep, dataset = self.getVarForCalculations('vpwp_ac', dataset)

        output_data = vpwp_ma + vpwp_ta + vpwp_ac

        return output_data, indep

    def calc_vpwp_pres(self, dataset_override=None):
        '''
        This function currently outputs an array of zeros since CLUBB does not have a comparable term
        to SAM.  In the future one could replace this with an appropriate pressure term.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vpwp_pr2, indep, dataset = self.getVarForCalculations('vpwp_pr2', dataset_override)

        output_data = vpwp_pr2 - vpwp_pr2

        return output_data, indep

    def calc_vpwp_aniz(self, dataset_override=None):
        '''
        In SAM, vpwp_aniz ends up being equal to (-w'*dp'/dy-v'*dp'/dz)-d(v'p')/dz.
        See around lines 990+ in statistics.f90 and the lines calculating presy.
        This is similar to the "scrambling" part of the pressure covariance,
        since it's the total pressure contribution minus a diffusion part.  (It's not 
        subtracting the total diffusion part, since there would also be d(w'p')/dz.)
        Hence I have set vpwp_aniz equal to the sum of all current pressure contributions,
        since they all parameterize the pressure-scrambling part of the equation.
        As of now CLUBB doesn't explicitly model the "diffusion" part of the pressure
        covariance vpwp_pres, so that part is set to zero.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vpwp_pr1, indep, dataset = self.getVarForCalculations('vpwp_pr1', dataset_override)
        vpwp_pr2, indep, dataset = self.getVarForCalculations('vpwp_pr2', dataset)
        vpwp_pr3, indep, dataset = self.getVarForCalculations('vpwp_pr3', dataset)
        vpwp_pr4, indep, dataset = self.getVarForCalculations('vpwp_pr4', dataset)

        output_data = vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4

        return output_data, indep

    def calc_vpwp_limiters(self, dataset_override=None):
        '''
        This term includes limiters, i.e. various types of clipping to prevent unwanted values.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vpwp_cl, indep, dataset = self.getVarForCalculations('vpwp_cl', dataset_override)
        vpwp_mfl, indep, dataset = self.getVarForCalculations('vpwp_mfl', dataset)

        output_data = vpwp_cl + vpwp_mfl

        return output_data, indep

    def getVpwpResidual(self, dataset_override=None):
        '''

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

        .. code-block:: python
            :linenos:

            vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_dp1 + vpwp_dp2 + vpwp_cl + vpwp_pd + vpwp_sf + vpwp_forcing)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        vpwp_cl, indep, dataset = self.getVarForCalculations('vpwp_cl', dataset_override)
        vpwp_tp, indep, dataset = self.getVarForCalculations('vpwp_tp', dataset)
        vpwp_ac, indep, dataset = self.getVarForCalculations('vpwp_ac', dataset)
        vpwp_bp, indep, dataset = self.getVarForCalculations('vpwp_bp', dataset)
        vpwp_dp1, indep, dataset = self.getVarForCalculations('vpwp_dp1', dataset)
        vpwp_bt, indep, dataset = self.getVarForCalculations('vpwp_bt', dataset)
        vpwp_ta, indep, dataset = self.getVarForCalculations('vpwp_ta', dataset)
        vpwp_pr1, indep, dataset = self.getVarForCalculations('vpwp_pr1', dataset)
        vpwp_pr2, indep, dataset = self.getVarForCalculations('vpwp_pr2', dataset)
        vpwp_pr3, indep, dataset = self.getVarForCalculations('vpwp_pr3', dataset)
        vpwp_pr4, indep, dataset = self.getVarForCalculations('vpwp_pr4', dataset)
        vpwp_mfl, indep, dataset = self.getVarForCalculations('vpwp_mfl', dataset)
        vpwp_ma, indep, dataset = self.getVarForCalculations('vpwp_ma', dataset)

        output_data = vpwp_bt - (
                vpwp_ma + vpwp_ta + vpwp_tp + vpwp_ac + vpwp_bp + vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4 +
                vpwp_dp1 + vpwp_mfl + vpwp_cl)

        return output_data, indep
