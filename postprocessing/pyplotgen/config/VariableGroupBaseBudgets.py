"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupBaseBudgets(VariableGroup):
    """

    """

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_dataset=None):

        self.name = "base variables budgets"
        
        thlm_budget_lines = [
            {'aliases': ['thlm_bt'], 'legend_label': 'thlm_bt'},
            {'aliases': ['thlm_ma'], 'legend_label': 'thlm_ma'},
            {'aliases': ['thlm_ta'], 'legend_label': 'thlm_ta'},
            {'aliases': ['thlm_mc'], 'legend_label': 'thlm_mc'},
            {'aliases': ['thlm_clipping'], 'legend_label': 'thlm_clipping', 'fallback_func': self.getThlmClipping},
            {'aliases': ['radht'], 'legend_label': 'radht'},
            {'aliases': ['lsforcing'], 'legend_label': 'lsforcing', 'fallback_func': self.getLsforcing},
            {'aliases': ['thlm_residual'], 'legend_label': 'thlm_residual', 'fallback_func': self.getThlmResidual},

        ]
        
        rtm_budget_lines = [
        {'aliases': ['rtm_bt'], 'legend_label': 'rtm_bt'},
            {'aliases': ['rtm_ma'], 'legend_label': 'rtm_ma'},
            {'aliases': ['rtm_ta'], 'legend_label': 'rtm_ta'},
            {'aliases': ['rtm_mc'], 'legend_label': 'rtm_mc'},
            {'aliases': ['rtm_clipping'], 'legend_label': 'rtm_bt', 'fallback_func': self.getRtmClipping},
            {'aliases': ['rtm_pd'], 'legend_label': 'rtm_pd'},
            {'aliases': ['rtm_forcing'], 'legend_label': 'rtm_forcing', 'fallback_func': self.getRtmForcing},
            {'aliases': ['rtm_residual'], 'legend_label': 'rtm_residual', 'fallback_func': self.getRtmResidual},

        ]

        wpthlp_budget_lines = [
            {'aliases': ['wpthlp_bt'], 'legend_label': 'wpthlp_bt'},
            {'aliases': ['wpthlp_ma'], 'legend_label': 'wpthlp_ma'},
            {'aliases': ['wpthlp_ta'], 'legend_label': 'wpthlp_ta'},
            {'aliases': ['wpthlp_tp'], 'legend_label': 'wpthlp_tp'},
            {'aliases': ['wpthlp_ac'], 'legend_label': 'wpthlp_ac'},
            {'aliases': ['wpthlp_bp'], 'legend_label': 'wpthlp_bp'},
            {'aliases': ['wpthlp_pr1'], 'legend_label': 'wpthlp_pr1'},
            {'aliases': ['wpthlp_pr2'], 'legend_label': 'wpthlp_pr2'},
            {'aliases': ['wpthlp_pr3'], 'legend_label': 'wpthlp_pr3'},
            {'aliases': ['wpthlp_dp1'], 'legend_label': 'wpthlp_dp1'},
            {'aliases': ['wpthlp_mfl'], 'legend_label': 'wpthlp_mfl'},
            {'aliases': ['wpthlp_cl'], 'legend_label': 'wpthlp_cl'},
            {'aliases': ['wpthlp_sicl'], 'legend_label': 'wpthlp_sicl'},
            {'aliases': ['wpthlp_forcing'], 'legend_label': 'wpthlp_forcing'},
            {'aliases': ['wpthlp_residual'], 'legend_label': 'wpthlp_residual', 'fallback_func': self.getWpthlpResidual},

        ]

        wprtp_budget_lines = [
            {'aliases': ['wprtp_bt'], 'legend_label': 'wprtp_bt'},
            {'aliases': ['wprtp_ma'], 'legend_label': 'wprtp_ma'},
            {'aliases': ['wprtp_ta'], 'legend_label': 'wprtp_ta'},
            {'aliases': ['wprtp_tp'], 'legend_label': 'wprtp_tp'},
            {'aliases': ['wprtp_ac'], 'legend_label': 'wprtp_ac'},
            {'aliases': ['wprtp_bp'], 'legend_label': 'wprtp_bp'},
            {'aliases': ['wprtp_pr1'], 'legend_label': 'wprtp_pr1'},
            {'aliases': ['wprtp_pr2'], 'legend_label': 'wprtp_pr2'},
            {'aliases': ['wprtp_pr3'], 'legend_label': 'wprtp_pr3'},
            {'aliases': ['wprtp_dp1'], 'legend_label': 'wprtp_dp1'},
            {'aliases': ['wprtp_mfl'], 'legend_label': 'wprtp_mfl'},
            {'aliases': ['wprtp_cl'], 'legend_label': 'wprtp_cl'},
            {'aliases': ['wprtp_sicl'], 'legend_label': 'wprtp_sicl'},
            {'aliases': ['wprtp_forcing'], 'legend_label': 'wprtp_forcing'},
            {'aliases': ['wprtp_residual'], 'legend_label': 'wprtp_residual', 'fallback_func': self.getWprtpResidual},
        ]

        wp2_budget_lines = [
            {'aliases': ['wp2_bt'], 'legend_label': 'wp2_bt'},
            {'aliases': ['wp2_ma'], 'legend_label': 'wp2_ma'},
            {'aliases': ['wp2_ta'], 'legend_label': 'wp2_ta'},
            {'aliases': ['wp2_ac'], 'legend_label': 'wp2_ac'},
            {'aliases': ['wp2_bp'], 'legend_label': 'wp2_bp'},
            {'aliases': ['wp2_pr1'], 'legend_label': 'wp2_pr1'},
            {'aliases': ['wp2_pr2'], 'legend_label': 'wp2_pr2'},
            {'aliases': ['wp2_pr3'], 'legend_label': 'wp2_pr3'},
            {'aliases': ['wp2_dp1'], 'legend_label': 'wp2_dp1'},
            {'aliases': ['wp2_dp2'], 'legend_label': 'wp2_dp2'},
            {'aliases': ['wp2_cl'], 'legend_label': 'wp2_cl'},
            {'aliases': ['wp2_pd'], 'legend_label': 'wp2_pd'},
            {'aliases': ['wp2_splat'], 'legend_label': 'wp2_splat'},
            {'aliases': ['wp2_sf'], 'legend_label': 'wp2_sf'},
            {'aliases': ['wp2_residual'], 'legend_label': 'wp2_residual', 'fallback_func': self.getWp2Residual},
        ]
        
        wp3_budget_lines = [
            {'aliases': ['wp3_bt'], 'legend_label': 'wp3_bt'},
            {'aliases': ['wp3_ma'], 'legend_label': 'wp3_ma'},
            {'aliases': ['wp3_ta'], 'legend_label': 'wp3_ta'},
            {'aliases': ['wp3_ac'], 'legend_label': 'wp3_ac'},
            {'aliases': ['wp3_pr1'], 'legend_label': 'wp3_pr1'},
            {'aliases': ['wp3_pr2'], 'legend_label': 'wp3_pr2'},
            {'aliases': ['wp3_pr3'], 'legend_label': 'wp3_pr3'},
            {'aliases': ['wp3_bp1'], 'legend_label': 'wp3_bp1'},
            {'aliases': ['wp3_bp2'], 'legend_label': 'wp3_bp2'},
            {'aliases': ['wp3_dp1'], 'legend_label': 'wp3_dp1'},
            {'aliases': ['wp3_tp'], 'legend_label': 'wp3_tp'},
            {'aliases': ['wp3_cl'], 'legend_label': 'wp3_cl'},
            {'aliases': ['wp3_splat'], 'legend_label': 'wp3_splat'},
            {'aliases': ['wp3_residual'], 'legend_label': 'wp3_residual', 'fallback_func': self.getWp3Residual},
        ]
        
        thlp2_budget_lines = [
            {'aliases': ['thlp2_bt'], 'legend_label': 'thlp2_bt'},
            {'aliases': ['thlp2_ma'], 'legend_label': 'thlp2_ma'},
            {'aliases': ['thlp2_ta'], 'legend_label': 'thlp2_ta'},
            {'aliases': ['thlp2_tp'], 'legend_label': 'thlp2_tp'},
            {'aliases': ['thlp2_dp1'], 'legend_label': 'thlp2_dp1'},
            {'aliases': ['thlp2_dp2'], 'legend_label': 'thlp2_dp2'},
            {'aliases': ['thlp2_cl'], 'legend_label': 'thlp2_cl'},
            {'aliases': ['thlp2_pd'], 'legend_label': 'thlp2_pd'},
            {'aliases': ['thlp2_sf'], 'legend_label': 'thlp2_sf'},
            {'aliases': ['thlp2_forcing'], 'legend_label': 'thlp2_forcing'},
            {'aliases': ['thlp2_residual'], 'legend_label': 'thlp2_residual', 'fallback_func': self.getThlp2Residual},
        ]
        
        rtp2_budget_lines = [
            {'aliases': ['rtp2_bt'], 'legend_label': 'rtp2_bt'},
            {'aliases': ['rtp2_ma'], 'legend_label': 'rtp2_ma'},
            {'aliases': ['rtp2_ta'], 'legend_label': 'rtp2_ta'},
            {'aliases': ['rtp2_tp'], 'legend_label': 'rtp2_tp'},
            {'aliases': ['rtp2_dp1'], 'legend_label': 'rtp2_dp1'},
            {'aliases': ['rtp2_dp2'], 'legend_label': 'rtp2_dp2'},
            {'aliases': ['rtp2_cl'], 'legend_label': 'rtp2_cl'},
            {'aliases': ['rtp2_pd'], 'legend_label': 'rtp2_pd'},
            {'aliases': ['rtp2_sf'], 'legend_label': 'rtp2_sf'},
            {'aliases': ['rtp2_forcing'], 'legend_label': 'rtp2_forcing'},
            {'aliases': ['rtp2_residual'], 'legend_label': 'rtp2_residual', 'fallback_func': self.getRtp2Residual},
        ]
        
        rtpthlp_budget_lines = [
            {'aliases': ['rtpthlp_bt'], 'legend_label': 'rtpthlp_bt'},
            {'aliases': ['rtpthlp_ma'], 'legend_label': 'rtpthlp_ma'},
            {'aliases': ['rtpthlp_ta'], 'legend_label': 'rtpthlp_ta'},
            {'aliases': ['rtpthlp_tp1'], 'legend_label': 'rtpthlp_tp1'},
            {'aliases': ['rtpthlp_dp1'], 'legend_label': 'rtpthlp_dp1'},
            {'aliases': ['rtpthlp_dp2'], 'legend_label': 'rtpthlp_dp2'},
            {'aliases': ['rtpthlp_cl'], 'legend_label': 'rtpthlp_cl'},
            {'aliases': ['rtpthlp_tp2'], 'legend_label': 'rtpthlp_tp2'},
            {'aliases': ['rtpthlp_sf'], 'legend_label': 'rtpthlp_sf'},
            {'aliases': ['rtpthlp_forcing'], 'legend_label': 'rtpthlp_forcing'},
            {'aliases': ['rtpthlp_residual'], 'legend_label': 'rtpthlp_residual', 'fallback_func': self.getRtpthlpResidual},
        ]
        
        upwp_budget_lines = [
            {'aliases': ['upwp_bt'], 'legend_label': 'upwp_bt'},
            {'aliases': ['upwp_ma'], 'legend_label': 'upwp_ma'},
            {'aliases': ['upwp_ta'], 'legend_label': 'upwp_ta'},
            {'aliases': ['upwp_tp'], 'legend_label': 'upwp_tp'},
            {'aliases': ['upwp_ac'], 'legend_label': 'upwp_ac'},
            {'aliases': ['upwp_bp'], 'legend_label': 'upwp_bp'},
            {'aliases': ['upwp_pr1'], 'legend_label': 'upwp_pr1'},
            {'aliases': ['upwp_pr2'], 'legend_label': 'upwp_pr2'},
            {'aliases': ['upwp_pr3'], 'legend_label': 'upwp_pr3'},
            {'aliases': ['upwp_pr4'], 'legend_label': 'upwp_pr4'},
            {'aliases': ['upwp_dp1'], 'legend_label': 'upwp_dp1'},
            {'aliases': ['upwp_cl'], 'legend_label': 'upwp_cl'},
            {'aliases': ['upwp_mfl'], 'legend_label': 'upwp_mfl'},
            {'aliases': ['upwp_residual'], 'legend_label': 'upwp_residual', 'fallback_func': self.getUpwpResidual},
        ]
        
        vpwp_budget_lines = [
            {'aliases': ['vpwp_bt'], 'legend_label': 'vpwp_bt'},
            {'aliases': ['vpwp_ma'], 'legend_label': 'vpwp_ma'},
            {'aliases': ['vpwp_ta'], 'legend_label': 'vpwp_ta'},
            {'aliases': ['vpwp_tp'], 'legend_label': 'vpwp_tp'},
            {'aliases': ['vpwp_ac'], 'legend_label': 'vpwp_ac'},
            {'aliases': ['vpwp_bp'], 'legend_label': 'vpwp_bp'},
            {'aliases': ['vpwp_pr1'], 'legend_label': 'vpwp_pr1'},
            {'aliases': ['vpwp_pr2'], 'legend_label': 'vpwp_pr2'},
            {'aliases': ['vpwp_pr3'], 'legend_label': 'vpwp_pr3'},
            {'aliases': ['vpwp_pr4'], 'legend_label': 'vpwp_pr4'},
            {'aliases': ['vpwp_dp1'], 'legend_label': 'vpwp_dp1'},
            {'aliases': ['vpwp_cl'], 'legend_label': 'vpwp_cl'},
            {'aliases': ['vpwp_mfl'], 'legend_label': 'vpwp_mfl'},
            {'aliases': ['vpwp_residual'], 'legend_label': 'vpwp_residual', 'fallback_func': self.getVpwpResidual},
        ]
        
        self.variable_definitions = [
            {'aliases': ['thlm'], 'lines': thlm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['rtm'], 'lines': rtm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['wpthlp'], 'lines': wpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['wprtp'], 'lines': wprtp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['wp2'], 'lines': wp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['wp3'], 'lines': wp3_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['thlp2'], 'lines': thlp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['rtp2'], 'lines': rtp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['rtpthlp'], 'lines': rtpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['upwp'], 'lines': upwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'aliases': ['vpwp'], 'lines': vpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},

        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_dataset = e3sm_dataset)

    def getThlmClipping(self, dataset_override = None):
        '''


        thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        :return:
        '''

        z = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset_override, fill_zeros=True)
        thlm_mfl = self.getVarForCalculations('thlm_mfl', dataset_override, fill_zeros=True)
        thlm_cl = self.getVarForCalculations('thlm_cl', dataset_override, fill_zeros=True)
        thlm_tacl = self.getVarForCalculations('thlm_tacl', dataset_override, fill_zeros=True)
        thlm_sdmp = self.getVarForCalculations('thlm_sdmp', dataset_override, fill_zeros=True)

        output_data = thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp

        return output_data, z

    def getLsforcing(self, dataset_override = None):
        '''


        thlm_forcing-radht-thlm_mc
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        thlm_forcing = self.getVarForCalculations('thlm_forcing', dataset_override, fill_zeros=True)
        radht = self.getVarForCalculations('radht', dataset_override, fill_zeros=True)
        thlm_mc = self.getVarForCalculations('thlm_mc', dataset_override, fill_zeros=True)

        output_data = thlm_forcing - radht - thlm_mc

        return output_data, z

    def getThlmResidual(self, dataset_override = None):
        '''


        thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        thlm_mfl = self.getVarForCalculations('thlm_mfl', dataset_override, fill_zeros=True)
        thlm_cl = self.getVarForCalculations('thlm_cl', dataset_override, fill_zeros=True)
        thlm_tacl = self.getVarForCalculations('thlm_tacl', dataset_override, fill_zeros=True)
        thlm_sdmp = self.getVarForCalculations('thlm_sdmp', dataset_override, fill_zeros=True)
        thlm_bt = self.getVarForCalculations('thlm_bt', dataset_override, fill_zeros=True)
        thlm_ta = self.getVarForCalculations('thlm_ta', dataset_override, fill_zeros=True)
        thlm_forcing = self.getVarForCalculations('thlm_forcing', dataset_override, fill_zeros=True)
        thlm_ma = self.getVarForCalculations('thlm_ma', dataset_override, fill_zeros=True)

        output_data = thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)

        return output_data, z

    def getRtmClipping(self, dataset_override = None):
        '''


        rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mfl = self.getVarForCalculations('rtm_mfl', dataset_override, fill_zeros=True)
        rtm_cl = self.getVarForCalculations('rtm_cl', dataset_override, fill_zeros=True)
        rtm_tacl = self.getVarForCalculations('rtm_tacl', dataset_override, fill_zeros=True)
        rtm_sdmp = self.getVarForCalculations('rtm_sdmp', dataset_override, fill_zeros=True)

        output_data = rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp

        return output_data, z

    def getRtmForcing(self, dataset_override = None):
        '''


        rtm_forcing - rtm_mc
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mc = self.getVarForCalculations('rtm_mc', dataset_override, fill_zeros=True)
        rtm_forcing = self.getVarForCalculations('rtm_forcing', dataset_override, fill_zeros=True)

        output_data = rtm_forcing - rtm_mc

        return output_data, z

    def getRtmResidual(self, dataset_override = None):
        '''


        rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mfl = self.getVarForCalculations('rtm_mfl', dataset_override, fill_zeros=True)
        rtm_cl = self.getVarForCalculations('rtm_cl', dataset_override, fill_zeros=True)
        rtm_tacl = self.getVarForCalculations('rtm_tacl', dataset_override, fill_zeros=True)
        rtm_sdmp = self.getVarForCalculations('rtm_sdmp', dataset_override, fill_zeros=True)
        rtm_bt = self.getVarForCalculations('rtm_bt', dataset_override, fill_zeros=True)
        rtm_ta = self.getVarForCalculations('rtm_ta', dataset_override, fill_zeros=True)
        rtm_forcing = self.getVarForCalculations('rtm_forcing', dataset_override, fill_zeros=True)
        rtm_pd = self.getVarForCalculations('rtm_pd', dataset_override, fill_zeros=True)
        rtm_ma = self.getVarForCalculations('rtm_ma', dataset_override, fill_zeros=True)

        output_data = rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)

        return output_data, z

    def getWpthlpResidual(self, dataset_override = None):
        '''


        wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wpthlp_mfl = self.getVarForCalculations('wpthlp_mfl', dataset_override, fill_zeros=True)
        wpthlp_cl = self.getVarForCalculations('wpthlp_cl', dataset_override, fill_zeros=True)
        wpthlp_tp = self.getVarForCalculations('wpthlp_tp', dataset_override, fill_zeros=True)
        wpthlp_ac = self.getVarForCalculations('wpthlp_ac', dataset_override, fill_zeros=True)
        wpthlp_pr1 = self.getVarForCalculations('wpthlp_pr1', dataset_override, fill_zeros=True)
        wpthlp_pr3 = self.getVarForCalculations('wpthlp_pr3', dataset_override, fill_zeros=True)
        wpthlp_pr2 = self.getVarForCalculations('wpthlp_pr2', dataset_override, fill_zeros=True)
        wpthlp_dp1 = self.getVarForCalculations('wpthlp_dp1', dataset_override, fill_zeros=True)
        wpthlp_sicl = self.getVarForCalculations('wpthlp_sicl', dataset_override, fill_zeros=True)
        wpthlp_bt = self.getVarForCalculations('wpthlp_bt', dataset_override, fill_zeros=True)
        wpthlp_ta = self.getVarForCalculations('wpthlp_ta', dataset_override, fill_zeros=True)
        wpthlp_forcing = self.getVarForCalculations('wpthlp_forcing', dataset_override, fill_zeros=True)
        wpthlp_bp = self.getVarForCalculations('wpthlp_bp', dataset_override, fill_zeros=True)
        wpthlp_ma = self.getVarForCalculations('wpthlp_ma', dataset_override, fill_zeros=True)

        output_data = wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)

        return output_data, z

    def getWprtpResidual(self, dataset_override = None):
        '''


        wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wprtp_mfl = self.getVarForCalculations('wprtp_mfl', dataset_override, fill_zeros=True)
        wprtp_cl = self.getVarForCalculations('wprtp_cl', dataset_override, fill_zeros=True)
        wprtp_tp = self.getVarForCalculations('wprtp_tp', dataset_override, fill_zeros=True)
        wprtp_ac = self.getVarForCalculations('wprtp_ac', dataset_override, fill_zeros=True)
        wprtp_pr1 = self.getVarForCalculations('wprtp_pr1', dataset_override, fill_zeros=True)
        wprtp_pr3 = self.getVarForCalculations('wprtp_pr3', dataset_override, fill_zeros=True)
        wprtp_pr2 = self.getVarForCalculations('wprtp_pr2', dataset_override, fill_zeros=True)
        wprtp_dp1 = self.getVarForCalculations('wprtp_dp1', dataset_override, fill_zeros=True)
        wprtp_sicl = self.getVarForCalculations('wprtp_sicl', dataset_override, fill_zeros=True)
        wprtp_bt = self.getVarForCalculations('wprtp_bt', dataset_override, fill_zeros=True)
        wprtp_ta = self.getVarForCalculations('wprtp_ta', dataset_override, fill_zeros=True)
        wprtp_forcing = self.getVarForCalculations('wprtp_forcing', dataset_override, fill_zeros=True)
        wprtp_bp = self.getVarForCalculations('wprtp_bp', dataset_override, fill_zeros=True)
        wprtp_ma = self.getVarForCalculations('wprtp_ma', dataset_override, fill_zeros=True)
        wprtp_pd = self.getVarForCalculations('wprtp_pd', dataset_override, fill_zeros=True)

        output_data = wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_pd + wprtp_forcing)

        return output_data, z

    def getWp2Residual(self, dataset_override = None):
        '''


        wp2_bt - (wp2_ma + wp2_ta + wp2_tp + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_mfl + wp2_cl + wp2_sicl + wp2_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wp2_sf = self.getVarForCalculations('wp2_sf', dataset_override, fill_zeros=True)
        wp2_cl = self.getVarForCalculations('wp2_cl', dataset_override, fill_zeros=True)
        wp2_ac = self.getVarForCalculations('wp2_ac', dataset_override, fill_zeros=True)
        wp2_pr1 = self.getVarForCalculations('wp2_pr1', dataset_override, fill_zeros=True)
        wp2_pr3 = self.getVarForCalculations('wp2_pr3', dataset_override, fill_zeros=True)
        wp2_pr2 = self.getVarForCalculations('wp2_pr2', dataset_override, fill_zeros=True)
        wp2_dp1 = self.getVarForCalculations('wp2_dp1', dataset_override, fill_zeros=True)
        wp2_dp2 = self.getVarForCalculations('wp2_dp2', dataset_override, fill_zeros=True)
        wp2_bt = self.getVarForCalculations('wp2_bt', dataset_override, fill_zeros=True)
        wp2_ta = self.getVarForCalculations('wp2_ta', dataset_override, fill_zeros=True)
        wp2_splat = self.getVarForCalculations('wp2_splat', dataset_override, fill_zeros=True)
        wp2_bp = self.getVarForCalculations('wp2_bp', dataset_override, fill_zeros=True)
        wp2_ma = self.getVarForCalculations('wp2_ma', dataset_override, fill_zeros=True)
        wp2_pd = self.getVarForCalculations('wp2_pd', dataset_override, fill_zeros=True)

        output_data = wp2_bt - (wp2_ma + wp2_ta + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_dp2 + wp2_cl + wp2_pd + wp2_sf + wp2_splat)

        return output_data, z

    def getWp3Residual(self, dataset_override = None):
        '''


        wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl+wp3_splat)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wp3_bp1 = self.getVarForCalculations('wp3_bp1', dataset_override, fill_zeros=True)
        wp3_bp2 = self.getVarForCalculations('wp3_bp2', dataset_override, fill_zeros=True)
        wp3_cl = self.getVarForCalculations('wp3_cl', dataset_override, fill_zeros=True)
        wp3_ac = self.getVarForCalculations('wp3_ac', dataset_override, fill_zeros=True)
        wp3_pr1 = self.getVarForCalculations('wp3_pr1', dataset_override, fill_zeros=True)
        wp3_pr3 = self.getVarForCalculations('wp3_pr3', dataset_override, fill_zeros=True)
        wp3_pr2 = self.getVarForCalculations('wp3_pr2', dataset_override, fill_zeros=True)
        wp3_dp1 = self.getVarForCalculations('wp3_dp1', dataset_override, fill_zeros=True)
        wp3_bt = self.getVarForCalculations('wp3_bt', dataset_override, fill_zeros=True)
        wp3_ta = self.getVarForCalculations('wp3_ta', dataset_override, fill_zeros=True)
        wp3_splat = self.getVarForCalculations('wp3_splat', dataset_override, fill_zeros=True)
        wp3_ma = self.getVarForCalculations('wp3_ma', dataset_override, fill_zeros=True)
        wp3_tp = self.getVarForCalculations('wp3_tp', dataset_override, fill_zeros=True)


        output_data = wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl+wp3_splat)

        return output_data, z

    def getThlp2Residual(self, dataset_override = None):
        '''


        thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        thlp2_cl = self.getVarForCalculations('thlp2_cl', dataset_override, fill_zeros=True)
        thlp2_dp2 = self.getVarForCalculations('thlp2_dp2', dataset_override, fill_zeros=True)
        thlp2_forcing = self.getVarForCalculations('thlp2_forcing', dataset_override, fill_zeros=True)
        thlp2_sf = self.getVarForCalculations('thlp2_sf', dataset_override, fill_zeros=True)
        thlp2_dp1 = self.getVarForCalculations('thlp2_dp1', dataset_override, fill_zeros=True)
        thlp2_bt = self.getVarForCalculations('thlp2_bt', dataset_override, fill_zeros=True)
        thlp2_ta = self.getVarForCalculations('thlp2_ta', dataset_override, fill_zeros=True)
        thlp2_pd = self.getVarForCalculations('thlp2_pd', dataset_override, fill_zeros=True)
        thlp2_ma = self.getVarForCalculations('thlp2_ma', dataset_override, fill_zeros=True)
        thlp2_tp = self.getVarForCalculations('thlp2_tp', dataset_override, fill_zeros=True)


        output_data = thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)

        return output_data, z

    def getRtp2Residual(self, dataset_override = None):
        '''


        rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtp2_cl = self.getVarForCalculations('rtp2_cl', dataset_override, fill_zeros=True)
        rtp2_dp2 = self.getVarForCalculations('rtp2_dp2', dataset_override, fill_zeros=True)
        rtp2_forcing = self.getVarForCalculations('rtp2_forcing', dataset_override, fill_zeros=True)
        rtp2_sf = self.getVarForCalculations('rtp2_sf', dataset_override, fill_zeros=True)
        rtp2_dp1 = self.getVarForCalculations('rtp2_dp1', dataset_override, fill_zeros=True)
        rtp2_bt = self.getVarForCalculations('rtp2_bt', dataset_override, fill_zeros=True)
        rtp2_ta = self.getVarForCalculations('rtp2_ta', dataset_override, fill_zeros=True)
        rtp2_pd = self.getVarForCalculations('rtp2_pd', dataset_override, fill_zeros=True)
        rtp2_ma = self.getVarForCalculations('rtp2_ma', dataset_override, fill_zeros=True)
        rtp2_tp = self.getVarForCalculations('rtp2_tp', dataset_override, fill_zeros=True)


        output_data = rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)

        return output_data, z

    def getRtpthlpResidual(self, dataset_override = None):
        '''


        rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_pd + rtpthlp_sf + rtpthlp_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtpthlp_cl = self.getVarForCalculations('rtpthlp_cl', dataset_override, fill_zeros=True)
        rtpthlp_dp2 = self.getVarForCalculations('rtpthlp_dp2', dataset_override, fill_zeros=True)
        rtpthlp_forcing = self.getVarForCalculations('rtpthlp_forcing', dataset_override, fill_zeros=True)
        rtpthlp_sf = self.getVarForCalculations('rtpthlp_sf', dataset_override, fill_zeros=True)
        rtpthlp_dp1 = self.getVarForCalculations('rtpthlp_dp1', dataset_override, fill_zeros=True)
        rtpthlp_bt = self.getVarForCalculations('rtpthlp_bt', dataset_override, fill_zeros=True)
        rtpthlp_ta = self.getVarForCalculations('rtpthlp_ta', dataset_override, fill_zeros=True)
        rtpthlp_tp2 = self.getVarForCalculations('rtpthlp_tp2', dataset_override, fill_zeros=True)
        rtpthlp_ma = self.getVarForCalculations('rtpthlp_ma', dataset_override, fill_zeros=True)
        rtpthlp_tp1 = self.getVarForCalculations('rtpthlp_tp1', dataset_override, fill_zeros=True)


        output_data = rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp1 + rtpthlp_tp2 + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_sf + rtpthlp_forcing)

        return output_data, z

    def getUpwpResidual(self, dataset_override = None):
        '''


        upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_dp1 + upwp_dp2 + upwp_cl + upwp_pd + upwp_sf + upwp_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        upwp_cl = self.getVarForCalculations('upwp_cl', dataset_override, fill_zeros=True)
        upwp_tp = self.getVarForCalculations('upwp_tp', dataset_override, fill_zeros=True)
        upwp_ac = self.getVarForCalculations('upwp_ac', dataset_override, fill_zeros=True)
        upwp_bp = self.getVarForCalculations('upwp_bp', dataset_override, fill_zeros=True)
        upwp_dp1 = self.getVarForCalculations('upwp_dp1', dataset_override, fill_zeros=True)
        upwp_bt = self.getVarForCalculations('upwp_bt', dataset_override, fill_zeros=True)
        upwp_ta = self.getVarForCalculations('upwp_ta', dataset_override, fill_zeros=True)
        upwp_pr1 = self.getVarForCalculations('upwp_pr1', dataset_override, fill_zeros=True)
        upwp_pr2 = self.getVarForCalculations('upwp_pr2', dataset_override, fill_zeros=True)
        upwp_pr3 = self.getVarForCalculations('upwp_pr3', dataset_override, fill_zeros=True)
        upwp_pr4 = self.getVarForCalculations('upwp_pr4', dataset_override, fill_zeros=True)
        upwp_mfl = self.getVarForCalculations('upwp_mfl', dataset_override, fill_zeros=True)
        upwp_ma = self.getVarForCalculations('upwp_ma', dataset_override, fill_zeros=True)


        output_data = upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_ac + upwp_bp + upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4 + upwp_dp1 + upwp_mfl + upwp_cl)

        return output_data, z

    def getVpwpResidual(self, dataset_override = None):
        '''


        vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_dp1 + vpwp_dp2 + vpwp_cl + vpwp_pd + vpwp_sf + vpwp_forcing)
        :return:
        '''
        z = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        vpwp_cl = self.getVarForCalculations('vpwp_cl', dataset_override, fill_zeros=True)
        vpwp_tp = self.getVarForCalculations('vpwp_tp', dataset_override, fill_zeros=True)
        vpwp_ac = self.getVarForCalculations('vpwp_ac', dataset_override, fill_zeros=True)
        vpwp_bp = self.getVarForCalculations('vpwp_bp', dataset_override, fill_zeros=True)
        vpwp_dp1 = self.getVarForCalculations('vpwp_dp1', dataset_override, fill_zeros=True)
        vpwp_bt = self.getVarForCalculations('vpwp_bt', dataset_override, fill_zeros=True)
        vpwp_ta = self.getVarForCalculations('vpwp_ta', dataset_override, fill_zeros=True)
        vpwp_pr1 = self.getVarForCalculations('vpwp_pr1', dataset_override, fill_zeros=True)
        vpwp_pr2 = self.getVarForCalculations('vpwp_pr2', dataset_override, fill_zeros=True)
        vpwp_pr3 = self.getVarForCalculations('vpwp_pr3', dataset_override, fill_zeros=True)
        vpwp_pr4 = self.getVarForCalculations('vpwp_pr4', dataset_override, fill_zeros=True)
        vpwp_mfl = self.getVarForCalculations('vpwp_mfl', dataset_override, fill_zeros=True)
        vpwp_ma = self.getVarForCalculations('vpwp_ma', dataset_override, fill_zeros=True)


        output_data = vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_ac + vpwp_bp + vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4 + vpwp_dp1 + vpwp_mfl + vpwp_cl)

        return output_data, z