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
                 e3sm_datasets=None):
        self.name = "base variables budgets"

        thlm_budget_lines = [
            {'var_names': ['thlm_bt'], 'legend_label': 'thlm_bt'},
            {'var_names': ['thlm_ma'], 'legend_label': 'thlm_ma'},
            {'var_names': ['thlm_ta'], 'legend_label': 'thlm_ta'},
            {'var_names': ['thlm_mc'], 'legend_label': 'thlm_mc'},
            {'var_names': ['thlm_clipping'], 'legend_label': 'thlm_clipping', 'fallback_func': self.getThlmClipping},
            {'var_names': ['radht'], 'legend_label': 'radht'},
            {'var_names': ['ls_forcing'], 'legend_label': 'thlm_ls_forcing', 'fallback_func': self.getLsforcing},
            {'var_names': ['thlm_residual'], 'legend_label': 'thlm_residual', 'fallback_func': self.getThlmResidual},

        ]

        rtm_budget_lines = [
            {'var_names': ['rtm_bt'], 'legend_label': 'rtm_bt'},
            {'var_names': ['rtm_ma'], 'legend_label': 'rtm_ma'},
            {'var_names': ['rtm_ta'], 'legend_label': 'rtm_ta'},
            {'var_names': ['rtm_mc'], 'legend_label': 'rtm_mc'},
            {'var_names': ['rtm_clipping'], 'legend_label': 'rtm_clipping', 'fallback_func': self.getRtmClipping},
            {'var_names': ['rtm_pd'], 'legend_label': 'rtm_pd'},
            {'var_names': ['ls_forcing'], 'legend_label': 'ls_forcing', 'fallback_func': self.getRtmForcing},
            {'var_names': ['rtm_residual'], 'legend_label': 'rtm_residual', 'fallback_func': self.getRtmResidual},

        ]

        wpthlp_budget_lines = [
            {'var_names': ['wpthlp_bt'], 'legend_label': 'wpthlp_bt'},
            {'var_names': ['wpthlp_ma'], 'legend_label': 'wpthlp_ma'},
            {'var_names': ['wpthlp_ta'], 'legend_label': 'wpthlp_ta'},
            {'var_names': ['wpthlp_tp'], 'legend_label': 'wpthlp_tp'},
            {'var_names': ['wpthlp_ac'], 'legend_label': 'wpthlp_ac'},
            {'var_names': ['wpthlp_bp'], 'legend_label': 'wpthlp_bp'},
            {'var_names': ['wpthlp_pr1'], 'legend_label': 'wpthlp_pr1'},
            {'var_names': ['wpthlp_pr2'], 'legend_label': 'wpthlp_pr2'},
            {'var_names': ['wpthlp_pr3'], 'legend_label': 'wpthlp_pr3'},
            {'var_names': ['wpthlp_dp1'], 'legend_label': 'wpthlp_dp1'},
            {'var_names': ['wpthlp_mfl'], 'legend_label': 'wpthlp_mfl'},
            {'var_names': ['wpthlp_cl'], 'legend_label': 'wpthlp_cl'},
            {'var_names': ['wpthlp_sicl'], 'legend_label': 'wpthlp_sicl'},
            {'var_names': ['wpthlp_forcing'], 'legend_label': 'wpthlp_forcing'},
            {'var_names': ['wpthlp_residual'], 'legend_label': 'wpthlp_residual',
             'fallback_func': self.getWpthlpResidual},

        ]

        wprtp_budget_lines = [
            {'var_names': ['wprtp_bt'], 'legend_label': 'wprtp_bt'},
            {'var_names': ['wprtp_ma'], 'legend_label': 'wprtp_ma'},
            {'var_names': ['wprtp_ta'], 'legend_label': 'wprtp_ta'},
            {'var_names': ['wprtp_tp'], 'legend_label': 'wprtp_tp'},
            {'var_names': ['wprtp_ac'], 'legend_label': 'wprtp_ac'},
            {'var_names': ['wprtp_bp'], 'legend_label': 'wprtp_bp'},
            {'var_names': ['wprtp_pr1'], 'legend_label': 'wprtp_pr1'},
            {'var_names': ['wprtp_pr2'], 'legend_label': 'wprtp_pr2'},
            {'var_names': ['wprtp_pr3'], 'legend_label': 'wprtp_pr3'},
            {'var_names': ['wprtp_dp1'], 'legend_label': 'wprtp_dp1'},
            {'var_names': ['wprtp_mfl'], 'legend_label': 'wprtp_mfl'},
            {'var_names': ['wprtp_cl'], 'legend_label': 'wprtp_cl'},
            {'var_names': ['wprtp_sicl'], 'legend_label': 'wprtp_sicl'},
            {'var_names': ['wprtp_pd'], 'legend_label': 'wprtp_pd'},
            {'var_names': ['wprtp_forcing'], 'legend_label': 'wprtp_forcing'},
            {'var_names': ['wprtp_residual'], 'legend_label': 'wprtp_residual', 'fallback_func': self.getWprtpResidual},
        ]

        wp2_budget_lines = [
            {'var_names': ['wp2_bt'], 'legend_label': 'wp2_bt'},
            {'var_names': ['wp2_ma'], 'legend_label': 'wp2_ma'},
            {'var_names': ['wp2_ta'], 'legend_label': 'wp2_ta'},
            {'var_names': ['wp2_ac'], 'legend_label': 'wp2_ac'},
            {'var_names': ['wp2_bp'], 'legend_label': 'wp2_bp'},
            {'var_names': ['wp2_pr1'], 'legend_label': 'wp2_pr1'},
            {'var_names': ['wp2_pr2'], 'legend_label': 'wp2_pr2'},
            {'var_names': ['wp2_pr3'], 'legend_label': 'wp2_pr3'},
            {'var_names': ['wp2_dp1'], 'legend_label': 'wp2_dp1'},
            {'var_names': ['wp2_dp2'], 'legend_label': 'wp2_dp2'},
            {'var_names': ['wp2_cl'], 'legend_label': 'wp2_cl'},
            {'var_names': ['wp2_pd'], 'legend_label': 'wp2_pd'},
            {'var_names': ['wp2_splat'], 'legend_label': 'wp2_splat'},
            {'var_names': ['wp2_sf'], 'legend_label': 'wp2_sf'},
            {'var_names': ['wp2_residual'], 'legend_label': 'wp2_residual', 'fallback_func': self.getWp2Residual},
        ]

        wp3_budget_lines = [
            {'var_names': ['wp3_bt'], 'legend_label': 'wp3_bt'},
            {'var_names': ['wp3_ma'], 'legend_label': 'wp3_ma'},
            {'var_names': ['wp3_ta'], 'legend_label': 'wp3_ta'},
            {'var_names': ['wp3_ac'], 'legend_label': 'wp3_ac'},
            {'var_names': ['wp3_pr1'], 'legend_label': 'wp3_pr1'},
            {'var_names': ['wp3_pr2'], 'legend_label': 'wp3_pr2'},
            {'var_names': ['wp3_pr3'], 'legend_label': 'wp3_pr3'},
            {'var_names': ['wp3_bp1'], 'legend_label': 'wp3_bp1'},
            {'var_names': ['wp3_bp2'], 'legend_label': 'wp3_bp2'},
            {'var_names': ['wp3_dp1'], 'legend_label': 'wp3_dp1'},
            {'var_names': ['wp3_tp'], 'legend_label': 'wp3_tp'},
            {'var_names': ['wp3_cl'], 'legend_label': 'wp3_cl'},
            {'var_names': ['wp3_splat'], 'legend_label': 'wp3_splat'},
            {'var_names': ['wp3_residual'], 'legend_label': 'wp3_residual', 'fallback_func': self.getWp3Residual},
        ]

        thlp2_budget_lines = [
            {'var_names': ['thlp2_bt'], 'legend_label': 'thlp2_bt'},
            {'var_names': ['thlp2_ma'], 'legend_label': 'thlp2_ma'},
            {'var_names': ['thlp2_ta'], 'legend_label': 'thlp2_ta'},
            {'var_names': ['thlp2_tp'], 'legend_label': 'thlp2_tp'},
            {'var_names': ['thlp2_dp1'], 'legend_label': 'thlp2_dp1'},
            {'var_names': ['thlp2_dp2'], 'legend_label': 'thlp2_dp2'},
            {'var_names': ['thlp2_cl'], 'legend_label': 'thlp2_cl'},
            {'var_names': ['thlp2_pd'], 'legend_label': 'thlp2_pd'},
            {'var_names': ['thlp2_sf'], 'legend_label': 'thlp2_sf'},
            {'var_names': ['thlp2_forcing'], 'legend_label': 'thlp2_forcing'},
            {'var_names': ['thlp2_residual'], 'legend_label': 'thlp2_residual', 'fallback_func': self.getThlp2Residual},
        ]

        rtp2_budget_lines = [
            {'var_names': ['rtp2_bt'], 'legend_label': 'rtp2_bt'},
            {'var_names': ['rtp2_ma'], 'legend_label': 'rtp2_ma'},
            {'var_names': ['rtp2_ta'], 'legend_label': 'rtp2_ta'},
            {'var_names': ['rtp2_tp'], 'legend_label': 'rtp2_tp'},
            {'var_names': ['rtp2_dp1'], 'legend_label': 'rtp2_dp1'},
            {'var_names': ['rtp2_dp2'], 'legend_label': 'rtp2_dp2'},
            {'var_names': ['rtp2_cl'], 'legend_label': 'rtp2_cl'},
            {'var_names': ['rtp2_pd'], 'legend_label': 'rtp2_pd'},
            {'var_names': ['rtp2_sf'], 'legend_label': 'rtp2_sf'},
            {'var_names': ['rtp2_forcing'], 'legend_label': 'rtp2_forcing'},
            {'var_names': ['rtp2_residual'], 'legend_label': 'rtp2_residual', 'fallback_func': self.getRtp2Residual},
        ]

        rtpthlp_budget_lines = [
            {'var_names': ['rtpthlp_bt'], 'legend_label': 'rtpthlp_bt'},
            {'var_names': ['rtpthlp_ma'], 'legend_label': 'rtpthlp_ma'},
            {'var_names': ['rtpthlp_ta'], 'legend_label': 'rtpthlp_ta'},
            {'var_names': ['rtpthlp_tp1'], 'legend_label': 'rtpthlp_tp1'},
            {'var_names': ['rtpthlp_dp1'], 'legend_label': 'rtpthlp_dp1'},
            {'var_names': ['rtpthlp_dp2'], 'legend_label': 'rtpthlp_dp2'},
            {'var_names': ['rtpthlp_cl'], 'legend_label': 'rtpthlp_cl'},
            {'var_names': ['rtpthlp_tp2'], 'legend_label': 'rtpthlp_tp2'},
            {'var_names': ['rtpthlp_sf'], 'legend_label': 'rtpthlp_sf'},
            {'var_names': ['rtpthlp_forcing'], 'legend_label': 'rtpthlp_forcing'},
            {'var_names': ['rtpthlp_residual'], 'legend_label': 'rtpthlp_residual',
             'fallback_func': self.getRtpthlpResidual},
        ]

        upwp_budget_lines = [
            {'var_names': ['upwp_bt'], 'legend_label': 'upwp_bt'},
            {'var_names': ['upwp_ma'], 'legend_label': 'upwp_ma'},
            {'var_names': ['upwp_ta'], 'legend_label': 'upwp_ta'},
            {'var_names': ['upwp_tp'], 'legend_label': 'upwp_tp'},
            {'var_names': ['upwp_ac'], 'legend_label': 'upwp_ac'},
            {'var_names': ['upwp_bp'], 'legend_label': 'upwp_bp'},
            {'var_names': ['upwp_pr1'], 'legend_label': 'upwp_pr1'},
            {'var_names': ['upwp_pr2'], 'legend_label': 'upwp_pr2'},
            {'var_names': ['upwp_pr3'], 'legend_label': 'upwp_pr3'},
            {'var_names': ['upwp_pr4'], 'legend_label': 'upwp_pr4'},
            {'var_names': ['upwp_dp1'], 'legend_label': 'upwp_dp1'},
            {'var_names': ['upwp_cl'], 'legend_label': 'upwp_cl'},
            {'var_names': ['upwp_mfl'], 'legend_label': 'upwp_mfl'},
            {'var_names': ['upwp_residual'], 'legend_label': 'upwp_residual', 'fallback_func': self.getUpwpResidual},
        ]

        vpwp_budget_lines = [
            {'var_names': ['vpwp_bt'], 'legend_label': 'vpwp_bt'},
            {'var_names': ['vpwp_ma'], 'legend_label': 'vpwp_ma'},
            {'var_names': ['vpwp_ta'], 'legend_label': 'vpwp_ta'},
            {'var_names': ['vpwp_tp'], 'legend_label': 'vpwp_tp'},
            {'var_names': ['vpwp_ac'], 'legend_label': 'vpwp_ac'},
            {'var_names': ['vpwp_bp'], 'legend_label': 'vpwp_bp'},
            {'var_names': ['vpwp_pr1'], 'legend_label': 'vpwp_pr1'},
            {'var_names': ['vpwp_pr2'], 'legend_label': 'vpwp_pr2'},
            {'var_names': ['vpwp_pr3'], 'legend_label': 'vpwp_pr3'},
            {'var_names': ['vpwp_pr4'], 'legend_label': 'vpwp_pr4'},
            {'var_names': ['vpwp_dp1'], 'legend_label': 'vpwp_dp1'},
            {'var_names': ['vpwp_cl'], 'legend_label': 'vpwp_cl'},
            {'var_names': ['vpwp_mfl'], 'legend_label': 'vpwp_mfl'},
            {'var_names': ['vpwp_residual'], 'legend_label': 'vpwp_residual', 'fallback_func': self.getVpwpResidual},
        ]

        um_budget_lines = [
            {'var_names': ['um_bt'], 'legend_label': 'um_bt'},
            {'var_names': ['um_ma'], 'legend_label': 'um_ma'},
            {'var_names': ['um_gf'], 'legend_label': 'um_gf'},
            {'var_names': ['um_cf'], 'legend_label': 'um_cf'},
            {'var_names': ['um_ta'], 'legend_label': 'um_ta'},
            {'var_names': ['um_f'], 'legend_label': 'um_f'},
            {'var_names': ['um_sdmp'], 'legend_label': 'um_sdmp'},
            {'var_names': ['um_ndg'], 'legend_label': 'um_ndg'},
            {'var_names': ['um_mfl'], 'legend_label': 'um_mfl'}
        ]

        vm_budget_lines = [
            {'var_names': ['vm_bt'], 'legend_label': 'vm_bt'},
            {'var_names': ['vm_ma'], 'legend_label': 'vm_ma'},
            {'var_names': ['vm_gf'], 'legend_label': 'vm_gf'},
            {'var_names': ['vm_cf'], 'legend_label': 'vm_cf'},
            {'var_names': ['vm_ta'], 'legend_label': 'vm_ta'},
            {'var_names': ['vm_f'], 'legend_label': 'vm_f'},
            {'var_names': ['vm_sdmp'], 'legend_label': 'vm_sdmp'},
            {'var_names': ['vm_ndg'], 'legend_label': 'vm_ndg'},
            {'var_names': ['vm_mfl'], 'legend_label': 'vm_mfl'}
        ]

        rrm_budget_lines = [
            {'var_names': ['rrm_bt'], 'legend_label': 'rrm_bt'},
            {'var_names': ['rrm_ma'], 'legend_label': 'rrm_ma'},
            {'var_names': ['rrm_sd'], 'legend_label': 'rrm_sd'},
            {'var_names': ['rrm_ta'], 'legend_label': 'rrm_ta'},
            {'var_names': ['rrm_ts'], 'legend_label': 'rrm_ts'},
            {'var_names': ['rrm_fill_clip'], 'legend_label': 'rrm_fill_clip', 'fallback_func': self.getRrmFillClip},
            {'var_names': ['rrm_mc'], 'legend_label': 'rrm_mc'},
            {'var_names': ['rrm_residual'], 'legend_label': 'rrm_residual', 'fallback_func': self.getRrmResidual},
        ]

        Nrm_budget_lines = [
            {'var_names': ['Nrm_bt'], 'legend_label': 'Nrm_bt'},
            {'var_names': ['Nrm_ma'], 'legend_label': 'Nrm_ma'},
            {'var_names': ['Nrm_sd'], 'legend_label': 'Nrm_sd'},
            {'var_names': ['Nrm_ta'], 'legend_label': 'Nrm_ta'},
            {'var_names': ['Nrm_ts'], 'legend_label': 'Nrm_ts'},
            {'var_names': ['Nrm_cl'], 'legend_label': 'Nrm_cl'},
            {'var_names': ['Nrm_mc'], 'legend_label': 'Nrm_mc'},
            {'var_names': ['Nrm_residual'], 'legend_label': 'Nrm_residual', 'fallback_func': self.getNrmResidual},
        ]

        self.variable_definitions = [
            {'var_names': {
                'clubb': ['thlm'],
                'sam': ['thlm'],
                'coamps': ['thlm'],
                'r408': ['thlm'],
                'hoc': ['thlm'],
                'e3sm': ['thlm']
            },
                'lines': thlm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['rtm'],
                'sam': ['rtm'],
                'coamps': ['rtm'],
                'r408': ['rtm'],
                'hoc': ['rtm'],
                'e3sm': ['rtm']
            },
                'lines': rtm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['wpthlp'],
                'sam': ['wpthlp'],
                'coamps': ['wpthlp'],
                'r408': ['wpthlp'],
                'hoc': ['wpthlp'],
                'e3sm': ['wpthlp']
            },
                'lines': wpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['wprtp'],
                'sam': ['wprtp'],
                'coamps': ['wprtp'],
                'r408': ['wprtp'],
                'hoc': ['wprtp'],
                'e3sm': ['wprtp']
            },
                'lines': wprtp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['wp2'],
                'sam': ['wp2'],
                'coamps': ['wp2'],
                'r408': ['wp2'],
                'hoc': ['wp2'],
                'e3sm': ['wp2']
            },
                'lines': wp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['wp3'],
                'sam': ['wp3'],
                'coamps': ['wp3'],
                'r408': ['wp3'],
                'hoc': ['wp3'],
                'e3sm': ['wp3']
            },
                'lines': wp3_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['thlp2'],
                'sam': ['thlp2'],
                'coamps': ['thlp2'],
                'r408': ['thlp2'],
                'hoc': ['thlp2'],
                'e3sm': ['thlp2']
            },
                'lines': thlp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['rtp2'],
                'sam': ['rtp2'],
                'coamps': ['rtp2'],
                'r408': ['rtp2'],
                'hoc': ['rtp2'],
                'e3sm': ['rtp2']
            },
                'lines': rtp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['rtpthlp'],
                'sam': ['rtpthlp'],
                'coamps': ['rtpthlp'],
                'r408': ['rtpthlp'],
                'hoc': ['rtpthlp'],
                'e3sm': ['rtpthlp']
            },
                'lines': rtpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['upwp'],
                'sam': ['upwp'],
                'coamps': ['upwp'],
                'r408': ['upwp'],
                'hoc': ['upwp'],
                'e3sm': ['upwp']
            },
                'lines': upwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['vpwp'],
                'sam': ['vpwp'],
                'coamps': ['vpwp'],
                'r408': ['vpwp'],
                'hoc': ['vpwp'],
                'e3sm': ['vpwp']
            },
                'lines': vpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['um'],
                'sam': ['um'],
                'coamps': ['um'],
                'r408': ['um'],
                'hoc': ['um'],
                'e3sm': ['um']
            },
                'lines': um_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['vm'],
                'sam': ['vm'],
                'coamps': ['vm'],
                'r408': ['vm'],
                'hoc': ['vm'],
                'e3sm': ['vm']
            },
                'lines': vm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},

            {'var_names': {
                'clubb': ['rrm'],
                'sam': ['rrm'],
                'coamps': ['rrm'],
                'r408': ['rrm'],
                'hoc': ['rrm'],
                'e3sm': ['rrm']
            },
                'lines': rrm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},

            {'var_names': {
                'clubb': ['Nrm'],
                'sam': ['Nrm'],
                'coamps': ['Nrm'],
                'r408': ['Nrm'],
                'hoc': ['Nrm'],
                'e3sm': ['Nrm']
            },
                'lines': Nrm_budget_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets)

    def getThlmClipping(self, dataset_override=None):
        '''


        thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        :return:
        '''

        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset_override, fill_zeros=True)
        thlm_mfl, z, dataset = self.getVarForCalculations('thlm_mfl', dataset_override, fill_zeros=True)
        thlm_cl, z, dataset = self.getVarForCalculations('thlm_cl', dataset_override, fill_zeros=True)
        thlm_tacl, z, dataset = self.getVarForCalculations('thlm_tacl', dataset_override, fill_zeros=True)
        thlm_sdmp, z, dataset = self.getVarForCalculations('thlm_sdmp', dataset_override, fill_zeros=True)

        output_data = thlm_mfl + thlm_cl + thlm_tacl + thlm_sdmp

        return output_data, z

    def getLsforcing(self, dataset_override=None):
        '''


        thlm_forcing-radht-thlm_mc
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset_override, fill_zeros=True)
        thlm_forcing, z, dataset = self.getVarForCalculations('thlm_forcing', dataset_override, fill_zeros=True)
        radht, z, dataset = self.getVarForCalculations('radht', dataset_override, fill_zeros=True)
        thlm_mc, z, dataset = self.getVarForCalculations('thlm_mc', dataset_override, fill_zeros=True)

        output_data = thlm_forcing - radht - thlm_mc

        return output_data, z

    def getThlmResidual(self, dataset_override=None):
        '''


        thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        thlm_mfl, z, dataset = self.getVarForCalculations('thlm_mfl', dataset_override, fill_zeros=True)
        thlm_cl, z, dataset = self.getVarForCalculations('thlm_cl', dataset_override, fill_zeros=True)
        thlm_tacl, z, dataset = self.getVarForCalculations('thlm_tacl', dataset_override, fill_zeros=True)
        thlm_sdmp, z, dataset = self.getVarForCalculations('thlm_sdmp', dataset_override, fill_zeros=True)
        thlm_bt, z, dataset = self.getVarForCalculations('thlm_bt', dataset_override, fill_zeros=True)
        thlm_ta, z, dataset = self.getVarForCalculations('thlm_ta', dataset_override, fill_zeros=True)
        thlm_forcing, z, dataset = self.getVarForCalculations('thlm_forcing', dataset_override, fill_zeros=True)
        thlm_ma, z, dataset = self.getVarForCalculations('thlm_ma', dataset_override, fill_zeros=True)

        output_data = thlm_bt - (thlm_ma + thlm_ta + thlm_mfl + thlm_cl + thlm_tacl + thlm_sdmp + thlm_forcing)

        return output_data, z

    def getRtmClipping(self, dataset_override=None):
        '''


        rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mfl, z, dataset = self.getVarForCalculations('rtm_mfl', dataset_override, fill_zeros=True)
        rtm_cl, z, dataset = self.getVarForCalculations('rtm_cl', dataset_override, fill_zeros=True)
        rtm_tacl, z, dataset = self.getVarForCalculations('rtm_tacl', dataset_override, fill_zeros=True)
        rtm_sdmp, z, dataset = self.getVarForCalculations('rtm_sdmp', dataset_override, fill_zeros=True)

        output_data = rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp

        return output_data, z

    def getRtmForcing(self, dataset_override=None):
        '''


        rtm_forcing - rtm_mc
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mc, z, dataset = self.getVarForCalculations('rtm_mc', dataset_override, fill_zeros=True)
        rtm_forcing, z, dataset = self.getVarForCalculations('rtm_forcing', dataset_override, fill_zeros=True)

        output_data = rtm_forcing - rtm_mc

        return output_data, z

    def getRtmResidual(self, dataset_override=None):
        '''


        rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtm_mfl, z, dataset = self.getVarForCalculations('rtm_mfl', dataset_override, fill_zeros=True)
        rtm_cl, z, dataset = self.getVarForCalculations('rtm_cl', dataset_override, fill_zeros=True)
        rtm_tacl, z, dataset = self.getVarForCalculations('rtm_tacl', dataset_override, fill_zeros=True)
        rtm_sdmp, z, dataset = self.getVarForCalculations('rtm_sdmp', dataset_override, fill_zeros=True)
        rtm_bt, z, dataset = self.getVarForCalculations('rtm_bt', dataset_override, fill_zeros=True)
        rtm_ta, z, dataset = self.getVarForCalculations('rtm_ta', dataset_override, fill_zeros=True)
        rtm_forcing, z, dataset = self.getVarForCalculations('rtm_forcing', dataset_override, fill_zeros=True)
        rtm_pd, z, dataset = self.getVarForCalculations('rtm_pd', dataset_override, fill_zeros=True)
        rtm_ma, z, dataset = self.getVarForCalculations('rtm_ma', dataset_override, fill_zeros=True)

        output_data = rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)

        return output_data, z

    def getWpthlpResidual(self, dataset_override=None):
        '''


        wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wpthlp_mfl, z, dataset = self.getVarForCalculations('wpthlp_mfl', dataset_override, fill_zeros=True)
        wpthlp_cl, z, dataset = self.getVarForCalculations('wpthlp_cl', dataset_override, fill_zeros=True)
        wpthlp_tp, z, dataset = self.getVarForCalculations('wpthlp_tp', dataset_override, fill_zeros=True)
        wpthlp_ac, z, dataset = self.getVarForCalculations('wpthlp_ac', dataset_override, fill_zeros=True)
        wpthlp_pr1, z, dataset = self.getVarForCalculations('wpthlp_pr1', dataset_override, fill_zeros=True)
        wpthlp_pr3, z, dataset = self.getVarForCalculations('wpthlp_pr3', dataset_override, fill_zeros=True)
        wpthlp_pr2, z, dataset = self.getVarForCalculations('wpthlp_pr2', dataset_override, fill_zeros=True)
        wpthlp_dp1, z, dataset = self.getVarForCalculations('wpthlp_dp1', dataset_override, fill_zeros=True)
        wpthlp_sicl, z, dataset = self.getVarForCalculations('wpthlp_sicl', dataset_override, fill_zeros=True)
        wpthlp_bt, z, dataset = self.getVarForCalculations('wpthlp_bt', dataset_override, fill_zeros=True)
        wpthlp_ta, z, dataset = self.getVarForCalculations('wpthlp_ta', dataset_override, fill_zeros=True)
        wpthlp_forcing, z, dataset = self.getVarForCalculations('wpthlp_forcing', dataset_override, fill_zeros=True)
        wpthlp_bp, z, dataset = self.getVarForCalculations('wpthlp_bp', dataset_override, fill_zeros=True)
        wpthlp_ma, z, dataset = self.getVarForCalculations('wpthlp_ma', dataset_override, fill_zeros=True)

        output_data = wpthlp_bt - (
                    wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)

        return output_data, z

    def getWprtpResidual(self, dataset_override=None):
        '''


        wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wprtp_mfl, z, dataset = self.getVarForCalculations('wprtp_mfl', dataset_override, fill_zeros=True)
        wprtp_cl, z, dataset = self.getVarForCalculations('wprtp_cl', dataset_override, fill_zeros=True)
        wprtp_tp, z, dataset = self.getVarForCalculations('wprtp_tp', dataset_override, fill_zeros=True)
        wprtp_ac, z, dataset = self.getVarForCalculations('wprtp_ac', dataset_override, fill_zeros=True)
        wprtp_pr1, z, dataset = self.getVarForCalculations('wprtp_pr1', dataset_override, fill_zeros=True)
        wprtp_pr3, z, dataset = self.getVarForCalculations('wprtp_pr3', dataset_override, fill_zeros=True)
        wprtp_pr2, z, dataset = self.getVarForCalculations('wprtp_pr2', dataset_override, fill_zeros=True)
        wprtp_dp1, z, dataset = self.getVarForCalculations('wprtp_dp1', dataset_override, fill_zeros=True)
        wprtp_sicl, z, dataset = self.getVarForCalculations('wprtp_sicl', dataset_override, fill_zeros=True)
        wprtp_bt, z, dataset = self.getVarForCalculations('wprtp_bt', dataset_override, fill_zeros=True)
        wprtp_ta, z, dataset = self.getVarForCalculations('wprtp_ta', dataset_override, fill_zeros=True)
        wprtp_forcing, z, dataset = self.getVarForCalculations('wprtp_forcing', dataset_override, fill_zeros=True)
        wprtp_bp, z, dataset = self.getVarForCalculations('wprtp_bp', dataset_override, fill_zeros=True)
        wprtp_ma, z, dataset = self.getVarForCalculations('wprtp_ma', dataset_override, fill_zeros=True)
        wprtp_pd, z, dataset = self.getVarForCalculations('wprtp_pd', dataset_override, fill_zeros=True)

        output_data = wprtp_bt - (
                    wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_pd + wprtp_forcing)

        return output_data, z

    def getWp2Residual(self, dataset_override=None):
        '''


        wp2_bt - (wp2_ma + wp2_ta + wp2_tp + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_mfl + wp2_cl + wp2_sicl + wp2_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wp2_sf, z, dataset = self.getVarForCalculations('wp2_sf', dataset_override, fill_zeros=True)
        wp2_cl, z, dataset = self.getVarForCalculations('wp2_cl', dataset_override, fill_zeros=True)
        wp2_ac, z, dataset = self.getVarForCalculations('wp2_ac', dataset_override, fill_zeros=True)
        wp2_pr1, z, dataset = self.getVarForCalculations('wp2_pr1', dataset_override, fill_zeros=True)
        wp2_pr3, z, dataset = self.getVarForCalculations('wp2_pr3', dataset_override, fill_zeros=True)
        wp2_pr2, z, dataset = self.getVarForCalculations('wp2_pr2', dataset_override, fill_zeros=True)
        wp2_dp1, z, dataset = self.getVarForCalculations('wp2_dp1', dataset_override, fill_zeros=True)
        wp2_dp2, z, dataset = self.getVarForCalculations('wp2_dp2', dataset_override, fill_zeros=True)
        wp2_bt, z, dataset = self.getVarForCalculations('wp2_bt', dataset_override, fill_zeros=True)
        wp2_ta, z, dataset = self.getVarForCalculations('wp2_ta', dataset_override, fill_zeros=True)
        wp2_splat, z, dataset = self.getVarForCalculations('wp2_splat', dataset_override, fill_zeros=True)
        wp2_bp, z, dataset = self.getVarForCalculations('wp2_bp', dataset_override, fill_zeros=True)
        wp2_ma, z, dataset = self.getVarForCalculations('wp2_ma', dataset_override, fill_zeros=True)
        wp2_pd, z, dataset = self.getVarForCalculations('wp2_pd', dataset_override, fill_zeros=True)

        output_data = wp2_bt - (
                    wp2_ma + wp2_ta + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_dp2 + wp2_cl + wp2_pd + wp2_sf + wp2_splat)

        return output_data, z

    def getWp3Residual(self, dataset_override=None):
        '''


        wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl+wp3_splat)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        wp3_bp1, z, dataset = self.getVarForCalculations('wp3_bp1', dataset_override, fill_zeros=True)
        wp3_bp2, z, dataset = self.getVarForCalculations('wp3_bp2', dataset_override, fill_zeros=True)
        wp3_cl, z, dataset = self.getVarForCalculations('wp3_cl', dataset_override, fill_zeros=True)
        wp3_ac, z, dataset = self.getVarForCalculations('wp3_ac', dataset_override, fill_zeros=True)
        wp3_pr1, z, dataset = self.getVarForCalculations('wp3_pr1', dataset_override, fill_zeros=True)
        wp3_pr3, z, dataset = self.getVarForCalculations('wp3_pr3', dataset_override, fill_zeros=True)
        wp3_pr2, z, dataset = self.getVarForCalculations('wp3_pr2', dataset_override, fill_zeros=True)
        wp3_dp1, z, dataset = self.getVarForCalculations('wp3_dp1', dataset_override, fill_zeros=True)
        wp3_bt, z, dataset = self.getVarForCalculations('wp3_bt', dataset_override, fill_zeros=True)
        wp3_ta, z, dataset = self.getVarForCalculations('wp3_ta', dataset_override, fill_zeros=True)
        wp3_splat, z, dataset = self.getVarForCalculations('wp3_splat', dataset_override, fill_zeros=True)
        wp3_ma, z, dataset = self.getVarForCalculations('wp3_ma', dataset_override, fill_zeros=True)
        wp3_tp, z, dataset = self.getVarForCalculations('wp3_tp', dataset_override, fill_zeros=True)

        output_data = wp3_bt - (
                    wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl + wp3_splat)

        return output_data, z

    def getThlp2Residual(self, dataset_override=None):
        '''


        thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        thlp2_cl, z, dataset = self.getVarForCalculations('thlp2_cl', dataset_override, fill_zeros=True)
        thlp2_dp2, z, dataset = self.getVarForCalculations('thlp2_dp2', dataset_override, fill_zeros=True)
        thlp2_forcing, z, dataset = self.getVarForCalculations('thlp2_forcing', dataset_override, fill_zeros=True)
        thlp2_sf, z, dataset = self.getVarForCalculations('thlp2_sf', dataset_override, fill_zeros=True)
        thlp2_dp1, z, dataset = self.getVarForCalculations('thlp2_dp1', dataset_override, fill_zeros=True)
        thlp2_bt, z, dataset = self.getVarForCalculations('thlp2_bt', dataset_override, fill_zeros=True)
        thlp2_ta, z, dataset = self.getVarForCalculations('thlp2_ta', dataset_override, fill_zeros=True)
        thlp2_pd, z, dataset = self.getVarForCalculations('thlp2_pd', dataset_override, fill_zeros=True)
        thlp2_ma, z, dataset = self.getVarForCalculations('thlp2_ma', dataset_override, fill_zeros=True)
        thlp2_tp, z, dataset = self.getVarForCalculations('thlp2_tp', dataset_override, fill_zeros=True)

        output_data = thlp2_bt - (
                    thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)

        return output_data, z

    def getRtp2Residual(self, dataset_override=None):
        '''


        rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtp2_cl, z, dataset = self.getVarForCalculations('rtp2_cl', dataset_override, fill_zeros=True)
        rtp2_dp2, z, dataset = self.getVarForCalculations('rtp2_dp2', dataset_override, fill_zeros=True)
        rtp2_forcing, z, dataset = self.getVarForCalculations('rtp2_forcing', dataset_override, fill_zeros=True)
        rtp2_sf, z, dataset = self.getVarForCalculations('rtp2_sf', dataset_override, fill_zeros=True)
        rtp2_dp1, z, dataset = self.getVarForCalculations('rtp2_dp1', dataset_override, fill_zeros=True)
        rtp2_bt, z, dataset = self.getVarForCalculations('rtp2_bt', dataset_override, fill_zeros=True)
        rtp2_ta, z, dataset = self.getVarForCalculations('rtp2_ta', dataset_override, fill_zeros=True)
        rtp2_pd, z, dataset = self.getVarForCalculations('rtp2_pd', dataset_override, fill_zeros=True)
        rtp2_ma, z, dataset = self.getVarForCalculations('rtp2_ma', dataset_override, fill_zeros=True)
        rtp2_tp, z, dataset = self.getVarForCalculations('rtp2_tp', dataset_override, fill_zeros=True)

        output_data = rtp2_bt - (
                    rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)

        return output_data, z

    def getRtpthlpResidual(self, dataset_override=None):
        '''


        rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_pd + rtpthlp_sf + rtpthlp_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rtpthlp_cl, z, dataset = self.getVarForCalculations('rtpthlp_cl', dataset_override, fill_zeros=True)
        rtpthlp_dp2, z, dataset = self.getVarForCalculations('rtpthlp_dp2', dataset_override, fill_zeros=True)
        rtpthlp_forcing, z, dataset = self.getVarForCalculations('rtpthlp_forcing', dataset_override, fill_zeros=True)
        rtpthlp_sf, z, dataset = self.getVarForCalculations('rtpthlp_sf', dataset_override, fill_zeros=True)
        rtpthlp_dp1, z, dataset = self.getVarForCalculations('rtpthlp_dp1', dataset_override, fill_zeros=True)
        rtpthlp_bt, z, dataset = self.getVarForCalculations('rtpthlp_bt', dataset_override, fill_zeros=True)
        rtpthlp_ta, z, dataset = self.getVarForCalculations('rtpthlp_ta', dataset_override, fill_zeros=True)
        rtpthlp_tp2, z, dataset = self.getVarForCalculations('rtpthlp_tp2', dataset_override, fill_zeros=True)
        rtpthlp_ma, z, dataset = self.getVarForCalculations('rtpthlp_ma', dataset_override, fill_zeros=True)
        rtpthlp_tp1, z, dataset = self.getVarForCalculations('rtpthlp_tp1', dataset_override, fill_zeros=True)

        output_data = rtpthlp_bt - (
                    rtpthlp_ma + rtpthlp_ta + rtpthlp_tp1 + rtpthlp_tp2 + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_sf + rtpthlp_forcing)

        return output_data, z

    def getUpwpResidual(self, dataset_override=None):
        '''


        upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_dp1 + upwp_dp2 + upwp_cl + upwp_pd + upwp_sf + upwp_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        upwp_cl, z, dataset = self.getVarForCalculations('upwp_cl', dataset_override, fill_zeros=True)
        upwp_tp, z, dataset = self.getVarForCalculations('upwp_tp', dataset_override, fill_zeros=True)
        upwp_ac, z, dataset = self.getVarForCalculations('upwp_ac', dataset_override, fill_zeros=True)
        upwp_bp, z, dataset = self.getVarForCalculations('upwp_bp', dataset_override, fill_zeros=True)
        upwp_dp1, z, dataset = self.getVarForCalculations('upwp_dp1', dataset_override, fill_zeros=True)
        upwp_bt, z, dataset = self.getVarForCalculations('upwp_bt', dataset_override, fill_zeros=True)
        upwp_ta, z, dataset = self.getVarForCalculations('upwp_ta', dataset_override, fill_zeros=True)
        upwp_pr1, z, dataset = self.getVarForCalculations('upwp_pr1', dataset_override, fill_zeros=True)
        upwp_pr2, z, dataset = self.getVarForCalculations('upwp_pr2', dataset_override, fill_zeros=True)
        upwp_pr3, z, dataset = self.getVarForCalculations('upwp_pr3', dataset_override, fill_zeros=True)
        upwp_pr4, z, dataset = self.getVarForCalculations('upwp_pr4', dataset_override, fill_zeros=True)
        upwp_mfl, z, dataset = self.getVarForCalculations('upwp_mfl', dataset_override, fill_zeros=True)
        upwp_ma, z, dataset = self.getVarForCalculations('upwp_ma', dataset_override, fill_zeros=True)

        output_data = upwp_bt - (
                    upwp_ma + upwp_ta + upwp_tp + upwp_ac + upwp_bp + upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4 + upwp_dp1 + upwp_mfl + upwp_cl)

        return output_data, z

    def getVpwpResidual(self, dataset_override=None):
        '''


        vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_dp1 + vpwp_dp2 + vpwp_cl + vpwp_pd + vpwp_sf + vpwp_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        vpwp_cl, z, dataset = self.getVarForCalculations('vpwp_cl', dataset_override, fill_zeros=True)
        vpwp_tp, z, dataset = self.getVarForCalculations('vpwp_tp', dataset_override, fill_zeros=True)
        vpwp_ac, z, dataset = self.getVarForCalculations('vpwp_ac', dataset_override, fill_zeros=True)
        vpwp_bp, z, dataset = self.getVarForCalculations('vpwp_bp', dataset_override, fill_zeros=True)
        vpwp_dp1, z, dataset = self.getVarForCalculations('vpwp_dp1', dataset_override, fill_zeros=True)
        vpwp_bt, z, dataset = self.getVarForCalculations('vpwp_bt', dataset_override, fill_zeros=True)
        vpwp_ta, z, dataset = self.getVarForCalculations('vpwp_ta', dataset_override, fill_zeros=True)
        vpwp_pr1, z, dataset = self.getVarForCalculations('vpwp_pr1', dataset_override, fill_zeros=True)
        vpwp_pr2, z, dataset = self.getVarForCalculations('vpwp_pr2', dataset_override, fill_zeros=True)
        vpwp_pr3, z, dataset = self.getVarForCalculations('vpwp_pr3', dataset_override, fill_zeros=True)
        vpwp_pr4, z, dataset = self.getVarForCalculations('vpwp_pr4', dataset_override, fill_zeros=True)
        vpwp_mfl, z, dataset = self.getVarForCalculations('vpwp_mfl', dataset_override, fill_zeros=True)
        vpwp_ma, z, dataset = self.getVarForCalculations('vpwp_ma', dataset_override, fill_zeros=True)

        output_data = vpwp_bt - (
                    vpwp_ma + vpwp_ta + vpwp_tp + vpwp_ac + vpwp_bp + vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4 + vpwp_dp1 + vpwp_mfl + vpwp_cl)

        return output_data, z


    def getRrmFillClip(self, dataset_override=None):
        '''


        rrm_hf + rrm_wvhf + rrm_cl
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rrm_hf, z, dataset = self.getVarForCalculations('rrm_hf', dataset_override, fill_zeros=True)
        rrm_wvhf, z, dataset = self.getVarForCalculations('rrm_wvhf', dataset_override, fill_zeros=True)
        rrm_cl, z, dataset = self.getVarForCalculations('rrm_cl', dataset_override, fill_zeros=True)

        output_data = rrm_hf + rrm_wvhf + rrm_cl

        return output_data, z

    def getRrmResidual(self, dataset_override=None):
        '''


        rrm_bt - (rrm_ma + rrm_sd + rrm_ta + rrm_ts + rrm_hf + rrm_wvhf + rrm_cl + rrm_mc)            :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        rrm_bt, z, dataset = self.getVarForCalculations('rrm_bt', dataset_override, fill_zeros=True)
        rrm_ma, z, dataset = self.getVarForCalculations('rrm_ma', dataset_override, fill_zeros=True)
        rrm_sd, z, dataset = self.getVarForCalculations('rrm_sd', dataset_override, fill_zeros=True)
        rrm_ta, z, dataset = self.getVarForCalculations('rrm_ta', dataset_override, fill_zeros=True)
        rrm_ts, z, dataset = self.getVarForCalculations('rrm_ts', dataset_override, fill_zeros=True)
        rrm_hf, z, dataset = self.getVarForCalculations('rrm_hf', dataset_override, fill_zeros=True)
        rrm_wvhf, z, dataset = self.getVarForCalculations('rrm_wvhf', dataset_override, fill_zeros=True)
        rrm_cl, z, dataset = self.getVarForCalculations('rrm_cl', dataset_override, fill_zeros=True)
        rrm_mc, z, dataset = self.getVarForCalculations('rrm_mc', dataset_override, fill_zeros=True)

        output_data = rrm_bt - (rrm_ma + rrm_sd + rrm_ta + rrm_ts + rrm_hf + rrm_wvhf + rrm_cl + rrm_mc)

        return output_data, z

    def getNrmResidual(self, dataset_override=None):
        '''


        Nrm_bt - (Nrm_ma + Nrm_sd + Nrm_ta + Nrm_ts + Nrm_cl + Nrm_mc)
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override, fill_zeros=True)
        Nrm_bt, z, dataset = self.getVarForCalculations('Nrm_bt', dataset_override, fill_zeros=True)
        Nrm_ma, z, dataset = self.getVarForCalculations('Nrm_ma', dataset_override, fill_zeros=True)
        Nrm_sd, z, dataset = self.getVarForCalculations('Nrm_sd', dataset_override, fill_zeros=True)
        Nrm_ta, z, dataset = self.getVarForCalculations('Nrm_ta', dataset_override, fill_zeros=True)
        Nrm_ts, z, dataset = self.getVarForCalculations('Nrm_ts', dataset_override, fill_zeros=True)
        Nrm_cl, z, dataset = self.getVarForCalculations('Nrm_cl', dataset_override, fill_zeros=True)
        Nrm_mc, z, dataset = self.getVarForCalculations('Nrm_mc', dataset_override, fill_zeros=True)

        output_data = Nrm_bt - (Nrm_ma + Nrm_sd + Nrm_ta + Nrm_ts + Nrm_cl + Nrm_mc)

        return output_data, z