"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupBaseBudgets(VariableGroup):
    """

    """

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, wrf_datasets=None):
        self.name = "base variables budgets"

        thlm_budget_lines = [
            {'var_names': ['thlm_bt'], 'legend_label': 'thlm_bt'},
            {'var_names': ['thlm_ma'], 'legend_label': 'thlm_ma'},
            {'var_names': ['thlm_ta'], 'legend_label': 'thlm_ta'},
            {'var_names': ['thlm_mc'], 'legend_label': 'thlm_mc'},
            {'var_names': ['thlm_clipping'], 'legend_label': 'thlm_clipping', 'clubb_calc': self.getThlmClipping,
             'e3sm_calc': self.getThlmClipping, 'wrf_calc': self.getThlmClipping},
            {'var_names': ['radht'], 'legend_label': 'radht'},
            {'var_names': ['ls_forcing'], 'legend_label': 'thlm_ls_forcing', 'clubb_calc': self.getThlmLsforcing,
             'e3sm_calc': self.getThlmLsforcing, 'wrf_calc': self.getThlmLsforcing},
            {'var_names': ['thlm_residual'], 'legend_label': 'thlm_residual', 'clubb_calc': self.getThlmResidual,
             'e3sm_calc': self.getThlmResidual, 'wrf_calc': self.getThlmResidual},

        ]

        rtm_budget_lines = [
            {'var_names': ['rtm_bt'], 'legend_label': 'rtm_bt'},
            {'var_names': ['rtm_ma'], 'legend_label': 'rtm_ma'},
            {'var_names': ['rtm_ta'], 'legend_label': 'rtm_ta'},
            {'var_names': ['rtm_mc'], 'legend_label': 'rtm_mc'},
            {'var_names': ['rtm_clipping'], 'legend_label': 'rtm_clipping', 'clubb_calc': self.getRtmClipping,
             'e3sm_calc': self.getRtmClipping, 'wrf_calc': self.getRtmClipping},
            {'var_names': ['rtm_pd'], 'legend_label': 'rtm_pd'},
            {'var_names': ['ls_forcing'], 'legend_label': 'ls_forcing', 'clubb_calc': self.getRtmForcing,
             'e3sm_calc': self.getRtmForcing, 'wrf_calc': self.getRtmForcing},
            {'var_names': ['rtm_residual'], 'legend_label': 'rtm_residual', 'clubb_calc': self.getRtmResidual,
             'e3sm_calc': self.getRtmResidual, 'wrf_calc': self.getRtmResidual},

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
             'clubb_calc': self.getWpthlpResidual, 'e3sm_calc': self.getWpthlpResidual,
             'wrf_calc': self.getWpthlpResidual},

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
            {'var_names': ['wprtp_residual'], 'legend_label': 'wprtp_residual', 'clubb_calc': self.getWprtpResidual,
             'e3sm_calc': self.getWprtpResidual, 'wrf_calc': self.getWprtpResidual},
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
            {'var_names': ['wp2_residual'], 'legend_label': 'wp2_residual', 'clubb_calc': self.getWp2Residual,
             'e3sm_calc': self.getWp2Residual, 'wrf_calc': self.getWp2Residual},
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
            {'var_names': ['wp3_residual'], 'legend_label': 'wp3_residual', 'clubb_calc': self.getWp3Residual,
             'e3sm_calc': self.getWp3Residual, 'wrf_calc': self.getWp3Residual},
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
            {'var_names': ['thlp2_residual'], 'legend_label': 'thlp2_residual', 'clubb_calc': self.getThlp2Residual,
             'e3sm_calc': self.getThlp2Residual, 'wrf_calc': self.getThlp2Residual},
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
            {'var_names': ['rtp2_residual'], 'legend_label': 'rtp2_residual', 'clubb_calc': self.getRtp2Residual,
             'e3sm_calc': self.getRtp2Residual, 'wrf_calc': self.getRtp2Residual},
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
             'clubb_calc': self.getRtpthlpResidual, 'e3sm_calc': self.getRtpthlpResidual,
             'wrf_calc': self.getRtpthlpResidual},
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
            {'var_names': ['upwp_residual'], 'legend_label': 'upwp_residual', 'clubb_calc': self.getUpwpResidual,
             'e3sm_calc': self.getUpwpResidual, 'wrf_calc': self.getUpwpResidual},
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
            {'var_names': ['vpwp_residual'], 'legend_label': 'vpwp_residual', 'clubb_calc': self.getVpwpResidual,
             'e3sm_calc': self.getVpwpResidual, 'wrf_calc': self.getVpwpResidual},
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

        up2_budget_lines = [
            {'var_names': ['up2_bt'], 'legend_label': 'up2_bt'},
            {'var_names': ['up2_ma'], 'legend_label': 'up2_ma'},
            {'var_names': ['up2_ta'], 'legend_label': 'up2_ta'},
            {'var_names': ['up2_tp'], 'legend_label': 'up2_tp'},
            {'var_names': ['up2_dp1'], 'legend_label': 'up2_dp1'},
            {'var_names': ['up2_dp2'], 'legend_label': 'up2_dp2'},
            {'var_names': ['up2_pr1'], 'legend_label': 'up2_pr1'},
            {'var_names': ['up2_pr2'], 'legend_label': 'up2_pr2'},
            {'var_names': ['up2_sdmp'], 'legend_label': 'up2_sdmp'},
            {'var_names': ['up2_cl'], 'legend_label': 'up2_cl'},
            {'var_names': ['up2_pd'], 'legend_label': 'up2_pd'},
            {'var_names': ['up2_sf'], 'legend_label': 'up2_sf'},
            {'var_names': ['up2_splat'], 'legend_label': 'up2_splat'}
        ]

        vp2_budget_lines = [
            {'var_names': ['vp2_bt'], 'legend_label': 'vp2_bt'},
            {'var_names': ['vp2_ma'], 'legend_label': 'vp2_ma'},
            {'var_names': ['vp2_ta'], 'legend_label': 'vp2_ta'},
            {'var_names': ['vp2_tp'], 'legend_label': 'vp2_tp'},
            {'var_names': ['vp2_dp1'], 'legend_label': 'vp2_dp1'},
            {'var_names': ['vp2_dp2'], 'legend_label': 'vp2_dp2'},
            {'var_names': ['vp2_pr1'], 'legend_label': 'vp2_pr1'},
            {'var_names': ['vp2_pr2'], 'legend_label': 'vp2_pr2'},
            {'var_names': ['vp2_sdmp'], 'legend_label': 'vp2_sdmp'},
            {'var_names': ['vp2_cl'], 'legend_label': 'vp2_cl'},
            {'var_names': ['vp2_pd'], 'legend_label': 'vp2_pd'},
            {'var_names': ['vp2_sf'], 'legend_label': 'vp2_sf'},
            {'var_names': ['vp2_splat'], 'legend_label': 'vp2_splat'}
        ]

        rrm_budget_lines = [
            {'var_names': ['rrm_bt'], 'legend_label': 'rrm_bt'},
            {'var_names': ['rrm_ma'], 'legend_label': 'rrm_ma'},
            {'var_names': ['rrm_sd'], 'legend_label': 'rrm_sd'},
            {'var_names': ['rrm_ta'], 'legend_label': 'rrm_ta'},
            {'var_names': ['rrm_ts'], 'legend_label': 'rrm_ts'},
            {'var_names': ['rrm_fill_clip'], 'legend_label': 'rrm_fill_clip', 'clubb_calc': self.getRrmFillClip,
             'e3sm_calc': self.getRrmFillClip, 'wrf_calc': self.getRrmFillClip},
            {'var_names': ['rrm_mc'], 'legend_label': 'rrm_mc'},
            {'var_names': ['rrm_residual'], 'legend_label': 'rrm_residual', 'clubb_calc': self.getRrmResidual,
             'e3sm_calc': self.getRrmResidual, 'wrf_calc': self.getRrmResidual},
        ]

        Nrm_budget_lines = [
            {'var_names': ['Nrm_bt'], 'legend_label': 'Nrm_bt'},
            {'var_names': ['Nrm_ma'], 'legend_label': 'Nrm_ma'},
            {'var_names': ['Nrm_sd'], 'legend_label': 'Nrm_sd'},
            {'var_names': ['Nrm_ta'], 'legend_label': 'Nrm_ta'},
            {'var_names': ['Nrm_ts'], 'legend_label': 'Nrm_ts'},
            {'var_names': ['Nrm_cl'], 'legend_label': 'Nrm_cl'},
            {'var_names': ['Nrm_mc'], 'legend_label': 'Nrm_mc'},
            {'var_names': ['Nrm_residual'], 'legend_label': 'Nrm_residual', 'clubb_calc': self.getNrmResidual,
             'e3sm_calc': self.getNrmResidual, 'wrf_calc': self.getNrmResidual},
        ]

        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['thlm'],
                'sam': ['thlm'],
                'coamps': ['thlm'],
                'r408': ['thlm'],
                'hoc': ['thlm'],
                'e3sm': ['thlm'],
                'cam': ['thlm'],
                'wrf': ['thlm'],
                },
             'lines': thlm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
            },
            {'var_names':
                {
                'clubb': ['rtm'],
                'sam': ['rtm'],
                'coamps': ['rtm'],
                'r408': ['rtm'],
                'hoc': ['rtm'],
                'e3sm': ['rtm'],
                'cam': ['rtm'],
                'wrf': ['rtm'],
                },
             'lines': rtm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
            },
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
            },
            {'var_names':
                {
                'clubb': ['um'],
                'sam': ['um'],
                'coamps': ['um'],
                'r408': ['um'],
                'hoc': ['um'],
                'e3sm': ['um'],
                'cam': ['um'],
                'wrf': ['um'],
                },
             'lines': um_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
            },
            {'var_names':
                {
                'clubb': ['vm'],
                'sam': ['vm'],
                'coamps': ['vm'],
                'r408': ['vm'],
                'hoc': ['vm'],
                'e3sm': ['vm'],
                'cam': ['vm'],
                'wrf': ['vm'],
                },
             'lines': vm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
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
	    },
            {'var_names':
                {
                'clubb': ['rrm'],
                'sam': ['rrm'],
                'coamps': ['rrm'],
                'r408': ['rrm'],
                'hoc': ['rrm'],
                'e3sm': ['rrm'],
                'cam': ['rrm'],
                'wrf': ['rrm'],
                },
             'lines': rrm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True,
             },
            {'var_names':
                {
                'clubb': ['Nrm'],
                'sam': ['Nrm'],
                'coamps': ['Nrm'],
                'r408': ['Nrm'],
                'hoc': ['Nrm'],
                'e3sm': ['Nrm'],
                'cam': ['Nrm'],
                'wrf': ['Nrm'],
                },
             'lines': Nrm_budget_lines, 'type': Panel.TYPE_BUDGET, 'centered': True
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets)

    def getThlmClipping(self, dataset_override=None):
        '''


        thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        :return:
        '''

        thlm_mfl, indep, dataset = self.getVarForCalculations('thlm_mfl', dataset_override)
        thlm_cl, indep, dataset = self.getVarForCalculations('thlm_cl', dataset)
        thlm_tacl, indep, dataset = self.getVarForCalculations('thlm_tacl', dataset)
        thlm_sdmp, indep, dataset = self.getVarForCalculations('thlm_sdmp', dataset)

        output_data = thlm_mfl + thlm_cl + thlm_tacl + thlm_sdmp

        return output_data, indep

    def getThlmLsforcing(self, dataset_override=None):
        '''


        thlm_forcing-radht-thlm_mc
        :return:
        '''
        thlm_forcing, indep, dataset = self.getVarForCalculations('thlm_forcing', dataset_override)
        radht, indep, dataset = self.getVarForCalculations('radht', dataset)
        thlm_mc, indep, dataset = self.getVarForCalculations('thlm_mc', dataset)

        output_data = thlm_forcing - radht - thlm_mc

        return output_data, indep

    def getThlmResidual(self, dataset_override=None):
        '''


        thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)
        :return:
        '''
        thlm_mfl, indep, dataset = self.getVarForCalculations('thlm_mfl', dataset_override)
        thlm_cl, indep, dataset = self.getVarForCalculations('thlm_cl', dataset)
        thlm_tacl, indep, dataset = self.getVarForCalculations('thlm_tacl', dataset)
        thlm_sdmp, indep, dataset = self.getVarForCalculations('thlm_sdmp', dataset)
        thlm_bt, indep, dataset = self.getVarForCalculations('thlm_bt', dataset)
        thlm_ta, indep, dataset = self.getVarForCalculations('thlm_ta', dataset)
        thlm_forcing, indep, dataset = self.getVarForCalculations('thlm_forcing', dataset)
        thlm_ma, indep, dataset = self.getVarForCalculations('thlm_ma', dataset)

        output_data = thlm_bt - (thlm_ma + thlm_ta + thlm_mfl + thlm_cl + thlm_tacl + thlm_sdmp + thlm_forcing)

        return output_data, indep

    def getRtmClipping(self, dataset_override=None):
        '''


        rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp
        :return:
        '''
        rtm_mfl, indep, dataset = self.getVarForCalculations('rtm_mfl', dataset_override)
        rtm_cl, indep, dataset = self.getVarForCalculations('rtm_cl', dataset)
        rtm_tacl, indep, dataset = self.getVarForCalculations('rtm_tacl', dataset)
        rtm_sdmp, indep, dataset = self.getVarForCalculations('rtm_sdmp', dataset)

        output_data = rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp

        return output_data, indep

    def getRtmForcing(self, dataset_override=None):
        '''


        rtm_forcing - rtm_mc
        :return:
        '''
        rtm_mc, indep, dataset = self.getVarForCalculations('rtm_mc', dataset_override)
        rtm_forcing, indep, dataset = self.getVarForCalculations('rtm_forcing', dataset)

        output_data = rtm_forcing - rtm_mc

        return output_data, indep

    def getRtmResidual(self, dataset_override=None):
        '''


        rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)
        :return:
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

    def getWpthlpResidual(self, dataset_override=None):
        '''


        wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 +
        wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)
        :return:
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

    def getWprtpResidual(self, dataset_override=None):
        '''


        wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 +
        wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_forcing)
        :return:
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

    def getWp2Residual(self, dataset_override=None):
        '''


        wp2_bt - (wp2_ma + wp2_ta + wp2_tp + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 +
        wp2_mfl + wp2_cl + wp2_sicl + wp2_forcing)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp2_sf, indep, dataset = self.getVarForCalculations('wp2_sf', dataset_override)
        wp2_cl, indep, dataset = self.getVarForCalculations('wp2_cl', dataset)
        wp2_ac, indep, dataset = self.getVarForCalculations('wp2_ac', dataset)
        wp2_pr1, indep, dataset = self.getVarForCalculations('wp2_pr1', dataset)
        wp2_pr3, indep, dataset = self.getVarForCalculations('wp2_pr3', dataset)
        wp2_pr2, indep, dataset = self.getVarForCalculations('wp2_pr2', dataset)
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
                wp2_cl + wp2_pd + wp2_sf + wp2_splat)

        return output_data, indep

    def getWp3Residual(self, dataset_override=None):
        '''


        wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 +
        wp3_dp1 + wp3_cl+wp3_splat)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        wp3_bp1, indep, dataset = self.getVarForCalculations('wp3_bp1', dataset_override)
        wp3_bp2, indep, dataset = self.getVarForCalculations('wp3_bp2', dataset)
        wp3_cl, indep, dataset = self.getVarForCalculations('wp3_cl', dataset)
        wp3_ac, indep, dataset = self.getVarForCalculations('wp3_ac', dataset)
        wp3_pr1, indep, dataset = self.getVarForCalculations('wp3_pr1', dataset)
        wp3_pr3, indep, dataset = self.getVarForCalculations('wp3_pr3', dataset)
        wp3_pr2, indep, dataset = self.getVarForCalculations('wp3_pr2', dataset)
        wp3_dp1, indep, dataset = self.getVarForCalculations('wp3_dp1', dataset)
        wp3_bt, indep, dataset = self.getVarForCalculations('wp3_bt', dataset)
        wp3_ta, indep, dataset = self.getVarForCalculations('wp3_ta', dataset)
        wp3_splat, indep, dataset = self.getVarForCalculations('wp3_splat', dataset)
        wp3_ma, indep, dataset = self.getVarForCalculations('wp3_ma', dataset)
        wp3_tp, indep, dataset = self.getVarForCalculations('wp3_tp', dataset)

        output_data = wp3_bt - (
                wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 +
                wp3_dp1 + wp3_cl + wp3_splat)

        return output_data, indep

    def getThlp2Residual(self, dataset_override=None):
        '''


        thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd +
        thlp2_sf + thlp2_forcing)
        :return:
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

    def getRtp2Residual(self, dataset_override=None):
        '''


        rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)
        :return:
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

    def getRtpthlpResidual(self, dataset_override=None):
        '''


        rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl +
        rtpthlp_pd + rtpthlp_sf + rtpthlp_forcing)
        :return:
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

    def getUpwpResidual(self, dataset_override=None):
        '''


        upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_dp1 + upwp_dp2 + upwp_cl + upwp_pd + upwp_sf + upwp_forcing)
        :return:
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

    def getVpwpResidual(self, dataset_override=None):
        '''


        vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_dp1 + vpwp_dp2 + vpwp_cl + vpwp_pd + vpwp_sf + vpwp_forcing)
        :return:
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

    def getRrmFillClip(self, dataset_override=None):
        '''


        rrm_hf + rrm_wvhf + rrm_cl
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rrm_hf, indep, dataset = self.getVarForCalculations('rrm_hf', dataset_override)
        rrm_wvhf, indep, dataset = self.getVarForCalculations('rrm_wvhf', dataset)
        rrm_cl, indep, dataset = self.getVarForCalculations('rrm_cl', dataset)

        output_data = rrm_hf + rrm_wvhf + rrm_cl

        return output_data, indep

    def getRrmResidual(self, dataset_override=None):
        '''


        rrm_bt - (rrm_ma + rrm_sd + rrm_ta + rrm_ts + rrm_hf + rrm_wvhf + rrm_cl + rrm_mc)
        :return:
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        rrm_bt, indep, dataset = self.getVarForCalculations('rrm_bt', dataset_override)
        rrm_ma, indep, dataset = self.getVarForCalculations('rrm_ma', dataset)
        rrm_sd, indep, dataset = self.getVarForCalculations('rrm_sd', dataset)
        rrm_ta, indep, dataset = self.getVarForCalculations('rrm_ta', dataset)
        rrm_ts, indep, dataset = self.getVarForCalculations('rrm_ts', dataset)
        rrm_hf, indep, dataset = self.getVarForCalculations('rrm_hf', dataset)
        rrm_wvhf, indep, dataset = self.getVarForCalculations('rrm_wvhf', dataset)
        rrm_cl, indep, dataset = self.getVarForCalculations('rrm_cl', dataset)
        rrm_mc, indep, dataset = self.getVarForCalculations('rrm_mc', dataset)

        output_data = rrm_bt - (rrm_ma + rrm_sd + rrm_ta + rrm_ts + rrm_hf + rrm_wvhf + rrm_cl + rrm_mc)

        return output_data, indep

    def getNrmResidual(self, dataset_override=None):
        '''


        Nrm_bt - (Nrm_ma + Nrm_sd + Nrm_ta + Nrm_ts + Nrm_cl + Nrm_mc)
        '''
        # z,z, dataset = self.getVarForCalculations('altitude', dataset_override)
        Nrm_bt, indep, dataset = self.getVarForCalculations('Nrm_bt', dataset_override)
        Nrm_ma, indep, dataset = self.getVarForCalculations('Nrm_ma', dataset)
        Nrm_sd, indep, dataset = self.getVarForCalculations('Nrm_sd', dataset)
        Nrm_ta, indep, dataset = self.getVarForCalculations('Nrm_ta', dataset)
        Nrm_ts, indep, dataset = self.getVarForCalculations('Nrm_ts', dataset)
        Nrm_cl, indep, dataset = self.getVarForCalculations('Nrm_cl', dataset)
        Nrm_mc, indep, dataset = self.getVarForCalculations('Nrm_mc', dataset)

        output_data = Nrm_bt - (Nrm_ma + Nrm_sd + Nrm_ta + Nrm_ts + Nrm_cl + Nrm_mc)

        return output_data, indep
