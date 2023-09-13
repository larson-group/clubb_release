'''
:author: Steffen Domke
:date: Mid 2019
TODO:   - Arrange lines so that styles match for different panels -> reduced momentum flux budgets
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupSamBudgets(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, verbose=False,
                 priority_vars=False, background_rcm=False, background_rcm_folder=None):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "sam budget variables"

        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)

        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000

        # ?? Which conversion factors? What is the variable description??
        hlp_budget_lines = [
            {'var_names': ['HLADV'], 'legend_label': 'HLADV',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLDFSN'], 'legend_label': 'HLDFSN',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLRAD'], 'legend_label': 'HLRAD',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLLAT'], 'legend_label': 'HLLAT',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['TTEND'], 'legend_label': 'TTEND',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['T_TNDCY'], 'legend_label': 'T_TNDCY'},
            {'var_names': ['HLSTOR'], 'legend_label': 'HLSTOR',
             'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HL_RES'], 'legend_label': 'HL_RES', 'sam_calc': self.getHlResidual},
        ]

        # Total water budget (same as rtm) (does not exist in _MICRO_M2005_UWM)
        qtp_budget_lines = [
            {'var_names': ['QTADV'], 'legend_label': 'QTADV',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTDFSN'], 'legend_label': 'QTDFSN',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSRC'], 'legend_label': 'QTSRC',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSINK'], 'legend_label': 'QTSINK',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTEND'], 'legend_label': 'QTEND',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSTOR'], 'legend_label': 'QTSTOR',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QV_TNDCY'], 'legend_label': 'QV_TNDCY',
             'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QT_RES'], 'legend_label': 'QT_RES', 'sam_calc': self.getQtResidual},
        ]

        # Vertical liquid water static energy flux
        tpwp_budget_lines = [
            {'var_names': ['TWGRAD'], 'legend_label': 'TWGRAD'},
            {'var_names': ['TWADV'], 'legend_label': 'TWADV'},
            {'var_names': ['TWDFSN'], 'legend_label': 'TWDFSN'},
            {'var_names': ['TWB+P'], 'legend_label': 'TWBUOY+TWPRES', 'sam_calc': self.getTwBuoyPlusPres},
            {'var_names': ['TWPREC'], 'legend_label': 'TWPREC'},
            {'var_names': ['TWRAD'], 'legend_label': 'TWRAD'},
            {'var_names': ['TWFORC'], 'legend_label': 'TWFORC'},
            {'var_names': ['TWBT'], 'legend_label': 'TWBT'},
            {'var_names': ['TW_RES'], 'legend_label': 'TW_RES', 'sam_calc': self.getTwResidual},
        ]

        tpwp_split_budget_lines = [
            {'var_names': ['TWGRAD'], 'legend_label': 'TWGRAD'},
            {'var_names': ['TWADV'], 'legend_label': 'TWADV'},
            {'var_names': ['TWDFSN'], 'legend_label': 'TWDFSN'},
            {'var_names': ['TWBUOY'], 'legend_label': 'TWBUOY'},
            {'var_names': ['TWPRES'], 'legend_label': 'TWPRES'},
            {'var_names': ['TWPREC'], 'legend_label': 'TWPREC'},
            {'var_names': ['TWRAD'], 'legend_label': 'TWRAD'},
            {'var_names': ['TWFORC'], 'legend_label': 'TWFORC'},
            {'var_names': ['TWBT'], 'legend_label': 'TWBT'},
            {'var_names': ['TW_RES'], 'legend_label': 'TWRES', 'sam_calc': self.getTwResidual},
        ]

        # Exists ? Vertical liquid water pot. temperature flux
        thlpwp_budget_lines = [
            {'var_names': ['THLWGRAD'], 'legend_label': 'THLWGRAD'},
            {'var_names': ['THLWADV'], 'legend_label': 'THLWADV'},
            {'var_names': ['THLWDFSN'], 'legend_label': 'THLWDFSN'},
            {'var_names': ['THLWB+P'], 'legend_label': 'THLWBUOY+PRES', 'sam_calc': self.getThlwBuoyPlusPres},
            {'var_names': ['THLWPREC'], 'legend_label': 'THLWPREC'},
            {'var_names': ['THLWRAD'], 'legend_label': 'THLWRAD'},
            {'var_names': ['THLWFORC'], 'legend_label': 'THLWFORC'},
            {'var_names': ['THLWBT'], 'legend_label': 'THLWBT'},
            {'var_names': ['THLW_RES'], 'legend_label': 'THLW_RES', 'sam_calc': self.getThlwResidual},
        ]

        thlpwp_split_budget_lines = [
            {'var_names': ['THLWGRAD'], 'legend_label': 'THLWGRAD'},
            {'var_names': ['THLWADV'], 'legend_label': 'THLWADV'},
            {'var_names': ['THLWDFSN'], 'legend_label': 'THLWDFSN'},
            {'var_names': ['THLWBUOY'], 'legend_label': 'THLWBUOY'},
            {'var_names': ['THLWPRES'], 'legend_label': 'THLWPRES'},
            {'var_names': ['THLWPREC'], 'legend_label': 'THLWPREC'},
            {'var_names': ['THLWRAD'], 'legend_label': 'THLWRAD'},
            {'var_names': ['THLWFORC'], 'legend_label': 'THLWFORC'},
            {'var_names': ['THLWBT'], 'legend_label': 'THLWBT'},
            {'var_names': ['THLW_RES'], 'legend_label': 'THLW_RES', 'sam_calc': self.getThlwResidual},
        ]

        # Vertical total water flux budget
        qpwp_budget_lines = [
            {'var_names': ['QWGRAD'], 'legend_label': 'QWGRAD'},
            {'var_names': ['QWADV'], 'legend_label': 'QWADV'},
            {'var_names': ['QWDFSN'], 'legend_label': 'QWDFSN'},
            {'var_names': ['QWB+P'], 'legend_label': 'QWBUOY+QWPRES', 'sam_calc': self.getQwBuoyPlusPres},
            {'var_names': ['QWPREC'], 'legend_label': 'QWPREC'},
            {'var_names': ['QWFORC'], 'legend_label': 'QWFORC'},
            {'var_names': ['QWBT'], 'legend_label': 'QWBT'},
            {'var_names': ['QW_RES'], 'legend_label': 'QW_RES', 'sam_calc': self.getQwResidual},
        ]

        qpwp_split_budget_lines = [
            {'var_names': ['QWGRAD'], 'legend_label': 'QWGRAD'},
            {'var_names': ['QWADV'], 'legend_label': 'QWADV'},
            {'var_names': ['QWDFSN'], 'legend_label': 'QWDFSN'},
            {'var_names': ['QWBUOY'], 'legend_label': 'QWBUOY'},
            {'var_names': ['QWPRES'], 'legend_label': 'QWPRES'},
            {'var_names': ['QWPREC'], 'legend_label': 'QWPREC'},
            {'var_names': ['QWFORC'], 'legend_label': 'QWFORC'},
            {'var_names': ['QWBT'], 'legend_label': 'QWBT'},
            {'var_names': ['QW_RES'], 'legend_label': 'QW_RES', 'sam_calc': self.getQwResidual},
        ]

        qtogpwp_budget_lines = [
            {'var_names': ['QTOGWGRAD'], 'legend_label': 'QTOGWGRAD'},
            {'var_names': ['QTOGWADV'], 'legend_label': 'QTOGWADV'},
            {'var_names': ['QTOGWDFSN'], 'legend_label': 'QTOGWDFSN'},
            {'var_names': ['QTOGWB+P'], 'legend_label': 'QTOGWBUOY+PRES',
             'sam_calc': self.getQtogwBuoyPlusPres},
            {'var_names': ['QTOGWPREC'], 'legend_label': 'QTOGWPREC'},
            {'var_names': ['QTOGWFORC'], 'legend_label': 'QTOGWFORC'},
            {'var_names': ['QTOGWBT'], 'legend_label': 'QTOGWBT'},
            {'var_names': ['QTOGW_RES'], 'legend_label': 'QTOGW_RES', 'sam_calc': self.getQtogwResidual},
        ]

        qtogpwp_split_budget_lines = [
            {'var_names': ['QTOGWGRAD'], 'legend_label': 'QTOGWGRAD'},
            {'var_names': ['QTOGWADV'], 'legend_label': 'QTOGWADV'},
            {'var_names': ['QTOGWDFSN'], 'legend_label': 'QTOGWDFSN'},
            {'var_names': ['QTOGWBUOY'], 'legend_label': 'QTOGWBUOY'},
            {'var_names': ['QTOGWPRES'], 'legend_label': 'QTOGWPRES'},
            {'var_names': ['QTOGWPREC'], 'legend_label': 'QTOGWPREC'},
            {'var_names': ['QTOGWFORC'], 'legend_label': 'QTOGWFORC'},
            {'var_names': ['QTOGWBT'], 'legend_label': 'QTOGWBT'},
            {'var_names': ['QTOGW_RES'], 'legend_label': 'QTOGW_RES', 'sam_calc': self.getQtogwResidual},
        ]

        tp2_budget_lines = [
            {'var_names': ['T2ADVTR'], 'legend_label': 'T2ADVTR'},
            {'var_names': ['T2GRAD'], 'legend_label': 'T2GRAD'},
            {'var_names': ['T2DISSIP'], 'legend_label': 'T2DISSIP'},
            {'var_names': ['T2DIFTR'], 'legend_label': 'T2DIFTR'},
            {'var_names': ['T2PREC'], 'legend_label': 'T2PREC'},
            {'var_names': ['T2RAD'], 'legend_label': 'T2RAD'},
            {'var_names': ['T2FORC'], 'legend_label': 'T2FORC'},
            {'var_names': ['T2BT'], 'legend_label': 'T2BT'},
            {'var_names': ['T2_RES'], 'legend_label': 'T2_RES', 'sam_calc': self.getT2Residual},
        ]

        # Exists?
        thlp2_budget_lines = [
            {'var_names': ['THL2ADVTR'], 'legend_label': 'THL2ADVTR'},
            {'var_names': ['THL2GRAD'], 'legend_label': 'THL2GRAD'},
            {'var_names': ['THL2DISSIP'], 'legend_label': 'THL2DISSIP'},
            {'var_names': ['THL2DIFTR'], 'legend_label': 'THL2DIFTR'},
            {'var_names': ['THL2PREC'], 'legend_label': 'THL2PREC'},
            {'var_names': ['THL2RAD'], 'legend_label': 'THL2RAD'},
            {'var_names': ['THL2FORC'], 'legend_label': 'THL2FORC'},
            {'var_names': ['THL2BT'], 'legend_label': 'THL2BT'},
            {'var_names': ['THL2_RES'], 'legend_label': 'THL2_RES', 'sam_calc': self.getThl2Residual},
        ]

        # Exists?
        rtp2_budget_lines = [
            {'var_names': ['Q2ADVTR'], 'legend_label': 'Q2ADVTR'},
            {'var_names': ['Q2GRAD'], 'legend_label': 'Q2GRAD'},
            {'var_names': ['Q2DISSIP'], 'legend_label': 'Q2DISSIP'},
            {'var_names': ['Q2DIFTR'], 'legend_label': 'Q2DIFTR'},
            {'var_names': ['Q2PREC'], 'legend_label': 'Q2PREC'},
            {'var_names': ['Q2FORC'], 'legend_label': 'Q2FORC'},
            {'var_names': ['Q2BT'], 'legend_label': 'Q2BT'},
            {'var_names': ['Q2_RES'], 'legend_label': 'Q2_RES', 'sam_calc': self.getQt2Residual},
        ]

        qtogp2_budget_lines = [
            {'var_names': ['QTOG2ADVTR'], 'legend_label': 'QTOG2ADVTR'},
            {'var_names': ['QTOG2GRAD'], 'legend_label': 'QTOG2GRAD'},
            {'var_names': ['QTOG2DISSIP'], 'legend_label': 'QTOG2DISSIP'},
            {'var_names': ['QTOG2DIFTR'], 'legend_label': 'QTOG2DIFTR'},
            {'var_names': ['QTOG2PREC'], 'legend_label': 'QTOG2PREC'},
            {'var_names': ['QTOG2FORC'], 'legend_label': 'QTOG2FORC'},
            {'var_names': ['QTOG2BT'], 'legend_label': 'QTOG2BT'},
            {'var_names': ['QTOG2_RES'], 'legend_label': 'QTOG2_RES',
             'sam_calc': self.getQtog2Residual},
        ]

        qpthlp_budget_lines = [
            {'var_names': ['QTHLADV'], 'legend_label': 'QTHLADV'},
            {'var_names': ['QTHLGRAD'], 'legend_label': 'QTHLGRAD'},
            {'var_names': ['QTHLDISSIP'], 'legend_label': 'QTHLDISSIP'},
            {'var_names': ['QTHLDIFTR'], 'legend_label': 'QTHLDIFTR'},
            {'var_names': ['QTHLPREC'], 'legend_label': 'QTHLPREC'},
            {'var_names': ['QTHLRAD'], 'legend_label': 'QTHLRAD'},
            {'var_names': ['QTHLFORC'], 'legend_label': 'QTHLFORC'},
            {'var_names': ['QTHLBT'], 'legend_label': 'QTHLBT'},
            {'var_names': ['QTHL_RES'], 'legend_label': 'QTHL_RES', 'sam_calc': self.getQThlResidual},
        ]

        ## Momentum flux budgets ##
        tke_budget_lines = [
            {'var_names': ['ADVTR'], 'legend_label': 'ADVTR'},
            {'var_names': ['SHEAR'], 'legend_label': 'SHEAR'},
            {'var_names': ['BUOYA'], 'legend_label': 'BUOYA'},
            {'var_names': ['PRESSTR'], 'legend_label': 'PRESSTR'},
            {'var_names': ['SDMP'], 'legend_label': 'SDMP'},
            {'var_names': ['DIFTR+DISS'], 'legend_label': 'DISSIP+DIFTR',
             'sam_calc': self.getTkeDissPlusDfsn},
            {'var_names': ['BT'], 'legend_label': 'BT'},
            {'var_names': ['TKE_RES'], 'legend_label': 'TKE_RES', 'sam_calc': self.getTkeResidual},
        ]

        tkes_budget_lines = [
            {'var_names': ['ADVTRS'], 'legend_label': 'ADVTRS'},
            {'var_names': ['SHEARS'], 'legend_label': 'SHEARS'},
            {'var_names': ['BUOYAS'], 'legend_label': 'BUOYAS'},
            {'var_names': ['DISSIPS'], 'legend_label': 'DISSIPS', 'sam_conv_factor': -1.},
            {'var_names': ['TKES_RES'], 'legend_label': 'TKES_RES', 'sam_calc': self.getTkesResidual},
        ]

        up2_vp2_budget_lines = [
            {'var_names': ['U2V2ADV'], 'legend_label': 'U2V2ADV', 'sam_calc': self.getU2V2Adv},
            {'var_names': ['U2V2BUOY'], 'legend_label': 'U2V2BUOY', 'sam_calc': self.getU2V2Buoy},
            {'var_names': ['U2V2PRESS'], 'legend_label': 'U2V2PRES', 'sam_calc': self.getU2V2Pres},
            {'var_names': ['W2REDIS'], 'legend_label': 'U2V2REDIS', 'sam_conv_factor': -1},
            {'var_names': ['U2V2DFSN'], 'legend_label': 'U2V2DFSN', 'sam_calc': self.getU2V2Dfsn},
            {'var_names': ['DISSIP'], 'legend_label': 'U2V2DISSIP', 'sam_conv_factor': 2},
            {'var_names': ['U2V2SDMP'], 'legend_label': 'U2V2SDMP', 'sam_calc': self.getU2V2Sdmp},
            {'var_names': ['SHEAR'], 'legend_label': 'U2V2SHEAR', 'sam_conv_factor': 2},
            {'var_names': ['U2V2BT'], 'legend_label': 'U2V2BT', 'sam_calc': self.getU2V2Bt},
            # {'var_names': ['U2V2_RES'], 'legend_label': '2 * TKE_RES - W2_RES'},
            {'var_names': ['U2V2_RES'], 'legend_label': 'U2V2_RES', 'sam_calc': self.getU2V2Residual},
        ]

        upwp_budget_lines = [
            {'var_names': ['WUDFSN'], 'legend_label': 'WUDFSN'},
            {'var_names': ['WU_RES'], 'legend_label': 'WU_RES', 'sam_calc': self.getUWResidual},
            {'var_names': ['WUADV'], 'legend_label': 'WUADV'},
            {'var_names': ['WUPRES'], 'legend_label': 'WUPRES'},
            {'var_names': ['WUANIZ'], 'legend_label': 'WUANIZ'},
            {'var_names': ['WUBUOY'], 'legend_label': 'WUBUOY'},
            {'var_names': ['WUSHEAR'], 'legend_label': 'WUSHEAR'},
            {'var_names': ['WUBT'], 'legend_label': 'WUBT'},
            # {'var_names': ['WUSDMP'], 'legend_label': 'WUSDMP'},
        ]

        vpwp_budget_lines = [
            {'var_names': ['WVDFSN'], 'legend_label': 'WVDFSN'},
            {'var_names': ['WV_RES'], 'legend_label': 'WV_RES', 'sam_calc': self.getVWResidual},
            {'var_names': ['WVADV'], 'legend_label': 'WVADV'},
            {'var_names': ['WVPRES'], 'legend_label': 'WVPRES'},
            {'var_names': ['WVANIZ'], 'legend_label': 'WVANIZ'},
            {'var_names': ['WVBUOY'], 'legend_label': 'WVBUOY'},
            {'var_names': ['WVSHEAR'], 'legend_label': 'WVSHEAR'},
            {'var_names': ['WVBT'], 'legend_label': 'WVBT'},
            # {'var_names': ['WVSDMP'], 'legend_label': 'WVSDMP'},
        ]

        up2_budget_lines = [
            {'var_names': ['U2ADV'], 'legend_label': 'U2ADV'},
            {'var_names': ['U2SHEAR'], 'legend_label': 'U2SHEAR'},
            {'var_names': ['U2REDIS'], 'legend_label': 'U2REDIS'},
            {'var_names': ['U2DFSN'], 'legend_label': 'U2DFSN'},
            {'var_names': ['U2BT'], 'legend_label': 'U2BT'},
            {'var_names': ['U2_RES'], 'legend_label': 'U2_RES', 'sam_calc': self.getU2Residual},
        ]

        vp2_budget_lines = [
            {'var_names': ['V2ADV'], 'legend_label': 'V2ADV'},
            {'var_names': ['V2SHEAR'], 'legend_label': 'V2SHEAR'},
            {'var_names': ['V2REDIS'], 'legend_label': 'V2REDIS'},
            {'var_names': ['V2DFSN'], 'legend_label': 'V2DFSN'},
            {'var_names': ['V2BT'], 'legend_label': 'V2BT'},
            {'var_names': ['V2_RES'], 'legend_label': 'V2_RES', 'sam_calc': self.getV2Residual},
        ]

        wp2_budget_lines = [
            {'var_names': ['W2ADV'], 'legend_label': 'W2ADV'},
            {'var_names': ['W2PRES'], 'legend_label': 'W2PRES'},
            {'var_names': ['W2REDIS'], 'legend_label': 'W2REDIS'},
            {'var_names': ['W2BUOY'], 'legend_label': 'W2BUOY'},
            {'var_names': ['W2DFSN'], 'legend_label': 'W2DFSN'},
            {'var_names': ['W2SDMP'], 'legend_label': 'W2SDMP'},
            {'var_names': ['W2BT'], 'legend_label': 'W2BT'},
            {'var_names': ['W2_RES'], 'legend_label': 'W2_RES', 'sam_calc': self.getW2Residual},
        ]

        wp3_budget_lines = [
            {'var_names': ['W3ADV'], 'legend_label': 'W3ADV'},
            {'var_names': ['W3BUOY'], 'legend_label': 'W3BUOY'},
            {'var_names': ['W3DFSN'], 'legend_label': 'W3DFSN'},
            {'var_names': ['W3BT'], 'legend_label': 'W3BT'},
            {'var_names': ['W3TP'], 'legend_label': 'W3TP'},
            {'var_names': ['W3PRESDFSN'], 'legend_label': 'W3PRDFSN'},
            {'var_names': ['W3PRESSCR'], 'legend_label': 'W3PRSCR'},
            {'var_names': ['W3_RES'], 'legend_label': 'W3_RES', 'sam_calc': self.getW3Residual},
           
        ]

        ## Tau plots showing C_2/tau and C_14/tau ##
        # No support for 2nd axis, not needed
        # up2vp2tau_budget_lines = [
        # {'var_names': [r'U2DFSN'], 'legend_label': 'U2DFSN'},
        # {'var_names': [r'V2DFSN'], 'legend_label': 'V2DFSN'},
        # {'var_names': [r'TKES'], 'legend_label': 'TKES', 1, 1],
        # {'var_names': [r'TKE'], 'legend_label': 'TKE', 1, 1],
        # {'var_names': ['U2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (u'^2)$",
        #                'sam_calc': self.getU2C14OverTau},
        # {'var_names': ['V2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (v'^2)$",
        #                'sam_calc': self.getV2C14OverTau},
        # ]

        # qtogp2tau_budget_lines = [
        # {'var_names': [r'QTOG2DISSIP'], 'legend_label': 'QTOG2DISSIP'},
        # {'var_names': [r'QTOG2DIFTR'], 'legend_label': 'QTOG2DIFTR'},
        # {'var_names': [r'QTOG2'], 'legend_label': 'QTOG2', 1, 1],
        # {'var_names': ['QTOG2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$",
        #  'sam_calc': self.getQtog2C2OverTau},
        # ]

        # qp2tau_budget_lines = [
        # {'var_names': [r'Q2DISSIP'], 'legend_label': 'Q2DISSIP'},
        # {'var_names': [r'Q2DIFTR'], 'legend_label': 'Q2DIFTR'},
        # {'var_names': [r'QT2'], 'legend_label': 'QT2', 1, 1],
        # {'var_names': ['QT2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q'^2)$",
        #  'sam_calc': self.getQ2C2OverTau},
        # ]

        # thlp2tau_budget_lines = [
        # {'var_names': ['THL2DISSIP'], 'legend_label': 'THL2DISSIP'},
        # {'var_names': ['THL2DIFTR'], 'legend_label': 'THL2DIFTR'},
        # {'var_names': ['THEL2'], 'legend_label': 'THEL2', 1, 1],
        # {'var_names': ['THL2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (\theta_l'^2)$",
        #  'sam_calc': self.getThl2C2OverTau},
        # ]

        # tau_comp_budget_lines = [
        # {'var_names': ['U2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (u'^2)$",
        #  'sam_calc': self.getU2C14OverTau},
        # {'var_names': ['V2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (v'^2)$",
        #  'sam_calc': self.getV2C14OverTau},
        # {'var_names': [r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$",
        #  'legend_label': '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2 * 1e-6, 1e-15 )'], 1, 1],#??
        # {'var_names': ['QTOG2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$",
        #  'sam_calc': self.getQtog2C2OverTau},
        # {'var_names': ['THL2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (\theta_l'^2)$",
        #  'sam_calc': self.getThl2C2OverTau},
        # {'var_names': ['QT2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q'^2)$",
        #  'sam_calc': self.getQ2C2OverTau},
        # ]

        ## Reduced second moment flux budgets ##TODO: How to sort lines and define colors?
        up2_reduced_budget_lines = [
            {'var_names': ['U2ADV'], 'legend_label': 'advection'},
            # {'var_names': ['U2BUOY'], 'legend_label': 'buoyancy', 'sam_calc': self.getNothing},
            {'var_names': ['U2DFSN'], 'legend_label': 'dissipation'},
            {'var_names': ['U2REDIS'], 'legend_label': 'isotropy'},
            # {'var_names': ['U2PRES'], 'legend_label': 'pressure', 'sam_calc': self.getNothing},
            {'var_names': ['U2SHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['U2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['U2_RES'], 'legend_label': 'residual', 'sam_calc': self.getU2Residual},
        ]

        vp2_reduced_budget_lines = [
            {'var_names': ['V2ADV'], 'legend_label': 'advection'},
            # {'var_names': ['V2BUOY'], 'legend_label': 'buoyancy', 'sam_calc': self.getNothing},
            {'var_names': ['V2DFSN'], 'legend_label': 'dissipation'},
            {'var_names': ['V2REDIS'], 'legend_label': 'isotropy'},
            # {'var_names': ['V2PRES'], 'legend_label': 'pressure', 'sam_calc': self.getNothing},
            {'var_names': ['V2SHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['V2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['V2_RES'], 'legend_label': 'residual', 'sam_calc': self.getV2Residual},
        ]

        wp2_reduced_budget_lines = [
            # {'var_names': ['W2SDMP'], 'legend_label': 'damping'},
            {'var_names': ['W2ADV'], 'legend_label': 'advection'},
            {'var_names': ['W2BUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['W2DFSN'], 'legend_label': 'dissipation'},
            {'var_names': ['W2REDIS+PRES'], 'legend_label': 'pressure', 'sam_calc': self.getW2RedisPlusPres},
            # {'var_names': ['W2SHEAR'], 'legend_label': 'turb. prod.', 'sam_calc': self.getNothing},
            {'var_names': ['W2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['W2_RES'], 'legend_label': 'residual', 'sam_calc': self.getW2Residual},
        ]

        upwp_reduced_budget_lines = [
            # {'var_names': ['WUSDMP'], 'legend_label': 'WUSDMP'},
            {'var_names': ['WUADV'], 'legend_label': 'advection'},
            {'var_names': ['WUBUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['WUDFSN'], 'legend_label': 'dissipation'},
            {'var_names': ['WUPRES+ANIZ'], 'legend_label': 'pressure', 'sam_calc': self.getUWPresPlusAniz},
            {'var_names': ['WUSHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['WUBT'], 'legend_label': 'time tndcy'},
            {'var_names': ['WU_RES'], 'legend_label': 'residual', 'sam_calc': self.getUWResidual},
        ]

        vpwp_reduced_budget_lines = [
            # {'var_names': ['WVSDMP'], 'legend_label': 'WUSDMP'},
            {'var_names': ['WVADV'], 'legend_label': 'advection'},
            {'var_names': ['WVBUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['WVDFSN'], 'legend_label': 'dissipation'},
            {'var_names': ['WVPRES+ANIZ'], 'legend_label': 'pressure', 'sam_calc': self.getVWPresPlusAniz},
            {'var_names': ['WVSHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['WVBT'], 'legend_label': 'time tndcy'},
            {'var_names': ['WV_RES'], 'legend_label': 'residual', 'sam_calc': self.getVWResidual},
        ]

        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['HL'],
                'sam': ['HL'],
                'coamps': ['HL'],
                'r408': ['HL'],
                'hoc': ['HL'],
                'e3sm': ['HL'],
                'cam': ['HL'],
                },
             'lines': hlp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': 'Liquid-Ice Static Energy Budget, LIMSE',
             'axis_title': r"LIMSE budget terms $\mathrm{\left[K\,s^{-1}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['QT'],
                'sam': ['QT'],
                'coamps': ['QT'],
                'r408': ['QT'],
                'hoc': ['QT'],
                'e3sm': ['QT'],
                'cam': ['QT'],
                },
             'lines': qtp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r'Total Water (No Rain/Snow Included) Budget, $\mathrm{r_t}$',
             'axis_title': r"$\mathrm{r_t}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,s^{-1}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['TW_B+P'],
                'sam': ['TW_B+P'],
                'coamps': ['TW_B+P'],
                'r408': ['TW_B+P'],
                'hoc': ['TW_B+P'],
                'e3sm': ['TW_B+P'],
                'cam': ['TW_B+P'],
                },
             'lines': tpwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Static Energy Flux Budget (Reduced), $\mathrm{\overline{s_v'w'}}$",
             'axis_title': r"$\mathrm{\overline{s_v'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['TW_COMPLETE'],
                'sam': ['TW_COMPLETE'],
                'coamps': ['TW_COMPLETE'],
                'r408': ['TW_COMPLETE'],
                'hoc': ['TW_COMPLETE'],
                'e3sm': ['TW_COMPLETE'],
                'cam': ['TW_COMPLETE'],
                },
             'lines': tpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Static Energy Flux Budget, $\mathrm{\overline{s_v'w'}}$",
             'axis_title': r"$\mathrm{\overline{s_v'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['THLW_B+P'],
                'sam': ['THLW_B+P'],
                'coamps': ['THLW_B+P'],
                'r408': ['THLW_B+P'],
                'hoc': ['THLW_B+P'],
                'e3sm': ['THLW_B+P'],
                'cam': ['THLW_B+P'],
                },
             'lines': thlpwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Potential Temperature Flux Budget (Reduced), "+
                      r"$\mathrm{\overline{\theta_l'w'}}$",
             'axis_title': r"$\mathrm{\overline{\theta_l'w'}}$ budget terms "+
                           r"$\mathrm{\left[m\,K\,s^{-2}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['THLW_COMPLETE'],
                'sam': ['THLW_COMPLETE'],
                'coamps': ['THLW_COMPLETE'],
                'r408': ['THLW_COMPLETE'],
                'hoc': ['THLW_COMPLETE'],
                'e3sm': ['THLW_COMPLETE'],
                'cam': ['THLW_COMPLETE'],
                },
             'lines': thlpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Potential Temperature Flux Budget, $\mathrm{\overline{\theta_l'w'}}$",
             'axis_title': r"$\mathrm{\overline{\theta_l'w'}}$ budget terms "+
                           r"$\mathrm{\left[m\,K\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['QW_B+P'],
                'sam': ['QW_B+P'],
                'coamps': ['QW_B+P'],
                'r408': ['QW_B+P'],
                'hoc': ['QW_B+P'],
                'e3sm': ['QW_B+P'],
                'cam': ['QW_B+P'],
                },
             'lines': qpwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Total Water Flux Budget (Reduced), $\mathrm{\overline{r_t'w'}}$",
             'axis_title': r"$\mathrm{\overline{r_t'w'}}$ budget terms "+
                           r"$\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['QW_COMPLETE'],
                'sam': ['QW_COMPLETE'],
                'coamps': ['QW_COMPLETE'],
                'r408': ['QW_COMPLETE'],
                'hoc': ['QW_COMPLETE'],
                'e3sm': ['QW_COMPLETE'],
                'cam': ['QW_COMPLETE'],
                },
             'lines': qpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Total Water Flux Budget, $\mathrm{\overline{r_t'w'}}$",
             'axis_title': r"$\mathrm{\overline{r_t'w'}}$ budget terms "+
                           r"$\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['QTOGW_B+P'],
                'sam': ['QTOGW_B+P'],
                'coamps': ['QTOGW_B+P'],
                'r408': ['QTOGW_B+P'],
                'hoc': ['QTOGW_B+P'],
                'e3sm': ['QTOGW_B+P'],
                'cam': ['QTOGW_B+P'],
                },
             'lines': qtogpwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Grand Total Water Flux Budget (Reduced), $\mathrm{\overline{q_{tog}'w'}}$",
             'axis_title': r"$\mathrm{\overline{q_{tog}'w'}}$ budget terms "+
                           r"$\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['QTOGW_COMPLETE'],
                'sam': ['QTOGW_COMPLETE'],
                'coamps': ['QTOGW_COMPLETE'],
                'r408': ['QTOGW_COMPLETE'],
                'hoc': ['QTOGW_COMPLETE'],
                'e3sm': ['QTOGW_COMPLETE'],
                'cam': ['QTOGW_COMPLETE'],
                },
             'lines': qtogpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Grand Total Water Flux Budget, $\mathrm{\overline{q_{tog}'w'}}$",
             'axis_title': r"$\mathrm{\overline{q_{tog}'w'}}$ budget terms "+
                           r"$\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['W2'],
                'sam': ['W2'],
                'coamps': ['W2'],
                'r408': ['W2'],
                'hoc': ['W2'],
                'e3sm': ['W2'],
                'cam': ['W2'],
                },
             'lines': wp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Vertical Momentum Variance Budget, $\mathrm{\overline{w'^2}}$",
             'axis_title': r"$\mathrm{\overline{w'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['W3'],
                'sam': ['W3'],
                'coamps': ['W3'],
                'r408': ['W3'],
                'hoc': ['W3'],
                'e3sm': ['W3'],
                'cam': ['W3'],
                },
             'lines': wp3_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Vertical Momentum Skewness Budget, $\mathrm{\overline{w'^3}}$",
             'axis_title': r"$\mathrm{\overline{w'^3}}$ budget terms $\mathrm{\left[m^3\,s^{-4}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['T2'],
                'sam': ['T2'],
                'coamps': ['T2'],
                'r408': ['T2'],
                'hoc': ['T2'],
                'e3sm': ['T2'],
                'cam': ['T2'],
                },
             'lines': tp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Static Energy Variance Budget, $\mathrm{\overline{s_v'^2}}$",
             'axis_title': r"$\mathrm{\overline{s_v'^2}}$ budget terms $\mathrm{\left[K^2\,s^{-1}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['THL2'],
                'sam': ['THL2'],
                'coamps': ['THL2'],
                'r408': ['THL2'],
                'hoc': ['THL2'],
                'e3sm': ['THL2'],
                'cam': ['THL2'],
                },
             'lines': thlp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Liquid Water Potential Temperature Variance Budget, $\mathrm{\overline{\theta_l'^2}}$",
             'axis_title': r"$\mathrm{\overline{\theta_l'^2}}$ budget terms $\mathrm{\left[K^2\,s^{-1}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['Q2'],
                'sam': ['Q2'],
                'coamps': ['Q2'],
                'r408': ['Q2'],
                'hoc': ['Q2'],
                'e3sm': ['Q2'],
                'cam': ['Q2'],
                },
             'lines': rtp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Total Water Variance Budget, $\mathrm{\overline{r_t'^2}}$",
             'axis_title': r"$\mathrm{\overline{r_t'^2}}$ budget terms "+
                           r"$\mathrm{\left[kg^2\,kg^{-2}\,s^{-1}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['QTOG2'],
                'sam': ['QTOG2'],
                'coamps': ['QTOG2'],
                'r408': ['QTOG2'],
                'hoc': ['QTOG2'],
                'e3sm': ['QTOG2'],
                'cam': ['QTOG2'],
                },
             'lines': qtogp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Grand Total Water Variance Budget, $\mathrm{\overline{q_{tog}'^2}}$",
             'axis_title': r"$\mathrm{\overline{q_{tog}'^2}}$ budget terms "+
                           r"$\mathrm{\left[kg^2\,kg^{-2}\,s^{-1}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['QTHL'],
                'sam': ['QTHL'],
                'coamps': ['QTHL'],
                'r408': ['QTHL'],
                'hoc': ['QTHL'],
                'e3sm': ['QTHL'],
                'cam': ['QTHL'],
                },
             'lines': qpthlp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Total Water, Liquid Water Pot. Temp. Covariance Budget, "+
                      r"$\mathrm{\overline{q_t'\theta_l'}}$",
             'axis_title': r"$\mathrm{\overline{q_t'\theta_l'}}$ budget terms "+
                           r"$\mathrm{\left[kg\,kg^{-1}\,K\,s^{-1}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['TKE'],
                'sam': ['TKE'],
                'coamps': ['TKE'],
                'r408': ['TKE'],
                'hoc': ['TKE'],
                'e3sm': ['TKE'],
                'cam': ['TKE'],
                },
             'lines': tke_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Turbulence Kinetic Energy Budget, TKE",
             'axis_title': r"TKE budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$", 'centered': True,
            },
            {'var_names':
                {
                'clubb': ['TKES'],
                'sam': ['TKES'],
                'coamps': ['TKES'],
                'r408': ['TKES'],
                'hoc': ['TKES'],
                'e3sm': ['TKES'],
                'cam': ['TKES'],
                },
             'lines': tkes_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Sub-Grid Turbulence Kinetic Energy Budget, $\mathrm{TKE_{SGS}}$",
             'axis_title': r"$\mathrm{TKE_{SGS}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['U2V2'],
                'sam': ['U2V2'],
                'coamps': ['U2V2'],
                'r408': ['U2V2'],
                'hoc': ['U2V2'],
                'e3sm': ['U2V2'],
                'cam': ['U2V2'],
                },
             'lines': up2_vp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Horizontal Turbulence Kinetic Energy Budget, $\mathrm{\overline{u'^2}+\overline{v'^2}}$",
             'axis_title': r"$\mathrm{\overline{u'^2}+\overline{v'^2}}$ budget terms "+
                           r"$\mathrm{\left[m^2\, s^{-3}\right]}$",
             'centered': True,
            },
            # TODO: Rename to UW?
            {'var_names':
                {
                'clubb': ['WU'],
                'sam': ['WU'],
                'coamps': ['WU'],
                'r408': ['WU'],
                'hoc': ['WU'],
                'e3sm': ['WU'],
                'cam': ['WU'],
                },
             'lines': upwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Eastward Momentum Flux Budget, $\mathrm{\overline{u'w'}}$",
             'axis_title': r"$\mathrm{\overline{u'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['WV'],
                'sam': ['WV'],
                'coamps': ['WV'],
                'r408': ['WV'],
                'hoc': ['WV'],
                'e3sm': ['WV'],
                'cam': ['WV'],
                },
             'lines': vpwp_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Northward Momentum Flux Budget, $\mathrm{\overline{v'w'}}$",
             'axis_title': r"$\mathrm{\overline{v'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['U2'],
                'sam': ['U2'],
                'coamps': ['U2'],
                'r408': ['U2'],
                'hoc': ['U2'],
                'e3sm': ['U2'],
                'cam': ['U2'],
                },
             'lines': up2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Eastward Momentum Variance Budget, $\mathrm{\overline{u'^2}}$",
             'axis_title': r"$\mathrm{\overline{u'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             'priority': False,
            },
            {'var_names':
                {
                'clubb': ['V2'],
                'sam': ['V2'],
                'coamps': ['V2'],
                'r408': ['V2'],
                'hoc': ['V2'],
                'e3sm': ['V2'],
                'cam': ['V2'],
                },
             'lines': vp2_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Northward Momentum Variance Budget, $\mathrm{\overline{v'^2}}$",
             'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             'priority': False,
            },
            #{'var_names':
            #    {
            #    'clubb': ['V2_COMP'],
            #    'sam': ['V2_COMP'],
            #    'coamps': ['V2_COMP'],
            #    'r408': ['V2_COMP'],
            #    'hoc': ['V2_COMP'],
            #    'e3sm': ['V2_COMP'],
            # 	 'cam': ['V2_COMP'],
            #    },
            # 'lines': vp2_comp_budget_lines, 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\ s^{-3}\right]}$",
            #},
            #{'var_names':
            #    {
            #    'clubb': ['U2_V2_TAU'],
            #    'sam': ['U2_V2_TAU'],
            #    'coamps': ['U2_V2_TAU'],
            #    'r408': ['U2_V2_TAU'],
            #    'hoc': ['U2_V2_TAU'],
            #    'e3sm': ['U2_V2_TAU'],
            #    'cam': ['U2_V2_TAU'],
            #    },
            # 'lines': up2vp2tau_budget_lines, 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\frac{C_{14}}{\tau}}\ \mathrm{\left[s^{-1}\right]}$",
            #},
            #{'var_names':
            #    {
            #    'clubb': ['QTOG2_TAU'],
            #    'sam': ['QTOG2_TAU'],
            #    'coamps': ['QTOG2_TAU'],
            #    'r408': ['QTOG2_TAU'],
            #    'hoc': ['QTOG2_TAU'],
            #    'e3sm': ['QTOG2_TAU'],
            # 	 'cam': ['QTOG2_TAU'],
            #    },
            # 'lines': qtogp2tau_budget_lines, 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$",
            #},
            #{'var_names':
            #    {
            #    'clubb': ['Q2_TAU'],
            #    'sam': ['Q2_TAU'],
            #    'coamps': ['Q2_TAU'],
            #    'r408': ['Q2_TAU'],
            #    'hoc': ['Q2_TAU'],
            #    'e3sm': ['Q2_TAU'],
            #    'cam': ['Q2_TAU'],
            #    },
            # 'lines': qp2tau_budget_lines, 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$",
            #},
            #{'var_names':
            #    {
            #    'clubb': ['THL2_TAU'],
            #    'sam': ['THL2_TAU'],
            #    'coamps': ['THL2_TAU'],
            #    'r408': ['THL2_TAU'],
            #    'hoc': ['THL2_TAU'],
            #    'e3sm': ['THL2_TAU'],
            #    'cam': ['THL2_TAU'],
            #    },
            # 'lines': qp2tau_budget_lines,
            # 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$",
            #},
            #{'var_names':
            #    {
            #    'clubb': ['C2_TAU_COMP'],
            #    'sam': ['C2_TAU_COMP'],
            #    'coamps': ['C2_TAU_COMP'],
            #    'r408': ['C2_TAU_COMP'],
            #    'hoc': ['C2_TAU_COMP'],
            #    'e3sm': ['C2_TAU_COMP'],
            #    'cam': ['C2_TAU_COMP'],
            #    },
            # 'lines': tau_comp_budget_lines,
            # 'type': Panel.TYPE_BUDGET,
            # 'axis_title': r"$\mathrm{\frac{C}{\tau}}\ \mathrm{\left[s^{-1}\right]}$",
            #},
            {'var_names':
                {
                'clubb': ['U2 REDUCED'],
                'sam': ['U2 REDUCED'],
                'coamps': ['U2 REDUCED'],
                'r408': ['U2 REDUCED'],
                'hoc': ['U2 REDUCED'],
                'e3sm': ['U2 REDUCED'],
                'cam': ['U2 REDUCED'],
                },
             'lines': up2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Eastward Momentum Variance Budget (Reduced), $\mathrm{\overline{u'^2}}$",
             'axis_title': r"$\mathrm{\overline{u'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['V2 REDUCED'],
                'sam': ['V2 REDUCED'],
                'coamps': ['V2 REDUCED'],
                'r408': ['V2 REDUCED'],
                'hoc': ['V2 REDUCED'],
                'e3sm': ['V2 REDUCED'],
                'cam': ['V2 REDUCED'],
                },
             'lines': vp2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Northward Momentum Variance Budget (Reduced), $\mathrm{\overline{v'^2}}$",
             'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['W2 REDUCED'],
                'sam': ['W2 REDUCED'],
                'coamps': ['W2 REDUCED'],
                'r408': ['W2 REDUCED'],
                'hoc': ['W2 REDUCED'],
                'e3sm': ['W2 REDUCED'],
                'cam': ['W2 REDUCED'],
                },
             'lines': wp2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Vertical Momentum Variance Budget (Reduced), $\mathrm{\overline{w'^2}}$",
             'axis_title': r"$\mathrm{\overline{w'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['WU REDUCED'],
                'sam': ['WU REDUCED'],
                'coamps': ['WU REDUCED'],
                'r408': ['WU REDUCED'],
                'hoc': ['WU REDUCED'],
                'e3sm': ['WU REDUCED'],
                'cam': ['WU REDUCED'],
                },
             'lines': upwp_reduced_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Eastward Momentum Flux Budget (Reduced), $\mathrm{\overline{u'w'}}$",
             'axis_title': r"$\mathrm{\overline{u'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
            },
            {'var_names':
                {
                'clubb': ['WV REDUCED'],
                'sam': ['WV REDUCED'],
                'coamps': ['WV REDUCED'],
                'r408': ['WV REDUCED'],
                'hoc': ['WV REDUCED'],
                'e3sm': ['WV REDUCED'],
                'cam': ['WV REDUCED'],
                },
             'lines': vpwp_reduced_budget_lines, 'type': Panel.TYPE_BUDGET,
             'title': r"Northward Momentum Flux Budget (Reduced), $\mathrm{\overline{v'w'}}$",
             'axis_title': r"$\mathrm{\overline{v'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$",
             'centered': True,
             },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)

    def getHlResidual(self, dataset_override=None):
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

        Calculates the residual for the HL budget using
        the following equation:
        ``(HLSTOR)+((-1)*(HLADV+HLDFSN+HLRAD+HLLAT+TTEND))*g_per_second_to_kg_per_day``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        HLADV, indep, dataset = self.getVarForCalculations('HLADV', dataset,
                                                           conversion_factor = self.kg_per_second_to_kg_per_day)

        HLDFSN, indep, dataset = self.getVarForCalculations('HLDFSN', dataset,
                                                        conversion_factor = self.kg_per_second_to_kg_per_day)

        HLLAT, indep, dataset = self.getVarForCalculations('HLLAT', dataset,
                                                       conversion_factor = self.kg_per_second_to_kg_per_day)

        HLRAD, indep, dataset = self.getVarForCalculations('HLRAD', dataset,
                                                       conversion_factor = self.kg_per_second_to_kg_per_day)

        HLSTOR, indep, dataset = self.getVarForCalculations('HLSTOR', dataset,
                                                        conversion_factor = self.kg_per_second_to_kg_per_day)

        TTEND, indep, dataset = self.getVarForCalculations('TTEND', dataset,
                                                       conversion_factor = self.kg_per_second_to_kg_per_day)

        # ALERT: Conversion factor correct?
        HL_RES = (HLSTOR - (HLADV + HLDFSN + HLLAT + HLRAD + TTEND)) * self.g_per_second_to_kg_per_day
        return HL_RES, indep

    def getQtResidual(self, dataset_override=None):
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

        Caclulates Nim from sam -> clubb using the equation

        .. code-block:: python
            :linenos:

            QTSTOR+(-1)*(QTADV+QTDFSN+QTSRC+QTSINK+QTEND)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTADV, indep, dataset = self.getVarForCalculations('QTADV', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QTDFSN, indep, dataset = self.getVarForCalculations('QTDFSN', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QTEND, indep, dataset = self.getVarForCalculations('QTEND', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QTSINK, indep, dataset = self.getVarForCalculations('QTSINK', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QTSRC, indep, dataset = self.getVarForCalculations('QTSRC', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QTSTOR, indep, dataset = self.getVarForCalculations('QTSTOR', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QV_TNDCY, indep, dataset = self.getVarForCalculations('QV_TNDCY', dataset,
                                                       conversion_factor  = self.g_per_second_to_kg_per_day)
        QT_RES = QTSTOR - (QTADV + QTDFSN + QTEND + QTSRC + QTSINK)
        return QT_RES, indep

    def getTwBuoyPlusPres(self, dataset_override=None):
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

        Calculates the combined BUOY and PRES variable value
        of the TW budget using the following equation:
        ``TWBUOY+TWPRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TWBUOY, indep, dataset = self.getVarForCalculations('TWBUOY', dataset)
        TWPRES, indep, dataset = self.getVarForCalculations('TWPRES', dataset)
        TW_BUOY_PRES = TWBUOY + TWPRES
        return TW_BUOY_PRES, indep

    def getTwResidual(self, dataset_override=None):
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

        Calculates the residual for the TW budget using
        the following equation:
        ``TWBT - (TWGRAD + TWADV + TWDFSN + TWBUOY + TWPRES + TWPREC + TWRAD + TWFORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TWADV, indep, dataset = self.getVarForCalculations('TWADV', dataset)
        TWBT, indep, dataset = self.getVarForCalculations('TWBT', dataset)
        TWBUOY, indep, dataset = self.getVarForCalculations('TWBUOY', dataset)
        TWDFSN, indep, dataset = self.getVarForCalculations('TWDFSN', dataset)
        TWFORC, indep, dataset = self.getVarForCalculations('TWFORC', dataset)
        TWGRAD, indep, dataset = self.getVarForCalculations('TWGRAD', dataset)
        TWPREC, indep, dataset = self.getVarForCalculations('TWPREC', dataset)
        TWPRES, indep, dataset = self.getVarForCalculations('TWPRES', dataset)
        TWRAD, indep, dataset = self.getVarForCalculations('TWRAD', dataset)
        TW_RES = TWBT - (TWADV + TWBUOY + TWDFSN + TWFORC + TWGRAD + TWPREC + TWPRES + TWRAD)
        return TW_RES, indep

    def getThlwBuoyPlusPres(self, dataset_override=None):
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

        Calculates the combined BUOY and PRES variable value
        of the THLW budget using the following equation:
        ``THLWBUOY + THLWPRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        THLWBUOY, indep, dataset = self.getVarForCalculations('THLWBUOY', dataset)
        THLWPRES, indep, dataset = self.getVarForCalculations('THLWPRES', dataset)
        THLW_BUOY_PRES = THLWBUOY + THLWPRES
        return THLW_BUOY_PRES, indep

    def getThlwResidual(self, dataset_override=None):
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

        Calculates the residual for the THLW budget using
        the following equation:
        ``THLWBT - (THLWGRAD + THLWADV + THLWDFSN + THLWBUOY + THLWPRES + THLWPREC + THLWRAD + THLWFORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        THLWADV, indep, dataset = self.getVarForCalculations('THLWADV', dataset)
        THLWBT, indep, dataset = self.getVarForCalculations('THLWBT', dataset)
        THLWBUOY, indep, dataset = self.getVarForCalculations('THLWBUOY', dataset)
        THLWDFSN, indep, dataset = self.getVarForCalculations('THLWDFSN', dataset)
        THLWFORC, indep, dataset = self.getVarForCalculations('THLWFORC', dataset)
        THLWGRAD, indep, dataset = self.getVarForCalculations('THLWGRAD', dataset)
        THLWPREC, indep, dataset = self.getVarForCalculations('THLWPREC', dataset)
        THLWPRES, indep, dataset = self.getVarForCalculations('THLWPRES', dataset)
        THLWRAD, indep, dataset = self.getVarForCalculations('THLWRAD', dataset)
        THLW_RES = THLWBT - (THLWADV + THLWBUOY + THLWDFSN + THLWFORC + THLWGRAD + THLWPREC + THLWPRES + THLWRAD)
        return THLW_RES, indep

    def getQwBuoyPlusPres(self, dataset_override=None):
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

        Calculates the combined BUOY and PRES variable value
        of the QW budget using the following equation:
        ``QWBUOY+QWPRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QWBUOY, indep, dataset = self.getVarForCalculations('QWBUOY', dataset)
        QWPRES, indep, dataset = self.getVarForCalculations('QWPRES', dataset)
        QW_BUOY_PRES = QWBUOY + QWPRES
        return QW_BUOY_PRES, indep

    def getQwResidual(self, dataset_override=None):
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

        Calculates the residual for the QW budget using
        the following equation:
        ``QWBT - (QWGRAD + QWADV + QWDFSN + QWBUOY + QWPRES + QWPREC + QWFORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QWADV, indep, dataset = self.getVarForCalculations('QWADV', dataset)
        QWBT, indep, dataset = self.getVarForCalculations('QWBT', dataset)
        QWBUOY, indep, dataset = self.getVarForCalculations('QWBUOY', dataset)
        QWDFSN, indep, dataset = self.getVarForCalculations('QWDFSN', dataset)
        QWFORC, indep, dataset = self.getVarForCalculations('QWFORC', dataset)
        QWGRAD, indep, dataset = self.getVarForCalculations('QWGRAD', dataset)
        QWPREC, indep, dataset = self.getVarForCalculations('QWPREC', dataset)
        QWPRES, indep, dataset = self.getVarForCalculations('QWPRES', dataset)
        QW_RES = QWBT - (QWGRAD + QWADV + QWDFSN + QWBUOY + QWPRES + QWPREC + QWFORC)
        return QW_RES, indep

    def getQtogwBuoyPlusPres(self, dataset_override=None):
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

        Calculates the combined BUOY and PRES variable value
        of the QTOGW budget using the following equation:
        ``QTOGWBUOY + QTOGWPRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTOGWBUOY, indep, dataset = self.getVarForCalculations('QTOGWBUOY', dataset)
        QTOGWPRES, indep, dataset = self.getVarForCalculations('QTOGWPRES', dataset)
        QTOGW_BUOY_PRES = QTOGWBUOY + QTOGWPRES
        return QTOGW_BUOY_PRES, indep

    def getQtogwResidual(self, dataset_override=None):
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

        Calculates the residual for the QTOGW budget using
        the following equation:
        ``QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDFSN + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTOGWADV, indep, dataset = self.getVarForCalculations('QTOGWADV', dataset)
        QTOGWBT, indep, dataset = self.getVarForCalculations('QTOGWBT', dataset)
        QTOGWBUOY, indep, dataset = self.getVarForCalculations('QTOGWBUOY', dataset)
        QTOGWDFSN, indep, dataset = self.getVarForCalculations('QTOGWDFSN', dataset)
        QTOGWFORC, indep, dataset = self.getVarForCalculations('QTOGWFORC', dataset)
        QTOGWGRAD, indep, dataset = self.getVarForCalculations('QTOGWGRAD', dataset)
        QTOGWPREC, indep, dataset = self.getVarForCalculations('QTOGWPREC', dataset)
        QTOGWPRES, indep, dataset = self.getVarForCalculations('QTOGWPRES', dataset)
        QTOGW_RES = QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDFSN + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)
        return QTOGW_RES, indep

    def getT2Residual(self, dataset_override=None):
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

        Calculates the residual for the T2 budget using
        the following equation:
        ``T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        T2ADVTR, indep, dataset = self.getVarForCalculations('T2ADVTR', dataset)
        T2BT, indep, dataset = self.getVarForCalculations('T2BT', dataset)
        T2DISSIP, indep, dataset = self.getVarForCalculations('T2DISSIP', dataset)
        T2DIFTR, indep, dataset = self.getVarForCalculations('T2DIFTR', dataset)
        T2FORC, indep, dataset = self.getVarForCalculations('T2FORC', dataset)
        T2GRAD, indep, dataset = self.getVarForCalculations('T2GRAD', dataset)
        T2PREC, indep, dataset = self.getVarForCalculations('T2PREC', dataset)
        T2RAD, indep, dataset = self.getVarForCalculations('T2RAD', dataset)
        T2_RES = T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)
        return T2_RES, indep

    def getThl2Residual(self, dataset_override=None):
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

        Calculates the residual for the THL2 budget using
        the following equation:
        ``THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        THL2ADVTR, indep, dataset = self.getVarForCalculations('THL2ADVTR', dataset)
        THL2BT, indep, dataset = self.getVarForCalculations('THL2BT', dataset)
        THL2DISSIP, indep, dataset = self.getVarForCalculations('THL2DISSIP', dataset)
        THL2DIFTR, indep, dataset = self.getVarForCalculations('THL2DIFTR', dataset)
        THL2FORC, indep, dataset = self.getVarForCalculations('THL2FORC', dataset)
        THL2GRAD, indep, dataset = self.getVarForCalculations('THL2GRAD', dataset)
        THL2PREC, indep, dataset = self.getVarForCalculations('THL2PREC', dataset)
        THL2RAD, indep, dataset = self.getVarForCalculations('THL2RAD', dataset)
        THL2_RES = THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)
        return THL2_RES, indep

    def getQt2Residual(self, dataset_override=None):
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

        Calculates the residual for the Q2 budget using
        the following equation:
        ``Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        Q2ADVTR, indep, dataset = self.getVarForCalculations('Q2ADVTR', dataset)
        Q2BT, indep, dataset = self.getVarForCalculations('Q2BT', dataset)
        Q2DISSIP, indep, dataset = self.getVarForCalculations('Q2DISSIP', dataset)
        Q2DIFTR, indep, dataset = self.getVarForCalculations('Q2DIFTR', dataset)
        Q2FORC, indep, dataset = self.getVarForCalculations('Q2FORC', dataset)
        Q2GRAD, indep, dataset = self.getVarForCalculations('Q2GRAD', dataset)
        Q2PREC, indep, dataset = self.getVarForCalculations('Q2PREC', dataset)
        Q2_RES = Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)
        return Q2_RES, indep

    def getQtog2Residual(self, dataset_override=None):
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

        Calculates the residual for the QTOG2 budget using
        the following equation:
        ``QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTOG2ADVTR, indep, dataset = self.getVarForCalculations('QTOG2ADVTR', dataset)
        QTOG2BT, indep, dataset = self.getVarForCalculations('QTOG2BT', dataset)
        QTOG2DIFTR, indep, dataset = self.getVarForCalculations('QTOG2DIFTR', dataset)
        QTOG2DISSIP, indep, dataset = self.getVarForCalculations('QTOG2DISSIP', dataset)
        QTOG2FORC, indep, dataset = self.getVarForCalculations('QTOG2FORC', dataset)
        QTOG2GRAD, indep, dataset = self.getVarForCalculations('QTOG2GRAD', dataset)
        QTOG2PREC, indep, dataset = self.getVarForCalculations('QTOG2PREC', dataset)
        QTOG2_RES = QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)
        return QTOG2_RES, indep

    def getQThlResidual(self, dataset_override=None):
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

        Calculates the residual for the QTHL budget using
        the following equation:
        ``QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTHLADV, indep, dataset = self.getVarForCalculations('QTHLADV', dataset)
        QTHLBT, indep, dataset = self.getVarForCalculations('QTHLBT', dataset)
        QTHLDIFTR, indep, dataset = self.getVarForCalculations('QTHLDIFTR', dataset)
        QTHLDISSIP, indep, dataset = self.getVarForCalculations('QTHLDISSIP', dataset)
        QTHLFORC, indep, dataset = self.getVarForCalculations('QTHLFORC', dataset)
        QTHLGRAD, indep, dataset = self.getVarForCalculations('QTHLGRAD', dataset)
        QTHLPREC, indep, dataset = self.getVarForCalculations('QTHLPREC', dataset)
        QTHLRAD, indep, dataset = self.getVarForCalculations('QTHLRAD', dataset)
        QTHLW_RES = QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)
        return QTHLW_RES, indep

    def getTkeDissPlusDfsn(self, dataset_override=None):
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

        Calculates the combined DISS and DIFTR variable value
        of the TKE budget using the following equation:
        ``DISSIP + DIFTR``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        DIFTR, indep, dataset = self.getVarForCalculations('DIFTR', dataset)
        DISSIP, indep, dataset = self.getVarForCalculations('DISSIP', dataset)
        TKE_DISS_DFSN = DIFTR + DISSIP
        return TKE_DISS_DFSN, indep

    def getTkeResidual(self, dataset_override=None):
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

        Calculates the residual for the TKE budget using
        the following equation:
        ``BT - (SHEAR + BUOYA + ADVTR + PRESSTR + DIFTR + SDMP + DISSIP)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        ADVTR, indep, dataset = self.getVarForCalculations('ADVTR', dataset)
        BT, indep, dataset = self.getVarForCalculations('BT', dataset)
        BUOYA, indep, dataset = self.getVarForCalculations('BUOYA', dataset)
        DIFTR, indep, dataset = self.getVarForCalculations('DIFTR', dataset)
        DISSIP, indep, dataset = self.getVarForCalculations('DISSIP', dataset)
        PRESSTR, indep, dataset = self.getVarForCalculations('PRESSTR', dataset)
        SDMP, indep, dataset = self.getVarForCalculations('SDMP', dataset)
        SHEAR, indep, dataset = self.getVarForCalculations('SHEAR', dataset)
        TKE_RES = BT - (SHEAR + BUOYA + ADVTR + PRESSTR + DIFTR + SDMP + DISSIP)
        return TKE_RES, indep

    def getTkesResidual(self, dataset_override=None):
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

        Calculates the residual for the TKES budget using
        the following equation:
        ``-(SHEARS + BUOYAS + ADVTRS + DISSIPS)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        ADVTRS, indep, dataset = self.getVarForCalculations('ADVTRS', dataset)
        BUOYAS, indep, dataset = self.getVarForCalculations('BUOYAS', dataset)
        DISSIPS, indep, dataset = self.getVarForCalculations('DISSIPS', dataset)
        SHEARS, indep, dataset = self.getVarForCalculations('SHEARS', dataset)
        TKES_RES = -(SHEARS + BUOYAS + ADVTRS + DISSIPS)
        return TKES_RES, indep

    def getU2Residual(self, dataset_override=None):
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

        Calculates the residual for the U2 budget using
        the following equation:
        ``U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DFSN)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        U2ADV, indep, dataset = self.getVarForCalculations('U2ADV', dataset)
        U2BT, indep, dataset = self.getVarForCalculations('U2BT', dataset)
        U2DFSN, indep, dataset = self.getVarForCalculations('U2DFSN', dataset)
        U2REDIS, indep, dataset = self.getVarForCalculations('U2REDIS', dataset)
        U2SHEAR, indep, dataset = self.getVarForCalculations('U2SHEAR', dataset)
        U2_RES = U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DFSN)
        return U2_RES, indep

    def getV2Residual(self, dataset_override=None):
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

        Calculates the residual for the V2 budget using
        the following equation:
        ``V2BT - (V2ADV + V2SHEAR + V2REDIS + V2DFSN)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        V2ADV, indep, dataset = self.getVarForCalculations('V2ADV', dataset)
        V2BT, indep, dataset = self.getVarForCalculations('V2BT', dataset)
        V2DFSN, indep, dataset = self.getVarForCalculations('V2DFSN', dataset)
        V2REDIS, indep, dataset = self.getVarForCalculations('V2REDIS', dataset)
        V2SHEAR, indep, dataset = self.getVarForCalculations('V2SHEAR', dataset)
        V2_RES = V2BT - (V2ADV + V2SHEAR + V2REDIS + V2DFSN)
        return V2_RES, indep

    def getW2RedisPlusPres(self, dataset_override=None):
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

        Calculates the combined REDIS and PRES variable value
        of the W2 budget using the following equation:
        ``W2REDIS + W2PRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W2PRES, indep, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, indep, dataset = self.getVarForCalculations('W2REDIS', dataset)
        W2_REDIS_PRES = W2REDIS + W2PRES
        return W2_REDIS_PRES, indep

    def getW2Residual(self, dataset_override=None):
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

        Calculates the residual for the W2 budget using
        the following equation:
        ``W2BT - (W2ADV + W2PRES + W2REDIS + W2BUOY + W2DFSN + W2SDMP)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W2ADV, indep, dataset = self.getVarForCalculations('W2ADV', dataset)
        W2BT, indep, dataset = self.getVarForCalculations('W2BT', dataset)
        W2BUOY, indep, dataset = self.getVarForCalculations('W2BUOY', dataset)
        W2DFSN, indep, dataset = self.getVarForCalculations('W2DFSN', dataset)
        W2PRES, indep, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, indep, dataset = self.getVarForCalculations('W2REDIS', dataset)
        W2SDMP, indep, dataset = self.getVarForCalculations('W2SDMP', dataset)
        W2_RES = W2BT - (W2ADV + W2PRES + W2REDIS + W2BUOY + W2DFSN + W2SDMP)
        return W2_RES, indep

    def getU2V2Adv(self, dataset_override):
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

        Calculates the ADV variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * ADVTR - W2ADV``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        ADVTR, indep, dataset = self.getVarForCalculations('ADVTR', dataset)
        W2ADV, indep, dataset = self.getVarForCalculations('W2ADV', dataset)
        U2V2_ADV = 2 * ADVTR + W2ADV
        return U2V2_ADV, indep

    def getU2V2Buoy(self, dataset_override):
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

        Calculates the BUOY variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * BUOYA - W2BUOY``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        BUOYA, indep, dataset = self.getVarForCalculations('BUOYA', dataset)
        W2BUOY, indep, dataset = self.getVarForCalculations('W2BUOY', dataset)
        U2V2_BUOY = 2 * BUOYA + W2BUOY
        return U2V2_BUOY, indep

    def getU2V2Pres(self, dataset_override):
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

        Calculates the PRES variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * PRESSTR - W2PRES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        PRESSTR, indep, dataset = self.getVarForCalculations('PRESSTR', dataset)
        W2PRES, indep, dataset = self.getVarForCalculations('W2PRES', dataset)
        U2V2_PRES = 2 * PRESSTR + W2PRES
        return U2V2_PRES, indep

    def getU2V2Dfsn(self, dataset_override):
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

        Calculates the DFSN variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * DIFTR - W2DFSN``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        DIFTR, indep, dataset = self.getVarForCalculations('DIFTR', dataset)
        W2DFSN, indep, dataset = self.getVarForCalculations('W2DFSN', dataset)
        U2V2_DFSN = 2 * DIFTR + W2DFSN
        return U2V2_DFSN, indep

    def getU2V2Sdmp(self, dataset_override):
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

        Calculates the SDMP variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * SDMP - W2SDMP``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        SDMP, indep, dataset = self.getVarForCalculations('SDMP', dataset)
        W2SDMP, indep, dataset = self.getVarForCalculations('W2SDMP', dataset)
        U2V2_SDMP = 2 * SDMP + W2SDMP
        return U2V2_SDMP, indep

    def getU2V2Bt(self, dataset_override):
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

        Calculates the BT variable value
        of the combined U2+V2 budget using the following equation:
        ``2 * BT - W2BT``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        BT, indep, dataset = self.getVarForCalculations('BT', dataset)
        W2BT, indep, dataset = self.getVarForCalculations('W2BT', dataset)
        U2V2_BT = 2 * BT + W2BT
        return U2V2_BT, indep

    def getU2V2Residual(self, dataset_override):
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

        Calculates the residual for the W2 budget using
        the following equation:
        ``2 * TKE_RES - W2_RES``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        ADVTR, indep, dataset = self.getVarForCalculations('ADVTR', dataset)
        W2ADV, indep, dataset = self.getVarForCalculations('W2ADV', dataset)
        BT, indep, dataset = self.getVarForCalculations('BT', dataset)
        W2BT, indep, dataset = self.getVarForCalculations('W2BT', dataset)
        BUOYA, indep, dataset = self.getVarForCalculations('BUOYA', dataset)
        W2BUOY, indep, dataset = self.getVarForCalculations('W2BUOY', dataset)
        DIFTR, indep, dataset = self.getVarForCalculations('DIFTR', dataset)
        W2DFSN, indep, dataset = self.getVarForCalculations('W2DFSN', dataset)
        DISSIP, indep, dataset = self.getVarForCalculations('DISSIP', dataset)
        PRESSTR, indep, dataset = self.getVarForCalculations('PRESSTR', dataset)
        W2PRES, indep, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, indep, dataset = self.getVarForCalculations('W2REDIS', dataset)
        SDMP, indep, dataset = self.getVarForCalculations('SDMP', dataset)
        W2SDMP, indep, dataset = self.getVarForCalculations('W2SDMP', dataset)
        SHEAR, indep, dataset = self.getVarForCalculations('SHEAR', dataset)
        U2V2_RES = 2. * (BT - (ADVTR + BUOYA + PRESSTR + DIFTR + DISSIP + SDMP + SHEAR)) - W2BT + (
                W2ADV + W2BUOY + W2PRES + W2DFSN + W2SDMP + W2REDIS)
        return U2V2_RES, indep

    def getW3PRESS(self, dataset_override=None):
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

        Calculates the correct value of W3PRES if an older version of SAM with he variable W3REDIS was used.
        In this case, the total effect of pressure on wp3 is W3PRES+W3REDIS
        ``W3PRES = W3PRES + W3REDIS``
     
        If W3REDIS is not present in the output, it will populate as zeros and this function will make no difference.

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W3PRES, indep, dataset = self.getVarForCalculations('W3PRES', dataset)
        W3REDIS, indep, dataset = self.getVarForCalculations('W3REDIS', dataset)
        W3PRESS = W3PRES + W3REDIS
        return W3PRESS, indep

    def getW3Residual(self, dataset_override=None):
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

        Calculates the residual for the W3 budget using
        the following equation:
        ``W3BT - (W3ADV + W3PRES + W3REDIS + W3BUOY + W3DFSN)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.

        This code does not extract W3PRESDFSN, W3PRESSCR, and W3TP directly, but these three variables sum to W3PRES which is extracted here.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W3ADV, indep, dataset = self.getVarForCalculations('W3ADV', dataset)
        W3BT, indep, dataset = self.getVarForCalculations('W3BT', dataset)
        W3BUOY, indep, dataset = self.getVarForCalculations('W3BUOY', dataset)
        W3DFSN, indep, dataset = self.getVarForCalculations('W3DFSN', dataset)
        W3PRES, indep, dataset = self.getVarForCalculations('W3PRES', dataset)
        W3_RES = W3BT - (W3ADV + W3PRES + W3BUOY + W3DFSN)
        return W3_RES, indep

    def getUWPresPlusAniz(self, dataset_override=None):
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

        Calculates the combined ANIZ and PRES variable value
        of the UW budget using the following equation:
        ``WUPRES + WUANIZ``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WUANIZ, indep, dataset = self.getVarForCalculations('WUANIZ', dataset)
        WUPRES, indep, dataset = self.getVarForCalculations('WUPRES', dataset)
        WU_ANIZ_PRES = WUANIZ + WUPRES
        return WU_ANIZ_PRES, indep

    def getUWResidual(self, dataset_override=None):
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

        Calculates the residual for the UW budget using
        the following equation:
        ``WUBT - (WUDFSN + WUSHEAR + WUADV + WUPRES + WUANIZ + WUBUOY + WUSDMP)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WUADV, indep, dataset = self.getVarForCalculations('WUADV', dataset)
        WUANIZ, indep, dataset = self.getVarForCalculations('WUANIZ', dataset)
        WUBT, indep, dataset = self.getVarForCalculations('WUBT', dataset)
        WUBUOY, indep, dataset = self.getVarForCalculations('WUBUOY', dataset)
        WUDFSN, indep, dataset = self.getVarForCalculations('WUDFSN', dataset)
        WUPRES, indep, dataset = self.getVarForCalculations('WUPRES', dataset)
        WUSHEAR, indep, dataset = self.getVarForCalculations('WUSHEAR', dataset)
        WUSDMP, indep, dataset = self.getVarForCalculations('WUSDMP', dataset)
        WU_RES = WUBT - (WUDFSN + WUSHEAR + WUADV + WUPRES + WUANIZ + WUBUOY + WUSDMP)
        return WU_RES, indep

    def getVWPresPlusAniz(self, dataset_override=None):
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

        Calculates the combined ANIZ and PRES variable value
        of the VW budget using the following equation:
        ``WVPRES + WVANIZ``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WVANIZ, indep, dataset = self.getVarForCalculations('WVANIZ', dataset)
        WVPRES, indep, dataset = self.getVarForCalculations('WVPRES', dataset)
        WV_ANIZ_PRES = WVANIZ + WVPRES
        return WV_ANIZ_PRES, indep

    def getVWResidual(self, dataset_override=None):
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

        Calculates the residual for the VW budget using
        the following equation:
        ``WVBT - (WVDFSN + WVSHEAR + WVADV + WVPRES + WVANIZ + WVBUOY + WVSDMP)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WVADV, indep, dataset = self.getVarForCalculations('WVADV', dataset)
        WVANIZ, indep, dataset = self.getVarForCalculations('WVANIZ', dataset)
        WVBT, indep, dataset = self.getVarForCalculations('WVBT', dataset)
        WVBUOY, indep, dataset = self.getVarForCalculations('WVBUOY', dataset)
        WVDFSN, indep, dataset = self.getVarForCalculations('WVDFSN', dataset)
        WVPRES, indep, dataset = self.getVarForCalculations('WVPRES', dataset)
        WVSHEAR, indep, dataset = self.getVarForCalculations('WVSHEAR', dataset)
        WVSDMP, indep, dataset = self.getVarForCalculations('WVSDMP', dataset)
        WV_RES = WVBT - (WVDFSN + WVSHEAR + WVADV + WVPRES + WVANIZ + WVBUOY + WVSDMP)
        return WV_RES, indep

    # NOTE: Not needed anymore?
    # def getU2C14OverTau(self, dataset_override = None):
    # formula = '- (3/2) * U2DFSN / np.maximum( TKE + TKES, 1e-6 )'

    # def getV2C14OverTau(self, dataset_override = None):
    # formula = ' - (3/2) * V2DFSN / np.maximum( TKE + TKES, 1e-6 )'

    # def getThl2C2OverTau(self, dataset_override = None):
    # formula = ' - ( THL2DISSIP + THL2DIFTR ) / np.maximum( THEL2, 1e-6 )'

    # def getQtog2C2OverTau(self, dataset_override = None):
    # formula = '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2*1e-6, 1e-15 )'

    # def getQ2C2OverTau(self, dataset_override = None):
    # formula = '- ( Q2DISSIP + Q2DIFTR ) / np.maximum( QT2 * 1e-6, 1e-15 )'
