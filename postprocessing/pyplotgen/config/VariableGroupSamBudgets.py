'''
:author: Steffen Domke
:date: Mid 2019
TODO:   - Should invalid variables produce zero lines or none at all?
        - Fallback <-> calc
        - Arrange lines so that styles match for different panels -> reduced momentum flux budgets
        - Does not work in the current implementation of pyplotgen, as SAM cannot be plotted alone
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupSamBudgets(VariableGroup):
    
    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None):
        """
        
        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "sam budget variables"
        
        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)
        
        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000
        
        # ?? Which conversion factors? What is the variable description??
        hlp_budget_lines = [
            {'var_names': ['HLADV'], 'legend_label': 'HLADV', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLDIFF'], 'legend_label': 'HLDIFF', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLRAD'], 'legend_label': 'HLRAD', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HLLAT'], 'legend_label': 'HLLAT', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['TTEND'], 'legend_label': 'TTEND', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['T_TNDCY'], 'legend_label': 'T_TNDCY'},
            {'var_names': ['HLSTOR'], 'legend_label': 'HLSTOR', 'sam_conv_factor': self.kg_per_second_to_kg_per_day},
            {'var_names': ['HL_RES'], 'legend_label': 'HL_RES', 'fallback_func': self.getHlResidual},
        ]
        
        # Total water budget (same as rtm) (does not exist in _MICRO_M2005_UWM)
        qtp_budget_lines = [
            {'var_names': ['QTADV'], 'legend_label': 'QTADV', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTDIFF'], 'legend_label': 'QTDIFF', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSRC'], 'legend_label': 'QTSRC', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSINK'], 'legend_label': 'QTSINK', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTEND'], 'legend_label': 'QTEND', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QTSTOR'], 'legend_label': 'QTSTOR', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QV_TNDCY'], 'legend_label': 'QV_TNDCY', 'sam_conv_factor': self.g_per_second_to_kg_per_day},
            {'var_names': ['QT_RES'], 'legend_label': 'QT_RES', 'fallback_func': self.getQtResidual},
        ]
        
        # Vertical liquid water static energy flux
        tpwp_budget_lines = [
            {'var_names': ['TWGRAD'], 'legend_label': 'TWGRAD'},
            {'var_names': ['TWADV'], 'legend_label': 'TWADV'},
            {'var_names': ['TWDIFF'], 'legend_label': 'TWDIFF'},
            {'var_names': ['TWB+P'], 'legend_label': 'TWBUOY+TWPRES', 'fallback_func': self.getTwBuoyPlusPres},
            {'var_names': ['TWPREC'], 'legend_label': 'TWPREC'},
            {'var_names': ['TWRAD'], 'legend_label': 'TWRAD'},
            {'var_names': ['TWFORC'], 'legend_label': 'TWFORC'},
            {'var_names': ['TWBT'], 'legend_label': 'TWBT'},
            {'var_names': ['TW_RES'], 'legend_label': 'TW_RES', 'fallback_func': self.getTwResidual},
        ]
        
        tpwp_split_budget_lines = [
            {'var_names': ['TWGRAD'], 'legend_label': 'TWGRAD'},
            {'var_names': ['TWADV'], 'legend_label': 'TWADV'},
            {'var_names': ['TWDIFF'], 'legend_label': 'TWDIFF'},
            {'var_names': ['TWBUOY'], 'legend_label': 'TWBUOY'},
            {'var_names': ['TWPRES'], 'legend_label': 'TWPRES'},
            {'var_names': ['TWPREC'], 'legend_label': 'TWPREC'},
            {'var_names': ['TWRAD'], 'legend_label': 'TWRAD'},
            {'var_names': ['TWFORC'], 'legend_label': 'TWFORC'},
            {'var_names': ['TWBT'], 'legend_label': 'TWBT'},
            {'var_names': ['TW_RES'], 'legend_label': 'TWRES', 'fallback_func': self.getTwResidual},
        ]
        
        # Exists ? Vertical liquid water pot. temperature flux
        thlpwp_budget_lines = [
            {'var_names': ['THLWGRAD'], 'legend_label': 'THLWGRAD'},
            {'var_names': ['THLWADV'], 'legend_label': 'THLWADV'},
            {'var_names': ['THLWDIFF'], 'legend_label': 'THLWDIFF'},
            {'var_names': ['THLWB+P'], 'legend_label': 'THLWBUOY+PRES', 'fallback_func': self.getThlwBuoyPlusPres},
            {'var_names': ['THLWPREC'], 'legend_label': 'THLWPREC'},
            {'var_names': ['THLWRAD'], 'legend_label': 'THLWRAD'},
            {'var_names': ['THLWFORC'], 'legend_label': 'THLWFORC'},
            {'var_names': ['THLWBT'], 'legend_label': 'THLWBT'},
            {'var_names': ['THLW_RES'], 'legend_label': 'THLW_RES', 'fallback_func': self.getThlwResidual},
        ]
        
        thlpwp_split_budget_lines = [
            {'var_names': ['THLWGRAD'], 'legend_label': 'THLWGRAD'},
            {'var_names': ['THLWADV'], 'legend_label': 'THLWADV'},
            {'var_names': ['THLWDIFF'], 'legend_label': 'THLWDIFF'},
            {'var_names': ['THLWBUOY'], 'legend_label': 'THLWBUOY'},
            {'var_names': ['THLWPRES'], 'legend_label': 'THLWPRES'},
            {'var_names': ['THLWPREC'], 'legend_label': 'THLWPREC'},
            {'var_names': ['THLWRAD'], 'legend_label': 'THLWRAD'},
            {'var_names': ['THLWFORC'], 'legend_label': 'THLWFORC'},
            {'var_names': ['THLWBT'], 'legend_label': 'THLWBT'},
            {'var_names': ['THLW_RES'], 'legend_label': 'THLW_RES', 'fallback_func': self.getThlwResidual},
        ]
        
        # Vertical total water flux budget
        qpwp_budget_lines = [
            {'var_names': ['QWGRAD'], 'legend_label': 'QWGRAD'},
            {'var_names': ['QWADV'], 'legend_label': 'QWADV'},
            {'var_names': ['QWDIFF'], 'legend_label': 'QWDIFF'},
            {'var_names': ['QWB+P'], 'legend_label': 'QWBUOY+QWPRES', 'fallback_func': self.getQwBuoyPlusPres},
            {'var_names': ['QWPREC'], 'legend_label': 'QWPREC'},
            {'var_names': ['QWFORC'], 'legend_label': 'QWFORC'},
            {'var_names': ['QWBT'], 'legend_label': 'QWBT'},
            {'var_names': ['QW_RES'], 'legend_label': 'QW_RES', 'fallback_func': self.getQwResidual},
        ]
        
        qpwp_split_budget_lines = [
            {'var_names': ['QWGRAD'], 'legend_label': 'QWGRAD'},
            {'var_names': ['QWADV'], 'legend_label': 'QWADV'},
            {'var_names': ['QWDIFF'], 'legend_label': 'QWDIFF'},
            {'var_names': ['QWBUOY'], 'legend_label': 'QWBUOY'},
            {'var_names': ['QWPRES'], 'legend_label': 'QWPRES'},
            {'var_names': ['QWPREC'], 'legend_label': 'QWPREC'},
            {'var_names': ['QWFORC'], 'legend_label': 'QWFORC'},
            {'var_names': ['QWBT'], 'legend_label': 'QWBT'},
            {'var_names': ['QW_RES'], 'legend_label': 'QW_RES', 'fallback_func': self.getQwResidual},
        ]
        
        qtogpwp_budget_lines = [
            {'var_names': ['QTOGWGRAD'], 'legend_label': 'QTOGWGRAD'},
            {'var_names': ['QTOGWADV'], 'legend_label': 'QTOGWADV'},
            {'var_names': ['QTOGWDIFF'], 'legend_label': 'QTOGWDIFF'},
            {'var_names': ['QTOGWB+P'], 'legend_label': 'QTOGWBUOY+PRES', 'fallback_func': self.getQtogwBuoyPlusPres},
            {'var_names': ['QTOGWPREC'], 'legend_label': 'QTOGWPREC'},
            {'var_names': ['QTOGWFORC'], 'legend_label': 'QTOGWFORC'},
            {'var_names': ['QTOGWBT'], 'legend_label': 'QTOGWBT'},
            {'var_names': ['QTOGW_RES'], 'legend_label': 'QTOGW_RES', 'fallback_func': self.getQtogwResidual},
        ]
        
        qtogpwp_split_budget_lines = [
            {'var_names': ['QTOGWGRAD'], 'legend_label': 'QTOGWGRAD'},
            {'var_names': ['QTOGWADV'], 'legend_label': 'QTOGWADV'},
            {'var_names': ['QTOGWDIFF'], 'legend_label': 'QTOGWDIFF'},
            {'var_names': ['QTOGWBUOY'], 'legend_label': 'QTOGWBUOY'},
            {'var_names': ['QTOGWPRES'], 'legend_label': 'QTOGWPRES'},
            {'var_names': ['QTOGWPREC'], 'legend_label': 'QTOGWPREC'},
            {'var_names': ['QTOGWFORC'], 'legend_label': 'QTOGWFORC'},
            {'var_names': ['QTOGWBT'], 'legend_label': 'QTOGWBT'},
            {'var_names': ['QTOGW_RES'], 'legend_label': 'QTOGW_RES', 'fallback_func': self.getQtogwResidual},
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
            {'var_names': ['T2_RES'], 'legend_label': 'T2_RES', 'fallback_func': self.getT2Residual},
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
            {'var_names': ['THL2_RES'], 'legend_label': 'THL2_RES', 'fallback_func': self.getThl2Residual},
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
            {'var_names': ['Q2_RES'], 'legend_label': 'Q2_RES', 'fallback_func': self.getQt2Residual},
        ]
        
        qtogp2_budget_lines = [
            {'var_names': ['QTOG2ADVTR'], 'legend_label': 'QTOG2ADVTR'},
            {'var_names': ['QTOG2GRAD'], 'legend_label': 'QTOG2GRAD'},
            {'var_names': ['QTOG2DISSIP'], 'legend_label': 'QTOG2DISSIP'},
            {'var_names': ['QTOG2DIFTR'], 'legend_label': 'QTOG2DIFTR'},
            {'var_names': ['QTOG2PREC'], 'legend_label': 'QTOG2PREC'},
            {'var_names': ['QTOG2FORC'], 'legend_label': 'QTOG2FORC'},
            {'var_names': ['QTOG2BT'], 'legend_label': 'QTOG2BT'},
            {'var_names': ['QTOG2_RES'], 'legend_label': 'QTOG2_RES', 'fallback_func': self.getQtog2Residual},
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
            {'var_names': ['QTHL_RES'], 'legend_label': 'QTHL_RES', 'fallback_func': self.getQThlResidual},
        ]
        
        ## Momentum flux budgets ##
        tke_budget_lines = [
            {'var_names': ['ADVTR'], 'legend_label': 'ADVTR'},
            {'var_names': ['SHEAR'], 'legend_label': 'SHEAR'},
            {'var_names': ['BUOYA'], 'legend_label': 'BUOYA'},
            {'var_names': ['PRESSTR'], 'legend_label': 'PRESSTR'},
            {'var_names': ['SDMP'], 'legend_label': 'SDMP'},
            {'var_names': ['DIFTR+DISS'], 'legend_label': 'DISSIP+DIFTR', 'fallback_func': self.getTkeDissPlusDiff},
            {'var_names': ['BT'], 'legend_label': 'BT'},
            {'var_names': ['TKE_RES'], 'legend_label': 'TKE_RES', 'fallback_func': self.getTkeResidual},
        ]
        
        tkes_budget_lines = [
            {'var_names': ['ADVTRS'], 'legend_label': 'ADVTRS'},
            {'var_names': ['SHEARS'], 'legend_label': 'SHEARS'},
            {'var_names': ['BUOYAS'], 'legend_label': 'BUOYAS'},
            {'var_names': ['DISSIPS'], 'legend_label': 'DISSIPS', 'sam_conv_factor': -1.},
            {'var_names': ['TKES_RES'], 'legend_label': 'TKES_RES', 'fallback_func': self.getTkesResidual},
        ]
        
        up2_vp2_budget_lines = [
            {'var_names': ['U2V2ADV'], 'legend_label': 'U2V2ADV', 'fallback_func': self.getU2V2Adv},
            {'var_names': ['U2V2BUOY'], 'legend_label': 'U2V2BUOY', 'fallback_func': self.getU2V2Buoy},
            {'var_names': ['U2V2PRESS'], 'legend_label': 'U2V2PRES', 'fallback_func': self.getU2V2Pres},
            {'var_names': ['W2REDIS'], 'legend_label': 'U2V2REDIS', 'sam_conv_factor': -1},
            {'var_names': ['U2V2DIFF'], 'legend_label': 'U2V2DIFF', 'fallback_func': self.getU2V2Diff},
            {'var_names': ['DISSIP'], 'legend_label': 'U2V2DISSIP', 'sam_conv_factor': 2},
            {'var_names': ['U2V2SDMP'], 'legend_label': 'U2V2SDMP', 'fallback_func': self.getU2V2Sdmp},
            {'var_names': ['SHEAR'], 'legend_label': 'U2V2SHEAR', 'sam_conv_factor': 2},
            {'var_names': ['U2V2BT'], 'legend_label': 'U2V2BT', 'fallback_func': self.getU2V2Bt},
            #{'var_names': ['U2V2_RES'], 'legend_label': '2 * TKE_RES - W2_RES'},
            {'var_names': ['U2V2_RES'], 'legend_label': 'U2V2_RES', 'fallback_func': self.getU2V2Residual},
        ]
        
        upwp_budget_lines = [
            {'var_names': ['WUDIFF'], 'legend_label': 'WUDIFF'},
            {'var_names': ['WU_RES'], 'legend_label': 'WU_RES', 'fallback_func': self.getUWResidual},
            {'var_names': ['WUADV'], 'legend_label': 'WUADV'},
            {'var_names': ['WUPRES'], 'legend_label': 'WUPRES'},
            {'var_names': ['WUANIZ'], 'legend_label': 'WUANIZ'},
            {'var_names': ['WUBUOY'], 'legend_label': 'WUBUOY'},
            {'var_names': ['WUSHEAR'], 'legend_label': 'WUSHEAR'},
            {'var_names': ['WUBT'], 'legend_label': 'WUBT'},
            #{'var_names': ['WUSDMP'], 'legend_label': 'WUSDMP'},
        ]
        
        vpwp_budget_lines = [
            {'var_names': ['WVDIFF'], 'legend_label': 'WVDIFF'},
            {'var_names': ['WV_RES'], 'legend_label': 'WV_RES', 'fallback_func': self.getVWResidual},
            {'var_names': ['WVADV'], 'legend_label': 'WVADV'},
            {'var_names': ['WVPRES'], 'legend_label': 'WVPRES'},
            {'var_names': ['WVANIZ'], 'legend_label': 'WVANIZ'},
            {'var_names': ['WVBUOY'], 'legend_label': 'WVBUOY'},
            {'var_names': ['WVSHEAR'], 'legend_label': 'WVSHEAR'},
            {'var_names': ['WVBT'], 'legend_label': 'WVBT'},
            #{'var_names': ['WVSDMP'], 'legend_label': 'WVSDMP'},
        ]
        
        up2_budget_lines = [
            {'var_names': ['U2ADV'], 'legend_label': 'U2ADV'},
            {'var_names': ['U2SHEAR'], 'legend_label': 'U2SHEAR'},
            {'var_names': ['U2REDIS'], 'legend_label': 'U2REDIS'},
            {'var_names': ['U2DIFF'], 'legend_label': 'U2DIFF'},
            {'var_names': ['U2BT'], 'legend_label': 'U2BT'},
            {'var_names': ['U2_RES'], 'legend_label': 'U2_RES', 'fallback_func': self.getU2Residual},
        ]
        
        vp2_budget_lines = [
            {'var_names': ['V2ADV'], 'legend_label': 'V2ADV'},
            {'var_names': ['V2SHEAR'], 'legend_label': 'V2SHEAR'},
            {'var_names': ['V2REDIS'], 'legend_label': 'V2REDIS'},
            {'var_names': ['V2DIFF'], 'legend_label': 'V2DIFF'},
            {'var_names': ['V2BT'], 'legend_label': 'V2BT'},
            {'var_names': ['V2_RES'], 'legend_label': 'V2_RES', 'fallback_func': self.getV2Residual},
        ]
        
        wp2_budget_lines = [
            {'var_names': ['W2ADV'], 'legend_label': 'W2ADV'},
            {'var_names': ['W2PRES'], 'legend_label': 'W2PRES'},
            {'var_names': ['W2REDIS'], 'legend_label': 'W2REDIS'},
            {'var_names': ['W2BUOY'], 'legend_label': 'W2BUOY'},
            {'var_names': ['W2DIFF'], 'legend_label': 'W2DIFF'},
            {'var_names': ['W2SDMP'], 'legend_label': 'W2SDMP'},
            {'var_names': ['W2BT'], 'legend_label': 'W2BT'},
            {'var_names': ['W2_RES'], 'legend_label': 'W2_RES', 'fallback_func': self.getW2Residual},
        ]
        
        wp3_budget_lines = [
            {'var_names': ['W3ADV'], 'legend_label': 'W3ADV'},
            {'var_names': ['W3PRES'], 'legend_label': 'W3PRES'},
            {'var_names': ['W3REDIS'], 'legend_label': 'W3REDIS'},
            {'var_names': ['W3BUOY'], 'legend_label': 'W3BUOY'},
            {'var_names': ['W3DIFF'], 'legend_label': 'W3DIFF'},
            {'var_names': ['W3BT'], 'legend_label': 'W3BT'},
            {'var_names': ['W3_RES'], 'legend_label': 'W3_RES', 'fallback_func': self.getW3Residual},
        ]
        
        ## Tau plots showing C_2/tau and C_14/tau ##
        # No support for 2nd axis, not needed
        #up2vp2tau_budget_lines = [
            #{'var_names': [r'U2DIFF'], 'legend_label': 'U2DIFF'},
            #{'var_names': [r'V2DIFF'], 'legend_label': 'V2DIFF'},
            #{'var_names': [r'TKES'], 'legend_label': 'TKES', 1, 1],
            #{'var_names': [r'TKE'], 'legend_label': 'TKE', 1, 1],
            #{'var_names': ['U2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (u'^2)$", 'fallback_func': self.getU2C14OverTau},
            #{'var_names': ['V2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (v'^2)$", 'fallback_func': self.getV2C14OverTau},
        #]
        
        #qtogp2tau_budget_lines = [
            #{'var_names': [r'QTOG2DISSIP'], 'legend_label': 'QTOG2DISSIP'},
            #{'var_names': [r'QTOG2DIFTR'], 'legend_label': 'QTOG2DIFTR'},
            #{'var_names': [r'QTOG2'], 'legend_label': 'QTOG2', 1, 1],
            #{'var_names': ['QTOG2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$", 'fallback_func': self.getQtog2C2OverTau},
        #]
        
        #qp2tau_budget_lines = [
            #{'var_names': [r'Q2DISSIP'], 'legend_label': 'Q2DISSIP'},
            #{'var_names': [r'Q2DIFTR'], 'legend_label': 'Q2DIFTR'},
            #{'var_names': [r'QT2'], 'legend_label': 'QT2', 1, 1],
            #{'var_names': ['QT2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q'^2)$", 'fallback_func': self.getQ2C2OverTau},
        #]
        
        #thlp2tau_budget_lines = [
            #{'var_names': ['THL2DISSIP'], 'legend_label': 'THL2DISSIP'},
            #{'var_names': ['THL2DIFTR'], 'legend_label': 'THL2DIFTR'},
            #{'var_names': ['THEL2'], 'legend_label': 'THEL2', 1, 1],
            #{'var_names': ['THL2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (\theta_l'^2)$", 'fallback_func': self.getThl2C2OverTau},
        #]
        
        #tau_comp_budget_lines = [
            #{'var_names': ['U2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (u'^2)$", 'fallback_func': self.getU2C14OverTau},
            #{'var_names': ['V2_C14_over_tau'], 'legend_label': r"$\frac{C_{14}}{\tau}\ (v'^2)$", 'fallback_func': self.getV2C14OverTau},
            #{'var_names': [r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$", 'legend_label': '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2 * 1e-6, 1e-15 )'], 1, 1],#??
            #{'var_names': ['QTOG2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q_{tog}'^2)$", 'fallback_func': self.getQtog2C2OverTau},
            #{'var_names': ['THL2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (\theta_l'^2)$", 'fallback_func': self.getThl2C2OverTau},
            #{'var_names': ['QT2_C2_over_tau'], 'legend_label': r"$\frac{C_2}{\tau}\ (q'^2)$", 'fallback_func': self.getQ2C2OverTau},
        #]
        
        ## Reduced second moment flux budgets ##TODO: How to sort lines and define colors?
        up2_reduced_budget_lines = [
            {'var_names': ['U2ADV'], 'legend_label': 'advection'},
            #{'var_names': ['U2BUOY'], 'legend_label': 'buoyancy', 'fallback_func': self.getNothing},
            {'var_names': ['U2DIFF'], 'legend_label': 'dissipation'},
            {'var_names': ['U2REDIS'], 'legend_label': 'isotropy'},
            #{'var_names': ['U2PRES'], 'legend_label': 'pressure', 'fallback_func': self.getNothing},
            {'var_names': ['U2SHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['U2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['U2_RES'], 'legend_label': 'residual', 'fallback_func': self.getU2Residual},
        ]

        vp2_reduced_budget_lines = [
            {'var_names': ['V2ADV'], 'legend_label': 'advection'},
            #{'var_names': ['V2BUOY'], 'legend_label': 'buoyancy', 'fallback_func': self.getNothing},
            {'var_names': ['V2DIFF'], 'legend_label': 'dissipation'},
            {'var_names': ['V2REDIS'], 'legend_label': 'isotropy'},
            #{'var_names': ['V2PRES'], 'legend_label': 'pressure', 'fallback_func': self.getNothing},
            {'var_names': ['V2SHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['V2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['V2_RES'], 'legend_label': 'residual', 'fallback_func': self.getV2Residual},
        ]
        
        wp2_reduced_budget_lines = [
            #{'var_names': ['W2SDMP'], 'legend_label': 'damping'},
            {'var_names': ['W2ADV'], 'legend_label': 'advection'},
            {'var_names': ['W2BUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['W2DIFF'], 'legend_label': 'dissipation'},
            {'var_names': ['W2REDIS+PRES'], 'legend_label': 'pressure', 'fallback_func': self.getW2RedisPlusPres},
            #{'var_names': ['W2SHEAR'], 'legend_label': 'turb. prod.', 'fallback_func': self.getNothing},
            {'var_names': ['W2BT'], 'legend_label': 'time tndcy'},
            {'var_names': ['W2_RES'], 'legend_label': 'residual', 'fallback_func': self.getW2Residual},
        ]

        upwp_reduced_budget_lines = [
            #{'var_names': ['WUSDMP'], 'legend_label': 'WUSDMP'},
            {'var_names': ['WUADV'], 'legend_label': 'advection'},
            {'var_names': ['WUBUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['WUDIFF'], 'legend_label': 'dissipation'},
            {'var_names': ['WUPRES+ANIZ'], 'legend_label': 'pressure', 'fallback_func': self.getUWPresPlusAniz},
            {'var_names': ['WUSHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['WUBT'], 'legend_label': 'time tndcy'},
            {'var_names': ['WU_RES'], 'legend_label': 'residual', 'fallback_func': self.getUWResidual},
        ]
        
        vpwp_reduced_budget_lines = [
            #{'var_names': ['WVSDMP'], 'legend_label': 'WUSDMP'},
            {'var_names': ['WVADV'], 'legend_label': 'advection'},
            {'var_names': ['WVBUOY'], 'legend_label': 'buoyancy'},
            {'var_names': ['WVDIFF'], 'legend_label': 'dissipation'},
            {'var_names': ['WVPRES+ANIZ'], 'legend_label': 'pressure', 'fallback_func': self.getVWPresPlusAniz},
            {'var_names': ['WVSHEAR'], 'legend_label': 'turb. prod.'},
            {'var_names': ['WVBT'], 'legend_label': 'time tndcy'},
            {'var_names': ['WV_RES'], 'legend_label': 'residual', 'fallback_func': self.getVWResidual},
        ]
        
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['HL'],
                'sam': ['HL'],
                'coamps': ['HL'],
                'r408': ['HL'],
                'hoc': ['HL'],
                'e3sm': ['HL'],
                },
            'lines': hlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': 'Liquid-Ice Static Energy Budget, LIMSE', 'axis_title': r"LIMSE budget terms $\mathrm{\left[K\,s^{-1}\right]}$",
            },
            
            {'var_names': {
                'clubb': ['QT'],
                'sam': ['QT'],
                'coamps': ['QT'],
                'r408': ['QT'],
                'hoc': ['QT'],
                'e3sm': ['QT'],
                },
            'lines': qtp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r'Total Water (No Rain/Snow Included) Budget, $\mathrm{r_t}$', 'axis_title': r"$\mathrm{r_t}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['TW_B+P'],
                'sam': ['TW_B+P'],
                'coamps': ['TW_B+P'],
                'r408': ['TW_B+P'],
                'hoc': ['TW_B+P'],
                'e3sm': ['TW_B+P'],
                },
            'lines': tpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Static Energy Flux Budget (Reduced), $\mathrm{\overline{s_v'w'}}$", 'axis_title': r"$\mathrm{\overline{s_v'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['TW_COMPLETE'],
                'sam': ['TW_COMPLETE'],
                'coamps': ['TW_COMPLETE'],
                'r408': ['TW_COMPLETE'],
                'hoc': ['TW_COMPLETE'],
                'e3sm': ['TW_COMPLETE'],
                },
            'lines': tpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Static Energy Flux Budget, $\mathrm{\overline{s_v'w'}}$", 'axis_title': r"$\mathrm{\overline{s_v'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['THLW_B+P'],
                'sam': ['THLW_B+P'],
                'coamps': ['THLW_B+P'],
                'r408': ['THLW_B+P'],
                'hoc': ['THLW_B+P'],
                'e3sm': ['THLW_B+P'],
                },
            'lines': thlpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Potential Temperature Flux Budget (Reduced), $\mathrm{\overline{\theta_l'w'}}$", 'axis_title': r"$\mathrm{\overline{\theta_l'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['THLW_COMPLETE'],
                'sam': ['THLW_COMPLETE'],
                'coamps': ['THLW_COMPLETE'],
                'r408': ['THLW_COMPLETE'],
                'hoc': ['THLW_COMPLETE'],
                'e3sm': ['THLW_COMPLETE'],
                },
            'lines': thlpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Potential Temperature Flux Budget, $\mathrm{\overline{\theta_l'w'}}$", 'axis_title': r"$\mathrm{\overline{\theta_l'w'}}$ budget terms $\mathrm{\left[m\,K\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['QW_B+P'],
                'sam': ['QW_B+P'],
                'coamps': ['QW_B+P'],
                'r408': ['QW_B+P'],
                'hoc': ['QW_B+P'],
                'e3sm': ['QW_B+P'],
                },
            'lines': qpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Total Water Flux Budget (Reduced), $\mathrm{\overline{r_t'w'}}$", 'axis_title': r"$\mathrm{\overline{r_t'w'}}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$"},
             
            {'var_names': {
                'clubb': ['QW_COMPLETE'],
                'sam': ['QW_COMPLETE'],
                'coamps': ['QW_COMPLETE'],
                'r408': ['QW_COMPLETE'],
                'hoc': ['QW_COMPLETE'],
                'e3sm': ['QW_COMPLETE'],
                },
            'lines': qpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Total Water Flux Budget, $\mathrm{\overline{r_t'w'}}$", 'axis_title': r"$\mathrm{\overline{r_t'w'}}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['QTOGW_B+P'],
                'sam': ['QTOGW_B+P'],
                'coamps': ['QTOGW_B+P'],
                'r408': ['QTOGW_B+P'],
                'hoc': ['QTOGW_B+P'],
                'e3sm': ['QTOGW_B+P'],
                },
            'lines': qtogpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Grand Total Water Flux Budget (Reduced), $\mathrm{\overline{q_{tog}'w'}}$", 'axis_title': r"$\mathrm{\overline{q_{tog}'w'}}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['QTOGW_COMPLETE'],
                'sam': ['QTOGW_COMPLETE'],
                'coamps': ['QTOGW_COMPLETE'],
                'r408': ['QTOGW_COMPLETE'],
                'hoc': ['QTOGW_COMPLETE'],
                'e3sm': ['QTOGW_COMPLETE'],
                },
            'lines': qtogpwp_split_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Grand Total Water Flux Budget, $\mathrm{\overline{q_{tog}'w'}}$", 'axis_title': r"$\mathrm{\overline{q_{tog}'w'}}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-2}\right]}$"},
            
            {'var_names': {
                'clubb': ['W2'],
                'sam': ['W2'],
                'coamps': ['W2'],
                'r408': ['W2'],
                'hoc': ['W2'],
                'e3sm': ['W2'],
                },
            'lines': wp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Vertical Momentum Variance Budget, $\mathrm{\overline{w'^2}}$", 'axis_title': r"$\mathrm{\overline{w'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['W3'],
                'sam': ['W3'],
                'coamps': ['W3'],
                'r408': ['W3'],
                'hoc': ['W3'],
                'e3sm': ['W3'],
                },
            'lines': wp3_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Vertical Momentum Skewness Budget, $\mathrm{\overline{w'^3}}$", 'axis_title': r"$\mathrm{\overline{w'^3}}$ budget terms $\mathrm{\left[m^3\,s^{-4}\right]}$"},
            
            {'var_names': {
                'clubb': ['T2'],
                'sam': ['T2'],
                'coamps': ['T2'],
                'r408': ['T2'],
                'hoc': ['T2'],
                'e3sm': ['T2'],
                },
            'lines': tp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Static Energy Variance Budget, $\mathrm{\overline{s_v'^2}}$", 'axis_title': r"$\mathrm{\overline{s_v'^2}}$ budget terms $\mathrm{\left[K^2\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['THL2'],
                'sam': ['THL2'],
                'coamps': ['THL2'],
                'r408': ['THL2'],
                'hoc': ['THL2'],
                'e3sm': ['THL2'],
                },
            'lines': thlp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Liquid Water Potential Temperature Variance Budget, $\mathrm{\overline{\theta_l'^2}}$", 'axis_title': r"$\mathrm{\overline{\theta_l'^2}}$ budget terms $\mathrm{\left[K^2\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['Q2'],
                'sam': ['Q2'],
                'coamps': ['Q2'],
                'r408': ['Q2'],
                'hoc': ['Q2'],
                'e3sm': ['Q2'],
                },
            'lines': rtp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Total Water Variance Budget, $\mathrm{\overline{r_t'^2}}$", 'axis_title': r"$\mathrm{\overline{r_t'^2}}$ budget terms $\mathrm{\left[kg^2\,kg^{-2}\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['QTOG2'],
                'sam': ['QTOG2'],
                'coamps': ['QTOG2'],
                'r408': ['QTOG2'],
                'hoc': ['QTOG2'],
                'e3sm': ['QTOG2'],
                },
            'lines': qtogp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Grand Total Water Variance Budget, $\mathrm{\overline{q_{tog}'^2}}$", 'axis_title': r"$\mathrm{\overline{q_{tog}'^2}}$ budget terms $\mathrm{\left[kg^2\,kg^{-2}\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['QTHL'],
                'sam': ['QTHL'],
                'coamps': ['QTHL'],
                'r408': ['QTHL'],
                'hoc': ['QTHL'],
                'e3sm': ['QTHL'],
                },
            'lines': qpthlp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Total Water, Liquid Water Pot. Temp. Covariance Budget, $\mathrm{\overline{q_t'\theta_l'}}$", 'axis_title': r"$\mathrm{\overline{q_t'\theta_l'}}$ budget terms $\mathrm{\left[kg\,kg^{-1}\,K\,s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['TKE'],
                'sam': ['TKE'],
                'coamps': ['TKE'],
                'r408': ['TKE'],
                'hoc': ['TKE'],
                'e3sm': ['TKE'],
                },
            'lines': tke_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Turbulence Kinetic Energy Budget, TKE", 'axis_title': r"TKE budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['TKES'],
                'sam': ['TKES'],
                'coamps': ['TKES'],
                'r408': ['TKES'],
                'hoc': ['TKES'],
                'e3sm': ['TKES'],
                },
            'lines': tkes_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Sub-Grid Turbulence Kinetic Energy Budget, $\mathrm{TKE_{SGS}}$", 'axis_title': r"$\mathrm{TKE_{SGS}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['U2V2'],
                'sam': ['U2V2'],
                'coamps': ['U2V2'],
                'r408': ['U2V2'],
                'hoc': ['U2V2'],
                'e3sm': ['U2V2'],
                },
            'lines': up2_vp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Horizontal Turbulence Kinetic Energy Budget, $\mathrm{\overline{u'^2}+\overline{v'^2}}$", 'axis_title': r"$\mathrm{\overline{u'^2}+\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\, s^{-3}\right]}$"},
            
            # TODO: Rename to UW?
            {'var_names': {
                'clubb': ['WU'],
                'sam': ['WU'],
                'coamps': ['WU'],
                'r408': ['WU'],
                'hoc': ['WU'],
                'e3sm': ['WU'],
                },
            'lines': upwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Eastward Momentum Flux Budget, $\mathrm{\overline{u'w'}}$", 'axis_title': r"$\mathrm{\overline{u'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['WV'],
                'sam': ['WV'],
                'coamps': ['WV'],
                'r408': ['WV'],
                'hoc': ['WV'],
                'e3sm': ['WV'],
                },
            'lines': vpwp_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Northward Momentum Flux Budget, $\mathrm{\overline{v'w'}}$", 'axis_title': r"$\mathrm{\overline{v'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['U2'],
                'sam': ['U2'],
                'coamps': ['U2'],
                'r408': ['U2'],
                'hoc': ['U2'],
                'e3sm': ['U2'],
                },
            'lines': up2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Eastward Momentum Variance Budget, $\mathrm{\overline{u'^2}}$", 'axis_title': r"$\mathrm{\overline{u'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['V2'],
                'sam': ['V2'],
                'coamps': ['V2'],
                'r408': ['V2'],
                'hoc': ['V2'],
                'e3sm': ['V2'],
                },
            'lines': vp2_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Northward Momentum Variance Budget, $\mathrm{\overline{v'^2}}$", 'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['V2_COMP'],
                #'sam': ['V2_COMP'],
                #'coamps': ['V2_COMP'],
                #'r408': ['V2_COMP'],
                #'hoc': ['V2_COMP'],
                #'e3sm': ['V2_COMP'],
                #},
            #'lines': vp2_comp_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\, s^{-3}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['U2_V2_TAU'],
                #'sam': ['U2_V2_TAU'],
                #'coamps': ['U2_V2_TAU'],
                #'r408': ['U2_V2_TAU'],
                #'hoc': ['U2_V2_TAU'],
                #'e3sm': ['U2_V2_TAU'],
                #},
            #'lines': up2vp2tau_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\frac{C_{14}}{\tau}}\ \mathrm{\left[s^{-1}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['QTOG2_TAU'],
                #'sam': ['QTOG2_TAU'],
                #'coamps': ['QTOG2_TAU'],
                #'r408': ['QTOG2_TAU'],
                #'hoc': ['QTOG2_TAU'],
                #'e3sm': ['QTOG2_TAU'],
                #},
            #'lines': qtogp2tau_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['Q2_TAU'],
                #'sam': ['Q2_TAU'],
                #'coamps': ['Q2_TAU'],
                #'r408': ['Q2_TAU'],
                #'hoc': ['Q2_TAU'],
                #'e3sm': ['Q2_TAU'],
                #},
            #'lines': qp2tau_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['THL2_TAU'],
                #'sam': ['THL2_TAU'],
                #'coamps': ['THL2_TAU'],
                #'r408': ['THL2_TAU'],
                #'hoc': ['THL2_TAU'],
                #'e3sm': ['THL2_TAU'],
                #},
            #'lines': qp2tau_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\frac{C_2}{\tau}}\ \mathrm{\left[s^{-1}\right]}$"},
            
            #{'var_names': {
                #'clubb': ['C2_TAU_COMP'],
                #'sam': ['C2_TAU_COMP'],
                #'coamps': ['C2_TAU_COMP'],
                #'r408': ['C2_TAU_COMP'],
                #'hoc': ['C2_TAU_COMP'],
                #'e3sm': ['C2_TAU_COMP'],
                #},
            #'lines': tau_comp_budget_lines, 'type': Panel.TYPE_BUDGET, 'axis_title': r"$\mathrm{\frac{C}{\tau}}\ \mathrm{\left[s^{-1}\right]}$"},
            
            {'var_names': {
                'clubb': ['U2 REDUCED'],
                'sam': ['U2 REDUCED'],
                'coamps': ['U2 REDUCED'],
                'r408': ['U2 REDUCED'],
                'hoc': ['U2 REDUCED'],
                'e3sm': ['U2 REDUCED'],
                },
            'lines': up2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Eastward Momentum Variance Budget (Reduced), $\mathrm{\overline{u'^2}}$", 'axis_title': r"$\mathrm{\overline{u'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['V2 REDUCED'],
                'sam': ['V2 REDUCED'],
                'coamps': ['V2 REDUCED'],
                'r408': ['V2 REDUCED'],
                'hoc': ['V2 REDUCED'],
                'e3sm': ['V2 REDUCED'],
                },
            'lines': vp2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Northward Momentum Variance Budget (Reduced), $\mathrm{\overline{v'^2}}$", 'axis_title': r"$\mathrm{\overline{v'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['W2 REDUCED'],
                'sam': ['W2 REDUCED'],
                'coamps': ['W2 REDUCED'],
                'r408': ['W2 REDUCED'],
                'hoc': ['W2 REDUCED'],
                'e3sm': ['W2 REDUCED'],
                },
            'lines': wp2_reduced_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Vertical Momentum Variance Budget, $\mathrm{\overline{w'^2}}$", 'axis_title': r"$\mathrm{\overline{w'^2}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['WU REDUCED'],
                'sam': ['WU REDUCED'],
                'coamps': ['WU REDUCED'],
                'r408': ['WU REDUCED'],
                'hoc': ['WU REDUCED'],
                'e3sm': ['WU REDUCED'],
                },
            'lines': upwp_reduced_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Eastward Momentum Flux Budget (Reduced), $\mathrm{\overline{u'w'}}$", 'axis_title': r"$\mathrm{\overline{u'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
            
            {'var_names': {
                'clubb': ['WV REDUCED'],
                'sam': ['WV REDUCED'],
                'coamps': ['WV REDUCED'],
                'r408': ['WV REDUCED'],
                'hoc': ['WV REDUCED'],
                'e3sm': ['WV REDUCED'],
                },
            'lines': vpwp_reduced_budget_lines, 'type': Panel.TYPE_BUDGET, 'title': r"Northward Momentum Flux Budget (Reduced), $\mathrm{\overline{v'w'}}$", 'axis_title': r"$\mathrm{\overline{v'w'}}$ budget terms $\mathrm{\left[m^2\,s^{-3}\right]}$"},
        ]
        
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets= e3sm_datasets)

    def getHlResidual(self, dataset_override = None):
        """
        Calculates the residual for the HL budget using
        the following equation:
        (HLSTOR)+((-1)*(HLADV+HLDIFF+HLRAD+HLLAT+TTEND))*g_per_second_to_kg_per_day
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        
        HLADV, z, dataset = self.getVarForCalculations('HLADV', dataset)
        HLADV *= self.kg_per_second_to_kg_per_day
        
        HLDIFF, z, dataset = self.getVarForCalculations('HLDIFF', dataset)
        HLDIFF *= self.kg_per_second_to_kg_per_day
        
        HLLAT, z, dataset = self.getVarForCalculations('HLLAT', dataset)
        HLLAT *= self.kg_per_second_to_kg_per_day
        
        HLRAD, z, dataset = self.getVarForCalculations('HLRAD', dataset)
        HLRAD *= self.kg_per_second_to_kg_per_day
        
        HLSTOR, z, dataset = self.getVarForCalculations('HLSTOR', dataset)
        HLSTOR *= self.kg_per_second_to_kg_per_day
        
        TTEND, z, dataset = self.getVarForCalculations('TTEND', dataset)
        TTEND *= self.kg_per_second_to_kg_per_day
        
        #ALERT: Conversion factor correct?
        HL_RES = (HLSTOR - (HLADV + HLDIFF + HLLAT + HLRAD + TTEND)) * self.g_per_second_to_kg_per_day
        return HL_RES, z
    
    def getQtResidual(self, dataset_override = None):
        """
        Calculates the residual for the QT budget using
        the following equation:
        QTSTOR+(-1)*(QTADV+QTDIFF+QTSRC+QTSINK+QTEND)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QTADV, z, dataset = self.getVarForCalculations('QTADV', dataset)
        QTADV *=  self.g_per_second_to_kg_per_day
        QTDIFF, z, dataset = self.getVarForCalculations('QTDIFF', dataset)
        QTDIFF *=  self.g_per_second_to_kg_per_day
        QTEND, z, dataset = self.getVarForCalculations('QTEND', dataset)
        QTEND*=  self.g_per_second_to_kg_per_day
        QTSINK, z, dataset = self.getVarForCalculations('QTSINK', dataset)
        QTSINK *=  self.g_per_second_to_kg_per_day
        QTSRC, z, dataset = self.getVarForCalculations('QTSRC', dataset)
        QTSRC *=  self.g_per_second_to_kg_per_day
        QTSTOR, z, dataset = self.getVarForCalculations('QTSTOR', dataset)
        QTSTOR *=  self.g_per_second_to_kg_per_day
        QV_TNDCY, z, dataset = self.getVarForCalculations('QV_TNDCY', dataset)
        QV_TNDCY*=  self.g_per_second_to_kg_per_day
        QT_RES = QTSTOR - (QTADV + QTDIFF + QTEND + QTSRC + QTSINK)
        return QT_RES, z

    def getTwBuoyPlusPres(self, dataset_override = None):
        """
        Calculates the combined BUOY and PRES variable value
        of the TW budget using the following equation:
        TWBUOY+TWPRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        TWBUOY, z, dataset = self.getVarForCalculations('TWBUOY', dataset)
        TWPRES, z, dataset = self.getVarForCalculations('TWPRES', dataset)
        TW_BUOY_PRES = TWBUOY + TWPRES
        return TW_BUOY_PRES, z
    
    def getTwResidual(self, dataset_override = None):
        """
        Calculates the residual for the TW budget using
        the following equation:
        TWBT - (TWGRAD + TWADV + TWDIFF + TWBUOY + TWPRES + TWPREC + TWRAD + TWFORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        TWADV, z, dataset = self.getVarForCalculations('TWADV', dataset)
        TWBT, z, dataset = self.getVarForCalculations('TWBT', dataset)
        TWBUOY, z, dataset = self.getVarForCalculations('TWBUOY', dataset)
        TWDIFF, z, dataset = self.getVarForCalculations('TWDIFF', dataset)
        TWFORC, z, dataset = self.getVarForCalculations('TWFORC', dataset)
        TWGRAD, z, dataset = self.getVarForCalculations('TWGRAD', dataset)
        TWPREC, z, dataset = self.getVarForCalculations('TWPREC', dataset)
        TWPRES, z, dataset = self.getVarForCalculations('TWPRES', dataset)
        TWRAD, z, dataset = self.getVarForCalculations('TWRAD', dataset)
        TW_RES = TWBT - (TWADV + TWBUOY + TWDIFF + TWFORC + TWGRAD + TWPREC + TWPRES + TWRAD)
        return TW_RES, z
    
    def getThlwBuoyPlusPres(self, dataset_override = None):
        """
        Calculates the combined BUOY and PRES variable value
        of the THLW budget using the following equation:
        THLWBUOY + THLWPRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        THLWBUOY, z, dataset = self.getVarForCalculations('THLWBUOY', dataset)
        THLWPRES, z, dataset = self.getVarForCalculations('THLWPRES', dataset)
        THLW_BUOY_PRES = THLWBUOY + THLWPRES
        return THLW_BUOY_PRES, z

    def getThlwResidual(self, dataset_override = None):
        """
        Calculates the residual for the THLW budget using
        the following equation:
        THLWBT - (THLWGRAD + THLWADV + THLWDIFF + THLWBUOY + THLWPRES + THLWPREC + THLWRAD + THLWFORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        THLWADV, z, dataset = self.getVarForCalculations('THLWADV', dataset);
        THLWBT, z, dataset = self.getVarForCalculations('THLWBT', dataset)
        THLWBUOY, z, dataset = self.getVarForCalculations('THLWBUOY', dataset)
        THLWDIFF, z, dataset = self.getVarForCalculations('THLWDIFF', dataset)
        THLWFORC, z, dataset = self.getVarForCalculations('THLWFORC', dataset)
        THLWGRAD, z, dataset = self.getVarForCalculations('THLWGRAD', dataset)
        THLWPREC, z, dataset = self.getVarForCalculations('THLWPREC', dataset)
        THLWPRES, z, dataset = self.getVarForCalculations('THLWPRES', dataset)
        THLWRAD, z, dataset = self.getVarForCalculations('THLWRAD', dataset)
        THLW_RES = THLWBT - (THLWADV + THLWBUOY + THLWDIFF + THLWFORC + THLWGRAD + THLWPREC + THLWPRES + THLWRAD)
        return THLW_RES, z

    def getQwBuoyPlusPres(self, dataset_override = None):
        """
        Calculates the combined BUOY and PRES variable value
        of the QW budget using the following equation:
        QWBUOY+QWPRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QWBUOY, z, dataset = self.getVarForCalculations('QWBUOY', dataset)
        QWPRES, z, dataset = self.getVarForCalculations('QWPRES', dataset)
        QW_BUOY_PRES = QWBUOY + QWPRES
        return QW_BUOY_PRES, z

    def getQwResidual(self, dataset_override = None):
        """
        Calculates the residual for the QW budget using
        the following equation:
        QWBT - (QWGRAD + QWADV + QWDIFF + QWBUOY + QWPRES + QWPREC + QWFORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QWADV, z, dataset = self.getVarForCalculations('QWADV', dataset)
        QWBT, z, dataset = self.getVarForCalculations('QWBT', dataset)
        QWBUOY, z, dataset = self.getVarForCalculations('QWBUOY', dataset)
        QWDIFF, z, dataset = self.getVarForCalculations('QWDIFF', dataset)
        QWFORC, z, dataset = self.getVarForCalculations('QWFORC', dataset)
        QWGRAD, z, dataset = self.getVarForCalculations('QWGRAD', dataset)
        QWPREC, z, dataset = self.getVarForCalculations('QWPREC', dataset)
        QWPRES, z, dataset = self.getVarForCalculations('QWPRES', dataset)
        QW_RES = QWBT - (QWGRAD + QWADV + QWDIFF + QWBUOY + QWPRES + QWPREC + QWFORC)
        return QW_RES, z

    def getQtogwBuoyPlusPres(self, dataset_override = None):
        """
        Calculates the combined BUOY and PRES variable value
        of the QTOGW budget using the following equation:
        QTOGWBUOY + QTOGWPRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QTOGWBUOY, z, dataset = self.getVarForCalculations('QTOGWBUOY', dataset)
        QTOGWPRES, z, dataset = self.getVarForCalculations('QTOGWPRES', dataset)
        QTOGW_BUOY_PRES = QTOGWBUOY + QTOGWPRES
        return QTOGW_BUOY_PRES, z

    def getQtogwResidual(self, dataset_override = None):
        """
        Calculates the residual for the QTOGW budget using
        the following equation:
        QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDIFF + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QTOGWADV, z, dataset = self.getVarForCalculations('QTOGWADV', dataset)
        QTOGWBT, z, dataset = self.getVarForCalculations('QTOGWBT', dataset)
        QTOGWBUOY, z, dataset = self.getVarForCalculations('QTOGWBUOY', dataset)
        QTOGWDIFF, z, dataset = self.getVarForCalculations('QTOGWDIFF', dataset)
        QTOGWFORC, z, dataset = self.getVarForCalculations('QTOGWFORC', dataset)
        QTOGWGRAD, z, dataset = self.getVarForCalculations('QTOGWGRAD', dataset)
        QTOGWPREC, z, dataset = self.getVarForCalculations('QTOGWPREC', dataset)
        QTOGWPRES, z, dataset = self.getVarForCalculations('QTOGWPRES', dataset)
        QTOGW_RES = QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDIFF + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)
        return QTOGW_RES, z

    def getT2Residual(self, dataset_override = None):
        """
        Calculates the residual for the T2 budget using
        the following equation:
        T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        T2ADVTR, z, dataset = self.getVarForCalculations('T2ADVTR', dataset)
        T2BT, z, dataset = self.getVarForCalculations('T2BT', dataset)
        T2DISSIP, z, dataset = self.getVarForCalculations('T2DISSIP', dataset)
        T2DIFTR, z, dataset = self.getVarForCalculations('T2DIFTR', dataset)
        T2FORC, z, dataset = self.getVarForCalculations('T2FORC', dataset)
        T2GRAD, z, dataset = self.getVarForCalculations('T2GRAD', dataset)
        T2PREC, z, dataset = self.getVarForCalculations('T2PREC', dataset)
        T2RAD, z, dataset = self.getVarForCalculations('T2RAD', dataset)
        T2_RES = T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)
        return T2_RES, z

    def getThl2Residual(self, dataset_override = None):
        """
        Calculates the residual for the THL2 budget using
        the following equation:
        THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        THL2ADVTR, z, dataset = self.getVarForCalculations('THL2ADVTR', dataset)
        THL2BT, z, dataset = self.getVarForCalculations('THL2BT', dataset)
        THL2DISSIP, z, dataset = self.getVarForCalculations('THL2DISSIP', dataset)
        THL2DIFTR, z, dataset = self.getVarForCalculations('THL2DIFTR', dataset)
        THL2FORC, z, dataset = self.getVarForCalculations('THL2FORC', dataset)
        THL2GRAD, z, dataset = self.getVarForCalculations('THL2GRAD', dataset)
        THL2PREC, z, dataset = self.getVarForCalculations('THL2PREC', dataset)
        THL2RAD, z, dataset = self.getVarForCalculations('THL2RAD', dataset)
        THL2_RES = THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)
        return THL2_RES, z

    def getQt2Residual(self, dataset_override = None):
        """
        Calculates the residual for the Q2 budget using
        the following equation:
        Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        Q2ADVTR, z, dataset = self.getVarForCalculations('Q2ADVTR', dataset)
        Q2BT, z, dataset = self.getVarForCalculations('Q2BT', dataset)
        Q2DISSIP, z, dataset = self.getVarForCalculations('Q2DISSIP', dataset)
        Q2DIFTR, z, dataset = self.getVarForCalculations('Q2DIFTR', dataset)
        Q2FORC, z, dataset = self.getVarForCalculations('Q2FORC', dataset)
        Q2GRAD, z, dataset = self.getVarForCalculations('Q2GRAD', dataset)
        Q2PREC, z, dataset = self.getVarForCalculations('Q2PREC', dataset)
        Q2_RES = Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)
        return Q2_RES, z

    def getQtog2Residual(self, dataset_override = None):
        """
        Calculates the residual for the QTOG2 budget using
        the following equation:
        QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QTOG2ADVTR, z, dataset = self.getVarForCalculations('QTOG2ADVTR', dataset)
        QTOG2BT, z, dataset = self.getVarForCalculations('QTOG2BT', dataset)
        QTOG2DIFTR, z, dataset = self.getVarForCalculations('QTOG2DIFTR', dataset)
        QTOG2DISSIP, z, dataset = self.getVarForCalculations('QTOG2DISSIP', dataset)
        QTOG2FORC, z, dataset = self.getVarForCalculations('QTOG2FORC', dataset)
        QTOG2GRAD, z, dataset = self.getVarForCalculations('QTOG2GRAD', dataset)
        QTOG2PREC, z, dataset = self.getVarForCalculations('QTOG2PREC', dataset)
        QTOG2_RES = QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)
        return QTOG2_RES, z

    def getQThlResidual(self, dataset_override = None):
        """
        Calculates the residual for the QTHL budget using
        the following equation:
        QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QTHLADV, z, dataset = self.getVarForCalculations('QTHLADV', dataset)
        QTHLBT, z, dataset = self.getVarForCalculations('QTHLBT', dataset)
        QTHLDIFTR, z, dataset = self.getVarForCalculations('QTHLDIFTR', dataset)
        QTHLDISSIP, z, dataset = self.getVarForCalculations('QTHLDISSIP', dataset)
        QTHLFORC, z, dataset = self.getVarForCalculations('QTHLFORC', dataset)
        QTHLGRAD, z, dataset = self.getVarForCalculations('QTHLGRAD', dataset)
        QTHLPREC, z, dataset = self.getVarForCalculations('QTHLPREC', dataset)
        QTHLRAD, z, dataset = self.getVarForCalculations('QTHLRAD', dataset)
        QTHLW_RES = QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)
        return QTHLW_RES, z

    def getTkeDissPlusDiff(self, dataset_override = None):
        """
        Calculates the combined DISS and DIFTR variable value
        of the TKE budget using the following equation:
        DISSIP + DIFTR
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        DIFTR, z, dataset = self.getVarForCalculations('DIFTR', dataset)
        DISSIP, z, dataset = self.getVarForCalculations('DISSIP', dataset)
        TKE_DISS_DIFF = DIFTR + DISSIP
        return TKE_DISS_DIFF, z

    def getTkeResidual(self, dataset_override = None):
        """
        Calculates the residual for the TKE budget using
        the following equation:
        BT - (SHEAR + BUOYA + ADVTR + PRESSTR + DIFTR + SDMP + DISSIP)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        ADVTR, z, dataset = self.getVarForCalculations('ADVTR', dataset)
        BT, z, dataset = self.getVarForCalculations('BT', dataset)
        BUOYA, z, dataset = self.getVarForCalculations('BUOYA', dataset)
        DIFTR, z, dataset = self.getVarForCalculations('DIFTR', dataset)
        DISSIP, z, dataset = self.getVarForCalculations('DISSIP', dataset)
        PRESSTR, z, dataset = self.getVarForCalculations('PRESSTR', dataset)
        SDMP, z, dataset = self.getVarForCalculations('SDMP', dataset)
        SHEAR, z, dataset = self.getVarForCalculations('SHEAR', dataset)
        TKE_RES = BT - (SHEAR + BUOYA + ADVTR + PRESSTR + DIFTR + SDMP + DISSIP)
        return TKE_RES, z
    
    def getTkesResidual(self, dataset_override = None):
        """
        Calculates the residual for the TKES budget using
        the following equation:
        -(SHEARS + BUOYAS + ADVTRS + DISSIPS)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        ADVTRS, z, dataset = self.getVarForCalculations('ADVTRS', dataset)
        BUOYAS, z, dataset = self.getVarForCalculations('BUOYAS', dataset)
        DISSIPS, z, dataset = self.getVarForCalculations('DISSIPS', dataset)
        SHEARS, z, dataset = self.getVarForCalculations('SHEARS', dataset)
        TKES_RES = -(SHEARS + BUOYAS + ADVTRS + DISSIPS)
        return TKES_RES, z
    
    def getU2Residual(self, dataset_override = None):
        """
        Calculates the residual for the U2 budget using
        the following equation:
        U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DIFF)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        U2ADV, z, dataset = self.getVarForCalculations('U2ADV', dataset)
        U2BT, z, dataset = self.getVarForCalculations('U2BT', dataset)
        U2DIFF, z, dataset = self.getVarForCalculations('U2DIFF', dataset)
        U2REDIS, z, dataset = self.getVarForCalculations('U2REDIS', dataset)
        U2SHEAR, z, dataset = self.getVarForCalculations('U2SHEAR', dataset)
        U2_RES = U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DIFF)
        return U2_RES, z
    
    def getV2Residual(self, dataset_override = None):
        """
        Calculates the residual for the V2 budget using
        the following equation:
        V2BT - (V2ADV + V2SHEAR + V2REDIS + V2DIFF)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        V2ADV, z, dataset = self.getVarForCalculations('V2ADV', dataset)
        V2BT, z, dataset = self.getVarForCalculations('V2BT', dataset)
        V2DIFF, z, dataset = self.getVarForCalculations('V2DIFF', dataset)
        V2REDIS, z, dataset = self.getVarForCalculations('V2REDIS', dataset)
        V2SHEAR, z, dataset = self.getVarForCalculations('V2SHEAR', dataset)
        V2_RES = V2BT - (V2ADV + V2SHEAR + V2REDIS + V2DIFF)
        return V2_RES, z
    
    def getW2RedisPlusPres(self, dataset_override = None):
        """
        Calculates the combined REDIS and PRES variable value
        of the W2 budget using the following equation:
        W2REDIS + W2PRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W2PRES, z, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, z, dataset = self.getVarForCalculations('W2REDIS', dataset)
        W2_REDIS_PRES = W2REDIS + W2PRES
        return W2_REDIS_PRES, z
    
    def getW2Residual(self, dataset_override = None):
        """
        Calculates the residual for the W2 budget using
        the following equation:
        W2BT - (W2ADV + W2PRES + W2REDIS + W2BUOY + W2DIFF + W2SDMP)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W2ADV, z, dataset = self.getVarForCalculations('W2ADV', dataset)
        W2BT, z, dataset = self.getVarForCalculations('W2BT', dataset)
        W2BUOY, z, dataset = self.getVarForCalculations('W2BUOY', dataset)
        W2DIFF, z, dataset = self.getVarForCalculations('W2DIFF', dataset)
        W2PRES, z, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, z, dataset = self.getVarForCalculations('W2REDIS', dataset)
        W2SDMP, z, dataset = self.getVarForCalculations('W2SDMP', dataset)
        W2_RES = W2BT - (W2ADV + W2PRES + W2REDIS + W2BUOY + W2DIFF + W2SDMP)
        return W2_RES, z
    
    def getU2V2Adv(self, dataset_override):
        """
        Calculates the ADV variable value
        of the combined U2+V2 budget using the following equation:
        2 * ADVTR - W2ADV
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        ADVTR, z, dataset = self.getVarForCalculations('ADVTR', dataset)
        W2ADV, z, dataset = self.getVarForCalculations('W2ADV', dataset)
        U2V2_ADV = 2 * ADVTR + W2ADV
        return U2V2_ADV, z
    
    def getU2V2Buoy(self, dataset_override):
        """
        Calculates the BUOY variable value
        of the combined U2+V2 budget using the following equation:
        2 * BUOYA - W2BUOY
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        BUOYA, z, dataset = self.getVarForCalculations('BUOYA', dataset)
        W2BUOY, z, dataset = self.getVarForCalculations('W2BUOY', dataset)
        U2V2_BUOY = 2 * BUOYA + W2BUOY
        return U2V2_BUOY, z
        
    
    def getU2V2Pres(self, dataset_override):
        """
        Calculates the PRES variable value
        of the combined U2+V2 budget using the following equation:
        2 * PRESSTR - W2PRES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        PRESSTR, z, dataset = self.getVarForCalculations('PRESSTR', dataset)
        W2PRES, z, dataset = self.getVarForCalculations('W2PRES', dataset)
        U2V2_PRES = 2 * PRESSTR + W2PRES
        return U2V2_PRES, z
    
    def getU2V2Diff(self, dataset_override):
        """
        Calculates the DIFF variable value
        of the combined U2+V2 budget using the following equation:
        2 * DIFTR - W2DIFF
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        DIFTR, z, dataset = self.getVarForCalculations('DIFTR', dataset)
        W2DIFF, z, dataset = self.getVarForCalculations('W2DIFF', dataset)
        U2V2_DIFF = 2 * DIFTR + W2DIFF
        return U2V2_DIFF, z
    
    def getU2V2Sdmp(self, dataset_override):
        """
        Calculates the SDMP variable value
        of the combined U2+V2 budget using the following equation:
        2 * SDMP - W2SDMP
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        SDMP, z, dataset = self.getVarForCalculations('SDMP', dataset)
        W2SDMP, z, dataset = self.getVarForCalculations('W2SDMP', dataset)
        U2V2_SDMP = 2 * SDMP + W2SDMP
        return U2V2_SDMP, z
    
    def getU2V2Bt(self, dataset_override):
        """
        Calculates the BT variable value
        of the combined U2+V2 budget using the following equation:
        2 * BT - W2BT
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        BT, z, dataset = self.getVarForCalculations('BT', dataset)
        W2BT, z, dataset = self.getVarForCalculations('W2BT', dataset)
        U2V2_BT = 2 * BT + W2BT
        return U2V2_BT, z
    
    def getU2V2Residual(self, dataset_override):
        """
        Calculates the residual for the W2 budget using
        the following equation:
        2 * TKE_RES - W2_RES
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        ADVTR, z, dataset = self.getVarForCalculations('ADVTR', dataset)
        W2ADV, z, dataset = self.getVarForCalculations('W2ADV', dataset)
        BT, z, dataset = self.getVarForCalculations('BT', dataset)
        W2BT, z, dataset = self.getVarForCalculations('W2BT', dataset)
        BUOYA, z, dataset = self.getVarForCalculations('BUOYA', dataset)
        W2BUOY, z, dataset = self.getVarForCalculations('W2BUOY', dataset)
        DIFTR, z, dataset = self.getVarForCalculations('DIFTR', dataset)
        W2DIFF, z, dataset = self.getVarForCalculations('W2DIFF', dataset)
        DISSIP, z, dataset = self.getVarForCalculations('DISSIP', dataset)
        PRESSTR, z, dataset = self.getVarForCalculations('PRESSTR', dataset)
        W2PRES, z, dataset = self.getVarForCalculations('W2PRES', dataset)
        W2REDIS, z, dataset = self.getVarForCalculations('W2REDIS', dataset)
        SDMP, z, dataset = self.getVarForCalculations('SDMP', dataset)
        W2SDMP, z, dataset = self.getVarForCalculations('W2SDMP', dataset)
        SHEAR, z, dataset = self.getVarForCalculations('SHEAR', dataset)
        U2V2_RES = 2. * (BT - (ADVTR + BUOYA + PRESSTR + DIFTR + DISSIP + SDMP + SHEAR)) - W2BT + (W2ADV + W2BUOY + W2PRES + W2DIFF + W2SDMP + W2REDIS)
        return U2V2_RES, z
    
    def getW3Residual(self, dataset_override = None):
        """
        Calculates the residual for the W3 budget using
        the following equation:
        W3BT - (W3ADV + W3PRES + W3REDIS + W3BUOY + W3DIFF)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W3ADV, z, dataset = self.getVarForCalculations('W3ADV', dataset)
        W3BT, z, dataset = self.getVarForCalculations('W3BT', dataset)
        W3BUOY, z, dataset = self.getVarForCalculations('W3BUOY', dataset)
        W3DIFF, z, dataset = self.getVarForCalculations('W3DIFF', dataset)
        W3PRES, z, dataset = self.getVarForCalculations('W3PRES', dataset)
        W3REDIS, z, dataset = self.getVarForCalculations('W3REDIS', dataset)
        W3_RES = W3BT - (W3ADV + W3PRES + W3REDIS + W3BUOY + W3DIFF)
        return W3_RES, z
    
    def getUWPresPlusAniz(self, dataset_override = None):
        """
        Calculates the combined ANIZ and PRES variable value
        of the UW budget using the following equation:
        WUPRES + WUANIZ
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        WUANIZ, z, dataset = self.getVarForCalculations('WUANIZ', dataset)
        WUPRES, z, dataset = self.getVarForCalculations('WUPRES', dataset)
        WU_ANIZ_PRES = WUANIZ + WUPRES
        return WU_ANIZ_PRES, z
    
    def getUWResidual(self, dataset_override = None):
        """
        Calculates the residual for the UW budget using
        the following equation:
        WUBT - (WUDIFF + WUSHEAR + WUADV + WUPRES + WUANIZ + WUBUOY + WUSDMP)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        WUADV, z, dataset = self.getVarForCalculations('WUADV', dataset)
        WUANIZ, z, dataset = self.getVarForCalculations('WUANIZ', dataset)
        WUBT, z, dataset = self.getVarForCalculations('WUBT', dataset)
        WUBUOY, z, dataset = self.getVarForCalculations('WUBUOY', dataset)
        WUDIFF, z, dataset = self.getVarForCalculations('WUDIFF', dataset)
        WUPRES, z, dataset = self.getVarForCalculations('WUPRES', dataset)
        WUSHEAR, z, dataset = self.getVarForCalculations('WUSHEAR', dataset)
        WUSDMP, z, dataset = self.getVarForCalculations('WUSDMP', dataset)
        WU_RES = WUBT - (WUDIFF + WUSHEAR + WUADV + WUPRES + WUANIZ + WUBUOY + WUSDMP)
        return WU_RES, z
    
    def getVWPresPlusAniz(self, dataset_override = None):
        """
        Calculates the combined ANIZ and PRES variable value
        of the VW budget using the following equation:
        WVPRES + WVANIZ
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        WVANIZ, z, dataset = self.getVarForCalculations('WVANIZ', dataset)
        WVPRES, z, dataset = self.getVarForCalculations('WVPRES', dataset)
        WV_ANIZ_PRES = WVANIZ + WVPRES
        return WV_ANIZ_PRES, z
    
    def getVWResidual(self, dataset_override = None):
        """
        Calculates the residual for the VW budget using
        the following equation:
        WVBT - (WVDIFF + WVSHEAR + WVADV + WVPRES + WVANIZ + WVBUOY + WVSDMP)
        :param dataset_override:
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        WVADV, z, dataset = self.getVarForCalculations('WVADV', dataset)
        WVANIZ, z, dataset = self.getVarForCalculations('WVANIZ', dataset)
        WVBT, z, dataset = self.getVarForCalculations('WVBT', dataset)
        WVBUOY, z, dataset = self.getVarForCalculations('WVBUOY', dataset)
        WVDIFF, z, dataset = self.getVarForCalculations('WVDIFF', dataset)
        WVPRES, z, dataset = self.getVarForCalculations('WVPRES', dataset)
        WVSHEAR, z, dataset = self.getVarForCalculations('WVSHEAR', dataset)
        WVSDMP, z, dataset = self.getVarForCalculations('WVSDMP', dataset)
        WV_RES = WVBT - (WVDIFF + WVSHEAR + WVADV + WVPRES + WVANIZ + WVBUOY + WVSDMP)
        return WV_RES, z
    
    #NOTE: Not needed anymore?
    #def getU2C14OverTau(self, dataset_override = None):
    #formula = '- (3/2) * U2DIFF / np.maximum( TKE + TKES, 1e-6 )'
    
    #def getV2C14OverTau(self, dataset_override = None):
    #formula = ' - (3/2) * V2DIFF / np.maximum( TKE + TKES, 1e-6 )'
    
    #def getThl2C2OverTau(self, dataset_override = None):
    #formula = ' - ( THL2DISSIP + THL2DIFTR ) / np.maximum( THEL2, 1e-6 )'
    
    #def getQtog2C2OverTau(self, dataset_override = None):
    #formula = '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2*1e-6, 1e-15 )'
    
    #def getQ2C2OverTau(self, dataset_override = None):
    #formula = '- ( Q2DISSIP + Q2DIFTR ) / np.maximum( QT2 * 1e-6, 1e-15 )'