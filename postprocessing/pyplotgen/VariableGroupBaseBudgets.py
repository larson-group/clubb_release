'''
:author: Nicolas Strike
:date: Mid 2019
'''

from Panel import Panel
from VariableGroup import VariableGroup
from Line import Line


class VariableGroupBaseBudgets(VariableGroup):
    '''

    '''

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "base variables budgets"
        
        thlm_lines = [
            {'clubb_name': 'thlm_bt', 'label': 'thlm_bt'},
            {'clubb_name': 'thlm_ma', 'label': 'thlm_ma'},
            {'clubb_name': 'thlm_ta', 'label': 'thlm_ta'},
            {'clubb_name': 'thlm_mc', 'label': 'thlm_mc'},
            {'clubb_name': 'thlm_clipping', 'label': 'thlm_bt', 'fallback_func': self.getThlmClipping},
            {'clubb_name': 'radht', 'label': 'radht'},
            {'clubb_name': 'lsforcing', 'label': 'lsforcing', 'fallback_func': self.getLsforcing},
            {'clubb_name': 'thlm_residual', 'label': 'thlm_residual', 'fallback_func': self.getThlmResidual},

        ]
        
        rtm_lines = [
        {'clubb_name': 'rtm_bt', 'label': 'rtm_bt'},
            {'clubb_name': 'rtm_ma', 'label': 'rtm_ma'},
            {'clubb_name': 'rtm_ta', 'label': 'rtm_ta'},
            {'clubb_name': 'rtm_mc', 'label': 'rtm_mc'},
            {'clubb_name': 'rtm_clipping', 'label': 'rtm_bt', 'fallback_func': self.getRtmClipping},
            {'clubb_name': 'rtm_pd', 'label': 'rtm_pd'},
            {'clubb_name': 'rtm_forcing', 'label': 'rtm_forcing', 'fallback_func': self.getRtmForcing},
            {'clubb_name': 'rtm_residual', 'label': 'rtm_residual', 'fallback_func': self.getRtmResidual},

        ]

        wpthlp_lines = [
            {'clubb_name': 'wpthlp_bt', 'label': 'wpthlp_bt'},
            {'clubb_name': 'wpthlp_ma', 'label': 'wpthlp_ma'},
            {'clubb_name': 'wpthlp_ta', 'label': 'wpthlp_ta'},
            {'clubb_name': 'wpthlp_tp', 'label': 'wpthlp_tp'},
            {'clubb_name': 'wpthlp_ac', 'label': 'wpthlp_ac'},
            {'clubb_name': 'wpthlp_bp', 'label': 'wpthlp_bp'},
            {'clubb_name': 'wpthlp_pr1', 'label': 'wpthlp_pr1'},
            {'clubb_name': 'wpthlp_pr2', 'label': 'wpthlp_pr2'},
            {'clubb_name': 'wpthlp_pr3', 'label': 'wpthlp_pr3'},
            {'clubb_name': 'wpthlp_dp1', 'label': 'wpthlp_dp1'},
            {'clubb_name': 'wpthlp_mfl', 'label': 'wpthlp_mfl'},
            {'clubb_name': 'wpthlp_cl', 'label': 'wpthlp_cl'},
            {'clubb_name': 'wpthlp_sicl', 'label': 'wpthlp_sicl'},
            {'clubb_name': 'wpthlp_forcing', 'label': 'wpthlp_forcing'},
            {'clubb_name': 'wpthlp_residual', 'label': 'wpthlp_residual', 'fallback_func': self.getWpthlpResidual},

        ]

        wprtp_lines = [
            {'clubb_name': 'wprtp_bt', 'label': 'wprtp_bt'},
            {'clubb_name': 'wprtp_ma', 'label': 'wprtp_ma'},
            {'clubb_name': 'wprtp_ta', 'label': 'wprtp_ta'},
            {'clubb_name': 'wprtp_tp', 'label': 'wprtp_tp'},
            {'clubb_name': 'wprtp_ac', 'label': 'wprtp_ac'},
            {'clubb_name': 'wprtp_bp', 'label': 'wprtp_bp'},
            {'clubb_name': 'wprtp_pr1', 'label': 'wprtp_pr1'},
            {'clubb_name': 'wprtp_pr2', 'label': 'wprtp_pr2'},
            {'clubb_name': 'wprtp_pr3', 'label': 'wprtp_pr3'},
            {'clubb_name': 'wprtp_dp1', 'label': 'wprtp_dp1'},
            {'clubb_name': 'wprtp_mfl', 'label': 'wprtp_mfl'},
            {'clubb_name': 'wprtp_cl', 'label': 'wprtp_cl'},
            {'clubb_name': 'wprtp_sicl', 'label': 'wprtp_sicl'},
            {'clubb_name': 'wprtp_forcing', 'label': 'wprtp_forcing'},
            {'clubb_name': 'wprtp_residual', 'label': 'wprtp_residual', 'fallback_func': self.getWprtpResidual},
        ]

        wp2_lines = [
            {'clubb_name': 'wp2_bt', 'label': 'wp2_bt'},
            {'clubb_name': 'wp2_ma', 'label': 'wp2_ma'},
            {'clubb_name': 'wp2_ta', 'label': 'wp2_ta'},
            {'clubb_name': 'wp2_ac', 'label': 'wp2_ac'},
            {'clubb_name': 'wp2_bp', 'label': 'wp2_bp'},
            {'clubb_name': 'wp2_pr1', 'label': 'wp2_pr1'},
            {'clubb_name': 'wp2_pr2', 'label': 'wp2_pr2'},
            {'clubb_name': 'wp2_pr3', 'label': 'wp2_pr3'},
            {'clubb_name': 'wp2_dp1', 'label': 'wp2_dp1'},
            {'clubb_name': 'wp2_dp2', 'label': 'wp2_dp2'},
            {'clubb_name': 'wp2_cl', 'label': 'wp2_cl'},
            {'clubb_name': 'wp2_pd', 'label': 'wp2_pd'},
            {'clubb_name': 'wp2_splat', 'label': 'wp2_splat'},
            {'clubb_name': 'wp2_sf', 'label': 'wp2_sf'},
            {'clubb_name': 'wp2_residual', 'label': 'wp2_residual', 'fallback_func': self.getWp2Residual},
        ]
        
        wp3_lines = [
            {'clubb_name': 'wp3_bt', 'label': 'wp3_bt'},
            {'clubb_name': 'wp3_ma', 'label': 'wp3_ma'},
            {'clubb_name': 'wp3_ta', 'label': 'wp3_ta'},
            {'clubb_name': 'wp3_ac', 'label': 'wp3_ac'},
            {'clubb_name': 'wp3_pr1', 'label': 'wp3_pr1'},
            {'clubb_name': 'wp3_pr2', 'label': 'wp3_pr2'},
            {'clubb_name': 'wp3_pr3', 'label': 'wp3_pr3'},
            {'clubb_name': 'wp3_bp1', 'label': 'wp3_bp1'},
            {'clubb_name': 'wp3_bp2', 'label': 'wp3_bp2'},
            {'clubb_name': 'wp3_dp1', 'label': 'wp3_dp1'},
            {'clubb_name': 'wp3_tp', 'label': 'wp3_tp'},
            {'clubb_name': 'wp3_cl', 'label': 'wp3_cl'},
            {'clubb_name': 'wp3_splat', 'label': 'wp3_splat'},
            {'clubb_name': 'wp3_residual', 'label': 'wp3_residual', 'fallback_func': self.getWp3Residual},
        ]
        
        thlp2_lines = [
            {'clubb_name': 'thlp2_bt', 'label': 'thlp2_bt'},
            {'clubb_name': 'thlp2_ma', 'label': 'thlp2_ma'},
            {'clubb_name': 'thlp2_ta', 'label': 'thlp2_ta'},
            {'clubb_name': 'thlp2_tp', 'label': 'thlp2_tp'},
            {'clubb_name': 'thlp2_dp1', 'label': 'thlp2_dp1'},
            {'clubb_name': 'thlp2_dp2', 'label': 'thlp2_dp2'},
            {'clubb_name': 'thlp2_cl', 'label': 'thlp2_cl'},
            {'clubb_name': 'thlp2_pd', 'label': 'thlp2_pd'},
            {'clubb_name': 'thlp2_sf', 'label': 'thlp2_sf'},
            {'clubb_name': 'thlp2_forcing', 'label': 'thlp2_forcing'},
            {'clubb_name': 'thlp2_residual', 'label': 'thlp2_residual', 'fallback_func': self.getThlp2Residual},
        ]
        
        rtp2_lines = [
            {'clubb_name': 'rtp2_bt', 'label': 'rtp2_bt'},
            {'clubb_name': 'rtp2_ma', 'label': 'rtp2_ma'},
            {'clubb_name': 'rtp2_ta', 'label': 'rtp2_ta'},
            {'clubb_name': 'rtp2_tp', 'label': 'rtp2_tp'},
            {'clubb_name': 'rtp2_dp1', 'label': 'rtp2_dp1'},
            {'clubb_name': 'rtp2_dp2', 'label': 'rtp2_dp2'},
            {'clubb_name': 'rtp2_cl', 'label': 'rtp2_cl'},
            {'clubb_name': 'rtp2_pd', 'label': 'rtp2_pd'},
            {'clubb_name': 'rtp2_sf', 'label': 'rtp2_sf'},
            {'clubb_name': 'rtp2_forcing', 'label': 'rtp2_forcing'},
            {'clubb_name': 'rtp2_residual', 'label': 'rtp2_residual', 'fallback_func': self.getRtp2Residual},
        ]
        
        rtpthlp_lines = [
            {'clubb_name': 'rtpthlp_bt', 'label': 'rtpthlp_bt'},
            {'clubb_name': 'rtpthlp_ma', 'label': 'rtpthlp_ma'},
            {'clubb_name': 'rtpthlp_ta', 'label': 'rtpthlp_ta'},
            {'clubb_name': 'rtpthlp_tp1', 'label': 'rtpthlp_tp1'},
            {'clubb_name': 'rtpthlp_dp1', 'label': 'rtpthlp_dp1'},
            {'clubb_name': 'rtpthlp_dp2', 'label': 'rtpthlp_dp2'},
            {'clubb_name': 'rtpthlp_cl', 'label': 'rtpthlp_cl'},
            {'clubb_name': 'rtpthlp_tp2', 'label': 'rtpthlp_tp2'},
            {'clubb_name': 'rtpthlp_sf', 'label': 'rtpthlp_sf'},
            {'clubb_name': 'rtpthlp_forcing', 'label': 'rtpthlp_forcing'},
            {'clubb_name': 'rtpthlp_residual', 'label': 'rtpthlp_residual', 'fallback_func': self.getRtpthlpResidual},
        ]
        
        upwp_lines = [
            {'clubb_name': 'upwp_bt', 'label': 'upwp_bt'},
            {'clubb_name': 'upwp_ma', 'label': 'upwp_ma'},
            {'clubb_name': 'upwp_ta', 'label': 'upwp_ta'},
            {'clubb_name': 'upwp_tp', 'label': 'upwp_tp'},
            {'clubb_name': 'upwp_ac', 'label': 'upwp_ac'},
            {'clubb_name': 'upwp_bp', 'label': 'upwp_bp'},
            {'clubb_name': 'upwp_pr1', 'label': 'upwp_pr1'},
            {'clubb_name': 'upwp_pr2', 'label': 'upwp_pr2'},
            {'clubb_name': 'upwp_pr3', 'label': 'upwp_pr3'},
            {'clubb_name': 'upwp_pr4', 'label': 'upwp_pr4'},
            {'clubb_name': 'upwp_dp1', 'label': 'upwp_dp1'},
            {'clubb_name': 'upwp_cl', 'label': 'upwp_cl'},
            {'clubb_name': 'upwp_mfl', 'label': 'upwp_mfl'},
            {'clubb_name': 'upwp_residual', 'label': 'upwp_residual', 'fallback_func': self.getUpwpResidual},
        ]
        
        vpwp_lines = [
            {'clubb_name': 'vpwp_bt', 'label': 'vpwp_bt'},
            {'clubb_name': 'vpwp_ma', 'label': 'vpwp_ma'},
            {'clubb_name': 'vpwp_ta', 'label': 'vpwp_ta'},
            {'clubb_name': 'vpwp_tp', 'label': 'vpwp_tp'},
            {'clubb_name': 'vpwp_ac', 'label': 'vpwp_ac'},
            {'clubb_name': 'vpwp_bp', 'label': 'vpwp_bp'},
            {'clubb_name': 'vpwp_pr1', 'label': 'vpwp_pr1'},
            {'clubb_name': 'vpwp_pr2', 'label': 'vpwp_pr2'},
            {'clubb_name': 'vpwp_pr3', 'label': 'vpwp_pr3'},
            {'clubb_name': 'vpwp_pr4', 'label': 'vpwp_pr4'},
            {'clubb_name': 'vpwp_dp1', 'label': 'vpwp_dp1'},
            {'clubb_name': 'vpwp_cl', 'label': 'vpwp_cl'},
            {'clubb_name': 'vpwp_mfl', 'label': 'vpwp_mfl'},
            {'clubb_name': 'vpwp_residual', 'label': 'vpwp_residual', 'fallback_func': self.getVpwpResidual},
        ]
        
        self.variable_definitions = [
            {'clubb_name': 'thlm', 'lines': thlm_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'rtm', 'lines': rtm_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'wpthlp', 'lines': wpthlp_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'wprtp', 'lines': wprtp_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'wp2', 'lines': wp2_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'wp3', 'lines': wp3_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'thlp2', 'lines': thlp2_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'rtp2', 'lines': rtp2_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'rtpthlp', 'lines': rtpthlp_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'upwp', 'lines': upwp_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},
            {'clubb_name': 'vpwp', 'lines': vpwp_lines, 'type': Panel.TYPE_BUDGET, 'fill_zeros': True},

        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file)

    def getThlmClipping(self):
        '''


        thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        :return:
        '''

        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        thlm_mfl = self.__getFallbackVar__('thlm_mfl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_cl = self.__getFallbackVar__('thlm_cl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_tacl = self.__getFallbackVar__('thlm_tacl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_sdmp = self.__getFallbackVar__('thlm_sdmp', self.ncdf_files['zt'], fill_zeros=True)

        output_data = thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        output_line = Line(output_data, z_ncdf, line_format="", label="thlm_clipping")

        return output_line

    def getLsforcing(self):
        '''


        thlm_forcing-radht-thlm_mc
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        thlm_forcing = self.__getFallbackVar__('thlm_forcing', self.ncdf_files['zt'], fill_zeros=True)
        radht = self.__getFallbackVar__('radht', self.ncdf_files['zt'], fill_zeros=True)
        thlm_mc = self.__getFallbackVar__('thlm_mc', self.ncdf_files['zt'], fill_zeros=True)

        output_data = thlm_forcing - radht - thlm_mc
        output_line = Line(output_data, z_ncdf, line_format="", label="lsforcing")

        return output_line
    
    def getThlmResidual(self):
        '''


        thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)
        :return: 
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        thlm_mfl = self.__getFallbackVar__('thlm_mfl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_cl = self.__getFallbackVar__('thlm_cl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_tacl = self.__getFallbackVar__('thlm_tacl', self.ncdf_files['zt'], fill_zeros=True)
        thlm_sdmp = self.__getFallbackVar__('thlm_sdmp', self.ncdf_files['zt'], fill_zeros=True)
        thlm_bt = self.__getFallbackVar__('thlm_bt', self.ncdf_files['zt'], fill_zeros=True)
        thlm_ta = self.__getFallbackVar__('thlm_ta', self.ncdf_files['zt'], fill_zeros=True)
        thlm_forcing = self.__getFallbackVar__('thlm_forcing', self.ncdf_files['zt'], fill_zeros=True)
        thlm_ma = self.__getFallbackVar__('thlm_ma', self.ncdf_files['zt'], fill_zeros=True)

        output_data = thlm_bt-(thlm_ma+thlm_ta+thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp+thlm_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="thlm_clipping")

        return output_line

    def getRtmClipping(self):
        '''


        rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        rtm_mfl = self.__getFallbackVar__('rtm_mfl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_cl = self.__getFallbackVar__('rtm_cl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_tacl = self.__getFallbackVar__('rtm_tacl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_sdmp = self.__getFallbackVar__('rtm_sdmp', self.ncdf_files['zt'], fill_zeros=True)

        output_data = rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp
        output_line = Line(output_data, z_ncdf, line_format="", label="rtm_clipping")

        return output_line

    def getRtmForcing(self):
        '''


        rtm_forcing - rtm_mc
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        rtm_mc = self.__getFallbackVar__('rtm_mc', self.ncdf_files['zt'], fill_zeros=True)
        rtm_forcing = self.__getFallbackVar__('rtm_forcing', self.ncdf_files['zt'], fill_zeros=True)

        output_data = rtm_forcing - rtm_mc
        output_line = Line(output_data, z_ncdf, line_format="", label="rtm_forcinng")

        return output_line

    def getRtmResidual(self):
        '''


        rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        rtm_mfl = self.__getFallbackVar__('rtm_mfl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_cl = self.__getFallbackVar__('rtm_cl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_tacl = self.__getFallbackVar__('rtm_tacl', self.ncdf_files['zt'], fill_zeros=True)
        rtm_sdmp = self.__getFallbackVar__('rtm_sdmp', self.ncdf_files['zt'], fill_zeros=True)
        rtm_bt = self.__getFallbackVar__('rtm_bt', self.ncdf_files['zt'], fill_zeros=True)
        rtm_ta = self.__getFallbackVar__('rtm_ta', self.ncdf_files['zt'], fill_zeros=True)
        rtm_forcing = self.__getFallbackVar__('rtm_forcing', self.ncdf_files['zt'], fill_zeros=True)
        rtm_pd = self.__getFallbackVar__('rtm_pd', self.ncdf_files['zt'], fill_zeros=True)
        rtm_ma = self.__getFallbackVar__('rtm_ma', self.ncdf_files['zt'], fill_zeros=True)

        output_data = rtm_bt - (rtm_ma + rtm_ta + rtm_mfl + rtm_cl + rtm_tacl + rtm_sdmp + rtm_forcing + rtm_pd)
        output_line = Line(output_data, z_ncdf, line_format="", label="rtm_residual")

        return output_line

    def getWpthlpResidual(self):
        '''


        wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_mfl = self.__getFallbackVar__('wpthlp_mfl', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_cl = self.__getFallbackVar__('wpthlp_cl', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_tp = self.__getFallbackVar__('wpthlp_tp', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_ac = self.__getFallbackVar__('wpthlp_ac', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_pr1 = self.__getFallbackVar__('wpthlp_pr1', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_pr3 = self.__getFallbackVar__('wpthlp_pr3', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_pr2 = self.__getFallbackVar__('wpthlp_pr2', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_dp1 = self.__getFallbackVar__('wpthlp_dp1', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_sicl = self.__getFallbackVar__('wpthlp_sicl', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_bt = self.__getFallbackVar__('wpthlp_bt', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_ta = self.__getFallbackVar__('wpthlp_ta', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_forcing = self.__getFallbackVar__('wpthlp_forcing', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_bp = self.__getFallbackVar__('wpthlp_bp', self.ncdf_files['zm'], fill_zeros=True)
        wpthlp_ma = self.__getFallbackVar__('wpthlp_ma', self.ncdf_files['zm'], fill_zeros=True)

        output_data = wpthlp_bt - (wpthlp_ma + wpthlp_ta + wpthlp_tp + wpthlp_ac + wpthlp_bp + wpthlp_pr1 + wpthlp_pr2 + wpthlp_pr3 + wpthlp_dp1 + wpthlp_mfl + wpthlp_cl + wpthlp_sicl + wpthlp_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="wpthlp_residual")

        return output_line
    
    def getWprtpResidual(self):
        '''


        wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_mfl = self.__getFallbackVar__('wprtp_mfl', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_cl = self.__getFallbackVar__('wprtp_cl', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_tp = self.__getFallbackVar__('wprtp_tp', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_ac = self.__getFallbackVar__('wprtp_ac', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_pr1 = self.__getFallbackVar__('wprtp_pr1', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_pr3 = self.__getFallbackVar__('wprtp_pr3', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_pr2 = self.__getFallbackVar__('wprtp_pr2', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_dp1 = self.__getFallbackVar__('wprtp_dp1', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_sicl = self.__getFallbackVar__('wprtp_sicl', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_bt = self.__getFallbackVar__('wprtp_bt', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_ta = self.__getFallbackVar__('wprtp_ta', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_forcing = self.__getFallbackVar__('wprtp_forcing', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_bp = self.__getFallbackVar__('wprtp_bp', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_ma = self.__getFallbackVar__('wprtp_ma', self.ncdf_files['zm'], fill_zeros=True)
        wprtp_pd = self.__getFallbackVar__('wprtp_pd', self.ncdf_files['zm'], fill_zeros=True)

        output_data = wprtp_bt - (wprtp_ma + wprtp_ta + wprtp_tp + wprtp_ac + wprtp_bp + wprtp_pr1 + wprtp_pr2 + wprtp_pr3 + wprtp_dp1 + wprtp_mfl + wprtp_cl + wprtp_sicl + wprtp_pd + wprtp_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="wprtp_residual")

        return output_line

    def getWp2Residual(self):
        '''


        wp2_bt - (wp2_ma + wp2_ta + wp2_tp + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_mfl + wp2_cl + wp2_sicl + wp2_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        wp2_sf = self.__getFallbackVar__('wp2_sf', self.ncdf_files['zm'], fill_zeros=True)
        wp2_cl = self.__getFallbackVar__('wp2_cl', self.ncdf_files['zm'], fill_zeros=True)
        wp2_ac = self.__getFallbackVar__('wp2_ac', self.ncdf_files['zm'], fill_zeros=True)
        wp2_pr1 = self.__getFallbackVar__('wp2_pr1', self.ncdf_files['zm'], fill_zeros=True)
        wp2_pr3 = self.__getFallbackVar__('wp2_pr3', self.ncdf_files['zm'], fill_zeros=True)
        wp2_pr2 = self.__getFallbackVar__('wp2_pr2', self.ncdf_files['zm'], fill_zeros=True)
        wp2_dp1 = self.__getFallbackVar__('wp2_dp1', self.ncdf_files['zm'], fill_zeros=True)
        wp2_dp2 = self.__getFallbackVar__('wp2_dp2', self.ncdf_files['zm'], fill_zeros=True)
        wp2_bt = self.__getFallbackVar__('wp2_bt', self.ncdf_files['zm'], fill_zeros=True)
        wp2_ta = self.__getFallbackVar__('wp2_ta', self.ncdf_files['zm'], fill_zeros=True)
        wp2_splat = self.__getFallbackVar__('wp2_splat', self.ncdf_files['zm'], fill_zeros=True)
        wp2_bp = self.__getFallbackVar__('wp2_bp', self.ncdf_files['zm'], fill_zeros=True)
        wp2_ma = self.__getFallbackVar__('wp2_ma', self.ncdf_files['zm'], fill_zeros=True)
        wp2_pd = self.__getFallbackVar__('wp2_pd', self.ncdf_files['zm'], fill_zeros=True)

        output_data = wp2_bt - (wp2_ma + wp2_ta + wp2_ac + wp2_bp + wp2_pr1 + wp2_pr2 + wp2_pr3 + wp2_dp1 + wp2_dp2 + wp2_cl + wp2_pd + wp2_sf + wp2_splat)
        output_line = Line(output_data, z_ncdf, line_format="", label="wp2_residual")

        return output_line

    def getWp3Residual(self):
        '''


        wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl+wp3_splat)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zt'], fill_zeros=True)
        wp3_bp1 = self.__getFallbackVar__('wp3_bp1', self.ncdf_files['zt'], fill_zeros=True)
        wp3_bp2 = self.__getFallbackVar__('wp3_bp2', self.ncdf_files['zt'], fill_zeros=True)
        wp3_cl = self.__getFallbackVar__('wp3_cl', self.ncdf_files['zt'], fill_zeros=True)
        wp3_ac = self.__getFallbackVar__('wp3_ac', self.ncdf_files['zt'], fill_zeros=True)
        wp3_pr1 = self.__getFallbackVar__('wp3_pr1', self.ncdf_files['zt'], fill_zeros=True)
        wp3_pr3 = self.__getFallbackVar__('wp3_pr3', self.ncdf_files['zt'], fill_zeros=True)
        wp3_pr2 = self.__getFallbackVar__('wp3_pr2', self.ncdf_files['zt'], fill_zeros=True)
        wp3_dp1 = self.__getFallbackVar__('wp3_dp1', self.ncdf_files['zt'], fill_zeros=True)
        wp3_bt = self.__getFallbackVar__('wp3_bt', self.ncdf_files['zt'], fill_zeros=True)
        wp3_ta = self.__getFallbackVar__('wp3_ta', self.ncdf_files['zt'], fill_zeros=True)
        wp3_splat = self.__getFallbackVar__('wp3_splat', self.ncdf_files['zt'], fill_zeros=True)
        wp3_ma = self.__getFallbackVar__('wp3_ma', self.ncdf_files['zt'], fill_zeros=True)
        wp3_tp = self.__getFallbackVar__('wp3_tp', self.ncdf_files['zt'], fill_zeros=True)


        output_data = wp3_bt - (wp3_ma + wp3_ta + wp3_tp + wp3_ac + wp3_bp1 + wp3_bp2 + wp3_pr1 + wp3_pr2 + wp3_pr3 + wp3_dp1 + wp3_cl+wp3_splat)
        output_line = Line(output_data, z_ncdf, line_format="", label="wp3_residual")

        return output_line

    def getThlp2Residual(self):
        '''


        thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_cl = self.__getFallbackVar__('thlp2_cl', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_dp2 = self.__getFallbackVar__('thlp2_dp2', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_forcing = self.__getFallbackVar__('thlp2_forcing', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_sf = self.__getFallbackVar__('thlp2_sf', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_dp1 = self.__getFallbackVar__('thlp2_dp1', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_bt = self.__getFallbackVar__('thlp2_bt', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_ta = self.__getFallbackVar__('thlp2_ta', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_pd = self.__getFallbackVar__('thlp2_pd', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_ma = self.__getFallbackVar__('thlp2_ma', self.ncdf_files['zm'], fill_zeros=True)
        thlp2_tp = self.__getFallbackVar__('thlp2_tp', self.ncdf_files['zm'], fill_zeros=True)


        output_data = thlp2_bt - (thlp2_ma + thlp2_ta + thlp2_tp + thlp2_dp1 + thlp2_dp2 + thlp2_cl + thlp2_pd + thlp2_sf + thlp2_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="thlp2_residual")

        return output_line

    def getRtp2Residual(self):
        '''


        rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_cl = self.__getFallbackVar__('rtp2_cl', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_dp2 = self.__getFallbackVar__('rtp2_dp2', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_forcing = self.__getFallbackVar__('rtp2_forcing', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_sf = self.__getFallbackVar__('rtp2_sf', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_dp1 = self.__getFallbackVar__('rtp2_dp1', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_bt = self.__getFallbackVar__('rtp2_bt', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_ta = self.__getFallbackVar__('rtp2_ta', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_pd = self.__getFallbackVar__('rtp2_pd', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_ma = self.__getFallbackVar__('rtp2_ma', self.ncdf_files['zm'], fill_zeros=True)
        rtp2_tp = self.__getFallbackVar__('rtp2_tp', self.ncdf_files['zm'], fill_zeros=True)


        output_data = rtp2_bt - (rtp2_ma + rtp2_ta + rtp2_tp + rtp2_dp1 + rtp2_dp2 + rtp2_cl + rtp2_pd + rtp2_sf + rtp2_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="rtp2_residual")

        return output_line

    def getRtpthlpResidual(self):
        '''


        rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_pd + rtpthlp_sf + rtpthlp_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_cl = self.__getFallbackVar__('rtpthlp_cl', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_dp2 = self.__getFallbackVar__('rtpthlp_dp2', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_forcing = self.__getFallbackVar__('rtpthlp_forcing', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_sf = self.__getFallbackVar__('rtpthlp_sf', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_dp1 = self.__getFallbackVar__('rtpthlp_dp1', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_bt = self.__getFallbackVar__('rtpthlp_bt', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_ta = self.__getFallbackVar__('rtpthlp_ta', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_tp2 = self.__getFallbackVar__('rtpthlp_tp2', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_ma = self.__getFallbackVar__('rtpthlp_ma', self.ncdf_files['zm'], fill_zeros=True)
        rtpthlp_tp1 = self.__getFallbackVar__('rtpthlp_tp1', self.ncdf_files['zm'], fill_zeros=True)


        output_data = rtpthlp_bt - (rtpthlp_ma + rtpthlp_ta + rtpthlp_tp1 + rtpthlp_tp2 + rtpthlp_dp1 + rtpthlp_dp2 + rtpthlp_cl + rtpthlp_sf + rtpthlp_forcing)
        output_line = Line(output_data, z_ncdf, line_format="", label="rtpthlp_residual")

        return output_line
    
    def getUpwpResidual(self):
        '''


        upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_dp1 + upwp_dp2 + upwp_cl + upwp_pd + upwp_sf + upwp_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        upwp_cl = self.__getFallbackVar__('upwp_cl', self.ncdf_files['zm'], fill_zeros=True)
        upwp_tp = self.__getFallbackVar__('upwp_tp', self.ncdf_files['zm'], fill_zeros=True)
        upwp_ac = self.__getFallbackVar__('upwp_ac', self.ncdf_files['zm'], fill_zeros=True)
        upwp_bp = self.__getFallbackVar__('upwp_bp', self.ncdf_files['zm'], fill_zeros=True)
        upwp_dp1 = self.__getFallbackVar__('upwp_dp1', self.ncdf_files['zm'], fill_zeros=True)
        upwp_bt = self.__getFallbackVar__('upwp_bt', self.ncdf_files['zm'], fill_zeros=True)
        upwp_ta = self.__getFallbackVar__('upwp_ta', self.ncdf_files['zm'], fill_zeros=True)
        upwp_pr1 = self.__getFallbackVar__('upwp_pr1', self.ncdf_files['zm'], fill_zeros=True)
        upwp_pr2 = self.__getFallbackVar__('upwp_pr2', self.ncdf_files['zm'], fill_zeros=True)
        upwp_pr3 = self.__getFallbackVar__('upwp_pr3', self.ncdf_files['zm'], fill_zeros=True)
        upwp_pr4 = self.__getFallbackVar__('upwp_pr4', self.ncdf_files['zm'], fill_zeros=True)
        upwp_mfl = self.__getFallbackVar__('upwp_mfl', self.ncdf_files['zm'], fill_zeros=True)
        upwp_ma = self.__getFallbackVar__('upwp_ma', self.ncdf_files['zm'], fill_zeros=True)


        output_data = upwp_bt - (upwp_ma + upwp_ta + upwp_tp + upwp_ac + upwp_bp + upwp_pr1 + upwp_pr2 + upwp_pr3 + upwp_pr4 + upwp_dp1 + upwp_mfl + upwp_cl)
        output_line = Line(output_data, z_ncdf, line_format="", label="upwp_residual")

        return output_line
    
    def getVpwpResidual(self):
        '''


        vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_dp1 + vpwp_dp2 + vpwp_cl + vpwp_pd + vpwp_sf + vpwp_forcing)
        :return:
        '''
        z_ncdf = self.__getFallbackVar__('altitude', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_cl = self.__getFallbackVar__('vpwp_cl', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_tp = self.__getFallbackVar__('vpwp_tp', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_ac = self.__getFallbackVar__('vpwp_ac', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_bp = self.__getFallbackVar__('vpwp_bp', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_dp1 = self.__getFallbackVar__('vpwp_dp1', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_bt = self.__getFallbackVar__('vpwp_bt', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_ta = self.__getFallbackVar__('vpwp_ta', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_pr1 = self.__getFallbackVar__('vpwp_pr1', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_pr2 = self.__getFallbackVar__('vpwp_pr2', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_pr3 = self.__getFallbackVar__('vpwp_pr3', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_pr4 = self.__getFallbackVar__('vpwp_pr4', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_mfl = self.__getFallbackVar__('vpwp_mfl', self.ncdf_files['zm'], fill_zeros=True)
        vpwp_ma = self.__getFallbackVar__('vpwp_ma', self.ncdf_files['zm'], fill_zeros=True)


        output_data = vpwp_bt - (vpwp_ma + vpwp_ta + vpwp_tp + vpwp_ac + vpwp_bp + vpwp_pr1 + vpwp_pr2 + vpwp_pr3 + vpwp_pr4 + vpwp_dp1 + vpwp_mfl + vpwp_cl)
        output_line = Line(output_data, z_ncdf, line_format="", label="vpwp_residual")

        return output_line