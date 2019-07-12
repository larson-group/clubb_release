'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Line import Line


class VariableGroupBase(VariableGroup):
    '''
    This is a panel group used for testing the functionality of pyplotgen.
    It contains a set of common panels being used for representing the majority
    of panels.
    '''

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "base variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'thlm', 'sam_calc': self.getThlmSamLine},
            {'clubb_name': 'rtm', 'sam_calc': self.getRtmSamLine},
            {'clubb_name': 'wpthlp', 'sam_name': 'WPTHLP', 'fallback_func': self.getWpthlpFallback},
            {'clubb_name': 'wprtp', 'sam_name': 'WPRTP', 'fallback_func': self.getWprtpFallback},
            {'clubb_name': 'rcm', 'sam_name': "QCL", 'sam_conv_factor': 1/1000},
            {'clubb_name': 'wp3', 'sam_name': 'W3'},
            {'clubb_name': 'thlp2', 'sam_name': 'THLP2', 'fallback_func': self.getThlp2Fallback},
            {'clubb_name': 'rtp2', 'sam_name': 'RTP2', 'fallback_func': self.getRtp2Fallback},
            {'clubb_name': 'wm', 'sam_name': 'WOBS'},
            {'clubb_name': 'um', 'sam_name': 'U'},
            {'clubb_name': 'vm', 'sam_name': 'V'},
            {'clubb_name': 'wp3', 'sam_name': 'W3'},
            {'clubb_name': 'wp2', 'sam_name': 'W2'},
            {'clubb_name': 'wp2_vert_avg', 'sam_name': 'CWP', 'type': Panel.TYPE_TIMESERIES},
            {'clubb_name': 'cloud_frac', 'sam_name': 'CLD'},
            {'clubb_name': 'upwp', 'sam_name': 'UW'},
            {'clubb_name': 'vpwp', 'sam_name': 'VW'},
            {'clubb_name': 'up2', 'sam_name': 'U2'},
            {'clubb_name': 'tau_zm'},
            {'clubb_name': 'vp2', 'sam_name': 'V2'},
            {'clubb_name': 'Lscale'},
            {'clubb_name': 'wpthvp', 'sam_name': 'WPTHVP', 'fallback_func': self.getWpthvpFallback},
            {'clubb_name': 'rtpthlp', 'sam_name': 'RTPTHLP', 'fallback_func': self.getRtpthlpFallback},
            # {'clubb_name': 'rtp3', 'sam_name': 'RTP3'},
            {'clubb_name': 'radht', 'sam_name': 'RADQR', 'sam_conv_factor': 1/86400},
            {'clubb_name': 'Skw_zt', 'sam_calc': self.getSkwZtSamLine},
            # {'clubb_name': 'thlp3', 'sam_name': 'THLP3'},
            # {'clubb_name': 'rtpthvp', 'sam_name': 'RTPTHVP'},
            {'clubb_name': 'Skrt_zt', 'sam_calc': self.getSkrtZtSamLine},
            {'clubb_name': 'Skthl_zt', 'sam_calc': self.getSkthlZtSamLine},
            {'clubb_name': 'corr_w_chi_1'},
            {'clubb_name': 'corr_chi_eta_1'},
            {'clubb_name': 'rcp2', 'sam_name': 'QC2', 'sam_conv_factor': 1/10**6},
            # {'clubb_name': 'thlpthvp', 'sam_name': 'THLPTHVP'},

        ]
        super().__init__(ncdf_datasets, case, sam_file)

    def getThlmSamLine(self):
        '''
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sec_per_min = 60
        sam_start_time = self.start_time / sec_per_min
        sam_end_time = self.end_time / sec_per_min

        thetal_ncdf = NetCdfVariable('THETAL', self.sam_file, 1, start_time=sam_start_time,
                                     end_time=sam_end_time)
        thetal = thetal_ncdf.data

        theta_ncdf = NetCdfVariable('THETA', self.sam_file, 1, start_time=sam_start_time,
                                    end_time=sam_end_time)
        theta = theta_ncdf.data

        tabs_ncdf = NetCdfVariable('TABS', self.sam_file, 1, start_time=sam_start_time,
                                   end_time=sam_end_time)
        tabs = tabs_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time, fill_zeros=True)
        qi = qi_ncdf.data

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        thlm = thlm[self.z_sam_min_idx:self.z_sam_max_idx]

        thlm_line = Line(thlm, self.z_sam.data, line_format="k-", label="LES output")
        return thlm_line

    def getRtmSamLine(self):
        '''
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        qt_ncdf = NetCdfVariable('QT', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        qt = qt_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time, fill_zeros=True)
        qi = qi_ncdf.data

        rtm = (qt - qi) / 1000
        rtm = rtm[self.z_sam_min_idx:self.z_sam_max_idx]

        rtm_line = Line(rtm, self.z_sam.data, line_format="k-", label="LES output")
        return rtm_line


    def getSkwZtSamLine(self):
        '''
        Calculates Skw_zt values from sam output using
        the following equation
        WP3 ./ (WP2 + 1.6e-3).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        wp3_ncdf = NetCdfVariable('WP3', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        wp3 = wp3_ncdf.data

        wp2_ncdf = NetCdfVariable('WP2', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        wp2 = wp2_ncdf.data

        skw_zt = wp3 / (wp2 + 1.6e-3 ) ** 1.5
        skw_zt = skw_zt[self.z_sam_min_idx:self.z_sam_max_idx]

        skw_zt_line = Line(skw_zt, self.z_sam.data, line_format="k-", label="LES output")
        return skw_zt_line

    def getSkrtZtSamLine(self):
        '''
        Calculates Skrt_zt values from sam output using
        the following equation
         # RTP3 ./ (RTP2 + 4e-16).^1.5  
         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        rtp3_ncdf = NetCdfVariable('RTP3', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rtp3 = rtp3_ncdf.data

        rtp2_ncdf = NetCdfVariable('RTP2', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rtp2 = rtp2_ncdf.data

        skrtp_zt = rtp3 / (rtp2 + 4e-16) ** 1.5
        skrtp_zt = skrtp_zt[self.z_sam_min_idx:self.z_sam_max_idx]

        skrtp_zt_line = Line(skrtp_zt, self.z_sam.data, line_format="k-", label="LES output")
        return skrtp_zt_line
    
    def getSkthlZtSamLine(self):
        '''
        Calculates Skthl_zt values from sam output using
        the following equation
        THLP3 ./ (THLP2 + 4e-4).^1.5
         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        thlp3_ncdf = NetCdfVariable('THLP3', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        thlp3 = thlp3_ncdf.data

        thlp2_ncdf = NetCdfVariable('THLP2', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        thlp2 = thlp2_ncdf.data

        skthl_zt = thlp3 / (thlp2 + 4e-16) ** 1.5
        skthl_zt = skthl_zt[self.z_sam_min_idx:self.z_sam_max_idx]

        skthl_zt_line = Line(skthl_zt, self.z_sam.data, line_format="k-", label="LES output")
        return skthl_zt_line

    def getWpthlpFallback(self):
        '''
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHLP = (TLFLUX) ./ (RHO * 1004)
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        tlflux_ncdf = NetCdfVariable('TLFLUX', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        tlflux = tlflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data

        wpthlp = tlflux / (rho * 1004)
        wpthlp = wpthlp[self.z_sam_min_idx:self.z_sam_max_idx]

        wpthlp = Line(wpthlp, self.z_sam.data, line_format="k-", label="LES output")
        return wpthlp

    def getWprtpFallback(self):
        '''
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPRTP = (QTFLUX) ./ (RHO * 2.5104e+6)
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        qtflux_ncdf = NetCdfVariable('QTFLUX', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        qtflux = qtflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data

        wprtp = qtflux / (rho * 2.5104e+6)
        wprtp = wprtp[self.z_sam_min_idx:self.z_sam_max_idx]

        wprtp = Line(wprtp, self.z_sam.data, line_format="k-", label="LES output")
        return wprtp

    def getWpthvpFallback(self):
        '''
        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHVP =  (TVFLUX) ./ ( RHO * 1004)
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        tvflux_ncdf = NetCdfVariable('TVFLUX', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        tvflux = tvflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data

        wpthvp = tvflux / (rho * 1004)
        wpthvp = wpthvp[self.z_sam_min_idx:self.z_sam_max_idx]

        wpthvp = Line(wpthvp, self.z_sam.data, line_format="k-", label="LES output")
        return wpthvp

    def getThlp2Fallback(self):
        '''
        This gets called if THLP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = TL2
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        tl2_ncdf = NetCdfVariable('TL2', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        tl2 = tl2_ncdf.data
        tl2 = tl2[self.z_sam_min_idx:self.z_sam_max_idx]

        tl2_line = Line(tl2, self.z_sam.data, line_format="k-", label="LES output")
        return tl2_line

    def getRtpthlpFallback(self):
        '''
        This gets called if Rtpthlp isn't outputted in an nc file as a backup way of gathering the data for plotting.
        Rtpthlp = TQ
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        tq_ncdf = NetCdfVariable('TQ', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        tq2 = tq_ncdf.data
        tq2 = tq2[self.z_sam_min_idx:self.z_sam_max_idx]

        thlp2 = Line(tq2, self.z_sam.data, line_format="k-", label="LES output")
        return thlp2

    def getRtp2Fallback(self):
        '''
        This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = QT2 / 1e+6
        :return:
        '''
        sam_start_time = self.start_time / 60
        sam_end_time = self.end_time / 60

        qt2_ncdf = NetCdfVariable('QT2', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        qt2 = qt2_ncdf.data

        rtp2 = qt2 / 1e6
        rtp2 = rtp2[self.z_sam_min_idx:self.z_sam_max_idx]

        rtp2_line = Line(rtp2, self.z_sam.data, line_format="k-", label="LES output")
        return rtp2_line

