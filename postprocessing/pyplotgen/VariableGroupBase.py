'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from Line import Line
from Panel import Panel
from VariableGroup import VariableGroup


class VariableGroupBase(VariableGroup):
    '''
    This is a panel group used for testing the functionality of pyplotgen.
    It contains a set of common panels being used for representing the majority
    of panels.
    '''

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "base variables"
        self.variable_definitions = [
            {'aliases': ['thlm'], 'sam_calc': self.getThlmSamLine},
            {'aliases': ['rtm', 'qtm'],	 'sam_calc': self.getRtmSamLine},
            {'aliases': ['wpthlp', 'WPTHLP'], 'fallback_func': self.getWpthlpFallback},
            {'aliases': ['wprtp', 'WPRTP', 'wpqtp'], 'fallback_func': self.getWprtpFallback},
            {'aliases': ['cloud_frac', 'cf', 'CLD']},
            {'aliases': ['rcm', 'QCL', 'qcm'], 'sam_conv_factor': 1 / 1000},
            {'aliases': ['wp2', 'W2', 'WP2']},
            {'aliases': ['wp3', 'W3', 'WP3'],	 'sam_name': 'W3'},
            {'aliases': ['thlp2', 'THLP2'], 'fallback_func': self.getThlp2Fallback},
            {'aliases': ['rtp2', 'RTP2', 'qtp2'], 'fallback_func': self.getRtp2Fallback},
            {'aliases': ['rtpthlp', 'RTPTHLP', 'qtpthlp'], 'fallback_func': self.getRtpthlpFallback},
            {'aliases': ['rtp3', 'RTP3', 'qtp3'], 'fallback_func': self.getRtp3Fallback},
            {'aliases': ['thlp3', 'THLP3']},
            {'aliases': ['Skw_zt'],	 'sam_calc': self.getSkwZtSamLine, 'fill_zeros':True},	# TODO SAM output (arm97) is wrong # TODO coamps eqn wp3 ./ (wp2 + 1.6e-3).^1.5
            {'aliases': ['Skrt_zt'],	 'sam_calc': self.getSkrtZtSamLine, 'fill_zeros': True},	 # TODO coamps eqn qtp3 ./ (qtp2 + 4e-16).^1.5
            {'aliases': ['Skthl_zt'],	 'sam_calc': self.getSkthlZtSamLine, 'fill_zeros': True},	 # TODO coamps eqn thlp3 ./ (thlp2 + 4e-4).^1.5
            {'aliases': ['wm', 'WOBS', 'wlsm']},
            {'aliases': ['um', 'U']},
            {'aliases': ['vm', 'V']},
            {'aliases': ['upwp', 'UW']}, # TODO coamps eqn wpup + wpup_sgs
            {'aliases': ['vpwp', 'VW']}, # TODO coamps eqn wpvp + wpvp_sgs
            {'aliases': ['up2', 'U2']},
            {'aliases': ['vp2', 'V2']},
            {'aliases': ['rcp2', 'QC2', 'qcp2'], 'sam_conv_factor': 1 / 10 ** 6},
            {'aliases': ['wp2_vert_avg', 'CWP'], 'type': Panel.TYPE_TIMESERIES,	 'fill_zeros': True},
            {'aliases': ['tau_zm'], 'fill_zeros': True},
            {'aliases': ['Lscale'], 'fill_zeros': True},
            {'aliases': ['wpthvp', 'WPTHVP'], 'fallback_func': self.getWpthvpFallback},
            {'aliases': ['radht', 'RADQR'], 'sam_conv_factor': 1 / 86400},
            {'aliases': ['rtpthvp', 'RTPTHVP', 'qtpthvp']},
            {'aliases': ['corr_w_chi_1'], 'fill_zeros': True},
            {'aliases': ['corr_chi_eta_1'], 'fill_zeros': True},
            {'aliases': ['thlpthvp']}, # TODO SAM output missing (arm)

            # TODO SAM output for these variables
            {'aliases': ['rc_coef_zm .* wprcp'],	 'fallback_func': self.get_rc_coef_zm_X_wprcp_clubb_line,
                'title': 'Contribution of Cloud Water Flux to wpthvp',	 'axis_title': 'rc_coef_zm * wprcp [K m/s]'},	 # TODO coamps eqn wpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
            {'aliases': ['rc_coef_zm .* thlprcp'],	 'fallback_func': self.get_rc_coef_zm_X_thlprcp_clubb_line,
                'title': 'Contribution of Cloud Water Flux to thlprcp',	 'axis_title': 'rc_coef_zm * thlprcp [K^2]'},	 # TODO coamps eqn thlpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
            {'aliases': ['rc_coef_zm .* rtprcp'],	 'fallback_func': self.get_rc_coef_zm_X_rtprcp_clubb_line,
                'title': 'Contribution of Cloud Water Flux to rtprcp',	 'axis_title': 'rc_coef_zm * rtprcp [kg/kg K]'},	 # TODO coamp eqn qtpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
            {'aliases': ['lwp', 'CWP'],	 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1/1000}

            # TODO rc_coev * wp2rcp


            # TODO corr chi 2's
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)

    def getThlmSamLine(self):
        '''
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''

        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        thetal_ncdf = NetCdfVariable('THETAL', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        thetal_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        thetal = thetal_ncdf.data

        theta_ncdf = NetCdfVariable('THETA', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        theta_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        theta = theta_ncdf.data

        tabs_ncdf = NetCdfVariable('TABS', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        tabs_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        tabs = tabs_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        qi_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qi = qi_ncdf.data

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        thlm_line = Line(thlm, z_ncdf.data, line_format="k-", label="SAM-LES")
        return thlm_line

    def getRtmSamLine(self):
        '''
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        qt_ncdf = NetCdfVariable('QT', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        qt_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qt = qt_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        qi_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qi = qi_ncdf.data

        rtm = (qt - qi) / 1000

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        rtm_line = Line(rtm, z_ncdf.data, line_format="k-", label="SAM-LES")
        return rtm_line

    def getSkwZtSamLine(self):
        '''
        Calculates Skw_zt values from sam output using
        the following equation
        WP3 ./ (WP2 + 1.6e-3).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        wp3_ncdf = NetCdfVariable('WP3', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True) # TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        wp3_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        wp3 = wp3_ncdf.data

        wp2_ncdf = NetCdfVariable('WP2', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)# TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        wp2_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        wp2 = wp2_ncdf.data

        skw_zt = wp3 / (wp2 + 1.6e-3) ** 1.5

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        skw_zt_line = Line(skw_zt, z_ncdf.data, line_format="k-", label="SAM-LES")
        return skw_zt_line

    def getSkrtZtSamLine(self):
        '''
        Calculates Skrt_zt values from sam output using
        the following equation
         # RTP3 ./ (RTP2 + 4e-16).^1.5  
         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        rtp3_ncdf = NetCdfVariable('RTP3', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)# TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        rtp3_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rtp3 = rtp3_ncdf.data

        rtp2_ncdf = NetCdfVariable('RTP2', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)# TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        rtp2_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rtp2 = rtp2_ncdf.data

        skrtp_zt = rtp3 / (rtp2 + 4e-16) ** 1.5

        z_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        skrtp_zt_line = Line(skrtp_zt, z_ncdf.data, line_format="k-", label="SAM-LES")
        return skrtp_zt_line

    def getSkthlZtSamLine(self):
        '''
        Calculates Skthl_zt values from sam output using
        the following equation
        THLP3 ./ (THLP2 + 4e-4).^1.5
         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        thlp3_ncdf = NetCdfVariable('THLP3', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True) # TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        thlp3_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        thlp3 = thlp3_ncdf.data

        thlp2_ncdf = NetCdfVariable('THLP2', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True) # TODO fill zeros tempfix until pyplotgen can handle multiple aliases
        thlp2_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        thlp2 = thlp2_ncdf.data

        skthl_zt = thlp3 / (thlp2 + 4e-16) ** 1.5

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        skthl_zt_line = Line(skthl_zt, z_ncdf.data, line_format="k-", label="SAM-LES")
        return skthl_zt_line

    def getWpthlpFallback(self, dataset_override = None):
        '''
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHLP = (TLFLUX) ./ (RHO * 1004)
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        tlflux_ncdf = NetCdfVariable('TLFLUX', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        tlflux_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        tlflux = tlflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        rho_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rho = rho_ncdf.data

        wpthlp = tlflux / (rho * 1004)

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        wpthlp = Line(wpthlp, z_ncdf.data, line_format="k-", label="SAM-LES")
        return wpthlp

    def getWprtpFallback(self, dataset_override = None):
        '''
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPRTP = (QTFLUX) ./ (RHO * 2.5104e+6)
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        qtflux_ncdf = NetCdfVariable('QTFLUX', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        qtflux_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qtflux = qtflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        rho_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rho = rho_ncdf.data

        wprtp = qtflux / (rho * 2.5104e+6)

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        wprtp = Line(wprtp, z_ncdf.data, line_format="k-", label="SAM-LES")
        return wprtp

    def getWpthvpFallback(self, dataset_override = None):
        '''
        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHVP =  (TVFLUX) ./ ( RHO * 1004)
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        tvflux_ncdf = NetCdfVariable('TVFLUX', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        tvflux_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        tvflux = tvflux_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        rho_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rho = rho_ncdf.data

        wpthvp = tvflux / (rho * 1004)

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        wpthvp = Line(wpthvp, z_ncdf.data, line_format="k-", label="SAM-LES")
        return wpthvp

    def getThlp2Fallback(self, dataset_override = None):
        '''
        This gets called if THLP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = TL2
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        tl2_ncdf = NetCdfVariable('TL2', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        tl2_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        tl2 = tl2_ncdf.data

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        tl2_line = Line(tl2, z_ncdf.data, line_format="k-", label="SAM-LES")
        return tl2_line

    def getRtpthlpFallback(self, dataset_override = None):
        '''
        This gets called if Rtpthlp isn't outputted in an nc file as a backup way of gathering the data for plotting.
        Rtpthlp = TQ
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        tq_ncdf = NetCdfVariable('TQ', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        tq_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        tq2 = tq_ncdf.data

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        thlp2 = Line(tq2, z_ncdf.data, line_format="k-", label="SAM-LES")
        return thlp2

    def getRtp2Fallback(self, dataset_override = None):
        '''
        This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = QT2 / 1e+6
        :return:
        '''
        self.start_time = self.start_time
        self.end_time = self.end_time

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        qt2_ncdf = NetCdfVariable('QT2', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        qt2_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qt2 = qt2_ncdf.data

        rtp2 = qt2 / 1e6

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        rtp2_line = Line(rtp2, z_ncdf.data, line_format="k-", label="SAM-LES")
        return rtp2_line

    def getRtp3Fallback(self, dataset_override=None):
        '''
        Caclulates Rtp3 output
        rc_coef_zm .* rtprcp

        :return:
        '''
        rtp3 = None
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        if 'rc_coef_zm' in dataset.variables.keys() and 'rtprcp' in dataset.variables.keys():
            rc_coef_zm = self.__getFallbackVar__('rc_coef_zm', dataset)
            rtprcp = self.__getFallbackVar__('rtprcp', dataset)
            rtp3 = rc_coef_zm * (rtprcp)

        elif 'QCFLUX' in dataset.variables.keys():
            QCFLUX = self.__getFallbackVar__('QCFLUX', dataset)
            RHO = self.__getFallbackVar__('RHO', dataset)
            PRES = self.__getFallbackVar__('PRES', dataset)
            THETAV = self.__getFallbackVar__('THETAV', dataset)
            rtp3 = ((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)
        return rtp3

    def get_rc_coef_zm_X_wprcp_clubb_line(self, dataset_override=None):
        '''
        Calculates the Contribution of Cloud Water Flux
        to wpthvp using the equation
        rc_coef_zm .* wprcp
        :return: Line representing rc_coef_zm .* wprcp
        '''
        z_ncdf = NetCdfVariable('altitude', self.ncdf_files['zm'], 1)

        rc_coef_zm_ncdf = NetCdfVariable('rc_coef_zm', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        rc_coef_zm_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rc_coef_zm = rc_coef_zm_ncdf.data

        wprcp_ncdf = NetCdfVariable('wprcp', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        wprcp_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        wprcp = wprcp_ncdf.data

        output = rc_coef_zm * wprcp

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        output = Line(output, z_ncdf.data, line_format='b-', label='current clubb')
        return output

    def get_rc_coef_zm_X_wprcp_sam_line(self, dataset_override = None):
        '''
        Calculates the Contribution of Cloud Water Flux
        to wpthvp for SAM using the equation

        wpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :return:
        '''
        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        wpgcp_ncdf = NetCdfVariable('wpgcp', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        wpgcp_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        wpgcp = wpgcp_ncdf.data

        thvm_ncdf = NetCdfVariable('thvm', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        thvm_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        thvm = thvm_ncdf.data

        ex0_ncdf = NetCdfVariable('ex0', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        ex0_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        ex0 = thvm_ncdf.data

        output = wpgcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        output = Line(output, z_ncdf.data, line_format='k-', label='SAM-LES')
        return output

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_clubb_line(self, dataset_override = None):
        '''
        Calculates the Contribution of Cloud Water Flux
        to thlprcp using the equation
        rc_coef_zm * thlprcp
        :return: Line representing rc_coef_zm .* thlprcp
        '''
        z_ncdf = NetCdfVariable('altitude', self.ncdf_files['zm'], 1)

        rc_coef_zm_ncdf = NetCdfVariable('rc_coef_zm', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        rc_coef_zm_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rc_coef_zm = rc_coef_zm_ncdf.data

        thlprcp_ncdf = NetCdfVariable('thlprcp', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        thlprcp_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        thlprcp = thlprcp_ncdf.data

        output = rc_coef_zm * thlprcp

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        output = Line(output, z_ncdf.data, line_format='b-', label='current clubb')
        return output
    
    def get_rc_coef_zm_X_rtprcp_clubb_line(self, dataset_override = None):
        '''
        Calculates the Contribution of Cloud Water Flux
        to rtprcp using the equation
        rc_coef_zm * rtprcp
        :return: Line representing rc_coef_zm .* rtprcp
        '''
        z_ncdf = NetCdfVariable('altitude', self.ncdf_files['zm'], 1)

        rc_coef_zm_ncdf = NetCdfVariable('rc_coef_zm', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        rc_coef_zm_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rc_coef_zm = rc_coef_zm_ncdf.data

        rtprcp_ncdf = NetCdfVariable('rtprcp', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        rtprcp_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rtprcp = rtprcp_ncdf.data

        output = rc_coef_zm * rtprcp

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        output = Line(output, z_ncdf.data, line_format='b-', label='current clubb')
        return output