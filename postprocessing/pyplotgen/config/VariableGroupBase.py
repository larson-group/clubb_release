"""
:author: Nicolas Strike
:date: Mid 2019
"""
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupBase(VariableGroup):
    """
    This is a panel group used for testing the functionality of pyplotgen.
    It contains a set of common panels being used for representing the majority
    of panels.
    """

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_dataset=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "base variables"
        self.variable_definitions = [
            {'aliases': ['thlm'], 'sam_calc': self.getThlmSamCalc},
            {'aliases': ['rtm', 'qtm'],	 'sam_calc': self.getRtmSamCalc},
            {'aliases': ['wpthlp', 'WPTHLP'], 'fallback_func': self.getWpthlpFallback},
            {'aliases': ['wprtp', 'WPRTP', 'wpqtp'], 'fallback_func': self.getWprtpFallback},
            {'aliases': ['cloud_frac', 'cf', 'CLD']},
            {'aliases': ['rcm', 'QCL', 'qcm'], 'sam_conv_factor': 1 / 1000},
            {'aliases': ['wp2', 'W2', 'WP2']},
            {'aliases': ['wp3', 'W3', 'WP3'],	 'sam_name': 'W3'},
            {'aliases': ['thlp2', 'THLP2', 'TL2']},
            {'aliases': ['rtp2', 'RTP2', 'qtp2'], 'fallback_func': self.getRtp2Fallback},
            {'aliases': ['rtpthlp', 'RTPTHLP', 'qtpthlp', 'TQ']},
            {'aliases': ['rtp3', 'RTP3', 'qtp3'], 'fallback_func': self.getRtp3Fallback},
            {'aliases': ['thlp3', 'THLP3']},
            {'aliases': ['Skw_zt'],	 'sam_calc': self.getSkwZtLesCalc, 'coamps_calc': self.getSkwZtLesCalc, 'fill_zeros':True},
            {'aliases': ['Skrt_zt'],	 'sam_calc': self.getSkrtZtLesCalc, 'coamps_calc': self.getSkrtZtLesCalc, 'fill_zeros': True},
            {'aliases': ['Skthl_zt'],	 'sam_calc': self.getSkthlZtLesCalc, 'coamps_calc': self.getSkthlZtLesCalc, 'fill_zeros': True},
            {'aliases': ['wm', 'WOBS', 'wlsm']},
            {'aliases': ['um', 'U']},
            {'aliases': ['vm', 'V']},
            {'aliases': ['upwp', 'UW'], 'coamps_calc': self.getUwCoampsData},
            {'aliases': ['vpwp', 'VW'], 'coamps_calc': self.getVwCoampsData},
            {'aliases': ['up2', 'U2']},
            {'aliases': ['vp2', 'V2']},
            {'aliases': ['rcp2', 'QC2', 'qcp2'], 'sam_conv_factor': 1 / 10 ** 6},
            {'aliases': ['lwp', 'CWP'],	 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1/1000},
            {'aliases': ['wp2_vert_avg', 'W2_VERT_AVG'], 'type': Panel.TYPE_TIMESERIES,	 'fill_zeros': True},
            {'aliases': ['tau_zm'], 'fill_zeros': True},
            {'aliases': ['Lscale'], 'fill_zeros': True},
            {'aliases': ['wpthvp', 'WPTHVP'], 'fallback_func': self.getWpthvpFallback},
            {'aliases': ['radht', 'RADQR'], 'sam_conv_factor': 1 / 86400},
            {'aliases': ['rtpthvp', 'RTPTHVP', 'qtpthvp']},
            {'aliases': ['corr_w_chi_1'], 'fill_zeros': True},
            {'aliases': ['corr_chi_eta_1'], 'fill_zeros': True},
            {'aliases': ['thlpthvp', 'THLPTHVP']},

            # TODO SAM output for these variables
            # TODO validate coamps output
            # TODO Fix output for these vars in some cases
            {'aliases': ['rc_coef_zm * wprcp'],
             'fallback_func': self.get_rc_coef_zm_X_wprcp_clubb_line,
             'sam_calc': self.get_rc_coef_zm_X_wprcp_sam_calc,
             'coamps_calc': self.get_rc_coef_zm_X_wprcp_coamps_calc,
             'title': 'Contribution of Cloud Water Flux to wpthvp',
             'axis_title': 'rc_coef_zm * wprcp [K m/s]'},

            {'aliases': ['rc_coef_zm * thlprcp'],
             'coamps_calc': self.get_rc_coef_zm_X_thlprcp_coamps_calc,
             'fallback_func': self.get_rc_coef_zm_X_thlprcp_clubb_fallback,
             'title': 'Contribution of Cloud Water Flux to thlprcp',
             'axis_title': 'rc_coef_zm * thlprcp [K^2]'},

            {'aliases': ['rc_coef_zm * rtprcp'],
             'coamps_calc': self.get_rc_coef_zm_X_rtprcp_coamps_calc,
             'fallback_func': self.get_rc_coef_zm_X_rtprcp_clubb_fallback,
             'title': 'Contribution of Cloud Water Flux to rtprcp',
             'axis_title': 'rc_coef_zm * rtprcp [kg/kg K]'}

            # TODO rc_coev * wp2rcp


            # TODO corr chi 2's
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_dataset = e3sm_dataset)

    def getThlmSamCalc(self, dataset_override = None):
        """
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        # z,z, dataset = self.getVarForCalculations('z', self.sam_file)
        thetal,z, dataset = self.getVarForCalculations('THETAL', self.sam_file)
        theta,z, dataset = self.getVarForCalculations('THETA', self.sam_file)
        tabs,z, dataset = self.getVarForCalculations('TABS', self.sam_file)
        qi,z, dataset = self.getVarForCalculations('QI', self.sam_file, fill_zeros=True)

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        return thlm, z

    # def getThlmE3smCalc(self, dataset_override = None):
    #     """
    #     Calculates thlm values from sam output using
    #     the following equation
    #     (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
    #     :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
    #     """
    #     # z,z, dataset = self.getVarForCalculations('Z3', self.e3sm_dataset)
    #     thlm,z, dataset = self.getVarForCalculations('THETAL', self.e3sm_dataset)
    #
    #     thlm = thlm - 650.0
    #     return thlm, z


    def getRtmSamCalc(self, dataset_override = None):
        """
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        # z,z, dataset = self.getVarForCalculations('z', self.sam_file)
        qt,z, dataset = self.getVarForCalculations('QT', self.sam_file)
        qi,z, dataset = self.getVarForCalculations('QI', self.sam_file, fill_zeros=True)

        rtm = (qt - qi) / 1000
        return rtm,z

    def getSkwZtLesCalc(self, dataset_override = None):
        """
        Calculates Skw_zt values from sam output using
        the following equation
        WP3 ./ (WP2 + 1.6e-3).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.sam_file is not None:
            dataset = self.sam_file

        if self.coamps_file is not None:
            dataset = self.coamps_file['sm']

        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wp3,z, dataset = self.getVarForCalculations(['WP3', 'W3', 'wp3'], dataset)
        wp2,z, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        skw_zt = wp3 / (wp2 + 1.6e-3) ** 1.5

        return skw_zt, z

    def getSkrtZtLesCalc(self, dataset_override = None):
        """
        Calculates Skrt_zt values from sam output using
        the following equation
         sam eqn RTP3 ./ (RTP2 + 4e-16).^1.5
         coamps eqn qtp3 ./ (qtp2 + 4e-16).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.sam_file is not None:
            dataset = self.sam_file

        if self.coamps_file is not None:
            dataset = self.coamps_file['sm']

        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rtp3,z, dataset = self.getVarForCalculations(['RTP3', 'qtp3'], dataset)
        rtp2,z, dataset = self.getVarForCalculations(['RTP2', 'qtp2'], dataset)
        skrtp_zt = rtp3 / (rtp2 + 4e-16) ** 1.5

        return skrtp_zt, z

    def getSkthlZtLesCalc(self, dataset_override = None):
        """
        Calculates Skthl_zt values from sam output using
        the following equation
        sam THLP3 ./ (THLP2 + 4e-4).^1.5
        coamps eqn thlp3 ./ (thlp2 + 4e-4).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.sam_file is not None:
            dataset = self.sam_file

        if self.coamps_file is not None:
            dataset = self.coamps_file['sm']

        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        thlp3,z, dataset = self.getVarForCalculations(['THLP3', 'thlp3'], dataset)
        thlp2,z, dataset = self.getVarForCalculations(['THLP2', 'thlp2'], dataset)

        skthl_zt = thlp3 / (thlp2 + 4e-16) ** 1.5
        return skthl_zt, z

    def getWpthlpFallback(self, dataset_override = None):
        """
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHLP = (TLFLUX) ./ (RHO * 1004)
        :return:
        """
        tlflux,z, dataset = self.getVarForCalculations(['TLFLUX'], self.sam_file)
        rho,z, dataset = self.getVarForCalculations(['RHO'], self.sam_file)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        wpthlp = tlflux / (rho * 1004)

        return wpthlp, z

    def getWprtpFallback(self, dataset_override = None):
        """
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPRTP = (QTFLUX) ./ (RHO * 2.5104e+6)
        :return:
        """
        qtflux,z, dataset = self.getVarForCalculations(['QTFLUX'], self.sam_file)
        rho,z, dataset = self.getVarForCalculations(['RHO'], self.sam_file)
        wprtp = qtflux / (rho * 2.5104e+6)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)
        return wprtp, z

    def getWpthvpFallback(self, dataset_override = None):
        """
        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHVP =  (TVFLUX) ./ ( RHO * 1004)
        :return:
        """
        tvflux,z, dataset = self.getVarForCalculations(['TVFLUX'], self.sam_file)
        rho,z, dataset = self.getVarForCalculations(['RHO'], self.sam_file)
        wpthvp = tvflux / (rho * 1004)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)
        return wpthvp, z

    def getRtp2Fallback(self, dataset_override = None):
        """
        This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = QT2 / 1e+6
        :return:
        """
        qt2,z, dataset = self.getVarForCalculations(['QT2'], self.sam_file)
        rtp2 = qt2 / 1e6
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)
        return rtp2, z

    def getRtp3Fallback(self, dataset_override=None):
        """
        Caclulates Rtp3 output
        rc_coef_zm .* rtprcp

        :return:
        """
        rtp3 = None
        z = None
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        for dataset in dataset.values():
            if 'rc_coef_zm' in dataset.variables.keys() and 'rtprcp' in dataset.variables.keys():
                rc_coef_zm,z, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
                rtprcp,z, dataset = self.getVarForCalculations('rtprcp', dataset)
                rtp3 = rc_coef_zm * (rtprcp)

            elif 'QCFLUX' in dataset.variables.keys():
                QCFLUX,z, dataset = self.getVarForCalculations('QCFLUX', dataset)
                RHO,z, dataset = self.getVarForCalculations('RHO', dataset)
                PRES,z, dataset = self.getVarForCalculations('PRES', dataset)
                THETAV,z, dataset = self.getVarForCalculations('THETAV', dataset)
                rtp3 = ((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)
            # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        return rtp3, z

    def get_rc_coef_zm_X_wprcp_clubb_line(self, dataset_override=None):
        """
        Calculates the Contribution of Cloud Water Flux
        to wpthvp using the equation
        rc_coef_zm .* wprcp
        :return: Line representing rc_coef_zm .* wprcp
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.ncdf_files['zm']
        rc_coef_zm,z, dataset = self.getVarForCalculations('rc_coef_zm', dataset, fill_zeros=True)
        wprcp,z, dataset = self.getVarForCalculations('wprcp', dataset, fill_zeros=True)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        output = rc_coef_zm * wprcp
        return output, z

    def get_rc_coef_zm_X_wprcp_sam_calc(self, dataset_override = None):
        """
        Calculates the Contribution of Cloud Water Flux
        to wpthvp for SAM using the equation

        sam eqn WPRCP * (2.5e6 / (1004.67*((PRES / 1000)^(287.04/1004.67))) - 1.61*THETAV)
        :return:
        """

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        WPRCP,z, dataset = self.getVarForCalculations('WPRCP', dataset, fill_zeros=True)
        PRES,z, dataset = self.getVarForCalculations('PRES', dataset, fill_zeros=True)
        THETAV,z, dataset = self.getVarForCalculations('THETAV', dataset, fill_zeros=True)

        output = WPRCP * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)
        return output, z

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_clubb_fallback(self, dataset_override = None):
        """
        Calculates the Contribution of Cloud Water Flux
        to thlprcp using the equation
        rc_coef_zm * thlprcp
        :return: Line representing rc_coef_zm .* thlprcp
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.ncdf_files['zm']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rc_coef_zm,z, dataset = self.getVarForCalculations('rc_coef_zm', dataset, fill_zeros=True)
        thlprcp,z, dataset = self.getVarForCalculations('thlprcp', dataset, fill_zeros=True)

        output = rc_coef_zm * thlprcp
        return output, z

    def get_rc_coef_zm_X_rtprcp_clubb_fallback(self, dataset_override = None):
        """
        Calculates the Contribution of Cloud Water Flux
        to rtprcp using the equation
        rc_coef_zm * rtprcp
        :return: Line representing rc_coef_zm .* rtprcp
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.ncdf_files['zm']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rc_coef_zm,z, dataset = self.getVarForCalculations('rc_coef_zm', dataset, fill_zeros=True)
        rtprcp,z, dataset = self.getVarForCalculations('rtprcp', dataset, fill_zeros=True)

        output = rc_coef_zm * rtprcp
        return output, z

    def getUwCoampsData(self, dataset_override = None):
        """
        coamps eqn upwp = wpup + wpup_sgs

         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_file['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpup,z, dataset = self.getVarForCalculations('wpup', dataset)
        wpup_sgs,z, dataset = self.getVarForCalculations('wpup_sgs', dataset)

        upwp = wpup + wpup_sgs
        return upwp, z

    def getVwCoampsData(self, dataset_override = None):
        """
        coamps eqn vpwp = wpvp + wpvp_sgs

         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_file['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpvp,z, dataset = self.getVarForCalculations('wpvp', dataset)
        wpvp_sgs,z, dataset = self.getVarForCalculations('wpvp_sgs', dataset)

        vpwp = wpvp + wpvp_sgs
        return vpwp, z

    def get_rc_coef_zm_X_wprcp_coamps_calc(self, dataset_override = None):
        """coamps eqn thlpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        coamps eqn wpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_file['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpqcp,z, dataset = self.getVarForCalculations('wqpcp', dataset, fill_zeros=True)
        ex0,z, dataset = self.getVarForCalculations('ex0', dataset, fill_zeros=True)
        thvm,z, dataset = self.getVarForCalculations('thvm', dataset, fill_zeros=True)

        output = wpqcp * (2.5e6 / (1004.67 * ex0) - 1.61*thvm)
        return output, z

    def get_rc_coef_zm_X_thlprcp_coamps_calc(self, dataset_override = None):
        """
        coamps eqn thlpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_file['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        thlpqcp,z, dataset = self.getVarForCalculations('thlpqcp', dataset, fill_zeros=True)
        ex0,z, dataset = self.getVarForCalculations('ex0', dataset, fill_zeros=True)
        thvm,z, dataset = self.getVarForCalculations('thvm', dataset, fill_zeros=True)

        output = thlpqcp * (2.5e6 / (1004.67*ex0) - 1.61*thvm)
        return output, z

    def get_rc_coef_zm_X_rtprcp_coamps_calc(self, dataset_override = None):
        """
        coamp eqn qtpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_file['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        qtpqcp,z, dataset = self.getVarForCalculations('qtpqcp', dataset, fill_zeros=True)
        ex0,z, dataset = self.getVarForCalculations('ex0', dataset, fill_zeros=True)
        thvm,z, dataset = self.getVarForCalculations('thvm', dataset, fill_zeros=True)

        output = qtpqcp * (2.5e6 / (1004.67*ex0) - 1.61*thvm)
        return output, z

