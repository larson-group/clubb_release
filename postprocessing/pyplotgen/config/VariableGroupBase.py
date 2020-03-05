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

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None, sam_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "base variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['thlm'],
                'sam': [],
                'coamps': ['thlm'],
                'r408': ['thlm'],
                'hoc': ['thlm'],
                'e3sm': ['thlm']
            },
                'sam_calc': self.getThlmSamCalc, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['rtm'],
                'sam': [],
                'coamps': ['qtm'],
                'r408': ['rtm'],
                'hoc': ['rtm'],
                'e3sm': ['rtm']
            },
                'sam_calc': self.getRtmSamCalc, 'sci_scale': -3},
            {'var_names': {
                'clubb': ['wpthlp'],
                'sam': ['WPTHLP'],
                'coamps': ['wpthlp'],
                'r408': ['wpthlp'],
                'hoc': ['wpthlp'],
                'e3sm': ['wpthlp']
            },
                'fallback_func': self.getWpthlpFallback, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['wprtp'],
                'sam': ['WPRTP'],
                'coamps': ['wpqtp'],
                'r408': ['wprtp'],
                'hoc': ['wprtp'],
                'e3sm': ['wprtp']
            },
                'fallback_func': self.getWprtpFallback, 'sci_scale': -4},
            {'var_names': {
                'clubb': ['cloud_frac'],
                'sam': ['CLD'],
                'coamps': ['cf'],
                'r408': ['cloud_frac', 'cf'],
                'hoc': ['cloud_frac', 'cf'],
                'e3sm': ['cloud_frac']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['rcm'],
                'sam': ['QCL'],
                'coamps': ['qcm'],
                'r408': ['rcm'],
                'hoc': ['rcm'],
                'e3sm': ['rcm']
            },
                'sam_conv_factor': 1 / 1000, 'sci_scale': -5},
            {'var_names': {
                'clubb': ['wp2', 'W2'],
                'sam': ['W2', 'WP2'],
                'coamps': ['wp2', 'W2'],
                'r408': ['wp2'],
                'hoc': ['wp2'],
                'e3sm': ['wp2']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['wp3'],
                'sam': ['wp3', 'W3', 'WP3'],
                'coamps': ['wp3', 'W3', 'WP3'],
                'r408': ['wp3'],
                'hoc': ['wp3'],
                'e3sm': ['wp3']
            },
                'sam_name': 'W3', 'sci_scale': 0, 'axis_title': "wp3"},
            {'var_names': {
                'clubb': ['thlp2'],
                'sam': ['THLP2', 'TL2'],
                'coamps': ['thlp2'],
                'r408': ['thlp2'],
                'hoc': ['thlp2'],
                'e3sm': ['thlp2']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['rtp2'],
                'sam': ['RTP2'],
                'coamps': ['qtp2'],
                'r408': ['rtp2'],
                'hoc': ['rtp2'],
                'e3sm': ['rtp2']
            },
                'fallback_func': self.getRtp2Fallback, 'sci_scale': -7},
            {'var_names': {
                'clubb': ['rtpthlp'],
                'sam': ['RTPTHLP', 'TQ'],
                'coamps': ['qtpthlp'],
                'r408': ['rtpthlp'],
                'hoc': ['rtpthlp'],
                'e3sm': ['rtpthlp']
            }, 'sci_scale': -4},

            {'var_names': {
                'clubb': ['rtp3'],
                'sam': ['RTP3'],
                'coamps': ['qtp3'],
                'r408': ['rtp3'],
                'hoc': ['rtp3'],
                'e3sm': ['rtp3']
            },
                'fallback_func': self.getRtp3Fallback, 'sci_scale': -9},
            {'var_names': {
                'clubb': ['thlp3'],
                'sam': ['THLP3'],
                'coamps': ['thlp3'],
                'r408': ['thlp3'],
                'hoc': ['thlp3'],
                'e3sm': ['thlp3']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['Skw_zt'],
                'sam': ['Skw_zt'],
                'coamps': ['Skw_zt'],
                'r408': ['Skw_zt'],
                'hoc': ['Skw_zt'],
                'e3sm': ['Skw_zt']
            },
                'sam_calc': self.getSkwZtLesCalc, 'coamps_calc': self.getSkwZtLesCalc, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['Skrt_zt'],
                'sam': ['Skrt_zt'],
                'coamps': ['Skrt_zt'],
                'r408': ['Skrt_zt'],
                'hoc': ['Skrt_zt'],
                'e3sm': ['Skrt_zt']
            },
                'sam_calc': self.getSkrtZtLesCalc, 'coamps_calc': self.getSkrtZtLesCalc, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['Skthl_zt'],
                'sam': ['Skthl_zt'],
                'coamps': ['Skthl_zt'],
                'r408': ['Skthl_zt'],
                'hoc': ['Skthl_zt'],
                'e3sm': ['Skthl_zt']
            },
                'sam_calc': self.getSkthlZtLesCalc, 'coamps_calc': self.getSkthlZtLesCalc, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['wm', 'wlsm'],
                'sam': ['wm', 'WOBS'],
                'coamps': ['wlsm'],
                'r408': ['wm'],
                'hoc': ['wm'],
                'e3sm': ['wm']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['um'],
                'sam': ['U'],
                'coamps': ['um'],
                'r408': ['um'],
                'hoc': ['um'],
                'e3sm': ['um']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['vm'],
                'sam': ['V'],
                'coamps': ['vm'],
                'r408': ['vm'],
                'hoc': ['vm'],
                'e3sm': ['vm']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['upwp'],
                'sam': ['UW'],
                'coamps': ['upwp'],
                'r408': ['upwp'],
                'hoc': ['upwp'],
                'e3sm': ['upwp']
            },
                'coamps_calc': self.getUwCoampsData, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['vpwp',],
                'sam': ['VW'],
                'coamps': ['vpwp'],
                'r408': ['vpwp'],
                'hoc': ['vpwp'],
                'e3sm': ['vpwp']
            },
                'coamps_calc': self.getVwCoampsData, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['up2'],
                'sam': ['U2'],
                'coamps': ['up2'],
                'r408': ['up2'],
                'hoc': ['up2'],
                'e3sm': ['up2']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['vp2'],
                'sam': ['V2'],
                'coamps': ['vp2'],
                'r408': ['vp2'],
                'hoc': ['vp2'],
                'e3sm': ['vp2']
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['rcp2'],
                'sam': ['QC2'],
                'coamps': ['qcp2'],
                'r408': ['rcp2'],
                'hoc': ['rcp2'],
                'e3sm': ['rcp2']
            },
                'sam_conv_factor': 1 / 10 ** 6, 'sci_scale': -8},
            {'var_names': {
                'clubb': ['lwp'],
                'sam': ['CWP'],
                'coamps': ['lwp'],
                'r408': ['lwp'],
                'hoc': ['lwp'],
                'e3sm': ['lwp']
            },
                'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'var_names': {
                'clubb': ['wp2_vert_avg'],
                'sam': ['W2_VERT_AVG'],
                'coamps': ['wp2_vert_avg'],
                'r408': ['wp2_vert_avg'],
                'hoc': ['wp2_vert_avg'],
                'e3sm': ['wp2_vert_avg']
            },
                'type': Panel.TYPE_TIMESERIES, },
            {'var_names': {
                'clubb': ['tau_zm'],
                'sam': ['tau_zm'],
                'coamps': ['tau_zm'],
                'r408': ['tau_zm'],
                'hoc': ['tau_zm'],
                'e3sm': ['tau_zm']
            },
                'sci_scale': 0},
            {'var_names': {
                'clubb': ['Lscale'],
                'sam': ['Lscale'],
                'coamps': ['Lscale'],
                'r408': ['Lscale'],
                'hoc': ['Lscale'],
                'e3sm': ['Lscale']
            },
                'sci_scale': 0},
            {'var_names': {
                'clubb': ['wpthvp'],
                'sam': ['WPTHVP'],
                'coamps': ['wpthvp'],
                'r408': ['wpthvp'],
                'hoc': ['wpthvp'],
                'e3sm': ['wpthvp']
            },
                'fallback_func': self.getWpthvpFallback},
            {'var_names': {
                'clubb': ['radht'],
                'sam': ['RADQR'],
                'coamps': ['radht'],
                'r408': ['radht'],
                'hoc': ['radht'],
                'e3sm': ['radht']
            },
                'sam_conv_factor': 1 / 86400},
            {'var_names': {
                'clubb': ['rtpthvp'],
                'sam': ['RTPTHVP'],
                'coamps': ['qtpthvp'],
                'r408': ['rtpthvp'],
                'hoc': ['rtpthvp'],
                'e3sm': ['rtpthvp']
            }, 'sci_scale': -5},

            {'var_names': {
                'clubb': ['corr_w_chi_1'],
                'sam': ['corr_w_chi_1'],
                'coamps': ['corr_w_chi_1'],
                'r408': ['corr_w_chi_1'],
                'hoc': ['corr_w_chi_1'],
                'e3sm': ['corr_w_chi_1']
            },
                'sci_scale': 0},
            {'var_names': {
                'clubb': ['corr_chi_eta_1'],
                'sam': ['corr_chi_eta_1'],
                'coamps': ['corr_chi_eta_1'],
                'r408': ['corr_chi_eta_1'],
                'hoc': ['corr_chi_eta_1'],
                'e3sm': ['corr_chi_eta_1']
            },
                'sci_scale': 0},
            {'var_names': {
                'clubb': ['thlpthvp'],
                'sam': ['thlpthvp', 'THLPTHVP'],
                'coamps': ['thlpthvp'],
                'r408': ['thlpthvp'],
                'hoc': ['thlpthvp'],
                'e3sm': ['thlpthvp']
            }},

            # TODO SAM output for these variables
            # TODO validate coamps output
            # TODO Fix output for these vars in some cases
            {'var_names': {
                'clubb': ['rc_coef_zm * wprcp'],
                'sam': ['rc_coef_zm * wprcp'],
                'coamps': ['rc_coef_zm * wprcp'],
                'r408': ['rc_coef_zm * wprcp'],
                'hoc': ['rc_coef_zm * wprcp'],
                'e3sm': ['rc_coef_zm * wprcp']
            },

                'fallback_func': self.get_rc_coef_zm_X_wprcp_clubb_line,
                'sam_calc': self.get_rc_coef_zm_X_wprcp_sam_calc,
                'coamps_calc': self.get_rc_coef_zm_X_wprcp_coamps_calc,
                'title': 'Contribution of Cloud Water Flux to wpthvp',
                'axis_title': 'rc_coef_zm * wprcp [K m/s]',
                'sci_scale': 0},

            {'var_names': {
                'clubb': ['rc_coef_zm * thlprcp'],
                'sam': ['rc_coef_zm * thlprcp'],
                'coamps': ['rc_coef_zm * thlprcp'],
                'r408': ['rc_coef_zm * thlprcp'],
                'hoc': ['rc_coef_zm * thlprcp'],
                'e3sm': ['rc_coef_zm * thlprcp']
            },

                'coamps_calc': self.get_rc_coef_zm_X_thlprcp_coamps_calc,
                'fallback_func': self.get_rc_coef_zm_X_thlprcp_clubb_fallback,
                'title': 'Contribution of Cloud Water Flux to thlprcp',
                'axis_title': 'rc_coef_zm * thlprcp [K^2]',
                'sci_scale': 0},

            {'var_names': {
                'clubb': ['rc_coef_zm * rtprcp'],
                'sam': ['rc_coef_zm * rtprcp'],
                'coamps': ['rc_coef_zm * rtprcp'],
                'r408': ['rc_coef_zm * rtprcp'],
                'hoc': ['rc_coef_zm * rtprcp'],
                'e3sm': ['rc_coef_zm * rtprcp']
            },

                'coamps_calc': self.get_rc_coef_zm_X_rtprcp_coamps_calc,
                'fallback_func': self.get_rc_coef_zm_X_rtprcp_clubb_fallback,
                'title': 'Contribution of Cloud Water Flux to rtprcp',
                'axis_title': 'rc_coef_zm * rtprcp [kg/kg K]',
                'sci_scale': 0}

            # TODO rc_coev * wp2rcp

            # TODO corr chi 2's
        ]
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset, r408_dataset=r408_dataset,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets)

    def getThlmSamCalc(self, dataset_override=None):
        """
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        # z,z, dataset = self.getVarForCalculations('z', self.les_dataset)
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        thetal, z, dataset = self.getVarForCalculations('THETAL', dataset)
        theta, z, dataset = self.getVarForCalculations('THETA', dataset)
        tabs, z, dataset = self.getVarForCalculations('TABS', dataset)
        qi, z, dataset = self.getVarForCalculations('QI', dataset)

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        return thlm, z

    # def getThlmE3smCalc(self, dataset_override = None):
    #     """
    #     Calculates thlm values from sam output using
    #     the following equation
    #     (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
    #     :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
    #     """
    #     # z,z, dataset = self.getVarForCalculations('Z3', self.e3sm_datasets)
    #     thlm,z, dataset = self.getVarForCalculations('THETAL', self.e3sm_datasets)
    #
    #     thlm = thlm - 650.0
    #     return thlm, z

    def getRtmSamCalc(self, dataset_override=None):
        """
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        # z,z, dataset = self.getVarForCalculations('z', self.les_dataset)
        qt, z, dataset = self.getVarForCalculations('QT', dataset)
        qi, z, dataset = self.getVarForCalculations('QI', dataset)

        rtm = (qt - qi) / 1000
        return rtm, z

    def getSkwZtLesCalc(self, dataset_override=None):
        """
        Calculates Skw_zt values from sam output using
        the following equation
        WP3 ./ (WP2 + 1.6e-3).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if self.coamps_dataset is not None:
            dataset = self.coamps_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override

            # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wp3, z, dataset = self.getVarForCalculations(['WP3', 'W3', 'wp3'], dataset)
        wp2, z, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        skw_zt = wp3 / (wp2 + 1.6e-3) ** 1.5

        return skw_zt, z

    def getSkrtZtLesCalc(self, dataset_override=None):
        """
        Calculates Skrt_zt values from sam output using
        the following equation
         sam eqn RTP3 ./ (RTP2 + 4e-16).^1.5
         coamps eqn qtp3 ./ (qtp2 + 4e-16).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset

        if self.coamps_dataset is not None:
            dataset = self.coamps_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override
            # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rtp3, z, dataset = self.getVarForCalculations(['RTP3', 'qtp3'], dataset)
        rtp2, z, dataset = self.getVarForCalculations(['RTP2', 'qtp2'], dataset)
        skrtp_zt = rtp3 / (rtp2 + 4e-16) ** 1.5

        return skrtp_zt, z

    def getSkthlZtLesCalc(self, dataset_override=None):
        """
        Calculates Skthl_zt values from sam output using
        the following equation
        sam THLP3 ./ (THLP2 + 4e-4).^1.5
        coamps eqn thlp3 ./ (thlp2 + 4e-4).^1.5
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset

        if self.coamps_dataset is not None:
            dataset = self.coamps_dataset['sm']

        if dataset_override is not None:
            dataset = dataset_override

            # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        thlp3, z, dataset = self.getVarForCalculations(['THLP3', 'thlp3'], dataset)
        thlp2, z, dataset = self.getVarForCalculations(['THLP2', 'thlp2'], dataset)

        skthl_zt = thlp3 / (thlp2 + 4e-16) ** 1.5
        return skthl_zt, z

    def getWpthlpFallback(self, dataset_override=None):
        """
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHLP = (TLFLUX) ./ (RHO * 1004)
        :return:
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        tlflux, z, dataset = self.getVarForCalculations(['TLFLUX'], dataset)
        rho, z, dataset = self.getVarForCalculations(['RHO'], dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)

        wpthlp = tlflux / (rho * 1004)

        return wpthlp, z

    def getWprtpFallback(self, dataset_override=None):
        """
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPRTP = (QTFLUX) ./ (RHO * 2.5104e+6)
        :return:
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qtflux, z, dataset = self.getVarForCalculations(['QTFLUX'], dataset)
        rho, z, dataset = self.getVarForCalculations(['RHO'], dataset)
        wprtp = qtflux / (rho * 2.5104e+6)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)
        return wprtp, z

    def getWpthvpFallback(self, dataset_override=None):
        """
        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHVP =  (TVFLUX) ./ ( RHO * 1004)
        :return:
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        tvflux, z, dataset = self.getVarForCalculations(['TVFLUX'], dataset)
        rho, z, dataset = self.getVarForCalculations(['RHO'], dataset)
        wpthvp = tvflux / (rho * 1004)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)
        return wpthvp, z

    def getRtp2Fallback(self, dataset_override=None):
        """
        This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the data for plotting.
        THLP2 = QT2 / 1e+6
        :return:
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qt2, z, dataset = self.getVarForCalculations(['QT2'], dataset)
        rtp2 = qt2 / 1e6
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)
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
            dataset = self.les_dataset
        for dataset in dataset.values():
            if 'rc_coef_zm' in dataset.variables.keys() and 'rtprcp' in dataset.variables.keys():
                rc_coef_zm, z, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
                rtprcp, z, dataset = self.getVarForCalculations('rtprcp', dataset)
                rtp3 = rc_coef_zm * (rtprcp)

            elif 'QCFLUX' in dataset.variables.keys():
                QCFLUX, z, dataset = self.getVarForCalculations('QCFLUX', dataset)
                RHO, z, dataset = self.getVarForCalculations('RHO', dataset)
                PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
                THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)
                rtp3 = ((QCFLUX) / (RHO * 2.5104e+6)) * (
                            2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
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
            dataset = self.clubb_datasets['zm']
        rc_coef_zm, z, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        wprcp, z, dataset = self.getVarForCalculations('wprcp', dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        output = rc_coef_zm * wprcp
        return output, z

    def get_rc_coef_zm_X_wprcp_sam_calc(self, dataset_override=None):
        """
        Calculates the Contribution of Cloud Water Flux
        to wpthvp for SAM using the equation

        sam eqn WPRCP * (2.5e6 / (1004.67*((PRES / 1000)^(287.04/1004.67))) - 1.61*THETAV)
        :return:
        """

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.les_dataset
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        WPRCP, z, dataset = self.getVarForCalculations('WPRCP', dataset)
        PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)

        output = WPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, z

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_clubb_fallback(self, dataset_override=None):
        """
        Calculates the Contribution of Cloud Water Flux
        to thlprcp using the equation
        rc_coef_zm * thlprcp
        :return: Line representing rc_coef_zm .* thlprcp
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rc_coef_zm, z, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        thlprcp, z, dataset = self.getVarForCalculations('thlprcp', dataset)

        output = rc_coef_zm * thlprcp
        return output, z

    def get_rc_coef_zm_X_rtprcp_clubb_fallback(self, dataset_override=None):
        """
        Calculates the Contribution of Cloud Water Flux
        to rtprcp using the equation
        rc_coef_zm * rtprcp
        :return: Line representing rc_coef_zm .* rtprcp
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        rc_coef_zm, z, dataset = self.getVarForCalculations('rc_coef_zm', dataset)
        rtprcp, z, dataset = self.getVarForCalculations('rtprcp', dataset)

        output = rc_coef_zm * rtprcp
        return output, z

    def getUwCoampsData(self, dataset_override=None):
        """
        coamps eqn upwp = wpup + wpup_sgs

         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_dataset['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpup, z, dataset = self.getVarForCalculations('wpup', dataset)
        wpup_sgs, z, dataset = self.getVarForCalculations('wpup_sgs', dataset)

        upwp = wpup + wpup_sgs
        return upwp, z

    def getVwCoampsData(self, dataset_override=None):
        """
        coamps eqn vpwp = wpvp + wpvp_sgs

         :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override['sw']
        else:
            dataset = self.coamps_dataset['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpvp, z, dataset = self.getVarForCalculations('wpvp', dataset)
        wpvp_sgs, z, dataset = self.getVarForCalculations('wpvp_sgs', dataset)

        vpwp = wpvp + wpvp_sgs
        return vpwp, z

    def get_rc_coef_zm_X_wprcp_coamps_calc(self, dataset_override=None):
        """coamps eqn thlpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        coamps eqn wpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        wpqcp, z, dataset = self.getVarForCalculations('wqpcp', dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        output = wpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        return output, z

    def get_rc_coef_zm_X_thlprcp_coamps_calc(self, dataset_override=None):
        """
        coamps eqn thlpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        thlpqcp, z, dataset = self.getVarForCalculations('thlpqcp', dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        output = thlpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        return output, z

    def get_rc_coef_zm_X_rtprcp_coamps_calc(self, dataset_override=None):
        """
        coamp eqn qtpqcp .* (2.5e6 ./ (1004.67*ex0) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sw']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        qtpqcp, z, dataset = self.getVarForCalculations('qtpqcp', dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        output = qtpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        return output, z
