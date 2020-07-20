"""
NOTE: Should be obsolete soon(TM) by Panel subclassing
:author: Steffen Domke
:date: Mid 2019
"""
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupTimeHeight(VariableGroup):
    """
    This is a panel group used for testing the functionality of pyplotgen's timeheight plots.
    It contains a set of common panels being used for representing the majority of panels.
    """

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None,
                 time_height=False, anim=False):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "Time-height plots"

        self.variable_definitions = [
            {'var_names':
                {
                    'clubb': ['thlm'],
                    'sam': ['THETAL', 'THETA'],
                    'coamps': ['thlm'],
                    'r408': ['thlm'],
                    'hoc': ['thlm'],
                    'e3sm': ['thlm'],
                    'cam': ['thlm'],
                    'wrf': ['thlm'],
                },
                'sam_calc': self.getThlmSamCalc, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': [],
                    'coamps': ['qtm', 'rtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                    'wrf': ['rtm'],
                },
                'sam_calc': self.getRtmSamCalc, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': ['WPTHLP'],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                    'cam': ['wpthlp'], # WPTHLP_CLUBB / (1 .* 1004)
                    'wrf': ['wpthlp'],
                },
                'sam_calc': self.getWpthlpSamCalc, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['wprtp'],
                    'sam': ['WPRTP'],
                    'coamps': ['wpqtp', 'wprtp'],
                    'r408': ['wprtp'],
                    'hoc': ['wprtp'],
                    'e3sm': ['wprtp'],
                    'cam': ['WPRTP_clubb', 'wprtp'],
                    'wrf': ['wprtp'],
                },
                'sam_calc': self.getWprtpSamCalc, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['cloud_frac'],
                    'sam': ['CLD'],
                    'coamps': ['cf'],
                    'r408': ['cloud_frac', 'cf'],
                    'hoc': ['cloud_frac', 'cf'],
                    'e3sm': ['cloud_frac'],
                    'cam': ['CLOUD', 'cloud_frac'],
                    'wrf': ['cloud_frac'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['rcm'],
                    'sam': ['QCL'],
                    'coamps': ['qcm', 'rcm'],
                    'r408': ['rcm'],
                    'hoc': ['rcm'],
                    'e3sm': ['rcm'],
                    'cam': ['CLDLIQ', 'rcm'],
                    'wrf': ['rcm'],
                },
                'sam_conv_factor': 1 / 1000, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['wp2', 'W2'],
                    'sam': ['W2', 'WP2'],
                    'coamps': ['wp2', 'W2'],
                    'r408': ['wp2'],
                    'hoc': ['wp2'],
                    'e3sm': ['wp2'],
                    'cam': ['WP2_CLUBB', 'wp2'],
                    'wrf': ['wp2'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['wp3'],
                    'sam': ['wp3', 'W3', 'WP3'],
                    'coamps': ['wp3', 'W3', 'WP3'],
                    'r408': ['wp3'],
                    'hoc': ['wp3'],
                    'e3sm': ['wp3'],
                    'cam': ['WP3_CLUBB', 'wp3'],
                    'wrf': ['wp3'],
                },
                'sam_name': 'W3', 'axis_title': "wp3", 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['thlp2'],
                    'sam': ['THLP2', 'TL2'],
                    'coamps': ['thlp2'],
                    'r408': ['thlp2'],
                    'hoc': ['thlp2'],
                    'e3sm': ['thlp2'],
                    'cam': ['THLP2_CLUBB', 'thlp2'],
                    'wrf': ['thlp2'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['rtp2'],
                    'sam': ['RTP2'],
                    'coamps': ['qtp2'],
                    'r408': ['rtp2'],
                    'hoc': ['rtp2'],
                    'e3sm': ['rtp2'],
                    'cam': ['RTP2_CLUBB', 'rtp2'],
                    'wrf': ['rtp2'],
                },
                'sam_calc': self.getRtp2SamCalc, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['rtpthlp'],
                    'sam': ['RTPTHLP', 'TQ'],
                    'coamps': ['qtpthlp', 'rtpthlp'],
                    'r408': ['rtpthlp'],
                    'hoc': ['rtpthlp'],
                    'e3sm': ['rtpthlp'],
                    'cam': ['rtpthlp'],
                    'wrf': ['rtpthlp'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['wm', 'wlsm'],
                    'sam': ['WOBS', 'WM'],
                    'coamps': ['wlsm', 'wm'],
                    'r408': ['wm'],
                    'hoc': ['wm'],
                    'e3sm': ['wm'],
                    'cam': ['wm'], # -OMEGA /(9.81.*1)
                    'wrf': ['wm', 'wlsm'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['um'],
                    'sam': ['U'],
                    'coamps': ['um'],
                    'r408': ['um'],
                    'hoc': ['um'],
                    'e3sm': ['um'],
                    'cam': ['U','um'],
                    'wrf': ['um'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['vm'],
                    'sam': ['V'],
                    'coamps': ['vm'],
                    'r408': ['vm'],
                    'hoc': ['vm'],
                    'e3sm': ['vm'],
                    'cam': ['V', 'vm'],
                    'wrf': ['vm'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['upwp'],
                    'sam': ['UW'],
                    'coamps': ['upwp'],
                    'r408': ['upwp'],
                    'hoc': ['upwp'],
                    'e3sm': ['upwp'],
                    'cam': ['UPWP_CLUBB', 'upwp'],
                    'wrf': ['upwp'],
                },
                'coamps_calc': self.getUwCoampsData, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['vpwp'],
                    'sam': ['VW'],
                    'coamps': ['vpwp'],
                    'r408': ['vpwp'],
                    'hoc': ['vpwp'],
                    'e3sm': ['vpwp'],
                    'cam': ['VPWP_CLUBB', 'vpwp'],
                    'wrf': ['vpwp'],
                },
                'coamps_calc': self.getVwCoampsData, 'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['up2'],
                    'sam': ['U2'],
                    'coamps': ['up2'],
                    'r408': ['up2'],
                    'hoc': ['up2'],
                    'e3sm': ['up2'],
                    'cam': ['UU', 'up2'],
                    'wrf': ['up2'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['vp2'],
                    'sam': ['V2'],
                    'coamps': ['vp2'],
                    'r408': ['vp2'],
                    'hoc': ['vp2'],
                    'e3sm': ['vp2'],
                    'cam': ['VV', 'vp2'],
                    'wrf': ['vp2'],
                },
                'type': Panel.TYPE_TIMEHEIGHT,
            },
            {'var_names':
                {
                    'clubb': ['rcp2'],
                    'sam': ['QC2'],
                    'coamps': ['qcp2', 'rcp2', 'rlp2'],
                    'r408': ['rcp2'],
                    'hoc': ['rcp2'],
                    'e3sm': ['rcp2'],
                    'cam': ['rcp2'],
                    'wrf': ['rcp2'],
                },
                'sam_conv_factor': 1 / 10 ** 6, 'type': Panel.TYPE_TIMEHEIGHT,
            },
        ]
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, les_dataset=les_dataset,
                         coamps_dataset=coamps_dataset, r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets)

    def getThlmSamCalc(self, dataset_override=None):
        """
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA/TABS).*(QI/1000))
        :return: requested variable dependent_data in the form of a list.
                Returned dependent_data is already cropped to the appropriate min,max indices
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

    def getRtmSamCalc(self, dataset_override=None):
        """
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) / 1000
        :return: requested variable dependent_data in the form of a list.
        Returned dependent_data is already cropped to the appropriate min,max indices
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
        WP3 / (WP2 + 1.6e-3)**1.5
        :return: requested variable dependent_data in the form of a list.
                Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset
        if self.coamps_dataset is not None:
            dataset = self.coamps_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override

        wp3, z, dataset = self.getVarForCalculations(['WP3', 'W3', 'wp3'], dataset)
        wp2, z, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        skw_zt = wp3 / (wp2 + 1.6e-3) ** 1.5

        return skw_zt, z

    def getSkrtZtLesCalc(self, dataset_override=None):
        """
        Calculates Skrt_zt values from sam output using
        the following equation
         sam eqn 
            RTP3 / (RTP2 + 4e-16)**1.5
         coamps eqn 
            qtp3 / (qtp2 + 4e-16)**1.5
            rtp3 / (rtp2 + 4e-16)**1.5
        :return: requested variable dependent_data in the form of a list.
                Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = None
        if self.les_dataset is not None:
            dataset = self.les_dataset

        if self.coamps_dataset is not None:
            dataset = self.coamps_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override
        rtp3, z, dataset = self.getVarForCalculations(['RTP3', 'qtp3', 'rtp3'], dataset)
        rtp2, z, dataset = self.getVarForCalculations(['RTP2', 'qtp2', 'rtp2', 'rlp2'], dataset)

        skrt_zt = rtp3 / (rtp2 + 4e-16) ** 1.5

        return skrt_zt, z

    def getSkthlZtLesCalc(self, dataset_override=None):
        """
        Calculates Skthl_zt values from sam output using
        the following equation
        sam THLP3 / (THLP2 + 4e-4)**1.5
        coamps eqn thlp3 / (thlp2 + 4e-4)**1.5
        :return: requested variable dependent_data in the form of a list.
                Returned dependent_data is already cropped to the appropriate min,max indices
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

        skthl_zt = thlp3 / (thlp2 + 4e-4)**1.5
        return skthl_zt, z

    def getWpthlpSamCalc(self, dataset_override=None):
        """
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.
        WPTHLP = (TLFLUX) / (RHO * 1004)
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

    def getWprtpSamCalc(self, dataset_override=None):
        """
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.
        WPRTP = (QTFLUX) / (RHO * 2.5104e+6)
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

    def getWpthvpSamCalc(self, dataset_override=None):
        """
        This gets called if WPTHVP isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.
        WPTHVP =  (TVFLUX) / ( RHO * 1004)
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

    def getRtp2SamCalc(self, dataset_override=None):
        """
        This gets called if RTP2 isn't outputted in an nc file as a backup way of gathering the dependent_data
        for plotting.
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

    def getRtp3SamCalc(self, dataset_override=None):
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
        sam eqn
        WPRCP                          * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)
        ((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)

        :return:
        """

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.les_dataset
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        WPRCP, z, dataset = self.getVarForCalculations('WPRCP', dataset)
        QCFLUX, z, dataset = self.getVarForCalculations('QCFLUX', dataset)
        RHO, z, dataset = self.getVarForCalculations('RHO', dataset)
        PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)


        output = WPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        wprcp_is_zeroes = min(WPRCP) == 0.0 and max(WPRCP) == 0.0
        if wprcp_is_zeroes:
            output = ((QCFLUX) / (RHO * 2.5104e+6)) * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)

        return output, z

    # rc_coef_zm. * thlprcp
    def get_rc_coef_zm_X_thlprcp_clubb_calc(self, dataset_override=None):
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

    def get_rc_coef_zm_X_rtprcp_clubb_calc(self, dataset_override=None):
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

         :return: requested variable dependent_data in the form of a list. Returned dependent_data is already
         cropped to the appropriate min,max indices
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

         :return: requested variable dependent_data in the form of a list. Returned dependent_data is already
         cropped to the appropriate min,max indices
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
        """
        coamps eqn's
        thlpqcp * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)
        wpqcp   * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)
        wprlp   * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)
        wprlp   * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)

        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sw']
        if dataset_override is not None:
            dataset = dataset_override['sw']

        wprlp, z, dataset = self.getVarForCalculations(['thlpqcp', 'wpqcp', 'wprlp'], dataset)
        ex0, z, dataset = self.getVarForCalculations(['ex0'], dataset)
        p, z, dataset = self.getVarForCalculations('p', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)
        output1 = wprlp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = wprlp * (2.5e6 / (1004.67 * ((p / 1.0e5) ** (287.04 / 1004.67))) - 1.61 * thvm)
        output = self.pickNonZeroOutput(output1, output2)
        return output, z

    def get_rc_coef_zm_X_thlprcp_sam_calc(self, dataset_override=None):
        """
        sam eqn
        THLPRCP .* (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)        :param dataset_override:
        :return:
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override
        THLPRCP, z, dataset = self.getVarForCalculations(['THLPRCP'], dataset)
        PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)

        output = THLPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, z

    def get_rc_coef_zm_X_thlprcp_coamps_calc(self, dataset_override=None):
        """
        coamps eqn
        thlpqcp * (2.5e6 / (1004.67 * ex0)                             - 1.61*thvm)
        thlprlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61*thvm)

        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sw']
        if dataset_override is not None:
            dataset = dataset_override
        thlpqcp, z, dataset = self.getVarForCalculations(['thlpqcp'], dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        thlprlp, z, dataset = self.getVarForCalculations(['thlprlp'], dataset)
        p, z, dataset = self.getVarForCalculations('p', dataset)

        output1 = thlpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = thlprlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61*thvm)
        output = self.pickNonZeroOutput(output1, output2)

        return output, z

    def get_rc_coef_zm_X_rtprcp_coamps_calc(self, dataset_override=None):
        """
        coamp eqn
        qtpqcp * (2.5e6 / (1004.67*ex0)                           - 1.61*thvm)
        qtpqcp * (2.5e6 / (1004.67*ex0)                           - 1.61*thvm)
        rtprlp * (2.5e6 / (1004.67*((p/1.0e5)**(287.04/1004.67))) - 1.61*thvm)
        :param dataset_override:
        :return:
        """
        dataset = self.coamps_dataset['sm']
        if dataset_override is not None:
            dataset = dataset_override['sm']
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        qtpqcp, z, dataset = self.getVarForCalculations(['qtpqcp', 'rtprcp'], dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        rtprlp, z, dataset = self.getVarForCalculations('rtprlp', dataset)
        p, z, dataset = self.getVarForCalculations('p', dataset)

        output1 = qtpqcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = rtprlp * (2.5e6 / (1004.67*((p/1.0e5)**(287.04/1004.67))) - 1.61*thvm)
        output = self.pickNonZeroOutput(output1, output2)

        return output, z

    def get_rc_coef_zm_X_rtprcp_sam_calc(self, dataset_override=None):
        """
        sam eqn
        RTPRCP * (2.5e6 / (1004.67*((PRES / 1000)**(287.04/1004.67))) - 1.61*THETAV)m)
        :param dataset_override:
        :return:
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], dataset)
        RTPRCP, z, dataset = self.getVarForCalculations('RTPRCP', dataset)
        PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)

        output = RTPRCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, z

    def get_rc_coef_X_wp2rcp_sam_calc(self, dataset_override=None):
        """
        WP2RCP * (2.5e6 / (1004.67*((PRES / 1000)^(287.04/1004.67))) - 1.61*THETAV)
        :param dataset_override:
        :return:
        """
        dataset = self.sam_datasets
        if dataset_override is not None:
            dataset = dataset_override['sam']
        WP2RCP, z, dataset = self.getVarForCalculations('WP2RCP', dataset)
        PRES, z, dataset = self.getVarForCalculations('PRES', dataset)
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)

        output = WP2RCP * (2.5e6 / (1004.67 * ((PRES / 1000) ** (287.04 / 1004.67))) - 1.61 * THETAV)
        return output, z

    def get_rc_coef_X_wp2rcp_clubb_calc(self, dataset_override=None):
        """

        :param dataset_override:
        :return:
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        rc_coef, z, dataset = self.getVarForCalculations('rc_coef', dataset)
        wp2rcp, z, dataset = self.getVarForCalculations('wp2rcp', dataset)

        output = rc_coef * wp2rcp
        return output, z

    def get_rc_coef_X_wp2rcp_coamps_calc(self, dataset_override=None):
        """
        wp2qcp * (2.5e6 / (1004.67 * ex0)                             - 1.61 * thvm)
        wp2rlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)
        :param dataset_override:
        :return:
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']
        wp2qcp, z, dataset = self.getVarForCalculations('wp2qcp', dataset)
        ex0, z, dataset = self.getVarForCalculations('ex0', dataset)
        thvm, z, dataset = self.getVarForCalculations('thvm', dataset)

        wp2rlp, z, dataset = self.getVarForCalculations('wp2rlp', dataset)
        p, z, dataset = self.getVarForCalculations('p', dataset)

        output1 = wp2qcp * (2.5e6 / (1004.67 * ex0) - 1.61 * thvm)
        output2 = wp2rlp * (2.5e6 / (1004.67 * (( p /1.0e5)**(287.04/1004.67))) - 1.61 * thvm)

        output = self.pickNonZeroOutput(output1, output2)

        return output, z
