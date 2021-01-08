"""


"""
import numpy as np
from netCDF4 import Dataset

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupNondimMoments(VariableGroup):
    """
    This class contains information for plotting normalized moments, 
    such as correlations or kurtosis.
    """

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "normalized moments"

        self.variable_definitions = [
            {'var_names':
                {
                    'clubb': [self.get_kurtosis_clubb],
                    'sam': [self.get_kurtosis_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'kurtosis, wp4/(wp2**2)',
                'axis_title': r'kurtosis    [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_wpthlp_corr_clubb],
                    'sam': [self.get_wpthlp_corr_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Correlation of w and thetal, wpthlp/sqrt(wp2*thlp2)',
                'axis_title': r'Correlation of w and thetal  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_wprtp_corr_clubb],
                    'sam': [self.get_wprtp_corr_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Correlation of w and rt, wprtp/sqrt(wp2*rtp2)',
                'axis_title': r'Correlation of w and rt  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_wprcp_corr_clubb],
                    'sam': [self.get_wprcp_corr_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Correlation of w and rc, wprcp/sqrt(wp2*rcp2)',
                'axis_title': r'Correlation of w and rc  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_upwp_corr_clubb],
                    'sam': [self.get_upwp_corr_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Correlation of u and w, upwp/sqrt(up2*wp2)',
                'axis_title': r'Correlation of u and w  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_vpwp_corr_clubb],
                    'sam': [self.get_vpwp_corr_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Correlation of v and w, vpwp/sqrt(vp2*wp2)',
                'axis_title': r'Correlation of v and w  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_nondim_wpthlp2_clubb],
                    'sam': [self.get_nondim_wpthlp2_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Nondimensionalized wpthlp2, wpthlp2/(sqrt(wp2)*thlp2)',
                'axis_title': r'Nondimensionalized wpthlp2  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_nondim_wprtp2_clubb],
                    'sam': [self.get_nondim_wprtp2_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Nondimensionalized wprtp2, wprtp2/(sqrt(wp2)*rtp2)',
                'axis_title': r'Nondimensionalized wprtp2  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_nondim_wp2thlp_clubb],
                    'sam': [self.get_nondim_wp2thlp_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Nondimensionalized wp2thlp, wp2thlp/(wp2*sqrt(thlp2))',
                'axis_title': r'Nondimensionalized wp2thlp  [-]',
            },
            {'var_names':
                {
                    'clubb': [self.get_nondim_wp2rtp_clubb],
                    'sam': [self.get_nondim_wp2rtp_sam],
                    'coamps': [''],
                    'r408': [''],
                    'hoc': [''],
                    'e3sm': [''],
                    'cam': [''],
                    'wrf': [''],
                },
                'sci_scale': 0,
                'priority': True,
                'title': r'Nondimensionalized wp2rtp, wp2rtp/(wp2*sqrt(rtp2))',
                'axis_title': r'Nondimensionalized wp2rtp  [-]',
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, les_dataset=les_dataset,
                         coamps_dataset=coamps_dataset, r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars)

    def get_kurtosis_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wp4, indep, dataset = self.getVarForCalculations('wp4', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)

        kurtosis = wp4 / ( wp2 * wp2 )
        return kurtosis, indep

    def get_kurtosis_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp4, indep, dataset = self.getVarForCalculations('WP4', dataset)
        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        kurtosis = wp4 / ( wp2 * wp2 )
        return kurtosis, indep

    def get_wpthlp_corr_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wpthlp, indep, dataset = self.getVarForCalculations('wpthlp', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        thlp2, indep, dataset = self.getVarForCalculations('thlp2', dataset)

        wpthlp_corr = wpthlp / np.sqrt( wp2 * thlp2 )
        return wpthlp_corr, indep

    def get_wpthlp_corr_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        TL2, indep, dataset = self.getVarForCalculations('TL2', dataset)
        THLP2_SGS, indep, dataset = self.getVarForCalculations('THLP2_SGS', dataset)
        thlp2 = TL2 + THLP2_SGS

        tlflux, indep, dataset = self.getVarForCalculations(['TLFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        wpthlp_sgs, indep, dataset = self.getVarForCalculations(['WPTHLP_SGS'], dataset)
        wpthlp = (tlflux / (rho * 1004))
        if not np.any(np.isnan(wpthlp_sgs)):
            wpthlp += wpthlp_sgs

        wpthlp_corr = wpthlp / np.sqrt( wp2 * thlp2 )
        return wpthlp_corr, indep

    def get_wprtp_corr_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wprtp, indep, dataset = self.getVarForCalculations('wprtp', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        rtp2, indep, dataset = self.getVarForCalculations('rtp2', dataset)

        wprtp_corr = wprtp / np.sqrt( wp2 * rtp2 )
        return wprtp_corr, indep

    def get_wprtp_corr_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        QT2, indep, dataset = self.getVarForCalculations(['QT2'], dataset)
        RTP2_SGS, indep, dataset = self.getVarForCalculations(['RTP2_SGS'], dataset)
        rtp2 = (QT2 / 1e+6) + RTP2_SGS

        qtflux, indep, dataset = self.getVarForCalculations(['QTFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        wprtp_sgs, indep, dataset = self.getVarForCalculations(['WPRTP_SGS'], dataset)
        wprtp = qtflux / (rho * 2.5104e+6)
        if not np.any(np.isnan(wprtp_sgs)):
            wprtp += wprtp_sgs

        wprtp_corr = wprtp / np.sqrt( wp2 * rtp2 )
        return wprtp_corr, indep

    def get_wprcp_corr_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wprcp, indep, dataset = self.getVarForCalculations('wprcp', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        rcp2, indep, dataset = self.getVarForCalculations('rcp2', dataset)

        wprcp_corr = wprcp / np.sqrt( wp2 * rcp2 )
        return wprcp_corr, indep

    def get_wprcp_corr_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        rcp2, indep, dataset = self.getVarForCalculations(['RCP2'], dataset)

        wprcp, indep, dataset = self.getVarForCalculations(['WPRCP'], dataset)

        wprcp_corr = wprcp / np.sqrt( wp2 * rcp2 )
        return wprcp_corr, indep

    def get_upwp_corr_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        upwp, indep, dataset = self.getVarForCalculations('upwp', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        up2, indep, dataset = self.getVarForCalculations('up2', dataset)

        upwp_corr = upwp / np.sqrt( wp2 * up2 )
        return upwp_corr, indep

    def get_upwp_corr_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        U2, z, dataset = self.getVarForCalculations('U2', dataset)
        UP2_SGS, z, dataset = self.getVarForCalculations('UP2_SGS', dataset)
        up2 = UP2_SGS + U2

        UW, z, dataset = self.getVarForCalculations('UW', dataset)
        UPWP_SGS, z, dataset = self.getVarForCalculations('UPWP_SGS', dataset)
        upwp = UW + UPWP_SGS

        upwp_corr = upwp / np.sqrt( wp2 * up2 )
        return upwp_corr, indep

    def get_vpwp_corr_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        vpwp, indep, dataset = self.getVarForCalculations('vpwp', dataset)
        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        vp2, indep, dataset = self.getVarForCalculations('vp2', dataset)

        vpwp_corr = vpwp / np.sqrt( wp2 * vp2 )
        return vpwp_corr, indep

    def get_vpwp_corr_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        V2, z, dataset = self.getVarForCalculations('V2', dataset)
        VP2_SGS, z, dataset = self.getVarForCalculations('VP2_SGS', dataset)
        vp2 = V2 + VP2_SGS 

        VW, z, dataset = self.getVarForCalculations('VW', dataset)
        VPWP_SGS, z, dataset = self.getVarForCalculations('VPWP_SGS', dataset)
        vpwp = VW + VPWP_SGS

        vpwp_corr = vpwp / np.sqrt( wp2 * vp2 )
        return vpwp_corr, indep

    def get_nondim_wpthlp2_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zt']

        wpthlp2, indep, dataset = self.getVarForCalculations('wpthlp2', dataset)


        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        thlp2, indep, dataset = self.getVarForCalculations('thlp2', dataset)

        nondim_wpthlp2 = wpthlp2 / ( np.sqrt( wp2 ) * thlp2 )
        return nondim_wpthlp2, indep

    def get_nondim_wpthlp2_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        TL2, indep, dataset = self.getVarForCalculations('TL2', dataset)
        THLP2_SGS, indep, dataset = self.getVarForCalculations('THLP2_SGS', dataset)
        thlp2 = TL2 + THLP2_SGS

        wpthlp2, indep, dataset = self.getVarForCalculations(['WPTHLP2'], dataset)

        nondim_wpthlp2 = wpthlp2 / ( np.sqrt( wp2 ) * thlp2 )
        return nondim_wpthlp2, indep

    def get_nondim_wprtp2_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zt']

        wprtp2, indep, dataset = self.getVarForCalculations('wprtp2', dataset)


        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        rtp2, indep, dataset = self.getVarForCalculations('rtp2', dataset)

        nondim_wprtp2 = wprtp2 / ( np.sqrt( wp2 ) * rtp2 )
        return nondim_wprtp2, indep

    def get_nondim_wprtp2_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        QT2, indep, dataset = self.getVarForCalculations(['QT2'], dataset)
        RTP2_SGS, indep, dataset = self.getVarForCalculations(['RTP2_SGS'], dataset)
        rtp2 = (QT2 / 1e+6) + RTP2_SGS

        wprtp2, indep, dataset = self.getVarForCalculations(['WPRTP2'], dataset)

        nondim_wprtp2 = wprtp2 / ( np.sqrt( wp2 ) * rtp2 )
        return nondim_wprtp2, indep

    def get_nondim_wp2thlp_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zt']

        wp2thlp, indep, dataset = self.getVarForCalculations('wp2thlp', dataset)


        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        thlp2, indep, dataset = self.getVarForCalculations('thlp2', dataset)

        nondim_wp2thlp = wp2thlp / ( wp2 * np.sqrt( thlp2 ) )
        return nondim_wp2thlp, indep

    def get_nondim_wp2thlp_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        TL2, indep, dataset = self.getVarForCalculations('TL2', dataset)
        THLP2_SGS, indep, dataset = self.getVarForCalculations('THLP2_SGS', dataset)
        thlp2 = TL2 + THLP2_SGS

        wp2thlp, indep, dataset = self.getVarForCalculations(['WP2THLP'], dataset)

        nondim_wp2thlp = wp2thlp / ( wp2 * np.sqrt( thlp2 ) )
        return nondim_wp2thlp, indep

    def get_nondim_wp2rtp_clubb(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zt']

        wp2rtp, indep, dataset = self.getVarForCalculations('wp2rtp', dataset)


        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.clubb_datasets['zm']

        wp2, indep, dataset = self.getVarForCalculations('wp2', dataset)
        rtp2, indep, dataset = self.getVarForCalculations('rtp2', dataset)

        nondim_wp2rtp = wp2rtp / ( wp2 * np.sqrt( rtp2 ) )
        return nondim_wp2rtp, indep

    def get_nondim_wp2rtp_sam(self, dataset_override=None):
        """

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model.
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being calculated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        wp2, indep, dataset = self.getVarForCalculations(['WP2', 'W2', 'wp2'], dataset)

        QT2, indep, dataset = self.getVarForCalculations(['QT2'], dataset)
        RTP2_SGS, indep, dataset = self.getVarForCalculations(['RTP2_SGS'], dataset)
        rtp2 = (QT2 / 1e+6) + RTP2_SGS

        wp2rtp, indep, dataset = self.getVarForCalculations(['WP2RTP'], dataset)

        nondim_wp2rtp = wp2rtp / ( wp2 * np.sqrt( rtp2 ) )
        return nondim_wp2rtp, indep

