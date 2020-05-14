"""
:author: Nicolas Strike
:date: Mid 2019
"""
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupIceMP(VariableGroup):

    def __init__(self, case, clubb_datasets=None, les_dataset=None, sam_datasets=None, coamps_dataset=None,
                 r408_dataset=None, hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, wrf_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "ice mp variables"
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['rim'],
                'sam': ['QI', 'QCI'],
                'coamps': ['rim', 'qim'],
                'r408': ['rim'],
                'hoc': ['rim'],
                'e3sm': ['rim'],
                'cam': ['rim'],
                'wrf': [],
                },
             'sam_conv_factor': 1 / 1000,
             #'fallback_func': self.getRimFallback,
            },
            {'var_names':
                {
                'clubb': ['Nim'],
                'sam': ['Nim'],
                'coamps': ['Nim'],
                'r408': ['Nim'],
                'hoc': ['Nim'],
                'e3sm': ['Nim'],
                'cam': ['Nim'],
                'wrf': [],
                },
             'sam_calc': self.getNimSamLine,
            },
            {'var_names':
                {
                'clubb': ['rsm'],
                'sam': ['QS', 'QPI'],
                'coamps': ['rsm', 'qsm'],
                'r408': ['rsm'],
                'hoc': ['rsm'],
                'e3sm': ['rsm'],
                'cam': ['rsm'],
                'wrf': [],
                },
             'sam_conv_factor': 1 / 1000,
             #'fallback_func': self.getRsmFallback,
            },
            {'var_names':
                {
                'clubb': ['Nsm'],
                'sam': ['Nsm'],
                'coamps': ['Nsm'],
                'r408': ['Nsm'],
                'hoc': ['Nsm'],
                'e3sm': ['Nsm'],
                'cam': ['Nsm'],
                'wrf': [],
                },
             'sam_calc': self.getNsmSamLine,
            },
            {'var_names':
                {
                'clubb': ['iwp'],
                'sam': ['IWP'],
                'coamps': ['iwp'],
                'r408': ['iwp'],
                'hoc': ['iwp'],
                'e3sm': ['iwp'],
                'cam': ['iwp'],
                'wrf': [],
                },
             'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000,
            },
            {'var_names':
                {
                'clubb': ['swp'],
                'sam': ['SWP'],
                'coamps': ['swp'],
                'r408': ['swp'],
                'hoc': ['swp'],
                'e3sm': ['swp'],
                'cam': ['swp'],
                'wrf': [],
                },
             'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000,
            },
            {'var_names':
                {
                'clubb': ['ice_supersat_frac'],
                'sam': ['ice_supersat_frac'],
                'coamps': ['ice_supersat_frac'],
                'r408': ['ice_supersat_frac'],
                'hoc': ['ice_supersat_frac'],
                'e3sm': ['ice_supersat_frac'],
                'cam': ['ice_supersat_frac'],
                'wrf': [],
                },
            },
            {'var_names':
                {
                'clubb': ['Ngm'],
                'sam': ['NG'],
                'coamps': ['Ngm'],
                'r408': ['Ngm'],
                'hoc': ['Ngm'],
                'e3sm': ['Ngm'],
                'cam': ['Ngm'],
                'wrf': [],
                },
             'sam_conv_factor': 10 ** 6,
            },
            {'var_names':
                {
                'clubb': ['rgm'],
                'sam': ['QG'],
                'coamps': ['rgm', 'qgm'],
                'r408': ['rgm'],
                'hoc': ['rgm'],
                'e3sm': ['rgm'],
                'cam': ['rgm'],
                'wrf': [],
                },
             'sam_conv_factor': 1 / 1000,
            },
            {'var_names':
                {
                'clubb': ['precip_rate_sfc'],
                'sam': ['PREC'],
                'coamps': ['precip_rate_sfc'],
                'r408': ['precip_rate_sfc'],
                'hoc': ['precip_rate_sfc'],
                'e3sm': ['precip_rate_sfc'],
                'cam': ['precip_rate_sfc'],
                'wrf': [],
                },
             'type': Panel.TYPE_TIMESERIES
            },
        ]
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets)

    def getNimSamLine(self, dataset_override=None):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NI * 1e+6) ./ RHO
        :return:
        """
        dataset = self.les_dataset
        if dataset_override != None:
            dataset = dataset_override
        ni, z, dataset = self.getVarForCalculations('NI', dataset)
        rho, z, dataset = self.getVarForCalculations('RHO', dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)

        nim = (ni * (10 ** 6) / rho)
        return nim, z

    def getNsmSamLine(self, dataset_override=None):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        """
        dataset = self.les_dataset
        if dataset_override != None:
            dataset = dataset_override
        ns, z, dataset = self.getVarForCalculations('NS', dataset)
        rho, z, dataset = self.getVarForCalculations('RHO', dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_dataset)

        nsm = (ns * (10 ** 6) / rho)
        return nsm, z
