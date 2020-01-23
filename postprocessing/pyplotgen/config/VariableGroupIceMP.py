"""
:author: Nicolas Strike
:date: Mid 2019
"""
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupIceMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "ice mp variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['rim'],
                'sam': ['QI', 'QCI'],
                'coamps': ['rim', 'qim'],
                'r408': ['rim'],
                'hoc': ['rim'],
                'e3sm': ['rim']
            },
                'sam_conv_factor': 1 / 1000},  # , 'fallback_func': self.getRimFallback},
            {'var_names': {
                'clubb': ['Nim'],
                'sam': ['Nim'],
                'coamps': ['Nim'],
                'r408': ['Nim'],
                'hoc': ['Nim'],
                'e3sm': ['Nim']
            },
                'sam_calc': self.getNimSamLine},
            {'var_names': {
                'clubb': ['rsm'],
                'sam': ['QS', 'QPI'],
                'coamps': ['rsm', 'qsm'],
                'r408': ['rsm'],
                'hoc': ['rsm'],
                'e3sm': ['rsm']
            },
                'sam_conv_factor': 1 / 1000},  # , 'fallback_func': self.getRsmFallback},
            {'var_names': {
                'clubb': ['Nsm'],
                'sam': ['Nsm'],
                'coamps': ['Nsm'],
                'r408': ['Nsm'],
                'hoc': ['Nsm'],
                'e3sm': ['Nsm']
            },
                'sam_calc': self.getNsmSamLine},
            {'var_names': {
                'clubb': ['iwp'],
                'sam': ['IWP'],
                'coamps': ['iwp'],
                'r408': ['iwp'],
                'hoc': ['iwp'],
                'e3sm': ['iwp']
            },
                'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'var_names': {
                'clubb': ['swp'],
                'sam': ['SWP'],
                'coamps': ['swp'],
                'r408': ['swp'],
                'hoc': ['swp'],
                'e3sm': ['swp']
            },
                'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'var_names': {
                'clubb': ['ice_supersat_frac'],
                'sam': ['ice_supersat_frac'],
                'coamps': ['ice_supersat_frac'],
                'r408': ['ice_supersat_frac'],
                'hoc': ['ice_supersat_frac'],
                'e3sm': ['ice_supersat_frac']
            }},

            {'var_names': {
                'clubb': ['Ngm'],
                'sam': ['NG'],
                'coamps': ['Ngm'],
                'r408': ['Ngm'],
                'hoc': ['Ngm'],
                'e3sm': ['Ngm']
            },
                'sam_conv_factor': 10 ** 6, 'fill_zeros': True},
            {'var_names': {
                'clubb': ['rgm'],
                'sam': ['QG'],
                'coamps': ['rgm', 'qgm'],
                'r408': ['rgm'],
                'hoc': ['rgm'],
                'e3sm': ['rgm']
            },
                'sam_conv_factor': 1 / 1000},
            {'var_names': {
                'clubb': ['precip_rate_sfc'],
                'sam': ['PREC'],
                'coamps': ['precip_rate_sfc'],
                'r408': ['precip_rate_sfc'],
                'hoc': ['precip_rate_sfc'],
                'e3sm': ['precip_rate_sfc']
            },
                'type': Panel.TYPE_TIMESERIES}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets)

    def getNimSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NI * 1e+6) ./ RHO
        :return:
        """
        ni, z, dataset = self.getVarForCalculations('NI', self.sam_file, fill_zeros=True)
        rho, z, dataset = self.getVarForCalculations('RHO', self.sam_file, fill_zeros=True)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        nim = (ni * (10 ** 6) / rho)
        return nim, z

    def getNsmSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        """
        ns, z, dataset = self.getVarForCalculations('NS', self.sam_file, fill_zeros=True)
        rho, z, dataset = self.getVarForCalculations('RHO', self.sam_file, fill_zeros=True)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        nsm = (ns * (10 ** 6) / rho)
        return nsm, z
