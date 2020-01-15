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
            {'var_names': ['rim', 'QI', 'QCI'], 'sam_conv_factor': 1 / 1000},#, 'fallback_func': self.getRimFallback},
            {'var_names': ['Nim'], 'sam_calc': self.getNimSamLine},
            {'var_names': ['rsm', 'QS', 'QPI'], 'sam_conv_factor': 1 / 1000},#, 'fallback_func': self.getRsmFallback},
            {'var_names': ['Nsm'], 'sam_calc': self.getNsmSamLine},
            {'var_names': ['iwp', 'IWP'], 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'var_names': ['swp', 'SWP'], 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'var_names': ['ice_supersat_frac']},
            {'var_names': ['Ngm', 'NG'], 'sam_conv_factor': 10 ** 6, 'fill_zeros': True},
            {'var_names': ['rgm', 'QG'], 'sam_conv_factor': 1 / 1000},
            {'var_names': ['precip_rate_sfc', 'PREC'], 'type': Panel.TYPE_TIMESERIES}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets= e3sm_datasets)

    def getNimSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NI * 1e+6) ./ RHO
        :return:
        """
        ni,z, dataset = self.getVarForCalculations('NI', self.sam_file, fill_zeros=True)
        rho,z, dataset = self.getVarForCalculations('RHO', self.sam_file, fill_zeros=True)
       # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        nim = (ni * (10 ** 6) / rho)
        return nim, z

    def getNsmSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        """
        ns,z, dataset = self.getVarForCalculations('NS', self.sam_file, fill_zeros=True)
        rho,z, dataset = self.getVarForCalculations('RHO', self.sam_file, fill_zeros=True)
       # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        nsm = (ns * (10 ** 6) / rho)
        return nsm, z

