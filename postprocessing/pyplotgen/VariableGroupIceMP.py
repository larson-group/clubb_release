'''
:author: Nicolas Strike
:date: Mid 2019
'''

from Line import Line
from Panel import Panel
from VariableGroup import VariableGroup


class VariableGroupIceMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "ice mp variables"
        self.variable_definitions = [
            {'aliases': ['rim', 'QI', 'QCI'], 'sam_conv_factor': 1 / 1000},#, 'fallback_func': self.getRimFallback},
            {'aliases': ['Nim'], 'sam_calc': self.getNimSamLine},
            {'aliases': ['rsm', 'QS', 'QPI'], 'sam_conv_factor': 1 / 1000},#, 'fallback_func': self.getRsmFallback},
            {'aliases': ['Nsm'], 'sam_calc': self.getNsmSamLine},
            {'aliases': ['iwp', 'IWP'], 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'aliases': ['swp', 'SWP'], 'type': Panel.TYPE_TIMESERIES, 'sam_conv_factor': 1 / 1000},
            {'aliases': ['ice_supersat_frac']},
            {'aliases': ['Ngm', 'NG'], 'sam_conv_factor': 10 ** 6, 'fill_zeros': True},
            {'aliases': ['rgm', 'QG'], 'sam_conv_factor': 1 / 1000},
            {'aliases': ['precip_rate_sfc', 'PREC'], 'type': Panel.TYPE_TIMESERIES}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)

    def getNimSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NI * 1e+6) ./ RHO
        :return:
        '''
        z = self.__getFallbackVar__('z', self.sam_file)
        ni = self.__getFallbackVar__('NI', self.sam_file, fill_zeros=True)
        rho = self.__getFallbackVar__('RHO', self.sam_file, fill_zeros=True)

        nim = (ni * (10 ** 6) / rho)
        nim_line = Line(nim, z, line_format='k-', label='SAM-LES')
        return nim_line

    def getNsmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        '''

        z = self.__getFallbackVar__('z', self.sam_file)
        ns = self.__getFallbackVar__('NS', self.sam_file, fill_zeros=True)
        rho = self.__getFallbackVar__('RHO', self.sam_file, fill_zeros=True)

        nsm = (ns * (10 ** 6) / rho)
        nsm_line = Line(nsm, z, line_format='k-', label='SAM-LES')
        return nsm_line

