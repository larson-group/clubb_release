'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Line import Line


class VariableGroupIceMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "ice mp variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'rim', 'sam_name': 'QI', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'Nim', 'sam_calc': self.getNimSamLine},
            {'clubb_name': 'rsm', 'sam_name': 'QS', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'Nsm', 'sam_calc': self.getNsmSamLine},
            {'clubb_name': 'iwp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'IWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'swp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'SWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'ice_supersat_frac'},
            {'clubb_name': 'Ngm', 'sam_name': 'NG', 'sam_conv_factor': 10 ** 6},
            {'clubb_name': 'rgm', 'sam_name': 'QG', 'sam_conv_factor': 1 / 1000}
            # {'clubb_name': 'morr_rain_rate'}
        ]
        super().__init__(ncdf_datasets, case, sam_file)


    def getNimSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NI * 1e+6) ./ RHO
        :return:
        '''
        sec_per_min = 60
        sam_start_time = self.start_time / sec_per_min
        sam_end_time = self.end_time / sec_per_min
        ni_ncdf = NetCdfVariable('NI', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        ni = ni_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data
        nim = (ni * (10 ** 6) / rho)
        nim = nim[self.z_sam_min_idx:self.z_sam_max_idx]
        nim_line = Line(nim, self.z_sam.data, line_format='k-', label='LES output')
        return nim_line

    def getNsmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        '''
        sec_per_min = 60
        sam_start_time = self.start_time / sec_per_min
        sam_end_time = self.end_time / sec_per_min
        ns_ncdf = NetCdfVariable('NS', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        ns = ns_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data
        nsm = (ns * (10 ** 6) / rho)
        nsm = nsm[self.z_sam_min_idx:self.z_sam_max_idx]
        nsm_line = Line(nsm, self.z_sam.data, line_format='k-', label='LES output')
        return nsm_line