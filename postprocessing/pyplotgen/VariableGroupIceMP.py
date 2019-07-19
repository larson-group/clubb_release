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
            {'clubb_name': 'rim', 'sam_name': 'QI', 'sam_conv_factor': 1 / 1000, 'fallback_func': self.getRimFallback},
            {'clubb_name': 'Nim', 'sam_calc': self.getNimSamLine},
            {'clubb_name': 'rsm', 'sam_name': 'QS', 'sam_conv_factor': 1 / 1000, 'fallback_func': self.getRsmFallback},
            {'clubb_name': 'Nsm', 'sam_calc': self.getNsmSamLine},
            {'clubb_name': 'iwp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'IWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'swp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'SWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'ice_supersat_frac'},
            {'clubb_name': 'Ngm', 'sam_name': 'NG', 'sam_conv_factor': 10 ** 6, 'fill_zeros': True},
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
        sam_start_time = self.start_time  # / sec_per_min
        sam_end_time = self.end_time  # / sec_per_min

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)
        z = z_ncdf.data
        start_idx, end_idx = self.__getStartEndIndex__(z, self.height_min_value, self.height_max_value)

        ni_ncdf = NetCdfVariable('NI', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        ni = ni_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data
        nim = (ni * (10 ** 6) / rho)
        nim = nim[start_idx:end_idx]
        nim_line = Line(nim, z, line_format='k-', label='LES output')
        return nim_line

    def getNsmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NS * 1e+6) ./ RHO
        :return:
        '''
        sec_per_min = 60
        sam_start_time = self.start_time  # / sec_per_min
        sam_end_time = self.end_time  # / sec_per_min

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)
        z = z_ncdf.data
        start_idx, end_idx = self.__getStartEndIndex__(z, self.height_min_value, self.height_max_value)

        ns_ncdf = NetCdfVariable('NS', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        ns = ns_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        rho = rho_ncdf.data
        nsm = (ns * (10 ** 6) / rho)
        nsm = nsm[start_idx:end_idx]
        nsm_line = Line(nsm, z, line_format='k-', label='LES output')
        return nsm_line

    def getRimFallback(self):
        '''
        This gets called if Rim isn't outputted in an nc file as a backup way of gathering the data for plotting.
        Rim = QCI / 1000
        :return:
        '''
        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        qci_ncdf = NetCdfVariable('QCI', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        qci_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qci = qci_ncdf.data

        rim = qci / 1000
        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        rim = Line(rim, z_ncdf.data, line_format="k-", label="LES output")
        return rim

    def getRsmFallback(self):
        '''
        This gets called if Rsm isn't outputted in an nc file as a backup way of gathering the data for plotting.
        Rsm = QPI / 1000
        :return:
        '''
        sam_start_time = self.start_time  # / 60
        sam_end_time = self.end_time  # / 60

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        qpi_ncdf = NetCdfVariable('QCI', self.sam_file, 1, start_time=sam_start_time, end_time=sam_end_time)
        qpi_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        qpi = qpi_ncdf.data

        rsm = qpi / 1000

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        rsm = Line(rsm, z_ncdf.data, line_format="k-", label="LES output")
        return rsm
