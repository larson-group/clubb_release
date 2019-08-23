'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from Line import Line
from Panel import Panel
from VariableGroup import VariableGroup


class VariableGroupLiquidMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "liquid mp variables"
        self.variable_definitions = [
            {'aliases': ['Ncm'], 'sam_calc': self.getNcmSamLine},
            {'aliases': ['Nc_in_cloud']},
            {'aliases': ['precip_frac']},
            {'aliases': ['rrm', 'QPL'], 'sam_conv_factor': 1 / 1000},
            {'aliases': ['Nrm'], 'sam_calc': self.getNrmSamLine},
            {'aliases': ['wprrp', 'WPRRP']},  # Not found in lba case file
            {'aliases': ['wpNrp', 'WPNRP']},  # Not found in lba case file
            {'aliases': ['rwp', 'RWP'], 'sam_conv_factor': 1 / 1000, 'type': Panel.TYPE_TIMESERIES},
            {'aliases': ['precip_rate_sfc'], 'type': Panel.TYPE_TIMESERIES}

        ]
        #rain_rate_sfc vs time

        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)


    def getNcmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NC * 1e+6) ./ RHO
        :return:
        '''

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        nc_ncdf = NetCdfVariable('NC', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        nc_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        nc = nc_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        rho_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rho = rho_ncdf.data

        ncm = (nc * (10 ** 6) / rho)

        z_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        ncm_line = Line(ncm, z_ncdf.data, line_format='k-', label='SAM-LES')
        return ncm_line

    def getNrmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NR * 1e+6) ./ RHO
        :return:
        '''

        z_ncdf = NetCdfVariable('z', self.sam_file, 1)

        nr_ncdf = NetCdfVariable('NR', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        nr_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        nr = nr_ncdf.data

        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, start_time=self.start_time, end_time=self.end_time)
        rho_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rho = rho_ncdf.data

        nrm = (nr * (10 ** 6) / rho)

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        nrm_line = Line(nrm, z_ncdf.data, line_format='k-', label='SAM-LES')
        return nrm_line