'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Line import Line


class VariableGroupLiquidMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "liquid mp variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'Ncm', 'sam_calc': self.getNcmSamLine},
            {'clubb_name': 'Nc_in_cloud'},
            {'clubb_name': 'precip_frac'},
            {'clubb_name': 'rrm', 'sam_name': 'QPL', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'Nrm', 'sam_calc': self.getNrmSamLine},
            # {'clubb_name': 'wprrp', 'sam_name': 'WPRRP'},  # Not found in lba case file
            # {'clubb_name': 'wpnrp', 'sam_name': 'WPNRP'},  # Not found in lba case file
            {'clubb_name': 'rwp', 'sam_name': 'RWP', 'sam_conv_factor': 1/1000, 'type': Panel.TYPE_TIMESERIES}

        ]
        super().__init__(ncdf_datasets, case, sam_file)


    def getNcmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NC * 1e+6) ./ RHO
        :return:
        '''
        sec_per_min = 60
        sam_start_time = self.averaging_start_time / sec_per_min
        sam_end_time = self.averaging_end_time / sec_per_min
        nc_ncdf = NetCdfVariable('NC', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time, fill_zeros=True)
        nc = nc_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time)
        rho = rho_ncdf.data
        nc = nc[self.z_sam_min_idx:self.z_sam_max_idx]
        rho = rho[self.z_sam_min_idx:self.z_sam_max_idx]
        ncm = (nc * (10 ** 6) / rho)
        ncm = ncm[self.z_sam_min_idx:self.z_sam_max_idx]
        ncm_line = Line(ncm, self.z_sam.data, line_format='k-', label='LES output')
        return ncm_line

    def getNrmSamLine(self):
        '''
        Caclulates Nim from sam -> clubb using the equation
        (NR * 1e+6) ./ RHO
        :return:
        '''
        sec_per_min = 60
        sam_start_time = self.averaging_start_time / sec_per_min
        sam_end_time = self.averaging_end_time / sec_per_min
        nr_ncdf = NetCdfVariable('NR', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time, fill_zeros=True)
        nr = nr_ncdf.data
        rho_ncdf = NetCdfVariable('RHO', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time)
        rho = rho_ncdf.data
        nr = nr[self.z_sam_min_idx:self.z_sam_max_idx]
        rho = rho[self.z_sam_min_idx:self.z_sam_max_idx]
        nrm = (nr * (10 ** 6) / rho)
        nrm = nrm[self.z_sam_min_idx:self.z_sam_max_idx]
        nrm_line = Line(nrm, self.z_sam.data, line_format='k-', label='LES output')
        return nrm_line