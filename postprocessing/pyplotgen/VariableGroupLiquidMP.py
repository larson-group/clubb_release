'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Lineplot import Lineplot


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
            {'clubb_name': 'Ncm'}, # LES
            {'clubb_name': 'Nc_in_cloud'},
            {'clubb_name': 'precip_frac'},
            {'clubb_name': 'rrm', 'sam_name': 'QPL', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'Nrm'},  # (NR * 1e+6) ./ RHO
            # {'clubb_name': 'wprrp', 'sam_name': 'WPRRP'},
            # {'clubb_name': 'wpnrp', 'sam_name': 'WPNRP'},
            {'clubb_name': 'rwp', 'sam_name': 'RWP', 'sam_conv_factor': 1/1000, 'type': Panel.TYPE_TIMESERIES}

        ]
        super().__init__(ncdf_datasets, case, sam_file)
