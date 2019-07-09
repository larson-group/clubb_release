'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Lineplot import Lineplot


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
            {'clubb_name': 'Nim'},  # (NI * 1e+6) ./ RHO
            {'clubb_name': 'rsm', 'sam_name': 'QS', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'Nsm'},  # (NS * 1e+6) ./ RHO
            {'clubb_name': 'iwp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'IWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'swp', 'type': Panel.TYPE_TIMESERIES, 'sam_name': 'SWP', 'sam_conv_factor': 1 / 1000},
            {'clubb_name': 'ice_supersat_frac'},
            {'clubb_name': 'Ngm'},  # NG * 1e+6
            {'clubb_name': 'rgm', 'sam_name': 'QG', 'sam_conv_factor': 1 / 1000}
            # {'clubb_name': 'morr_rain_rate'}
        ]
        super().__init__(ncdf_datasets, case, sam_file)
