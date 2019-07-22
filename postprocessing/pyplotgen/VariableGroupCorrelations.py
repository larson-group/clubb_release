'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Line import Line


class VariableGroupCorrelations(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "w variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'corr_w_rr_1'},
            {'clubb_name': 'corr_w_Nr_1'},
            {'clubb_name': 'corr_w_Ncn_1'},
            {'clubb_name': 'corr_chi_rr_1'},
            {'clubb_name': 'corr_chi_Nr_1'},
            {'clubb_name': 'corr_chi_Ncn_1'},
            {'clubb_name': 'corr_rr_Nr_1'},

        ]
        super().__init__(ncdf_datasets, case, sam_file)
