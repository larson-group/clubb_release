'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from Panel import Panel
from VariableGroup import VariableGroup
from Line import Line


class VariableGroupKKMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "kk mp variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'rrm_cond', 'sam_name': 'EVAPM'},
            {'clubb_name': 'rrm_accr', 'sam_name': 'ACCRM'},
            {'clubb_name': 'rrm_auto', 'sam_name': 'AUTOM'}
        ]
        super().__init__(ncdf_datasets, case, sam_file)
