'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from VariableGroup import VariableGroup


class VariableGroupKKMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "kk mp variables"
        self.variable_definitions = [
            {'clubb_name': 'rrm_cond', 'sam_name': 'EVAPM'},
            {'clubb_name': 'rrm_accr', 'sam_name': 'ACCRM'},
            {'clubb_name': 'rrm_auto', 'sam_name': 'AUTOM'}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file)
