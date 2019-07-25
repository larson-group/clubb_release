'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Line import Line


class VariableGroupBaseBudgets(VariableGroup):
    '''

    '''

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "base variables budgets"
        
        thlm_overplots = [
            {'clubb_name': 'thlm_bt', 'label': 'thlm_bt'},
            {'clubb_name': 'thlm_ma', 'label': 'thlm_ma'},
            {'clubb_name': 'thlm_ta', 'label': 'thlm_ta'},
            {'clubb_name': 'thlm_mc', 'label': 'thlm_mc'},
            # {'clubb_name': 'thlm_clipping', 'label': 'thlm_bt', 'fallback_func': self.getThlmClipping},
            {'clubb_name': 'radht', 'label': 'radht'},
            # {'clubb_name': 'lsforcing', 'label': 'lsforcing', 'fallback_func': self.getThlmLsforcing},
            # {'clubb_name': 'thlm_residual', 'label': 'thlm_residual', 'fallback_func': self.getThlmResidual},

        ]
        self.variable_definitions = [
            {'clubb_name': 'thlm', 'lines': thlm_overplots, 'type': Panel.TYPE_BUDGET},

        ]
        super().__init__(ncdf_datasets, case, sam_file)

    def getThlmClipping(self):
        '''


        thlm_mfl+thlm_cl+thlm_tacl+thlm_sdmp
        :return:
        '''

    def getLsforcing(self):
        '''


        thlm_forcing-radht-thlm_mc
        :return:
        '''