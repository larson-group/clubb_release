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
            {'clubb_name': 'corr_w_rr_1'},#, 'sam_name': 'WP4'},
            {'clubb_name': 'corr_w_Nr_1'},#, 'sam_name': 'WP2THLP'},
            {'clubb_name': 'corr_w_Ncn_1'},#, 'sam_name': 'WP2RTP'},
            {'clubb_name': 'corr_chi_rr_1'},#, 'sam_name': 'WPTHLP2'},
            {'clubb_name': 'corr_chi_Nr_1'},#, 'sam_name': 'WPRTP2'},
            {'clubb_name': 'corr_chi_Ncn_1'},#, 'sam_name': 'WPRTPTHLP'},
            {'clubb_name': 'corr_rr_Nr_1'},#, 'sam_name': 'WP2THVP'},
            # {'clubb_name': 'corr_w_Ncn_1'},#, 'sam_name': 'WPTHLP2'},
            # {'clubb_name': 'corr_chi_Ncn_1'},#, 'sam_name': 'WPRTP2'},
            # {'clubb_name': 'corr_rr_Nr_1'},#, 'sam_name': 'WPRTPTHLP'}

        ]
        super().__init__(ncdf_datasets, case, sam_file)

