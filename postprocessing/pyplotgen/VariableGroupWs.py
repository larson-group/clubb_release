'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Lineplot import Lineplot


class VariableGroupWs(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "w variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'wp4', 'sam_name': 'WP4'},
            {'clubb_name': 'wp2thlp', 'sam_name': 'WP2THLP'},
            {'clubb_name': 'wp2rtp', 'sam_name': 'WP2RTP'},
            {'clubb_name': 'wpthlp2', 'sam_name': 'WPTHLP2'},
            {'clubb_name': 'wprtp2', 'sam_name': 'WPRTP2'},
            {'clubb_name': 'wprtpthlp', 'sam_name': 'WPRTPTHLP'},
            {'clubb_name': 'wp2thvp', 'sam_name': 'WP2THVP'} # LES

        ]
        super().__init__(ncdf_datasets, case, sam_file)

