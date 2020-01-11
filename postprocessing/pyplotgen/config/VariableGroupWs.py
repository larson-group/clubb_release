"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupWs(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_dataset=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "w variables"
        self.variable_definitions = [
            {'aliases': ['wp4', 'WP4']},
            {'aliases': ['wp2thlp', 'WP2THLP'], 'fill_zeros':'True'},
            {'aliases': ['wp2rtp', 'WP2RTP', 'wp2qtp']},
            {'aliases': ['wpthlp2', 'WPTHLP2']},
            {'aliases': ['wprtp2', 'WPRTP2', 'wpqtp2']},
            {'aliases': ['wprtpthlp', 'WPRTPTHLP', 'wpqtpthlp']},
            {'aliases': ['wp2thvp','WP2THVP']},  # TODO LES
        ]
        super().__init__(ncdf_datasets, case, sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_dataset = e3sm_dataset)

