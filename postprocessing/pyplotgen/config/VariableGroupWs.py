"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupWs(VariableGroup):

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "w variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['wp4'],
                'sam': ['WP4'],
                'coamps': ['wp4'],
                'r408': ['wp4'],
                'hoc': ['wp4'],
                'e3sm': ['wp4'],
                'cam': ['wp4'],
                'wrf': []
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['wp2thlp'],
                'sam': ['WP2THLP'],
                'coamps': ['wp2thlp'],
                'r408': ['wp2thlp'],
                'hoc': ['wp2thlp'],
                'e3sm': ['wp2thlp'],
                'cam': ['wp2thlp'],
                'wrf': []
            }, 'sci_scale': 0},
            {'var_names': {
                'clubb': ['wp2rtp'],
                'sam': ['WP2RTP'],
                'coamps': ['wp2qtp'],
                'r408': ['wp2rtp'],
                'hoc': ['wp2rtp'],
                'e3sm': ['wp2rtp'],
                'cam': ['wp2rtp'],
                'wrf': []
            }, 'sci_scale': -4},

            {'var_names': {
                'clubb': ['wpthlp2'],
                'sam': ['WPTHLP2'],
                'coamps': ['wpthlp2'],
                'r408': ['wpthlp2'],
                'hoc': ['wpthlp2'],
                'e3sm': ['wpthlp2'],
                'cam': ['wpthlp2'],
                'wrf': []
            }, 'sci_scale': 0},

            {'var_names': {
                'clubb': ['wprtp2'],
                'sam': ['WPRTP2'],
                'coamps': ['wpqtp2'],
                'r408': ['wprtp2'],
                'hoc': ['wprtp2'],
                'e3sm': ['wprtp2'],
                'cam': ['wprtp2'],
                'wrf': []
            }, 'sci_scale': -7},

            {'var_names': {
                'clubb': ['wprtpthlp'],
                'sam': ['WPRTPTHLP'],
                'coamps': ['wpqtpthlp'],
                'r408': ['wprtpthlp'],
                'hoc': ['wprtpthlp'],
                'e3sm': ['wprtpthlp'],
                'cam': ['wprtpthlp'],
                'wrf': []
            }, 'sci_scale': -4},

            {'var_names': {
                'clubb': ['wp2thvp'],
                'sam': ['WP2THVP'],
                'coamps': ['wp2thvp'],
                'r408': ['wp2thvp'],
                'hoc': ['wp2thvp'],
                'e3sm': ['wp2thvp'],
                'cam': ['wp2thvp'],
                'wrf': []
            }, 'sci_scale': 0},
            # TODO LES
        ]
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets)
