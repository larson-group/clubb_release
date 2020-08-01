"""
:author: Nicolas Strike
:date: Mid 2020
"""

from src.VariableGroup import VariableGroup


class VariableGroupScalars(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "scalar variables"
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['sclr1m'],
                'sam': [],
                'coamps': ['sclr1m'],
                'r408': ['sclr1m'],
                'hoc': ['sclr1m'],
                'e3sm': ['sclr1m'],
                'cam': ['sclr1m'],
                'wrf': ['sclr1m'],
                },
            },
            {'var_names':
                {
                    'clubb': ['sclr2m'],
                    'sam': [],
                    'coamps': ['sclr2m'],
                    'r408': ['sclr2m'],
                    'hoc': ['sclr2m'],
                    'e3sm': ['sclr2m'],
                    'cam': ['sclr2m'],
                    'wrf': ['sclr2m'],
                },
            },
            {'var_names':
                {
                    'clubb': ['sclr1p2'],
                    'sam': [],
                    'coamps': ['sclr1p2'],
                    'r408': ['sclr1p2'],
                    'hoc': ['sclr1p2'],
                    'e3sm': ['sclr1p2'],
                    'cam': ['sclr1p2'],
                    'wrf': ['sclr1p2'],
                },
            },
            {'var_names':
                {
                    'clubb': ['sclr2p2'],
                    'sam': [],
                    'coamps': ['sclr2p2'],
                    'r408': ['sclr2p2'],
                    'hoc': ['sclr2p2'],
                    'e3sm': ['sclr2p2'],
                    'cam': ['sclr2p2'],
                    'wrf': ['sclr2p2'],
                },
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets)
