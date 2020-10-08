'''
:author: Nicolas Strike
:date: Mid 2019
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupSubcolumns(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """

        """
        self.name = "Subcolumn variables"
        self.variable_definitions = [
            {'var_names':
                {
                    'clubb': ['w'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'w'
            },
            {'var_names':
                {
                    'clubb': ['chi'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'chi'
            },
            {'var_names':
                {
                    'clubb': ['eta'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'eta'
            },
            {'var_names':
                {
                    'clubb': ['Ncn'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'Ncn'
            },

            {'var_names':
                {
                    'clubb': ['rr'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'rr'
            },
            {'var_names':
                {
                    'clubb': ['Nr'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'Nr'
            },
            {'var_names':
                {
                    'clubb': ['ri'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'ri'
            },
            {'var_names':
                {
                    'clubb': ['Ni'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'Ni'
            },
            {'var_names':
                {
                    'clubb': ['rg'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'rg'
            },
            {'var_names':
                {
                    'clubb': ['Ng'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'Ng'
            },
            {'var_names':
                {
                    'clubb': ['rs'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'rs'
            },
            {'var_names':
                {
                    'clubb': ['Ns'],
                    'sam': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                    'wrf': [],
                },
                'type': Panel.TYPE_SUBCOLUMN, 'legend_label': 'Ns'
            },

        ]

        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets)
