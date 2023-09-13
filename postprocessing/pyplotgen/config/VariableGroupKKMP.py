"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupKKMP(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None,
                 coamps_benchmark_dataset=None, r408_dataset=None, wrf_benchmark_dataset=None, 
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False,
                 background_rcm=False, background_rcm_folder=None):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "kk mp variables"
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': ['rrm_evap'],
                'sam': ['EVAPM'],
                'coamps': ['rrm_evap'],
                'r408': ['rrm_evap'],
                'hoc': ['rrm_evap'],
                'e3sm': ['rrm_evap'],
                'cam': ['rrm_evap'],
                'wrf': [],
                },
            },
            {'var_names':
                {
                'clubb': ['rrm_accr'],
                'sam': ['ACCRM'],
                'coamps': ['rrm_accr'],
                'r408': ['rrm_accr'],
                'hoc': ['rrm_accr'],
                'e3sm': ['rrm_accr'],
                'cam': ['rrm_accr'],
                'wrf': [],
                },
            },
            {'var_names':
                {
                'clubb': ['rrm_auto'],
                'sam': ['AUTOM'],
                'coamps': ['rrm_auto'],
                'r408': ['rrm_auto'],
                'hoc': ['rrm_auto'],
                'e3sm': ['rrm_auto'],
                'cam': ['rrm_auto'],
                'wrf': [],
                },
            },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset, coamps_benchmark_dataset=coamps_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         wrf_benchmark_dataset=wrf_benchmark_dataset, 
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)
