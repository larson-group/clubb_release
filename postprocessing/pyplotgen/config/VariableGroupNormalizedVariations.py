"""
:author: Jan Gruenenwald
:date: 04/2022
"""

from src.VariableGroup import VariableGroup
from src.Panel import Panel


class VariableGroupNormalizedVariations(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "normalized variation variables"
        self.variable_definitions = [
            {'var_names':
                 {
                 'clubb': ['tot_vartn_normlzd_rtm'],
                 'sam': [''],
                 'coamps': [''],
                 'r408': [''],
                 'hoc': [''],
                 'e3sm': [''],
                 'cam': [''],
                 'wrf': [''],
                 },
                 'type': Panel.TYPE_TIMESERIES,
                 'priority': False,
             },
             {'var_names':
                  {
                  'clubb': ['tot_vartn_normlzd_thlm'],
                  'sam': [''],
                  'coamps': [''],
                  'r408': [''],
                  'hoc': [''],
                  'e3sm': [''],
                  'cam': [''],
                  'wrf': [''],
                  },
                  'type': Panel.TYPE_TIMESERIES,
                  'priority': False,
              },
              {'var_names':
                   {
                   'clubb': ['tot_vartn_normlzd_wprtp'],
                   'sam': [''],
                   'coamps': [''],
                   'r408': [''],
                   'hoc': [''],
                   'e3sm': [''],
                   'cam': [''],
                   'wrf': [''],
                   },
                   'type': Panel.TYPE_TIMESERIES,
                   'priority': False,
               }
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars)
 
