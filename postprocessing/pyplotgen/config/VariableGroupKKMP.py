"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupKKMP(VariableGroup):

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "kk mp variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['rrm_evap'],
                'sam': ['EVAPM'],
                'coamps': ['rrm_evap'],
                'r408': ['rrm_evap'],
                'hoc': ['rrm_evap'],
                'e3sm': ['rrm_evap'],
                'wrf': []
            }},

            {'var_names': {
                'clubb': ['rrm_accr'],
                'sam': ['ACCRM'],
                'coamps': ['rrm_accr'],
                'r408': ['rrm_accr'],
                'hoc': ['rrm_accr'],
                'e3sm': ['rrm_accr'],
                'wrf': []
            }},

            {'var_names': {
                'clubb': ['rrm_auto'],
                'sam': ['AUTOM'],
                'coamps': ['rrm_auto'],
                'r408': ['rrm_auto'],
                'hoc': ['rrm_auto'],
                'e3sm': ['rrm_auto'],
                'wrf': []
            }}

        ]
        super().__init__(case, clubb_datasets=clubb_datasets, les_dataset=les_dataset, coamps_dataset=coamps_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         sam_datasets=sam_datasets, wrf_datasets=wrf_datasets)
