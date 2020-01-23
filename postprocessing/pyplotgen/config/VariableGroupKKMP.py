"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupKKMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "kk mp variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['rrm_evap'],
                'sam': ['EVAPM'],
                'coamps': ['rrm_evap'],
                'r408': ['rrm_evap'],
                'hoc': ['rrm_evap'],
                'e3sm': ['rrm_evap']
            }},

            {'var_names': {
                'clubb': ['rrm_accr'],
                'sam': ['ACCRM'],
                'coamps': ['rrm_accr'],
                'r408': ['rrm_accr'],
                'hoc': ['rrm_accr'],
                'e3sm': ['rrm_accr']
            }},

            {'var_names': {
                'clubb': ['rrm_auto'],
                'sam': ['AUTOM'],
                'coamps': ['rrm_auto'],
                'r408': ['rrm_auto'],
                'hoc': ['rrm_auto'],
                'e3sm': ['rrm_auto']
            }}

        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets)
