"""
:author: Nicolas Strike
:date: Mid 2019
"""

from VariableGroup import VariableGroup


class VariableGroupKKMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "kk mp variables"
        self.variable_definitions = [
            {'aliases': ['rrm_cond', 'EVAPM']},
            {'aliases': ['rrm_accr', 'ACCRM']},
            {'aliases': ['rrm_auto', 'AUTOM']}
        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)
