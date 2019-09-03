"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupCorrelations(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):

        self.name = "w variables"
        self.variable_definitions = [
            {'aliases': ['corr_w_rr_1']},
            {'aliases': ['corr_w_Nr_1']},
            {'aliases': ['corr_w_Ncn_1']},
            {'aliases': ['corr_chi_rr_1']},
            {'aliases': ['corr_chi_Nr_1']},
            {'aliases': ['corr_chi_Ncn_1']},
            {'aliases': ['corr_rr_Nr_1']},

        ]
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)