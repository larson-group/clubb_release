'''
:author: Nicolas Strike
:date: Mid 2019
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupLiquidMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, les_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None, sam_datasets=None):
        """

        :param ncdf_datasets:
        :param case:
        :param les_file:
        """
        self.name = "liquid mp variables"
        self.variable_definitions = [
            {'var_names': {
                'clubb': ['Ncm'],
                'sam': ['Ncm'],
                'coamps': ['Ncm'],
                'r408': ['Ncm'],
                'hoc': ['Ncm'],
                'e3sm': ['Ncm']
            },
                'sam_calc': self.getNcmSamLine},
            {'var_names': {
                'clubb': ['Nc_in_cloud'],
                'sam': ['Nc_in_cloud'],
                'coamps': ['Nc_in_cloud'],
                'r408': ['Nc_in_cloud'],
                'hoc': ['Nc_in_cloud'],
                'e3sm': ['Nc_in_cloud']
            }},

            {'var_names': {
                'clubb': ['precip_frac'],
                'sam': ['precip_frac'],
                'coamps': ['precip_frac'],
                'r408': ['precip_frac'],
                'hoc': ['precip_frac'],
                'e3sm': ['precip_frac']
            }},

            {'var_names': {
                'clubb': ['rrm'],
                'sam': ['QPL'],
                'coamps': ['rrm'],
                'r408': ['rrm'],
                'hoc': ['rrm'],
                'e3sm': ['rrm']
            },
                'sam_conv_factor': 1 / 1000},
            {'var_names': {
                'clubb': ['Nrm'],
                'sam': ['Nrm'],
                'coamps': ['Nrm'],
                'r408': ['Nrm'],
                'hoc': ['Nrm'],
                'e3sm': ['Nrm']
            },
                'sam_calc': self.getNrmSamLine},
            {'var_names': {
                'clubb': ['wprrp'],
                'sam': ['WPRRP'],
                'coamps': ['wprrp'],
                'r408': ['wprrp'],
                'hoc': ['wprrp'],
                'e3sm': ['wprrp']
            }},
            # Not found in lba case file
            {'var_names': {
                'clubb': ['wpNrp'],
                'sam': ['WPNRP'],
                'coamps': ['wpNrp'],
                'r408': ['wpNrp'],
                'hoc': ['wpNrp'],
                'e3sm': ['wpNrp']
            }},
            # Not found in lba case file
            {'var_names': {
                'clubb': ['rwp'],
                'sam': ['RWP'],
                'coamps': ['rwp'],
                'r408': ['rwp'],
                'hoc': ['rwp'],
                'e3sm': ['rwp']
            },
                'sam_conv_factor': 1 / 1000, 'type': Panel.TYPE_TIMESERIES},
            {'var_names': {
                'clubb': ['precip_rate_sfc'],
                'sam': ['precip_rate_sfc'],
                'coamps': ['precip_rate_sfc'],
                'r408': ['precip_rate_sfc'],
                'hoc': ['precip_rate_sfc'],
                'e3sm': ['precip_rate_sfc']
            },
                'type': Panel.TYPE_TIMESERIES}

        ]
        # rain_rate_sfc vs time

        super().__init__(ncdf_datasets, case, les_file=les_file, coamps_file=coamps_file, r408_dataset=r408_dataset,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, sam_datasets=sam_datasets)

    def getNcmSamLine(self, dataset_override=None):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NC * 1e+6) ./ RHO
        :return:
        """

        dataset = self.les_file
        if dataset_override != None:
            dataset = dataset_override
        nc, z, dataset = self.getVarForCalculations('NC', dataset)
        rho, z, dataset = self.getVarForCalculations('RHO', dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_file)

        ncm = (nc * (10 ** 6) / rho)
        return ncm, z

    def getNrmSamLine(self, dataset_override=None):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NR * 1e+6) ./ RHO
        :return:
        """

        dataset = self.les_file
        if dataset_override != None:
            dataset = dataset_override

        nr, z, dataset = self.getVarForCalculations('NR', dataset)
        rho, z, dataset = self.getVarForCalculations('RHO', dataset)
        # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.les_file)

        nrm = (nr * (10 ** 6) / rho)

        return nrm, z
