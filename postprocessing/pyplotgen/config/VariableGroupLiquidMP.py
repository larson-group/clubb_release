'''
:author: Nicolas Strike
:date: Mid 2019
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupLiquidMP(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None, hoc_dataset=None,
                 e3sm_datasets=None):
        """

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "liquid mp variables"
        self.variable_definitions = [
            {'var_names': ['Ncm'], 'sam_calc': self.getNcmSamLine},
            {'var_names': ['Nc_in_cloud']},
            {'var_names': ['precip_frac']},
            {'var_names': ['rrm', 'QPL'], 'sam_conv_factor': 1 / 1000},
            {'var_names': ['Nrm'], 'sam_calc': self.getNrmSamLine},
            {'var_names': ['wprrp', 'WPRRP']},  # Not found in lba case file
            {'var_names': ['wpNrp', 'WPNRP']},  # Not found in lba case file
            {'var_names': ['rwp', 'RWP'], 'sam_conv_factor': 1 / 1000, 'type': Panel.TYPE_TIMESERIES},
            {'var_names': ['precip_rate_sfc'], 'type': Panel.TYPE_TIMESERIES}

        ]
        #rain_rate_sfc vs time

        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets= e3sm_datasets)


    def getNcmSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NC * 1e+6) ./ RHO
        :return:
        """
        nc,z, dataset = self.getVarForCalculations('NC', self.sam_file, fill_zeros=True)
        rho,z, dataset = self.getVarForCalculations('RHO', self.sam_file)
       # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        ncm = (nc * (10 ** 6) / rho)
        return ncm, z

    def getNrmSamLine(self):
        """
        Caclulates Nim from sam -> clubb using the equation
        (NR * 1e+6) ./ RHO
        :return:
        """
        nr,z, dataset = self.getVarForCalculations('NR', self.sam_file, fill_zeros=True)
        rho,z, dataset = self.getVarForCalculations('RHO', self.sam_file)
       # z,z, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_file)

        nrm = (nr * (10 ** 6) / rho)

        return nrm, z