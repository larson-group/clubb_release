'''
:author: Nicolas Strike
:date: Mid 2019
'''
from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupLiquidMP(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None,priority_vars=False):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "liquid mp variables"
        self.variable_definitions = [
            {'var_names':
                {
                    'clubb': ['Ncm', 'ncm'],
                    'sam': [self.getNcmSamLine, 'Ncm', 'ncm'],
                    'coamps': ['Ncm', 'ncm'],
                    'r408': ['Ncm', 'ncm'],
                    'hoc': ['Ncm', 'ncm'],
                    'e3sm': ['Ncm', 'ncm'],
                    'cam': ['AWNC', 'Ncm', 'ncm'],
                    'wrf': ['Ncm', 'ncm'],
                },
                # 'sam_calc': self.getNcmSamLine
            },
            {'var_names':
                {
                    'clubb': ['Nc_in_cloud'],
                    'sam': [self.getNcInCloudSamLine, 'Nc_in_cloud'],
                    'coamps': ['Nc_in_cloud'],
                    'r408': ['Nc_in_cloud'],
                    'hoc': ['Nc_in_cloud'],
                    'e3sm': ['Nc_in_cloud'],
                    'cam': ['Nc_in_cloud'],
                    'wrf': ['Nc_in_cloud'],
                },
                # 'sam_calc': self.getNcInCloudSamLine
            },
            {'var_names':
                {
                    'clubb': ['precip_frac'],
                    'sam': ['precip_frac'],
                    'coamps': ['precip_frac'],
                    'r408': ['precip_frac'],
                    'hoc': ['precip_frac'],
                    'e3sm': ['precip_frac'],
                    'cam': ['precip_frac'],
                    'wrf': ['precip_frac'],
                },
            },
            {'var_names':
                {
                    'clubb': ['rrm'],
                    'sam': ['QPL'],
                    'coamps': ['rrm'],
                    'r408': ['rrm'],
                    'hoc': ['rrm'],
                    'e3sm': ['rrm'],
                    'cam': ['AQRAIN', 'rrm'],
                    'wrf': ['rrm'],
                },
                'sam_conv_factor': 1 / 1000,
                'sci_scale': -6
            },
            {'var_names':
                {
                    'clubb': ['Nrm'],
                    'sam': [self.getNrmSamLine, 'Nrm'],
                    'coamps': ['Nrm'],
                    'r408': ['Nrm'],
                    'hoc': ['Nrm'],
                    'e3sm': ['Nrm'],
                    'cam': ['ANRAIN', 'Nrm'],
                    'wrf': ['Nrm'],
                },
                # 'sam_calc': self.getNrmSamLine,
            },
            {'var_names':
                {
                    'clubb': ['wprrp'],
                    'sam': ['WPRRP'],
                    'coamps': ['wprrp'],
                    'r408': ['wprrp'],
                    'hoc': ['wprrp'],
                    'e3sm': ['wprrp'],
                    'cam': ['wprrp'],
                    'wrf': ['wprrp'],
                },
            },
            # Not found in lba case file
            {'var_names':
                {
                    'clubb': ['wpNrp'],
                    'sam': ['WPNRP'],
                    'coamps': ['wpNrp'],
                    'r408': ['wpNrp'],
                    'hoc': ['wpNrp'],
                    'e3sm': ['wpNrp'],
                    'cam': ['wpNrp'],
                    'wrf': ['wpNrp'],
                },
            },
            # Not found in lba case file
            {'var_names':
                {
                    'clubb': ['rwp'],
                    'sam': ['RWP'],
                    'coamps': ['rwp'],
                    'r408': ['rwp'],
                    'hoc': ['rwp'],
                    'e3sm': ['rwp'],
                    'cam': ['rwp'],
                    'wrf': ['rwp'],
                },
                'sam_conv_factor': 1 / 1000, 'type': Panel.TYPE_TIMESERIES,
            },
            {'var_names':
                {
                    'clubb': ['precip_rate_sfc'],
                    'sam': ['PREC','precip_rate_sfc'],
                    'coamps': ['precip_rate_sfc'],
                    'r408': ['precip_rate_sfc'],
                    'hoc': ['precip_rate_sfc'],
                    'e3sm': ['precip_rate_sfc'],
                    'cam': ['precip_rate_sfc'], # PRECT.*1000.*86400
                    'wrf': ['precip_rate_sfc'],
                },
                'type': Panel.TYPE_TIMESERIES
            },
        ]
        # rain_rate_sfc vs time

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset,wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars)

    def getNcmSamLine(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Caclulates Nim from sam -> clubb using the equation

        .. code-block:: python
            :linenos:

            (NC * 1e+6) ./ RHO
            (GCSSNC * 1e+6) ./ RHO
            CLD .* (NC * 1e+6) ./ RHO

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """

        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        nc, indep, dataset = self.getVarForCalculations(['NC'], dataset)
        gcssnc, indep, dataset = self.getVarForCalculations(['GCSSNC'], dataset)
        rho, indep, dataset = self.getVarForCalculations('RHO', dataset)
        cld, indep, dataset = self.getVarForCalculations('CLD', dataset)

        output1 = cld * (nc * (10 ** 6) / rho) #ncm
        output2 = (gcssnc * (10 ** 6) / rho) #nc in cloud

        output = self.pickNonZeroOutput(output1, output2)
        return output, indep

    def getNrmSamLine(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Caclulates Nim from sam -> clubb using the equation

        .. code-block:: python
            :linenos:

            (NR * 1e+6) ./ RHO
            (CONP * 1e+6) ./ RHO

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """

        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override

        nr, indep, dataset = self.getVarForCalculations(['NR','CONP'], dataset)
        rho, indep, dataset = self.getVarForCalculations('RHO', dataset)
        # indep, indep, dataset = self.getVarForCalculations(['z', 'lev', 'altitude'], self.sam_benchmark_dataset)

        nrm = (nr * (10 ** 6) / rho)

        return nrm, indep

    def getNcInCloudSamLine(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        ``(NC * 1e+6) ./ RHO``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        nc, indep, dataset = self.getVarForCalculations(['NC'], dataset)
        rho, indep, dataset = self.getVarForCalculations('RHO', dataset)

        output = (nc * (10 ** 6) / rho)

        return output, indep
