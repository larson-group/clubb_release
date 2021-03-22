# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 14:57:19 2021

"""

def main():

    import numpy as np
    import pdb
    from analyze_sensitivity_matrix \
        import analyzeSensMatrix#, plotNormlzdSensMatri

    # Metrics are observed quantities that we want to match.
    metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

    # Parameters are tunable model parameters.
    #paramsNames = np.array(['C5','C8'])
    paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

    transformedParams = np.array(['clubb_c8']) #np.array(['clubb_c8'])

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    # Global mean, NOTE PRECT is in the unit of m/s
    obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

    # Create fake testing data with simple values
    write_test_netcdf_files(obsMetricValsDict)

    # Filenames of netcdf files containing fake testing data
    sensNcFilenames = np.array(['sens1.nc', 'sens2.nc'])
    defaultNcFilename = 'default.nc'

    # Calculate changes in parameter values needed to match metrics.
    sensMatrix, normlzdSensMatrix, svdInvrsNormlzd, dparamsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParams,
                      sensNcFilenames, defaultNcFilename,
                       obsMetricValsDict)

    # Check whether the expected answer for the fake data
    #    has been calculated correctly.
    if np.all( np.isclose( dparamsSoln, np.array([[0.48658305], [2./3.]]), \
                              rtol=1e-5, atol=1e-5 ) ):
        print("\nPassed test.")
    else:
        print("\ndparamsSoln = ")
        print(dparamsSoln)
        print("\nError: dparamsSoln should equal [0.48658305 2/3], but it does not.")

    # Create a heatmap matrix that displays the normalized sensitivity matrix
#    plotNormlzdSensMatrix(normlzdSensMatrix, metricsNames, paramsNames)

    print("\nReached the end of main in test harness.")

    return


def write_test_netcdf_files(obsMetricValsDict):

    import netCDF4

    SWCF_obs = obsMetricValsDict['SWCF']
    LWCF_obs = obsMetricValsDict['LWCF']
    PRECT_obs = obsMetricValsDict['PRECT']

    f_default = netCDF4.Dataset('default.nc', 'w')
    f_default.createDimension('scalar', 1)
    SWCF = f_default.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 1.0
    LWCF = f_default.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs - 1.0
    PRECT = f_default.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs - 1.0
    clubb_c8 = f_default.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = 0.1812692 #0.2
    clubb_c_invrs_tau_n2 = f_default.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 3.
    f_default.close()

    f_sens1 = netCDF4.Dataset('sens1.nc', 'w')
    f_sens1.createDimension('scalar', 1)
    SWCF = f_sens1.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 0.5
    LWCF = f_sens1.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs - 1.
    PRECT = f_sens1.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs - 0.5
    clubb_c8 = f_sens1.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = 0.50341469 #0.7
    clubb_c_invrs_tau_n2 = f_sens1.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 3.
    f_sens1.close()

    f_sens2 = netCDF4.Dataset('sens2.nc', 'w')
    f_sens2.createDimension('scalar', 1)
    SWCF = f_sens2.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 1.
    LWCF = f_sens2.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs + 1.
    PRECT = f_sens2.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs + 1.
    clubb_c8 = f_sens2.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = 0.1812692 #0.2
    clubb_c_invrs_tau_n2 = f_sens2.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 5.
    f_sens2.close()

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()

