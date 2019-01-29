import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


class Plotter:

    def __init__(self):
        '''
        Initialize a plotter object TODO
        '''



    def plot_figure(self, case, variable):
        '''
        Plot a single figure/panel. This is 1 graph but
        may include multiple sets of data (e.g. multiple lines)
        :param case:
        :param data:
        :return:
        '''
        # title = data.
        figure = plt.figure(1)
        plt.subplot(111)
        plt.plot(thlm_time[11])


    def plot_case(self, case, data):
        '''
        Plots a case onto a set of graphs (panels)
        :param case: The case to be plotted
        :param data: The data to plot with the case
        :return:
        '''
        # print("\n\nVariables:\n")
        # for key in data.variables.keys():
        #     print("\t" + key)

        #*** From Stackoverflow here https://stackoverflow.com/questions/45582344/extracting-data-from-netcdf-by-python
        thlm_values = self.get_budgets_from_nc(data, "thlm", 1.0,1,1)
        height_values = self.get_budgets_from_nc(data, "altitude", 1.0,1,1)


        # start_value = 181
        # end_value = 240
        thlm_values = self.mean_profiles(thlm_values,0,0,0,0)
        start_index, end_idx = self.getStartEndIndex(thlm_values, 281, 310)

        thlm_values = thlm_values[start_index:end_idx + 1]
        height_values = height_values[start_index: end_idx + 1]

        # thlm_values = self.mean_profiles(thlm_values, 181, 240, start_index, end_idx)
        # lon = data.variables['longitude']
        # lat = data.variables['latitude']
        # time = data.variables['time']
        # altitudes = data.variables['altitude']

        # lon_array = lon[:]
        # lat_array = lat[:]
        # time_array = time[:]
        # altitude_array = altitudes[:]

        # thlm = data.variables["thlm"]

        # i = np.abs(lon_array - 10).argmin()
        # j = np.abs(lat_array - 30).argmin()
        # altitude_data = np.abs(altitude_array - 20).argmin()
        # t = np.abs(time_array).argmin()

        # thlm_values = thlm[:,:,i,j]
        #print("\n\nthlm_values\n" + str(thlm_values))
        #************************************************************************************************************
        #print("\nTest data: " + str(data.variables["thlm"]))

        plt.figure(1)
        plt.subplot(111)
        plt.title(self.get_long_name(data, 'thlm'))
        plt.plot(thlm_values[:], height_values[:], 'r--')
        plt.ylabel("Height[m]")
        plt.xlabel("thlm [K]")
        plt.show()

    def getStartEndIndex(self, data, start_value, end_value):
        '''
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending
        index representing the start and end time passed into the function
        '''
        start_idx = 0
        end_idx = len(data) -1
        for i in range(0,len(data)):
                # Check for start index
                test_value = data[i]
                if test_value <= start_value and test_value > data[start_idx]:
                    start_idx = i
                # Check for end index
                if test_value >= end_value and test_value < data[end_idx]:
                    end_idx = i

        return start_idx, end_idx


    def mean_profiles(self, var, idx_t0, idx_t1, idx_z0, idx_z1):
        # logger.info('mean_profiles')
        """
        Input:
          var    -- time x height array of some property
          idx_t0 -- Index corrosponding to the beginning of the averaging interval
          idx_t1 -- Index corrosponding to the end of the averaging interval
          idx_z0 -- Index corrosponding to the lowest model level of the averaging interval
          idx_z1 -- Index corrosponding to the highest model level of the averaging interval

        Output:
          var    -- time averaged vertical profile of the specified variable
        """
        #TODO Nic changed nanmean() to mean()
        var_average = np.mean(var[:,:],axis=0)
        return var_average

    def get_units(self, nc, varname):
        # logger.info('get_units:%s', varname)
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          unit as string
        """

        keys = nc.variables.keys()
        if varname in keys:
            # logger.debug('%s is in keys', varname)
            unit = nc.variables[varname].units
        else:
            # logger.debug('%s is not in keys', varname)
            unit = "nm"

        return unit

    def get_long_name(self, nc, varname):
        # logger.info('get_long_name:%s', varname)
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          long_name as string
        """

        keys = nc.variables.keys()
        if varname in keys:
            # logger.debug('%s is in keys', varname)
            long_name = nc.variables[varname].long_name
        else:
            # logger.debug('%s is not in keys', varname)
            long_name = "nm"

        return long_name

    def get_budgets_from_nc(self, nc, varname, conversion, n, t):
        # logger.info('get_budgets_from_nc:%s', varname)
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
          conversion --  Conversion factor
          n          --  amount of level
          t          --  amount of timesteps
          n and t are used, if the variable cannot be found

        Output:
          time x height array of the specified variable, scaled by conversion factor
        """

        keys = nc.variables.keys()
        if varname in keys:
            # logger.debug('%s is in keys', varname)
            var = nc.variables[varname]
            var = np.squeeze(var)
            var = var*conversion
        else:
            # logger.debug('%s is not in keys', varname)
            var = np.zeros(shape=(n,t)) - 1000.
        return var