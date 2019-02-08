'''
This file is used to load data in from netcdf and grads
formatted files and provide a Plot class instance to
a Plotter instance for graphing.

:author: Nicolas Strike
:date: Early 2019
'''

import os
from netCDF4 import Dataset
import numpy as np
from collections import OrderedDict
from collections import namedtuple

class DataReader():
    '''
    This class is responsible for handling input files. Given a
    input file, it determines what format it is in (grads vs
    netcdf) and then uses the correct helper file (e.g.
    NetcdfReader.py) to load the data into python.

    :author: Nicolas Strike
    :data: January 2019
    '''


    def __init__(self):
        '''
        This is the object constructer. Initalizes class variables and data
        :author: Nicolas Strike
        '''
        self.nc_filenames = []
        self.grads_dat_filenames = []
        self.grads_ctl_filenames = []
        self.nc_datasets = []
        self.PlotData = namedtuple("PlotData", "x_values y_values title x_title y_title")

    def cleanup(self):
        '''
        This is the cleanup method. This is called on the instance's destruction
        to unallocate resources that may be held (e.g. dataset files).
        :return:
        :author: Nicolas Strike
        '''
        for dataset in self.nc_datasets:
            dataset.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''
        Calls the cleanup and cleanly closes out the object isntance
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        :author: Nicolas Strike
        '''
        self.cleanup()

    def load_nc_file(self, filename):
        '''
        Load the given NetCDF file
        :param filename: the netcdf file to be loaded
        :return: a netcdf dataset containing the data from the given file
        '''
        dataset = Dataset(filename, "r+", format="NETCDF4")
        # print("\n\nVariables:\n")
        # for key in dataset.variables.keys():
        #     print("\t" + key)
        #print(filename + ": " + str(dataset))
        return dataset


    def load_folder(self, folder_path):
        '''
        Finds all dataset files in a given folder and loads
        them using the appropriet helper class.
        :param folder_path: The path of the folder to be loaded
        :return: An array of datasets
        :author: Nicolas Strike
        '''
        for root, dirs, files in os.walk(folder_path):
            for filename in files:
                abs_filename = os.path.abspath(os.path.join(root, filename))
                file_ext = os.path.splitext(filename)[1]
                if file_ext == ".nc":
                    self.nc_filenames.append(abs_filename)
                    self.nc_datasets.append(self.load_nc_file(abs_filename))
                elif file_ext == ".dat":
                    self.grads_dat_filenames.append(abs_filename)
                elif file_ext == ".ctl":
                    self.grads_ctl_filenames.append(abs_filename)
                else:
                    print("Filetype " + file_ext + " is not supported. Attempted to load " + abs_filename)
        return self.nc_datasets
        # print("Files loaded:\n\n----nc files----\n" + str(self.nc_filenames) + "\n\n----dat files----\n" + str(self.grads_dat_filenames) +
        #       "\n\n----ctl files----\n" + str(self.grads_ctl_filenames))


    def mean_profiles(self, var):
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
        # Nic changed nanmean() to mean()
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

    def get_values_from_nc(self, nc, varname, conversion, level_amount, num_timesteps):
        # logger.info('get_budgets_from_nc:%s', varname)
        """
        Get data values out of a netcdf object, returning them as an array

        :param nc: Netcdf file object
        :param varname: Variable name string
        :param conversion: Conversion factor
        :param level_amount: amount of level
        :param  num_timesteps: amount of timesteps level_amount and num_timesteps are used, if the variable cannot be found

        :return: time x height array of the specified variable, scaled by conversion factor
        """

        keys = nc.variables.keys()
        if varname in keys:
            # logger.debug('%s is in keys', varname)
            var = nc.variables[varname]
            var = np.squeeze(var)
            var = var*conversion
        else:
            # logger.debug('%s is not in keys', varname)
            var = np.zeros(shape=(level_amount, num_timesteps)) - 1000.
        return var

    def getStartEndIndex(self, data, start_value, end_value):
        '''
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending index representing the start and end time passed into the function
        :author: Nicolas Strike
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


    def getPlotsData(self, netcdf_data, case):
        '''
        Create a plot tuple containing the data needed to
        create a graph.

        Plot tuple definition:
        Plot(x_values=<array with x values>, y_values=<array with y values>, title='Name of Graph',
            x_title='X Axis Title', y_title='Y Axis Title')

        Example creation:
            myPlot = Plot(x_values=arrayOfXValues, y_values=arrayOfYValues, title='Example Plot', x_title='My X Axis', y_title='My Y Axis')

        :param netcdf_data: The NetCDF data object containing the desired data, e.g. variables to be plotted
        :param case: The case to be plotted (e.g. contains title, data labels, start/end values, etc
        :return: A plot struct containing the data elements listed above
        :author: Nicolas Strike
        '''

        # TODO load these, don't hardcode them
        x_variable_name = "thlm"
        x_conversion_factor = 1.0
        x_level_amount = 1
        x_num_timesteps = 1
        # TODO load these, don't hardcode them
        y_variable_name = "altitude"
        y_conversion_factor = 1.0
        y_level_amount = 1
        y_num_timesteps = 1
        # TODO load these, don't hardcode them
        start_x_value = 281 # Used to determine what x value to begin the graph at
        end_x_value = 310 # Used to determine what x value to end the graph at
        # TODO load these, don't hardcode them
        title = self.get_long_name(netcdf_data, x_variable_name)
        x_axis_title = x_variable_name + "[K]"
        y_axis_title = "Height [m]"


        # Process inspired by Stackoverflow: https://stackoverflow.com/questions/45582344/extracting-data-from-netcdf-by-python
        x_axis_values = self.get_values_from_nc(netcdf_data, x_variable_name, x_conversion_factor, x_level_amount, x_num_timesteps)
        y_axis_values = self.get_values_from_nc(netcdf_data, y_variable_name, y_conversion_factor, y_level_amount, y_num_timesteps)

        x_axis_values = self.mean_profiles(x_axis_values)

        start_index, end_idx = self.getStartEndIndex(x_axis_values, start_x_value, end_x_value) # Get the index values that correspond to the desired start/end x values

        x_axis_values = x_axis_values[start_index:end_idx + 1] # we use end_idx + 1 to ensure python uses the last value (off by one)
        y_axis_values = y_axis_values[start_index: end_idx + 1] # we use end_idx + 1 to ensure python uses the last value (off by one)

        plot_data = self.PlotData(x_values=x_axis_values, y_values=y_axis_values, title=title, x_title=x_axis_title, y_title=y_axis_title)

        return plot_data