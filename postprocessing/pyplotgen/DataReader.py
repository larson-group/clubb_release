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
import pathlib as pathlib

class NetCdfVariable:
    '''
    Class used for conviniently storying the information about a given netcdf variable
    '''
    def __init__(self, name, ncdf_file, conversion_factor = 1, avging_start_time = -1, avging_end_time = -1, timestep = -1, avg_axis = 0):
        '''

        :param name:
        :param ncdf_file:
        :param conversion_factor:
        :param avging_start_time:
        :param avging_end_time:
        :param timestep:
        :param avg_axis: The axis to avg data over. 0 for time-avg, 1 for height avg
        '''
        data_reader = DataReader()
        self.name = name
        self.start_time = avging_start_time
        self.end_time = avging_end_time
        self.timestep = timestep
        self.conv_factor = conversion_factor
        self.avg_axis = avg_axis
        self.data = data_reader.getVarData(ncdf_file, self)

    def constrain(self, min_value, max_value):
        '''
        Remove everything in the data from before min_value
        and after max_value

        :param min_value:
        :param max_value:
        :return:
        '''
        start_idx, end_idx = self.__getStartEndIndex__(self.data, min_value, max_value)
        self.data = self.data[start_idx:end_idx]

    def __getStartEndIndex__(self, data, start_value, end_value):
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
        self.root_dir =  pathlib.Path(__file__).parent
        self.panels_dir = self.root_dir.as_uri() + "/cases/panels/"

    def cleanup(self):
        '''
        This is the cleanup method. This is called on the instance's destruction
        to unallocate resources that may be held (e.g. dataset files).
        :return:
        :author: Nicolas Strike
        '''
        for dataset in self.nc_datasets:
            dataset.close()

    def loadFolder(self, folder_path, ignore_git = True):
        '''
        Finds all dataset files in a given folder and loads
        them using the appropriate helper class.
        :param folder_path: The path of the folder to be loaded
        :param ignore_git: Ignore files and paths that contain '.git' in their name
        :return: An array of datasets
        :author: Nicolas Strike
        '''
        for root, dirs, files in os.walk(folder_path):
            for filename in files:
                abs_filename = os.path.abspath(os.path.join(root, filename))
                file_ext = os.path.splitext(filename)[1]
                if ignore_git and '.git' in abs_filename:
                    continue
                if file_ext == ".nc":
                    self.nc_filenames.append(abs_filename)
                    self.nc_datasets.append(self.__loadNcFile__(abs_filename))
                elif file_ext == ".dat":
                    self.grads_dat_filenames.append(abs_filename)
                elif file_ext == ".ctl":
                    self.grads_ctl_filenames.append(abs_filename)
                else:
                    print("Filetype " + file_ext + " is not supported. Attempted to load " + abs_filename)
        return self.nc_datasets

    def getVarData(self, netcdf_data, variable: NetCdfVariable):
        '''

        :author: Nicolas Strike
        '''

        # TODO load these, don't hardcode them
        start_time_value = variable.start_time
        end_time_value = variable.end_time
        variable_name = variable.name
        conv_factor = variable.conv_factor
        avg_axis = variable.avg_axis
        level_amount = 1
        num_timesteps = 1
        time_values = self.__getValuesFromNc__(netcdf_data, "time", 1, 1, 1) #TODO conversion shouldn't be only 1 value
        (start_avging_index, end_avging_idx) = self.__getStartEndIndex__(time_values, start_time_value, end_time_value) # Get the index values that correspond to the desired start/end x values

        values = self.__getValuesFromNc__(netcdf_data, variable_name, conv_factor, 1, 1) #TODO conversion shouldn't be only 1 value\
        if values.ndim > 1:#not variable.one_dimentional:
            values = self.__meanProfiles__(values, start_avging_index, end_avging_idx, avg_axis=avg_axis)

        return values


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

    def __loadNcFile__(self, filename):
        '''
        Load the given NetCDF file
        :param filename: the netcdf file to be loaded
        :return: a netcdf dataset containing the data from the given file
        '''
        dataset = Dataset(filename, "r+", format="NETCDF4")
        return dataset

    def __meanProfiles__(self, var, idx_t0 = 0, idx_t1 = -1, idx_z0 = 0, idx_z1 = -1, avg_axis=0):
        # logger.info('mean_profiles')
        """
        Input:
          var    -- time x height array of some property
          idx_t0 -- Index corrosponding to the beginning of the averaging interval
          idx_t1 -- Index corrosponding to the end of the averaging interval
          idx_z0 -- Index corrosponding to the lowest model height of the averaging interval
          idx_z1 -- Index corrosponding to the highest model height of the averaging interval

        Output:
          var    -- time averaged vertical profile of the specified variable
        """
        if idx_t1 is -1:
            idx_t1 = len(var)

        # Nic changed nanmean() to mean()
        if avg_axis is 0: # if time-averaged
            var_average = np.mean(var[idx_t0:idx_t1,:],axis=avg_axis)
        else: # if height averaged
            var_average = np.mean(var[:,idx_z0:idx_z1],axis=avg_axis)
        return var_average

    def __getUnits__(self, nc, varname):
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

    def __getLongName__(self, nc, varname):
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

    def __getValuesFromNc__(self, nc, varname, conversion, level_amount, num_timesteps):
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
        #varname = varname.upper()
        if varname in keys:
            # logger.debug('%s is in keys', varname)
            var = nc.variables[varname]
            var = np.squeeze(var)
            var = var*conversion
        else:
            # logger.debug('%s is not in keys', varname)
            # var = np.zeros(shape=(level_amount, num_timesteps)) - 1.
            raise ValueError("Variable " + varname + " does not exist in nc file")
        return var

    def __getStartEndIndex__(self, data, start_value, end_value):
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
