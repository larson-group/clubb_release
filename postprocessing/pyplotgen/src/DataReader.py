"""
:author: Nicolas Strike
:date: Early 2019
"""

import os
import pathlib as pathlib
from warnings import warn

import numpy as np
from netCDF4 import Dataset


class NetCdfVariable:
    """
    Class used for conveniently storing the information about a given netcdf variable
    """

    def __init__(self, name, ncdf_data, conversion_factor=1, start_time=0, end_time=-1, avg_axis=0, fill_zeros=False):
        """

        :param name: the name of the variable as defined in the ncdf_data
        :param ncdf_data: Accepts either a Dataset object to pull data from, or a dict of Datasets which will automatically be searched for the variable
        :param conversion_factor: A multiplication factor multiplied to every value in the dataset. Defaults to 1
        :param start_time: The time value to being the averaging period, e.g. 181 minutes. Defaults to 0.
        :param end_time: The time value to stop the averaging period, e.g. 240 minutes. Defaults to -1.
        :param avg_axis: The axis to avg data over. 0 for time-avg, 1 for height avg
        :param fill_zeros: If a variable is not found, setting fill_zeros to True allows pyplotgen to create a fake variable containing 0 for its values.
        """

        # If argument name is a list of aliases, find the correct one and assign it
        # if isinstance(ncdf_data, Dataset):
        #     ncdf_data = {'temp': ncdf_data}
        if isinstance(name, list):
            for tempname in name:
                if tempname in ncdf_data.variables.keys():
                    name = tempname
                    break
        if isinstance(name, list):
            raise NameError("None of the values " + str(name) + " were found in the dataset " + str(ncdf_data))
        data_reader = DataReader()
        self.ncdf_data = ncdf_data
        self.name = name
        self.start_time = start_time
        self.end_time = end_time
        self.conv_factor = conversion_factor
        self.avg_axis = avg_axis
        if isinstance(ncdf_data, dict):
            self.filter_datasets()
        self.data = data_reader.getVarData(self.ncdf_data, self, fill_zeros=fill_zeros)

    def __getStartEndIndex__(self, data, start_value, end_value):
        """
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending index representing the start and end time passed into the function
        :author: Nicolas Strike
        """
        start_idx = 0
        end_idx = len(data) - 1
        for i in range(0, len(data)):
            # Check for start index
            test_value = data[i]
            if test_value <= start_value and test_value > data[start_idx]:
                start_idx = i
            # Check for end index
            if test_value >= end_value and test_value < data[end_idx]:
                end_idx = i

        return start_idx, end_idx + 1

    def constrain(self, start_value, end_value, data=None):
        """
        Remove all data elements from the variable that are not between the start and end value. Assumes
        the data is always increasing. If the optional data parameter is used, it will restrict to the indices
        of where the start/end values are in that dataset rather than the variables data itself.

        :param start_value: The smallest possible value
        :param end_value: The largest possible value
        :param data: A list of data values to find start/ending indicies for constraining rather than using the NetCdfVariable's data.
        :return: None, operates in place
        """
        if data is None:
            data = self.data
        start_idx, end_idx = self.__getStartEndIndex__(data, start_value, end_value)
        self.data = self.data[start_idx:end_idx]


    def filter_datasets(self):
        """
        Looks through the input files for a case and finds the Dataset
        containing the requested variable and assigns it to self.ncdf_data
        Assumes self.ncdf_data is a dictionary.

        :return: None
        """
        for subdataset in self.ncdf_data.values():
            if self.name in subdataset.variables.keys():
                self.ncdf_data = subdataset
                break

    def __len__(self):
        return self.data.size


class DataReader():
    """
    This class is responsible for handling input files. Given input
    files it reads them and returns NetCDF Dataset objects to python. It is
    also responsible for performing the time-averaging calculation at load time.

    :author: Nicolas Strike
    :data: January 2019
    """

    def __init__(self):
        """
        This is the object constructor. Initializes class variables and data
        :author: Nicolas Strike
        """
        self.nc_filenames = {}
        self.nc_datasets = {}
        self.root_dir = pathlib.Path(__file__).parent
        self.panels_dir = self.root_dir.as_uri() + "/cases/panels/"

    def cleanup(self):
        """
        This is the cleanup method. This is called on the instance's destruction
        to deallocate resources that may be held (e.g. dataset files).
        :return:
        :author: Nicolas Strike
        """
        for dataset in self.nc_datasets:
            dataset.close()

    def loadFolder(self, folder_path, ignore_git=True):
        """
        Finds all dataset files in a given folder and loads
        them using the appropriate helper class.
        :param folder_path: The path of the folder to be loaded
        :param ignore_git: Ignore files and paths that contain '.git' in their name
        :return: A 2 dimentional dictionary where a case_key (e.g. gabls3_rad)
        contains a dictionary defining filenames behind a filetype key (e.g. zm)
        For example: to access the zm data for gabls3_rad, simply call nc_datasets['gabls3']['zm']
        :author: Nicolas Strike
        """
        for root, dirs, files in os.walk(folder_path):
            for filename in files:
                abs_filename = os.path.abspath(os.path.join(root, filename))
                file_ext = os.path.splitext(filename)[1]
                if ignore_git and '.git' in abs_filename or file_ext != '.nc':
                    continue
                ext_offset = filename.rindex(
                    '_')  # Find offset to eliminate trailing chars like "_zt.nc", "_zm.nc", and "_sfc.nc
                file_type = filename[ext_offset + 1:-3]  # Type of current file (zm, zt, or sfc)
                case_key = filename[:ext_offset]
                if case_key in self.nc_filenames.keys():
                    self.nc_filenames[case_key][file_type] = abs_filename
                    self.nc_datasets[case_key][file_type] = self.__loadNcFile__(abs_filename)
                else:
                    self.nc_filenames[case_key] = {file_type: abs_filename}
                    self.nc_datasets[case_key] = {file_type: self.__loadNcFile__(abs_filename)}
        return self.nc_datasets

    def getVarData(self, netcdf_dataset, ncdf_variable, fill_zeros=False):
        """
        Given a Dataset and NetCdfVariable object, this function returns the numerical
        data for the given variable.

        :param netcdf_dataset: Dataset containing the given variable
        :param ncdf_variable: NetCdfVariable object to get data values for
        :param fill_zeros: If the variable is not found in the Dataset, setting fill_zeros = True will allow
            pyplotgen to autofill 0 for each datapoint.
        :return: The list of numeric values for the given variable
        """
        start_time_value = ncdf_variable.start_time
        end_time_value = ncdf_variable.end_time
        variable_name = ncdf_variable.name
        conv_factor = ncdf_variable.conv_factor
        avg_axis = ncdf_variable.avg_axis
        time_conv_factor = 1

        time_values = self.__getValuesFromNc__(netcdf_dataset, "time", time_conv_factor)
        if end_time_value == -1:
            if variable_name != 'time' and variable_name != 'z' and variable_name != 'altitude' and variable_name != 'lev':
                warn("End time value was not specified (or was set to -1) for variable "+variable_name+". Automatically using last time in dataset.")
            end_time_value = time_values[-1]
        # Get the index values that correspond to the desired start/end x values
        start_avg_index, end_avg_idx = self.__getStartEndIndex__(time_values, start_time_value, end_time_value)

        try:
            values = self.__getValuesFromNc__(netcdf_dataset, variable_name, conv_factor)
        except ValueError:
            # Auto fill values with zeros if allowed, otherwise propagate the error
            if fill_zeros:
                size = 0
                for dimention in netcdf_dataset.dimensions.values():
                    if dimention.size > size:
                        size = dimention.size
                values = np.zeros(size)
            else:
                raise
        if values.ndim > 1:  # not ncdf_variable.one_dimensional:
            values = self.__meanProfiles__(values, start_avg_index, end_avg_idx + 1, avg_axis=avg_axis)

        return values

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Calls the cleanup and cleanly closes out the object isntance
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        :author: Nicolas Strike
        """
        self.cleanup()

    def __loadNcFile__(self, filename):
        """
        Load the given NetCDF file
        :param filename: the netcdf file to be loaded
        :return: a netcdf dataset containing the data from the given file
        """
        dataset = Dataset(filename, "r+", format="NETCDF4")
        return dataset

    def __meanProfiles__(self, var, idx_t0=0, idx_t1=-1, idx_z0=0, idx_z1=-1, avg_axis=0):
        """
        Input:
          var    -- time x height array of some property
          idx_t0 -- Index corresponding to the beginning of the averaging interval
          idx_t1 -- Index corresponding to the end of the averaging interval, this value is exclusive.
                    If you want the data at this index to be included in averaging do your index + 1
          idx_z0 -- Index corresponding to the lowest model height of the averaging interval
          idx_z1 -- Index corresponding to the highest model height of the averaging interval

        Output:
          var    -- time averaged vertical profile of the specified variable
        """
        idx_t1 = idx_t1
        if idx_t1 is -1:
            idx_t1 = len(var)
            warn("An end index for the time averaging interval was not specified. Automatically using the last index.")
        if idx_t1 - idx_t0 <= 1:
            warn("Time averaging interval is less than or equal to 1 (idx_t0 = " + str(idx_t0) + ", idx_t1 = " + str(
                idx_t1) + ").")
        if avg_axis is 0:  # if time-averaged
            var_average = np.nanmean(var[idx_t0:idx_t1, :], axis=avg_axis)
        else:  # if height averaged
            var_average = np.nanmean(var[:, idx_z0:idx_z1], axis=avg_axis)
            warn("Using height averaging. If this is not desirable, change your averaging axis from 1 to 0 for time "
                 "averaging instead.")
        return var_average

    def __getUnits__(self, datasets, varname):
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          unit as string
        """
        if not isinstance(datasets, dict):
            datasets = {'auto dict': datasets}
        for dataset in datasets.values():
            keys = dataset.variables.keys()
            if varname in keys:
                unit = dataset.variables[varname].units
                break
            else:
                unit = "n/a"
        return unit

    def getLongName(self, ncdf_datasets, varname):
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          long_name as string, often used as a title
        """

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'auto': ncdf_datasets}

        long_name = "longname not found"
        for dataset in ncdf_datasets.values():
            keys = dataset.variables.keys()
            if varname in keys:
                long_name = dataset.variables[varname].long_name
                break
        return long_name

    def getAxisTitle(self, ncdf_datasets, varname):
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          long_name as string, often used as a title
        """
        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'auto': ncdf_datasets}

        axis_title = "axis title not found"
        for dataset in ncdf_datasets.values():
            keys = dataset.variables.keys()
            if varname in keys:
                imported_name = dataset.variables[varname].name
                units = dataset.variables[varname].units
                axis_title = imported_name + ' ' + '[' + units + ']'
                break
        return axis_title

    def __getValuesFromNc__(self, ncdf_data, varname, conversion):
        """
        Get data values out of a netcdf object, returning them as an array

        :param ncdf_data: Netcdf file object
        :param varname: Variable name string
        :param conversion: Conversion factor

        :return: time x height array of the specified variable, scaled by conversion factor
        """
        if ncdf_data is None:
            raise ValueError("ncdf_data was passed as None into __getValuesFromNc__ while looking for variable " + varname)
        var_values = None
        keys = ncdf_data.variables.keys()
        if varname in keys:
            var_values = ncdf_data.variables[varname]
            var_values = np.squeeze(var_values)
            var_values = var_values * conversion
            # Variables with 0-1 data points/values return a float after being 'squeeze()'ed, this converts it back to an array
            if isinstance(var_values, float):
                var_values = np.array([var_values])
            # break
        if var_values is None:
            raise ValueError("Variable " + varname + " does not exist in ncdf_data file. If this is expected,"
                                                     " try passing fill_zeros=True when you create the "
                                                     "NetCdfVariable for " + varname + ".\nVariables found in dataset: " + str(ncdf_data))
        if varname == 'time' and 'hour' in ncdf_data.variables[varname].units:
            var_values = var_values[:] * 60
            var_values = var_values[:] - var_values[0] + 1

        if varname == 'time' and 'sec' in ncdf_data.variables[varname].units:
            var_values = var_values[:] / 60
            var_values = var_values[:] - var_values[0] + 1

        elif self.getNcdfSourceModel(ncdf_data) == 'unknown-model' and 'sfc' not in ncdf_data.history:
            warn("Warning, unknown model detected. PyPlotgen doesn't know where this netcdf data is from." + str(ncdf_data))
            if varname == 'time':
                warn("Attempting to autoshift time values")
                var_values = var_values[:] - var_values[0] + 1

        if varname == 'time' and var_values[0] != 1:
            warn("First time value is " + str(var_values[0]) + " instead of 1. Are these time values supposed to be scaled to minutes?")

        return var_values

    def __getStartEndIndex__(self, data, start_value, end_value):
        """
        Get the list floor index that contains the value to start graphing at and the
        ceiling index that contains the end value to stop graphing at

        If neither are found, returns the entire array back
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx   which contains the starting and ending index representing the start and end time passed into the function
        :author: Nicolas Strike
        """
        start_idx = 0
        end_idx = len(data) - 1
        for i in range(0, len(data)):
            # Check for start index
            test_value = data[i]
            if test_value <= start_value and test_value > data[start_idx]:
                start_idx = i
            # Check for end index
            if test_value >= end_value and test_value < data[end_idx]:
                end_idx = i
        if end_idx == -1:
            raise ValueError("Event end_value " + str(
                end_value) + " is too large. There is no subset of data [start:end] that contains " + str(end_value))

        return start_idx, end_idx

    def getNcdfSourceModel(self, ncdf_dataset):
        """
        Guesses the model that outputted a given ncdf file. Currently does this by investigating which type of
        elevation is outputted (e.g. altitude, z, lev). If there is not elevation parameter found, and the filename
        contains 'sfc' it will return a generic 'sfc calculations' value.

        :param ncdf_dataset: NetCDF Dataset object imported from output file from some supported model (e.g. CLUBB, SAM)
        :return: the (estimated) name of the model that outputted the input file
        """
        if isinstance(ncdf_dataset, Dataset):
            ncdf_dataset = {'temp': ncdf_dataset}
        for dataset in ncdf_dataset.values():
            if 'altitude' in dataset.variables.keys():
                return 'clubb'
            elif 'z' in dataset.variables.keys():
                return 'sam'
            elif 'lev' in dataset.variables.keys():
                return 'coamps'
            elif 'sfc' in dataset.history:
                return 'sfc calculations'
            else:
                return 'unknown-model'
