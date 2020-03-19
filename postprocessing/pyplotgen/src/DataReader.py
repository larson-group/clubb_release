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

    def __init__(self, names, ncdf_data, independant_var_names=None, conversion_factor=1, start_time=0, end_time=-1, avg_axis=0):
        """

        :param names: the name of the variable as defined in the ncdf_data
        :param ncdf_data: Accepts either a Dataset object to pull dependent_data from, or a dict of Datasets which will automatically be searched for the variable
        :param conversion_factor: A multiplication factor multiplied to every value in the dataset. Defaults to 1
        :param start_time: The time value to being the averaging period, e.g. 181 minutes. Defaults to 0.
        :param end_time: The time value to stop the averaging period, e.g. 240 minutes. Defaults to -1.
        :param avg_axis: The axis to avg dependent_data over. 0 for time-avg, 1 for height avg
        """
        dataset_with_var = None
        varname = ""
        var_found_in_dataset = False
        # If argument names is a list of aliases, find the correct one and assign it
        if isinstance(ncdf_data, Dataset):
            ncdf_data = {'temp': ncdf_data}
        if not isinstance(names, list):
            names = [names]
        if not isinstance(independant_var_names, list) and independant_var_names is not None:
            independant_var_names = [independant_var_names]
        for dataset in ncdf_data.values():
            for tempname in names:
                if tempname in dataset.variables.keys():
                    dataset_with_var = dataset
                    varname = tempname
                    var_found_in_dataset = True
                    break
            if var_found_in_dataset:
                break
        if not var_found_in_dataset:
            # changed to warn from 'raise NameError'
            dataset_with_var = next(iter(ncdf_data.values()))
            warn("None of the values " + str(names) + " were found in the dataset " + str(dataset_with_var.filepath()))
            varname = names[0]
        indep_var_name_in_dataset = None
        for tempname in independant_var_names:
            if tempname in dataset_with_var.variables.keys():
                indep_var_name_in_dataset = tempname
                break

        self.ncdf_data = dataset_with_var
        self.name = varname
        self.independent_var_name = indep_var_name_in_dataset
        self.start_time = start_time
        self.end_time = end_time
        self.conv_factor = conversion_factor
        self.avg_axis = avg_axis
        data_reader = DataReader()
        self.dependent_data, self.independent_data = data_reader.getVarData(self.ncdf_data, self)


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

        # If dependent_data is a an array with 1 element, return 0's
        if(len(data) is 1):
            return 0, 1

        start_idx = 0
        end_idx = len(data) - 1
        ascending_data = data[0] < data[1]
        if ascending_data:
            for i in range(0, len(data)):
                # Check for start index
                test_value = data[i]
                if test_value <= start_value and test_value > data[start_idx]:
                    start_idx = i
                # Check for end index
                if test_value >= end_value and test_value < data[end_idx]:
                    end_idx = i
        else: # if dependent_data is descending
            for i in range(0, len(data)):
                test_value = data[i]
                if test_value >= end_value and test_value < data[start_idx]:
                    start_idx = i
                if test_value <= start_value and test_value > data[end_idx]:
                    end_idx = i

        return start_idx, end_idx

    def constrain(self, start_value, end_value, data=None):
        """
        Remove all dependent_data elements from the variable that are not between the start and end value. Assumes
        the dependent_data is always increasing. If the optional dependent_data parameter is used, it will restrict to the indices
        of where the start/end values are in that dataset rather than the variables dependent_data itself.

        :param start_value: The smallest possible value
        :param end_value: The largest possible value
        :param data: A list of dependent_data values to find start/ending indicies for constraining rather than using the NetCdfVariable's dependent_data.
        :return: None, operates in place
        """
        if data is None:
            data = self.dependent_data
        start_idx, end_idx = self.__getStartEndIndex__(data, start_value, end_value)
        self.dependent_data = self.dependent_data[start_idx:end_idx]
        if self.independent_data is not None:
            self.independent_data = self.independent_data[start_idx:end_idx]


    def filter_datasets(self):
        """
        Looks through the clubb files for a case and finds the Dataset
        containing the requested variable and assigns it to self.ncdf_data
        Assumes self.ncdf_data is a dictionary.

        :return: None
        """
        for subdataset in self.ncdf_data.values():
            if self.name in subdataset.variables.keys():
                self.ncdf_data = subdataset
                break

    def __len__(self):
        return self.dependent_data.size


class DataReader():
    """
    This class is responsible for handling clubb files. Given clubb
    files it reads them and returns NetCDF Dataset objects to python. It is
    also responsible for performing the time-averaging calculation at load time.

    :author: Nicolas Strike
    :dependent_data: January 2019
    """

    def __init__(self):
        """
        This is the object constructor. Initializes class variables and dependent_data
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
        For example: to access the zm dependent_data for gabls3_rad, simply call clubb_datasets['gabls3']['zm']
        :author: Nicolas Strike
        """
        for sub_folder in folder_path:
            for root, dirs, files in os.walk(sub_folder):
                for filename in files:
                    abs_filename = os.path.abspath(os.path.join(root, filename))
                    file_ext = os.path.splitext(filename)[1]
                    if ignore_git and '.git' in abs_filename or file_ext != '.nc':
                        continue
                    ext_offset = filename.rindex('_')  # Find offset to eliminate trailing chars like "_zt.nc", "_zm.nc", and "_sfc.nc
                    file_type = filename[ext_offset + 1:-3]  # Type of current file (zm, zt, or sfc)
                    case_key = filename[:ext_offset]
                    if case_key in self.nc_filenames.keys() and sub_folder in self.nc_filenames[case_key].keys():
                        self.nc_filenames[case_key][sub_folder][file_type] = abs_filename
                        self.nc_datasets[case_key][sub_folder][file_type] = self.__loadNcFile__(abs_filename)
                    elif case_key in self.nc_filenames.keys() and sub_folder not in self.nc_filenames[case_key].keys():
                        self.nc_filenames[case_key][sub_folder] = {file_type: abs_filename}
                        self.nc_datasets[case_key][sub_folder] = {file_type: self.__loadNcFile__(abs_filename)}
                    else:
                        self.nc_filenames[case_key] = {sub_folder: {file_type: abs_filename}}
                        self.nc_datasets[case_key] = {sub_folder: {file_type: self.__loadNcFile__(abs_filename)} }
        return self.nc_datasets

    def getVarData(self, netcdf_dataset, ncdf_variable):
        """
        Given a Dataset and NetCdfVariable object, this function returns the numerical
        dependent_data for the given variable.

        :param netcdf_dataset: Dataset containing the given variable
        :param ncdf_variable: NetCdfVariable object to get dependent_data values for
        :return: The list of numeric values for the given variable
        """
        start_time_value = ncdf_variable.start_time
        end_time_value = ncdf_variable.end_time
        variable_name = ncdf_variable.name
        conv_factor = ncdf_variable.conv_factor
        avg_axis = ncdf_variable.avg_axis
        independent_var_name = ncdf_variable.independent_var_name
        time_conv_factor = 1
        time_values = self.__getValuesFromNc__(netcdf_dataset, "time", time_conv_factor)
        if end_time_value == -1:
            if variable_name != 'time' and variable_name != 'z' and variable_name != 'altitude' and variable_name != 'lev' and variable_name != 'Z3':
                warn("End time value was not specified (or was set to -1) for variable "+variable_name+". Automatically using last time in dataset.")
            end_time_value = time_values[-1]
        # Get the index values that correspond to the desired start/end x values
        start_avg_index, end_avg_idx = self.__getStartEndIndex__(time_values, start_time_value, end_time_value)
        independent_values = self.__getValuesFromNc__(netcdf_dataset, independent_var_name, 1)

        try:
            dependent_values = self.__getValuesFromNc__(netcdf_dataset, variable_name, conv_factor)
        except ValueError:
            dependent_values = np.zeros(len(independent_values))

        if dependent_values.ndim > 1:  # not ncdf_variable.one_dimensional:
            dependent_values = self.__meanProfiles__(dependent_values, start_avg_index, end_avg_idx + 1, avg_axis=avg_axis)

        # E3SM outputs Z3 as it's height variable, which may also contain an offset
        # (e.g. e3sm height = clubb height + 650 for the dycoms2_rf01 case). This eliminates that offset.
        if variable_name == 'Z3':
            dependent_values = dependent_values - dependent_values[-1]

        # Change all remaining NaN values to 0
        if 'SAM version' in netcdf_dataset.ncattrs():
            dependent_values = np.where(np.isnan(dependent_values), 0, dependent_values)

        return dependent_values, independent_values

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Calls the cleanup and cleanly closes out the object instance
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
        :return: a netcdf dataset containing the dependent_data from the given file
        """
        dataset = Dataset(filename, "r+", format="NETCDF4")
        return dataset

    def __meanProfiles__(self, var, idx_t0=0, idx_t1=-1, idx_z0=0, idx_z1=-1, avg_axis=0):
        """
        Input:
          var    -- time x height array of some property
          idx_t0 -- Index corresponding to the beginning of the averaging interval
          idx_t1 -- Index corresponding to the end of the averaging interval, this value is exclusive.
                    If you want the dependent_data at this index to be included in averaging do your index + 1
          idx_z0 -- Index corresponding to the lowest model height of the averaging interval
          idx_z1 -- Index corresponding to the highest model height of the averaging interval

        Output:
          var    -- time averaged vertical profile of the specified variable
        """
        idx_t1 = idx_t1
        if idx_t1 is -1:
            idx_t1 = len(var)
            warn("An end index for the time averaging interval was not specified. Automatically using the last index.")
        if idx_t1 - idx_t0 <= 10:
            warn("Time averaging interval is small (less than or equal to 10): " + str(idx_t1 - idx_t0) + " | (idx_t0 = " + str(idx_t0) + ", idx_t1 = " + str(
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
        units = "units n/a"
        for dataset in datasets:
            keys = dataset.variables.keys()
            if varname in keys and hasattr(dataset.variables[varname], 'units'):
                raw_units = dataset.variables[varname].units
                raw_units = self.__remove_invalid_unit_chars__(raw_units)
                if len(raw_units) > 0:
                    units = "$" + raw_units + "$" # $'s are used to format equations
                    break
        return units

    def getLongName(self, ncdf_datasets, varnames):
        """
        Input:
          nc         --  Netcdf file object
          varname    --  Variable name string
        Output:
          long_name as string, often used as a title
        """

        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'auto': ncdf_datasets}

        long_name = "title not found"
        for varname in varnames:
            for dataset in ncdf_datasets:
                keys = dataset.variables.keys()
                if varname in keys:
                    long_name = dataset.variables[varname].long_name
                    return long_name
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
        for dataset in ncdf_datasets:
            keys = dataset.variables.keys()
            if varname in keys:
                imported_name = dataset.variables[varname].name
                units = self.__getUnits__([dataset], varname)#dataset.variables[varname].units
                axis_title = imported_name
                axis_title += ' ' + '[' + units + ']' # $'s are used to format equations
                break
        return axis_title

    def __getValuesFromNc__(self, ncdf_data, varname, conversion):
        """
        Get dependent_data values out of a netcdf object, returning them as an array

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
            # Check if data comes from SAM and convert -9999 values to NaN
            if 'SAM version' in ncdf_data.ncattrs():
                var_values = np.where(np.isclose(var_values, -9999), np.nan, var_values)
            var_values = var_values * conversion
            # Variables with 0-1 dependent_data points/values return a float after being 'squeeze()'ed, this converts it back to an array
            if isinstance(var_values, float):
                var_values = np.array([var_values])
            # break
        if var_values is None:
            raise ValueError("Variable " + str(varname) + " does not exist in ncdf_data file.\nVariables found in dataset: " + str(ncdf_data))


        hrs_in_day = 24
        min_in_hr = 60
        sec_in_min = 60
        if varname == 'time' and 'day' in ncdf_data.variables[varname].units:
            var_values = var_values[:] * hrs_in_day * min_in_hr
            var_values = var_values[:] - var_values[0] + 1
            # var_values = [i for i in range(1,len(var_values) + 1)]

        if varname == 'time' and 'hour' in ncdf_data.variables[varname].units:
            var_values = var_values[:] * min_in_hr
            var_values = var_values[:] - var_values[0] + 1

        if varname == 'time' and 'sec' in ncdf_data.variables[varname].units:
            var_values = var_values[:] / sec_in_min
            var_values = var_values[:] - var_values[0] + 1

        elif self.getNcdfSourceModel(ncdf_data) == 'unknown-model' and 'sfc' not in ncdf_data.history:
            warn("Warning, unknown model detected. PyPlotgen doesn't know where this netcdf dependent_data is from." + str(ncdf_data))
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

        # If dependent_data is a an array with 1 element, return 0's
        if(len(data) is 1):
            return 0, 0

        start_idx = 0
        end_idx = len(data) - 1
        ascending_data = data[0] < data[1]
        if ascending_data:
            for i in range(0, len(data)):
                # Check for start index
                test_value = data[i]
                if test_value <= start_value and test_value > data[start_idx]:
                    start_idx = i
                # Check for end index
                if test_value >= end_value and test_value < data[end_idx]:
                    end_idx = i
        else:
            for i in range(0, len(data)):
                test_value = data[i]
                if test_value <= end_value and test_value :
                    start_idx = i
                pass
        return start_idx, end_idx

    def getNcdfSourceModel(self, ncdf_dataset):
        """
        Guesses the model that outputted a given ncdf file. Currently does this by investigating which type of
        elevation is outputted (e.g. altitude, z, lev). If there is not elevation parameter found, and the filename
        contains 'sfc' it will return a generic 'sfc calculations' value.

        :param ncdf_dataset: NetCDF Dataset object imported from output file from some supported model (e.g. CLUBB, SAM)
        :return: the (estimated) name of the model that outputted the clubb file
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
            elif 'Z3' in dataset.history:
                return 'e3sm'
            else:
                return 'unknown-model'


    def __remove_invalid_unit_chars__(self, units):
        """
        Removes characters from a string that are not valid for mathtext
        """
        units = units.replace('#', '')
        return units
