"""
TODO: Change independent <-> dependent, since independent data are the dimension variable
        and dependent are all the other variables

:author: Nicolas Strike
:date: Early 2019
"""
import os
import pathlib as pathlib
from collections.abc import Iterable
from warnings import warn

import numpy as np
from netCDF4 import Dataset

from config import Case_definitions


class NetCdfVariable:
    """
    Class used for conveniently storing the information about a given netcdf variable
    """

    def __init__(self, names, ncdf_data, independent_var_names=None, conversion_factor=1, start_time=0, end_time=-1,
                 min_height=0, max_height=-1, avg_axis=0):
        """
        Initialize a new instance of this class and return it.

        TODO: Test with multidimensional variables!

        :param names: The name of the variable as defined in the ncdf_data
        :param ncdf_data: Accepts either a Dataset object to pull dependent_data from,
            or a dict of Datasets which will automatically be searched for the variable
        :param independent_var_names: List of or single name of independent variable to be used for this object
        :param conversion_factor: A multiplication factor applied to every value in the dataset. Defaults to 1
        :param start_time: The time value to begin the averaging period, e.g. 181 minutes. Defaults to 0.
        :param end_time: The time value to stop the averaging period, e.g. 240 minutes. Defaults to -1.
        :param avg_axis: The axis to average dependent_data over. 0 for time-avg, 1 for height avg, 2 for no averaging
        """
        dataset_with_var = None
        varname = ""
        var_found_in_dataset = False
        independent_var_name_in_dataset = {}
        independent_keys = ['time', 'height']
        # If ncdf_data is an netCDF4.Dataset, wrap it into a dict
        if isinstance(ncdf_data, Dataset):
            ncdf_data = {'temp': ncdf_data}
        # If names and independent_var_name are not lists, wrap them into lists
        if not isinstance(names, list):
            names = [names]

        # TODO: Figure out the best way to implement independent_var_names for time-height plots
        if isinstance(independent_var_names, dict):
            # if avg_axis != 2:
            #     raise ValueError('Error in parameter independent_var_names: Dict type only valid if avg_axis==2.')
            if any([key not in independent_var_names.keys() for key in independent_keys]):
                raise KeyError('Error in parameter independent_var_names: Dict keys must include "time" and "height".')
            # time_vars = independent_var_names['time']
            # height_vars = independent_var_names['height']
        elif isinstance(independent_var_names, Iterable):
                # If independent_var_names is iterable, split it up into time and height varnames
                time_vars = list(set(independent_var_names).intersection(set(Case_definitions.TIME_VAR_NAMES)))
                height_vars = list(set(independent_var_names).intersection(set(Case_definitions.HEIGHT_VAR_NAMES)))
                independent_var_names = {'time': time_vars, 'height': height_vars}
        else:
            # Asssume passed value is a string and either a valid time or height varname
            if independent_var_names in Case_definitions.TIME_VAR_NAMES:
                time_vars = [independent_var_names]
                height_vars = []
            else:
                time_vars = []
                height_vars = [independent_var_names]
            independent_var_names = {'time': time_vars, 'height': height_vars}

        # Find a variable with a name listed in <<names>> in any of the datasets
        # and save both its variable name and the dataset
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
            if len(ncdf_data.values()) == 0:
                # TODO patch this bug
                warn("Some model is missing files for " + str(names) +". Try either including these files or "
                     "removing their names from the Case_definitions.py config file. "
                     "This warning is a temporary notice until a bug related to missing filenames is patched")
            # changed to warn from 'raise NameError'
            # If variable was not found, select a valid dataset by looking for a matching independent variable name
            # Extension for multiple dimensions:
            # In each dataset, find an independent varname for each dimension
            # Check for success depending on averaging mode
            for dataset in ncdf_data.values():
                keys = dataset.variables.keys()
                for independent_key in independent_var_names:
                    varnames = independent_var_names[independent_key]
                    for temp_independent_var in varnames:
                        if temp_independent_var in keys:
                            if temp_independent_var == "lev" and "Z3" in keys:
                                continue
                            dataset_with_var = dataset
                            independent_var_name_in_dataset[independent_key] = temp_independent_var
                            break
                # Check if for each independent dimension a variable was found
                key_test = [key not in independent_var_name_in_dataset.keys() for key in independent_keys]
                if avg_axis == 2: # non-averaged case (time-height plots)
                    # For both time and height an independent variable must be found, otherwise try the next dataset
                    if any(key_test):
                        dataset_with_var = None
                        independent_var_name_in_dataset = {}
                    else:
                        break
                else: # averaged case
                    # TODO: differentiate between time and height averaged variables!
                    # (Depending on avg_axis, only one of either time or height independent variables must be found)
                    if all(key_test):
                        dataset_with_var = None
                        independent_var_name_in_dataset = {}
                    else:
                        break
            if dataset_with_var is None:
                raise KeyError("Could not find both ", names, " and ", independent_var_names, " in it in datasets ",
                               ncdf_data)
            # dataset_with_var = next(iter(ncdf_data.values()))
            warn("None of the values " + str(names) + " were found in the dataset " + str(dataset_with_var.filepath()))
            varname = names[0]
        # If not already found, find a matching independent variable in the same dataset
        # TODO: Accomodate finding time AND height variables
        if not independent_var_name_in_dataset: # Check if independent_var_name_in_dataset is still empty
            # If empty, try and find independent variable(s) in dataset_with_var
            for independent_key in independent_var_names:
                varnames = independent_var_names[independent_key]
                for tempname in varnames:
                    if tempname in dataset_with_var.variables.keys():
                        # Attempt to determine whether or not the height var is actually a height var (e.g. cam "lev"
                        # var is not)
                        if independent_key == "height" and hasattr(dataset_with_var.variables[tempname], 'long_name') :
                            # and if that title contains "height", then it's a height var
                            if "height" in dataset_with_var.variables[tempname].long_name.lower():
                                independent_var_name_in_dataset[independent_key] = tempname
                            # skip varname matches that aren't height vars if we know they're not height vars
                            else:
                                continue
                        # assume it actually is a height var if we can't prove it's not or a time var
                        else:
                            independent_var_name_in_dataset[independent_key] = tempname
                        break
            # Check if independent variables were found
            key_test = [key not in independent_var_name_in_dataset.keys() for key in independent_keys]
            if avg_axis == 2: # non-averaged case (time-height plots)
                # For both time and height an independent variable must be found, otherwise try the next dataset
                if any(key_test):
                    independent_var_name_in_dataset = None
            else: # averaged case
                # Depending on avg_axis, only one of either time or height independent variables must be found
                if avg_axis == 0 and 'height' not in independent_var_name_in_dataset.keys() \
                or avg_axis == 1 and 'time' not in independent_var_name_in_dataset.keys():
                    independent_var_name_in_dataset = None
        # Store the dataset in which the variable was found
        self.ncdf_data = dataset_with_var
        self.varname = varname
        self.independent_var_name = independent_var_name_in_dataset
        self.start_time = start_time
        self.end_time = end_time
        self.min_height = min_height
        self.max_height = max_height
        self.conversion_factor = conversion_factor
        self.avg_axis = avg_axis
        data_reader = DataReader()
        # Get dependent and independent data from the chosen dataset
        self.dependent_data, self.independent_data = data_reader.getVarData(self.ncdf_data, self)

    def trimArray(self, start_value, end_value, data=None, axis=0):
        """
        Remove all dependent_data elements from the variable that are not between the start and end value. Assumes
        the dependent_data is sorted. If the optional data parameter is used, it will restrict to
        the indices of where the start/end values are in that dataset rather than the variables dependent_data itself.

        TODO: Test with multiple dimensions!

        :param start_value: The smallest possible value
        :param end_value: The largest possible value
        :param data: A list of dependent_data values to find start/ending indices for constraining rather than using
            the NetCdfVariable's dependent_data. Must be specified if called on a multidimensional variable.
        :param axis: The index for the axis to be trimmed (0=time, 1=height). Only needed for multidimensional data.
        :return: None, operates in place
        """
        if data is None:
            data = self.dependent_data
            if len(data.shape)>1:
                warn('Warning! trimArray was called on a multidimensional array without specifying data.')
                return None

        start_idx, end_idx = DataReader.__getStartEndIndex__(data, start_value, end_value)
        if len(self.dependent_data.shape)>1: # multiple dimensions
            if axis==0:
                # Trim 1st coordinate -> time dimension
                self.dependent_data = self.dependent_data[start_idx:end_idx]
                if self.independent_data is not None:
                    self.independent_data['time'] = self.independent_data['time'][start_idx:end_idx]
            else:
                # Trim 2nd coordinate -> height dimension
                self.dependent_data = self.dependent_data[:,start_idx:end_idx]
                if self.independent_data is not None:
                    self.independent_data['height'] = self.independent_data['height'][start_idx:end_idx]
        else:
            self.dependent_data = self.dependent_data[start_idx:end_idx]
            if self.independent_data is not None:
                self.independent_data = self.independent_data[start_idx:end_idx]

    def filterDatasets(self):
        """
        Looks through the nc files for a case and finds the Dataset
        containing the requested variable and assigns it to self.ncdf_data.
        Assumes self.ncdf_data is a dictionary.

        :return: None
        """
        for subdataset in self.ncdf_data.values():
            if self.varname in subdataset.variables.keys():
                self.ncdf_data = subdataset
                break

    def __len__(self):
        return self.dependent_data.size


class DataReader():
    """
    This class is responsible for handling nc files.
    Given nc files it reads them and returns NetCDF Dataset objects to python.
    It is also responsible for performing the time-averaging calculation at load time.

    :author: Nicolas Strike
    :date: January 2019
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

        :return: None
        :author: Nicolas Strike
        """
        # self.nc_datasets is a 3-level dictionary with the following key levels:
        # 1. case
        # 2. folder
        # 3. filetype
        # So in order to close the Dataset objects we need to cycle through all three levels
        for case in self.nc_datasets:
            casedict = self.nc_datasets[case]
            for folder in casedict:
                folderdict = casedict[folder]
                for filetype in folderdict:
                    dataset = folderdict[filetype]
                    dataset.close()

    def loadFolder(self, folder_path, ignore_git=True):
        """
        Finds all dataset files in a given folder and loads them using the appropriate helper class.

        :param folder_path: A list of paths to folders containing nc files
        :param ignore_git: Ignore files and paths that contain '.git' in their name
        :return: A 2 dimensional dictionary where a case_key (e.g. gabls3_rad)
            contains a dictionary defining filenames behind a filetype key (e.g. zm)
            For example: to access the zm dependent_data for gabls3_rad, simply call clubb_datasets['gabls3']['zm']
        :author: Nicolas Strike
        """
        for sub_folder in folder_path:
            # Recursively generate tuples for each folder that is found in the directory tree below sub_folder
            # containining the name of the folder, its subdirectories, and its files
            for root, dirs, files in os.walk(sub_folder):
                # Check each file in root folder
                for filename in files:
                    abs_filename = os.path.abspath(os.path.join(root, filename))
                    file_ext = os.path.splitext(filename)[1]
                    # Only process nc files (not containing '.git')
                    if ignore_git and '.git' in abs_filename or file_ext != '.nc':
                        continue
                    # Find offset to eliminate trailing chars like "_zt.nc", "_zm.nc", and "_sfc.nc
                    ext_offset = filename.rindex('_')
                    # Get type of current file (zm, zt, or sfc)
                    file_type = filename[ext_offset + 1:-3]
                    # Remaining prefix of filename is the case's name to which the file data belongs
                    case_key = filename[:ext_offset]

                    if case_key in self.nc_filenames.keys() and sub_folder in self.nc_filenames[case_key].keys():
                        # Another file from the same folder for the same case was loaded before
                        self.nc_filenames[case_key][sub_folder][file_type] = abs_filename
                        self.nc_datasets[case_key][sub_folder][file_type] = self.__loadNcFile__(abs_filename)
                    elif case_key in self.nc_filenames.keys() and sub_folder not in self.nc_filenames[case_key].keys():
                        # Another file for the same case but from a different folder was loaded before
                        self.nc_filenames[case_key][sub_folder] = {file_type: abs_filename}
                        self.nc_datasets[case_key][sub_folder] = {file_type: self.__loadNcFile__(abs_filename)}
                    else:
                        # This is the first file for this specific case
                        self.nc_filenames[case_key] = {sub_folder: {file_type: abs_filename}}
                        self.nc_datasets[case_key] = {sub_folder: {file_type: self.__loadNcFile__(abs_filename)}}
        return self.nc_datasets

    def getVarData(self, netcdf_dataset, ncdf_variable):
        """
        Given a Dataset and NetCdfVariable object, this function returns the numerical
        dependent_data and independent_data for the given variable.

        TODO: Split up reading and averaging into different functions
            since time-height plots do not need averaging

        :param netcdf_dataset: Dataset containing the given variable
        :param ncdf_variable: NetCdfVariable object to get dependent_data values for
        :return: The list of numeric values for the given variable
        """
        start_time_value = ncdf_variable.start_time
        end_time_value = ncdf_variable.end_time
        # min_height = ncdf_variable.min_height
        # max_height = ncdf_variable.max_height
        variable_name = ncdf_variable.varname
        conv_factor = ncdf_variable.conversion_factor
        avg_axis = ncdf_variable.avg_axis
        independent_var_name = ncdf_variable.independent_var_name
        time_conv_factor = 1
        # height_conv_factor = 1
        time_values = None
        # height_values = None
        # Get time dimension from netcdf_dataset
        for time_var in Case_definitions.TIME_VAR_NAMES:
            if time_var in netcdf_dataset.variables.keys():
                time_values = self.__getValuesFromNc__(netcdf_dataset, time_var, time_conv_factor)
                # np.savetxt("time.csv", time_values, delimiter=',', fmt='%f') # occasionally used when debugging

        if time_values is None:
            raise NameError("None of the time variables " + str(Case_definitions.TIME_VAR_NAMES) +
                            " could be found in the dataset " + str(netcdf_dataset.filepath()))

        # If end_time_value was set to its default value, set to last time value
        if end_time_value == -1:
            if variable_name not in Case_definitions.HEIGHT_VAR_NAMES \
                    and variable_name not in Case_definitions.TIME_VAR_NAMES:
                warn("End time value was not specified (or was set to -1) for variable " + variable_name +
                     ". Automatically using last time in dataset.")
            end_time_value = time_values[-1]

        # Get the index values that correspond to the desired start/end x values
        start_avg_idx, end_avg_idx = 0,-1
        if ncdf_variable.avg_axis==0:
            start_avg_idx, end_avg_idx = self.__getStartEndIndex__(time_values, start_time_value, end_time_value)
        elif ncdf_variable.avg_axis==1:
            warn('Warning! Averaging over heights is not yet implemented.')
            start_avg_idx, end_avg_idx = self.__getStartEndIndex__(time_values, start_time_value, end_time_value)
        # Get independent values from nc file
        if independent_var_name is None:
            raise NameError("None of the independent variables " + str(independent_var_name) +
                            " could be found in the dataset " + str(ncdf_variable.ncdf_data.filepath()))
        else:
            independent_values = dict()
            if ncdf_variable.avg_axis == 0:
                independent_values = self.__getValuesFromNc__(netcdf_dataset, independent_var_name['height'],1)
                if independent_var_name['height'] == 'Z3':
                    independent_values = self.__averageData__(independent_values,idx_t0=start_avg_idx, idx_t1=end_avg_idx, avg_axis=0)


            elif ncdf_variable.avg_axis == 1:
                independent_values = self.__getValuesFromNc__(netcdf_dataset, independent_var_name['time'], 1)
            elif ncdf_variable.avg_axis == 2:
                independent_values['height'] = self.__getValuesFromNc__(netcdf_dataset, independent_var_name['height'],
                                                                        1)
                independent_values['time'] = self.__getValuesFromNc__(netcdf_dataset, independent_var_name['time'], 1)
            # np.savetxt("" + independent_var_name + ".csv", independent_values,  delimiter=',', fmt='%f')  # occasionally used when debugging



        # Try and get dependent_data from nc file
        try:
            dependent_values = self.__getValuesFromNc__(netcdf_dataset, variable_name, conv_factor)
            # np.savetxt("" + variable_name + ".csv", dependent_values,  delimiter=',', fmt='%f')  # occasionally used when debugging

        except ValueError:
            dependent_values = np.zeros(len(independent_values))

        # If dependent_data is more than 1-dimensional, average over avg_axis if needed
        if dependent_values.ndim > 1 and avg_axis in [0,1]:  # not ncdf_variable.one_dimensional:
            dependent_values = self.__averageData__(dependent_values, start_avg_idx, end_avg_idx, avg_axis=avg_axis)
        # SAM data may contain NaNs at this point. Change those to 0
        if 'SAM version' in netcdf_dataset.ncattrs():
            dependent_values = np.where(np.isnan(dependent_values), 0, dependent_values)

        if dependent_values.ndim > 1:
            raise ValueError("Number of dependent dimentions is greater than 1")
        if independent_values.ndim > 1:
            raise ValueError("Number of independent dimentions is greater than 1")
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

    def __del__(self):
        """
        Destructor for DataReader class.
        Calls the cleanup and cleanly closes out the object instance

        :author: Nicolas Strike
        """
        self.cleanup()

    def __loadNcFile__(self, filename):
        """
        Load the given NetCDF file

        :param filename: The netcdf file to be loaded
        :return: A netCDF4 Dataset object containing the data from the given file
        """
        dataset = Dataset(filename, "r", format="NETCDF4")
        return dataset

    def __averageData__(self, var, idx_t0=0, idx_t1=-1, idx_z0=0, idx_z1=-1, avg_axis=0):
        """
        Averages 2d data with one time and one height dimension and returns the averaged data

        TODO: Change avg_axis to receive string instead of int -> 'time'/'height'/...?

        :param var: 2d data array with one time and one height dimension
        :param idx_t0: Index corresponding to the beginning of the averaging interval
        :param idx_t1: Index corresponding to the end of the averaging interval, this value is exclusive.
            If you want the dependent_data at this index to be included in averaging do your index + 1
        :param idx_z0: Index corresponding to the lowest model height of the averaging interval
        :param idx_z1: Index corresponding to the highest model height of the averaging interval
        :param avg_axis: If 0, the data will be averaged over axis 0, the time axis;
            If 1, the data will be averaged over axis 1, the height axis
        :return: 1d data, averaged over the interval [idx_t0:idx_t1] or [idx_z0:idx_z1],
            depending on the averaging axis
        """
        if avg_axis not in [0,1]:
            raise ValueError('An invalid value for avg_axis was specified. '+
                             'Only 0 (time) and 1 (height) are valid values.')
        if idx_t1 - idx_t0 <= 10:
            warn("Time averaging interval is small (less than or equal to 10): " + str(idx_t1 - idx_t0) +
                 " | (idx_t0 = " + str(idx_t0) + ", idx_t1 = " + str(
                idx_t1) + "). Note, start index is inclusive, end index is exclusive.")
        if avg_axis == 0:  # if time-averaged
            if idx_t1 == -1:
                idx_t1 = len(var)
                warn("An end index for the time averaging interval was not specified. Automatically using the last index.")
            var_average = np.nanmean(var[idx_t0:idx_t1, :], axis=avg_axis)
        elif avg_axis == 1:  # if height averaged
            if idx_z1 == -1:
                idx_z1 = len(var)
                warn("An end index for the height averaging interval was not specified. Automatically using the last index.")
            var_average = np.nanmean(var[:, idx_z0:idx_z1], axis=avg_axis)
            warn("Using height averaging. If this is not desirable, change your averaging axis from 1 to 0 for time "
                 "averaging instead.")
        return var_average

    def __getUnits__(self, datasets, varnames):
        """
        Find netcdf variable in datasets matching any name in varnames and get its units attribute

        :param datasets: List of datasets
        :param varnames: List of variable names
        :return: String containing the physical units of the variable that was found
        """
        units = "Units n/a"
        if not isinstance(varnames, list):
            varnames = [varnames]

        # Find dataset and variable with matching name and 'units' attribute,
        # Get the unit string, prepare it for latex text rendering and return it
        for name in varnames:
            var_found = False
            for dataset in datasets:
                keys = dataset.variables.keys()
                # If variable exists and has a "units" attribute
                if name in keys and hasattr(dataset.variables[name], 'units'):
                    var_found = True
                    raw_units = dataset.variables[name].units
                    raw_units = self.__remove_invalid_unit_chars__(raw_units)
                    # Encase unit string in '$' characters, signaling matplotlib
                    # that this string should be rendered with latex
                    if len(raw_units) > 0:
                        units = "$" + raw_units + "$"  # $'s are used to format equations
                        break
            if var_found:
                break
        if units == "Units n/a":
            warn("Failed to find units for variables " + str(varnames))
        return units

    def getLongName(self, ncdf_datasets, varnames):
        """
        Find netcdf variable in ncdf_datasets matching any name in varnames and get its long_name attribute

        :param ncdf_datasets: List of datasets
        :param varnames: List of variable names
        :return: String containing the long name of the variable that was found, often used as a title
        """
        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = [ncdf_datasets]

        long_name = "Title not found"
        # Find dataset and with matching variable name,
        # Get the string stored as long_name attribute and return it
        for varname in varnames:
            for dataset in ncdf_datasets:
                keys = dataset.variables.keys()
                if varname in keys:
                    long_name = dataset.variables[varname].long_name
                    return long_name
        return long_name

    def getAxisTitle(self, ncdf_datasets, varname):
        """
        Find netcdf variable in ncdf_datasets matching varname and generate the axis title from its attributes

        :param ncdf_datasets: List of datasets
        :param varname: Variable name
        :return: String containing units and variable name to be used as axis title
        """
        if isinstance(ncdf_datasets, Dataset):
            ncdf_datasets = {'auto': ncdf_datasets}

        axis_title = "Axis title not found"
        # Find dataset and with matching variable name,
        # Get name and unit attributes, generate an axis title string from those and return it
        for dataset in ncdf_datasets:
            keys = dataset.variables.keys()
            if varname in keys:
                imported_name = dataset.variables[varname].name
                units = self.__getUnits__([dataset], varname)  # dataset.variables[varname].units
                axis_title = imported_name
                axis_title += ' ' + '[' + units + ']'  # $'s are used to format equations
                break
        return axis_title

    def __getValuesFromNc__(self, ncdf_data, varname, conversion):
        """
        Get dependent_data values out of a netcdf object, returning them as an array

        :param ncdf_data: Netcdf file object
        :param varname: Variable name string
        :param conversion: Conversion factor
        :return: Data array of the specified variable, scaled by conversion factor
        """
        if ncdf_data is None:
            raise ValueError("ncdf_data was passed as None into __getValuesFromNc__ while looking for variable " +
                             varname)

        # TODO this model detection method is old and can no longer be trusted
        src_model = self.getNcdfSourceModel(ncdf_data)
        if src_model == 'unknown-model':
            warn("Warning, unknown model detected. PyPlotgen doesn't know where this netcdf dependent_data is from. "
                 + str(ncdf_data))

        var_values = None
        keys = ncdf_data.variables.keys()
        if varname in keys:
            var_values = ncdf_data.variables[varname]
            var_values = np.squeeze(var_values)
            # Check if data comes from SAM and convert -9999 values to NaN
            if 'SAM version' in ncdf_data.ncattrs():
                var_values = np.where(np.isclose(var_values, -9999), np.nan, var_values)
            var_values = var_values * conversion
            # Variables with 0-1 dependent_data points/values return a float after being 'squeeze()'ed,
            # this converts it back to an array
            if isinstance(var_values, float):
                var_values = np.array([var_values])
            # break
        if var_values is None:
            raise ValueError("Variable " + str(varname) +
                             " does not exist in ncdf_data file.\nVariables found in dataset: " + str(ncdf_data))

        # If variable contains time data, convert to minutes
        if varname in Case_definitions.TIME_VAR_NAMES:
            hrs_in_day = 24
            min_in_hr = 60
            sec_in_min = 60

            if 'day' in ncdf_data.variables[varname].units:
                var_values = var_values[:] * hrs_in_day * min_in_hr

            if 'hour' in ncdf_data.variables[varname].units:
                var_values = var_values[:] * min_in_hr

            if 'sec' in ncdf_data.variables[varname].units:
                var_values = var_values[:] / sec_in_min

            delta_t = 1
            if len(var_values) > 1:
                delta_t = var_values[1] - var_values[0]

            # In a lot of cases this loop has no effect, but for some cases (e.g. r408 lines on atex case)
            # it corrects time data.
            for i in range(len(var_values)):
                var_values[i] = delta_t * (i + 1)

            if var_values[0] > 1:
                warn("First time value is " + str(var_values[0]) +
                     " instead of 0-1. Are these time values supposed to be scaled to minutes?")

        return var_values

    @staticmethod
    def __getStartEndIndex__(data, start_value, end_value):
        """
        Get indices for values from data array corresponding to start_value and end_value.
        This function is intended to be used to get the indices for array slicing.
        E.g. given a bottom and top altitude, return the indices of the values in the array corresponding to those altitudes.
        This function can be used for height and time data.
        The data array MUST be presorted and start_value <= end_value for pyplotgen to work correctly!
        If neither are found, returns 0 and array size - 1.

        :param data: Array of numerical values. Must be presorted!
        :param start_value: The first value to be graphed (may return indexes to values smaller than this)
        :param end_value: The last value that needs to be graphed (may return indexes to values larger than this)
        :return: (tuple) start_idx, end_idx which contains the starting and ending index representing the start and
            end time passed into the function
        :author: Nicolas Strike
        """
        # If dependent_data is a an array with 1 element, return 0's
        if (len(data) == 1):
            return 0, 1

        start_idx = 0
        end_idx = len(data)
        ascending_data = data[0] < data[1]

        start_idx_found = False
        if ascending_data:
            # dependent_data is ascending
            for i in range(0, len(data)):
                # Check for start index
                test_value = data[i]
                if test_value >= start_value and not start_idx_found:
                    start_idx_found = True
                    start_idx = i
                if test_value == end_value:
                    end_idx = i + 1
                # Check for end index
                # end_idx -1 is in place because this check is inclusive, but the index is compatible with exclusivity
                if end_value < test_value < data[end_idx - 1]:
                    end_idx = i
        else:
            # Check for start index
            start_idx_found = False
            for i in range(0, len(data)):
                test_value = data[i]
                if test_value <= end_value and not start_idx_found:
                    start_idx_found = True
                    start_idx = i
                if test_value == start_value:
                    end_idx = i + 1
                # Check for end index
                # end_idx -1 is in place because this check is inclusive, but the index is compatible with exclusivity
                if start_value > test_value > data[end_idx - 1]:
                    end_idx = i

        return start_idx, end_idx

    def getNcdfSourceModel(self, ncdf_dataset):
        """
        Guesses the model that outputted a given ncdf file. Currently does this by investigating which type of
        elevation is outputted (e.g. altitude, z, lev). If there is not elevation parameter found, and the filename
        contains 'sfc' it will return a generic 'sfc calculations' value.

        :param ncdf_dataset: A netCDF4 Dataset object containing data from some supported model (e.g. CLUBB, SAM)
        :return: The (estimated) name of the model that outputted the nc file
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
            elif hasattr(dataset, 'TITLE') and 'WRF' in dataset.TITLE:
                return 'wrf'
            elif hasattr(dataset, 'history') and 'sfc' in dataset.history:
                return 'sfc calculations'
            elif hasattr(dataset, 'history') and 'Z3' in dataset.history:
                return 'e3sm'
            else:
                return 'unknown-model'

    def __remove_invalid_unit_chars__(self, units):
        """
        Removes characters from a string that are not valid for mathtext

        :param units: String to be used as unit label in plot
        :return: Cleaned string in which any invalid character has been removed or replaced
        """
        units = units.replace('#', '')
        return units
