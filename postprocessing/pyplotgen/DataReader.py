'''
This file is used to load data in from netcdf and grads
formatted files and provide a Plot class instance to
a Plotter instance for graphing.

:author: Nicolas Strike
:date: Early 2019
'''

import os
import re
from collections import namedtuple
from netCDF4 import Dataset
import configparser
import pandas as pd

import numpy as np


class DataReader():
    '''
    This class is responsible for handling input files. Given a
    input file, it determines what format it is in (grads vs
    netcdf) and then uses the correct helper file (e.g.
    NetcdfReader.py) to load the data into python.

    :author: Nicolas Strike
    :data: January 2019
    '''

    panel_filename = '/home/strike/clubb/postprocessing/pyplotgen/cases/base/Panel_thlm.ini'
    # panel_type_filename = '/home/strike/clubb/postprocessing/pyplotgen/cases/PanelType_Profile.ini'


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

    def __loadNcFile__(self, filename):
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


    def loadFolder(self, folder_path):
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
                    self.nc_datasets.append(self.__loadNcFile__(abs_filename))
                elif file_ext == ".dat":
                    self.grads_dat_filenames.append(abs_filename)
                elif file_ext == ".ctl":
                    self.grads_ctl_filenames.append(abs_filename)
                else:
                    print("Filetype " + file_ext + " is not supported. Attempted to load " + abs_filename)
        return self.nc_datasets
        # print("Files loaded:\n\n----nc files----\n" + str(self.nc_filenames) + "\n\n----dat files----\n" + str(self.grads_dat_filenames) +
        #       "\n\n----ctl files----\n" + str(self.grads_ctl_filenames))


    def __meanProfiles__(self, var):
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
            var = np.zeros(shape=(level_amount, num_timesteps)) - 1.
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

    def __getArrayFromString__(self, data_string, rows = -1, columns = 1):
        '''
        Some elements of our config ini files are stored as 2d arrays
        but are read by config parser as a string. This method takes in
        the string and returns it as an array. By default it interprets
        it as a 2d array with n rows and 2 columns. In numpy, -1 denotes
        'autofill', filling that dimension until it runs out of data.

        The string given must be a series of comma seperated values.
        Do not use quote markers to denote strings, as it will blindly
        match whatever is between commas. Any newline characters in
        the string will be removed.

        :param data_string: The array string to be interpreted
        :param rows: The number of rows in the data. It is reccomended to keep this at -1 to allow numpy to autofit to the amount of data
        :param columns: The number of columns in the dataset
        :return: a numpy array containing the data from the data_string formated into the specified number of rows and columns
        '''
        data_array = re.compile("\s*,\s*").split(data_string) #np.array(re.findall("\".*?\"",data_string, re.MULTILINE)) # Turn the string into a 1d array
        # for i in data_array
        #data_array = data_array.reshape(rows,columns) # Turn the 1d array into a rows X columns array
        return data_array


    def __parseAxisData__(self, axis_data_str):
        '''
        This function takes in a string defined in a ini config
        file that defines plot variables and an equation to graph.

        Example of an input string:
            variables = THETAL, False, THETAL, 1., 0,
                THETA, False, THETA, 1., 0,
                TABS, False, TABS, 1., 0,
                QI, False, QI, 1./KG, 0,
                THETAL, True, THETAL + 2500.4 * (THETA/TABS) * QI, 1., 0
        :param axis_data_str:
        :return: Array of numeric values given from equation given in the string
        '''
        variable_lines = axis_data_str.strip(' ').split('\n', re.MULTILINE)
        lines_data = []
        for line in variable_lines:
            data = line.split(',')
            if "true" in data[1].lower():
                data[1] = True
            else:
                data[1] = False
            data[3] = float(data[3])
            data[4] = float(data[4])

            lines_data.append(data)
        return lines_data

    def isFunction(self, value):
        '''
        Given a string, this function
        determines if it represents a
        mathmatical function. In particular
        it will return true if the string
        contains any of +-*/.

        :return: True if the string is a math function, False otherwise
        '''

        # logger.info('isFunction')
        isFunc = False
        if '+' in value:
            isFunc = True
        elif '-' in value:
            isFunc = True
        elif '*' in value:
            isFunc = True
        elif '/' in value:
            isFunc = True
        return isFunc

    def __getLineToPlotIndex__(self, array2d):
        '''
        Given a 2d array containing Axis Data, aka
        variable_lines data, return the row that
        contains True in the second column. This
        is used to find the index of the line that
        contains the variable to be plotted to a graph
        (which also likely contains an equation in column
        3).

        :param array: a 2d array containing parsed axis data (see parseAxisData())
        :return: the row aka dim 0 index containing a value of True in the second column
        :author: Nicolas Strike
        '''
        rowIndex = -1
        for i in range(0, len(array2d)):
            if array2d[i][1] == True:
                rowIndex = i
                break
        return rowIndex


    def __calcVariableValues__(self, netcdf_data, plot_vars, equation, conversion_factors):
        '''
        This function takes a line config such as the one used by
        __parseAxisData__() and calculates the new data array from
        the equation listed in the variable data. If the equation
        is simply a variable name, the original data will be
        returned back.

        Example input:
            [[THETAL, False, THETAL, 1., 0],
            [THETA, False, THETA, 1., 0],
            [TABS, False, TABS, 1., 0],
            [QI, False, QI, 1., 0],
            [THETAL, True, THETAL + 2500.4 * (THETA/TABS) * QI, 1., 0]]

        :param plot_vars: A 2d array containing the variable names, conversion factors, and equation needed to calculate the new values
        :return: A 1d array containing the calculated values
        '''
        variable_values = {}
        # equation_row = self.__getLineToPlotIndex__(plot_vars)
        # equation = plot_vars[equation_row][2]
        variables_evaled = []

        #pad the equation to help with regex
        equation = " " + equation + " "
        i = 0
        for variable in plot_vars:
            #each variable should be of the form: variable name within python, is this the variable w/ equation to be plotted, variable name in SAM output, model conversion rate
            level_amount = 1 # TODO stop hardcoding this
            num_timesteps = 1 # TODO stop hardcoding this

            variable_values[variable] = self.__getValuesFromNc__(netcdf_data, variable, conversion_factors[i], level_amount, num_timesteps)

            #replace each variable name with a reference to python's equivalent variable,
            # but avoid overwriting variables with the same beginning
            # (e.g. overwrite THETAL variable when THETA get's processed)
            var_replacement = " variable_values['"+variable+"'] "
            #equation = equation.replace(var_to_replace, var_replacement)
            if not variable in variables_evaled:
                equation = re.sub('(\W)' + variable + '(\W)',r'\1' + var_replacement + r'\2', equation)
            variables_evaled.append(variable)
            i += 1

        new_values = eval(equation)
        return new_values

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
        panel_config = configparser.ConfigParser()
        panel_config.read(self.panel_filename)
        panel_type_config = configparser.RawConfigParser()
        panel_parent_dir = '/home/strike/clubb/postprocessing/pyplotgen/cases/' # TODO stop hardcoding this
        panel_type_filename = panel_config['defaults']['panel-type_file']
        panel_type_config.read(panel_parent_dir + panel_type_filename)

        all_plots_data = []
        panel_vars_str = panel_config['defaults']['model_vars']
        panel_vars = self.__getArrayFromString__(panel_vars_str)
        # num_lines = panel_vars.__len__()


        # TODO load these, don't hardcode them
        y_variable_name = panel_config['defaults']['y_axis_data']#"z"
        conversion_factors_str = panel_config['defaults']['conversion_factors']
        conversion_factors = self.__getArrayFromString__(conversion_factors_str)
        for i in range(len(conversion_factors)):
            conversion_factors[i] = float(conversion_factors[i])
        y_level_amount = 1
        y_num_timesteps = 1
        y_axis_values = self.__getValuesFromNc__(netcdf_data, y_variable_name, conversion_factors[0], y_level_amount, y_num_timesteps) #TODO conversion shouldn't be only 1 value

        # for i in range(0, num_lines):
            #TODO TEMP DEV TEST
            # if panel_vars[i] != "thetal":
            #     break
            # line_config = panel_type_config['default'][panel_vars[i]]
            # line_data = self.__parseAxisData__(line_config)
            # plot_index = self.__getLineToPlotIndex__(line_data)
            # TODO load these, don't hardcode them
            # nc_var_name = panel_type_config['default'][panel_vars[i]]
            # var_name_raw_str = panel_type_config[nc_var_name._name]['nc_name']
            # x_variable_name = line_data[plot_index][0] #var_name_raw_str.split('\n', re.MULTILINE)[0]
            # x_conversion_factor = 1.0
            # x_level_amount = 1
            # x_num_timesteps = 1

            # TODO load these, don't hardcode them
        start_x_value = int(panel_config['defaults']['x_min']) # Used to determine what x value to begin the graph at
        end_x_value = int(panel_config['defaults']['x_max']) # Used to determine what x value to end the graph at
        # TODO load these, don't hardcode them
        # plot_names = self.__getArrayFromString__(panel_type_config['default']['plotNames'])

        title = panel_config['defaults']['title'] #self.get_long_name(netcdf_data, x_variable_name)
        x_axis_title = panel_config['defaults']['x_label'] # x_variable_name + "[K]"
        y_axis_title = panel_config['defaults']['y_label']
        equation = panel_config['defaults']['vars_relationship']
        # Process inspired by Stackoverflow: https://stackoverflow.com/questions/45582344/extracting-data-from-netcdf-by-python
        x_axis_values = self.__calcVariableValues__(netcdf_data, panel_vars, equation, conversion_factors)# self.__getValuesFromNc__(netcdf_data, x_variable_name, x_conversion_factor, x_level_amount, x_num_timesteps)
        x_axis_values = self.__meanProfiles__(x_axis_values)
        x_axis_values = x_axis_values.reshape((-1,))
        (start_index, end_idx) = self.getStartEndIndex(x_axis_values, start_x_value, end_x_value) # Get the index values that correspond to the desired start/end x values

        x_axis_values = x_axis_values[start_index:end_idx + 1] # we use end_idx + 1 to ensure python uses the last value (off by one)
        y_axis_values = y_axis_values[start_index: end_idx + 1] # we use end_idx + 1 to ensure python uses the last value (off by one)

        plot_data = self.PlotData(x_values=x_axis_values, y_values=y_axis_values, title=title, x_title=x_axis_title, y_title=y_axis_title)
        all_plots_data.append(plot_data)

        return all_plots_data
