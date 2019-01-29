import os
from netCDF4 import Dataset
import numpy as np
from collections import OrderedDict

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
        '''
        self.nc_filenames = []
        self.grads_dat_filenames = []
        self.grads_ctl_filenames = []
        self.nc_datasets = []

    def cleanup(self):
        '''
        This is the cleanup method. This is called on the instance's destruction
        to unallocate resources that may be held (e.g. dataset files).
        :return:
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


    # def isFunction(value):
    #     '''
    #     Determine if an element of a budget is a function
    #     :return True if it is a function, false if not:
    #     '''
    #     # logger.info('isFunction')
    #     isFunc = False
    #     if '+' in value:
    #         isFunc = True
    #     elif '-' in value:
    #         isFunc = True
    #     elif '*' in value:
    #         isFunc = True
    #     elif '/' in value:
    #         isFunc = True
    #     else:
    #         isFunc = False
    #     return isFunc


    # def seperateDataAndFunctions(self, budget):
    #     '''
    #     collect and separate data variables and function variables
    #     calculate mean of data over time -> array containing means for all height levels z0<=h<=z1
    #     :param budget: The budget to seperate values out of
    #     :return:
    #     '''
    #     budgets_data = []
    #     functions = []
    #
    #     for j in range(len(budget)):
    #         if not self.isFunction(budget[j][2]):
    #             # grap data of each variable that is not a function
    #             # logger.info("Grap data of: %s", budget[j][0])
    #             value = self.mean_profiles(self.get_budgets_from_nc(nc, budget[j][2], budget[j][3], n, t), idx_t0, idx_t1, idx_z0, idx_z1)
    #             if np.any(np.isnan(value)) or np.any(value <= -1000):
    #                 # if there are no values for the variable
    #                 value = np.zeros(n)
    #                 # logger.warning("Could not find the variable %s of %s", budget[j][0], bv.sortPlots[i])
    #             budgets_data.append([budget[j][0], budget[j][1], budget[j][2], budget[j][4], value])
    #         else:
    #             # save a function for an evaluation
    #             #functions.append([budget[j][2]])
    #             functions.append([budget[j][0], budget[j][1], budget[j][2], budget[j][4]])
    #             #func_names.append(budget[j][0])
    #     return budgets_data, functions

    def load_folder(self, folder_path):
        '''
        Finds all dataset files in a given folder and loads
        them using the appropriet helper class.
        :param folder_path: The path of the folder to be loaded
        :return: An array of datasets
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

    # def mean_profiles(var, idx_t0, idx_t1, idx_z0, idx_z1):
    #     # logger.info('mean_profiles')
    #     """
    #     Input:
    #       var    -- time x height array of some property
    #       idx_t0 -- Index corrosponding to the beginning of the averaging interval
    #       idx_t1 -- Index corrosponding to the end of the averaging interval
    #       idx_z0 -- Index corrosponding to the lowest model level of the averaging interval
    #       idx_z1 -- Index corrosponding to the highest model level of the averaging interval
    #
    #     Output:
    #       var    -- time averaged vertical profile of the specified variable
    #     """
    #
    #     return np.nanmean(var[idx_t0:idx_t1,idx_z0:idx_z1],axis=0)


    #
