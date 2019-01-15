import os
from netCDF4 import Dataset

class DataReader():
    '''
    This class is responsible for handling input files. Given a
    input file, it determines what format it is in (grads vs
    netcdf) and then uses the correct helper file (e.g.
    NetcdfReader.py) to load the data into python.
    '''


    def __init__(self):
        '''
        This is the object constructer. Initalizes class variables and data
        '''
        print("DataReader not implemented")
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
        print("Variables: " + str(dataset.variables))
        #print(filename + ": " + str(dataset))
        return dataset


    def load_folder(self, folder_path):
        '''
        Finds all dataset files in a given folder and loads
        them using the appropriet helper class.
        :param folder_path: The path of the folder to be loaded
        :return: None
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
        # print("Files loaded:\n\n----nc files----\n" + str(self.nc_filenames) + "\n\n----dat files----\n" + str(self.grads_dat_filenames) +
        #       "\n\n----ctl files----\n" + str(self.grads_ctl_filenames))

