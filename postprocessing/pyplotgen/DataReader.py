import os
class DataReader():
    '''
    This class is responsible for handling input files. Given a
    input file, it determines what format it is in (grads vs
    netcdf) and then uses the correct helper file (e.g.
    NetcdfReader.py) to load the data into python.
    '''
    def __init__(self):
        print("DataReader not implemented")

    def load_folder(self, folder_path):
        '''
        Finds all dataset files in a given folder and loads
        them using the appropriet helper class.
        :param folder_path: The path of the folder to be loaded
        :return: None
        '''
        data_files = []
        for root, dirs, files in os.walk(folder_path):
            for filename in files:
                abs_filename = os.path.abspath(os.path.join(root, filename))
                data_files.append(abs_filename)