'''
:author: Nicolas Strike
:date: Mid 2019
'''


from pyplotgen.DataReader import DataReader

class VariableGroup:
    '''
    This is the parent PanelGroup class. All other panel groups
    should be created as a subclass of this one. Properties and
    methods common to each PanelGroup can be found here.


    A PanelGroup child defines Lineplots, Panels, and is responsible for
    calculating any 'calculated' variables from netcdf
    '''
    def __init__(self, ncdf_file):
        self.ncdf_file = ncdf_file
        self.panels = []

    def plot(self):
        '''
        Plots every panel in this group to the output folder specified
        in the pyplotgen launch parameters, unless the panel is blacklisted.
        
        :return: n/a
        '''
        for panel in self.panels:
            if not panel.blacklisted:
                panel.plot(self.ncdf_file)

    def get_var_from_ncdf(self, varname):
        '''
        Retrieve numerical data from netcdf

        :param varname:
        :return:
        '''
        data_reader = DataReader()
        return data_reader.getVarData(self.ncdf_file, varname)


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
