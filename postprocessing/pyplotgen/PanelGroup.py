from pyplotgen.DataReader import DataReader


class PanelGroup:
    '''

    '''
    def __init__(self, ncdf_file):
        '''

        '''
        # self.panel_type = panel_type
        # self.blacklisted_panels = blacklisted_panels
        self.ncdf_file = ncdf_file
        self.panels = []

    def plot(self, netcdf_data):
        '''

        :param plotter:
        :return:
        '''
        for panel in self.panels:
            if not panel.blacklisted:
                panel.plot(netcdf_data)

    def get_var_from_ncdf(self, varname):
        '''

        :param varname:
        :return:
        '''
        data_reader = DataReader()
        return data_reader.getVarData(self.ncdf_file, varname, 1.0, 1, 1) #TODO hardcoded


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

    def __constrain__(self, min_index, max_index):
        '''
        Restrict all panels to the min_idx->max_idx set of data for each lineplot

        :param min_index:
        :param max_index:
        :return:
        '''
        pass

    # def get_xy_data_from_ncdf(self, x_name, y_name): # 99999 is an arbitrarily chosen large number, replace it
    #     '''
    #
    #     :param varname:
    #     :return:
    #     '''
    #     data_reader = DataReader()
    #     return data_reader.getPlotsData(self.ncdf_file, x_name, y_name, [1], 0, 1200)