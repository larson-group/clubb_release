from pyplotgen.DataReader import NetCdfVariable


class Line:
    '''
    This class holds information commonly passed around when using plotAll data.
    Data included is anything necessary to graph a line using pyplot, including the
    numerical values, line format style, and a label for the plot
    '''
    def __init__(self, x_data, y_data, line_format = "", label ="Unlabeled plot"):
        '''

        '''
        if len(x_data) != len(y_data):
            raise ValueError("The size of x is not the same as the size of Y. " + str(len(x_data)) + " (x) vs " + str(len(y_data)) + "(y)")
        if isinstance(x_data, NetCdfVariable):
            x_data = x_data.data
        if isinstance(y_data, NetCdfVariable):
            y_data = y_data.data
        self.x = x_data
        self.y = y_data
        self.line_format = line_format
        self.label = label

