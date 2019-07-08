

class Lineplot:
    '''
    This class holds information commonly passed around when using plotAll data.
    Data included is anything necessary to graph a line using pyplot, including the
    numerical values, line format style, and a label for the plot
    '''
    def __init__(self, x_data, y_data, line_format = "", label ="Unlabeled plot"):
        '''

        '''
        self.x = x_data
        self.y = y_data
        self.line_format = line_format
        self.label = label

