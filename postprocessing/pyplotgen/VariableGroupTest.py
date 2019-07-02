'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Lineplot import Lineplot


class VariableGroupBase(VariableGroup):
    '''
    This is a panel group used for testing the functionality of pyplotgen.
    It contains a set of common panels being used for representing the majority
    of panels.
    '''

    def __init__(self, ncdf_files, sam_file=None):
        super(VariableGroupBase, self).__init__(ncdf_files, sam_file)

        self.name = "base variables"

        ### Panel config ###
        sec_per_min = 60
        self.averaging_start_time = 181 * sec_per_min  # minutes * seconds
        self.averaging_end_time = 240 * sec_per_min  # minutes * seconds
        self.timeseries_start_time = 0 * sec_per_min  # minutes * seconds
        self.timeseries_end_time = 240 * sec_per_min  # minutes * seconds
        self.height_min_value = 0
        self.height_max_value = 1200

        ### Initialize Height ###
        self.z = NetCdfVariable('altitude', ncdf_files['zm'], avging_start_time=self.averaging_start_time,
                                avging_end_time=self.averaging_end_time)
        self.z_min_idx, self.z_max_idx = self.__getStartEndIndex__(self.z.data, self.height_min_value,
                                                                   self.height_max_value)
        self.z.data = self.z.data[self.z_min_idx:self.z_max_idx]

        ### Initialize Time ###
        sec_to_min = 1 / sec_per_min
        self.time = NetCdfVariable('time', ncdf_files, conversion_factor=sec_to_min)

        ### Initialize Sam Height ###
        if sam_file != None:
            self.z_sam = NetCdfVariable('z', sam_file, avging_start_time=self.averaging_start_time,
                                        avging_end_time=self.averaging_end_time)
            self.z_sam_min_idx, self.z_sam_max_idx = self.__getStartEndIndex__(self.z_sam.data, self.height_min_value,
                                                                               self.height_max_value)
            self.z_sam.data = self.z_sam.data[self.z_sam_min_idx:self.z_sam_max_idx]

        variables = [
            {'clubb_name': 'thlm', 'title': 'Liquid Water Potential Temperature, Theta l', 'axis_title':'thlm [K]', 'sam_calc':self.getThlmSamPlot},
            {'clubb_name': 'rtm', 'title': 'Turbulent Flux of Theta l', 'axis_title': 'wpthlp [K m/s]', 'sam_calc':self.getRtmSamPlot},
            {'clubb_name': 'wpthlp', 'title': 'Turbulent Flux of Theta l', 'axis_title': 'wpthlp [K m/s]'},
            {'clubb_name': 'wprtp', 'title': 'Turbulent Flux of rt', 'axis_title': 'wprtp [(kg/kg) m/s]'},
            {'clubb_name': 'wp2', 'title': 'Variance of w', 'axis_title': 'wp2 [m^2/s^2]'},
            {'clubb_name': 'wp2_vert_avg', 'title': 'Density-Weighted Vertically Averaged wp2', 'axis_title': 'wp2 [m^2/s^2]', 'type':Panel.TYPE_TIMESERIES},
            {'clubb_name': 'cloud_frac', 'title': 'Cloud Frac', 'axis_title': 'cloud_frac[%]'}

        ]

        for variable in variables:
            super().addClubbVariable(variable)
        super().generatePanels()


    def plot(self):
        '''
        Plots every panel in this group to the output folder specified
        in the pyplotgen launch parameters, unless the panel is blacklisted.

        Uses the PanelGroup.plotAll() function
        :return: n/a
        '''
        super().plot()

    def getThlmSamPlot(self):
        '''
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sec_per_min = 60
        sam_start_time = self.averaging_start_time / sec_per_min
        sam_end_time = self.averaging_end_time / sec_per_min

        thetal_ncdf = NetCdfVariable('THETAL', self.sam_file, 1, avging_start_time=sam_start_time,
                                     avging_end_time=sam_end_time)
        thetal = thetal_ncdf.data

        theta_ncdf = NetCdfVariable('THETA', self.sam_file, 1, avging_start_time=sam_start_time,
                                    avging_end_time=sam_end_time)
        theta = theta_ncdf.data

        tabs_ncdf = NetCdfVariable('TABS', self.sam_file, 1, avging_start_time=sam_start_time,
                                   avging_end_time=sam_end_time)
        tabs = tabs_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time, fill_zeros=True)
        qi = qi_ncdf.data

        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        thlm = thlm[self.z_sam_min_idx:self.z_sam_max_idx]

        thlm_lineplot = Lineplot(thlm, self.z_sam.data, line_format="k-", label="LES output")
        return thlm_lineplot

    def getRtmSamPlot(self):
        '''
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        '''
        sam_start_time = self.averaging_start_time / 60
        sam_end_time = self.averaging_end_time / 60

        qt_ncdf = NetCdfVariable('QT', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time)
        qt = qt_ncdf.data

        qi_ncdf = NetCdfVariable('QI', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time, fill_zeros=True)
        qi = qi_ncdf.data

        rtm = (qt - qi) / 1000
        rtm = rtm[self.z_sam_min_idx:self.z_sam_max_idx]

        rtm_lineplot = Lineplot(rtm, self.z_sam.data, line_format="k-", label="LES output")
        return rtm_lineplot
