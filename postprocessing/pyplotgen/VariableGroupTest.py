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
        self.z = NetCdfVariable('altitude', ncdf_files, avging_start_time=self.averaging_start_time,
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

        ### Initialize plots in order to make panels later ###
        # Step 1) Generate current clubb plots using super().getVarLinePlots()
        # Step 2) If plotting SAM and variable needs to be calculated, create and call a calulating function
        # Step 3) Append the SAM plotAll from step 2 into the list generated from step 1
        thlm_plots = super().getVarLinePlots('thlm', ncdf_files, averaging_start_time=self.averaging_start_time,
                                             averaging_end_time=self.averaging_end_time, label="current clubb",
                                             line_format='r--')
        if sam_file != None:
            thlm_sam = self.getThlmSamPlot()
            thlm_plots.append(thlm_sam)

        rtm_plots = super().getVarLinePlots('rtm', ncdf_files, averaging_start_time=self.averaging_start_time,
                                            averaging_end_time=self.averaging_end_time, label="current clubb",
                                            line_format='r--')
        if sam_file != None:
            rtm_sam = self.getRtmSamPlot()
            rtm_plots.append(rtm_sam)

        # SAM variables that don't require calculation beyond a conversion factor can be plotted by defining the
        # variable name conversion in the file VarnameConversions.py
        wpthlp_plots = super().getVarLinePlots('wpthlp', ncdf_files, averaging_start_time=self.averaging_start_time,
                                               averaging_end_time=self.averaging_end_time, sam_file=sam_file,
                                               label="current clubb", line_format='r--')
        wprtp_plots = super().getVarLinePlots('wprtp', ncdf_files, averaging_start_time=self.averaging_start_time,
                                              averaging_end_time=self.averaging_end_time, sam_file=sam_file,
                                              label="current clubb", line_format='r--')
        wp2_plots = super().getVarLinePlots('wp2', ncdf_files, averaging_start_time=self.averaging_start_time,
                                            averaging_end_time=self.averaging_end_time, sam_file=sam_file,
                                            label="current clubb", line_format='r--')
        den_wgt_avg_plots = super().getVarLinePlots('wp2_vert_avg', ncdf_files['sfc'],
                                                    averaging_start_time=self.averaging_start_time,
                                                    averaging_end_time=self.averaging_end_time,
                                                    override_panel_type=Panel.TYPE_TIMESERIES, label="current clubb",
                                                    sam_conv_factor=1 / 1000, line_format='r--')
        cloud_frac_plots = super().getVarLinePlots('cloud_frac', ncdf_files, sam_conv_factor=100, conversion_factor=100,
                                                   averaging_start_time=self.averaging_start_time,
                                                   averaging_end_time=self.averaging_end_time, sam_file=sam_file,
                                                   label="current clubb", line_format='r--')

        ### Initialize Panels using plots from above ###
        liq_pot_temp = Panel(thlm_plots, title="Liquid Water Potential Temperature, Theta l", dependant_title="thlm [K]")
        turbulent_flux__theta = Panel(wpthlp_plots, title="Turbulent Flux of Theta l", dependant_title='wpthlp [K m/s]')
        turbulent_flux_r = Panel(wprtp_plots, title="Turbulent Flux of rt", dependant_title='wprtp [(kg/kg) m/s]')
        wp2 = Panel(wp2_plots, title="Variance of w", dependant_title='wp2 [m^2/s^2]')
        den_wgt_avg = Panel(den_wgt_avg_plots, title="Density-Weighted Vertically Averaged wp2",
                            dependant_title="wp2 [m^2/s^2]", panel_type=Panel.TYPE_TIMESERIES)
        rtm_panel = Panel(rtm_plots, title="Total Water Mixing Ratio, _t", dependant_title="rtm [kg/kg]")
        cld_panel = Panel(cloud_frac_plots, title="Cloud Fraction", dependant_title="cloud_frac[%]")

        ### Build panels list ###
        self.panels = [turbulent_flux__theta, turbulent_flux_r, wp2, den_wgt_avg, liq_pot_temp, rtm_panel, cld_panel]

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
