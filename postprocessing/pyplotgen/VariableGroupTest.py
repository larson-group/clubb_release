'''
:author: Nicolas Strike
:date: Mid 2019
'''

from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.VariableGroup import VariableGroup
from pyplotgen.Lineplot import Lineplot
from pyplotgen.VarnameConversions import CLUBB_TO_SAM


class VariableGroupBase(VariableGroup):
    '''
    This is a panel group used for testing the functionality of pyplotgen.
    '''
    def __init__(self, ncdf_files, sam_file = None):
        super(VariableGroupBase, self).__init__(ncdf_files)

        self.sam_file = sam_file
        self.ncdf_files = ncdf_files
        ### Panel config ###
        y_title = "Height [m]"
        self.averaging_start_time = 181 * 60 # minutes * seconds
        self.averaging_end_time = 240 * 60 # minutes * seconds
        self.timeseries_start_time = 0 * 60 # minutes * seconds
        self.timeseries_end_time  = 240 * 60 # minutes * seconds
        self.height_min_value = 0
        self.height_max_value = 1200

        self.z = NetCdfVariable('altitude', ncdf_files, avging_start_time=self.averaging_start_time, avging_end_time=self.averaging_end_time)
        self.z_min_idx, self.z_max_idx = self.__getStartEndIndex__(self.z.data, self.height_min_value, self.height_max_value)
        self.z.data = self.z.data[self.z_min_idx:self.z_max_idx]
        self.time = NetCdfVariable('time', ncdf_files)

        self.z_sam = NetCdfVariable('z', [sam_file], avging_start_time=self.averaging_start_time, avging_end_time=self.averaging_end_time)
        self.z_sam_min_idx, self.z_sam_max_idx = self.__getStartEndIndex__(self.z_sam.data, self.height_min_value, self.height_max_value)
        self.z_sam.data = self.z_sam.data[self.z_sam_min_idx:self.z_sam_max_idx]

        thlm_plots = super().getLinePlots('thlm', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, label="current clubb")
        thlm_sam = self.getThlmSamPlot()
        thlm_plots.append(thlm_sam)

        rtm_plots = super().getLinePlots('rtm', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, label="current clubb")
        rtm_sam = self.getRtmSamPlot()
        rtm_plots.append(rtm_sam)

        wpthlp_plots = super().getLinePlots('wpthlp', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, sam_file=sam_file, label="current clubb")
        wprtp_plots = super().getLinePlots('wprtp', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, sam_file=sam_file, label="current clubb")
        wp2_plots = super().getLinePlots('wp2', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, sam_file=sam_file, label="current clubb")
        den_wgt_avg_plots = super().getLinePlots('wp2', ncdf_files, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, sam_file=sam_file, override_panel_type=Panel.TYPE_TIMESERIES, label="current clubb")
        cloud_frac_plots = super().getLinePlots('cloud_frac', ncdf_files,sam_conv_factor=100, conversion_factor=100, averaging_start_time=self.averaging_start_time, averaging_end_time=self.averaging_end_time, sam_file=sam_file, label="current clubb")


        ### Plot / Panel pairs
        liq_pot_temp = Panel(thlm_plots, title="Liquid Water Potential Temperature, Theta l", dependant_title="thlm [K]")
        turbulent_flux__theta = Panel(wpthlp_plots, title="Turbulent Flux of Theta l", dependant_title='wpthlp [K m/s]')
        turbulent_flux_r = Panel(wprtp_plots, title="Turbulent Flux of rt", dependant_title='wprtp [(kg/kg) m/s]')
        wp2 = Panel(wp2_plots, title="Variance of w", dependant_title='wp2 [m^2/s^2]')
        den_wgt_avg = Panel(den_wgt_avg_plots, title="Density-Weighted Vertically Averaged wp2", dependant_title="wp2 [m^2/s^2]")
        rtm_panel = Panel(rtm_plots, title="Total Water Mixing Ratio, _t", dependant_title="rtm [kg/kg]")
        cld_panel = Panel(cloud_frac_plots, title="Cloud Fraction", dependant_title="cloud_frac[%]")

        ### Build panels ###
        self.panels = [ turbulent_flux__theta, turbulent_flux_r, wp2, den_wgt_avg, liq_pot_temp, rtm_panel, cld_panel]

    def plot(self):
        '''
        Plots every panel in this group to the output folder specified
        in the pyplotgen launch parameters, unless the panel is blacklisted.

        Uses the PanelGroup.plot() function
        :return: n/a
        '''
        super().plot()


    def getThlmSamPlot(self):
        '''
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return:
        '''
        sam_start_time = self.averaging_start_time / 60
        sam_end_time = self.averaging_end_time / 60

        thetal_ncdf = NetCdfVariable('THETAL', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        thetal = self.get_var_from_ncdf(thetal_ncdf)

        theta_ncdf = NetCdfVariable('THETA', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        theta = self.get_var_from_ncdf(theta_ncdf)

        tabs_ncdf = NetCdfVariable('TABS', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        tabs = self.get_var_from_ncdf(tabs_ncdf)

        qi_ncdf = NetCdfVariable('QI', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        qi = self.get_var_from_ncdf(qi_ncdf)

        thlm = thetal + (2500.4 * (theta/tabs) * (qi/1000))
        thlm = thlm[self.z_sam_min_idx:self.z_sam_max_idx]

        thlm_lineplot = Lineplot(thlm, self.z_sam.data, line_format="k-", label="LES output")
        return thlm_lineplot

    def getRtmSamPlot(self):
        '''
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return:
        '''
        sam_start_time = self.averaging_start_time / 60
        sam_end_time = self.averaging_end_time / 60

        qt_ncdf = NetCdfVariable('QT', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        qt = self.get_var_from_ncdf(qt_ncdf)

        qi_ncdf = NetCdfVariable('QI', [self.sam_file], 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time,timestep=1)
        qi = self.get_var_from_ncdf(qi_ncdf)

        rtm = (qt - qi) / 1000
        rtm = rtm[self.z_sam_min_idx:self.z_sam_max_idx]

        rtm_lineplot = Lineplot(rtm, self.z_sam.data, line_format="k-", label="LES output")
        return rtm_lineplot
