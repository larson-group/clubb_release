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
    '''
    def __init__(self, ncdf_file):
        super(VariableGroupBase, self).__init__(ncdf_file)

        ### Panel config ###
        y_title = "Height [m]"
        averaging_start_time = 181
        averaging_end_time = 240
        timeseries_start_time = 0
        timeseries_end_time  = 240
        height_min_value = 0
        height_max_value = 1200


        ### Plot / Panel pairs
        liq_pot_temp_panel = Panel('THLM', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Liquid Water Potential Temperature, Theta l", dependant_title="thlm [K]", independent_min_value=height_min_value, independent_max_value=height_max_value)
        tot_water_mix_panel = Panel('RTM', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Total Water Mixing Ratio, ", dependant_title="rtm [kg/kg]", independent_min_value=height_min_value, independent_max_value=height_max_value)
        turbulent_flux__theta_panel = Panel('WPTHLP', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Turbulent Flux of Theta l", dependant_title='wpthlp [K m/s]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        turbulent_flux_r_panel = Panel('WPRTP', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Turbulent Flux of rt", dependant_title='wprtp [(kg/kg) m/s]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        wp3_panel = Panel('WP3', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Third order moment of w", dependant_title='wp3 [m^3/s^3]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        wp2_panel = Panel('WP2', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Variance of w", dependant_title='wp2 [m^2/s^2]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        # thlp2_panel = Panel(thlp2_plot, title="Variance of theta l", x_title='thlp2 [K^2]', y_title=y_title, height_min_value=height_min_value, height_max_value=height_max_value)
        rtp2_panel = Panel('RTP2', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Variance of r t", dependant_title='rtp2 [(kg/kg)^2] x 10^-6', independent_min_value=height_min_value, independent_max_value=height_max_value)
        wm_panel = Panel('WM', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Vertical Wind Component, w (subsidence)", dependant_title='wm [m/s]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        # um_panel = Panel('UM', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Zonal Wind Component, u", dependant_title='um [m/s]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        # vm_panel = Panel('VM', ncdf_file, averaging_start_time=averaging_start_time, averaging_end_time=averaging_end_time, title="Zonal Wind Component, u", dependant_title='um [m/s]', independent_min_value=height_min_value, independent_max_value=height_max_value)
        #
        #
        liq_water_path = Panel('CWP', ncdf_file, panel_type=Panel.TYPE_TIMESERIES, title="Liquid Water Path", dependant_title="lwp [kg/m^2]", independent_min_value=timeseries_start_time, independent_max_value=timeseries_end_time)
        den_wgt_avg = Panel('W2', ncdf_file, averaging_start_time=averaging_start_time, panel_type=Panel.TYPE_TIMESERIES, averaging_end_time=averaging_end_time, title="Density-Weighted Vertically Averaged wp2", dependant_title="wp2 [m^2/s^2]", independent_min_value=timeseries_start_time, independent_max_value=timeseries_end_time)


        ### Build panels ###
        self.panels = [liq_pot_temp_panel, tot_water_mix_panel, turbulent_flux__theta_panel, turbulent_flux_r_panel,
                       wp2_panel, wp3_panel, liq_water_path, den_wgt_avg]


    # def get_wp2(self):'''
    #
    #     :return:
    #     '''

    def plot(self):
        '''
        Plots every panel in this group to the output folder specified
        in the pyplotgen launch parameters, unless the panel is blacklisted.

        Uses the PanelGroup.plot() function
        :return: n/a
        '''
        super().plot()


