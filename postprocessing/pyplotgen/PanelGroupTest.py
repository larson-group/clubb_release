# from pyplotgen.NetCdfVariable import NetCdfVariable
from pyplotgen.DataReader import NetCdfVariable
from pyplotgen.Panel import Panel
from pyplotgen.PanelGroup import PanelGroup
from pyplotgen.Lineplot import Lineplot


class PanelGroupTest(PanelGroup):
    '''

    '''
    def __init__(self, ncdf_file):
        '''

        '''
        super(PanelGroupTest, self).__init__(ncdf_file)

        ### Panel config ###
        y_title = "Height [m]"
        start_time = 181
        end_time = 240
        independent_min_value = 0
        independent_max_value = 1200

        ### Independant Variable ###
        self.z = NetCdfVariable('z', ncdf_file, one_dimentional=True, start_time=start_time, end_time=end_time)
        independent_min_idx, independent_max_idx = super().__getStartEndIndex__(self.z.data, independent_min_value, independent_max_value)

        ### Dependant Variables ###
        self.thlm = NetCdfVariable('THLM', ncdf_file, start_time=start_time, end_time=end_time)
        self.rtm = NetCdfVariable('RTM', ncdf_file, start_time=start_time, end_time=end_time)

        self.__constrain__(independent_min_idx, independent_max_idx)

        ### Plot / Panel pairs
        liq_pot_temp_plot =  Lineplot(self.thlm, self.z, plot_name="clubb current")
        liq_pot_temp_panel = Panel(liq_pot_temp_plot, title="Liquid Water Potential Temperature, θl", x_title="thlm [K]", y_title=y_title)

        tot_water_mix_plot =  Lineplot(self.rtm, self.z, plot_name="clubb current")
        tot_water_mix_panel = Panel(tot_water_mix_plot, title="Total Water Mixing Ratio, @", x_title="rtm [kg/kg]", y_title=y_title)

        ### Build panels ###
        self.panels = [liq_pot_temp_panel, tot_water_mix_panel]

    def plot(self, netcdf_data):
        '''

        :param plotter:
        :param netcdf_data:
        :return:
        '''
        super().plot(netcdf_data)

    def __constrain__(self, min_index, max_index):
        '''
        Restrict all panels to the min_idx->max_idx set of data for each lineplot

        :param min_index:
        :param max_index:
        :return:
        '''
        self.z.data = self.z.data[min_index:max_index]
        self.rtm.data = self.rtm.data[min_index:max_index]
        self.thlm.data = self.thlm.data[min_index:max_index]
