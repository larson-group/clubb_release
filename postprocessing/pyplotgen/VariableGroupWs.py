'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from VariableGroup import VariableGroup
from Line import Line


class VariableGroupWs(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "w variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            {'clubb_name': 'wp4', 'sam_name': 'WP4'},
            {'clubb_name': 'wp2thlp', 'sam_name': 'WP2THLP'},
            {'clubb_name': 'wp2rtp', 'sam_name': 'WP2RTP'},
            {'clubb_name': 'wpthlp2', 'sam_name': 'WPTHLP2'},
            {'clubb_name': 'wprtp2', 'sam_name': 'WPRTP2'},
            {'clubb_name': 'wprtpthlp', 'sam_name': 'WPRTPTHLP'},
            {'clubb_name': 'wp2thvp', 'sam_name': 'WP2THVP'},  # LES
            {'clubb_name': 'rc_coef * wprcp', 'fallback_func': self.get_rc_coef_zm_X_wprcp_clubb_line,
                'title': 'Contribution of Cloud Water Flux to wprcp', 'axis_title': 'rc_coef_zm * wprcp [K m/s]'}

        ]
        super().__init__(ncdf_datasets, case, sam_file)



    def get_rc_coef_zm_X_wprcp_clubb_line(self):
        '''
        Calculates the Contribution of Cloud Water Flux
        to wprcp using the equation
        rc_coef_zm * wprcp
        :return: Line representing rc_coef_zm .* wprcp
        '''
        z_ncdf = NetCdfVariable('altitude', self.ncdf_files['zm'], 1)

        rc_coef_zm_ncdf = NetCdfVariable('rc_coef_zm', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        rc_coef_zm_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        rc_coef_zm = rc_coef_zm_ncdf.data

        wprcp_ncdf = NetCdfVariable('wprcp', self.ncdf_files['zm'], 1, start_time=self.start_time, end_time=self.end_time, fill_zeros=True)
        wprcp_ncdf.constrain(self.height_min_value, self.height_max_value, data=z_ncdf.data)
        wprcp = wprcp_ncdf.data

        output = rc_coef_zm * wprcp

        z_ncdf.constrain(self.height_min_value, self.height_max_value)
        output = Line(output, z_ncdf.data, line_format='b-', label='current clubb')
        return output