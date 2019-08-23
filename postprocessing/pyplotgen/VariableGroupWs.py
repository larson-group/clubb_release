'''
:author: Nicolas Strike
:date: Mid 2019
'''

from DataReader import NetCdfVariable
from Line import Line
from VariableGroup import VariableGroup


class VariableGroupWs(VariableGroup):

    def __init__(self, ncdf_datasets, case, sam_file=None, coamps_file=None, r408_dataset=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "w variables"
        self.variable_definitions = [
            {'aliases': ['wp4', 'WP4']},
            {'aliases': ['wp2thlp', 'WP2THLP'], 'fill_zeros':'True'},
            {'aliases': ['wp2rtp', 'WP2RTP', 'wp2qtp']},
            {'aliases': ['wpthlp2', 'WPTHLP2']},
            {'aliases': ['wprtp2', 'WPRTP2', 'wpqtp2']},
            {'aliases': ['wprtpthlp', 'WPRTPTHLP', 'wpqtpthlp']},
            {'aliases': ['wp2thvp','WP2THVP']},  # LES
            {'aliases': ['rc_coef * wprcp'], 'fallback_func': self.get_rc_coef_zm_X_wprcp_clubb_line,
                'title': 'Contribution of Cloud Water Flux to wprcp', 'axis_title': 'rc_coef_zm * wprcp [K m/s]'}

        ]
        super().__init__(ncdf_datasets, case, sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset)



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