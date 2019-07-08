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

    def __init__(self, ncdf_datasets, case, sam_file=None):
        '''

        :param ncdf_datasets:
        :param case:
        :param sam_file:
        '''
        self.name = "base variables"
        # TODO Support fill_zeros
        self.variable_definitions = [
            # {'clubb_name': 'thlm', 'sam_calc': self.getThlmSamPlot},
            # {'clubb_name': 'rtm', 'sam_calc': self.getRtmSamPlot},
            # # {'clubb_name': 'wpthlp', 'sam_name': 'WPTHLP'},
            # {'clubb_name': 'wpthlp', 'sam_calc': self.getWpthlpSamPlot},
            # {'clubb_name': 'wprtp', 'sam_name': 'WPRTP'},
            # {'clubb_name': 'rcm', 'sam_name': "QCL", 'sam_conv_factor': 1/1000},
            # {'clubb_name': 'wp3', 'sam_name': 'W3'},
            # {'clubb_name': 'thlp2', 'sam_name': 'THLP2'},
            # {'clubb_name': 'rtp2', 'sam_name': 'RTP2'},
            # {'clubb_name': 'wm', 'sam_name': 'WOBS'},
            # {'clubb_name': 'um', 'sam_name': 'U'},
            # {'clubb_name': 'vm', 'sam_name': 'V'},
            # {'clubb_name': 'wp3', 'sam_name': 'W3'},
            # {'clubb_name': 'wp2', 'sam_name': 'W2'},
            {'clubb_name': 'wp2_vert_avg', 'sam_name': 'CWP', 'type': Panel.TYPE_TIMESERIES}
            # {'clubb_name': 'cloud_frac', 'sam_name': 'CLD'},
            # {'clubb_name': 'upwp', 'sam_name': 'UW'},
            # {'clubb_name': 'vpwp', 'sam_name': 'VW'},
            # {'clubb_name': 'up2', 'sam_name': 'U2'},
            # {'clubb_name': 'tau_zm'},
            # {'clubb_name': 'vp2', 'sam_name': 'V2'},
            # {'clubb_name': 'Lscale'},
            # {'clubb_name': 'wpthvp'},
            # {'clubb_name': 'rtpthlp', 'sam_name': 'RTPTHLP'},
            # {'clubb_name': 'rtp3', 'sam_name': 'RTP3'},
            # {'clubb_name': 'radht'},
            # {'clubb_name': 'Skw_zt'},
            # {'clubb_name': 'thlp3', 'sam_name': 'THLP3'},
            # {'clubb_name': 'rtpthvp'},
            # {'clubb_name': 'Skrt_zt'},
            # {'clubb_name': 'Skthl_zt'},
            # {'clubb_name': 'corr_w_chi_2'},
            # {'clubb_name': 'corr_chi_eta_2'},
            # {'clubb_name': 'rcp2'},
            # {'clubb_name': 'thlpthvp'},

        ]
        super().__init__(ncdf_datasets, case, sam_file)

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

    def getWpthlpSamPlot(self):
        '''

        :return:
        '''
        sam_start_time = self.averaging_start_time / 60
        sam_end_time = self.averaging_end_time / 60

        wpthlp_ncdf = NetCdfVariable('WPTHLP', self.sam_file, 1, avging_start_time=sam_start_time, avging_end_time=sam_end_time)
        wpthlp = wpthlp_ncdf.data
        wpthlp = wpthlp[self.z_sam_min_idx:self.z_sam_max_idx]


        wpthlp_lineplot = Lineplot(wpthlp, self.z_sam.data, line_format='k-', label='LES output')
        return wpthlp_lineplot