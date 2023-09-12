"""
:author: Nicolas Strike
:date: Mid 2019
"""

from src.VariableGroup import VariableGroup


class VariableGroupWs(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None, r408_dataset=None,
                 wrf_benchmark_dataset=None, 
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False, background_rcm=False):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "w variables"
        self.variable_definitions = [
             {'var_names':
                 {
                 'clubb': ['wpup2'],
                 'sam': [''],
                 'coamps': ['wpup2'],
                 'r408': ['wpup2'],
                 'hoc': ['wpup2'],
                 'e3sm': ['wpup2'],
                 'cam': ['wpup2'],
                 'wrf': ['wpup2'],
                 },
             },
             {'var_names':
                 {
                 'clubb': ['wpvp2'],
                 'sam': [''],
                 'coamps': ['wpvp2'],
                 'r408': ['wpvp2'],
                 'hoc': ['wpvp2'],
                 'e3sm': ['wpvp2'],
                 'cam': ['wpvp2'],
                 'wrf': ['wpvp2'],
                 },
             },
             {'var_names':
                 {
                 'clubb': [self.get_wpuiui],
                 'sam': [''],
                 'coamps': [''],
                 'r408': [''],
                 'hoc': [''],
                 'e3sm': [''],
                 'cam': [''],
                 'wrf': [''],
                 },
              'title': "3rd-order moment $\overline{w'u_i'u_i'}$",
              'axis_title': "$\overline{w'u_i'u_i'}$ [m3/s3]",
             },
             {'var_names':
                 {
                 'clubb': ['wp2up2'],
                 'sam': ['WP2UP2'],
                 'coamps': ['wp2up2'],
                 'r408': ['wp2up2'],
                 'hoc': ['wp2up2'],
                 'e3sm': ['wp2up2'],
                 'cam': ['wp2up2'],
                 'wrf': ['wp2up2'],
                 },
             },
             {'var_names':
                 {
                 'clubb': ['wp2vp2'],
                 'sam': ['WP2VP2'],
                 'coamps': ['wp2vp2'],
                 'r408': ['wp2vp2'],
                 'hoc': ['wp2vp2'],
                 'e3sm': ['wp2vp2'],
                 'cam': ['wp2vp2'],
                 'wrf': ['wp2vp2'],
                 },
             },
             {'var_names':
                 {
                 'clubb': [self.get_wp2uiui],
                 'sam': ['WP2TKE'],
                 'coamps': ['wp2vp2'],
                 'r408': ['wp2vp2'],
                 'hoc': ['wp2vp2'],
                 'e3sm': ['wp2vp2'],
                 'cam': ['wp2vp2'],
                 'wrf': ['wp2vp2'],
                 },
              'title': "4th-order moment $0.5\overline{w'^2u_i'u_i'}$",
              'axis_title': "$0.5\overline{w'^2u_i'u_i'}$ [m4/s4]",
             },
             {'var_names':
                {
                'clubb': ['wp4'],
                'sam': ['WP4'],
                'coamps': ['wp4'],
                'r408': ['wp4'],
                'hoc': ['wp4'],
                'e3sm': ['wp4'],
                'cam': ['wp4'],
                'wrf': ['wp4'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['wp2thlp'],
                'sam': ['WP2THLP'],
                'coamps': ['wp2thlp'],
                'r408': ['wp2thlp'],
                'hoc': ['wp2thlp'],
                'e3sm': ['wp2thlp'],
                'cam': ['wp2thlp'],
                'wrf': ['wp2thlp'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['wp2rtp'],
                'sam': ['WP2RTP'],
                'coamps': ['wp2qtp', 'wp2rtp'],
                'r408': ['wp2rtp'],
                'hoc': ['wp2rtp'],
                'e3sm': ['wp2rtp'],
                'cam': ['wp2rtp'],
                'wrf': ['wp2rtp'],
                },
             'sci_scale': -4,
            },
            {'var_names':
                {
                'clubb': ['wpthlp2'],
                'sam': ['WPTHLP2'],
                'coamps': ['wpthlp2'],
                'r408': ['wpthlp2'],
                'hoc': ['wpthlp2'],
                'e3sm': ['wpthlp2'],
                'cam': ['wpthlp2'],
                'wrf': ['wpthlp2'],
                },
             'sci_scale': 0,
            },
            {'var_names':
                {
                'clubb': ['wprtp2'],
                'sam': ['WPRTP2'],
                'coamps': ['wpqtp2', 'wprtp2'],
                'r408': ['wprtp2'],
                'hoc': ['wprtp2'],
                'e3sm': ['wprtp2'],
                'cam': ['wprtp2'],
                'wrf': ['wprtp2'],
                },
             'sci_scale': -7,
            },
            {'var_names':
                {
                'clubb': ['wprtpthlp'],
                'sam': ['WPRTPTHLP'],
                'coamps': ['wpqtpthlp','wprtpthlp'],
                'r408': ['wprtpthlp'],
                'hoc': ['wprtpthlp'],
                'e3sm': ['wprtpthlp'],
                'cam': ['wprtpthlp'],
                'wrf': ['wprtpthlp'],
                },
             'sci_scale': -4,
            },
            {'var_names':
                 {
                 'clubb': ['w_up_in_cloud'],
                 'sam': ['WSUP'],
                 'coamps': [''],
                 'r408': [''],
                 'hoc': [''],
                 'e3sm': [''],
                 'cam': [''],
                 'wrf': [''],
                 },
             },
            {'var_names':
                 {
                 'clubb': ['w_down_in_cloud'],
                 'sam': [''],
                 'coamps': [''],
                 'r408': [''],
                 'hoc': [''],
                 'e3sm': [''],
                 'cam': [''],
                 'wrf': [''],
                 },
             },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_benchmark_dataset=sam_benchmark_dataset, 
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets,
                         cam_datasets=cam_datasets, sam_datasets=sam_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm)

    def get_wpuiui(self, dataset_override=None):

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        wpup2, z, dataset = self.getVarForCalculations('wpup2', dataset)
        wpvp2, z, dataset = self.getVarForCalculations('wpvp2', dataset)
        wp3, z, dataset = self.getVarForCalculations('wp3', dataset)

        output = wpup2 + wpvp2 + wp3

        return output, z
    
    def get_wp2uiui(self, dataset_override=None):

        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_datasets
        wp2up2, z, dataset = self.getVarForCalculations('wp2up2', dataset)
        wp2vp2, z, dataset = self.getVarForCalculations('wp2vp2', dataset)
        wp4, z, dataset = self.getVarForCalculations('wp4', dataset)

        output = 0.5 * ( wp2up2 + wp2vp2 + wp4 )

        return output, z
