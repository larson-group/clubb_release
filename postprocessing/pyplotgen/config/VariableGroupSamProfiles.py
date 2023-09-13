'''
:author: Steffen Domke
:date: Mid 2019
TODO:   - Figure out how to include standard plots in VariableGroupBase (labels etc.)
'''
import numpy as np

from src.VariableGroup import VariableGroup


class VariableGroupSamProfiles(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, sam_benchmark_dataset=None, coamps_benchmark_dataset=None,
                 wrf_benchmark_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None, priority_vars=False,
                 background_rcm=False, background_rcm_folder=None):
        """

        :param clubb_datasets:
        :param case:
        :param sam_benchmark_dataset:
        """
        self.name = "sam profile variables"
        
        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)
        
        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000
        
        self.variable_definitions = [
            {'var_names':
                {
                'clubb': [''],
                'sam': ['WPPP'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"Covariance of w' and p'/$\rho$",
             'axis_title': r"$\mathrm{\overline{w'p'/\rho}}$ $\mathrm{\left[m^3\,s^{-3}\right]}$",
             },
            {'var_names':
                {
                'clubb': [''],
                'sam': ['WP2PP'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"Covariance of w'^2 and p'/$\rho$",
             'axis_title': r"$\mathrm{\overline{w'^2p'/\rho}}$ $\mathrm{\left[m^4\,s^{-4}\right]}$",
             },
            {'var_names':
                {
                'clubb': [''],
                'sam': ['WP2UP2'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"4th-order moment, w'^2u'^2",
             'axis_title': r"$\mathrm{\overline{w'^2u'^2}}$ $\mathrm{\left[m^4\,s^{-4}\right]}$",
             },
            {'var_names':
                {
                'clubb': [''],
                'sam': ['WP2VP2'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"4th-order moment, w'^2v'^2",
             'axis_title': r"$\mathrm{\overline{w'^2v'^2}}$ $\mathrm{\left[m^4\,s^{-4}\right]}$",
},
            {'var_names':
                {
                'clubb': [''],
                'sam': ['WP2TKE'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"4th-order moment, w'^2(0.5*(u'^2+v'^2+w'^2))",
             'axis_title': r"$\mathrm{\overline{w'^2\;0.5u_i'u_i'}}$ $\mathrm{\left[m^4\,s^{-4}\right]}$",
             },
            {'var_names':
                {
                'clubb': [''],
                'sam': ['PPTKE'],
                'coamps': [''],
                'r408': [''],
                'hoc': [''],
                'e3sm': [''],
                'cam': [''],
                },
             'title': r"3rd-order moment, p'/$\rho$(0.5*(u'^2+v'^2+w'^2))",
             'axis_title': r"$\mathrm{\overline{p'/\rho\;0.5u_i'u_i'}}$ $\mathrm{\left[m^4\,s^{-4}\right]}$",
             },
            # CORR(W, THETA_L)
            {'var_names':
                {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_calc': self.getCorrWpThlpCalc, 'sam_conv_factor': 1,
             'title': r"Corr(w,$\mathrm{\theta_l}$)",
             'axis_title': r"Correlation "+
                           r"$\mathrm{\overline{w'\theta_l'} / \sqrt{\overline{w'^2}\;\overline{\theta_l'^2}}}$"+
                           r" $\mathrm{\left[-\right]}$",
             'legend_label': r"$\mathrm{Corr(w,\theta_l)}$",
            },
            # CORR(W, R_T)
            {'var_names':
                {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_calc': self.getCorrWpRtpCalc, 'sam_conv_factor': 1,
             'title': r"$\mathrm{Corr(w,r_t)}$",
             'axis_title': r"Correlation "+
                           r"$\mathrm{\overline{w'r_t'} / \sqrt{\overline{w'^2}\;\overline{r_t'^2}}}$"+
                           r" $\mathrm{\left[-\right]}$",
             'legend_label': r"$\mathrm{Corr(w,r_t)}$",
            },
            # cloudliq_frac_em6
            {'var_names':
                {
                'clubb': [],
                'sam': ['cloudliq_frac_em6'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "Cloud Liquid Fraction",
             'axis_title': r"Cloud liquid fraction $\mathrm{\left[\frac{\%}{100}\right]}$",
             'legend_label': 'cloudliq_frac_em6',
            },
            # WCLD unweighted
            {'var_names':
                {
                'clubb': [],
                'sam': ['WCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"In-cloud mean wind, $\overline{W}^{cld}$",
             'axis_title': r"In-cloud mean wind, "+
                           r"$\overline{w}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
            },
            # UCLD unweighted
            {'var_names':
                {
                'clubb': [],
                'sam': ['UCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"In-cloud mean wind, $\overline{u}^{cld}$",
             'axis_title': r"In-cloud mean wind, "+
                           r"$\overline{u}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
            },
            # VCLD unweighted
            {'var_names':
                {
                'clubb': [],
                'sam': ['VCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"In-cloud mean wind, $\overline{v}^{cld}$",
             'axis_title': r"In-cloud mean wind, "+
                           r"$\overline{v}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
            },
            # CORR(U, W)
            {'var_names':
                {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_calc': self.getUpWpCorrCalc, 'sam_conv_factor': 1,
             'title': r"Corr(u,w)",
             'axis_title': r"Correlation "+
                           r"$\mathrm{\overline{u'w'} / \sqrt{\overline{u'^2}\;\overline{w'^2}}\ \left[-\right]}$",
             'legend_label': r"Corr(u,w)",
            },
            # CORR(V, W)
            {'var_names':
                {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_calc': self.getVpWpCorrCalc, 'sam_conv_factor': 1,
             'title': r"Corr(v,w)",
             'axis_title': r"Correlation "+
                           r"$\mathrm{\overline{v'w'} / \sqrt{\overline{v'^2}\;\overline{w'^2}}\ \left[-\right]}$",
             'legend_label': r"Corr(v,w)",
            },
            ### PRECIPITATION
            
            ## Rain Water Mixing Ratio (only MICRO_DRIZZLE, MICRO_2MOM, and MICRO_SAM1MOM)
            # QR
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': ['QR'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_conv_factor': 1/1000,
            # 'title': "Rain Water Mixing Ratio",
            # 'axis_title': r"$\mathrm{\overline{q_r}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': r"QR",
            #},
            ## QR_IP (only MICRO_DRIZZLE, MICRO_M2005_UWM)
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': ['qrainm_ip'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_conv_factor': 1,
            # 'title': "Rain Water Mixing Ratio in Rain",
            # 'axis_title': r"$\mathrm{q_r^{ip}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'QR_IP',
            #},
            ## QR'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': ['qrainp2'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Rain Water Mixing Ratio",
            # 'axis_title': r"$\mathrm{\overline{q_r'^2}}$ $\mathrm{\left[kg^2\,kg^{-2}\right]}$",
            # 'legend_label': 'qrainp2',
            #},
            ## QR'^2_IP / QR_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': self.getQRP2_QRIP, 'sam_conv_factor': 1,
            # 'title': "Within-rain Variance of Rain Water Mixing Ratio",
            # 'axis_title': r"$\mathrm{\overline{q_r'^2}^{ip} / \overline{q_r^{ip}}^2}$ [-]",
            # 'legend_label': 'QRP2_QRIP',
            #},
            ### Rain Drop Number Concentration
            ## NR
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Rain Drop Concentration",
            # 'axis_title': r"Nrm $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nrm',
            #},
            ## NR_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Rain Drop Concentration in Rain",
            # 'axis_title': r"Nrm_ip $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nrm_IP',
            #},
            ## NR'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance\nof Rain Drop Concentration",
            # 'axis_title': r"Nrp2 $\mathrm{\left[num^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Nrp2',
            #},
            ## NR'^2_IP / NR_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-rain Variance of Rain Drop Concentration",
            # 'axis_title': r"Nrp2_ip / Nrm_ip^2 [-]",
            # 'legend_label': 'Nrp2_NrmIP',
            #},
            ### Cloud Droplet Number Concentration
            ## NC
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Droplet Number Concentration",
            # 'axis_title': r"Ncm $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Ncm',
            #},
            ## NC_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Droplet Number Concentration in Cloud",
            # 'axis_title': r"Ncm_ip $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Ncm_IP',
            #},
            ## NC'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Cloud Droplet Number Concentration",
            # 'axis_title': r"Ncp2 $\mathrm{\left[\#^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Ncp2',
            #},
            ## NC'^2_IP / NC_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-cloud Variance of Cloud Droplet Number Concentration",
            # 'axis_title': r"Ncp2_ip / Ncm_ip^2 [-]",
            # 'legend_label': 'Ncp2_NcmIP',
            #},
            ### Graupel Number Concentration
            ## NG
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Graupel Number Concentration",
            # 'axis_title': r"Ngm $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Ngm',
            #},
            ## NG_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Graupel Number Concentration in Graupel",
            # 'axis_title': r"Ngm_ip $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Ngm_IP',
            #},
            ## NG'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Graupel Number Concentration",
            # 'axis_title': r"Ngp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Ngp2',
            #},
            ## NG'^2_IP / NG_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-graupel Variance\nof Graupel Number Concentration",
            # 'axis_title': r"Ngp2_ip / Ngm_ip^2 [-]",
            # 'legend_label': 'Ngp2_NgmIP',
            #},
            ### Graupel Mixing Ratio
            ## QG
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Graupel Mixing Ratio",
            # 'axis_title': r"qgm $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qgm',
            #},
            ## QG_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Graupel Mixing Ratio in Graupel",
            # 'axis_title': r"qgm_ip $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qgm_IP',
            #},
            ## QG'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance\nof Graupel Mixing Ratio",
            # 'axis_title': r"qgp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Qgp2',
            #},
            ## QG'^2_IP / QG_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-graupel Variance of Graupel Mixing Ratio",
            # 'axis_title': r"qgp2_ip / qgm_ip^2 [-]",
            # 'legend_label': 'Qgp2_QgmIP',
            #},
            ### Cloud Ice Mixing Ratio
            ## QI
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Ice Mixing Ratio",
            # 'axis_title': r"qim $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qim',
            #},
            ## QI_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Ice Mixing Ratio in Cloud Ice",
            # 'axis_title': r"qim_ip $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qim_IP',
            #},
            ## QI'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Cloud Ice Mixing Ratio",
            # 'axis_title': r"qip2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Qip2',
            #},
            ## QI'^2_IP / QI_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-cloud-ice Variance of Cloud Ice  Mixing Ratio",
            # 'axis_title': r"qip2_ip / qim_ip^2 [-]",
            # 'legend_label': 'Qip2_QimIP',
            #},
            ### Cloud Ice Number Concentration
            ## NI
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Ice Concentration",
            # 'axis_title': r"Nim $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nim',
            #},
            ## NI_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Cloud Ice Number Concentration in Cloud Ice",
            # 'axis_title': r"Ni_ip $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nim_IP',
            #},
            ## NI'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Cloud Ice Number Concentration",
            # 'axis_title': r"Nip2 $\mathrm{\left[num^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Nip2',
            #},
            ## NI'^2_IP / NI_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-cloud-ice Variance of Cloud Ice Number Concentration",
            # 'axis_title': r"Nip2_ip / Nim_ip^2 [-]",
            # 'legend_label': 'Nip2_NimIP',
            #},
            ## Snow Mixing Ratio
            ## QS
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Snow Mixing Ratio",
            # 'axis_title': r"qsm $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qsm',
            #},
            ## QS_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Snow Mixing Ratio in Snow",
            # 'axis_title': r"qsm_ip $\mathrm{\left[kg\,kg^{-1}\right]}$",
            # 'legend_label': 'Qsm_IP',
            #},
            ## QS'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Snow Mixing Ratio",
            # 'axis_title': r"qsp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Qsp2',
            #},
            ## QS'^2_IP / QS_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-snow Variance of Snow Mixing Ratio",
            # 'axis_title': r"qsp2_ip / qsm_ip^2 [-]",
            # 'legend_label': 'Qsp2_QsmIP',
            #},
            ## Snow Number Concentration
            ## NS
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Snow Number Concentration",
            # 'axis_title': r"Nsm $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nsm',
            #},
            ## NS_IP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Snow Number Concentration in Snow",
            # 'axis_title': r"Nsm_ip $\mathrm{\left[num\,kg^{-1}\right]}$",
            # 'legend_label': 'Nsm_IP',
            #},
            ## NS'^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Domain-wide Variance of Snow Number Concentration",
            # 'axis_title': r"Nsp2 $\mathrm{\left[\#^2\,kg^{-2}\right]}$",
            # 'legend_label': 'Nsp2',
            #},
            ## NS'^2_IP / NS_IP^2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': "Within-snow Variance of Snow Number Concentration",
            # 'axis_title': r"Nsp2_ip / Nsm_ip^2 [-]",
            # 'legend_label': 'Nsp2_NsmIP',
            #},
            ## MICROFRACTIONS
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': 'MicroFractions',
            # 'axis_title': "Micro Fractions", r"[%/100]",
            # 'legend_label': 'MicroFractions',
            #},
            ## W'THETA_V'
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': r"Buoyancy flux, $\overline{w'\theta_v'}$",
            # 'axis_title': r"wpthvp / tlflux $\mathrm{\left[K\,m\,s^{-1}\right]}$",
            # 'legend_label': r"$\overline{w'\theta_v'}$",
            #},
            ## LWP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1/1000,
            # 'title': r"CWP",
            # 'axis_title': r"CWP$\mathrm{\left[bla\right]}$",
            # 'legend_label': 'CWP',
            #},
            ## PREC
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': r"PREC",
            # 'axis_title': r"PREC $\mathrm{\left[bla\right]}$",
            # 'legend_label': r"PREC",
            #},
            ## WP2_W2
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1,
            # 'title': r"NS",
            # 'axis_title': r"NS $\mathrm{\left[bla\right]}$",
            # 'legend_label': r"NS",
            #},
            ## IWP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1/1000,
            # 'title': r"IWP",
            # 'axis_title': r"IWP $\mathrm{\left[bla\right]}$",
            # 'legend_label': r"IWP",
            #},
            ## SWP
            #{'var_names':
                #{
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
				#'cam': [],
                #},
            # 'sam_calc': , 'sam_conv_factor': 1/1000,
            # 'title': r"SWP",
            # 'axis_title': r"SWP $\mathrm{\left[bla\right]}$",
            # 'legend_label': r"SWP",
            #},
            # U'R_C'
            {'var_names':
                {
                'clubb': [],
                'sam': ['UPRCP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{u'r_c'}$",
             'axis_title': r"Liquid water flux, "+
                           r"$\overline{u'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
             'legend_label': r"$\overline{u'r_c'}$",
            },
            # U'R_T'
            {'var_names':
                {
                'clubb': [],
                'sam': ['UPRTP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{u'r_t'}$",
             'axis_title': r"Total water flux, "+
                           r"$\overline{u'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
             'legend_label': r"$\overline{u'r_t'}$",
            },
            # U'THETA_L'
            {'var_names':
                {
                'clubb': [],
                'sam': ['UPTHLP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{u'\theta_l'}$",
             'axis_title': r"Liq. water pot. temp. flux, "+
                           r"$\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{u'\theta_l'}$",
            },
            # U'THETA_V'
            {'var_names':
                {
                'clubb': [],
                'sam': ['UPTHVP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{u'\theta_v'}$",
             'axis_title': r"Virtual pot. temp. flux, "+
                           r"$\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{u'\theta_v'}$",
            },
            # V'R_C'
            {'var_names':
                {
                'clubb': [],
                'sam': ['VPRCP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{v'r_c'}$",
             'axis_title': r"Liquid water flux, "+
                           r"$\overline{v'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
             'legend_label': r"$\overline{v'r_c'}$",
            },
            # V'R_t'
            {'var_names':
                {
                'clubb': [],
                'sam': ['VPRTP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{v'r_t'}$",
             'axis_title': r"Total water flux, "+
                           r"$\overline{v'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
             'legend_label': r"$\overline{v'r_t'}$",
            },
            # V'THETA_L'
            {'var_names':
                {
                'clubb': [],
                'sam': ['VPTHLP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{v'\theta_l'}$",
             'axis_title': r"Liq. water pot. temp. flux, "+
                           r"$\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{v'\theta_l'}$",
            },
            # V'THETA_V'
            {'var_names':
                {
                'clubb': [],
                'sam': ['VPTHVP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': r"$\overline{v'\theta_v'}$",
             'axis_title': r"Virtual pot.temp. flux, "+
                           r"$\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
             'legend_label': r"$\overline{v'\theta_v'}$",
            },
            ## Tracer variables
            
            # TR01
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "Tracer 01",
             'axis_title': "Tracer 01 [TR]",
             'legend_label': 'TR01',
            },
            # TR01ADV
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01ADV'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "TR01 vert. adv. tend.",
             'axis_title': r"TR01 tendency due to vert. adv. $\mathrm{\left[TR/day\right]}$",
             'legend_label': 'TR01ADV',
            },
            # TR01DIFF
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01DIFF'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "TR01 vert. SGS transp. tend.",
             'axis_title': r"TR01 tendency due to vert. SGS transport $\mathrm{\left[TR/day\right]}$",
             'legend_label': 'TR01DIFF',
            },
            # TR01FLX
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01FLX'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "Total TR01 flux",
             'axis_title': r"Total flux of TR01  $\mathrm{\left[TR kg m^{-2} s^{-1}\right]}$",
             'legend_label': 'TR01FLX',
            },
            # TR01FLXS
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01FLXS'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "TR01 SGS flux",
             'axis_title': r"SGS flux of TR01  $\mathrm{\left[TR kg m^{-2} s^{-1}\right]}$",
             'legend_label': 'TR01FLXS',
            },
            # TR01PHYS
            {'var_names':
                {
                'clubb': [],
                'sam': ['TR01PHYS'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "TR01 phys. tend.",
             'axis_title': r"TR01 tendency due to physics $\mathrm{\left[TR/day\right]}$",
             'legend_label': 'TR01PHYS',
            },
            # RHO
            {'var_names':
                {
                'clubb': [],
                'sam': ['RHO'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1,
             'title': "Air density",
             'axis_title': r"Air density $\rho\ \mathrm{\left[kgm^{-3}\right]}$",
             'legend_label': r'$\rho$',
             },
        ]

        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, sam_benchmark_dataset=sam_benchmark_dataset,
                         coamps_benchmark_dataset=coamps_benchmark_dataset, wrf_benchmark_dataset=wrf_benchmark_dataset,
                         r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets= e3sm_datasets, wrf_datasets=wrf_datasets,
                         priority_vars=priority_vars, background_rcm=background_rcm,
                         background_rcm_folder=background_rcm_folder)
            
    def getThlmSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates thlm values from sam output using
        the following equation
        ``(THETAL + 2500.4.*(THETA./TABS).*(QI./1000))``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        thetal, indep, dataset = self.getVarForCalculations('THETAL', dataset)
        theta, indep, dataset = self.getVarForCalculations('THETA', dataset)
        tabs, indep, dataset = self.getVarForCalculations('TABS', dataset)
        qi, indep, dataset = self.getVarForCalculations('QI', dataset)
        
        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        return thlm, indep
    
    def getRtmSamCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates rtm values from sam output using
        the following equation
        ``(QT-QI) ./ 1000``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qt, indep, dataset = self.getVarForCalculations('QT', dataset)
        qi, indep, dataset = self.getVarForCalculations('QI', dataset)
        
        rtm = (qt - qi) / 1000
        return rtm, indep
    
    def getWpthlpCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        This gets called if WPTHLP isn't outputted in an nc file
        as a backup way of gathering the dependent_data for plotting.
        ``WPTHLP = (TLFLUX) ./ (RHO * 1004) + WPTHLP_SGS (CLUBB variable)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        tlflux, indep, dataset = self.getVarForCalculations(['TLFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        WPTHLP_SGS, indep, dataset = self.getVarForCalculations('WPTHLP_SGS', dataset)
        
        wpthlp = tlflux / (rho * 1004) + WPTHLP_SGS
        
        return wpthlp, indep
    
    def getCorrWpThlpCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the correlation of W and THETAL from SAM output
        using the following equation:
        ``(((TLFLUX) / (RHO * 1004.)) + WPTHLP_SGS)/np.sqrt(W2*TL2 + 1e-4)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TLFLUX, indep, dataset = self.getVarForCalculations('TLFLUX', dataset)
        RHO, indep, dataset = self.getVarForCalculations('RHO', dataset)
        WPTHLP_SGS, indep, dataset = self.getVarForCalculations('WPTHLP_SGS', dataset)
        W2, indep, dataset = self.getVarForCalculations('W2', dataset)
        TL2, indep, dataset = self.getVarForCalculations('TL2', dataset)
        
        CorrWpThlp = ( TLFLUX / (RHO * 1004) + WPTHLP_SGS ) / np.sqrt(W2 * TL2 + 1e-4)
        return CorrWpThlp, indep
    
    def getWprtpCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        This gets called if WPRTP isn't outputted in an nc file
        as a backup way of gathering the dependent_data for plotting.
        ``WPRTP = (QTFLUX) / (RHO * 2.5104e+6) + WPRTP_SGS (CLUBB variable)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qtflux, indep, dataset = self.getVarForCalculations(['QTFLUX'], dataset)
        rho, indep, dataset = self.getVarForCalculations(['RHO'], dataset)
        WPRTP_SGS, indep, dataset = self.getVarForCalculations('WPRTP_SGS', dataset)
        
        wprtp = qtflux / (rho * 2.5104e+6) + WPRTP_SGS
        return wprtp, indep
    
    def getCorrWpRtpCalc(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the correlation of W and THETAL from SAM output
        using the following equation:
        ``WPRTP/(np.sqrt(W2*QT2*1e-6)+1e-8)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WPRTP, indep, dataset = self.getVarForCalculations('WPRTP', dataset)
        WPRTP_SGS, indep, dataset = self.getVarForCalculations('WPRTP_SGS', dataset)
        W2, indep, dataset = self.getVarForCalculations('W2', dataset)
        QT2, indep, dataset = self.getVarForCalculations('QT2', dataset)
        
        CorrWpRtp = WPRTP / (np.sqrt(W2*QT2*1e-6)+1e-8)
        return CorrWpRtp, indep
    
    def getWp2Calc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total vertical momentum variance W'^2 from SAM output
        using the following equation:
        ``W2 = W2 + WP2_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W2, indep, dataset = self.getVarForCalculations('W2', dataset)
        WP2_SGS, indep, dataset = self.getVarForCalculations('WP2_SGS', dataset)
        WP2 = W2 + WP2_SGS
        return WP2, indep
    
    def getWp3Calc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total vertical momentum skewness W'^3 from SAM output
        using the following equation:
        ``W3 = W3 + WP3_SGS (CLUBB variable)``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        W3, indep, dataset = self.getVarForCalculations('W3', dataset)
        WP3_SGS, indep, dataset = self.getVarForCalculations('WP3_SGS', dataset)
        WP3 = W3 + WP3_SGS
        return WP3, indep
    
    def getThetalVarCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total variance of THETAL from SAM output
        using the following equation:
        ``THETALVAR = THL2 + THLP2_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TL2, indep, dataset = self.getVarForCalculations('TL2', dataset)
        THLP2_SGS, indep, dataset = self.getVarForCalculations('THLP2_SGS', dataset)
        THETALVAR = TL2 + THLP2_SGS
        return THETALVAR, indep
    
    def getRtVarCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total variance of RT from SAM output
        using the following equation:
        ``RTVAR = (QT2 * 1e-6) + RTP2_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QT2, indep, dataset = self.getVarForCalculations('QT2', dataset)
        RTP2_SGS, indep, dataset = self.getVarForCalculations('RTP2_SGS', dataset)
        RTVAR = (QT2 * 1e-6) + RTP2_SGS
        return RTVAR, indep
    
    def getUpWpCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total covariance of U and W from SAM output
        using the following equation:
        ``UPWP = UW + UPWP_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UW, indep, dataset = self.getVarForCalculations('UW', dataset)
        UPWP_SGS, indep, dataset = self.getVarForCalculations('UPWP_SGS', dataset)
        UPWP = UW + UPWP_SGS
        return UPWP, indep
    
    def getVpWpCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total covariance of V and W from SAM output
        using the following equation:
        ``VPWP = VW + VPWP_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VW, indep, dataset = self.getVarForCalculations('VW', dataset)
        VPWP_SGS, indep, dataset = self.getVarForCalculations('VPWP_SGS', dataset)
        VPWP = VW + VPWP_SGS
        return VPWP, indep
    
    def getUp2Calc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total variance of U from SAM output
        using the following equation:
        ``UVAR = U2 + UP2_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        U2, indep, dataset = self.getVarForCalculations('U2', dataset)
        UP2_SGS, indep, dataset = self.getVarForCalculations('UP2_SGS', dataset)
        UVAR = U2 + UP2_SGS
        return UVAR, indep
    
    def getVp2Calc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total variance of V from SAM output
        using the following equation:
        ``VVAR = V2 + VP2_SGS`` (CLUBB variable)

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        V2, indep, dataset = self.getVarForCalculations('V2', dataset)
        VP2_SGS, indep, dataset = self.getVarForCalculations('VP2_SGS', dataset)
        VVAR = V2 + VP2_SGS
        return VVAR, indep

    def getUpWpCorrCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total correlation of U and W from SAM output
        using the following equation:
        ``(UW+UPWP_SGS)/(np.sqrt((U2+UP2_SGS)*(W2+WP2_SGS)+1e-4))``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UW, indep, dataset = self.getVarForCalculations('UW', dataset)
        UPWP_SGS, indep, dataset = self.getVarForCalculations('UPWP_SGS', dataset)
        UPWP = UW + UPWP_SGS
        return UPWP, indep
    
    def getVpWpCorrCalc(self, dataset_override = None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the total correlation of V and W from SAM output
        using the following equation:
        ``(VW+VPWP_SGS)/(np.sqrt((V2+VP2_SGS)*(W2+WP2_SGS)+1e-4))``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VW, indep, dataset = self.getVarForCalculations('VW', dataset)
        VPWP_SGS, indep, dataset = self.getVarForCalculations('VPWP_SGS', dataset)
        VPWP = VW + VPWP_SGS
        return VPWP, indep
    
    def getQRP2_QRIP(self, dataset_override=None):
        """
        This is a "calculate function". Calculate functions are intended to be written by the user in the event that
        they need a variable that is not output by their atmospheric model. The general format for these functions
        is:
            1. Get the proper dataset. This is either passed in as dataset_override, or some benchmark dataset
            2. Get the equations needed variables from the dataset using ``self.getVarForCalculations()``
            3. Calculate the new variable
            4. (optional) If there are multiple valid equations, pick the one that worked using
               ``self.pickNonZeroOutput()``
            5. Return the data as (dependent,independent)

        For more information on calculate functions, see the "Creating a new calculated function (for calculated
        variables)" section of the README.md

        Calculates the Within-Rain Variance of Rain Water Mixing Ratio from SAM output
        using the following equation:
        ``(qrainp2_ip / (np.maximum(qrainm_ip, 1e-5)**2))``

        :param dataset_override: If passed, this netcdf dataset will be used to gather the data needed to calculate the
          given variable. if not passed, this function should attempt to find the best source for the data, e.g.
          the benchmark data for the given model
        :return: tuple of numeric lists of the form (dependent_data, independent_data) for the given variable being caluclated.
          Lists will be filled with NaN's if the variable could not be calculated.
        """
        dataset = self.sam_benchmark_dataset
        if dataset_override is not None:
            dataset = dataset_override
        qrainp2_ip, indep, dataset = self.getVarForCalculations('qrainp2_ip', dataset)
        qrainm_ip, indep, dataset = self.getVarForCalculations('qrainm_ip', dataset)
        QRP2_QRIP = qrainp2_ip / (np.maximum(qrainm_ip, 1e-5)**2)
        return QRP2_QRIP, indep
