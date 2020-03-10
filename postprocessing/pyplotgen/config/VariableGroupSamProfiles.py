'''
:author: Steffen Domke
:date: Mid 2019
TODO:   - Should invalid variables produce zero lines or none at all?
        - fallback <-> calc
        - Update axis labels at TASK
        - Implement/copy sam calcs
        - Figure out how to include standard plots in VariableGroupBase (labels etc.)
'''
import numpy as np

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupSamProfiles(VariableGroup):
    
    def __init__(self, case, clubb_datasets=None, sam_file=None, coamps_dataset=None, r408_dataset=None, hoc_dataset=None, e3sm_datasets=None):
        """
        
        :param clubb_datasets:
        :param case:
        :param sam_file:
        """
        self.name = "sam profile variables"
        
        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)
        
        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000
        
        u_cld_lines = [
            {'var_names': ['U'], 'legend_label': r"$\overline{u}$"},
            {'var_names': ['UCLD'], 'legend_label': r"$\overline{u}^\mathrm{{cld}}$"}
            ]
        
        v_cld_lines = [
            {'var_names': ['V'], 'legend_label': r"$\overline{v}$"},
            {'var_names': ['VCLD'], 'legend_label': r"$\overline{v}^\mathrm{{cld}}$"}
            ]
        
        w_cld_lines = [
            {'var_names': ['WM'], 'legend_label': r"$\overline{w}$"},
            {'var_names': ['WCLD'], 'legend_label': r"$\overline{w}^\mathrm{{cld}}$"}
            ]
        
        u_cond_lines = [
            {'var_names': ['U'], 'legend_label': r"$\overline{u}$"},
            {'var_names': ['UCLD'], 'legend_label': '$\overline{u}^\mathrm{cld}$'},
            {'var_names': ['UENV'], 'legend_label': r"$\overline{u}^\mathrm{{env}}$", 'fallback_func': self.getUEnvUnweighted}
            ]
        
        u_weight_lines = [
            {'var_names': ['U'], 'legend_label': r"$\overline{u}$"},
            {'var_names': ['UCLD'], 'legend_label': r"$\overline{u}^\mathrm{{cld}}$", 'fallback_func': self.getUCldWeighted},
            {'var_names': ['UENV'], 'legend_label': r"$\overline{u}^\mathrm{{env}}$", 'fallback_func': self.getUEnvWeighted}
            ]
        
        v_cond_lines = [
            {'var_names': ['V'], 'legend_label': r"$\overline{v}$"},
            {'var_names': ['VCLD'], 'legend_label': '$\overline{v}^\mathrm{cld}$'},
            {'var_names': ['VENV'], 'legend_label': r"$\overline{v}^\mathrm{{env}}$", 'fallback_func': self.getVEnvUnweighted}
            ]
        
        v_weight_lines = [
            {'var_names': ['V'], 'legend_label': r"$\overline{v}$"},
            {'var_names': ['VCLD'], 'legend_label': r"$\overline{v}^\mathrm{{cld}}$", 'fallback_func': self.getVCldWeighted},
            {'var_names': ['VENV'], 'legend_label': r"$\overline{v}^\mathrm{{env}}$", 'fallback_func': self.getWEnvWeighted}
            ]
        
        w_cond_lines = [
            {'var_names': ['WM'], 'legend_label': r"$\overline{w}$"},
            {'var_names': ['WCLD'], 'legend_label': '$\overline{w}^\mathrm{cld}$'},
            {'var_names': ['WENV'], 'legend_label': r"$\overline{w}^\mathrm{{env}}$", 'fallback_func': self.getWEnvUnweighted}
            ]
        
        w_weight_lines = [
            {'var_names': ['WM'], 'legend_label': r"$\overline{w}$"},
            {'var_names': ['WCLD'], 'legend_label': r"$\overline{w}^\mathrm{{cld}}$", 'fallback_func': self.getWCldWeighted},
            {'var_names': ['WENV'], 'legend_label': r"$\overline{w}^\mathrm{{env}}$", 'fallback_func': self.getWEnvWeighted}
            ]
        
        uw_cond_lines = [
            {'var_names': ['UW'], 'legend_label': r"$\overline{u'w'}$"},
            {'var_names': ['UWCLD'], 'legend_label': r"$\overline{u'w'}^\mathrm{cld}$"},
            {'var_names': ['UWENV'], 'legend_label': r"$\overline{u'w'}^\mathrm{env}$", 'fallback_func': self.getUWEnvUnweighted},
            ]
        
        uw_weight_lines = [
            {'var_names': ['UW'], 'legend_label': r"$\overline{u'w'}$"},
            {'var_names': ['UWCLD'], 'legend_label': r"$\overline{u'w'}^\mathrm{cld}$", 'fallback_func': self.getUWCldWeighted},
            {'var_names': ['UWENV'], 'legend_label': r"$\overline{u'w'}^\mathrm{env}$", 'fallback_func': self.getUWEnvWeighted},
            ]
        
        vw_cond_lines = [
            {'var_names': ['VW'], 'legend_label': r"$\overline{v'w'}$"},
            {'var_names': ['VWCLD'], 'legend_label': r"$\overline{v'w'}^\mathrm{cld}$"},
            {'var_names': ['VWENV'], 'legend_label': r"$\overline{v'w'}^\mathrm{env}$", 'fallback_func': self.getVWEnvUnweighted},
            ]
        
        vw_weight_lines = [
            {'var_names': ['VW'], 'legend_label': r"$\overline{v'w'}$"},
            {'var_names': ['VWCLD'], 'legend_label': r"$\overline{v'w'}^\mathrm{cld}$", 'fallback_func': self.getVWCldWeighted},
            {'var_names': ['VWENV'], 'legend_label': r"$\overline{v'w'}^\mathrm{env}$", 'fallback_func': self.getVWEnvWeighted},
            ]
        
        thv_cond_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$"},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$", 'fallback_func': self.getTHVEnvUnweighted},
            ]
        
        thv_weight_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$", 'fallback_func': self.getTHVCldWeighted},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$", 'fallback_func': self.getTHVEnvWeighted},
            ]
        
        thl_weight_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$"},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$", 'fallback_func': self.getTHVEnvUnweighted},
            ]
        
        qtw_weight_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$", 'fallback_func': self.getTHVCldWeighted},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$", 'fallback_func': self.getTHVEnvWeighted},
            ]
        
        qt_cond_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$", 'fallback_func': self.getTHVCldWeighted},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$", 'fallback_func': self.getTHVEnvWeighted},
            ]
        
        uw_lines = [
            {'var_names': ['UW'], 'legend_label': r"$\overline{u'w'}$"},
            {'var_names': ['U2'], 'legend_label': r"$\overline{u'^2}$"},
            {'var_names': ['W2'], 'legend_label': r"$\overline{w'^2}$"},
            ]
        
        vw_lines = [
            {'var_names': ['VW'], 'legend_label': r"$\overline{v'w'}$"},
            {'var_names': ['V2'], 'legend_label': r"$\overline{v'^2}$"},
            {'var_names': ['W2'], 'legend_label': r"$\overline{w'^2}$"},
            ]
        
        self.variable_definitions = [
            # THETA_L
            {'var_names': {
                'clubb': ['thlm'],
                'sam': ['THETAL'],
                'coamps': ['thlm'],
                'r408': ['thlm'],
                'hoc': ['thlm'],
                'e3sm': ['thlm'],
                },
            'sam_calc': self.getThlmSamCalc, 'sam_conv_factor': 1, 'title': r"Liquid Water Potential Temperature, $\mathrm{\theta_l}$", 'axis_title': r"$\mathrm{\theta_l}$ $\mathrm{\left[K\right]}$", 'legend_label': r"$\mathrm{\theta_l}$",
            },
            
            # R_T
            {'var_names': {
                'clubb': ['rtm'],
                'sam': ['QT'],
                'coamps': ['qtm'],
                'r408': ['rtm'],
                'hoc': ['rtm'],
                'e3sm': ['rtm'],
                },
            'sam_calc': self.getRtmSamCalc, 'sam_conv_factor': 1, 'title': r"Total Water Mixing Ratio, $\mathrm{r_t}$", 'axis_title': r"$\mathrm{r_t}$ $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': r'$\mathrm{r_t}$',
            },
            
            # W'THETA_L'
            {'var_names': {
                'clubb': ['wpthlp'],
                'sam': ['TLFLUX'],
                'coamps': ['wpthlp'],
                'r408': ['wpthlp'],
                'hoc': ['wpthlp'],
                'e3sm': ['wpthlp']
                },
            'sam_calc': self.getWpthlpCalc, 'sam_conv_factor': 1, 'title': r"Turbulent Flux of $\mathrm{\theta_l}$", 'axis_title': r"$\mathrm{\overline{w'\theta_l'}}$ / thflux(s) $\mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\mathrm{\overline{w'\theta_l'}}$",
            },
            
            # CORR(W, THETA_L)
            {'var_names': {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getCorrWpThlpCalc, 'sam_conv_factor': 1, 'title': r"Corr(w,$\mathrm{\theta_l}$)", 'axis_title': r"Correlation $\mathrm{\overline{w'\theta_l'} / \sqrt{\overline{w'^2}\;\overline{\theta_l'^2}}}$ $\mathrm{\left[-\right]}$", 'legend_label': r"$\mathrm{Corr(w,\theta_l)}$",
            },
            
            # W'R_T'
            {'var_names': {
                'clubb': ['wprtp'],
                'sam': ['QTFLUX'],
                'coamps': ['wpqtp'],
                'r408': ['wprtp'],
                'hoc': ['wprtp'],
                'e3sm': ['wprtp']
                },
            'sam_calc': self.getWprtpCalc, 'sam_conv_factor': 1, 'title': r"Turbulent Flux of $\mathrm{r_t}$", 'axis_title': r"$\mathrm{\overline{w'r_t'}}$ / qtflux(s) $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$", 'legend_label': r"$\mathrm{\overline{w'r_t'}}$",
            },
            
            # CORR(W, R_T)
            {'var_names': {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getCorrWpRtpCalc, 'sam_conv_factor': 1, 'title': r"$\mathrm{Corr(w,r_t)}$", 'axis_title': r"Correlation, $\mathrm{\overline{w'r_t'} / \sqrt{\overline{w'^2}\;\overline{r_t'^2}}}$ $\mathrm{\left[-\right]}$", 'legend_label': r"$\mathrm{Corr(w,r_t)}$",
            },
            
            # cloudliq_frac_em6
            {'var_names': {
                'clubb': [],
                'sam': ['cloudliq_frac_em6'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "Cloud Liquid Fraction", 'axis_title': r"Cloud liquid fraction $\mathrm{\left[\frac{\%}{100}\right]}$", 'legend_label': 'cloudliq_frac_em6',
            },
            
            # R_C
            {'var_names': {
                'clubb': [],
                'sam': ['QCL'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1e-3, 'title': r"Cloud Water Mixing Ratio, $\mathrm{r_c}$", 'axis_title': r"$\mathrm{r_c}$ / $\mathrm{q_{cl}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': r'$\mathrm{r_c}$',
            },
            
            # W'^2
            {'var_names': {
                'clubb': [],
                'sam': ['W2'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getWp2Calc, 'sam_conv_factor': 1, 'title': r"Vertical Momentum Variance, $\mathrm{\overline{w'^2}}$", 'axis_title': r"Momentum variance, $\mathrm{\overline{w'^2}}$ $\mathrm{\left[m^2\,s^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{w'^2}}$",
            },
            
            # W'^3
            {'var_names': {
                'clubb': [],
                'sam': ['W3'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getWp3Calc, 'sam_conv_factor': 1, 'title': r"Vertical Momentum Skewness, $\mathrm{\overline{w'^3}}$", 'axis_title': r"Momentum Skewness, $\mathrm{\overline{w'^3}}$ $\mathrm{\left[m^3\,s^{-3}\right]}$", 'legend_label': r"$\mathrm{\overline{w'^3}}$",
            },
            
            #THETA_L'^2
            {'var_names': {
                'clubb': [],
                'sam': ['TL2'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getThetalVarCalc, 'sam_conv_factor': 1, 'title': r"Variance of Liquid Water Potential Temperature, $\mathrm{\overline{\theta_l'^2}}$", 'axis_title': r"$\mathrm{\overline{\theta_l'^2}}$ $\mathrm{\left[K^2\right]}$", 'legend_label': r"$\mathrm{\overline{\theta_l'^2}}$",
            },
            
            # R_T'^2
            {'var_names': {
                'clubb': [],
                'sam': ['QT2'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getRtVarCalc, 'sam_conv_factor': 1, 'title': r"Variance of Total Water Mixing Ratio, $\mathrm{\overline{r_t'^2}}$", 'axis_title': r"$\mathrm{\overline{r_t'^2}}$ / $\mathrm{\overline{q_t'^2}}$ $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{r_t'^2}}$",
            },
            
            # R_T'THETA_L'
            {'var_names': {
                'clubb': [],
                'sam': ['RTPTHLP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"Covariance of $\mathrm{r_t}$ & $\mathrm{\theta_l}$", 'axis_title': r"$\mathrm{\overline{r_t'\theta_l'}}$ $\mathrm{\left[kg\,kg^{-1}\,K\right]}$", 'legend_label': r"$\mathrm{\overline{r_t'\theta_l'}}$",
            },
            
            # W_OBS
            {'var_names': {
                'clubb': [],
                'sam': ['WOBS'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\mathrm{w_{obs}}$", 'axis_title': r"Observed wind, $\mathrm{w_{obs}}\ \mathrm{\left[m\,s^{-1}\right]}$", 'legend_label': r"$\mathrm{w_{obs}}$",
            },
            
            # U
            {'var_names': {
                'clubb': [],
                'sam': ['U'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"Zonal Wind Component, $\mathrm{\overline{u}}$", 'axis_title': r"$\mathrm{\overline{u}}\ \mathrm{\left[m\,s^{-1}\right]}$", 'legend_label': r'\mathrm{\overline{u}}',
            },
            
            # V
            {'var_names': {
                'clubb': [],
                'sam': ['V'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"Meridonal Wind Component, $\mathrm{\overline{v}}$", 'axis_title': r"$\mathrm{\overline{v}}\ \mathrm{\left[m\,s^{-1}\right]}$", 'legend_label': r"$\mathrm{\overline{v}}$",
            },
            
            # U'W'
            {'var_names': {
                'clubb': [],
                'sam': ['UW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getUpWpCalc, 'sam_conv_factor': 1, 'title': r"$\mathrm{\overline{u'w'}}$", 'axis_title': r"Momentum flux, $\mathrm{\overline{u'w'}\ \left[m^2\,s^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{u'w'}}$",
            },
            
            # V'W'
            {'var_names': {
                'clubb': [],
                'sam': ['VW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getVpWpCalc, 'sam_conv_factor': 1, 'title': r"$\mathrm{\overline{v'w'}}$", 'axis_title': r"Momentum flux, $\mathrm{\overline{v'w'}\ \left[m^2\,s^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{v'w'}}$",
            },
            
            # U'^2
            {'var_names': {
                'clubb': [],
                'sam': ['U2'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getUp2Calc, 'sam_conv_factor': 1, 'title': r"$\mathrm{\overline{u'^2}}$", 'axis_title': r"Momentum variance, $\mathrm{\overline{u'^2}\ \left[m^2\,s^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{u'^2}}$",
            },
            
            # V'^2
            {'var_names': {
                'clubb': [],
                'sam': ['V2'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getVp2Calc, 'sam_conv_factor': 1, 'title': r"$\mathrm{\overline{v'^2}}$", 'axis_title': r"Momentum variance, $\mathrm{\overline{v'^2}\ \left[m^2\,s^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{v'^2}}$",
            },
            
            # CORR(U, W)
            {'var_names': {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getUpWpCorrCalc, 'sam_conv_factor': 1, 'title': r"Corr(u,w)", 'axis_title': r"Correlation, $\mathrm{\overline{u'w'} / \sqrt{\overline{u'^2}\;\overline{w'^2}}\ \left[-\right]}$", 'legend_label': r"Corr(u,w)",
            },
            
            # CORR(V, W)
            {'var_names': {
                'clubb': [],
                'sam': [],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': self.getVpWpCorrCalc, 'sam_conv_factor': 1, 'title': r"Corr(v,w)", 'axis_title': r"Correlation, $\mathrm{\overline{v'w'} / \sqrt{\overline{v'^2}\;\overline{w'^2}}\ \left[-\right]}$", 'legend_label': r"Corr(v,w)",
            },
            
            ### PRECIPITATION
            
            ## Rain Water Mixing Ratio (only MICRO_DRIZZLE, MICRO_2MOM, and MICRO_SAM1MOM)
            # QR
            #{'var_names': {
                #'clubb': [],
                #'sam': ['QR'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_conv_factor': 1/1000, 'title': "Rain Water Mixing Ratio", 'axis_title': r"qrm $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': r"QR",
            #},
            ##'sam_conv_factor': 1/1000, 'title': "Rain Water Mixing Ratio", 'axis_title': r"$\mathrm{\overline{q_r}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': r"$\mathrm{\overline{q_r}}$",
            ##},
            ## QR_IP (only MICRO_DRIZZLE, MICRO_M2005_UWM)
            #{'var_names': {
                #'clubb': [],
                #'sam': ['qrainm_ip'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_conv_factor': 1, 'title': "Rain Water Mixing Ratio in Rain", 'axis_title': r"qrm_ip $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'QR_IP',
            ##'sam_conv_factor': 1, 'title': "Rain Water Mixing Ratio in Rain", 'axis_title': r"$\mathrm{q_r^{ip}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': r"$\mathrm{q_r^{ip}}$",
            ##},
            #},
            ## QR'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': ['qrainp2'],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_conv_factor': 1, 'title': "Domain-wide Variance of Rain Water Mixing Ratio", 'axis_title': r"qrp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': 'qrainp2',
            #},
            ##'sam_conv_factor': 1, 'title': "Domain-wide Variance of Rain Water Mixing Ratio", 'axis_title': r"$\mathrm{\overline{q_r'^2}}$ $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': r"$\mathrm{\overline{q_r'^2}}$",
            ##},
            ## QR'^2_IP / QR_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': self.getQRP2_QRIP, 'sam_conv_factor': 1, 'title': "Within-rain Variance of Rain Water Mixing Ratio", 'axis_title': r"qrp2_ip / qrm_ip^2 [-]", 'legend_label': 'QRP2_QRIP',
            #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-rain Variance of Rain Water Mixing Ratio", 'axis_title': r"$\mathrm{\overline{q_r'^2}^{ip} / \overline{q_r^{ip}}^2}$ [-]", 'legend_label': 'QRP2_QRIP',
            #},
            
            ### Rain Drop Number Concentration
            ## NR
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Rain Drop Concentration", 'axis_title': r"Nrm $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nrm',
            #},
            ## NR_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Rain Drop Concentration in Rain", 'axis_title': r"Nrm_ip $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nrm_IP',
            #},
            ## NR'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance\nof Rain Drop Concentration", 'axis_title': r"Nrp2 $\mathrm{\left[num^2\,kg^{-2}\right]}$", 'legend_label': 'Nrp2',
            #},
            ## NR'^2_IP / NR_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-rain Variance of Rain Drop Concentration", 'axis_title': r"Nrp2_ip / Nrm_ip^2 [-]", 'legend_label': 'Nrp2_NrmIP',
            #},
            
            ### Cloud Droplet Number Concentration
            ## NC
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Droplet Number Concentration", 'axis_title': r"Ncm $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Ncm',
            #},
            ## NC_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Droplet Number Concentration in Cloud", 'axis_title': r"Ncm_ip $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Ncm_IP',
            #},
            ## NC'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Cloud Droplet Number Concentration", 'axis_title': r"Ncp2 $\mathrm{\left[\#^2\,kg^{-2}\right]}$", 'legend_label': 'Ncp2',
            #},
            ## NC'^2_IP / NC_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-cloud Variance of Cloud Droplet Number Concentration", 'axis_title': r"Ncp2_ip / Ncm_ip^2 [-]", 'legend_label': 'Ncp2_NcmIP',
            #},
            
            ### Graupel Number Concentration
            ## NG
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Graupel Number Concentration", 'axis_title': r"Ngm $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Ngm',
            #},
            ## NG_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Graupel Number Concentration in Graupel", 'axis_title': r"Ngm_ip $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Ngm_IP',
            #},
            ## NG'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Graupel Number Concentration", 'axis_title': r"Ngp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': 'Ngp2',
            #},
            ## NG'^2_IP / NG_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-graupel Variance\nof Graupel Number Concentration", 'axis_title': r"Ngp2_ip / Ngm_ip^2 [-]", 'legend_label': 'Ngp2_NgmIP',
            #},
            
            ### Graupel Mixing Ratio
            ## QG
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Graupel Mixing Ratio", 'axis_title': r"qgm $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qgm',
            #},
            ## QG_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Graupel Mixing Ratio in Graupel", 'axis_title': r"qgm_ip $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qgm_IP',
            #},
            ## QG'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance\nof Graupel Mixing Ratio", 'axis_title': r"qgp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': 'Qgp2',
            #},
            ## QG'^2_IP / QG_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-graupel Variance of Graupel Mixing Ratio", 'axis_title': r"qgp2_ip / qgm_ip^2 [-]", 'legend_label': 'Qgp2_QgmIP',
            #},
            
            ### Cloud Ice Mixing Ratio
            ## QI
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Ice Mixing Ratio", 'axis_title': r"qim $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qim',
            #},
            ## QI_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Ice Mixing Ratio in Cloud Ice", 'axis_title': r"qim_ip $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qim_IP',
            #},
            ## QI'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Cloud Ice Mixing Ratio", 'axis_title': r"qip2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': 'Qip2',
            #},
            ## QI'^2_IP / QI_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-cloud-ice Variance of Cloud Ice  Mixing Ratio", 'axis_title': r"qip2_ip / qim_ip^2 [-]", 'legend_label': 'Qip2_QimIP',
            #},
            
            ### Cloud Ice Number Concentration
            ## NI
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Ice Concentration", 'axis_title': r"Nim $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nim',
            #},
            ## NI_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Cloud Ice Number Concentration in Cloud Ice", 'axis_title': r"Ni_ip $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nim_IP',
            #},
            ## NI'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Cloud Ice Number Concentration", 'axis_title': r"Nip2 $\mathrm{\left[num^2\,kg^{-2}\right]}$", 'legend_label': 'Nip2',
            #},
            ## NI'^2_IP / NI_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-cloud-ice Variance of Cloud Ice Number Concentration", 'axis_title': r"Nip2_ip / Nim_ip^2 [-]", 'legend_label': 'Nip2_NimIP',
            #},
            
            ## Snow Mixing Ratio
            ## QS
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Snow Mixing Ratio", 'axis_title': r"qsm $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qsm',
            #},
            ## QS_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Snow Mixing Ratio in Snow", 'axis_title': r"qsm_ip $\mathrm{\left[kg\,kg^{-1}\right]}$", 'legend_label': 'Qsm_IP',
            #},
            ## QS'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Snow Mixing Ratio", 'axis_title': r"qsp2 $\mathrm{\left[kg^2\,kg^{-2}\right]}$", 'legend_label': 'Qsp2',
            #},
            ## QS'^2_IP / QS_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-snow Variance of Snow Mixing Ratio", 'axis_title': r"qsp2_ip / qsm_ip^2 [-]", 'legend_label': 'Qsp2_QsmIP',
            #},
            
            ## Snow Number Concentration
            ## NS
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Snow Number Concentration", 'axis_title': r"Nsm $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nsm',
            #},
            ## NS_IP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Snow Number Concentration in Snow", 'axis_title': r"Nsm_ip $\mathrm{\left[num\,kg^{-1}\right]}$", 'legend_label': 'Nsm_IP',
            #},
            ## NS'^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Domain-wide Variance of Snow Number Concentration", 'axis_title': r"Nsp2 $\mathrm{\left[\#^2\,kg^{-2}\right]}$", 'legend_label': 'Nsp2',
            #},
            ## NS'^2_IP / NS_IP^2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': "Within-snow Variance of Snow Number Concentration", 'axis_title': r"Nsp2_ip / Nsm_ip^2 [-]", 'legend_label': 'Nsp2_NsmIP',
            #},
            ## MICROFRACTIONS
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': 'MicroFractions', 'axis_title': "Micro Fractions", r"[%/100]", 'legend_label': 'MicroFractions',
            #},
            ## W'THETA_V'
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': r"Buoyancy flux, $\overline{w'\theta_v'}$", 'axis_title': r"wpthvp / tlflux $\mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\overline{w'\theta_v'}$",
            #},
            
            ## LWP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1/1000, 'title': r"CWP", 'axis_title': r"CWP$\mathrm{\left[bla\right]}$", 'legend_label': 'CWP',
            #},
            
            ## PREC
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': r"PREC", 'axis_title': r"PREC $\mathrm{\left[bla\right]}$", 'legend_label': r"PREC",
            #},
            
            ## WP2_W2
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1, 'title': r"NS", 'axis_title': r"NS $\mathrm{\left[bla\right]}$", 'legend_label': r"NS",
            #},
            
            ## IWP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1/1000, 'title': r"IWP", 'axis_title': r"IWP $\mathrm{\left[bla\right]}$", 'legend_label': r"IWP",
            #},
            
            ## SWP
            #{'var_names': {
                #'clubb': [],
                #'sam': [],
                #'coamps': [],
                #'r408': [],
                #'hoc': [],
                #'e3sm': [],
                #},
            #'sam_calc': , 'sam_conv_factor': 1/1000, 'title': r"SWP", 'axis_title': r"SWP $\mathrm{\left[bla\right]}$", 'legend_label': r"SWP",
            #},
            
            # U'R_C'
            {'var_names': {
                'clubb': [],
                'sam': ['URCP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{u'r_c'}$", 'axis_title': r"Liquid water flux, $\overline{u'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$", 'legend_label': r"$\overline{u'r_c'}$",
            },
            
            # U'R_T'
            {'var_names': {
                'clubb': [],
                'sam': ['UPRTP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{u'r_t'}$", 'axis_title': r"Total water flux, $\overline{u'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$", 'legend_label': r"$\overline{u'r_t'}$",
            }, 
            
            # U'THETA_L'
            {'var_names': {
                'clubb': [],
                'sam': ['UPTHLP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{u'\theta_l'}$", 'axis_title': r"Liq. water pot. temp. flux, $\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\overline{u'\theta_l'}$",
            },
            
            # U'THETA_V'
            {'var_names': {
                'clubb': [],
                'sam': ['UPTHVP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{u'\theta_v'}$", 'axis_title': r"Virtual pot. temp. flux, $\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\overline{u'\theta_v'}$",
            },
            
            # V'R_C'
            {'var_names': {
                'clubb': [],
                'sam': ['VPRCP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{v'r_c'}$", 'axis_title': r"Liquid water flux, $\overline{v'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$", 'legend_label': r"$\overline{v'r_c'}$",
            },
            
            # V'R_t'
            {'var_names': {
                'clubb': [],
                'sam': ['VPRTP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{v'r_t'}$", 'axis_title': r"Total water flux, $\overline{v'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$", 'legend_label': r"$\overline{v'r_t'}$",
            },
            
            # V'THETA_L'
            {'var_names': {
                'clubb': [],
                'sam': ['VPTHLP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{v'\theta_l'}$", 'axis_title': r"Liq. water pot. temp. flux, $\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\overline{v'\theta_l'}$",
            },
            
            # V'THETA_V'
            {'var_names': {
                'clubb': [],
                'sam': ['VPTHVP'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': r"$\overline{v'\theta_v'}$", 'axis_title': r"Virtual pot.temp. flux, $\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$", 'legend_label': r"$\overline{v'\theta_v'}$",
            },
            
            ## Cloud conditional plots
            
            # UCLD
            {'var_names': {
                'clubb': [],
                'sam': ['UCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': u_cld_lines, 'title': r"Conditional mean wind, $\overline{u}^\mathrm{{cld}}$", 'axis_title': r"Conditional mean wind, $\overline{u}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
            },
            
            # VCLD
            {'var_names': {
                'clubb': [],
                'sam': ['VCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': v_cld_lines, 'title': r"Conditional mean wind, $\overline{v}^\mathrm{{cld}}$", 'axis_title': r"Conditional mean wind, $\overline{v}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
            },
            
            # WCLD
            {'var_names': {
                'clubb': [],
                'sam': ['WCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': w_cld_lines, 'title': r"Conditional mean wind, $\overline{w}^\mathrm{{cld}}$", 'axis_title': r"Conditional mean wind, $\overline{w}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
            },
            
            # UCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['UCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': u_cond_lines, 'title': r"Conditional mean wind, $\overline{u}$", 'axis_title': r"Conditional mean wind, $\overline{w}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
            },
            
            # UCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['UCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': u_weight_lines, 'title': r"Weighted mean wind, $\overline{u}$", 'axis_title': r"Weighted mean wind, $\overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
            },
            
            # VCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['VCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': v_cond_lines, 'title': r"Conditional mean wind, $\overline{v}$", 'axis_title': r"Conditional mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
            },
            
            # VCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['VCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': v_weight_lines, 'title': r"Weighted mean wind, $\overline{v}$", 'axis_title': r"Weighted mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
            },
            
            # WCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['WCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': w_cond_lines, 'title': r"Conditional mean wind, $\overline{w}$", 'axis_title': r"Conditional mean wind, $\overline{w}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
            },
            
            # WCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['WCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': w_weight_lines, 'title': r"Weighted mean wind, $\overline{w}$", 'axis_title': r"Weighted mean wind, $\overline{w}\ \mathrm{\left[m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
            },
            
            # UWCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['UWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_cond_lines, 'title': r"Cloud-conditonal $\overline{u'w'}$", 'axis_title': r"Conditional flux, $\overline{u'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            #'legend_label': r"$\overline{u'w'}^\mathrm{cld}$",
            },
            
            # UWCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['UWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_weight_lines, 'title': r"Cloud-weighted $\overline{u'w'}$", 'axis_title': r"Weighted flux, $\overline{u'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            #'legend_label': r"$\overline{u'w'}^\mathrm{cld}$",
            },
            
            # VWCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['VWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_cond_lines, 'title': r"Cloud-conditonal $\overline{v'w'}$", 'axis_title': r"Conditional flux, $\overline{v'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            #'legend_label': r"$\overline{v'w'}^\mathrm{cld}$",
            },
            
            # VWCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['VWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_weight_lines, 'title': r"Cloud-weighted $\overline{v'w'}$", 'axis_title': r"Weighted flux, $\overline{v'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            #'legend_label': r"$\overline{v'w'}^\mathrm{cld}$",
            },
            
            # TVCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['TVCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': thv_cond_lines, 'title': r"Conditional virt. pot. temp., $\overline{\theta}_v$", 'axis_title': r"Conditional virt. pot. temp., $\overline{\theta}_v\ \mathrm{\left[K\right]}$",
            #'legend_label': r"$\overline{\theta}_v^\mathrm{cld}",
            },
            
            # TVCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['TVCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': thv_weight_lines, 'title': r"Weighted virt. pot. temp., $\overline{\theta}_v$", 'axis_title': r"Weighted virt. pot. temp., $\overline{\theta}_v\ \mathrm{\left[K\right]}$",
            #'legend_label': r"$\overline{\theta}_v^\mathrm{cld}",
            },
            
            # TLCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['TLCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': thl_weight_lines, 'title': r"Weighted flux, $\overline{w's_L'}$", 'axis_title': r"Weighted flux, $\overline{w's_L'}\ \mathrm{\left[K\,m\, s^{-1}\right]}$",
            #'legend_label': r"$\overline{w's_L'}^\mathrm{cld}$",
            },
            
            # QTWCLD weighted
            {'var_names': {
                'clubb': [],
                'sam': ['QTWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': qtw_weight_lines, 'title': r"Weighted flux, $\overline{w'r_t'}$", 'axis_title': r"Weighted flux, $\overline{w'r_t'}\ \mathrm{\left[kg\,kg^{-1}\, m\,s^{-1}\right]}$",
            #'legend_label': r"$\overline{w'r_t'}^\mathrm{cld}$",
            },
            
            # QTCLD unweighted
            {'var_names': {
                'clubb': [],
                'sam': ['QTCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': qt_cond_lines, 'title': r"Conditional total water mixing ratio, $\mathrm{r_t}$", 'axis_title': r"Conditional total water mixing ratio, $\mathrm{r_t}\ \mathrm{\left[g\,kg^{-1}\right]}$",
            #'legend_label': r"$\mathrm{r_t^{cld}}$",
            },
            
            # U 2ND MOMENTS
            {'var_names': {
                'clubb': [],
                'sam': ['UW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_lines, 'title': "Eastward 2nd-moments", 'axis_title': r"2nd moments $\mathrm{\left[m^2\,s^{-2}\right]}$",
            #'legend_label': '',
            },
            
            # V 2ND MOMENTS
            {'var_names': {
                'clubb': [],
                'sam': ['VW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_lines, 'title': "Northward 2nd-moments", 'axis_title': r"2nd moments $\mathrm{\left[m^2\,s^{-2}\right]}$",
            #'legend_label': '',
            },
            
            # CLD
            {'var_names': {
                'clubb': [],
                'sam': ['CLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "Cloud fraction", 'axis_title': "Cloud fraction [-]", 'legend_label': 'CLD',
            },
            
            ## Tracer variables
            
            # TR01
            {'var_names': {
                'clubb': [],
                'sam': ['TR01'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "Tracer 01", 'axis_title': "Tracer 01 [TR]", 'legend_label': 'TR01',
            },
            
            # TR01ADV
            {'var_names': {
                'clubb': [],
                'sam': ['TR01ADV'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "TR01 vert. adv. tend.", 'axis_title': r"TR01 tendency due to vert. adv. $\mathrm{\left[TR/day\right]}$", 'legend_label': 'TR01ADV',
            },
            
            # TR01DIFF
            {'var_names': {
                'clubb': [],
                'sam': ['TR01DIFF'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "TR01 vert. SGS transp. tend.", 'axis_title': r"TR01 tendency due to vert. SGS transport $\mathrm{\left[TR/day\right]}$", 'legend_label': 'TR01DIFF',
            },
            
            # TR01FLX
            {'var_names': {
                'clubb': [],
                'sam': ['TR01FLX'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "Total TR01 flux", 'axis_title': r"Total flux of TR01  $\mathrm{\left[TR kg m^{-2} s^{-1}\right]}$", 'legend_label': 'TR01FLX',
            },
            
            # TR01FLXS
            {'var_names': {
                'clubb': [],
                'sam': ['TR01FLXS'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "TR01 SGS flux", 'axis_title': r"SGS flux of TR01  $\mathrm{\left[TR kg m^{-2} s^{-1}\right]}$", 'legend_label': 'TR01FLXS',
            },
            
            # TR01PHYS
            {'var_names': {
                'clubb': [],
                'sam': ['TR01PHYS'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_calc': , 'sam_conv_factor': 1, 'title': "TR01 phys. tend.", 'axis_title': r"TR01 tendency due to physics $\mathrm{\left[TR/day\right]}$", 'legend_label': 'TR01PHYS',
            },
            
            # RHO
            {'var_names': {
                'clubb': [],
                'sam': ['RHO'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
                },
            'sam_conv_factor': 1, 'title': "Air density", 'axis_title': r"Air density $\rho\ \mathrm{\left[kgm^{-3}\right]}$", 'legend_label': r'$\rho$',
            },
            ]
        
        # Call ctor of parent class
        super().__init__(ncdf_datasets, case, sam_file=sam_file, coamps_file=coamps_file, r408_dataset=r408_dataset, hoc_dataset=hoc_dataset, e3sm_datasets= e3sm_datasets)
            
    def getThlmSamCalc(self, dataset_override=None):
        """
        Calculates thlm values from sam output using
        the following equation
        (THETAL + 2500.4.*(THETA./TABS).*(QI./1000))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        thetal, z, dataset = self.getVarForCalculations('THETAL', self.sam_file)
        theta, z, dataset = self.getVarForCalculations('THETA', self.sam_file)
        tabs, z, dataset = self.getVarForCalculations('TABS', self.sam_file)
        qi, z, dataset = self.getVarForCalculations('QI', self.sam_file)
        
        thlm = thetal + (2500.4 * (theta / tabs) * (qi / 1000))
        return thlm, z
    
    def getRtmSamCalc(self, dataset_override=None):
        """
        Calculates rtm values from sam output using
        the following equation
        (QT-QI) ./ 1000
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        qt, z, dataset = self.getVarForCalculations('QT', self.sam_file)
        qi, z, dataset = self.getVarForCalculations('QI', self.sam_file)
        
        rtm = (qt - qi) / 1000
        return rtm, z
    
    def getWpthlpCalc(self, dataset_override=None):
        """
        This gets called if WPTHLP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPTHLP = (TLFLUX) ./ (RHO * 1004) + WPTHLP_SGS (CLUBB variable)
        :return:
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        tlflux, z, dataset = self.getVarForCalculations(['TLFLUX'], self.sam_file)
        rho, z, dataset = self.getVarForCalculations(['RHO'], self.sam_file)
        WPTHLP_SGS, z, dataset = self.getVarForCalculations('WPTHLP_SGS', self.sam_file)
        
        wpthlp = tlflux / (rho * 1004) + WPTHLP_SGS
        
        return wpthlp, z
    
    def getCorrWpThlpCalc(self, dataset_override=None):
        """
        Calculates the correlation of W and THETAL from SAM output
        using the following equation:
        (((TLFLUX) / (RHO * 1004.)) + WPTHLP_SGS)/np.sqrt(W2*TL2 + 1e-4)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        TLFLUX, z, dataset = self.getVarForCalculations('TLFLUX', self.sam_file)
        RHO, z, dataset = self.getVarForCalculations('RHO', self.sam_file)
        WPTHLP_SGS, z, dataset = self.getVarForCalculations('WPTHLP_SGS', self.sam_file)
        W2, z, dataset = self.getVarForCalculations('W2', self.sam_file)
        TL2, z, dataset = self.getVarForCalculations('TL2', self.sam_file)
        
        CorrWpThlp = ( TLFLUX / (RHO * 1004) + WPTHLP_SGS ) / np.sqrt(W2 * TL2 + 1e-4)
        return CorrWpThlp, z
    
    def getWprtpCalc(self, dataset_override=None):
        """
        This gets called if WPRTP isn't outputted in an nc file as a backup way of gathering the data for plotting.
        WPRTP = (QTFLUX) / (RHO * 2.5104e+6) + WPRTP_SGS (CLUBB variable)
        :return:
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        qtflux, z, dataset = self.getVarForCalculations(['QTFLUX'], self.sam_file)
        rho, z, dataset = self.getVarForCalculations(['RHO'], self.sam_file)
        WPRTP_SGS, z, dataset = self.getVarForCalculations('WPRTP_SGS', self.sam_file)
        
        wprtp = qtflux / (rho * 2.5104e+6) + WPRTP_SGS
        return wprtp, z
    
    def getCorrWpRtpCalc(self, dataset_override=None):
        """
        Calculates the correlation of W and THETAL from SAM output
        using the following equation:
        WPRTP/(np.sqrt(W2*QT2*1e-6)+1e-8)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        WPRTP, z, dataset = self.getVarForCalculations('WPRTP', self.sam_file)
        #RHO, z, dataset = self.getVarForCalculations('RHO', self.sam_file)
        WPRTP_SGS, z, dataset = self.getVarForCalculations('WPRTP_SGS', self.sam_file)
        W2, z, dataset = self.getVarForCalculations('W2', self.sam_file)
        QT2, z, dataset = self.getVarForCalculations('QT2', self.sam_file)
        
        CorrWpRtp = WPRTP / (np.sqrt(W2*QT2*1e-6)+1e-8)
        return CorrWpRtp, z
    
    def getWp2Calc(self, dataset_override = None):
        """
        Calculates the total vertical momentum variance W'^2 from SAM output
        using the following equation:
        W2 = W2 + WP2_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W2, z, dataset = self.getVarForCalculations('W2', self.sam_file)
        WP2_SGS, z, dataset = self.getVarForCalculations('WP2_SGS', self.sam_file)
        WP2 = W2 + WP2_SGS
        return WP2, z
    
    def getWp3Calc(self, dataset_override = None):
        """
        Calculates the total vertical momentum skewness W'^3 from SAM output
        using the following equation:
        W3 = W3 + WP3_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W3, z, dataset = self.getVarForCalculations('W3', self.sam_file)
        WP3_SGS, z, dataset = self.getVarForCalculations('WP3_SGS', self.sam_file)
        WP3 = W3 + WP3_SGS
        return WP3, z
    
    def getThetalVarCalc(self, dataset_override = None):
        """
        Calculates the total variance of THETAL from SAM output
        using the following equation:
        THETALVAR = THL2 + THLP2_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        TL2, z, dataset = self.getVarForCalculations('TL2', self.sam_file)
        THLP2_SGS, z, dataset = self.getVarForCalculations('THLP2_SGS', self.sam_file)
        THETALVAR = TL2 + THLP2_SGS
        return THETALVAR, z
    
    def getRtVarCalc(self, dataset_override = None):
        """
        TODO
        Calculates the total variance of RT from SAM output
        using the following equation:
        RTVAR = (QT2 * 1e-6) + RTP2_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        QT2, z, dataset = self.getVarForCalculations('QT2', self.sam_file)
        RTP2_SGS, z, dataset = self.getVarForCalculations('RTP2_SGS', self.sam_file)
        RTVAR = (QT2 * 1e-6) + RTP2_SGS
        return RTVAR, z
    
    def getUpWpCalc(self, dataset_override = None):
        """
        Calculates the total covariance of U and W from SAM output
        using the following equation:
        UPWP = UW + UPWP_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        UW, z, dataset = self.getVarForCalculations('UW', self.sam_file)
        UPWP_SGS, z, dataset = self.getVarForCalculations('UPWP_SGS', self.sam_file)
        UPWP = UW + UPWP_SGS
        return UPWP, z
    
    def getVpWpCalc(self, dataset_override = None):
        """
        Calculates the total covariance of V and W from SAM output
        using the following equation:
        VPWP = VW + VPWP_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        VW, z, dataset = self.getVarForCalculations('VW', self.sam_file)
        VPWP_SGS, z, dataset = self.getVarForCalculations('VPWP_SGS', self.sam_file)
        VPWP = VW + VPWP_SGS
        return VPWP, z
    
    def getUp2Calc(self, dataset_override = None):
        """
        Calculates the total variance of U from SAM output
        using the following equation:
        UVAR = U2 + UP2_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        U2, z, dataset = self.getVarForCalculations('U2', self.sam_file)
        UP2_SGS, z, dataset = self.getVarForCalculations('UP2_SGS', self.sam_file)
        UVAR = U2 + UP2_SGS
        return UVAR, z
    
    def getVp2Calc(self, dataset_override = None):
        """
        Calculates the total variance of V from SAM output
        using the following equation:
        VVAR = V2 + VP2_SGS (CLUBB variable)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        V2, z, dataset = self.getVarForCalculations('V2', self.sam_file)
        VP2_SGS, z, dataset = self.getVarForCalculations('VP2_SGS', self.sam_file)
        VVAR = V2 + VP2_SGS
        return VVAR, z
    
    def getUpWpCorrCalc(self, dataset_override = None):
        """
        Calculates the total correlation of U and W from SAM output
        using the following equation:
        (UW+UPWP_SGS)/(np.sqrt((U2+UP2_SGS)*(W2+WP2_SGS)+1e-4))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        UW, z, dataset = self.getVarForCalculations('UW', self.sam_file)
        UPWP_SGS, z, dataset = self.getVarForCalculations('UPWP_SGS', self.sam_file)
        UPWP = UW + UPWP_SGS
        return UPWP, z
    
    def getVpWpCorrCalc(self, dataset_override = None):
        """
        Calculates the total correlation of V and W from SAM output
        using the following equation:
        (VW+VPWP_SGS)/(np.sqrt((V2+VP2_SGS)*(W2+WP2_SGS)+1e-4))
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        VW, z, dataset = self.getVarForCalculations('VW', self.sam_file)
        VPWP_SGS, z, dataset = self.getVarForCalculations('VPWP_SGS', self.sam_file)
        VPWP = VW + VPWP_SGS
        return VPWP, z
    
    def getQRP2_QRIP(self, dataset_override=None):
        """
        Calculates the Within-Rain Variance of Rain Water Mixing Ratio from SAM output
        using the following equation:
        '(qrainp2_ip / (np.maximum(np.full(n,1e-5),qrainm_ip)**2))'
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        ALERT: How does this work in the budgets plotter?
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        qrainp2_ip, z, dataset = self.getVarForCalculations('qrainp2_ip', self.sam_file)
        qrainm_ip, z, dataset = self.getVarForCalculations('qrainm_ip', self.sam_file)
        QRP2_QRIP = qrainp2_ip / (np.maximum(np.full(n,1e-5),qrainm_ip)**2)
        return QRP2_QRIP, z
    
    ## Conditional average fallback functions
    def getUEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of U from SAM output
        using the following equation:
        (U - CLD * UCLD)/(1 - CLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        U, z, dataset = self.getVarForCalculations('U', self.sam_file)
        UCLD, z, dataset = self.getVarForCalculations('UCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UENV = (U - CLD * UCLD)/(1 - CLD)
        return UENV, z
    
    def getUEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of U from SAM output
        using the following equation:
        (U - CLD * UCLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        U, z, dataset = self.getVarForCalculations('U', self.sam_file)
        UCLD, z, dataset = self.getVarForCalculations('UCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UENV = (U - CLD * UCLD)
        return UENV, z
    
    def getUCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of U from SAM output
        using the following equation:
        CLD * UCLD
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        U, z, dataset = self.getVarForCalculations('U', self.sam_file)
        UCLD, z, dataset = self.getVarForCalculations('UCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UCLDW = CLD * UCLD
        return UCLDW, z
    
    def getVEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of V from SAM output
        using the following equation:
        (V - CLD * VCLD)/(1 - CLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        V, z, dataset = self.getVarForCalculations('V', self.sam_file)
        VCLD, z, dataset = self.getVarForCalculations('VCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VENV = (V - CLD * VCLD)/(1 - CLD)
        return VENV, z
    
    def getVEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of V from SAM output
        using the following equation:
        (V - CLD * VCLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        V, z, dataset = self.getVarForCalculations('V', self.sam_file)
        VCLD, z, dataset = self.getVarForCalculations('VCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VENV = (V - CLD * VCLD)
        return VENV, z
    
    def getVCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of V from SAM output
        using the following equation:
        CLD * VCLD
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        V, z, dataset = self.getVarForCalculations('V', self.sam_file)
        VCLD, z, dataset = self.getVarForCalculations('VCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VCLDW = CLD * VCLD
        return VCLDW, z
    
    def getWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of W from SAM output
        using the following equation:
        (W - CLD * WCLD)/(1 - CLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W, z, dataset = self.getVarForCalculations('W', self.sam_file)
        WCLD, z, dataset = self.getVarForCalculations('WCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        WENV = (W - CLD * WCLD)/(1 - CLD)
        return WENV, z
    
    def getWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of W from SAM output
        using the following equation:
        (W - CLD * WCLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W, z, dataset = self.getVarForCalculations('W', self.sam_file)
        WCLD, z, dataset = self.getVarForCalculations('WCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        WENV = (W - CLD * WCLD)
        return WENV, z

    def getWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of W from SAM output
        using the following equation:
        CLD * WCLD
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        W, z, dataset = self.getVarForCalculations('W', self.sam_file)
        WCLD, z, dataset = self.getVarForCalculations('WCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        WCLDW = CLD * WCLD
        return WCLDW, z
    
    
    def getUWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of UW from SAM output
        using the following equation:
        (UW - CLD * UWCLD)/(1 - CLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        UW, z, dataset = self.getVarForCalculations('UW', self.sam_file)
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UWENV = (UW - CLD * UWCLD)/(1 - CLD)
        return UWENV, z
    
    def getUWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of UW from SAM output
        using the following equation:
        (UW - CLD * UWCLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        UW, z, dataset = self.getVarForCalculations('UW', self.sam_file)
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UWENV = (UW - CLD * UWCLD)
        return UWENV, z
    
    def getUWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of UW from SAM output
        using the following equation:
        CLD * UWCLD
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        UW, z, dataset = self.getVarForCalculations('UW', self.sam_file)
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        UWCLDW = CLD * UWCLD
        return UWCLDW, z
    
    
    def getVWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of VW from SAM output
        using the following equation:
        (VW - CLD * VWCLD)/(1 - CLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        VW, z, dataset = self.getVarForCalculations('VW', self.sam_file)
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VWENV = (VW - CLD * VWCLD)/(1 - CLD)
        return VWENV, z
    
    def getVWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of VW from SAM output
        using the following equation:
        (VW - CLD * VWCLD)
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        VW, z, dataset = self.getVarForCalculations('VW', self.sam_file)
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VWENV = (VW - CLD * VWCLD)
        return VWENV, z
    
    def getVWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of VW from SAM output
        using the following equation:
        CLD * VWCLD
        :return: requested variable data in the form of a list. Returned data is already cropped to the appropriate min,max indices
        """
        if dataset_override is not None:
            dataset = dataset_override
        else:
            dataset = self.sam_file
        VW, z, dataset = self.getVarForCalculations('VW', self.sam_file)
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', self.sam_file)
        CLD, z, dataset = self.getVarForCalculations('CLD', self.sam_file)
        VWCLDW = CLD * VWCLD
        return VWCLDW, z