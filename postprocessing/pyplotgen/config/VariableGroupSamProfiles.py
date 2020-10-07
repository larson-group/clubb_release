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

    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None, silhs_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """

        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "sam profile variables"

        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)

        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000

        self.variable_definitions = [
            # THETA_L
            {'var_names':
                {
                    'clubb': ['thlm'],
                    'sam': ['THETAL'],
                    'silhs': [],
                    'coamps': ['thlm'],
                    'r408': ['thlm'],
                    'hoc': ['thlm'],
                    'e3sm': ['thlm'],
                    'cam': ['thlm'],
                },
                'sam_calc': self.getThlmSamCalc, 'sam_conv_factor': 1,
                'title': r"Liquid Water Potential Temperature, $\mathrm{\theta_l}$",
                'axis_title': r"$\mathrm{\theta_l}$ $\mathrm{\left[K\right]}$",
                'legend_label': r"$\mathrm{\theta_l}$",
            },
            # R_T
            {'var_names':
                {
                    'clubb': ['rtm'],
                    'sam': ['QT'],
                    'silhs': [],
                    'coamps': ['qtm'],
                    'r408': ['rtm'],
                    'hoc': ['rtm'],
                    'e3sm': ['rtm'],
                    'cam': ['rtm'],
                },
                'sam_calc': self.getRtmSamCalc, 'sam_conv_factor': 1,
                'title': r"Total Water Mixing Ratio, $\mathrm{r_t}$",
                'axis_title': r"$\mathrm{r_t}$ $\mathrm{\left[kg\,kg^{-1}\right]}$",
                'legend_label': r'$\mathrm{r_t}$',
            },
            # W'THETA_L'
            {'var_names':
                {
                    'clubb': ['wpthlp'],
                    'sam': ['TLFLUX'],
                    'silhs': [],
                    'coamps': ['wpthlp'],
                    'r408': ['wpthlp'],
                    'hoc': ['wpthlp'],
                    'e3sm': ['wpthlp'],
                },
                'sam_calc': self.getWpthlpCalc, 'sam_conv_factor': 1,
                'title': r"Turbulent Flux of $\mathrm{\theta_l}$",
                'axis_title': r"$\mathrm{\overline{w'\theta_l'}}$ / thflux(s) $\mathrm{\left[K\,m\,s^{-1}\right]}$",
                'legend_label': r"$\mathrm{\overline{w'\theta_l'}}$",
            },
            # CORR(W, THETA_L)
            {'var_names':
                {
                    'clubb': [],
                    'sam': [],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getCorrWpThlpCalc, 'sam_conv_factor': 1,
                'title': r"Corr(w,$\mathrm{\theta_l}$)",
                'axis_title': r"Correlation " +
                              r"$\mathrm{\overline{w'\theta_l'} / \sqrt{\overline{w'^2}\;\overline{\theta_l'^2}}}$" +
                              r" $\mathrm{\left[-\right]}$",
                'legend_label': r"$\mathrm{Corr(w,\theta_l)}$",
            },
            # W'R_T'
            {'var_names':
                {
                    'clubb': ['wprtp'],
                    'sam': ['QTFLUX'],
                    'silhs': [],
                    'coamps': ['wpqtp'],
                    'r408': ['wprtp'],
                    'hoc': ['wprtp'],
                    'e3sm': ['wprtp'],
                },
                'sam_calc': self.getWprtpCalc, 'sam_conv_factor': 1,
                'title': r"Turbulent Flux of $\mathrm{r_t}$",
                'axis_title': r"$\mathrm{\overline{w'r_t'}}$ (QTFLUX) $\mathrm{\left[kg\,kg^{-1}\,m\,s^{-1}\right]}$",
                'legend_label': r"$\mathrm{\overline{w'r_t'}}$",
            },
            # CORR(W, R_T)
            {'var_names':
                {
                    'clubb': [],
                    'sam': [],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getCorrWpRtpCalc, 'sam_conv_factor': 1,
                'title': r"$\mathrm{Corr(w,r_t)}$",
                'axis_title': r"Correlation " +
                              r"$\mathrm{\overline{w'r_t'} / \sqrt{\overline{w'^2}\;\overline{r_t'^2}}}$" +
                              r" $\mathrm{\left[-\right]}$",
                'legend_label': r"$\mathrm{Corr(w,r_t)}$",
            },
            # cloudliq_frac_em6
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['cloudliq_frac_em6'],
                    'silhs': [],
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
            # R_C
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['QCL'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1e-3,
                'title': r"Cloud Water Mixing Ratio, $\mathrm{r_c}$",
                'axis_title': r"$\mathrm{r_c}$ / $\mathrm{q_{cl}}$ $\mathrm{\left[kg\,kg^{-1}\right]}$",
                'legend_label': r'$\mathrm{r_c}$',
            },
            # W'^2
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['W2'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getWp2Calc, 'sam_conv_factor': 1,
                'title': r"Vertical Momentum Variance, $\mathrm{\overline{w'^2}}$",
                'axis_title': r"Momentum variance, $\mathrm{\overline{w'^2}}$ $\mathrm{\left[m^2\,s^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{w'^2}}$",
            },
            # W'^3
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['W3'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getWp3Calc, 'sam_conv_factor': 1,
                'title': r"Vertical Momentum Skewness, $\mathrm{\overline{w'^3}}$",
                'axis_title': r"Momentum Skewness, $\mathrm{\overline{w'^3}}$ $\mathrm{\left[m^3\,s^{-3}\right]}$",
                'legend_label': r"$\mathrm{\overline{w'^3}}$",
            },
            # THETA_L'^2
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['TL2'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getThetalVarCalc, 'sam_conv_factor': 1,
                'title': r"Variance of Liquid Water Potential Temperature, $\mathrm{\overline{\theta_l'^2}}$",
                'axis_title': r"$\mathrm{\overline{\theta_l'^2}}$ $\mathrm{\left[K^2\right]}$",
                'legend_label': r"$\mathrm{\overline{\theta_l'^2}}$",
            },
            # R_T'^2
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['QT2'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getRtVarCalc, 'sam_conv_factor': 1,
                'title': r"Variance of Total Water Mixing Ratio, $\mathrm{\overline{r_t'^2}}$",
                'axis_title': r"$\mathrm{\overline{r_t'^2}}$ / $\mathrm{\overline{q_t'^2}}$ " +
                              r"$\mathrm{\left[kg^2\,kg^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{r_t'^2}}$",
            },
            # R_T'THETA_L'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['RTPTHLP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"Covariance of $\mathrm{r_t}$ & $\mathrm{\theta_l}$",
                'axis_title': r"$\mathrm{\overline{r_t'\theta_l'}}$ $\mathrm{\left[kg\,kg^{-1}\,K\right]}$",
                'legend_label': r"$\mathrm{\overline{r_t'\theta_l'}}$",
            },
            # W_OBS
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['WOBS'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\mathrm{w_{obs}}$",
                'axis_title': r"Observed wind, $\mathrm{w_{obs}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r"$\mathrm{w_{obs}}$",
            },
            # WCLD unweighted
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['WCLD'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"In-cloud mean wind, $\overline{W}^{cld}$",
                'axis_title': r"In-cloud mean wind, " +
                              r"$\overline{w}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
            },
            # U
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['U'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"Zonal Wind Component, $\mathrm{\overline{u}}$",
                'axis_title': r"$\mathrm{\overline{u}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r'\mathrm{\overline{u}}',
            },
            # UCLD unweighted
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UCLD'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"In-cloud mean wind, $\overline{u}^{cld}$",
                'axis_title': r"In-cloud mean wind, " +
                              r"$\overline{u}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
            },
            # V
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['V'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"Meridonal Wind Component, $\mathrm{\overline{v}}$",
                'axis_title': r"$\mathrm{\overline{v}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r"$\mathrm{\overline{v}}$",
            },
            # VCLD unweighted
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VCLD'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"In-cloud mean wind, $\overline{v}^{cld}$",
                'axis_title': r"In-cloud mean wind, " +
                              r"$\overline{v}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
            },
            # U'W'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UW'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getUpWpCalc, 'sam_conv_factor': 1,
                'title': r"$\mathrm{\overline{u'w'}}$",
                'axis_title': r"Momentum flux, $\mathrm{\overline{u'w'}\ \left[m^2\,s^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{u'w'}}$",
            },
            # V'W'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VW'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getVpWpCalc, 'sam_conv_factor': 1,
                'title': r"$\mathrm{\overline{v'w'}}$",
                'axis_title': r"Momentum flux, $\mathrm{\overline{v'w'}\ \left[m^2\,s^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{v'w'}}$",
            },
            # U'^2
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['U2'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getUp2Calc, 'sam_conv_factor': 1,
                'title': r"$\mathrm{\overline{u'^2}}$",
                'axis_title': r"Momentum variance, $\mathrm{\overline{u'^2}\ \left[m^2\,s^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{u'^2}}$",
            },
            # V'^2
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['V2'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getVp2Calc, 'sam_conv_factor': 1,
                'title': r"$\mathrm{\overline{v'^2}}$",
                'axis_title': r"Momentum variance, $\mathrm{\overline{v'^2}\ \left[m^2\,s^{-2}\right]}$",
                'legend_label': r"$\mathrm{\overline{v'^2}}$",
            },
            # CORR(U, W)
            {'var_names':
                {
                    'clubb': [],
                    'sam': [],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getUpWpCorrCalc, 'sam_conv_factor': 1,
                'title': r"Corr(u,w)",
                'axis_title': r"Correlation " +
                              r"$\mathrm{\overline{u'w'} / \sqrt{\overline{u'^2}\;\overline{w'^2}}\ \left[-\right]}$",
                'legend_label': r"Corr(u,w)",
            },
            # CORR(V, W)
            {'var_names':
                {
                    'clubb': [],
                    'sam': [],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_calc': self.getVpWpCorrCalc, 'sam_conv_factor': 1,
                'title': r"Corr(v,w)",
                'axis_title': r"Correlation " +
                              r"$\mathrm{\overline{v'w'} / \sqrt{\overline{v'^2}\;\overline{w'^2}}\ \left[-\right]}$",
                'legend_label': r"Corr(v,w)",
            },
            # U'R_C'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UPRCP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{u'r_c'}$",
                'axis_title': r"Liquid water flux, " +
                              r"$\overline{u'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
                'legend_label': r"$\overline{u'r_c'}$",
            },
            # U'R_T'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UPRTP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{u'r_t'}$",
                'axis_title': r"Total water flux, " +
                              r"$\overline{u'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
                'legend_label': r"$\overline{u'r_t'}$",
            },
            # U'THETA_L'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UPTHLP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{u'\theta_l'}$",
                'axis_title': r"Liq. water pot. temp. flux, " +
                              r"$\overline{u'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{u'\theta_l'}$",
            },
            # U'THETA_V'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['UPTHVP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{u'\theta_v'}$",
                'axis_title': r"Virtual pot. temp. flux, " +
                              r"$\overline{u'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{u'\theta_v'}$",
            },
            # V'R_C'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VPRCP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{v'r_c'}$",
                'axis_title': r"Liquid water flux, " +
                              r"$\overline{v'r_c'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
                'legend_label': r"$\overline{v'r_c'}$",
            },
            # V'R_t'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VPRTP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{v'r_t'}$",
                'axis_title': r"Total water flux, " +
                              r"$\overline{v'r_t'}\ \mathrm{\left[m\,s^{-1}\,kg\,kg^{-1}\right]}$",
                'legend_label': r"$\overline{v'r_t'}$",
            },
            # V'THETA_L'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VPTHLP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{v'\theta_l'}$",
                'axis_title': r"Liq. water pot. temp. flux, " +
                              r"$\overline{v'\theta_l'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{v'\theta_l'}$",
            },
            # V'THETA_V'
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['VPTHVP'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': r"$\overline{v'\theta_v'}$",
                'axis_title': r"Virtual pot.temp. flux, " +
                              r"$\overline{v'\theta_v'}\ \mathrm{\left[K\,m\,s^{-1}\right]}$",
                'legend_label': r"$\overline{v'\theta_v'}$",
            },

            # CLD
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['CLD'],
                    'silhs': [],
                    'coamps': [],
                    'r408': [],
                    'hoc': [],
                    'e3sm': [],
                    'cam': [],
                },
                'sam_conv_factor': 1,
                'title': "Cloud fraction",
                'axis_title': "Cloud fraction [-]",
                'legend_label': 'CLD',
            },

            ## Tracer variables

            # TR01
            {'var_names':
                {
                    'clubb': [],
                    'sam': ['TR01'],
                    'silhs': [],
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
                    'silhs': [],
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
                    'silhs': [],
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
                    'silhs': [],
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
                    'silhs': [],
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
                    'silhs': [],
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
                    'silhs': [],
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
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, les_dataset=les_dataset,
                         coamps_dataset=coamps_dataset, r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets,
                         silhs_datasets=silhs_datasets)

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

        CorrWpThlp = (TLFLUX / (RHO * 1004) + WPTHLP_SGS) / np.sqrt(W2 * TL2 + 1e-4)
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

        CorrWpRtp = WPRTP / (np.sqrt(W2 * QT2 * 1e-6) + 1e-8)
        return CorrWpRtp, indep

    def getWp2Calc(self, dataset_override=None):
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

    def getWp3Calc(self, dataset_override=None):
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

    def getThetalVarCalc(self, dataset_override=None):
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

    def getRtVarCalc(self, dataset_override=None):
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

    def getUpWpCalc(self, dataset_override=None):
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

    def getVpWpCalc(self, dataset_override=None):
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

    def getUp2Calc(self, dataset_override=None):
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

    def getVp2Calc(self, dataset_override=None):
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

    def getUpWpCorrCalc(self, dataset_override=None):
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

    def getVpWpCorrCalc(self, dataset_override=None):
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
        QRP2_QRIP = qrainp2_ip / (np.maximum(qrainm_ip, 1e-5) ** 2)
        return QRP2_QRIP, indep
