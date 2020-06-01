'''
:author: Steffen Domke
:date: Mid 2019
TODO:   - Should invalid variables produce zero lines or none at all?
        - fallback <-> calc
        - Update axis labels at TASK
        - Implement/copy sam calcs
        - Figure out how to include standard plots in VariableGroupBase (labels etc.)
'''

from src.Panel import Panel
from src.VariableGroup import VariableGroup


class VariableGroupSamMultilineProfiles(VariableGroup):
    """

    """
    def __init__(self, case, clubb_datasets=None, les_dataset=None, coamps_dataset=None, r408_dataset=None,
                 hoc_dataset=None, cam_datasets=None,
                 e3sm_datasets=None, sam_datasets=None, wrf_datasets=None):
        """
        
        :param clubb_datasets:
        :param case:
        :param les_dataset:
        """
        self.name = "sam multiline profile variables"
        
        self.kg_per_second_to_kg_per_day = 1. / (24 * 3600)
        
        self.g_per_second_to_kg_per_day = self.kg_per_second_to_kg_per_day / 1000
        
        u_cond_lines = [
            {'var_names': ['U'], 'legend_label': r"$\overline{u}$"},
            {'var_names': ['UCLD'], 'legend_label': '$\overline{u}^\mathrm{cld}$'},
            {'var_names': ['UENV'], 'legend_label': r"$\overline{u}^\mathrm{{env}}$",
             'sam_calc': self.getUEnvUnweighted},
            ]
        
        u_weight_lines = [
            {'var_names': ['U'], 'legend_label': r"$\overline{u}$"},
            {'var_names': ['UCLD'], 'legend_label': r"$\overline{u}^\mathrm{{cld}}$",
             'sam_calc': self.getUCldWeighted},
            {'var_names': ['UENV'], 'legend_label': r"$\overline{u}^\mathrm{{env}}$",
             'sam_calc': self.getUEnvWeighted},
            ]
        
        v_cond_lines = [
            {'var_names': ['V'], 'legend_label': r"$\overline{v}$"},
            {'var_names': ['VCLD'], 'legend_label': '$\overline{v}^\mathrm{cld}$'},
            {'var_names': ['VENV'], 'legend_label': r"$\overline{v}^\mathrm{{env}}$",
             'sam_calc': self.getVEnvUnweighted},
            ]
        
        v_weight_lines = [
            {'var_names': ['V'], 'legend_label': r"$\overline{v}$"},
            {'var_names': ['VCLD'], 'legend_label': r"$\overline{v}^\mathrm{{cld}}$",
             'sam_calc': self.getVCldWeighted},
            {'var_names': ['VENV'], 'legend_label': r"$\overline{v}^\mathrm{{env}}$",
             'sam_calc': self.getVEnvWeighted},
            ]
        
        w_cond_lines = [
            {'var_names': ['WM'], 'legend_label': r"$\overline{w}$"},
            {'var_names': ['WCLD'], 'legend_label': '$\overline{w}^\mathrm{cld}$'},
            {'var_names': ['WENV'], 'legend_label': r"$\overline{w}^\mathrm{{env}}$",
             'sam_calc': self.getWEnvUnweighted},
            ]
        
        w_weight_lines = [
            {'var_names': ['WM'], 'legend_label': r"$\overline{w}$"},
            {'var_names': ['WCLD'], 'legend_label': r"$\overline{w}^\mathrm{{cld}}$",
             'sam_calc': self.getWCldWeighted},
            {'var_names': ['WENV'], 'legend_label': r"$\overline{w}^\mathrm{{env}}$",
             'sam_calc': self.getWEnvWeighted},
            ]
        
        uw_cond_lines = [
            {'var_names': ['UW'], 'legend_label': r"$\overline{u'w'}$"},
            {'var_names': ['UWCLD'], 'legend_label': r"$\overline{u'w'}^\mathrm{cld}$"},
            {'var_names': ['UWENV'], 'legend_label': r"$\overline{u'w'}^\mathrm{env}$",
             'sam_calc': self.getUWEnvUnweighted},
            ]
        
        uw_weight_lines = [
            {'var_names': ['UW'], 'legend_label': r"$\overline{u'w'}$"},
            {'var_names': ['UWCLD'], 'legend_label': r"$\overline{u'w'}^\mathrm{cld}$",
             'sam_calc': self.getUWCldWeighted},
            {'var_names': ['UWENV'], 'legend_label': r"$\overline{u'w'}^\mathrm{env}$",
             'sam_calc': self.getUWEnvWeighted},
            ]
        
        vw_cond_lines = [
            {'var_names': ['VW'], 'legend_label': r"$\overline{v'w'}$"},
            {'var_names': ['VWCLD'], 'legend_label': r"$\overline{v'w'}^\mathrm{cld}$"},
            {'var_names': ['VWENV'], 'legend_label': r"$\overline{v'w'}^\mathrm{env}$",
             'sam_calc': self.getVWEnvUnweighted},
            ]
        
        vw_weight_lines = [
            {'var_names': ['VW'], 'legend_label': r"$\overline{v'w'}$"},
            {'var_names': ['VWCLD'], 'legend_label': r"$\overline{v'w'}^\mathrm{cld}$",
             'sam_calc': self.getVWCldWeighted},
            {'var_names': ['VWENV'], 'legend_label': r"$\overline{v'w'}^\mathrm{env}$",
             'sam_calc': self.getVWEnvWeighted},
            ]
        
        thv_cond_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$"},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
             'sam_calc': self.getTHVEnvUnweighted},
            ]
        
        thv_weight_lines = [
            {'var_names': ['THETAV'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TVCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$",
             'sam_calc': self.getTHVCldWeighted},
            {'var_names': ['TVENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
             'sam_calc': self.getTHVEnvWeighted},
            ]
        
        #thl_weight_lines = [
        ## Cloud conditional data does not exist. TLCLD is cld cond. liquid water static energy
            #{'var_names': ['THETAL'], 'legend_label': r"$\overline{\theta_v}$"},
            #{'var_names': ['TLCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$",
            # 'sam_calc': self.getTHLCldWeighted},
            #{'var_names': ['TLENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
            # 'sam_calc': self.getTHLEnvWeighted},
            #]
        
        qt_cond_lines = [
            {'var_names': ['QT'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['QTCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$"},
            {'var_names': ['QTENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
             'sam_calc': self.getQTEnvUnweighted},
            ]
        
        qtw_weight_lines = [
             # QTFLUX is only defined in microphysics packages
            {'var_names': ['QTFLUX'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['QTWCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$",
             'sam_calc': self.getQTWCldWeighted},
            {'var_names': ['QTWENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
             'sam_calc': self.getQTWEnvWeighted},
            ]
        
        tlw_weight_lines = [
            {'var_names': ['TLFLUX'], 'legend_label': r"$\overline{\theta_v}$"},
            {'var_names': ['TLWCLD'], 'legend_label': r"$\overline{\theta_v}^\mathrm{cld}$",
             'sam_calc': self.getTLWCldWeighted},
            {'var_names': ['TLWENV'], 'legend_label': r"$\overline{\theta_v}^\mathrm{env}$",
             'sam_calc': self.getTLWEnvWeighted},
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
            ## Cloud conditional plots
            
            # UCLD unweighted
            # Display U, unweighted UCLD, and unweighted calculated UENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': u_cond_lines,
             'title': r"Conditional mean wind, $\overline{u}$",
             'axis_title': r"Conditional mean wind, "+
                           r"$\overline{u}^\mathrm{{cld}}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # UCLD weighted
            # Display U, weighted calculated UCLD, and weighted calculated UENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': u_weight_lines,
             'title': r"Weighted mean wind, $\overline{u}$",
             'axis_title': r"Weighted mean wind, $\overline{u}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # VCLD unweighted
            # Display V, unweighted VCLD, and unweighted calculated VENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': v_cond_lines,
             'title': r"Conditional mean wind, $\overline{v}$",
             'axis_title': r"Conditional mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # VCLD weighted
            # Display V, weighted calculated VCLD, and weighted calculated VENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': v_weight_lines,
             'title': r"Weighted mean wind, $\overline{v}$",
             'axis_title': r"Weighted mean wind, $\overline{v}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # WCLD unweighted
            # Display W, unweighted WCLD, and unweighted calculated WENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': w_cond_lines,
             'title': r"Conditional mean wind, $\overline{w}$",
             'axis_title': r"Conditional mean wind, $\overline{w}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # WCLD weighted
            # Display W, weighted calculated WCLD, and weighted calculated WENV
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
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': w_weight_lines,
             'title': r"Weighted mean wind, $\overline{w}$",
             'axis_title': r"Weighted mean wind, $\overline{w}\ \mathrm{\left[m\,s^{-1}\right]}$",
            },
            # UWCLD unweighted
            # Display UW, unweighted UWCLD, and unweighted calculated UWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['UWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_cond_lines,
             'title': r"Cloud-conditonal $\overline{u'w'}$",
             'axis_title': r"Conditional flux, $\overline{u'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            },
            # UWCLD weighted
            # Display UW, weighted calculated UWCLD, and weighted calculated UWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['UWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_weight_lines,
             'title': r"Cloud-weighted $\overline{u'w'}$",
             'axis_title': r"Weighted flux, $\overline{u'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            },
            # VWCLD unweighted
            # Display VW, unweighted VWCLD, and unweighted calculated VWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['VWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_cond_lines,
             'title': r"Cloud-conditonal $\overline{v'w'}$",
             'axis_title': r"Conditional flux, $\overline{v'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            },
            # VWCLD weighted
            # Display VW, weighted calculated VWCLD, and weighted calculated VWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['VWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_weight_lines,
             'title': r"Cloud-weighted $\overline{v'w'}$",
             'axis_title': r"Weighted flux, $\overline{v'w'}\ \mathrm{\left[m^2\, s^{-2}\right]}$",
            },
            # TVCLD unweighted
            # Display THETAV, unweighted TVCLD, and unweighted calculated UWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['TVCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': thv_cond_lines,
             'title': r"Conditional virt. pot. temp., $\overline{\theta}_v$",
             'axis_title': r"Conditional virt. pot. temp., $\overline{\theta}_v\ \mathrm{\left[K\right]}$",
            },
            # TVCLD weighted
            # Display THETAV, weighted calculatedTVCLD, and weighted calculated TVENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['TVCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': thv_weight_lines,
             'title': r"Weighted virt. pot. temp., $\overline{\theta}_v$",
             'axis_title': r"Weighted virt. pot. temp., $\overline{\theta}_v\ \mathrm{\left[K\right]}$",
            },
            # TLWCLD weighted
            # Display TLFLUX, weighted calculated TLWCLD, and weighted calculated TLWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['TLWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': tlw_weight_lines,
             'title': r"Weighted flux, $\overline{w's_L'}$",
             'axis_title': r"Weighted flux, $\overline{w's_L'}\ \mathrm{\left[K\,m\, s^{-1}\right]}$",
            },
            # QTWCLD weighted
            # Display QTFLUX, weighted calculated QTWCLD, and weighted calculated QTWENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['QTWCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': qtw_weight_lines,
             'title': r"Weighted flux, $\overline{w'r_t'}$",
             'axis_title': r"Weighted flux, $\overline{w'r_t'}\ \mathrm{\left[kg\,kg^{-1}\, m\,s^{-1}\right]}$",
            },
            # QTCLD unweighted
            # Display QT, unweighted QTCLD, and unweighted QTENV
            {'var_names':
                {
                'clubb': [],
                'sam': ['QTCLD'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': qt_cond_lines,
             'title': r"Conditional total water mixing ratio, $\mathrm{r_t}$",
             'axis_title': r"Conditional total water mixing ratio, "+
                           r"$\mathrm{r_t}\ \mathrm{\left[g\,kg^{-1}\right]}$",
            },
            # U 2ND MOMENTS
            # Display UW, U2, W2
            {'var_names':
                {
                'clubb': [],
                'sam': ['UW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': uw_lines,
             'title': "Eastward 2nd-moments",
             'axis_title': r"2nd moments $\mathrm{\left[m^2\,s^{-2}\right]}$",
            },
            # V 2ND MOMENTS
            # Display VW, V2, W2
            {'var_names':
                {
                'clubb': [],
                'sam': ['VW'],
                'coamps': [],
                'r408': [],
                'hoc': [],
                'e3sm': [],
				'cam': [],
                },
             'sam_conv_factor': 1, 'type': Panel.TYPE_BUDGET, 'lines': vw_lines,
             'title': "Northward 2nd-moments",
             'axis_title': r"2nd moments $\mathrm{\left[m^2\,s^{-2}\right]}$",
            },
        ]
        
        # Call ctor of parent class
        super().__init__(case, clubb_datasets=clubb_datasets, sam_datasets=sam_datasets, les_dataset=les_dataset,
                         coamps_dataset=coamps_dataset, r408_dataset=r408_dataset, cam_datasets=cam_datasets,
                         hoc_dataset=hoc_dataset, e3sm_datasets=e3sm_datasets, wrf_datasets=wrf_datasets)

    def getUpWpCalc(self, dataset_override = None):
        """
        Calculates the total covariance of U and W from SAM output
        using the following equation:
        UPWP = UW + UPWP_SGS (CLUBB variable)
        :return: requested variable dependent_data in the form of a list.
                 Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UW, z, dataset = self.getVarForCalculations('UW', dataset)
        UPWP_SGS, z, dataset = self.getVarForCalculations('UPWP_SGS', dataset)
        UPWP = UW + UPWP_SGS
        return UPWP, z
    
    def getVpWpCalc(self, dataset_override = None):
        """
        Calculates the total covariance of V and W from SAM output
        using the following equation:
        VPWP = VW + VPWP_SGS (CLUBB variable)
        :return: requested variable dependent_data in the form of a list.
                 Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VW, z, dataset = self.getVarForCalculations('VW', dataset)
        VPWP_SGS, z, dataset = self.getVarForCalculations('VPWP_SGS', dataset)
        VPWP = VW + VPWP_SGS
        return VPWP, z
    
    def getUp2Calc(self, dataset_override = None):
        """
        Calculates the total variance of U from SAM output
        using the following equation:
        UVAR = U2 + UP2_SGS (CLUBB variable)
        :return: requested variable dependent_data in the form of a list.
                 Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        U2, z, dataset = self.getVarForCalculations('U2', dataset)
        UP2_SGS, z, dataset = self.getVarForCalculations('UP2_SGS', dataset)
        UVAR = U2 + UP2_SGS
        return UVAR, z
    
    def getVp2Calc(self, dataset_override = None):
        """
        Calculates the total variance of V from SAM output
        using the following equation:
        VVAR = V2 + VP2_SGS (CLUBB variable)
        :return: requested variable dependent_data in the form of a list.
                 Returned dependent_data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        V2, z, dataset = self.getVarForCalculations('V2', dataset)
        VP2_SGS, z, dataset = self.getVarForCalculations('VP2_SGS', dataset)
        VVAR = V2 + VP2_SGS
        return VVAR, z

    ## Conditional average fallback functions
    # UCLD
    def getUEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of U from SAM output
        using the following equation:
        (U - CLD * UCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        U, z, dataset = self.getVarForCalculations('U', dataset)
        UCLD, z, dataset = self.getVarForCalculations('UCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UENV = (U - CLD * UCLD) / (1 - CLD)
        return UENV, z
    
    def getUEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of U from SAM output
        using the following equation:
        (U - CLD * UCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        U, z, dataset = self.getVarForCalculations('U', dataset)
        UCLD, z, dataset = self.getVarForCalculations('UCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UENV = (U - CLD * UCLD)
        return UENV, z
    
    def getUCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of U from SAM output
        using the following equation:
        CLD * UCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UCLD, z, dataset = self.getVarForCalculations('UCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UCLDW = CLD * UCLD
        return UCLDW, z
    
    # VCLD
    def getVEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of V from SAM output
        using the following equation:
        (V - CLD * VCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        V, z, dataset = self.getVarForCalculations('V', dataset)
        VCLD, z, dataset = self.getVarForCalculations('VCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VENV = (V - CLD * VCLD) / (1 - CLD)
        return VENV, z
    
    def getVEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of V from SAM output
        using the following equation:
        (V - CLD * VCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        V, z, dataset = self.getVarForCalculations('V', dataset)
        VCLD, z, dataset = self.getVarForCalculations('VCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VENV = (V - CLD * VCLD)
        return VENV, z
    
    def getVCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of V from SAM output
        using the following equation:
        CLD * VCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VCLD, z, dataset = self.getVarForCalculations('VCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VCLDW = CLD * VCLD
        return VCLDW, z
    
    # WCLD
    def getWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of W from SAM output
        using the following equation:
        (WM - CLD * WCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WM, z, dataset = self.getVarForCalculations('WM', dataset)
        WCLD, z, dataset = self.getVarForCalculations('WCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        WENV = (WM - CLD * WCLD) / (1 - CLD)
        return WENV, z
    
    def getWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of W from SAM output
        using the following equation:
        (WM - CLD * WCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WM, z, dataset = self.getVarForCalculations('WM', dataset)
        WCLD, z, dataset = self.getVarForCalculations('WCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        WENV = (WM - CLD * WCLD)
        return WENV, z

    def getWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of W from SAM output
        using the following equation:
        CLD * WCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WCLD, z, dataset = self.getVarForCalculations('WCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        WCLDW = CLD * WCLD
        return WCLDW, z
    
    # UWCLD
    def getUWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of UW from SAM output
        using the following equation:
        (UW - CLD * UWCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UW, z, dataset = self.getVarForCalculations('UW', dataset)
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UWENV = (UW - CLD * UWCLD) / (1 - CLD)
        return UWENV, z
    
    def getUWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of UW from SAM output
        using the following equation:
        (UW - CLD * UWCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UW, z, dataset = self.getVarForCalculations('UW', dataset)
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UWENV = (UW - CLD * UWCLD)
        return UWENV, z
    
    def getUWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of UW from SAM output
        using the following equation:
        CLD * UWCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        UWCLD, z, dataset = self.getVarForCalculations('UWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        UWCLDW = CLD * UWCLD
        return UWCLDW, z
    
    # VWCLD
    def getVWEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of VW from SAM output
        using the following equation:
        (VW - CLD * VWCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VW, z, dataset = self.getVarForCalculations('VW', dataset)
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VWENV = (VW - CLD * VWCLD) / (1 - CLD)
        return VWENV, z
    
    def getVWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of VW from SAM output
        using the following equation:
        (VW - CLD * VWCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VW, z, dataset = self.getVarForCalculations('VW', dataset)
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VWENV = (VW - CLD * VWCLD)
        return VWENV, z
    
    def getVWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of VW from SAM output
        using the following equation:
        CLD * VWCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        VWCLD, z, dataset = self.getVarForCalculations('VWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        VWCLDW = CLD * VWCLD
        return VWCLDW, z
    
    # TVCLD
    def getTHVEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of THETAV from SAM output
        using the following equation:
        (THETAV - CLD * TVCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)
        TVCLD, z, dataset = self.getVarForCalculations('TVCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        TVENV = (THETAV - CLD * TVCLD) / (1 - CLD)
        return TVENV, z
    
    def getTHVEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of THETAV from SAM output
        using the following equation:
        (THETAV - CLD * TVCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        THETAV, z, dataset = self.getVarForCalculations('THETAV', dataset)
        TVCLD, z, dataset = self.getVarForCalculations('TVCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        TVENV = (THETAV - CLD * TVCLD)
        return TVENV, z
    
    def getTHVCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of THETAV from SAM output
        using the following equation:
        CLD * TVCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TVCLD, z, dataset = self.getVarForCalculations('TVCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        TVCLDW = CLD * TVCLD
        return TVCLDW, z
    
    # QTCLD
    # ALERT: Check equations/units of SAM QT variables
    def getQTEnvUnweighted(self, dataset_override=None):
        """
        Calculates the unweighted environment-conditional average of QT from SAM output
        using the following equation:
        (QT - CLD * QTCLD) / (1 - CLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QT, z, dataset = self.getVarForCalculations('QT', dataset)
        QTCLD, z, dataset = self.getVarForCalculations('QTCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        QTENV = (QT - CLD * QTCLD) / (1 - CLD)
        return QTENV, z
    
    # QTWCLD
    # ALERT: Check equations/units of SAM QTW variables
    def getQTWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of QTW from SAM output
        using the following equation:
        (WPRTP - CLD * QTWCLD / 1000)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        WPRTP, z, dataset = self.getVarForCalculations('WPRTP', dataset)
        QTWCLD, z, dataset = self.getVarForCalculations('QTWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        QTWENV = (WPRTP - CLD * QTWCLD / 1000)
        return QTWENV, z
    
    def getQTWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of QTW from SAM output
        using the following equation:
        CLD * QTWCLD / 1000
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        QTWCLD, z, dataset = self.getVarForCalculations('QTWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        QTWCLDW = CLD * QTWCLD
        return QTWCLDW, z
    
    # TLWCLD
    # ALERT: Check equations/units of SAM TLW variables
    def getTLWEnvWeighted(self, dataset_override=None):
        """
        Calculates the weighted environment-conditional average of TL from SAM output
        using the following equation:
        (TLFLUX / RHO / 1004 - CDL * TLWCLD)
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TLFLUX, z, dataset = self.getVarForCalculations('TLFLUX', dataset)
        TLWCLD, z, dataset = self.getVarForCalculations('TLWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        RHO, z, dataset = self.getVarForCalculations('RHO', dataset)
        TLWENV = (TLFLUX / RHO / 1004 - CLD * TLWCLD)
        return TLWENV, z
    
    def getTLWCldWeighted(self, dataset_override=None):
        """
        Calculates the weighted cloud-conditional average of TL from SAM output
        using the following equation:
        CLD * TLWCLD
        :return: requested variable data in the form of a list.
                 Returned data is already cropped to the appropriate min,max indices
        """
        dataset = self.les_dataset
        if dataset_override is not None:
            dataset = dataset_override
        TLWCLD, z, dataset = self.getVarForCalculations('TLWCLD', dataset)
        CLD, z, dataset = self.getVarForCalculations('CLD', dataset)
        TLWCLDW = CLD * TLWCLD
        return TLWCLDW, z
