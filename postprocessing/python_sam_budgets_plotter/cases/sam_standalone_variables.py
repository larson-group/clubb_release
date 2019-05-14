"""
-------------------------------------------------------------------------------
   G E N E R A L   I N F O R M A T I O N
-------------------------------------------------------------------------------
This file contains general constants and information about the SAM variables saved
in the netCDF file needed for plotgen.py.

The list variables sortPlots, plotNames and lines are sorted identically in
order to relate the individual variables.
"""

#-------------------------------------------------------------------------------
#   I M P O R T S
#-------------------------------------------------------------------------------
from numpy import nan

#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24                                                    # 1d  = 24h
HOUR = 3600                                                 # 1h  = 3600s
KG = 1000.                                                  # 1kg = 1000g
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)
filler = nan                                                # Define the fill value which should replace invalid values in the data
startLevel = 0                                              # Set the lower height level at which the plots should begin. For example, startLevel=2 would cut off the lowest 2 data points for each line.
header = 'SAM standalone profiles'
name = 'sam_standalone'                                     # String used as part of the output file name
nc_files = ['sam']                                          # NetCDF files needed for plots, paths are defined

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
# Names of the variables
sortPlots = ['theta_l', 'r_t', 'theta_l_flux', 'r_t_flux', 'cloudliq_frac', 'r_c', 'w_var', 'w3', 'theta_l_var', 'r_t_var', 'covar_thetal_rt', 'wobs', 'U', 'V', 'covar_uw', 'covar_vw', 'u_var', 'v_var', 'corr_uw', 'corr_vw',\
                'QR', 'QR_IP', 'QRP2', 'QRP2_QRIP', \
                'Nrm', 'Nrm_IP', 'Nrp2', 'Nrp2_NrmIP', \
                'Ncm', 'Ncm_IP', 'Ncp2', 'Ncp2_NcmIP', \
                'Ngm', 'Ngm_IP', 'Ngp2', 'Ngp2_NgmIP', \
                'Qgm', 'Qgm_IP', 'Qgp2', 'Qgp2_QgmIP', \
                'Qim', 'Qim_IP', 'Qip2', 'Qip2_QimIP', \
                'Nim', 'Nim_IP', 'Nip2', 'Nip2_NimIP', \
                'Qsm', 'Qsm_IP', 'Qsp2', 'Qsp2_QsmIP', \
                'Nsm', 'Nsm_IP', 'Nsp2', 'Nsp2_NsmIP', \
                'MicroFractions', 'Buoy_Flux', \
                'uprcp', 'uprtp', 'upthlp', 'upthvp', \
                'vprcp', 'vprtp', 'vpthlp', 'vpthvp', \
                'ucld' , 'vcld' , 'wcld',\
                'ucomp', 'vcomp', 'wcomp', 'uwcomp', 'vwcomp', 'tvcomp', 'tlcomp', 'qtwcomp',\
                ]
                
# settings of each plot:
# plot number, plot title, axis label
plotNames = [\
                ['Liquid Water Potential Temperature, theta_l', 'thetal [K]'],\
                ['Total Water Mixing Ratio, r_t', 'rtm / qt [kg/kg]'],\
                ['Turbulent Flux of theta_l', 'wpthlp / thflux(s) [K m/s]'],\
                ['Turbulent Flux of r_t', 'wprtp / qtflux(s) [(kg/kg) m/s]'],\
                ['Cloud Liquid Fraction', ' [%/100]'],\
                ['Cloud Water Mixing Ratio, r_c', 'rcm / qcl [kg/kg]'],\
                [r"$\mathrm{\overline{w'^2}}$", r"Momentum variance $\mathrm{\overline{w'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{w'^3}}$", r"Third-order Moment $\mathrm{\overline{w'^3}\ \left[\frac{m^3}{s^3}\right]}$"],\
                ['Variance of theta_l', 'thlp2 / tl2 [K^2]'],\
                ['Variance of r_t', 'rtp2 / qtp2 [(kg/kg)^2]'],\
                ['Covariance of r_t & theta_l', 'rtpthlp [(kg/kg) K]'],\
                [r"$\mathrm{w_{obs}}$", r"Observed wind $\mathrm{w_{obs}\ \left[\frac{m}{s}\right]}$"],\
                [r"Zonal Wind Component, u", r"$\mathrm{\bar{u}\ \left[\frac{m}{s}\right]}$"],\
                [r"Meridonal Wind Component, v", r"$\mathrm{\bar{v}\ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\overline{u'w'}}$", r"Momentum flux $\mathrm{\overline{u'w'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'w'}}$", r"Momentum flux $\mathrm{\overline{v'w'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{u'^2}}$", r"Momentum variance $\mathrm{\overline{u'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'^2}}$", r"Momentum variance $\mathrm{\overline{v'^2}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"Corr(u,w)", r"Flux correlation $\mathrm{\overline{u'w'} / \sqrt{\overline{u'^2}\overline{w'^2}}\ [-]}$"],\
                [r"Corr(v,w)", r"Flux correlation $\mathrm{\overline{v'w'} / \sqrt{\overline{v'^2}\overline{w'^2}}\ [-]}$"],\
                # Rain Water Mixing Ratio
                ['Rain Water Mixing Ratio', 'qrm [kg/kg]'],\
                ['Rain Water Mixing Ratio in Rain', 'qrm_ip [kg/kg]'],\
                ['Domain-wide Variance\nof Rain Water Mixing Ratio', 'qrp2 [(kg/kg)^2]'],\
                ['Within-rain Variance\nof Rain Water Mixing Ratio', 'qrp2_ip / qrm_ip^2 [-]'],\
                #Rain Drop Number Concentration
                ['Rain Drop Concentration', 'Nrm [num/kg]'],\
                ['Rain Drop Concentration in Rain', 'Nrm_ip [num/kg]'],\
                ['Domain-wide Variance\nof Rain Drop Concentration', 'Nrp2 [(num/kg)^2]'],\
                ['Within-rain Variance\nof Rain Drop Concentration', 'Nrp2_ip / Nrm_ip^2 [-]'],\
                #Cloud Droplet Number Concentration
                ['Cloud Droplet Number Concentration', 'Ncm [num/kg]'],\
                ['Cloud Droplet Number Concentration in Cloud', 'Ncm_ip [num/kg]'],\
                ['Domain-wide Variance\nof Cloud Droplet Number Concentration', 'Ncp2 [(#/kg)^2]'],\
                ['Within-cloud Variance\nof Cloud Droplet Number Concentration', 'Ncp2_ip / Ncm_ip^2 [-]'],\
                #Graupel Number Concentration
                ['Graupel Number Concentration', 'Ngm [kg/kg]'],\
                ['Graupel Number Concentration in Graupel', 'Ngm_ip [num/kg]'],\
                ['Domain-wide Variance\nof Graupel Number Concentration', 'Ngp2 [(kg/kg)^2]'],\
                ['Within-graupel Variance\nof Graupel Number Concentration', 'Ngp2_ip / Ngm_ip^2 [-]'],\
                #Graupel Mixing Ratio
                ['Graupel Mixing Ratio', 'qgm [kg/kg]'],\
                ['Graupel Mixing Ratio in Graupel', 'qgm_ip [kg/kg]'],\
                ['Domain-wide Variance\nof Graupel Mixing Ratio', 'qgp2 [(kg/kg)^2]'],\
                ['Within-graupel Variance\nof Graupel Mixing Ratio', 'qgp2_ip / qgm_ip^2 [-]'],\
                #Cloud Ice Mixing Ratio
                ['Cloud Ice Mixing Ratio', 'qim [kg/kg]'],\
                ['Cloud Ice Mixing Ratio in Cloud Ice', 'qim_ip [kg/kg]'],\
                ['Domain-wide Variance\nof Cloud Ice Mixing Ratio', 'qip2 [(kg/kg)^2]'],\
                ['Within-cloud-ice Variance\nof Cloud Ice  Mixing Ratio', 'qip2_ip / qim_ip^2 [-]'],\
                #Cloud Ice Number Concentration
                ['Cloud Ice Concentration', 'Nim [num/kg]'],\
                ['Cloud Ice Number Concentration in Cloud Ice', 'Ni_ip [num/kg]'],\
                ['Domain-wide Variance\nof Cloud Ice Number Concentration', 'Nip2 [(num/kg)^2]'],\
                ['Within-cloud-ice Variance\nof Cloud Ice Number Concentration', 'Nip2_ip / Nim_ip^2 [-]'],\
                #Snow Mixing Ratio
                ['Snow Mixing Ratio ', 'qsm [kg/kg]'],\
                ['Snow Mixing Ratio in Snow', 'qsm_ip [kg/kg]'],\
                ['Domain-wide Variance\nof Snow Mixing Ratio', 'qsp2 [(kg/kg)^2]'],\
                ['Within-snow Variance\nof Snow Mixing Ratio ', 'qsp2_ip / qsm_ip^2 [-]'],\
                #Snow Number Concentration
                ['Snow Number Concentration', 'Nsm [num/kg]'],\
                ['Snow Number Concentration in Snow', 'Nsm_ip [num/kg]'],\
                ['Domain-wide Variance\nof Snow Number Concentration', 'Nsp2 [(#/kg)^2]'],\
                ['Within-snow Variance\nof Snow Number Concentration', 'Nsp2_ip / Nsm_ip^2 [-]'],\
                ['Micro Fractions', '[%/100]'],\
                ['Buoyancy flux', 'wpthvp / tlflux [K m/s]'],\
                #['Liquid Water Path', 'lwp [kg/m^2]'],\
                #['Surface rainfall rate', 'rain_rate_sfc[mm/day]'],\
                #['Density-Weighted Vertically Averaged wp2', 'wp2 / w2 [m^2/s^2]'],\
                #['Cloud Ice Water Path', 'iwp [kg/m^2]'],\
                #['Snow Water Path', 'swp [kg/m^2]'],\
                # buoyancy sub-terms for parameterization in upwp budget
                [r"$\mathrm{\overline{u'r_c'}}$", r"Liquid water flux $\mathrm{\overline{u'r_c'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{u'r_t'}}$", r"Total water flux $\mathrm{\overline{u'r_t'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{u'\theta_l'}}$", r"Liq. water pot. temp. flux $\mathrm{\overline{u'\theta_l'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{u'\theta_v'}}$", r"Virtual pot. temp. flux $\mathrm{\overline{u'\theta_v'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                # buoyancy sub-terms for parameterization in vpwp budget
                [r"$\mathrm{\overline{v'r_c'}}$", r"Liquid water flux $\mathrm{\overline{v'r_c'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'r_t'}}$", r"Total water flux $\mathrm{\overline{v'r_t'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'\theta_l'}}$", r"Liq. water pot. temp. flux $\mathrm{\overline{v'\theta_l'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'\theta_v'}}$", r"Virtual pot.temp. flux $\mathrm{\overline{v'\theta_v'}\ \left[\frac{m^2}{s^2}\right]}$"],\
                # Conditional mean wind speeds in clouds
                [r"$\mathrm{\bar{u}_{cld}}$", r"(In-cloud) mean wind $\mathrm{u_{cld}\ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\bar{v}_{cld}}$", r"(In-cloud) mean wind $\mathrm{v_{cld}\ \ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\bar{w}_{cld}}$", r"(In-cloud) mean wind $\mathrm{w_{cld}\ \ \left[\frac{m}{s}\right]}$"],\
                # Conditional comparison plots
                [r"$\mathrm{\bar{u}_{cld}, \bar{u}_{env}}$", r"Cloud conditional mean wind $\mathrm{\bar{u}*\ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\bar{v}_{cld}, \bar{v}_{env}}$", r"Cloud conditional mean wind $\mathrm{\bar{v}*\ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\bar{w}_{cld}, \bar{w}_{env}}$", r"Cloud conditional mean wind $\mathrm{\bar{w}*\ \left[\frac{m}{s}\right]}$"],\
                [r"$\mathrm{\overline{u'w'}_{cld}, \overline{u'w'}_{env}}$", r"Cloud conditional flux $\mathrm{\overline{u'w'}*\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\overline{v'w'}_{cld}, \overline{v'w'}_{env}}$", r"Cloud conditional flux $\mathrm{\overline{v'w'}*\ \left[\frac{m^2}{s^2}\right]}$"],\
                [r"$\mathrm{\bar{\theta}_v^{cld}, \bar{\theta}_v^{env}}$", r"Cloud cond. virt. pot. temp. $\mathrm{\bar{\theta}_v*\ \left[K\right]}$"],\
                [r"$\mathrm{\overline{s_L'w'}_{cld}, \overline{s_L'w'}_{env}}$", r"Cloud conditional flux $\mathrm{\overline{s_L'w'}*\ \left[\frac{K\,m}{s}\right]}$"],\
                [r"$\mathrm{\overline{r_t'w'}_{cld}, \overline{r_t'w'}_{env}}$", r"Cloud conditional flux $\mathrm{\overline{r_t'w'}*\ \left[\frac{g\,m}{kg\,s}\right]}$"],\
            ]

# lines of each plot:
# variable name within python, shall this variable be plotted?, variable name in SAM output, conversion

thetal = [\
         # variables of thetal
         ['THETAL', False, 'THETAL', 1., 0],\
         ['THETA', False, 'THETA', 1., 0],\
         ['TABS', False, 'TABS', 1., 0],\
         ['QI', False, 'QI', 1./KG, 0],\
         ['THETAL', True, 'THETAL + 2500.4 * (THETA/TABS) * QI', 1., 0],\
        ]        

rt = [\
         # variables of rt
         ['QI', False, 'QI', 1., 0],\
         ['QT', False, 'QT', 1., 0],\
         ['RT', True, '(QT-QI)', 1./KG, 0],\
        ]
        
thetalflux = [\
         # variables of thetalflux
         ['TLFLUX', False, 'TLFLUX', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['WPTHLP_SGS', False, 'WPTHLP_SGS', 1., 0],\
         ['THETALFLUX', True, '((TLFLUX) / (RHO * 1004.)) + WPTHLP_SGS', 1., 0],\
        ]
        
rtflux = [\
         # variables of rtflux
         ['QTFLUX', False, 'TLFLUX', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['WPRTP_SGS', False, 'WPRTP_SGS', 1., 0],\
         ['RTFLUX', True, '(TLFLUX / (RHO * 2.5104e+6)) + WPRTP_SGS', 1., 0],\
        ]
    
cloudliqfrac = [\
        # variables of cloudliqfrac
        ['cloudliq_frac_em6', True, 'cloudliq_frac_em6', 1., 0],\
        ]
        
qcl = [\
         # variables of qcl
         ['QCL', True, 'QCL', 1./KG, 0],\
        ]
        
wVar = [\
         # variables of wVar
         ['WP2_SGS', False, 'WP2_SGS', 1., 0],\
         ['W2', False, 'W2', 1., 0],\
         [r"$\mathrm{\overline{w'^2}}$", True, 'WP2_SGS + W2', 1., 0],\
        ]
        
w3 = [\
         # variables of wVar
         ['WP3_SGS', False, 'WP3_SGS', 1., 0],\
         ['W3', False, 'W3', 1., 0],\
         [r"$\mathrm{\overline{w'^3}}$", True, 'WP3_SGS + W3', 1., 0],\
        ]
        
thetalVar = [\
         # variables of thetalVar
         ['THLP2_SGS', False, 'THLP2_SGS', 1., 0],\
         ['TL2', False, 'TL2', 1., 0],\
         ['THETALVAR', True, 'THLP2_SGS + TL2', 1., 0],\
        ]
        
rtVar = [\
         # variables of rtVar
         ['RTP2_SGS', False, 'RTP2_SGS', 1., 0],\
         ['QT2', False, 'QT2', 1., 0],\
         ['RTVAR', True, '(QT2 / 1e+6) + RTP2_SGS', 1., 0],\
        ]
        
covarThetalRt = [\
         # variables of covarThetalRt
         ['CovarThetaLRT', True, 'RTPTHLP_SGS', 1., 0],\
        ]
        
wobs = [\
         # variables of wobs
         ['WOBS', True, 'WOBS', 1., 0],\
        ]
        
U = [\
         # variables of U
        [r"$\mathrm{\bar{u}}$", True, 'U', 1., 0],\
        ]
        
V = [\
         # variables of V
        [r"$\mathrm{\bar{v}}$", True, 'V', 1., 0],\
        ]
        
covarUW = [\
         # variables of covarUW (standard SAM run has no SGS output, only SB)
         ['UPWP_SGS', False, 'UPWP_SGS', 1., 0],\
         [r"$\mathrm{\overline{u'w'}}$ (subgrid)", True, 'UWSB', 1., 0],\
         ['UW', False, 'UW', 1., 0],\
         [r"$\mathrm{\overline{u'w'}}$ (resolved)", True,'UW-UWSB+UPWP_SGS',1,0],\
         [r"$\mathrm{\overline{u'w'}}$ (total)", True, 'UW+UPWP_SGS', 1., 0],\
        ]
        
covarVW = [\
         # variables of covarVW (standard SAM run has no SGS output, only SB)
         ['VPWP_SGS', False, 'VPWP_SGS', 1., 0],\
         [r"$\mathrm{\overline{v'w'}}$ (subgrid)", True, 'VWSB', 1., 0],\
         ['VW', False, 'VW', 1., 0],\
         [r"$\mathrm{\overline{v'w'}}$ (resolved)", True,'VW-VWSB',1,0],\
         [r"$\mathrm{\overline{v'w'}}$ (total)", True, 'VW + VPWP_SGS', 1., 0],\
        ]

uVar = [\
         # variables of uVar
         ['UP2_SGS', False, 'UP2_SGS', 1., 0],\
         ['U2', False, 'U2', 1., 0],\
         [r"$\mathrm{\overline{u'^2}}$", True, 'UP2_SGS + U2', 1., 0],\
        ]

vVar = [\
         # variables of vVar
         ['VP2_SGS', False, 'VP2_SGS', 1., 0],\
         ['V2', False, 'V2', 1., 0],\
         [r"$\mathrm{\overline{v'^2}}$", True, 'VP2_SGS + V2', 1., 0],\
        ]

corrUW = [\
        # variables of uVar
        ['UP2_SGS', False, 'UP2_SGS', 1., 0],\
        ['U2', False, 'U2', 1., 0],\
        # variables of wVar
        ['WP2_SGS', False, 'WP2_SGS', 1., 0],\
        ['W2', False, 'W2', 1., 0],\
        # variables of covarUW
        ['UPWP_SGS', False, 'UPWP_SGS', 1., 0],\
        ['UW', False, 'UW', 1., 0],\
        ['corr(u,w)', True,'(UW+UPWP_SGS)/(np.sqrt((U2+UP2_SGS)*(W2+WP2_SGS)+1e-4))', 1., 0],\
        ]
    
corrVW = [\
        # variables of vVar
        ['VP2_SGS', False, 'VP2_SGS', 1., 0],\
        ['V2', False, 'V2', 1., 0],\
        # variables of wVar
        ['WP2_SGS', False, 'WP2_SGS', 1., 0],\
        ['W2', False, 'W2', 1., 0],\
        # variables of covarUW
        ['VPWP_SGS', False, 'VPWP_SGS', 1., 0],\
        ['VW', False, 'VW', 1., 0],\
        ['corr(v,w)', True,'(VW+VPWP_SGS)/(np.sqrt((V2+VP2_SGS)*(W2+WP2_SGS)+1e-4))', 1., 0],\
        ]

# Rain Water Mixing Ratio

QR = [\
         # variables of QR
         ['QR', True, 'QR', 1./KG, 0],\
        ]
        
QRIP = [\
         # variables of QRIP
        ['qrainm_ip', True, 'qrainm_ip', 1., 0],\
        ]
          
QRP2 = [\
         # variables of QRP2
        ['qrainp2', True, 'qrainp2', 1., 0],\
        ]
        
QRP2_QRIP = [\
         # variables of QRP2_QRIP
         ['qrainp2_ip', False, 'qrainp2_ip', 1., 0],\
         ['qrainm_ip', False, 'qrainm_ip', 1., 0],\
         ['QRP2_QRIP', True, '(qrainp2_ip / (np.maximum(np.full(n,1e-5),qrainm_ip)**2))', 1., 0],\
        ]
        
#Rain Drop Number Concentration

Nrm = [\
         # variables of Nrm
         ['NR', False, 'NR', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['NRM', True, '(NR * 1e+6) / RHO', 1., 0],\
        ]
             
Nrm_IP = [\
         # variables of Nrm_IP
        ['nrainm_ip', True, 'nrainm_ip', 1., 0],\
        ]
        
Nrp2 = [\
         # variables of Nrp2
        ['nrainp2', True, 'nrainp2', 1., 0],\
        ]
        
Nrp2_NrmIP = [\
         # variables of Nrp2_NrmIP
         ['nrainp2_ip', False, 'nrainp2_ip', 1., 0],\
         ['nrainm_ip', False, 'nrainm_ip', 1., 0],\
         ['Nrp2_NrmIP', True, '(nrainp2_ip / (np.maximum(np.full(n,1e-5),nrainm_ip)**2))', 1., 0],\
        ]
        
#Cloud Droplet Number Concentration

Ncm = [\
         # variables of Ncm
         ['NC', False, 'NC', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['NCM', True, '(NC * 1e+6) / RHO', 1., 0],\
        ]
             
Ncm_IP = [\
         # variables of Ncm_IP
         ['ncloudliqm_ip', True, 'ncloudliqm_ip', 1., 0],\
        ]
        
Ncp2 = [\
         # variables of Ncp2
         ['Ncp2', True, 'Ncp2', 1., 0],\
        ]
        
Ncp2_NcmIP = [\
         # variables of Ncp2_NcmIP
         ['ncloudliqp2_ip', False, 'ncloudliqp2_ip', 1., 0],\
         ['ncloudliqm_ip', False, 'ncloudliqm_ip', 1., 0],\
         ['Ncp2_NcmIP', True, '(ncloudliqp2_ip / (np.maximum(np.full(n,1e-5),ncloudliqm_ip)**2))', 1., 0],\
        ]
        
#Graupel Number Concentration

Ngm = [\
         # variables of Ngm
         ['NG', False, 'NG', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['NGM', True, '(NG * 1e+6) / RHO', 1., 0],\
        ]
             
Ngm_IP = [\
         # variables of Ngm_IP
         ['ngraupelm_ip', True, 'ngraupelm_ip', 1., 0],\
        ]
        
Ngp2 = [\
         # variables of Ngp2
         ['ngraupelp2', True, 'ngraupelp2', 1., 0],\
        ]
        
Ngp2_NgmIP = [\
         # variables of Ngp2_NgmIP
         ['ngraupelp2_ip', False, 'ngraupelp2_ip', 1., 0],\
         ['ngraupelm_ip', False, 'ngraupelm_ip', 1., 0],\
         ['Ngp2_NgmIP', True, '(ngraupelp2_ip / (np.maximum(np.full(n,1e-5),ngraupelm_ip)**2))', 1., 0],\
        ]
        
#Graupel Mixing Ratio

Qgm = [\
         # variables of Qgm
         ['QG', True, 'QG', 1./KG, 0],\
        ]
             
Qgm_IP = [\
         # variables of Qgm_IP
         ['qgraupelm_ip', True, 'qgraupelm_ip', 1., 0],\
        ]
        
Qgp2 = [\
         # variables of Qgp2
        ['qgraupelp2', True, 'qgraupelp2', 1., 0],\
        ]
        
Qgp2_QgmIP = [\
         # variables of Qgp2_QgmIP
         ['qgraupelp2_ip', False, 'qgraupelp2_ip', 1., 0],\
         ['qgraupelm_ip', False, 'qgraupelm_ip', 1., 0],\
         ['Qgp2_QgmIP', True, '(qgraupelp2_ip / (np.maximum(np.full(n,1e-5),qgraupelm_ip)**2))', 1., 0],\
        ]
        
#Cloud Ice Mixing Ratio

# Note: redundant, could not find variable in sam
Qim = [\
         # variables of Qim
         ['QG', True, 'QG', 1./KG, 0],\
        ]
             
Qim_IP = [\
         # variables of Qim_IP
        ['qcloudicem_ip', True, 'qcloudicem_ip', 1., 0],\
        ]
        
Qip2 = [\
         # variables of Qip2
        ['qcloudicep2', True, 'qcloudicep2', 1., 0],\
        ]
        
Qip2_QimIP = [\
         # variables of Qip2_QimIP
         ['qcloudicep2_ip', False, 'qcloudicep2_ip', 1., 0],\
         ['qcloudicem_ip', False, 'qcloudicem_ip', 1., 0],\
         ['Qip2_QimIP', True, '(qcloudicep2_ip / (np.maximum(np.full(n,1e-5),qcloudicem_ip)**2))', 1., 0],\
        ]
        
#Cloud Ice Number Concentration

Nim = [\
         # variables of Nim
         ['NI', False, 'NI', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['NIM', True, '(NI * 1e+6) / RHO', 1., 0],\
        ]
             
Nim_IP = [\
         # variables of Nim_IP
        ['ncloudicem_ip', True, 'ncloudicem_ip', 1., 0],\
        ]
        
Nip2 = [\
         # variables of Nip2
        ['ncloudicep2', True, 'ncloudicep2', 1., 0],\
        ]
        
Nip2_NimIP = [\
         # variables of Nip2_NimIP
         ['ncloudicep2_ip', False, 'ncloudicep2_ip', 1., 0],\
         ['ncloudicem_ip', False, 'ncloudicem_ip', 1., 0],\
         ['Nip2_NimIP', True, '(ncloudicep2_ip / (np.maximum(np.full(n,1e-5),ncloudicem_ip)**2))', 1., 0],\
        ]

#Snow Mixing Ratio

Qsm = [\
         # variables of Qsm
         ['QSM', True, 'QS', 1./KG, 0],\
        ]
             
Qsm_IP = [\
         # variables of Qsm_IP
        ['qsnowm_ip', True, 'qsnowm_ip', 1., 0],\
        ]
        
Qsp2 = [\
         # variables of Qsp2
        ['qsnowp2', True, 'qsnowp2', 1., 0],\
        ]
        
Qsp2_QsmIP = [\
         # variables of Qsp2_QsmIP
         ['qsnowp2_ip', False, 'qsnowp2_ip', 1., 0],\
         ['qsnowm', False, 'qsnowm', 1., 0],\
         ['Qsp2_QsmIP', True, '(qsnowp2_ip / (np.maximum(np.full(n,1e-5),qsnowm)**2))', 1., 0],\
        ]
        
#Snow Number Concentration

Nsm = [\
         # variables of Nsm
         ['NS', False, 'NS', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['NSM', True, '(NS * 1e+6) / RHO', 1., 0],\
        ]
             
Nsm_IP = [\
         # variables of Nsm_IP
        ['nsnowm_ip', True, 'nsnowm_ip', 1., 0],\
        ]
        
Nsp2 = [\
         # variables of Nsp2
        ['nsnowp2', True, 'nsnowp2', 1., 0],\
        ]
        
Nsp2_NsmIP = [\
         # variables of Nsp2_NsmIP
         ['nsnowp2_ip', False, 'nsnowp2_ip', 1., 0],\
         ['nsnowm_ip', False, 'nsnowm_ip', 1., 0],\
         ['Nsp2_NsmIP', True, '(nsnowp2_ip / (np.maximum(np.full(n,1e-5),nsnowm_ip)**2))', 1., 0],\
        ]
        

MicroFractions = [\
         # variables of MicroFractions
         ['Cloud_liq', True, 'cloudliq_frac_em6', 1., 0],\
         ['Rain', True, 'rain_frac_em6', 1., 0],\
         ['Cloud_ice', True, 'cloudice_frac_em6', 1., 0],\
         ['Snow', True, 'snow_frac_em6', 1., 0],\
         ['Graupel', True, 'graupel_frac_em6', 1., 0],\
        ]
        
Buoy_Flux = [\
         # variables of Buoy_Flux
         ['TVFLUX', False, 'TVFLUX', 1., 0],\
         ['RHO', False, 'RHO', 1., 0],\
         ['Buoy_Flux', True, 'TVFLUX / (RHO * 1004)', 1., 0],\
        ]
         
lwp = [\
         # variables of lwp
         ['CWP', True, 'CWP', 1./KG, 0],\
        ]
        
PREC = [\
         # variables of lwp
         ['PREC', True, 'PREC', 1., 0],\
        ]
           
WP2_W2 = [\
         # variables of WP2_W2
         ['NS', True, 'NS', 0., 0],\
        ]

IWP = [\
         # variables of IWP
         ['IWP', True, 'IWP', 1./KG, 0],\
        ]
        
SWP = [\
         # variables of SWP
         ['SWP', True, 'SWP', 1./KG, 0],\
        ]

uprcp = [\
        [r"$\mathrm{\overline{u'r_c'}}$", True, 'UPRCP', 1., 0],\
        ]

uprtp = [\
        [r"$\mathrm{\overline{u'r_t'}}$", True, 'UPRTP', 1., 0],\
        ]

upthlp = [\
        [r"$\mathrm{\overline{u'\theta_l'}}$", True, 'UPTHLP', 1., 0],\
        ]
    
upthvp = [\
        [r"$\mathrm{\overline{u'\theta_v'}}$", True, 'UPTHVP', 1., 0],\
        ]
    
vprcp = [\
        [r"$\mathrm{\overline{v'r_c'}}$", True, 'VPRCP', 1., 0],\
        ]

vprtp = [\
        [r"$\mathrm{\overline{v'r_t'}}$", True, 'VPRTP', 1., 0],\
        ]

vpthlp = [\
        [r"$\mathrm{\overline{v'\theta_l'}}$", True, 'VPTHLP', 1., 0],\
        ]

vpthvp = [\
        [r"$\mathrm{\overline{v'\theta_v'}}$", True, 'VPTHVP', 1., 0],\
        ]

# Conditional plots
ucld = [\
    ['In-cloud mean of u', True, 'UCLD', 1., 0 ],\
    ['Total mean of u', True, 'U', 1., 0 ],\
    ]


vcld = [\
    ['In-cloud mean of v', True, 'VCLD', 1., 0 ],\
    ['Total mean of v', True, 'V', 1., 0 ],\
    ]

wcld = [\
    ['In-cloud mean of w', True, 'WCLD', 1., 0 ],\
    ['Total mean of w', True, 'WM', 1., 0 ],\
    ]

# Cnditional comparison plots

ucomp = [\
    ['In-cloud mean of u', True, 'UCLD', 1., 0 ],\
    ['Total mean of u', False, 'U', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    ['Environmental mean of u', True, '(U - CLD*UCLD)/(1-CLD)', 1., 0],\
    ]

vcomp = [\
    ['In-cloud mean of v', True, 'VCLD', 1., 0 ],\
    ['V', False, 'V', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    ['Environmental mean of v', True, '(V - CLD*VCLD)/(1-CLD)', 1., 0],\
    ]

wcomp = [\
    ['In-cloud mean of w', True, 'WCLD', 1., 0 ],\
    ['WM', False, 'WM', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    ['Environmental mean of w', True, '(WM - CLD*WCLD)/(1-CLD)', 1., 0],\
    ]

uwcomp = [\
    [r"In-cloud $\mathrm{\overline{u'w'}}$", True, 'UWCLD', 1., 0 ],\
    ['UW', False, 'UW', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    [r"Environmental $\mathrm{\overline{u'w'}}$", True, '(UW - CLD*UWCLD)/(1-CLD)', 1., 0],\
    ]

vwcomp = [\
    [r"In-cloud $\mathrm{\overline{v'w'}}$", True, 'VWCLD', 1., 0 ],\
    ['VW', False, 'VW', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    [r"Environmental $\mathrm{\overline{v'w'}}$", True, '(VW - CLD*VWCLD)/(1-CLD)', 1., 0],\
    ]

tvcomp = [\
    [r'In-cloud $\mathrm{\bar{\theta}_v}$', True, 'TVCLD', 1., 0 ],\
    ['THETAV', False, 'THETAV', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    [r'Environmental $\mathrm{\bar{\theta}_v}$', True, '(THETAV - CLD*TVCLD)/(1-CLD)', 1., 0],\
    ]

tlcomp = [\
    [r"In-cloud $\mathrm{\overline{s'_Lw'}}$", True, 'TLWCLD', 1., 0 ],\
    ['TLFLUX', False, 'TLFLUX', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    [r"Environmental $\mathrm{\overline{s'_Lw'}}$", True, '(TLFLUX - CLD*TLWCLD)/(1-CLD)', 1., 0],\
    ]

qtwcomp = [\
    [r"In-cloud $\mathrm{\overline{r_t'w'}}$", True, 'QTWCLD', 1., 0 ],\
    ['WPRTP', False, 'WPRTP', 1., 0 ],\
    ['CLD', False, 'CLD', 1., 0],\
    [r"Environmental $\mathrm{\overline{r_t'w'}}$", True, '(WPRTP - CLD*QTWCLD)/(1-CLD)', 1., 0],\
    ]

lines = [thetal, rt, thetalflux, rtflux, cloudliqfrac, qcl, \
        wVar, w3, thetalVar, rtVar, covarThetalRt, wobs, U, V, covarUW, covarVW, uVar, vVar, corrUW, corrVW, \
        QR, QRIP, QRP2, QRP2_QRIP, \
        Nrm, Nrm_IP, Nrp2, Nrp2_NrmIP, \
        Ncm, Ncm_IP, Ncp2, Ncp2_NcmIP, \
        Ngm, Ngm_IP, Ngp2, Ngp2_NgmIP, \
        Qgm, Qgm_IP, Qgp2, Qgp2_QgmIP, \
        Qim, Qim_IP, Qip2, Qip2_QimIP, \
        Nim, Nim_IP, Nip2, Nip2_NimIP, \
        Qsm, Qsm_IP, Qsp2, Qsp2_QsmIP, \
        Nsm, Nsm_IP, Nsp2, Nsp2_NsmIP, \
        MicroFractions, Buoy_Flux, \
        uprcp, uprtp, upthlp, upthvp, \
        vprcp, vprtp, vpthlp, vpthvp,\
        ucld, vcld, wcld,\
        ucomp, vcomp, wcomp, uwcomp, vwcomp, tvcomp, tlcomp, qtwcomp,\
        ]
