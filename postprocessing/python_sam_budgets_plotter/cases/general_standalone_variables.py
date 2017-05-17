#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000.
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
sortPlots = ['theta_l', 'r_t', 'theta_l_flux', 'r_t_flux', 'cloudliq_frac', 'r_c', 'w_var', 'w3', 'theta_l_var', 'r_t_var', 'covar_thetal_rt', 'wobs', 'U', 'V', 'covar_uw', 'covar_vw', 'u_var', 'v_var',\
                'QR', 'QR_IP', 'QRP2', 'QRP2_QRIP', \
                'Nrm', 'Nrm_IP', 'Nrp2', 'Nrp2_NrmIP', \
                'Ncm', 'Ncm_IP', 'Ncp2', 'Ncp2_NcmIP', \
                'Ngm', 'Ngm_IP', 'Ngp2', 'Ngp2_NgmIP', \
                'Qgm', 'Qgm_IP', 'Qgp2', 'Qgp2_QgmIP', \
                'Qim', 'Qim_IP', 'Qip2', 'Qip2_QimIP', \
                'Nim', 'Nim_IP', 'Nip2', 'Nip2_NimIP', \
                'Qsm', 'Qsm_IP', 'Qsp2', 'Qsp2_QsmIP', \
                'Nsm', 'Nsm_IP', 'Nsp2', 'Nsp2_NsmIP', \
                'MicroFractions', 'Buoy_Flux' \
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
                ['Variance of w', 'wp2 / w2 [m^2/s^2]'],\
                ['Third-order Moment of w', 'wp3 / w3 [m^3/s^3]'],\
                ['Variance of theta_l', 'thlp2 / tl2 [K^2]'],\
                ['Variance of r_t', 'rtp2 / qtp2 [(kg/kg)^2]'],\
                ['Covariance of r_t & theta_l', 'rtpthlp [(kg/kg) K]'],\
                ['Vertical Wind Component, w (subsidence)', 'wobs [m/s]'],\
                ['Zonal Wind Component, u', 'um / u [m/s]'],\
                ['Meridonal Wind Component, v', 'vm / v [m/s]'],\
                ['Covariance of u & w', 'upwp / uw [m^2/s^2]'],\
                ['Covariance of v & w', 'vpwp / vw [m^2/s^2]'],\
                ['Variance of u wind', 'up2 / u2 [m^2/s^2]'],\
                ['Variance of v wind', 'vp2 / v2 [m^2/s^2]'],\
                # Rain Water Mixing Ratio
                ['Rain Water Mixing Ratio', 'qrm [kg/kg]'],\
                ['Rain Water Mixing Ratio in Rain', 'qrm_ip [kg/kg]'],\
                ['Domain-wide Variance of Rain Water Mixing Ratio', 'qrp2 [(kg/kg)^2]'],\
                ['Within-rain Variance of Rain Water Mixing Ratio', 'qrp2_ip / qrm_ip^2 [-]'],\
                #Rain Drop Number Concentration
                ['Rain Drop Concentration', 'Nrm [num/kg]'],\
                ['Rain Drop Concentration in Rain', 'Nrm_ip [num/kg]'],\
                ['Domain-wide Variance of Rain Drop Concentration', 'Nrp2 [(num/kg)^2]'],\
                ['Within-rain Variance of Rain Drop Concentration', 'Nrp2_ip / Nrm_ip^2 [-]'],\
                #Cloud Droplet Number Concentration
                ['Cloud Droplet Number Concentration', 'Ncm [num/kg]'],\
                ['Cloud Droplet Number Concentration in Cloud', 'Ncm_ip [num/kg]'],\
                ['Domain-wide Variance of Cloud Droplet Number Concentration', 'Ncp2 [(#/kg)^2]'],\
                ['Within-cloud Variance of Cloud Droplet Number Concentration', 'Ncp2_ip / Ncm_ip^2 [-]'],\
                #Graupel Number Concentration
                ['Graupel Number Concentration', 'Ngm [kg/kg]'],\
                ['Graupel Number Concentration in Graupel', 'Ngm_ip [num/kg]'],\
                ['Domain-wide Variance of Graupel Number ConcentratMicroFractionsion', 'Ngp2 [(kg/kg)^2]'],\
                ['Within-graupel Variance of Graupel Number Concentration', 'Ngp2_ip / Ngm_ip^2 [-]'],\
                #Graupel Mixing Ratio
                ['Graupel Mixing Ratio', 'qgm [kg/kg]'],\
                ['Graupel Mixing Ratio in Graupel', 'qgm_ip [kg/kg]'],\
                ['Domain-wide Variance of Graupel Mixing Ratio', 'qgp2 [(kg/kg)^2]'],\
                ['Within-graupel Variance of Graupel Mixing Ratio', 'qgp2_ip / qgm_ip^2 [-]'],\
                #Cloud Ice Mixing Ratio
                ['Cloud Ice Mixing Ratio', 'qim [kg/kg]'],\
                ['Cloud Ice Mixing Ratio in Cloud Ice', 'qim_ip [kg/kg]'],\
                ['Domain-wide Variance of Cloud Ice Mixing Ratio', 'qip2 [(kg/kg)^2]'],\
                ['Within-cloud-ice Variance of Cloud Ice  Mixing Ratio', 'qip2_ip / qim_ip^2 [-]'],\
                #Cloud Ice Number Concentration
                ['Cloud Ice Concentration', 'Nim [num/kg]'],\
                ['Cloud Ice Number Concentration in Cloud Ice', 'Ni_ip [num/kg]'],\
                ['Domain-wide Variance of Cloud Ice Number Concentration', 'Nip2 [(num/kg)^2]'],\
                ['Within-cloud-ice Variance of Cloud Ice Number Concentration', 'Nip2_ip / Nim_ip^2 [-]'],\
                #Snow Mixing Ratio
                ['Snow Mixing Ratio ', 'qsm [kg/kg]'],\
                ['Snow Mixing Ratio in Snow', 'qsm_ip [kg/kg]'],\
                ['Domain-wide Variance of Snow Mixing Ratio', 'qsp2 [(kg/kg)^2]'],\
                ['Within-snow Variance of Snow Mixing Ratio ', 'qsp2_ip / qsm_ip^2 [-]'],\
                #Snow Number Concentration
                ['Snow Number Concentration', 'Nsm [num/kg]'],\
                ['Snow Number Concentration in Snow', 'Nsm_ip [num/kg]'],\
                ['Domain-wide Variance of Snow Number Concentration', 'Nsp2 [(#/kg)^2]'],\
                ['Within-snow Variance of Snow Number Concentration', 'Nsp2_ip / Nsm_ip^2 [-]'],\
                ['Micro Fractions', '[%/100]'],\
                ['Buoyancy flux', 'wpthvp / tlflux [K m/s]'],\
                #['Liquid Water Path', 'lwp [kg/m^2]'],\
                #['Surface rainfall rate', 'rain_rate_sfc[mm/day]'],\
                #['Density-Weighted Vertically Averaged wp2', 'wp2 / w2 [m^2/s^2]'],\
                #['Cloud Ice Water Path', 'iwp [kg/m^2]'],\
                #['Snow Water Path', 'swp [kg/m^2]'],\
            ]

# lines of each plot:
# variable name within python, shall this variable be plotted?, variable name in SAM output, conversion

thetal = [\
         # variables of thetal
         ['THETAL', False, 'THETAL', 1.],\
         ['THETA', False, 'THETA', 1.],\
         ['TABS', False, 'TABS', 1.],\
         ['QI', False, 'QI', 1./KG],\
         ['', True, 'THETAL + 2500.4 * (THETA/TABS) * QI', 1.],\
        ]        

rt = [\
         # variables of rt
         ['QI', False, 'QI', 1.],\
         ['QT', False, 'QT', 1.],\
         ['', True, '(QT-QI)', 1./KG],\
        ]
        
thetalflux = [\
         # variables of thetalflux
         ['TLFLUX', False, 'TLFLUX', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['WPTHLP_SGS', False, 'WPTHLP_SGS', 1.],\
         ['', True, '((TLFLUX) / (RHO * 1004.)) + WPTHLP_SGS', 1.],\
        ]
        
rtflux = [\
         # variables of rtflux
         ['QTFLUX', False, 'TLFLUX', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['WPRTP_SGS', False, 'WPRTP_SGS', 1.],\
         ['', True, '(QTFLUX / (RHO * 2.5104e+6)) + WPRTP_SGS', 1.],\
        ]
    
cloudliqfrac = [\
         # variables of cloudliqfrac
         ['', True, 'cloudliq_frac_em6', 1.],\
        ]
        
qcl = [\
         # variables of qcl
         ['', True, 'QCL', 1./KG],\
        ]
        
wVar = [\
         # variables of wVar
         ['WP2_SGS', False, 'WP2_SGS', 1.],\
         ['W2', False, 'W2', 1.],\
         ['', True, 'WP2_SGS + W2', 1.],\
        ]
        
w3 = [\
         # variables of wVar
         ['WP3_SGS', False, 'WP3_SGS', 1.],\
         ['W3', False, 'W3', 1.],\
         ['', True, 'WP3_SGS + W3', 1.],\
        ]
        
thetalVar = [\
         # variables of thetalVar
         ['THLP2_SGS', False, 'THLP2_SGS', 1.],\
         ['TL2', False, 'TL2', 1.],\
         ['', True, 'THLP2_SGS + TL2', 1.],\
        ]
        
rtVar = [\
         # variables of rtVar
         ['RTP2_SGS', False, 'RTP2_SGS', 1.],\
         ['QT2', False, 'QT2', 1.],\
         ['', True, '(QT2 / 1e+6) + RTP2_SGS', 1.],\
        ]
        
covarThetalRt = [\
         # variables of covarThetalRt
         ['', True, 'RTPTHLP_SGS', 1.],\
        ]
        
wobs = [\
         # variables of wobs
         ['', True, 'WOBS', 1.],\
        ]
        
U = [\
         # variables of U
         ['', True, 'U', 1.],\
        ]
        
V = [\
         # variables of V
         ['', True, 'V', 1.],\
        ]
        
covarUV = [\
         # variables of covarUV
         ['UPWP_SGS', False, 'UPWP_SGS', 1.],\
         ['UW', False, 'UW', 1.],\
         ['', True, 'UW + UPWP_SGS', 1.],\
        ]
        
covarVW = [\
         # variables of covarVW
         ['VPWP_SGS', False, 'VPWP_SGS', 1.],\
         ['VW', False, 'VW', 1.],\
         ['', True, 'VW + VPWP_SGS', 1.],\
        ]
        
uVar = [\
         # variables of uVar
         ['UP2_SGS', False, 'UP2_SGS', 1.],\
         ['U2', False, 'U2', 1.],\
         ['', True, 'UP2_SGS + U2', 1.],\
        ]
        
vVar = [\
         # variables of vVar
         ['VP2_SGS', False, 'VP2_SGS', 1.],\
         ['V2', False, 'V2', 1.],\
         ['', True, 'VP2_SGS + V2', 1.],\
        ]
        
# Rain Water Mixing Ratio

QR = [\
         # variables of QR
         ['', True, 'QR', 1./KG],\
        ]
        
QRIP = [\
         # variables of QRIP
         ['', True, 'qrainm_ip', 1.],\
        ]
          
QRP2 = [\
         # variables of QRP2
         ['', True, 'qrainp2', 1.],\
        ]
        
QRP2_QRIP = [\
         # variables of QRP2_QRIP
         ['qrainp2_ip', False, 'qrainp2_ip', 1.],\
         ['qrainm_ip', False, 'qrainm_ip', 1.],\
         ['', True, '(qrainp2_ip / (np.maximum(np.full(n,1e-5),qrainm_ip)**2))', 1.],\
        ]
        
#Rain Drop Number Concentration

Nrm = [\
         # variables of Nrm
         ['NR', False, 'NR', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, '(NR * 1e+6) / RHO', 1.],\
        ]
             
Nrm_IP = [\
         # variables of Nrm_IP
         ['', True, 'nrainm_ip', 1.],\
        ]
        
Nrp2 = [\
         # variables of Nrp2
         ['', True, 'nrainp2', 1.],\
        ]
        
Nrp2_NrmIP = [\
         # variables of Nrp2_NrmIP
         ['nrainp2_ip', False, 'nrainp2_ip', 1.],\
         ['nrainm_ip', False, 'nrainm_ip', 1.],\
         ['', True, '(nrainp2_ip / (np.maximum(np.full(n,1e-5),nrainm_ip)**2))', 1.],\
        ]
        
#Cloud Droplet Number Concentration

Ncm = [\
         # variables of Ncm
         ['NC', False, 'NC', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, '(NC * 1e+6) / RHO', 1.],\
        ]
             
Ncm_IP = [\
         # variables of Ncm_IP
         ['', True, 'ncloudliqm_ip', 1.],\
        ]
        
Ncp2 = [\
         # variables of Ncp2
         ['', True, 'nrainp2', 1.],\
        ]
        
Ncp2_NcmIP = [\
         # variables of Ncp2_NcmIP
         ['ncloudliqp2_ip', False, 'ncloudliqp2_ip', 1.],\
         ['ncloudliqm_ip', False, 'ncloudliqm_ip', 1.],\
         ['', True, '(ncloudliqp2_ip / (np.maximum(np.full(n,1e-5),ncloudliqm_ip)**2))', 1.],\
        ]
        
#Graupel Number Concentration

Ngm = [\
         # variables of Ngm
         ['NG', False, 'NG', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, '(NG * 1e+6) / RHO', 1.],\
        ]
             
Ngm_IP = [\
         # variables of Ngm_IP
         ['', True, 'ngraupelm_ip', 1.],\
        ]
        
Ngp2 = [\
         # variables of Ngp2
         ['', True, 'ngraupelp2', 1.],\
        ]
        
Ngp2_NgmIP = [\
         # variables of Ngp2_NgmIP
         ['ngraupelp2_ip', False, 'ngraupelp2_ip', 1.],\
         ['ngraupelm_ip', False, 'ngraupelm_ip', 1.],\
         ['', True, '(ngraupelp2_ip / (np.maximum(np.full(n,1e-5),ngraupelm_ip)**2))', 1.],\
        ]
        
#Graupel Mixing Ratio

Qgm = [\
         # variables of Qgm
         ['', True, 'QG', 1./KG],\
        ]
             
Qgm_IP = [\
         # variables of Qgm_IP
         ['', True, 'qgraupelm_ip', 1.],\
        ]
        
Qgp2 = [\
         # variables of Qgp2
         ['', True, 'qgraupelp2', 1.],\
        ]
        
Qgp2_QgmIP = [\
         # variables of Qgp2_QgmIP
         ['qgraupelp2_ip', False, 'qgraupelp2_ip', 1.],\
         ['qgraupelm_ip', False, 'qgraupelm_ip', 1.],\
         ['', True, '(qgraupelp2_ip / (np.maximum(np.full(n,1e-5),qgraupelm_ip)**2))', 1.],\
        ]
        
#Cloud Ice Mixing Ratio

Qim = [\
         # variables of Qim
         ['', True, 'QG', 1./KG],\
        ]
             
Qim_IP = [\
         # variables of Qim_IP
         ['', True, 'qcloudicem_ip', 1.],\
        ]
        
Qip2 = [\
         # variables of Qip2
         ['', True, 'qcloudicep2', 1.],\
        ]
        
Qip2_QimIP = [\
         # variables of Qip2_QimIP
         ['qcloudicep2_ip', False, 'qcloudicep2_ip', 1.],\
         ['qcloudicem_ip', False, 'qcloudicem_ip', 1.],\
         ['', True, '(qcloudicep2_ip / (np.maximum(np.full(n,1e-5),qcloudicem_ip)**2))', 1.],\
        ]
        
#Cloud Ice Number Concentration

Nim = [\
         # variables of Nim
         ['NI', False, 'NI', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, '(NI * 1e+6) / RHO', 1.],\
        ]
             
Nim_IP = [\
         # variables of Nim_IP
         ['', True, 'ncloudicem_ip', 1.],\
        ]
        
Nip2 = [\
         # variables of Nip2
         ['', True, 'ncloudicep2', 1.],\
        ]
        
Nip2_NimIP = [\
         # variables of Nip2_NimIP
         ['ncloudicep2_ip', False, 'ncloudicep2_ip', 1.],\
         ['ncloudicem_ip', False, 'ncloudicem_ip', 1.],\
         ['', True, '(ncloudicep2_ip / (np.maximum(np.full(n,1e-5),ncloudicem_ip)**2))', 1.],\
        ]

#Snow Mixing Ratio

Qsm = [\
         # variables of Qsm
         ['', True, 'QS', 1./KG],\
        ]
             
Qsm_IP = [\
         # variables of Qsm_IP
         ['', True, 'qsnowm_ip', 1.],\
        ]
        
Qsp2 = [\
         # variables of Qsp2
         ['', True, 'qsnowp2', 1.],\
        ]
        
Qsp2_QsmIP = [\
         # variables of Qsp2_QsmIP
         ['qsnowp2_ip', False, 'qsnowp2_ip', 1.],\
         ['qsnowm', False, 'qsnowm', 1.],\
         ['', True, '(qsnowp2_ip / (np.maximum(np.full(n,1e-5),qsnowm)**2))', 1.],\
        ]
        
#Snow Number Concentration

Nsm = [\
         # variables of Nsm
         ['NS', False, 'NS', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, '(NS * 1e+6) / RHO', 1.],\
        ]
             
Nsm_IP = [\
         # variables of Nsm_IP
         ['', True, 'nsnowm_ip', 1.],\
        ]
        
Nsp2 = [\
         # variables of Nsp2
         ['', True, 'nsnowp2', 1.],\
        ]
        
Nsp2_NsmIP = [\
         # variables of Nsp2_NsmIP
         ['nsnowp2_ip', False, 'nsnowp2_ip', 1.],\
         ['nsnowm_ip', False, 'nsnowm_ip', 1.],\
         ['', True, '(nsnowp2_ip / (np.maximum(np.full(n,1e-5),nsnowm_ip)**2))', 1.],\
        ]
        

MicroFractions = [\
         # variables of MicroFractions
         ['Cloud_liq', True, 'cloudliq_frac_em6', 1.],\
         ['Rain', True, 'rain_frac_em6', 1.],\
         ['Cloud_ice', True, 'cloudice_frac_em6', 1.],\
         ['Snow', True, 'snow_frac_em6', 1.],\
         ['Graupel', True, 'graupel_frac_em6', 1.],\
        ]
        
Buoy_Flux = [\
         # variables of Buoy_Flux
         ['TVFLUX', False, 'TVFLUX', 1.],\
         ['RHO', False, 'RHO', 1.],\
         ['', True, 'TVFLUX / (RHO * 1004)', 1.],\
        ]
         
lwp = [\
         # variables of lwp
         ['', True, 'CWP', 1./KG],\
        ]
        
PREC = [\
         # variables of lwp
         ['', True, 'PREC', 1.],\
        ]
           
WP2_W2 = [\
         # variables of WP2_W2
         ['', True, 'NS', 0.],\
        ]

IWP = [\
         # variables of IWP
         ['', True, 'IWP', 1./KG],\
        ]
        
SWP = [\
         # variables of SWP
         ['', True, 'SWP', 1./KG],\
        ]


lines = [thetal, rt, thetalflux, rtflux, cloudliqfrac, qcl, wVar, w3, thetalVar, rtVar, covarThetalRt, wobs, U, V, covarUV, covarVW, uVar, vVar,\
        QR, QRIP, QRP2, QRP2_QRIP, \
        Nrm, Nrm_IP, Nrp2, Nrp2_NrmIP, \
        Ncm, Ncm_IP, Ncp2, Ncp2_NcmIP, \
        Ngm, Ngm_IP, Ngp2, Ngp2_NgmIP, \
        Qgm, Qgm_IP, Qgp2, Qgp2_QgmIP, \
        Qim, Qim_IP, Qip2, Qip2_QimIP, \
        Nim, Nim_IP, Nip2, Nip2_NimIP, \
        Qsm, Qsm_IP, Qsp2, Qsp2_QsmIP, \
        Nsm, Nsm_IP, Nsp2, Nsp2_NsmIP, \
        MicroFractions, Buoy_Flux \
        ]