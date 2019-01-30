#-------------------------------------------------------------------------------
#   C O N S T A N T S
#-------------------------------------------------------------------------------
DAY = 24
HOUR = 3600
KG = 1000
g_per_second_to_kg_per_day = 1. / (DAY * HOUR * KG)
kg_per_second_to_kg_per_day = 1. / (DAY * HOUR)

#-------------------------------------------------------------------------------
# P L O T S
#-------------------------------------------------------------------------------
sortPlots = ['HL', 'QT', 'TW', 'THLW', 'QW', 'QTOGW', 'W2', 'W3', 'T2', 'THL2', 'Q2', 'QTOG2', 'QTHL', 'TKE', 'TKES', 'U2V2', 'WU','WV','U2', 'V2', 'V2_COMP', 'U2_V2_TAU', 'QTOG2_TAU', 'Q2_TAU', 'THL2_TAU', 'C2_TAU_COMP']
# settings of each plot:
# plot number, plot title, axis label
plotNames = [\
                ['HL 3D SAM Benchmark Budgets', 'HL [K/s]'],\
                ['QT 3D SAM Benchmark Budgets', 'QTO [kg/kg/s]'],\
                ['TW SAM 3D Benchmark Budgets', 'TW [(m K)/s^2]'],\
                ['THLW SAM 3D Benchmark Budgets', 'THLW [(m K)/s^2]'],\
                ['QW SAM 3D Benchmark Budgets', 'QW [(kg/kg)m/s^2]'],\
                ['QTOGW SAM 3D Benchmark Budgets', 'QTOGW [(kg/kg)m/s^2]'],\
                ['W2 SAM 3D Benchmark Budgets', 'W2 [m^2/s^3]'],\
                ['W3 SAM 3D Benchmark Budgets', 'W3 [m^3/s^4]'],\
                ['T2 3D SAM Benchmark Budgets', 'T2 [K^2/s]'],\
                ['THL2 3D SAM Benchmark Budgets', 'THL2 [K^2/s]'],\
                ['Q2 3D SAM Benchmark Budgets', 'Q2 [(kg/kg)^2/s]'],\
                ['QTOG2 3D SAM Benchmark Budgets', 'QTOG2 [(kg/kg)^2/s]'],\
                ['QTHL 3D SAM Benchmark Budgets', 'QTHL [(kg/kg)K/s]'],\
                ['TKE 3D SAM Benchmark Budgets', 'TKE [m^2/s^3]'],\
                ['TKE (SGS) 3D SAM Benchmark Budgets', 'TKE [m^2/s^3]'],\
                ['U2+V2 3D SAM Benchmark Budgets', 'U2+V2 [m^2/s^2]'],\
                ['WU 3D SAM Benchmark Budgets', 'WU [m^2/s^2]'],\
                ['WV 3D SAM Benchmark Budgets', 'WV [m^2/s^2]'],\
                ['U2 3D SAM Benchmark Budgets', 'U2 [m^2/s^2]'],\
                ['V2 3D SAM Benchmark Budgets', 'V2 [m^2/s^2]'],\
                ['V2 Comparison from U2V2-U2', 'V2 [m^2/s^2]'],\
                ['C14/TAU Calculation based on\nU2 and V2 Budgets', 'C14/TAU [1/s]'],\
                ['C2/TAU Calculation based on QTOG2 Budget', 'C2/TAU [1/s]'],\
                ['C2/TAU Calculation based on QT2 Budget', 'C2/TAU [1/s]'],\
                ['C2/TAU Calculation based on THL2 Budget', 'C2/TAU [1/s]'],\
                ['C2/TAU Comparison', 'C2/TAU [1/s]'],\
            ]

# lines of each plot:
# variable name within python, shall this variable be plotted?, variable name in SAM output, conversion

HL = [\
         # variables of HL
         ['HLADV', True, 'HLADV', kg_per_second_to_kg_per_day, 0],\
         ['HLDIFF', True, 'HLDIFF', kg_per_second_to_kg_per_day, 0],\
         ['HLRAD', True, 'HLRAD', kg_per_second_to_kg_per_day, 0],\
         ['HLLAT', True, 'HLLAT', kg_per_second_to_kg_per_day, 0],\
         ['TTEND', True, 'TTEND', kg_per_second_to_kg_per_day, 0],\
         ['T_TNDCY', True, 'T_TNDCY', 1.0, 0],\
         ['HLSTOR', True, 'HLSTOR', kg_per_second_to_kg_per_day, 0],\
         ['HL_RES', True, '((HLSTOR)+((-1)*(HLADV+HLDIFF+HLRAD+HLLAT+TTEND)))', g_per_second_to_kg_per_day, 0],\
        ]
         
QT = [\
         # variables of QTO
         ['QTADV', True, 'QTADV', g_per_second_to_kg_per_day, 0],\
         ['QTDIFF', True, 'QTDIFF', g_per_second_to_kg_per_day, 0],\
         ['QTSRC', True, 'QTSRC', g_per_second_to_kg_per_day, 0],\
         ['QTSINK', True, 'QTSINK', g_per_second_to_kg_per_day, 0],\
         ['QTEND', True, 'QTEND', g_per_second_to_kg_per_day, 0],\
         ['QTSTOR', True, 'QTSTOR', g_per_second_to_kg_per_day, 0],\
         ['QV_TNDCY', True, 'QV_TNDCY', g_per_second_to_kg_per_day, 0],\
         ['QC_TNDCY', True, 'QC_TNDCY', g_per_second_to_kg_per_day, 0],\
         ['QT_RES', True, 'QTSTOR+(-1)*(QTADV+QTDIFF+QTSRC+QTSINK+QTEND)', 1, 0],\
        ]
         
TW = [\
         # variables of TW
         ['TWGRAD', True, 'TWGRAD', 1, 0],\
         ['TWADV', True, 'TWADV', 1, 0],\
         ['TWDIFF', True, 'TWDIFF', 1, 0],\
         ['TWBUOY', False, 'TWBUOY', 1, 0],\
         ['TWPRES', False, 'TWPRES', 1, 0],\
         ['TWB+P', True, 'TWBUOY+TWPRES', 1, 0],\
         ['TWPREC', True, 'TWPREC', 1, 0],\
         ['TWRAD', True, 'TWRAD', 1, 0],\
         ['TWFORC', True, 'TWFORC', 1, 0],\
         ['TWBT', True, 'TWBT', 1, 0],\
         ['TW_RES', True, 'TWBT - (TWGRAD + TWADV + TWDIFF + TWBUOY + TWPRES + TWPREC + TWRAD + TWFORC)', 1, 0],\
        ]
         
THLW = [\
         # variables of THLW
         ['THLWGRAD', True, 'THLWGRAD', 1, 0],\
         ['THLWADV', True, 'THLWADV', 1, 0],\
         ['THLWDIFF', True, 'THLWDIFF', 1, 0],\
         ['THLWB+P', True, 'THLWBUOY+THLWPRES', 1, 0],\
         ['THLWBUOY', False, 'THLWBUOY', 1, 0],\
         ['THLWPRES', False, 'THLWPRES', 1, 0],\
         ['THLWPREC', True, 'THLWPREC', 1, 0],\
         ['THLWRAD', True, 'THLWRAD', 1, 0],\
         ['THLWFORC', True, 'THLWFORC', 1, 0],\
         ['THLWBT', True, 'THLWBT', 1, 0],\
         ['THLW_Res', True, 'THLWBT - (THLWGRAD + THLWADV + THLWDIFF + THLWBUOY + THLWPRES + THLWPREC + THLWRAD + THLWFORC)', 1, 0],\
        ]
         
QW = [\
         # variables of QW
         ['QWGRAD', True, 'QWGRAD', 1, 0],\
         ['QWADV', True, 'QWADV', 1, 0],\
         ['QWDIFF', True, 'QWDIFF', 1, 0],\
         ['QWB+P', True, 'QWBUOY+QWPRES', 1, 0],\
         ['QWBUOY', False, 'QWBUOY', 1, 0],\
         ['QWPRES', False, 'QWPRES', 1, 0],\
         ['QWPREC', True, 'QWPREC', 1, 0],\
         ['QWFORC', True, 'QWFORC', 1, 0],\
         ['QWBT', True, 'QWBT', 1, 0],\
         ['QW_Res', True, 'QWBT - (QWGRAD + QWADV + QWDIFF + QWBUOY + QWPRES + QWPREC + QWFORC)', 1, 0],\
        ]
         
QTOGW = [\
         # variables of QTOGW
         ['QTOGWGRAD', True, 'QTOGWGRAD', 1, 0],\
         ['QTOGWADV', True, 'QTOGWADV', 1, 0],\
         ['QTOGWDIFF', True, 'QTOGWDIFF', 1, 0],\
         ['QTOGWB+P', True, 'QTOGWBUOY+QTOGWPRES', 1, 0],\
         ['QTOGWBUOY', False, 'QTOGWBUOY', 1, 0],\
         ['QTOGWPRES', False, 'QTOGWPRES', 1, 0],\
         ['QTOGWPREC', True, 'QTOGWPREC', 1, 0],\
         ['QTOGWFORC', True, 'QTOGWFORC', 1, 0],\
         ['QTOGWBT', True, 'QTOGWBT', 1, 0],\
         ['QTOGW_Res', True, 'QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDIFF + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)', 1, 0],\
        ]
         
W2 = [\
         # variables of W2
         ['W2ADV', True, 'W2ADV', 1, 0],\
         ['W2PRES', True, 'W2PRES', 1, 0],\
         ['W2REDIS', True, 'W2REDIS', 1, 0],\
         ['W2BUOY', True, 'W2BUOY', 1, 0],\
         ['W2DIFF', True, 'W2DIFF', 1, 0],\
         ['W2SDMP', True, 'W2SDMP',1, 0],\
         ['W2BT', True, 'W2BT', 1, 0],\
         ['W2_Res', True, 'W2BT - (W2ADV + W2PRES + W2REDIS + W2BUOY + W2DIFF + W2SDMP)', 1, 0],\
        ]

W3 = [\
        # variables of W3
         ['W3ADV', True, 'W3ADV', 1, 0],\
         ['W3PRES', True, 'W3PRES', 1, 0],\
         ['W3REDIS', True, 'W3REDIS', 1, 0],\
         ['W3BUOY', True, 'W3BUOY', 1, 0],\
         ['W3DIFF', True, 'W3DIFF', 1, 0],\
         ['W3BT', True, 'W3BT', 1, 0],\
         ['W3_Res', True, 'W3BT - (W3ADV + W3PRES + W3REDIS + W3BUOY + W3DIFF)', 1, 0],\
       ]
                
T2 = [\
         # variables of T2
         ['T2ADVTR', True, 'T2ADVTR', 1, 0],\
         ['T2GRAD', True, 'T2GRAD', 1, 0],\
         ['T2DISSIP', True, 'T2DISSIP', 1, 0],\
         ['T2DIFTR', True, 'T2DIFTR', 1, 0],\
         ['T2PREC', True, 'T2PREC', 1, 0],\
         ['T2RAD', True, 'T2RAD', 1, 0],\
         ['T2FORC', True, 'T2FORC', 1, 0],\
         ['T2BT', True, 'T2BT', 1, 0],\
         ['T2_Res', True, 'T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)', 1, 0],\
        ]
         
THL2 = [\
         # variables of THL2
         ['THL2ADVTR', True, 'THL2ADVTR', 1, 0],\
         ['THL2GRAD', True, 'THL2GRAD', 1, 0],\
         ['THL2DISSIP', True, 'THL2DISSIP', 1, 0],\
         ['THL2DIFTR', True, 'THL2DIFTR', 1, 0],\
         ['THL2PREC', True, 'THL2PREC', 1, 0],\
         ['THL2RAD', True, 'THL2RAD', 1, 0],\
         ['THL2FORC', True, 'THL2FORC', 1, 0],\
         ['THL2BT', True, 'THL2BT', 1, 0],\
         ['THL2_Res', True, 'THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)', 1, 0],\
        ]
         
Q2 = [\
         # variables of Q2
         ['Q2ADVTR', True, 'Q2ADVTR', 1, 0],\
         ['Q2GRAD', True, 'Q2GRAD', 1, 0],\
         ['Q2DISSIP', True, 'Q2DISSIP', 1, 0],\
         ['Q2DIFTR', True, 'Q2DIFTR', 1, 0],\
         ['Q2PREC', True, 'Q2PREC', 1, 0],\
         ['Q2FORC', True, 'Q2FORC', 1, 0],\
         ['Q2BT', True, 'Q2BT', 1, 0],\
         ['Q2_Res', True, 'Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)', 1, 0],\
        ]
         
QTOG2 = [\
         # variables of QTOG2
         ['QTOG2ADVTR', True, 'QTOG2ADVTR', 1, 0],\
         ['QTOG2GRAD', True, 'QTOG2GRAD', 1, 0],\
         ['QTOG2DISSIP', True, 'QTOG2DISSIP', 1, 0],\
         ['QTOG2DIFTR', True, 'QTOG2DIFTR', 1, 0],\
         ['QTOG2PREC', True, 'QTOG2PREC', 1, 0],\
         ['QTOG2FORC', True, 'QTOG2FORC', 1, 0],\
         ['QTOG2BT', True, 'QTOG2BT', 1, 0],\
         ['QTOG2_Res', True, 'QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)', 1, 0],\
        ]
         
QTHL = [\
         # variables of QTHL
         ['QTHLADV', True, 'QTHLADV', 1, 0],\
         ['QTHLGRAD', True, 'QTHLGRAD', 1, 0],\
         ['QTHLDISSIP', True, 'QTHLDISSIP', 1, 0],\
         ['QTHLDIFTR', True, 'QTHLDIFTR', 1, 0],\
         ['QTHLPREC', True, 'QTHLPREC', 1, 0],\
         ['QTHLRAD', True, 'QTHLRAD', 1, 0],\
         ['QTHLFORC', True, 'QTHLFORC', 1, 0],\
         ['QTHLBT', True, 'QTHLBT', 1, 0],\
         ['QTHL_Res', True, 'QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)', 1, 0],\
        ]
        
TKE = [\
        # variables of TKE resolved
        ['ADVTR', True, 'ADVTR', 1, 0],\
        ['SHEAR', True, 'SHEAR', 1, 0],\
        ['BUOYA', True, 'BUOYA', 1, 0],\
        ['PRESSTR', True, 'PRESSTR', 1, 0],\
        ['DIFTR', False, 'DIFTR', 1, 0],\
        ['SDMP', True, 'SDMP', 1, 0],\
        ['DISSIP', False, 'DISSIP', 1, 0],\
        ['DIFTR+DISS', True, 'DISSIP + DIFTR', 1, 0],\
        ['TKEBT', True, 'BT', 1., 0],\
        ['TKE_RES', True, 'TKEBT - (SHEAR + BUOYA + ADVTR + PRESSTR + DIFTR + SDMP + DISSIP)', 1, 0],\
       ]
         
TKES = [\
        # variables of TKE SGS
        ['ADVTRS', True, 'ADVTRS', 1, 0],\
        ['SHEARS', True, 'SHEARS', 1, 0],\
        ['BUOYAS', True, 'BUOYAS', 1, 0],\
        ['DISSIPS', True, 'DISSIPS', -1., 0],\
        ['TKE_RES', True, '-(SHEARS + BUOYAS + ADVTRS + DISSIPS)', 1, 0],\
       ]
       
       # NOTE(Steffen Domke): TKESHEAR is duplicate
U2V2 = [\
        # variables of U2 + V2 = 2 * TKE - W2
        ['W2ADV', False, 'W2ADV', 1, 0],\
        ['W2PRES', False, 'W2PRES', 1, 0],\
        ['W2BUOY', False, 'W2BUOY', 1, 0],\
        ['W2DIFF', False, 'W2DIFF', 1, 0],\
        ['W2SDMP', False, 'W2SDMP', 1, 0],\
        ['W2REDIS', False, 'W2REDIS', 1, 0],\
        ['W2BT', False, 'W2BT', 1, 0],\
        ['TKEADVTR', False, 'ADVTR', 1, 0],\
        ['TKESHEAR', False, 'SHEAR', 1, 0],\
        ['TKEBUOYA', False, 'BUOYA', 1, 0],\
        ['TKEPRESSTR', False, 'PRESSTR', 1, 0],\
        ['TKEDIFTR', False, 'DIFTR', 1, 0],\
        ['TKEDISSIP', False, 'DISSIP', 1, 0],\
        ['TKESDMP', False, 'SDMP', 1, 0],\
        ['TKESHEAR', False, 'SHEAR', 1, 0],\
        ['TKEBT', False, 'BT', 1., 0],\
        ['U2V2ADV', True, '2. * TKEADVTR - W2ADV', 1., 0],\
        ['U2V2BUOY', True, '2. * TKEBUOYA - W2BUOY', 1., 0],\
        ['U2V2PRESS', True, '2. * TKEPRESSTR - W2PRES', 1., 0],\
        ['U2V2REDIS', True, '- W2REDIS', 1, 0],\
        ['U2V2DIFF', True, '2. * TKEDIFTR - W2DIFF', 1., 0],\
        ['U2V2DISSIP', True, '2. * TKEDISSIP', 1, 0],\
        ['U2V2SDMP', True, '2. * TKESDMP - W2SDMP', 1., 0],\
        ['U2V2SHEAR', True, '2. * TKESHEAR', 1, 0],\
        ['U2V2BT', True, '2. * TKEBT - W2BT', 1, 0],\
        ['U2V2_RES', True, 'U2V2BT - (U2V2ADV + U2V2BUOY + U2V2PRESS + U2V2DIFF + U2V2DISSIP + U2V2SDMP + U2V2REDIS + U2V2SHEAR)', 1., 0],\
       ]
  
WU = [\
    # variables of WV
        ['WUDIFF', True, 'WUDIFF', 1, 0],\
        ['WU_RES', True, 'WUBT - (WUDIFF + WUSHEAR + WUADV + WUPRES + WUANIZ + WUBUOY + WUSDMP)', 1, 0],\
        ['WUADV', True, 'WUADV', 1, 0],\
        ['WUPRES', True, 'WUPRES', 1, 0],\
        ['WUANIZ', True, 'WUANIZ', 1, 0],\
        ['WUBUOY', True, 'WUBUOY', 1, 0],\
        ['WUSHEAR', True, 'WUSHEAR', 1, 0],\
        ['WUBT', True, 'WUBT', 1, 0],\
        ['WUSDMP', True, 'WUSDMP', 1, 0],\
      ]

WV = [\
    # variables of WV
        ['WVDIFF', True, 'WVDIFF', 1, 0],\
        ['WV_RES', True, 'WVBT - (WVDIFF + WVSHEAR + WVADV + WVPRES + WVANIZ + WVBUOY + WVSDMP)', 1, 0],\
        ['WVADV', True, 'WVADV', 1, 0],\
        ['WVPRES', True, 'WVPRES', 1, 0],\
        ['WVANIZ', True, 'WVANIZ', 1, 0],\
        ['WVBUOY', True, 'WVBUOY', 1, 0],\
        ['WVSHEAR', True, 'WVSHEAR', 1, 0],\
        ['WVBT', True, 'WVBT', 1, 0],\
        ['WVSDMP', True, 'WVSDMP', 1, 0],\
     ]


U2 = [\
    # variables of U2
        ['U2ADV', True, 'U2ADV', 1, 0],\
        ['U2SHEAR', True, 'U2SHEAR', 1, 0],\
        ['U2REDIS', True, 'U2REDIS', 1, 0],\
        ['U2DIFF', True, 'U2DIFF', 1, 0],\
        ['U2BT', True, 'U2BT', 1, 0],\
        ['U2_RES', True, 'U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DIFF)', 1, 0],\
    ]

V2 = [\
    # variables of V2
        ['V2ADV', True, 'V2ADV', 1, 0],\
        ['V2SHEAR', True, 'V2SHEAR', 1, 0],\
        ['V2REDIS', True, 'V2REDIS', 1, 0],\
        ['V2DIFF', True, 'V2DIFF', 1, 0],\
        ['V2BT', True, 'V2BT', 1, 0],\
        ['V2_RES', True, 'V2BT - (V2ADV + V2SHEAR + V2REDIS + V2DIFF)', 1, 0],\
    ]

V2_COMP = [\
    # Use U2V2 to calculate V2
        ['U2ADV', False, 'U2ADV', 1, 0],\
        ['U2SHEAR', False, 'U2SHEAR', 1, 0],\
        ['U2REDIS', False, 'U2REDIS', 1, 0],\
        ['U2DIFF', False, 'U2DIFF', 1, 0],\
        ['U2BT', True, 'U2BT', 1, 0],\
        ['U2_RES', True, 'U2BT - (U2ADV + U2SHEAR + U2REDIS + U2DIFF)', 1, 0],\
        ['W2ADV', False, 'W2ADV', 1, 0],\
        ['W2PRES', False, 'W2PRES', 1, 0],\
        ['W2BUOY', False, 'W2BUOY', 1, 0],\
        ['W2DIFF', False, 'W2DIFF', 1, 0],\
        ['W2SDMP', False, 'W2SDMP', 1, 0],\
        ['W2REDIS', False, 'W2REDIS', 1, 0],\
        ['W2BT', False, 'W2BT', 1, 0],\
        ['TKEADVTR', False, 'ADVTR', 1, 0],\
        ['TKESHEAR', False, 'SHEAR', 1, 0],\
        ['TKEBUOYA', False, 'BUOYA', 1, 0],\
        ['TKEPRESSTR', False, 'PRESSTR', 1, 0],\
        ['TKEDIFTR', False, 'DIFTR', 1, 0],\
        ['TKEDISSIP', False, 'DISSIP', 1, 0],\
        ['TKESDMP', False, 'SDMP', 1, 0],\
        ['TKESHEAR', False, 'SHEAR', 1, 0],\
        ['TKEBT', False, 'BT', 1., 0],\
        ['U2V2ADV', False, '2. * TKEADVTR - W2ADV', 1., 0],\
        ['U2V2BUOY', False, '2. * TKEBUOYA - W2BUOY', 1., 0],\
        ['U2V2PRESS', False, '2. * TKEPRESSTR - W2PRES', 1., 0],\
        ['U2V2REDIS', False, '- W2REDIS', 1, 0],\
        ['U2V2DIFF', False, '2. * TKEDIFTR - W2DIFF', 1., 0],\
        ['U2V2DISSIP', False, '2. * TKEDISSIP', 1, 0],\
        ['U2V2SDMP', False, '2. * TKESDMP - W2SDMP', 1., 0],\
        ['U2V2SHEAR', False, '2. * TKESHEAR', 1, 0],\
        ['U2V2BT', True, '2. * TKEBT - W2BT', 1, 0],\
        ['U2V2_RES', True, 'U2V2BT - (U2V2ADV + U2V2BUOY + U2V2PRESS + U2V2DIFF + U2V2DISSIP + U2V2SDMP + U2V2REDIS + U2V2SHEAR)', 1., 0],\
        ['V2BT', True, 'U2V2BT - U2BT', 1., 0],\
        ['V2_RES', True, 'U2V2_RES - U2_RES', 1., 0],\
    ]

U2_V2_TAU = [\
    ['U2DIFF', True, 'U2DIFF', 1, 0],\
    ['V2DIFF', True, 'V2DIFF', 1, 0],\
    ['TKES', True, 'TKES', 1, 1],\
    ['TKE', True, 'TKE', 1, 1],\
    ['C14/tau_U2', True, '- (3/2) * U2DIFF / np.maximum( TKE + TKES, 1e-6 )', 1, 0],\
    ['C14/tau_V2', True, ' - (3/2) * V2DIFF / np.maximum( TKE + TKES, 1e-6 )', 1, 0],\
    ]


QTOG2_TAU = [\
    ['QTOG2DISSIP', True, 'QTOG2DISSIP', 1, 0],\
    ['QTOG2DIFTR', True, 'QTOG2DIFTR', 1, 0],\
    ['QTOG2', True, 'QTOG2', 1, 1],\
    ['C2/tau_QTOG2', True, '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2*1e-6, 1e-12 )', 1, 0],\
    ]

Q2_TAU = [\
    ['Q2DISSIP', True, 'Q2DISSIP', 1, 0],\
    ['Q2DIFTR', True, 'Q2DIFTR', 1, 0],\
    ['QT2', True, 'QT2', 1, 1],\
    ['C2/tau_Q2', True, '- ( Q2DISSIP + Q2DIFTR ) / np.maximum( QT2 * 1e-6, 1e-12 )', 1, 0],\
    ]

THL2_TAU = [\
    ['THL2DISSIP', True, 'THL2DISSIP', 1, 0],\
    ['THL2DIFTR', True, 'THL2DIFTR', 1, 0],\
    ['THEL2', True, 'THEL2', 1, 1],\
    ['C2/tau_THL2', True, ' - ( THL2DISSIP + THL2DIFTR ) / np.maximum( THEL2, 1e-6 )', 1, 0],\
    ]

C2_TAU_COMP = [\
    ['QTOG2DISSIP', False, 'QTOG2DISSIP', 1, 0],\
    ['QTOG2DIFTR', False, 'QTOG2DIFTR', 1, 0],\
    ['QTOG2', False, 'QTOG2', 1, 0],\
    ['THL2DISSIP', False, 'THL2DISSIP', 1, 0],\
    ['THL2DIFTR', False, 'THL2DIFTR', 1, 0],\
    ['THEL2', False, 'THEL2', 1, 0],\
    ['Q2DISSIP', False, 'Q2DISSIP', 1, 0],\
    ['Q2DIFTR', False, 'Q2DIFTR', 1, 0],\
    ['QT2', False, 'QT2', 1, 0],\
    ['U2DIFF', False, 'U2DIFF', 1, 0],\
    ['V2DIFF', False, 'V2DIFF', 1, 0],\
    ['TKES', False, 'TKES', 1, 0],\
    ['TKE', False, 'TKE', 1, 0],\
    ['C14/tau_U2', True, '- (3/2) * U2DIFF / np.maximum( TKE + TKES, 1e-6 )', 1, 0],\
    ['C14/tau_V2', True, ' - (3/2) * V2DIFF / np.maximum( TKE + TKES, 1e-6 )', 1, 0],\
    ['C2/tau_QTOG2', True, '- ( QTOG2DISSIP + QTOG2DIFTR ) / np.maximum( QTOG2 * 1e-6, 1e-12 )', 1, 0],\
    ['C2/tau_THL2', True, ' - ( THL2DISSIP + THL2DIFTR ) / np.maximum( THEL2, 1e-6 )', 1, 0],\
    ['C2/tau_Q2', True, '- ( Q2DISSIP + Q2DIFTR ) / np.maximum( QT2 * 1e-6, 1e-12 )', 1, 0],\
    ]

lines = [HL, QT, TW, THLW, QW, QTOGW, W2, W3, T2, THL2, Q2, QTOG2, QTHL, TKE, TKES, U2V2, WU, WV, U2, V2, V2_COMP, U2_V2_TAU, QTOG2_TAU, Q2_TAU, THL2_TAU, C2_TAU_COMP]