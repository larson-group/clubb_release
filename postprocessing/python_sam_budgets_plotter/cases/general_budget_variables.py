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
sortPlots = ['HL', 'QT', 'TW', 'THLW', 'QW', 'QTOGW', 'W2', 'T2', 'THL2', 'Q2', 'QTOG2', 'QTHL']
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
                ['T2 3D SAM Benchmark Budgets', 'T2 [K^2/s]'],\
                ['THL2 3D SAM Benchmark Budgets', 'THL2 [K^2/s]'],\
                ['Q2 3D SAM Benchmark Budgets', 'Q2 [(kg/kg)^2/s]'],\
                ['QTOG2 3D SAM Benchmark Budgets', 'QTOG2 [(kg/kg)^2/s]'],\
                ['QTHL 3D SAM Benchmark Budgets', 'QTHL [(kg/kg)K/s]']\
            ]

# lines of each plot:
# plot number, variable name, plot?, computation, conversion

HL = [\
         # variables of HL
         ['HLADV', True, 'HLADV', kg_per_second_to_kg_per_day],\
         ['HLDIFF', True, 'HLDIFF', kg_per_second_to_kg_per_day],\
         ['HLRAD', True, 'HLRAD', kg_per_second_to_kg_per_day],\
         ['HLLAT', True, 'HLLAT', kg_per_second_to_kg_per_day],\
         ['TTEND', True, 'TTEND', kg_per_second_to_kg_per_day],\
         ['T_TNDCY', True, 'T_TNDCY', 1.0],\
         ['HLSTOR', True, 'HLSTOR', kg_per_second_to_kg_per_day],\
         ['HL_RES', True, '((HLSTOR)+((-1)*(HLADV+HLDIFF+HLRAD+HLLAT+TTEND)))', g_per_second_to_kg_per_day],\
        ]
         
QT = [\
         # variables of QTO
         ['QTADV', True, 'QTADV', g_per_second_to_kg_per_day],\
         ['QTDIFF', True, 'QTDIFF', g_per_second_to_kg_per_day],\
         #['QTLSADV', True, 'QTLSADV', g_per_second_to_kg_per_day],\
         ['QTSRC', True, 'QTSRC', g_per_second_to_kg_per_day],\
         ['QTSINK', True, 'QTSINK', g_per_second_to_kg_per_day],\
         ['QTEND', True, 'QTEND', g_per_second_to_kg_per_day],\
         ['QTSTOR', True, 'QTSTOR', g_per_second_to_kg_per_day],\
         ['QV_TNDCY', True, 'QV_TNDCY', g_per_second_to_kg_per_day],\
         ['QC_TNDCY', True, 'QC_TNDCY', g_per_second_to_kg_per_day],\
         ['QT_RES', True, 'QTSTOR+(-1)*(QTADV+QTDIFF+QTSRC+QTSINK+QTEND)', 1],\
        ]
         
TW = [\
         # variables of TW
         ['TWGRAD', True, 'TWGRAD', 1],\
         ['TWADV', True, 'TWADV', 1],\
         ['TWDIFF', True, 'TWDIFF', 1],\
         ['TWBUOY', False, 'TWBUOY', 1],\
         ['TWPRES', False, 'TWPRES', 1],\
         ['TWB+P', True, 'TWBUOY+TWPRES', 1],\
         ['TWPREC', True, 'TWPREC', 1],\
         ['TWRAD', True, 'TWRAD', 1],\
         ['TWFORC', True, 'TWFORC', 1],\
         ['TWBT', True, 'TWBT', 1],\
         ['TW_RES', True, 'TWBT - (TWGRAD + TWADV + TWDIFF + TWBUOY + TWPRES + TWPREC + TWRAD + TWFORC)', 1],\
        ]
         
THLW = [\
         # variables of THLW
         ['THLWGRAD', True, 'THLWGRAD', 1],\
         ['THLWADV', True, 'THLWADV', 1],\
         ['THLWDIFF', True, 'THLWDIFF', 1],\
         ['THLWB+P', True, 'THLWBUOY+THLWPRES', 1],\
         ['THLWBUOY', False, 'THLWBUOY', 1],\
         ['THLWPRES', False, 'THLWPRES', 1],\
         ['THLWPREC', True, 'THLWPREC', 1],\
         ['THLWRAD', True, 'THLWRAD', 1],\
         ['THLWFORC', True, 'THLWFORC', 1],\
         ['THLWBT', True, 'THLWBT', 1],\
         ['THLW_Res', True, 'THLWBT - (THLWGRAD + THLWADV + THLWDIFF + THLWBUOY + THLWPRES + THLWPREC + THLWRAD + THLWFORC)', 1],\
        ]
         
QW = [\
         # variables of QW
         ['QWGRAD', True, 'QWGRAD', 1],\
         ['QWADV', True, 'QWADV', 1],\
         ['QWDIFF', True, 'QWDIFF', 1],\
         ['QWB+P', True, 'QWBUOY+QWPRES', 1],\
         ['QWBUOY', False, 'QWBUOY', 1],\
         ['QWPRES', False, 'QWPRES', 1],\
         ['QWPREC', True, 'QWPREC', 1],\
         ['QWFORC', True, 'QWFORC', 1],\
         ['QWBT', True, 'QWBT', 1],\
         ['QW_Res', True, 'QWBT - (QWGRAD + QWADV + QWDIFF + QWBUOY + QWPRES + QWPREC + QWFORC)', 1],\
        ]
         
QTOGW = [\
         # variables of QTOGW
         ['QTOGWGRAD', True, 'QTOGWGRAD', 1],\
         ['QTOGWADV', True, 'QTOGWADV', 1],\
         ['QTOGWDIFF', True, 'QTOGWDIFF', 1],\
         ['QTOGWB+P', True, 'QTOGWBUOY+QTOGWPRES', 1],\
         ['QTOGWBUOY', False, 'QTOGWBUOY', 1],\
         ['QTOGWPRES', False, 'QTOGWPRES', 1],\
         ['QTOGWPREC', True, 'QTOGWPREC', 1],\
         ['QTOGWFORC', True, 'QTOGWFORC', 1],\
         ['QTOGWBT', True, 'QTOGWBT', 1],\
         ['QTOGW_Res', True, 'QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDIFF + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)', 1],\
        ]
         
W2 = [\
         # variables of W2
         ['W2ADV', True, 'W2ADV', 1],\
         ['W2PRES', True, 'W2PRES', 1],\
         ['W2REDIS', True, 'W2REDIS', 1],\
         ['W2BUOY', True, 'W2BUOY', 1],\
         ['W2DIFF', True, 'W2DIFF', 1],\
         ['W2_Res', True, '-1*(W2ADV + W2PRES + W2REDIS + W2BUOY + W2DIFF)', 1],\
        ]
         
T2 = [\
         # variables of T2
         ['T2ADVTR', True, 'T2ADVTR', 1],\
         ['T2GRAD', True, 'T2GRAD', 1],\
         ['T2DISSIP', True, 'T2DISSIP', 1],\
         ['T2DIFTR', True, 'T2DIFTR', 1],\
         ['T2PREC', True, 'T2PREC', 1],\
         ['T2RAD', True, 'T2RAD', 1],\
         ['T2FORC', True, 'T2FORC', 1],\
         ['T2BT', True, 'T2BT', 1],\
         ['T2_Res', True, 'T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)', 1],\
        ]
         
THL2 = [\
         # variables of THL2
         ['THL2ADVTR', True, 'THL2ADVTR', 1],\
         ['THL2GRAD', True, 'THL2GRAD', 1],\
         ['THL2DISSIP', True, 'THL2DISSIP', 1],\
         ['THL2DIFTR', True, 'THL2DIFTR', 1],\
         ['THL2PREC', True, 'THL2PREC', 1],\
         ['THL2RAD', True, 'THL2RAD', 1],\
         ['THL2FORC', True, 'THL2FORC', 1],\
         ['THL2BT', True, 'THL2BT', 1],\
         ['THL2_Res', True, 'THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)', 1],\
        ]
         
Q2 = [\
         # variables of Q2
         ['Q2ADVTR', True, 'Q2ADVTR', 1],\
         ['Q2GRAD', True, 'Q2GRAD', 1],\
         ['Q2DISSIP', True, 'Q2DISSIP', 1],\
         ['Q2DIFTR', True, 'Q2DIFTR', 1],\
         ['Q2PREC', True, 'Q2PREC', 1],\
         ['Q2FORC', True, 'Q2FORC', 1],\
         ['Q2BT', True, 'Q2BT', 1],\
         ['Q2_Res', True, 'Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)', 1],\
        ]
         
QTOG2 = [\
         # variables of QTOG2
         ['QTOG2ADVTR', True, 'QTOG2ADVTR', 1],\
         ['QTOG2GRAD', True, 'QTOG2GRAD', 1],\
         ['QTOG2DISSIP', True, 'QTOG2DISSIP', 1],\
         ['QTOG2DIFTR', True, 'QTOG2DIFTR', 1],\
         ['QTOG2PREC', True, 'QTOG2PREC', 1],\
         ['QTOG2FORC', True, 'QTOG2FORC', 1],\
         ['QTOG2BT', True, 'QTOG2BT', 1],\
         ['QTOG2_Res', True, 'QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)', 1],\
        ]
         
QTHL = [\
         # variables of QTHL
         ['QTHLADV', True, 'QTHLADV', 1],\
         ['QTHLGRAD', True, 'QTHLGRAD', 1],\
         ['QTHLDISSIP', True, 'QTHLDISSIP', 1],\
         ['QTHLDIFTR', True, 'QTHLDIFTR', 1],\
         ['QTHLPREC', True, 'QTHLPREC', 1],\
         ['QTHLRAD', True, 'QTHLRAD', 1],\
         ['QTHLFORC', True, 'QTHLFORC', 1],\
         ['QTHLBT', True, 'QTHLBT', 1],\
         ['QTHL_Res', True, 'QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)', 1],\
        ]

lines = [HL, QT, TW, THLW, QW, QTOGW, W2, T2, THL2, Q2, QTOG2, QTHL]