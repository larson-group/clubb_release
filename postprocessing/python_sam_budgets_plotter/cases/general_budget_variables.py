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

# settings of each plot:
# plot number, plot title, axis label
plotNames = [['01', 'HL 3D SAM Benchmark Budgets', 'HL [K/s]'],\
                ['02', 'QTO 3D SAM Benchmark Budgets', 'QTO [kg/kg/s]'],\
                ['03', 'TW SAM 3D Benchmark Budgets', 'TW [(m K)/s^2]'],\
                ['04', 'THLW SAM 3D Benchmark Budgets', 'THLW [(m K)/s^2]'],\
                ['05', 'QW SAM 3D Benchmark Budgets', 'QW [(kg/kg)m/s^2]'],\
                ['06', 'QTOGW SAM 3D Benchmark Budgets', 'QTOGW [(kg/kg)m/s^2]'],\
                ['07', 'W2 SAM 3D Benchmark Budgets', 'W2 [m^2/s^3]'],\
                ['08', 'T2 3D SAM Benchmark Budgets', 'T2 [K^2/s]'],\
                ['09', 'THL2 3D SAM Benchmark Budgets', 'THL2 [K^2/s]'],\
                ['10', 'Q2 3D SAM Benchmark Budgets', 'Q2 [(kg/kg)^2/s]'],\
                ['11', 'QTOG2 3D SAM Benchmark Budgets', 'QTOG2 [(kg/kg)^2/s]'],\
                ['12', 'QTHL 3D SAM Benchmark Budgets', 'QTHL [(kg/kg)K/s]']\
               ]

# lines of each plot:
# plot number, variable name, plot?, computation, conversion

lines = [\
         # variables of HL
         ['01', 'HLADV', True, 'HLADV', kg_per_second_to_kg_per_day],\
         ['01', 'HLDIFF', True, 'HLDIFF', kg_per_second_to_kg_per_day],\
         ['01', 'HLRAD', True, 'HLRAD', kg_per_second_to_kg_per_day],\
         ['01', 'HLLAT', True, 'HLLAT', kg_per_second_to_kg_per_day],\
         ['01', 'TTEND', True, 'TTEND', kg_per_second_to_kg_per_day],\
         ['01', 'T_TNDCY', True, 'T_TNDCY', 1.0],\
         ['01', 'HLSTOR', True, 'HLSTOR', kg_per_second_to_kg_per_day],\
         ['01', 'HL_RES', True, '((HLSTOR)+((-1)*(HLADV+HLDIFF+HLRAD+HLLAT+TTEND)))', g_per_second_to_kg_per_day],\
         
         # variables of QTO
         ['02', 'QTADV', True, 'QTADV', g_per_second_to_kg_per_day],\
         ['02', 'QTDIFF', True, 'QTDIFF', g_per_second_to_kg_per_day],\
         #['02', 'QTLSADV', True, 'QTLSADV', g_per_second_to_kg_per_day],\
         ['02', 'QTSRC', True, 'QTSRC', g_per_second_to_kg_per_day],\
         ['02', 'QTSINK', True, 'QTSINK', g_per_second_to_kg_per_day],\
         ['02', 'QTEND', True, 'QTEND', g_per_second_to_kg_per_day],\
         ['02', 'QTSTOR', True, 'QTSTOR', g_per_second_to_kg_per_day],\
         ['02', 'QV_TNDCY', True, 'QV_TNDCY', g_per_second_to_kg_per_day],\
         ['02', 'QC_TNDCY', True, 'QC_TNDCY', g_per_second_to_kg_per_day],\
         ['02', 'QTO_RES', True, 'QTSTOR+(-1)*(QTADV+QTDIFF+QTSRC+QTSINK+QTEND)', 1],\
         
         # variables of TW
         ['03', 'TWGRAD', True, 'TWGRAD', 1],\
         ['03', 'TWADV', True, 'TWADV', 1],\
         ['03', 'TWDIFF', True, 'TWDIFF', 1],\
         ['03', 'TWBUOY', False, 'TWBUOY', 1],\
         ['03', 'TWPRES', False, 'TWPRES', 1],\
         ['03', 'TWB+P', True, 'TWBUOY+TWPRES', 1],\
         ['03', 'TWPREC', True, 'TWPREC', 1],\
         ['03', 'TWRAD', True, 'TWRAD', 1],\
         ['03', 'TWFORC', True, 'TWFORC', 1],\
         ['03', 'TWBT', True, 'TWBT', 1],\
         ['03', 'TW_RES', True, 'TWBT - (TWGRAD + TWADV + TWDIFF + TWBUOY + TWPRES + TWPREC + TWRAD + TWFORC)', 1],\
         
         # variables of THLW
         ['04', 'THLWGRAD', True, 'THLWGRAD', 1],\
         ['04', 'THLWADV', True, 'THLWADV', 1],\
         ['04', 'THLWDIFF', True, 'THLWDIFF', 1],\
         ['04', 'THLWB+P', True, 'THLWBUOY+THLWPRES', 1],\
         ['04', 'THLWBUOY', False, 'THLWBUOY', 1],\
         ['04', 'THLWPRES', False, 'THLWPRES', 1],\
         ['04', 'THLWPREC', True, 'THLWPREC', 1],\
         ['04', 'THLWRAD', True, 'THLWRAD', 1],\
         ['04', 'THLWFORC', True, 'THLWFORC', 1],\
         ['04', 'THLWBT', True, 'THLWBT', 1],\
         ['04', 'THLW_Res', True, 'THLWBT - (THLWGRAD + THLWADV + THLWDIFF + THLWBUOY + THLWPRES + THLWPREC + THLWRAD + THLWFORC)', 1],\
         
         # variables of QW
         ['05', 'QWGRAD', True, 'QWGRAD', 1],\
         ['05', 'QWADV', True, 'QWADV', 1],\
         ['05', 'QWDIFF', True, 'QWDIFF', 1],\
         ['05', 'QWB+P', True, 'QWBUOY+QWPRES', 1],\
         ['05', 'QWBUOY', False, 'QWBUOY', 1],\
         ['05', 'QWPRES', False, 'QWPRES', 1],\
         ['05', 'QWPREC', True, 'QWPREC', 1],\
         ['05', 'QWFORC', True, 'QWFORC', 1],\
         ['05', 'QWBT', True, 'QWBT', 1],\
         ['05', 'QW_Res', True, 'QWBT - (QWGRAD + QWADV + QWDIFF + QWBUOY + QWPRES + QWPREC + QWFORC)', 1],\
         
         # variables of QTOGW
         ['06', 'QTOGWGRAD', True, 'QTOGWGRAD', 1],\
         ['06', 'QTOGWADV', True, 'QTOGWADV', 1],\
         ['06', 'QTOGWDIFF', True, 'QTOGWDIFF', 1],\
         ['06', 'QTOGWB+P', True, 'QTOGWBUOY+QTOGWPRES', 1],\
         ['06', 'QTOGWBUOY', False, 'QTOGWBUOY', 1],\
         ['06', 'QTOGWPRES', False, 'QTOGWPRES', 1],\
         ['06', 'QTOGWPREC', True, 'QTOGWPREC', 1],\
         ['06', 'QTOGWFORC', True, 'QTOGWFORC', 1],\
         ['06', 'QTOGWBT', True, 'QTOGWBT', 1],\
         ['06', 'QTOGW_Res', True, 'QTOGWBT - (QTOGWGRAD + QTOGWADV + QTOGWDIFF + QTOGWBUOY + QTOGWPRES + QTOGWPREC + QTOGWFORC)', 1],\
         
         # variables of W2
         ['07', 'W2ADV', True, 'W2ADV', 1],\
         ['07', 'W2PRES', True, 'W2PRES', 1],\
         ['07', 'W2REDIS', True, 'W2REDIS', 1],\
         ['07', 'W2BUOY', True, 'W2BUOY', 1],\
         ['07', 'W2DIFF', True, 'W2DIFF', 1],\
         ['07', 'W2_Res', True, '-1*(W2ADV + W2PRES + W2REDIS + W2BUOY + W2DIFF)', 1],\
         
         # variables of T2
         ['08', 'T2ADVTR', True, 'T2ADVTR', 1],\
         ['08', 'T2GRAD', True, 'T2GRAD', 1],\
         ['08', 'T2DISSIP', True, 'T2DISSIP', 1],\
         ['08', 'T2DIFTR', True, 'T2DIFTR', 1],\
         ['08', 'T2PREC', True, 'T2PREC', 1],\
         ['08', 'T2RAD', True, 'T2RAD', 1],\
         ['08', 'T2FORC', True, 'T2FORC', 1],\
         ['08', 'T2BT', True, 'T2BT', 1],\
         ['08', 'T2_Res', True, 'T2BT - (T2ADVTR + T2GRAD + T2DISSIP + T2DIFTR + T2PREC + T2RAD + T2FORC)', 1],\
         
         # variables of THL2
         ['09', 'THL2ADVTR', True, 'THL2ADVTR', 1],\
         ['09', 'THL2GRAD', True, 'THL2GRAD', 1],\
         ['09', 'THL2DISSIP', True, 'THL2DISSIP', 1],\
         ['09', 'THL2DIFTR', True, 'THL2DIFTR', 1],\
         ['09', 'THL2PREC', True, 'THL2PREC', 1],\
         ['09', 'THL2RAD', True, 'THL2RAD', 1],\
         ['09', 'THL2FORC', True, 'THL2FORC', 1],\
         ['09', 'THL2BT', True, 'THL2BT', 1],\
         ['09', 'THL2_Res', True, 'THL2BT - (THL2ADVTR + THL2GRAD + THL2DISSIP + THL2DIFTR + THL2PREC + THL2RAD + THL2FORC)', 1],\
         
         # variables of Q2
         ['10', 'Q2ADVTR', True, 'Q2ADVTR', 1],\
         ['10', 'Q2GRAD', True, 'Q2GRAD', 1],\
         ['10', 'Q2DISSIP', True, 'Q2DISSIP', 1],\
         ['10', 'Q2DIFTR', True, 'Q2DIFTR', 1],\
         ['10', 'Q2PREC', True, 'Q2PREC', 1],\
         ['10', 'Q2FORC', True, 'Q2FORC', 1],\
         ['10', 'Q2BT', True, 'Q2BT', 1],\
         ['10', 'Q2_Res', True, 'Q2BT - (Q2ADVTR + Q2GRAD + Q2DISSIP + Q2DIFTR + Q2PREC + Q2FORC)', 1],\
         
         # variables of QTOG2
         ['11', 'QTOG2ADVTR', True, 'QTOG2ADVTR', 1],\
         ['11', 'QTOG2GRAD', True, 'QTOG2GRAD', 1],\
         ['11', 'QTOG2DISSIP', True, 'QTOG2DISSIP', 1],\
         ['11', 'QTOG2DIFTR', True, 'QTOG2DIFTR', 1],\
         ['11', 'QTOG2PREC', True, 'QTOG2PREC', 1],\
         ['11', 'QTOG2FORC', True, 'QTOG2FORC', 1],\
         ['11', 'QTOG2BT', True, 'QTOG2BT', 1],\
         ['11', 'QTOG2_Res', True, 'QTOG2BT - (QTOG2ADVTR + QTOG2GRAD + QTOG2DISSIP + QTOG2DIFTR + QTOG2PREC + QTOG2FORC)', 1],\
         
         # variables of QTHL
         ['12', 'QTHLADV', True, 'QTHLADV', 1],\
         ['12', 'QTHLGRAD', True, 'QTHLGRAD', 1],\
         ['12', 'QTHLDISSIP', True, 'QTHLDISSIP', 1],\
         ['12', 'QTHLDIFTR', True, 'QTHLDIFTR', 1],\
         ['12', 'QTHLPREC', True, 'QTHLPREC', 1],\
         ['12', 'QTHLRAD', True, 'QTHLRAD', 1],\
         ['12', 'QTHLFORC', True, 'QTHLFORC', 1],\
         ['12', 'QTHLBT', True, 'QTHLBT', 1],\
         ['12', 'QTHL_Res', True, 'QTHLBT - (QTHLADV + QTHLGRAD + QTHLDISSIP + QTHLDIFTR + QTHLPREC + QTHLRAD + QTHLFORC)', 1],\
        ]