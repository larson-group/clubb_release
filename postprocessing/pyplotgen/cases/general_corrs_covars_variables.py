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
sortPlots = ['corr_chi_w', 'covarnce_chi_w', 'corr_chi_rr', 'covarnce_chi_rr', 'corr_chi_Nr', 'covarnce_chi_Nr', 'corr_chi_ri', 'covarnce_chi_ri', 'corr_chi_Ni', 'covarnce_chi_Ni', 'corr_chi_rs', 'covarnce_chi_rs', 'corr_chi_Ns', 'covarnce_chi_Ns', 'corr_chi_rg', 'covarnce_chi_rg', 'corr_chi_Ng', 'covarnce_chi_Ng',\
            'corr_w_rr', 'covarnce_w_rr', 'corr_w_Nr', 'covarnce_w_Nr', 'corr_w_ri', 'covarnce_w_ri', 'corr_w_Ni', 'covarnce_w_Ni', 'corr_w_rs', 'covarnce_w_rs', 'corr_w_Ns', 'covarnce_w_Ns', 'corr_w_rg', 'covarnce_w_rg', 'corr_w_Ng', 'covarnce_w_Ng',\
            'corr_rr_Nr', 'covarnce_rr_Nr', 'corr_rr_ri', 'covarnce_rr_ri', 'corr_rr_Ni', 'covarnce_rr_Ni', 'corr_rr_rs', 'covarnce_rr_rs', 'corr_rr_Ns', 'covarnce_rr_Ns', 'corr_rr_rg', 'covarnce_rr_rg', 'corr_rr_Ng', 'covarnce_rr_Ng',\
            'corr_Nr_ri', 'covarnce_Nr_ri', 'corr_Nr_Ni', 'covarnce_Nr_Ni', 'corr_Nr_rs', 'covarnce_Nr_rs', 'corr_Nr_Ns', 'covarnce_Nr_Ns', 'corr_Nr_rg', 'covarnce_Nr_rg', 'corr_Nr_Ng', 'covarnce_Nr_Ng',\
            'corr_ri_Ni', 'covarnce_ri_Ni', 'corr_ri_rs', 'covarnce_ri_rs', 'corr_ri_Ns', 'covarnce_ri_Ns', 'corr_ri_rg', 'covarnce_ri_rg', 'corr_ri_Ng', 'covarnce_ri_Ng',\
            'corr_Ni_rs', 'covarnce_Ni_rs', 'corr_Ni_Ns', 'covarnce_Ni_Ns', 'corr_Ni_rg', 'covarnce_Ni_rg', 'corr_Ni_Ng', 'covarnce_Ni_Ng',\
            'corr_rs_Ns', 'covarnce_rs_Ns', 'corr_rs_rg', 'covarnce_rs_rg', 'corr_rs_Ng', 'covarnce_rs_Ng',\
            'corr_Ns_rg', 'covarnce_Ns_rg', 'corr_Ns_Ng', 'covarnce_Ns_Ng',\
            'corr_rg_Ng', 'covarnce_rg_Ng']
# settings of each plot:
# plot number, plot title, axis label
plotNames = [\
                ['Correlation of CHI and W', '[-]'],\
                ['Covariance of CHI and W', '[ (kg kg^-1) (m s^-1)]'],\
                ['Correlation of CHI and RR', '[-]'],\
                ['Covariance of CHI and RR', '[(kg kg^-1) (kg/kg^-1)]'],\
                ['Correlation of CHI and NR', '[-]'],\
                ['Covariance of CHI and NR', '[(kg kg^-1) (#/kg^-1) ]'],\
                ['Correlation of CHI and RI', '[-]'],\
                ['Covariance of CHI and RI', '[(kg kg^-1) (kg/kg^-1)]'],\
                ['Correlation of CHI and NI', '[-]'],\
                ['Covariance of CHI and NI', '[(kg kg^-1) (#/kg^-1) ]'],\
                ['Correlation of CHI and RS', '[-]'],\
                ['Covariance of CHI and RS', '[(kg kg^-1) (kg/kg^-1)]'],\
                ['Correlation of CHI and NS', '[-]'],\
                ['Covariance of CHI and NS', '[(kg kg^-1) (#/kg^-1) ]'],\
                ['Correlation of CHI and RG', '[-]'],\
                ['Covariance of CHI and RG', '[(kg kg^-1) (kg/kg^-1)'],\
                ['Correlation of CHI and NG', '[-]'],\
                ['Covariance of CHI and NG', '[(kg kg^-1) (#/kg^-1) ]'],\
                ['Correlation of W and RR', '[-]'],\
                ['Covariance of W and RR', '[(m s^-1) (kg kg^-1)]'],\
                ['Correlation of W and NR', '[-]'],\
                ['Covariance of W and NR', '[(m s^-1) (# kg^-1)]'],\
                ['Correlation of W and RI', '[-]'],\
                ['Covariance of W and RI', '[(m s^-1) (kg kg^-1)]'],\
                ['Correlation of W and NI', '[-]'],\
                ['Covariance of W and NI', '[(m s^-1) (# kg^-1)]'],\
                ['Correlation of W and RS', '[-]'],\
                ['Covariance of W and RS', '[(m s^-1) (kg kg^-1)]'],\
                ['Correlation of W and NS', '[-]'],\
                ['Covariance of W and NS', '[(m s^-1) (# kg^-1)]'],\
                ['Correlation of W and RG', '[-]'],\
                ['Covariance of W and RG', '[(m s^-1) (kg kg^-1)]'],\
                ['Correlation of W and NG', '[-]'],\
                ['Covariance of W and NG', '[(m s^-1) (# kg^-1)]'],\
                ['Correlation of RR and NR', '[-]'],\
                ['Covariance of RR and NR', '[ (kg kg^-1) (# kg^-1) ]'],\
                ['Correlation of RR and RI', '[-]'],\
                ['Covariance of RR and RI', '[(kg kg^-1) (kg kg^-1) ]'],\
                ['Correlation of RR and NI', '[-]'],\
                ['Covariance of RR and NI', '[ (kg kg^-1) (# kg^-1) ]'],\
                ['Correlation of RR and RS', '[-]'],\
                ['Covariance of RR and RS', '[(kg kg^-1) (kg kg^-1) ]'],\
                ['Correlation of RR and NS', '[-]'],\
                ['Covariance of RR and NS', '[ (kg kg^-1) (# kg^-1) ]'],\
                ['Correlation of RR and RG', '[-]'],\
                ['Covariance of RR and RG', '[(kg kg^-1) (kg kg^-1) ]'],\
                ['Correlation of RR and NG', '[-]'],\
                ['Covariance of RR and NG', '[ (kg kg^-1) (# kg^-1) ]'],\
                ['Correlation of NR and RI', '[-]'],\
                ['Covariance of NR and RI', '[ (# kg^-1) (kg kg^-1) ]'],\
                ['Correlation of NR and NI', '[-]'],\
                ['Covariance of NR and NI', '[(# kg^-1) (# kg^-1) ]'],\
                ['Correlation of NR and RS', '[-]'],\
                ['Covariance of NR and RS', '[ (# kg^-1) (kg kg^-1) ]'],\
                ['Correlation of NR and NS', '[-]'],\
                ['Covariance of NR and NS', '[(# kg^-1) (# kg^-1) ]'],\
                ['Correlation of NR and RG', '[-]'],\
                ['Covariance of NR and RG', '[ (# kg^-1) (kg kg^-1) ]'],\
                ['Correlation of NR and NG', '[-]'],\
                ['Covariance of NR and NG', '[(# kg^-1) (# kg^-1) ]'],\
                ['Correlation of RI and NI', '[-]'],\
                ['Covariance of RI and NI', '[(kg kg^-1) (# kg^-1)]'],\
                ['Correlation of RI and RS', '[-]'],\
                ['Covariance of RI and RS', '[(kg kg^-1) (kg kg^-1) ]'],\
                ['Correlation of RI and NS', '[-]'],\
                ['Covariance of RI and NS', '[(kg kg^-1) (# kg^-1)]'],\
                ['Correlation of RI and RG', '[-]'],\
                ['Covariance of RI and RG', '[(kg kg^-1) (kg kg^-1) ]'],\
                ['Correlation of RI and NG', '[-]'],\
                ['Covariance of RI and NG', '[(kg kg^-1) (# kg^-1)]'],\
                ['Correlation of NI and RS', '[-]'],\
                ['Covariance of NI and RS', '[(# kg^-1) (kg kg^-1)]'],\
                ['Correlation of NI and NS', '[-]'],\
                ['Covariance of NI and NS', '[(# kg^-1) (# kg^-1)]'],\
                ['Correlation of NI and RG', '[-]'],\
                ['Covariance of NI and RG', '[(# kg^-1) (kg kg^-1)]'],\
                ['Correlation of NI and NG', '[-]'],\
                ['Covariance of NI and NG', '[(# kg^-1) (# kg^-1)]'],\
                ['Correlation of RS and NS', '[-]'],\
                ['Covariance of RS and NS', '[(kg kg^-1) (# kg^-1)]'],\
                ['Correlation of RS and RG', '[-]'],\
                ['Covariance of RS and RG', '[(kg kg^-1) (kg kg^-1)]'],\
                ['Correlation of RS and NG', '[-]'],\
                ['Covariance of RS and NG', '[(kg kg^-1) (# kg^-1)]'],\
                ['Correlation of NS and RG', '[-]'],\
                ['Covariance of NS and RG', '[(# kg^-1) (kg kg^-1)]'],\
                ['Correlation of NS and NG', '[-]'],\
                ['Covariance of NS and NG', '[(# kg^-1) (# kg^-1)]'],\
                ['Correlation of RR and NG', '[-]'],\
                ['Covariance of RR and NG', '[(kg kg^-1) (# kg^-1)]'],\
            ]

# lines of each plot:
# variable name within python, shall this variable be plotted?, variable name in SAM output, conversion

corr_chi_w = [ ['', True, 'corr_chi_w', 1.] ]
        
covarnce_chi_w = [ ['', True, 'covarnce_chi_w', 1.] ]
        
corr_chi_rr = [ ['', True, 'corr_chi_rr', 1.] ]
        
covarnce_chi_rr = [ ['', True, 'covarnce_chi_rr', 1.] ]
        
corr_chi_Nr = [ ['', True, 'corr_chi_Nr', 1.] ]
        
covarnce_chi_Nr = [ ['', True, 'covarnce_chi_Nr', 1.] ]
         
corr_chi_ri = [ ['', True, 'corr_chi_ri', 1.] ]
        
covarnce_chi_ri = [ ['', True, 'covarnce_chi_ri', 1.] ]
        
corr_chi_Ni = [ ['', True, 'corr_chi_Ni', 1.] ]
        
covarnce_chi_Ni = [ ['', True, 'covarnce_chi_Ni', 1.] ]

corr_chi_rs = [ ['', True, 'corr_chi_rs', 1.] ]
        
covarnce_chi_rs = [ ['', True, 'covarnce_chi_rs', 1.] ]

corr_chi_Ns = [ ['', True, 'corr_chi_Ns', 1.] ]
        
covarnce_chi_Ns = [ ['', True, 'covarnce_chi_Ns', 1.] ]

corr_chi_rg = [ ['', True, 'corr_chi_rg', 1.] ]
        
covarnce_chi_rg = [ ['', True, 'covarnce_chi_rg', 1.] ]

corr_chi_Ng = [ ['', True, 'corr_chi_Ng', 1.] ]
        
covarnce_chi_Ng = [ ['', True, 'covarnce_chi_Ng', 1.] ]
    
corr_w_rr = [ ['', True, 'corr_w_rr', 1.] ]
        
covarnce_w_rr = [ ['', True, 'covarnce_w_rr', 1.] ]
    
corr_w_Nr = [ ['', True, 'corr_w_Nr', 1.] ]
        
covarnce_w_Nr = [ ['', True, 'covarnce_w_Nr', 1.] ]
         
corr_w_ri = [ ['', True, 'corr_w_ri', 1.] ]
        
covarnce_w_ri = [ ['', True, 'covarnce_w_ri', 1.] ]
        
corr_w_Ni = [ ['', True, 'corr_w_Ni', 1.] ]
        
covarnce_w_Ni = [ ['', True, 'covarnce_w_Ni', 1.] ]

corr_w_rs = [ ['', True, 'corr_w_rs', 1.] ]
        
covarnce_w_rs = [ ['', True, 'covarnce_w_rs', 1.] ]

corr_w_Ns = [ ['', True, 'corr_w_Ns', 1.] ]
        
covarnce_w_Ns = [ ['', True, 'covarnce_w_Ns', 1.] ]

corr_w_rg = [ ['', True, 'corr_w_rg', 1.] ]
        
covarnce_w_rg = [ ['', True, 'covarnce_w_rg', 1.] ]

corr_w_Ng = [ ['', True, 'corr_w_Ng', 1.] ]
        
covarnce_w_Ng = [ ['', True, 'covarnce_w_Ng', 1.] ]

corr_rr_Nr = [ ['', True, 'corr_rr_Nr', 1.] ]
        
covarnce_rr_Nr = [ ['', True, 'covarnce_rr_Nr', 1.] ]
         
corr_rr_ri = [ ['', True, 'corr_rr_ri', 1.] ]
        
covarnce_rr_ri = [ ['', True, 'covarnce_rr_ri', 1.] ]
        
corr_rr_Ni = [ ['', True, 'corr_rr_Ni', 1.] ]
        
covarnce_rr_Ni = [ ['', True, 'covarnce_rr_Ni', 1.] ]

corr_rr_rs = [ ['', True, 'corr_rr_rs', 1.] ]
        
covarnce_rr_rs = [ ['', True, 'covarnce_rr_rs', 1.] ]

corr_rr_Ns = [ ['', True, 'corr_rr_Ns', 1.] ]
        
covarnce_rr_Ns = [ ['', True, 'covarnce_rr_Ns', 1.] ]

corr_rr_rg = [ ['', True, 'corr_rr_rg', 1.] ]
        
covarnce_rr_rg = [ ['', True, 'covarnce_rr_rg', 1.] ]

corr_rr_Ng = [ ['', True, 'corr_rr_Ng', 1.] ]
        
covarnce_rr_Ng = [ ['', True, 'covarnce_rr_Ng', 1.] ]
        
corr_Nr_ri = [ ['', True, 'corr_Nr_ri', 1.] ]
        
covarnce_Nr_ri = [ ['', True, 'covarnce_Nr_ri', 1.] ]
        
corr_Nr_Ni = [ ['', True, 'corr_Nr_Ni', 1.] ]
        
covarnce_Nr_Ni = [ ['', True, 'covarnce_Nr_Ni', 1.] ]

corr_Nr_rs = [ ['', True, 'corr_Nr_rs', 1.] ]
        
covarnce_Nr_rs = [ ['', True, 'covarnce_Nr_rs', 1.] ]

corr_Nr_Ns = [ ['', True, 'corr_Nr_Ns', 1.] ]
        
covarnce_Nr_Ns = [ ['', True, 'covarnce_Nr_Ns', 1.] ]

corr_Nr_rg = [ ['', True, 'corr_Nr_rg', 1.] ]
        
covarnce_Nr_rg = [ ['', True, 'covarnce_Nr_rg', 1.] ]

corr_Nr_Ng = [ ['', True, 'corr_Nr_Ng', 1.] ]
        
covarnce_Nr_Ng = [ ['', True, 'covarnce_Nr_Ng', 1.] ]
        
corr_ri_Ni = [ ['', True, 'corr_ri_Ni', 1.] ]
        
covarnce_ri_Ni = [ ['', True, 'covarnce_ri_Ni', 1.] ]

corr_ri_rs = [ ['', True, 'corr_ri_rs', 1.] ]
        
covarnce_ri_rs = [ ['', True, 'covarnce_ri_rs', 1.] ]

corr_ri_Ns = [ ['', True, 'corr_ri_Ns', 1.] ]
        
covarnce_ri_Ns = [ ['', True, 'covarnce_ri_Ns', 1.] ]

corr_ri_rg = [ ['', True, 'corr_ri_rg', 1.] ]
        
covarnce_ri_rg = [ ['', True, 'covarnce_ri_rg', 1.] ]

corr_ri_Ng = [ ['', True, 'corr_ri_Ng', 1.] ]
        
covarnce_ri_Ng = [ ['', True, 'covarnce_ri_Ng', 1.] ]

corr_Ni_rs = [ ['', True, 'corr_Ni_rs', 1.] ]
        
covarnce_Ni_rs = [ ['', True, 'covarnce_Ni_rs', 1.] ]

corr_Ni_Ns = [ ['', True, 'corr_Ni_Ns', 1.] ]
        
covarnce_Ni_Ns = [ ['', True, 'covarnce_Ni_Ns', 1.] ]

corr_Ni_rg = [ ['', True, 'corr_Ni_rg', 1.] ]
        
covarnce_Ni_rg = [ ['', True, 'covarnce_Ni_rg', 1.] ]

corr_Ni_Ng = [ ['', True, 'corr_Ni_Ng', 1.] ]
        
covarnce_Ni_Ng = [ ['', True, 'covarnce_Ni_Ng', 1.] ]

corr_rs_Ns = [ ['', True, 'corr_rs_Ns', 1.] ]
        
covarnce_rs_Ns = [ ['', True, 'covarnce_rs_Ns', 1.] ]

corr_rs_rg = [ ['', True, 'corr_rs_rg', 1.] ]
        
covarnce_rs_rg = [ ['', True, 'covarnce_rs_rg', 1.] ]

corr_rs_Ng = [ ['', True, 'corr_rs_Ng', 1.] ]
        
covarnce_rs_Ng = [ ['', True, 'covarnce_rs_Ng', 1.] ]

corr_Ns_rg = [ ['', True, 'corr_Ns_rg', 1.] ]
        
covarnce_Ns_rg = [ ['', True, 'covarnce_Ns_rg', 1.] ]

corr_Ns_Ng = [ ['', True, 'corr_Ns_Ng', 1.] ]
        
covarnce_Ns_Ng = [ ['', True, 'covarnce_Ns_Ng', 1.] ]

corr_rg_Ng = [ ['', True, 'corr_rg_Ng', 1.] ]
        
covarnce_rg_Ng = [ ['', True, 'covarnce_rg_Ng', 1.] ]


lines = [corr_chi_w, covarnce_chi_w, corr_chi_rr, covarnce_chi_rr, corr_chi_Nr, covarnce_chi_Nr, corr_chi_ri, covarnce_chi_ri, corr_chi_Ni, covarnce_chi_Ni, corr_chi_rs, covarnce_chi_rs, corr_chi_Ns, covarnce_chi_Ns, corr_chi_rg, covarnce_chi_rg, corr_chi_Ng, covarnce_chi_Ng,\
        corr_w_rr, covarnce_w_rr, corr_w_Nr, covarnce_w_Nr, corr_w_ri, covarnce_w_ri, corr_w_Ni, covarnce_w_Ni, corr_w_rs, covarnce_w_rs, corr_w_Ns, covarnce_w_Ns, corr_w_rg, covarnce_w_rg, corr_w_Ng, covarnce_w_Ng,\
        corr_rr_Nr, covarnce_rr_Nr, corr_rr_ri, covarnce_rr_ri, corr_rr_Ni, covarnce_rr_Ni, corr_rr_rs, covarnce_rr_rs, corr_rr_Ns, covarnce_rr_Ns, corr_rr_rg, covarnce_rr_rg, corr_rr_Ng, covarnce_rr_Ng,\
        corr_Nr_ri, covarnce_Nr_ri, corr_Nr_Ni, covarnce_Nr_Ni, corr_Nr_rs, covarnce_Nr_rs, corr_Nr_Ns, covarnce_Nr_Ns, corr_Nr_rg, covarnce_Nr_rg, corr_Nr_Ng, covarnce_Nr_Ng,\
        corr_ri_Ni, covarnce_ri_Ni, corr_ri_rs, covarnce_ri_rs, corr_ri_Ns, covarnce_ri_Ns, corr_ri_rg, covarnce_ri_rg, corr_ri_Ng, covarnce_ri_Ng,\
        corr_Ni_rs, covarnce_Ni_rs, corr_Ni_Ns, covarnce_Ni_Ns, corr_Ni_rg, covarnce_Ni_rg, corr_Ni_Ng, covarnce_Ni_Ng,\
        corr_rs_Ns, covarnce_rs_Ns, corr_rs_rg, covarnce_rs_rg, corr_rs_Ng, covarnce_rs_Ng,\
        corr_Ns_rg, covarnce_Ns_rg, corr_Ns_Ng, covarnce_Ns_Ng,\
        corr_rg_Ng, covarnce_rg_Ng]
