import netCDF4
import numpy as np

case_name = 'rico_lh'
num_categories = 8

nc        = netCDF4.Dataset(case_name + '_zt.nc')
lh_nc     = netCDF4.Dataset(case_name + '_lh_zt.nc')
lh_sfc_nc = netCDF4.Dataset(case_name + '_lh_sfc.nc')
u_2D      = netCDF4.Dataset(case_name + '_u_lh_sample_points_2D.nc')
nl_2D     = netCDF4.Dataset(case_name + '_nl_lh_sample_points_2D.nc')

def compute_category_real_probs(t,k):
  cloud_frac_1  = nc.variables['cloud_frac_1']
  cloud_frac_2  = nc.variables['cloud_frac_2']
  precip_frac_1 = nc.variables['precip_frac_1']
  precip_frac_2 = nc.variables['precip_frac_2']
  mixt_frac     = nc.variables['mixt_frac']

  cf1 = cloud_frac_1 [t,k,0,0]
  cf2 = cloud_frac_2 [t,k,0,0]
  pf1 = precip_frac_1[t,k,0,0]
  pf2 = precip_frac_2[t,k,0,0]
  mf1 = mixt_frac    [t,k,0,0]

  category_real_probs = np.empty(num_categories)
  category_real_probs[0] = cf1 * pf1 * mf1
  category_real_probs[1] = cf2 * pf2 * (1-mf1)
  category_real_probs[2] = (1-cf1) * pf1 * mf1
  category_real_probs[3] = (1-cf2) * pf2 * (1-mf1)
  category_real_probs[4] = cf1 * (1-pf1) * mf1
  category_real_probs[5] = cf2 * (1-pf2) * (1-mf1)
  category_real_probs[6] = (1-cf1) * (1-pf1) * mf1
  category_real_probs[7] = (1-cf2) * (1-pf2) * (1-mf1)

  return category_real_probs

def which_category(t,k,i):
  chi = nl_2D.variables['chi']        [t,k,i,0]
  X_m =  u_2D.variables['X_mixt_comp'][t,k,i,0]
  rrl = nl_2D.variables['rr']         [t,k,i,0]

  if chi < 0.0:
    l_in_cloud = False
  else:
    l_in_cloud = True

  if rrl > 0.0:
    l_in_precip = True
  else:
    l_in_precip = False

  if X_m == 1:
    l_in_component_1 = True
  else:
    l_in_component_1 = False

  if l_in_cloud and l_in_precip and l_in_component_1:
    category = 1
  elif l_in_cloud and l_in_precip and not l_in_component_1:
    category = 2
  elif not l_in_cloud and l_in_precip and l_in_component_1:
    category = 3
  elif not l_in_cloud and l_in_precip and not l_in_component_1:
    category = 4
  elif l_in_cloud and not l_in_precip and l_in_component_1:
    category = 5
  elif l_in_cloud and not l_in_precip and not l_in_component_1:
    category = 6
  elif not l_in_cloud and not l_in_precip and l_in_component_1:
    category = 7
  elif not l_in_cloud and not l_in_precip and not l_in_component_1:
    category = 8
  else:
    print("Error determining category.")
    exit(1)

  return category - 1

# Choose a height level and timestep to perform analysis on
t = 3000
k = int(lh_sfc_nc.variables['k_lh_start'][t,0,0,0]) - 1

# Determine number of sample points
num_samples = len(u_2D.dimensions['latitude'])

# Some variables of interest:
rho = nc.variables['rho'][t,k,0,0]

# Mean in each category
mean_cat  = np.zeros(num_categories)
# Loop over sample points
for i in range(0,num_samples):
    # Calculate value of the variable
    sample_weight = u_2D.variables['lh_sample_point_weights'][t,k,i,0]
    chi = nl_2D.variables['chi'][t,k,i,0]
    Ncn = nl_2D.variables['Ncn'][t,k,i,0]

#   if chi > 0.0:
#       var_value = 1350.0 * ((rho / 1.0e6) ** -1.79) * \
#           (chi ** 2.47) * (Ncn ** -1.79)
#   else:
#       var_value = 0.0
    var_value = chi

    icat = which_category(t,k,i)

    mean_cat[icat]  = mean_cat[icat]  + sample_weight * (var_value ** 2)

mean_cat = np.divide(mean_cat, float(num_samples))

category_real_probs = compute_category_real_probs(t,k)
mean_cat = np.sqrt(np.multiply(mean_cat, category_real_probs))

variance = mean_cat

#################################################################
# Determine Fortran variance as sum as variance in each category
fortran_variance = np.empty(num_categories)

for i in range(0,num_categories):
    variance_var = lh_nc.variables['silhs_var_cat_'+str(i+1)]
    fortran_variance[i] = variance_var[t,k,0,0]

# Finally, compare calculated to what is contained in the file
print('Fortran: ' + str(fortran_variance))
print('Python:  ' + str(variance))
