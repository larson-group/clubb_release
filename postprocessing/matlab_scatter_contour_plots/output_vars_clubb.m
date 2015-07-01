% $Id$
function [ varname_clubb, units_corrector_type, ...
           num_var_clubb, num_tot_var_clubb ] = output_vars_clubb

idx_count = 0;

%============= Add Model Field Variables Here =============================

% w_1 (m/s).
idx_count = idx_count + 1;
global idx_w_1
idx_w_1 = idx_count;
units_corrector_type(idx_w_1) = 0;
varname_clubb(idx_w_1,1:3) = 'w_1';

% w_2 (m/s).
idx_count = idx_count + 1;
global idx_w_2
idx_w_2 = idx_count;
units_corrector_type(idx_w_2) = 0;
varname_clubb(idx_w_2,1:3) = 'w_2';

% rt_1 (kg/kg).
idx_count = idx_count + 1;
global idx_rt_1
idx_rt_1 = idx_count;
units_corrector_type(idx_rt_1) = 0;
varname_clubb(idx_rt_1,1:4) = 'rt_1';

% rt_2 (kg/kg).
idx_count = idx_count + 1;
global idx_rt_2
idx_rt_2 = idx_count;
units_corrector_type(idx_rt_2) = 0;
varname_clubb(idx_rt_2,1:4) = 'rt_2';

% thl_1 (K).
idx_count = idx_count + 1;
global idx_thl_1
idx_thl_1 = idx_count;
units_corrector_type(idx_thl_1) = 0;
varname_clubb(idx_thl_1,1:5) = 'thl_1';

% thl_2 (K).
idx_count = idx_count + 1;
global idx_thl_2
idx_thl_2 = idx_count;
units_corrector_type(idx_thl_2) = 0;
varname_clubb(idx_thl_2,1:5) = 'thl_2';

% chi_1 (kg/kg).
idx_count = idx_count + 1;
global idx_chi_1
idx_chi_1 = idx_count;
units_corrector_type(idx_chi_1) = 0;
varname_clubb(idx_chi_1,1:5) = 'chi_1';

% chi_2 (kg/kg).
idx_count = idx_count + 1;
global idx_chi_2
idx_chi_2 = idx_count;
units_corrector_type(idx_chi_2) = 0;
varname_clubb(idx_chi_2,1:5) = 'chi_2';

% mu_rr_1_n (ln(kg/kg)).
idx_count = idx_count + 1;
global idx_mu_rr_1_n
idx_mu_rr_1_n = idx_count;
units_corrector_type(idx_mu_rr_1_n) = 0;
varname_clubb(idx_mu_rr_1_n,1:9) = 'mu_rr_1_n';

% mu_rr_2_n (ln(kg/kg)).
idx_count = idx_count + 1;
global idx_mu_rr_2_n
idx_mu_rr_2_n = idx_count;
units_corrector_type(idx_mu_rr_2_n) = 0;
varname_clubb(idx_mu_rr_2_n,1:9) = 'mu_rr_2_n';

% mu_Nr_1_n (ln(num/kg)).
idx_count = idx_count + 1;
global idx_mu_Nr_1_n
idx_mu_Nr_1_n = idx_count;
units_corrector_type(idx_mu_Nr_1_n) = 0;
varname_clubb(idx_mu_Nr_1_n,1:9) = 'mu_Nr_1_n';

% mu_Nr_2_n (ln(num/kg)).
idx_count = idx_count + 1;
global idx_mu_Nr_2_n
idx_mu_Nr_2_n = idx_count;
units_corrector_type(idx_mu_Nr_2_n) = 0;
varname_clubb(idx_mu_Nr_2_n,1:9) = 'mu_Nr_2_n';

% varnce_w_1 (m^2/s^2).
idx_count = idx_count + 1;
global idx_varnce_w_1
idx_varnce_w_1 = idx_count;
units_corrector_type(idx_varnce_w_1) = 0;
varname_clubb(idx_varnce_w_1,1:10) = 'varnce_w_1';

% varnce_w_2 (m^2/s^2).
idx_count = idx_count + 1;
global idx_varnce_w_2
idx_varnce_w_2 = idx_count;
units_corrector_type(idx_varnce_w_2) = 0;
varname_clubb(idx_varnce_w_2,1:10) = 'varnce_w_2';

% varnce_rt_1 (kg^2/kg^2).
idx_count = idx_count + 1;
global idx_varnce_rt_1
idx_varnce_rt_1 = idx_count;
units_corrector_type(idx_varnce_rt_1) = 0;
varname_clubb(idx_varnce_rt_1,1:11) = 'varnce_rt_1';

% varnce_rt_2 (kg^2/kg^2).
idx_count = idx_count + 1;
global idx_varnce_rt_2
idx_varnce_rt_2 = idx_count;
units_corrector_type(idx_varnce_rt_2) = 0;
varname_clubb(idx_varnce_rt_2,1:11) = 'varnce_rt_2';

% varnce_thl_1 (K^2).
idx_count = idx_count + 1;
global idx_varnce_thl_1
idx_varnce_thl_1 = idx_count;
units_corrector_type(idx_varnce_thl_1) = 0;
varname_clubb(idx_varnce_thl_1,1:12) = 'varnce_thl_1';

% varnce_thl_2 (K^2).
idx_count = idx_count + 1;
global idx_varnce_thl_2
idx_varnce_thl_2 = idx_count;
units_corrector_type(idx_varnce_thl_2) = 0;
varname_clubb(idx_varnce_thl_2,1:12) = 'varnce_thl_2';

% stdev_chi_1 (kg/kg).
idx_count = idx_count + 1;
global idx_stdev_chi_1
idx_stdev_chi_1 = idx_count;
units_corrector_type(idx_stdev_chi_1) = 0;
varname_clubb(idx_stdev_chi_1,1:11) = 'stdev_chi_1';

% stdev_chi_2 (kg/kg).
idx_count = idx_count + 1;
global idx_stdev_chi_2
idx_stdev_chi_2 = idx_count;
units_corrector_type(idx_stdev_chi_2) = 0;
varname_clubb(idx_stdev_chi_2,1:11) = 'stdev_chi_2';

% stdev_eta_1 (kg/kg).
idx_count = idx_count + 1;
global idx_stdev_eta_1
idx_stdev_eta_1 = idx_count;
units_corrector_type(idx_stdev_eta_1) = 0;
varname_clubb(idx_stdev_eta_1,1:11) = 'stdev_eta_1';

% stdev_eta_2 (kg/kg).
idx_count = idx_count + 1;
global idx_stdev_eta_2
idx_stdev_eta_2 = idx_count;
units_corrector_type(idx_stdev_eta_2) = 0;
varname_clubb(idx_stdev_eta_2,1:11) = 'stdev_eta_2';

% sigma_rr_1_n (-).
idx_count = idx_count + 1;
global idx_sigma_rr_1_n
idx_sigma_rr_1_n = idx_count;
units_corrector_type(idx_sigma_rr_1_n) = 0;
varname_clubb(idx_sigma_rr_1_n,1:12) = 'sigma_rr_1_n';

% sigma_rr_2_n (-).
idx_count = idx_count + 1;
global idx_sigma_rr_2_n
idx_sigma_rr_2_n = idx_count;
units_corrector_type(idx_sigma_rr_2_n) = 0;
varname_clubb(idx_sigma_rr_2_n,1:12) = 'sigma_rr_2_n';

% sigma_Nr_1_n (-).
idx_count = idx_count + 1;
global idx_sigma_Nr_1_n
idx_sigma_Nr_1_n = idx_count;
units_corrector_type(idx_sigma_Nr_1_n) = 0;
varname_clubb(idx_sigma_Nr_1_n,1:12) = 'sigma_Nr_1_n';

% sigma_Nr_2_n (-).
idx_count = idx_count + 1;
global idx_sigma_Nr_2_n
idx_sigma_Nr_2_n = idx_count;
units_corrector_type(idx_sigma_Nr_2_n) = 0;
varname_clubb(idx_sigma_Nr_2_n,1:12) = 'sigma_Nr_2_n';

% Correlation (within-component) of rt and thl (both PDF comps.), rrtthl (-).
idx_count = idx_count + 1;
global idx_corr_rt_thl
idx_corr_rt_thl = idx_count;
units_corrector_type(idx_corr_rt_thl) = 0;
varname_clubb(idx_corr_rt_thl,1:6) = 'rrtthl';

% corr_chi_eta_1_ca (-).
idx_count = idx_count + 1;
global idx_corr_chi_eta_1_ca
idx_corr_chi_eta_1_ca = idx_count;
units_corrector_type(idx_corr_chi_eta_1_ca) = 0;
varname_clubb(idx_corr_chi_eta_1_ca,1:17) = 'corr_chi_eta_1_ca';

% corr_chi_eta_2_ca (-).
idx_count = idx_count + 1;
global idx_corr_chi_eta_2_ca
idx_corr_chi_eta_2_ca = idx_count;
units_corrector_type(idx_corr_chi_eta_2_ca) = 0;
varname_clubb(idx_corr_chi_eta_2_ca,1:17) = 'corr_chi_eta_2_ca';

% corr_w_rr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_w_rr_1_n
idx_corr_w_rr_1_n = idx_count;
units_corrector_type(idx_corr_w_rr_1_n) = 0;
varname_clubb(idx_corr_w_rr_1_n,1:13) = 'corr_w_rr_1_n';

% corr_w_rr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_w_rr_2_n
idx_corr_w_rr_2_n = idx_count;
units_corrector_type(idx_corr_w_rr_2_n) = 0;
varname_clubb(idx_corr_w_rr_2_n,1:13) = 'corr_w_rr_2_n';

% corr_w_Nr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_w_Nr_1_n
idx_corr_w_Nr_1_n = idx_count;
units_corrector_type(idx_corr_w_Nr_1_n) = 0;
varname_clubb(idx_corr_w_Nr_1_n,1:13) = 'corr_w_Nr_1_n';

% corr_w_Nr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_w_Nr_2_n
idx_corr_w_Nr_2_n = idx_count;
units_corrector_type(idx_corr_w_Nr_2_n) = 0;
varname_clubb(idx_corr_w_Nr_2_n,1:13) = 'corr_w_Nr_2_n';

% corr_chi_rr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_chi_rr_1_n
idx_corr_chi_rr_1_n = idx_count;
units_corrector_type(idx_corr_chi_rr_1_n) = 0;
varname_clubb(idx_corr_chi_rr_1_n,1:15) = 'corr_chi_rr_1_n';

% corr_chi_rr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_chi_rr_2_n
idx_corr_chi_rr_2_n = idx_count;
units_corrector_type(idx_corr_chi_rr_2_n) = 0;
varname_clubb(idx_corr_chi_rr_2_n,1:15) = 'corr_chi_rr_2_n';

% corr_chi_Nr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_chi_Nr_1_n
idx_corr_chi_Nr_1_n = idx_count;
units_corrector_type(idx_corr_chi_Nr_1_n) = 0;
varname_clubb(idx_corr_chi_Nr_1_n,1:15) = 'corr_chi_Nr_1_n';

% corr_chi_Nr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_chi_Nr_2_n
idx_corr_chi_Nr_2_n = idx_count;
units_corrector_type(idx_corr_chi_Nr_2_n) = 0;
varname_clubb(idx_corr_chi_Nr_2_n,1:15) = 'corr_chi_Nr_2_n';

% corr_eta_rr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_eta_rr_1_n
idx_corr_eta_rr_1_n = idx_count;
units_corrector_type(idx_corr_eta_rr_1_n) = 0;
varname_clubb(idx_corr_eta_rr_1_n,1:15) = 'corr_eta_rr_1_n';

% corr_eta_rr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_eta_rr_2_n
idx_corr_eta_rr_2_n = idx_count;
units_corrector_type(idx_corr_eta_rr_2_n) = 0;
varname_clubb(idx_corr_eta_rr_2_n,1:15) = 'corr_eta_rr_2_n';

% corr_eta_Nr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_eta_Nr_1_n
idx_corr_eta_Nr_1_n = idx_count;
units_corrector_type(idx_corr_eta_Nr_1_n) = 0;
varname_clubb(idx_corr_eta_Nr_1_n,1:15) = 'corr_eta_Nr_1_n';

% corr_eta_Nr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_eta_Nr_2_n
idx_corr_eta_Nr_2_n = idx_count;
units_corrector_type(idx_corr_eta_Nr_2_n) = 0;
varname_clubb(idx_corr_eta_Nr_2_n,1:15) = 'corr_eta_Nr_2_n';

% corr_rr_Nr_1_n (-).
idx_count = idx_count + 1;
global idx_corr_rr_Nr_1_n
idx_corr_rr_Nr_1_n = idx_count;
units_corrector_type(idx_corr_rr_Nr_1_n) = 0;
varname_clubb(idx_corr_rr_Nr_1_n,1:14) = 'corr_rr_Nr_1_n';

% corr_rr_Nr_2_n (-).
idx_count = idx_count + 1;
global idx_corr_rr_Nr_2_n
idx_corr_rr_Nr_2_n = idx_count;
units_corrector_type(idx_corr_rr_Nr_2_n) = 0;
varname_clubb(idx_corr_rr_Nr_2_n,1:14) = 'corr_rr_Nr_2_n';

% mixt_frac (-).
idx_count = idx_count + 1;
global idx_mixt_frac
idx_mixt_frac = idx_count;
units_corrector_type(idx_mixt_frac) = 0;
varname_clubb(idx_mixt_frac,1:9) = 'mixt_frac';

% precip_frac_1 (-).
idx_count = idx_count + 1;
global idx_precip_frac_1
idx_precip_frac_1 = idx_count;
units_corrector_type(idx_precip_frac_1) = 0;
varname_clubb(idx_precip_frac_1,1:13) = 'precip_frac_1';

% precip_frac_2 (-).
idx_count = idx_count + 1;
global idx_precip_frac_2
idx_precip_frac_2 = idx_count;
units_corrector_type(idx_precip_frac_2) = 0;
varname_clubb(idx_precip_frac_2,1:13) = 'precip_frac_2';

% sigma_sqd_w (-).
idx_count = idx_count + 1;
global idx_sigma_sqd_w
idx_sigma_sqd_w = idx_count;
units_corrector_type(idx_sigma_sqd_w) = 0;
varname_clubb(idx_sigma_sqd_w,1:14) = 'sigma_sqd_w_zt';

num_var_clubb = idx_count;

%============= Add Model Coordinates (Height and Time) Here ===============

% Altitude (meters)
idx_count = idx_count + 1;
global idx_z
idx_z = idx_count;
units_corrector_type(idx_z) = 0;
varname_clubb(idx_z,1:8) = 'altitude';

% Elapsed time (minutes)
idx_count = idx_count + 1;
global idx_time
idx_time = idx_count;
units_corrector_type(idx_time) = 0;
varname_clubb(idx_time,1:4) = 'time';

num_tot_var_clubb = idx_count;
