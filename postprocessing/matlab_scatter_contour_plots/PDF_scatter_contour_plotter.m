% $Id$
function PDF_scatter_contour_plotter( input_file_sam, input_file_clubb )

% SAM LES 3D file variables
global idx_3D_w
global idx_3D_rr
% CLUBB variables
global idx_w_1
global idx_w_2
global idx_mu_rr_1_n
global idx_mu_rr_2_n
global idx_varnce_w_1
global idx_varnce_w_2
global idx_sigma_rr_1_n
global idx_sigma_rr_2_n
global idx_corr_w_rr_1_n
global idx_corr_w_rr_2_n
global idx_mixt_frac
global idx_precip_frac_1
global idx_precip_frac_2

% SAM LES 3D NetCDF filename.
if ( strcmp( input_file_sam, 'default' ) )
   filename_sam = '/home/griffinb/hoc_v2.2_tuner/output/LES_output/RICO_128x128x100_drizzle_128_0000255600.nc';
else
   filename_sam = input_file_sam;
end

% CLUBB zt NetCDF filename.
if ( strcmp( input_file_clubb, 'default' ) )
   filename_clubb = '/home/griffinb/hoc_v2.2_tuner/output/rico_zt.nc';
else
   filename_clubb = input_file_clubb;
end

% CLUBB's time output at index 852 is at a time of 4260 minutes.
clubb_time_idx = 852;

% Read SAM NetCDF file and obtain variables.
[ z_sam, time_sam, var_sam, units_corrector_type_sam, ...
  nx_sam, ny_sam, nz_sam, num_t_sam, num_var_sam ] ...
= read_SAM_3D_file( filename_sam );

% Read CLUBB zt NetCDF file and obtain variables.
[ z_clubb, time_clubb, var_clubb, units_corrector_type_clubb, ...
  nz_clubb, num_t_clubb, num_var_clubb ] ...
= read_CLUBB_file( filename_clubb );

% Use appropriate units (SI units).
[ var_sam ] ...
   = unit_corrector( num_var_sam, var_sam, units_corrector_type_sam, -1 );
[ var_clubb ] ...
   = unit_corrector( num_var_clubb, var_clubb, ...
                     units_corrector_type_clubb, -1 );

% SAM output is taken from 1500 m.  CLUBB output is taken from 1524 m.
sam_height_idx = 38;
clubb_height_idx = 19;

% Place the SAM variables from the same location in the same 1-D index for use
% in a scatterplot.
var1 = zeros( nx_sam*ny_sam, 1 );
var2 = zeros( nx_sam*ny_sam, 1 );
for i = 1:1:nx_sam
   for j = 1:1:ny_sam
      var1((i-1)*nx_sam+j) = var_sam(idx_3D_w,i,j,sam_height_idx,1);
      var2((i-1)*nx_sam+j) = var_sam(idx_3D_rr,i,j,sam_height_idx,1);
   end
end

% Unpack CLUBB variables (PDF parameters).
mu_w_1 = var_clubb( idx_w_1, 1, 1, clubb_height_idx, clubb_time_idx );
mu_w_2 = var_clubb( idx_w_2, 1, 1, clubb_height_idx, clubb_time_idx );
mu_rr_1_n = var_clubb( idx_mu_rr_1_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
mu_rr_2_n = var_clubb( idx_mu_rr_2_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
sigma_w_1 = sqrt( var_clubb( idx_varnce_w_1, 1, 1, ...
                             clubb_height_idx, clubb_time_idx ) );
sigma_w_2 = sqrt( var_clubb( idx_varnce_w_2, 1, 1, ...
                             clubb_height_idx, clubb_time_idx ) );
sigma_rr_1_n = var_clubb( idx_sigma_rr_1_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );
sigma_rr_2_n = var_clubb( idx_sigma_rr_2_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );
corr_w_rr_1_n = var_clubb( idx_corr_w_rr_1_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_w_rr_2_n = var_clubb( idx_corr_w_rr_2_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
mixt_frac = var_clubb( idx_mixt_frac, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
precip_frac_1 = var_clubb( idx_precip_frac_1, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
precip_frac_2 = var_clubb( idx_precip_frac_2, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
                       
% Plot the CLUBB PDF and LES points.
plot_CLUBB_PDF_LES_pts_NL( var1, var2, nx_sam, ...
                           ny_sam, 100, 100, ...
                           mu_w_1, mu_w_2, mu_rr_1_n, mu_rr_2_n, ...
                           sigma_w_1, sigma_w_2, sigma_rr_1_n, ...
                           sigma_rr_2_n, corr_w_rr_1_n, ...
                           corr_w_rr_2_n, precip_frac_1, ...
                           precip_frac_2, mixt_frac, ...
                           'w    [m/s]', 'r_{r}    [kg/kg]', ...
                           '\bf RICO', 'w vs. r_{r}', ...
                           'Time = 4260 minutes', ...
                           'Altitude = 1524 meters', ...
                           'Input fields (predictive fields)' )
