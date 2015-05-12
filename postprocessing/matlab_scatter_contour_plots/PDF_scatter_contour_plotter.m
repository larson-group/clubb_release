% $Id$
function PDF_scatter_contour_plotter( input_file_sam, input_file_clubb )

% SAM LES 3D NetCDF filename.
if ( strcmp( input_file_sam, 'default' ) )
   filename_sam = '../../output/LES_output/RICO_128x128x100_drizzle_128_0000255600_micro.nc';
else
   filename_sam = input_file_sam;
end

% CLUBB zt NetCDF filename.
if ( strcmp( input_file_clubb, 'default' ) )
   filename_clubb = '../../output/rico_zt.nc';
else
   filename_clubb = input_file_clubb;
end

% Declare the CLUBB vertical level index.
% CLUBB contour and line plots will be generated from PDF parameters from
% this level.  SAM LES data will be interpolated to the altitude of this
% level and then displayed in scatterplots and histograms.
clubb_height_idx = 19;
%clubb_height_idx = 59;

% Information to be printed on the plots.
casename = 'RICO';
%casename = 'RF02';
print_note = 'Input fields (predictive fields)';

% Select the plots that will be plotted.
plot_rt_thl  = true;
plot_chi_eta = true;
plot_w_rr    = true;
plot_w_Nr    = true;
plot_chi_rr  = true;
plot_chi_Nr  = true;
plot_eta_rr  = true;
plot_eta_Nr  = true;
plot_rr_Nr   = true;
plot_w       = true;
plot_ln_rr   = true;
plot_ln_Nr   = true;

%==========================================================================

% SAM LES 3D file variables
global idx_3D_w
global idx_3D_rt
global idx_3D_thl
global idx_3D_chi
global idx_3D_eta
global idx_3D_rr
global idx_3D_Nr
% CLUBB variables
global idx_w_1
global idx_w_2
global idx_rt_1
global idx_rt_2
global idx_thl_1
global idx_thl_2
global idx_chi_1
global idx_chi_2
global idx_mu_rr_1_n
global idx_mu_rr_2_n
global idx_mu_Nr_1_n
global idx_mu_Nr_2_n
global idx_varnce_w_1
global idx_varnce_w_2
global idx_varnce_rt_1
global idx_varnce_rt_2
global idx_varnce_thl_1
global idx_varnce_thl_2
global idx_stdev_chi_1
global idx_stdev_chi_2
global idx_stdev_eta_1
global idx_stdev_eta_2
global idx_sigma_rr_1_n
global idx_sigma_rr_2_n
global idx_sigma_Nr_1_n
global idx_sigma_Nr_2_n
global idx_corr_rt_thl
global idx_corr_chi_eta_1_ca
global idx_corr_chi_eta_2_ca
global idx_corr_w_rr_1_n
global idx_corr_w_rr_2_n
global idx_corr_w_Nr_1_n
global idx_corr_w_Nr_2_n
global idx_corr_chi_rr_1_n
global idx_corr_chi_rr_2_n
global idx_corr_chi_Nr_1_n
global idx_corr_chi_Nr_2_n
global idx_corr_eta_rr_1_n
global idx_corr_eta_rr_2_n
global idx_corr_eta_Nr_1_n
global idx_corr_eta_Nr_2_n
global idx_corr_rr_Nr_1_n
global idx_corr_rr_Nr_2_n
global idx_mixt_frac
global idx_precip_frac_1
global idx_precip_frac_2
global idx_sigma_sqd_w

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

% Find the time in the CLUBB zt output file that is equal (or closest) to
% the SAM LES output time.
time_sam_sec = time_sam * 86400.0;
time_diff_clubb_sam_sec = abs( time_clubb - time_sam_sec );
% Initialize the minimum time difference and its (CLUBB) index.
idx_min_time_diff = 1;
min_time_diff = time_diff_clubb_sam_sec(idx_min_time_diff);
% Find the index of the minimum time difference between CLUBB output time
% and the requested SAM 3D file output time.
for iter = 2:1:num_t_clubb
   if ( time_diff_clubb_sam_sec(iter) < min_time_diff )
      min_time_diff = time_diff_clubb_sam_sec(iter);
      idx_min_time_diff = iter;
   end
end
% The CLUBB output index is the index that corresponds with the minimum
% time difference between the CLUBB output time and the SAM LES 3D file
% output time.
clubb_time_idx = idx_min_time_diff;

% Print the time of the SAM LES output.
fprintf( 'Time of SAM LES output (seconds): %g\n', time_sam_sec );

% Print the CLUBB output time index and the associated time.
fprintf( [ 'Time index of CLUBB output: %g;', ...
           ' Time of CLUBB output (seconds): %g\n' ], ...
           clubb_time_idx, time_clubb(clubb_time_idx) );

% Print the altitude at the CLUBB vertical level index.
fprintf( 'Altitude of CLUBB zt grid level (meters): %g\n', ...
         z_clubb(clubb_height_idx) );

% Place the SAM variables from the same location in the same 1-D index for
% use in a scatterplot.
sam_var_lev = zeros( num_var_sam, nx_sam*ny_sam );

if ( z_clubb(clubb_height_idx) > z_sam(nz_sam) )

   % The height of the CLUBB grid level is above the highest SAM LES grid
   % level.  Use SAM values from the highest SAM LES grid level.
   fprintf( [ 'The altitude of the CLUBB zt grid level is higher ', ...
              'than the highest SAM LES grid level.  The highest ', ...
              'SAM LES grid level will be used.\n' ] );
   fprintf( 'Altitude of SAM LES grid level (meters): %g\n', ...
            z_sam(nz_sam) );

   for i = 1:1:nx_sam
      for j = 1:1:ny_sam
         for idx_var = 1:1:num_var_sam
            sam_var_lev(idx_var,(i-1)*ny_sam+j) ...
            = var_sam(idx_var,i,j,nz_sam,1);
         end % idx_var = 1:1:num_var_sam
      end % j = 1:1:ny_sam
   end % i = 1:1:nx_sam

elseif ( z_clubb(clubb_height_idx) < z_sam(1) )
         
   % The height of the CLUBB grid level is below the lowest SAM LES grid
   % level.  Use SAM values from the lowest SAM LES grid level.
   fprintf( [ 'The altitude of the CLUBB zt grid level is lower ', ...
              'than the lowest SAM LES grid level.  The lowest ', ...
              'SAM LES grid level will be used.\n' ] );
   fprintf( 'Altitude of SAM LES grid level (meters):  %g\n', ...
            z_sam(1) );

   for i = 1:1:nx_sam
      for j = 1:1:ny_sam
         for idx_var = 1:1:num_var_sam
            sam_var_lev(idx_var,(i-1)*ny_sam+j) = var_sam(idx_var,i,j,1,1);
         end % idx_var = 1:1:num_var_sam
      end % j = 1:1:ny_sam
   end % i = 1:1:nx_sam

else % z_sam(1) <= z_clubb(clubb_height_idx) <= z_sam(nz_sam)

   % The height of the CLUBB grid level is found within the SAM LES
   % vertical domain.
   exact_lev_idx = -1;
   lower_lev_idx = -1;
   upper_lev_idx = -1;
   for k = 1:1:nz_sam

      if ( z_sam(k) == z_clubb(clubb_height_idx) )

         % The SAM LES grid level is at the exact same altitude as the
         % requested CLUBB grid level.
         exact_lev_idx = k;
         break

      elseif ( z_sam(k) < z_clubb(clubb_height_idx) )

         % The SAM LES grid level is below the requested CLUBB grid level.
         lower_lev_idx = k;

      else % z_sam(k) > z_clubb(clubb_height_idx)

         % The SAM LES grid level is above the requested CLUBB grid level.
         upper_lev_idx = k;

      end

      if ( upper_lev_idx == lower_lev_idx + 1 )
         break
      end

   end % k = 1:1:nz_sam

   if ( exact_lev_idx > 0 )

      fprintf( [ 'The altitude of the SAM LES grid level is the same ', ...
                 'as the CLUBB zt grid level.\n' ] );

      for i = 1:1:nx_sam
         for j = 1:1:ny_sam
            for idx_var = 1:1:num_var_sam
               sam_var_lev(idx_var,(i-1)*ny_sam+j) ...
               = var_sam(idx_var,i,j,exact_lev_idx,1);
            end % idx_var = 1:1:num_var_sam
         end % j = 1:1:ny_sam
      end % i = 1:1:nx_sam

   else % interpolate between two levels.
 
      interp_weight ...
      = ( z_clubb(clubb_height_idx) - z_sam(lower_lev_idx) ) ...
        / ( z_sam(upper_lev_idx) - z_sam(lower_lev_idx) );

      fprintf( [ 'The altitude of the CLUBB zt grid level is between ', ...
                 'two SAM LES grid levels.\n' ] );
      fprintf( [ 'Altitude of the SAM LES grid level above the ', ...
                 'CLUBB zt grid level (meters):  %g\n' ], ...
                 z_sam(upper_lev_idx) );
      fprintf( [ 'Altitude of the SAM LES grid level below the ', ...
                 'CLUBB zt grid level (meters):  %g\n' ], ...
                 z_sam(lower_lev_idx) );

      for i = 1:1:nx_sam
         for j = 1:1:ny_sam
            for idx_var = 1:1:num_var_sam

               sam_var_lev(idx_var,(i-1)*ny_sam+j) ...
               = interp_weight ...
                 * var_sam(idx_var,i,j,upper_lev_idx,1) ...
                 + ( 1.0 - interp_weight ) ...
                   * var_sam(idx_var,i,j,lower_lev_idx,1);

            end % idx_var = 1:1:num_var_sam
         end % j = 1:1:ny_sam
      end % i = 1:1:nx_sam

   end % exact_lev_idx > 0

end % z_clubb(clubb_height_idx)

%==========================================================================

% Unpack CLUBB variables (PDF parameters).

% PDF component means.
mu_w_1 = var_clubb( idx_w_1, 1, 1, clubb_height_idx, clubb_time_idx );
mu_w_2 = var_clubb( idx_w_2, 1, 1, clubb_height_idx, clubb_time_idx );
mu_rt_1 = var_clubb( idx_rt_1, 1, 1, clubb_height_idx, clubb_time_idx );
mu_rt_2 = var_clubb( idx_rt_2, 1, 1, clubb_height_idx, clubb_time_idx );
mu_thl_1 = var_clubb( idx_thl_1, 1, 1, clubb_height_idx, clubb_time_idx );
mu_thl_2 = var_clubb( idx_thl_2, 1, 1, clubb_height_idx, clubb_time_idx );
mu_chi_1 = var_clubb( idx_chi_1, 1, 1, clubb_height_idx, clubb_time_idx );
mu_chi_2 = var_clubb( idx_chi_2, 1, 1, clubb_height_idx, clubb_time_idx );
mu_eta_1 = 0.0; % The component mean of eta is always defined as 0.
mu_eta_2 = 0.0; % The component mean of eta is always defined as 0.
mu_rr_1_n = var_clubb( idx_mu_rr_1_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
mu_rr_2_n = var_clubb( idx_mu_rr_2_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
mu_Nr_1_n = var_clubb( idx_mu_Nr_1_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
mu_Nr_2_n = var_clubb( idx_mu_Nr_2_n, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );

% PDF component standard deviations.
sigma_w_1 = sqrt( var_clubb( idx_varnce_w_1, 1, 1, ...
                             clubb_height_idx, clubb_time_idx ) );
sigma_w_2 = sqrt( var_clubb( idx_varnce_w_2, 1, 1, ...
                             clubb_height_idx, clubb_time_idx ) );
sigma_rt_1 = sqrt( var_clubb( idx_varnce_rt_1, 1, 1, ...
                              clubb_height_idx, clubb_time_idx ) );
sigma_rt_2 = sqrt( var_clubb( idx_varnce_rt_2, 1, 1, ...
                              clubb_height_idx, clubb_time_idx ) );
sigma_thl_1 = sqrt( var_clubb( idx_varnce_thl_1, 1, 1, ...
                               clubb_height_idx, clubb_time_idx ) );
sigma_thl_2 = sqrt( var_clubb( idx_varnce_thl_2, 1, 1, ...
                               clubb_height_idx, clubb_time_idx ) );
sigma_chi_1 = var_clubb( idx_stdev_chi_1, 1, 1, ...
                         clubb_height_idx, clubb_time_idx );
sigma_chi_2 = var_clubb( idx_stdev_chi_2, 1, 1, ...
                         clubb_height_idx, clubb_time_idx );
sigma_eta_1 = var_clubb( idx_stdev_eta_1, 1, 1, ...
                         clubb_height_idx, clubb_time_idx );
sigma_eta_2 = var_clubb( idx_stdev_eta_2, 1, 1, ...
                         clubb_height_idx, clubb_time_idx );
sigma_rr_1_n = var_clubb( idx_sigma_rr_1_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );
sigma_rr_2_n = var_clubb( idx_sigma_rr_2_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );
sigma_Nr_1_n = var_clubb( idx_sigma_Nr_1_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );
sigma_Nr_2_n = var_clubb( idx_sigma_Nr_2_n, 1, 1, ...
                          clubb_height_idx, clubb_time_idx );

% PDF component correlations.
corr_rt_thl_1 = var_clubb( idx_corr_rt_thl, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_rt_thl_2 = corr_rt_thl_1; % CLUBB sets corr_rt_thl_1 = corr_rt_thl_2.
corr_chi_eta_1 = var_clubb( idx_corr_chi_eta_1_ca, 1, 1, ...
                            clubb_height_idx, clubb_time_idx );
corr_chi_eta_2 = var_clubb( idx_corr_chi_eta_2_ca, 1, 1, ...
                            clubb_height_idx, clubb_time_idx );
corr_w_rr_1_n = var_clubb( idx_corr_w_rr_1_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_w_rr_2_n = var_clubb( idx_corr_w_rr_2_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_w_Nr_1_n = var_clubb( idx_corr_w_Nr_1_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_w_Nr_2_n = var_clubb( idx_corr_w_Nr_2_n, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
corr_chi_rr_1_n = var_clubb( idx_corr_chi_rr_1_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_chi_rr_2_n = var_clubb( idx_corr_chi_rr_2_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_chi_Nr_1_n = var_clubb( idx_corr_chi_Nr_1_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_chi_Nr_2_n = var_clubb( idx_corr_chi_Nr_2_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_eta_rr_1_n = var_clubb( idx_corr_eta_rr_1_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_eta_rr_2_n = var_clubb( idx_corr_eta_rr_2_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_eta_Nr_1_n = var_clubb( idx_corr_eta_Nr_1_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_eta_Nr_2_n = var_clubb( idx_corr_eta_Nr_2_n, 1, 1, ...
                             clubb_height_idx, clubb_time_idx );
corr_rr_Nr_1_n = var_clubb( idx_corr_rr_Nr_1_n, 1, 1, ...
                            clubb_height_idx, clubb_time_idx );
corr_rr_Nr_2_n = var_clubb( idx_corr_rr_Nr_2_n, 1, 1, ...
                            clubb_height_idx, clubb_time_idx );

% Other variables involved in the PDF.
mixt_frac = var_clubb( idx_mixt_frac, 1, 1, ...
                       clubb_height_idx, clubb_time_idx );
precip_frac_1 = var_clubb( idx_precip_frac_1, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );
precip_frac_2 = var_clubb( idx_precip_frac_2, 1, 1, ...
                           clubb_height_idx, clubb_time_idx );

sigma_sqd_w = var_clubb( idx_sigma_sqd_w, 1, 1, ...
                         clubb_height_idx, clubb_time_idx );

%==========================================================================

% Calculations of distribution properties (skewness, etc.) from both SAM
% and CLUBB.

% Calculations based on SAM LES 3D data.
mean_w_sam = mean( sam_var_lev(idx_3D_w,:) );
mean_rt_sam = mean( sam_var_lev(idx_3D_rt,:) );
mean_thetal_sam = mean( sam_var_lev(idx_3D_thl,:) );
sumw2 = 0.0;
sumw3 = 0.0;
sumrt2 = 0.0;
sumrt3 = 0.0;
sumthl2 = 0.0;
sumthl3 = 0.0;
sumwrt = 0.0;
sumwthl = 0.0;
sumrtthl = 0.0;
for idx = 1:1:nx_sam*ny_sam
   sumw2   = sumw2   + ( sam_var_lev(idx_3D_w,idx)   - mean_w_sam )^2;
   sumw3   = sumw3   + ( sam_var_lev(idx_3D_w,idx)   - mean_w_sam )^3;
   sumrt2  = sumrt2  + ( sam_var_lev(idx_3D_rt,idx)  - mean_rt_sam )^2;
   sumrt3  = sumrt3  + ( sam_var_lev(idx_3D_rt,idx)  - mean_rt_sam )^3;
   sumthl2 = sumthl2 + ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam )^2;
   sumthl3 = sumthl3 + ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam )^3;
   sumwrt   = sumwrt ...
              + ( sam_var_lev(idx_3D_w,idx) - mean_w_sam ) ...
                * ( sam_var_lev(idx_3D_rt,idx) - mean_rt_sam );
   sumwthl  = sumwthl ...
              + ( sam_var_lev(idx_3D_w,idx) - mean_w_sam ) ...
                * ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam );
   sumrtthl = sumrtthl ...
              + ( sam_var_lev(idx_3D_rt,idx) - mean_rt_sam ) ...
                * ( sam_var_lev(idx_3D_thl,idx) - mean_thetal_sam );
end
wp2_sam   = sumw2   / ( nx_sam * ny_sam );
wp3_sam   = sumw3   / ( nx_sam * ny_sam );
rtp2_sam  = sumrt2  / ( nx_sam * ny_sam );
rtp3_sam  = sumrt3  / ( nx_sam * ny_sam );
thlp2_sam = sumthl2 / ( nx_sam * ny_sam );
thlp3_sam = sumthl3 / ( nx_sam * ny_sam );
wprtp_sam   = sumwrt   / ( nx_sam * ny_sam );
wpthlp_sam  = sumwthl  / ( nx_sam * ny_sam );
rtpthlp_sam = sumrtthl / ( nx_sam * ny_sam );

Skw_sam   = wp3_sam   / wp2_sam^1.5;
Skrt_sam  = rtp3_sam  / rtp2_sam^1.5;
Skthl_sam = thlp3_sam / thlp2_sam^1.5;
corr_w_rt_ov_sam ...
   = wprtp_sam   / ( sqrt( wp2_sam )  * sqrt( rtp2_sam ) );
corr_w_thl_ov_sam ...
   = wpthlp_sam  / ( sqrt( wp2_sam )  * sqrt( thlp2_sam ) );
corr_rt_thl_ov_sam ...
   = rtpthlp_sam / ( sqrt( rtp2_sam ) * sqrt( thlp2_sam ) );

% Calculations based on SAM LES 3D data from this level.
fprintf( '\n' )
fprintf( 'Calculations based on SAM LES 3D data at this level:\n' )
fprintf( 'SAM wp2 = %g\n', wp2_sam );
fprintf( 'SAM wp3 = %g\n', wp3_sam );
fprintf( 'SAM rtp2 = %g\n', rtp2_sam );
fprintf( 'SAM rtp3 = %g\n', rtp3_sam );
fprintf( 'SAM thlp2 = %g\n', thlp2_sam );
fprintf( 'SAM thlp3 = %g\n', thlp3_sam );
fprintf( 'SAM wprtp = %g\n', wprtp_sam );
fprintf( 'SAM wpthlp = %g\n', wpthlp_sam );
fprintf( 'SAM rtpthlp = %g\n', rtpthlp_sam );
fprintf( 'SAM skewness of w = %g\n', Skw_sam );
fprintf( 'SAM skewness of rt = %g\n', Skrt_sam );
fprintf( 'SAM skewness of thl = %g\n', Skthl_sam );
fprintf( 'SAM overall correlation of w and rt = %g\n', ...
         corr_w_rt_ov_sam );
fprintf( 'SAM overall correlation of w and thl = %g\n', ...
         corr_w_thl_ov_sam );
fprintf( 'SAM overall correlation of rt and thl = %g\n', ...
         corr_rt_thl_ov_sam );

% Calculations based on CLUBB PDF parameters.
wm   = mixt_frac * mu_w_1   + ( 1.0 - mixt_frac ) * mu_w_2;
rtm  = mixt_frac * mu_rt_1  + ( 1.0 - mixt_frac ) * mu_rt_2;
thlm = mixt_frac * mu_thl_1 + ( 1.0 - mixt_frac ) * mu_thl_2;
wp2_clubb_pdf ...
= mixt_frac * ( ( mu_w_1 - wm )^2 + sigma_w_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( ( mu_w_2 - wm )^2 + sigma_w_2^2 );
rtp2_clubb_pdf ...
= mixt_frac * ( ( mu_rt_1 - rtm )^2 + sigma_rt_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( ( mu_rt_2 - rtm )^2 + sigma_rt_2^2 );
thlp2_clubb_pdf ...
= mixt_frac * ( ( mu_thl_1 - thlm )^2 + sigma_thl_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( ( mu_thl_2 - thlm )^2 + sigma_thl_2^2 );
wp3_clubb_pdf ...
= mixt_frac * ( mu_w_1 - wm ) ...
            * ( ( mu_w_1 - wm )^2 + 3.0*sigma_w_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( mu_w_2 - wm ) ...
                        * ( ( mu_w_2 - wm )^2 + 3.0*sigma_w_2^2 );
rtp3_clubb_pdf ...
= mixt_frac * ( mu_rt_1 - rtm ) ...
            * ( ( mu_rt_1 - rtm )^2 + 3.0*sigma_rt_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( mu_rt_2 - rtm ) ...
                        * ( ( mu_rt_2 - rtm )^2 + 3.0*sigma_rt_2^2 );
thlp3_clubb_pdf ...
= mixt_frac * ( mu_thl_1 - thlm ) ...
            * ( ( mu_thl_1 - thlm )^2 + 3.0*sigma_thl_1^2 ) ...
  + ( 1.0 - mixt_frac ) * ( mu_thl_2 - thlm ) ...
                        * ( ( mu_thl_2 - thlm )^2 + 3.0*sigma_thl_2^2 );

Skw_clubb_pdf   = wp3_clubb_pdf   / wp2_clubb_pdf^1.5;
Skrt_clubb_pdf  = rtp3_clubb_pdf  / rtp2_clubb_pdf^1.5;
Skthl_clubb_pdf = thlp3_clubb_pdf / thlp2_clubb_pdf^1.5;

% Calculations based on CLUBB PDF parameters from this level.
fprintf( '\n' )
fprintf( 'Calculations based on CLUBB PDF parameters at this level:\n' )
fprintf( 'CLUBB PDF wp2 = %g\n', wp2_clubb_pdf );
fprintf( 'CLUBB PDF wp3 = %g\n', wp3_clubb_pdf );
fprintf( 'CLUBB PDF rtp2 = %g\n', rtp2_clubb_pdf );
fprintf( 'CLUBB PDF rtp3 = %g\n', rtp3_clubb_pdf );
fprintf( 'CLUBB PDF thlp2 = %g\n', thlp2_clubb_pdf );
fprintf( 'CLUBB PDF thlp3 = %g\n', thlp3_clubb_pdf );
fprintf( 'CLUBB PDF skewness of w = %g\n', Skw_clubb_pdf );
fprintf( 'CLUBB PDF skewness of rt = %g\n', Skrt_clubb_pdf );
fprintf( 'CLUBB PDF skewness of thl = %g\n', Skthl_clubb_pdf );

% Backsolve for the appropriate values of beta (for each of rt and thl)
% using SAM LES data (and CLUBB output for sigma_sqd_w) based on Larson
% and Golaz (2005), Eq. 33.
%Skw_hat = Skw_sam / ( 1.0 - sigma_sqd_w )^1.5;
%corr_w_rt_ov_hat  = corr_w_rt_ov_sam  / sqrt( 1.0 - sigma_sqd_w );
%corr_w_thl_ov_hat = corr_w_thl_ov_sam / sqrt( 1.0 - sigma_sqd_w );
%
%beta_rt ...
%= ( Skrt_sam / ( Skw_hat * corr_w_rt_ov_hat ) - corr_w_rt_ov_hat^2 ) ...
%  / ( 1.0 - corr_w_rt_ov_hat^2 );
%
%beta_thl ...
%= ( Skthl_sam / ( Skw_hat * corr_w_thl_ov_hat ) - corr_w_thl_ov_hat^2 ) ...
%  / ( 1.0 - corr_w_thl_ov_hat^2 );
%
% Calculations based on CLUBB PDF parameters from this level.
%fprintf( '\n' )
%fprintf( 'Backsolve for the parameter beta based largely on SAM LES:\n' )
%fprintf( 'Ideal beta for skewness of rt = %g\n', beta_rt );
%fprintf( 'Ideal beta for skewness of thl = %g\n', beta_thl );

% Calculations for hydrometeor fields based on SAM LES 3D data.
mean_rr_sam = mean( sam_var_lev(idx_3D_rr,:) );
mean_Nr_sam = mean( sam_var_lev(idx_3D_Nr,:) );
sumrr2 = 0.0;
sumrr3 = 0.0;
sumNr2 = 0.0;
sumNr3 = 0.0;
for idx = 1:1:nx_sam*ny_sam
   sumrr2 = sumrr2 + ( sam_var_lev(idx_3D_rr,idx) - mean_rr_sam )^2;
   sumrr3 = sumrr3 + ( sam_var_lev(idx_3D_rr,idx) - mean_rr_sam )^3;
   sumNr2 = sumNr2 + ( sam_var_lev(idx_3D_Nr,idx) - mean_Nr_sam )^2;
   sumNr3 = sumNr3 + ( sam_var_lev(idx_3D_Nr,idx) - mean_Nr_sam )^3;
end
rrp2_sam = sumrr2 / ( nx_sam * ny_sam );
rrp3_sam = sumrr3 / ( nx_sam * ny_sam );
Nrp2_sam = sumNr2 / ( nx_sam * ny_sam );
Nrp3_sam = sumNr3 / ( nx_sam * ny_sam );

Skrr_sam = rrp3_sam / rrp2_sam^1.5;
SkNr_sam = Nrp3_sam / Nrp2_sam^1.5;

precipitating = zeros(nx_sam*ny_sam,1);
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_rr,idx) > 0.0 ...
        || sam_var_lev(idx_3D_Nr,idx) > 0.0 )
      precipitating(idx) = 1.0;
   else
      precipitating(idx) = 0.0;
   end
end
precip_frac_sam = sum( precipitating ) / ( nx_sam * ny_sam );

if ( mean_rr_sam > 0.0 )
   % Rain water mixing ratio is found at this level in SAM.
   [ mean_rr_ip_sam, rrp2_ip_sam, rrp3_ip_sam, Skrr_ip_sam ] ...
   = calc_in_precip_values( mean_rr_sam, rrp2_sam, ...
                            rrp3_sam, precip_frac_sam );
   % Calculate in-precip variance-over-mean-squared for rr.
   voms_ip_rr_sam = rrp2_ip_sam / mean_rr_ip_sam^2;
else
   % Set values to 0.
   mean_rr_ip_sam = 0.0;
   rrp2_ip_sam = 0.0;
   rrp2_ip_sam = 0.0;
   Skrr_ip_sam = 0.0;
   voms_ip_rr_sam = 0.0;
end

if ( mean_Nr_sam > 0.0 )
   % Rain drop concentration is found at this level in SAM.
   [ mean_Nr_ip_sam, Nrp2_ip_sam, Nrp3_ip_sam, SkNr_ip_sam ] ...
   = calc_in_precip_values( mean_Nr_sam, Nrp2_sam, ...
                            Nrp3_sam, precip_frac_sam );
   % Calculate in-precip variance-over-mean-squared for Nr.
   voms_ip_Nr_sam = Nrp2_ip_sam / mean_Nr_ip_sam^2;
else
   % Set values to 0.
   mean_Nr_ip_sam = 0.0;
   Nrp2_ip_sam = 0.0;
   Nrp2_ip_sam = 0.0;
   SkNr_ip_sam = 0.0;
   voms_ip_Nr_sam = 0.0;
end

% Calculate the mean, variance, 3rd-order central moment, and skewness of
% ln rr for all rr in precip. (where rr > 0) for SAM LES.
count_lnrr = 0;
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_rr,idx) > 0.0 )
      count_lnrr = count_lnrr + 1;
      sam_lnrr(count_lnrr) = log( sam_var_lev(idx_3D_rr,idx) );
   end
end
if ( count_lnrr > 0 )
   mean_lnrr_sam = mean( sam_lnrr );
   sum_lnrr2 = 0.0;
   sum_lnrr3 = 0.0;
   for idx = 1:1:count_lnrr
      sum_lnrr2 = sum_lnrr2 + ( sam_lnrr(idx) - mean_lnrr_sam )^2;
      sum_lnrr3 = sum_lnrr3 + ( sam_lnrr(idx) - mean_lnrr_sam )^3;
   end
   lnrrp2_sam = sum_lnrr2 / count_lnrr;
   lnrrp3_sam = sum_lnrr3 / count_lnrr;
   Sk_lnrr_sam = lnrrp3_sam / lnrrp2_sam^1.5;
end % count_lnrr > 0
count_lnNr = 0;
for idx = 1:1:nx_sam*ny_sam
   if ( sam_var_lev(idx_3D_Nr,idx) > 0.0 )
      count_lnNr = count_lnNr + 1;
      sam_lnNr(count_lnNr) = log( sam_var_lev(idx_3D_Nr,idx) );
   end
end
% Calculate the mean, variance, 3rd-order central moment, and skewness of
% ln Nr for all Nr in precip. (where Nr > 0) for SAM LES.
if ( count_lnNr > 0 )
   mean_lnNr_sam = mean( sam_lnNr );
   sum_lnNr2 = 0.0;
   sum_lnNr3 = 0.0;
   for idx = 1:1:count_lnNr
      sum_lnNr2 = sum_lnNr2 + ( sam_lnNr(idx) - mean_lnNr_sam )^2;
      sum_lnNr3 = sum_lnNr3 + ( sam_lnNr(idx) - mean_lnNr_sam )^3;
   end
   lnNrp2_sam = sum_lnNr2 / count_lnNr;
   lnNrp3_sam = sum_lnNr3 / count_lnNr;
   Sk_lnNr_sam = lnNrp3_sam / lnNrp2_sam^1.5;
end % count_lnNr > 0

% Calculations for hydrometeor fields based on CLUBB PDF parameters.
rrm_clubb_pdf = mixt_frac * precip_frac_1 ...
                * exp( mu_rr_1_n + 0.5 * sigma_rr_1_n^2 ) ...
                + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                  * exp( mu_rr_2_n + 0.5 * sigma_rr_2_n^2 );
Nrm_clubb_pdf = mixt_frac * precip_frac_1 ...
                * exp( mu_Nr_1_n + 0.5 * sigma_Nr_1_n^2 ) ...
                + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                  * exp( mu_Nr_2_n + 0.5 * sigma_Nr_2_n^2 );
rrp2_clubb_pdf = mixt_frac * precip_frac_1 ...
                 * exp( 2.0 * mu_rr_1_n + 2.0 * sigma_rr_1_n^2 ) ...
                 + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                   * exp( 2.0 * mu_rr_2_n + 2.0 * sigma_rr_2_n^2 ) ...
                 - rrm_clubb_pdf^2;
Nrp2_clubb_pdf = mixt_frac * precip_frac_1 ...
                 * exp( 2.0 * mu_Nr_1_n + 2.0 * sigma_Nr_1_n^2 ) ...
                 + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                   * exp( 2.0 * mu_Nr_2_n + 2.0 * sigma_Nr_2_n^2 ) ...
                 - Nrm_clubb_pdf^2;
rrp3_clubb_pdf = mixt_frac * precip_frac_1 ...
                 * exp( 3.0 * mu_rr_1_n + 4.5 * sigma_rr_1_n^2 ) ...
                 + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                   * exp( 3.0 * mu_rr_2_n + 4.5 * sigma_rr_2_n^2 ) ...
                 - 3.0 * rrp2_clubb_pdf * rrm_clubb_pdf ...
                 - rrm_clubb_pdf^3;
Nrp3_clubb_pdf = mixt_frac * precip_frac_1 ...
                 * exp( 3.0 * mu_Nr_1_n + 4.5 * sigma_Nr_1_n^2 ) ...
                 + ( 1.0 - mixt_frac ) * precip_frac_2 ...
                   * exp( 3.0 * mu_Nr_2_n + 4.5 * sigma_Nr_2_n^2 ) ...
                 - 3.0 * Nrp2_clubb_pdf * Nrm_clubb_pdf ...
                 - Nrm_clubb_pdf^3;

Skrr_clubb_pdf = rrp3_clubb_pdf / rrp2_clubb_pdf^1.5;
SkNr_clubb_pdf = Nrp3_clubb_pdf / Nrp2_clubb_pdf^1.5;

precip_frac_clubb = mixt_frac * precip_frac_1 ...
                    + ( 1.0 - mixt_frac ) * precip_frac_2;

if ( rrm_clubb_pdf > 0.0 )
   % Rain water mixing ratio is found at this level in CLUBB.
   [ rrm_ip_clubb_pdf, rrp2_ip_clubb_pdf, ...
     rrp3_ip_clubb_pdf, Skrr_ip_clubb_pdf ] ...
   = calc_in_precip_values( rrm_clubb_pdf, rrp2_clubb_pdf, ...
                            rrp3_clubb_pdf, precip_frac_clubb );
   % Calculate in-precip variance-over-mean-squared for rr.
   voms_ip_rr_clubb_pdf = rrp2_ip_clubb_pdf / rrm_ip_clubb_pdf^2;
else
   % Set values to 0.
   rrm_ip_clubb_pdf = 0.0;
   rrp2_ip_clubb_pdf = 0.0;
   rrp2_ip_clubb_pdf = 0.0;
   Skrr_ip_clubb_pdf = 0.0;
   voms_ip_rr_clubb_pdf = 0.0;
end

if ( Nrm_clubb_pdf > 0.0 )
   % Rain drop concentration is found at this level in CLUBB.
   [ Nrm_ip_clubb_pdf, Nrp2_ip_clubb_pdf, ...
     Nrp3_ip_clubb_pdf, SkNr_ip_clubb_pdf ] ...
   = calc_in_precip_values( Nrm_clubb_pdf, Nrp2_clubb_pdf, ...
                            Nrp3_clubb_pdf, precip_frac_clubb );
   % Calculate in-precip variance-over-mean-squared for Nr.
   voms_ip_Nr_clubb_pdf = Nrp2_ip_clubb_pdf / Nrm_ip_clubb_pdf^2;
else
   % Set values to 0.
   Nrm_ip_clubb_pdf = 0.0;
   Nrp2_ip_clubb_pdf = 0.0;
   Nrp2_ip_clubb_pdf = 0.0;
   SkNr_ip_clubb_pdf = 0.0;
   voms_ip_Nr_clubb_pdf = 0.0;
end

% Calculate the mean, variance, 3rd-order central moment, and skewness of
% ln rr for all rr in precip. (where rr > 0) for CLUBB.
if ( mu_rr_1_n > -realmax('single') && mu_rr_2_n > -realmax('single') )

   mean_lnrr_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) * mu_rr_1_n ...
     + ( 1.0 - mixt_frac ) ...
       * ( precip_frac_2 / precip_frac_clubb ) * mu_rr_2_n;

   lnrrp2_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_rr_1_n - mean_lnrr_clubb_pdf )^2 + sigma_rr_1_n^2 ) ...
     + ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
       * ( ( mu_rr_2_n - mean_lnrr_clubb_pdf )^2 + sigma_rr_2_n^2 );

   lnrrp3_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_rr_1_n - mean_lnrr_clubb_pdf )^3 ...
         + 3.0 * ( mu_rr_1_n - mean_lnrr_clubb_pdf ) * sigma_rr_1_n^2 ) ...
     + ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
       * ( ( mu_rr_2_n - mean_lnrr_clubb_pdf )^3 ...
           + 3.0 * ( mu_rr_2_n - mean_lnrr_clubb_pdf ) * sigma_rr_2_n^2 );

   Sk_lnrr_clubb_pdf = lnrrp3_clubb_pdf / lnrrp2_clubb_pdf^1.5;

elseif ( mu_rr_1_n > -realmax('single') )

   mean_lnrr_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) * mu_rr_1_n;

   lnrrp2_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_rr_1_n - mean_lnrr_clubb_pdf )^2 + sigma_rr_1_n^2 );

   lnrrp3_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_rr_1_n - mean_lnrr_clubb_pdf )^3 ...
         + 3.0 * ( mu_rr_1_n - mean_lnrr_clubb_pdf ) * sigma_rr_1_n^2 );

   Sk_lnrr_clubb_pdf = lnrrp3_clubb_pdf / lnrrp2_clubb_pdf^1.5;

elseif ( mu_rr_2_n > -realmax('single') )

   mean_lnrr_clubb_pdf ...
   = ( 1.0 - mixt_frac ) ...
     * ( precip_frac_2 / precip_frac_clubb ) * mu_rr_2_n;

   lnrrp2_clubb_pdf ...
   = ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
     * ( ( mu_rr_2_n - mean_lnrr_clubb_pdf )^2 + sigma_rr_2_n^2 );

   lnrrp3_clubb_pdf ...
   = ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
     * ( ( mu_rr_2_n - mean_lnrr_clubb_pdf )^3 ...
         + 3.0 * ( mu_rr_2_n - mean_lnrr_clubb_pdf ) * sigma_rr_2_n^2 );

   Sk_lnrr_clubb_pdf = lnrrp3_clubb_pdf / lnrrp2_clubb_pdf^1.5;

else

   mean_lnrr_clubb_pdf = -realmax('single');
   lnrrp2_clubb_pdf = 0.0;
   lnrrp3_clubb_pdf = 0.0;
   Sk_lnrr_clubb_pdf = 0.0;

end

if ( mu_Nr_1_n > -realmax('single') && mu_Nr_2_n > -realmax('single') )

   mean_lnNr_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) * mu_Nr_1_n ...
     + ( 1.0 - mixt_frac ) ...
       * ( precip_frac_2 / precip_frac_clubb ) * mu_Nr_2_n;

   lnNrp2_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_Nr_1_n - mean_lnNr_clubb_pdf )^2 + sigma_Nr_1_n^2 ) ...
     + ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
       * ( ( mu_Nr_2_n - mean_lnNr_clubb_pdf )^2 + sigma_Nr_2_n^2 );

   lnNrp3_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_Nr_1_n - mean_lnNr_clubb_pdf )^3 ...
         + 3.0 * ( mu_Nr_1_n - mean_lnNr_clubb_pdf ) * sigma_Nr_1_n^2 ) ...
     + ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
       * ( ( mu_Nr_2_n - mean_lnNr_clubb_pdf )^3 ...
           + 3.0 * ( mu_Nr_2_n - mean_lnNr_clubb_pdf ) * sigma_Nr_2_n^2 );

   Sk_lnNr_clubb_pdf = lnNrp3_clubb_pdf / lnNrp2_clubb_pdf^1.5;

elseif ( mu_Nr_1_n > -realmax('single') )

   mean_lnNr_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) * mu_Nr_1_n;

   lnNrp2_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_Nr_1_n - mean_lnNr_clubb_pdf )^2 + sigma_Nr_1_n^2 );

   lnNrp3_clubb_pdf ...
   = mixt_frac * ( precip_frac_1 / precip_frac_clubb ) ...
     * ( ( mu_Nr_1_n - mean_lnNr_clubb_pdf )^3 ...
         + 3.0 * ( mu_Nr_1_n - mean_lnNr_clubb_pdf ) * sigma_Nr_1_n^2 );

   Sk_lnNr_clubb_pdf = lnNrp3_clubb_pdf / lnNrp2_clubb_pdf^1.5;

elseif ( mu_Nr_2_n > -realmax('single') )

   mean_lnNr_clubb_pdf ...
   = ( 1.0 - mixt_frac ) ...
     * ( precip_frac_2 / precip_frac_clubb ) * mu_Nr_2_n;

   lnNrp2_clubb_pdf ...
   = ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
     * ( ( mu_Nr_2_n - mean_lnNr_clubb_pdf )^2 + sigma_Nr_2_n^2 );

   lnNrp3_clubb_pdf ...
   = ( 1.0 - mixt_frac ) * ( precip_frac_2 / precip_frac_clubb ) ...
     * ( ( mu_Nr_2_n - mean_lnNr_clubb_pdf )^3 ...
         + 3.0 * ( mu_Nr_2_n - mean_lnNr_clubb_pdf ) * sigma_Nr_2_n^2 );

   Sk_lnNr_clubb_pdf = lnNrp3_clubb_pdf / lnNrp2_clubb_pdf^1.5;

else

   mean_lnNr_clubb_pdf = -realmax('single');
   lnNrp2_clubb_pdf = 0.0;
   lnNrp3_clubb_pdf = 0.0;
   Sk_lnNr_clubb_pdf = 0.0;

end


% Print the values of hydrometeor statistics.
% Note:  v.o.m.s. stands for variance over mean-squared.
fprintf( '\n' )
fprintf( 'Calculations for hydrometeors at this level:\n' )
fprintf( '\n' )
fprintf( 'SAM precipitation fraction = %g\n', precip_frac_sam );
fprintf( 'CLUBB precipitation fraction = %g\n', precip_frac_clubb );
fprintf( '\n' )
fprintf( 'Rain water mixing ratio, rr\n' )
fprintf( 'Overall values:\n' )
fprintf( 'SAM rrm = %g\n', mean_rr_sam );
fprintf( 'SAM rrp2 = %g\n', rrp2_sam );
fprintf( 'SAM rrp3 = %g\n', rrp3_sam );
fprintf( 'SAM skewness of rr = %g\n', Skrr_sam );
fprintf( 'CLUBB PDF rrm = %g\n', rrm_clubb_pdf );
fprintf( 'CLUBB PDF rrp2 = %g\n', rrp2_clubb_pdf );
fprintf( 'CLUBB PDF rrp3 = %g\n', rrp3_clubb_pdf );
fprintf( 'CLUBB PDF skewness of rr = %g\n', Skrr_clubb_pdf );
fprintf( 'In-precip. (i.p.) values:\n' )
fprintf( 'SAM mean (i.p.) of rr = %g\n', mean_rr_ip_sam );
fprintf( 'SAM variance (i.p.) of rr = %g\n', rrp2_ip_sam );
fprintf( 'SAM 3rd moment (i.p.) of rr = %g\n', rrp3_ip_sam );
fprintf( 'SAM skewness (i.p.) of rr = %g\n', Skrr_ip_sam );
fprintf( 'SAM v.o.m.s (i.p.) of rr = %g\n', voms_ip_rr_sam );
fprintf( 'SAM mean (i.p.) of ln rr = %g\n', mean_lnrr_sam );
fprintf( 'SAM variance (i.p.) of ln rr = %g\n', lnrrp2_sam );
fprintf( 'SAM 3rd moment (i.p.) of ln rr = %g\n', lnrrp3_sam );
fprintf( 'SAM skewness (i.p.) of ln rr = %g\n', Sk_lnrr_sam );
fprintf( 'CLUBB PDF mean (i.p.) of rr = %g\n', rrm_ip_clubb_pdf );
fprintf( 'CLUBB PDF variance (i.p.) of rr = %g\n', rrp2_ip_clubb_pdf );
fprintf( 'CLUBB PDF 3rd moment (i.p.) of rr = %g\n', rrp3_ip_clubb_pdf );
fprintf( 'CLUBB PDF skewness (i.p.) of rr = %g\n', Skrr_ip_clubb_pdf );
fprintf( 'CLUBB PDF v.o.m.s (i.p.) of rr = %g\n', voms_ip_rr_clubb_pdf );
fprintf( 'CLUBB PDF mean (i.p.) of ln rr = %g\n', mean_lnrr_clubb_pdf );
fprintf( 'CLUBB PDF variance (i.p.) of ln rr = %g\n', lnrrp2_clubb_pdf );
fprintf( 'CLUBB PDF 3rd moment (i.p.) of ln rr = %g\n', lnrrp3_clubb_pdf );
fprintf( 'CLUBB PDF skewness (i.p.) of ln rr = %g\n', Sk_lnrr_clubb_pdf );
fprintf( '\n' )
fprintf( 'Rain drop concentration, Nr\n' )
fprintf( 'Overall values:\n' )
fprintf( 'SAM Nrm = %g\n', mean_Nr_sam );
fprintf( 'SAM Nrp2 = %g\n', Nrp2_sam );
fprintf( 'SAM Nrp3 = %g\n', Nrp3_sam );
fprintf( 'SAM skewness of Nr = %g\n', SkNr_sam );
fprintf( 'CLUBB PDF Nrm = %g\n', Nrm_clubb_pdf );
fprintf( 'CLUBB PDF Nrp2 = %g\n', Nrp2_clubb_pdf );
fprintf( 'CLUBB PDF Nrp3 = %g\n', Nrp3_clubb_pdf );
fprintf( 'CLUBB PDF skewness of Nr = %g\n', SkNr_clubb_pdf );
fprintf( 'In-precip. (i.p.) values:\n' )
fprintf( 'SAM mean (i.p.) of Nr = %g\n', mean_Nr_ip_sam );
fprintf( 'SAM variance (i.p.) of Nr = %g\n', Nrp2_ip_sam );
fprintf( 'SAM 3rd moment (i.p.) of Nr = %g\n', Nrp3_ip_sam );
fprintf( 'SAM skewness (i.p.) of Nr = %g\n', SkNr_ip_sam );
fprintf( 'SAM v.o.m.s (i.p.) of Nr = %g\n', voms_ip_Nr_sam );
fprintf( 'SAM mean (i.p.) of ln Nr = %g\n', mean_lnNr_sam );
fprintf( 'SAM variance (i.p.) of ln Nr = %g\n', lnNrp2_sam );
fprintf( 'SAM 3rd moment (i.p.) of ln Nr = %g\n', lnNrp3_sam );
fprintf( 'SAM skewness (i.p.) of ln Nr = %g\n', Sk_lnNr_sam );
fprintf( 'CLUBB PDF mean (i.p.) of Nr = %g\n', Nrm_ip_clubb_pdf );
fprintf( 'CLUBB PDF variance (i.p.) of Nr = %g\n', Nrp2_ip_clubb_pdf );
fprintf( 'CLUBB PDF 3rd moment (i.p.) of Nr = %g\n', Nrp3_ip_clubb_pdf );
fprintf( 'CLUBB PDF skewness (i.p.) of Nr = %g\n', SkNr_ip_clubb_pdf );
fprintf( 'CLUBB PDF v.o.m.s (i.p.) of Nr = %g\n', voms_ip_Nr_clubb_pdf );
fprintf( 'CLUBB PDF mean (i.p.) of ln Nr = %g\n', mean_lnNr_clubb_pdf );
fprintf( 'CLUBB PDF variance (i.p.) of ln Nr = %g\n', lnNrp2_clubb_pdf );
fprintf( 'CLUBB PDF 3rd moment (i.p.) of ln Nr = %g\n', lnNrp3_clubb_pdf );
fprintf( 'CLUBB PDF skewness (i.p.) of ln Nr = %g\n', Sk_lnNr_clubb_pdf );

% Print a blank line after the calculations have been printed.
fprintf( '\n' )

%==========================================================================

% Information to be printed on the plots.
print_alt = int2str( round( z_clubb(clubb_height_idx) ) );
print_time = int2str( round( time_clubb(clubb_time_idx) / 60.0 ) );

% When the SAM LES data for rr is 0 everywhere at the level, the plots
% involving rr will fail, causing an exit with an error.  Since these
% plots aren't interesting anyway when there's no rr, simply turn them
% of to avoid the error.
if ( all( sam_var_lev(idx_3D_rr,:) == 0.0 ) )

   fprintf( [ 'The SAM LES values of rr are 0 everywhere at this ', ...
              'level.  Any plots involving rr will be disabled.\n' ] )

   plot_w_rr   = false;
   plot_chi_rr = false;
   plot_eta_rr = false;
   plot_rr_Nr  = false;
   plot_ln_rr  = false;

end

% When the SAM LES data for Nr is 0 everywhere at the level, the plots
% involving Nr will fail, causing an exit with an error.  Since these
% plots aren't interesting anyway when there's no Nr, simply turn them
% of to avoid the error.
if ( all( sam_var_lev(idx_3D_Nr,:) == 0.0 ) )

   fprintf( [ 'The SAM LES values of Nr are 0 everywhere at this ', ...
              'level.  Any plots involving Nr will be disabled.\n' ] )

   plot_w_Nr   = false;
   plot_chi_Nr = false;
   plot_eta_Nr = false;
   plot_rr_Nr  = false;
   plot_ln_Nr  = false;

end
 
% Plot the CLUBB PDF and LES points for rt and thl.
if ( plot_rt_thl )

   fprintf( 'Plotting scatter/contour plot for rt and thl\n' );

   % The number of rt bins/points to plot.
   num_rt_pts = 100;

   % The number of thl bins/points to plot.
   num_thl_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 100;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NN( sam_var_lev(idx_3D_rt,:), ...
                              sam_var_lev(idx_3D_thl,:), ...
                              nx_sam, ny_sam, ...
                              num_rt_pts, num_thl_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, ...
                              sigma_rt_1, sigma_rt_2, sigma_thl_1, ...
                              sigma_thl_2, corr_rt_thl_1, ...
                              corr_rt_thl_2, mixt_frac, ...
                              'r_{t}    [kg/kg]', '\theta_{l}    [K]', ...
                              [ '\bf ', casename ], 'r_{t} vs. \theta_{l}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_rt_thl_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_rt_thl

% Plot the CLUBB PDF and LES points for chi and eta.
if ( plot_chi_eta )

   fprintf( 'Plotting scatter/contour plot for chi and eta\n' );

   % The number of chi bins/points to plot.
   num_chi_pts = 100;

   % The number of eta bins/points to plot.
   num_eta_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 100;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NN( sam_var_lev(idx_3D_chi,:), ...
                              sam_var_lev(idx_3D_eta,:), ...
                              nx_sam, ny_sam, ...
                              num_chi_pts, num_eta_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_chi_1, mu_chi_2, mu_eta_1, mu_eta_2, ...
                              sigma_chi_1, sigma_chi_2, sigma_eta_1, ...
                              sigma_eta_2, corr_chi_eta_1, ...
                              corr_chi_eta_2, mixt_frac, ...
                              '\chi    [kg/kg]', '\eta    [kg/kg]', ...
                              [ '\bf ', casename ], '\chi vs. \eta', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_eta_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_chi_eta

% Plot the CLUBB PDF and LES points for w and rr.
if ( plot_w_rr )

   fprintf( 'Plotting scatter/contour plot for w and rr\n' );

   % The number of w bins/points to plot.
   num_w_pts = 100;

   % The number of rr bins/points to plot.
   num_rr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_w,:), ...
                              sam_var_lev(idx_3D_rr,:), ...
                              nx_sam, ny_sam, ...
                              num_w_pts, num_rr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_w_1, mu_w_2, mu_rr_1_n, mu_rr_2_n, ...
                              sigma_w_1, sigma_w_2, sigma_rr_1_n, ...
                              sigma_rr_2_n, corr_w_rr_1_n, ...
                              corr_w_rr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              'w    [m/s]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], 'w vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_w_rr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_w_rr

% Plot the CLUBB PDF and LES points for w and Nr.
if ( plot_w_Nr )

   fprintf( 'Plotting scatter/contour plot for w and Nr\n' );

   % The number of w bins/points to plot.
   num_w_pts = 100;

   % The number of Nr bins/points to plot.
   num_Nr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_w,:), ...
                              sam_var_lev(idx_3D_Nr,:), ...
                              nx_sam, ny_sam, ...
                              num_w_pts, num_Nr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_w_1, mu_w_2, mu_Nr_1_n, mu_Nr_2_n, ...
                              sigma_w_1, sigma_w_2, sigma_Nr_1_n, ...
                              sigma_Nr_2_n, corr_w_Nr_1_n, ...
                              corr_w_Nr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              'w    [m/s]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], 'w vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_w_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_w_Nr

% Plot the CLUBB PDF and LES points for chi and rr.
if ( plot_chi_rr )

   fprintf( 'Plotting scatter/contour plot for chi and rr\n' );

   % The number of chi bins/points to plot.
   num_chi_pts = 100;

   % The number of rr bins/points to plot.
   num_rr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_chi,:), ...
                              sam_var_lev(idx_3D_rr,:), ...
                              nx_sam, ny_sam, ...
                              num_chi_pts, num_rr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_chi_1, mu_chi_2, mu_rr_1_n, mu_rr_2_n, ...
                              sigma_chi_1, sigma_chi_2, sigma_rr_1_n, ...
                              sigma_rr_2_n, corr_chi_rr_1_n, ...
                              corr_chi_rr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              '\chi    [kg/kg]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], '\chi vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_rr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_chi_rr

% Plot the CLUBB PDF and LES points for chi and Nr.
if ( plot_chi_Nr )

   fprintf( 'Plotting scatter/contour plot for chi and Nr\n' );

   % The number of chi bins/points to plot.
   num_chi_pts = 100;

   % The number of Nr bins/points to plot.
   num_Nr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_chi,:), ...
                              sam_var_lev(idx_3D_Nr,:), ...
                              nx_sam, ny_sam, ...
                              num_chi_pts, num_Nr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_chi_1, mu_chi_2, mu_Nr_1_n, mu_Nr_2_n, ...
                              sigma_chi_1, sigma_chi_2, sigma_Nr_1_n, ...
                              sigma_Nr_2_n, corr_chi_Nr_1_n, ...
                              corr_chi_Nr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              '\chi    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], '\chi vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_chi_Nr

% Plot the CLUBB PDF and LES points for eta and rr.
if ( plot_eta_rr)

   fprintf( 'Plotting scatter/contour plot for eta and rr\n' );

   % The number of eta bins/points to plot.
   num_eta_pts = 100;

   % The number of rr bins/points to plot.
   num_rr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_eta,:), ...
                              sam_var_lev(idx_3D_rr,:), ...
                              nx_sam, ny_sam, ...
                              num_eta_pts, num_rr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_eta_1, mu_eta_2, mu_rr_1_n, mu_rr_2_n, ...
                              sigma_eta_1, sigma_eta_2, sigma_rr_1_n, ...
                              sigma_rr_2_n, corr_eta_rr_1_n, ...
                              corr_eta_rr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              '\eta    [kg/kg]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], '\eta vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_eta_rr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_eta_rr

% Plot the CLUBB PDF and LES points for eta and Nr.
if ( plot_eta_Nr )

   fprintf( 'Plotting scatter/contour plot for eta and Nr\n' );

   % The number of eta bins/points to plot.
   num_eta_pts = 100;

   % The number of Nr bins/points to plot.
   num_Nr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 50;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_NL( sam_var_lev(idx_3D_eta,:), ...
                              sam_var_lev(idx_3D_Nr,:), ...
                              nx_sam, ny_sam, ...
                              num_eta_pts, num_Nr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_eta_1, mu_eta_2, mu_Nr_1_n, mu_Nr_2_n, ...
                              sigma_eta_1, sigma_eta_2, sigma_Nr_1_n, ...
                              sigma_Nr_2_n, corr_eta_Nr_1_n, ...
                              corr_eta_Nr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              '\eta    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], '\eta vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_eta_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_eta_Nr

% Plot the CLUBB PDF and LES points for rr and Nr.
if ( plot_rr_Nr )

   fprintf( 'Plotting scatter/contour plot for rr and Nr\n' );

   % The number of rr bins/points to plot.
   num_rr_pts = 100;

   % The number of Nr bins/points to plot.
   num_Nr_pts = 100;

   % The number of contours on the scatter/contour plot.
   num_contours = 100;

   % The number of standard deviations used to help find a minimum contour
   % for the scatter/contour plot.
   num_std_devs_min_contour = 2.0;

   plot_CLUBB_PDF_LES_pts_LL( sam_var_lev(idx_3D_rr,:), ...
                              sam_var_lev(idx_3D_Nr,:), ...
                              nx_sam, ny_sam, ...
                              num_rr_pts, num_Nr_pts, num_contours, ...
                              num_std_devs_min_contour, ...
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, ...
                              sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, ...
                              sigma_Nr_2_n, corr_rr_Nr_1_n, ...
                              corr_rr_Nr_2_n, precip_frac_1, ...
                              precip_frac_2, mixt_frac, ...
                              'r_{r}    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], 'r_{r} vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_rr_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_rr_Nr

% Plot the CLUBB PDF and LES points for w.
if ( plot_w )

   fprintf( 'Plotting marginal plot for w\n' );

   % The number of ln rr bins/points to plot.
   num_w_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_w,:), ...
                             nx_sam, ny_sam, ...
                             num_w_pts, ...
                             mu_w_1, mu_w_2, ...
                             sigma_w_1, sigma_w_2, ...
                             mixt_frac, 'w    [m/s]', ...
                             [ '\bf ', casename ], 'w', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_w_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_w

% Plot the CLUBB PDF and LES points for ln rr.
if ( plot_ln_rr )

   fprintf( 'Plotting marginal plot for ln rr\n' );

   % The number of ln rr bins/points to plot.
   num_rr_pts = 100;

   % Absolute minimum value in ln rr plot.
   ln_rr_abs_min = log( 1.0e-8 );

   plot_CLUBB_PDF_LES_pts_ln_L( sam_var_lev(idx_3D_rr,:), ...
                                nx_sam, ny_sam, ...
                                num_rr_pts, ln_rr_abs_min, ...
                                mu_rr_1_n, mu_rr_2_n, ...
                                sigma_rr_1_n, sigma_rr_2_n, ...
                                precip_frac_1, precip_frac_2, ...
                                mixt_frac, 'ln r_{r}    [ln(kg/kg)]', ...
                                [ '\bf ', casename ], 'ln r_{r}', ...
                                [ 'Time = ', print_time, ' minutes' ], ...
                                [ 'Altitude = ', print_alt, ' meters' ], ...
                                print_note )

   output_filename = [ 'output/', casename, '_ln_rr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_ln_rr

% Plot the CLUBB PDF and LES points for ln Nr.
if ( plot_ln_Nr )

   fprintf( 'Plotting marginal plot for ln Nr\n' );

   % The number of ln Nr bins/points to plot.
   num_Nr_pts = 100;

   % Absolute minimum value in ln Nr plot.
   ln_Nr_abs_min = log( 1.0 );

   plot_CLUBB_PDF_LES_pts_ln_L( sam_var_lev(idx_3D_Nr,:), ...
                                nx_sam, ny_sam, ...
                                num_Nr_pts, ln_Nr_abs_min, ...
                                mu_Nr_1_n, mu_Nr_2_n, ...
                                sigma_Nr_1_n, sigma_Nr_2_n, ...
                                precip_frac_1, precip_frac_2, ...
                                mixt_frac, 'ln N_{r}    [ln(num/kg)]', ...
                                [ '\bf ', casename ], 'ln N_{r}', ...
                                [ 'Time = ', print_time, ' minutes' ], ...
                                [ 'Altitude = ', print_alt, ' meters' ], ...
                                print_note )

   output_filename = [ 'output/', casename, '_ln_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( '-dpng', output_filename );

end % plot_ln_Nr
