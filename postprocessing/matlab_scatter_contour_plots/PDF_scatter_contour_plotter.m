% $Id$
function PDF_scatter_contour_plotter( input_file_sam, input_file_clubb, ...
                                      num_clubb_files )

% Input Variables:
%
% 1) input_file_sam:  A string that contains the path-and-filename of the
%                     SAM LES 3D NetCDF file.  It can also have a value of
%                     'default', in which case the SAM file is set to the
%                     default path-and-filename used below.
%
% 2) input_file_clubb:  A string that contains the path(s)-and-filename(s)
%                       of the CLUBB zt NetCDF file(s).  When more than one
%                       CLUBB file is entered, the files need to be listed
%                       in a single string with the files delimited by a
%                       single space.  This is done this way because this
%                       is the way that the bash script passes arrays.
%                       This string is split up into individual filenames
%                       in the code below.  Additionally, the value of
%                       input_file_clubb can be set to 'default', in which
%                       case the CLUBB file is set to the default
%                       path-and-filename used below.
%
% 3) num_clubb_files:  The number of CLUBB zt NetCDF files found in
%                      input_file_clubb.  When run with the bash script
%                      (run_scatter_contour_plots.bash), this is
%                      automatically calculated.  Otherwise, this number
%                      needs to be passed into this function.

% SAM LES 3D NetCDF filename.
if ( strcmp( input_file_sam, 'default' ) )
   filename_sam = '../../output/LES_output/RICO_128x128x100_drizzle_128_0000255600_micro.nc';
else
   filename_sam = input_file_sam;
end

% CLUBB zt NetCDF filename.
if ( strcmp( input_file_clubb, 'default' ) )

   % The value in input_file_clubb is 'default'.
   filename_clubb = '../../output/rico_zt.nc';

else

   % Filenames (and paths) are found in input_file_clubb.

   if ( num_clubb_files > 1 )

      % There multiple CLUBB zt files named in input_file_clubb.

      % Loop over all characters in the CLUBB input file string to find
      % spaces, which are the filename delimiters.
      string_length = zeros( num_clubb_files, 1 );
      str_len_idx = 0;
      delim_idx_prev = 0;
      for char_idx = 1:1:size( input_file_clubb, 2 )-1
         if ( isspace( input_file_clubb(1,char_idx) ) )
            str_len_idx = str_len_idx + 1;
            string_length(str_len_idx) = char_idx - delim_idx_prev - 1;
            delim_idx_prev = char_idx;
         end
      end % char_idx = 1:1:size( input_file_clubb, 2 )-1
      string_length(num_clubb_files) ...
      = size( input_file_clubb, 2 ) - delim_idx_prev;

      % Setup size for filename_clubb by longest string.
      char_max_len = blanks( max( string_length ) );

      % Initialize filename_clubb.
      filename_clubb(1,:) = char_max_len;

      % Enter filenames into filename_clubb.
      remaining_string = input_file_clubb;
      for clubb_idx = 1:1:num_clubb_files-1
         [ filename_clubb(clubb_idx,1:string_length(clubb_idx)) ...
           remain ] = strtok( remaining_string, ' ' );
         remaining_string = remain;
      end % clubb_idx = 1:1:num_clubb_files
      filename_clubb(num_clubb_files,1:string_length(num_clubb_files)) ...
      = strtok( remaining_string, ' ' );

   else % num_clubb_files = 1

      % There is only one CLUBB zt file named in input_file_clubb.
      filename_clubb = input_file_clubb;

   end % num_clubb_files > 1

end % strcmp( input_file_clubb, 'default' )

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
plot_rt      = true;
plot_thl     = true;
plot_chi     = true;
plot_eta     = true;
plot_rr      = true;
plot_Nr      = true;
plot_ln_rr   = true;
plot_ln_Nr   = true;

% Flag to plot marginal plots for w, rt, thl, chi, eta, rr, and Nr with a
% log-scaled y-axis (the value of P(x)).  This flag is for the aforementioned
% single-variable plots only.
log_Px_plot = false;

% Select the type of output file for MATLAB.
output_type = '-dpng';

% Text for legends of plots.
legend_text(1,1:3) = 'LES';
% Flag to automatically generate legend text for CLUBB files.  CLUBB files
% are listed as a generic 'CLUBB 1', 'CLUBB 2', etc., when there are
% multiple CLUBB files and as 'CLUBB' when there is a single CLUBB file.
% Otherwise, legend text must be manually entered below.
auto_legend_text = true;
if ( auto_legend_text )
   if ( num_clubb_files > 1 )
      for clubb_idx = 1:1:num_clubb_files
         legend_text(clubb_idx+1,1:6+size(int2str(clubb_idx),2)) ...
         = [ 'CLUBB ', int2str(clubb_idx) ];
      end % clubb_idx = 1:1:num_clubb_files
   else % num_clubb_files = 1
      legend_text(2,1:5) = 'CLUBB';
   end % num_clubb_files > 1
else
   % Customize legend text for CLUBB files.
   % The auto_legend_text flag must be turned off to do this.
   % The first CLUBB file has an index of 2.
   legend_text(2,1:9) = 'CLUBB DDL';
   legend_text(3,1:8) = 'CLUBB DL';
   legend_text(4,1:8) = 'CLUBB SL';
end % auto_legend_text

% Assign line colors, line styles, and line widths for CLUBB lines in
% plots.  The LES uses solid red lines that have a width of 3 (when the
% histograms are not in use).
for clubb_idx = 1:1:num_clubb_files
   % CLUBB line colors.
   if ( mod( clubb_idx, 5 ) == 1 )
      clubb_linecolor(clubb_idx) = 'b';
   elseif ( mod( clubb_idx, 5 ) == 2 )
      clubb_linecolor(clubb_idx) = 'g';
   elseif ( mod( clubb_idx, 5 ) == 3 )
      clubb_linecolor(clubb_idx) = 'm';
   elseif ( mod( clubb_idx, 5 ) == 4 )
      clubb_linecolor(clubb_idx) = 'c';
   elseif ( mod( clubb_idx, 5 ) == 0 )
      clubb_linecolor(clubb_idx) = 'k';
   end
   % CLUBB line styles.
   if ( mod( clubb_idx, 4 ) == 1 )
      clubb_linestyle(clubb_idx,1:1) = '-';
   elseif ( mod( clubb_idx, 4 ) == 2 )
      clubb_linestyle(clubb_idx,1:2) = '--';
   elseif ( mod( clubb_idx, 4 ) == 3 )
      clubb_linestyle(clubb_idx,1:2) = '-.';
   elseif ( mod( clubb_idx, 4 ) == 0 )
      clubb_linestyle(clubb_idx,1:1) = ':';
   end
   % CLUBB line widths.
   if ( clubb_idx == 1 )
      clubb_linewidth(clubb_idx) = 2;
   elseif ( clubb_idx == 2 )
      clubb_linewidth(clubb_idx) = 1.5;
   else
      clubb_linewidth(clubb_idx) = 1;
   end
end % clubb_idx = 1:1:num_clubb_files

%==========================================================================

% SAM LES 3D file variable indices
global idx_3D_w
global idx_3D_rt
global idx_3D_thl
global idx_3D_chi
global idx_3D_eta
global idx_3D_rr
global idx_3D_Nr

% Read SAM NetCDF file and obtain variables.
[ z_sam, time_sam, var_sam, units_corrector_type_sam, ...
  nx_sam, ny_sam, nz_sam, num_t_sam, num_var_sam ] ...
= read_SAM_3D_file( filename_sam );

for clubb_idx = 1:1:num_clubb_files

   % Read CLUBB zt NetCDF files and obtain variables.
   [ z_clubb_out, time_clubb_out, var_clubb_out,...
     units_corrector_type_clubb_out, nz_clubb_out, ...
     num_t_clubb_out, num_var_clubb_out, varid_clubb_out ] ...
   = read_CLUBB_file( strtrim( filename_clubb(clubb_idx,:) ) );

   % Enter the number of vertical grid levels into an array.
   nz_clubb_file(clubb_idx) = nz_clubb_out;
   % Check that all CLUBB files have the same number of vertical grid
   % levels.
   if ( num_clubb_files > 1 )
      if ( nz_clubb_file(clubb_idx) ~= nz_clubb_file(1) )
         fprintf( [ 'Error:  CLUBB files do not have the same number', ...
                    ' of vertical grid levels.\n', ...
                    'File:  %s;  Number of grid levels:  %d\n', ...
                    'File:  %s;  Number of grid levels:  %d\n' ], ...
                  strtrim( filename_clubb(1,:) ), nz_clubb_file(1), ...
                  strtrim( filename_clubb(clubb_idx,:) ), ...
                  nz_clubb_file(clubb_idx) );
         fprintf( '\n' );
         quit force
      end % nz_clubb_file(clubb_idx) ~= nz_clubb_file(1)
   end % num_clubb_files > 1

   % Enter the values of altitudes into an array.
   z_clubb_file(clubb_idx,:) = z_clubb_out;
   % Check that all CLUBB files have the same values of altitude.
   if ( num_clubb_files > 1 )
      if ( any( z_clubb_file(clubb_idx,:) ~= z_clubb_file(1,:) ) )
         fprintf( [ 'Error:  CLUBB files do not have the same values', ...
                    ' of altitude.\n' ] );
         fprintf( '\n' );
         quit force
      end % z_clubb_file(clubb_idx,:) ~= z_clubb_file(1,:)
   end % num_clubb_files > 1

   % Enter the number of statistical output timesteps into an array.
   num_t_clubb_file(clubb_idx) = num_t_clubb_out;
   % Check that all CLUBB files have the same number of output timesteps.
   if ( num_clubb_files > 1 )
      if ( num_t_clubb_file(clubb_idx) ~= num_t_clubb_file(1) )
         fprintf( [ 'Error:  CLUBB files do not have the same number', ...
                    ' of output timesteps.\n', ...
                    'File:  %s;  Number of output timesteps:  %d\n', ...
                    'File:  %s;  Number of output timesteps:  %d\n' ], ...
                  strtrim( filename_clubb(1,:) ), num_t_clubb_file(1), ...
                  strtrim( filename_clubb(clubb_idx,:) ), ...
                  num_t_clubb_file(clubb_idx) );
         fprintf( '\n' );
         quit force
      end % num_t_clubb_file(clubb_idx) ~= num_t_clubb_file(1)
   end % num_clubb_files > 1

   % Enter the values of statistical output times into an array.
   time_clubb_file(clubb_idx,:) = time_clubb_out;
   % Check that all CLUBB files have the same statistical output times.
   if ( num_clubb_files > 1 )
      if ( any( time_clubb_file(clubb_idx,:) ~= time_clubb_file(1,:) ) )
         fprintf( [ 'Error:  CLUBB files do not have the same', ...
                    ' statistical output times.\n' ] );
         fprintf( '\n' );
         quit force
      end % time_clubb_file(clubb_idx,:) ~= time_clubb_file(1,:)
   end % num_clubb_files > 1

   % Since we define the number of relevant CLUBB variales (num_var_clubb)
   % in output_vars_clubb, num_var_clubb will be the same for all CLUBB
   % output files.  However, we can compare variables I.D.s (varid_clubb)
   % between files to check that the relevant variables have the same
   % indices.
   varid_clubb_file(clubb_idx,:) = varid_clubb_out;
   % Check that all CLUBB files have the same relevant variables listed
   % in the same order in the output array.
   if ( num_clubb_files > 1 )
      if ( any( varid_clubb_file(clubb_idx,:) ~= varid_clubb_file(1,:) ) )
         fprintf( [ 'Error:  CLUBB files do not have the same', ...
                    ' variables listed in the same order.\n' ] );
         fprintf( '\n' );
         quit force
      end % varid_clubb_file(clubb_idx,:) ~= varid_clubb_file(1,:)
   end % num_clubb_files > 1

   % Enter the number of relevant clubb variables into an array.
   % This will be the same between all CLUBB files.
   num_var_clubb_file(clubb_idx) = num_var_clubb_out;

   % Enter the values of units_corrector_type into an array.
   % This will be the same between all CLUBB files.
   units_corrector_type_clubb_file(clubb_idx,:) ...
      = units_corrector_type_clubb_out;

   % Enter the values of relevant CLUBB variables into an array.
   var_clubb(clubb_idx,:,:,:,:,:) = var_clubb_out;

end % clubb_idx = 1:1:num_clubb_files

% The following CLUBB variables are the same between all the CLUBB files.
% Set them equal to the values found in CLUBB output file 1.
nz_clubb = nz_clubb_file(1);
z_clubb = z_clubb_file(1,:);
num_t_clubb = num_t_clubb_file(1);
time_clubb = time_clubb_file(1,:);
num_var_clubb = num_var_clubb_file(1);
units_corrector_type_clubb = units_corrector_type_clubb_file(1,:);

% Use appropriate units (SI units).
[ var_sam ] ...
   = unit_corrector( num_var_sam, var_sam, units_corrector_type_sam, -1 );

for clubb_idx = 1:1:num_clubb_files

   var_clubb_inout = reshape( var_clubb(clubb_idx,:,:,:,:,:), ...
                              num_var_clubb, 1, 1, nz_clubb, num_t_clubb );

   [ var_clubb_inout ] ...
   = unit_corrector( num_var_clubb, var_clubb_inout, ...
                     units_corrector_type_clubb, -1 );

   var_clubb(clubb_idx,:,:,:,:,:) = var_clubb_inout;

end % clubb_idx = 1:1:num_clubb_files

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
for clubb_idx = 1:1:num_clubb_files

   var_clubb_in = reshape( var_clubb(clubb_idx,:,:,:,:,:), ...
                           num_var_clubb, 1, 1, nz_clubb, num_t_clubb );

   [ mu_w_1(clubb_idx), mu_w_2(clubb_idx), mu_rt_1(clubb_idx), ...
     mu_rt_2(clubb_idx), mu_thl_1(clubb_idx), mu_thl_2(clubb_idx), ...
     mu_chi_1(clubb_idx), mu_chi_2(clubb_idx), mu_eta_1(clubb_idx), ...
     mu_eta_2(clubb_idx), mu_rr_1_n(clubb_idx), mu_rr_2_n(clubb_idx), ...
     mu_Nr_1_n(clubb_idx), mu_Nr_2_n(clubb_idx), ...
     sigma_w_1(clubb_idx), sigma_w_2(clubb_idx), ...
     sigma_rt_1(clubb_idx), sigma_rt_2(clubb_idx), ...
     sigma_thl_1(clubb_idx), sigma_thl_2(clubb_idx), ...
     sigma_chi_1(clubb_idx), sigma_chi_2(clubb_idx), ...
     sigma_eta_1(clubb_idx), sigma_eta_2(clubb_idx), ...
     sigma_rr_1_n(clubb_idx), sigma_rr_2_n(clubb_idx), ...
     sigma_Nr_1_n(clubb_idx), sigma_Nr_2_n(clubb_idx), ...
     corr_rt_thl_1(clubb_idx), corr_rt_thl_2(clubb_idx), ...
     corr_chi_eta_1(clubb_idx), corr_chi_eta_2(clubb_idx), ...
     corr_w_rr_1_n(clubb_idx), corr_w_rr_2_n(clubb_idx), ...
     corr_w_Nr_1_n(clubb_idx), corr_w_Nr_2_n(clubb_idx), ...
     corr_chi_rr_1_n(clubb_idx), corr_chi_rr_2_n(clubb_idx), ...
     corr_chi_Nr_1_n(clubb_idx), corr_chi_Nr_2_n(clubb_idx), ...
     corr_eta_rr_1_n(clubb_idx), corr_eta_rr_2_n(clubb_idx), ...
     corr_eta_Nr_1_n(clubb_idx), corr_eta_Nr_2_n(clubb_idx), ...
     corr_rr_Nr_1_n(clubb_idx), corr_rr_Nr_2_n(clubb_idx), ...
     mixt_frac(clubb_idx), precip_frac_1(clubb_idx), ...
     precip_frac_2(clubb_idx), sigma_sqd_w(clubb_idx) ] ...
   = unpack_CLUBB_vars( var_clubb_in, clubb_height_idx, clubb_time_idx );
   
end % clubb_idx = 1:1:num_clubb_files

%==========================================================================

% Calculations of distribution properties (skewness, etc.) from both SAM
% and CLUBB.

% Calculate SAM higher-order moments and other statistics by analyzing
% the 3D output data.
[ wp2_sam, wp3_sam, rtp2_sam, rtp3_sam, thlp2_sam, thlp3_sam, ...
  wprtp_sam, wpthlp_sam, rtpthlp_sam, ...
  Skw_sam, Skrt_sam, Skthl_sam, ...
  corr_w_rt_ov_sam, corr_w_thl_ov_sam, corr_rt_thl_ov_sam, ...
  precip_frac_sam, mean_rr_sam, rrp2_sam, rrp3_sam, Skrr_sam, ...
  mean_rr_ip_sam, rrp2_ip_sam, rrp3_ip_sam, ...
  Skrr_ip_sam, voms_ip_rr_sam, mean_lnrr_sam, ...
  lnrrp2_sam, lnrrp3_sam, Sk_lnrr_sam, ...
  mean_Nr_sam, Nrp2_sam, Nrp3_sam, SkNr_sam, ...
  mean_Nr_ip_sam, Nrp2_ip_sam, Nrp3_ip_sam, ...
  SkNr_ip_sam, voms_ip_Nr_sam, mean_lnNr_sam, ...
  lnNrp2_sam, lnNrp3_sam, Sk_lnNr_sam ] ...
= calc_stats_SAM( sam_var_lev, nx_sam, ny_sam );

% Calculate CLUBB higher-order moments and other statistics by integration
% over the PDF for each CLUBB data set (file).
for clubb_idx = 1:1:num_clubb_files

   [ wp2_clubb_pdf(clubb_idx), wp3_clubb_pdf(clubb_idx), ...
     rtp2_clubb_pdf(clubb_idx), rtp3_clubb_pdf(clubb_idx), ...
     thlp2_clubb_pdf(clubb_idx), thlp3_clubb_pdf(clubb_idx), ...
     Skw_clubb_pdf(clubb_idx), Skrt_clubb_pdf(clubb_idx), ...
     Skthl_clubb_pdf(clubb_idx), precip_frac_clubb(clubb_idx), ...
     rrm_clubb_pdf(clubb_idx), rrp2_clubb_pdf(clubb_idx), ...
     rrp3_clubb_pdf(clubb_idx), Skrr_clubb_pdf(clubb_idx), ...
     rrm_ip_clubb_pdf(clubb_idx), rrp2_ip_clubb_pdf(clubb_idx), ...
     rrp3_ip_clubb_pdf(clubb_idx), Skrr_ip_clubb_pdf(clubb_idx), ...
     voms_ip_rr_clubb_pdf(clubb_idx), mean_lnrr_clubb_pdf(clubb_idx), ...
     lnrrp2_clubb_pdf(clubb_idx), lnrrp3_clubb_pdf(clubb_idx), ...
     Sk_lnrr_clubb_pdf(clubb_idx), ...
     Nrm_clubb_pdf(clubb_idx), Nrp2_clubb_pdf(clubb_idx), ...
     Nrp3_clubb_pdf(clubb_idx), SkNr_clubb_pdf(clubb_idx), ...
     Nrm_ip_clubb_pdf(clubb_idx), Nrp2_ip_clubb_pdf(clubb_idx), ...
     Nrp3_ip_clubb_pdf(clubb_idx), SkNr_ip_clubb_pdf(clubb_idx), ...
     voms_ip_Nr_clubb_pdf(clubb_idx), mean_lnNr_clubb_pdf(clubb_idx), ...
     lnNrp2_clubb_pdf(clubb_idx), lnNrp3_clubb_pdf(clubb_idx), ...
     Sk_lnNr_clubb_pdf(clubb_idx) ] ...
   = calc_stats_CLUBB( mu_w_1(clubb_idx), mu_w_2(clubb_idx), ...
                       mu_rt_1(clubb_idx), mu_rt_2(clubb_idx), ...
                       mu_thl_1(clubb_idx), mu_thl_2(clubb_idx), ...
                       mu_rr_1_n(clubb_idx), mu_rr_2_n(clubb_idx), ...
                       mu_Nr_1_n(clubb_idx), mu_Nr_2_n(clubb_idx), ...
                       sigma_w_1(clubb_idx), sigma_w_2(clubb_idx), ...
                       sigma_rt_1(clubb_idx), sigma_rt_2(clubb_idx), ...
                       sigma_thl_1(clubb_idx), sigma_thl_2(clubb_idx), ...
                       sigma_rr_1_n(clubb_idx), sigma_rr_2_n(clubb_idx),...
                       sigma_Nr_1_n(clubb_idx), sigma_Nr_2_n(clubb_idx),...
                       mixt_frac(clubb_idx), precip_frac_1(clubb_idx), ...
                       precip_frac_2(clubb_idx) );

end % clubb_idx = 1:1:num_clubb_files

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

% Calculations based on CLUBB PDF parameters from this level.
fprintf( '\n' )
fprintf( 'Calculations based on CLUBB PDF parameters at this level:\n' )
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB PDF wp2 = %g\n', wp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF wp3 = %g\n', wp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF rtp2 = %g\n', rtp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF rtp3 = %g\n', rtp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF thlp2 = %g\n', thlp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF thlp3 = %g\n', thlp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness of w = %g\n', Skw_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness of rt = %g\n', Skrt_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness of thl = %g\n', ...
            Skthl_clubb_pdf(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files

% Backsolve for the appropriate values of beta (for each of rt and thl)
% using SAM LES data (and CLUBB output for sigma_sqd_w) based on Larson
% and Golaz (2005), Eq. 33.
%fprintf( '\n' )
%fprintf( 'Backsolve for the parameter beta based largely on SAM LES:\n' )
%for clubb_idx = 1:1:num_clubb_files
%
%   if ( num_clubb_files > 1 )
%      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
%               strtrim( filename_clubb(clubb_idx,:) ) )
%   end % num_clubb_files > 1
%
%   Skw_hat = Skw_sam / ( 1.0 - sigma_sqd_w(clubb_idx) )^1.5;
%   corr_w_rt_ov_hat ...
%   = corr_w_rt_ov_sam / sqrt( 1.0 - sigma_sqd_w(clubb_idx) );
%   corr_w_thl_ov_hat ...
%   = corr_w_thl_ov_sam / sqrt( 1.0 - sigma_sqd_w(clubb_idx) );
%
%   beta_rt ...
%   = ( Skrt_sam / ( Skw_hat * corr_w_rt_ov_hat ) ...
%       - corr_w_rt_ov_hat^2 ) ...
%     / ( 1.0 - corr_w_rt_ov_hat^2 );
%
%   beta_thl ...
%   = ( Skthl_sam / ( Skw_hat * corr_w_thl_ov_hat ) ...
%       - corr_w_thl_ov_hat^2 ) ...
%     / ( 1.0 - corr_w_thl_ov_hat^2 );
%
%   % Calculations based on CLUBB PDF parameters from this level.
%   fprintf( 'Ideal beta for skewness of rt = %g\n', beta_rt );
%   fprintf( 'Ideal beta for skewness of thl = %g\n', beta_thl );
%
%end % clubb_idx = 1:1:num_clubb_files

% Print the values of hydrometeor statistics.
% Note:  v.o.m.s. stands for variance over mean-squared.
fprintf( '\n' )
fprintf( 'Calculations for hydrometeors at this level:\n' )
fprintf( '\n' )
fprintf( 'SAM precipitation fraction = %g\n', precip_frac_sam );
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB precipitation fraction = %g\n', ...
            precip_frac_clubb(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files
fprintf( '\n' )
fprintf( 'Rain water mixing ratio, rr\n' )
fprintf( 'Overall values:\n' )
fprintf( 'SAM rrm = %g\n', mean_rr_sam );
fprintf( 'SAM rrp2 = %g\n', rrp2_sam );
fprintf( 'SAM rrp3 = %g\n', rrp3_sam );
fprintf( 'SAM skewness of rr = %g\n', Skrr_sam );
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB PDF rrm = %g\n', rrm_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF rrp2 = %g\n', rrp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF rrp3 = %g\n', rrp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness of rr = %g\n', Skrr_clubb_pdf(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files
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
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB PDF mean (i.p.) of rr = %g\n', ...
            rrm_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF variance (i.p.) of rr = %g\n', ...
            rrp2_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF 3rd moment (i.p.) of rr = %g\n', ...
            rrp3_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness (i.p.) of rr = %g\n', ...
            Skrr_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF v.o.m.s (i.p.) of rr = %g\n', ...
            voms_ip_rr_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF mean (i.p.) of ln rr = %g\n', ...
            mean_lnrr_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF variance (i.p.) of ln rr = %g\n', ...
            lnrrp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF 3rd moment (i.p.) of ln rr = %g\n', ...
            lnrrp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness (i.p.) of ln rr = %g\n', ...
            Sk_lnrr_clubb_pdf(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files
fprintf( '\n' )
fprintf( 'Rain drop concentration, Nr\n' )
fprintf( 'Overall values:\n' )
fprintf( 'SAM Nrm = %g\n', mean_Nr_sam );
fprintf( 'SAM Nrp2 = %g\n', Nrp2_sam );
fprintf( 'SAM Nrp3 = %g\n', Nrp3_sam );
fprintf( 'SAM skewness of Nr = %g\n', SkNr_sam );
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB PDF Nrm = %g\n', Nrm_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF Nrp2 = %g\n', Nrp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF Nrp3 = %g\n', Nrp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness of Nr = %g\n', SkNr_clubb_pdf(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files
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
for clubb_idx = 1:1:num_clubb_files
   if ( num_clubb_files > 1 )
      fprintf( '*** CLUBB file %d:  %s\n', clubb_idx, ...
               strtrim( filename_clubb(clubb_idx,:) ) )
   end % num_clubb_files > 1
   fprintf( 'CLUBB PDF mean (i.p.) of Nr = %g\n', ...
            Nrm_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF variance (i.p.) of Nr = %g\n', ...
            Nrp2_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF 3rd moment (i.p.) of Nr = %g\n', ...
            Nrp3_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness (i.p.) of Nr = %g\n', ...
            SkNr_ip_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF v.o.m.s (i.p.) of Nr = %g\n', ...
            voms_ip_Nr_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF mean (i.p.) of ln Nr = %g\n', ...
            mean_lnNr_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF variance (i.p.) of ln Nr = %g\n', ...
            lnNrp2_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF 3rd moment (i.p.) of ln Nr = %g\n', ...
            lnNrp3_clubb_pdf(clubb_idx) );
   fprintf( 'CLUBB PDF skewness (i.p.) of ln Nr = %g\n', ...
            Sk_lnNr_clubb_pdf(clubb_idx) );
end % clubb_idx = 1:1:num_clubb_files

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
   plot_rr     = false;
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
   plot_Nr     = false;
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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              'r_{t}    [kg/kg]', '\theta_{l}    [K]', ...
                              [ '\bf ', casename ], 'r_{t} vs. \theta_{l}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_rt_thl_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              '\chi    [kg/kg]', '\eta    [kg/kg]', ...
                              [ '\bf ', casename ], '\chi vs. \eta', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_eta_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              'w    [m/s]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], 'w vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_w_rr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              'w    [m/s]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], 'w vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_w_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              '\chi    [kg/kg]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], '\chi vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_rr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              '\chi    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], '\chi vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_chi_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              '\eta    [kg/kg]', 'r_{r}    [kg/kg]', ...
                              [ '\bf ', casename ], '\eta vs. r_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_eta_rr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              '\eta    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], '\eta vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_eta_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

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
                              num_clubb_files, ...
                              clubb_linecolor, clubb_linestyle, ...
                              clubb_linewidth, legend_text, ...
                              'r_{r}    [kg/kg]', 'N_{r}    [num/kg]', ...
                              [ '\bf ', casename ], 'r_{r} vs. N_{r}', ...
                              [ 'Time = ', print_time, ' minutes' ], ...
                              [ 'Altitude = ', print_alt, ' meters' ], ...
                              print_note )

   output_filename = [ 'output/', casename, '_rr_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_rr_Nr

% Plot the CLUBB PDF and LES points for w.
if ( plot_w )

   fprintf( 'Plotting marginal plot for w\n' );

   % The number of w bins/points to plot.
   num_w_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_w,:), ...
                             nx_sam, ny_sam, ...
                             num_w_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_w_1, mu_w_2, ...
                             sigma_w_1, sigma_w_2, ...
                             mixt_frac, 'w    [m/s]', ...
                             [ '\bf ', casename ], 'w', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_w_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_w

% Plot the CLUBB PDF and LES points for rt.
if ( plot_rt )

   fprintf( 'Plotting marginal plot for rt\n' );

   % The number of rt bins/points to plot.
   num_rt_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_rt,:), ...
                             nx_sam, ny_sam, ...
                             num_rt_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_rt_1, mu_rt_2, ...
                             sigma_rt_1, sigma_rt_2, ...
                             mixt_frac, 'r_{t}    [kg/kg]', ...
                             [ '\bf ', casename ], 'r_{t}', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_rt_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_rt

% Plot the CLUBB PDF and LES points for theta-l.
if ( plot_thl )

   fprintf( 'Plotting marginal plot for thl\n' );

   % The number of thl bins/points to plot.
   num_thl_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_thl,:), ...
                             nx_sam, ny_sam, ...
                             num_thl_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_thl_1, mu_thl_2, ...
                             sigma_thl_1, sigma_thl_2, ...
                             mixt_frac, '\theta_{l}    [K]', ...
                             [ '\bf ', casename ], '\theta_{l}', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_thl_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_thl

% Plot the CLUBB PDF and LES points for chi.
if ( plot_chi )

   fprintf( 'Plotting marginal plot for chi\n' );

   % The number of chi bins/points to plot.
   num_chi_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_chi,:), ...
                             nx_sam, ny_sam, ...
                             num_chi_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_chi_1, mu_chi_2, ...
                             sigma_chi_1, sigma_chi_2, ...
                             mixt_frac, '\chi    [kg/kg]', ...
                             [ '\bf ', casename ], '\chi', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_chi_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_chi

% Plot the CLUBB PDF and LES points for eta.
if ( plot_eta )

   fprintf( 'Plotting marginal plot for eta\n' );

   % The number of eta bins/points to plot.
   num_eta_pts = 100;

   plot_CLUBB_PDF_LES_pts_N( sam_var_lev(idx_3D_eta,:), ...
                             nx_sam, ny_sam, ...
                             num_eta_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_eta_1, mu_eta_2, ...
                             sigma_eta_1, sigma_eta_2, ...
                             mixt_frac, '\eta    [kg/kg]', ...
                             [ '\bf ', casename ], '\eta', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_eta_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_eta

% Plot the CLUBB PDF and LES points for rr.
if ( plot_rr )

   fprintf( 'Plotting marginal plot for rr\n' );

   % The number of rr bins/points to plot.
   num_rr_pts = 100;

   plot_CLUBB_PDF_LES_pts_L( sam_var_lev(idx_3D_rr,:), ...
                             nx_sam, ny_sam, ...
                             num_rr_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_rr_1_n, mu_rr_2_n, ...
                             sigma_rr_1_n, sigma_rr_2_n, ...
                             precip_frac_1, precip_frac_2, ...
                             mixt_frac, 'r_{r}    [kg/kg]', ...
                             [ '\bf ', casename ], 'r_{r}', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_rr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_rr

% Plot the CLUBB PDF and LES points for Nr.
if ( plot_Nr )

   fprintf( 'Plotting marginal plot for Nr\n' );

   % The number of Nr bins/points to plot.
   num_Nr_pts = 100;

   plot_CLUBB_PDF_LES_pts_L( sam_var_lev(idx_3D_Nr,:), ...
                             nx_sam, ny_sam, ...
                             num_Nr_pts, log_Px_plot, ...
                             num_clubb_files, ...
                             clubb_linecolor, clubb_linestyle, ...
                             clubb_linewidth, legend_text, ...
                             mu_Nr_1_n, mu_Nr_2_n, ...
                             sigma_Nr_1_n, sigma_Nr_2_n, ...
                             precip_frac_1, precip_frac_2, ...
                             mixt_frac, 'N_{r}    [num/kg]', ...
                             [ '\bf ', casename ], 'N_{r}', ...
                             [ 'Time = ', print_time, ' minutes' ], ...
                             [ 'Altitude = ', print_alt, ' meters' ], ...
                             print_note )

   output_filename = [ 'output/', casename, '_Nr_z', ...
                       print_alt, '_t', print_time ];

   print( output_type, output_filename );

end % plot_Nr

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
                                num_clubb_files, ...
                                clubb_linecolor, clubb_linestyle, ...
                                clubb_linewidth, legend_text, ...
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

   print( output_type, output_filename );

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
                                num_clubb_files, ...
                                clubb_linecolor, clubb_linestyle, ...
                                clubb_linewidth, legend_text, ...
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

   print( output_type, output_filename );

end % plot_ln_Nr
