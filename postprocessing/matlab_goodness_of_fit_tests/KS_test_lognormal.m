% $Id$
function [ KS_number_x ] ...
= KS_test_lognormal( sam_var_x, mu_x_1_n, mu_x_2_n, sigma_x_1_n, ...
                     sigma_x_2_n, mixt_frac, precip_frac_1, ...
                     precip_frac_2, flag_ip_only, num_clubb_files )

% Function that calculates the Kolmogorov-Smirnov (K-S) statistic for
% variable x comparing CLUBB to SAM LES 3D results.  This function is for
% variables that CLUBB uses a sum of two delta-logormal distributions to
% estimate the distribution.

% Sort SAM LES 3D data in ascending order.
sam_var_x_sort = sort( sam_var_x );

% Find the number of SAM LES 3D data points.
num_sam_pts = size( sam_var_x, 2 );

% Find the index in sam_var_x_sort where x is greater than 0.
for idx = 1:1:num_sam_pts
   if ( sam_var_x_sort(idx) > 0.0 )
      first_sam_x_pt_gt_0_idx = idx;
      break
   elseif ( idx == num_sam_pts )
      % The value of sam_var_x is 0 everywhere.
      first_sam_x_pt_gt_0_idx = num_sam_pts + 1;
   end % sam_var_x_sort(idx) > 0
end % idx = 1:1:num_sam_pts

% Calculate CLUBB overall precipitation fraction.
for clubb_idx = 1:1:num_clubb_files
   precip_frac_clubb(clubb_idx) ...
   = mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
     + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx);
end % clubb_idx = 1:1:num_clubb_files

if ( flag_ip_only )

   % Don't include x = 0 in CLUBB's CDF.
   for clubb_idx = 1:1:num_clubb_files
      clubb_cdf_at_0(clubb_idx) = 0.0;
   end % clubb_idx = 1:1:num_clubb_files

else % !~flag_ip_only

   % Calculate the CLUBB CDF at x = 0.
   for clubb_idx = 1:1:num_clubb_files
      clubb_cdf_at_0(clubb_idx) = 1.0 - precip_frac_clubb(clubb_idx);
   end % clubb_idx = 1:1:num_clubb_files

   % Calculate SAM precipitation fraction.
   precip_frac_sam = ( num_sam_pts - first_sam_x_pt_gt_0_idx + 1 ) ...
                     / num_sam_pts;

   % Calculate the K-S statistic at x = 0.
   for clubb_idx = 1:1:num_clubb_files
      KS_number_x_at_0(clubb_idx) ...
      = abs( precip_frac_clubb(clubb_idx) - precip_frac_sam );
   end % clubb_idx = 1:1:num_clubb_files

end % flag_ip_only

% Loop over all CLUBB data sets (files).
for clubb_idx = 1:1:num_clubb_files
   for idx = first_sam_x_pt_gt_0_idx:1:num_sam_pts

      % Find the value of x at the index location.
      var_x = sam_var_x_sort(idx);

      % Calculate the value of the CLUBB CDF at x.
      clubb_cdf ...
      = clubb_cdf_at_0(clubb_idx) ...
        + mixt_frac(clubb_idx) * precip_frac_1(clubb_idx) ...
          * CDF_comp_Lognormal( var_x, mu_x_1_n(clubb_idx), ...
                                sigma_x_1_n(clubb_idx) ) ...
          + ( 1.0 - mixt_frac(clubb_idx) ) * precip_frac_2(clubb_idx) ...
            * CDF_comp_Lognormal( var_x, mu_x_2_n(clubb_idx), ...
                                  sigma_x_2_n(clubb_idx) );

      if ( flag_ip_only )
         clubb_cdf = clubb_cdf / precip_frac_clubb(clubb_idx);
      end % flag_ip_only

      % Find the difference between the CLUBB CDF and the SAM LES CDF at
      % the index point.
      if ( flag_ip_only )
         diff(clubb_idx,idx-first_sam_x_pt_gt_0_idx+1) ...
         = max( abs( clubb_cdf ...
                     - (idx-first_sam_x_pt_gt_0_idx+1) ...
                       /(num_sam_pts-first_sam_x_pt_gt_0_idx+1) ), ...
                abs( clubb_cdf ...
                     - (idx-first_sam_x_pt_gt_0_idx) ...
                       /(num_sam_pts-first_sam_x_pt_gt_0_idx+1) ) );
      else % !~flag_ip_only
         diff(clubb_idx,idx-first_sam_x_pt_gt_0_idx+1) ...
         = max( abs( clubb_cdf - idx/num_sam_pts ), ...
                abs( clubb_cdf - (idx-1)/num_sam_pts ) );
      end % flag_ip_only

   end % idx = 1:1:num_sam_pts

   % The K-S statisitic is the maximum difference between the CLUBB CDF and
   % the SAM LES CDF.
   if ( flag_ip_only )
      KS_number_x(clubb_idx) = max( diff(clubb_idx,:) );
   else % !~flag_ip_only
      KS_number_x(clubb_idx) = max( max( diff(clubb_idx,:) ), ...
                                    KS_number_x_at_0(clubb_idx) );
   end % flag_ip_only

end % clubb_idx = 1:1:num_clubb_files
