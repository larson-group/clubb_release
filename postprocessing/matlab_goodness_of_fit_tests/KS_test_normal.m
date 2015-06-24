% $Id$
function [ KS_number_x ] ...
= KS_test_normal( sam_var_x, mu_x_1, mu_x_2, sigma_x_1, ...
                  sigma_x_2, mixt_frac, num_clubb_files )

% Function that calculates the Kolmogorov-Smirnov (K-S) statistic for
% variable x comparing CLUBB to SAM LES 3D results.  This function is for
% variables that CLUBB uses a sum of two normal distributions to estimate
% the distribution.

% Sort SAM LES 3D data in ascending order.
sam_var_x_sort = sort( sam_var_x );

% Find the number of SAM LES 3D data points.
num_sam_pts = size( sam_var_x, 2 );

% If there are an insufficient number of sample points, set the test value
% to -1 and return.
if ( num_sam_pts < 50 )
   KS_number_x(1:num_clubb_files) = -1.0;
   return
end % num_sam_pts < 50

% Loop over all CLUBB data sets (files).
for clubb_idx = 1:1:num_clubb_files
   for idx = 1:1:num_sam_pts

      % Find the value of x at the index location.
      var_x = sam_var_x_sort(idx);

      % Calculate the value of the CLUBB CDF at x.
      clubb_cdf = mixt_frac(clubb_idx) ...
                  * CDF_comp_Normal( var_x, mu_x_1(clubb_idx), ...
                                     sigma_x_1(clubb_idx) ) ...
                  + ( 1.0 - mixt_frac(clubb_idx) ) ...
                    * CDF_comp_Normal( var_x, mu_x_2(clubb_idx), ...
                                       sigma_x_2(clubb_idx) );

      % Find the difference between the CLUBB CDF and the SAM LES CDF at
      % the index point.
      diff(clubb_idx,idx) = max( abs( clubb_cdf - idx/num_sam_pts ), ...
                                 abs( clubb_cdf - (idx-1)/num_sam_pts ) );

   end % idx = 1:1:num_sam_pts

   % The K-S statisitic is the maximum difference between the CLUBB CDF and
   % the SAM LES CDF.
   KS_number_x(clubb_idx) = max( diff(clubb_idx,:) );

end % clubb_idx = 1:1:num_clubb_files
