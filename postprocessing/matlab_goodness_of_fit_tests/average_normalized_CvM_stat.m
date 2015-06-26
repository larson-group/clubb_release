% $Id$
function [ avg_nrmlzed_CvM_number_x ] ...
= average_normalized_CvM_stat( CvM_number_x, num_pts_test_x, ...
                               num_clubb_files )

% Find the number of statistical output times involved.
num_times = size( CvM_number_x, 2 );

% Find the number of statistical output altitudes involved.
num_heights = size( CvM_number_x, 3 );

% Find the average normalized C-vM number for x for each CLUBB run.
% The normalized C-vM number is CvM_number_x / num_pts_test_x.
for clubb_idx = 1:1:num_clubb_files

   % Calculate the average normalized C-vM number of x.
   sum_nrmlzed_CvM = 0.0;
   num_nrmlzed_CvM = 0;
   for time_idx = 1:1:num_times
      for height_idx = 1:1:num_heights
         if ( CvM_number_x(clubb_idx,time_idx,height_idx) ~= -1.0 )
            % Include the normalized C-vM number for x from this particular
            % level and time in the average.
            sum_nrmlzed_CvM ...
            = sum_nrmlzed_CvM ...
              + CvM_number_x(clubb_idx,time_idx,height_idx) ...
                / num_pts_test_x(time_idx,height_idx);
            num_nrmlzed_CvM = num_nrmlzed_CvM + 1;
         end % CvM_number_x(clubb_idx,time_idx,height_idx) ~= -1.0
      end % height_idx = 1:1:num_heights
   end % time_idx = 1:1:num_times

   if ( num_nrmlzed_CvM > 0 )
      avg_nrmlzed_CvM_number_x(clubb_idx) ...
      = sum_nrmlzed_CvM / num_nrmlzed_CvM;
   else % num_nrmlzed_CvM
      avg_nrmlzed_CvM_number_x(clubb_idx) = -1.0;
   end % num_nrmlzed_CvM > 0
   
end % clubb_idx = 1:1:num_clubb_files
