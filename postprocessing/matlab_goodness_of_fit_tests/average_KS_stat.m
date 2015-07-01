% $Id$
function [ avg_KS_number_x ] ...
= average_KS_stat( KS_number_x, num_clubb_files )

% Find the number of statistical output times involved.
num_times = size( KS_number_x, 2 );

% Find the number of statistical output altitudes involved.
num_heights = size( KS_number_x, 3 );

% Find the average K-S number for x for each CLUBB run.
for clubb_idx = 1:1:num_clubb_files

   % Calculate the average K-S number of x.
   sum_KS = 0.0;
   num_KS = 0;
   for time_idx = 1:1:num_times
      for height_idx = 1:1:num_heights
         if ( KS_number_x(clubb_idx,time_idx,height_idx) ~= -1.0 )
            % Include the K-S number for x from this particular level and
            % time in the average.
            sum_KS = sum_KS + KS_number_x(clubb_idx,time_idx,height_idx);
            num_KS = num_KS + 1;
         end % KS_number_x(clubb_idx,time_idx,height_idx) ~= -1.0
      end % height_idx = 1:1:num_heights
   end % time_idx = 1:1:num_times

   if ( num_KS > 0 )
      avg_KS_number_x(clubb_idx) = sum_KS / num_KS;
   else % num_KS
      avg_KS_number_x(clubb_idx) = -1.0;
   end % num_KS > 0
   
end % clubb_idx = 1:1:num_clubb_files
