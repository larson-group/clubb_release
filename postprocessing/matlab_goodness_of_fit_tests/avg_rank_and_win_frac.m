% $Id$
function [ avg_rank_stat_x win_frac_stat_x ] ...
= avg_rank_and_win_frac( stat_number_x, num_clubb_files )

% Find the number of statistical output times involved.
num_times = size( stat_number_x, 2 );

% Find the number of statistical output altitudes involved.
num_heights = size( stat_number_x, 3 );

% Initialize storage arrays.
sort_stat_number_x = -ones( num_clubb_files, 1 );
sort_stat_index = -ones( num_clubb_files, 1 );
rank = -ones( num_clubb_files, num_times, num_heights );
win = -ones( num_clubb_files, num_times, num_heights );
avg_rank_stat_x = -ones( num_clubb_files, 1 );
win_frac_stat_x = -ones( num_clubb_files, 1 );

% Loop over each statistical output time and altitude.
for time_idx = 1:1:num_times
   for height_idx = 1:1:num_heights

      % Sort the CLUBB files in rank order based on the stat number for x.
      [ sort_stat_number_x, sort_stat_index ] ...
      = sort( stat_number_x(:,time_idx,height_idx) );

      if ( all( sort_stat_number_x == sort_stat_number_x(1) ) )

         % All CLUBB files have the same value.
         rank(:,time_idx,height_idx) = -1;
         win(:,time_idx,height_idx)  = -1;

      else

         % Store the rank for each CLUBB file.
         for clubb_idx = 1:1:num_clubb_files
            for rank_idx = 1:1:num_clubb_files
               if ( clubb_idx == sort_stat_index(rank_idx) )
                  rank(clubb_idx,time_idx,height_idx) = rank_idx;
                  break
               end % clubb_idx == sort_stat_index(rank_idx)
            end % rank_idx = 1:1:num_clubb_files
         end % clubb_idx = 1:1:num_clubb_files

         % Set the value of win to 1 if the CLUBB file is ranked 1st.
         % Otherwise, set the value of win to 0.
         for clubb_idx = 1:1:num_clubb_files
            if ( rank(clubb_idx,time_idx,height_idx) == 1 )
               win(clubb_idx,time_idx,height_idx) = 1.0;
            else
               win(clubb_idx,time_idx,height_idx) = 0.0;
            end
         end % clubb_idx = 1:1:num_clubb_files
      end % all( sort_stat_number_x == sort_stat_number_x(1) )

   end % height_idx = 1:1:num_heights
end % time_idx = 1:1:num_times

% Find the average rank and winning fraction for all CLUBB runs for the
% stat number for x.
for clubb_idx = 1:1:num_clubb_files

   % Calculate the average rank and winning fraction for each CLUBB file.
   sum_rank = 0.0;
   sum_win  = 0.0;
   num_rank = 0;
   num_win  = 0;
   for time_idx = 1:1:num_times
      for height_idx = 1:1:num_heights

         if ( rank(clubb_idx,time_idx,height_idx) ~= -1.0 )
            % Include the rank value from this particular level and time in
            % the average for this CLUBB file.
            sum_rank = sum_rank + rank(clubb_idx,time_idx,height_idx);
            num_rank = num_rank + 1;
         end % rank(clubb_idx,time_idx,height_idx) ~= -1.0

         if ( win(clubb_idx,time_idx,height_idx) ~= -1.0 )
            % Include the win value from this particular level and time in
            % the average for this CLUBB file.
            sum_win = sum_win + win(clubb_idx,time_idx,height_idx);
            num_win = num_win + 1;
         end % rank(clubb_idx,time_idx,height_idx) ~= -1.0

      end % height_idx = 1:1:num_heights
   end % time_idx = 1:1:num_times

   if ( num_rank > 0 )
      avg_rank_stat_x(clubb_idx) = sum_rank / num_rank;
   else % num_rank
      avg_rank_stat_x(clubb_idx) = -1.0;
   end % num_rank > 0

   if ( num_win > 0 )
      win_frac_stat_x(clubb_idx) = sum_win / num_win;
   else % num_win
      win_frac_stat_x(clubb_idx) = -1.0;
   end % num_win > 0

end % clubb_idx = 1:1:num_clubb_files
