% $Id$
function [ clubb_time_idx ] ...
= get_clubb_time_index( time_sam, time_clubb, num_t_clubb )

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
