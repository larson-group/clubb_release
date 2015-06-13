% $Id$
function [ sam_var_lev ] ...
= interp_SAM_3D_points( z_clubb, clubb_height_idx, ...
                        z_sam, nx_sam, ny_sam, nz_sam, ...
                        var_sam, num_var_sam );

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
