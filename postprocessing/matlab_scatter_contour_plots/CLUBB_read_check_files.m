% $Id$
function [ nz_clubb, z_clubb, num_t_clubb, time_clubb, num_var_clubb, ...
           units_corrector_type_clubb, var_clubb ] ...
= CLUBB_read_check_files( filename_clubb, num_clubb_files )

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
