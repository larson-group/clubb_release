% $Id$
function [ z_clubb, time_clubb, var_clubb, units_corrector_type, ...
           nz_clubb, num_t_clubb, num_var_clubb, varid_clubb ] ...
= read_CLUBB_file( filename_clubb )

% Paths to MEXNC and NetCDF Toolbox utilities.
addpath ( '/usr/share/mexcdf/mexnc' )
addpath ( '/usr/local/netcdf_toolbox/netcdf', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/nctype', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/ncutility', '-end' )

global idx_z
global idx_time

% Open the CLUBB NetCDF file for reading.
[ncid_clubb, status] = mexnc ( 'open', filename_clubb, nc_nowrite_mode );

% Get the file attributes from the CLUBB NetCDF file.
[ dimid_clubb_x,  dimid_clubb_y,  dimid_clubb_z,  dimid_clubb_t, ...
  name_clubb_x,   name_clubb_y,   name_clubb_z,   name_clubb_t, ...
  nx_clubb, ny_clubb, nz_clubb, num_t_clubb ] ...
= get_file_attributes( ncid_clubb, 'longitude', 'latitude', 'altitude', ...
                       'time' );

% Read in CLUBB NetCDF output variable names and assign indices.
[ varname_clubb, units_corrector_type, ...
  num_var_clubb, num_tot_var_clubb ] = output_vars_clubb;

% Get the variable IDs for the variable in question.
varid_clubb = zeros(1,num_tot_var_clubb);
for i = 1:1:num_tot_var_clubb
   [ varid_clubb(i), status ] ...
   = nc_variable_id( ncid_clubb, varname_clubb(i,:) );
   % Print a warning message when the status of the variable id function
   % is not equal to 0, meaning that there is a problem.
   % This is most likely due to a requested variable missing from the
   % CLUBB statistical file.
   if ( status ~= 0 )
      fprintf( 'Error:  variable %s is not found in %s\n', ...
               strtrim( varname_clubb(i,:) ), filename_clubb );
      fprintf( '\n' );
   end
end

% Read CLUBB variables.

% Altitiude (nz_clubb x 1 array)
[ z_clubb, status ] ...
   = mexnc ( 'get_var_double', ncid_clubb, varid_clubb(idx_z) );

% Time (num_t_clubb x 1 array)
[ time_clubb, status ] ...
   = mexnc ( 'get_var_double', ncid_clubb, varid_clubb(idx_time) );

% Variables (nx_clubb x ny_clubb x nz_clubb x num_t_clubb array for each)

% Initialize (we can place variables all variables into a
% num_var_clubb x nx_clubb x ny_clubb x nz_clubb x num_t_clubb array). 
var_clubb = zeros( num_var_clubb, nx_clubb, ny_clubb, nz_clubb, num_t_clubb );

% Read in variables
for i = 1:1:num_var_clubb
   [ var_clubb(i,:,:,:,:), status ] ...
      = mexnc ( 'get_var_double', ncid_clubb, varid_clubb(i) );
end

% Close CLUBB NetCDF file.
status = mexnc ( 'close', ncid_clubb );
