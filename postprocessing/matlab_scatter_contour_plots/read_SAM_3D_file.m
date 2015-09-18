% $Id$
function [ z_sam, time_sam, var_sam, units_corrector_type, ...
           nx_sam, ny_sam, nz_sam, num_t_sam, num_var_sam ] ...
= read_SAM_3D_file( filename_sam, casename )

% Paths to MEXNC and NetCDF Toolbox utilities.
addpath ( '/usr/share/mexcdf/mexnc' )
addpath ( '/usr/local/netcdf_toolbox/netcdf', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/nctype', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/ncutility', '-end' )

global idx_3D_x
global idx_3D_y
global idx_3D_z
global idx_3D_time

% Open the SAM NetCDF file for reading.
[ncid_sam, status] = mexnc ( 'open', filename_sam, nc_nowrite_mode );

% Get the file attributes from the SAM NetCDF file.
[ dimid_sam_x,  dimid_sam_y,  dimid_sam_z,  dimid_sam_t, ...
  name_sam_x,   name_sam_y,   name_sam_z,   name_sam_t, ...
  nx_sam, ny_sam, nz_sam, num_t_sam ] ...
= get_file_attributes( ncid_sam, 'x', 'y', 'z', 'time' );

% Read in SAM NetCDF output variable names and assign indices.
[ varname_sam, units_corrector_type, ...
  num_var_sam, num_tot_var_sam ] = output_vars_3D( casename );

% Get the variable IDs for the variable in question.
varid_sam = zeros(1,num_tot_var_sam);
for i = 1:1:num_tot_var_sam
   [ varid_sam(i), status ] = nc_variable_id( ncid_sam, varname_sam(i,:) );
   % Print a warning message when the status of the variable id function
   % is not equal to 0, meaning that there is a problem.
   % This is most likely due to a requested variable missing from the
   % SAM LES 3D statistical file.
   if ( status ~= 0 )
      fprintf( 'Error:  variable %s is not found in %s\n', ...
               strtrim( varname_sam(i,:) ), filename_sam );
      fprintf( '\n' );
   end
end

% Read SAM variables.

% West-East Distance (nx_sam x 1 array)
[ x_sam, status ] ...
   = mexnc ( 'get_var_double', ncid_sam, varid_sam(idx_3D_x) );

% South-North Distance (ny_sam x 1 array)
[ y_sam, status ] ...
   = mexnc ( 'get_var_double', ncid_sam, varid_sam(idx_3D_y) );

% Altitiude (nz_sam x 1 array)
[ z_sam, status ] ...
   = mexnc ( 'get_var_double', ncid_sam, varid_sam(idx_3D_z) );

% Time (num_t_sam x 1 array)
[ time_sam, status ] ...
   = mexnc ( 'get_var_double', ncid_sam, varid_sam(idx_3D_time) );

% Variables (nx_sam x ny_sam x nz_sam x num_t_sam array for each)

% Initialize (we can place variables all variables into a
% num_var_sam x nx_sam x ny_sam x nz_sam x num_t_sam array). 
var_sam = zeros( num_var_sam, nx_sam, ny_sam, nz_sam, num_t_sam );

% Read in variables
for i = 1:1:num_var_sam
   [ var_sam(i,:,:,:,:), status ] ...
      = mexnc ( 'get_var_double', ncid_sam, varid_sam(i) );
end

% Close SAM NetCDF file.
status = mexnc ( 'close', ncid_sam );
