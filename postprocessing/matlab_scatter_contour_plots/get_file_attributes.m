% $Id$
function [ dimid_x,  dimid_y,  dimid_z,  dimid_t, ...
           name_x,   name_y,   name_z,   name_t, ...
           length_x, length_y, length_z, length_t ] ...
  = get_file_attributes( ncid, xdim_name, ydim_name, zdim_name, tdim_name )

% Paths to MEXNC and NetCDF Toolbox utilities.
addpath ( '/usr/share/mexcdf/mexnc' )
addpath ( '/usr/local/netcdf_toolbox/netcdf', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/nctype', '-end' )
addpath ( '/usr/local/netcdf_toolbox/netcdf/ncutility', '-end' )

% Get the Dimension IDs from the NetCDF file.
[dimid_x, status] = mexnc ( 'inq_dimid', ncid, xdim_name );
[dimid_y, status] = mexnc ( 'inq_dimid', ncid, ydim_name );
[dimid_z, status] = mexnc ( 'inq_dimid', ncid, zdim_name );
[dimid_t, status] = mexnc ( 'inq_dimid', ncid, tdim_name );
% Find the length of each of the dimension entries.
[name_x, length_x, status] = mexnc ( 'inq_dim', ncid, dimid_x );
[name_y, length_y, status] = mexnc ( 'inq_dim', ncid, dimid_y );
[name_z, length_z, status] = mexnc ( 'inq_dim', ncid, dimid_z );
[name_t, length_t, status] = mexnc ( 'inq_dim', ncid, dimid_t );

